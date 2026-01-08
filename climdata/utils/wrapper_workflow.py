import os
import json
import pandas as pd
import xarray as xr
from pathlib import Path
import warnings
warnings.filterwarnings("default")
from dataclasses import dataclass
from typing import Optional, List, Dict
from functools import wraps

import xclim
import climdata
from climdata.utils.config import _ensure_local_conf
from climdata.utils.utils_download import get_output_filename
from climdata.extremes.indices import extreme_index
from climdata.impute.impute_xarray import Imputer
    
from hydra import initialize, compose
from hydra.core.global_hydra import GlobalHydra
from omegaconf import DictConfig
from shapely.geometry import shape, Polygon, Point
import logging

logger = logging.getLogger(__name__)

# ----------------------------
# Dataclass for workflow result
# ----------------------------
@dataclass
class WorkflowResult:
    cfg: DictConfig
    dataset: Optional[xr.Dataset] = None
    dataframe: Optional[pd.DataFrame] = None
    filename: Optional[str] = None
    index_ds: Optional[xr.Dataset] = None
    index_filename: Optional[str] = None
    impute_ds: Optional[xr.Dataset] = None
    impute_filename: Optional[str] = None

    def keys(self):
        return [k for k, v in self.__dict__.items() if v is not None]

def update_ds(attr_name=None):
    def decorator(func):
        @wraps(func)
        def wrapper(self, *args, **kwargs):
            ds = func(self, *args, **kwargs)
            if ds is not None:
                self.current_ds = ds
                if attr_name:
                    setattr(self, attr_name, ds)
                # log update
                log = getattr(self, "logger", logger)
                log.debug(f"Dataset updated by {func.__name__}; attr_name={attr_name}")
                # Ensure filenames are generated after the dataset update so
                # that filename templates use the newly produced dataset (e.g. index datasets).
                try:
                    if hasattr(self, "_gen_fn_cfg") and callable(getattr(self, "_gen_fn_cfg")):
                        self._gen_fn_cfg()
                except Exception:
                    log.exception("Generating filenames after %s failed", func.__name__)
            return ds
        return wrapper
    return decorator
def update_df(attr_name=None):
    """Decorator to update ``self.current_df`` with the result of a method.

    Args:
        attr_name (str, optional): If provided, also store the result on ``self`` under this attribute name.
    """
    def decorator(func):
        @wraps(func)
        def wrapper(self, *args, **kwargs):
            df = func(self, *args, **kwargs)
            if df is not None:
                self.current_df = df
                self.df = df
                if attr_name:
                    setattr(self, attr_name, df)
                log = getattr(self, "logger", logger)
                log.debug(f"DataFrame updated by {func.__name__}; attr_name={attr_name}")
            return df
        return wrapper
    return decorator
# ----------------------------
# ClimateExtractor Class
# ----------------------------
class ClimateExtractor:
    """Climate data extraction and extreme-index workflow manager.

    Provides a high-level API for:
        - loading/configuring dataset providers via Hydra config
        - uploading NetCDF/CSV content into xarray Datasets
        - extracting data from supported providers (CMIP, DWD, MSWX, HYRAS, POWER)
        - computing extreme indices using configured xclim indices
        - converting datasets to long-form DataFrames and saving results

    Attributes:
        cfg (DictConfig): Hydra configuration object describing dataset, region/time/variables, outputs.
        current_ds (xr.Dataset): The most recently loaded or extracted dataset.
        current_df (pd.DataFrame): The most recently produced long-form DataFrame.
        filename_csv/filename_nc/filename_zarr (str): Generated output filename templates/paths.

    Example:
        extractor = ClimateExtractor(overrides=['dataset=cmip', 'region=europe'])
        extractor.extract()
        idx_ds = extractor.calc_index()
        df = extractor.to_dataframe(idx_ds)
        extractor.to_csv(df)
    """

    def __init__(self, cfg_name="config", conf_path=None, overrides: Optional[List[str]] = None):
        """Initialize the workflow manager and load configuration.

        Args:
            cfg_name (str): Name of the Hydra configuration (default: "config").
            conf_path (str, optional): Optional config path override.
            overrides (list[str], optional): Hydra overrides to apply to the configuration.
        """
        self.cfg_name = cfg_name
        self.conf_path = conf_path
        self.cfg: Optional[DictConfig] = None

        # Stage datasets
        self.ds = None
        self.current_ds = None
        self.index_ds = None
        self.impute_ds = None
        self.bias_corrected_ds = None

        # Stage DataFrames
        self.raw_df = None
        self.current_df = None
        self.index_df = None
        self.impute_df = None
        self.bias_corrected_df = None
        self.df = None  # alias for current_df

        # filenames
        self.filename = None
        self.filetype = None

        # Automatically load config on init
        self.load_config(overrides)
        self.cfg = self.preprocess_aoi(self.cfg)

        # instance logger for this extractor
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
    def _gen_fn(self, ds: xr.Dataset):
        """Create filenames (csv, nc, zarr) using config templates and dataset metadata.

        The method automatically handles coordinate aliases such as lat/latitude, lon/longitude and time/date.
        """

        # ------------------------
        # Helper: find coord alias
        # ------------------------
        def find_coord(ds, names):
            """Return ds coordinate if any alias in names exists."""
            for name in names:
                if name in ds.coords:
                    return ds[name]
            return None

        # ------------------------
        # Extract coordinates
        # ------------------------
        lat = find_coord(ds, ["lat", "latitude"])
        lon = find_coord(ds, ["lon", "longitude"])
        time = find_coord(ds, ["time", "date"])

        # ------------------------
        # Provider
        # ------------------------
        provider = ds.attrs.get("source", "unknown")

        # ------------------------
        # Parameter(s)
        # ------------------------
        vars_list = list(ds.data_vars)
        parameter = vars_list[0] if len(vars_list) == 1 else "_".join(vars_list)

        # ------------------------
        # Latitude range
        # ------------------------
        if lat is None:
            lat_str = "unknown"
            lat_range = "unknown"
        else:
            # Flatten in case of 2D lat/lon grid
            lat_vals = lat.values.reshape(-1)
            lat_min = float(lat_vals.min())
            lat_max = float(lat_vals.max())

            if lat_min == lat_max:
                lat_str = f"{lat_min}"
                lat_range = f"{lat_min}"
            else:
                lat_str = f"{lat_min}_{lat_max}"
                lat_range = f"{lat_min}-{lat_max}"

        # ------------------------
        # Longitude range
        # ------------------------
        if lon is None:
            lon_str = "unknown"
            lon_range = "unknown"
        else:
            lon_vals = lon.values.reshape(-1)
            lon_min = float(lon_vals.min())
            lon_max = float(lon_vals.max())

            if lon_min == lon_max:
                lon_str = f"{lon_min}"
                lon_range = f"{lon_min}"
            else:
                lon_str = f"{lon_min}_{lon_max}"
                lon_range = f"{lon_min}-{lon_max}"

        # ------------------------
        # Time range
        # ------------------------
        if time is None:
            start = end = "unknown"
        else:
            tvals = pd.to_datetime(time.values)
            start = tvals.min().strftime("%Y-%m-%d")
            end = tvals.max().strftime("%Y-%m-%d")

        # ------------------------
        # Build filenames
        # ------------------------
        outdir = Path(self.cfg.output.out_dir)

        def build(fn_template):
            return fn_template.format(
                provider=provider,
                parameter=parameter,
                lat=lat_str,
                lon=lon_str,
                start=start,
                end=end,
                lat_range=lat_range,
                lon_range=lon_range,
            )

        self.filename_csv = str(outdir / build(self.cfg.output.filename_csv))
        self.filename_nc = str(outdir / build(self.cfg.output.filename_nc))
        self.filename_zarr = str(outdir / build(self.cfg.output.filename_zarr))

    
    def _gen_fn_cfg(self):
        """Generate output filenames using configuration and extracted dataset metadata.

        Uses settings from ``self.cfg`` and ``self.current_ds`` (if available) to build filename templates.
        """

        cfg = self.cfg
        out = cfg.output
        provider = cfg.dataset.lower()
        if self.current_ds:
            if len(self.current_ds.data_vars) == 0:
                parameter = "unknown"
            elif len(self.current_ds.data_vars) == 1:
                parameter = next(iter(self.current_ds.data_vars))
            else:
                parameter = "_".join(self.current_ds.data_vars)
        else:
            parameter = "_".join(self.cfg.variables)
        # --------------------------------
        # Determine lat/lon values
        # --------------------------------
        if cfg.lat is not None and cfg.lon is not None:
            lat_range = lon_range = None   # single point
            lat_str = str(cfg.lat)
            lon_str = str(cfg.lon)
        else:
            b = cfg.bounds[cfg.region]
            lat_min, lat_max = b["lat_min"], b["lat_max"]
            lon_min, lon_max = b["lon_min"], b["lon_max"]

            lat_str = f"{lat_min}_{lat_max}"
            lon_str = f"{lon_min}_{lon_max}"
            lat_range = f"{lat_min}-{lat_max}"
            lon_range = f"{lon_min}-{lon_max}"

        # --------------------------------
        # Time range from cfg
        # --------------------------------
        start = pd.to_datetime(cfg.time_range.start_date).strftime("%Y-%m-%d")
        end = pd.to_datetime(cfg.time_range.end_date).strftime("%Y-%m-%d")

        # --------------------------------
        # Format filenames
        # --------------------------------
        def format_template(template):
            return template.format(
                provider=provider,
                parameter=parameter,
                lat=lat_str,
                lon=lon_str,
                start=start,
                end=end,
                lat_range=lat_range or lat_str,
                lon_range=lon_range or lon_str,
            )

        out_dir = Path('./')
        out_dir.mkdir(parents=True, exist_ok=True)

        self.filename_csv = str(out_dir / format_template(out.filename_csv))
        self.filename_nc = str(out_dir / format_template(out.filename_nc))
        self.filename_zarr = str(out_dir / format_template(out.filename_zarr))

    # ----------------------------
    # Hydra config
    # ----------------------------
    def load_config(self, overrides: Optional[List[str]] = None) -> DictConfig:
        """Load and compose the Hydra configuration.

        Args:
            overrides (list[str], optional): Hydra overrides to apply when composing the configuration.

        Returns:
            DictConfig: Composed Hydra configuration object stored on ``self.cfg``.
        """
        overrides = overrides or []
        conf_dir = _ensure_local_conf()
        rel_conf_dir = os.path.relpath(conf_dir, os.path.dirname(__file__))

        if not GlobalHydra.instance().is_initialized():
            hydra_ctx = initialize(config_path=rel_conf_dir, version_base=None)
        else:
            hydra_ctx = None

        if hydra_ctx:
            with hydra_ctx:
                self.cfg = compose(config_name=self.cfg_name, overrides=overrides)
        else:
            self.cfg = compose(config_name=self.cfg_name, overrides=overrides)
        return self.cfg

    # ----------------------------
    # AOI preprocessing
    # ----------------------------
    def preprocess_aoi(self, cfg: DictConfig) -> DictConfig:
        """Process an 'aoi' specification in the configuration.

        Supports GeoJSON strings or dictionaries for FeatureCollection, Feature, or simple geometry objects (Point/Polygon).

        Args:
            cfg (DictConfig): Configuration object with optional ``aoi`` entry.

        Returns:
            DictConfig: The modified configuration. When a Point is provided, ``cfg.lat`` and ``cfg.lon`` are set; when a Polygon is provided, ``cfg.bounds`` is set and ``cfg.region`` is set to "custom".
        """
        if not hasattr(cfg, "aoi") or cfg.aoi is None:
            return cfg

        if isinstance(cfg.aoi, str):
            try:
                cfg.aoi = json.loads(cfg.aoi)
            except json.JSONDecodeError:
                raise ValueError("Invalid AOI JSON string")

        aoi = cfg.aoi

        if aoi.get("type") == "FeatureCollection":
            geom = shape(aoi["features"][0]["geometry"])
        elif aoi.get("type") == "Feature":
            geom = shape(aoi["geometry"])
        elif "type" in aoi:
            geom = shape(aoi)
        else:
            raise ValueError(f"Unsupported AOI format: {aoi}")

        if isinstance(geom, Point):
            cfg.lat = geom.y
            cfg.lon = geom.x
            cfg.bounds = None
        elif isinstance(geom, Polygon):
            minx, miny, maxx, maxy = geom.bounds
            cfg.bounds = {"custom": {"lat_min": miny, "lat_max": maxy,
                                     "lon_min": minx, "lon_max": maxx}}
            cfg.region = "custom"
            cfg.lat = None
            cfg.lon = None
        else:
            raise ValueError(f"Unknown geometry type {geom.geom_type}")

        return cfg

    # ----------------------------
    # Upload NetCDF
    # ----------------------------
    @update_ds(attr_name='ds')
    def upload_netcdf(self, nc_file: str) -> xr.Dataset:
        """Load a NetCDF file into an xarray.Dataset and update file metadata.

        Args:
            nc_file (str): Path to the NetCDF file to open.

        Returns:
            xr.Dataset: The loaded dataset (also sets ``self.current_ds``).
        """
        if not os.path.exists(nc_file):
            raise FileNotFoundError(f"{nc_file} does not exist")

        ds = xr.open_dataset(nc_file)

        # Update cfg variables & varinfo
        if not hasattr(self.cfg, "variables") or not self.cfg.variables:
            self.cfg.variables = list(ds.data_vars)
        if not hasattr(self.cfg, "varinfo") or not self.cfg.varinfo:
            self.cfg.varinfo = {v: {"units": ds[v].attrs.get("units", "unknown")}
                                for v in ds.data_vars}
        self._gen_fn(ds)
        return ds

    # ----------------------------
    # Upload CSV → xarray.Dataset
    # ----------------------------
    @update_ds(attr_name='ds')
    def upload_csv(self, csv_file: str) -> xr.Dataset:
        """Load a long-form CSV into an xarray.Dataset.

        The CSV must contain ``time`` and ``lat``/``latitude``, ``lon``/``longitude``, ``variable``, ``value``. Units may be supplied in a ``units`` column and an optional ``source`` column is recognized.

        Args:
            csv_file (str): Path to the CSV file to load.

        Returns:
            xr.Dataset: The converted dataset (also sets ``self.current_ds``).
        """
        if not os.path.exists(csv_file):
            raise FileNotFoundError(f"{csv_file} does not exist")

        df = pd.read_csv(csv_file, parse_dates=["time"])

        lat_col = next((c for c in ["lat", "latitude"] if c in df.columns), None)
        lon_col = next((c for c in ["lon", "longitude"] if c in df.columns), None)
        if lat_col is None or lon_col is None:
            raise ValueError("CSV must have 'lat'/'latitude' and 'lon'/'longitude' columns")

        id_vars = ["time", lat_col, lon_col]
        df_wide = df.pivot_table(index=id_vars, columns="variable", values="value").reset_index()
        ds = df_wide.set_index(id_vars).to_xarray()

        # Attach units from CSV
        for var in ds.data_vars:
            units_series = df[df["variable"] == var]["units"]
            ds[var].attrs["units"] = units_series.iloc[0] if not units_series.empty else "unknown"

        # Global source attribute
        if "source" in df.columns:
            source_series = df["source"].dropna().unique()
            if len(source_series) > 0:
                ds.attrs["source"] = source_series[0]

        # Update cfg variables & varinfo
        if not hasattr(self.cfg, "variables") or not self.cfg.variables:
            self.cfg.variables = list(ds.data_vars)
        if not hasattr(self.cfg, "varinfo") or not self.cfg.varinfo:
            self.cfg.varinfo = {v: {"units": ds[v].attrs.get("units", "unknown")} for v in ds.data_vars}
        self._gen_fn(ds)
        return ds

    # ----------------------------
    # Extract data from datasets like CMIP, DWD, etc.
    # ----------------------------
    @update_ds(attr_name='ds')
    def extract(self) -> xr.Dataset:
        """Extract data from the configured provider using ``self.cfg``.

        Uses provider-specific classes (e.g., ``CMIP``, ``DWD``, ``MSWX``, ``HYRAS``, ``POWER``)
        to fetch, load and extract datasets. When extraction completes, units are converted to those declared in ``cfg.varinfo``, the dataset is computed, and filenames are generated from the configuration.

        Returns:
            xr.Dataset: The extracted and computed dataset (also sets ``self.current_ds``).
        """
        cfg = self.cfg
        extract_kwargs = {}

        if cfg.lat is not None and cfg.lon is not None:
            extract_kwargs["point"] = (cfg.lon, cfg.lat)
            if cfg.dataset == "dwd":
                extract_kwargs["buffer_km"] = 30
        elif cfg.region is not None:
            extract_kwargs["box"] = cfg.bounds[cfg.region]
        elif cfg.shapefile is not None:
            extract_kwargs["shapefile"] = cfg.shapefile

        ds = None
        dataset_upper = cfg.dataset.upper()

        if dataset_upper == "MSWX":
            ds_vars = []
            for var in cfg.variables:
                mswx = climdata.MSWX(cfg)
                mswx.extract(**extract_kwargs)
                mswx.load(var)
                ds_vars.append(mswx.dataset)
            ds = xr.merge(ds_vars)

        elif dataset_upper == "CMIP":
            cmip = climdata.CMIP(cfg)
            cmip.fetch()
            cmip.load()
            cmip.extract(**extract_kwargs)
            ds = cmip.ds

        elif dataset_upper == "POWER":
            power = climdata.POWER(cfg)
            power.fetch()
            power.load()
            ds = power.ds

        elif dataset_upper == "DWD":
            ds_vars = []
            for var in cfg.variables:
                dwd = climdata.DWD(cfg)
                ds_var = dwd.extract(variable=var, **extract_kwargs)
                ds_vars.append(ds_var)
            ds = xr.merge(ds_vars)

        elif dataset_upper == "HYRAS":
            hyras = climdata.HYRAS(cfg)
            ds_vars = []
            for var in cfg.variables:
                hyras.extract(**extract_kwargs)
                ds_vars.append(hyras.load(var)[[var]])
            ds = xr.merge(ds_vars, compat="override")

        for var in ds.data_vars:
            ds[var] = xclim.core.units.convert_units_to(ds[var], cfg.varinfo[var].units)

        ds = ds.compute()

        return ds
    # ----------------------------
    # Compute extreme index
    # ----------------------------
    @update_ds(attr_name='index_ds')
    def calc_index(self, ds: xr.Dataset = None) -> xr.Dataset:
        """Calculate the configured extreme index using xclim indices.

        Args:
            ds (xr.Dataset, optional): Dataset to operate on. If ``None``, ``self.current_ds`` is used.

        Returns:
            xr.Dataset: The computed index as an xarray Dataset (also sets ``self.index_ds``).
        """
        cfg = self.cfg

        # Use provided ds or fallback
        ds = ds or self.current_ds
        if ds is None:
            raise ValueError("No dataset provided and no current_ds is available.")

        if cfg.index is None:
            self.logger.info("No index selected.")
            return None

        if "time" in ds.coords:
            years = pd.to_datetime(ds.time.values).year
            n_years = len(pd.unique(years))
            if n_years < 30:
                warnings.warn(f"Index {cfg.index} usually requires ≥30 years, got {n_years}", UserWarning)

        indices = extreme_index(cfg, ds)
        index_ds = indices.calculate(cfg.index).compute()
        index_ds = index_ds.to_dataset(name=cfg.index)

        return index_ds
    # ----------------------------
    # Dataset → Long-form DataFrame
    # ----------------------------
    @update_df()
    def to_dataframe(self, ds: xr.Dataset = None) -> pd.DataFrame:
        """Convert a dataset to a long-form pandas DataFrame.

        The output contains columns: time, lat, lon (or latitude/longitude), variable, value, units, source.

        Args:
            ds (xr.Dataset, optional): Dataset to convert. If ``None``, uses ``self.current_ds``.

        Returns:
            pd.DataFrame: Long-form DataFrame (also sets ``self.current_df``).
        """
        ds = ds or self.current_ds
        if ds is None:
            raise ValueError("No dataset provided and no current_ds is available.")
        
        df = ds.to_dataframe().reset_index()
        
        id_vars = [c for c in ("time", "lat", "lon", "latitude", "longitude") if c in df]
        value_vars = [v for v in ds.data_vars if v in df.columns]
        
        if not value_vars:
            raise ValueError("No variables in dataset available to melt into long format")
        
        df_long = df.melt(
            id_vars=id_vars,
            value_vars=value_vars,
            var_name="variable",
            value_name="value"
        )
        
        df_long["units"] = df_long["variable"].apply(
            lambda v: ds[v].attrs.get("units", "unknown")
        )
        df_long["source"] = getattr(self.cfg, "dataset", ds.attrs.get("source", "unknown"))
        self._gen_fn_cfg()
        return df_long

    # ----------------------------
    # Save CSV
    # ----------------------------
    def to_csv(self, df: Optional[pd.DataFrame] = None, filename: Optional[str] = None) -> str:
        """Save a DataFrame to CSV.

        Args:
            df (pd.DataFrame, optional): DataFrame to save. Defaults to ``self.current_df``.
            filename (str, optional): Output filename. Defaults to ``self.filename_csv``.

        Returns:
            str: The path of the written CSV file.
        """
        df = df if df is not None else self.current_df

        filename = filename or getattr(self, "filename_csv", None)
        if filename is None:
            raise ValueError("No filename provided and filename_csv is not set")

        path = Path(filename)
        path.parent.mkdir(parents=True, exist_ok=True)

        df.to_csv(filename, index=False)
        self.filename_csv = str(path)
        self.current_filename = str(path)
        
        # print(f"DataFrame saved to CSV file: {self.current_filename}")
        self.logger.info(f"DataFrame saved to CSV file: {self.current_filename}")

        return filename
    
    def to_nc(self, ds: Optional[xr.Dataset] = None, filename: Optional[str] = None) -> str:
        """Save an xarray Dataset to NetCDF.

        Notes:
            - If ``ds`` is ``None``: save ``current_ds``.
            - If ``filename`` is ``None``: use ``self.filename_nc``.
            - Creates directories if needed and updates ``self.filename_nc`` and ``self.current_filename``.

        Args:
            ds (xr.Dataset, optional): Dataset to save. If ``None``, uses ``self.current_ds``.
            filename (str, optional): Output filename. Defaults to ``self.filename_nc``.

        Returns:
            str: The path of the written NetCDF file.
        """

        # -------------------------------
        # 1. Determine dataset to save
        # -------------------------------
        ds = ds or getattr(self, "current_ds", None)
        if ds is None:
            raise ValueError("No dataset available to save")

        # -------------------------------
        # 2. Determine filename
        # -------------------------------
        filename = filename or getattr(self, "filename_nc", None)
        if filename is None:
            raise ValueError("No filename provided and filename_nc is not set")

        path = Path(filename)
        path.parent.mkdir(parents=True, exist_ok=True)

        # -------------------------------
        # 3. Save to NetCDF
        # -------------------------------
        ds.to_netcdf(path)

        # -------------------------------
        # 4. Track filenames
        # -------------------------------
        self.filename_nc = str(path)
        self.current_filename = str(path)
        
        # print(f"Dataset saved to NetCDF file: {self.current_filename}")
        self.logger.info(f"Dataset saved to NetCDF file: {self.current_filename}")

        return str(path)

    # ----------------------------
    # Unified workflow
    # ----------------------------
    def run_workflow(self, overrides: Optional[List[str]] = None,
                     actions: Optional[List[str]] = None,
                     file: Optional[str] = None) -> WorkflowResult:
        """Execute a sequence of workflow actions.

        Args:
            overrides (list[str], optional): Hydra overrides to apply (not all actions will use these).
            actions (list[str], optional): Ordered list of actions to perform. Supported actions include: 'upload_netcdf', 'upload_csv', 'extract', 'calc_index', 'to_dataframe', 'to_csv', 'to_nc'.
            file (str, optional): File path used for upload actions when required.

        Returns:
            WorkflowResult: Named result container with populated fields for dataframe/dataset/filenames.
        """
        actions = actions or ["extract", "calc_index", "to_csv", "to_nc"]
        result = WorkflowResult(cfg=self.cfg)
        for action in actions:
            self.logger.info("Starting action: %s", action)
            try:
                if action == "upload_netcdf":
                    if file is None:
                        raise ValueError(
                            "Action 'upload_netcdf' requires argument 'netcdf_file', "
                            "but none was provided."
                        )
                        # Validate extension
                    valid_nc_ext = (".nc", ".nc4", ".nc.gz")
                    if not any(str(file).lower().endswith(ext) for ext in valid_nc_ext):
                        raise ValueError(
                            f"Invalid file format for upload_netcdf: '{file}'. "
                            f"Expected one of: {valid_nc_ext}"
                        )
                    self.upload_netcdf(file)
                    result.dataset = self.current_ds

                elif action == "upload_csv":
                    if file is None:
                        raise ValueError(
                            "Action 'upload_csv' requires argument 'csv_file', "
                            "but none was provided."
                        )

                    # Validate CSV extension
                    valid_csv_ext = (".csv", ".csv.gz")
                    if not any(str(file).lower().endswith(ext) for ext in valid_csv_ext):
                        raise ValueError(
                            f"Invalid file format for upload_csv: '{file}'. "
                            f"Expected one of: {valid_csv_ext}"
                        )

                    self.upload_csv(file)
                    result.dataset = self.current_ds

                elif action == "extract":
                    if self.cfg.dataset is None:
                        raise ValueError(
                            "Action 'extract' cannot run because no dataset provider is set "
                            "(cfg.dataset is None)."
                        )
                    self.extract()
                    result.dataset = self.current_ds

                elif action == "calc_index":
                    if self.current_ds is None:
                        raise ValueError(
                            "Action 'calc_index' requires a dataset, but no dataset is available. "
                            "Upload or extract a dataset before computing an index."
                        )
                    self.calc_index()
                    result.index_ds = self.current_ds

                elif action == "to_csv":
                    if self.current_ds is None:
                        raise ValueError(
                            "Action 'to_dataframe' requires a dataset, but no dataset is available. "
                            "Upload or extract a dataset before converting to a DataFrame."
                        )
                    self.to_dataframe()
                    result.dataframe = self.current_df
                    result.filename = self.to_csv()

                elif action == "to_nc":
                    if self.current_ds is None:
                        raise ValueError(
                            "Action 'to_nc' requires a dataset, but no dataset is available. "
                            "Upload or extract a dataset before saving to NetCDF."
                        )
                    result.filename = self.to_nc()

                elif action == "impute":
                    if self.current_ds is None:
                        raise ValueError("Action 'impute' requires a dataset, but no dataset is available.")
                    self.impute()
                    result.dataset = self.current_ds
                    result.impute_ds = getattr(self, "impute_ds", None)

                else:
                    raise ValueError(f"Unknown action '{action}'")
                self.logger.info("Completed action: %s", action)
            except Exception:
                self.logger.exception("Action '%s' failed", action)
                raise

        return result

    # ----------------------------
    # Exploration helpers using cfg.dsinfo
    # ----------------------------
    def get_datasets(self) -> List[str]:
        """Return the list of dataset provider names available in configuration.

        Returns:
            List[str]: Names of available dataset providers from ``cfg.dsinfo``.
        """
        if not self.cfg or not hasattr(self.cfg, "dsinfo"):
            raise ValueError("Configuration or dsinfo not loaded")
        return list(self.cfg.dsinfo.keys())

    def get_variables(self, dataset: Optional[str] = None) -> List[str]:
        """Return the list of variables available for a dataset.

        Args:
            dataset (str, optional): Dataset name to query. Defaults to ``cfg.dataset``.

        Returns:
            List[str]: List of variable names.
        """
        if not self.cfg or not hasattr(self.cfg, "dsinfo"):
            raise ValueError("Configuration or dsinfo not loaded")

        dataset_name = dataset or getattr(self.cfg, "dataset", None)
        if dataset_name is None:
            raise ValueError("Dataset not specified and cfg.dataset is None")

        dsinfo = self.cfg.dsinfo.get(dataset_name)
        if not dsinfo or "variables" not in dsinfo:
            raise ValueError(f"No variable info available for dataset '{dataset_name}'")

        return list(dsinfo["variables"].keys())

    def get_varinfo(self, var: str) -> dict:
        """Get metadata for a variable from varinfo.

        Args:
            var (str): Name of the variable, e.g., 'tas', 'tasmax', 'pr'.

        Returns:
            dict: Metadata dictionary containing cf_name, long_name, units, etc.

        Raises:
            ValueError: If varinfo is not loaded or variable not found.
        """
        if not self.cfg or not hasattr(self.cfg, "varinfo") or not self.cfg.varinfo:
            raise ValueError("Configuration or varinfo not loaded")

        if var not in self.cfg.varinfo:
            raise ValueError(f"Variable '{var}' not found in varinfo")

        return self.cfg.varinfo[var]

    
    def get_actions(self) -> dict:
        """Return a dictionary of workflow actions with their outputs and descriptions.

        Supports ``actionsinfo`` in mapping style or list style and returns a consistent mapping of action name to description/output.

        Returns:
            dict: Mapping action name -> {'output': ..., 'description': ...}
        """
        if not self.cfg or not hasattr(self.cfg, "actionsinfo"):
            raise ValueError("Configuration or actionsinfo not loaded")

        actions_map = getattr(self.cfg, "actionsinfo")
        
        # If 'actions' key exists, fallback to list style
        if "actions" in actions_map:
            actions_map = {a["name"]: {"output": a["output"], "description": a["description"]}
                        for a in actions_map["actions"]}
        
        return actions_map
    def get_indices(self, variables: List[str], require_all: bool = True) -> Dict[str, dict]:
        """Fetch climate extreme indices from ``cfg.extinfo`` that involve the given variables.

        Args:
            variables (list[str]): Variables to filter indices by (if ``None``, uses ``cfg.variables``).
            require_all (bool): If True, return indices that require all provided variables; otherwise return indices if any variable matches.

        Returns:
            dict: Mapping index_name -> index_definition.
        """
        cfg = self.cfg
        variables = variables or cfg.variables 
        if not hasattr(cfg, "extinfo") or not cfg.extinfo:
            raise ValueError("cfg.extinfo is not defined or empty")

        indices_def = cfg.extinfo.get("indices", {})
        if not indices_def:
            return {}

        matched_indices = {}
        for idx_name, idx_info in indices_def.items():
            idx_vars = idx_info.get("variables", [])
            if require_all:
                if all(var in variables for var in idx_vars):
                    matched_indices[idx_name] = idx_info
            else:
                if any(var in variables for var in idx_vars):
                    matched_indices[idx_name] = idx_info

        return matched_indices

    # ----------------------------
    # Imputation
    # ----------------------------
    @update_ds(attr_name='impute_ds')
    def impute(self, ds: xr.Dataset = None) -> xr.Dataset:
        """Impute missing values using the configured imputation method.

        Args:
            ds (xr.Dataset, optional): Dataset to impute. If None, uses
                ``self.current_ds``.

        Returns:
            xr.Dataset | None: The imputed dataset (also sets
                ``self.current_ds`` and ``self.impute_ds``). Returns ``None``
                if no imputation method is configured.

        Raises:
            ValueError: If ``ds`` is ``None`` and ``self.current_ds`` is not set.
        """
        cfg = self.cfg
        impute_cfg = cfg.imputeinfo
        ds = ds or self.current_ds
        if ds is None:
            raise ValueError("No dataset provided and no current_ds is available.")

        if cfg.impute is None:
            self.logger.warning("No imputation method selected.")
            return None
        # select variables (optional)
        # variables = cfg.get("variables", None)
        # if variables:
        #     missing = [v for v in variables if v not in self.current_ds.data_vars]
        #     if missing:
        #         raise ValueError(f"Variables not present in dataset: {missing}")
        #     ds_in = self.current_ds[variables]
        # else:
        #     ds_in = self.current_ds

        method = cfg.impute
        normalize = impute_cfg[method].get("normalize", True)
        time_dim = cfg.dsinfo[cfg.dataset].get("time_dim", "time")
        lat_dim = cfg.dsinfo[cfg.dataset].get("lat_dim", "lat")
        lon_dim = cfg.dsinfo[cfg.dataset].get("lon_dim", "lon")
        # epochs = impute_cfg[method].get("epochs", 300)

        # run imputer (Imputer expects dims (time, lat, lon))
        imputer = Imputer(
            ds,
            time_dim=time_dim,
            lat_dim=lat_dim,
            lon_dim=lon_dim,
            method=method,
            normalize=normalize,
        )
        recovered = imputer.impute()

        # merge imputed variables back into original dataset if we operated on a subset

        ds_out = recovered

        # Return dataset (decorator will set current_ds and impute_ds and generate filenames)
        return ds_out

    def get_impute_methods(self) -> Dict[str, dict]:
        """Return mapping of available imputation methods from config.

        Returns:
            Dict[str, dict]: Mapping of method name -> config (empty dict if none configured).
        """
        if not hasattr(self.cfg, "imputeinfo") or not self.cfg.imputeinfo:
            return {}
        return dict(self.cfg.imputeinfo)
    
    def configure_logging(self, level=logging.INFO, handler: logging.Handler = None):
        """Configure logging for this extractor instance.

        Args:
            level (int, optional): Logging level (default: ``logging.INFO``).
            handler (logging.Handler, optional): Handler to add; if ``None``, a default StreamHandler is created.
        """
        if handler is None:
            handler = logging.StreamHandler()
            handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
        # Avoid adding duplicate handlers
        if not any(isinstance(h, handler.__class__) for h in self.logger.handlers):
            self.logger.addHandler(handler)
        self.logger.setLevel(level)
        # also set module logger level
        logger.setLevel(level)