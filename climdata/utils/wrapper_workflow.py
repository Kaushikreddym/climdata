"""The workflow layer: one class that drives a whole climdata run.

:class:`ClimateExtractor` (exported as ``climdata.ClimData``) composes the
provider classes, the extreme-index engine, the imputer and the regridder behind
a single Hydra configuration, so a run is described by overrides rather than by
code. Individual steps can be called directly, or chained through
:meth:`ClimateExtractor.run_workflow`.

State is threaded between steps by two decorators, :func:`update_ds` and
:func:`update_df`, which capture each method's return value into
``self.current_ds`` / ``self.current_df``. That is what lets ``calc_index()``
follow ``extract()`` with no argument passing.

Example:
    >>> import climdata as cd
    >>> ex = cd.ClimData(overrides=["dataset=mswx", "lat=52.5", "lon=13.4"])   # doctest: +SKIP
    >>> result = ex.run_workflow(actions=["extract", "calc_index", "to_nc"])   # doctest: +SKIP
"""

import os
import json
import pandas as pd
import numpy as np
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
from climdata.extremes.indices import extreme_index
from climdata.impute.impute_xarray import Imputer
    
from hydra import initialize_config_dir, compose
from hydra.core.global_hydra import GlobalHydra
from omegaconf import DictConfig
from shapely.geometry import shape, Polygon, Point
import logging

logger = logging.getLogger(__name__)

# ----------------------------
# CF to DWD Variable Mapping
# ----------------------------
# Maps Climate and Forecast (CF) standard names to DWD (Deutscher Wetterdienst) names
CF_TO_DWD_NAMES = {
    # Temperature variables
    'tas': 'TempMean',           # Daily mean temperature
    'tasmin': 'TempMin',         # Daily minimum temperature
    'tasmax': 'TempMax',         # Daily maximum temperature
    # Precipitation
    'pr': 'Precipitation',       # Daily precipitation
    # Solar radiation
    'rsds': 'Radiation',         # Shortwave downwelling radiation
    'rlds': 'LongwaveRadiation', # Longwave downwelling radiation
    # Wind
    'sfcWind': 'Windspeed',      # Surface wind speed
    # Humidity
    'hurs': 'RelHum',            # Relative humidity
    'huss': 'SpecHum',           # Specific humidity
    # Other
    'ps': 'SurfPressure',        # Surface pressure
}

#: Actions ``ClimateExtractor.run_workflow()`` dispatches. Kept in step with
#: ``conf/mappings/actions.yaml``, which is what ``get_actions()`` reports.
_WORKFLOW_ACTIONS = (
    "upload_netcdf", "upload_csv", "extract", "calc_index",
    "impute", "reproject", "to_csv", "to_nc", "to_fair",
)

#: Providers ``ClimateExtractor.extract()`` knows how to drive, for error messages.
_EXTRACTABLE = (
    "MSWX", "CMIP", "POWER", "DWD", "HYRAS", "HOSTRADA",
    "W5E5", "CMIP_W5E5", "NEXGDDP", "AGRI_ISIMIP",
)

# ----------------------------
# Dataclass for workflow result
# ----------------------------
@dataclass
class WorkflowResult:
    """What :meth:`ClimateExtractor.run_workflow` produces.

    Every field but ``cfg`` is optional and stays ``None`` unless an action in
    the sequence populated it, so which attributes are set records which actions
    ran. :meth:`keys` names them.

    Attributes:
        cfg (DictConfig): The configuration the workflow ran with.
        dataset (xr.Dataset | None): Extracted or uploaded data.
        dataframe (pd.DataFrame | None): Long-form conversion of the dataset.
        filename (str | None): Path written by the last file-writing action.
        index_ds (xr.Dataset | None): Computed extreme index.
        index_filename (str | None): Path of the written index.
        impute_ds (xr.Dataset | None): Gap-filled dataset.
        impute_filename (str | None): Path of the written imputed dataset.
        reprojected_ds (xr.Dataset | None): Regridded dataset.
        reprojected_filename (str | None): Path of the written regridded dataset.
    """

    cfg: DictConfig
    dataset: Optional[xr.Dataset] = None
    dataframe: Optional[pd.DataFrame] = None
    filename: Optional[str] = None
    index_ds: Optional[xr.Dataset] = None
    index_filename: Optional[str] = None
    impute_ds: Optional[xr.Dataset] = None
    impute_filename: Optional[str] = None
    reprojected_ds: Optional[xr.Dataset] = None
    reprojected_filename: Optional[str] = None

    def keys(self):
        """Name the fields that were actually populated.

        Returns:
            list[str]: Attribute names whose value is not ``None`` — in effect,
            a record of which workflow actions produced output.
        """
        return [k for k, v in self.__dict__.items() if v is not None]

def update_ds(attr_name=None):
    """Capture a method's returned Dataset as the extractor's current state.

    Stores the result on ``self.current_ds`` (and optionally a second attribute),
    then regenerates the output filenames so they describe what was just
    produced rather than the previous stage — an index dataset gets an
    index-shaped filename, not the raw extraction's.

    A method returning ``None`` — as :meth:`ClimateExtractor.calc_index` does
    when no index is configured — leaves the current state untouched. Filename
    regeneration is logged and swallowed on failure, so a template that does not
    fit an unusual dataset costs the filenames, not the data.

    Args:
        attr_name (str, optional): Second attribute to store the result under,
            e.g. ``"index_ds"``, giving each stage a stable handle.

    Returns:
        Callable: A decorator for methods that return an ``xr.Dataset``.
    """
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
    """Capture a method's returned DataFrame as the extractor's current state.

    The frame counterpart of :func:`update_ds`; it stores the result on both
    ``self.current_df`` and ``self.df`` but does not touch filenames.

    Args:
        attr_name (str, optional): Second attribute to store the result under,
            e.g. ``"index_df"``.

    Returns:
        Callable: A decorator for methods that return a ``pd.DataFrame``.
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
        self.reprojected_ds = None
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

        # Optional shared Dask distributed client (created lazily by extract())
        self._dask_client = None

        # Automatically load config on init
        self.load_config(overrides)
        self.cfg = self.preprocess_aoi(self.cfg)

        # instance logger for this extractor
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
    def _gen_fn(self, *, ds: xr.Dataset = None, df: pd.DataFrame = None):
        """Build output filenames from *data* metadata rather than from config.

        The counterpart to :meth:`_gen_fn_cfg`, used after uploading a file: the
        provider, variables and extent are read from the data itself, because an
        uploaded dataset need not match what the configuration describes.

        Args:
            ds (xr.Dataset, optional): Dataset to describe. Keyword-only.
            df (pd.DataFrame, optional): Long-form frame with ``lat``/``lon``,
                ``time`` or ``date``, ``variable``, ``value`` and optionally
                ``source`` columns. Keyword-only.

        Returns:
            ClimateExtractor: ``self``, with ``filename_csv``, ``filename_nc``
            and ``filename_zarr`` set.

        Raises:
            ValueError: If both or neither of ``ds`` and ``df`` are given.
        """

        # ------------------------
        # Validate inputs
        # ------------------------
        if (ds is None) == (df is None):
            raise ValueError("Provide exactly one of `ds` or `df` as a keyword argument.")

        # ------------------------
        # Helper: find coord alias (xarray)
        # ------------------------
        def find_coord(ds, names):
            for name in names:
                if name in ds.coords:
                    return ds[name]
            return None

        # ------------------------
        # Case 1: xarray.Dataset
        # ------------------------
        if ds is not None:
            lat = find_coord(ds, ["lat", "latitude"])
            lon = find_coord(ds, ["lon", "longitude"])
            time = find_coord(ds, ["time", "date"])

            provider = ds.attrs.get("source", "unknown")
            vars_list = list(ds.data_vars)
            parameter = vars_list[0] if len(vars_list) == 1 else "_".join(vars_list)

            # Latitude range
            if lat is not None:
                lat_vals = lat.values.reshape(-1)
                lat_min, lat_max = float(lat_vals.min()), float(lat_vals.max())
            else:
                lat_min = lat_max = None

            # Longitude range
            if lon is not None:
                lon_vals = lon.values.reshape(-1)
                lon_min, lon_max = float(lon_vals.min()), float(lon_vals.max())
            else:
                lon_min = lon_max = None

            # Time range
            if time is not None:
                tvals = pd.to_datetime(time.values)
                start, end = tvals.min().strftime("%Y-%m-%d"), tvals.max().strftime("%Y-%m-%d")
            else:
                start = end = "unknown"

        # ------------------------
        # Case 2: pandas.DataFrame (long form)
        # ------------------------
        else:
            cols = df.columns.astype(str)

            # Identify coordinate columns
            lat_cols = [c for c in cols if "lat" in c.lower()]
            lon_cols = [c for c in cols if "lon" in c.lower()]
            time_cols = [c for c in cols if "time" in c.lower() or "date" in c.lower()]

            # Provider from 'source' column
            if "source" in df.columns:
                unique_sources = df["source"].dropna().unique()
                provider = unique_sources[0] if len(unique_sources) == 1 else "_".join(map(str, unique_sources))
            else:
                provider = "unknown"

            # Unique parameters from 'variable' column
            if "variable" in df.columns:
                unique_parameters = sorted(df["variable"].dropna().unique())
                parameter = unique_parameters[0] if len(unique_parameters) == 1 else "_".join(unique_parameters)
            else:
                parameter = "unknown"

            # Latitude range
            if lat_cols:
                lat_vals = pd.to_numeric(df[lat_cols[0]], errors="coerce")
                lat_min, lat_max = float(lat_vals.min()), float(lat_vals.max())
            else:
                lat_min = lat_max = None

            # Longitude range
            if lon_cols:
                lon_vals = pd.to_numeric(df[lon_cols[0]], errors="coerce")
                lon_min, lon_max = float(lon_vals.min()), float(lon_vals.max())
            else:
                lon_min = lon_max = None

            # Time range
            if time_cols:
                tvals = pd.to_datetime(df[time_cols[0]], errors="coerce")
                start = tvals.min().strftime("%Y-%m-%d")
                end = tvals.max().strftime("%Y-%m-%d")
            else:
                start = end = "unknown"

        # ------------------------
        # Format lat/lon strings
        # ------------------------
        if lat_min is None:
            lat_str = lat_range = "unknown"
        elif lat_min == lat_max:
            lat_str = lat_range = f"{lat_min}"
        else:
            lat_str = f"{lat_min}_{lat_max}"
            lat_range = f"{lat_min}-{lat_max}"

        if lon_min is None:
            lon_str = lon_range = "unknown"
        elif lon_min == lon_max:
            lon_str = lon_range = f"{lon_min}"
        else:
            lon_str = f"{lon_min}_{lon_max}"
            lon_range = f"{lon_min}-{lon_max}"

        # ------------------------
        # Build filenames
        # ------------------------
        outdir = Path(self.cfg.output.out_dir)
        outdir.mkdir(parents=True, exist_ok=True)
        def build(fn_template):
            # The same placeholder set as _gen_fn_cfg, so one template works for
            # both the configured and the data-derived path. `lat_or_lat_range`
            # is what conf/config.yaml actually uses.
            return fn_template.format(
                provider=provider,
                parameter=parameter,
                lat=lat_str,
                lon=lon_str,
                start=start,
                end=end,
                lat_range=lat_range,
                lon_range=lon_range,
                lat_or_lat_range=lat_range,
                lon_or_lon_range=lon_range,
            )

        self.filename_csv = str(outdir / build(self.cfg.output.filename_csv))
        self.filename_nc = str(outdir / build(self.cfg.output.filename_nc))
        self.filename_zarr = str(outdir / build(self.cfg.output.filename_zarr))
        return self
    def _gen_fn_cfg(self):
        """Build output filenames from the configuration and the current dataset.

        The usual path, called automatically by :func:`update_ds` after every
        stage. The provider name, region and period come from ``cfg``; the
        variable part comes from ``current_ds`` when one exists, so the filename
        follows the data through the workflow — after ``calc_index()`` it names
        the index rather than the input variables.

        Returns:
            None: ``filename_csv``, ``filename_nc`` and ``filename_zarr`` are set
            on the instance, and ``cfg.output.out_dir`` is created.
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
        elif cfg.region is not None:
            b = cfg.bounds[cfg.region]
            lat_min, lat_max = b["lat_min"], b["lat_max"]
            lon_min, lon_max = b["lon_min"], b["lon_max"]

            lat_str = f"{lat_min}_{lat_max}"
            lon_str = f"{lon_min}_{lon_max}"
            lat_range = f"{lat_min}-{lat_max}"
            lon_range = f"{lon_min}-{lon_max}"
        else:
            # Fallback for undefined regions
            lat_str = "unknown"
            lon_str = "unknown"
            lat_range = None
            lon_range = None

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
                lat_or_lat_range=lat_range or lat_str,
                lon_or_lon_range=lon_range or lon_str,
                start=start,
                end=end,
            )

        out_dir = Path(self.cfg.output.out_dir)
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
        conf_dir_abs = os.path.abspath(_ensure_local_conf())

        GlobalHydra.instance().clear()
        with initialize_config_dir(config_dir=conf_dir_abs, version_base=None):
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
        self._gen_fn(ds=ds)
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

        # Attach units from CSV. A blank cell reads back as NaN, which is not a
        # units string — it must become "unknown" like an absent column, or
        # xclim fails downstream on a float where it expects a unit.
        for var in ds.data_vars:
            units_series = (
                df[df["variable"] == var]["units"].dropna()
                if "units" in df.columns
                else pd.Series(dtype=object)
            )
            ds[var].attrs["units"] = (
                str(units_series.iloc[0]) if not units_series.empty else "unknown"
            )

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
        self._gen_fn(ds=ds)
        return ds

    # ----------------------------
    # Extract data from datasets like CMIP, DWD, etc.
    # ----------------------------
    def _ensure_dask_client(self):
        """Optionally start a shared Dask distributed cluster for parallel opens/compute.

        Controlled by ``cfg.load.dask``. Opt-in (``enabled: false`` by default) so
        other entry points are unaffected. Idempotent: reuses a client this
        extractor already created, or one the user started elsewhere in the
        process, instead of spinning up a second cluster.

        Once a distributed client exists in the process, ``open_mfdataset(parallel=True)``
        routes its per-file opens through it automatically — no other code changes.

        Returns:
            distributed.Client | None: The active client, or None when disabled/unavailable.
        """
        from omegaconf import OmegaConf

        if not OmegaConf.select(self.cfg, "load.dask.enabled", default=False):
            return None

        # Reuse a client we already created for this extractor.
        if self._dask_client is not None:
            return self._dask_client

        try:
            from dask.distributed import Client, LocalCluster, get_client
        except ImportError:
            self.logger.warning(
                "cfg.load.dask.enabled=true but dask.distributed is not installed; "
                "falling back to the default scheduler."
            )
            return None

        # Reuse an externally-created client if the user already started one.
        try:
            self._dask_client = get_client()
            self.logger.info("Reusing existing Dask client (dashboard: %s)", self._dask_client.dashboard_link)
            return self._dask_client
        except ValueError:
            pass  # no active client -> create our own below

        cluster = LocalCluster(
            n_workers=int(OmegaConf.select(self.cfg, "load.dask.n_workers", default=4)),
            threads_per_worker=int(OmegaConf.select(self.cfg, "load.dask.threads_per_worker", default=2)),
            memory_limit=OmegaConf.select(self.cfg, "load.dask.memory_limit", default="auto"),
        )
        self._dask_client = Client(cluster)
        self.logger.info("Started Dask cluster (dashboard: %s)", self._dask_client.dashboard_link)
        return self._dask_client

    def close_dask(self):
        """Shut down the Dask cluster this extractor started, if any.

        Safe to call unconditionally: it is a no-op when no cluster was started.
        A cluster the extractor merely reused, rather than created, is left
        running for whoever owns it.

        Returns:
            None
        """
        client = self._dask_client
        if client is None:
            return
        cluster = getattr(client, "cluster", None)
        try:
            client.close()
            if cluster is not None:
                cluster.close()
        finally:
            self._dask_client = None

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

        # Optionally start a shared Dask cluster so per-file opens (and later
        # compute) run in parallel across worker processes.
        self._ensure_dask_client()

        if cfg.lat is not None and cfg.lon is not None:
            extract_kwargs["point"] = (cfg.lon, cfg.lat)
            if cfg.dataset.upper() == "DWD":
                extract_kwargs["buffer_km"] = 30
        elif cfg.box.lat_min is not None:
            extract_kwargs["box"] = {
                "lat_min": cfg.box.lat_min,
                "lat_max": cfg.box.lat_max,
                "lon_min": cfg.box.lon_min,
                "lon_max": cfg.box.lon_max,
            }
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
            self.dataset_class = mswx
        elif dataset_upper == "CMIP":
            cmip = climdata.CMIP(cfg)
            cmip.fetch()
            cmip.load()
            cmip.extract(**extract_kwargs)
            ds = cmip.ds
            self.dataset_class = cmip
        elif dataset_upper == "POWER":
            power = climdata.POWER(cfg)
            power.fetch()
            power.load()
            ds = power.ds
            self.dataset_class = power
        elif dataset_upper == "DWD":
            ds_vars = []
            for var in cfg.variables:
                dwd = climdata.DWD(cfg)
                ds_var = dwd.extract(variable=var, **extract_kwargs)
                ds_vars.append(ds_var)
            ds = xr.merge(ds_vars)
            self.dataset_class = dwd
        elif dataset_upper == "HYRAS":
            hyras = climdata.HYRAS(cfg)
            ds_vars = []
            for var in cfg.variables:
                hyras.extract(**extract_kwargs)
                ds_vars.append(hyras.load(var, chunking={'time':"auto"})[[var]])
            ds = xr.merge(ds_vars, compat="override")
            self.dataset_class = hyras
        elif dataset_upper == "HOSTRADA":
            hostrada = climdata.HOSTRADA(cfg)
            hostrada.extract(**extract_kwargs)
            ds_vars = []
            for var in cfg.variables:
                ds_vars.append(hostrada.load(var, chunking={'time': 'auto'})[[var]])
            ds = xr.merge(ds_vars, compat="override")
            self.dataset_class = hostrada
        elif dataset_upper == "W5E5":
            w5e5 = climdata.W5E5(cfg)
            w5e5.fetch()  # Download from ISIMIP
            w5e5.load()   # Load into xarray
            w5e5.extract(**extract_kwargs)
            ds = w5e5.ds
            self.dataset_class = w5e5
        elif dataset_upper == "CMIP_W5E5":
            cmip_w5e5 = climdata.CMIPW5E5(cfg)
            cmip_w5e5.fetch()  # Download CMIP6 data from ISIMIP
            cmip_w5e5.load()   # Load into xarray
            cmip_w5e5.extract(**extract_kwargs)
            ds = cmip_w5e5.ds
            self.dataset_class = cmip_w5e5
        elif dataset_upper == "NEXGDDP":
            nexgddp = climdata.NEXGDDP(cfg)
            nexgddp.fetch()
            nexgddp.extract(**extract_kwargs)  # set params BEFORE load so preprocess clips per-file
            nexgddp.load()                     # spatial clip applied inside open_mfdataset preprocess
            ds = nexgddp.ds
            self.dataset_class = nexgddp
        elif dataset_upper == "AGRI_ISIMIP":
            agri = climdata.AGRI_ISIMIP(cfg)
            agri.fetch()  # Download agricultural data from ISIMIP
            ds_vars = []
            for var in cfg.variables:
                ds_vars.append(agri.load(var, chunking={'time': 'auto'})[[var]])
            ds = xr.merge(ds_vars, compat="override")
            agri.extract(**extract_kwargs)
            self.dataset_class = agri
        else:
            # Without this the unhandled provider leaves ds as None and fails
            # below with a bare `TypeError: 'NoneType' is not subscriptable`,
            # which says nothing about the actual mistake.
            raise ValueError(
                f"Dataset provider '{cfg.dataset}' has no extraction path in "
                f"ClimateExtractor.extract(). Supported: {', '.join(_EXTRACTABLE)}.\n"
                f"   ERA5 is available as climdata.ERA5, but it mirrors whole "
                f"months to Zarr rather than extracting, so it is driven directly "
                f"rather than through this workflow."
            )
        for var in cfg.variables:
            ds[var] = xclim.core.units.convert_units_to(ds[var], cfg.varinfo[var].units)

        # ds = ds.compute()

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
            # Handle both numpy datetime64 and cftime datetime objects
            time_values = ds.time.values
            try:
                # Try pandas conversion first (works for datetime64)
                years = pd.to_datetime(time_values).year
            except (TypeError, ValueError):
                # Fall back to cftime handling
                years = np.array([t.year for t in time_values])
            
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
    def plot(self, variable=None, ds: xr.Dataset = None, **kwargs):
        """Plot a variable on a Cartopy map with coastlines and country borders.

        Thin wrapper around :func:`climdata.viz.plot_map` operating on
        ``self.current_ds`` by default.

        Args:
            variable (str, optional): Variable to plot. Defaults to the sole
                data variable if the dataset has exactly one.
            ds (xr.Dataset, optional): Dataset to plot. If ``None``, uses ``self.current_ds``.
            **kwargs: Forwarded to :func:`climdata.viz.plot_map`
                (e.g. ``time``, ``reduce``, ``cmap``, ``ax``).

        Returns:
            The Cartopy ``GeoAxes`` the field was drawn on.
        """
        from climdata.viz import plot_map

        ds = ds if ds is not None else self.current_ds
        if ds is None:
            raise ValueError("No dataset provided and no current_ds is available.")
        return plot_map(ds, variable=variable, **kwargs)

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
        if getattr(self.cfg, "dataset") == 'cmip':
            df_long["source_id"] = getattr(self.cfg, "source_id")
        df_long["source"] = getattr(self.cfg, "dataset", ds.attrs.get("source", "unknown"))
        df_long = df_long.drop_duplicates()
        self._gen_fn_cfg()
        return df_long

    # ----------------------------
    # Save CSV
    # ----------------------------
    def to_csv(self, df: Optional[pd.DataFrame] = None, filename: Optional[str] = None, format: str = "default") -> str:
        """Save a DataFrame to CSV with optional format specification.

        Args:
            df (pd.DataFrame, optional): DataFrame to save. Defaults to ``self.current_df``.
            filename (str, optional): Output filename/directory. Defaults to ``self.filename_csv``.
            format (str): Output format. Options:
                - 'default': Long-form (single file)
                - 'simplace': SIMPLACE format (splits by location, tab-separated)
                - 'monica': MONICA format (splits by location, tab-separated)
                Defaults to 'default'.

        Returns:
            str: The path of the written CSV file(s). For SIMPLACE/MONICA, returns base directory.
            
        Raises:
            ValueError: If format is not supported or required data is missing.
        """
        df = df if df is not None else self.current_df
        format_lower = format.lower()

        # Apply format conversion if requested
        if format_lower in ("simplace", "monica"):
            return self._convert_to_gridded_format(df, filename, format=format_lower)
        elif format_lower == "default":
            filename = filename or getattr(self, "filename_csv", None)
            if filename is None:
                raise ValueError("No filename provided and filename_csv is not set")

            path = Path(filename)
            path.parent.mkdir(parents=True, exist_ok=True)

            df.to_csv(filename, index=False, sep='\t')
            self.filename_csv = str(path)
            self.current_filename = str(path)
            
            self.logger.info(f"DataFrame saved to CSV file ({format} format): {self.current_filename}")

            return filename
        else:
            raise ValueError(f"Unsupported format: {format}. Supported formats: 'default', 'simplace', 'monica'")

    def _convert_to_gridded_format(self, df: pd.DataFrame, filename: Optional[str] = None, format: str = "simplace") -> str:
        """Convert long-form DataFrame to gridded SIMPLACE/MONICA format by row/column indices.

        Creates a folder structure: <base_dir>/col_<col_number>/<format>_variables_row_<r>_col_<c>.csv

        Args:
            df (pd.DataFrame): Long-form DataFrame with columns: date/time, lat/x, lon/y, variable, value
            filename (str, optional): Base output directory. If None, creates 'simplace_output' or 'monica_output'
            format (str): Format type - 'simplace' or 'monica'

        Returns:
            str: Base directory path where all files were created

        Raises:
            ValueError: If DataFrame format is incompatible or required columns missing.
        """
        # Validate input DataFrame structure
        time_col = None
        for col in df.columns:
            if col.lower() in ('date', 'time'):
                time_col = col
                break
        if time_col is None:
            raise ValueError("DataFrame must contain 'date' or 'time' column")

        var_col = None
        for col in df.columns:
            if col.lower() == 'variable':
                var_col = col
                break
        if var_col is None:
            raise ValueError("DataFrame must contain 'variable' column")

        val_col = None
        for col in df.columns:
            if col.lower() in ('value', 'data'):
                val_col = col
                break
        if val_col is None:
            raise ValueError("DataFrame must contain 'value' or 'data' column")

        # Handle both lat/lon and x/y dimensions
        lat_col = next((c for c in df.columns if c.lower() in ('lat', 'latitude', 'y')), None)
        lon_col = next((c for c in df.columns if c.lower() in ('lon', 'longitude', 'x')), None)

        if lat_col is None or lon_col is None:
            raise ValueError("DataFrame must contain 'lat'/'latitude'/'y' and 'lon'/'longitude'/'x' columns")

        # Get dimensions from xarray Dataset if available
        if hasattr(self, 'ds') and self.ds is not None:
            ds = self.ds
            # Determine dimension names (handle lat/lon or y/x)
            lat_dim = next((d for d in ds.dims if d.lower() in ('lat', 'latitude', 'y')), None)
            lon_dim = next((d for d in ds.dims if d.lower() in ('lon', 'longitude', 'x')), None)
            
            if lat_dim and lon_dim:
                n_rows = ds.sizes[lat_dim]
                n_cols = ds.sizes[lon_dim]
                self.logger.info(f"Dataset dimensions: {n_rows} rows (lat/y) × {n_cols} columns (lon/x)")
            else:
                n_rows = df[lat_col].nunique()
                n_cols = df[lon_col].nunique()
                self.logger.info(f"Inferred dimensions: {n_rows} rows × {n_cols} columns")
        else:
            n_rows = df[lat_col].nunique()
            n_cols = df[lon_col].nunique()
            self.logger.info(f"DataFrame dimensions: {n_rows} rows × {n_cols} columns")

        # Create base output directory
        if filename is None:
            filename = f"{format}_{self.cfg.dataset}"
        base_dir = Path(filename)
        base_dir.mkdir(parents=True, exist_ok=True)

        # Get unique lat/lon pairs and assign row/column indices
        unique_locations = df.drop_duplicates(subset=[lat_col, lon_col]).sort_values([lat_col, lon_col]).reset_index(drop=True)
        
        # Create mapping from lat/lon to row/col indices
        lat_sorted = sorted(df[lat_col].unique(), reverse=True)  # rows decrease from north to south
        lon_sorted = sorted(df[lon_col].unique())  # columns increase from west to east
        
        lat_to_row = {lat: i for i, lat in enumerate(lat_sorted)}
        lon_to_col = {lon: j for j, lon in enumerate(lon_sorted)}
        
        n_locations = len(unique_locations)
        self.logger.info(f"Converting to {format.upper()} format with {n_locations} unique locations")

        created_files = []

        # Process each location
        for idx, (_, loc_row) in enumerate(unique_locations.iterrows(), 1):
            lat = loc_row[lat_col]
            lon = loc_row[lon_col]
            
            row_num = lat_to_row[lat]
            col_num = lon_to_col[lon]

            # Filter data for this location
            df_loc = df[(df[lat_col] == lat) & (df[lon_col] == lon)].copy()

            # Pivot: time x variables
            loc_df = df_loc.pivot_table(
                index=time_col,
                columns=var_col,
                values=val_col,
                aggfunc='first'
            ).reset_index()

            # Rename Date column
            loc_df.rename(columns={time_col: 'Date'}, inplace=True)

            # Rename variables from CF names to DWD/standard names
            rename_mapping = {}
            for col in loc_df.columns:
                if col in CF_TO_DWD_NAMES:
                    rename_mapping[col] = CF_TO_DWD_NAMES[col]
            loc_df.rename(columns=rename_mapping, inplace=True)

            # Generate filename
            # Get dataset name from config
            dataset_name = getattr(self.cfg, 'dataset', 'unknown')
            
            # Create column-based folder: col_<col_number> (no zero-padding)
            col_folder = base_dir / f"{col_num}"
            col_folder.mkdir(parents=True, exist_ok=True)

            # Save file: <format>_<variables>_<dataset>.csv
            output_filename = col_folder / f"{dataset_name.upper()}_Daily_C{col_num}R{row_num}.csv"
            loc_df.to_csv(output_filename, index=False, sep='\t')

            created_files.append(str(output_filename))
            self.logger.debug(f"[{idx}/{n_locations}] Saved col_{col_num} (row {row_num}) -> {output_filename}")

        self.logger.info(f"Successfully created {len(created_files)} {format.upper()} format files in {base_dir}")
        self.logger.info(f"Grid dimensions: {n_rows} rows × {n_cols} columns")
        
        self.filename_csv = str(base_dir)
        self.current_filename = str(base_dir)

        return str(base_dir)

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
    # FAIR data object (RO-Crate)
    # ----------------------------
    def to_fair(
        self,
        ds: Optional[xr.Dataset] = None,
        output_dir: Optional[str] = None,
        title: Optional[str] = None,
        description: Optional[str] = None,
        license_url: str = "https://creativecommons.org/licenses/by/4.0/",
        creator: Optional[str] = None,
        zip_crate: bool = False,
    ) -> str:
        """Export the dataset as a FAIR Research Object Crate (RO-Crate).

        Packages the dataset into a self-describing directory that satisfies
        the FAIR principles (Findable, Accessible, Interoperable, Reusable):

        * **Findable** – unique dataset identifier & rich JSON-LD metadata
        * **Accessible** – standard CF-compliant NetCDF + open JSON-LD manifest
        * **Interoperable** – CF conventions in NetCDF; Schema.org / RO-Crate vocab
        * **Reusable** – license, provenance, variable descriptions embedded

        The output folder contains:

        .. code-block:: text

            <output_dir>/
                data.nc                  # CF-compliant NetCDF dataset
                ro-crate-metadata.json   # JSON-LD manifest (RO-Crate 1.1)
                ro-crate-preview.html    # human-readable summary (optional)

        Uses the ``rocrate`` Python package when available; falls back to a
        hand-written JSON-LD manifest so the method always works.

        Args:
            ds (xr.Dataset, optional): Dataset to package. Defaults to
                ``self.current_ds``.
            output_dir (str, optional): Destination folder. Defaults to
                ``<filename_nc_stem>_fair_crate``.
            title (str, optional): Human-readable dataset title.
            description (str, optional): Free-text description of the dataset.
            license_url (str): SPDX or Creative-Commons URL for the data
                licence. Defaults to CC-BY-4.0.
            creator (str, optional): Name of the person / organisation that
                produced the data.
            zip_crate (bool): If *True*, zip the crate directory and return
                the ``.zip`` path. Defaults to *False*.

        Returns:
            str: Absolute path to the RO-Crate directory (or ``.zip`` file).

        Raises:
            ValueError: If no dataset is available.

        Example:
            >>> extractor = ClimData(overrides=["dataset=mswx", ...])
            >>> ds = extractor.extract()
            >>> crate_path = extractor.to_fair(title="MSWX Germany 2014")
            >>> print(crate_path)      # …/mswx_fair_crate/
        """
        import datetime
        import zipfile
        import uuid

        ds = ds if ds is not None else self.current_ds
        if ds is None:
            raise ValueError(
                "No dataset available. Run extract() or upload_netcdf() first."
            )

        cfg = self.cfg

        # ── 1. Resolve output directory ─────────────────────────────────────
        if output_dir is None:
            base = getattr(self, "filename_nc", None)
            stem = Path(base).stem if base else getattr(cfg, "dataset", "climdata")
            output_dir = str(Path(base).parent / f"{stem}_fair_crate") if base else f"{stem}_fair_crate"

        crate_dir = Path(output_dir)
        crate_dir.mkdir(parents=True, exist_ok=True)

        # ── 2. Write the NetCDF payload ──────────────────────────────────────
        nc_path = crate_dir / "data.nc"

        # Embed CF global attributes before saving
        ds_out = ds.copy()
        ds_out.attrs.setdefault("Conventions", "CF-1.8")
        ds_out.attrs.setdefault("institution", creator or "climdata")
        ds_out.attrs.setdefault(
            "source",
            getattr(cfg, "dataset", "unknown"),
        )
        ds_out.attrs.setdefault(
            "history",
            f"{datetime.datetime.utcnow().isoformat()}Z climdata to_fair()",
        )
        ds_out.to_netcdf(nc_path)
        self.logger.info("FAIR payload written to %s", nc_path)

        # ── 3. Collect metadata from cfg ────────────────────────────────────
        dataset_name = getattr(cfg, "dataset", "unknown")
        variables = list(getattr(cfg, "variables", list(ds.data_vars)))

        t_start = getattr(getattr(cfg, "time_range", None), "start_date", None)
        t_end = getattr(getattr(cfg, "time_range", None), "end_date", None)
        if t_start is None and "time" in ds.coords:
            t_start = str(ds.time.values[0])[:10]
        if t_end is None and "time" in ds.coords:
            t_end = str(ds.time.values[-1])[:10]

        # spatial extent
        lat_min = lat_max = lon_min = lon_max = None
        try:
            box = getattr(cfg, "box", None)
            if box and getattr(box, "lat_min", None) is not None:
                lat_min = float(box.lat_min)
                lat_max = float(box.lat_max)
                lon_min = float(box.lon_min)
                lon_max = float(box.lon_max)
            elif getattr(cfg, "lat", None) is not None:
                lat_min = lat_max = float(cfg.lat)
                lon_min = lon_max = float(cfg.lon)
            elif "lat" in ds.coords:
                lat_min = float(ds.lat.values.min())
                lat_max = float(ds.lat.values.max())
                lon_min = float(ds.lon.values.min())
                lon_max = float(ds.lon.values.max())
        except Exception:
            pass

        # variable-level metadata
        var_descriptions = []
        varinfo = getattr(cfg, "varinfo", {})
        for v in variables:
            info = varinfo.get(v, {}) if varinfo else {}
            units = getattr(info, "units", None) or ds[v].attrs.get("units", "1") if v in ds else "1"
            long_name = ds[v].attrs.get("long_name", v) if v in ds else v
            var_descriptions.append({"@type": "PropertyValue", "name": v,
                                      "description": long_name, "unitText": units})

        title = title or f"{dataset_name.upper()} climate data extract"
        description = description or (
            f"Climate dataset extracted by climdata v{climdata.__version__} "
            f"from {dataset_name.upper()}. "
            f"Variables: {', '.join(variables)}. "
            f"Period: {t_start} to {t_end}."
        )
        crate_id = str(uuid.uuid4())
        now_iso = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")

        # ── 4. Try rocrate library; fall back to hand-crafted JSON-LD ────────
        metadata_path = crate_dir / "ro-crate-metadata.json"
        try:
            from rocrate.rocrate import ROCrate  # type: ignore
            from rocrate.model import ContextEntity  # type: ignore

            crate = ROCrate()
            crate.name = title
            crate.description = description

            crate.add_file(
                str(nc_path),
                dest_path="data.nc",
                properties={
                    "name": "data.nc",
                    "description": "CF-compliant NetCDF dataset",
                    "encodingFormat": "application/x-netcdf",
                    "variableMeasured": var_descriptions,
                },
            )

            if t_start:
                crate.root_dataset["temporalCoverage"] = f"{t_start}/{t_end}"
            if lat_min is not None:
                crate.root_dataset["spatialCoverage"] = {
                    "@type": "Place",
                    "geo": {
                        "@type": "GeoShape",
                        "box": f"{lat_min} {lon_min} {lat_max} {lon_max}",
                    },
                }
            if creator:
                person = crate.add(
                    ContextEntity(crate, f"#{creator.replace(' ', '_')}",
                                  properties={"@type": "Person", "name": creator})
                )
                crate.root_dataset["creator"] = person
            crate.root_dataset["license"] = license_url
            crate.root_dataset["datePublished"] = now_iso
            crate.root_dataset["identifier"] = crate_id
            crate.root_dataset["keywords"] = variables
            crate.root_dataset["producer"] = {
                "@type": "SoftwareApplication",
                "name": "climdata",
                "version": climdata.__version__,
                "url": "https://github.com/climdata/climdata",
            }
            crate.write(str(crate_dir))
            self.logger.info("RO-Crate written via rocrate library.")

        except ImportError:
            # ── fallback: hand-craft the JSON-LD ────────────────────────────
            self.logger.warning(
                "'rocrate' package not found. Writing plain JSON-LD metadata. "
                "Install with: pip install rocrate"
            )
            geo_shape: Optional[Dict] = None
            if lat_min is not None:
                geo_shape = {
                    "@type": "GeoShape",
                    "box": f"{lat_min} {lon_min} {lat_max} {lon_max}",
                }

            metadata = {
                "@context": [
                    "https://w3id.org/ro/crate/1.1/context",
                    {"climdata": "https://climdata.readthedocs.io/"},
                ],
                "@graph": [
                    # ── RO-Crate metadata descriptor ────────────────────
                    {
                        "@id": "ro-crate-metadata.json",
                        "@type": "CreativeWork",
                        "conformsTo": {"@id": "https://w3id.org/ro/crate/1.1"},
                        "about": {"@id": "./"},
                    },
                    # ── Root data entity (the crate itself) ─────────────
                    {
                        "@id": "./",
                        "@type": "Dataset",
                        "identifier": crate_id,
                        "name": title,
                        "description": description,
                        "datePublished": now_iso,
                        "license": {"@id": license_url},
                        "keywords": variables,
                        **({
                            "creator": {
                                "@type": "Person",
                                "name": creator,
                            }
                        } if creator else {}),
                        **({
                            "temporalCoverage": f"{t_start}/{t_end}"
                        } if t_start else {}),
                        **({
                            "spatialCoverage": {
                                "@type": "Place",
                                "geo": geo_shape,
                            }
                        } if geo_shape else {}),
                        "producer": {
                            "@type": "SoftwareApplication",
                            "name": "climdata",
                            "version": climdata.__version__,
                            "url": "https://github.com/Kaushikreddym/climdata",
                        },
                        "hasPart": [{"@id": "data.nc"}],
                    },
                    # ── File entity (the NetCDF) ─────────────────────────
                    {
                        "@id": "data.nc",
                        "@type": "File",
                        "name": "data.nc",
                        "description": "CF-compliant NetCDF dataset",
                        "encodingFormat": "application/x-netcdf",
                        "variableMeasured": var_descriptions,
                        **({
                            "temporalCoverage": f"{t_start}/{t_end}"
                        } if t_start else {}),
                    },
                    # ── License entity ───────────────────────────────────
                    {
                        "@id": license_url,
                        "@type": "CreativeWork",
                        "name": "CC-BY 4.0" if "creativecommons" in license_url else license_url,
                        "url": license_url,
                    },
                ],
            }
            with open(metadata_path, "w", encoding="utf-8") as f:
                json.dump(metadata, f, indent=2, ensure_ascii=False)

        # ── 5. Write a simple HTML preview ──────────────────────────────────
        html_path = crate_dir / "ro-crate-preview.html"
        var_rows = "".join(
            f"<tr><td><code>{v['name']}</code></td>"
            f"<td>{v['description']}</td>"
            f"<td>{v['unitText']}</td></tr>"
            for v in var_descriptions
        )
        spatial_str = (
            f"{lat_min:.3f}°N – {lat_max:.3f}°N, {lon_min:.3f}°E – {lon_max:.3f}°E"
            if lat_min is not None
            else "(not specified)"
        )
        html_content = f"""<!DOCTYPE html>
<html lang="en">
<head><meta charset="UTF-8">
<title>{title}</title>
<style>
  body{{font-family:sans-serif;max-width:900px;margin:2em auto;padding:0 1em}}
  table{{border-collapse:collapse;width:100%}}
  th,td{{border:1px solid #ccc;padding:.4em .8em;text-align:left}}
  th{{background:#f0f0f0}}
  .badge{{background:#0078d4;color:#fff;padding:.2em .6em;border-radius:4px;font-size:.85em}}
</style>
</head>
<body>
<h1>{title} <span class="badge">RO-Crate</span></h1>
<p>{description}</p>
<table>
  <tr><th>Field</th><th>Value</th></tr>
  <tr><td>Dataset</td><td>{dataset_name.upper()}</td></tr>
  <tr><td>Period</td><td>{t_start} &#8594; {t_end}</td></tr>
  <tr><td>Spatial coverage</td><td>{spatial_str}</td></tr>
  <tr><td>License</td><td><a href="{license_url}">{license_url}</a></td></tr>
  <tr><td>Produced by</td><td>climdata v{climdata.__version__}</td></tr>
  <tr><td>Crate ID</td><td><code>{crate_id}</code></td></tr>
</table>
<h2>Variables</h2>
<table>
  <tr><th>CF name</th><th>Long name</th><th>Units</th></tr>
  {var_rows}
</table>
<h2>Files</h2>
<ul>
  <li><a href="data.nc">data.nc</a> – CF-compliant NetCDF</li>
  <li><a href="ro-crate-metadata.json">ro-crate-metadata.json</a> – JSON-LD metadata</li>
</ul>
</body>
</html>
"""
        html_path.write_text(html_content, encoding="utf-8")

        # ── 6. Optionally zip ────────────────────────────────────────────────
        if zip_crate:
            zip_path = crate_dir.with_suffix(".zip")
            with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
                for file in crate_dir.rglob("*"):
                    zf.write(file, file.relative_to(crate_dir.parent))
            self.logger.info("RO-Crate zipped to %s", zip_path)
            return str(zip_path)

        self.logger.info("FAIR RO-Crate ready at %s", crate_dir)
        return str(crate_dir)

    # ----------------------------
    # Unified workflow
    # ----------------------------
    def run_workflow(self, overrides: Optional[List[str]] = None,
                     actions: Optional[List[str]] = None,
                     file: Optional[str] = None) -> WorkflowResult:
        """Execute a sequence of workflow actions.

        Args:
            overrides (list[str], optional): Hydra overrides to apply (not all actions will use these).
            actions (list[str], optional): Ordered list of actions to perform. Supported actions include: 'upload_netcdf', 'upload_csv', 'extract', 'calc_index', 'to_dataframe', 'to_csv', 'to_nc', 'to_fair'.
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

                elif action == "reproject":
                    if self.current_ds is None:
                        raise ValueError(
                            "Action 'reproject' requires a dataset, but no dataset is available. "
                            "Upload or extract a dataset before regridding."
                        )
                    self.reproject()
                    result.dataset = self.current_ds
                    result.reprojected_ds = getattr(self, "reprojected_ds", None)

                elif action == "to_fair":
                    if self.current_ds is None:
                        raise ValueError(
                            "Action 'to_fair' requires a dataset, but no dataset is available. "
                            "Upload or extract a dataset before exporting a FAIR crate."
                        )
                    result.filename = self.to_fair()

                else:
                    raise ValueError(
                        f"Unknown action '{action}'. Supported: "
                        + ", ".join(sorted(_WORKFLOW_ACTIONS))
                    )
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
        time_dim = cfg.dsinfo[cfg.dataset.upper()].get("time_dim", "time")
        lat_dim = cfg.dsinfo[cfg.dataset.upper()].get("lat_dim", "lat")
        lon_dim = cfg.dsinfo[cfg.dataset.upper()].get("lon_dim", "lon")
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

    @update_ds("reprojected_ds")
    def reproject(self, ds: xr.Dataset = None, **overrides) -> xr.Dataset:
        """Reproject / resample the dataset onto the configured target grid.

        Reads ``target_projection`` and ``target_resolution`` from the config, plus
        the secondary knobs under ``regrid``. Returns ``None`` (leaving the current
        dataset untouched) when neither target is configured.

        The units of ``target_resolution`` must match the axis units of
        ``target_projection``: a metric resolution such as ``"10 km"`` has no fixed
        angular size, so combining it with a geographic CRS raises
        :class:`~climdata.grid.units.ResolutionCRSMismatch` rather than silently
        approximating. The check runs before any data is touched.

        Args:
            ds (xr.Dataset, optional): Dataset to transform. If ``None``, uses
                ``self.current_ds``.
            **overrides: Any argument of :func:`climdata.grid.reproject`, overriding
                the configured value for this call (e.g. ``method="average"``).

        Returns:
            xr.Dataset | None: The reprojected dataset (also sets ``self.current_ds``
            and ``self.reprojected_ds``), or ``None`` if no target grid is configured.

        Raises:
            ValueError: If ``ds`` is ``None`` and ``self.current_ds`` is not set.
            ResolutionCRSMismatch: If the resolution units contradict the target CRS.

        Example:
            extractor = ClimData(overrides=[
                "target_projection=EPSG:3035", "target_resolution=10 km",
            ])
            extractor.extract()
            grid_ds = extractor.reproject()
        """
        from climdata.grid import reproject as _reproject

        cfg = self.cfg
        ds = ds if ds is not None else self.current_ds
        if ds is None:
            raise ValueError("No dataset provided and no current_ds is available.")

        target_projection = overrides.pop("target_projection", cfg.get("target_projection"))
        target_resolution = overrides.pop("target_resolution", cfg.get("target_resolution"))
        if target_projection is None and target_resolution is None and "like" not in overrides:
            self.logger.info(
                "No target_projection or target_resolution configured; skipping reprojection."
            )
            return None

        regrid_cfg = cfg.get("regrid") or {}
        kwargs = {
            "method": regrid_cfg.get("method", "bilinear"),
            "align": regrid_cfg.get("align", True),
            "bounds": regrid_cfg.get("bounds"),
            "engine": regrid_cfg.get("engine", "rasterio"),
        }
        kwargs.update(overrides)

        # Dimension names are already declared per dataset in parameters.yaml.
        try:
            dsinfo = cfg.dsinfo[cfg.dataset.upper()]
            kwargs.setdefault(
                "dsinfo_dims", (dsinfo.get("lon_dim", "lon"), dsinfo.get("lat_dim", "lat"))
            )
        except Exception:
            pass

        if kwargs.get("bounds") is not None:
            kwargs["bounds"] = [float(v) for v in kwargs["bounds"]]
        if isinstance(kwargs.get("method"), DictConfig):
            kwargs["method"] = dict(kwargs["method"])

        self.logger.info(
            "Reprojecting to %s at %s", target_projection or "source CRS", target_resolution
        )
        return _reproject(ds, target_projection, target_resolution, **kwargs)

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