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
import time
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
from climdata.utils import provenance as prov
    
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
    "impute", "reproject", "to_csv", "to_nc",
)

#: ``output.fair`` defaults, mirrored from ``conf/config.yaml`` for configs
#: copied before the key existed.
_FAIR_DEFAULTS = {
    "enabled": True, "title": None, "description": None,
    "creator": None, "license": None, "zip": False,
}

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
        crate_dir (str | None): FAIR crate folder the last write went into.
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
    crate_dir: Optional[str] = None

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
def track_step(func):
    """Record each call of a processing method in ``self.provenance``.

    Logs the step name, UTC start/end time, duration, status, the parameters
    the method declared via :meth:`ClimateExtractor._note`, any input files
    declared via :meth:`ClimateExtractor._note_input`, and a structural summary
    of the return value. A method returning ``None`` is logged as ``skipped``;
    one that raises is logged as ``failed`` and the exception re-raised.

    Apply it beneath :func:`update_ds` / :func:`update_df` so the summary
    describes what the method itself returned.

    Args:
        func (Callable): A ``ClimateExtractor`` method.

    Returns:
        Callable: The wrapped method.
    """
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        if getattr(self, "provenance", None) is None:
            # Instance built without __init__ (e.g. in tests): nothing to log to.
            return func(self, *args, **kwargs)
        outer = self._step
        self._step = {"action": func.__name__, "started": prov.now_iso(),
                      "t0": time.perf_counter(), "parameters": {}, "inputs": [], "results": []}
        try:
            out = func(self, *args, **kwargs)
        except Exception as exc:
            self._log_step("failed", error=f"{type(exc).__name__}: {exc}")
            raise
        finally:
            step, self._step = self._step, outer
        self._log_step("skipped" if out is None else "completed", output=out, step=step)
        return out
    return wrapper
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

        # Processing record: one dict per tracked step (see track_step).
        self.provenance: List[dict] = []
        self._step = None
        self._crates: Dict[str, dict] = {}     # crate folder -> id, created, files
        self.crate_dir: Optional[str] = None

        # Automatically load config on init
        self.load_config(overrides)
        self.cfg = self.preprocess_aoi(self.cfg)

        # instance logger for this extractor
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")

    # ----------------------------
    # Provenance
    # ----------------------------
    def _note(self, **params):
        """Attach parameters to the step currently being recorded.

        Called from inside a :func:`track_step` method; a no-op outside one.
        """
        if self._step is not None:
            self._step["parameters"].update(prov.jsonable(params))

    def _note_input(self, path):
        """Record a file the current step read, with its size and checksum."""
        if self._step is not None:
            self._step["inputs"].append(prov.file_record(path))

    def _make_entry(self, status="completed", output=None, error=None, step=None):
        """Build a provenance entry for ``step`` (default: the step now running).

        Returns ``None`` outside a tracked step.
        """
        step = step if step is not None else self._step
        if step is None:
            return None
        entry = {
            "step": len(self.provenance) + 1,
            "action": step["action"],
            "status": status,
            "started": step["started"],
            "ended": prov.now_iso(),
            "duration_s": round(time.perf_counter() - step["t0"], 3),
            "parameters": dict(step["parameters"]),
        }
        if step["inputs"]:
            entry["inputs"] = list(step["inputs"])
        if step["results"]:
            entry["results"] = list(step["results"])
        if output is not None:
            entry["output"] = prov.summarize_output(output)
        if error:
            entry["error"] = error
        return entry

    def _log_step(self, status, output=None, error=None, step=None):
        """Append an entry for ``step`` (default: the step now running) to ``self.provenance``."""
        entry = self._make_entry(status, output=output, error=error, step=step)
        self.provenance.append(entry)
        return entry

    def provenance_record(self) -> List[dict]:
        """The processing steps recorded so far, oldest first.

        Each entry holds ``step``, ``action``, ``status`` (``completed``,
        ``skipped`` or ``failed``), ``started``/``ended`` UTC timestamps,
        ``duration_s``, the ``parameters`` that step used, ``inputs`` (files
        read, with SHA-256), ``results`` (files written into a crate) and an
        ``output`` summary. :meth:`to_nc` and :meth:`to_csv` write this list to
        the crate's ``meta.json``.

        Returns:
            list[dict]: A copy of ``self.provenance``.
        """
        return [dict(e) for e in self.provenance]

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
        elif (cfg.get("box") or {}).get("lat_min") is not None:
            b = cfg.box
            lat_str = f"{b.lat_min}_{b.lat_max}"
            lon_str = f"{b.lon_min}_{b.lon_max}"
            lat_range = f"{b.lat_min}-{b.lat_max}"
            lon_range = f"{b.lon_min}-{b.lon_max}"
        elif self.current_ds is not None and prov.spatial_coverage(self.current_ds):
            # No spatial selection configured (e.g. an uploaded file): describe the data.
            sc = prov.spatial_coverage(self.current_ds)
            lat_min, lat_max = f"{sc['lat_min']:g}", f"{sc['lat_max']:g}"
            lon_min, lon_max = f"{sc['lon_min']:g}", f"{sc['lon_max']:g}"
            lat_str, lon_str = f"{lat_min}_{lat_max}", f"{lon_min}_{lon_max}"
            lat_range, lon_range = f"{lat_min}-{lat_max}", f"{lon_min}-{lon_max}"
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

        # A conf/ copied into the working directory by an older climdata has no
        # output.fair; fill in the defaults so FAIR output behaves the same.
        from omegaconf import OmegaConf, open_dict
        if self.cfg.get("output") is not None:
            with open_dict(self.cfg):
                self.cfg.output.fair = OmegaConf.merge(
                    _FAIR_DEFAULTS, self.cfg.output.get("fair") or {}
                )
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
    @track_step
    def upload_netcdf(self, nc_file: str) -> xr.Dataset:
        """Load a NetCDF file into an xarray.Dataset and update file metadata.

        Args:
            nc_file (str): Path to the NetCDF file to open.

        Returns:
            xr.Dataset: The loaded dataset (also sets ``self.current_ds``).
        """
        if not os.path.exists(nc_file):
            raise FileNotFoundError(f"{nc_file} does not exist")

        self._note(file=os.path.abspath(nc_file))
        self._note_input(nc_file)
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
    @track_step
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

        self._note(file=os.path.abspath(csv_file))
        self._note_input(csv_file)
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
    @track_step
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
        units = {var: cfg.varinfo[var].units for var in cfg.variables}
        for var in cfg.variables:
            ds[var] = xclim.core.units.convert_units_to(ds[var], units[var])

        self._note(
            dataset=cfg.dataset,
            provider_class=type(self.dataset_class).__name__,
            variables=list(cfg.variables),
            start_date=str(cfg.time_range.start_date),
            end_date=str(cfg.time_range.end_date),
            selection=extract_kwargs,
            region=cfg.region,
            converted_units=units,
            data_dir=os.path.abspath(cfg.data_dir),
        )
        if dataset_upper in ("CMIP", "CMIP_W5E5", "NEXGDDP"):
            self._note(**{k: cfg.get(k) for k in ("source_id", "experiment_id", "table_id")
                          if cfg.get(k) is not None})

        # ds = ds.compute()

        return ds
    # ----------------------------
    # Compute extreme index
    # ----------------------------
    @update_ds(attr_name='index_ds')
    @track_step
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

        idx_def = cfg.extinfo.get("indices", {}).get(cfg.index, {})
        self._note(index=cfg.index, function=idx_def.get("function"),
                   input_variables=idx_def.get("variables"), args=idx_def.get("args"))
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
    @track_step
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
    @track_step
    def to_csv(self, df: Optional[pd.DataFrame] = None, filename: Optional[str] = None,
               format: str = "default", fair: Optional[bool] = None) -> str:
        """Save a DataFrame to CSV, by default inside a FAIR crate.

        With FAIR output on (the default, ``output.fair.enabled``), the file is
        written into a crate folder named after it — ``out/x.csv`` becomes
        ``out/x_crate/x.csv`` — next to ``meta.json`` (every processing step,
        the configuration, software versions, checksums),
        ``ro-crate-metadata.json`` and ``ro-crate-preview.html``. The NetCDF
        from :meth:`to_nc` shares the stem, so both go into the same crate.

        Args:
            df (pd.DataFrame, optional): DataFrame to save. Defaults to ``self.current_df``.
            filename (str, optional): Output filename/directory. Defaults to ``self.filename_csv``.
            format (str): Output format. Options:
                - 'default': Long-form (single file)
                - 'simplace': SIMPLACE format (splits by location, tab-separated)
                - 'monica': MONICA format (splits by location, tab-separated)
                Defaults to 'default'.
            fair (bool, optional): Write into a FAIR crate. Defaults to
                ``output.fair.enabled`` (``True``); ``False`` writes the bare
                file at ``filename``.

        Returns:
            str: The path of the written CSV file. For SIMPLACE/MONICA, the base directory.
            The crate folder, when one is written, is in ``self.crate_dir``.

        Raises:
            ValueError: If format is not supported or required data is missing.
        """
        df = df if df is not None else self.current_df
        format_lower = format.lower()
        fair = self._fair_enabled(fair)
        self._note(format=format_lower, fair=fair)
        if format_lower not in ("default", "simplace", "monica"):
            raise ValueError(f"Unsupported format: {format}. Supported formats: 'default', 'simplace', 'monica'")

        if format_lower in ("simplace", "monica"):
            crate = None
            if fair:
                if filename is not None:
                    crate, filename = self._crate_paths(filename)
                else:
                    base = getattr(self, "filename_csv", None) or f"{format_lower}_{self.cfg.dataset}"
                    crate, _ = self._crate_paths(base)
                    filename = crate / format_lower
            out = self._convert_to_gridded_format(df, None if filename is None else str(filename), format=format_lower)
        else:
            filename = filename or getattr(self, "filename_csv", None)
            if filename is None:
                raise ValueError("No filename provided and filename_csv is not set")
            crate, path = self._crate_paths(filename) if fair else (None, Path(filename))
            path.parent.mkdir(parents=True, exist_ok=True)

            df.to_csv(path, index=False, sep='\t')
            self.filename_csv = str(path)
            self.current_filename = str(path)
            self.logger.info(f"DataFrame saved to CSV file ({format} format): {self.current_filename}")
            out = str(path)

        if crate is not None:
            variables = sorted(map(str, df["variable"].dropna().unique())) if "variable" in df else []
            self._finalize_crate(crate, Path(out), variables, ds=self.current_ds, df=df)
        return out

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

    @track_step
    def to_nc(self, ds: Optional[xr.Dataset] = None, filename: Optional[str] = None,
              fair: Optional[bool] = None) -> str:
        """Save an xarray Dataset to NetCDF, by default inside a FAIR crate.

        With FAIR output on (the default, ``output.fair.enabled``), the file is
        written into a crate folder named after it:

        .. code-block:: text

            <provider>_<variables>_LAT…_LON…_<start>_<end>_crate/
                <provider>_<variables>_….nc   # CF `history` lists every step
                meta.json                     # full provenance record
                ro-crate-metadata.json        # RO-Crate 1.1, one CreateAction per step
                ro-crate-preview.html         # human-readable summary

        ``meta.json`` holds every step recorded in :attr:`provenance` with its
        parameters, inputs, timings and output; the configuration the run used;
        software versions; and checksums of the data files. Title, creator,
        licence and zipping come from ``output.fair`` in the configuration.

        Notes:
            - If ``ds`` is ``None``: save ``current_ds``.
            - If ``filename`` is ``None``: use ``self.filename_nc``.
            - Creates directories if needed and updates ``self.filename_nc`` and ``self.current_filename``.

        Args:
            ds (xr.Dataset, optional): Dataset to save. If ``None``, uses ``self.current_ds``.
            filename (str, optional): Output filename. Defaults to ``self.filename_nc``.
            fair (bool, optional): Write into a FAIR crate. Defaults to
                ``output.fair.enabled`` (``True``); ``False`` writes the bare
                file at ``filename``.

        Returns:
            str: The path of the written NetCDF file. The crate folder, when one
            is written, is in ``self.crate_dir``.
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

        fair = self._fair_enabled(fair)
        self._note(fair=fair)
        crate, path = self._crate_paths(filename) if fair else (None, Path(filename))
        path.parent.mkdir(parents=True, exist_ok=True)

        # -------------------------------
        # 3. Save to NetCDF
        # -------------------------------
        if crate is not None:
            variables = list(ds.data_vars)
            fields = self._fair_fields(variables)
            self._nc_with_provenance(ds, self._crate_entry(crate)["id"], fields).to_netcdf(path)
            self._finalize_crate(crate, path, variables, ds=ds)
        else:
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
    # FAIR crate (written by to_nc / to_csv)
    # ----------------------------
    def _fair_enabled(self, fair: Optional[bool]) -> bool:
        """Whether a write goes into a FAIR crate: the call's ``fair``, else ``output.fair.enabled``."""
        if fair is not None:
            return bool(fair)
        from omegaconf import OmegaConf
        return bool(OmegaConf.select(self.cfg, "output.fair.enabled", default=True))

    @staticmethod
    def _crate_paths(path):
        """Map an output path to ``(crate_dir, path inside the crate)``.

        ``out/x.nc`` becomes ``out/x_crate/x.nc``. A path already inside its
        crate is returned unchanged, so writing twice never nests crates. The
        ``.nc`` and ``.csv`` filenames share a stem, so both land in one crate.
        """
        p = Path(path)
        if p.parent.name == f"{p.stem}_crate":
            return p.parent, p
        crate = p.parent / f"{p.stem}_crate"
        return crate, crate / p.name

    def _crate_entry(self, crate_dir: Path) -> dict:
        """The crate's identity for this extractor: ID, creation time, files written."""
        import uuid
        key = str(Path(crate_dir).resolve())
        if key not in self._crates:
            self._crates[key] = {"id": str(uuid.uuid4()), "created": prov.now_iso(), "files": []}
        return self._crates[key]

    def _fair_fields(self, variables: List[str]) -> dict:
        """Title, description, creator and licence for a crate, from ``output.fair``."""
        from omegaconf import OmegaConf
        cfg = self.cfg
        dataset_name = str(getattr(cfg, "dataset", None) or "unknown")
        try:
            dsinfo = dict(cfg.dsinfo[dataset_name.upper()])
        except Exception:
            dsinfo = {}
        fair = lambda key: OmegaConf.select(cfg, f"output.fair.{key}", default=None)
        tc = prov.temporal_coverage(self.current_ds) if self.current_ds is not None else None
        title = fair("title") or f"{dataset_name.upper()} {', '.join(variables)}"
        description = fair("description") or (
            f"{', '.join(variables)} derived from {dataset_name.upper()} with climdata"
            + (f", {tc['start']} to {tc['end']}" if tc else "")
            + ". See meta.json for every processing step."
        )
        return {
            "title": title,
            "description": description,
            "creator": fair("creator"),
            # Derived data inherits the source's terms unless configured otherwise.
            "license": fair("license") or dsinfo.get("license")
                       or "https://creativecommons.org/licenses/by/4.0/",
            "dataset_name": dataset_name,
            "dsinfo": {k: v for k, v in dsinfo.items() if k not in ("variables", "params")},
        }

    def _nc_with_provenance(self, ds: xr.Dataset, crate_id: str, fields: dict) -> xr.Dataset:
        """Copy of ``ds`` with CF global attributes and the processing history."""
        version = getattr(climdata, "__version__", "unknown")
        steps = self.provenance_record()
        inflight = self._make_entry()
        if inflight:
            steps.append(inflight)
        out = ds.copy()
        out.attrs.setdefault("Conventions", "CF-1.8")
        out.attrs.setdefault("institution", fields["creator"] or "climdata")
        out.attrs.setdefault("source", fields["dataset_name"])
        out.attrs["title"] = fields["title"]
        out.attrs["license"] = fields["license"]
        out.attrs["climdata_version"] = version
        out.attrs["climdata_crate_id"] = crate_id
        out.attrs["provenance"] = "Full processing record in meta.json in the same folder."
        # CF convention: newest entry first, prepended to any existing history.
        new = list(reversed(prov.history_lines(steps, version)))
        old = ds.attrs.get("history")
        out.attrs["history"] = "\n".join(new + ([old] if old else []))
        return out

    def _finalize_crate(self, crate_dir: Path, written: Path, variables: List[str],
                        ds: Optional[xr.Dataset] = None, df: Optional[pd.DataFrame] = None) -> str:
        """Register ``written`` in its crate and (re)write the crate metadata.

        Writes ``meta.json`` (the full record), ``ro-crate-metadata.json`` and
        ``ro-crate-preview.html``, all from the provenance log including the
        step currently running. Called by :meth:`to_nc` and :meth:`to_csv`
        after they write their data; a second write into the same crate adds
        its file and refreshes the metadata.

        Returns:
            str: The crate folder.
        """
        import zipfile
        from omegaconf import OmegaConf

        crate_dir = Path(crate_dir)
        entry = self._crate_entry(crate_dir)
        new = [written] if written.is_file() else sorted(f for f in written.rglob("*") if f.is_file())
        for f in new:
            rel = str(f.relative_to(crate_dir))
            if rel not in entry["files"]:
                entry["files"].append(rel)
            if self._step is not None and rel not in self._step["results"]:
                self._step["results"].append(rel)

        fields = self._fair_fields(variables)
        steps = self.provenance_record()
        inflight = self._make_entry()
        if inflight:
            steps.append(inflight)

        if ds is not None:
            data = prov.summarize_dataset(ds)
        else:
            data = {"type": "pandas.DataFrame",
                    "variables": [{"name": v} for v in variables],
                    "rows": int(len(df)) if df is not None else None}

        kinds = {".nc": ("application/x-netcdf", "NetCDF dataset (CF-1.8)"),
                 ".csv": ("text/tab-separated-values", "Tab-separated table")}
        files = []
        for rel in entry["files"]:
            fmt, desc = kinds.get(Path(rel).suffix, ("application/octet-stream", "Data file"))
            files.append(dict(prov.file_record(crate_dir / rel, name=rel),
                              role="data", encodingFormat=fmt, description=desc))
        files += [
            {"path": "meta.json", "role": "metadata", "description": "Full processing record"},
            {"path": "ro-crate-metadata.json", "role": "metadata", "description": "RO-Crate 1.1 JSON-LD manifest"},
            {"path": "ro-crate-preview.html", "role": "metadata", "description": "Human-readable summary"},
        ]

        meta = {
            "schema_version": prov.META_SCHEMA_VERSION,
            "crate_id": entry["id"],
            "title": fields["title"],
            "description": fields["description"],
            "created": entry["created"],
            "updated": prov.now_iso(),
            "creator": fields["creator"],
            "license": fields["license"],
            "source_dataset": fields["dataset_name"].upper(),
            "source_dataset_info": prov.jsonable(fields["dsinfo"]),
            "data": data,
            "files": files,
            "processing": steps,
            "config": prov.config_snapshot(self.cfg),
            "software": prov.software_environment(),
        }
        (crate_dir / "meta.json").write_text(
            json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
        (crate_dir / "ro-crate-metadata.json").write_text(
            json.dumps(prov.build_ro_crate(meta), indent=2, ensure_ascii=False), encoding="utf-8")
        (crate_dir / "ro-crate-preview.html").write_text(prov.render_preview(meta), encoding="utf-8")

        if OmegaConf.select(self.cfg, "output.fair.zip", default=False):
            zip_path = crate_dir.with_suffix(".zip")
            with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
                for f in crate_dir.rglob("*"):
                    zf.write(f, f.relative_to(crate_dir.parent))
            self.logger.info("Crate zipped to %s", zip_path)

        self.crate_dir = str(crate_dir)
        self.logger.info("FAIR crate updated at %s", crate_dir)
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
            actions (list[str], optional): Ordered list of actions to perform. Supported actions include: 'upload_netcdf', 'upload_csv', 'extract', 'calc_index', 'impute', 'reproject', 'to_csv', 'to_nc'. 'to_csv' and 'to_nc' write into a FAIR crate unless ``output.fair.enabled=false``.
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
                    result.crate_dir = getattr(self, "crate_dir", None)

                elif action == "to_nc":
                    if self.current_ds is None:
                        raise ValueError(
                            "Action 'to_nc' requires a dataset, but no dataset is available. "
                            "Upload or extract a dataset before saving to NetCDF."
                        )
                    result.filename = self.to_nc()
                    result.crate_dir = getattr(self, "crate_dir", None)

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
    @track_step
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
        self._note(method=method, normalize=normalize,
                   dims=[time_dim, lat_dim, lon_dim], method_config=impute_cfg[method])
        recovered = imputer.impute()

        # merge imputed variables back into original dataset if we operated on a subset

        ds_out = recovered

        # Return dataset (decorator will set current_ds and impute_ds and generate filenames)
        return ds_out

    @update_ds("reprojected_ds")
    @track_step
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

        self._note(target_projection=target_projection, target_resolution=target_resolution,
                   **{k: v for k, v in kwargs.items() if k != "like"})
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