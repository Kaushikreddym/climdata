"""
HOSTRADA – Hourly Station-based and Raster-based Terrestrial Atmospheric Dataset
DWD Open Data server:
  https://opendata.dwd.de/climate_environment/CDC/grids_germany/hourly/hostrada/

File naming convention (monthly files):
  {var}_1hr_HOSTRADA-v1-0_BE_gn_{YYYYMMDDHH}-{YYYYMMDDHH}.nc

Currently implemented variables
--------------------------------
  rsds  – Surface downwelling shortwave radiation (W m⁻²)
          subfolder: radiation_downwelling/
"""

from __future__ import annotations

import calendar
import os
from datetime import datetime
from typing import Dict, List, Optional, Tuple

import numpy as np
import requests
import xarray as xr
from omegaconf import DictConfig


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _last_hour_of_month(year: int, month: int) -> str:
    """Return the last time-stamp string (YYYYMMDDHH) for a monthly file."""
    last_day = calendar.monthrange(year, month)[1]
    return f"{year}{month:02d}{last_day:02d}23"


def _first_hour_of_month(year: int, month: int) -> str:
    return f"{year}{month:02d}0100"


def fetch_hostrada(cfg: DictConfig, variable: str) -> None:
    """Download HOSTRADA NetCDF files for *variable* over the configured time range.

    Files are stored as::
        <data_dir>/hostrada/<VARIABLE>/<filename>.nc
    """
    provider = cfg.dataset.lower()          # "hostrada"
    param_info = cfg.dsinfo[provider]["variables"][variable]
    base_url: str = param_info["base_url"]
    prefix: str = param_info["prefix"]      # e.g. "rsds"
    version: str = param_info["version"]    # e.g. "v1-0"

    start_dt = datetime.fromisoformat(cfg.time_range.start_date)
    end_dt = datetime.fromisoformat(cfg.time_range.end_date)

    var_dir = os.path.join(cfg.data_dir, provider, variable.upper())
    os.makedirs(var_dir, exist_ok=True)

    year = start_dt.year
    month = start_dt.month

    while (year, month) <= (end_dt.year, end_dt.month):
        t_start = _first_hour_of_month(year, month)
        t_end = _last_hour_of_month(year, month)
        filename = f"{prefix}_1hr_HOSTRADA-{version}_BE_gn_{t_start}-{t_end}.nc"
        file_url = f"{base_url}{filename}"
        local_path = os.path.join(var_dir, filename)

        if os.path.exists(local_path):
            print(f"✔️  Already exists: {local_path}")
        else:
            print(f"⬇️  Downloading: {file_url}")
            head = requests.head(file_url, timeout=10)
            if head.status_code != 200:
                print(f"⚠️  Not found on server: {file_url} (HTTP {head.status_code})")
            else:
                try:
                    resp = requests.get(file_url, stream=True, timeout=60)
                    resp.raise_for_status()
                    with open(local_path, "wb") as fh:
                        for chunk in resp.iter_content(chunk_size=8192):
                            fh.write(chunk)
                    print(f"✅ Saved: {local_path}")
                except requests.HTTPError as exc:
                    print(f"⚠️  Download failed: {file_url} — {exc}")

        # advance to next month
        if month == 12:
            year += 1
            month = 1
        else:
            month += 1


# ---------------------------------------------------------------------------
# Dataset class
# ---------------------------------------------------------------------------

class HOSTRADAmirror:
    """Load and subset HOSTRADA gridded hourly data.

    The data uses a regular lat/lon grid (EPSG:4326) so spatial subsetting is
    straightforward index slicing.

    Parameters
    ----------
    cfg : DictConfig
        Hydra/OmegaConf config produced by ``load_config``.

    Usage
    -----
    ::

        from climdata.datasets.HOSTRADA import HOSTRADAmirror

        mirror = HOSTRADAmirror(cfg)
        mirror.extract(box={"lat_min": 47, "lat_max": 55,
                             "lon_min": 6, "lon_max": 15})
        ds = mirror.load("rsds")
    """

    def __init__(self, cfg: DictConfig) -> None:
        self.cfg = cfg
        self.dataset: Optional[xr.Dataset] = None
        self.variables: List[str] = list(cfg.variables)

        self._extract_mode: Optional[str] = None
        self._extract_params = None

    # ------------------------------------------------------------------
    # File discovery / download
    # ------------------------------------------------------------------

    def fetch(self, variable: str) -> List[str]:
        """Download files for *variable* and return list of local paths."""
        fetch_hostrada(self.cfg, variable)
        return self._find_files(variable)

    def _find_files(self, variable: str) -> List[str]:
        """Return sorted list of existing local NetCDF files for *variable*."""
        import glob

        provider = self.cfg.dataset.lower()
        var_dir = os.path.join(self.cfg.data_dir, provider, variable.upper())

        start_dt = datetime.fromisoformat(self.cfg.time_range.start_date)
        end_dt = datetime.fromisoformat(self.cfg.time_range.end_date)

        files: List[str] = []
        year, month = start_dt.year, start_dt.month

        while (year, month) <= (end_dt.year, end_dt.month):
            pattern = os.path.join(
                var_dir,
                f"{variable}_1hr_HOSTRADA-*_{year}{month:02d}*.nc",
            )
            found = sorted(glob.glob(pattern))
            files.extend(found)

            if month == 12:
                year += 1
                month = 1
            else:
                month += 1

        return files

    # ------------------------------------------------------------------
    # Extraction intent setter (mirrors HYRAS API)
    # ------------------------------------------------------------------

    def extract(
        self,
        *,
        point: Optional[Tuple[float, float]] = None,
        box: Optional[Dict[str, float]] = None,
        shapefile=None,
        buffer_km: float = 0.0,
    ) -> "HOSTRADAmirror":
        if point is not None:
            self._extract_mode = "point"
            self._extract_params = point          # (lon, lat)
        elif box is not None:
            self._extract_mode = "box"
            self._extract_params = box
        elif shapefile is not None:
            self._extract_mode = "shapefile"
            self._extract_params = shapefile
        else:
            self._extract_mode = None
            self._extract_params = None
        return self

    # ------------------------------------------------------------------
    # Load
    # ------------------------------------------------------------------

    def load(
        self,
        variable: str,
        chunking: Optional[Dict] = None,
        use_dask: bool = True,
    ) -> xr.Dataset:
        """Fetch (if needed) and open the dataset for *variable*.

        HOSTRADA files use a Lambert Conformal Conic projection (EPSG:3034)
        with projected-metre dimensions ``Y`` / ``X`` and 2-D auxiliary
        coordinates ``lat(Y, X)`` / ``lon(Y, X)``.  Spatial subsetting is
        therefore done via a boolean mask over those 2-D arrays and then
        converting to a tight rectangular ``isel`` slice.

        Parameters
        ----------
        variable : str
            CF variable name, e.g. ``"rsds"``.
        chunking : dict, optional
            Chunk specification passed to ``dset.chunk()``.
        use_dask : bool
            Whether to use dask-backed lazy loading.

        Returns
        -------
        xr.Dataset
        """
        files = self.fetch(variable)

        if not files:
            raise FileNotFoundError(
                f"No HOSTRADA files found for variable '{variable}' "
                f"in the requested time range "
                f"({self.cfg.time_range.start_date} – {self.cfg.time_range.end_date}). "
                f"Check data_dir and that the files have been downloaded."
            )

        mode = self._extract_mode

        # Compute Y/X index bounds once from the first file so we can use
        # fast isel in the preprocess callbacks.
        y_slice, x_slice = self._compute_yx_slices(files[0], mode)

        def _preprocess(ds: xr.Dataset) -> xr.Dataset:
            # Apply Y/X rectangular slice (works whether mode is box/shapefile/point)
            if y_slice is not None and x_slice is not None:
                ds = ds.isel(Y=y_slice, X=x_slice)
            return ds

        dset = xr.open_mfdataset(
            files,
            combine="nested",
            concat_dim="time",
            preprocess=_preprocess,
            engine="netcdf4",
            parallel=False,
        )

        if use_dask and chunking:
            dset = dset.chunk(chunking)

        dset = self._apply_time_subset(dset)
        self.dataset = dset
        return dset

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _compute_yx_slices(self, sample_file: str, mode: Optional[str]):
        """Return (y_slice, x_slice) integer index slices for the requested
        spatial subset, computed from the 2-D ``lat``/``lon`` auxiliary
        coordinates in *sample_file*.

        Returns ``(None, None)`` when no spatial subsetting is requested.
        """
        import numpy as np

        if mode is None:
            return None, None

        ds_sample = xr.open_dataset(sample_file, engine="netcdf4")
        lat2d = ds_sample["lat"].values   # shape (Y, X)
        lon2d = ds_sample["lon"].values   # shape (Y, X)
        ds_sample.close()

        if mode == "point":
            lon_pt, lat_pt = self._extract_params
            dist = (lat2d - lat_pt) ** 2 + (lon2d - lon_pt) ** 2
            iy, ix = np.unravel_index(int(np.argmin(dist)), dist.shape)
            return slice(iy, iy + 1), slice(ix, ix + 1)

        elif mode in ("box", "shapefile"):
            if mode == "box":
                box = self._extract_params
                lat_min, lat_max = box["lat_min"], box["lat_max"]
                lon_min, lon_max = box["lon_min"], box["lon_max"]
            else:
                import geopandas as gpd
                if isinstance(self._extract_params, str):
                    gdf = gpd.read_file(self._extract_params)
                else:
                    gdf = self._extract_params
                lon_min, lat_min, lon_max, lat_max = gdf.total_bounds

            mask = (
                (lat2d >= lat_min) & (lat2d <= lat_max) &
                (lon2d >= lon_min) & (lon2d <= lon_max)
            )
            if not mask.any():
                raise ValueError(
                    f"Bounding box lat=[{lat_min},{lat_max}] "
                    f"lon=[{lon_min},{lon_max}] does not intersect the HOSTRADA grid."
                )
            y_idx, x_idx = np.where(mask)
            y_slice = slice(int(y_idx.min()), int(y_idx.max()) + 1)
            x_slice = slice(int(x_idx.min()), int(x_idx.max()) + 1)
            return y_slice, x_slice

        return None, None

    def _apply_time_subset(self, ds: xr.Dataset) -> xr.Dataset:
        start = getattr(self.cfg.time_range, "start_date", None)
        end = getattr(self.cfg.time_range, "end_date", None)
        if start or end:
            try:
                ds = ds.sel(time=slice(start, end))
            except Exception as exc:
                print(f"⚠️  Time subsetting failed: {exc}")
        return ds
