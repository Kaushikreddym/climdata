"""HOSTRADA hourly gridded data for Germany, from the DWD Open Data server.

HOSTRADA (Hourly Station-based and Raster-based Terrestrial Atmospheric Dataset)
is published as one NetCDF per variable per month at
https://opendata.dwd.de/climate_environment/CDC/grids_germany/hourly/hostrada/,
named ``{var}_1hr_HOSTRADA-v1-0_BE_gn_{YYYYMMDDHH}-{YYYYMMDDHH}.nc``.

The grid is Lambert Conformal Conic (EPSG:3034) with projected-metre ``Y``/``X``
dimensions and two-dimensional ``lat(Y, X)`` / ``lon(Y, X)`` auxiliary
coordinates. Geographic subsetting therefore cannot be a coordinate slice; it
goes through a mask over those 2-D arrays, reduced to an integer ``isel``
rectangle (see :meth:`HOSTRADAmirror._compute_yx_slices`). The rectangle is the
bounding box of the mask, so a request that is not axis-aligned in the
projection returns somewhat more area than asked for.

Which variables are available depends on ``cfg.dsinfo.HOSTRADA.variables``;
``rsds`` (surface downwelling shortwave radiation) is the reference entry.

Example:
    >>> mirror = HOSTRADAmirror(cfg)                             # doctest: +SKIP
    >>> mirror.extract(box={"lat_min": 47, "lat_max": 55,        # doctest: +SKIP
    ...                     "lon_min": 6, "lon_max": 15})
    >>> ds = mirror.load("rsds")                                 # doctest: +SKIP
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
    """Build the closing timestamp of a monthly HOSTRADA filename.

    Args:
        year (int): Calendar year.
        month (int): Calendar month, 1-12.

    Returns:
        str: ``YYYYMMDDHH`` for hour 23 of the month's last day, e.g.
        ``"2020013123"``.
    """
    last_day = calendar.monthrange(year, month)[1]
    return f"{year}{month:02d}{last_day:02d}23"


def _first_hour_of_month(year: int, month: int) -> str:
    """Build the opening timestamp of a monthly HOSTRADA filename.

    Args:
        year (int): Calendar year.
        month (int): Calendar month, 1-12.

    Returns:
        str: ``YYYYMMDDHH`` for hour 00 of the first day, e.g. ``"2020010100"``.
    """
    return f"{year}{month:02d}0100"


def fetch_hostrada(cfg: DictConfig, variable: str) -> None:
    """Download the monthly HOSTRADA files covering the configured period.

    Walks the range one month at a time, storing each file under
    ``<data_dir>/HOSTRADA/<VARIABLE>/``. Files already on disk are skipped, so an
    interrupted run resumes. A month the server does not have is reported and
    skipped rather than raised — the caller finds out from the shorter file list
    that :meth:`HOSTRADAmirror.fetch` returns.

    Args:
        cfg (DictConfig): Configuration with ``dataset``, ``data_dir``,
            ``time_range`` and ``dsinfo.HOSTRADA.variables``.
        variable (str): CF variable name, e.g. ``"rsds"``.

    Returns:
        None: Files are written to disk.

    Raises:
        KeyError: If ``variable`` is not declared in
            ``cfg.dsinfo.HOSTRADA.variables``.
    """
    provider = cfg.dataset.upper()          # "HOSTRADA"
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
    """Download, open and subset HOSTRADA hourly grids.

    Follows the same call order as :class:`~climdata.datasets.MSWX.MSWXmirror`:
    :meth:`extract` records the region of interest, then :meth:`load` downloads
    what is missing and applies the subset per file as it opens.

    The grid is projected with 2-D latitude/longitude coordinates, so the subset
    is an integer ``Y``/``X`` rectangle rather than a coordinate slice — see the
    module docstring for what that implies. Output keeps the projected ``Y``/``X``
    dimensions; pass it through :func:`climdata.reproject` for a geographic grid.

    Attributes:
        cfg (DictConfig): Hydra configuration.
        variables (list[str]): ``cfg.variables``.
        dataset (xr.Dataset | None): The loaded dataset, set by :meth:`load`.

    Example:
        >>> mirror = HOSTRADAmirror(cfg)                          # doctest: +SKIP
        >>> mirror.extract(box={"lat_min": 47, "lat_max": 55,     # doctest: +SKIP
        ...                     "lon_min": 6, "lon_max": 15})
        >>> ds = mirror.load("rsds", chunking={"time": "auto"})   # doctest: +SKIP
    """

    def __init__(self, cfg: DictConfig) -> None:
        """Bind a configuration.

        Args:
            cfg (DictConfig): Configuration with ``dataset``, ``data_dir``,
                ``variables``, ``time_range`` and ``dsinfo.HOSTRADA``.
        """
        self.cfg = cfg
        self.dataset: Optional[xr.Dataset] = None
        self.variables: List[str] = list(cfg.variables)

        self._extract_mode: Optional[str] = None
        self._extract_params = None

    # ------------------------------------------------------------------
    # File discovery / download
    # ------------------------------------------------------------------

    def fetch(self, variable: str) -> List[str]:
        """Download the files for a variable and list what is now on disk.

        Args:
            variable (str): CF variable name, e.g. ``"rsds"``.

        Returns:
            list[str]: Sorted local paths covering the configured period. Shorter
            than the period if the server was missing months.
        """
        fetch_hostrada(self.cfg, variable)
        return self._find_files(variable)

    def _find_files(self, variable: str) -> List[str]:
        """List the local files for a variable across the configured period.

        Globs one month at a time so the result is in chronological order, which
        matters because :meth:`load` concatenates along ``time`` positionally.

        Args:
            variable (str): CF variable name.

        Returns:
            list[str]: Existing local paths, in month order.
        """
        import glob

        provider = self.cfg.dataset.upper()
        var_dir = os.path.join(self.cfg.data_dir, provider.upper(), variable.upper())

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
        """Record the region of interest for :meth:`load` to apply.

        Nothing is read here. The three modes are mutually exclusive; the first
        supplied wins, in the order ``point``, ``box``, ``shapefile``. Calling
        with no geometry clears any previous request, so the next :meth:`load`
        returns the full German domain.

        A ``shapefile`` is reduced to its bounding box, not clipped to the
        polygon — cells outside the geometry but inside its envelope are kept.

        Args:
            point (tuple[float, float], optional): ``(lon, lat)`` in degrees,
                longitude first. Selects the single nearest grid cell.
            box (dict, optional): Bounding box with keys ``lat_min``, ``lat_max``,
                ``lon_min``, ``lon_max``.
            shapefile (str | geopandas.GeoDataFrame, optional): Geometry whose
                bounding box is used.
            buffer_km (float): Accepted for signature compatibility with the
                other providers; ignored.

        Returns:
            HOSTRADAmirror: ``self``, so the call can be chained into :meth:`load`.
        """
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
        """Download what is missing, then open every monthly file as one dataset.

        The ``Y``/``X`` index rectangle is computed once from the first file and
        reused for all of them, so the per-file ``preprocess`` hook is a cheap
        ``isel`` rather than a repeated mask evaluation. Files are opened
        serially (``parallel=False``): HOSTRADA's NetCDF-4 files are not reliably
        thread-safe to open concurrently. Finally the result is clipped to the
        configured time range.

        Args:
            variable (str): CF variable name, e.g. ``"rsds"``.
            chunking (dict, optional): Chunk sizes applied after opening, e.g.
                ``{"time": "auto"}``. Only used when ``use_dask`` is true.
            use_dask (bool): Whether ``chunking`` is applied. Defaults to ``True``.

        Returns:
            xr.Dataset: The concatenated, subset dataset, also stored on
            :attr:`dataset`. Dimensions are ``(time, Y, X)``.

        Raises:
            FileNotFoundError: If no files exist for the variable and period.
            ValueError: If the requested box does not intersect the German domain.
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
        """Reduce the requested geometry to integer ``Y``/``X`` index slices.

        Reads the 2-D ``lat``/``lon`` auxiliary coordinates from one file and
        masks them against the request. For ``point`` this is the cell minimising
        squared distance in degrees — an approximation, since a degree of
        longitude is shorter than one of latitude, but harmless at HOSTRADA's
        1 km spacing. For ``box`` and ``shapefile`` it is the bounding rectangle
        of every matching cell, so the result is a superset of the request.

        Args:
            sample_file (str): Path to one HOSTRADA file, used for its grid.
            mode (str | None): ``"point"``, ``"box"``, ``"shapefile"``, or
                ``None`` for no subsetting.

        Returns:
            tuple[slice, slice] | tuple[None, None]: Index slices for ``Y`` and
            ``X``, or ``(None, None)`` when ``mode`` is ``None``.

        Raises:
            ValueError: If the requested area does not intersect the grid.
        """
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
        """Clip a dataset to the configured time range.

        Monthly files always overrun the requested range at both ends. A failure
        here is reported and swallowed, returning the unclipped dataset, so a
        surprising time encoding costs extra timesteps rather than the whole load.

        Args:
            ds (xr.Dataset): Dataset with a ``time`` dimension.

        Returns:
            xr.Dataset: The clipped dataset, or ``ds`` unchanged on failure.
        """
        start = getattr(self.cfg.time_range, "start_date", None)
        end = getattr(self.cfg.time_range, "end_date", None)
        if start or end:
            try:
                ds = ds.sel(time=slice(start, end))
            except Exception as exc:
                print(f"⚠️  Time subsetting failed: {exc}")
        return ds
