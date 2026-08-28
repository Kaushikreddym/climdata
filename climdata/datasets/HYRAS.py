"""HYRAS 1 km gridded observations for Germany, from the DWD Open Data server.

HYRAS is published in two shapes that this module has to reconcile. The older
variables (``pr``, ``tas``, ``tasmin``, ``tasmax``, ``hurs``) ship as one NetCDF
per year on a curvilinear grid: projected ``x``/``y`` dimensions in Gauss-Krueger
zone 3 (EPSG:31467) with 2-D ``lat(y, x)`` / ``lon(y, x)`` coordinates. The newer
soil and evaporation variables (``evpot``, ``soilMoist``, ``soilTemp``) ship as
monthly ``.tgz`` archives of daily ESRI ASCII rasters, which are read straight
out of the archive without ever being extracted to disk.

Because the geographic coordinates are two-dimensional, a lat/lon request cannot
be a coordinate slice; it is resolved to integer indices through a k-d tree
(:func:`find_nearest_xy`) or a rasterised mask, computed once from a sample file
and reused for every file in the request.

Example:
    >>> hyras = HYRASmirror(cfg)                                  # doctest: +SKIP
    >>> hyras.extract(point=(13.4, 52.5))                         # doctest: +SKIP
    >>> ds = hyras.load("pr", chunking={"time": "auto"})          # doctest: +SKIP
"""

import os
import tarfile
from datetime import datetime
from typing import Dict, Optional, Tuple, List

import numpy as np
import xarray as xr
import geopandas as gpd
import pandas as pd
import requests
from scipy.spatial import cKDTree
import rasterio.features as rfeatures
from rasterio.transform import from_bounds
from shapely.geometry import mapping

from omegaconf import DictConfig

def fetch_dwd(var_cfg, var):
    """Download the HYRAS files covering the configured period for one variable.

    Picks the layout from the variable's ``prefix``: monthly ``.tgz`` archives
    for the soil and evaporation variables, yearly ``.nc`` files for the rest.
    Files already on disk are skipped, so an interrupted run resumes.

    The two layouts differ in how they treat a missing file, deliberately. The
    ``.tgz`` variables are still being published, so a month the server does not
    yet have is reported and skipped. The ``.nc`` variables are a closed
    archive, so a missing year means the request is wrong and it raises.

    Args:
        var_cfg (DictConfig): Configuration with ``dataset``, ``data_dir``,
            ``time_range`` and ``dsinfo.HYRAS.variables``.
        var (str): CF variable name, e.g. ``"pr"``.

    Returns:
        None: Files are written under ``<data_dir>/HYRAS/<VARIABLE>/``.

    Raises:
        FileNotFoundError: If a required yearly ``.nc`` file is absent from the server.
        RuntimeError: If a yearly ``.nc`` download fails partway through.
        KeyError: If ``var`` is not declared in ``dsinfo.HYRAS.variables``.
    """
    param_mapping = var_cfg.dsinfo
    provider = var_cfg.dataset.upper()
    parameter_key = var
    # Validate provider and parameter

    param_info = param_mapping[provider]['variables'][parameter_key]
    base_url = param_info["base_url"]
    prefix = param_info["prefix"]
    version = param_info["version"]

    start_date = var_cfg.time_range.start_date
    end_date = var_cfg.time_range.end_date

    # Parse dates & extract unique years and months
    start_year = datetime.fromisoformat(start_date).year
    start_month = datetime.fromisoformat(start_date).month
    end_year = datetime.fromisoformat(end_date).year
    end_month = datetime.fromisoformat(end_date).month

    # output_file = cfg.output.filename
    # os.makedirs(parameter_key, exist_ok=True)

    # Determine if this is a newer dataset (uses .tgz with YYYYMM format) or older (uses .nc with YYYY_version_de.nc)
    is_tgz_format = prefix in ["grids_germany_daily_evapo_p", "grids_germany_daily_soil_moist", "grids_germany_daily_soil_temperature_5cm"]

    if is_tgz_format:
        # Handle new format: grids_germany_daily_*.tgz with year-month format
        for year in range(start_year, end_year + 1):
            start_m = start_month if year == start_year else 1
            end_m = end_month if year == end_year else 12
            
            for month in range(start_m, end_m + 1):
                date_str = f"{year}{month:02d}"
                file_name = f"{prefix}_{date_str}.tgz"
                file_url = f"{base_url}{file_name}"
                local_path = os.path.join(var_cfg.data_dir, provider, parameter_key.upper(), file_name)
                os.makedirs(os.path.dirname(local_path), exist_ok=True)
                print(f"⬇️  Checking: {file_url}")

                if os.path.exists(local_path):
                    print(f"✔️  Exists locally: {local_path}")
                    continue

                # Check if file exists on server first (HEAD request)
                head = requests.head(file_url)
                if head.status_code != 200:
                    print(f"⚠️  Not found on server: {file_url} (HTTP {head.status_code})")
                    continue

                print(f"⬇️  Downloading: {file_url}")
                try:
                    response = requests.get(file_url, stream=True)
                    response.raise_for_status()
                    with open(local_path, "wb") as f:
                        for chunk in response.iter_content(chunk_size=8192):
                            f.write(chunk)
                    print(f"✅ Saved: {local_path}")
                except requests.HTTPError as e:
                    print(f"⚠️  Failed download: {file_url} — {e}")
    else:
        # Handle old format: *_{year}_{version}_de.nc or *_{year}_{version}.nc (for ET0/evaporation_fao)
        # ET0 files don't have '_de' suffix
        has_de_suffix = prefix != "grids_germany_daily_evaporation_fao"
        
        for year in range(start_year, end_year + 1):
            # Construct filename based on whether it has _de suffix
            if has_de_suffix:
                file_name = f"{prefix}_{year}_{version}_de.nc"
            else:
                file_name = f"{prefix}_{year}_{version}.nc"
            
            file_url = f"{base_url}{file_name}"
            local_path = os.path.join(var_cfg.data_dir,provider,parameter_key.upper(), file_name)
            os.makedirs(os.path.dirname(local_path), exist_ok=True)
            print(f"⬇️  Checking: {file_url}")

            if os.path.exists(local_path):
                print(f"✔️  Exists locally: {local_path}")
                continue

            # Check if file exists on server first (HEAD request)
            head = requests.head(file_url)
            if head.status_code != 200:
                raise FileNotFoundError(f"❌ Not found on server: {file_url} (HTTP {head.status_code})")

            print(f"⬇️  Downloading: {file_url}")
            try:
                response = requests.get(file_url, stream=True)
                response.raise_for_status()
                with open(local_path, "wb") as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)
                print(f"✅ Saved: {local_path}")
            except requests.HTTPError as e:
                raise RuntimeError(f"❌ Failed download: {file_url} — {e}")

def read_asc_file(asc_file: str, varname: str = None, units: str = None) -> xr.DataArray:
    """Read one ESRI ASCII raster into a DataArray on the Gauss-Krueger grid.

    The six-line header (``ncols``, ``nrows``, ``xllcorner``, ``yllcorner``,
    ``cellsize``, ``NODATA_value``) is parsed, the nodata sentinel replaced with
    NaN, and the grid given descending ``y`` coordinates to match the
    north-to-south row order the format stores.

    DWD publishes ``evpot`` and ``soilTemp`` scaled by ten; those are divided
    back, while ``soilMoist`` is left alone. Getting ``varname`` wrong therefore
    yields values off by a factor of ten rather than an error.

    Args:
        asc_file (str | os.PathLike | IO): Path to a ``.asc`` file, or an open
            binary/text file object — which is how members are read directly out
            of a ``.tgz`` archive.
        varname (str, optional): CF variable name, used only to decide scaling.
        units (str, optional): Units string to record in the attributes.

    Returns:
        xr.DataArray: Dimensions ``(y, x)`` in EPSG:31467 metres, with
        ``xllcorner``, ``yllcorner``, ``cellsize`` and ``crs_grid`` attributes.
    """
    import io
    if isinstance(asc_file, (str, os.PathLike)):
        ctx = open(asc_file, 'r')
    else:
        # file-like object (e.g. from tarfile.extractfile) — wrap bytes as text
        raw = asc_file.read()
        ctx = io.StringIO(raw.decode('latin-1') if isinstance(raw, bytes) else raw)
    with ctx as f:
        header = {}
        for _ in range(6):
            line = f.readline().strip().split()
            header[line[0].lower()] = float(line[1])
        
        ncols = int(header['ncols'])
        nrows = int(header['nrows'])
        xllcorner = header['xllcorner']
        yllcorner = header['yllcorner']
        cellsize = header['cellsize']
        nodata = header['nodata_value']
        
        # Read data
        # loadtxt squeezes a single-row or single-column grid to 1-D, which does
        # not match the ('y', 'x') dims below. The header already declares the
        # shape, so reshape to it — that also catches a truncated file.
        data = np.loadtxt(f).reshape(nrows, ncols)

        # Handle NODATA values
        data[data == nodata] = np.nan
        
        # Conditional scaling based on variable type
        # Only evpot and soilTemp are divided by 10; soilMoist is not scaled
        if varname in ['evpot', 'soilTemp']:
            data = data / 10.0
        
        # Create coordinate arrays (Gauss Krüger 3 grid coordinates)
        # Note: ASCII raster data is stored top-to-bottom (north to south)
        # so we need y coordinates in descending order to match data order
        x = xllcorner + np.arange(ncols) * cellsize
        y = yllcorner + (np.arange(nrows))[::-1] * cellsize  # Reverse to match data top-to-bottom order
        
        # Create DataArray with grid attributes
        da = xr.DataArray(
            data,
            dims=('y', 'x'),
            coords={'y': y, 'x': x}
        )
        
        # Store grid reference information for coordinate transformation
        da.attrs['xllcorner'] = xllcorner
        da.attrs['yllcorner'] = yllcorner
        da.attrs['cellsize'] = cellsize
        da.attrs['crs_grid'] = 'EPSG:31467'  # Gauss Krüger 3 (standard for German grids)
        if units:
            da.attrs['units'] = units
        
        return da

def read_asc_timeseries(asc_files: List[str], varname: str = 'value', units: str = None) -> xr.Dataset:
    """Stack per-day ASCII rasters into one dataset with geographic coordinates.

    Files are sorted by name, so the ``YYYYMMDD`` stamp they carry puts them in
    chronological order. The projected ``x``/``y`` grid is transformed to 2-D
    ``lat``/``lon`` auxiliary coordinates, matching the layout of the NetCDF
    HYRAS variables so both branches of :meth:`HYRASmirror.load` return the same
    shape.

    Unlike :func:`read_asc_timeseries_from_tgz`, every file is read eagerly.

    Args:
        asc_files (list[str]): Paths to ``.asc`` files, one per day.
        varname (str): Name for the resulting data variable. Defaults to ``"value"``.
        units (str, optional): Units string to record on the variable.

    Returns:
        xr.Dataset: Dimensions ``(time, y, x)``, with ``lat``/``lon`` coordinates
        and ``crs_grid`` / ``crs_geographic`` attributes.
    """
    import re
    
    data_arrays = []
    times = []
    
    # Sort files to ensure chronological order
    asc_files = sorted(asc_files)
    
    for asc_file in asc_files:
        # Try to extract date from filename (format: prefix_YYYYMMDD.asc)
        basename = os.path.basename(asc_file)
        # Extract YYYYMMDD from filename
        match = re.search(r'(\d{8})', basename)
        if match:
            date_str = match.group(1)
            time = pd.to_datetime(date_str, format='%Y%m%d')
        else:
            time = pd.to_datetime(len(data_arrays), unit='D')
        
        times.append(time)
        da = read_asc_file(asc_file, varname=varname, units=units)
        data_arrays.append(da.values)
    
    # Stack into 3D array (time, y, x)
    data_3d = np.stack(data_arrays, axis=0)
    
    # Get coordinates from first file
    da_first = read_asc_file(asc_files[0], varname=varname, units=units)
    y_grid = da_first.coords['y'].values  # Grid Y (Gauss Krüger 3)
    x_grid = da_first.coords['x'].values  # Grid X (Gauss Krüger 3)
    
    # Transform from Gauss Krüger 3 (EPSG:31467) to WGS84 (EPSG:4326) lat/lon
    try:
        from pyproj import Transformer
        # Create transformer from Gauss Krüger 3 to WGS84
        transformer = Transformer.from_crs(31467, 4326, always_xy=True)
        
        # Create meshgrid of Gauss Krüger coordinates
        x_grid_2d, y_grid_2d = np.meshgrid(x_grid, y_grid)
        
        # Transform all points at once
        lon_2d, lat_2d = transformer.transform(x_grid_2d, y_grid_2d)
    except ImportError:
        # Fallback: use x/y directly as lon/lat (not ideal but works for visualization)
        lon_2d, lat_2d = np.meshgrid(x_grid, y_grid)
    
    # Create dataset with grid coordinates (x/y), geographic coordinates (lat/lon), and time
    ds = xr.Dataset(
        {varname: (['time', 'y', 'x'], data_3d)},
        coords={
            'time': times, 
            'y': y_grid, 
            'x': x_grid,
            'lat': (['y', 'x'], lat_2d),
            'lon': (['y', 'x'], lon_2d),
        }
    )
    
    # Add units to the data variable (critical for xclim and unit conversions)
    if units:
        ds[varname].attrs['units'] = units
    
    # Add CRS information at dataset level
    ds.attrs['crs_grid'] = 'EPSG:31467'  # Grid CRS: Gauss Krüger 3
    ds.attrs['crs_geographic'] = 'EPSG:4326'  # Geographic CRS: WGS84
    ds.attrs['description'] = 'HYRAS 1km gridded climate data (Germany)'
    if units:
        ds.attrs['units'] = units
    
    # Add coordinate reference system using rioxarray if available
    try:
        import rasterio.crs as rio_crs
        crs = rio_crs.CRS.from_epsg(4326)
        ds.rio.write_crs(crs, inplace=True)
    except (ImportError, AttributeError):
        # Fallback if rasterio/rioxarray not available
        pass
    
    return ds

def _iter_tgz_members(tgz_path: str):
    """Iterate the raster members of a ``.tgz`` archive without extracting it.

    Args:
        tgz_path (str): Path to the archive.

    Yields:
        tuple[str | None, str, io.BytesIO]: ``(date_str, member_name, file_like)``
        for each ``.asc`` or ``.nc`` member, in name order. ``date_str`` is the
        ``YYYYMMDD`` stamp found in the member name, or ``None`` if it has none.
        The file object holds the member in memory.
    """
    import io, re
    with tarfile.open(tgz_path, 'r:gz') as tar:
        for member in sorted(tar.getmembers(), key=lambda m: m.name):
            name = os.path.basename(member.name)
            if not (name.endswith('.asc') or name.endswith('.nc')):
                continue
            m = re.search(r'(\d{8})', name)
            date_str = m.group(1) if m else None
            raw = tar.extractfile(member)
            if raw is None:
                continue
            yield date_str, name, io.BytesIO(raw.read())


def _read_one_asc_member(tgz_path: str, member_name: str,
                          varname: str, units: str) -> np.ndarray:
    """Read one ``.asc`` member out of an archive as a plain array.

    Reopens the archive on every call rather than holding a handle, because this
    runs inside a ``dask.delayed`` task: a shared ``tarfile`` object is neither
    picklable nor safe to read from several workers at once.

    Args:
        tgz_path (str): Path to the archive.
        member_name (str): Name of the member to read.
        varname (str): CF variable name, passed to :func:`read_asc_file` for scaling.
        units (str): Units string, passed through to :func:`read_asc_file`.

    Returns:
        numpy.ndarray: The 2-D grid as ``float64``.

    Raises:
        RuntimeError: If the member cannot be extracted.
    """
    import io, tarfile as _tarfile
    with _tarfile.open(tgz_path, 'r:gz') as tar:
        raw = tar.extractfile(member_name)
        if raw is None:
            raise RuntimeError(f"Cannot extract {member_name} from {tgz_path}")
        fileobj = io.BytesIO(raw.read())
    da = read_asc_file(fileobj, varname=varname, units=units)
    return da.values.astype(np.float64)


def read_asc_timeseries_from_tgz(tgz_files: List[str], varname: str = 'value',
                                  units: str = None) -> xr.Dataset:
    """Build a lazy dataset over the daily rasters inside several ``.tgz`` archives.

    Scans every archive once to index its members and read the grid geometry from
    the first, then wraps each daily raster in a ``dask.delayed`` call. Nothing
    but that first raster is decompressed until the result is computed, so a
    multi-year request costs an index scan rather than a full decompression.

    Args:
        tgz_files (list[str]): Paths to monthly archives. Sorted internally.
        varname (str): Name for the resulting data variable. Defaults to ``"value"``.
        units (str, optional): Units string to record on the variable.

    Returns:
        xr.Dataset: Dask-backed, dimensions ``(time, y, x)``, with ``lat``/``lon``
        coordinates and the same attributes as :func:`read_asc_timeseries`.
    """
    import dask
    import dask.array as da

    # ------------------------------------------------------------------ #
    # 1) Scan all archives to collect (tgz_path, member_name, date) tuples
    #    and determine the 2-D grid shape from the very first .asc member.
    # ------------------------------------------------------------------ #
    index: List[tuple] = []          # (tgz_path, member_name, pd.Timestamp)
    first_y = first_x = None
    grid_shape = None

    for tgz_path in sorted(tgz_files):
        for date_str, name, fileobj in _iter_tgz_members(tgz_path):
            if not name.endswith('.asc'):
                continue
            time_stamp = pd.to_datetime(date_str, format='%Y%m%d') if date_str else None
            if time_stamp is None:
                continue
            if grid_shape is None:
                # Read the very first slice eagerly to get coords + shape
                ref_da = read_asc_file(fileobj, varname=varname, units=units)
                first_y = ref_da.coords['y'].values
                first_x = ref_da.coords['x'].values
                grid_shape = ref_da.shape          # (ny, nx)
            # Store the *member path inside the archive* (str) for the worker
            index.append((tgz_path, name, time_stamp))

    if not index:
        raise FileNotFoundError(f"No .asc members found in: {tgz_files}")

    # ------------------------------------------------------------------ #
    # 2) Build one dask.delayed per time step and stack into a dask array
    # ------------------------------------------------------------------ #
    ny, nx = grid_shape
    delayed_slices = [
        dask.delayed(_read_one_asc_member)(tgz, member, varname, units)
        for tgz, member, _ in index
    ]
    dask_slices = [
        da.from_delayed(s, shape=(ny, nx), dtype=np.float64)
        for s in delayed_slices
    ]
    data_dask = da.stack(dask_slices, axis=0)          # (time, y, x)
    times = [t for _, _, t in index]

    # ------------------------------------------------------------------ #
    # 3) Build lat/lon auxiliary coordinates (tiny eager computation)
    # ------------------------------------------------------------------ #
    try:
        from pyproj import Transformer
        transformer = Transformer.from_crs(31467, 4326, always_xy=True)
        x_2d, y_2d = np.meshgrid(first_x, first_y)
        lon_2d, lat_2d = transformer.transform(x_2d, y_2d)
    except ImportError:
        lon_2d, lat_2d = np.meshgrid(first_x, first_y)

    # ------------------------------------------------------------------ #
    # 4) Wrap in xarray Dataset
    # ------------------------------------------------------------------ #
    ds = xr.Dataset(
        {varname: xr.Variable(['time', 'y', 'x'], data_dask)},
        coords={
            'time': times,
            'y': first_y,
            'x': first_x,
            'lat': (['y', 'x'], lat_2d),
            'lon': (['y', 'x'], lon_2d),
        }
    )
    if units:
        ds[varname].attrs['units'] = units
    ds.attrs['crs_grid'] = 'EPSG:31467'
    ds.attrs['crs_geographic'] = 'EPSG:4326'
    ds.attrs['description'] = 'HYRAS 1km gridded climate data (Germany)'
    return ds


def find_nearest_xy(ds, target_lat, target_lon):
    """Find the grid indices nearest a geographic point on a curvilinear grid.

    Builds a k-d tree over the flattened 2-D ``lat``/``lon`` coordinates and
    queries it in degree space — treating latitude and longitude as if they were
    Cartesian. That overweights longitude at German latitudes, but at HYRAS's
    1 km spacing the answer is the same cell.

    The tree is rebuilt on every call, so this belongs outside per-file loops.

    Args:
        ds (xr.Dataset): Dataset with 2-D ``lat`` and ``lon`` coordinates.
        target_lat (float): Latitude in degrees.
        target_lon (float): Longitude in degrees.

    Returns:
        tuple[int, int]: ``(iy, ix)`` indices into the ``y`` and ``x`` dimensions.
    """
    lat = ds['lat'].values  # shape (y,x) or (x,y)
    lon = ds['lon'].values

    # Flatten to 1D for k-d tree
    lat_flat = lat.flatten()
    lon_flat = lon.flatten()

    tree = cKDTree(np.column_stack((lat_flat, lon_flat)))
    _, idx = tree.query([target_lat, target_lon])
    iy, ix = np.unravel_index(idx, lat.shape)

    return iy, ix


class HYRASmirror:
    """Download, open and subset HYRAS grids across both of its file layouts.

    Call :meth:`extract` to record the region of interest, then :meth:`load`,
    which downloads what is missing and applies the subset. Each of the three
    subsetting modes takes a different route, chosen for cost:

    * **point** — resolved per file inside ``open_mfdataset``'s ``preprocess``
      hook, so only one cell of each file is ever materialised.
    * **box** — index bounds computed once from a sample file, then applied to
      each file opened individually and concatenated. Files are opened one at a
      time rather than through ``open_mfdataset`` because the index bounds must
      be known before the first file is read.
    * **shapefile** — a mask rasterised once, then applied with ``where``. Unlike
      the box path this keeps the full grid extent and blanks cells outside the
      geometry, so the result has the source's spatial shape.

    Loading also normalises two unit strings DWD writes in a form xclim does not
    accept: ``pr`` in ``"mm"`` becomes ``"mm/day"``, and ``hurs`` in
    ``"Percent"`` becomes ``"%"``.

    Attributes:
        cfg (DictConfig): Hydra configuration.
        variables (list[str]): ``cfg.variables``.
        files (list[str]): Local paths resolved by the last :meth:`fetch`.
        dataset (xr.Dataset | None): The loaded dataset, set by :meth:`load`.

    Example:
        >>> hyras = HYRASmirror(cfg)                                 # doctest: +SKIP
        >>> hyras.extract(box={"lat_min": 47, "lat_max": 55,         # doctest: +SKIP
        ...                    "lon_min": 6, "lon_max": 15})
        >>> ds = hyras.load("pr")                                    # doctest: +SKIP
    """

    def __init__(self, cfg: DictConfig):
        """Bind a configuration.

        Args:
            cfg (DictConfig): Configuration with ``dataset``, ``data_dir``,
                ``variables``, ``time_range`` and ``dsinfo.HYRAS``.
        """
        self.cfg = cfg
        self.dataset: Optional[xr.Dataset] = None
        self.variables = cfg.variables
        self.files: List[str] = []

        # extraction state
        self._extract_mode: Optional[str] = None
        self._extract_params = None

        # cached grid helpers (computed from first file when needed)
        self._cached_box_idx = None  # (y0,y1,x0,x1)
        self._cached_mask = None     # 2D boolean mask for shapefile
        self._cached_lonlat_info = None  # (lon_min, lon_max, lat_min, lat_max, nx, ny)

    # --------------------------
    # File discovery / fetch
    # --------------------------
    def fetch(self, variable: str) -> List[str]:
        """Download a variable's files and list the local paths now available.

        Args:
            variable (str): CF variable name, e.g. ``"pr"``.

        Returns:
            list[str]: Local paths in chronological order — monthly ``.tgz``
            archives or yearly ``.nc`` files depending on the variable. Also
            stored on :attr:`files`.
        """
        # keep your fetch behavior (calls fetch_dwd)
        fetch_dwd(self.cfg, variable)

        provider = self.cfg.dataset.upper()
        param_info = self.cfg.dsinfo[provider]['variables'][variable]
        prefix = param_info["prefix"]
        version = param_info["version"]

        start_date = datetime.fromisoformat(self.cfg.time_range.start_date)
        end_date = datetime.fromisoformat(self.cfg.time_range.end_date)
        start_year = start_date.year
        start_month = start_date.month
        end_year = end_date.year
        end_month = end_date.month

        files = []
        var_dir = os.path.join(self.cfg.data_dir, provider.upper(), variable.upper())
        
        # Determine if this is a newer dataset (uses .tgz/.asc with YYYYMM naming)
        is_new_format = prefix in ["grids_germany_daily_evapo_p", "grids_germany_daily_soil_moist", "grids_germany_daily_soil_temperature_5cm"]

        if is_new_format:
            # Return .tgz paths directly — members are read in-memory at load time
            for year in range(start_year, end_year + 1):
                start_m = start_month if year == start_year else 1
                end_m = end_month if year == end_year else 12
                
                for month in range(start_m, end_m + 1):
                    tgz_name = f"{prefix}_{year}{month:02d}.tgz"
                    tgz_path = os.path.join(var_dir, tgz_name)
                    if os.path.exists(tgz_path):
                        files.append(tgz_path)
        else:
            # Handle old format: *_{year}_{version}_de.nc
            for year in range(start_year, end_year + 1):
                file_name = f"{prefix}_{year}_{version}_de.nc"
                file_path = os.path.join(var_dir, file_name)
                if os.path.exists(file_path):
                    files.append(file_path)

        self.files = files
        
        return files

    # --------------------------
    # Public extract setter
    # --------------------------
    def extract(self, *, point: Tuple[float, float] = None, box: Dict[str, float] = None,
                shapefile: str = None, buffer_km: float = 0.0):
        """Record the region of interest for :meth:`load` to apply.

        Nothing is read here, and any cached grid indices from a previous call
        are discarded. The three modes are mutually exclusive; the first supplied
        wins, in the order ``point``, ``box``, ``shapefile``.

        Args:
            point (tuple[float, float], optional): ``(lon, lat)`` in degrees,
                longitude first. Selects the single nearest grid cell.
            box (dict, optional): Bounding box; all four of ``lat_min``,
                ``lat_max``, ``lon_min``, ``lon_max`` are required.
            shapefile (str | geopandas.GeoDataFrame, optional): Polygon(s) to
                mask to, as a path or an in-memory frame.
            buffer_km (float): Dilation applied to ``shapefile`` geometries in
                Web Mercator. Ignored for ``point`` and ``box``. Defaults to ``0.0``.

        Returns:
            HYRASmirror: ``self``, so the call can be chained into :meth:`load`.

        Raises:
            ValueError: If no geometry was given, or if ``box`` is missing a key.
        """
        if point is not None:
            lon, lat = point
            self._extract_mode = "point"
            self._extract_params = (lon, lat)

        elif box is not None:
            # expect keys lat_min, lat_max, lon_min, lon_max
            for k in ("lat_min", "lat_max", "lon_min", "lon_max"):
                if k not in box:
                    raise ValueError(f"box missing key {k}")
            self._extract_mode = "box"
            self._extract_params = box

        elif shapefile is not None:
            gdf = gpd.read_file(shapefile) if isinstance(shapefile, str) else shapefile
            if buffer_km > 0:
                gdf = gdf.to_crs(epsg=3857)
                gdf["geometry"] = gdf.buffer(buffer_km * 1000)
                gdf = gdf.to_crs(epsg=4326)
            self._extract_mode = "shapefile"
            self._extract_params = gdf

        else:
            raise ValueError("Must provide point, box, or shapefile.")

        # Clear cached helpers when extraction changes
        self._cached_box_idx = None
        self._cached_mask = None
        self._cached_lonlat_info = None

        return self

    # --------------------------
    # Helpers to compute indices/mask from a sample file
    # --------------------------
    def _load_sample_grid(self, sample_file: str, varname: Optional[str] = None):
        """Read just the coordinate arrays from one file.

        Opens with ``decode_times=False`` and touches only the coordinates, so
        no data array is materialised.

        Args:
            sample_file (str): Path to a HYRAS NetCDF file.
            varname (str, optional): Unused; accepted for call-site symmetry.

        Returns:
            tuple[numpy.ndarray, numpy.ndarray]: ``(lat, lon)``, each 1-D or 2-D
            depending on the file.
        """
        ds = xr.open_dataset(sample_file, engine="netcdf4", decode_times=False)
        # Try to access coordinates in common names; adapt if your files differ.
        # Accept either 1D 'lat','lon' or 2D 'lat','lon' on dims (y,x).
        if ("lat" in ds.coords) and ("lon" in ds.coords):
            lat = ds["lat"].values
            lon = ds["lon"].values
        else:
            # fallback: if coordinates stored as variables
            lat = ds["lat"].values
            lon = ds["lon"].values

        ds.close()
        return lat, lon

    def _compute_box_indices(self, sample_file: str):
        """Resolve the requested box to cached integer index bounds.

        Maps the box's two opposite corners to grid indices through
        :func:`find_nearest_xy` and takes the enclosing rectangle. Because the
        grid is curvilinear, that rectangle is an approximation of the
        geographic box, and a rotated domain yields slightly more area than
        asked for. The result is cached until :meth:`extract` is called again.

        Args:
            sample_file (str): Path to a file whose grid represents the request.

        Returns:
            tuple[int, int, int, int]: ``(y0, y1, x0, x1)`` as Python slice
            endpoints, so ``y1`` and ``x1`` are exclusive.
        """
        if self._cached_box_idx is not None:
            return self._cached_box_idx

        lat, lon = self._load_sample_grid(sample_file)
        box = self._extract_params
        # find_nearest_xy expects ds-like input; we'll open a tiny ds for indices
        ds_sample = xr.open_dataset(sample_file, engine="netcdf4", decode_times=False)
        iy_min, ix_min = find_nearest_xy(ds_sample, box["lat_min"], box["lon_min"])
        iy_max, ix_max = find_nearest_xy(ds_sample, box["lat_max"], box["lon_max"])
        ds_sample.close()

        y0, y1 = sorted([iy_min, iy_max])
        x0, x1 = sorted([ix_min, ix_max])
        self._cached_box_idx = (y0, y1 + 1, x0, x1 + 1)  # python slice endpoints
        return self._cached_box_idx

    def _compute_shapefile_mask(self, sample_file: str):
        """Rasterise the requested geometry onto the grid, and cache the mask.

        The transform is derived from the coordinate bounding box and the grid
        shape, which treats the curvilinear grid as regular in lat/lon. Over
        Germany the distortion is under a cell; over a larger or more rotated
        domain it would not be.

        Args:
            sample_file (str): Path to a file whose grid represents the request.

        Returns:
            numpy.ndarray: Boolean mask of shape ``(ny, nx)``, ``True`` inside
            the geometry.

        Raises:
            RuntimeError: If the extraction mode is not ``"shapefile"``, or the
                coordinate arrays are neither 1-D nor 2-D.
        """
        if self._cached_mask is not None:
            return self._cached_mask

        if self._extract_mode != "shapefile":
            raise RuntimeError("shapefile mask requested but extract mode is not 'shapefile'")

        gdf = self._extract_params
        lat, lon = self._load_sample_grid(sample_file)

        # handle 1D or 2D lat/lon
        if lat.ndim == 1 and lon.ndim == 1:
            ny = lat.size
            nx = lon.size
            lon_min, lon_max = float(lon.min()), float(lon.max())
            lat_min, lat_max = float(lat.min()), float(lat.max())
        elif lat.ndim == 2 and lon.ndim == 2:
            ny, nx = lat.shape
            lon_min, lon_max = float(lon.min()), float(lon.max())
            lat_min, lat_max = float(lat.min()), float(lat.max())
        else:
            raise RuntimeError("Unsupported lat/lon shapes for rasterization")

        transform = from_bounds(lon_min, lat_min, lon_max, lat_max, nx, ny)

        shapes = ((mapping(geom), 1) for geom in gdf.geometry)

        mask = rfeatures.rasterize(
            shapes=shapes,
            out_shape=(ny, nx),
            transform=transform,
            fill=0,
            default_value=1,
            dtype="uint8",
        ).astype(bool)

        self._cached_mask = mask
        # also cache lon/lat info (useful if needed)
        self._cached_lonlat_info = (lon_min, lon_max, lat_min, lat_max, nx, ny)
        return mask

    # --------------------------
    # Core loading logic
    # --------------------------

    def _apply_time_subset(self, ds):
        """Clip a dataset to the configured time range.

        Yearly and monthly files always overrun the request at both ends. If the
        time coordinate has not been decoded, ``decode_cf`` is applied and the
        slice retried; if that also fails the dataset is returned unclipped, so a
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
            except Exception:
                # If time coord not decoded yet, try forcing decode then slice
                try:
                    ds = xr.decode_cf(ds)
                    ds = ds.sel(time=slice(start, end))
                except Exception:
                    # fallback: return ds unchanged
                    pass
        return ds
    def load(self, variable: str, use_dask: bool = True, chunking: dict = None):
        """Download what is missing, then open the variable with the subset applied.

        Dispatches on file layout first (``.tgz`` archives versus yearly NetCDF)
        and then on extraction mode; see the class docstring for why each mode
        takes the route it does. Unit strings are normalised and the result is
        clipped to the configured time range.

        Args:
            variable (str): CF variable name, e.g. ``"pr"``.
            use_dask (bool): Whether ``chunking`` is applied. Defaults to ``True``.
            chunking (dict, optional): Chunk sizes, e.g. ``{"time": "auto"}``.
                Only used when ``use_dask`` is true.

        Returns:
            xr.Dataset: The subset dataset, also stored on :attr:`dataset`.
            Dimensions are ``(time, y, x)``, or ``(time,)`` in point mode.

        Raises:
            FileNotFoundError: If no files exist for the variable and period.
        """

        files = self.fetch(variable)
        if not files:
            raise FileNotFoundError(f"No files found for variable {variable}")

        mode = self._extract_mode

        # Check if files are .tgz archives (read members in-memory)
        is_tgz_format = any(f.endswith('.tgz') for f in files)
        if is_tgz_format:
            tgz_files = [f for f in files if f.endswith('.tgz')]
            provider = self.cfg.dataset.upper()
            units = None
            try:
                param_info = self.cfg.dsinfo[provider]['variables'].get(variable, {})
                units = param_info.get('units')
            except (AttributeError, KeyError, TypeError):
                pass
            dset = read_asc_timeseries_from_tgz(tgz_files, varname=variable, units=units)
            # Apply same spatial extraction as the .asc branch below
            original_attrs = dset.attrs.copy()
            if mode == "point":
                lon_pt, lat_pt = self._extract_params
                ix = np.abs(dset['x'].values - lon_pt).argmin()
                iy = np.abs(dset['y'].values - lat_pt).argmin()
                dset = dset.isel(x=ix, y=iy)
            elif mode == "box":
                box = self._extract_params
                lat_2d = dset['lat'].values
                lon_2d = dset['lon'].values
                mask = (
                    (lat_2d >= box['lat_min']) & (lat_2d <= box['lat_max']) &
                    (lon_2d >= box['lon_min']) & (lon_2d <= box['lon_max'])
                )
                mask_da = xr.DataArray(mask, dims=('y', 'x'),
                                       coords={'y': dset['y'], 'x': dset['x']})
                dset = dset.where(mask_da, drop=True)
            dset.attrs.update(original_attrs)
            start_date = getattr(self.cfg.time_range, "start_date", None)
            end_date = getattr(self.cfg.time_range, "end_date", None)
            if start_date or end_date:
                try:
                    dset = dset.sel(time=slice(start_date, end_date))
                except Exception as e:
                    print(f"⚠️  Time subsetting failed: {e}")
            self.dataset = dset
            return dset

        # -------------------------
        # POINT mode: preprocess per-file and use open_mfdataset
        # -------------------------
        if mode == "point":
            lon, lat = self._extract_params
            def preprocess_point(ds):
                iy, ix = find_nearest_xy(ds, lat, lon)
                # ensure dimension order if present
                if variable in ds:
                    try:
                        ds[variable] = ds[variable].transpose("time", "y", "x")
                        ds["time"] = ds["time"].dt.floor("D")
                    except Exception:
                        pass
                # point selection via nearest index (fast)
                return ds.isel(y=iy, x=ix)

            dset = xr.open_mfdataset(
                files,
                combine="nested",
                concat_dim="time",
                preprocess=preprocess_point,
                engine="netcdf4",
                parallel=False,
            )
            if use_dask and chunking:
                dset = dset.chunk(chunking)

            # normalize pr units
            if "pr" in dset:
                if dset["pr"].attrs.get("units") == "mm":
                    dset["pr"].attrs["units"] = "mm/day"
            if "hurs" in dset:
                if dset["hurs"].attrs.get("units") == "Percent":
                    dset["hurs"].attrs["units"] = "%"

            self.dataset = dset
            return dset

        # -------------------------
        # BOX or SHAPEFILE mode: compute indices/mask once and apply per-file
        # -------------------------
        elif mode in ("box", "shapefile"):
            sample_file = files[0]
            datasets = []

            if mode == "box":
                y0, y1, x0, x1 = self._compute_box_indices(sample_file)
            else:  # shapefile
                mask = self._compute_shapefile_mask(sample_file)

            for f in files:
                # open each file lightly
                ds = xr.open_dataset(f, engine="netcdf4", decode_times=True)

                # ensure dims and variable layout
                if variable in ds:
                    try:
                        ds[variable] = ds[variable].transpose("time", "y", "x")
                        ds["time"] = ds["time"].dt.floor("D")
                    except Exception:
                        pass

                # apply slice or mask
                if mode == "box":
                    sub = ds.isel(y=slice(y0, y1), x=slice(x0, x1))
                else:  # shapefile
                    mask_da = xr.DataArray(mask, dims=("y", "x"))
                    sub = ds.where(mask_da, drop=False)

                # optionally chunk lazily
                if use_dask and chunking:
                    sub = sub.chunk(chunking)

                datasets.append(sub)

            # concatenate along time
            dset = xr.concat(datasets, dim="time", combine_attrs="override")

            # normalize pr units
            if "pr" in dset:
                if dset["pr"].attrs.get("units") == "mm":
                    dset["pr"].attrs["units"] = "mm/day"
            if "hurs" in dset:
                if dset["hurs"].attrs.get("units") == "Percent":
                    dset["hurs"].attrs["units"] = "%"
            dset = self._apply_time_subset(dset)
            self.dataset = dset
            return dset

        else:
            # no extraction mode -> just open normally (light transpose)
            def preprocess_identity(ds):
                if variable in ds:
                    try:
                        ds[variable] = ds[variable].transpose("time", "y", "x")
                        ds["time"] = ds["time"].dt.floor("D")
                    except Exception:
                        pass
                return ds

            dset = xr.open_mfdataset(
                files,
                combine="nested",
                concat_dim="time",
                preprocess=preprocess_identity,
                engine="netcdf4",
                parallel=False,
            )

            if use_dask and chunking:
                dset = dset.chunk(chunking)

            if "pr" in dset:
                if dset["pr"].attrs.get("units") == "mm":
                    dset["pr"].attrs["units"] = "mm/day"

            if "hurs" in dset:
                if dset["hurs"].attrs.get("units") == "Percent":
                    dset["hurs"].attrs["units"] = "%"

            self.dataset = dset
            return dset

    # --------------------------
    # Utility: save current dataset to CSV
    # --------------------------
    def save_csv(self, filename: str, df: pd.DataFrame = None):
        """Write a frame, or the loaded dataset, to CSV.

        A gridded HYRAS dataset flattens to one row per cell per timestep, so
        this can be very large for anything but a point or a small box.

        Args:
            filename (str): Destination path. The parent directory must exist.
            df (pd.DataFrame, optional): Frame to write. Defaults to the loaded
                dataset converted with ``to_dataframe()``.

        Returns:
            None

        Raises:
            ValueError: If no frame is given and no dataset is loaded.
        """
        if df is None:
            if self.dataset is None:
                raise ValueError("No dataset loaded")
            # convert to dataframe (may be large)
            df = self.dataset.to_dataframe().reset_index()
        df.to_csv(filename, index=False)
        print(f"Saved CSV to {filename}")
