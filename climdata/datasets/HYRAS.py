
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

def fetch_dwd(var_cfg,var):
    """Download HYRAS data for one variable and a list of years. Handles both .nc and .tgz formats."""
    param_mapping = var_cfg.dsinfo
    provider = var_cfg.dataset.lower()
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
                    
                    # Extract .tgz file
                    extract_path = os.path.dirname(local_path)
                    print(f"📦 Extracting: {file_name}")
                    with tarfile.open(local_path, 'r:gz') as tar:
                        tar.extractall(path=extract_path)
                    print(f"✅ Extracted: {file_name}")
                except requests.HTTPError as e:
                    print(f"⚠️  Failed download: {file_url} — {e}")
                except tarfile.TarError as e:
                    print(f"⚠️  Failed extraction: {file_name} — {e}")
    else:
        # Handle old format: *_{year}_{version}_de.nc
        for year in range(start_year, end_year + 1):
            file_name = f"{prefix}_{year}_{version}_de.nc"
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
    """
    Read an ASCII raster file (ESRI .asc format) and return as xarray DataArray.
    
    ASCII raster format header:
        ncols         NUMBER
        nrows         NUMBER
        xllcorner     NUMBER
        yllcorner     NUMBER
        cellsize      NUMBER
        NODATA_value  NUMBER
        
    Followed by the data grid.
    
    Parameters:
    -----------
    asc_file : str
        Path to ASCII raster file
    varname : str, optional
        Variable name to determine scaling rules (evpot, soilTemp -> divide by 10; soilMoist -> no scaling)
    units : str, optional
        Units string to include in attributes
    """
    with open(asc_file, 'r') as f:
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
        data = np.loadtxt(f)
        
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
    """
    Read a list of ASCII raster files (one per time step) and stack into a 3D xarray Dataset.
    Transforms grid coordinates (Gauss Krüger 3) to geographic lat/lon and includes CRS information.
    """
    from dateutil import parser
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

def find_nearest_xy(ds, target_lat, target_lon):
    """
    Given a dataset with curvilinear grid, find the nearest x,y index.
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
    """
    Optimized HYRAS mirror loader.

    - Point extraction: done per-file inside preprocess (open_mfdataset).
    - Box / shapefile extraction: done outside open_mfdataset:
        * compute indices / mask from a sample file once
        * apply indices / mask to each file opened individually
        * concat along time
    - Optional dask chunking via use_dask flag.
    """

    def __init__(self, cfg: DictConfig):
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
        """Download HYRAS data for a given variable and time range and return file list. Handles both .nc, .asc, and .tgz formats."""
        # keep your fetch behavior (calls fetch_dwd)
        fetch_dwd(self.cfg, variable)

        provider = self.cfg.dataset.lower()
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
        var_dir = os.path.join(self.cfg.data_dir, provider, variable.upper())
        
        # Determine if this is a newer dataset (uses .tgz/.asc with YYYYMM naming)
        is_new_format = prefix in ["grids_germany_daily_evapo_p", "grids_germany_daily_soil_moist", "grids_germany_daily_soil_temperature_5cm"]

        if is_new_format:
            # Handle new format: .asc files or extracted .nc files from .tgz archives with YYYYMM naming
            import glob
            
            # First, try to find .asc files directly (daily ASCII raster format)
            for year in range(start_year, end_year + 1):
                start_m = start_month if year == start_year else 1
                end_m = end_month if year == end_year else 12
                
                for month in range(start_m, end_m + 1):
                    # Look for .asc files (format: prefix_YYYYMMDD.asc)
                    pattern = f"{prefix}_{year}{month:02d}*.asc"
                    asc_files = sorted(glob.glob(os.path.join(var_dir, pattern)))
                    files.extend(asc_files)
                    
                    # Also look for extracted .nc files from .tgz archives
                    pattern_nc = f"{prefix}_{year}{month:02d}*.nc"
                    nc_files = sorted(glob.glob(os.path.join(var_dir, pattern_nc)))
                    files.extend(nc_files)
            
            # If no direct files found, look in subdirectories (common in .tgz extractions)
            if not files:
                asc_files = sorted(glob.glob(os.path.join(var_dir, "**/*.asc"), recursive=True))
                nc_files = sorted(glob.glob(os.path.join(var_dir, "**/*.nc"), recursive=True))
                files.extend(asc_files)
                files.extend(nc_files)
        else:
            # Handle old format: *_{year}_{version}_de.nc
            import glob
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
        """
        Specify extraction intent.
        - point: (lon, lat)
        - box: dict(lat_min, lat_max, lon_min, lon_max)
        - shapefile: path or GeoDataFrame (if str -> read file)
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
        """
        Open one file (lightweight) and return lat/lon arrays and shape.
        We don't load big data arrays here; just coordinates.
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
        """Compute nearest-array indices for the box on sample grid and cache them."""
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
        """Rasterize shapefile on the sample grid and cache mask (y,x boolean)."""
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

    def compute_point_indices(self, sample_file, lon, lat):
        ds = xr.open_dataset(sample_file)
        x = ds['lon'].values
        y = ds['lat'].values

        ix = np.abs(x - lon).argmin()
        iy = np.abs(y - lat).argmin()
        ds.close()
        return ix, iy
    def _apply_time_subset(self, ds):
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
        """
        Load variable with extraction applied.
        Handles both NetCDF (.nc) and ASCII raster (.asc) files.

        - If extraction mode == 'point' -> uses open_mfdataset with preprocess (fast)
        - If extraction mode in ('box','shapefile') -> open per-file, apply index/mask, concat
        - If files are .asc -> converts to xarray via read_asc_timeseries
        """

        files = self.fetch(variable)
        if not files:
            raise FileNotFoundError(f"No files found for variable {variable}")

        mode = self._extract_mode
        
        # Check if files are ASCII raster format (.asc)
        is_asc_format = any(f.endswith('.asc') for f in files)
        
        # -------------------------
        # Handle ASCII raster (.asc) files
        # -------------------------
        if is_asc_format:
            # Filter to only .asc files (exclude any .nc files in the list)
            asc_files = [f for f in files if f.endswith('.asc')]
            
            # Get units from config if available
            provider = self.cfg.dataset.lower()
            units = None
            try:
                param_info = self.cfg.dsinfo[provider]['variables'].get(variable, {})
                units = param_info.get('units')
            except (AttributeError, KeyError, TypeError):
                pass
            
            # Load all .asc files as a time series
            dset = read_asc_timeseries(asc_files, varname=variable, units=units)
            
            # Preserve dataset attributes before extraction (xarray indexing removes them)
            original_attrs = dset.attrs.copy()
            
            # Apply extraction mode
            if mode == "point":
                lon, lat = self._extract_params
                # Find nearest x, y indices
                ix = np.abs(dset['x'].values - lon).argmin()
                iy = np.abs(dset['y'].values - lat).argmin()
                # Extract point
                dset = dset.isel(x=ix, y=iy)
                
            elif mode == "box":
                box = self._extract_params
                lat_min, lat_max = box['lat_min'], box['lat_max']
                lon_min, lon_max = box['lon_min'], box['lon_max']
                # For ASCII data with curvilinear lat/lon, use 2D coordinate indexing
                # Find y,x indices based on geographic lat/lon coordinates
                lat_2d = dset['lat'].values
                lon_2d = dset['lon'].values
                mask = (lat_2d >= lat_min) & (lat_2d <= lat_max) & (lon_2d >= lon_min) & (lon_2d <= lon_max)
                # Convert mask to DataArray with proper dimensions
                mask_da = xr.DataArray(mask, dims=('y', 'x'), coords={'y': dset['y'], 'x': dset['x']})
                dset = dset.where(mask_da, drop=True)
                
            elif mode == "shapefile":
                # For shapefile mode with .asc, would need to convert to lat/lon mask
                pass
            
            # Restore dataset attributes after extraction
            dset.attrs.update(original_attrs)
            
            # Apply time subsetting from config
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
                parallel=False,  # point preproc is tiny; parallel could be True on dask cluster
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
                    # mask may be (ny,nx) and ds dims are (y,x)
                    # create DataArray mask aligned to y,x
                    mask_da = xr.DataArray(mask, dims=("y", "x"))
                    # where keeps coords; drop=False keeps dims even if all-NaN
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
        if df is None:
            if self.dataset is None:
                raise ValueError("No dataset loaded")
            # convert to dataframe (may be large)
            df = self.dataset.to_dataframe().reset_index()
        df.to_csv(filename, index=False)
        print(f"Saved CSV to {filename}")
