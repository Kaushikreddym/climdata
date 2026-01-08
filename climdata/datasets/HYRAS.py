# import os
# import pandas as pd
# import xarray as xr
# from datetime import datetime
# from omegaconf import DictConfig
# from climdata.utils.utils_download import find_nearest_xy, fetch_dwd
# import geopandas as gpd

# class HYRASmirror:
#     def __init__(self, cfg: DictConfig):
#         self.cfg = cfg
#         self.dataset = None
#         self.variables = cfg.variables
#         self.files = []
#         self._extract_mode = None
#         self._extract_params = None

#     def fetch(self, variable: str):
#         """
#         Download HYRAS NetCDF files for a given variable and time range.
#         """
#         fetch_dwd(self.cfg,variable)
#         # Build file list for the variable and time range
#         param_mapping = self.cfg.dsinfo
#         provider = self.cfg.dataset.lower()
#         parameter_key = variable
#         param_info = param_mapping[provider]['variables'][parameter_key]
#         prefix = param_info["prefix"]
#         version = param_info["version"]
#         start_year = datetime.fromisoformat(self.cfg.time_range.start_date).year
#         end_year = datetime.fromisoformat(self.cfg.time_range.end_date).year
#         files = []
#         for year in range(start_year, end_year + 1):
#             file_name = f"{prefix}_{year}_{version}_de.nc"
#             files.append(os.path.join(self.cfg.data_dir, provider, parameter_key.upper(), file_name))
#         self.files = files
#         return files

#     def load(self, variable: str):
#         files = self.fetch(variable)

#         def preprocess(ds):
#             # force transpose first
#             if variable in ds:
#                 ds[variable] = ds[variable].transpose("time", "y", "x")

#             # apply spatial extraction here
#             ds = self._extract_preprocess(ds)

#             return ds

#         # Open files with preprocess
#         dset = xr.open_mfdataset(
#             files,
#             combine="nested",
#             concat_dim="time",
#             preprocess=preprocess,
#             engine="netcdf4",
#             parallel=False,
#         )
#         if 'pr' in dset: 
#             dset['pr'].attrs['units'] = "mm/day"
#         self.dataset = dset
#         return dset


#     def extract(self, *, point=None, box=None, shapefile=None, buffer_km=0.0):
#         """Store extraction instructions; extraction happens per-file in preprocess()."""

#         if point is not None:
#             lon, lat = point
#             self._extract_mode = "point"
#             self._extract_params = (lon, lat)

#         elif box is not None:
#             self._extract_mode = "box"
#             self._extract_params = box

#         elif shapefile is not None:
#             gdf = gpd.read_file(shapefile) if isinstance(shapefile, str) else shapefile

#             if buffer_km > 0:
#                 gdf = gdf.to_crs(epsg=3857)
#                 gdf["geometry"] = gdf.buffer(buffer_km * 1000)
#                 gdf = gdf.to_crs(epsg=4326)

#             self._extract_mode = "shapefile"
#             self._extract_params = gdf

#         else:
#             raise ValueError("Must provide point, box, or shapefile.")

#         return self
#     def _extract_preprocess(self, ds):
#         """Apply point/box/shapefile extraction to a single HYRAS file."""

#         mode = self._extract_mode
#         params = self._extract_params

#         if mode is None:
#             return ds   # no extraction requested

#         # ---- point extraction ----
#         if mode == "point":
#             lon, lat = params
#             iy, ix = find_nearest_xy(ds, lat, lon)
#             return ds.isel(x=ix, y=iy)

#         # ---- box extraction ----
#         elif mode == "box":
#             box = params
#             iy_min, ix_min = find_nearest_xy(ds, box["lat_min"], box["lon_min"])
#             iy_max, ix_max = find_nearest_xy(ds, box["lat_max"], box["lon_max"])
#             y0, y1 = sorted([iy_min, iy_max])
#             x0, x1 = sorted([ix_min, ix_max])
#             return ds.isel(y=slice(y0, y1 + 1), x=slice(x0, x1 + 1))

#         # ---- shapefile extraction ----
#         elif mode == "shapefile":
#             gdf = params
#             # flatten coords
#             latv = ds["lat"].values
#             lonv = ds["lon"].values
#             mask = np.zeros_like(latv, dtype=bool)

#             for geom in gdf.geometry:
#                 inside = np.array([
#                     geom.contains(Point(lon, lat))
#                     for lon, lat in zip(lonv.ravel(), latv.ravel())
#                 ])
#                 mask |= inside.reshape(latv.shape)

#             return ds.where(mask)

#         return ds

#     def save_csv(self, filename, df=None):
#         """
#         Save the extracted time series to CSV.
#         """
#         if df is None:
#             if self.dataset is None:
#                 raise ValueError("No dataset loaded or extracted.")
#             # If dataset is a DataArray, convert to DataFrame
#             if isinstance(self.dataset, xr.Dataset):
#                 df = self.dataset.to_dataframe().reset_index()
#             else:
#                 raise ValueError("Please provide a DataFrame or extract a point first.")
#         df.to_csv(filename, index=False)
#         print(f"Saved CSV to {filename}")

import os
from datetime import datetime
from typing import Dict, Optional, Tuple, List

import numpy as np
import xarray as xr
import geopandas as gpd
import pandas as pd
import rasterio.features as rfeatures
from rasterio.transform import from_bounds
from shapely.geometry import mapping

from omegaconf import DictConfig
from climdata.utils.utils_download import find_nearest_xy, fetch_dwd


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
        """Download HYRAS NetCDF files for a given variable and time range and return file list."""
        # keep your fetch behavior (calls fetch_dwd)
        fetch_dwd(self.cfg, variable)

        provider = self.cfg.dataset.lower()
        param_info = self.cfg.dsinfo[provider]['variables'][variable]
        prefix = param_info["prefix"]
        version = param_info["version"]

        start_year = datetime.fromisoformat(self.cfg.time_range.start_date).year
        end_year = datetime.fromisoformat(self.cfg.time_range.end_date).year

        files = []
        for year in range(start_year, end_year + 1):
            file_name = f"{prefix}_{year}_{version}_de.nc"
            files.append(os.path.join(self.cfg.data_dir, provider, variable.upper(), file_name))

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

        - If extraction mode == 'point' -> uses open_mfdataset with preprocess (fast)
        - If extraction mode in ('box','shapefile') -> open per-file, apply index/mask, concat
        """

        files = self.fetch(variable)
        if not files:
            raise FileNotFoundError(f"No files found for variable {variable}")

        mode = self._extract_mode

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
