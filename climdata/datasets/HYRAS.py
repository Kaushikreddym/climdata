import os
import pandas as pd
import xarray as xr
from datetime import datetime
from omegaconf import DictConfig
from climdata.utils.utils_download import find_nearest_xy, fetch_dwd
import geopandas as gpd

class HYRASmirror:
    def __init__(self, cfg: DictConfig):
        self.cfg = cfg
        self.dataset = None
        self.variables = cfg.variables
        self.files = []

    def fetch(self, variable: str):
        """
        Download HYRAS NetCDF files for a given variable and time range.
        """
        fetch_dwd(self.cfg,variable)
        # Build file list for the variable and time range
        param_mapping = self.cfg.dsinfo
        provider = self.cfg.dataset.lower()
        parameter_key = variable
        param_info = param_mapping[provider]['variables'][parameter_key]
        prefix = param_info["prefix"]
        version = param_info["version"]
        start_year = datetime.fromisoformat(self.cfg.time_range.start_date).year
        end_year = datetime.fromisoformat(self.cfg.time_range.end_date).year
        files = []
        for year in range(start_year, end_year + 1):
            file_name = f"{prefix}_{year}_{version}_de.nc"
            files.append(os.path.join(self.cfg.data_dir, provider, parameter_key.upper(), file_name))
        self.files = files
        return files

    def load(self, variable: str):
        """
        Load HYRAS NetCDFs for a given variable into a single xarray Dataset.
        """
        files = self.fetch(variable)
        datasets = []
        for f in files:
            if not os.path.exists(f):
                print(f"File not found: {f}")
                continue
            try:
                ds = xr.open_dataset(f)
                datasets.append(ds)
            except Exception as e:
                print(f"Skipping file {f} due to error: {e}")
        if not datasets:
            raise RuntimeError(f"No datasets could be loaded for {variable}.")
        dset = xr.concat(datasets, dim="time")
        dset[variable] = dset[variable].transpose("time", "y", "x")
        self.dataset = dset
        return self.dataset

    def extract(self, *, point=None, box=None, shapefile=None, buffer_km=0.0):
        """
        Extract data from the loaded HYRAS dataset.

        Parameters
        ----------
        point : tuple (lon, lat), optional
            Extracts a time series at the nearest grid point.
        box : dict with lat/lon bounds, optional
            Example: {"lat_min": 47, "lat_max": 49, "lon_min": 10, "lon_max": 12}
        shapefile : str, optional
            Path to a shapefile to clip the dataset spatially.
        buffer_km : float, optional
            Buffer distance (in kilometers) applied to the shapefile before clipping.
        """
        if self.dataset is None:
            raise ValueError("No dataset loaded. Call `load()` first.")
        ds = self.dataset

        # Point extraction
        if point is not None:
            lat, lon = point[1], point[0]
            iy, ix = find_nearest_xy(ds, lat, lon)
            print(f"📌 Nearest grid point at (y,x)=({iy},{ix})")
            ts = ds.isel(x=ix, y=iy)
            self.dataset = ts
            return ts

        # Box extraction
        elif box is not None:
            # Validate box keys
            if not all(k in box for k in ["lat_min", "lat_max", "lon_min", "lon_max"]):
                raise ValueError("Box must contain lat_min, lat_max, lon_min, lon_max.")

            # Find nearest indices for box boundaries
            iy_min, ix_min = find_nearest_xy(ds, box["lat_min"], box["lon_min"])
            iy_max, ix_max = find_nearest_xy(ds, box["lat_max"], box["lon_max"])
            # print(iy_min,ix_min,iy_max,ix_max)
            # Ensure proper ordering
            y_start, y_end = sorted([iy_min, iy_max])
            x_start, x_end = sorted([ix_min, ix_max])

            # Extract subset using indices
            dset_box = ds.isel(y=slice(y_start, y_end + 1), x=slice(x_start, x_end + 1))

            print(f"📦 Extracted curvilinear box with shape: {dset_box.dims}")
            self.dataset = dset_box
            return dset_box

        # Shapefile extraction
        elif shapefile is not None:
            """
            Clip a curvilinear xarray dataset using a shapefile with optional buffer in km.
            Works for 2D lat/lon coordinates.
            """
            # Read shapefile
            gdf = gpd.read_file(shapefile)
            
            # Apply buffer if needed
            if buffer_km > 0:
                gdf = gdf.to_crs(epsg=3857)
                gdf["geometry"] = gdf.buffer(buffer_km * 1000)
                gdf = gdf.to_crs(epsg=4326)
            
            # Flatten 2D lat/lon arrays
            lat_vals = ds['lat'].values
            lon_vals = ds['lon'].values
            
            # Create mask: True inside shapefile
            mask = np.zeros_like(lat_vals, dtype=bool)
            for polygon in gdf.geometry:
                # vectorized check for points inside polygon
                inside = np.array([polygon.contains(Point(lon, lat)) 
                                for lon, lat in zip(lon_vals.ravel(), lat_vals.ravel())])
                mask |= inside.reshape(lat_vals.shape)
            
            # Apply mask
            ds_clipped = ds.where(mask)
            return ds_clipped

        else:
            raise NotImplementedError("Must provide either point, box, or shapefile.")
        
    def save_csv(self, filename, df=None):
        """
        Save the extracted time series to CSV.
        """
        if df is None:
            if self.dataset is None:
                raise ValueError("No dataset loaded or extracted.")
            # If dataset is a DataArray, convert to DataFrame
            if isinstance(self.dataset, xr.Dataset):
                df = self.dataset.to_dataframe().reset_index()
            else:
                raise ValueError("Please provide a DataFrame or extract a point first.")
        df.to_csv(filename, index=False)
        print(f"Saved CSV to {filename}")