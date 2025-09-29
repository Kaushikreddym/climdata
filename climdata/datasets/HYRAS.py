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
        param_mapping = self.cfg.mappings
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
            if not all(k in box for k in ["lat_min", "lat_max", "lon_min", "lon_max"]):
                raise ValueError("Box must contain lat_min, lat_max, lon_min, lon_max.")
            dset_box = ds.sel(
                y=slice(box["lat_max"], box["lat_min"]),  # y usually decreasing (north -> south)
                x=slice(box["lon_min"], box["lon_max"])
            )
            print(f"📦 Extracted box with shape: {dset_box.dims}")
            self.dataset = dset_box
            return dset_box

        # Shapefile extraction
        elif shapefile is not None:
            gdf = gpd.read_file(shapefile)

            if buffer_km > 0:
                gdf = gdf.to_crs(epsg=3857)  # project to meters
                gdf["geometry"] = gdf.buffer(buffer_km * 1000)  # buffer in meters
                gdf = gdf.to_crs(epsg=4326)  # back to lat/lon

            # Ensure dataset has CRS info for clipping
            if not ds.rio.crs:
                ds = ds.rio.write_crs("EPSG:4326")

            dset_clipped = ds.rio.clip(gdf.geometry, gdf.crs, drop=True)
            print(f"🗺️ Extracted shapefile area with dims: {dset_clipped.dims}")
            self.dataset = dset_clipped
            return dset_clipped

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