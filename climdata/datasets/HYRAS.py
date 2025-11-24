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
        self._extract_mode = None
        self._extract_params = None

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
        files = self.fetch(variable)

        def preprocess(ds):
            # force transpose first
            if variable in ds:
                ds[variable] = ds[variable].transpose("time", "y", "x")

            # apply spatial extraction here
            ds = self._extract_preprocess(ds)

            return ds

        # Open files with preprocess
        dset = xr.open_mfdataset(
            files,
            combine="nested",
            concat_dim="time",
            preprocess=preprocess,
            engine="netcdf4",
            parallel=False,
        )

        self.dataset = dset
        return dset


    def extract(self, *, point=None, box=None, shapefile=None, buffer_km=0.0):
        """Store extraction instructions; extraction happens per-file in preprocess()."""

        if point is not None:
            lon, lat = point
            self._extract_mode = "point"
            self._extract_params = (lon, lat)

        elif box is not None:
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

        return self
    def _extract_preprocess(self, ds):
        """Apply point/box/shapefile extraction to a single HYRAS file."""

        mode = self._extract_mode
        params = self._extract_params

        if mode is None:
            return ds   # no extraction requested

        # ---- point extraction ----
        if mode == "point":
            lon, lat = params
            iy, ix = find_nearest_xy(ds, lat, lon)
            return ds.isel(x=ix, y=iy)

        # ---- box extraction ----
        elif mode == "box":
            box = params
            iy_min, ix_min = find_nearest_xy(ds, box["lat_min"], box["lon_min"])
            iy_max, ix_max = find_nearest_xy(ds, box["lat_max"], box["lon_max"])
            y0, y1 = sorted([iy_min, iy_max])
            x0, x1 = sorted([ix_min, ix_max])
            return ds.isel(y=slice(y0, y1 + 1), x=slice(x0, x1 + 1))

        # ---- shapefile extraction ----
        elif mode == "shapefile":
            gdf = params
            # flatten coords
            latv = ds["lat"].values
            lonv = ds["lon"].values
            mask = np.zeros_like(latv, dtype=bool)

            for geom in gdf.geometry:
                inside = np.array([
                    geom.contains(Point(lon, lat))
                    for lon, lat in zip(lonv.ravel(), latv.ravel())
                ])
                mask |= inside.reshape(latv.shape)

            return ds.where(mask)

        return ds

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