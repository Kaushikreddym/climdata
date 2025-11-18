import pandas as pd
import geopandas as gpd
import os
from tqdm import tqdm
import warnings
from datetime import datetime, timedelta
import xarray as xr
from omegaconf import DictConfig

from google.oauth2 import service_account
from googleapiclient.discovery import build

from climdata.utils.utils_download import list_drive_files, download_drive_file
from shapely.geometry import mapping
import cf_xarray

warnings.filterwarnings("ignore", category=Warning)


class MSWXmirror:
    def __init__(self, cfg: DictConfig):
        self.cfg = cfg
        self.dataset = None
        self.variables = cfg.variables
        self.files = []

    def _fix_coords(self, ds: xr.Dataset | xr.DataArray):
        """Ensure latitude is ascending and longitude is in the range [0, 360]."""
        ds = ds.cf.sortby("latitude")
        lon_name = ds.cf["longitude"].name
        ds = ds.assign_coords({lon_name: ds.cf["longitude"] % 360})
        return ds.sortby(lon_name)

    def fetch(self, folder_id: str, variable: str):
        """
        Fetch MSWX files from Google Drive for a given variable.
        """
        start = datetime.fromisoformat(self.cfg.time_range.start_date)
        end = datetime.fromisoformat(self.cfg.time_range.end_date)

        expected_files = []
        current = start
        while current <= end:
            doy = current.timetuple().tm_yday
            basename = f"{current.year}{doy:03d}.nc"
            expected_files.append(basename)
            current += timedelta(days=1)

        output_dir = self.cfg.data_dir
        local_files, missing_files = [], []

        for basename in expected_files:
            local_path = os.path.join(output_dir,self.cfg.dataset.lower(), variable, basename)
            if os.path.exists(local_path):
                local_files.append(basename)
            else:
                missing_files.append(basename)
        
        if not missing_files:
            print(f"✅ All {len(expected_files)} {variable} files already exist locally.")
            return local_files

        print(f"📂 {len(local_files)} exist, {len(missing_files)} missing — fetching {variable} from Drive...")

        SCOPES = ['https://www.googleapis.com/auth/drive.readonly']
        creds = service_account.Credentials.from_service_account_file(
            self.cfg.dsinfo.mswx.params.google_service_account, scopes=SCOPES
        )
        service = build('drive', 'v3', credentials=creds)

        drive_files = list_drive_files(folder_id, service)
        valid_filenames = set(missing_files)
        files_to_download = [f for f in drive_files if f['name'] in valid_filenames]

        if not files_to_download:
            print(f"⚠️ No {variable} files found in Drive for requested dates.")
            return local_files

        for file in files_to_download:
            filename = file['name']
            local_path = os.path.join(output_dir, self.cfg.dataset, variable, filename)
            os.makedirs(os.path.dirname(local_path), exist_ok=True)
            print(f"⬇️ Downloading {filename} ...")
            download_drive_file(file['id'], local_path, service)
            local_files.append(filename)

        return local_files

    def load(self, variable: str):
        """
        Load MSWX NetCDFs for a given variable into a single xarray Dataset.
        """
        folder_id = self.cfg.dsinfo["mswx"]["variables"][variable]["folder_id"]
        files = self.fetch(folder_id, variable)
        datasets = []

        for f in files:
            local_path = os.path.join(self.cfg.data_dir, self.cfg.dataset.lower(), variable, f)
            try:
                ds = xr.open_dataset(local_path, chunks="auto", engine="netcdf4")[self.cfg.dsinfo[self.cfg.dataset].variables[variable].name]
                ds = ds.rename(variable)
                datasets.append(ds)
            except Exception as e:
                print(f"Skipping file {f} due to error: {e}")

        if not datasets:
            raise RuntimeError(f"No datasets could be loaded for {variable}.")

        dset = xr.concat(datasets, dim="time")
        dset = dset.transpose("time", "lat", "lon")
        dset = self._fix_coords(dset)

        self.dataset = dset
        return self.dataset

    def to_zarr(self, zarr_filename: str):
        if self.dataset is None:
            raise ValueError("No dataset loaded. Call `load()` first.")

        var_name = self.dataset.name
        if var_name == 'pr':
            self.dataset.attrs['units'] = 'mm/day'
        elif var_name in ['tas', 'tasmax', 'tasmin']:
            self.dataset.attrs['units'] = 'degC'

        zarr_path = os.path.join("data/MSWX", zarr_filename)
        os.makedirs(os.path.dirname(zarr_path), exist_ok=True)

        print(f"💾 Saving {var_name} to Zarr: {zarr_path}")
        self.dataset.to_zarr(zarr_path, mode="w")

    def extract(self, *, point=None, box=None, shapefile=None, buffer_km=0.0):
        if self.dataset is None:
            raise ValueError("No dataset loaded. Call `load()` first.")

        ds = self.dataset

        # Ensure CRS and spatial dimensions
        if "x" not in ds.dims or "y" not in ds.dims:
            ds = ds.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=False)
        if not ds.rio.crs:
            ds = ds.rio.write_crs("EPSG:4326", inplace=False)

        # ---- Point extraction ----
        if point is not None:
            lon, lat = point
            if buffer_km > 0:
                buffer_deg = buffer_km / 111
                ds_subset = ds.sel(
                    lon=slice(lon-buffer_deg, lon+buffer_deg),
                    lat=slice(lat-buffer_deg, lat+buffer_deg),
                )
            else:
                ds_subset = ds.sel(lon=lon, lat=lat, method="nearest")

        # ---- Box extraction ----
        elif box is not None:
            ds_subset = ds.sel(
                lon=slice(box["lon_min"], box["lon_max"]),
                lat=slice(box["lat_min"], box["lat_max"]),
            )

        # ---- Shapefile extraction ----
        elif shapefile is not None:
            # Read shapefile if path provided
            if isinstance(shapefile, str):
                gdf = gpd.read_file(shapefile)
            else:
                gdf = shapefile

            # Optional buffer in km
            if buffer_km > 0:
                gdf = gdf.to_crs(epsg=3857)
                gdf["geometry"] = gdf.buffer(buffer_km * 1000)
                gdf = gdf.to_crs(epsg=4326)

            # Create a new dimension for each geometry
            clipped_list = []
            for i, geom in enumerate(gdf.geometry):
                clipped = ds.rio.clip([mapping(geom)], gdf.crs, drop=True)
                clipped = clipped.expand_dims(geom_id=[i])
                if "geometry_name" in gdf.columns:
                    clipped = clipped.assign_coords(
                        geom_name=("geom_id", [gdf.loc[i, "geometry_name"]])
                    )
                clipped_list.append(clipped)

            # Concatenate along new "geom_id" dimension
            ds_subset = xr.concat(clipped_list, dim="geom_id")

        else:
            raise ValueError("Must provide either point, box, or shapefile.")

        self.dataset = ds_subset.to_dataset()
        return ds_subset


    # def to_dataframe(self, ds=None):
    #     if ds is None:
    #         if self.dataset is None:
    #             raise ValueError("No dataset loaded. Call `load()` first or pass `ds`.")
    #         ds = self.dataset

    #     if isinstance(ds, xr.Dataset):
    #         if len(ds.data_vars) != 1:
    #             raise ValueError("Dataset has multiple variables. Please select one.")
    #         ds = ds[list(ds.data_vars)[0]]

    #     df = ds.to_dataframe().reset_index()
    #     df = df[["time", "lat", "lon", ds.name]]
    #     df = df.rename(columns={"lat": "latitude", "lon": "longitude", ds.name: "value"})
    #     return df

    def _format(self, df):
        """Format dataframe for standardized output."""
        value_vars = [v for v in self.variables if v in df.columns]
        id_vars = [c for c in df.columns if c not in value_vars]

        df_long = df.melt(
            id_vars=id_vars,
            value_vars=value_vars,
            var_name="variable",
            value_name="value",
        )

        df_long["units"] = df_long["variable"].map(
            lambda v: self.dataset[v].attrs.get("units", "unknown")
            if v in self.dataset.data_vars
            else "unknown"
        )

        df_long["source"] = self.cfg.dataset

        cols = [
            "source",
            "table",
            "time",
            "lat",
            "lon",
            "variable",
            "value",
            "units",
        ]
        df_long = df_long[[c for c in cols if c in df_long.columns]]

        return df_long

    def save_csv(self, filename):
        if self.dataset is not None:
            df = self.dataset.to_dataframe().reset_index()
            df = self._format(df)
            df.to_csv(filename, index=False)
            print(f"Saved CSV to {filename}")
