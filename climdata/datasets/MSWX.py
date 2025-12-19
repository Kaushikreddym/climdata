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
        self._extract_mode = None
        self._extract_params = None

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
    def _extract_preprocess(self, ds):
        """Apply extraction to a single-daily dataset during preprocessing."""
        
        # Fix coords first
        ds = self._fix_coords(ds)

        # ---- Point extraction ----
        if self._extract_mode == "point":
            lon, lat, buffer_deg = self._extract_params
            if buffer_deg > 0:
                ds = ds.sel(
                    lon=slice(lon-buffer_deg, lon+buffer_deg),
                    lat=slice(lat-buffer_deg, lat+buffer_deg),
                ).mean(["lat", "lon"])
            else:
                ds = ds.sel(lon=lon, lat=lat, method="nearest")

        # ---- Box extraction ----
        elif self._extract_mode == "box":
            box = self._extract_params
            ds = ds.sel(
                lon=slice(box["lon_min"], box["lon_max"]),
                lat=slice(box["lat_min"], box["lat_max"]),
            )

        # ---- Shapefile extraction ----
        elif self._extract_mode == "shapefile":
            gdf = self._extract_params
            
            # Suppose your dataset uses 'lon' and 'lat' as coordinates
            ds = ds.rio.set_spatial_dims(x_dim="lon", y_dim="lat")

            # Also ensure CRS is set
            ds = ds.rio.write_crs("EPSG:4326", inplace=True)
            
            clipped_list = []
            for geom in gdf.geometry:
                clipped = ds.rio.clip([mapping(geom)], gdf.crs, drop=True)
                clipped_list.append(clipped)

            ds = xr.concat(clipped_list, dim="geom_id")

        return ds
    def extract(self, *, point=None, box=None, shapefile=None, buffer_km=0.0):
        """Store extraction instructions; the actual extraction happens during load()."""

        if point is not None:
            lon, lat = point
            buffer_deg = buffer_km / 111.0
            self._extract_mode = "point"
            self._extract_params = (lon, lat, buffer_deg)

        elif box is not None:
            self._extract_mode = "box"
            self._extract_params = box

        elif shapefile is not None:
            if isinstance(shapefile, str):
                gdf = gpd.read_file(shapefile)
            else:
                gdf = shapefile
            
            if buffer_km > 0:
                gdf = gdf.to_crs(epsg=3857)
                gdf["geometry"] = gdf.buffer(buffer_km * 1000)
                gdf = gdf.to_crs(epsg=4326)

            self._extract_mode = "shapefile"
            self._extract_params = gdf

        else:
            raise ValueError("Must provide point, box, or shapefile.")

        return self

    def load(self, variable: str):
        """
        Load MSWX NetCDF files for a given variable into a single xarray Dataset using open_mfdataset.
        This method supports lazy loading, parallel processing, and large numbers of files efficiently.

        Args:
            variable (str): Variable name as defined in cfg.variables.

        Returns:
            xr.Dataset: Concatenated dataset along the 'time' dimension with fixed coordinates.
        """
        # Get folder ID and list of files
        folder_id = self.cfg.dsinfo["mswx"]["variables"][variable]["folder_id"]
        files = self.fetch(folder_id, variable)
        if not files:
            raise RuntimeError(f"No files found for variable '{variable}' in Drive or local directory.")

        # Full paths
        file_paths = [
            os.path.join(self.cfg.data_dir, self.cfg.dataset.lower(), variable, f)
            for f in files
        ]

        # MSWX internal variable name
        varname = self.cfg.dsinfo[self.cfg.dataset].variables[variable].name

        # Optional: preprocess each file (e.g., rename variable)
        def preprocess(ds):
            ds = self._extract_preprocess(ds)
            return ds[[varname]].rename({varname: variable})
        import dask
        # Open all files as a single dataset
        try:
            @dask.delayed
            def open_point(f):
                ds = preprocess(xr.open_dataset(f, engine="h5netcdf"))
                return ds
            batch_size = 500  # process 500 files at a time
            dsets = []

            for i in range(0, len(file_paths), batch_size):
                batch_files = file_paths[i:i+batch_size]
                delayed_batch = [dask.delayed(open_point)(f) for f in batch_files]
                batch_ds = list(dask.compute(*delayed_batch))
                dsets.extend(batch_ds)

            dset = xr.concat(dsets, dim="time")

            # dset = xr.open_mfdataset(
            #     file_paths,
            #     combine="nested",
            #     concat_dim="time",
            #     parallel=True,            # uses Dask for parallel reading
            #     engine="h5netcdf",       # faster than netcdf4
            #     # chunks = {"time": 90, "lat": 200, "lon": 200},  # quarterly chunks
            #     preprocess=preprocess
            # )
        except Exception as e:
            raise RuntimeError(f"Failed to load dataset for variable '{variable}': {e}")

        # Ensure consistent dimension order
        if self._extract_mode is not "point":
            dset = dset.transpose("time", "lat", "lon")

        # Store in the class
        self.dataset = dset
        return dset


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
            # print(f"Saved CSV to {filename}")
