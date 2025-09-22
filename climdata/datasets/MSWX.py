import pandas as pd
import numpy as np
from wetterdienst import Settings
from wetterdienst.provider.dwd.observation import DwdObservationRequest
import geemap
import ee
import ipdb
import geopandas as gpd
from omegaconf import DictConfig
import os
import yaml
import time
from tqdm import tqdm
import warnings
from datetime import datetime, timedelta
import xarray as xr
import hydra
from omegaconf import DictConfig
import pint
import pint_pandas

from google.oauth2 import service_account
from googleapiclient.discovery import build
from googleapiclient.http import MediaIoBaseDownload

from climdata.utils.utils_download import list_drive_files, download_drive_file

import io
import requests
from scipy.spatial import cKDTree
import argparse
import re

import requests
from bs4 import BeautifulSoup
import concurrent.futures

import gzip
# from utils.utils import *
# from datasets.datasets import *
import rioxarray
from shapely.geometry import mapping

warnings.filterwarnings("ignore", category=Warning)

import cf_xarray

class MSWXmirror:
    def __init__(self, var_cfg: DictConfig):
        self.var_cfg = var_cfg
        self.files = []
        self.dataset = None

    def _fix_coords(self, ds: xr.Dataset | xr.DataArray) -> xr.Dataset | xr.DataArray:
        """
        Ensure latitude is ascending and longitude is in the range [0, 360].

        Parameters
        ----------
        ds : xr.Dataset or xr.DataArray
            Input dataset or dataarray with latitude and longitude coordinates.

        Returns
        -------
        xr.Dataset or xr.DataArray
            Dataset with latitude ascending and longitude wrapped to [0, 360].
        """
        # Flip latitude to ascending
        ds = ds.cf.sortby("latitude")

        # Wrap longitude into [0, 360]
        lon_name = ds.cf["longitude"].name
        ds = ds.assign_coords({lon_name: ds.cf["longitude"] % 360})

        # Sort by longitude
        ds = ds.sortby(lon_name)

        return ds


    def fetch(self):
        param_mapping = self.var_cfg.mappings
        provider = self.var_cfg.dataset.lower()
        parameter_key = self.var_cfg.weather.parameter

        param_info = param_mapping[provider]['variables'][parameter_key]
        folder_id = param_info["folder_id"]

        start_date = self.var_cfg.time_range.start_date
        end_date = self.var_cfg.time_range.end_date

        start = datetime.fromisoformat(start_date)
        end = datetime.fromisoformat(end_date)

        expected_files = []
        current = start
        while current <= end:
            doy = current.timetuple().tm_yday
            basename = f"{current.year}{doy:03d}.nc"
            expected_files.append(basename)
            current += timedelta(days=1)

        output_dir = self.var_cfg.data_dir
        provider = self.var_cfg.dataset.lower()
        parameter_key = self.var_cfg.weather.parameter
        local_files = []
        missing_files = []

        for basename in expected_files:
            local_path = os.path.join(output_dir, provider, parameter_key, basename)
            if os.path.exists(local_path):
                local_files.append(basename)
            else:
                missing_files.append(basename)

        if not missing_files:
            print(f"✅ All {len(expected_files)} files already exist locally. No download needed.")
            self.files = local_files
            return local_files

        print(f"📂 {len(local_files)} exist, {len(missing_files)} missing — fetching from Drive...")

        SCOPES = ['https://www.googleapis.com/auth/drive.readonly']
        creds = service_account.Credentials.from_service_account_file(
            param_mapping[provider].params.google_service_account, scopes=SCOPES
        )
        service = build('drive', 'v3', credentials=creds)

        drive_files = list_drive_files(folder_id, service)
        valid_filenames = set(missing_files)
        files_to_download = [f for f in drive_files if f['name'] in valid_filenames]

        if not files_to_download:
            print(f"⚠️ None of the missing files found in Drive. Check folder & date range.")
            self.files = local_files
            return local_files

        for file in files_to_download:
            filename = file['name']
            local_path = os.path.join(output_dir, provider, parameter_key, filename)
            print(f"⬇️ Downloading {filename} ...")
            download_drive_file(file['id'], local_path, service)
            local_files.append(filename)

        self.files = local_files
        return local_files

    def load(self):
        param_mapping = self.var_cfg.mappings
        provider = self.var_cfg.dataset.lower()
        parameter_key = self.var_cfg.weather.parameter
        region = self.var_cfg.region
        bounds = self.var_cfg.bounds[region]

        param_info = param_mapping[provider]['variables'][parameter_key]
        output_dir = self.var_cfg.data_dir
        valid_dsets = []

        for f in self.files:
            local_path = os.path.join(output_dir, provider, parameter_key, f)
            try:
                ds = xr.open_dataset(local_path, chunks='auto', engine='netcdf4')[param_info.name]
                valid_dsets.append(ds)
            except Exception as e:
                print(f"Skipping file due to error: {f}\n{e}")

        dset = xr.concat(valid_dsets, dim='time')
        dset = dset.transpose('time', 'lat', 'lon')
        self.dataset = self._fix_coords(dset)
        return dset

    def to_zarr(self, zarr_filename):
        if self.dataset is None:
            raise ValueError("No dataset loaded. Call `load()` before `to_zarr()`.")

        var_name = self.var_cfg.weather.parameter
        dataset_name = self.var_cfg.dataset
        region = self.var_cfg.region

        # Add standard units metadata
        if var_name == 'pr':
            self.dataset.attrs['units'] = 'mm/day'
        elif var_name in ['tas', 'tasmax', 'tasmin']:
            self.dataset.attrs['units'] = 'degC'

        zarr_path = os.path.join("data/MSWX/", zarr_filename)
        os.makedirs(os.path.dirname(zarr_path), exist_ok=True)

        print(f"💾 Saving {var_name} to Zarr: {zarr_path}")
        self.dataset.to_zarr(zarr_path, mode="w")

    def extract(self, *, point=None, box=None, shapefile=None, buffer_km=0.0):
        """
        Extract a subset of the dataset by point, bounding box, or shapefile.

        Parameters
        ----------
        point : tuple(float, float), optional
            (lon, lat) coordinates for a single point.
        box : tuple(float, float, float, float), optional
            (min_lon, min_lat, max_lon, max_lat) bounding box.
        shapefile : str or geopandas.GeoDataFrame, optional
            Path to shapefile or a GeoDataFrame.
        buffer_km : float, optional
            Buffer distance in kilometers (for point or shapefile).
        
        Returns
        -------
        xarray.Dataset or xarray.DataArray
            Subset of the dataset.
        """
        if self.dataset is None:
            raise ValueError("No dataset loaded. Call `load()` first.")

        ds = self.dataset.rio.write_crs("EPSG:4326", inplace=False)

        if point is not None:
            lon, lat = point
            if buffer_km > 0:
                # buffer around point
                buffer_deg = buffer_km / 111  # rough conversion km→degrees
                ds_subset = ds.sel(
                    lon=slice(lon-buffer_deg, lon+buffer_deg),
                    lat=slice(lat-buffer_deg, lat+buffer_deg)
                )
            else:
                ds_subset = ds.sel(lon=lon, lat=lat, method="nearest")

        elif box is not None:
            min_lon, min_lat, max_lon, max_lat = box
            ds_subset = ds.sel(
                lon=slice(min_lon, max_lon),
                lat=slice(min_lat, max_lat)
            )

        elif shapefile is not None:
            if isinstance(shapefile, str):
                gdf = gpd.read_file(shapefile)
            else:
                gdf = shapefile
            if buffer_km > 0:
                gdf = gdf.to_crs(epsg=3857)  # project to meters
                gdf["geometry"] = gdf.buffer(buffer_km * 1000)
                gdf = gdf.to_crs(epsg=4326)

            geom = [mapping(g) for g in gdf.geometry]
            ds_subset = ds.rio.clip(geom, gdf.crs, drop=True)

        else:
            raise ValueError("Must provide either point, box, or shapefile.")

        return ds_subset
    
    def to_dataframe(self, ds=None):
        """
        Convert extracted xarray dataset to a tidy dataframe.

        Parameters
        ----------
        ds : xr.DataArray or xr.Dataset, optional
            Dataset to convert. If None, use self.dataset.

        Returns
        -------
        pd.DataFrame
        """
        if ds is None:
            if self.dataset is None:
                raise ValueError("No dataset loaded. Call `load()` first or pass `ds`.")
            ds = self.dataset

        # If Dataset, pick first variable
        if isinstance(ds, xr.Dataset):
            if len(ds.data_vars) != 1:
                raise ValueError("Dataset has multiple variables. Please select one.")
            ds = ds[list(ds.data_vars)[0]]

        df = ds.to_dataframe().reset_index()

        # Keep only relevant cols
        df = df[["time", "lat", "lon", ds.name]]

        # Rename
        df = df.rename(columns={
            "lat": "latitude",
            "lon": "longitude",
            ds.name: "value"
        })
        return df

    def format(self, df):
        """
        Format dataframe into standard schema.
        """
        df = df.copy()
        df["variable"] = self.var_cfg.weather.parameter
        df["source"] = self.var_cfg.dataset.upper()
        df["units"] = self.dataset.attrs.get("units", "unknown")

        df = df[["latitude", "longitude", "time", "source", "variable", "value", "units"]]
        return df
    
