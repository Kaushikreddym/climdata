import pandas as pd
import numpy as np
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


from google.oauth2 import service_account
from googleapiclient.discovery import build
from googleapiclient.http import MediaIoBaseDownload

import io
import requests
from scipy.spatial import cKDTree
import argparse
import re

import requests
from bs4 import BeautifulSoup
import concurrent.futures

warnings.filterwarnings("ignore", category=Warning)

def list_drive_files(folder_id, service):
    """
    List all files in a Google Drive folder, handling pagination.
    """
    files = []
    page_token = None

    while True:
        results = service.files().list(
            q=f"'{folder_id}' in parents and trashed = false",
            fields="files(id, name), nextPageToken",
            pageToken=page_token
        ).execute()

        files.extend(results.get("files", []))
        page_token = results.get("nextPageToken", None)

        if not page_token:
            break

    return files
def download_drive_file(file_id, local_path, service):
    """
    Download a single file from Drive to a local path.
    """
    request = service.files().get_media(fileId=file_id)
    os.makedirs(os.path.dirname(local_path), exist_ok=True)

    with io.FileIO(local_path, 'wb') as fh:
        downloader = MediaIoBaseDownload(fh, request)

        done = False
        while not done:
            status, done = downloader.next_chunk()
            print(f"   → Download {int(status.progress() * 100)}% complete")

def fetch_dwd(var_cfg,var):
    """Download HYRAS data for one variable and a list of years."""
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

    # Parse dates & extract unique years
    start_year = datetime.fromisoformat(start_date).year
    end_year = datetime.fromisoformat(end_date).year
    years = list(range(start_year, end_year + 1))

    # output_file = cfg.output.filename
    os.makedirs(parameter_key, exist_ok=True)

    for year in years:
        file_name = f"{prefix}_{year}_{version}_de.nc"
        file_url = f"{base_url}{file_name}"
        local_path = os.path.join(var_cfg.data_dir,provider,parameter_key.upper(), file_name)
        os.makedirs(os.path.dirname(local_path), exist_ok=True)
        print(f"⬇️  Checking: {file_url}")

        # Check if file exists on server first (HEAD request)
        head = requests.head(file_url)
        if head.status_code != 200:
            raise FileNotFoundError(f"❌ Not found on server: {file_url} (HTTP {head.status_code})")

        if os.path.exists(local_path):
            print(f"✔️  Exists locally: {local_path}")
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
            raise RuntimeError(f"❌ Failed download: {file_url} — {e}")

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

def extract_ts_dwd(cfg: DictConfig):
    param_mapping = cfg.mappings
    provider = cfg.dataset.lower()
    parameter_key = cfg.weather.parameter
    # Validate provider and parameter

    param_info = param_mapping[provider]['variables'][parameter_key]
    prefix = param_info["prefix"]
    version = param_info["version"]

    start_date = cfg.time_range.start_date
    end_date = cfg.time_range.end_date

    # Parse dates & extract unique years
    start_year = datetime.fromisoformat(start_date).year
    end_year = datetime.fromisoformat(end_date).year
    years = list(range(start_year, end_year + 1))
    files=[]
    for year in years:
        file_name = f"{prefix}_{year}_{version}_de.nc"
        files.append(os.path.join(cfg.data_dir,provider,parameter_key.upper(), file_name))

    if not files:
        raise FileNotFoundError(f"No NetCDF files found for {parameter_key}")

    target_lat = cfg.location.lat
    target_lon = cfg.location.lon

    ts_list = []

    for f in files:
        print(f"📂 Opening: {f}")
        ds = xr.open_dataset(f)

        # Dimensions: (time, y, x) or (time, x, y)
        # lat/lon: 2D
        time_name = [x for x in ds.coords if "time" in x.lower()][0]
        
        iy, ix = find_nearest_xy(ds, target_lat, target_lon)

        print(f"📌 Nearest grid point at (y,x)=({iy},{ix})")
        
        ts = ds[parameter_key].isel(x=ix, y=iy)  # watch order: dims must match

        df = ts.to_dataframe().reset_index()[[time_name, parameter_key]]
        ts_list.append(df)

    # Combine all time series
    ts_all = pd.concat(ts_list).sort_values(by=time_name).reset_index(drop=True)
    
    # Slice on combined DataFrame
    ts_all[time_name] = pd.to_datetime(ts_all[time_name])
    mask = (ts_all[time_name] >= start_date) & (ts_all[time_name] <= end_date)
    ts_all = ts_all.loc[mask].reset_index(drop=True)

    out_dir = hydra.utils.to_absolute_path(cfg.output.out_dir)
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, cfg.output.filename)
    
    ts_all["variable"] = param_info['name']
    ts_all["latitude"] = target_lat
    ts_all["longitude"] = target_lon
    ts_all['source'] = provider.upper()
    ts_all['units'] = ts.attrs['units']
    ts_all.rename(columns={param_info['name']: 'value'}, inplace=True)
    ts_all = ts_all[["latitude", "longitude", "time", "source", "variable", "value",'units']]
    ts_all.to_csv(out_path, index=False)
    print(f"✅ Saved time series to: {out_path}")

    return ts_all

import os
from omegaconf import DictConfig

def get_output_filename(cfg, output_type="nc", lat=None, lon=None, shp_name = None):
    """
    Generate output filename based on config, output type, and extraction mode.
    output_type: "nc", "csv", or "zarr"
    """
    if output_type == "csv":
        template = cfg.output.filename_csv
    elif output_type == "zarr":
        template = cfg.output.filename_zarr
    else:
        template = cfg.output.filename_nc

    # If lat/lon are provided, use point template
    if lat is not None and lon is not None:
        filename = template.format(
            provider=cfg.dataset,
            parameter="surface",
            lat=f"{lat}",
            lon=f"{lon}",
            start=cfg.time_range.start_date.replace("-", ""),
            end=cfg.time_range.end_date.replace("-", ""),
        )
    elif shp_name is not None:
        filename = template.format(
            provider=cfg.dataset,
            parameter="surface",
            lat_range=f"{shp_name}",
            lon_range=f"{shp_name}",
            start=cfg.time_range.start_date.replace("-", ""),
            end=cfg.time_range.end_date.replace("-", ""),
        )
    else:
        # Use region bounds
        region_bounds = cfg.bounds[cfg.region]
        filename = template.format(
            provider=cfg.dataset,
            parameter="surface",
            lat_range=f"{region_bounds['lat_min']}-{region_bounds['lat_max']}",
            lon_range=f"{region_bounds['lon_min']}-{region_bounds['lon_max']}",
            start=cfg.time_range.start_date.replace("-", ""),
            end=cfg.time_range.end_date.replace("-", ""),
        )
    return filename