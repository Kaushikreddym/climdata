"""
W5E5 dataset access via ISIMIP client

W5E5 (WFDE5 over land merged with ERA5 over the ocean) is a global meteorological 
forcing dataset available through ISIMIP (Inter-Sectoral Impact Model Intercomparison Project).
It provides daily climate data at 0.5° resolution from 1979 onwards.

This module uses the isimip-client library to search and download W5E5 data from the 
ISIMIP data repository.
"""

import os
import xarray as xr
import pandas as pd
from pathlib import Path
from datetime import datetime
from omegaconf import DictConfig
from typing import Optional, Tuple, Dict, List
import warnings

warnings.filterwarnings("ignore", category=Warning)

class W5E5:
    """
    A class to download and process W5E5 climate data from ISIMIP repository.
    
    W5E5 is available through ISIMIP3a (historical observations) and is used as 
    bias adjustment reference for ISIMIP3b climate projections.
    
    Attributes
    ----------
    cfg : DictConfig
        Configuration containing lat, lon, variables, time_range, etc.
    ds : xr.Dataset
        Loaded xarray dataset
    client : ISIMIPClient
        ISIMIP API client for data access
    """
    
    def __init__(self, cfg: DictConfig):
        self.cfg = cfg
        self.ds = None
        self.client = None
        self.downloaded_files = []
        self._extract_mode = None
        self._extract_params = None
        
        # Initialize ISIMIP client
        try:
            from isimip_client.client import ISIMIPClient
            self.client = ISIMIPClient()
        except ImportError:
            raise ImportError(
                "isimip-client is required for W5E5 data access. "
                "Install it with: pip install isimip-client"
            )
    
    def fetch(self):
        """
        Search and download W5E5 files from ISIMIP repository for the requested
        variables and time range.
        
        Uses ISIMIP3a simulation round for W5E5 observational data.
        """
        print("🔍 Searching for W5E5 datasets in ISIMIP repository...")
        
        start_date = datetime.fromisoformat(self.cfg.time_range.start_date)
        end_date = datetime.fromisoformat(self.cfg.time_range.end_date)
        
        output_dir = Path(self.cfg.data_dir) / self.cfg.dataset.lower()
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Search for each variable separately
        for var in self.cfg.variables:
            print(f"\n📥 Fetching {var}...")
            
            # Map variable names to W5E5 names if needed
            w5e5_var = self._map_variable_name(var)
            
            # Search ISIMIP repository for W5E5 data
            # W5E5 is available in ISIMIP3a as observational input data
            try:
                response = self.client.datasets(
                    simulation_round='ISIMIP3a',
                    product='InputData',
                    climate_forcing='20crv3-w5e5',  # W5E5 version 2.0
                    climate_scenario='obsclim',  # Observed climate
                    climate_variable=w5e5_var
                )
                
                if not response.get('results'):
                    print(f"⚠️ No W5E5 datasets found for {var}")
                    continue
                
                # Get the first matching dataset
                dataset = response['results'][0]
                print(f"✅ Found dataset: {dataset.get('name', 'unnamed')}")
                
                # Filter files by date range
                for file_info in dataset.get('files', []):
                    file_path = file_info['path']
                    file_name = file_info['name']
                    
                    # Parse date from filename (W5E5 files typically contain year ranges)
                    # Example: w5e5v2.0_obsclim_tas_global_daily_1979_1989.nc
                    if self._is_file_in_date_range(file_name, start_date, end_date):
                        local_path = output_dir / var / file_name
                        local_path.parent.mkdir(parents=True, exist_ok=True)
                        
                        if local_path.exists():
                            print(f"  ✓ Already exists: {file_name}")
                            self.downloaded_files.append(str(local_path))
                        else:
                            print(f"  ⬇️ Downloading: {file_name}")
                            # Download directly using the file URL
                            self.client.download(
                                file_info['file_url'],
                                path=str(local_path.parent),
                                validate=False
                            )
                            self.downloaded_files.append(str(local_path))
                            
            except Exception as e:
                print(f"❌ Error fetching {var}: {str(e)}")
                continue
        
        print(f"\n✅ Downloaded {len(self.downloaded_files)} files")
    
    def load(self):
        """
        Load the downloaded W5E5 netCDF files into an xarray Dataset.
        Combines multiple files if necessary and selects the requested time range.
        """
        if not self.downloaded_files:
            raise ValueError("No files to load. Run fetch() first.")
        
        print(f"📂 Loading {len(self.downloaded_files)} W5E5 files...")
        
        # Group files by variable
        files_by_var = {}
        for fpath in self.downloaded_files:
            # Determine which variable this file contains
            for var in self.cfg.variables:
                if f"/{var}/" in fpath or f"_{self._map_variable_name(var)}_" in fpath:
                    if var not in files_by_var:
                        files_by_var[var] = []
                    files_by_var[var].append(fpath)
                    break
        
        # Load each variable separately and merge
        datasets = []
        for var, file_list in files_by_var.items():
            print(f"  Loading {var} from {len(file_list)} file(s)...")
            
            if len(file_list) == 1:
                ds_var = xr.open_dataset(file_list[0])
            else:
                # Multiple files - open as multi-file dataset
                ds_var = xr.open_mfdataset(
                    file_list,
                    combine='by_coords',
                    parallel=True
                )
            
            datasets.append(ds_var)
        
        # Merge all variables into one dataset
        if len(datasets) == 1:
            self.ds = datasets[0]
        else:
            self.ds = xr.merge(datasets)
        
        # Subset to requested time range
        start = self.cfg.time_range.start_date
        end = self.cfg.time_range.end_date
        self.ds = self.ds.sel(time=slice(start, end))
        
        # Add metadata
        self.ds.attrs['source'] = 'W5E5 via ISIMIP'
        self.ds.attrs['dataset'] = 'W5E5v2.0'
        self.ds.attrs['description'] = 'WFDE5 over land merged with ERA5 over ocean'
        
        print(f"✅ Loaded dataset with {len(self.ds.data_vars)} variables")
    
    def extract(self, *, point: Optional[Tuple[float, float]] = None, 
                box: Optional[Dict] = None, 
                shapefile: Optional[str] = None, 
                buffer_km: float = 0.0):
        """
        Store extraction instructions to be applied during or after load.
        
        Parameters
        ----------
        point : tuple of (lon, lat), optional
            Extract data for a specific point location
        box : dict, optional
            Extract data for a bounding box with keys: lon_min, lon_max, lat_min, lat_max
        shapefile : str or GeoDataFrame, optional
            Extract data for a shapefile region
        buffer_km : float, default=0.0
            Buffer distance in kilometers around point (converted to degrees)
        """
        if point is not None:
            lon, lat = point
            buffer_deg = buffer_km / 111.0
            self._extract_mode = "point"
            self._extract_params = (lon, lat, buffer_deg)
            
            # Apply extraction if dataset is already loaded
            if self.ds is not None:
                self._apply_extraction()
        
        elif box is not None:
            self._extract_mode = "box"
            self._extract_params = box
            
            if self.ds is not None:
                self._apply_extraction()
        
        elif shapefile is not None:
            import geopandas as gpd
            if isinstance(shapefile, str):
                gdf = gpd.read_file(shapefile)
            else:
                gdf = shapefile
            
            self._extract_mode = "shapefile"
            self._extract_params = gdf
            
            if self.ds is not None:
                self._apply_extraction()
    
    def _apply_extraction(self):
        """Apply the stored extraction instructions to the dataset."""
        if self._extract_mode == "point":
            lon, lat, buffer_deg = self._extract_params
            
            if buffer_deg > 0:
                self.ds = self.ds.sel(
                    lon=slice(lon - buffer_deg, lon + buffer_deg),
                    lat=slice(lat - buffer_deg, lat + buffer_deg)
                ).mean(["lat", "lon"])
            else:
                self.ds = self.ds.sel(lon=lon, lat=lat, method="nearest")
        
        elif self._extract_mode == "box":
            box = self._extract_params
            self.ds = self.ds.sel(
                lon=slice(box["lon_min"], box["lon_max"]),
                lat=slice(box["lat_min"], box["lat_max"])
            )
        
        elif self._extract_mode == "shapefile":
            import rioxarray
            from shapely.geometry import mapping
            
            gdf = self._extract_params
            self.ds = self.ds.rio.set_spatial_dims(x_dim="lon", y_dim="lat")
            self.ds = self.ds.rio.write_crs("EPSG:4326", inplace=True)
            
            clipped_list = []
            for geom in gdf.geometry:
                clipped = self.ds.rio.clip([mapping(geom)], gdf.crs, drop=True)
                clipped_list.append(clipped)
            
            self.ds = xr.concat(clipped_list, dim="geom_id")
    
    def save_netcdf(self, filename: str):
        """Save the dataset to a NetCDF file."""
        if self.ds is None:
            raise ValueError("No dataset loaded. Run load() first.")
        
        output_path = Path(filename)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        self.ds.to_netcdf(output_path)
        print(f"💾 Saved to: {output_path}")
    
    def save_csv(self, filename: str):
        """Save the dataset to a CSV file."""
        if self.ds is None:
            raise ValueError("No dataset loaded. Run load() first.")
        
        output_path = Path(filename)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        df = self.ds.to_dataframe()
        df.to_csv(output_path)
        print(f"💾 Saved to: {output_path}")
    
    def _map_variable_name(self, var: str) -> str:
        """
        Map standard variable names to W5E5 variable names.
        
        W5E5 uses standard CMIP variable names.
        """
        # W5E5 uses standard names, so most map directly
        variable_map = {
            'tas': 'tas',        # Near-surface air temperature
            'tasmax': 'tasmax',  # Maximum temperature
            'tasmin': 'tasmin',  # Minimum temperature
            'pr': 'pr',          # Precipitation
            'rsds': 'rsds',      # Surface downwelling shortwave radiation
            'hurs': 'hurs',      # Near-surface relative humidity
            'sfcWind': 'sfcwind', # Near-surface wind speed (note: lowercase 'w')
            'ps': 'ps',          # Surface air pressure
            'rlds': 'rlds',      # Surface downwelling longwave radiation
        }
        
        return variable_map.get(var, var)
    
    def _is_file_in_date_range(self, filename: str, start_date: datetime, end_date: datetime) -> bool:
        """
        Check if a file covers the requested date range.
        
        W5E5 files typically have year ranges in their names like:
        w5e5v2.0_obsclim_tas_global_daily_1979_1989.nc
        """
        import re
        
        # Extract year range from filename
        match = re.search(r'_(\d{4})_(\d{4})\.nc', filename)
        if match:
            file_start_year = int(match.group(1))
            file_end_year = int(match.group(2))
            
            # Check if there's any overlap
            return not (file_end_year < start_date.year or file_start_year > end_date.year)
        
        # If we can't parse the date, include the file to be safe
        return True


class W5E5Mirror(W5E5):
    """
    Alias for W5E5 class to maintain consistent naming with other datasets.
    """
    pass
