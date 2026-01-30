"""
NEX-GDDP-CMIP6 dataset access module

This module provides access to NASA Earth Exchange Global Daily Downscaled Projections (NEX-GDDP-CMIP6) data.
NEX-GDDP-CMIP6 provides downscaled climate projections from CMIP6 at 0.25° resolution (~25km) globally.

Data is accessed via NASA's THREDDS Data Server using the HTTP file server for direct file downloads.
The dataset includes daily climate data for various CMIP6 models and scenarios.

More info: https://www.nccs.nasa.gov/services/data-collections/land-based-products/nex-gddp-cmip6
"""

import os
import xarray as xr
import pandas as pd
from pathlib import Path
from datetime import datetime
from omegaconf import DictConfig
from typing import Optional, Tuple, Dict, List
import warnings
import requests
from tqdm import tqdm

warnings.filterwarnings("ignore", category=Warning)


class NEXGDDP:
    """
    A class to download and process NEX-GDDP-CMIP6 climate data from NASA's THREDDS server.
    
    NEX-GDDP-CMIP6 provides statistically downscaled CMIP6 climate projections at 0.25° resolution.
    Data is available for multiple climate models, scenarios, and variables at daily temporal resolution.
    
    Attributes
    ----------
    cfg : DictConfig
        Configuration containing lat, lon, variables, time_range, experiment_id, source_id, etc.
    ds : xr.Dataset
        Loaded xarray dataset
    experiment_id : str
        CMIP6 experiment identifier (e.g., 'historical', 'ssp126', 'ssp585')
    source_id : str
        CMIP6 model identifier (e.g., 'GFDL-ESM4', 'UKESM1-0-LL', 'MRI-ESM2-0')
    member_id : str
        CMIP6 realization identifier (e.g., 'r1i1p1f1')
    base_url : str
        Base URL for NASA THREDDS HTTP file server
    """
    
    # Available models in NEX-GDDP-CMIP6
    AVAILABLE_MODELS = [
        "ACCESS-CM2",
        "ACCESS-ESM1-5",
        "BCC-CSM2-MR",
        "CESM2",
        "CESM2-WACCM",
        "CMCC-CM2-SR5",
        "CMCC-ESM2",
        "CNRM-CM6-1",
        "CNRM-ESM2-1",
        "CanESM5",
        "EC-Earth3",
        "EC-Earth3-Veg-LR",
        "FGOALS-g3",
        "GFDL-CM4",
        "GFDL-ESM4",
        "GISS-E2-1-G",
        "HadGEM3-GC31-LL",
        "HadGEM3-GC31-MM",
        "IITM-ESM",
        "INM-CM4-8",
        "INM-CM5-0",
        "IPSL-CM6A-LR",
        "KACE-1-0-G",
        "KIOST-ESM",
        "MIROC-ES2L",
        "MIROC6",
        "MPI-ESM1-2-HR",
        "MPI-ESM1-2-LR",
        "MRI-ESM2-0",
        "NESM3",
        "NorESM2-LM",
        "NorESM2-MM",
        "TaiESM1",
        "UKESM1-0-LL",
    ]
    
    # Available experiments
    AVAILABLE_EXPERIMENTS = [
        "historical",
        "ssp126",
        "ssp245",
        "ssp370",
        "ssp585",
    ]
    
    # Available variables
    AVAILABLE_VARIABLES = {
        'tas': 'Near-Surface Air Temperature',
        'tasmax': 'Daily Maximum Near-Surface Air Temperature',
        'tasmin': 'Daily Minimum Near-Surface Air Temperature',
        'pr': 'Precipitation',
        'hurs': 'Near-Surface Relative Humidity',
        'huss': 'Near-Surface Specific Humidity',
        'sfcWind': 'Near-Surface Wind Speed',
        'rsds': 'Surface Downwelling Shortwave Radiation',
        'rlds': 'Surface Downwelling Longwave Radiation',
    }
    
    def __init__(self, cfg: DictConfig):
        """
        Initialize NEX-GDDP data accessor.
        
        Parameters
        ----------
        cfg : DictConfig
            Configuration with required fields:
            - variables: list of variable names
            - time_range: dict with start_date and end_date
            - data_dir: directory to save downloaded files
            Optional fields:
            - experiment_id: experiment/scenario (default: 'historical')
            - source_id: model name (default: 'GFDL-ESM4')
            - member_id: realization (default: 'r1i1p1f1')
        """
        self.cfg = cfg
        self.ds = None
        self.downloaded_files = []
        self._extract_mode = None
        self._extract_params = None
        
        # Extract NEX-GDDP-specific parameters
        self.experiment_id = cfg.get('experiment_id', 'historical')
        self.source_id = cfg.get('source_id', 'GFDL-ESM4')
        self.member_id = cfg.get('member_id', 'r1i1p1f1')
        
        # Base URL for NASA THREDDS HTTP File Server
        self.base_url = "https://ds.nccs.nasa.gov/thredds/fileServer/AMES/NEX/GDDP-CMIP6"
        
        # Validate inputs
        self._validate_inputs()
        self._validate_time_range()
    
    def _validate_inputs(self):
        """Validate model, experiment, and variable selections."""
        # Normalize model name to uppercase with hyphens
        self.source_id = self.source_id.upper().replace('_', '-')
        
        if self.source_id not in self.AVAILABLE_MODELS:
            print(f"⚠️  Warning: Model '{self.source_id}' may not be available.")
            print(f"   Available models: {', '.join(self.AVAILABLE_MODELS[:5])}...")
        
        if self.experiment_id not in self.AVAILABLE_EXPERIMENTS:
            print(f"⚠️  Warning: Experiment '{self.experiment_id}' may not be available.")
            print(f"   Available experiments: {', '.join(self.AVAILABLE_EXPERIMENTS)}")
        
        for var in self.cfg.variables:
            if var not in self.AVAILABLE_VARIABLES:
                print(f"⚠️  Warning: Variable '{var}' may not be available.")
                print(f"   Available variables: {', '.join(self.AVAILABLE_VARIABLES.keys())}")
    
    def _validate_time_range(self):
        """
        Validate that the requested time range is appropriate for the experiment.
        
        Historical runs: 1950-2014
        SSP scenarios: 2015-2100
        
        Raises
        ------
        ValueError
            If the time range doesn't match the experiment period
        """
        start_date = datetime.fromisoformat(self.cfg.time_range.start_date)
        end_date = datetime.fromisoformat(self.cfg.time_range.end_date)
        
        start_year = start_date.year
        end_year = end_date.year
        
        # Define valid periods for each experiment type
        if self.experiment_id == 'historical':
            valid_start = 1950
            valid_end = 2014
            period_name = "Historical"
        elif self.experiment_id.startswith('ssp'):
            valid_start = 2015
            valid_end = 2100
            period_name = f"SSP scenario ({self.experiment_id})"
        else:
            # Unknown experiment, skip validation
            return
        
        # Check if requested period is outside valid range
        if end_year < valid_start or start_year > valid_end:
            raise ValueError(
                f"❌ Time range mismatch for experiment '{self.experiment_id}'!\n"
                f"   Requested: {start_year}-{end_year}\n"
                f"   Valid period for {period_name}: {valid_start}-{valid_end}\n"
                f"   \n"
                f"   Hint: Use 'historical' for years 1950-2014, and SSP scenarios for 2015-2100."
            )
        
        # Warn if requested period extends beyond valid range
        if start_year < valid_start or end_year > valid_end:
            print(f"⚠️  Warning: Requested time range {start_year}-{end_year} extends beyond")
            print(f"   the typical {period_name} period ({valid_start}-{valid_end}).")
            print(f"   Data availability may be limited.")
    
    def get_experiment_ids(self) -> List[str]:
        """
        Get available NEX-GDDP-CMIP6 experiment IDs.
        
        Returns
        -------
        List[str]
            List of available experiment IDs
        """
        return self.AVAILABLE_EXPERIMENTS.copy()
    
    def get_source_ids(self, experiment_id: Optional[str] = None) -> List[str]:
        """
        Get available NEX-GDDP-CMIP6 model (source) IDs.
        
        Parameters
        ----------
        experiment_id : str, optional
            Experiment ID (not used for NEX-GDDP as models are consistent across experiments)
            
        Returns
        -------
        List[str]
            List of available model IDs
        """
        return self.AVAILABLE_MODELS.copy()
    
    def get_variables(self) -> Dict[str, str]:
        """
        Get available variables with descriptions.
        
        Returns
        -------
        Dict[str, str]
            Dictionary mapping variable names to descriptions
        """
        return self.AVAILABLE_VARIABLES.copy()
    
    def _construct_download_url(self, variable: str, year: int) -> str:
        """
        Construct THREDDS HTTP file server URL for downloading NEX-GDDP data.
        
        Parameters
        ----------
        variable : str
            Variable name (e.g., 'tasmax', 'pr')
        year : int
            Year to download
            
        Returns
        -------
        tuple
            (url, filename) - Complete HTTP URL and filename
        """
        # Filename format: {var}_day_{model}_{experiment}_{member}_gn_{year}_v2.0.nc
        filename = f"{variable}_day_{self.source_id}_{self.experiment_id}_{self.member_id}_gn_{year}_v2.0.nc"
        
        # Construct full URL using HTTP file server
        url = f"{self.base_url}/{self.source_id}/{self.experiment_id}/{self.member_id}/{variable}/{filename}"
        
        return url, filename
    
    def fetch(self):
        """
        Download NEX-GDDP-CMIP6 files from NASA THREDDS server
        for the requested variables, time range, experiment, and model.
        """
        print(f"🔍 Downloading NEX-GDDP-CMIP6 data from NASA THREDDS...")
        print(f"   Model: {self.source_id}, Experiment: {self.experiment_id}")
        
        start_date = datetime.fromisoformat(self.cfg.time_range.start_date)
        end_date = datetime.fromisoformat(self.cfg.time_range.end_date)
        
        # Create directory structure: nexgddp/{MODEL}/{experiment}/{variable}/
        base_dir = Path(self.cfg.data_dir) / "nexgddp" / self.source_id / self.experiment_id
        base_dir.mkdir(parents=True, exist_ok=True)
        
        # Generate list of years to download
        years = range(start_date.year, end_date.year + 1)
        
        # Download each variable for each year
        for var in self.cfg.variables:
            print(f"\n📥 Fetching {var} ({self.AVAILABLE_VARIABLES.get(var, var)})...")
            
            var_dir = base_dir / var
            var_dir.mkdir(parents=True, exist_ok=True)
            
            for year in tqdm(list(years), desc=f"  Downloading {var}"):
                url, filename = self._construct_download_url(var, year)
                local_path = var_dir / filename
                
                # Skip if file already exists
                if local_path.exists():
                    # print(f"  ✓ Already exists: {filename}")
                    self.downloaded_files.append(str(local_path))
                    continue
                
                # Download file
                try:
                    response = requests.get(url, stream=True, timeout=120)
                    response.raise_for_status()
                    
                    # Save file
                    with open(local_path, 'wb') as f:
                        for chunk in response.iter_content(chunk_size=8192):
                            f.write(chunk)
                    
                    self.downloaded_files.append(str(local_path))
                    # print(f"  ✓ Downloaded: {filename}")
                    
                except requests.exceptions.RequestException as e:
                    print(f"  ❌ Failed to download {filename}: {str(e)}")
                    continue
        
        print(f"\n✅ Downloaded {len(self.downloaded_files)} files")
    
    def load(self):
        """
        Load the downloaded NEX-GDDP netCDF files into an xarray Dataset.
        Combines multiple files if necessary and selects the requested time range.
        """
        if not self.downloaded_files:
            raise ValueError("No files to load. Run fetch() first.")
        
        print(f"📂 Loading {len(self.downloaded_files)} NEX-GDDP files...")
        
        # Group files by variable
        files_by_var = {}
        for fpath in self.downloaded_files:
            # Determine which variable this file contains
            for var in self.cfg.variables:
                if f"/{var}/" in fpath or f"_{var}_day_" in fpath:
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
                    sorted(file_list),
                    combine='by_coords',
                    parallel=True,
                    engine='netcdf4'
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
        self.ds.attrs['source'] = f'NEX-GDDP-CMIP6'
        self.ds.attrs['source_url'] = 'https://www.nccs.nasa.gov/services/data-collections/land-based-products/nex-gddp-cmip6'
        self.ds.attrs['experiment_id'] = self.experiment_id
        self.ds.attrs['source_id'] = self.source_id
        self.ds.attrs['member_id'] = self.member_id
        self.ds.attrs['resolution'] = '0.25 degrees'
        self.ds.attrs['description'] = f'NEX-GDDP-CMIP6 downscaled {self.experiment_id} from {self.source_id}'
        
        print(f"✅ Loaded dataset with {len(self.ds.data_vars)} variables")
        
        # Apply extraction if it was set before loading
        if self._extract_mode is not None:
            self._apply_extraction()
    
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
            
            # Find coordinate names
            lat_name = None
            lon_name = None
            for coord in self.ds.coords:
                if coord.lower() in ['lat', 'latitude']:
                    lat_name = coord
                elif coord.lower() in ['lon', 'longitude']:
                    lon_name = coord
            
            if lat_name is None or lon_name is None:
                raise ValueError("Could not find latitude/longitude coordinates in dataset")
            
            if buffer_deg > 0:
                self.ds = self.ds.sel(
                    {lon_name: slice(lon - buffer_deg, lon + buffer_deg),
                     lat_name: slice(lat - buffer_deg, lat + buffer_deg)}
                ).mean([lat_name, lon_name])
            else:
                self.ds = self.ds.sel({lon_name: lon, lat_name: lat}, method="nearest")
        
        elif self._extract_mode == "box":
            box = self._extract_params
            
            # Find coordinate names
            lat_name = None
            lon_name = None
            for coord in self.ds.coords:
                if coord.lower() in ['lat', 'latitude']:
                    lat_name = coord
                elif coord.lower() in ['lon', 'longitude']:
                    lon_name = coord
            
            if lat_name is None or lon_name is None:
                raise ValueError("Could not find latitude/longitude coordinates in dataset")
            
            self.ds = self.ds.sel(
                {lon_name: slice(box["lon_min"], box["lon_max"]),
                 lat_name: slice(box["lat_min"], box["lat_max"])}
            )
        
        elif self._extract_mode == "shapefile":
            import rioxarray
            from shapely.geometry import mapping
            
            gdf = self._extract_params
            
            # Find coordinate names
            lat_name = None
            lon_name = None
            for coord in self.ds.coords:
                if coord.lower() in ['lat', 'latitude']:
                    lat_name = coord
                elif coord.lower() in ['lon', 'longitude']:
                    lon_name = coord
            
            if lat_name is None or lon_name is None:
                raise ValueError("Could not find latitude/longitude coordinates in dataset")
            
            self.ds = self.ds.rio.set_spatial_dims(x_dim=lon_name, y_dim=lat_name)
            self.ds = self.ds.rio.write_crs("EPSG:4326", inplace=True)
            
            clipped_list = []
            for geom in gdf.geometry:
                clipped = self.ds.rio.clip([mapping(geom)], gdf.crs, drop=True)
                clipped_list.append(clipped)
            
            self.ds = xr.concat(clipped_list, dim="geom_id")
    
    def save_netcdf(self, filename: str):
        """
        Save the dataset to a NetCDF file.
        
        Parameters
        ----------
        filename : str
            Output filename or path
        """
        if self.ds is None:
            raise ValueError("No dataset loaded. Run load() first.")
        
        output_path = Path(filename)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        self.ds.to_netcdf(output_path)
        print(f"💾 Saved to: {output_path}")
    
    def save_csv(self, filename: str):
        """
        Save the dataset to a CSV file.
        
        Parameters
        ----------
        filename : str
            Output filename or path
        """
        if self.ds is None:
            raise ValueError("No dataset loaded. Run load() first.")
        
        output_path = Path(filename)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        df = self.ds.to_dataframe()
        df.to_csv(output_path)
        print(f"💾 Saved to: {output_path}")
    
    def to_dataframe(self) -> pd.DataFrame:
        """
        Convert the dataset to a pandas DataFrame.
        
        Returns
        -------
        pd.DataFrame
            Dataset as DataFrame
        """
        if self.ds is None:
            raise ValueError("No dataset loaded. Run load() first.")
        
        return self.ds.to_dataframe()


# Alias for consistency
NEXGDDPMirror = NEXGDDP
