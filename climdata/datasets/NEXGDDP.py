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
import time

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
        self.catalog_base_url = "https://ds.nccs.nasa.gov/thredds/catalog/AMES/NEX/GDDP-CMIP6"
        
        # Validate inputs
        self._validate_inputs()
        self._validate_time_range()
        
        # Auto-discover member_id and grid_label if not explicitly set
        self.grid_label = cfg.get('grid_label', 'gn')
        self._auto_discover_metadata()
    
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
    
    def _auto_discover_metadata(self):
        """
        Auto-discover member_id and grid_label by querying THREDDS catalog.
        If member_id was explicitly provided in config, use that instead.
        """
        # Only auto-discover if member_id was not explicitly set
        if 'member_id' in self.cfg:
            print(f"ℹ️  Using configured member_id: {self.member_id}")
            return
        
        print(f"🔍 Auto-discovering metadata for {self.source_id}/{self.experiment_id}...")
        
        try:
            # Get available member_ids from THREDDS catalog
            member_ids = self._get_available_member_ids()
            
            if not member_ids:
                print(f"⚠️  Could not auto-discover member_ids, using default: {self.member_id}")
                return
            
            # Use the first available member_id
            discovered_member_id = member_ids[0]
            
            # Get grid label from a sample file
            grid_label = self._get_grid_label(discovered_member_id)
            
            if grid_label:
                self.grid_label = grid_label
                print(f"✓ Discovered grid_label: {self.grid_label}")
            
            if discovered_member_id != self.member_id:
                self.member_id = discovered_member_id
                print(f"✓ Discovered member_id: {self.member_id}")
            
        except Exception as e:
            print(f"⚠️  Auto-discovery failed: {str(e)}")
            print(f"   Using defaults: member_id={self.member_id}, grid_label={self.grid_label}")
    
    def _get_available_member_ids(self) -> List[str]:
        """
        Get available member_ids (realizations) for the current model and experiment
        by probing the server with common member_id patterns.
        
        Returns
        -------
        List[str]
            List of available member_ids (e.g., ['r1i1p1f1', 'r2i1p1f1'])
        """
        # Common member_id patterns to check
        common_members = [
            'r1i1p1f1', 'r1i1p1f2', 'r1i1p1f3',
            'r2i1p1f1', 'r2i1p1f2', 'r3i1p1f1',
            'r4i1p1f1', 'r5i1p1f1'
        ]
        
        available_members = []
        
        # Use the first variable to check availability
        if not self.cfg.variables:
            return []
        
        test_var = self.cfg.variables[0]
        
        # Determine a test year based on experiment
        if self.experiment_id == 'historical':
            test_year = 2000  # Use a year that's more likely to exist
        else:
            test_year = 2050
        
        print(f"   Checking available realizations...")
        
        for member_id in common_members:
            # Try common grid labels (gr1 is used by many models)
            for grid_label in ['gr1', 'gn', 'gr']:
                url = self._construct_test_url(test_var, test_year, member_id, grid_label)
                
                try:
                    # Use HEAD request to check if file exists without downloading
                    response = requests.head(url, timeout=10, allow_redirects=True)
                    
                    if response.status_code == 200:
                        available_members.append(member_id)
                        print(f"   ✓ Found: {member_id} (grid: {grid_label})")
                        break  # Found this member_id, move to next
                    
                except requests.exceptions.RequestException:
                    continue
        
        return sorted(set(available_members))
    
    def _construct_test_url(self, variable: str, year: int, member_id: str, grid_label: str) -> str:
        """
        Construct a test URL for probing file availability.
        
        Parameters
        ----------
        variable : str
            Variable name
        year : int
            Year
        member_id : str
            Member/realization ID
        grid_label : str
            Grid label (e.g., 'gn', 'gr')
            
        Returns
        -------
        str
            Complete URL for testing
        """
        filename = f"{variable}_day_{self.source_id}_{self.experiment_id}_{member_id}_{grid_label}_{year}_v2.0.nc"
        url = f"{self.base_url}/{self.source_id}/{self.experiment_id}/{member_id}/{variable}/{filename}"
        return url
    
    def _get_grid_label(self, member_id: str) -> Optional[str]:
        """
        Get grid label by probing with common grid labels.
        
        Parameters
        ----------
        member_id : str
            The member_id to check
            
        Returns
        -------
        str or None
            Grid label (e.g., 'gn', 'gr') or None if not found
        """
        # Use the first available variable
        if not self.cfg.variables:
            return None
        
        sample_var = self.cfg.variables[0]
        
        # Determine a test year based on experiment
        if self.experiment_id == 'historical':
            test_year = 2000  # Use a year that's more likely to exist
        else:
            test_year = 2050
        
        # Try common grid labels (gr1 is common for NEX-GDDP)
        for grid_label in ['gr1', 'gn', 'gr']:
            url = self._construct_test_url(sample_var, test_year, member_id, grid_label)
            
            try:
                response = requests.head(url, timeout=10, allow_redirects=True)
                if response.status_code == 200:
                    return grid_label
            except requests.exceptions.RequestException:
                continue
        
        return None
    
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
    
    def get_member_ids(self) -> List[str]:
        """
        Get available member IDs (realizations) for the current model and experiment
        by querying the THREDDS server.
        
        Returns
        -------
        List[str]
            List of available member_ids (e.g., ['r1i1p1f1', 'r2i1p1f1'])
        """
        return self._get_available_member_ids()
    
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
        # Filename format: {var}_day_{model}_{experiment}_{member}_{grid}_{year}_v2.0.nc
        filename = f"{variable}_day_{self.source_id}_{self.experiment_id}_{self.member_id}_{self.grid_label}_{year}_v2.0.nc"
        
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
                
                # Skip if file already exists and is complete
                if local_path.exists() and self._verify_file_complete(local_path, url):
                    self.downloaded_files.append(str(local_path))
                    continue
                
                # Download file with retry logic
                success = self._download_with_retry(url, local_path, max_retries=5)
                
                if success:
                    self.downloaded_files.append(str(local_path))
                else:
                    print(f"  ⚠️  Failed to download {filename} after multiple retries")
        
        print(f"\n✅ Downloaded {len(self.downloaded_files)} files")
    
    def _verify_file_complete(self, local_path: Path, url: str) -> bool:
        """
        Verify if a local file is complete by comparing size with server.
        
        Parameters
        ----------
        local_path : Path
            Path to local file
        url : str
            URL of the file on server
            
        Returns
        -------
        bool
            True if file is complete, False otherwise
        """
        try:
            local_size = local_path.stat().st_size
            
            # Get expected size from server
            response = requests.head(url, timeout=30, allow_redirects=True)
            if response.status_code == 200:
                expected_size = int(response.headers.get('Content-Length', 0))
                if expected_size > 0 and local_size == expected_size:
                    return True
            
            return False
        except Exception:
            return False
    
    def _download_with_retry(self, url: str, local_path: Path, max_retries: int = 5) -> bool:
        """
        Download a file with retry logic and resume capability.
        
        Parameters
        ----------
        url : str
            URL to download from
        local_path : Path
            Local path to save to
        max_retries : int
            Maximum number of retry attempts
            
        Returns
        -------
        bool
            True if download succeeded, False otherwise
        """
        for attempt in range(max_retries):
            try:
                # Get expected file size
                head_response = requests.head(url, timeout=30, allow_redirects=True)
                head_response.raise_for_status()
                expected_size = int(head_response.headers.get('Content-Length', 0))
                
                # Check if we have a partial download
                existing_size = 0
                if local_path.exists():
                    existing_size = local_path.stat().st_size
                    
                    # If file is complete, we're done
                    if existing_size == expected_size:
                        return True
                    
                    # If partial, try to resume
                    if existing_size > 0 and existing_size < expected_size:
                        headers = {'Range': f'bytes={existing_size}-'}
                        mode = 'ab'  # Append mode
                    else:
                        # File exists but wrong size, start over
                        local_path.unlink()
                        headers = {}
                        mode = 'wb'
                        existing_size = 0
                else:
                    headers = {}
                    mode = 'wb'
                
                # Download with streaming
                response = requests.get(
                    url, 
                    headers=headers,
                    stream=True, 
                    timeout=120,
                    allow_redirects=True
                )
                response.raise_for_status()
                
                # Save file
                with open(local_path, mode) as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                
                # Verify download completed
                final_size = local_path.stat().st_size
                if final_size == expected_size or (existing_size > 0 and final_size == expected_size):
                    return True
                else:
                    if attempt < max_retries - 1:
                        print(f"  ↻ Incomplete download (attempt {attempt + 1}/{max_retries}), retrying...")
                        time.sleep(2 ** attempt)  # Exponential backoff
                        continue
                    else:
                        return False
                
            except (requests.exceptions.RequestException, IOError) as e:
                if attempt < max_retries - 1:
                    wait_time = 2 ** attempt
                    print(f"  ↻ Download error (attempt {attempt + 1}/{max_retries}): {type(e).__name__}")
                    print(f"     Retrying in {wait_time} seconds...")
                    time.sleep(wait_time)
                    continue
                else:
                    print(f"  ❌ Failed after {max_retries} attempts: {str(e)}")
                    # Clean up partial file on final failure
                    if local_path.exists():
                        local_path.unlink()
                    return False
        
        return False
    
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
