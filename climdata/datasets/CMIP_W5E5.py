"""CMIP6 projections bias-adjusted to W5E5, via the ISIMIP client.

ISIMIP3b publishes CMIP6 model output already bias-adjusted and downscaled to
the W5E5 0.5 degree grid. That makes it directly comparable to the W5E5
reference - see :class:`~climdata.datasets.W5E5.W5E5` - without a
bias-correction step of your own, which is the main reason to prefer it over raw
CMIP6 from :class:`~climdata.datasets.CMIPCloud.CMIPCloud`.

ISIMIP names models in lowercase (``gfdl-esm4``), unlike the CMIP6 catalogue's
uppercase (``GFDL-ESM4``).

Example:
    >>> cmip = CMIPW5E5(cfg)                                  # doctest: +SKIP
    >>> cmip.fetch(); cmip.load()                             # doctest: +SKIP
    >>> cmip.extract(point=(13.4, 52.5))                      # doctest: +SKIP
"""

import xarray as xr
from pathlib import Path
from datetime import datetime
from omegaconf import DictConfig
from typing import Optional, Tuple, Dict, List
import warnings
import re

warnings.filterwarnings("ignore", category=Warning)

class CMIPW5E5:
    """Download and assemble ISIMIP3b bias-adjusted CMIP6 projections.

    The lifecycle is ``fetch`` -> ``load`` -> ``extract``. :meth:`extract` may be
    called before :meth:`load`, in which case the subset is stored and applied as
    soon as data arrives.

    Attributes:
        cfg (DictConfig): Configuration supplying ``variables``, ``time_range``,
            ``data_dir`` and ``dataset``.
        ds (xr.Dataset | None): The loaded dataset, set by :meth:`load`.
        client (ISIMIPClient): ISIMIP API client, created in :meth:`__init__`.
        experiment_id (str): Experiment, e.g. ``"historical"``, ``"ssp585"``.
        source_id (str): Model, lowercase, e.g. ``"gfdl-esm4"``.
        member_id (str): Realisation, e.g. ``"r1i1p1f1"``.
        downloaded_files (list[str]): Local paths accumulated by :meth:`fetch`.

    Example:
        >>> cmip = CMIPW5E5(cfg)                                  # doctest: +SKIP
        >>> cmip.get_source_ids("ssp585")                         # doctest: +SKIP
        ['gfdl-esm4', 'ipsl-cm6a-lr', 'mpi-esm1-2-hr', ...]
    """

    def __init__(self, cfg: DictConfig):
        """Bind a configuration, open an ISIMIP client and validate the period.

        Args:
            cfg (DictConfig): Configuration with ``variables``, ``time_range``,
                ``data_dir``, ``dataset``, and optionally ``experiment_id``
                (default ``"historical"``), ``source_id`` (default
                ``"gfdl-esm4"``) and ``member_id`` (default ``"r1i1p1f1"``).

        Raises:
            ImportError: If ``isimip-client`` is not installed.
            ValueError: If the time range cannot belong to the experiment - see
                :meth:`_validate_time_range`.
        """
        self.cfg = cfg
        self.ds = None
        self.client = None
        self.downloaded_files = []
        self._extract_mode = None
        self._extract_params = None
        
        # Extract CMIP6-specific parameters
        self.experiment_id = cfg.get('experiment_id', 'historical')
        self.source_id = cfg.get('source_id', 'gfdl-esm4')
        self.member_id = cfg.get('member_id', 'r1i1p1f1')
        
        # Initialize ISIMIP client
        try:
            from isimip_client.client import ISIMIPClient
            self.client = ISIMIPClient()
        except ImportError:
            raise ImportError(
                "isimip-client is required for CMIP-W5E5 data access. "
                "Install it with: pip install isimip-client"
            )
        
        # Validate time range for experiment
        self._validate_time_range()
    
    def _fix_coords(self, ds):
        """Put latitude in ascending order and longitude on the 0-360 convention.

        Args:
            ds (xr.Dataset | xr.DataArray): Object with CF-identifiable latitude
                and longitude coordinates.

        Returns:
            xr.Dataset | xr.DataArray: Same type as the input, reoriented.
        """
        ds = ds.cf.sortby("latitude")
        lon_name = ds.cf["longitude"].name
        ds = ds.assign_coords({lon_name: ds.cf["longitude"] % 360})
        return ds.sortby(lon_name)
    
    def _validate_time_range(self):
        """Validate that the requested time range is appropriate for the experiment.

        Historical runs: 1850-2014
        SSP scenarios: 2015-2100

        Raises:
            ValueError: If the time range doesn't match the experiment period
        """
        start_date = datetime.fromisoformat(self.cfg.time_range.start_date)
        end_date = datetime.fromisoformat(self.cfg.time_range.end_date)
        
        start_year = start_date.year
        end_year = end_date.year
        
        # Define valid periods for each experiment type
        if self.experiment_id == 'historical':
            valid_start = 1850
            valid_end = 2014
            period_name = "Historical"
        elif self.experiment_id.startswith('ssp'):
            valid_start = 2015
            valid_end = 2100
            period_name = f"SSP scenario ({self.experiment_id})"
        elif self.experiment_id == 'picontrol':
            # Pre-industrial control - typically long runs, less strict
            return
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
                f"   Hint: Use 'historical' for years 1850-2014, and SSP scenarios (ssp126, ssp370, ssp585) for 2015-2100."
            )
        
        # Warn if requested period extends beyond valid range
        if start_year < valid_start or end_year > valid_end:
            print(f"⚠️  Warning: Requested time range {start_year}-{end_year} extends beyond")
            print(f"   the typical {period_name} period ({valid_start}-{valid_end}).")
            print(f"   Data availability may be limited.")
    
    def get_experiment_ids(self) -> List[str]:
        """Get available CMIP6 experiment IDs from ISIMIP repository.

        Returns:
            List[str]: List of available experiment IDs (e.g., ['historical', 'ssp126',
                'ssp245', 'ssp370', 'ssp585'])
        """
        print("🔍 Fetching available experiment IDs from ISIMIP...")
        
        try:
            # Query ISIMIP3b for CMIP6 scenarios
            response = self.client.datasets(
                simulation_round='ISIMIP3b',
                product='InputData'
            )

            results = self._normalize_dataset_results(response)
            if not results:
                print("⚠️ No datasets found")
                return []

            # Extract unique experiment IDs from results
            experiment_ids = set()
            for dataset in results:
                # Parse experiment ID from dataset metadata or path
                climate_scenario = dataset.get('climate_scenario', '')
                if climate_scenario and climate_scenario not in ['obsclim', 'counterclim']:
                    experiment_ids.add(climate_scenario)
            
            # Common CMIP6 experiments in ISIMIP3b
            common_experiments = ['historical', 'ssp126', 'ssp370', 'ssp585']
            available = sorted([exp for exp in common_experiments if exp in experiment_ids or exp == 'historical'])
            
            print(f"✅ Found {len(available)} experiment IDs: {available}")
            return available
            
        except Exception as e:
            print(f"❌ Error fetching experiment IDs: {str(e)}")
            # Return common CMIP6 experiments as fallback
            return ['historical', 'ssp126', 'ssp370', 'ssp585']
    
    def get_source_ids(self, experiment_id: Optional[str] = None) -> List[str]:
        """Get available CMIP6 model (source) IDs for a given experiment from ISIMIP repository.

        Args:
            experiment_id (str, optional): CMIP6 experiment ID. If None, uses
                self.experiment_id

        Returns:
            List[str]: List of available model IDs (e.g., ['gfdl-esm4', 'ipsl-cm6a-lr',
                'mpi-esm1-2-hr', 'mri-esm2-0', 'ukesm1-0-ll'])
        """
        if experiment_id is None:
            experiment_id = self.experiment_id
        
        print(f"🔍 Fetching available source IDs for experiment '{experiment_id}'...")
        
        try:
            # Map experiment_id to climate_scenario for ISIMIP
            climate_scenario = self._map_experiment_to_scenario(experiment_id)
            
            # Query ISIMIP3b for available models
            response = self.client.datasets(
                simulation_round='ISIMIP3b',
                product='InputData',
                climate_scenario=climate_scenario
            )
            
            results = self._normalize_dataset_results(response)
            if not results:
                print(f"⚠️ No datasets found for experiment '{experiment_id}'")
                return []
            
            # Extract unique source IDs from results
            source_ids = set()
            for dataset in results:
                # Parse source ID from dataset metadata or filename
                climate_forcing = dataset.get('climate_forcing', '')
                if climate_forcing:
                    # ISIMIP format typically: gfdl-esm4, ipsl-cm6a-lr, etc.
                    # Extract model name from climate forcing string
                    model_match = re.search(r'(gfdl-esm4|ipsl-cm6a-lr|mpi-esm1-2-hr|mri-esm2-0|ukesm1-0-ll)', 
                                          climate_forcing, re.IGNORECASE)
                    if model_match:
                        source_ids.add(model_match.group(1).lower())
            
            # Common ISIMIP3b CMIP6 models
            common_models = ['gfdl-esm4', 'ipsl-cm6a-lr', 'mpi-esm1-2-hr', 'mri-esm2-0', 'ukesm1-0-ll']
            available = sorted([model for model in common_models if model in source_ids])
            
            if not available:
                available = common_models  # Fallback to common models
            
            print(f"✅ Found {len(available)} source IDs: {available}")
            return available
            
        except Exception as e:
            print(f"❌ Error fetching source IDs: {str(e)}")
            # Return common ISIMIP3b models as fallback
            return ['gfdl-esm4', 'ipsl-cm6a-lr', 'mpi-esm1-2-hr', 'mri-esm2-0', 'ukesm1-0-ll']
    
    def _map_experiment_to_scenario(self, experiment_id: str) -> str:
        """Map CMIP6 experiment ID to ISIMIP climate scenario name.

        Args:
            experiment_id (str): CMIP6 experiment ID

        Returns:
            str: ISIMIP climate scenario name
        """
        scenario_map = {
            'historical': 'historical',
            'ssp126': 'ssp126',
            'ssp245': 'ssp245',
            'ssp370': 'ssp370',
            'ssp585': 'ssp585',
            'picontrol': 'picontrol',
            'ssp119': 'ssp119',
            'ssp434': 'ssp434',
            'ssp460': 'ssp460',
        }
        return scenario_map.get(experiment_id, experiment_id)
    
    def fetch(self):
        """Search ISIMIP3b and download the files covering the request.

        Queries once per configured variable for the configured model and
        experiment, keeps the files whose filename year-range overlaps the period
        (see :meth:`_is_file_in_date_range`), and downloads any not already under
        ``<data_dir>/<DATASET>/<variable>/``.

        Per-variable failures are reported and skipped rather than raised, so one
        unavailable variable does not abort a multi-variable request - check
        :attr:`downloaded_files` to see what actually arrived.

        Returns:
            None: Paths are appended to :attr:`downloaded_files`.
        """
        print(f"🔍 Searching for CMIP6 datasets in ISIMIP repository...")
        print(f"   Model: {self.source_id}, Experiment: {self.experiment_id}")
        
        start_date = datetime.fromisoformat(self.cfg.time_range.start_date)
        end_date = datetime.fromisoformat(self.cfg.time_range.end_date)
        
        # Create directory structure: cmip_w5e5/global/daily/{experiment_id}/{MODEL}/
        base_dir = Path(self.cfg.data_dir) / "cmip_w5e5" / "global" / "daily" / self.experiment_id / self.source_id.upper()
        base_dir.mkdir(parents=True, exist_ok=True)
        
        # Map experiment_id to ISIMIP climate scenario
        climate_scenario = self._map_experiment_to_scenario(self.experiment_id)
        
        # Search for each variable separately
        for var in self.cfg.variables:
            print(f"\n📥 Fetching {var}...")
            
            # Map variable names to CMIP6/W5E5 names if needed
            cmip_var = self._map_variable_name(var)
            
            # Search ISIMIP repository for CMIP6 data
            # CMIP6 is available in ISIMIP3b as secondary input data
            try:
                response = self.client.datasets(
                    simulation_round='ISIMIP3b',
                    product='InputData',
                    climate_forcing=self.source_id.lower(),  # CMIP6 model
                    climate_scenario=climate_scenario,
                    climate_variable=cmip_var
                )
                
                results = self._normalize_dataset_results(response)
                if not results:
                    print(f"⚠️ No CMIP6 datasets found for {var} with model {self.source_id}")
                    continue
                
                # Get the first matching dataset
                dataset = results[0]
                print(f"✅ Found dataset: {dataset.get('name', 'unnamed')}")
                
                # Filter files by date range
                for file_info in self._extract_files_list(dataset):
                    file_name = file_info['name']
                    
                    # Parse date from filename
                    # Example: gfdl-esm4_r1i1p1f1_ssp585_tas_global_daily_2015_2024.nc
                    if self._is_file_in_date_range(file_name, start_date, end_date):
                        # Add variable to path: cmip_w5e5/global/daily/{experiment_id}/{MODEL}/{variable}/
                        var_dir = base_dir / var
                        var_dir.mkdir(parents=True, exist_ok=True)
                        local_path = var_dir / file_name
                        
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

    def _normalize_dataset_results(self, response) -> List[Dict]:
        """Coerce an isimip-client response into a list of dataset records.

        The client returns a bare list from some endpoints and a paginated
        ``{"results": [...]}`` envelope from others.

        Args:
            response: The value returned by ``client.datasets(...)``.

        Returns:
            list[dict]: Dataset records, empty if the response has neither shape.
        """
        if response is None:
            return []
        if isinstance(response, list):
            return response
        if isinstance(response, dict):
            results = response.get('results')
            if isinstance(results, list):
                return results
        return []

    def _extract_files_list(self, dataset: Dict) -> List[Dict]:
        """Pull the file records out of one ISIMIP dataset record.

        Args:
            dataset (dict): A single dataset record.

        Returns:
            list[dict]: File records, each with ``path``, ``name`` and
            ``file_url``. Empty if the record carries neither ``files`` nor a
            single ``file``.
        """
        files = dataset.get('files')
        if isinstance(files, list):
            return files
        single_file = dataset.get('file')
        if isinstance(single_file, dict):
            return [single_file]
        return []
    
    def load(self):
        """Open the downloaded files and merge them into one dataset.

        Files are grouped by variable, concatenated along time where a variable
        spans several decades, merged across variables, and clipped to the
        configured time range.

        Returns:
            None: The dataset is stored on :attr:`ds`.

        Raises:
            ValueError: If :meth:`fetch` has not run, or produced no files.
        """
        if not self.downloaded_files:
            raise ValueError("No files to load. Run fetch() first.")
        
        print(f"📂 Loading {len(self.downloaded_files)} CMIP6 files...")
        
        # Group files by variable
        files_by_var = {}
        for fpath in self.downloaded_files:
            # Determine which variable this file contains
            for var in self.cfg.variables:
                cmip_var = self._map_variable_name(var)
                if f"/{var}/" in fpath or f"_{cmip_var}_" in fpath:
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
            self.ds = self._fix_coords(datasets[0])
        else:
            self.ds = self._fix_coords(xr.merge(datasets))
        
        # Subset to requested time range
        start = self.cfg.time_range.start_date
        end = self.cfg.time_range.end_date
        self.ds = self.ds.sel(time=slice(start, end))
        
        # Add metadata
        self.ds.attrs['source'] = f'CMIP6 {self.source_id} via ISIMIP'
        self.ds.attrs['dataset'] = f'CMIP6-W5E5'
        self.ds.attrs['experiment_id'] = self.experiment_id
        self.ds.attrs['source_id'] = self.source_id
        self.ds.attrs['member_id'] = self.member_id
        self.ds.attrs['description'] = f'CMIP6 {self.experiment_id} scenario from {self.source_id} in W5E5 format'
        
        print(f"✅ Loaded dataset with {len(self.ds.data_vars)} variables")
    
    def extract(self, *, point: Optional[Tuple[float, float]] = None, 
                box: Optional[Dict] = None, 
                shapefile: Optional[str] = None, 
                buffer_km: float = 0.0):
        """Store extraction instructions to be applied during or after load.

        Args:
            point (tuple of (lon, lat), optional): Extract data for a specific point
                location
            box (dict, optional): Extract data for a bounding box with keys: lon_min,
                lon_max, lat_min, lat_max
            shapefile (str or GeoDataFrame, optional): Extract data for a shapefile
                region
            buffer_km (float, default=0.0): Buffer distance in kilometers around point
                (converted to degrees)
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
        """Run the subset recorded by :meth:`extract` against :attr:`ds`.

        Point requests either take the nearest cell or area-average a buffered
        box; box requests slice; shapefile requests clip each geometry through
        rioxarray and stack the results along a new ``geom_id`` dimension.

        Returns:
            None: :attr:`ds` is replaced in place.
        """
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
            import rioxarray  # noqa: F401 — registers the .rio accessor
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
        """Write the dataset to NetCDF, creating parent directories.

        Args:
            filename (str): Destination path.

        Returns:
            None

        Raises:
            ValueError: If no dataset has been loaded.
        """
        if self.ds is None:
            raise ValueError("No dataset loaded. Run load() first.")
        
        output_path = Path(filename)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        self.ds.to_netcdf(output_path)
        print(f"💾 Saved to: {output_path}")
    
    def save_csv(self, filename: str):
        """Write the dataset to CSV, creating parent directories.

        The frame is written in xarray's wide form - one column per variable,
        indexed by the dataset's dimensions - not the long form produced by
        :meth:`~climdata.utils.wrapper_workflow.ClimateExtractor.to_dataframe`.

        Args:
            filename (str): Destination path.

        Returns:
            None

        Raises:
            ValueError: If no dataset has been loaded.
        """
        if self.ds is None:
            raise ValueError("No dataset loaded. Run load() first.")
        
        output_path = Path(filename)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        df = self.ds.to_dataframe()
        df.to_csv(output_path)
        print(f"💾 Saved to: {output_path}")
    
    def _map_variable_name(self, var: str) -> str:
        """Translate a CF variable name to the name ISIMIP3b files use.

        ISIMIP follows CMIP naming, so nearly all names pass through unchanged;
        ``sfcWind`` is the exception, stored lowercase as ``sfcwind``.

        Args:
            var (str): CF variable name, e.g. ``"tasmax"``.

        Returns:
            str: The ISIMIP name, or ``var`` unchanged if it is not in the table.
        """
        # CMIP6 uses standard names
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
            'psl': 'psl',        # Sea level pressure
        }
        
        return variable_map.get(var, var)
    
    def _is_file_in_date_range(self, filename: str, start_date: datetime, end_date: datetime) -> bool:
        """Test whether a file's year range overlaps the requested period.

        ISIMIP encodes coverage in the filename, e.g.
        ``gfdl-esm4_r1i1p1f1_ssp585_tas_global_daily_2015_2024.nc``. Comparison
        is at year granularity, so a file may be kept for a request that only
        touches part of it.

        Args:
            filename (str): Basename of the candidate file.
            start_date (datetime): Start of the requested period.
            end_date (datetime): End of the requested period.

        Returns:
            bool: ``True`` if the ranges overlap. Also ``True`` when the filename
            carries no parseable year range - an unrecognised name is downloaded
            rather than silently dropped.
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
