"""Agricultural impact-model output from ISIMIP.

Unlike the other providers, which serve climate variables, this serves the
*output* of crop models — yields, biomass, planting dates — produced by models
such as DSSAT, EPIC, GEPIC, LPJmL and PROMET driven by CMIP6 projections.

That difference shapes the API. A file is identified by a five-part combination
(impact model, climate forcing, period, crop, irrigation regime), and which
combinations actually exist cannot be guessed. So this class is the only one in
climdata with a discovery step: :meth:`AGRI_ISIMIP.query_available_metadata`
asks the API what exists for a configuration before anything is requested, and
the ``AVAILABLE_*`` properties raise until it has run.

More info: https://www.isimip.org/

Example:
    >>> agri = AGRI_ISIMIP(cfg)                                  # doctest: +SKIP
    >>> agri.explore()             # what exists for this config  # doctest: +SKIP
    >>> agri.fetch()                                             # doctest: +SKIP
    >>> ds = agri.load("yield-mai-noirr")                        # doctest: +SKIP
"""

import pandas as pd
import xarray as xr
from pathlib import Path
from omegaconf import DictConfig
from typing import Optional, Tuple, Dict, List, Set
import warnings
import re

warnings.filterwarnings("ignore", category=Warning)


class AGRI_ISIMIP:
    """Discover, download and subset ISIMIP agricultural model output.

    The lifecycle has an extra first step compared to the climate providers:
    :meth:`query_available_metadata` (or :meth:`explore`, which loops it over
    every configured combination) must run before the ``AVAILABLE_*`` properties
    return anything, because what a given impact model published for a given
    forcing and period can only be learned from the API.

    After that the usual ``fetch`` → ``load`` → ``extract`` applies. Requests are
    described by ``cfg.query_configs``, a list of ``impact_model`` /
    ``climate_forcing`` / ``period`` combinations.

    Attributes:
        cfg (DictConfig): Configuration with ``query_configs`` and ``data_dir``.
        client (ISIMIPClient): ISIMIP API client, created in :meth:`__init__`.
        ds (xr.Dataset | None): The loaded dataset, set by :meth:`load`.
        downloaded_files (list[str]): Local paths accumulated by :meth:`fetch`.

    Example:
        >>> agri = AGRI_ISIMIP(cfg)                                    # doctest: +SKIP
        >>> agri.query_available_metadata("LPJmL", "gfdl-esm4", "future")
        >>> sorted(agri.AVAILABLE_CROPS)[:3]                           # doctest: +SKIP
        ['mai', 'ri1', 'soy']
    """

    def __init__(self, cfg: DictConfig):
        """Bind a configuration and open an ISIMIP client.

        No metadata is queried here; the ``AVAILABLE_*`` caches start empty.

        Args:
            cfg (DictConfig): Configuration with ``query_configs``, ``data_dir``
                and ``variables``.

        Raises:
            ImportError: If ``isimip-client`` is not installed.
        """
        self.cfg = cfg
        self.ds = None
        self.downloaded_files = []
        self._extract_mode = None
        self._extract_params = None

        # Initialize ISIMIP client
        try:
            from isimip_client.client import ISIMIPClient
            self.client = ISIMIPClient()
        except ImportError:
            raise ImportError(
                "isimip-client is required for AGRI-ISIMIP data access. "
                "Install it with: pip install isimip-client"
            )

        # Cached metadata from API queries
        self._available_crops: Optional[Set[str]] = None
        self._available_variables: Optional[Set[str]] = None
        self._available_irrigations: Optional[Set[str]] = None
        self._available_scenarios: Optional[Set[str]] = None

    @property
    def AVAILABLE_CROPS(self) -> Set[str]:
        """Crops the API reported for the queried configuration.

        Returns:
            set[str]: Crop codes, e.g. ``{"mai", "ri1", "soy", "swh", "wwh"}``.

        Raises:
            ValueError: If :meth:`query_available_metadata` has not run.
        """
        if self._available_crops is None:
            raise ValueError("No crops loaded. Call query_available_metadata() first.")
        return self._available_crops

    @property
    def AVAILABLE_VARIABLES(self) -> Set[str]:
        """Variables the API reported for the queried configuration.

        Returns:
            set[str]: Variable names, e.g. ``{"yield", "biom", "plantyear"}``.

        Raises:
            ValueError: If :meth:`query_available_metadata` has not run.
        """
        if self._available_variables is None:
            raise ValueError("No variables loaded. Call query_available_metadata() first.")
        return self._available_variables

    @property
    def AVAILABLE_IRRIGATIONS(self) -> Set[str]:
        """Irrigation regimes the API reported for the queried configuration.

        Returns:
            set[str]: Regime codes — typically ``{"noirr", "firr"}`` for rainfed
            and fully irrigated.

        Raises:
            ValueError: If :meth:`query_available_metadata` has not run.
        """
        if self._available_irrigations is None:
            raise ValueError("No irrigations loaded. Call query_available_metadata() first.")
        return self._available_irrigations

    @property
    def AVAILABLE_SCENARIOS(self) -> Set[str]:
        """Climate scenarios the API reported for the queried configuration.

        Returns:
            set[str]: Scenario names, e.g. ``{"historical", "ssp126", "ssp585"}``.

        Raises:
            ValueError: If :meth:`query_available_metadata` has not run.
        """
        if self._available_scenarios is None:
            raise ValueError("No scenarios loaded. Call query_available_metadata() first.")
        return self._available_scenarios

    def query_available_metadata(self, impact_model: str, climate_forcing: str, period: str) -> None:
        """Query ISIMIP API to get all available metadata for a specific configuration.
        Extracts crops, variables, irrigation types, and scenarios.

        Args:
            impact_model (str): Impact model name (e.g., 'SIMPLACE-LINTUL5')
            climate_forcing (str): Climate forcing/model (e.g., 'mpi-esm1-2-hr')
            period (str): Time period (e.g., 'historical')

        Raises:
            Exception: If API query fails or no datasets found
        """
        path = f'ISIMIP3b/OutputData/agriculture/{impact_model}/{climate_forcing}/{period}'
        
        crops: Set[str] = set()
        variables: Set[str] = set()
        irrigation: Set[str] = set()
        scenarios: Set[str] = set()
        
        page = 1
        total_fetched = 0
        
        while True:
            response = self.client.datasets(
                path=path,
                page=page,
                page_size=100
            )
            
            results = response.get('results', [])
            
            if not results:
                break
            
            total_fetched += len(results)
            
            for ds in results:
                spec = ds.get('specifiers', {})
                
                crop = spec.get('crop')
                var = spec.get('variable')
                irr = spec.get('irrigation')
                scen = spec.get('climate_scenario')
                
                if crop:
                    crops.add(crop)
                if var:
                    variables.add(var)
                if irr:
                    irrigation.add(irr)
                if scen:
                    scenarios.add(scen)
            
            # Check pagination
            if not response.get('next'):
                break
            
            page += 1
        
        self._available_crops = crops
        self._available_variables = variables
        self._available_irrigations = irrigation
        self._available_scenarios = scenarios
        
        print(f"✅ Fetched {total_fetched} datasets from ISIMIP API")
        print(f"   Crops: {sorted(crops)}")
        print(f"   Variables: {sorted(variables)}")
        print(f"   Irrigation: {sorted(irrigation)}")
        print(f"   Scenarios: {sorted(scenarios)}")

    def explore(self):
        """Print what ISIMIP holds for every configured query combination.

        Runs :meth:`query_available_metadata` over each entry of
        ``cfg.query_configs`` and prints the crops, variables, irrigation
        regimes and scenarios found. The intended first call on a new
        configuration, before committing to a download.

        Note that the ``AVAILABLE_*`` caches are overwritten by each combination
        in turn, so after this they reflect the last one only.

        Returns:
            None: The report is written to stdout.
        """
        print("\n" + "=" * 70)
        print("🌾 AGRI-ISIMIP Dataset Explorer")
        print("=" * 70)
        
        query_configs = self.cfg.get('query_configs', [])
        
        if not query_configs:
            print("\n⚠️  No query_configs found in configuration")
            return
        
        for i, config in enumerate(query_configs, 1):
            print(f"\n[{i}/{len(query_configs)}] Configuration:")
            print(f"  Impact Model: {config['impact_model']}")
            print(f"  Climate Forcing: {config['climate_forcing']}")
            print(f"  Period: {config['period']}")
            print("-" * 70)
            
            try:
                self.query_available_metadata(
                    impact_model=config['impact_model'],
                    climate_forcing=config['climate_forcing'],
                    period=config['period']
                )
            except Exception as e:
                print(f"❌ Error: {str(e)}")

    def list_files(self, impact_model: str, climate_forcing: str, period: str,
                   crop: Optional[str] = None, irrigation: Optional[str] = None,
                   variable: Optional[str] = None) -> List[Dict]:
        """List available files matching the given criteria.

        Args:
            impact_model (str): Impact model name (e.g., 'SIMPLACE-LINTUL5')
            climate_forcing (str): Climate forcing (e.g., 'mpi-esm1-2-hr')
            period (str): Time period (e.g., 'historical')
            crop (str, optional): Filter by crop name
            irrigation (str, optional): Filter by irrigation type ('firr' or 'noirr')
            variable (str, optional): Filter by variable name

        Returns:
            List[Dict]: List of file metadata dictionaries
        """
        path = f'ISIMIP3b/OutputData/agriculture/{impact_model}/{climate_forcing}/{period}'
        
        print(f"\n🔍 Listing AGRI-ISIMIP files:")
        print(f"   Path: {path}")
        if crop:
            print(f"   Crop: {crop}")
        if irrigation:
            print(f"   Irrigation: {irrigation}")
        if variable:
            print(f"   Variable: {variable}")

        files_list = []
        page = 1
        total_fetched = 0

        try:
            while True:
                response = self.client.datasets(
                    path=path,
                    page=page,
                    page_size=100
                )
                
                results = response.get('results', [])
                
                if not results:
                    break
                
                total_fetched += len(results)
                
                for dataset in results:
                    spec = dataset.get('specifiers', {})
                    file_crop = spec.get('crop')
                    file_irr = spec.get('irrigation')
                    file_var = spec.get('variable')
                    
                    # Apply filters
                    if crop and file_crop != crop:
                        continue
                    if irrigation and file_irr != irrigation:
                        continue
                    if variable and file_var != variable:
                        continue
                    
                    files_list.append({
                        'name': dataset.get('name', ''),
                        'path': dataset.get('path', ''),
                        'size': dataset.get('size', 'N/A'),
                        'crop': file_crop,
                        'irrigation': file_irr,
                        'variable': file_var,
                        'dataset': dataset,
                    })
                
                # Check pagination
                if not response.get('next'):
                    break
                
                page += 1
            
            print(f"✅ Found {len(files_list)} matching files (from {total_fetched} total)\n")
            
            # Display first 20 results
            if files_list:
                print("Files found:")
                print("-" * 100)
                for i, f in enumerate(files_list[:20], 1):
                    print(f"{i:3d}. {f['name']}")
                    print(f"     Path: {f['path']}")
                    print(f"     Crop: {f['crop']}, Irrigation: {f['irrigation']}, Variable: {f['variable']}")
                    if f['size'] != 'N/A':
                        size_gb = f['size'] / (1024**3) if isinstance(f['size'], (int, float)) else 'N/A'
                        if size_gb != 'N/A':
                            print(f"     Size: {size_gb:.2f} GB")
                if len(files_list) > 20:
                    print(f"\n     ... and {len(files_list) - 20} more files")
            
            return files_list

        except Exception as e:
            print(f"❌ Error listing files: {str(e)}")
            return []

    def fetch(self, impact_model: str, climate_forcing: str, period: str,
              crop: Optional[str] = None, irrigation: Optional[str] = None,
              variable: Optional[str] = None):
        """Download AGRI-ISIMIP files matching the criteria using isimip_client.download().

        Files are stored as::
        <data_dir>/AGRI_ISIMIP/<impact_model>/<climate_forcing>/<period>/<crop>/<variable>/<irrigation>/

        Args:
            impact_model (str): Impact model name
            climate_forcing (str): Climate forcing
            period (str): Time period
            crop (str, optional): Specific crop to download
            irrigation (str, optional): Specific irrigation type to download
            variable (str, optional): Specific variable to download
        """
        files_list = self.list_files(
            impact_model=impact_model,
            climate_forcing=climate_forcing,
            period=period,
            crop=crop,
            irrigation=irrigation,
            variable=variable
        )
        
        if not files_list:
            print("⚠️  No files found to download")
            return
        
        # Create directory structure
        base_dir = (
            Path(self.cfg.data_dir) /
            "AGRI_ISIMIP" /
            impact_model /
            climate_forcing /
            period
        )
        base_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"\n📥 Downloading {len(files_list)} files to {base_dir}...\n")
        
        for i, file_info in enumerate(files_list, 1):
            try:
                dataset_info = file_info['dataset']
                
                # Create subdirectory based on crop/variable/irrigation
                crop_dir = base_dir / (file_info['crop'] or 'unknown')
                var_dir = crop_dir / (file_info['variable'] or 'unknown')
                irr_dir = var_dir / (file_info['irrigation'] or 'unknown')
                irr_dir.mkdir(parents=True, exist_ok=True)
                
                # Get files from dataset
                files = dataset_info.get('files', [])
                
                if not files:
                    print(f"[{i}/{len(files_list)}] ⚠️  No files in dataset")
                    continue
                
                # Download each file using client.download()
                for file_dict in files:
                    file_url = file_dict.get('file_url')
                    file_name = file_dict.get('path', '').split('/')[-1] or 'file'
                    
                    if not file_url:
                        print(f"[{i}/{len(files_list)}] ⚠️  No file_url for {file_name}")
                        continue
                    
                    out_file = irr_dir / file_name
                    
                    # NEW: Check if file already exists
                    local_file = irr_dir / file_dict['name']

                    if local_file.exists():
                        # Check file size matches
                        local_size = local_file.stat().st_size
                        remote_size = file_dict.get('size', 0)
                        
                        if local_size == remote_size:
                            print(f"   ⏭️  Skipping (exists): {file_dict['name']}")
                            self.downloaded_files.append(str(local_file))
                            continue
                        else:
                            print(f"   ⚠️  Incomplete file, re-downloading: {file_dict['name']}")

                    # NEW: Only download if file doesn't exist or is incomplete
                    print(f"   ⬇️  Downloading: {file_dict['name']}")
                    self.client.download(
                        file_url,
                        path=str(irr_dir),
                        validate=True
                    )
                    
                    print(f"     ✅ Downloaded {file_name}")
                    self.downloaded_files.append(str(out_file))
                    
            except Exception as e:
                print(f"     ❌ Error: {str(e)}")
        
        print(f"\n✅ Downloaded {len(self.downloaded_files)} files")

    def load(self, file_path: str, use_dask: bool = True, decode_times: bool = False) -> xr.Dataset:
        """Load a NetCDF file into an xarray Dataset with custom time handling for ISIMIP data.

        Args:
            file_path (str): Path to the NetCDF file
            use_dask (bool, optional): Whether to use dask for lazy loading (default:
                True)
            decode_times (bool, optional): Whether to decode times (default: False for
                ISIMIP data with custom calendars) Set to True if your data uses
                standard time encoding

        Returns:
            xr.Dataset: Loaded xarray dataset with properly decoded time coordinates

        Note:
            ISIMIP data often uses custom time encoding like "growing seasons since 1601-01-01"
            with custom calendars. This function automatically detects and converts such times
            to proper datetime coordinates by:
            1. Loading with decode_times=False
            2. Parsing the time unit string (e.g., "since 1601-01-01")
            3. Interpreting time values as years offset from the reference year
            4. Creating proper datetime coordinates
        """
        print(f"📂 Loading {file_path}...")
        
        try:
            # Load with decode_times=False to avoid calendar errors
            if use_dask:
                self.ds = xr.open_dataset(file_path, chunks='auto', decode_times=False)
            else:
                self.ds = xr.open_dataset(file_path, decode_times=False)
            
            # Check if time dimension exists and needs custom decoding
            if 'time' in self.ds.dims and 'time' in self.ds.data_vars or 'time' in self.ds.coords:
                try:
                    # Try to decode ISIMIP custom time format
                    time_attrs = self.ds.time.attrs.get('units', '')
                    
                    # Check if it's a custom ISIMIP time format (e.g., "growing seasons since ...")
                    if 'since' in time_attrs:
                        print(f"   📅 Detected custom time format: {time_attrs}")
                        
                        # Extract start year from units string
                        match = re.search(r"since\s(.+)", time_attrs)
                        if match:
                            start_date_str = match.group(1).strip()
                            start = pd.Timestamp(start_date_str)
                            
                            # Interpret time values as years offset from start year
                            years = start.year + self.ds['time'].values.astype(int)
                            
                            # Assign proper datetime coordinates
                            self.ds = self.ds.assign_coords(
                                time=pd.to_datetime(years, format="%Y")
                            )
                            print(f"   ✅ Converted time coordinates (reference year: {start.year})")
                except Exception as time_error:
                    print(f"   ⚠️  Could not decode custom time format: {str(time_error)}")
                    print(f"   Keeping original time values")
            
            print(f"✅ Loaded dataset with dimensions: {dict(self.ds.dims)}")
            return self.ds
            
        except Exception as e:
            print(f"❌ Error loading file: {str(e)}")
            if 'decode_times' in str(e) or 'calendar' in str(e):
                print("💡 Hint: ISIMIP data uses custom calendars. Load is called with decode_times=False")
            return None

    def extract(self, *, point: Tuple[float, float] = None, box: Dict[str, float] = None,
                shapefile: str = None, buffer_km: float = 0.0):
        """Specify spatial extraction intent.

        Args:
            point (tuple, optional): (lon, lat) for point extraction
            box (dict, optional): dict(lat_min, lat_max, lon_min, lon_max) for box
                extraction
            shapefile (str, optional): Path to shapefile for shapefile extraction
            buffer_km (float, optional): Buffer in km for shapefile extraction
        """
        if point is not None:
            lon, lat = point
            self._extract_mode = "point"
            self._extract_params = (lon, lat)
            print(f"📍 Set point extraction mode: ({lon}, {lat})")

        elif box is not None:
            for k in ("lat_min", "lat_max", "lon_min", "lon_max"):
                if k not in box:
                    raise ValueError(f"box missing key {k}")
            self._extract_mode = "box"
            self._extract_params = box
            print(f"📦 Set box extraction mode: {box}")

        elif shapefile is not None:
            import geopandas as gpd
            gdf = gpd.read_file(shapefile) if isinstance(shapefile, str) else shapefile
            if buffer_km > 0:
                gdf = gdf.to_crs(epsg=3857)
                gdf["geometry"] = gdf.buffer(buffer_km * 1000)
                gdf = gdf.to_crs(epsg=4326)
            self._extract_mode = "shapefile"
            self._extract_params = gdf
            print(f"🗺️  Set shapefile extraction mode with buffer {buffer_km} km")

        else:
            raise ValueError("Must provide point, box, or shapefile.")

        return self

    def get_extracted_data(self) -> xr.Dataset:
        """Extract data based on spatial mode set via extract().

        Returns:
            xr.Dataset: Extracted dataset at specified location(s)
        """
        if self.ds is None:
            raise ValueError("No dataset loaded. Call load() first.")
        
        if self._extract_mode == "point":
            lon, lat = self._extract_params
            extracted = self.ds.sel(lon=lon, lat=lat, method='nearest')
            print(f"✅ Extracted point data at ({lon}, {lat})")
            return extracted
        
        elif self._extract_mode == "box":
            box = self._extract_params
            extracted = self.ds.sel(
                lon=slice(box['lon_min'], box['lon_max']),
                lat=slice(box['lat_min'], box['lat_max'])
            )
            print(f"✅ Extracted box data")
            return extracted
        
        elif self._extract_mode == "shapefile":
            # TODO: Implement shapefile extraction
            raise NotImplementedError("Shapefile extraction not yet implemented")
        
        else:
            raise ValueError("No extraction mode set. Call extract() first.")

    def save_csv(self, filename: str):
        """Write the loaded dataset to CSV, one row per cell per timestep.

        Args:
            filename (str): Destination path. The parent directory must exist.

        Returns:
            None

        Raises:
            ValueError: If no dataset has been loaded.
        """
        if self.ds is None:
            raise ValueError("No dataset loaded. Call load() first.")
        df = self.ds.to_dataframe().reset_index()
        df.to_csv(filename, index=False)
        print(f"✅ Saved CSV to {filename}")

    def save_netcdf(self, filename: str):
        """Write the loaded dataset to NetCDF.

        Args:
            filename (str): Destination path. The parent directory must exist.

        Returns:
            None

        Raises:
            ValueError: If no dataset has been loaded.
        """
        if self.ds is None:
            raise ValueError("No dataset loaded. Call load() first.")
        self.ds.to_netcdf(filename)
        print(f"✅ Saved NetCDF to {filename}")
