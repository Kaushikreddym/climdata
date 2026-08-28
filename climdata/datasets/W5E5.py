"""
W5E5 dataset access via ISIMIP client

W5E5 (WFDE5 over land merged with ERA5 over the ocean) is a global meteorological 
forcing dataset available through ISIMIP (Inter-Sectoral Impact Model Intercomparison Project).
It provides daily climate data at 0.5° resolution from 1979 onwards.

This module uses the isimip-client library to search and download W5E5 data from the 
ISIMIP data repository.
"""

import xarray as xr
from pathlib import Path
from datetime import datetime
from omegaconf import DictConfig
from typing import Optional, Tuple, Dict, List
import warnings

warnings.filterwarnings("ignore", category=Warning)

class W5E5:
    """Download and assemble W5E5 climate data from the ISIMIP repository.

    W5E5 is distributed through ISIMIP3a as observational input data and serves
    as the bias-adjustment reference for ISIMIP3b projections — pair it with
    :class:`~climdata.datasets.CMIP_W5E5.CMIPW5E5` to compare a projection
    against the reference it was adjusted to.

    The lifecycle is ``fetch`` → ``load`` → ``extract``. :meth:`extract` may be
    called before :meth:`load`, in which case the subsetting instructions are
    stored and applied as soon as data arrives.

    Attributes:
        cfg (DictConfig): Configuration supplying ``variables``, ``time_range``,
            ``data_dir`` and ``dataset``.
        ds (xr.Dataset | None): The loaded dataset, set by :meth:`load`.
        client (ISIMIPClient): ISIMIP API client, created in :meth:`__init__`.
        downloaded_files (list[str]): Local paths accumulated by :meth:`fetch`.

    Example:
        >>> w5e5 = W5E5(cfg)                                  # doctest: +SKIP
        >>> w5e5.fetch()                                      # doctest: +SKIP
        >>> w5e5.load()                                       # doctest: +SKIP
        >>> w5e5.extract(point=(13.4, 52.5))                  # doctest: +SKIP
    """
    
    def __init__(self, cfg: DictConfig):
        """Bind a configuration and open an ISIMIP client.

        Args:
            cfg (DictConfig): Configuration with ``variables``,
                ``time_range.start_date`` / ``.end_date``, ``data_dir`` and
                ``dataset``.

        Raises:
            ImportError: If ``isimip-client`` is not installed.
        """
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
        """Search ISIMIP and download the W5E5 files covering the request.

        Queries the ISIMIP3a ``obsclim`` catalogue once per configured variable,
        keeps the files whose filename year-range overlaps the requested period
        (see :meth:`_is_file_in_date_range`), and downloads any that are not
        already on disk under ``<data_dir>/<DATASET>/<variable>/``.

        Per-variable failures are reported and skipped rather than raised, so one
        unavailable variable does not abort a multi-variable request — check
        :attr:`downloaded_files` to see what actually arrived.

        Returns:
            None: Paths are appended to :attr:`downloaded_files`.
        """
        print("🔍 Searching for W5E5 datasets in ISIMIP repository...")
        
        start_date = datetime.fromisoformat(self.cfg.time_range.start_date)
        end_date = datetime.fromisoformat(self.cfg.time_range.end_date)
        
        output_dir = Path(self.cfg.data_dir) / self.cfg.dataset.upper()
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

                results = self._normalize_dataset_results(response)
                if not results:
                    print(f"⚠️ No W5E5 datasets found for {var}")
                    continue

                # Get the first matching dataset
                dataset = results[0]
                print(f"✅ Found dataset: {dataset.get('name', 'unnamed')}")
                
                # Filter files by date range
                for file_info in self._extract_files_list(dataset):
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
        spans several decades, merged across variables, and finally clipped to
        the configured time range.

        Returns:
            None: The dataset is stored on :attr:`ds`.

        Raises:
            ValueError: If :meth:`fetch` has not run, or produced no files.
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
        """Record a spatial subset, and apply it if data is already loaded.

        The instruction is stored rather than executed, so this may be called
        before :meth:`load`; :meth:`_apply_extraction` runs it once :attr:`ds`
        exists. The three modes are mutually exclusive — the first one supplied
        wins, in the order ``point``, ``box``, ``shapefile``. Calling with no
        argument clears nothing and does nothing.

        Args:
            point (tuple[float, float], optional): ``(lon, lat)`` in degrees.
                Note the order: longitude first.
            box (dict, optional): Bounding box with keys ``lon_min``, ``lon_max``,
                ``lat_min``, ``lat_max``.
            shapefile (str | geopandas.GeoDataFrame, optional): Polygon(s) to clip
                to, as a path or an in-memory frame.
            buffer_km (float): Half-width of a box averaged around ``point``,
                converted at a flat 111 km per degree. ``0.0``, the default,
                selects the single nearest grid cell instead.

        Returns:
            None
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

        The frame is written in xarray's wide form — one column per variable,
        indexed by the dataset's dimensions — not the long form produced by
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
        """Translate a CF variable name to the name W5E5 files use.

        W5E5 follows CMIP naming, so nearly all names pass through unchanged;
        ``sfcWind`` is the exception, stored lowercase as ``sfcwind``.

        Args:
            var (str): CF variable name, e.g. ``"tasmax"``.

        Returns:
            str: The W5E5 name, or ``var`` unchanged if it is not in the table.
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
        """Test whether a file's year range overlaps the requested period.

        W5E5 encodes coverage in the filename, e.g.
        ``w5e5v2.0_obsclim_tas_global_daily_1979_1989.nc``. Comparison is at
        year granularity, so a file may be kept for a request that only touches
        part of it.

        Args:
            filename (str): Basename of the candidate file.
            start_date (datetime): Start of the requested period.
            end_date (datetime): End of the requested period.

        Returns:
            bool: ``True`` if the ranges overlap. Also ``True`` when the filename
            carries no parseable year range — an unrecognised name is downloaded
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
