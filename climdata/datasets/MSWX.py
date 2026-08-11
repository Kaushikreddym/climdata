import pandas as pd
import geopandas as gpd
import os
from tqdm import tqdm
import warnings
from datetime import datetime, timedelta
import xarray as xr
from omegaconf import DictConfig, OmegaConf

try:
    from google.oauth2 import service_account
    from googleapiclient.discovery import build
    _GOOGLE_AVAILABLE = True
except ImportError:
    _GOOGLE_AVAILABLE = False

from climdata.utils.utils_download import list_drive_files, download_drive_file
from shapely.geometry import mapping
import cf_xarray

warnings.filterwarnings("ignore", category=Warning)

class MSWXmirror:
    def __init__(self, cfg: DictConfig):
        if not _GOOGLE_AVAILABLE:
            raise ImportError(
                "MSWX requires google-api-python-client and google-auth. "
                "Install with: pip install google-api-python-client google-auth"
            )
        self.cfg = cfg
        self.dataset = None
        self.variables = cfg.variables
        self.files = []
        self._extract_mode = None
        self._extract_params = None

    def _load_opts(self):
        """Read the optional ``cfg.load`` block, falling back to sensible defaults.

        Keeping this defensive means configs/overrides without a ``load:`` section
        continue to work unchanged.
        """
        defaults = {
            "engine": "h5netcdf",
            "parallel": True,
            # One chunk per global frame. Leaving lat/lon unset makes xarray
            # inherit the file's native ~(8, 32) on-disk chunking, i.e. ~25k
            # tiny chunks per frame (~46M tasks/global-year) that stall dask's
            # graph optimisation. Pin the spatial dims to the whole extent.
            "chunks": {"time": 1, "lat": -1, "lon": -1},
            "fix_coords": True,
        }
        opts = dict(defaults)

        chunks = OmegaConf.select(self.cfg, "load.chunks", default="__missing__")
        if chunks != "__missing__":
            # Allow either a mapping of dim -> size, or the string "auto".
            if isinstance(chunks, str):
                opts["chunks"] = chunks
            else:
                chunks = OmegaConf.to_container(chunks, resolve=True)
                # Any spatial dim the caller didn't pin defaults to the whole
                # extent (-1) rather than the tiny native on-disk chunking.
                for dim in ("lat", "lon"):
                    chunks.setdefault(dim, -1)
                opts["chunks"] = chunks

        for key in ("engine", "parallel", "fix_coords"):
            val = OmegaConf.select(self.cfg, f"load.{key}", default="__missing__")
            if val != "__missing__":
                opts[key] = val

        # When the Dask cluster is disabled, treat extraction as eager: the
        # dataset is materialised into memory (NumPy) at the end of load() so
        # nothing downstream stays lazily dask-backed. `load.dask.enabled`
        # defaults to False (opt-in cluster), matching the wrapper's default.
        opts["eager"] = not bool(
            OmegaConf.select(self.cfg, "load.dask.enabled", default=False)
        )

        return opts

    def _fix_coords(self, ds: xr.Dataset | xr.DataArray):
        """Ensure latitude is ascending and longitude is in the range [-180, 180).

        Latitude is flipped with a plain slice reversal instead of ``sortby``
        when it is monotonically decreasing: ``sortby`` builds an argsort +
        fancy-index reindex layer *per file*, which — over ~1800 daily global
        files — explodes the dask graph and stalls the client in graph
        optimisation. A ``[::-1]`` reversal is exact for a monotone axis and
        adds essentially no graph.
        """
        lat_name = ds.cf["latitude"].name
        lat = ds[lat_name]
        # A single .values read of the (tiny) 1-D coord; cheaper than argsort.
        lat_vals = lat.values
        if lat_vals[0] > lat_vals[-1]:               # descending -> flip
            ds = ds.isel({lat_name: slice(None, None, -1)})
        elif not (lat_vals[1:] >= lat_vals[:-1]).all():
            ds = ds.sortby(lat_name)                 # non-monotone: fall back

        lon_name = ds.cf["longitude"].name
        lon = ds[lon_name]

        # Only convert if longitudes are in the [0, 360] convention
        if lon.max().item() > 180:
            ds = ds.assign_coords(
                {lon_name: ((lon + 180) % 360) - 180}
            )
            ds = ds.sortby(lon_name)

        return ds
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
            local_path = os.path.join(output_dir, self.cfg.dataset.upper(), variable.upper(), basename)
            if os.path.exists(local_path):
                local_files.append(basename)
            else:
                missing_files.append(basename)
        
        if not missing_files:
            print(f"✅ All {len(expected_files)} {variable} files already exist locally.")
            return local_files

        print(f"📂 {len(local_files)} exist, {len(missing_files)} missing — fetching {variable} from Drive...")

        SCOPES = ['https://www.googleapis.com/auth/drive.readonly']
        key_file = self.cfg.dsinfo.MSWX.params.google_service_account
        # Hydra/YAML can hand back the literal string "None" as well as a real null.
        if key_file in (None, "None", ""):
            raise ValueError(
                "MSWX downloads need a Google service-account key. Add the override "
                "'dsinfo.MSWX.params.google_service_account=/path/to/service.json' "
                "(note the key is case-sensitive and must not be prefixed with '+')."
            )
        key_file = os.path.expanduser(str(key_file))
        if not os.path.exists(key_file):
            raise FileNotFoundError(f"MSWX service-account key not found: {key_file}")

        creds = service_account.Credentials.from_service_account_file(
            key_file, scopes=SCOPES
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
            local_path = os.path.join(output_dir, self.cfg.dataset.upper(), variable.upper(), filename)
            os.makedirs(os.path.dirname(local_path), exist_ok=True)
            print(f"⬇️ Downloading {filename} ...")
            download_drive_file(file['id'], local_path, service)
            local_files.append(filename)

        return local_files
    def _extract_preprocess(self, ds):
        """Apply extraction to a single-daily dataset during preprocessing."""

        # Fix coords first (can be disabled via cfg.load.fix_coords for speed
        # when the source grid is already ascending-lat / [0,360]-lon).
        if getattr(self, "_fix_coords_enabled", True):
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
        folder_id = self.cfg.dsinfo[self.cfg.dataset.upper()]["variables"][variable]["folder_id"]
        files = self.fetch(folder_id, variable)
        if not files:
            raise RuntimeError(f"No files found for variable '{variable}' in Drive or local directory.")

        # Full paths
        file_paths = [
            os.path.join(self.cfg.data_dir, self.cfg.dataset.upper(), variable.upper(), f)
            for f in files
        ]

        # MSWX internal variable name
        varname = self.cfg.dsinfo[self.cfg.dataset.upper()].variables[variable].name

        # Resolve open/preprocess options from cfg.load (with defaults).
        opts = self._load_opts()
        self._fix_coords_enabled = opts["fix_coords"]

        # Optional: preprocess each file (e.g., rename variable)
        def preprocess(ds):
            ds = self._extract_preprocess(ds)
            return ds[[varname]].rename({varname: variable})

        # Open all files as a single, lazily-loaded (Dask-backed) dataset.
        # `chunks` triggers xarray's Dask backend so data stays on disk until
        # explicitly computed; `parallel=True` opens files concurrently.
        try:
            dset = xr.open_mfdataset(
                file_paths,
                combine="nested",
                concat_dim="time",
                parallel=opts["parallel"],   # uses Dask for parallel file opening
                engine=opts["engine"],       # e.g. h5netcdf, faster than netcdf4
                chunks=opts["chunks"],       # Dask chunks -> lazy access
                preprocess=preprocess,
            )
        except Exception as e:
            raise RuntimeError(f"Failed to load dataset for variable '{variable}': {e}")

        # Ensure consistent dimension order
        if self._extract_mode != "point":
            dset = dset.transpose("time", "lat", "lon")

        # Eager mode (Dask cluster disabled): pull the data into memory now so
        # the returned dataset is NumPy-backed, not lazily Dask-backed. The
        # bounded chunking above keeps this a straightforward read rather than a
        # graph-optimisation stall.
        if opts["eager"]:
            dset = dset.load()

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
