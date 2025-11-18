import intake
import xarray as xr
import pandas as pd
from omegaconf import DictConfig
import intake
import xarray as xr
import pandas as pd
from omegaconf import DictConfig


class CMIPCloud:
    def __init__(self, cfg: DictConfig):
        # Directly read from flat config
        self.experiment_id = cfg.experiment_id
        self.source_id = cfg.source_id
        self.table_id = cfg.table_id
        self.variables = cfg.variables
        self.start_date = cfg.time_range.start_date
        self.end_date = cfg.time_range.end_date

        self.col_subsets = []
        self.ds = None
        self.col = None

    def fetch(self):
        """Collect intake catalog subsets for each variable."""
        col = intake.open_esm_datastore(
            "https://storage.googleapis.com/cmip6/pangeo-cmip6.json"
        )
        self.col_subsets = []

        for var in self.variables:
            query = dict(
                experiment_id=[self.experiment_id],
                source_id=self.source_id,
                table_id=self.table_id,
                variable_id=var,
            )
            col_subset = col.search(require_all_on=["source_id"], **query)

            if len(col_subset.df) == 0:
                continue

            self.col_subsets.append(col_subset)
            self.col = col

        # ✔ raise error if nothing found
        if not self.col_subsets:
            raise ValueError(
                f"No matching CMIP6 data found for: "
                f"experiment_id={self.experiment_id}, "
                f"source_id={self.source_id}, "
                f"table_id={self.table_id}, "
                f"variables={self.variables}"
            )

        return self.col_subsets

    def load(self):
        """Load and merge datasets from collected col_subsets."""
        datasets = []
        for col_subset in self.col_subsets:
            zstore_path = col_subset.df.zstore.values[0].replace(
                "gs:/", "https://storage.googleapis.com"
            )
            ds_var = xr.open_zarr(zstore_path)
            datasets.append(ds_var)
        if datasets:
            self.ds = xr.merge(datasets)
        else:
            self.ds = None
        return self.ds

    def extract(self, *, point=None, box=None, shapefile=None, buffer_km=0.0):
        """
        Extract a subset of the dataset by point, bounding box (dict), or shapefile.
        """
        import geopandas as gpd
        from shapely.geometry import mapping

        if self.ds is None:
            raise ValueError("No dataset loaded. Call `load()` first.")
        
        self._subset_time(self.start_date, self.end_date) 
        
        ds = self.ds
        if point is not None:
            lon, lat = point
            if buffer_km > 0:
                buffer_deg = buffer_km / 111
                ds_subset = ds.sel(
                    lon=slice(lon - buffer_deg, lon + buffer_deg),
                    lat=slice(lat - buffer_deg, lat + buffer_deg),
                )
            else:
                ds_subset = ds.sel(lon=lon, lat=lat, method="nearest")

        elif box is not None:
            ds_subset = ds.sel(
                lon=slice(box["lon_min"], box["lon_max"]),
                lat=slice(box["lat_min"], box["lat_max"]),
            )

        elif shapefile is not None:
            if isinstance(shapefile, str):
                gdf = gpd.read_file(shapefile)
            else:
                gdf = shapefile
            if buffer_km > 0:
                gdf = gdf.to_crs(epsg=3857)
                gdf["geometry"] = gdf.buffer(buffer_km * 1000)
                gdf = gdf.to_crs(epsg=4326)
            geom = [mapping(g) for g in gdf.geometry]
            import rioxarray

            ds = ds.rio.write_crs("EPSG:4326", inplace=False)
            ds_subset = ds.rio.clip(geom, gdf.crs, drop=True)

        else:
            raise ValueError("Must provide either point, box, or shapefile.")
        self.ds = ds_subset
        return ds_subset

    def _subset_time(self, start_date, end_date):
        """Subset the dataset by time range."""
        if self.ds is None:
            return None
        ds_time = self.ds.sel(time=slice(start_date, end_date))
        self.ds = ds_time
        return ds_time

    def save_netcdf(self, filename):
        if self.ds is not None:
            if "time" in self.ds.variables:
                self.ds["time"].encoding.clear()
            self.ds.to_netcdf(filename)
            print(f"Saved NetCDF to {filename}")

    def save_zarr(self, store_path):
        if self.ds is not None:
            self.ds.to_zarr(store_path, mode="w")
            print(f"Saved Zarr to {store_path}")

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
            lambda v: self.ds[v].attrs.get("units", "unknown")
            if v in self.ds.data_vars
            else "unknown"
        )

        df_long["source"] = self.source_id
        df_long["experiment"] = self.experiment_id
        df_long["table"] = self.table_id

        cols = [
            "source",
            "experiment",
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
        if self.ds is not None:
            df = self.ds.to_dataframe().reset_index()
            df = self._format(df)
            df.to_csv(filename, index=False)
            print(f"Saved CSV to {filename}")