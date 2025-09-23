import intake
import xarray as xr
import pandas as pd

class CMIPCloud:
    def __init__(self, experiment_id, source_id, table_id, variables, region_bounds=None):
        self.experiment_id = experiment_id
        self.source_id = source_id
        self.table_id = table_id
        self.variables = variables
        self.region_bounds = region_bounds
        self.col_subsets = []
        self.ds = None

    def fetch(self):
        """Collect intake catalog subsets for each variable."""
        col = intake.open_esm_datastore("https://storage.googleapis.com/cmip6/pangeo-cmip6.json")
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
        return self.col_subsets

    def load(self):
        """Load and merge datasets from collected col_subsets."""
        datasets = []
        for col_subset in self.col_subsets:
            zstore_path = col_subset.df.zstore.values[0].replace('gs:/', "https://storage.googleapis.com")
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

        ds = self.ds

        if point is not None:
            lon, lat = point
            if buffer_km > 0:
                buffer_deg = buffer_km / 111
                ds_subset = ds.sel(
                    lon=slice(lon-buffer_deg, lon+buffer_deg),
                    lat=slice(lat-buffer_deg, lat+buffer_deg)
                )
            else:
                ds_subset = ds.sel(lon=lon, lat=lat, method="nearest")

        elif box is not None:
            # Accept dict: {'lat_min': ..., 'lat_max': ..., 'lon_min': ..., 'lon_max': ...}
            ds_subset = ds.sel(
                lon=slice(box['lon_min'], box['lon_max']),
                lat=slice(box['lat_min'], box['lat_max'])
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
        """
        Subset the dataset by time range.
        Dates should be strings in 'YYYY-MM-DD' format.
        """
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

    def save_csv(self, filename):
        if self.ds is not None:
            df = self.ds.to_dataframe().reset_index()
            df.to_csv(filename, index=False)
            print(f"Saved CSV to {filename}")