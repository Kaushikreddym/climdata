import os
import pandas as pd
import geopandas as gpd
import hydra
try:
    from wetterdienst import Settings
    from wetterdienst.provider.dwd.observation import DwdObservationRequest
    _WETTERDIENST_AVAILABLE = True
except ImportError:
    _WETTERDIENST_AVAILABLE = False

class DWDmirror:
    def __init__(self, cfg):
        if not _WETTERDIENST_AVAILABLE:
            raise ImportError(
                "DWD requires wetterdienst. "
                "Install with: pip install wetterdienst"
            )
        self.cfg = cfg
        self.param_mapping = cfg.dsinfo
        self.start_date = cfg.time_range.start_date
        self.end_date = cfg.time_range.end_date
        self.df = None
    def get_stations(self, variable='pr'):
        """
        Load DWD station metadata for the chosen parameter using the updated wetterdienst API.
        """
        # Lookup info from your mapping
        param_info = self.param_mapping.DWD.variables[variable]
        resolution = param_info["resolution"]  # e.g., "daily"
        dataset = param_info["dataset"]        # e.g., "climate_summary"

        # Create request with updated API
        request = DwdObservationRequest(
            parameters=[(resolution, dataset)]
        )

        # Get station metadata as pandas DataFrame
        stations_df = request.all().df.to_pandas()

        # Cache per dataset so different variables with different datasets
        # (e.g. climate_summary vs solar) don't share the same station list
        if not hasattr(self, "_stations_cache"):
            self._stations_cache = {}
        self._stations_cache[dataset] = stations_df
        self.stations = stations_df
        return stations_df
    def load(self, variable, lat_loc, lon_loc, buffer_km = 50):
        param_info = self.param_mapping.DWD.variables[variable]
        resolution = param_info["resolution"]
        dataset = param_info["dataset"]
        variable_name = param_info["name"]

        settings = Settings(ts_shape="long", ts_humanize=True)
        request = DwdObservationRequest(
            parameters=(resolution, dataset, variable_name),
            start_date=self.start_date,
            end_date=self.end_date,
            settings=settings
        ).filter_by_distance(
            latlon=(lat_loc, lon_loc),
            distance=buffer_km,
            unit="km"
        )

        df = request.values.all().df.to_pandas()
        self.df = df
        return self.df
    def extract(self, *,variable, point=None, box=None, shapefile=None, buffer_km=25.0):
        param_info = self.param_mapping.DWD.variables[variable]
        resolution = param_info["resolution"]
        dataset = param_info["dataset"]
        variable_name = param_info["name"]

        # Fetch station list for this variable's dataset (cached per dataset)
        cache = getattr(self, "_stations_cache", {})
        if dataset not in cache:
            self.get_stations(variable=variable)
        stations_df = self._stations_cache[dataset]

        # ---- Point extraction ----
        if point is not None:
            lon, lat = point
            if buffer_km > 0:
                buffer_deg = buffer_km / 111
                subset = stations_df[
                    (stations_df.longitude.between(lon - buffer_deg, lon + buffer_deg)) &
                    (stations_df.latitude.between(lat - buffer_deg, lat + buffer_deg))
                ]
            else:
                # Find nearest station
                subset = stations_df.copy()
                subset["distance"] = ((subset.longitude - lon)**2 + (subset.latitude - lat)**2)**0.5
                subset = subset.nsmallest(1, "distance")
            
        # ---- Box extraction ----
        elif box is not None:
            subset = stations_df[
                (stations_df.longitude.between(box["lon_min"], box["lon_max"])) &
                (stations_df.latitude.between(box["lat_min"], box["lat_max"]))
            ]

        # ---- Shapefile extraction ----
        elif shapefile is not None:
            if isinstance(shapefile, str):
                gdf = gpd.read_file(shapefile)
            else:
                gdf = shapefile

            # Optional buffer
            if buffer_km > 0:
                gdf = gdf.to_crs(epsg=3857)
                gdf["geometry"] = gdf.buffer(buffer_km * 1000)
                gdf = gdf.to_crs(epsg=4326)

            points = gpd.GeoDataFrame(
                stations_df, geometry=gpd.points_from_xy(stations_df.longitude, stations_df.latitude), crs="EPSG:4326"
            )
            # Keep stations inside any of the geometries
            mask = points.geometry.apply(lambda p: any(g.contains(p) for g in gdf.geometry))
            subset = stations_df[mask]

        else:
            raise ValueError("Must provide either point, box, or shapefile.")

        # ---- Download data from selected stations ----
        station_ids = subset.index.tolist()
        if not station_ids:
            raise ValueError("No stations found in selection.")

        request = DwdObservationRequest(
            parameters=(resolution, dataset, variable_name),
            start_date=self.start_date,
            end_date=self.end_date,
        ).filter_by_station_id(station_id=station_ids)
        print(point)
        data = request.values.all().df.to_pandas()  # pandas DataFrame

        # Normalise column names across wetterdienst versions
        col_map = {}
        for old, new in [("datetime", "date"), ("stationid", "station_id"), ("station", "station_id")]:
            if old in data.columns and new not in data.columns:
                col_map[old] = new
        if col_map:
            data = data.rename(columns=col_map)

        # If station_id is the index (some versions), promote it to a column
        if data.index.name == "station_id":
            data = data.reset_index()

        # Strip timezone and cast to nanosecond precision so xarray won't warn
        import pandas as pd
        if isinstance(data["date"].dtype, pd.DatetimeTZDtype):
            data["date"] = data["date"].dt.tz_localize(None)
        data["date"] = data["date"].astype("datetime64[ns]")

        # Convert to xarray
        ds = data.set_index(["station_id", "date"]).to_xarray()
        ds = ds.rename(
            {
                "date": "time",
                "value": variable,
                "quality": f"quality_{variable}"
            }
        )
        attrs = {}
        for attr,n in zip(["resolution", "dataset", "parameter"],["resolution", "dataset", "long_name"]):
            attrs[n] = ds[attr].values[0, 0]

        ds = ds.assign_attrs(attrs)
        vars_to_drop = [v for v in ["index", "resolution", "dataset", "parameter"] if v in ds]
        ds = ds.drop_vars(vars_to_drop)
        # Assign 'units' attribute so xclim can read/convert units
        unit = self.param_mapping.DWD.variables[variable].unit
        ds[variable].attrs["units"] = unit
        self.dataset = ds
        return ds
    def format(self, variable, lat_loc, lon_loc):
        self.df['date'] = pd.to_datetime(self.df['date'])
        self.df = self.df.groupby(['date']).agg({
            'value': 'mean',
            'station_id': lambda x: x.mode().iloc[0] if not x.mode().empty else None,
            'resolution': lambda x: x.mode().iloc[0] if not x.mode().empty else None,
            'dataset': lambda x: x.mode().iloc[0] if not x.mode().empty else None,
            'parameter': lambda x: x.mode().iloc[0] if not x.mode().empty else None,
            'quality': lambda x: x.mode().iloc[0] if not x.mode().empty else None,
        }).reset_index()

        self.df = self.df.rename(columns={
            "date": "time",
            "value": "value",
            "station_id": "frequent_station",
        })
        self.df["variable"] = variable
        self.df["lat"] = lat_loc
        self.df["lon"] = lon_loc
        self.df['source'] = 'DWD'
        self.df['units'] = self.param_mapping.DWD.variables[variable].unit
        self.df = self.df[["lat", "lon", "time", "source", "variable", "value", "units"]]
        # self.df = df
        return self.df

    def save_csv(self,filename):
        self.df.to_csv(filename, index=False)
        print(f"✅ Saved time series to: {filename}")
        return filename
    