import os
import pandas as pd
import hydra
from wetterdienst import Settings
from wetterdienst.provider.dwd.observation import DwdObservationRequest

class DWDmirror:
    def __init__(self, cfg):
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
        param_info = self.param_mapping.dwd.variables[variable]
        resolution = param_info["resolution"]  # e.g., "daily"
        dataset = param_info["dataset"]        # e.g., "temperature_air_200"
        
        # Create request with updated API
        request = DwdObservationRequest(
            parameters=(resolution, dataset)   # <-- new API
        )
        
        # Get station metadata as pandas DataFrame
        stations_df = request.all().df.to_pandas()
        
        # Store stations in the object
        self.stations = stations_df
        return stations_df
    def load(self, variable, lat_loc, lon_loc, buffer_km = 50):
        param_info = self.param_mapping.dwd.variables[variable]
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
        param_info = self.param_mapping.dwd.variables[variable]
        resolution = param_info["resolution"]
        dataset = param_info["dataset"]
        variable_name = param_info["name"]
        if not hasattr(self, "stations"):
            self.get_stations()

        stations_df = self.stations

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

        data = request.values.all().df.to_pandas()  # pandas DataFrame

        # Convert to xarray
        ds = data.reset_index().set_index(["station_id", "date"]).to_xarray()
        ds = ds.assign_coords(date=ds["date"].to_pandas().tz_localize(None))
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
        ds = ds.drop_vars(["index", "resolution", "dataset", "parameter"])
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
        self.df['units'] = self.param_mapping.dwd.variables[variable].unit
        self.df = self.df[["lat", "lon", "time", "source", "variable", "value", "units"]]
        # self.df = df
        return self.df

    def save_csv(self,filename):
        self.df.to_csv(filename, index=False)
        print(f"✅ Saved time series to: {filename}")
        return filename
    