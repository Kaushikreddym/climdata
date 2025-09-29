import os
import pandas as pd
import hydra
from wetterdienst import Settings
from wetterdienst.provider.dwd.observation import DwdObservationRequest

class DWDmirror:
    def __init__(self, cfg):
        self.cfg = cfg
        self.param_mapping = cfg.mappings
        self.start_date = cfg.time_range.start_date
        self.end_date = cfg.time_range.end_date
        self.df = None
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
    