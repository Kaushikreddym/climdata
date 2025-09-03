import os
import pandas as pd
import hydra
from wetterdienst import Settings
from wetterdienst.provider.dwd.observation import DwdObservationRequest
from climdata.utils.utils_download import build_output_filename

class DWDmirror:
    def __init__(self, cfg):
        self.cfg = cfg
        self.param_mapping = cfg.mappings
        self.provider = cfg.dataset.lower()
        self.parameter_key = cfg.weather.parameter
        self.lat = cfg.location.lat
        self.lon = cfg.location.lon
        self.distance = cfg.location.buffer_km
        self.start_date = cfg.time_range.start_date
        self.end_date = cfg.time_range.end_date
        self.units = self.param_mapping[self.provider]['variables'][self.parameter_key].get("unit", None)
        self.df = None
    def fetch(self):
        param_info = self.param_mapping[self.provider]['variables'][self.parameter_key]
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
            latlon=(self.lat, self.lon),
            distance=self.distance,
            unit="km"
        )

        df = request.values.all().df.to_pandas()
        self.df = df
        return self.df

    def format(self):
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
        self.df["variable"] = self.parameter_key
        self.df["latitude"] = self.lat
        self.df["longitude"] = self.lon
        self.df['source'] = 'DWD'
        self.df['units'] = self.units
        self.df = self.df[["latitude", "longitude", "time", "source", "variable", "value", "units"]]
        # self.df = df
        return self.df

    def save(self):
        filename = build_output_filename(self.cfg)
        self.df.to_csv(self.cfg.output.out_dir+filename, index=False)
        print(f"✅ Saved time series to: {filename}")
        return filename
    