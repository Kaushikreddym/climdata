import requests
import pandas as pd
import xarray as xr
from datetime import datetime

class POWER:
    def __init__(self, cfg):
        self.cfg = cfg
        self.raw = None
        self.ds = None

    # ---------------------------------------------------------
    # Step 1: Fetch data from NASA POWER
    # ---------------------------------------------------------
    def fetch(self):
        lat = self.cfg.lat
        lon = self.cfg.lon

        params = self.cfg.variables
        api_vars =  ",".join([
            self.cfg.dsinfo.power.variables[v].api_id
            for v in params
        ])
        start = self.cfg.time_range.start_date.replace("-", "")
        end   = self.cfg.time_range.end_date.replace("-", "")

        url = (
            "https://power.larc.nasa.gov/api/temporal/daily/point"
            f"?parameters={api_vars}"
            f"&community=AG"
            f"&latitude={lat}"
            f"&longitude={lon}"
            f"&start={start}"
            f"&end={end}"
            f"&format=JSON"
        )

        r = requests.get(url)
        r.raise_for_status()
        self.raw = r.json()

    # ---------------------------------------------------------
    # Step 2: Load JSON into xarray Dataset
    # ---------------------------------------------------------
    def load(self):
        data = self.raw["properties"]["parameter"]

        df = pd.DataFrame(data)
        df.index = pd.to_datetime(df.index, format="%Y%m%d")
        var_map = {
            v.api_id: cmip_id
            for cmip_id, v in self.cfg.dsinfo.power.variables.items()
        }
        df = df.rename(columns=var_map)
        self.ds = xr.Dataset.from_dataframe(df)

        # Add coords
        self.ds = self.ds.assign_coords(
            latitude=self.cfg.lat,
            longitude=self.cfg.lon
        )
        self.ds = self.ds.rename({"index":"time"})
        self.ds.latitude.attrs["units"] = "degrees_north"
        self.ds.longitude.attrs["units"] = "degrees_east"

        for cmip_id in self.ds.data_vars:
            if cmip_id in self.cfg.dsinfo.power.variables:
                vinfo = self.cfg.dsinfo.power.variables[cmip_id]

                self.ds[cmip_id].attrs.update({
                    "long_name": vinfo.long_name,
                    "units": vinfo.units,
                    "source": "NASA POWER",
                    "api_id": vinfo.api_id,
                })

    # ---------------------------------------------------------
    # Step 3: Extract (temporal subsetting etc.)
    # ---------------------------------------------------------
    def extract(self, start=None, end=None):
        if start or end:
            self.ds = self.ds.sel(time=slice(start, end))

    # ---------------------------------------------------------
    # Step 4: Save methods (same API as CMIP)
    # ---------------------------------------------------------
    def save_netcdf(self, filename):
        self.ds.to_netcdf(filename)

    def save_csv(self, filename):
        self.ds.to_dataframe().to_csv(filename)
