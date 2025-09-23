import os
import glob
import pandas as pd
import xarray as xr
from datetime import datetime
from typing import Optional, Dict, Union
from omegaconf import DictConfig
import warnings
from pathlib import Path
from tqdm.notebook import tqdm
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from xclim.core import units
warnings.filterwarnings("ignore", category=Warning)


class CMIPmirror:
    def __init__(self, var_cfg: DictConfig, experiments):
        self.var_cfg = var_cfg
        self.files = []
        self.dataset = None
        self.experiments = experiments

    def _subset_by_bounds(self, ds, bounds, lat_name='lat', lon_name='lon'):
        return ds.sel(
            **{
                lat_name: slice(bounds['lat_min'], bounds['lat_max']),
                lon_name: slice(bounds['lon_min'], bounds['lon_max'])
            }
        )

    def _check_lat_lon(self, ds: xr.Dataset) -> xr.Dataset:
        # Fix latitude ascending order
        if "lat" in ds.coords:
            lat = ds["lat"]
            if lat.values[0] > lat.values[-1]:  # descending
                ds = ds.sortby("lat")

        # Fix longitude range to -180 to 180
        if "lon" in ds.coords:
            lon = ds["lon"]
            lon_vals = lon.values
            if lon_vals.max() > 180:
                lon_fixed = ((lon_vals + 180) % 360) - 180
                ds = ds.assign_coords(lon=lon_fixed)
                ds = ds.sortby("lon")
        return ds

    def fetch(self, base_dir,tbl_id):
        nc_files = [
            f
            for exp in self.experiments
            for f in glob.glob(
                os.path.join(base_dir, "*/*/*", exp, f"*/{tbl_id}/*/*/*/*.nc"),
                recursive=True
            )
        ]
        rows = []
        for file_path in tqdm(nc_files, desc="Indexing CMIP6 files"):
            parts = file_path.split(os.sep)
            try:
                activity_id   = parts[6]
                institution_id = parts[7]
                source_id      = parts[8]
                experiment_id  = parts[9]
                member_id      = parts[10]
                table_id       = parts[11]
                variable_id    = parts[12]
                grid_label     = parts[13]
                version        = parts[14]
            except IndexError:
                continue

            # Extract start and end date from filename
            fname = os.path.basename(file_path)
            # Example: pr_day_MIROC6_ssp245-nat_r8i1p1f1_gn_20210101-20301231.nc
            date_part = fname.split("_")[-1].replace(".nc", "")
            start_str, end_str = date_part.split("-")
            
            if tbl_id == 'Amon':
                start_date = pd.to_datetime(start_str, format="%Y%m")
                end_date   = pd.to_datetime(end_str, format="%Y%m")
            elif tbl_id == 'day':
                start_date = pd.to_datetime(start_str, format="%Y%m%d")
                end_date   = pd.to_datetime(end_str, format="%Y%m%d")
            rows.append({
                "path": file_path,
                "activity_id": activity_id,
                "institution_id": institution_id,
                "source_id": source_id,
                "experiment_id": experiment_id,
                "member_id": member_id,
                "table_id": table_id,
                "variable_id": variable_id,
                "grid_label": grid_label,
                "version": version,
                "start_date": start_date,
                "end_date": end_date
            })

        df = pd.DataFrame(rows)
        # import ipdb; ipdb.set_trace()
        # keep only experiments that match all requested
        grouped = df.groupby(["institution_id", "source_id"])["experiment_id"].unique()
        valid_pairs = grouped[grouped.apply(lambda exps: set(self.experiments).issubset(set(exps)))].index
        df = df[df.set_index(["institution_id", "source_id"]).index.isin(valid_pairs)]

        # keep only versions with "v"
        df = df[df['version'].str.contains('v')]

        # compute file-level duration
        df["years"] = (df["end_date"] - df["start_date"]).dt.days / 365.25

        # compute total duration per dataset
        coverage = df.groupby(
            ["institution_id", "source_id", "experiment_id", "member_id", "variable_id", "grid_label"]
        ).agg(
            total_years=("years", "sum"),
            start=("start_date", "min"),
            end=("end_date", "max"),
            nfiles=("path", "count")
        ).reset_index()

        # keep only groups with ≥ 60 years
        valid_groups = coverage[coverage["total_years"] >= 60]

        # filter original dataframe
        df_filtered = df.merge(
            valid_groups,
            on=["institution_id", "source_id", "experiment_id", "member_id", "variable_id", "grid_label"],
            how="inner"
        )

        return df_filtered

    def _process_var_model(self, var, model, df_filtered,subset_experiments):
        ds_list = []
        for exp in subset_experiments:
            df_filtered_sub = df_filtered[
            (df_filtered['variable_id'] == var) &
            (df_filtered['source_id'] == model) &
            (df_filtered['experiment_id'] == exp)
            ]
            members = df_filtered_sub['member_id'].unique()
            for i,member in enumerate(members[:3]):
                df_filt = df_filtered_sub[
                    (df_filtered_sub['experiment_id'] == exp) &
                    (df_filtered_sub['member_id'] == member)
                ]
                if df_filt.empty:
                    continue

                paths = df_filt['path'].values
                ds = xr.open_mfdataset(paths, combine="by_coords", chunks={"time": 365})
                if var == "pr":
                    ds[var] = units.convert_units_to(ds[var], "mm d-1")
                elif var in ["tas", "tasmax", "tasmin"]:
                    ds[var] = units.convert_units_to(ds[var], "degC")
                ds = self._check_lat_lon(ds)
                ds_europe = self._subset_by_bounds(
                    ds,
                    self.var_cfg.bounds[self.var_cfg.region]
                )
                ds_list.append(ds_europe.expand_dims({
                    "experiment": [exp],
                    "member": [i]
                }))

        if ds_list:
            ds_list = xr.align(*ds_list, join="inner", exclude=["experiment", "member"])
            combined_ds = xr.combine_by_coords(ds_list, combine_attrs="override")
            return (var, model, combined_ds)
        else:
            return (var, model, None)

    def load(self, df_filtered, vars_of_interest, subset_experiments = ["historical", "hist-aer", "hist-GHG"]):
        data_dict = defaultdict(dict)
        var_model_pairs = list(
            df_filtered[df_filtered['variable_id'].isin(vars_of_interest)]
            [['variable_id', 'source_id']]
            .drop_duplicates()
            .itertuples(index=False, name=None)
        )

        with ProcessPoolExecutor(max_workers=4) as executor:
            futures = [
                executor.submit(self._process_var_model, var, model, df_filtered, subset_experiments)
                for var, model in var_model_pairs
            ]
            for f in futures:
                var, model, ds = f.result()
                if ds is not None:
                    data_dict[model][var] = ds.chunk({'lat': 10, 'lon': 10, 'time': -1})[var]
        self.dataset = data_dict
        return data_dict

    def to_zarr(self,dataset):
        if self.dataset is None:
            raise ValueError("No dataset loaded. Call `load()` before `to_zarr()`.")
        for var_name in self.dataset.keys():
            for mod_name in self.dataset[var_name].keys():
                ds_model = self.dataset[var_name][mod_name]
            
                dataset_name = mod_name
                region = self.var_cfg.region

                if var_name == 'pr':
                    self.dataset.attrs['units'] = 'kg m-2 s-1'
                elif var_name in ['tas', 'tasmax', 'tasmin']:
                    self.dataset.attrs['units'] = 'degC'
        
                zarr_filename = self.var_cfg.output.filename.format(
                    index=var_name,
                    dataset=dataset_name,
                    region=region,
                    start=self.var_cfg.time_range.start_date,
                    end=self.var_cfg.time_range.end_date,
                    freq='1D',
                )
                zarr_path = os.path.join(f"data/{mod_name}/", zarr_filename)
                os.makedirs(os.path.dirname(zarr_path), exist_ok=True)
        
                print(f"💾 Saving {var_name} to Zarr: {zarr_path}")
                self.dataset.to_zarr(zarr_path, mode="w")
