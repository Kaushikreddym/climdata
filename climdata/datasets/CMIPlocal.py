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

    def fetch(self, base_dir, tbl_id):
        # Collect all matching NetCDF files
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
            path_parts = Path(file_path).parts

            # Find the base_dir in path_parts to determine the offset
            base_parts = Path(base_dir).parts
            try:
                base_index = len(base_parts)
            except:
                continue

            # Ensure path has enough parts after base_dir
            if len(path_parts) < base_index + 9:
                continue

            activity_id, institution_id, source_id, experiment_id, member_id, table_id, variable_id, grid_label, version = (
                path_parts[base_index:base_index + 9]
            )

            # Extract start and end dates from filename
            fname = Path(file_path).name
            date_part = fname.split("_")[-1].replace(".nc", "")
            start_str, end_str = date_part.split("-")

            if tbl_id.lower() == 'amon':
                date_fmt = "%Y%m"
            elif tbl_id.lower() == 'day':
                date_fmt = "%Y%m%d"
            else:
                raise ValueError(f"Unknown table_id: {tbl_id}")

            start_date = pd.to_datetime(start_str, format=date_fmt)
            end_date = pd.to_datetime(end_str, format=date_fmt)

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

        # Create DataFrame
        df = pd.DataFrame(rows)
        if df.empty:
            return df

        # Keep only requested experiments per institution/source pair
        grouped = df.groupby(["institution_id", "source_id"])["experiment_id"].unique()
        valid_pairs = grouped[grouped.apply(lambda exps: set(self.experiments).issubset(exps))].index
        df = df[df.set_index(["institution_id", "source_id"]).index.isin(valid_pairs)]

        # Keep only versions containing "v"
        df = df[df['version'].str.contains('v')]

        # Compute file-level duration
        df["years"] = (df["end_date"] - df["start_date"]).dt.days / 365.25
        self.df = df
        # Compute total duration per dataset and filter for ≥60 years
        coverage = (
            df.groupby(
                ["institution_id", "source_id", "experiment_id", "member_id", "variable_id", "grid_label"]
            )
            .agg(
                total_years=("years", "sum"),
                start=("start_date", "min"),
                end=("end_date", "max"),
                nfiles=("path", "count")
            )
            .reset_index()
        )

        valid_groups = coverage[coverage["total_years"] >= 60]

        # Filter original df to keep only valid groups
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

    def to_zarr(self, dataset=None):
        if self.dataset is None:
            raise ValueError("No dataset loaded. Call `load()` before `to_zarr()`.")
        for mod_name in self.dataset.keys():
            for var_name in self.dataset[mod_name].keys():
                ds_model = self.dataset[mod_name][var_name]
            
                dataset_name = mod_name
                region = self.var_cfg.region

                if var_name == 'pr':
                    ds_model.attrs['units'] = 'kg m-2 s-1'
                elif var_name in ['tas', 'tasmax', 'tasmin']:
                    ds_model.attrs['units'] = 'degC'
        
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
                ds_model.to_zarr(zarr_path, mode="w")
