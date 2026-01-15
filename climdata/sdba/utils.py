import os
import xarray as xr
import numpy as np
from xsdba.adjustment import DetrendedQuantileMapping as DQM
from xsdba.processing import adapt_freq 
from scipy.fft import fft, ifft
import xesmf as xe
# from dask.distributed import Client, LocalCluster
from glob import glob
import ipdb
import xarray as xr

import warnings
warnings.filterwarnings("ignore")

import argparse
from multiprocessing import Pool
from functools import partial

import re
import xclim.core.units as xunits
parser = argparse.ArgumentParser()
parser.add_argument('--variable', required=True, help='Climate variable to process (e.g. tasmax, pr)')
args = parser.parse_args()
variable = args.variable

def compute_daily_climatology(data):
    doy_clim = data.groupby("time.dayofyear").mean("time")
    smoothed = xr.apply_ufunc(
        smooth_fft, doy_clim,
        input_core_dims=[["dayofyear"]],
        output_core_dims=[["dayofyear"]],
        vectorize=True, dask="parallelized", output_dtypes=[float],
        dask_gufunc_kwargs={"allow_rechunk": True},
    )
    return smoothed

def smooth_fft(x, n_harmonics=3):
    N = x.size
    f = fft(x)
    f[n_harmonics+1:-n_harmonics] = 0
    return np.real(ifft(f))

def preserve_attrs(source, target):
    target.attrs = source.attrs
    return target

def get_realization_path(base_path, variable, gcm, scenario, year):
    search_pattern = os.path.join(
        base_path, gcm, scenario, variable,
        f'{variable}_day_{gcm}_{scenario}_*_gn_{year}.nc'
    )
    files = glob(search_pattern)
    if not files:
        raise FileNotFoundError(f"No file found for {gcm} {scenario} {year}")
    
    return files[0]  # or choose the first match, or apply filtering
def process_future_year(variable, gcm, year, scenario, data_paths):
    dqm = GLOBALS['dqm']
    regridder = GLOBALS['regridder']
    try:
        fut_file = get_realization_path(data_paths["NEX"], variable, gcm, scenario, year)
        print(f"  Processing {scenario} {year} from {fut_file}...")

        fut_raw = xr.open_dataset(fut_file)[variable]
        fut_regr = preserve_attrs(fut_raw, regridder(fut_raw)).chunk({'time': -1, 'x': 50, 'y': 50})
        fut_bc = dqm.adjust(fut_regr)

        realization = re.search(rf"{scenario}_(.*?)_gn_{year}", os.path.basename(fut_file)).group(1)
        out_file = f'{variable}_{year}_{gcm}_{scenario}_{realization}_BA.nc'
        out_path = os.path.join('./data/DQM/', gcm, scenario, variable, out_file)
        os.makedirs(os.path.dirname(out_path), exist_ok=True)

        print(f"Saving bias-corrected data to {out_path}...")
        fut_bc.to_netcdf(out_path)

    except Exception as e:
        print(f"Failed to process {scenario} {year} for {gcm}: {e}")
GLOBALS = {}

def init_worker(DQM_, REGRIDDER_):
    GLOBALS['dqm'] = DQM_
    GLOBALS['regridder'] = REGRIDDER_
