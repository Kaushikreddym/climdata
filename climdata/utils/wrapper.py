# configuration imports
from hydra import initialize, compose
from hydra.core.global_hydra import GlobalHydra

from omegaconf import DictConfig
# climdata imports
import climdata
from climdata.utils.config import _ensure_local_conf
from climdata.utils.utils_download import get_output_filename
#data processing
import xarray as xr
import xclim
import pandas as pd

# system imports
import os

# overrides=[
#     "dataset=CMIP",
#     "lat=50",
#     "lon=10",
#     "variables=[tasmax,tasmin]"
# ]
def extract_data(cfg_name: str = "config", overrides: list = None) -> str:
    overrides = overrides or []

    # 1. Ensure local configs are available
    conf_dir = _ensure_local_conf()  # copies conf/ to cwd
    rel_conf_dir = os.path.relpath(conf_dir, os.path.dirname(__file__))
    print(rel_conf_dir)
    # 2. Initialize Hydra only if not already initialized
    if not GlobalHydra.instance().is_initialized():
        hydra_context = initialize(config_path=rel_conf_dir, version_base=None)
    else:
        # If already initialized, just set context to None for clarity
        hydra_context = None

    # Use compose within context manager if newly initialized
    if hydra_context is not None:
        with hydra_context:
            cfg: DictConfig = compose(config_name=cfg_name, overrides=overrides)
    else:
        # Already initialized: compose directly
        cfg: DictConfig = compose(config_name=cfg_name, overrides=overrides)
    extract_kwargs = {}
    filename = None
    # Determine extraction type
    if cfg.lat is not None and cfg.lon is not None:
        extract_kwargs["point"] = (cfg.lon, cfg.lat)
        filename = get_output_filename(cfg, output_type="csv", lat=cfg.lat, lon=cfg.lon)
    elif cfg.region is not None:
        extract_kwargs["box"] = cfg.bounds[cfg.region]
        filename = get_output_filename(cfg, output_type="nc")
    elif cfg.shapefile is not None:
        extract_kwargs["shapefile"] = cfg.shapefile
        filename = get_output_filename(cfg, output_type="nc", shp_name=cfg.shp_name)

    dataset_upper = cfg.dataset.upper()

    if dataset_upper == "MSWX":
        ds_vars = []
        for var in cfg.variables:
            mswx = climdata.MSWX(cfg)
            mswx.load(var)
            mswx.extract(**extract_kwargs)
            ds_vars.append(mswx.dataset)
        ds = xr.merge(ds_vars)
        for var in ds.data_vars:
            ds[var] = xclim.core.units.convert_units_to(ds[var], cfg.varinfo[var].units)
        if filename.endswith(".nc"):
            ds.to_netcdf(filename)
        else:
            mswx.dataset = ds
            mswx.save_csv(filename)

    elif dataset_upper == "CMIP":
        ds_vars = []
        print(cfg.variables)
        cmip = climdata.CMIP(cfg)
        cmip.fetch()
        cmip.load()
        cmip.extract(**extract_kwargs)
        ds = cmip.ds
        for var in ds.data_vars:
            ds[var] = xclim.core.units.convert_units_to(ds[var], cfg.varinfo[var].units)
        if filename.endswith(".nc"):
            cmip.save_netcdf(filename)
        else:
            cmip.save_csv(filename)
    elif dataset_upper == "DWD":
        # if "box" in extract_kwargs:
        #     raise ValueError("Region extraction is not supported for DWD. Please provide lat and lon.")
        ds_vars = []
        for var in cfg.variables:
            dwd = climdata.DWD(cfg)
            lat_val, lon_val = lat, lon
            extract_kwargs['variable'] = var
            ds = dwd.extract(**extract_kwargs)
            # dwd.format(var, lat_val, lon_val)
            ds_vars.append(ds)
        ds = xr.merge(ds_vars)
        # dwd.df = df
        # dwd.save_csv(filename)

    elif dataset_upper == "HYRAS":
        hyras = climdata.HYRAS(cfg)
        ds_vars = []
        for var in cfg.variables:
            hyras.load(var)
            ds = hyras.extract(**extract_kwargs)
            ds_vars.append(ds[[var]])
        ds = xr.merge(ds_vars, compat="override")
        hyras.dataset = ds
        if filename.endswith(".nc"):
            ds.to_netcdf(filename)
        else:
            hyras.save_csv(filename)

    print(f"✅ Saved output to {filename}")
    
    return cfg, filename, ds