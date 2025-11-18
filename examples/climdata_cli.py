import climdata
import xarray as xr
import xclim
import pandas as pd
import hydra
from omegaconf import DictConfig
from climdata.utils.utils_download import get_output_filename

## uncomment the below snippet for parallel processing 
# import dask
# from dask.distributed import Client

# # Configure Dask
# client = Client(
#     n_workers=20,        # or match number of physical cores
#     threads_per_worker=2,
#     memory_limit="10GB"  # per worker (8 * 10GB = 80GB total)
# )
# from multiprocessing import freeze_support

@hydra.main(config_path="../climdata/conf", config_name="config", version_base=None)
def main(cfg: DictConfig) -> None:
    
    # Extraction box or point
    extract_kwargs = {}
    filename = None
    # import ipdb; ipdb.set_trace()
    # Use config values (CLI overrides work automatically with Hydra)
    if cfg.lat is not None and cfg.lon is not None:
        extract_kwargs["point"] = (cfg.lon, cfg.lat)
        filename = get_output_filename(cfg, output_type="csv", lat=cfg.lat, lon=cfg.lon)
    elif cfg.region is not None:
        extract_kwargs["box"] = cfg.bounds[cfg.region]
        filename = get_output_filename(cfg, output_type="nc")
    elif cfg.shapefile is not None:
        extract_kwargs["shapefile"] = cfg.shapefile
        filename = get_output_filename(cfg, output_type="nc",shp_name=cfg.shp_name)

    # MSWX
    if cfg.dataset.upper() == "MSWX":
        ds_vars = []
        for var in cfg.variables:
            mswx = climdata.MSWX(cfg)
            mswx.load(var)
            mswx.extract(**extract_kwargs)
            ds_vars.append(mswx.dataset)
        ds = xr.merge(ds_vars)
        for var in ds.data_vars:
            ds[var] = xclim.core.units.convert_units_to(
                ds[var], cfg.varinfo[var].units
            )
        if filename.endswith(".nc"):
            ds.to_netcdf(filename)
        else:
            mswx.dataset = ds
            mswx.save_csv(filename)

    # CMIP
    elif cfg.dataset.upper() == "CMIP":
        ds_vars = []
        for var in cfg.variables:
            cmip = climdata.CMIP(cfg)
            cmip.fetch()
            cmip.load()
            cmip.extract(**extract_kwargs)
            ds_vars.append(cmip.ds)
        ds_merged = xr.merge(ds_vars)
        for var in ds_merged.data_vars:
            ds_merged[var] = xclim.core.units.convert_units_to(
                ds_merged[var], cfg.mappings["info"][var].units
            )
        cmip.ds = ds_merged
        if filename.endswith(".nc"):
            cmip.save_netcdf(filename)
        else:
            cmip.save_csv(filename)

    # DWD
    elif cfg.dataset.upper() == "DWD":
        if "box" in extract_kwargs:
            raise ValueError(
                "Region extraction is not supported for DWD. Please provide lat and lon for point extraction."
            )
        df_vars = []
        for var in cfg.variables:
            dwd = climdata.DWD(cfg)
            lat, lon = cfg.lat, cfg.lon
            dwd.load(var, lat, lon)
            dwd.format(var, lat, lon)
            df_vars.append(dwd.df)
        df = pd.concat(df_vars, axis=0)
        dwd.df = df
        dwd.save_csv(filename)

    # HYRAS
    elif cfg.dataset.upper() == "HYRAS":
        hyras = climdata.HYRAS(cfg)
        ds_vars = []
        for var in cfg.variables:
            hyras.load(var)
            ds = hyras.extract(**extract_kwargs)
            ds_vars.append(ds[[var]])
        ds = xr.merge(ds_vars,compat="override")
        hyras.dataset = ds  # assign for saving
        
        if filename.endswith(".nc"):
            ds.to_netcdf(filename)
        else:
            hyras.save_csv(filename)
    print(f"✅ Saved output to {filename}")

if __name__ == "__main__":
    main()
