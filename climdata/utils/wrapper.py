# configuration imports
from hydra import initialize, compose
from hydra.core.global_hydra import GlobalHydra

from omegaconf import DictConfig
# climdata imports
import climdata
from climdata.utils.config import _ensure_local_conf
from climdata.utils.utils_download import get_output_filename
from climdata.extremes.indices import extreme_index
from pathlib import Path
from warnings import warn

#data processing
import xarray as xr
import xclim
import pandas as pd

# system imports
import os

import json
from shapely.geometry import shape, Polygon, Point

def preprocess_aoi(cfg):
    """
    Normalize AOI in cfg:
    - Converts string AOI → dict
    - Extracts shapely geometry from FeatureCollection, Feature, or raw geometry
    - Detects type: point, bbox, or polygon
    - Updates cfg with:
        cfg.aoi_type: "point" | "bbox" | "polygon"
        cfg.lat, cfg.lon: for point
        cfg.bounds: [minx, miny, maxx, maxy] for bbox/polygon
        cfg.geometry: shapely geometry object
    """
    # -----------------------------
    # 1) Load AOI value from string
    # -----------------------------
    if not hasattr(cfg, "aoi") or cfg.aoi is None:
        return cfg

    if isinstance(cfg.aoi, str):
        try:
            cfg.aoi = json.loads(cfg.aoi)
        except json.JSONDecodeError:
            raise ValueError(f"AOI string is not valid JSON: {cfg.aoi}")

    aoi = cfg.aoi

    # -----------------------------
    # 2) Extract shapely geometry
    # -----------------------------
    if aoi.get("type") == "FeatureCollection":
        if not aoi.get("features"):
            raise ValueError("FeatureCollection contains no features")
        geom = shape(aoi["features"][0]["geometry"])
    elif aoi.get("type") == "Feature":
        geom = shape(aoi["geometry"])
    elif "type" in aoi and "coordinates" in aoi:
        geom = shape(aoi)
    else:
        raise ValueError(f"Unsupported AOI format: {aoi}")
    # -----------------------------
    # 3) Determine AOI type
    # -----------------------------
    if isinstance(geom, Point):
        # cfg.aoi_type = "point"
        cfg.lat = geom.y
        cfg.lon = geom.x
        cfg.bounds = None
    elif isinstance(geom, Polygon):
        # Check if axis-aligned bbox
        coords = list(geom.exterior.coords)
        is_bbox = len(coords) == 5 and len({c[0] for c in coords}) == 2 and len({c[1] for c in coords}) == 2

        # if is_bbox:
        #     cfg.aoi_type = "bbox"
        # else:
        #     cfg.aoi_type = "polygon"

        minx, miny, maxx, maxy = geom.bounds
        cfg.bounds['custom'] = {
                "lat_min": miny,
                "lat_max": maxy,
                "lon_min": minx,
                "lon_max": maxx
            }
        cfg.region='custom'
        cfg.lat = None
        cfg.lon = None
    else:
        raise ValueError(f"Unsupported geometry type: {geom.geom_type}")

    # -----------------------------
    # 4) Store geometry itself
    # -----------------------------
    # cfg.shapefile = geom

    return cfg



def extract_data(cfg_name: str = "config", overrides: list = None, save_to_file = True) -> str:
    overrides = overrides or []
    
    # 1. Ensure local configs are available
    conf_dir = _ensure_local_conf()  # copies conf/ to cwd
    rel_conf_dir = os.path.relpath(conf_dir, os.path.dirname(__file__))

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
    cfg = preprocess_aoi(cfg)
    extract_kwargs = {}
    filename = None
    # Determine extraction type
    if cfg.lat is not None and cfg.lon is not None:
        extract_kwargs["point"] = (cfg.lon, cfg.lat)
        if cfg.dataset=="dwd":
            extract_kwargs["buffer_km"] = 30
        filename = get_output_filename(cfg, output_type="csv", lat=cfg.lat, lon=cfg.lon)
        filename_index = get_output_filename(cfg, output_type="csv", lat=cfg.lat, lon=cfg.lon, param=cfg.index)
    elif cfg.region is not None:
        extract_kwargs["box"] = cfg.bounds[cfg.region]
        filename = get_output_filename(cfg, output_type="nc")
        filename_index = get_output_filename(cfg, output_type="nc", param=cfg.index)
    elif cfg.shapefile is not None:
        extract_kwargs["shapefile"] = cfg.shapefile
        filename = get_output_filename(cfg, output_type="nc", shp_name=cfg.shp_name)
        filename_index = get_output_filename(cfg, output_type="nc", shp_name=cfg.shp_name, param=cfg.index)

    dataset_upper = cfg.dataset.upper()

    if dataset_upper == "MSWX":
        ds_vars = []
        for var in cfg.variables:
            mswx = climdata.MSWX(cfg)
            mswx.extract(**extract_kwargs)
            mswx.load(var) 
            ds_vars.append(mswx.dataset)
        ds = xr.merge(ds_vars)
        for var in ds.data_vars:
            ds[var] = xclim.core.units.convert_units_to(ds[var], cfg.varinfo[var].units)
        if save_to_file:
            if filename.endswith(".nc"):
                ds.to_netcdf(filename)
            else:
                mswx.dataset = ds
                mswx.save_csv(filename)

    elif dataset_upper == "CMIP":
        ds_vars = []
        cmip = climdata.CMIP(cfg)
        cmip.fetch()
        cmip.load()
        cmip.extract(**extract_kwargs)
        ds = cmip.ds
        for var in ds.data_vars:
            ds[var] = xclim.core.units.convert_units_to(ds[var], cfg.varinfo[var].units)
        if save_to_file:  
            if filename.endswith(".nc"):
                cmip.save_netcdf(filename)
            else:
                cmip.save_csv(filename)

    elif dataset_upper == 'POWER':
        if "box" in extract_kwargs:
            raise ValueError("Region extraction is not supported for DWD. Please provide point.")
        
        power = climdata.POWER(cfg)
        power.fetch()
        power.load()

        ds = power.ds

        for var in ds.data_vars:
            ds[var] = xclim.core.units.convert_units_to(
                ds[var],
                cfg.varinfo[var].units
            )
        if save_to_file:
            if filename.endswith(".nc"):
                ds.to_netcdf(filename)
            else:
                df = ds.to_dataframe().reset_index()
                id_vars = [c for c in ("time", "latitude", "longitude") if c in df.columns]
                value_vars = [v for v in cfg.variables if v in df.columns]
                if not value_vars:
                    # fallback to any non-id columns as variables
                    value_vars = [c for c in df.columns if c not in id_vars]

                df_long = df.melt(
                    id_vars=id_vars,
                    value_vars=value_vars,
                    var_name="variable",
                    value_name="value",
                )

                # attach units from dataset attrs when available
                units_map = {
                    v: (ds[v].attrs.get("units", "unknown") if v in getattr(ds, "data_vars", {}) else "unknown")
                    for v in value_vars
                }
                df_long["units"] = df_long["variable"].map(units_map)
                df_long["source"] = "NASA-POWER"

                # reorder columns
                cols = ["source", "time", "latitude", "longitude", "variable", "value", "units"]
                df_long = df_long[[c for c in cols if c in df_long.columns]]

                df_long.to_csv(filename, index=False)
    elif dataset_upper == "DWD":
        # if "box" in extract_kwargs:
            # raise ValueError("Region extraction is not supported for DWD. Please provide point.")
        ds_vars = []
        for var in cfg.variables:
            dwd = climdata.DWD(cfg)
            extract_kwargs['variable'] = var
            ds = dwd.extract(**extract_kwargs)
            # dwd.format(var, lat_val, lon_val)
            ds_vars.append(ds)
        ds = xr.merge(ds_vars)
        if save_to_file:
            dwd.df = ds.mean("station_id").to_dataframe()
            dwd.save_csv(filename)

    elif dataset_upper == "HYRAS":
        hyras = climdata.HYRAS(cfg)
        ds_vars = []
        for var in cfg.variables:
            hyras.extract(**extract_kwargs)
            ds = hyras.load(var)
            ds_vars.append(ds[[var]])
        ds = xr.merge(ds_vars, compat="override")
        hyras.dataset = ds
        if save_to_file:
            if filename.endswith(".nc"):
                ds.to_netcdf(filename)
            else:
                df = ds.to_dataframe().reset_index()
                id_vars = [c for c in ("time", "lat", "lon") if c in df.columns]
                value_vars = [v for v in cfg.variables if v in df.columns]
                if not value_vars:
                    # fallback to any non-id columns as variables
                    value_vars = [c for c in df.columns if c not in id_vars]

                df_long = df.melt(
                    id_vars=id_vars,
                    value_vars=value_vars,
                    var_name="variable",
                    value_name="value",
                )

                # attach units from dataset attrs when available
                units_map = {
                    v: (ds[v].attrs.get("units", "unknown") if v in getattr(ds, "data_vars", {}) else "unknown")
                    for v in value_vars
                }
                df_long["units"] = df_long["variable"].map(units_map)
                df_long["source"] = "HYRAS"

                # reorder columns
                cols = ["source", "time", "lat", "lon", "variable", "value", "units"]
                df_long = df_long[[c for c in cols if c in df_long.columns]]

                df_long.to_csv(filename, index=False)
    ds = ds.compute()
    index = None
    if cfg.index is not None:
        try:
            # ----------------------------------------
            # CHECK FOR MINIMUM 30 YEARS OF DATA
            # ----------------------------------------
            if "time" not in ds.coords:
                warn(
                    "Dataset has no 'time' coordinate; cannot determine number of years.",
                    UserWarning
                )
            else:
                time_index = pd.to_datetime(ds.time.values)
                n_years = len(pd.unique(time_index.year))

                if n_years < 30:
                    warn(
                        f"Index '{cfg.index}' typically requires ≥30 years of data, "
                        f"but dataset contains only {n_years} years. Proceeding anyway.",
                        UserWarning
                    )

            # ----------------------------------------
            # Compute the index
            # ----------------------------------------
            indices = extreme_index(cfg, ds)

            print(f"Calculating index: {cfg.index}")
            index = indices.calculate(cfg.index).compute()

            if index is None:
                raise ValueError(f"Index calculation returned None for '{cfg.index}'")

            # Prepare output path
            out_path = Path(f"{filename_index}")
            out_path.parent.mkdir(parents=True, exist_ok=True)

            # Save to file
            if save_to_file:
                if filename_index.endswith(".nc"):
                    index.to_netcdf(out_path)
                else:
                    df = index.to_dataframe().reset_index()
                    df.to_csv(out_path)

                print(f"Saved index to: {out_path}")

        except Exception as e:
            print(f"❌ Failed to compute or save index '{cfg.index}': {e}")

    else:
        print("ℹ️ No index selected (cfg.index is None). Skipping index computation.")
    if save_to_file:
        print(f"✅ Saved output to {filename}")
        return {'cfg':cfg, 'filename': filename, 'dataset':ds, 'index_filename': filename_index, 'index': index}
    else:
        return {'cfg':cfg, 'dataset':ds, 'index': index}
    

import pandas as pd
import xarray as xr
from hydra import initialize, compose
from dask.distributed import Client
import climdata

climdata.utils.config._ensure_local_conf()

def extract_index(
    csv_path: str,
    cfg_name: str = "config",
    extra_overrides: list = None,       # Only additional overrides (e.g., "index=tx90p")
    save_to_file: bool = True,
    output_path: str = "output_index.csv",
):
    """
    Load climate CSV → generate dataset + time_range overrides automatically →
    convert to xarray → compute extreme index.
    """

    # ---- Load CSV ----
    df = pd.read_csv(csv_path, parse_dates=["time"])

    # ---- Always extract dataset + time_range from CSV ----
    dataset = df["source"].unique()[0]
    time_range_start = df["time"].min().strftime("%Y-%m-%d")
    time_range_end   = df["time"].max().strftime("%Y-%m-%d")

    auto_overrides = [
        f"dataset={dataset}",
        f"time_range.start_date={time_range_start}",
        f"time_range.end_date={time_range_end}",
    ]

    # ---- Merge with any user-provided overrides ----
    # (e.g., "index=tn10p")
    if extra_overrides:
        overrides = auto_overrides + extra_overrides
    else:
        overrides = auto_overrides + ["index=tn10p"]  # default index

    print("Using Hydra overrides:", overrides)
    
    # 1. Ensure local configs are available
    conf_dir = _ensure_local_conf()  # copies conf/ to cwd
    rel_conf_dir = os.path.relpath(conf_dir, os.path.dirname(__file__))

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
    # ---- Convert CSV → xarray ----
    df_pivot = df.pivot_table(
        index=["time", "lat", "lon"],
        columns="variable",
        values="value"
    ).reset_index()

    ds = df_pivot.set_index(["time", "lat", "lon"]).to_xarray()

    # ---- Attach units ----
    for var in df["variable"].unique():
        units = df[df["variable"] == var]["units"].iloc[0]
        ds[var].attrs["units"] = units

    print(ds)

    # ---- Compute extreme index ----
    indices = climdata.extreme_index(cfg, ds)
    print(f"Calculating index: {cfg.index}")

    result = indices.calculate(cfg.index).compute()

    # ---- Save result ----
    if save_to_file:
        result.to_dataframe().reset_index().to_csv(output_path)
        print(f"Saved result → {output_path}")

    return result



