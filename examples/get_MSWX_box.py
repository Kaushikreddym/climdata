import argparse
import climdata
from hydra import initialize, compose
import xarray as xr
import xclim

def main():
    parser = argparse.ArgumentParser(
        description="Fetch MSWX data for multiple variables over a box and save merged NetCDF."
    )
    parser.add_argument("--region", type=str, default="europe", help="Region name as in config (e.g. europe)")
    parser.add_argument(
        "--variables",
        type=str,
        default="tasmax,tasmin,pr",
        help="Comma-separated list of weather variables (e.g. tasmax,tasmin,pr).",
    )
    parser.add_argument(
        "--start-date", type=str, default="2000-08-01", help="Start date (YYYY-MM-DD)."
    )
    parser.add_argument(
        "--end-date", type=str, default="2000-08-31", help="End date (YYYY-MM-DD)."
    )
    parser.add_argument("--output", type=str, default="./data/mswx_box.nc", help="Output NetCDF file name.")
    parser.add_argument("--service_account", type=str, default="./conf/service_client.json", help="Path to Google service account JSON file.")

    args = parser.parse_args()
    variables = [v.strip() for v in args.variables.split(",")]
    ds_var = []
    for var in variables:
        with initialize(config_path="../../conf", version_base=None):
            cfg = compose(
                config_name="config",
                overrides=[
                    f"region={args.region}",
                    "dataset=mswx",
                    f"weather.parameter={var}",
                    f"time_range.start_date={args.start_date}",
                    f"time_range.end_date={args.end_date}",
                    f"mappings.mswx.params.google_service_account={args.service_account}",
                ],
            )

        mswx = climdata.MSWX(cfg)
        mswx.fetch()      # gets file lists
        mswx.load()       # loads and merges datasets

        mswx.extract(box=cfg.bounds[cfg.region])
        ds_var.append(mswx.dataset)
        
    # import ipdb;ipdb.set_trace()
    mswx.dataset = xr.merge(ds_var)
    mswx.dataset = xclim.core.units.convert_units_to(
        mswx.dataset, cfg.mappings["info"][cfg.weather.parameter].units
        )
    mswx.save_netcdf(args.output)
    print(f"Saved MSWX subset to {args.output}")

if __name__ == "__main__":
    main()