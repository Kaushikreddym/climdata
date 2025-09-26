import argparse
import climdata
from hydra import initialize, compose
import pandas as pd
import xclim

def main():
    parser = argparse.ArgumentParser(
        description="Fetch CMIP data for multiple variables and save merged NetCDF."
    )
    parser.add_argument("--lat", type=float, default=52.0, help="Latitude of the point of interest.")
    parser.add_argument("--lon", type=float, default=15.0, help="Longitude of the point of interest.")
    parser.add_argument(
        "--variables",
        type=str,
        default="tasmax,tasmin,pr",
        help="Comma-separated list of weather variables (e.g. tasmax,tasmin,pr).",
    )
    parser.add_argument(
        "--start-date", type=str, default="2000-01-01", help="Start date (YYYY-MM-DD)."
    )
    parser.add_argument(
        "--end-date", type=str, default="2000-01-31", help="End date (YYYY-MM-DD)."
    )
    parser.add_argument("--experiment-id", type=str, default="historical", help="CMIP experiment ID.")
    parser.add_argument("--source-id", type=str, default="MIROC6", help="CMIP source/model ID.")
    parser.add_argument("--table-id", type=str, default="day", help="CMIP table ID.")
    parser.add_argument("--output", type=str, default="./data/cmip_loc.csv", help="Output NetCDF file name.")

    args = parser.parse_args()
    variables = [v.strip() for v in args.variables.split(",")]

    with initialize(config_path="../../conf", version_base=None):
                cfg = compose(
                    config_name="config",
                    overrides=[
                        "dataset=cmip6",
                        f"time_range.start_date={args.start_date}",
                        f"time_range.end_date={args.end_date}",
                    ],
                )

    cmip = climdata.CMIP(
        experiment_id=args.experiment_id,
        source_id=args.source_id,
        table_id=args.table_id,
        variables=variables,
    )
    cmip.fetch()      # gets file lists
    cmip.load()       # loads and merges datasets
    cmip.extract(point=(args.lat, args.lon))
    ds = cmip._subset_time(cfg.time_range.start_date, cfg.time_range.end_date)
    for var in ds.data_vars:
        ds[var] =  xclim.core.units.convert_units_to(
            ds[var], cfg.mappings["info"][var].units
        )
    cmip.save_csv(args.output)
    print(f"Saved CMIP subset to {args.output}")

if __name__ == "__main__":
    main()
