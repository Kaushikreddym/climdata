import argparse
import climdata
from hydra import initialize, compose

def main():
    parser = argparse.ArgumentParser(
        description="Fetch CMIP data for multiple variables over a box and save merged NetCDF."
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
    parser.add_argument("--experiment-id", type=str, default="historical", help="CMIP experiment ID.")
    parser.add_argument("--source-id", type=str, default="MIROC6", help="CMIP source/model ID.")
    parser.add_argument("--table-id", type=str, default="day", help="CMIP table ID.")
    parser.add_argument("--output", type=str, default="./data/cmip_box.nc", help="Output NetCDF file name.")

    args = parser.parse_args()
    variables = [v.strip() for v in args.variables.split(",")]

    with initialize(config_path="../../conf", version_base=None):
        cfg = compose(
            config_name="config",
            overrides=[
                f"region={args.region}",
                f"time_range.start_date={args.start_date}",
                f"time_range.end_date={args.end_date}",
            ],
        )

    cmip = climdata.CMIP(
        experiment_id=args.experiment_id,
        source_id=args.source_id,
        table_id=args.table_id,
        variables=variables,
        region_bounds=cfg.bounds[cfg.region]
    )
    cmip.fetch()      # gets file lists
    cmip.load()       # loads and merges datasets

    cmip.extract(box=cmip.region_bounds)
    ds = cmip._subset_time(cfg.time_range.start_date, cfg.time_range.end_date)
    
    cmip.save_netcdf(args.output)
    print(f"Saved CMIP subset to {args.output}")

if __name__ == "__main__":
    main()