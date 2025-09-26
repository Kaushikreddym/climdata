import argparse
import climdata
from hydra import initialize, compose
import xclim
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Fetch MSWX data for multiple variables and save merged CSV."
    )
    parser.add_argument("--lat", type=float, default=52.0, help="Latitude of the point of interest.")
    parser.add_argument("--lon", type=float, default=20.0, help="Longitude of the point of interest.")
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
    parser.add_argument("--output", type=str, default="./data/mswx_loc.csv", help="Output CSV file name.")
    parser.add_argument("--service_account", type=str, default="./conf/service_client.json", help="Path to Google service account JSON file.")

    args = parser.parse_args()
    variables = [v.strip() for v in args.variables.split(",")]

    all_dfs = []

    for var in variables:
        with initialize(config_path="../../conf", version_base=None):
            cfg = compose(
                config_name="config",
                overrides=[
                    "dataset=mswx",
                    f"weather.parameter={var}",
                    f"time_range.start_date={args.start_date}",
                    f"time_range.end_date={args.end_date}",
                    f"mappings.mswx.params.google_service_account={args.service_account}",
                ],
            )
        mswx_ob = climdata.MSWX(cfg)
        mswx_ob.fetch()
        dset = mswx_ob.load()

        dset_subset = mswx_ob.extract(point=(args.lat, args.lon)).compute()
        dset_subset = xclim.core.units.convert_units_to(
            dset_subset, cfg.mappings["info"][cfg.weather.parameter].units
        )

        df = mswx_ob.to_dataframe(dset_subset)
        df = mswx_ob.format(df)

        # Rename column to variable name to avoid conflicts
        # df = df.rename(columns={df.columns[-1]: var})
        all_dfs.append(df)

    # Merge all dataframes on index (time)
    merged_df = pd.concat(all_dfs, axis=0)

    # Save to CSV
    merged_df.to_csv(args.output)
    print(f"Saved merged dataframe with variables {variables} to {args.output}")


if __name__ == "__main__":
    main()