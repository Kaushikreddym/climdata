#!/usr/bin/env python
"""
hyras_obs_to_csv.py
───────────────────
Validation companion to nexgddp_hyras_bc_zarr.py.

Writes the HYRAS REC2D *observation* reference to per-location CSV files using
exactly the same grid indexing, file naming and output schema as the
bias-corrected CMIP output — but with NO bias correction applied.  This gives a
ground-truth set of CSVs that can be compared cell-for-cell against the
bias-corrected model runs.

One CSV is written per grid cell:
  {OUTPUT_DIR}/{col_num}/obs_{start}_{end}_C{col_num}R{row_num}.csv

Row index  = decreasing-latitude order (north → south, 0-based)
Col index  = increasing-longitude order (west  → east,  0-based)

HYRAS obs is already in degC / mm / W m-2 / % — the only conversion applied is
rsds (W m-2 → kJ m-2, ×86.4), identical to the BC output, so the columns are
directly comparable.

Output CSV columns match the BC script exactly (see OUTPUT_COLUMNS); columns
without an observation source are written as -999.0.

Usage
-----
    python hyras_obs_to_csv.py
    python hyras_obs_to_csv.py --start 1981-01-01 --end 2010-12-31
    python hyras_obs_to_csv.py --variables pr tas tasmax tasmin --debug
"""

import os
import json
import argparse
import warnings
import time
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import xarray as xr

# Reuse the exact schema / helpers from the BC script so the two stay in sync.
from nexgddp_hyras_bc_zarr import (
    OUTPUT_COLUMNS,
    obs_zarr_path,
    open_zarr,
    _build_output_df,
)

print("✓ All imports successful")


SPATIAL_CHUNK = 10  # matches zarr lat/lon chunk size


def write_obs_to_csv(
    obs: xr.Dataset,
    output_dir: str,
    variables: list,
    start_str: str,
    end_str: str,
    label: str = "obs",
    debug: bool = False,
) -> None:
    """
    Iterate over SPATIAL_CHUNK×SPATIAL_CHUNK tiles of the domain, load each tile
    eagerly, and write one CSV per grid cell.  No bias correction — the raw obs
    time series (rsds converted to kJ m-2) is written through the shared
    _build_output_df so the columns match the BC output exactly.
    """
    lat_name = 'lat' if 'lat' in obs.coords else 'latitude'
    lon_name = 'lon' if 'lon' in obs.coords else 'longitude'
    lats = obs[lat_name].values
    lons = obs[lon_name].values
    n_lat, n_lon = len(lats), len(lons)

    lat_to_row = {float(lat): i for i, lat in enumerate(np.sort(lats)[::-1])}
    lon_to_col = {float(lon): j for j, lon in enumerate(np.sort(lons))}

    os.makedirs(output_dir, exist_ok=True)
    json_path = Path(output_dir) / "rowcol_to_latlon.json"
    if not json_path.exists():
        row_to_lat = {v: k for k, v in lat_to_row.items()}
        col_to_lon = {v: k for k, v in lon_to_col.items()}
        mapping = {
            f"{row},{col}": [row_to_lat[row], col_to_lon[col]]
            for row in range(n_lat)
            for col in range(n_lon)
        }
        with open(json_path, "w") as fh:
            json.dump(mapping, fh)
        print(f"   Saved {json_path}")

    if debug:
        d0 = (500 // SPATIAL_CHUNK) * SPATIAL_CHUNK
        obs = obs.isel({lat_name: slice(d0, d0 + SPATIAL_CHUNK),
                        lon_name: slice(d0, d0 + SPATIAL_CHUNK)})
        lats = obs[lat_name].values
        lons = obs[lon_name].values
        n_lat, n_lon = len(lats), len(lons)
        print(f"   [DEBUG] tile ({d0}:{d0+SPATIAL_CHUNK}, {d0}:{d0+SPATIAL_CHUNK})")

    time_dates = [str(d)[:10] for d in obs.time.values]

    n_chunks = ((n_lat + SPATIAL_CHUNK - 1) // SPATIAL_CHUNK) * \
               ((n_lon + SPATIAL_CHUNK - 1) // SPATIAL_CHUNK)
    chunk_idx = 0
    n_done = n_skipped = 0
    t0 = time.time()

    for i0 in range(0, n_lat, SPATIAL_CHUNK):
        i1 = min(i0 + SPATIAL_CHUNK, n_lat)

        for j0 in range(0, n_lon, SPATIAL_CHUNK):
            j1 = min(j0 + SPATIAL_CHUNK, n_lon)
            chunk_idx += 1
            chunk_lats = lats[i0:i1]
            chunk_lons = lons[j0:j1]

            # Collect output paths; skip whole tile if every CSV exists
            cells = []
            for ci, lat in enumerate(chunk_lats):
                for cj, lon in enumerate(chunk_lons):
                    row_num = lat_to_row[float(lat)]
                    col_num = lon_to_col[float(lon)]
                    csv_path = (
                        Path(output_dir) / str(col_num)
                        / f"{label}_{start_str}_{end_str}_C{col_num}R{row_num}.csv"
                    )
                    cells.append((ci, cj, col_num, row_num, csv_path))

            pending = [c for c in cells if not c[4].exists()]
            if not pending:
                n_skipped += len(cells)
                continue

            # Load this tile eagerly (reads only this tile's zarr chunks)
            tile = obs.isel({lat_name: slice(i0, i1),
                             lon_name: slice(j0, j1)}).load()

            for ci, cj, col_num, row_num, csv_path in pending:
                col_data = {
                    var: tile[var].isel({lat_name: ci, lon_name: cj}).values
                    for var in variables if var in tile.data_vars
                }

                df = _build_output_df(col_data, time_dates)
                if df is None:
                    continue  # all-NaN / off-domain cell

                csv_path.parent.mkdir(parents=True, exist_ok=True)
                df.to_csv(csv_path, index=False)
                n_done += 1

            elapsed = time.time() - t0
            rate = n_done / elapsed if elapsed > 0 else 0
            print(
                f"   tile {chunk_idx}/{n_chunks} [{i0}:{i1},{j0}:{j1}]"
                f" | done={n_done:,} skipped={n_skipped:,} | {rate:.2f} loc/s"
            )

    print(f"   ✅ {n_done:,} new CSVs | {n_skipped:,} skipped → {output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Write HYRAS observation reference to per-location CSV "
                    "(validation companion to the NEX-GDDP BC run)"
    )
    parser.add_argument("--variables", type=str, nargs="+",
                        default=["pr", "tas", "tasmax", "tasmin", "rsds", "hurs"])
    parser.add_argument("--start", type=str, default="1951-01-01")
    parser.add_argument("--end",   type=str, default="2014-12-31")
    parser.add_argument("--label", type=str, default="obs",
                        help="Filename prefix for the obs CSVs (default: obs)")
    parser.add_argument("--debug", action="store_true",
                        help="Process only the tile containing iloc (500, 500)")
    parser.add_argument("--dask-workers", type=int, default=None,
                        help="Dask scheduler thread count for xarray/zarr I/O")
    args = parser.parse_args()

    if args.dask_workers:
        import dask
        dask.config.set(num_workers=args.dask_workers)
        print(f"  Dask workers set to {args.dask_workers}")

    VARIABLES = args.variables

    ZARR_DIR   = "/data01/FDS/muduchuru/Atmos/zarr_stores"
    OUTPUT_DIR = "/data01/FDS/muduchuru/Atmos/NEXGDDP_HYRAS_BC_CSV/OBS"

    start_str = args.start.replace('-', '')
    end_str   = args.end.replace('-', '')

    print("\n" + "=" * 70)
    print(f"  HYRAS observation reference → CSV (no bias correction)")
    print(f"  Variables : {VARIABLES}")
    print(f"  Period    : {args.start} → {args.end}")
    print(f"  Zarr dir  : {ZARR_DIR}")
    print(f"  Output    : {OUTPUT_DIR}")
    print("=" * 70 + "\n")

    t_total = time.time()

    print("── Open HYRAS obs Zarr ──")
    obs = open_zarr(obs_zarr_path(ZARR_DIR), args.start, args.end, VARIABLES)
    # HYRAS is already in degC / mm / W m-2 / % — no unit conversion before output;
    # rsds W m-2 → kJ m-2 is applied inside _build_output_df, matching the BC run.
    print(f"  obs : {dict(obs.dims)}\n")

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    write_obs_to_csv(
        obs, OUTPUT_DIR, VARIABLES, start_str, end_str,
        label=args.label, debug=args.debug,
    )

    print(f"\n{'=' * 70}")
    print(f"✅ Done in {(time.time() - t_total) / 60:.1f} min")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
