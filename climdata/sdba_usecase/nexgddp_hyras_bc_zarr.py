#!/usr/bin/env python
"""
nexgddp_hyras_bc_zarr.py
────────────────────────
Per-location multivariate bias correction of NEX-GDDP-CMIP6 against HYRAS REC2D,
reading all inputs lazily from pre-built Zarr stores (lat=10, lon=10, time=-1 chunks).

One CSV is written per grid cell:
  {OUTPUT_DIR}/{model}/{scenario}/{col_num}/{model}_{scenario}_{start}_{end}_C{col_num}R{row_num}.csv

Row index  = decreasing-latitude order (north → south, 0-based)
Col index  = increasing-longitude order (west  → east,  0-based)

Units normalised before BC:
  NEX-GDDP tas/tasmax/tasmin : K         → degC   (sub 273.15)
  NEX-GDDP pr                : kg m-2 s-1 → mm d-1 (× 86400)
  HYRAS obs                  : already in degC / mm / W m-2 / %

Output CSV columns (in order):
  time, Precipitation (mm/day), TempMin (°C), TempMean (°C), TempMax (°C),
  Radiation (kJ m-2), SunshineDuration, SoilMoisture, SoilTemperature,
  Windspeed, RefETcalc, RefETdwd, RelHumCalc (%)
  Columns without a data source are written as NaN.

Usage
-----
    python nexgddp_hyras_bc_zarr.py --model ACCESS-CM2
    python nexgddp_hyras_bc_zarr.py --model GFDL-ESM4 --scenarios ssp370 ssp126
    python nexgddp_hyras_bc_zarr.py --model CanESM5 --n-iterations 20
"""

import os
import json
import shutil
import argparse
import warnings
import time
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import xarray as xr

from climdata.sdba.bcsd import BiasCorrection

print("✓ All imports successful")


# ── Output schema ─────────────────────────────────────────────────────────────

VAR_RENAME = {
    'pr':     'Precipitation',
    'tas':    'TempMean',
    'tasmax': 'TempMax',
    'tasmin': 'TempMin',
    'rsds':   'Radiation',
    'hurs':   'RelHumCalc',
}

OUTPUT_COLUMNS = [
    'Precipitation',
    'TempMin',
    'TempMean',
    'TempMax',
    'Radiation',
    'SunshineDuration',
    'SoilMoisture',
    'SoilTemperature',
    'Windspeed',
    'RefETcalc',
    'RefETdwd',
    'RelHumCalc',
]

# Columns that have no BC source — always NaN after reindex, filled with -999.0.
DUMMY_COLUMNS = [c for c in OUTPUT_COLUMNS if c not in set(VAR_RENAME.values())]

# Fill sentinel written by adjust_bias_one_location when all values are missing.
VENDOR_FILL = 1.e20


# ── Unit conversion ───────────────────────────────────────────────────────────

def to_standard_units(ds: xr.Dataset) -> xr.Dataset:
    """
    Normalise simulation units to match HYRAS obs and ISIMIP3BASD thresholds:
      tas/tasmax/tasmin : K         → degC      (additive; threshold-agnostic)
      pr               : kg m-2 s-1 → mm d-1    (lower_threshold=0.1 mm d-1)
      rsds, hurs       : W m-2 / %  (unchanged)

    Only converts when the 'units' attribute matches the source unit.
    """
    ds = ds.copy()

    for var in ('tas', 'tasmax', 'tasmin'):
        if var not in ds:
            continue
        u = str(ds[var].attrs.get('units', '')).strip().upper()
        if u in ('K', 'KELVIN'):
            ds[var] = (ds[var] - 273.15).assign_attrs({**ds[var].attrs, 'units': 'degC'})
            print(f"   [{var}] K → degC")

    if 'pr' in ds:
        u = str(ds['pr'].attrs.get('units', '')).strip().lower()
        if 'kg' in u:
            ds['pr'] = (ds['pr'] * 86400.0).assign_attrs({**ds['pr'].attrs, 'units': 'mm d-1'})
            print("   [pr] kg m-2 s-1 → mm d-1")

    return ds


def rename_and_convert(row_data: dict) -> dict:
    """Rename CF names and convert units for CSV output.
    rsds: W m-2 → kJ m-2 (multiply by 86.4 = 86400 s/day / 1000)
    """
    out = {}
    for var, arr in row_data.items():
        if var == 'rsds':
            arr = arr * 86.4
        out[VAR_RENAME.get(var, var)] = arr
    return out


def _build_output_df(col_data: dict, time_dates: list):
    """
    Build the output DataFrame for one grid cell.
    Returns None if any core BC variable contains NaN or the vendor fill value
    (1e+20, written when all input values are missing) — caller should skip the cell.
    Check is done on raw CF arrays before unit conversion so the rsds ×86.4
    factor does not inflate the sentinel into a different magnitude.
    Dummy columns (no BC source) are written as -999.0 instead of NaN.
    """
    for arr in col_data.values():
        if np.any(np.isnan(arr)) or np.any(arr >= VENDOR_FILL):
            return None
    df = (
        pd.DataFrame({'time': time_dates, **rename_and_convert(col_data)})
        .reindex(columns=['time'] + OUTPUT_COLUMNS)
    )
    df[DUMMY_COLUMNS] = df[DUMMY_COLUMNS].fillna(-999.0)
    return df


# ── Zarr store helpers ────────────────────────────────────────────────────────

def obs_zarr_path(zarr_dir: str) -> str:
    return os.path.join(zarr_dir, "HYRAS_REC2D_obs_1951-01-01_2014-12-31.zarr")


def sim_zarr_path(zarr_dir: str, model: str, scenario: str) -> str:
    if scenario == "historical":
        return os.path.join(zarr_dir, f"NEX_GDDP_{model}_historical_1951_2014.zarr")
    return os.path.join(zarr_dir, f"NEX_GDDP_{model}_{scenario}_2015_2100.zarr")


def open_zarr(store_path: str, start_date: str, end_date: str, variables: list) -> xr.Dataset:
    """Open a Zarr store lazily, slice time, keep requested variables."""
    if not os.path.exists(store_path):
        raise FileNotFoundError(f"Zarr store not found: {store_path}")
    print(f"  📂 {os.path.basename(store_path)}")
    ds = xr.open_zarr(store_path, consolidated=True, chunks="auto")
    keep    = [v for v in variables if v in ds.data_vars]
    missing = set(variables) - set(keep)
    if missing:
        print(f"  ⚠️  Variables missing in store (skipped): {sorted(missing)}")
    ds = ds[keep].sel(time=slice(start_date, end_date))
    print(f"     dims : {dict(ds.dims)}")
    print(f"     vars : {list(ds.data_vars)}")
    return ds


# ── Per-location bias correction → CSV ───────────────────────────────────────

SPATIAL_CHUNK = 10  # matches zarr lat/lon chunk size


def run_bc_to_csv(
    bc: BiasCorrection,
    obs_ref: xr.Dataset,
    sim_ref: xr.Dataset,
    sim_fut: xr.Dataset,
    output_dir: str,
    model: str,
    scenario: str,
    n_processes: int = 1,
    debug: bool = False,
) -> None:
    """
    Iterate over SPATIAL_CHUNK×SPATIAL_CHUNK tiles of the domain.  Each tile is
    bias-corrected independently via ba.adjust_bias (vendor-parallel), CSVs are
    written, then the tile is discarded.  Memory footprint is one tile at a time
    rather than the full domain.
    """
    from climdata._vendor.isimip3basd import bias_adjustment as ba
    from climdata._vendor.isimip3basd import utility_functions as uf

    lat_name = 'lat' if 'lat' in obs_ref.coords else 'latitude'
    lon_name = 'lon' if 'lon' in obs_ref.coords else 'longitude'
    lats = obs_ref[lat_name].values
    lons = obs_ref[lon_name].values
    n_lat, n_lon = len(lats), len(lons)

    lat_to_row = {float(lat): i for i, lat in enumerate(np.sort(lats)[::-1])}
    lon_to_col = {float(lon): j for j, lon in enumerate(np.sort(lons))}

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

    start_str = str(sim_fut.time.values[0])[:10].replace('-', '')
    end_str   = str(sim_fut.time.values[-1])[:10].replace('-', '')

    if debug:
        n_processes = 1
        # Restrict to the single 10×10 tile that contains (500, 500)
        d0 = (500 // SPATIAL_CHUNK) * SPATIAL_CHUNK
        obs_ref = obs_ref.isel({lat_name: slice(d0, d0 + SPATIAL_CHUNK),
                                lon_name: slice(d0, d0 + SPATIAL_CHUNK)})
        sim_ref = sim_ref.isel({lat_name: slice(d0, d0 + SPATIAL_CHUNK),
                                lon_name: slice(d0, d0 + SPATIAL_CHUNK)})
        sim_fut = sim_fut.isel({lat_name: slice(d0, d0 + SPATIAL_CHUNK),
                                lon_name: slice(d0, d0 + SPATIAL_CHUNK)})
        lats = obs_ref[lat_name].values
        lons = obs_ref[lon_name].values
        n_lat, n_lon = len(lats), len(lons)
        print(f"   [DEBUG] tile ({d0}:{d0+SPATIAL_CHUNK}, {d0}:{d0+SPATIAL_CHUNK})")

    # Pre-build ba kwargs — same for every tile
    vkw = bc._ba_vendor_kwargs()
    ba_kwargs = dict(
        n_processes=n_processes,
        n_iterations=bc.n_iterations,
        randomization_seed=bc.randomization_seed,
        n_quantiles=bc.n_quantiles,
        detrend=vkw.pop('detrend'),
        halfwin_upper_bound_climatology=vkw.pop('halfwin_upper_bound_climatology'),
        **vkw,
    )

    n_chunks  = ((n_lat + SPATIAL_CHUNK - 1) // SPATIAL_CHUNK) * \
                ((n_lon + SPATIAL_CHUNK - 1) // SPATIAL_CHUNK)
    chunk_idx = 0
    n_done = n_skipped = 0
    time_dates = None  # extracted from first tile's iris cube, reused for all
    t0 = time.time()

    tmpdir = Path(output_dir) / '_tmp_npy'
    tmpdir.mkdir(parents=True, exist_ok=True)

    try:
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
                            / f"{model}_{scenario}_{start_str}_{end_str}"
                              f"_C{col_num}R{row_num}.csv"
                        )
                        cells.append((ci, cj, col_num, row_num, csv_path))

                pending = [c for c in cells if not c[4].exists()]
                if not pending:
                    n_skipped += len(cells)
                    continue

                # Slice zarr stores to this tile (lazy — no data loaded yet)
                obs_c  = obs_ref.isel({lat_name: slice(i0, i1), lon_name: slice(j0, j1)})
                sref_c = sim_ref.isel({lat_name: slice(i0, i1), lon_name: slice(j0, j1)})
                sfut_c = sim_fut.isel({lat_name: slice(i0, i1), lon_name: slice(j0, j1)})

                # Per-tile mask: reads only this tile's zarr chunks.
                # .all(dim=time) ensures cells with any NaN in any variable are excluded —
                # required because rsds has halfwin_upper_bound_climatology=15, and partial
                # NaN in the input causes NaN in the running-window climatology → assertion.
                _tdims = [d for d in obs_c.dims if d not in (lat_name, lon_name)]
                tile_mask = (
                    obs_c.notnull()
                    .all(dim=_tdims)
                    .to_array(dim='variable')
                    .all('variable')
                    .compute()
                )
                if not bool(tile_mask.any()):
                    n_skipped += len(cells)
                    continue

                obs_c  = obs_c.where(tile_mask)
                sref_c = sref_c.where(tile_mask)
                sfut_c = sfut_c.where(tile_mask)

                obs_cubes, sh_cubes, sf_cubes = bc._to_iris_cubes(obs_c, sref_c, sfut_c)

                if time_dates is None:
                    tc = sf_cubes[0].coord('time')
                    time_dates = [str(d)[:10] for d in tc.units.num2date(tc.points)]

                tile_shape = sf_cubes[0].shape[1:]
                n_times    = sf_cubes[0].shape[0]

                chunk_tmpdir = tmpdir / f"c_{i0}_{j0}"
                chunk_tmpdir.mkdir(exist_ok=True)
                npy_paths = [(chunk_tmpdir / f"ba_{v}.nc").resolve() for v in bc.variables]
                for p in npy_paths:
                    uf.setup_npy_stack(str(p), (n_times,) + tile_shape)

                ba.adjust_bias(
                    obs_hist=obs_cubes,
                    sim_hist=sh_cubes,
                    sim_fut=sf_cubes,
                    sim_fut_ba_path=[str(p) for p in npy_paths],
                    **ba_kwargs,
                )

                for ci, cj, col_num, row_num, csv_path in pending:
                    i_1d = np.ravel_multi_index((ci, cj), tile_shape)
                    npy_file_0 = uf.npy_stack_dir(str(npy_paths[0])) + f'{i_1d}.npy'
                    if not os.path.exists(npy_file_0):
                        continue  # vendor skipped — all-missing

                    col_data = {}
                    for var, p in zip(bc.variables, npy_paths):
                        npy_file = uf.npy_stack_dir(str(p)) + f'{i_1d}.npy'
                        col_data[var] = np.load(npy_file).squeeze()
                        os.remove(npy_file)

                    df = _build_output_df(col_data, time_dates)
                    if df is None:
                        continue

                    csv_path.parent.mkdir(parents=True, exist_ok=True)
                    df.to_csv(csv_path, index=False)
                    n_done += 1

                shutil.rmtree(str(chunk_tmpdir), ignore_errors=True)

                elapsed = time.time() - t0
                rate = n_done / elapsed if elapsed > 0 else 0
                print(
                    f"   tile {chunk_idx}/{n_chunks} [{i0}:{i1},{j0}:{j1}]"
                    f" | done={n_done:,} skipped={n_skipped:,} | {rate:.2f} loc/s"
                )

    finally:
        shutil.rmtree(str(tmpdir), ignore_errors=True)

    print(f"   ✅ {n_done:,} new CSVs | {n_skipped:,} skipped → {output_dir}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Per-location BC of NEX-GDDP vs HYRAS — Zarr input → CSV output"
    )
    parser.add_argument("--model",        type=str, default="ACCESS-CM2")
    parser.add_argument("--member",       type=str, default="r1i1p1f1")
    parser.add_argument("--variables",    type=str, nargs="+",
                        default=["pr", "tas", "tasmax", "tasmin", "rsds", "hurs"])
    parser.add_argument("--scenarios",    type=str, nargs="+",
                        default=["historical", "ssp370", "ssp126"])
    parser.add_argument("--n-iterations", type=int, default=20)
    parser.add_argument("--n-processes",  type=int, default=1,
                        help="Worker processes for the per-location BC loop (default: 1 = serial)")
    parser.add_argument("--debug",        action="store_true",
                        help="Process only iloc (500, 500) with n_processes=1")
    parser.add_argument("--dask-workers", type=int, default=None,
                        help="Dask scheduler thread count for xarray/zarr I/O (default: dask default)")
    args = parser.parse_args()

    if args.dask_workers:
        import dask
        dask.config.set(num_workers=args.dask_workers)
        print(f"  Dask workers set to {args.dask_workers}")

    MODEL     = args.model
    MEMBER    = args.member
    VARIABLES = args.variables
    SCENARIOS = args.scenarios

    ZARR_DIR   = "/data01/FDS/muduchuru/Atmos/zarr_stores"
    OUTPUT_DIR = "/data01/FDS/muduchuru/Atmos/NEXGDDP_HYRAS_BC_CSV"

    REF_START  = "1951-01-01"
    REF_END    = "1980-12-31"
    HIST_START = "1951-01-01"
    HIST_END   = "2014-12-31"
    FUT_START  = "2015-01-01"
    FUT_END    = "2100-12-31"

    print("\n" + "=" * 70)
    print(f"  Model     : {MODEL}  ({MEMBER})")
    print(f"  Variables : {VARIABLES}")
    print(f"  Scenarios : {SCENARIOS}")
    print(f"  Ref period: {REF_START} → {REF_END}")
    print(f"  Zarr dir  : {ZARR_DIR}")
    print(f"  Output    : {OUTPUT_DIR}")
    print("=" * 70 + "\n")

    t_total = time.time()

    # ── 1. Open HYRAS obs Zarr (training period only) ─────────────────────────
    print("── Step 1: Open HYRAS obs Zarr ──")
    obs_raw = open_zarr(obs_zarr_path(ZARR_DIR), HIST_START, HIST_END, VARIABLES)
    # HYRAS is already in degC / mm / W m-2 / % — no conversion needed
    obs_ref = obs_raw.sel(time=slice(REF_START, REF_END))
    print(f"  obs_ref : {dict(obs_ref.dims)}\n")
    
    # ── 2. Open sim historical Zarr (for training period) ────────────────────
    print("── Step 2: Open sim historical Zarr ──")
    sim_hist_full = open_zarr(
        sim_zarr_path(ZARR_DIR, MODEL, "historical"), HIST_START, HIST_END, VARIABLES)
    sim_ref = to_standard_units(sim_hist_full.sel(time=slice(REF_START, REF_END)))
    print(f"  sim_ref : {dict(sim_ref.dims)}\n")

    # ── 3. Per-scenario bias correction ───────────────────────────────────────
    for scenario in SCENARIOS:
        print(f"\n{'=' * 70}")
        print(f"  Scenario : {scenario}")
        print(f"{'=' * 70}")

        if scenario == "historical":
            sim_fut = to_standard_units(sim_hist_full)
        else:
            print(f"── Opening sim {scenario} Zarr ──")
            sim_fut = to_standard_units(
                open_zarr(
                    sim_zarr_path(ZARR_DIR, MODEL, scenario),
                    FUT_START, FUT_END, VARIABLES,
                )
            )
        print(f"  sim_fut : {dict(sim_fut.dims)}")
        out_dir = os.path.join(OUTPUT_DIR, MODEL, scenario)
        os.makedirs(out_dir, exist_ok=True)

        bc = BiasCorrection(
            variable           = VARIABLES,
            n_iterations       = args.n_iterations,
            n_processes        = args.n_processes,
            randomization_seed = 42,
        )

        t0 = time.time()
        
        # import ipdb; ipdb.set_trace()
        
        print(f"\n── Step 3: Bias Correction → CSV ──")
        run_bc_to_csv(bc, obs_ref, sim_ref, sim_fut, out_dir, MODEL, scenario,
                      n_processes=args.n_processes, debug=args.debug)
        print(f"  ✅ {scenario} done in {(time.time() - t0) / 60:.1f} min")

    print(f"\n{'=' * 70}")
    print(f"✅ All scenarios done in {(time.time() - t_total) / 60:.1f} min")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
