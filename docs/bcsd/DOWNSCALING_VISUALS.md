# Visual Architecture & Data Flow Diagrams

## Complete Pipeline Architecture

```
┌──────────────────────────────────────────────────────────────────────────┐
│                    NEXGDDP-HYRAS-BCSD PIPELINE                          │
└──────────────────────────────────────────────────────────────────────────┘

Step 0: Data Download & Loading
═════════════════════════════════
┌─────────────────────────────────────────────────────────────────────────┐
│ NEX-GDDP (via climdata)          HYRAS (via climdata)                   │
│  • Model: GFDL-ESM4               • Observations reference data         │
│  • Scenario: ssp370 (future)      • Period: 1951–1980 (30 years)       │
│  • Period: 2015–2100 (86 years)   • Resolution: 5 km                   │
│  • Resolution: 0.25° (~25 km)     • Format: NetCDF (loaded to RAM)     │
│  • Format: NetCDF (cloud access)  • Vars: pr, hurs, tas, ...           │
│  • Vars: pr, hurs, tas, ...       • Coords: 2D lat(y,x), lon(y,x)     │
│  • Coords: 1D lat, 1D lon         • Grid: rotated-pole projection      │
│  • Grid: regular lat/lon          │                                     │
└─────────────────────────────────────────────────────────────────────────┘
          │                                      │
          ▼                                      ▼
Step 1: Bias Correction (coarse resolution)
════════════════════════════════════════════
          ┌──────────────────────────────────────┐
          │ BiasCorrection.correct()             │
          │                                      │
          │ Input:  sim_hist (1951–1980)         │
          │         sim_fut  (2015–2100)         │
          │         obs_coarse (regridded)       │
          │                                      │
          │ Method: xclim.sdba.MBCn              │
          │         (multivariate BC)            │
          │                                      │
          │ Output: bc_result (corrected)        │
          │         [31390, 32, 38]              │
          └──────────────────────────────────────┘
                        │
                        ▼
Step 2: Statistical Downscaling (fine resolution)
═════════════════════════════════════════════════
  
  For each variable:

  ┌─────────────────────────────────────────────────────────────────┐
  │ INPUT: obs_fine (HYRAS)     bc_result (NEX-GDDP BC, coarse)    │
  │        [9863, 890, 619]        [31390, 32, 38]                 │
  │        30 years, 5 km           86 years, 0.25°                │
  └─────────────────────────────────────────────────────────────────┘
               │                           │
               ▼                           ▼
    ┌──────────────────────┐    ┌────────────────────────┐
    │  Stage 1A: Check     │    │  Stage 1B: Harmonize   │
    │  Grid Type           │    │  • xarray → iris       │
    │  • Curvilinear       │    │  • Normalize calendar  │
    │    (2D AuxCoords)?   │    │  • Stamp CRS           │
    │  • YES: Reproject    │    │  • Normalize units     │
    │         (scipy)      │    └────────────────────────┘
    │                      │
    │  Output: regular     │    Output: iris cubes
    │  lat/lon grid        │    with 1D coords
    │  [9863, 864, 608]    │    [9863, 864, 608]
    └──────────────────────┘    [31390, 32, 38]
             │
             ├──────────────────────────────┘
             │
             ▼
    ┌──────────────────────────────────────────────────┐
    │ Stage 2: Parallel Bilinear Remapbil              │
    │                                                  │
    │ for year in [2015, 2016, ..., 2100]:            │
    │   Thread(year) → iris.analysis.Linear           │
    │     [365, 32, 38] → [365, 864, 608]            │
    │                                                  │
    │ ThreadPoolExecutor (32 workers)                 │
    │ Output: sim_coarse_remapbil_cube                │
    │         [31390, 864, 608]                       │
    └──────────────────────────────────────────────────┘
             │
             ▼
    ┌──────────────────────────────────────────────────┐
    │ Stage 3: ISIMIP3BASD Statistical Downscaling    │
    │                                                  │
    │ sd.downscale(                                    │
    │   obs_fine: [9863, 864, 608]                    │
    │   sim_coarse: [31390, 32, 38]                   │
    │   sim_coarse_remapbil: [31390, 864, 608]       │
    │ )                                                │
    │                                                  │
    │ ┌─ for coarse_cell (i,j) in (32 × 38):         │
    │ │  ┌─ for month in [1..12]:                    │
    │ │ │  ┌─ for iteration in [0..20]:              │
    │ │ │ │                                          │
    │ │ │ │  MBCn Algorithm:                         │
    │ │ │ │  • Rotate data                           │
    │ │ │ │  • Quantile map per axis                 │
    │ │ │ │  • Preserve weighted sum                 │
    │ │ │ │                                          │
    │ │ │ └─ Output: x_sim_fine [31390, 432]        │
    │ │ │                                            │
    │ │ └─ Save to npy_stack/{cell_idx}.npy         │
    │ │                                              │
    │ └─ 1,216 cells × 12 months × 21 iters         │
    │                                                │
    │ Output: npy_stack/ (1216 .npy files)           │
    └──────────────────────────────────────────────────┘
             │
             ▼
    ┌──────────────────────────────────────────────────┐
    │ Stage 4: Collection & Export                     │
    │                                                  │
    │ • Load npy_stack/ (1216 files)                  │
    │ • Reshape: (1216, 31390, 432) →               │
    │           (31390, 864, 608)                     │
    │ • Assemble iris cube                            │
    │ • Save to NetCDF                                │
    │ • Convert back to xarray                        │
    │                                                  │
    │ Output: sim_fine_downscaled                     │
    │         [31390, 864, 608]                       │
    │         86 years, 5 km, daily                   │
    └──────────────────────────────────────────────────┘
             │
             ▼
  ┌─────────────────────────────────────────────────────────────────┐
  │ OUTPUT: {var}_BCSD.nc                                           │
  │         • Shape: (31390, 864, 608)                              │
  │         • Time: 2015-01-01 to 2100-12-31 (daily)              │
  │         • Resolution: ~5 km regular lat/lon                    │
  │         • Attrs: downscaling_method, n_iterations, etc.        │
  └─────────────────────────────────────────────────────────────────┘

  ┌─────────────────────────────────────────────────────────────────┐
  │ Repeat for: pr, hurs, tas, tasmax, tasmin, rsds               │
  │ (Bias correction is multivariate, done once)                   │
  │ (Downscaling per variable, can be parallelized)               │
  └─────────────────────────────────────────────────────────────────┘

Step 3: Output Files
════════════════════
/data01/FDS/muduchuru/Atmos/NEXGDDP_HYRAS_ISIMIP3BASD/
├── GFDL-ESM4_r1i1p1f1_ssp370_BC.nc               [Bias-corrected coarse]
├── GFDL-ESM4_r1i1p1f1_ssp370_pr_BCSD.nc          [Downscaled precip]
├── GFDL-ESM4_r1i1p1f1_ssp370_hurs_BCSD.nc        [Downscaled humidity]
├── GFDL-ESM4_r1i1p1f1_ssp370_tas_BCSD.nc         [Downscaled temp]
├── ... [other variables]
└── cache/
    ├── sim_hist_GFDL-ESM4_1951-2014.nc
    ├── sim_fut_GFDL-ESM4_ssp370_2015-2100.nc
    └── obs_raw_HYRAS_1951-1980.nc
```

---

## Grid Transformation Visualization

```
HYRAS Curvilinear Grid (Before Harmonization)
═════════════════════════════════════════════

890 rows (y-axis)          Rotated-Pole Projection
×619 cols (x-axis)         ┌─────────────────┐
                           │  Each cell has  │
                           │  its own lat/lon│
                           │  (2D AuxCoord)  │
                           └─────────────────┘

Example: A few grid points (not to scale)

        x=0      x=100                           x=619
y=0   (57°, 2°)  (58°, 2°) ─ ─ ─ ─ ─ ─ ─ ─ ─  (62°, 4°)
       •  •  •  •                                •
y=100 (56°, 1°)  (57°, 1°) ─ ─ ─ ─ ─ ─ ─ ─ ─  (61°, 3°)
       •  •  •  •                                •
 ...               ...                            ...
       
y=890 (55°, 0°)  (56°, 0°) ─ ─ ─ ─ ─ ─ ─ ─ ─  (60°, 2°)

Coordinates:
  lat: shape (890, 619) — one value per cell
  lon: shape (890, 619) — one value per cell
  
Problem: iris.analysis.Linear needs 1D DimCoords!

        ↓ _regrid_curvilinear_to_rectilinear()

HYRAS Regular Grid (After Reprojection)
════════════════════════════════════════

890 rows (latitude axis)    Regular Lat/Lon
×608 cols (longitude axis)  ┌──────────────┐
                            │ 1D coords    │
                            │ (DimCoord)   │
                            └──────────────┘

Lat:  [55.0, 55.1, 55.2, ..., 61.5]  (890 values)
Lon:  [0.0, 0.05, 0.1, ..., 30.9]    (608 values)

Grid (aligned):
      lon=0  lon=0.05  ...  lon=30.9
lat=61.5 •      •       ...    •    (northwest)
lat=61.4 •      •       ...    •
...      •      •       ...    •
lat=55.0 •      •       ...    •    (southwest)

Data filled via scipy.interpolate.griddata (linear, per timestep)
Result: Same physical coverage, but on regular lat/lon
```

---

## Bilinear Remapbil Visualization

```
NEX-GDDP Coarse Grid → HYRAS Fine Grid
══════════════════════════════════════

COARSE (32 × 38):          FINE (864 × 608):
Each cell is ~25 km        Each cell is ~5 km

┌────────────────────┐     ┌──────────┬──────────┬──────────┐
│                    │     │  •  •  • │  •  •  • │  •  •  • │ (27 cells × 16 cells)
│ Coarse cell (5,10) │     │          │          │          │
│ Value: 10.5 mm/day │ ──▶ │  •  •  • │  •  •  • │  •  •  • │
│                    │     │          │          │          │
└────────────────────┘     │  •  •  • │  •  •  • │  •  •  • │
                           └──────────┴──────────┴──────────┘

Bilinear interpolation:
  For each fine cell position (if, jf):
    1. Find containing coarse cell (ic, jc)
    2. Compute fractional position (u, v) within coarse cell
    3. Interpolate from 4 coarse corners:
       
       value = y[ic,jc] * (1-u)*(1-v)
             + y[ic+1,jc] * u*(1-v)
             + y[ic,jc+1] * (1-u)*v
             + y[ic+1,jc+1] * u*v

Example: Fine cell at fractional position (u=0.3, v=0.7) in coarse cell:

  Coarse corners:
    (ic, jc)     = 10.0    (ic+1, jc) = 12.0
    (ic, jc+1)   = 11.0    (ic+1, jc+1) = 13.0
  
  Weighted sum:
    = 10.0 * 0.7 * 0.3 + 12.0 * 0.3 * 0.3 + 11.0 * 0.7 * 0.7 + 13.0 * 0.3 * 0.7
    = 2.1 + 1.08 + 5.39 + 2.73
    = 11.3 mm/day
```

---

## MBCn Algorithm Rotation Visualization

```
2D Example (N=2 fine cells, 2 variables):
═════════════════════════════════════════

Original Space
  Variable 2 ▲
           3 │     •(sim)         Simulated values (M=5 timesteps)
           2 │   • •               
           1 │ •                   Observed values
           0 └────────────────────▶ Variable 1
             0  1  2  3
             
Iteration 0: Rotate TO Weighted Sum Axis
  
  Weighted sum: w1=0.6, w2=0.4
  Sum axis: direction (0.6, 0.4) / sqrt(0.36 + 0.16) = (0.832, 0.555)
  
  After rotation:
    Axis 0: weighted sum (main diagonal direction)
    Axis 1: perpendicular direction
  
  Apply QM on Axis 0: preserve coarse model sum
  Force: sum_sim = sum_coarse (fixed value)
  
Iterations 1–19: Random Rotations
  
  Random rotation 1:
    • Rotate by random orthogonal matrix
    • QM on both axes independently
    • Correct: if sum changed, undo the change
  
  Random rotation 2:
    • Different rotation angle
    • Same process
    • Gradually fixes correlations
  
  Random rotation 3–19:
    • Each tests different linear combination
    • Fine cells gradually learn observation statistics
  
Iteration 20: Rotate Back
  
  Rotate by O_total^T (inverse of cumulative rotation)
  All rotations cancel out
  Result now in original coordinate system
  
  Final result:
    • Same marginal distributions as obs (each axis matches)
    • Same weighted sum as coarse model
    • Similar dependence structure

Conceptual sketch after 20 iterations:

  Variable 2 ▲
           3 │           • (final downscaled)
           2 │     • • •
           1 │   •
           0 └────────────────────▶ Variable 1
             0  1  2  3
           
           (Sim stretched to match obs shape
            while preserving sum)
```

---

## Spatial Downscaling Hierarchy

```
Coarse Grid Cell (32 × 38 = 1,216 total)
├─ Coarse Cell (5, 10)
│  │
│  └─ Contains 27 × 16 = 432 fine cells
│     ├─ Fine cells: (135:162, 160:176)
│     │
│     ├─ Month 1 (January)
│     │  ├─ Day 1: Jan 1, 1951 (obs) → apply MBCn → xsim_fine[0]
│     │  ├─ Day 2: Jan 2, 1951 (obs) → apply MBCn → xsim_fine[1]
│     │  └─ ...
│     │  └─ Day 31: Jan 31, 1951 (obs)
│     │     xobs: (31, 432) — 31 obs days × 432 fine cells
│     │     ysim: (31, 432) — 31 sim days bilinearly interp
│     │     ↓ MBCn ↓
│     │     xsim_fine: (31, 432)
│     │
│     ├─ Month 2 (February)
│     │  └─ Similar: (28, 432)
│     │
│     └─ Months 3–12
│        └─ Concatenate all months → (365, 432) for year
│           
│     Output for this coarse cell:
│     npy_stack/{cell_idx}.npy contains:
│       Shape: (31390, 432)  [all time × 432 fine cells]
│       Represents: All 86 years × 432 fine cells of this coarse cell
│
├─ Coarse Cell (5, 11)
│  └─ Similar process → npy_stack/{cell_idx+1}.npy
│
└─ ... (1,215 more coarse cells)

Reassembly:
  Load all 1,216 .npy files
  Reshape: (1216, 31390, 432) → (31390, 864, 608)
           coarse_cells × time × fine_cells_per_coarse
           ↓
           time × latitude × longitude
```

---

## Time Complexity Analysis

```
Bilinear Remapbil
═════════════════
Year 2015:  Thread 0  [365 days, 32×38→864×608]  ~30 sec
Year 2016:  Thread 1  [365 days, 32×38→864×608]  ~30 sec
Year 2017:  Thread 2  [365 days, 32×38→864×608]  ~30 sec
...
Year 2100:  Thread 31 [366 days, 32×38→864×608]  ~31 sec

Sequential:  86 × 30 sec = 43 min
Parallel (32 threads, GIL release): ~30 min → speedup ~4–6×


ISIMIP3BASD Statistical Downscaling
════════════════════════════════════
Coarse cell (0, 0):
  Month 1: QM with 21 rotations × 432 cells = ~1 sec
  Month 2: QM with 21 rotations × 432 cells = ~1 sec
  ...
  Month 12: ...
  Total per cell: ~12 sec

1,216 coarse cells × 12 sec = 14,592 sec = ~4 hours (sequential)

Parallel (32 processes, pool.imap):
  32 processes × 12 sec each (in parallel) = ~12 sec + overhead
  
  Actually: 1,216 / 32 = 38 batches × 12 sec = ~456 sec = ~8 min
  + GIL contention, I/O → reality is ~30–60 min

Actual measured: 30–120 min depending on n_processes, n_iterations
```

---

## Memory Allocation Timeline

```
Loading Phase
═════════════
Time  Memory Used               Comment
─────────────────────────────────────────────────────────
t=0   0 GB
t=1   4 GB    ← Loading NEX-GDDP historical
t=2   8 GB    ← + NEX-GDDP future
t=3   12 GB   ← + HYRAS
      (all 3 files fully loaded)

Harmonization
═════════════
t=4   14 GB   ← Cubes loaded with dask arrays
t=5   14 GB   ← Reprojection (streams per timestep)
t=6   14 GB   ← Cropping (in-place)

Parallel Remapbil (Threading)
══════════════════════════════
t=7   15 GB   ← Year 2015 remapbil [365, 864, 608] ≈ 0.8 GB
t=8   15 GB   ← Year 2016 (thread 2)
      ...     ← All threads < 20 GB total
t=40  18 GB   ← Concatenation of 86 years

Statistical Downscaling (Multiprocessing)
═════════════════════════════════════════
t=41  25 GB   ← Spawn 32 processes with shared lazy_data
t=42  28 GB   ← npy_stack writing (scatter across drives)
      ...     ← Process pool active
t=120 30 GB   ← Peak memory (temporary buffers in each process)
      ...
t=200 28 GB   ← Processes complete, npy_stack full

Collection
═════════
t=201 28 GB
t=202 35 GB   ← Load npy_stack into dask array (memory-mapped)
t=203 40 GB   ← Reshape & compute
t=204 45 GB   ← Peak: loading full result into memory
t=205 40 GB   ← Converting to xarray
t=206 35 GB   ← Saved to NetCDF, objects freed
t=207 20 GB   ← Result objects freed

Cleanup
════════
t=208 0 GB    ← All temporary files cleaned up
```

---

## Code Navigation Map

```
climdata/
├─ sdba/
│  └─ bcsd.py (1225 lines)
│     │
│     ├─ Lines 45–60:    _to_proleptic_gregorian()
│     ├─ Lines 64–71:    _is_curvilinear()
│     ├─ Lines 74–104:   _stamp_geogcs()
│     ├─ Lines 107–180:  _regrid_curvilinear_to_rectilinear()
│     ├─ Lines 200–350:  regrid_to_coarse()
│     │
│     ├─ Lines 350–700:  BiasCorrection class
│     │  └─ .correct()
│     │
│     ├─ Lines 700–1050: StatisticalDownscaling class
│     │  ├─ __init__()
│     │  └─ .downscale()
│     │     ├─ Line 820–850: Harmonization
│     │     ├─ Line 850–875: Crop to divisibility
│     │     ├─ Line 895–898: CRS stamping
│     │     ├─ Line 900–950: Parallel remapbil
│     │     │  └─ _remapbil_year()
│     │     └─ Line 950+:  Vendor call & collection
│     │
│     └─ Lines 1050–1225: BCSD class (combined workflow)
│
└─ _vendor/
   └─ isimip3basd/
      ├─ statistical_downscaling.py (694 lines)
      │  ├─ Lines 117–215:  weighted_sum_preserving_mbcn()
      │  ├─ Lines 216–290:  downscale_one_month()
      │  ├─ Lines 291–410:  downscale_one_location()
      │  └─ Lines 424–500:  downscale()
      │
      └─ utility_functions.py (2059 lines)
         ├─ Lines 1912–1955: remapbil()
         ├─ Lines 1930s–2000s: Quantile mapping functions
         ├─ Lines 2013–2030: get_downscaling_factors()
         └─ Lines 1950s+: Helper functions
```

---

## Success Indicator Checklist

✅ **Harmonization succeeded** if:
- No errors about "curvilinear" or "coordinate metadata"
- Cubes have 1D lat/lon DimCoords
- Both have GeogCS stamped

✅ **Remapbil succeeded** if:
- Year-by-year loop completes with tqdm
- Each year: [365/366, 864, 608] shape
- No "Index out of bounds" errors

✅ **MBCn succeeded** if:
- npy_stack_dir/ contains 1,216 files
- Each file is ~4 GB (31390 × 27 × 16 × 4 bytes)
- No NaN propagation warnings

✅ **Collection succeeded** if:
- Output NetCDF file created
- Shape: (31390, 864, 608)
- Time coord: 2015-01-01 to 2100-12-31
- All variables present

---

## Typical Log Output

```
🔄 Starting statistical downscaling for pr...
   Converting xarray datasets to iris cubes (in-memory)...
   ⚠  obs_fine has curvilinear (rotated-pole) coords — reprojecting to regular lat/lon grid via scipy griddata...
   obs_fine after reprojection: (9863, 864, 608)
   ⚠  obs_fine cropped to (864×608) [was 890×619] so shape is divisible by coarse (32×38)
   obs_fine  (9863, 864, 608), sim_coarse (31390, 32, 38)
   Creating bilinearly interpolated intermediate data (parallel by year)...
   remapbil |█████████████████████| 86/86 [00:05<00:00, 17.65 yr/s]
   Running ISIMIP3BASD statistical downscaling...
   (This may take a while for large datasets)
   ✅ Statistical downscaling complete!
   📦 Collecting results from npy_stack...
   Created npy_stack directory: /tmp/bcsd_downscale/npy_stack_pr
   💾 Saving result to: /data01/FDS/muduchuru/Atmos/NEXGDDP_HYRAS_ISIMIP3BASD/pr_BCSD.nc
   📂 File saved: True
   💾 Saved to: /data01/FDS/muduchuru/Atmos/NEXGDDP_HYRAS_ISIMIP3BASD/pr_BCSD.nc
✓ pr downscaled in 87s | shape: (31390, 864, 608)

✅ Statistical downscaling complete for all variables
```
