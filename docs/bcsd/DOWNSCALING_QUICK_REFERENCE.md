# Downscale Method: Quick Reference

## Method Signature

```python
result = sd.downscale(
    obs_fine=obs_raw,                   # xr.Dataset
    sim_coarse=bc_result[[var]],        # xr.Dataset
    output_path="output.nc",            # str or Path
)
```

**Returns**: `xr.Dataset` with downscaled variable at fine resolution

---

## What Gets Called (Call Stack)

```
StatisticalDownscaling.downscale()          [bcsd.py line 789]
├─ _to_proleptic_gregorian()                [normalise calendar]
├─ _is_curvilinear()                        [detect rotated-pole grid]
├─ _regrid_curvilinear_to_rectilinear()     [scipy.griddata reprojection]
├─ _stamp_geogcs()                          [attach CRS + normalise units]
├─ _remapbil_year() × 86                    [ThreadPoolExecutor parallel]
│  └─ uf.remapbil()                         [iris.analysis.Linear regrid]
│     ├─ coarse.regrid(fine, iris.analysis.Linear)
│     └─ coarse.regrid(fine, iris.analysis.Nearest) [fallback]
├─ iris.cube.CubeList().concatenate_cube()  [merge years]
└─ sd.downscale()                           [VENDOR: _vendor/isimip3basd]
   ├─ initializer(lazy_data, month_numbers)  [setup multiprocessing globals]
   ├─ downscale_one_location() × 1216        [Pool or sequential]
   │  ├─ get_fine_location_indices()
   │  └─ downscale_one_month() × 12
   │     └─ weighted_sum_preserving_mbcn()   [core MBCn algorithm]
   │        ├─ sample_invalid_values()       [NaN handling]
   │        ├─ randomize_censored_values()   [threshold handling]
   │        └─ map_quantiles_non_parametric_with_constant_extrapolation() × (n_iterations+2)
   ├─ Load results from npy_stack/
   └─ Reshape & return
```

---

## Key Transformations

### 1. xarray → iris
```
xr.Dataset[var]
  ↓ .to_iris()
iris.Cube (with dask arrays)
  ↓ _to_proleptic_gregorian()
iris.Cube (calendar updated)
```

### 2. Curvilinear → Rectilinear (if needed)
```
HYRAS rotated-pole:
  lat: AuxCoord (890, 619)  [2D, one value per grid cell]
  lon: AuxCoord (890, 619)

    ↓ scipy.interpolate.griddata × 9863 timesteps

Regular lat/lon:
  lat: DimCoord (890,)      [1D, values only]
  lon: DimCoord (619,)
```

### 3. Bilinear Remapbil (per year)
```
sim_coarse: [365, 32, 38]        [NEX-GDDP for 1 year, coarse]
obs_fine:   [1, 890, 619]        [HYRAS spatial template]

    ↓ iris.analysis.Linear.regrid()

sim_coarse_remapbil: [365, 890, 619]  [NEX-GDDP on HYRAS grid]
```

### 4. ISIMIP3BASD Statistical Downscaling
```
Training:
  obs_fine:             [9863, 890, 619]     (1951–1980, daily)
  sim_coarse:           [31390, 32, 38]      (2015–2100, daily)
  sim_coarse_remapbil:  [31390, 890, 619]    (derived above)

For each coarse grid cell (0,0)→(31,37), each month:
  x_obs_fine:           [N_obs, 432]         (all obs for month in cell)
  x_sim_coarse:         [N_sim,]             (scalar coarse model value)
  x_sim_coarse_remapbil:[N_sim, 432]         (bilinearly interpolated)
  
  ↓ weighted_sum_preserving_mbcn() with 20 random rotations
  
  x_sim_fine:           [N_sim, 432]         (downscaled per fine cell)

    ↓ Save: npy_stack/{coarse_cell_idx}.npy

Reshape & collect:
  All coarse cells → npy_stack/ (1216 files, each [31390, 432])
  ↓ dask.array.from_npy_stack()
  ↓ reshape([31390, 890, 619])
  
Result: sim_fine [31390, 890, 619]  (downscaled 2015–2100 daily)
```

---

## Input/Output Shapes

| Stage | obs_fine | sim_coarse | Output |
|-------|----------|-----------|--------|
| **Input** | (9863, 890, 619) | (31390, 32, 38) | — |
| After reprojection | (9863, 864, 608) | (31390, 32, 38) | — |
| After remapbil | (9863, 864, 608) | (31390, 32, 38) | (31390, 864, 608) |
| After ISIMIP3BASD | — | — | (31390, 864, 608) |

**Final output**: `(86 years × 365 days, 864 lat, 608 lon) = (31390, 864, 608)`

---

## Parallelism Strategy

### Bilinear Remapbil
```
Year 2015 → _remapbil_year(2015) → iris regrid [365 timesteps]
Year 2016 → _remapbil_year(2016) → iris regrid [365 timesteps]
...
Year 2100 → _remapbil_year(2100) → iris regrid [366 timesteps]

ThreadPoolExecutor(max_workers=32)
├─ Thread 0: year 2015 (iris.analysis.Linear releases GIL)
├─ Thread 1: year 2016
├─ Thread 2: year 2017
├─ ...
└─ Thread 31: year 2046

Speedup: ~4–6× vs sequential (scipy.RegularGridInterpolator releases GIL)
```

### ISIMIP3BASD Downscaling
```
Coarse cell (0,0)  → downscale_one_location() [32 processes]
Coarse cell (0,1)  → downscale_one_location()
Coarse cell (0,2)  → downscale_one_location()
...
Coarse cell (31,37) → downscale_one_location()

1216 coarse cells total, processed in parallel with n_processes=32

Each process:
  For month in [1–12]:
    For each fine cell in coarse cell (27×16 = 432 cells):
      weighted_sum_preserving_mbcn(20 rotations)
      ↓
      [31390 timesteps] → [31390 timesteps]
```

---

## Critical Parameters

### From `StatisticalDownscaling` class
```python
sd = StatisticalDownscaling(
    variable='pr',              # Climate variable
    n_iterations=20,            # MBCn rotation iterations (default: 20)
    n_processes=32,             # Parallel processes for downscaling
    randomization_seed=42,      # For reproducible copula rotations
    
    # Variable-specific thresholds (defaults shown for 'pr'):
    lower_bound=0,              # Physical lower limit
    lower_threshold=0.1,        # Below: randomize before MBCn
    upper_bound=None,
    upper_threshold=None,
)
```

### Variable-Specific Bounds

**Precipitation (pr)**:
```python
lower_bound=0,           # Can't be negative
lower_threshold=0.1,     # Drizzle events (< 0.1 mm/day) randomised
```

**Humidity (hurs)**:
```python
lower_bound=0,           # Can't be < 0%
lower_threshold=0.01,
upper_bound=100,         # Can't be > 100%
upper_threshold=99.99,   # Super-saturation events clamped
```

---

## Troubleshooting: Common Errors

### ❌ ValueError: Cube must contain single 1D x coordinate
**Cause**: HYRAS has 2D AuxCoords (rotated-pole), iris.analysis.Linear needs 1D
**Fix**: Auto-handled by `_is_curvilinear()` → `_regrid_curvilinear_to_rectilinear()`

### ❌ ValueError: coordinates... do not have matching coordinate metadata
**Cause**: xarray→iris produces `units='degrees_north'`, scipy.griddata produces `units='degrees'`
**Fix**: `_stamp_geogcs()` attaches identical CRS + normalises units

### ❌ ValueError: Unsupported units for coordinate system. Expected 'degrees' got Unit('m')
**Cause**: iris asserts `units=='degrees'` strictly once a CRS is attached
**Fix**: Part of `_stamp_geogcs()` — normalises before stamping CRS

### ❌ IndexError: Index is out of bounds for axis 0 with size 10958
**Cause**: `obs_fine_cube[_idx]` where `_idx` > 10958 (obs_fine only has 30 years, sim_coarse has 86)
**Fix**: Use `obs_fine_cube[0]` (remapbil only needs spatial template, not temporal data)

### ❌ Assertion Error: downscaling_factors not exact divisor
**Cause**: 890 % 32 ≠ 0 (not evenly divisible)
**Fix**: Crop to divisible shape before downscaling (864 % 32 = 0 ✓)

---

## Output File Structure

```
{OUTPUT_DIR}/
├── {SOURCE_ID}_{MEMBER_ID}_{SCENARIO}_BC.nc          [Bias-corrected coarse]
├── {SOURCE_ID}_{MEMBER_ID}_{SCENARIO}_pr_BCSD.nc     [Downscaled precip]
├── {SOURCE_ID}_{MEMBER_ID}_{SCENARIO}_hurs_BCSD.nc   [Downscaled humidity]
├── {SOURCE_ID}_{MEMBER_ID}_{SCENARIO}_tas_BCSD.nc    [Downscaled temp]
├── ... [other variables]
└── cache/
    ├── sim_hist_GFDL-ESM4_1951-2014.nc
    ├── sim_fut_GFDL-ESM4_ssp370_2015-2100.nc
    └── obs_raw_HYRAS_1951-1980.nc
```

Each `_BCSD.nc` file:
- **Dimensions**: `(31390, 864, 608)` — time, lat, lon
- **Coords**: 
  - `time`: 2015-01-01 to 2100-12-31 (daily)
  - `latitude`: regular 1D, ~890 values
  - `longitude`: regular 1D, ~608 values
- **Data variables**: single variable (pr, hurs, tas, etc.)
- **Attributes**:
  - `downscaling_method`: "ISIMIP3BASD modified MBCn"
  - `n_iterations`: 20
  - Original variable metadata (units, standard_name, etc.)

---

## References in Code

| Topic | File | Lines |
|-------|------|-------|
| High-level wrapper | `bcsd.py` | 789–1050 |
| Helper: proleptic_gregorian | `bcsd.py` | 45–60 |
| Helper: curvilinear detection | `bcsd.py` | 64–71 |
| Helper: CRS stamping | `bcsd.py` | 74–104 |
| Helper: reprojection | `bcsd.py` | 107–180 |
| Parallel remapbil | `bcsd.py` | 900–950 |
| Vendor: downscale entry | `_vendor/isimip3basd/statistical_downscaling.py` | 424–500 |
| Vendor: per-location loop | `_vendor/isimip3basd/statistical_downscaling.py` | 291–410 |
| Vendor: MBCn core | `_vendor/isimip3basd/statistical_downscaling.py` | 117–215 |
| Vendor: remapbil | `_vendor/isimip3basd/utility_functions.py` | 1912–1955 |
| Vendor: downscaling factors | `_vendor/isimip3basd/utility_functions.py` | 2013–2030 |
