# Statistical Downscaling Pipeline: Detailed Walkthrough

## Overview

The downscaling pipeline in `climdata.sdba.bcsd.StatisticalDownscaling` wraps the **ISIMIP3BASD** (Integrated Scenarios for Mitigation and Adaptation Policy) bias adjustment and statistical downscaling package.

**Goal**: Take coarse-resolution bias-corrected climate model data (e.g., NEX-GDDP at 0.25°) and downscale it to fine-resolution observations (e.g., HYRAS at 5 km) using a modified multivariate MBCn algorithm.

---

## Pipeline Stages

### Stage 1: Input Preparation & Grid Harmonization

**Location**: `bcsd.py` lines ~800–900

#### Step 1.1: Data Format Conversion (xarray → iris)

```python
# Convert xarray variables to iris cubes
_of_cube = _to_proleptic_gregorian(obs_fine[self.variable].to_iris())
_sc_cube = _to_proleptic_gregorian(sim_coarse[self.variable].to_iris())
```

**Why iris?**
- ISIMIP3BASD is built on iris (UK MetOffice climate data library)
- iris handles complex coordinate systems (rotated poles, irregular grids)
- ISIMIP3BASD's `remapbil()` function uses `iris.analysis.Linear` for spatial interpolation

**Calendar normalization**:
- Reinterprets time coordinate calendar to `proleptic_gregorian` (required by vendor code)
- Safe because modern climate datasets use compatible calendars after 1582

#### Step 1.2: Curvilinear Grid Detection

```python
if _is_curvilinear(obs_fine_cube):
    obs_fine_cube = _regrid_curvilinear_to_rectilinear(obs_fine_cube)
```

**Why check?**
- HYRAS uses a **rotated-pole** projection → 2D latitude/longitude AuxCoords (890 y × 619 x)
- `iris.analysis.Linear` (RectilinearRegridder) requires **1D DimCoords** (rectilinear grids)
- Curvilinear grids must be reprojected first

**Reprojection method** (`_regrid_curvilinear_to_rectilinear`):
```
Input:  HYRAS 890×619 on rotated-pole → 2D lat(y,x), lon(y,x) AuxCoords
        └─ Each of ~551K grid cells has a different lat/lon value

↓ scipy.interpolate.griddata (per timestep)

Output: Regular lat/lon grid 890×619 with:
        └─ 1D latitude DimCoord (890 values)
        └─ 1D longitude DimCoord (619 values)
```

**Algorithm**:
1. Extract 2D lat/lon from AuxCoords
2. Build regular target grid: `lat = linspace(min, max, 890)`, `lon = linspace(min, max, 619)`
3. For each timestep: scatter interpolate from 2D curvilinear → regular grid using `scipy.interpolate.griddata(..., method='linear')`
4. Create iris Cube with 1D DimCoords

#### Step 1.3: Spatial Divisibility Check & Crop

```python
nlat_coarse, nlon_coarse = len(lat_coarse), len(lon_coarse)  # NEX: 32×38
nlat_fine, nlon_fine = obs_fine_cube.shape[1], obs_fine_cube.shape[2]  # HYRAS: 890×619

nlat_crop = (nlat_fine // nlat_coarse) * nlat_coarse  # (890 // 32) * 32 = 864
nlon_crop = (nlon_fine // nlon_coarse) * nlon_coarse  # (619 // 38) * 38 = 608

if nlat_crop != nlat_fine or nlon_crop != nlon_fine:
    # Trim edges symmetrically
    lat_trim = (nlat_fine - nlat_crop) // 2
    lon_trim = (nlon_fine - nlon_crop) // 2
    obs_fine_cube = obs_fine_cube[:, lat_trim:lat_trim+nlat_crop, lon_trim:lon_trim+nlon_crop]
```

**Why?**
- ISIMIP3BASD's `get_downscaling_factors()` asserts `shape_fine % shape_coarse == 0`
- Downscaling factors = (864÷32, 608÷38) = (27, 16)
- Each coarse cell is represented by a 27×16 block of fine cells

#### Step 1.4: CRS/Units Harmonization

```python
_stamp_geogcs(obs_fine_cube)
_stamp_geogcs(sim_coarse_cube)
```

**What it does** (`_stamp_geogcs` function):
```python
def _stamp_geogcs(cube):
    import iris.coord_systems as ics
    import cf_units
    geog_cs = ics.GeogCS(6371229.0)
    deg = cf_units.Unit('degrees')
    for coord in cube.coords(dim_coords=True):
        axis = iris.util.guess_coord_axis(coord)
        if axis in ('X', 'Y'):
            coord.coord_system = geog_cs      # Attach CRS
            coord.units = deg                 # Normalize to 'degrees'
```

**Why two fixes?**

| Problem | Cause | Solution |
|---------|-------|----------|
| **Metadata string mismatch** | xarray→iris: `units='degrees_north'`; scipy griddata: `units='degrees'` → iris compares metadata strings when no CRS → ValueError | Stamp **identical GeogCS** on both cubes |
| **Unit strictness with CRS** | Once a coord has GeogCS, iris asserts `units=='degrees'` exactly (not `'degrees_north'`) | Normalize to `'degrees'` explicitly |

---

### Stage 2: Bilinear Interpolation (Coarse → Fine Grid)

**Location**: `bcsd.py` lines ~900–950

#### Overview: Why remapbil?

We need to go from **coarse grid** (NEX 32×38) → **fine grid** (HYRAS 890×619).

Standard xarray `.interp()` cannot handle rotated-pole grids (HYRAS). Solution: use vendor's `remapbil()`.

#### Step 2.1: Parallel Year-by-Year Processing

```python
from concurrent.futures import ThreadPoolExecutor, as_completed

def _remapbil_year(yr):
    _idx = np.where(_years == yr)[0]  # All timesteps in year yr
    return yr, uf.remapbil(sim_coarse_cube[_idx], obs_fine_cube[0])

_unique_yrs = np.unique(_years)  # [2015, 2016, ..., 2100] for 86 years
_results_map = {}

with ThreadPoolExecutor(max_workers=32) as _pool:
    _futures = {_pool.submit(_remapbil_year, yr): yr for yr in _unique_yrs}
    _iter = as_completed(_futures)
    if tqdm is not None:
        _iter = tqdm(_iter, total=len(_unique_yrs), desc="remapbil", unit="yr")
    for _fut in _iter:
        _yr, _slice = _fut.result()
        _results_map[_yr] = _slice

_remap_slices = [_results_map[yr] for yr in _unique_yrs]
sim_coarse_remapbil_cube = iris.cube.CubeList(_remap_slices).concatenate_cube()
```

**Why parallel?**
- Each year's interpolation is independent
- `iris.analysis.Linear` calls `scipy.RegularGridInterpolator` which **releases the GIL**
- ThreadPoolExecutor effective without multiprocessing overhead (~4–6× speedup for 86 years)

**Why `obs_fine_cube[0]`?**
- `obs_fine_cube[0]` = single spatial grid (1 timestep, 890×619)
- `remapbil()` only uses the **target cube's coordinate system**, not its data or time axis
- `obs_fine_cube[_idx]` would fail: obs_fine only has ~10,958 timesteps (1951–1980), but `_idx` from sim_coarse has values > 31,000 (2015–2100)

#### Step 2.2: The remapbil Function

**Location**: `_vendor/isimip3basd/utility_functions.py` line 1912

```python
def remapbil(coarse, fine):
    """
    Bilinearly interpolates coarse iris cube to grid of fine iris cube.
    Falls back to nearest-neighbour where bilinear interpolation fails.
    """
    # Mark longitude circular for proper wraparound
    coarse.coord('longitude').circular = True
    
    # Bilinear interpolation
    bil = coarse.regrid(fine, iris.analysis.Linear(extrapolation_mode='mask'))
    
    # Where bilinear failed, use nearest-neighbour
    if np.ma.is_masked(bil.data):
        nn = coarse.regrid(fine, iris.analysis.Nearest(extrapolation_mode='extrapolate'))
        mask = np.logical_and(bil.data.mask, np.logical_not(nn.data.mask))
        bil.data.data[mask] = nn.data.data[mask]
        bil.data.mask = nn.data.mask
    
    return bil
```

**What happens**:
```
Input:  coarse cube [365, 32, 38] (1 year of NEX-GDDP)
        fine cube   [  1, 890, 619] (spatial template from HYRAS)

↓ iris.analysis.Linear.regrid()

Output: [365, 890, 619]
        └─ Each coarse grid cell linearly interpolated to 27×16 fine cells
        └─ Bilinear everywhere possible
        └─ Nearest-neighbour fallback at masked edges
```

---

### Stage 3: Statistical Downscaling (ISIMIP3BASD Core Algorithm)

**Location**: `_vendor/isimip3basd/statistical_downscaling.py` line 424

#### Overview: Modified MBCn Algorithm

**MBCn** = Multivariate Bias Correction in N dimensions

"Modified" = adds spatial downscaling via weighted sum preservation

#### Step 3.1: Vendor downscale() Initialization

```python
sd.downscale(
    obs_fine=obs_fine_cube,           # [9863, 890, 619] (1951–1980, HYRAS fine)
    sim_coarse=sim_coarse_cube,       # [31390, 32, 38] (2015–2100, NEX coarse)
    sim_coarse_remapbil=sim_coarse_remapbil_cube,  # [31390, 890, 619] (NEX on HYRAS grid)
    sim_fine_path=output_path,        # Where to save results
    n_processes=32,
    n_iterations=20,
    randomization_seed=42,
    lower_bound=..., upper_threshold=...  # Variable-specific bounds
)
```

**Key inputs**:

| Variable | Role | Shape |
|----------|------|-------|
| `obs_fine` | Training data: fine observations at fine resolution | (9863, 890, 619) |
| `sim_coarse` | Input to correct: bias-corrected model at coarse resolution | (31390, 32, 38) |
| `sim_coarse_remapbil` | Intermediate: sim_coarse bilinearly interpolated to fine grid | (31390, 890, 619) |

#### Step 3.2: Preprocessing in Vendor

```python
# Extract dask lazy arrays
lazy_data = {
    'obs_fine': obs_fine_cube.core_data(),
    'sim_coarse': sim_coarse_cube.core_data(),
    'sim_coarse_remapbil': sim_coarse_remapbil_cube.core_data(),
}

# Extract month numbers for each cube
month_numbers = {
    key: convert_datetimes(cube.time, 'month_number') for key, cube in cubes.items()
}

# Compute downscaling factors
downscaling_factors = get_downscaling_factors(
    obs_fine.shape[1:],      # (890, 619)
    sim_coarse.shape[1:]     # (32, 38)
)  # → (27, 16)

# Compute area weights for cosine latitude (HYRAS is on regular lat/lon)
sum_weights = iris.analysis.cartography.cosine_latitude_weights(obs_fine[0])
# → (890, 619) array: lower values at poles, ~1 at equator

# Generate random rotation matrices for copula (to preserve inter-variable dependence)
rotation_matrices = [uf.generateCREmatrix(np.prod(downscaling_factors))  # 27*16=432
                     for i in range(n_iterations)]  # 20 random 432×432 orthogonal matrices
```

#### Step 3.3: Core Algorithm: `downscale_one_location()`

**Called for each coarse grid cell**: (32 × 38 = 1,216 calls)

```python
for i_loc_coarse in np.ndindex(sim_coarse.shape[1:]):  # (0,0), (0,1), ..., (31,37)
    downscale_one_location(
        i_loc_coarse=i_loc_coarse,
        months=[1,2,...,12],
        downscaling_factors=(27, 16),
        ...
    )
```

**What downscale_one_location does**:

1. **Get fine grid cells** corresponding to this coarse cell:
   ```python
   i_locs_fine = get_fine_location_indices(i_loc_coarse, (27, 16))
   # e.g., for coarse (5,10) → fine cells [(135,160), (135,161), ..., (161,175)]
   ```

2. **For each calendar month**:
   ```python
   for month in months:
       # Extract timeseries for this month across all years
       x_obs_fine = obs_fine[month_numbers==month, i_locs_fine]  # (N_obs, 432)
       x_sim_coarse = sim_coarse[month_numbers==month, i_loc_coarse]  # (N_sim,)
       x_sim_coarse_remapbil = sim_coarse_remapbil[month_numbers==month, i_locs_fine]  # (N_sim, 432)
   ```

3. **Apply modified MBCn** (see next section):
   ```python
   x_sim_fine = weighted_sum_preserving_mbcn(
       x_obs_fine,              # (N_obs, 432) — observation training timeseries
       x_sim_coarse,            # (N_sim,) — scalar coarse model data
       x_sim_coarse_remapbil,   # (N_sim, 432) — bilinearly interpolated model
       sum_weights=area_weights,  # (432,) — cosine latitude weights for this coarse cell
       rotation_matrices=rotation_matrices  # 20 random 432×432 orthogonal matrices
   )  # → (N_sim, 432) downscaled timeseries
   ```

4. **Save result**:
   ```python
   np.save(f"{npy_stack_dir}/loc_1d.npy", x_sim_fine)
   ```

#### Step 3.4: The weighted_sum_preserving_mbcn() Core

**Location**: `_vendor/isimip3basd/statistical_downscaling.py` line 117

**Purpose**: Transform coarse model output to fine-resolution using observations

**Algorithm** (pseudo-code):
```
Initialize: o_total = identity matrix (432×432)

For i in range(n_iterations + 2):
    1. Compute rotation matrix o
       - Iteration 0: Rotate TO sum axis (weighted sum direction)
       - Iteration 1–19: Rotate BY random orthogonal matrix
       - Iteration 20: Rotate BACK to original axes

    2. Rotate all data
       x_obs_fine   ← x_obs_fine @ o
       x_sim        ← x_sim @ o
       sum_weights  ← sum_weights @ o

    3. Apply quantile mapping (non-parametric)
       if i==0:
           # Special handling: preserve coarse model value on sum axis
           x_sim[:, 0] ← x_sim_coarse (keep coarse value)
           # Map obs to match coarse model distribution
           x_obs_fine[:, 0] ← quantile_map(x_obs_fine[:, 0], x_sim_coarse)
       else:
           # Univariate quantile mapping for each fine cell
           for j in range(432):
               x_sim[:, j] ← quantile_map(x_sim[:, j], x_obs_fine[:, j])
           # Preserve weighted sum
           diff = x_sim - x_sim_prev
           x_sim ← x_sim - outer(diff @ sum_weights, sum_weights)

Return: x_sim (rotated back implicitly by next iteration's o_total^T)
```

**Key insight**: By rotating to the **sum axis first**, we preserve the coarse model's sum (total precipitation, energy balance, etc.) while allowing fine cells to adjust to match observation statistics.

---

### Stage 4: Result Collection & Export

**Location**: `bcsd.py` lines ~950–1050

#### Step 4.1: Load from npy_stack

The vendor code saves individual coarse grid cells' results to `.npy` files:

```
npy_stack_dir/
├── 0.npy          # Result for coarse cell (0,0), shape (31390, 432)
├── 1.npy          # Result for coarse cell (0,1), shape (31390, 432)
├── ...
└── 1215.npy       # Result for coarse cell (31,37), shape (31390, 432)
```

Collection:
```python
from dask.array import from_npy_stack

d = da.from_npy_stack(npy_stack_path, mmap_mode=None)
# → (1216, 31390, 432) dask array (using memory-mapping to avoid loading all at once)

d_reshaped = d.reshape(sim_fine_cube.shape)
# → (31390, 890, 619) — back to spatial shape [time, lat, lon]
```

#### Step 4.2: Assemble iris Cube & Save

```python
sim_fine_cube = obs_fine_cube.copy()  # Use HYRAS metadata as template
sim_fine_cube.data = d_reshaped.compute()  # Load into memory

# Save to netCDF
iris.save(sim_fine_cube, output_path,
          saver=iris.fileformats.netcdf.save,
          unlimited_dimensions=['time'],
          zlib=True, complevel=1)
```

#### Step 4.3: Convert Back to xarray

```python
result_cube = iris.load_cube(output_path)
result = xr.DataArray.from_iris(result_cube).to_dataset(name=self.variable)
result.attrs['downscaling_method'] = 'ISIMIP3BASD modified MBCn'
result.attrs['n_iterations'] = 20
result.to_netcdf(output_path)
```

**Output shape**: `(31390, 890, 619)` — 86 years of daily data at HYRAS resolution

---

## Data Flow Diagram

```
┌─────────────────────────────────────────────────────────────────────┐
│ INPUT: obs_fine (HYRAS)           sim_coarse (NEX-GDDP BC)         │
│   Shape: (9863, 890, 619)           Shape: (31390, 32, 38)          │
│   Period: 1951–1980 (30 yrs)        Period: 2015–2100 (86 yrs)      │
│   Grid: 5 km rotated-pole           Grid: 0.25° regular lat/lon     │
└─────────────────────────────────────────────────────────────────────┘
                                 ↓
        ┌────────────────────────────────────────────┐
        │ STAGE 1: Harmonization                     │
        ├────────────────────────────────────────────┤
        │ • xarray → iris conversion                 │
        │ • Calendar: → proleptic_gregorian         │
        │ • Curvilinear check: HYRAS rotated-pole   │
        │ • Reproject to regular lat/lon (scipy)    │
        │ • Crop to divisible shape (864×608)       │
        │ • Stamp CRS on both cubes                 │
        └────────────────────────────────────────────┘
                                 ↓
        obs_fine_cube: (9863, 864, 608)
        sim_coarse_cube: (31390, 32, 38)
                                 ↓
        ┌────────────────────────────────────────────┐
        │ STAGE 2: Parallel Bilinear Remapbil        │
        ├────────────────────────────────────────────┤
        │ • Year-by-year (ThreadPoolExecutor)        │
        │ • iris.analysis.Linear regrid              │
        │ • Fallback to Nearest for masked edges     │
        │ • 32×38 → 864×608 (27×16 factor)          │
        └────────────────────────────────────────────┘
                                 ↓
        sim_coarse_remapbil_cube: (31390, 864, 608)
                                 ↓
        ┌────────────────────────────────────────────┐
        │ STAGE 3: ISIMIP3BASD Statistical DS        │
        ├────────────────────────────────────────────┤
        │ • Per coarse cell (1,216 iterations)       │
        │ • Per calendar month (12 sub-iterations)   │
        │ • Modified MBCn with 20 rotation matrices  │
        │ • Non-parametric quantile mapping          │
        │ • Weighted sum preservation                │
        │ • Results → npy_stack/                     │
        └────────────────────────────────────────────┘
                                 ↓
        ┌────────────────────────────────────────────┐
        │ STAGE 4: Result Collection & Export        │
        ├────────────────────────────────────────────┤
        │ • Load from npy_stack (1,216 files)        │
        │ • Reshape (1216, 31390, 432) → (31390, 864, 608) │
        │ • Assemble iris cube                       │
        │ • Save to NetCDF                           │
        │ • Convert back to xarray                   │
        └────────────────────────────────────────────┘
                                 ↓
        ┌─────────────────────────────────────────────────────────────┐
        │ OUTPUT: sim_fine_downscaled                                 │
        │   Shape: (31390, 864, 608)                                  │
        │   Period: 2015–2100 (86 yrs daily)                         │
        │   Grid: ~5 km regular lat/lon                              │
        │   File: {SOURCE_ID}_{MEMBER_ID}_{SCENARIO}_{var}_BCSD.nc  │
        └─────────────────────────────────────────────────────────────┘
```

---

## Key Design Decisions & Rationale

### 1. Why Parallel Year-by-Year Bilinear Remapbil?

| Alternative | Problem |
|---|---|
| Sequential (all years at once) | Allocates (31390, 32, 38) → (31390, 864, 608) ≈ 8 GB peak memory |
| Parallel year-by-year (365, 32, 38) → (365, 864, 608) ≈ 20 MB per year | Chunked allocation, ~160× smaller |
| Multiprocessing | Pickling/unpickling iris cubes adds overhead |
| **Threading** ✓ | scipy.RegularGridInterpolator releases GIL, effective speedup |

### 2. Why Not Use xarray.interp()?

xarray `.interp()` assumes **Cartesian coordinates**. HYRAS after reprojection has:
- Irregular spatial spacing (scipy griddata preserves cell count but changes spacing)
- iris-compatible coordinate metadata required by RectilinearRegridder

Answer: Use vendor's `remapbil()` which handles CRS and irregular grids.

### 3. Why Crop Before Downscaling?

```python
if nlat_crop != nlat_fine or nlon_crop != nlon_fine:
    obs_fine_cube = obs_fine_cube[:, lat_trim:lat_trim+nlat_crop, ...]
```

ISIMIP3BASD assertion:
```python
downscaling_factors = shape_fine // shape_coarse
assert np.all(downscaling_factors * shape_coarse == shape_fine)  # Must be exact
```

- Without crop: 890÷32 = 27.8... (not divisible) → crash
- With crop: 864÷32 = 27 exactly ✓

### 4. Why obs_fine_cube[0] Instead of obs_fine_cube[_idx]?

```python
# WRONG:
return yr, uf.remapbil(sim_coarse_cube[_idx], obs_fine_cube[_idx])

# CORRECT:
return yr, uf.remapbil(sim_coarse_cube[_idx], obs_fine_cube[0])
```

- `obs_fine_cube.shape = (9863, 864, 608)` — only 30 years of data
- `sim_coarse_cube.shape = (31390, 32, 38)` — 86 years of data
- For year 2050, `_idx` might be `[12775, 12776, ..., 13139]` — all > 9863 → IndexError

**The fix**: `remapbil()` only uses **target cube's coordinate system**, not its time dimension. A single spatial snapshot is sufficient.

### 5. Why Modified MBCn with Rotation Matrices?

Standard univariate quantile mapping independently corrects each grid cell → loses inter-variable dependencies (e.g., temperature–humidity correlation).

MBCn with random rotations:
- Decorrelates variables via copula transforms
- 20 rotations → 20 chances to correct correlations at different angles
- Weighted sum preservation → preserves physically important totals

---

## Variable-Specific Parameters

From vendor documentation:

| Variable | Lower Bound | Lower Threshold | Upper Bound | Upper Threshold |
|---|---|---|---|---|
| `pr` (mm/day) | 0 | 0.1 | — | — |
| `hurs` (%) | 0 | 0.01 | 100 | 99.99 |
| `tas` (K) | — | — | — | — |
| `rsds` (W/m²) | 0 | 0.01 | — | — |
| `sfcWind` (m/s) | 0 | 0.01 | — | — |

**Applied in `downscale_one_month()`**:
```python
# Remove invalid values (NaN → interpolate from long-term mean)
x[key] = sample_invalid_values(d, seed, long_term_mean[key])[0]

# Replace censored values with random numbers
# e.g., pr < 0.1 mm/day → replace with random(0, 0.1)
x[key] = randomize_censored_values(x[key], 
    lower_bound=0, lower_threshold=0.1, upper_bound=None, upper_threshold=None)

# Run MBCn...

# De-randomize: replace random → bound
de_randomize_censored_values(x_sim_fine, ...)
```

---

## Performance Characteristics

**Single Variable Downscaling (1 NEX-GDDP model, 1 scenario)**:

| Step | Time (est.) | Memory |
|---|---|---|
| Load data | 2–5 min | 8–12 GB |
| Bias correction | 15–30 min | 12 GB |
| Bilinear remapbil (32 threads) | 3–5 min | 20 MB peak (per year) |
| ISIMIP3BASD statistical DS | 30–120 min* | 2–4 GB (npy_stack + processing) |
| Result collection | 2–5 min | 10 GB peak |
| **Total** | **60–180 min** | **~12 GB** |

*Depends on n_processes (32 → ~30 min) and n_iterations (20).

**Full Workflow** (all 6 variables in sequence):
- Bias correction: **once** (multivariate, 30 min)
- Per-variable downscaling: **6× in parallel** (~2 hours total)
- **Total: ~2.5–3 hours**

---

## Error Handling & Edge Cases

### Curvilinear Grid Error
```
ValueError: Cube must contain single 1D x coordinate
```
**Solution**: Auto-detect `_is_curvilinear()` → reproject via scipy.griddata

### Coordinate Metadata Mismatch
```
ValueError: The rectilinear grid coordinates... do not have matching coordinate metadata
```
**Solutions**:
1. Stamp identical GeogCS on both cubes
2. Normalize units to 'degrees'

### Index Out of Bounds
```
IndexError: Index is out of bounds for axis 0 with size 10958
```
**Solution**: Use `obs_fine_cube[0]` (spatial template) instead of time-indexed `obs_fine_cube[_idx]`

### Missing Values & Masked Arrays
- `dask.array.masked_invalid()` → mark NaN as masked
- Vendor `randomize_censored_values()` → impute censored values before MBCn
- `de_randomize_censored_values()` → clip back to bounds after MBCn

---

## References

- **ISIMIP3BASD Paper**: Lange, S. (2019). Trend-preserving bias adjustment and statistical downscaling with ISIMIP3BASD (v1.0). Geoscientific Model Development, 12(7), 3055–3070. https://doi.org/10.5194/gmd-12-3055-2019

- **ISIMIP3BASD Repo**: https://github.com/ISI-MIP/isimip3basd

- **iris Regridding**: https://scitools-iris.readthedocs.io/en/latest/userguide/regridding.html

- **scipy.interpolate.griddata**: https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html
