# ISIMIP3BASD Vendor Code Deep Dive

## Overview

The **ISIMIP3BASD** package (vendored in `climdata/_vendor/isimip3basd/`) provides the core statistical downscaling algorithm. It's authored by Stefan Lange at the Potsdam Institute for Climate Impact Research (PIK).

**Key Papers**:
- Lange, S. (2019). Trend-preserving bias adjustment and statistical downscaling with ISIMIP3BASD (v1.0). GMD, 12(7), 3055–3070. https://doi.org/10.5194/gmd-12-3055-2019

---

## Vendor Module Structure

```
climdata/_vendor/isimip3basd/
├── __init__.py                    [exports main functions]
├── bias_adjustment.py             [not used in this pipeline]
├── statistical_downscaling.py    [main: downscale() function]
├── utility_functions.py           [helpers: remapbil, quantile mapping, etc.]
├── README.md
├── LICENSE.txt
└── CHANGELOG.txt
```

---

## Key Vendor Functions

### 1. `statistical_downscaling.downscale()`

**Location**: `statistical_downscaling.py` line 424

**Signature**:
```python
def downscale(obs_fine, sim_coarse, sim_coarse_remapbil, sim_fine_path,
              n_processes=1, n_iterations=20, randomization_seed=None, **kwargs):
    """
    Applies modified MBCn algorithm for statistical downscaling 
    calendar month by calendar month and coarse grid cell by coarse grid cell.
    """
```

**Inputs**:

| Parameter | Type | Shape | Purpose |
|-----------|------|-------|---------|
| `obs_fine` | iris.Cube | (9863, 864, 608) | Observations at fine resolution (1951–1980 HYRAS) |
| `sim_coarse` | iris.Cube | (31390, 32, 38) | Bias-corrected simulations at coarse resolution (2015–2100 NEX) |
| `sim_coarse_remapbil` | iris.Cube | (31390, 864, 608) | sim_coarse bilinearly interpolated to fine grid |
| `sim_fine_path` | str | — | Where to save output (npy_stack created here) |
| `n_processes` | int | — | Number of parallel processes (1 = serial) |
| `n_iterations` | int | — | Number of MBCn rotation iterations (default 20) |
| `randomization_seed` | int | — | For reproducibility of random rotations |
| `**kwargs` | dict | — | Passed to `downscale_one_location()` (variable bounds, fill_value, etc.) |

**Algorithm**:

```python
# 1. Extract lazy dask arrays from cubes
lazy_data = {
    'obs_fine': obs_fine_cube.core_data(),           # dask array
    'sim_coarse': sim_coarse_cube.core_data(),
    'sim_coarse_remapbil': sim_coarse_remapbil_cube.core_data(),
}

# 2. Extract month numbers for each cube
month_numbers = {
    'obs_fine': convert_datetimes(obs_fine.coord('time'), 'month_number'),
    'sim_coarse': convert_datetimes(sim_coarse.coord('time'), 'month_number'),
    'sim_coarse_remapbil': convert_datetimes(sim_coarse_remapbil.coord('time'), 'month_number'),
}

# 3. Compute downscaling factors
downscaling_factors = get_downscaling_factors(
    obs_fine.shape[1:],      # (864, 608)
    sim_coarse.shape[1:]     # (32, 38)
)  # → (27, 16)

# 4. Generate random rotation matrices for MBCn copula
rotation_matrices = [uf.generateCREmatrix(np.prod(downscaling_factors))  # 432×432
                     for i in range(n_iterations)]  # 20 matrices

# 5. Compute area weights (cosine latitude)
sum_weights = iris.analysis.cartography.cosine_latitude_weights(obs_fine[0])
# → (864, 608) — cosine weighting for fair area representation

# 6. Create multiprocessing pool and send shared resources
if n_processes > 1:
    pool = mp.Pool(n_processes, 
                   initializer=initializer, 
                   initargs=(lazy_data, month_numbers))
    # ↑ All workers have access to lazy_data and month_numbers without copying
else:
    initializer(lazy_data, month_numbers)

# 7. Process each coarse grid cell in parallel
sdol = partial(downscale_one_location, 
               sim_fine_path=sim_fine_path,
               downscaling_factors=(27, 16),
               sum_weights=sum_weights,
               rotation_matrices=rotation_matrices,
               randomization_seed=randomization_seed,
               **kwargs)

i_locations_coarse = np.ndindex(sim_coarse.shape[1:])  # (0,0), (0,1), ..., (31,37)
foo = list(pool.imap(sdol, i_locations_coarse))  # 1216 calls total
pool.close()
pool.join()
pool.terminate()
```

**No explicit return** — results written to `npy_stack_dir(sim_fine_path)/`

---

### 2. `downscale_one_location()`

**Location**: `statistical_downscaling.py` line 291

**Signature**:
```python
def downscale_one_location(i_loc_coarse, sim_fine_path, downscaling_factors, sum_weights,
                          months=[1,2,3,4,5,6,7,8,9,10,11,12], fill_value=1.e20,
                          lower_bound=None, lower_threshold=None,
                          upper_bound=None, upper_threshold=None,
                          if_all_invalid_use=np.nan, resume_job=False, **kwargs):
    """
    Applies modified MBCn for statistical downscaling month by month 
    within one coarse grid cell.
    """
```

**Called once per coarse grid cell**: 1,216 times total (32 × 38)

**Example: Coarse cell (5, 10)**

```
┌─────────────────────────────────────────────────┐
│ Coarse grid cell (5, 10)                        │
│ Shape on coarse grid: (1, 1)                    │
│ Corresponding fine cells: (135:162, 160:176)    │
│ Number of fine cells: 27 × 16 = 432             │
└─────────────────────────────────────────────────┘

For each calendar month (12 iterations):
│
├─ Month 1 (January)
│  ├─ Extract observations: obs_fine[month==1, i_locs_fine]
│  │  Shape: (N_obs_jan, 432) — all Jan days in obs, all fine cells
│  │
│  ├─ Extract simulation (coarse): sim_coarse[month==1, (5,10)]
│  │  Shape: (N_sim_jan,) — scalar, just this coarse cell
│  │
│  ├─ Extract simulation (remapbil): sim_coarse_remapbil[month==1, i_locs_fine]
│  │  Shape: (N_sim_jan, 432) — coarse sim interpolated to fine cells
│  │
│  └─ Apply MBCn → downscaled_jan: (N_sim_jan, 432)
│
├─ Month 2 (February)
│  └─ [same as January]
│
└─ ...Month 12
   └─ [same as January]

Concatenate all months → final result for coarse cell (5,10)
Shape: (31390, 432) — all time, all 432 fine cells

Save: npy_stack_dir / (5*38 + 10).npy  # 1D index = 200
```

**Implementation**:

```python
def downscale_one_location(i_loc_coarse, ...):
    # Get fine location indices
    i_locs_fine = uf.get_fine_location_indices(i_loc_coarse, downscaling_factors)
    # e.g., (5, 10) + (27, 16) → [(135, 160), (135, 161), ..., (161, 175)]
    
    # Convert to 1D index for saving
    i_loc_1d = np.ravel_multi_index(i_loc_coarse, (32, 38))
    i_loc_path = lambda p, i: uf.npy_stack_dir(p) + '%i.npy' % i
    
    # Check if already done (resume_job)
    if resume_job and os.path.isfile(i_loc_path(sim_fine_path, i_loc_1d)):
        return None  # Already computed
    
    # For each calendar month
    x_sim_fine_all_months = []
    for month in months:
        # Get month mask for each dataset (from global_month_numbers set by initializer)
        month_mask_obs = global_month_numbers['obs_fine'] == month
        month_mask_sim = global_month_numbers['sim_coarse'] == month
        
        # Get data for this month (from global_lazy_data)
        x_obs_fine = global_lazy_data['obs_fine'][month_mask_obs][:, i_locs_fine]  # (N_obs_month, 432)
        x_sim_coarse = global_lazy_data['sim_coarse'][month_mask_sim, i_loc_coarse]  # (N_sim_month,)
        x_sim_coarse_remapbil = global_lazy_data['sim_coarse_remapbil'][month_mask_sim][:, i_locs_fine]  # (N_sim_month, 432)
        
        # Compute long-term means for this month
        long_term_mean = {
            'x_obs_fine': uf.mean_ignoring_nans(x_obs_fine),  # (432,)
            'x_sim_coarse': uf.mean_ignoring_nans(x_sim_coarse),  # scalar
            'x_sim_coarse_remapbil': uf.mean_ignoring_nans(x_sim_coarse_remapbil),  # (432,)
        }
        
        # Downscale this month
        data = {
            'x_obs_fine': x_obs_fine,
            'x_sim_coarse': x_sim_coarse,
            'sim_coarse_remapbil': x_sim_coarse_remapbil,
        }
        x_sim_fine_month = downscale_one_month(
            data, long_term_mean,
            lower_bound=lower_bound,
            lower_threshold=lower_threshold,
            upper_bound=upper_bound,
            upper_threshold=upper_threshold,
            randomization_seed=randomization_seed,
            **kwargs
        )  # → (N_sim_month, 432)
        
        x_sim_fine_all_months.append(x_sim_fine_month)
    
    # Concatenate all months in order
    x_sim_fine = np.concatenate(x_sim_fine_all_months, axis=0)  # (31390, 432)
    
    # Reshape: (31390, 432) → (31390, 27, 16) to match spatial structure
    x_sim_fine = x_sim_fine.reshape((-1,) + downscaling_factors)  # (31390, 27, 16)
    
    # Save to disk
    np.save(i_loc_path(sim_fine_path, i_loc_1d), x_sim_fine)
    
    print(i_loc_coarse, 'done')
    sys.stdout.flush()
```

**Outputs**:
- Creates `npy_stack_dir(sim_fine_path)/{i_loc_1d}.npy` files (one per coarse cell)
- Each file: `(31390, 27, 16)` array

---

### 3. `downscale_one_month()`

**Location**: `statistical_downscaling.py` line 216

**Core algorithm** (before MBCn):

```python
def downscale_one_month(data, long_term_mean, lower_bound=None, 
                       lower_threshold=None, upper_bound=None, 
                       upper_threshold=None, randomization_seed=None, **kwargs):
    """
    1. Replace invalid values (NaN)
    2. Replace censored values (< lower_threshold or > upper_threshold) with random
    3. Apply modified MBCn algorithm
    4. De-randomize: clip to bounds
    """
    
    # STEP 1: Replace NaN with long-term mean or random
    x = {}
    for key in ['x_obs_fine', 'x_sim_coarse', 'sim_coarse_remapbil']:
        x[key], n_sample_invalid = uf.sample_invalid_values(
            data[key], randomization_seed, long_term_mean[key]
        )
    
    # STEP 2: Randomize censored values
    for key in x:
        x[key] = uf.randomize_censored_values(
            x[key],
            lower_bound, lower_threshold, upper_bound, upper_threshold,
            False,  # not_randomize_censored
            False,  # not_randomize_invalid
            randomization_seed,
            10.,   # randomization_exponent (creates values close to bounds)
            10.
        )
    
    # STEP 3: CORE MBCn ALGORITHM (next section)
    x_sim_coarse_remapbil_copy = x['sim_coarse_remapbil'].copy()
    x_sim_fine = weighted_sum_preserving_mbcn(
        x['x_obs_fine'],
        x['x_sim_coarse'],
        x['sim_coarse_remapbil'],
        sum_weights=area_weights,  # passed via kwargs
        rotation_matrices=rotation_matrices,
        **kwargs
    )  # → (N_sim_month, 432)
    
    # STEP 4: De-randomize censored values (clip to bounds)
    uf.randomize_censored_values(
        x_sim_fine,
        lower_bound, lower_threshold, upper_bound, upper_threshold,
        True,  # not_randomize_censored ← DE-randomize
        True   # not_randomize_invalid ← Replace NaN back
    )
    
    # Verify no infs/nans remain
    uf.assert_no_infs_or_nans(x_sim_coarse_remapbil_copy, x_sim_fine)
    
    return x_sim_fine
```

---

### 4. `weighted_sum_preserving_mbcn()` — THE CORE

**Location**: `statistical_downscaling.py` line 117

**Most important function!**

```python
def weighted_sum_preserving_mbcn(
        x_obs, x_sim_coarse, x_sim,
        sum_weights, rotation_matrices=[], n_quantiles=50):
    """
    Applies modified MBCn algorithm.
    
    Parameters:
    -----------
    x_obs : (M, N) ndarray
        M observations, N fine grid cells (e.g., 2500 days × 432 cells)
    
    x_sim_coarse : (M,) array
        M simulated scalar values at coarse cell (e.g., 2500 days)
    
    x_sim : (M, N) ndarray
        M simulations bilinearly interpolated to N fine cells (e.g., 2500 × 432)
    
    sum_weights : (N,) array
        Grid cell area weights, normalized. Represents spatial importance
        of each fine cell within the coarse cell.
    
    rotation_matrices : list of (N, N) ndarrays
        Pre-computed random orthogonal matrices (one per MBCn iteration).
        Each rotates the (M, N) data to a new coordinate system.
    
    Returns:
    --------
    x_sim : (M, N) ndarray
        Downscaled simulations at fine resolution
    """
```

**Algorithm**:

```python
# Initialize total rotation matrix (starts as identity)
n_variables = sum_weights.size  # N = 432
o_total = np.diag(np.ones(n_variables))

# Quantile levels for non-parametric mapping
p = np.linspace(0., 1., n_quantiles+1)  # 0, 0.02, 0.04, ..., 1.0

# Normalize sum_weights to unit vector (ensures coarse model value preserved)
sum_weights = sum_weights / np.sqrt(np.sum(np.square(sum_weights)))

# Rescale x_sim_coarse for initial step (ensures sum is preserved)
x_sim_coarse = x_sim_coarse * np.sum(sum_weights)

# MAIN LOOP: n_iterations + 2 rotations total
n_loops = len(rotation_matrices) + 2  # 20 + 2 = 22
for i in range(n_loops):
    print(f"  MBCn iteration {i+1}/{n_loops}")
    
    # STEP 1: Choose rotation
    if i == 0:
        # First: rotate TO the weighted sum direction
        o = uf.generate_rotation_matrix_fixed_first_axis(sum_weights)
        # After this rotation, axis 0 is the weighted sum of all cells
    
    elif i == n_loops - 1:
        # Last: rotate BACK to original axes
        o = o_total.T  # Transpose of cumulative rotation
    
    else:
        # Middle iterations: random rotations
        o = rotation_matrices[i - 1]  # Use pre-computed random orthogonal matrix
    
    # Accumulate total rotation
    o_total = np.dot(o_total, o)
    
    # STEP 2: Rotate all data to new coordinate system
    x_sim = np.dot(x_sim, o)       # (M, 432) → rotated
    x_obs = np.dot(x_obs, o)
    sum_weights = np.dot(sum_weights, o)  # Weight vector also rotates
    
    # STEP 3: Apply quantile mapping on rotated axes
    if i == 0:
        # Iteration 0: special handling for the sum axis (axis 0)
        # ────────────────────────────────────────────────────
        # Preserve the coarse model's value on the sum axis
        x_sim[:, 0] = x_sim_coarse  # Keep original coarse values
        
        # Quantile-map observations to match coarse model distribution
        q_sim = uf.percentile1d(x_sim_coarse, p)        # (51,) quantiles of coarse
        q_obs = uf.percentile1d(x_obs[:, 0], p)         # (51,) quantiles of obs sum
        
        # Map obs values to sim distribution while preserving order
        x_obs[:, 0] = uf.map_quantiles_non_parametric_with_constant_extrapolation(
            x_obs[:, 0], q_obs, q_sim
        )
    
    else:
        # Iterations 1–21: standard univariate quantile mapping per cell
        # ──────────────────────────────────────────────────────────────
        x_sim_previous = x_sim.copy()
        
        for j in range(n_variables):
            # Get quantiles of simulated and observed distributions
            q_sim = uf.percentile1d(x_sim[:, j], p)     # (51,) quantiles of sim cell j
            q_obs = uf.percentile1d(x_obs[:, j], p)     # (51,) quantiles of obs cell j
            
            # Map simulated values to observed quantiles
            x_sim[:, j] = uf.map_quantiles_non_parametric_with_constant_extrapolation(
                x_sim[:, j], q_sim, q_obs
            )
            # Now x_sim[:, j] has the same marginal distribution as x_obs[:, j]
        
        # Preserve weighted sum of original variables
        # If the sum changed due to quantile mapping, adjust back
        if i < n_loops - 1:  # Don't preserve on last iteration
            diff = x_sim - x_sim_previous  # How much changed
            x_sim -= np.outer(np.dot(diff, sum_weights), sum_weights)
            # Subtract correction proportional to weights
    
    print(f"    x_sim shape: {x_sim.shape}, mean={x_sim.mean():.4f}")

return x_sim
```

**What this does**:

1. **Iteration 0**:
   - Rotates to the "sum axis" (weighted average direction)
   - Forces x_sim to preserve coarse model's sum value
   - Adjusts x_obs to match the coarse model's distribution

2. **Iterations 1–19**:
   - Random rotations (different combination of fine cells at each axis)
   - Univariate quantile mapping: each axis matches obs distribution
   - Correction step: if sum changed, undo it (weighted sum preservation)
   - Effectively: decorrelates variables and fixes them independently

3. **Iteration 20**:
   - Rotate back to original space
   - All rotations cancel out → result is in original coordinate system

**Result**: xsim now has:
- ✅ Same marginal distributions as observations (univariate)
- ✅ Same weighted sum as original coarse model (energy/mass conservation)
- ✅ Copula structure similar to observations (multivariate dependence)

---

### 5. `remapbil()` — Bilinear Regridding

**Location**: `utility_functions.py` line 1912

**Signature**:
```python
def remapbil(coarse, fine):
    """Bilinearly interpolates coarse iris cube to grid of fine iris cube."""
```

**Implementation**:

```python
def remapbil(coarse, fine):
    # Mark longitude as circular for proper wraparound
    try:
        coarse.coord('longitude').circular = True
    except iris.exceptions.CoordinateNotFoundError:
        warnings.warn('could not find longitude coordinate')
    
    # ─ Bilinear interpolation ─
    # iris.analysis.Linear uses scipy.RegularGridInterpolator
    # Requires both cubes to have the same CRS (or compatible coordinates)
    bil = coarse.regrid(fine, iris.analysis.Linear(extrapolation_mode='mask'))
    # Shape: coarse [365, 32, 38] → bil [365, 864, 608]
    
    # ─ Handle masked values ─
    # Where bilinear interpolation failed (e.g., at masked coarse cells),
    # use nearest-neighbour as fallback
    if np.ma.is_masked(bil.data):
        nn = coarse.regrid(fine, iris.analysis.Nearest(
            extrapolation_mode='extrapolate'
        ))
        # nn has same coarse values mapped to fine grid
        
        # Find where bil is masked but nn is not
        mask = np.logical_and(
            bil.data.mask, 
            np.logical_not(nn.data.mask)
        )
        # Fill those locations with nearest-neighbour values
        bil.data.data[mask] = nn.data.data[mask]
        bil.data.mask = nn.data.mask
    
    return bil
```

**Why two-step approach?**

| Method | Pros | Cons |
|--------|------|------|
| **Bilinear only** | Smooth, accurate | Fails at masked edges |
| **Nearest-only** | Always works | Blocky, pixelated |
| **Bilinear + NN fallback** ✓ | Smooth interior, reliable edges | Slight artifacts at boundaries |

---

## Utility Helper Functions

### `get_downscaling_factors(shape_fine, shape_coarse)`

```python
def get_downscaling_factors(shape_fine, shape_coarse):
    """
    Returns downscaling factors (fine / coarse for each dimension).
    Asserts that fine is exactly divisible by coarse.
    """
    factors = np.array(shape_fine) // np.array(shape_coarse)
    assert np.all((factors * np.array(shape_coarse)) == np.array(shape_fine)), \
        f"Shape {shape_fine} not divisible by {shape_coarse}"
    return factors
```

**Example**:
```python
get_downscaling_factors((864, 608), (32, 38))
# Checks: (864 // 32 * 32 == 864) ✓ and (608 // 38 * 38 == 608) ✓
# Returns: array([27, 16])
```

### `get_fine_location_indices(i_loc_coarse, downscaling_factors)`

```python
def get_fine_location_indices(i_loc_coarse, downscaling_factors):
    """Returns list of fine grid indices contained in this coarse cell."""
    i_loc_fine = []
    for dim, (i_coarse, factor) in enumerate(zip(i_loc_coarse, downscaling_factors)):
        i_fine_min = i_coarse * factor
        i_fine_max = i_fine_min + factor
        i_loc_fine.append(slice(i_fine_min, i_fine_max))
    return tuple(i_loc_fine)
```

**Example**:
```python
get_fine_location_indices((5, 10), (27, 16))
# Coarse (5, 10) contains fine cells:
#   latitude: 5*27 to 5*27+27 = 135:162
#   longitude: 10*16 to 10*16+16 = 160:176
# Returns: (slice(135, 162), slice(160, 176))
```

### `map_quantiles_non_parametric_with_constant_extrapolation(x, q_x, q_y)`

```python
def map_quantiles_non_parametric_with_constant_extrapolation(x, q_x, q_y):
    """
    Maps values from one distribution to another using quantile matching.
    
    Parameters:
    -----------
    x : array
        Values to map (e.g., simulated data)
    
    q_x : array
        Quantiles of x's original distribution (e.g., [min, p25, p50, p75, max])
    
    q_y : array
        Quantiles of target distribution (e.g., [min, p25, p50, p75, max])
    
    Returns:
    --------
    x_mapped : array
        x values mapped to y's distribution
    """
    # For each value in x, find where it falls in q_x, then replace with q_y value
    # Monotonic: preserves order (smallest stays smallest)
    # Constant extrapolation: x < min(q_x) → min(q_y)
    return np.interp(x, q_x, q_y, left=q_y[0], right=q_y[-1])
```

**Example**:
```
x = [−2, 0, 5, 10, 15]          # Simulated data
q_x = [−5, −2, 5, 12, 20]       # Its quantiles
q_y = [−10, 0, 0, 5, 10]        # Observed quantiles

result = map_quantiles_non_parametric_with_constant_extrapolation(x, q_x, q_y)
# −2 is at p25 in x → gets p25 from y → 0
# 0 is between p25 and p50 → linear interpolation → 0
# 5 is at p50 in x → gets p50 from y → 0
# 10 is between p50 and p75 → linear interpolation → 3
# 15 is between p75 and max → linear interpolation → 8

result ≈ [0, 0, 0, 3, 8]
```

### `generateCREmatrix(n)`

```python
def generateCREmatrix(n):
    """
    Returns random orthogonal n × n matrix from circular real ensemble.
    Reference: Mezzadri (2007) http://arxiv.org/abs/math-ph/0609050v2
    """
    # Compute QR decomposition of random normal matrix
    z = np.random.randn(n, n)
    q, r = scipy.linalg.qr(z)
    
    # Adjust signs so diagonal of r is positive
    d = np.diagonal(r)
    return q * (d / np.abs(d))
```

**Why this?**
- Uniformly distributed random orthogonal matrix
- Used for random rotations in MBCn iterations 1–19
- Each rotation tests a different "slice" through the data space
- 20 rotations → 20 chances to fix variable correlations

---

## Data Flow Through Vendor Code

```
Input Cubes
├─ obs_fine[t, lat, lon]        (9863, 864, 608)
├─ sim_coarse[t, lat, lon]      (31390, 32, 38)
└─ sim_coarse_remapbil[t, lat, lon]  (31390, 864, 608)

        ↓

downscale()
├─ Extract lazy arrays → lazy_data dict
├─ Extract month numbers → month_numbers dict
├─ Compute downscaling_factors = (27, 16)
├─ Compute sum_weights = cosine_latitude_weights [864, 608]
├─ Generate 20 random CRE matrices [432, 432 each]
│
└─ For each coarse cell (i, j) ∈ {0..31} × {0..37} [1216 total]:
   
   downscale_one_location(i, j)
   │
   ├─ Get fine cell indices: (slice(i*27, (i+1)*27), slice(j*16, (j+1)*16))
   │  [27×16 = 432 fine cells per coarse cell]
   │
   └─ For each month m ∈ {1..12}:
      
      downscale_one_month()
      │
      ├─ Extract obs_fine[month==m, fine_cells]        (N_m, 432)
      ├─ Extract sim_coarse[month==m, (i,j)]           (N_m,)
      ├─ Extract sim_coarse_remapbil[month==m, fine_cells]  (N_m, 432)
      │
      ├─ Sample invalid values → impute NaN
      ├─ Randomize censored values (< lower_threshold or > upper_threshold)
      │
      ├─ weighted_sum_preserving_mbcn()
      │  │
      │  ├─ Iteration 0: Rotate to sum axis
      │  │              Preserve coarse model sum
      │  │              Map obs to coarse distribution
      │  │
      │  ├─ Iterations 1–19: Random rotations
      │  │                   Univariate quantile mapping per axis
      │  │                   Restore weighted sum
      │  │
      │  └─ Iteration 20: Rotate back to original axes
      │     → x_sim_fine (N_m, 432)
      │
      ├─ De-randomize censored values → clip to bounds
      └─ Return: x_sim_fine (N_m, 432)
      
      All months concatenated → (31390, 432)

   Save to npy_stack_dir/{1D_index}.npy

        ↓

Wrapper collects results:
└─ Load all 1216 .npy files → dask array (1216, 31390, 432)
   └─ Reshape → (31390, 864, 608)
   └─ Assemble iris cube
   └─ Save to NetCDF
   └─ Convert back to xarray

        ↓

Output: sim_fine_downscaled (31390, 864, 608)
```

---

## Key Takeaways

1. **Modified MBCn is clever**: By rotating to the sum axis first, it preserves totals while allowing fine-scale adjustment.

2. **Parallel efficiency**: Workers store lazy dask arrays + month masks in shared memory → no data copying across processes.

3. **Month-by-month processing**: Allows different statistics for winter vs. summer; climatological consistency.

4. **Two-step regridding**: Bilinear for accuracy, nearest-neighbour fallback for reliability.

5. **Non-parametric quantile mapping**: Preserves order; handles extreme values via constant extrapolation.

---

## References

- **Lange (2019)**: https://doi.org/10.5194/gmd-12-3055-2019
- **ISIMIP3BASD GitHub**: https://github.com/ISI-MIP/isimip3basd
- **CRE Matrix**: Mezzadri, F. (2007). How to generate random matrices from the classical compact groups. arXiv preprint math-ph/0609050.
