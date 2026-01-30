# Bias Correction and Statistical Downscaling (BCSD)

This module provides tools for bias correction and statistical downscaling of climate model data using the ISIMIP3BASD methodology.

## Overview

The BCSD workflow has two main steps:

1. **Bias Correction (BC)**: Adjusts systematic biases in climate model outputs using quantile mapping
2. **Statistical Downscaling (SD)**: Increases spatial resolution using fine-scale observations

## Installation

### Step 1: Install climdata dependencies

The BCSD module requires additional dependencies:

```bash
pip install iris dask scipy
```

### Step 2: Download ISIMIP3BASD code

The original ISIMIP3BASD code is not included due to licensing. Download it from:
https://github.com/ISI-MIP/isimip3basd

Place these files in `climdata/_vendor/isimip3basd/`:
- `bias_adjustment.py`
- `statistical_downscaling.py`
- `utility_functions.py`

## Quick Start

```python
from climdata.sdba import BCSD

# Initialize BCSD for temperature
bcsd = BCSD(
    variable='tas',
    regridding_tool='xesmf',  # or 'cdo'
    regridding_method='conservative'
)

# Run complete workflow - obs_hist_coarse automatically derived!
result = bcsd.run(
    obs_fine=obs_025deg,             # Historical obs at fine resolution
    sim_hist_coarse=gcm_hist_1deg,  # Historical GCM
    sim_fut_coarse=gcm_fut_1deg,    # Future GCM
    output_path='./outputs/tas_bcsd.nc'
)
```

**Key Feature**: No need to manually create coarse observations! The module automatically regrids `obs_fine` to match the GCM grid using your choice of xESMF (Python) or CDO (command-line).

## Usage Examples

### Example 1: Temperature Bias Correction and Downscaling

```python
from climdata.sdba import BCSD

# Automatic configuration for temperature
bcsd = BCSD(
    variable='tas',
    regridding_tool='xesmf',  # Use Python-based xESMF
    regridding_method='conservative',  # Conservative for aggregation
    bias_correction_kwargs={
        'n_processes': 4  # Use 4 CPU cores
    }
)

# Only need to provide 3 datasets (fine obs auto-regridded to coarse!)
result = bcsd.run(
    obs_fine=obs_fine,
    sim_hist_coarse=sim_hist_coarse,
    sim_fut_coarse=sim_fut_coarse
)
```

### Example 2: Using CDO for Regridding

```python
from climdata.sdba import BCSD

# Use CDO instead of xESMF
bcsd = BCSD(
    variable='tas',
    regridding_tool='cdo',        # Use CDO command-line tool
    cdo_method='remapcon',        # Conservative remapping
    cdo_env='cdo_stable',         # Your conda env with CDO
    weights_dir='./weights'       # Save weights for reuse
)

result = bcsd.run(
    obs_fine=obs_fine,
    sim_hist_coarse=sim_hist_coarse,
    sim_fut_4: Providing Your Own Coarse Observations

If you already have native coarse-resolution observations (not derived from fine obs):

```python
# Load your own coarse observations
obs_hist_coarse = xr.open_dataset('./data/obs_coarse_native.nc')

bcsd = BCSD(variable='tas')

result = bcsd.run(
    obs_fine=obs_fine,
    sim_hist_coarse=sim_hist_coarse,
    sim_fut_coarse=sim_fut_coarse,
    obs_hist_coarse=obs_hist_coarse  # Provide your own!
)
```

### Example 5oarse=sim_fut_coarse
)
```

### Example 3: Precipitation with Custom Settings

```python
from climdata.sdba import BCSD

bcsd = BCSD(
    variable='pr',
    bias_correction_kwargs={
        'distribution': 'gamma',
        'trend_preservation': 'mixed',
        'lower_bound': 0,
        'lower_threshold': 0.1,
        'n_quantiles': 100
    }
)

result = bcsd.run(...)
```

### Example 3: Bias Correction Only
6
```python
from climdata.sdba import BiasCorrection

bc = BiasCorrection(
    variable='tas',
    distribution='normal',
    trend_preservation='additive',
    detrend=True
)

corrected = bc.correct(
    obs_hist=obs,
    sim_hist=gcm_hist,
    sim_fut=gcm_fut
)
```

### Example 4: Downscaling Only

```python
from climdata.sdba import StatisticalDownscaling

sd = StatisticalDownscaling(
    variable='tas',
    n_iterations=20
)

downscaled = sd.downscale(
    obs_fine=obs_fine,
    sim_coarse=sim_coarse_already_corrected
)
```Automatic Coarse Observation Derivation

When `obs_hist_coarse` is not provided, the module automatically regrids `obs_fine` to match the GCM grid:

**xESMF (Python):**
- Conservative remapping (area-weighted, mass-preserving)
- Bilinear interpolation (faster, less accurate)
- Nearest neighbor
- Weights cached and reused for efficiency

**CDO (Command-line):**
- `remapcon`: Conservative remapping (recommended)
- `remapbil`: Bilinear interpolation
- `remapdis`: Distance-weighted averaging
- `remapnn`: Nearest neighbor

Both methods ensure the coarse observations properly match the GCM grid structure.

### 

## Variable-Specific Configurations

The module includes optimal default configurations for common variables:

| Variable | Distribution | Trend Preservation | Bounds | Detrend |
|----------|--------------|-------------------|--------|---------|
| `tas` (Temperature) | Normal | Additive | None | Yes |
| `pr` (Precipitation) | Gamma | Mixed | [0, ∞) | No |
| `rsds` (Solar radiation) | Beta | Bounded | [0, max] | No |
| `hurs` (Humidity) | Beta | Bounded | [0, 100] | No |
| `sfcWind` (Wind) | Weibull | Mixed | [0, ∞) | No |

These defaults are based on Lange (2019) and can be overridden.

## Method Details

### Bias Correction

Uses trend-preserving quantile mapping:

1. **Build CDFs** for each day-of-year (±15 day window)
2. **Map quantiles** from model to observations
3. **Preserve trends** using additive/multiplicative/mixed methods
4. **Optional detrending** for temperature variables

Trend preservation methods:
- **Additive**: `corrected = obs_hist_quantile + (sim_fut - sim_hist)`
- **Multiplicative**: `corrected = obs_hist_quantile × (sim_fut / sim_hist)`
- **Mixed**: Combination based on threshold exceedance
- **Bounded**: Special handling for variables with physical bounds

### Statistical Downscaling

Uses modified MBCn (Multivariate Bias Correction) algorithm:

1. **Bilinear interpolation** to fine resolution (initial guess)
2. **Iterative quantile mapping** (default 20 iterations)
3. **Preserve spatial sums** within each coarse grid cell
4. **Match fine-scale patterns** from observations

## Complete Workflow with climdata

```python
from omegaconf import DictConfig
from climdata.datasets.CMIP_W5E5 import CMIPW5E5
from climdata.datasets.W5E5 import W5E5
from climdata.sdba import BCSD

# 1. Download CMIP6 historical (coarse)
cfg_hist = DictConfig({
    'variables': ['tas'],
    'time_range': {'start_date': '1980-01-01', 'end_date': '2014-12-31'},
    'experiment_id': 'historical',
    'source_id': 'gfdl-esm4',
    'data_dir': './data'
})
cmip_hist = CMIPW5E5(cfg_hist)
cmip_hist.fetch()
cmip_hist.load()

# 2. Download CMIP6 future (coarse)
cfg_fut = DictConfig({
    'variables': ['tas'],
    'time_range': {'start_date': '2015-01-01', 'end_date': '2050-12-31'},
    'experiment_id': 'ssp585',
    'source_id': 'gfdl-esm4',
    'data_dir': './data'
})
cmip_fut = CMIPW5E5(cfg_fut)
cmip_fut.fetch()
cmip_futfine=w5e5.ds,                # Fine obs (auto-regridded to coarse)
    sim_hist_coarse=cmip_hist.ds,
    sim_fut_coarse=cmip_fut.ds,
    output_path='./outputs/tas_ssp585_bcsd.nc'
)
```

## Regridding Options

### xESMF (Python-based)

**ProCache regridding weights**: Use `weights_dir` parameter
3. **Use conservative regridding**: Best for obs aggregation (mass-preserving)
4. **Process by region**: Extract smaller regions before BCSD
5. **Save intermediate results**: Use `save_intermediate=True`
6. **Chunk datasets**: Ensure data is well-chunked in time dimension
7. **Choose regridding tool wisely**:
   - xESMF: Better for repeated operations (cached weights)
   - CDO: Better for one-time large operations
- Good integration with xarray

**Cons:**
- Requires xESMF installation: `pip install xesmf`
- May need additional dependencies (ESMF library)

**Usage:**
```python
bcsd = BCSD(
    variable='tas',
    regridding_tool='xesmf',
    regridding_method='conservative',  # or 'bilinear', 'nearest'
    weights_dir='./weights'  # Cache weights here
)
```

### CDO (Command-line)

**Pros:**
- Very mature and reliable
- Excellent for large datasets
- Multiple specialized methods

**Cons:**
- Requires CDO installation
- Subprocess overhead
- May be slower for many small operations

**Usage:**
```python
bcsd = BCSD(
    variable='tas',
    regridding_tool='cdo',
    cdo_method='remapcon',  # or 'remapbil', 'remapdis', 'remapnn'
    cdo_env='cdo_stable'    # Conda env with CDO
)
```

**Install CDO:**
```bash
conda create -n cdo_stable -c conda-forge cdo   'time_range': {'start_date': '1980-01-01', 'end_date': '2014-12-31'},
    'data_dir': './data'
})
w5e5 = W5E5(cfg_obs)
w5e5.fetch()
w5e5.load()

# 4. Create coarse observations (aggregate fine obs)
obs_coarse = w5e5.ds.coarsen(lat=2, lon=2, boundary='trim').mean()

# 5. Run BCSD
bcsd = BCSD(variable='tas')
result = bcsd.run(
    obs_hist_coarse=obs_coarse,
    obs_fine=w5e5.ds,
    sim_hist_coarse=cmip_hist.ds,
    sim_fut_coarse=cmip_fut.ds,
    output_path='./outputs/tas_ssp585_bcsd.nc'
)
```

## Performance Tips

1. **Use parallel processing**: Set `n_processes=4` or higher
2. **Process by region**: Extract smaller regions before BCSD
3. **Save intermediate results**: Use `save_intermediate=True`
4. **Chunk datasets**: Ensure data is well-chunked in time dimension

## Reference

Lange, S. (2019). Trend-preserving bias adjustment and statistical downscaling with ISIMIP3BASD (v1.0). *Geoscientific Model Development*, 12(7), 3055-3070. https://doi.org/10.5194/gmd-12-3055-2019

## License

The ISIMIP3BASD code is licensed under GNU Affero General Public License v3.0.
See `climdata/_vendor/isimip3basd/README.md` for details.
