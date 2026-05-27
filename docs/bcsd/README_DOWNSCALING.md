# Downscaling Documentation Index

Welcome! This directory contains comprehensive documentation of the statistical downscaling pipeline used in the NEXGDDP-HYRAS-BCSD workflow.

## Documents

### 1. **[DOWNSCALING_PIPELINE_EXPLAINED.md](DOWNSCALING_PIPELINE_EXPLAINED.md)** — Start Here
**Best for**: Understanding the complete workflow, step-by-step

- 📊 Full pipeline architecture with data flow diagrams
- 🔄 Four stages: Harmonization → Remapping → Statistical DS → Collection
- 🛠️ Helper functions (`_stamp_geogcs()`, `_regrid_curvilinear_to_rectilinear()`, etc.)
- ⚙️ Key design decisions and their rationale
- 🐛 Error handling guide
- 📈 Performance characteristics

**Read this first** if you're new to the downscaling process.

---

### 2. **[DOWNSCALING_QUICK_REFERENCE.md](DOWNSCALING_QUICK_REFERENCE.md)** — Cheat Sheet
**Best for**: Quick lookup, troubleshooting, code navigation

- 📋 Method signature and parameters
- 🔗 Complete call stack with line numbers
- 🎯 Input/output shapes at each stage
- ⚡ Parallelism strategy
- 🔧 Variable-specific parameter tables
- ❌ Common errors and fixes
- 📁 Output file structure

**Keep this open** while reading code or debugging.

---

### 3. **[VENDOR_CODE_DEEP_DIVE.md](VENDOR_CODE_DEEP_DIVE.md)** — Technical Details
**Best for**: Understanding ISIMIP3BASD internals, algorithm details

- 📦 Vendor module structure (`statistical_downscaling.py`, `utility_functions.py`)
- 🔍 Line-by-line breakdown of:
  - `downscale()` — main entry point
  - `downscale_one_location()` — per-coarse-cell loop
  - `downscale_one_month()` — monthly subdivision
  - `weighted_sum_preserving_mbcn()` — **the core algorithm**
  - `remapbil()` — bilinear regridding
- 🧮 Helper functions (quantile mapping, random rotation generation, etc.)
- 🔄 Complete data flow through vendor code with pseudocode

**Study this** if you need to understand MBCn algorithm internals or modify vendor behavior.

---

### 4. **[DOWNSCALING_MATHEMATICS.md](DOWNSCALING_MATHEMATICS.md)** — Equations & Theory
**Best for**: Mathematical formulation, algorithm validation, paper writing

- 📐 Formal problem statement with constraints
- 🔄 Bilinear interpolation equations
- 📊 Modified MBCn algorithm with full mathematical notation
- 📈 Quantile mapping theory
- 🎯 Weighted sum preservation proof sketch
- 📊 Complexity analysis (time, memory)
- 🔗 Notation table and academic references

**Reference this** for papers, talks, or rigorous understanding.

---

## Quick Navigation

### "How do I...?"

| Question | Document | Section |
|----------|----------|---------|
| Understand what downscaling does | PIPELINE | Overview |
| Debug an error | QUICK_REF | Troubleshooting |
| Find line numbers in code | QUICK_REF | References in Code |
| Understand MBCn algorithm | VENDOR_CODE | Section "weighted_sum_preserving_mbcn()" |
| Know input/output shapes | QUICK_REF | Input/Output Shapes |
| Learn mathematical details | MATHEMATICS | All sections |
| See performance estimates | PIPELINE | Performance Characteristics |
| Understand parallelisation | QUICK_REF | Parallelism Strategy |

### By Role

**Climate scientist** → [PIPELINE](DOWNSCALING_PIPELINE_EXPLAINED.md) + [QUICK_REF](DOWNSCALING_QUICK_REFERENCE.md)

**Software engineer** → [QUICK_REF](DOWNSCALING_QUICK_REFERENCE.md) + [VENDOR_CODE](VENDOR_CODE_DEEP_DIVE.md)

**Researcher / PhD** → [PIPELINE](DOWNSCALING_PIPELINE_EXPLAINED.md) + [MATHEMATICS](DOWNSCALING_MATHEMATICS.md) + [VENDOR_CODE](VENDOR_CODE_DEEP_DIVE.md)

**Algorithm developer** → [VENDOR_CODE](VENDOR_CODE_DEEP_DIVE.md) + [MATHEMATICS](DOWNSCALING_MATHEMATICS.md)

---

## Key Concepts at a Glance

### What is Statistical Downscaling?

```
Coarse-resolution model data (0.25°, ~25 km)
↓
Use fine-resolution observations (5 km) as training data
↓
Learn fine-scale spatial patterns
↓
Apply to coarse model → fine-resolution output (5 km)
```

**Not interpolation** — statistically corrects and enriches simulations.

### The Modified MBCn Algorithm

**MBCn** = Multivariate Bias Correction in N dimensions

**Modified** = adds spatial downscaling + weighted sum preservation

**Key steps**:
1. Rotate data to **weighted sum axis** → force conservation of totals
2. Apply **20 random rotations** → decorrelate variables, fix each independently
3. Use **non-parametric quantile mapping** → match observation statistics
4. Rotate **back to original space** → maintain correlations

**Result**: Data matches observations' marginal distributions + coarse model's physical totals + inter-variable dependence.

### Why It Works

| Issue | Solution |
|-------|----------|
| Univariate QM loses correlations | Rotate through multiple axes, fix each, rotate back |
| QM breaks energy/mass conservation | First axis is weighted sum; preserve it; let others adjust |
| Can't handle extreme values | Pre/post-randomization keeps values near physically realistic bounds |
| Computationally expensive | Month-by-month, coarse-cell-by-coarse-cell parallelisation |

---

## Pipeline at a Glance

```
Input                              Output
├─ obs_fine (HYRAS)               └─ sim_fine (BCSD)
│  9863×864×608 (30 years daily)      31390×864×608 (86 years daily)
│  2D rotated-pole coords
│
├─ sim_coarse (NEX-GDDP BC)
│  31390×32×38 (86 years daily)
│  1D regular lat/lon
│
└─ [Variables: pr, hurs, tas, etc.]

Stage 1: Harmonization (5 min)
├─ xarray → iris conversion
├─ Curvilinear → rectilinear reprojection (if HYRAS detected)
├─ Crop to divisible shape
├─ CRS/units stamping

Stage 2: Bilinear Remapbil (5 min, parallel)
├─ Year-by-year processing
├─ ThreadPoolExecutor (32 workers)
├─ iris.analysis.Linear regridding
└─ → sim_coarse_remapbil_cube [31390, 864, 608]

Stage 3: ISIMIP3BASD Statistical Downscaling (60–120 min)
├─ Per coarse cell (1,216)
├─ Per calendar month (12)
├─ Modified MBCn with 20 rotations
├─ Save to npy_stack/

Stage 4: Collection & Export (5 min)
├─ Load from npy_stack/
├─ Reshape to [31390, 864, 608]
├─ Assemble iris cube
└─ Save NetCDF

Total time per variable: 60–180 min (depending on n_processes, n_iterations)
```

---

## Terminology

| Term | Meaning | Example |
|------|---------|---------|
| **Downscaling** | Increasing spatial resolution | 25 km → 5 km |
| **Bias correction** | Removing systematic model errors | Applied before downscaling |
| **Remapping** | Regridding to a different grid | Bilinear interpolation |
| **Curvilinear grid** | Lat/lon varies per cell (rotated-pole) | HYRAS |
| **Rectilinear grid** | 1D lat/lon DimCoords | NEX-GDDP |
| **MBCn** | Multivariate quantile mapping algorithm | Core statistical method |
| **Quantile mapping** | Match CDF of two distributions | Non-parametric, via percentiles |
| **Copula** | Dependence structure between variables | Preserved via rotations |
| **CRE matrix** | Random orthogonal matrix | Used for random rotations |

---

## Code Files

| File | Role | Key Functions |
|------|------|---|
| `bcsd.py` | Wrapper & helpers | `StatisticalDownscaling.downscale()` |
| `bcsd.py` | Helper functions | `_stamp_geogcs()`, `_regrid_curvilinear_to_rectilinear()`, `_is_curvilinear()` |
| `statistical_downscaling.py` | Vendor core | `downscale()`, `downscale_one_location()`, `weighted_sum_preserving_mbcn()` |
| `utility_functions.py` | Vendor utils | `remapbil()`, `get_downscaling_factors()`, quantile mapping functions |

---

## Performance Summary

| Task | Time | Memory |
|------|------|--------|
| Load data | 2–5 min | 8–12 GB |
| Bias correction | 15–30 min | 12 GB |
| Downscale 1 var (32 threads) | 60–120 min | 2–4 GB (active) |
| Downscale 6 vars (sequential) | ~10 hrs | ~12 GB |
| **Total workflow** | **~12 hours** | **~12 GB** |

Parallel multi-variable processing reduces to ~3–4 hours.

---

## Troubleshooting Workflow

1. **Error during harmonization?** → [QUICK_REF Troubleshooting](DOWNSCALING_QUICK_REFERENCE.md#troubleshooting-common-errors)
2. **Need to find code location?** → [QUICK_REF References](DOWNSCALING_QUICK_REFERENCE.md#references-in-code)
3. **Want to understand algorithm?** → [VENDOR_CODE](VENDOR_CODE_DEEP_DIVE.md)
4. **Performance issues?** → [PIPELINE Performance](DOWNSCALING_PIPELINE_EXPLAINED.md#performance-characteristics)
5. **Modifying the code?** → [MATHEMATICS](DOWNSCALING_MATHEMATICS.md) for constraints/guarantees

---

## References & Links

- **ISIMIP3BASD Paper**: Lange, S. (2019). Trend-preserving bias adjustment and statistical downscaling with ISIMIP3BASD (v1.0). GMD, 12(7), 3055–3070. https://doi.org/10.5194/gmd-12-3055-2019

- **ISIMIP3BASD GitHub**: https://github.com/ISI-MIP/isimip3basd

- **iris Documentation**: https://scitools-iris.readthedocs.io/

- **xarray Documentation**: https://xarray.dev/

- **NEX-GDDP-CMIP6**: https://www.ncei.noaa.gov/products/nex-gddp-cmip6-downscaled-climate-simulations

- **HYRAS**: https://www.dwd.de/DE/leistungen/met_verfahren/mosmix/mosmix_stationskatalog.cfg?view=nasPublication

---

## Document Maintenance

These documents were created to explain:
- **Wrapper code**: `climdata/sdba/bcsd.py` (`StatisticalDownscaling` class)
- **Vendor code**: `climdata/_vendor/isimip3basd/` (statistical downscaling algorithm)
- **Notebook**: `usecase/nexgddp_hyras_bcsd.ipynb` (application example)

Last updated: 2024-12-27

---

## Getting Help

- **"What does this function do?"** → Read docstring + [QUICK_REF](DOWNSCALING_QUICK_REFERENCE.md)
- **"Why did downscaling fail?"** → [QUICK_REF Troubleshooting](DOWNSCALING_QUICK_REFERENCE.md#troubleshooting-common-errors)
- **"How do I modify X?"** → Read relevant section in [VENDOR_CODE](VENDOR_CODE_DEEP_DIVE.md)
- **"What's the theory?"** → [MATHEMATICS](DOWNSCALING_MATHEMATICS.md)

---

**Happy downscaling!** 🌍📈
