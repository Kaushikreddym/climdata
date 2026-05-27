# Downscaling: Mathematical & Algorithmic Reference

## Problem Statement

**Goal**: Given
- $\mathbf{x}_{\text{obs}}^{\text{fine}}$ — fine-resolution observations (training data)
- $\mathbf{y}_{\text{sim}}^{\text{coarse}}$ — coarse-resolution model simulations (to correct)
- $\mathbf{y}_{\text{sim,\,remap}}^{\text{fine}}$ — coarse simulations bilinearly interpolated to fine grid

Find:
- $\hat{\mathbf{y}}_{\text{sim}}^{\text{fine}}$ — model data at fine resolution matching observational statistics

**Constraints**:
- Marginal distributions: $\text{CDF}_{\hat{\mathbf{y}}_i} \approx \text{CDF}_{\mathbf{x}_{\text{obs}, i}}$ (per fine cell $i$)
- Weighted sum: $\sum_i w_i \hat{y}_i \approx \sum_j y_{\text{coarse}, j}$ (energy/mass conservation)
- Multivariate structure: correlations between fine cells preserved

---

## Stage 1: Bilinear Remapping

### Problem
Need to go from coarse grid $(32 \times 38)$ → fine grid $(864 \times 608)$

$$\mathbf{y}_{\text{sim}}^{\text{coarse}} \in \mathbb{R}^{T \times 32 \times 38} \to \mathbf{y}_{\text{sim,\,remap}}^{\text{fine}} \in \mathbb{R}^{T \times 864 \times 608}$$

### Solution: iris.analysis.Linear (RectilinearRegridder)

For each timestep $t$, for each fine cell $(i_f, j_f)$:

1. **Find containing coarse cell** $(i_c, j_c)$

2. **Bilinear interpolation** using four surrounding coarse values:

$$\hat{y}(i_f, j_f) = y(i_c, j_c) \cdot (1-u)(1-v) + y(i_c+1, j_c) \cdot u(1-v)$$
$$+ y(i_c, j_c+1) \cdot (1-u)v + y(i_c+1, j_c+1) \cdot uv$$

where $u, v \in [0, 1]$ are the fractional position within the coarse cell.

3. **Fallback**: Where bilinear fails (masked coarse cells), use nearest-neighbour.

### Result
$$\mathbf{y}_{\text{sim,\,remap}}^{\text{fine}}[t] = \text{Regrid}_{\text{bilinear}}(\mathbf{y}_{\text{sim}}^{\text{coarse}}[t], \text{grid}_{\text{fine}})$$

---

## Stage 2: Modified MBCn Algorithm

### Core Problem
Transform $\mathbf{y}_{\text{sim,\,remap}}^{\text{fine}}$ to match observational statistics while preserving coarse model's physical total.

### Algorithm Overview

**Inputs** (for one coarse cell, one month):
- $\mathbf{x} \in \mathbb{R}^{M \times N}$ — observations at $N$ fine cells over $M$ days
- $y_{\text{coarse}} \in \mathbb{R}^{M}$ — scalar coarse model value
- $\mathbf{y}_{\text{remap}} \in \mathbb{R}^{M \times N}$ — bilinearly interpolated model at fine cells
- $\mathbf{w} \in \mathbb{R}^{N}$ — area weights (normalized: $\|\mathbf{w}\|_2 = 1$)

**Output**:
- $\hat{\mathbf{y}} \in \mathbb{R}^{M \times N}$ — downscaled model at fine cells

### The Weighted Sum Axis

Define the **weighted sum** (linear combination) of all fine cells:

$$s_i = \sum_{j=1}^{N} w_j y_{i,j}$$

where $i$ ranges over timesteps.

**Key insight**: If we rotate to a coordinate system where the first axis is this weighted sum direction, we can:
1. Force simulations to preserve the coarse model's sum
2. Independently adjust other axes to match observations

### Rotation Matrix for Sum Axis

Given normalized weights $\mathbf{w}$, compute orthogonal matrix $\mathbf{O}_0$ such that:

$$\mathbf{e}_1^T \mathbf{O}_0^T = \mathbf{w}^T$$

i.e., the first row of $\mathbf{O}_0$ is $\mathbf{w}$.

This can be constructed via QR decomposition or Householder reflections.

### MBCn Main Loop

```
Initialize: O_total = I_N (identity)
For iteration i = 0 to n_iterations + 1:
    1. Choose rotation matrix O_i
    2. Rotate: X ← X O_i, Y ← Y O_i, w ← w O_i
    3. Apply quantile mapping on rotated axes
    4. Accumulate: O_total ← O_total @ O_i
```

#### Iteration 0: To Sum Axis

$$\mathbf{O}_0 = \begin{bmatrix} \mathbf{w}^T \\ \vdots \\ \text{(other rows orthogonal to } \mathbf{w}) \end{bmatrix}$$

After rotation:
- Axis 0 is the weighted sum direction
- $y_{i,0}^{\text{rot}} = \sum_j w_j y_{i,j}$ for all rows

**Apply quantile mapping on axis 0 only**:

$$y_{i,0}^{\text{mapped}} = Q_{\text{coarse}}^{-1}(F_{\text{obs, axis 0}}(x_{i,0}^{\text{rot}}))$$

where:
- $F_{\text{obs, axis 0}}$ = empirical CDF of observed axis-0 values
- $Q_{\text{coarse}}^{-1}$ = inverse quantile of coarse model

**Force preservation**:
$$y_{i,0}^{\text{mapped}} := \text{sum}(y_{\text{coarse}, i} \cdot \mathbf{w})$$

i.e., keep the coarse model's sum exactly.

#### Iterations 1 to $n_{\text{iter}}-1$: Random Rotations

Choose random orthogonal matrix $\mathbf{O}_i$ (from circular real ensemble):

$$\mathbf{O}_i \sim \text{CRE}(N)$$

**Apply univariate quantile mapping per axis**:

For axis $j = 0, 1, \ldots, N-1$:

$$y_{:,j}^{\text{mapped}} = Q_{\text{obs}, j}^{-1}(F_{\text{sim}, j}(y_{:,j}^{\text{rot}}))$$

**Preserve weighted sum**:

Compute how much the sum changed:

$$\Delta s_i = \sum_j w_j (y_{i,j}^{\text{mapped}} - y_{i,j}^{\text{rot}})$$

Distribute correction back proportionally:

$$y_{i,j}^{\text{corrected}} = y_{i,j}^{\text{mapped}} - w_j \Delta s_i$$

#### Iteration $n_{\text{iter}}$: Back to Original Axes

$$\mathbf{O}_{n_{\text{iter}}} = \mathbf{O}_{\text{total}}^T$$

This cancels out all previous rotations.

### Non-Parametric Quantile Mapping

For axis $j$, map $\mathbf{y}$ from its distribution to match $\mathbf{x}$:

1. **Compute empirical quantiles**:
   - $\{q^{(y)}_k = \text{quantile}(\mathbf{y}_j, k/(n_q+1))\}_{k=0}^{n_q+1}$
   - $\{q^{(x)}_k = \text{quantile}(\mathbf{x}_j, k/(n_q+1))\}_{k=0}^{n_q+1}$
   
   (typically $n_q = 50$, giving 51 quantile pairs)

2. **Map each value**:
   $$y_i^{\text{mapped}} = \text{interp}(y_i, q^{(y)}, q^{(x)}, \text{left}=q^{(x)}_0, \text{right}=q^{(x)}_{n_q+1})$$
   
   i.e., linear interpolation with constant extrapolation

**Properties**:
- ✅ Monotonic: order is preserved
- ✅ Non-parametric: no distribution assumption
- ✅ Handles extremes: constant extrapolation outside quantile range

### Mathematical Properties

**After full MBCn loop**:

1. **Univariate matching** (via quantile mapping):
   $$\text{CDF}_{\hat{y}_j}(x) \approx \text{CDF}_{\mathbf{x}_j}(x) \quad \forall j$$

2. **Weighted sum preservation** (via correction step):
   $$\sum_j w_j \hat{y}_{i,j} = \sum_j w_j y_{\text{remap}, i,j} \quad \forall i$$

3. **Copula structure** (via multiple rotations):
   - Each rotation tests a different linear combination
   - Correlations between axes gradually adjusted
   - 20 rotations → rich dependence structure

---

## Stage 3: Biased Correction Handling

### Pre-Processing: Invalid Value Imputation

For each invalid (NaN) value in $\mathbf{y}$ or $\mathbf{x}$:

$$x_{i,j}^{\text{imputed}} = \text{sample from } N(\bar{x}_j, \sigma_j^2)$$

where $\bar{x}_j, \sigma_j$ are computed from valid values.

### Censoring & Randomization

Some variables have physical bounds, e.g.:
- Precipitation: $p \geq 0$
- Humidity: $0 \leq h \leq 100\%$

**Before MBCn**: Replace censored values (too close to bounds) with random noise:

$$y_i^{\text{randomized}} = \begin{cases}
\text{Uniform}(y_{\min}, y_{\text{threshold}}) & \text{if } y_i < y_{\text{threshold}} \\
y_i & \text{otherwise}
\end{cases}$$

**Why?** MBCn's quantile mapping can amplify values at quantile extremes. Pre-randomization avoids unrealistic concentrations at bounds.

### Post-Processing: De-Randomization & Clipping

After MBCn, clip back to physical bounds:

$$\hat{y}_i^{\text{final}} = \begin{cases}
y_{\max} & \text{if } \hat{y}_i > y_{\max} \\
y_{\min} & \text{if } \hat{y}_i < y_{\min} \\
\hat{y}_i & \text{otherwise}
\end{cases}$$

---

## Downscaling Factors

### Definition

Given:
- Fine grid shape: $(N_{\text{lat,fine}}, N_{\text{lon,fine}}) = (864, 608)$
- Coarse grid shape: $(N_{\text{lat,coarse}}, N_{\text{lon,coarse}}) = (32, 38)$

Downscaling factors:

$$f_{\text{lat}} = N_{\text{lat,fine}} / N_{\text{lat,coarse}} = 864 / 32 = 27$$
$$f_{\text{lon}} = N_{\text{lon,fine}} / N_{\text{lon,coarse}} = 608 / 38 = 16$$

**Constraint**: Must divide evenly (no remainder).

### Fine Cell Mapping

For coarse cell $(i_c, j_c)$, corresponding fine cells are:

$$\{(i_f, j_f) : i_f \in [i_c \cdot f_{\text{lat}}, (i_c+1) \cdot f_{\text{lat}}), j_f \in [j_c \cdot f_{\text{lon}}, (j_c+1) \cdot f_{\text{lon}})\}$$

**Total fine cells per coarse cell**:

$$N_{\text{fine|coarse}} = f_{\text{lat}} \times f_{\text{lon}} = 27 \times 16 = 432$$

---

## Area-Weighted Statistics

### Cosine Latitude Weighting

For a regular lat/lon grid, true grid cell areas vary by latitude due to meridian convergence:

$$A(\phi) \propto \cos(\phi)$$

Normalized weights:

$$w_i = \frac{\cos(\phi_i)}{\sum_j \cos(\phi_j)}$$

**Used for**:
- Computing weighted averages (spatial mean)
- Preserving total physical quantity (e.g., total precipitation)
- Fair contribution of polar vs. equatorial cells

### Example Weights (HYRAS 864 cells in latitude)

```
Latitude     Cos(φ)    Weight
90° (pole)     0       0.0
60°          0.5      ~0.012
30°        0.866      ~0.020
0° (equator)   1       ~0.024
```

Equator cells get ~2× more weight than 60° cells, reflecting their larger true area.

---

## Computational Complexity

### Time Complexity

| Step | Operations | Time |
|------|-----------|------|
| Bilinear remapbil | $T \times N_{\text{coarse}} \times$ (local interp) | $O(T \cdot N_c)$ |
| Per coarse cell | 1,216 iterations | — |
| Per month | 12 iterations | — |
| Per MBCn iteration | $n_{\text{iter}}+2$ rotations, $N^2$ matrix mult | $O(N^2)$ |
| Quantile mapping | $M \cdot N \cdot \log(M)$ | $O(M N \log M)$ |
| **Total for 1 variable** | | **~30–120 min** |

### Memory Complexity

| Array | Shape | Size (bytes) |
|-------|-------|---|
| obs_fine | (9863, 864, 608) | ~41 GB |
| sim_coarse | (31390, 32, 38) | ~36 MB |
| sim_coarse_remapbil | (31390, 864, 608) | ~120 GB |
| Rotation matrix | (432, 432) | ~1.5 MB × 20 rotations |
| **Peak** | | **~161 GB (all at once)** |
| **Per year** (threading) | (365, 864, 608) | **~0.9 GB** |

Threading reduces peak to ~10–20 GB by processing years sequentially.

---

## Notation Summary

| Symbol | Meaning | Shape |
|--------|---------|-------|
| $\mathbf{x}$ | Observations (training) | $(M_{\text{obs}}, N_{\text{fine}})$ |
| $\mathbf{y}_{\text{sim}}^{\text{coarse}}$ | Simulations at coarse resolution | $(M_{\text{sim}}, N_{\text{coarse}})$ |
| $\mathbf{y}_{\text{remap}}^{\text{fine}}$ | Sim bilinearly interpolated to fine | $(M_{\text{sim}}, N_{\text{fine}})$ |
| $\hat{\mathbf{y}}^{\text{fine}}$ | Downscaled output | $(M_{\text{sim}}, N_{\text{fine}})$ |
| $\mathbf{w}$ | Area-weighted cosine lat weights | $(N_{\text{fine}},)$ |
| $\mathbf{O}$ | Orthogonal rotation matrix | $(N_{\text{fine}}, N_{\text{fine}})$ |
| $Q$ | Quantile function (inverse CDF) | — |
| $F$ | Empirical CDF | — |

---

## References

- **Quantile mapping**: Themeßl et al. (2012). Empirical-statistical downscaling and error correction of regional climate models and its impact on the climate change signal. Climatic Change, 112(2), 449–468.

- **MBCn algorithm**: Cannon et al. (2015). Multivariate quantile mapping bias correction. Journal of Hydrology, 529, 753–769.

- **ISIMIP3BASD**: Lange, S. (2019). Trend-preserving bias adjustment and statistical downscaling with ISIMIP3BASD (v1.0). GMD, 12(7), 3055–3070.

- **CRE matrices**: Mezzadri, F. (2007). How to generate random matrices from the classical compact groups. arXiv:math-ph/0609050.
