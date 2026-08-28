# Implementation Plan: Resolution & Projection Transformation (`climdata.grid`)

**Status:** Implemented (PRs 1, 2 and 4 of the sequence in §12; PR 3 (dask) and PR 5
(dataset dedup) outstanding)
**Author:** Kaushik Muduchuru
**Date:** 2026-08-27

> **Implementation notes.** `climdata/grid/{units,crs,reproject}.py` are in place and
> `tests/test_grid.py` passes 66/66 against the `sdba` env. Two deviations from the plan
> as drafted, both found while building:
>
> 1. `str2pint` also accepts a *unit with no magnitude* -- `"km"` parses as 1 kilometre
>    and `""` as dimensionless 1. Neither is a resolution, so `_parse_scalar` now
>    requires the normalised string to start with a number.
> 2. The mismatch error echoes the caller's original string (`"10 km"`) rather than the
>    canonicalised form (`"10000 m"`), which reads far better in a traceback.
> 3. **HYRAS is EPSG:3035, not EPSG:31467**, and its CRS is not where the plan assumed.
>    See §7.1 -- found by running the real extractor, and it made `infer_src_crs` fail
>    outright on genuine HYRAS data until the inference chain was extended.
>
> One pre-existing issue surfaced during integration and is **not** fixed here:
> `_ensure_local_conf` (`climdata/utils/config.py`) copies the packaged `conf/` into the
> working directory only when no `conf/` exists, so a local copy never picks up new
> config keys after an upgrade. New options such as `target_projection` are invisible to
> Hydra until the stale `conf/` is refreshed. See §13.4.
**Scope:** Add first-class resolution and projection transformation to `climdata`, by
**composing existing libraries** — xclim/pint for units, rasterio/rioxarray for warping
and grid alignment, pyproj for CRS — and writing only the thin layer none of them provide.

---

## 1. Goals

1. Accept human-readable resolution strings — `"10 km"`, `"0.1 deg"`, `"30 arcsec"` —
   and parse units cleanly.
2. Accept any rasterio/pyproj-compatible projection — `"EPSG:4326"`, WKT, PROJ strings,
   integer EPSG codes, `pyproj.CRS`, `rasterio.crs.CRS`.
3. Wrap rasterio/rioxarray for reprojection and reproducible grid alignment.
4. Enforce the projection/resolution compatibility rules explicitly, and fail loudly and
   early when a requested combination is ill-defined.

### Non-goals

- Replacing the existing xESMF-based `regrid_to_coarse` (`climdata/sdba/bcsd.py:241`).
  It stays, and becomes reachable through the new API as an alternate engine.
- True area-weighted conservative remapping in the rasterio path. See §9.3.
- **Writing anything an existing dependency already does.** See §3.

---

## 2. Current state of the codebase

| Observation | Location |
|---|---|
| `regrid_to_coarse` exists — xESMF/CDO, template-driven, no CRS awareness, no resolution parsing | `climdata/sdba/bcsd.py:241` |
| rioxarray used ad-hoc in 5 datasets, each repeating `set_spatial_dims` + `write_crs("EPSG:4326")` inline, only for clipping | `climdata/datasets/W5E5.py:277`, `NEXGDDP.py:708`, `CMIPCloud.py:277`, `HYRAS.py`, `HOSTRADA.py` |
| HYRAS is the only non-4326 source. `HYRAS.py:263` writes `crs_grid = 'EPSG:31467'`, but the v6-1 product actually delivered is **EPSG:3035** -- see §7.1 | `climdata/datasets/HYRAS.py:263` |
| Per-dataset dim names already declared as `lat_dim`/`lon_dim` | `climdata/conf/mappings/parameters.yaml` |
| Workflow stages datasets via the `@update_ds(attr)` decorator | `climdata/utils/wrapper_workflow.py:70` |
| CMIP / NEX-GDDP arrive on 0–360 longitude | `climdata/datasets/CMIPCloud.py`, `NEXGDDP.py` |

### Verified environment

Probed in the `sdba` conda env (`/home/muduchuru/miniforge3/envs/sdba`), Python 3.10.17:

| Package | Version | Already a `climdata` dependency? |
|---|---|---|
| xclim | 0.56.0 | yes — `requirements.txt` |
| xsdba | 0.4.0 | yes |
| pint | 0.24.4 | yes (via xclim, and `pint-pandas`) |
| rioxarray | 0.18.2 | yes |
| rasterio | 1.4.3 | transitive (rioxarray) |
| pyproj | 3.7.1 | transitive (rioxarray/geopandas) |
| xesmf | 0.8.8 | yes |
| cf-xarray | 0.10.5 | transitive (xclim) |
| affine | 2.4.0 | transitive (rasterio) |
| odc-geo | **absent** | no — PyPI has 0.5.3 |

The base env lacks rasterio/rioxarray; `dev_env` has rioxarray 0.19.0 / rasterio 1.4.3 /
pyproj 3.7.2. Tests must therefore tolerate their absence (§12).

---

## 3. Build vs. reuse

An earlier draft of this plan hand-rolled a unit registry, a resolution parser, and grid
alignment arithmetic. Probing the installed stack showed almost all of that already
exists. **What survives as new code is small and specific.**

| Capability | Verdict | Provided by |
|---|---|---|
| Parse `"10 km"`, `"0.1 deg"`, `"30 arcsec"`, `"0.5 mi"`, `"30 ft"` | **reuse** | `xclim.core.units.str2pint` |
| Linear vs. angular discrimination | **reuse** | pint `Quantity.dimensionality` — `[length]` vs dimensionless |
| Unit conversion km→m, arcsec→deg, m→US survey foot | **reuse** | pint `Quantity.to()` |
| Parse `"EPSG:4326"` / WKT / PROJ / int / CRS objects | **reuse** | `pyproj.CRS.from_user_input` |
| CRS axis units & conversion factor | **reuse** | `pyproj.CRS.axis_info[i].unit_conversion_factor` |
| Destination transform + shape from a resolution | **reuse** | `rasterio.warp.calculate_default_transform(..., resolution=)` |
| **Grid alignment / origin snapping** | **reuse** | `rasterio.warp.aligned_target` |
| The warp itself, with nodata and resampling | **reuse** | `DataArray.rio.reproject(dst_crs, resolution=, transform=, shape=, resampling=, nodata=)` |
| Match an existing grid exactly | **reuse** | `rio.reproject_match` |
| Resampling method by name | **reuse** | `rasterio.enums.Resampling[name]` |
| Conservative (area-weighted) remapping | **reuse** | xESMF, via existing `regrid_to_coarse` |
| Dask-parallel reprojection | **reuse** | odc-geo (new optional dep) |
| — | — | — |
| **Projection ↔ resolution compatibility gate** | **BUILD** | nothing does this (§5) |
| **Unspaced/symbol input normalisation** (`"10km"`, `"0.1°"`) | **BUILD** — ~6 lines | `str2pint` rejects these (§6.1) |
| **Bare-number interception** | **BUILD** — critical | pint mis-converts it (§6.2) |
| **climdata-specific source-CRS inference** (`crs_grid` attr, dsinfo dims) | **BUILD** | repo-specific (§7) |
| **0–360 longitude normalisation** | **BUILD** — ~5 lines | repo-specific need (§7) |
| **Workflow / config / capability integration** | **BUILD** | repo-specific (§10) |

Note that **xclim itself does no reprojection and has no CRS handling.** It contributes
the *units* half only. The *projection* half is entirely rasterio/rioxarray/pyproj.

### Verified composition

The recommended call chain was executed end-to-end on a synthetic 4326 field over Germany,
reprojected to EPSG:3035 at 10 km:

```
calculate_default_transform -> 78 x 93, origin (3930307.6, 3569366.8)
aligned_target              -> 79 x 94, origin (3930000.0, 3570000.0)   # snapped
rio.reproject               -> EPSG:3035, resolution (10000.0, -10000.0)
```

Alignment invariant checked directly: the same field cropped differently, put through the
same chain, shared **72 of 72** cell centres. Stock `aligned_target` delivers the
reproducibility guarantee — no bespoke snapping arithmetic required.

---

## 4. Module layout

`climdata/grid/`, deliberately separate from `sdba` (which stays
statistical-downscaling-only). Smaller than the first draft, because §3 removed most of it:

| File | Responsibility | Approx. size |
|---|---|---|
| `climdata/grid/units.py` | `parse_resolution` (normalise → `str2pint` → classify), `ResolutionCRSMismatch`, the compatibility gate | ~120 lines |
| `climdata/grid/crs.py` | `parse_crs`, `infer_src_crs`, `infer_spatial_dims`, `normalize_longitude` | ~120 lines |
| `climdata/grid/reproject.py` | `reproject` — validate, compose `calculate_default_transform` + `aligned_target` + `rio.reproject` | ~150 lines |
| `climdata/grid/__init__.py` | Exports | — |

The separate `spec.py` / `GridSpec` from the first draft is **dropped**. Its entire job was
transform-and-shape arithmetic that `calculate_default_transform` and `aligned_target`
already do. A `NamedTuple` of `(crs, transform, width, height)` is enough to pass between
functions.

---

## 5. Projection ↔ resolution dependencies

This is the one piece of real logic that no dependency provides, so it gets an explicit
name and a single enforcement point.

### 5.1 The dependency chain

Resolution validity is **not** a property of `target_resolution` alone. It depends on the
*effective* target CRS, which is itself derived:

```
target_projection ─┐
                   ├─► effective_crs ─► axis unit kind ─┐
source CRS (fallback if target_projection is None)      ├─► COMPATIBILITY GATE
                                                        │
target_resolution ─────────► resolution kind ───────────┘
                                                        │
                                                        └─► calculate_default_transform
                                                            + aligned_target
                                                                 ▲
              regrid.bounds ───────────────────────────────────────┘
                            (must be expressed in effective_crs units)
```

Two consequences that must be **coded**, not merely documented:

1. **The gate runs on the effective CRS, after defaulting.** `target_projection=None` +
   `target_resolution="10 km"` on a 4326 source must fail identically to spelling out
   `target_projection="EPSG:4326"`. Validating the user's literal input would let this through.
2. **`like=` supersedes both.** A template carries CRS *and* resolution. Passing `like=`
   together with `target_resolution`/`target_projection` is contradictory and raises
   `ValueError`, rather than silently letting one win.

### 5.2 Compatibility matrix

| `target_resolution` | `target_projection` | Result |
|---|---|---|
| angular — `"0.1 deg"`, `"30 arcsec"` | geographic — 4326, 4258 | OK — pint `.to("degree")` |
| linear — `"10 km"`, `"30 ft"` | projected — 3035, 31467, 32632, 2263 | OK — pint `.to("metre")`, then ÷ `unit_conversion_factor` |
| **linear** | **geographic** | **ERROR — `ResolutionCRSMismatch`** |
| angular | projected | ERROR — `ResolutionCRSMismatch` (symmetric case) |
| unitless — `"0.1"`, `0.1` | either | OK — adopts CRS axis units verbatim (§6.2) |

The angular + projected case errors too. It is the same defect with the arguments swapped;
allowing one direction but not the other is an asymmetry that generates bug reports.

### 5.3 Why linear + geographic must be an error

A metric resolution has no fixed angular size. Measured with `pyproj.Geod(ellps='WGS84')`:

| 10 km at | Δlon | Δlat | anisotropy |
|---|---|---|---|
| 0°N | 0.0898° | 0.0904° | 1.01× |
| 50°N | 0.1395° | 0.0899° | 1.55× |
| 70°N | 0.2619° | 0.0896° | 2.92× |

A single `"10 km"` cannot describe that grid. There is no correct value to pick, so the
library must not pick one.

### 5.4 Unit derivation, positionally safe

```python
def crs_axis_unit(crs) -> tuple[str, float]:
    """-> ("angular"|"linear", conversion_factor_to_deg_or_m).
    Raises on a CRS whose horizontal axes carry different units."""
```

- Kind comes from `crs.is_geographic` / `crs.is_projected` — **never** from `axis_info[0]`.
  Verified trap: **EPSG:4326's `axis_info[0]` is latitude, not longitude.** Positional
  indexing silently mis-assigns anisotropic resolutions.
- Factor from `axis_info[i].unit_conversion_factor`. Verified:

  | CRS | geographic | axis units | factor |
  |---|---|---|---|
  | EPSG:4326 / 4258 | True | degree | 0.0174532925 |
  | EPSG:3035 / 3857 / 31467 / 32632 | False | metre | 1.0 |
  | EPSG:2263 | False | US survey foot | 0.3048006096 |

- Assert both horizontal axes share a unit; raise on the rare mixed-unit CRS.
- Compound/3D CRSs: reduce to the horizontal sub-CRS first.

### 5.5 Documented non-goal: Web Mercator

EPSG:3857 passes the gate (projected, metres) but is only true-scale at the equator —
`"10 km"` there is ~10 km of ground distance at 0° and ~3.4 km at 70°N. The gate cannot
catch this, so `reproject` emits a `UserWarning` when the target is Web Mercator and the
domain spans more than 20° of latitude. Warning, not error: 3857 is legitimate for tiling,
just not for area analysis.

---

## 6. Resolution parsing — `climdata/grid/units.py`

**Delegates to `xclim.core.units.str2pint`.** No hand-rolled unit tables, no new dependency —
xclim is already in `requirements.txt`. Verified working:

| Input | `str2pint` result | dimensionality |
|---|---|---|
| `"10 km"` | `Quantity(10.0, 'kilometer')` | `[length]` |
| `"1000 m"` | `Quantity(1000.0, 'meter')` | `[length]` |
| `"0.5 mi"` | `Quantity(0.5, 'mile')` | `[length]` |
| `"30 ft"` | `Quantity(30.0, 'foot')` | `[length]` |
| `"0.1 deg"` | `Quantity(0.1, 'degree')` | dimensionless |
| `"0.25 degree"` | `Quantity(0.25, 'degree')` | dimensionless |
| `"30 arcsec"` | `Quantity(30.0, 'arcsecond')` | dimensionless |
| `"5 arcmin"` | `Quantity(5.0, 'arcminute')` | dimensionless |

### 6.1 Gotcha — `str2pint` rejects unspaced input

Verified failures: `"10km"`, `"0.1°"`, `"0.25degree"` all raise
`ValueError: Unit expression cannot have a scaling factor.` Since `"10km"` is exactly what
gets typed on a Hydra CLI (where a space needs quoting), a normalisation shim runs first:

```python
_NUM = r'^([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)\s*'

def _normalize(s: str) -> str:
    s = s.strip().replace('°', ' deg')
    s = re.sub(r"'$", ' arcmin', s)      # anchored: trailing symbol only
    s = re.sub(r'"$', ' arcsec', s)
    return re.sub(_NUM, r'\1 ', s).strip()
```

Verified round-trips: `"10km"` → `10 km` → `Quantity(10.0, 'kilometer')`;
`"0.1°"` → `0.1 deg` → `Quantity(0.1, 'degree')`; `"30arcsec"` → `Quantity(30.0, 'arcsecond')`;
`"1e3 m"` → `Quantity(1000.0, 'meter')`.

The `'`/`"` substitutions are **anchored to end-of-string**. Unanchored replacement corrupts
input — verified: `"'x'"` became `arcminx arcmin`.

### 6.2 Gotcha — pint treats a bare number as radians

**This is the sharp edge.** In pint, degrees are *dimensionless*, so a bare number and an
angle share a dimensionality. Verified:

```
str2pint("0.1")   -> Quantity(0.1, 'dimensionless')  ->  .to("degree") == 5.7295779   # radians!
str2pint("0.1 deg") -> Quantity(0.1, 'degree')       ->  .to("degree") == 0.1
```

A bare `0.1` silently becomes 5.73° if passed through conversion. So classification must
check **units, not just dimensionality**, and unitless values must be intercepted *before*
any `.to()` call and adopted as CRS axis units verbatim:

```python
if q.dimensionality == "[length]":        kind = "linear"
elif str(q.units) in ("dimensionless", "1"): kind = "unitless"   # intercept — never .to()
else:                                     kind = "angular"       # deg, arcsec, arcmin, rad
```

### 6.3 The compatibility gate

```python
def resolution_in_crs_units(value, crs) -> tuple[float, float]:
    """Parse `value` and express it in the axis units of `crs`.

    Raises ResolutionCRSMismatch if the kinds are incompatible.
    """
```

This is the **single enforcement point**, called by `reproject`, so every entry point —
direct call, `ClimateExtractor.reproject`, Hydra override, LLM tool call — inherits the
check. It runs **before** any data is read or warped, so a bad config fails in milliseconds
rather than after a download.

Accepted input: `str`, `float`/`int`, a pint `Quantity`, or a 2-tuple of any of those for
anisotropic `(x, y)` resolution.

### 6.4 The exception

```python
class ResolutionCRSMismatch(ValueError):
    """target_resolution units are incompatible with target_projection axis units."""
```

Subclassing `ValueError` keeps existing `except ValueError` call sites working, including
the LLM execution layer in `climdata/LLM/execution.py`.

Message, with numbers computed at the domain's actual mid-latitude where a dataset is
available (falling back to the 50°N/70°N illustration otherwise):

```
ResolutionCRSMismatch: target_resolution="10 km" is a linear resolution, but
target_projection="EPSG:4326" is a geographic CRS with axis units of degree.

A metric resolution has no fixed angular size: 10 km spans 0.1395° of longitude
but 0.0899° of latitude at 50.0°N, and 0.2619° × 0.0896° at 70°N. There is no
single correct value for this grid.

Choose one:
  • angular resolution   target_resolution="0.1 deg"
  • projected CRS        target_projection="EPSG:3035"   (ETRS89-LAEA, metres — Europe)
                         target_projection="EPSG:32632"  (UTM 32N, metres)
```

### 6.5 Escape hatch — Python only, by design

```python
from climdata.grid import to_angular
res = to_angular("10 km", latitude=51.2)     # explicit, no error
reproject(ds, target_projection="EPSG:4326", target_resolution=res)
```

`to_angular()` requires an explicit `latitude` — no default — so the approximation cannot be
made accidentally. Implemented with `pyproj.Geod.fwd()`: walk `x` metres at azimuth 90° for
Δlon, `y` metres at azimuth 0° for Δlat. Exact at that latitude, and only there.

There is deliberately **no** `allow_approximate: true` config flag. A YAML toggle is exactly
how this would creep back into production runs unnoticed. The default path stays a hard error.

---

## 7. CRS handling — `climdata/grid/crs.py`

```python
def parse_crs(value) -> pyproj.CRS                     # thin wrapper, better errors
def infer_spatial_dims(obj, x_dim=None, y_dim=None) -> tuple[str, str]
def infer_src_crs(obj, default="EPSG:4326") -> pyproj.CRS
def normalize_longitude(obj, x_dim)
```

`parse_crs` delegates entirely to `pyproj.CRS.from_user_input`, which already covers
`"EPSG:4326"`, `4326`, WKT1/WKT2, PROJ strings, `"+proj=..."`, CF `grid_mapping` dicts,
`rasterio.crs.CRS`, and `pyproj.CRS`. The only value added is catching `CRSError` and
re-raising it naming the offending input.

`infer_src_crs` is genuinely climdata-specific. Resolution order:

1. `obj.rio.crs` (already set)
2. dataset-level attributes (`epsg_code`, `crs_wkt`, `spatial_ref`, `crs_grid`, ...)
3. a CF `grid_mapping` variable, **when it survived subsetting**
4. CRS attributes carried on the data variables themselves — notably `esri_pe_string`
5. a `spatial_ref` coordinate
6. heuristic: dims named lat/lon with values within ±90/±360 → EPSG:4326, with a warning
7. otherwise `ValueError` instructing the caller to pass `src_crs=`

### 7.1 What real HYRAS actually looks like

The plan originally assumed HYRAS was EPSG:31467 with its CRS in `attrs['crs_grid']`,
following `climdata/datasets/HYRAS.py:263`. Extracting January 2014 for real showed
otherwise. From `tas_hyras_1_2014_v6-1_de.nc`, the CF `crs` variable reads:

```
grid_mapping_name : lambert_azimuthal_equal_area
false_easting     : 4321000.0        latitude_of_projection_origin : 52.0
false_northing    : 3210000.0        longitude_of_central_meridian : 10.0
epsg_code         : EPSG:3035
long_name         : DWD HYRAS ETRS89 LAEA grid with 665 columns and 890 rows
```

The product is **ETRS89-LAEA (EPSG:3035)**, an equal-area projection in metres —
convenient, since that is also the natural target for European analysis.

Worse for inference: `ClimData.extract()` **drops the `crs` variable** during subsetting
while leaving `tas.attrs['grid_mapping'] = 'crs'` pointing at nothing. The only surviving
record of the projection is `tas.attrs['esri_pe_string']`, the ESRI WKT. Steps 3 and 4
above exist precisely for this: step 3 follows the pointer while the variable is still
there, step 4 reads the WKT off the variable when it is not.

The `crs_grid = 'EPSG:31467'` line at `HYRAS.py:263` is inconsistent with the delivered
data. It is left untouched here — it may serve an older or regional product — but it
warrants a separate look: anything trusting it will place HYRAS roughly 800 km from where
it belongs.

`infer_spatial_dims` tries `obj.rio.x_dim` first, then a name table
(`x, lon, longitude, rlon, nav_lon` / `y, lat, latitude, rlat, nav_lat`), then falls back to
the `cfg.dsinfo[DATASET].lat_dim/lon_dim` values the caller passes in.

`normalize_longitude` is **not optional polish.** CMIP and NEX-GDDP arrive on 0–360 and
rasterio will reproject that into garbage without raising. It rolls to −180…180 and
re-sorts; called unconditionally when the source CRS is geographic and `x.max() > 180`.

---

## 8. Reprojection wrapper — `climdata/grid/reproject.py`

```python
def reproject(
    obj,                        # xr.Dataset | xr.DataArray
    target_projection=None,     # "EPSG:4326" | WKT | PROJ | int | pyproj.CRS | rasterio CRS
    target_resolution=None,     # "10 km" | "0.1 deg" | float | pint Quantity | (x, y)
    *,
    like=None,                  # template Dataset/DataArray
    method="bilinear",          # str, or {var: str}
    bounds=None,
    align=True,
    src_crs=None,
    x_dim=None, y_dim=None,
    nodata=np.nan,
    engine="rasterio",          # or "xesmf" for true conservative remapping
)
```

Both primary options are leading and positional, in the order one would say them aloud.

### 8.1 Execution order

1. Coerce `target_projection` via `parse_crs`; infer `src_crs`; infer dims.
2. Resolve the **effective** target CRS (falls back to `src_crs`).
3. Reject contradictory `like=` + explicit target (`ValueError`).
4. `normalize_longitude` if geographic and needed.
5. `obj.rio.set_spatial_dims(...).rio.write_crs(src_crs)`.
6. `rio.write_nodata(np.nan, encoded=False)` on float vars — **without this rasterio fills
   gaps with 0**, silently corrupting masked domains such as HYRAS-over-Germany.
7. `resolution_in_crs_units(target_resolution, effective_crs)` —
   **this is where `ResolutionCRSMismatch` is raised.**
8. `calculate_default_transform(src_crs, dst_crs, width, height, *bounds, resolution=res)`.
9. If `align`: `aligned_target(transform, width, height, res)`.
10. Dispatch:
    - `like` given → `rio.reproject_match(template, resampling=...)`
    - otherwise → `rio.reproject(dst_crs, transform=..., shape=(h, w), resampling=..., nodata=...)`
11. Restore per-variable `attrs`/`units`/`cell_methods` (rioxarray drops some), write the CF
    `grid_mapping`, and rename output coords to `x`/`y` for projected CRSs and `lon`/`lat`
    for geographic ones so downstream `to_dataframe`
    (`climdata/utils/wrapper_workflow.py:791`) keeps working.

Steps 8–10 are library calls. The new code is steps 1–7 and 11.

### 8.2 Resampling methods

`Resampling[name]` — no hand-written mapping. Verified members: `nearest, bilinear, cubic,
cubic_spline, lanczos, average, mode, gauss, max, min, med, q1, q3, sum, rms`. Invalid names
raise `KeyError`, re-raised as `ValueError` listing valid options.

`method` also accepts a per-variable mapping, e.g. `{"pr": "average", "tas": "bilinear"}`.
Proposed: add an optional `regrid_method:` key per variable in
`climdata/conf/mappings/variables.yaml`, so the sensible default (fluxes → `average`,
states → `bilinear`, categorical → `nearest`) comes from config rather than each call site.

### 8.3 Conservative remapping — an honest limitation

rasterio's `average` is **not** area-weighted conservative remapping. For coarsening
precipitation it is a reasonable approximation on near-equal-area grids and wrong on a
lat/lon grid spanning many degrees of latitude.

Therefore `engine="xesmf"` delegates to the existing `regrid_to_coarse`
(`climdata/sdba/bcsd.py:241`), building the target template from the computed transform.
Keeping both engines behind one entry point is better than pretending rasterio covers the
conservative case — and reuses code that is already written and in production here.

### 8.4 Dask — use odc-geo, do not hand-roll

rioxarray's own `reproject` docstring says it plainly:

> To re-project with dask, see [odc-geo](https://odc-geo.readthedocs.io/) & pyresample.

So the first draft's bespoke `map_blocks` scheme is dropped. Plan: for dask-backed inputs,
delegate to `odc.geo.xr.xr_reproject`, added as an **optional** dependency
(`climdata[grid]`), with a clear error if a chunked array is passed without odc-geo
installed. Non-dask inputs continue through plain `rio.reproject`.

odc-geo is not currently installed in any local env (PyPI has 0.5.3), so this path needs an
install and a benchmark before PR 3 lands. It remains the one item to verify empirically.

---

## 9. Integration

### 9.1 Public API — `climdata/__init__.py`

```python
from .grid import reproject, parse_crs, parse_resolution, to_angular, ResolutionCRSMismatch
```

### 9.2 Workflow method — `climdata/utils/wrapper_workflow.py`

New `ClimateExtractor.reproject()`, decorated `@update_ds("reprojected_ds")`, mirroring the
shape of `impute()` (`climdata/utils/wrapper_workflow.py:1748`): read config, fall back to
`self.current_ds`, return `None` with a log line when no target grid is configured. Add
`reprojected_ds` / `reprojected_filename` to `WorkflowResult`
(`climdata/utils/wrapper_workflow.py:57`).

### 9.3 Config — `climdata/conf/config.yaml`

The two primary options sit at **top level**, matching the existing flat `index:` / `impute:`
convention:

```yaml
target_projection: null      # e.g. "EPSG:3035"
target_resolution: null      # e.g. "10 km" (projected) or "0.1 deg" (geographic)

regrid:                      # secondary knobs; mirrors extinfo / imputeinfo
  method: bilinear
  align: true
  bounds: null
  engine: rasterio
```

All null, so existing runs are untouched. Overrides read naturally:

```python
ClimData(overrides=["target_projection=EPSG:3035", "target_resolution=10 km"])
```

**Hydra detail:** `"10 km"` contains a space, so CLI use needs quoting —
`climdata 'target_resolution=10 km'`. The §6.1 shim also accepts `"10km"` unspaced, which is
the reason that shim exists. Docs should show the unspaced form for CLI, spaced for Python.

### 9.4 Registration

- An entry in `climdata/conf/mappings/actions.yaml` (`type: primary`, `output: dataset`).
- A `Capability` in `climdata/LLM/capabilities.py` alongside the existing `regrid_to_coarse`
  entry (`climdata/LLM/capabilities.py:344`), under a new `CapabilityKind.REGRIDDING` rather
  than reusing `BIAS_CORRECTION`, which is a misfile in the current code. `ParamSpec` names
  use `target_projection` / `target_resolution` so the planner emits working overrides.

### 9.5 Dataset cleanup (separate PR)

Replace the five inline `set_spatial_dims` / `write_crs` blocks with the
`climdata/grid/crs.py` helpers. Pure deduplication, no behaviour change.

---

## 10. Dependencies

**No new required dependency.** Everything in §3's reuse column is already installed:
xclim, pint, rioxarray, rasterio, pyproj, xesmf.

Changes to `requirements.txt`: name `rasterio`, `pyproj`, and `affine` explicitly, since the
new module imports them directly rather than relying on them transitively.

New optional extra in `pyproject.toml`:

```toml
[project.optional-dependencies]
grid = ["odc-geo>=0.4"]     # dask-parallel reprojection only
```

---

## 11. Tests — `tests/test_grid.py`

**Test our code, not our dependencies.** pint's unit table, pyproj's CRS parsing and
rasterio's warp arithmetic are upstream-tested; duplicating that here buys nothing. The
suite targets the shim, the gate, the inference chain, and the integration.

rasterio-dependent tests are guarded with `pytest.importorskip("rioxarray")` so the suite
still runs in the base env, where those libraries are absent.

| # | Test | Asserts |
|---|---|---|
| 1 | Normalisation shim (§6.1) | `"10km"`, `"0.1°"`, `"30arcsec"`, `"0.25degree"`, `"1e3 m"` parse; anchored `'`/`"` do not corrupt `"'x'"` |
| 2 | Kind classification (§6.2) | linear / angular / unitless assigned correctly; **`"0.1"` is unitless, never 5.73°** |
| 3 | Compatibility gate (§5.2) | All kind×kind combinations; `ResolutionCRSMismatch` on exactly the two invalid ones |
| 4 | Effective-CRS path | `target_projection=None` on a 4326 source with `"10 km"` raises |
| 5 | Error message quality | Contains the resolved mid-latitude and both suggested fixes |
| 6 | Escape hatch | `to_angular("10 km", latitude=50)` → 0.1395° / 0.0899° within 1e-4 of a direct `Geod` call; missing `latitude` raises `TypeError` |
| 7 | Non-metre linear CRS | `"10 km"` → EPSG:2263 yields 32808.4 US survey feet |
| 8 | Contradictory inputs | `like=` plus `target_resolution` raises `ValueError` |
| 9 | Axis-order safety (§5.4) | Anisotropic resolution into EPSG:4326 (lat-first axis order) maps to the correct axes |
| 10 | Known-answer reprojection | Synthetic 4326 → EPSG:3035 at 10 km gives 79×94, origin (3930000, 3570000) — the values verified in §3 |
| 11 | **Alignment invariant** | Two differently-cropped inputs → identical shared cell centres (verified achievable: 72/72) |
| 12 | Idempotence | Reproject to native CRS + native resolution ≈ input within tolerance |
| 13 | Nodata | Masked input does not acquire zeros in the masked region |
| 14 | Longitude wrap | 0–360 input produces correctly ordered −180…180 output |
| 15 | Conservation smoke test | Constant field coarsened with `average` stays that constant |
| 16 | Integration | Two HYRAS-shaped fixtures: the real v6-1 shape (EPSG:3035, dangling `grid_mapping`, CRS only in `esri_pe_string`) and the `crs_grid` variant, both → EPSG:4326, exercising the full inference chain |
| 17 | grid_mapping pointer | A surviving CF `crs` variable is followed to its `epsg_code` |

Tests 10 and 11 assert values already confirmed against the real library stack (§3), so they
are regression tests from day one rather than aspirations.

---

## 12. Sequencing

| PR | Content | Depends on |
|---|---|---|
| 1 | `units.py` + `crs.py` + tests 1–7, 14 | — |
| 2 | `reproject.py` (rasterio path) + tests 8–13, 15, 16; requirements | 1 |
| 3 | odc-geo dask path + benchmark + `climdata[grid]` extra | 2 |
| 4 | `ClimateExtractor.reproject`, config, actions.yaml, capabilities | 2 |
| 5 | Dataset dedup, docs page, `regrid_to_coarse` delegation | 4 |

PR 1 is pure and runs without a geo stack installed. Merging `spec.py` into `reproject.py`
(§4) collapses the first draft's five PRs' worth of scope into a smaller PR 2.

---

## 13. Open decisions

1. **`regrid_to_coarse`.** It is exported from `climdata/sdba/__init__.py` and registered as
   an LLM capability, so it is public surface. Preference: keep it, and have `engine="xesmf"`
   call into it (§8.3) — reuse rather than deprecate.
2. **Angular resolution + projected CRS.** Currently planned as an error, for symmetry with
   the linear + geographic case. Could be made permissive if a real use case appears.
3. **odc-geo as an optional extra vs. required.** Optional keeps the base install light;
   required makes the dask path always available. Needs the §8.4 benchmark first.
4. **Stale local `conf/`.** `_ensure_local_conf` never refreshes an existing working-
   directory `conf/`, so config keys added by an upgrade are invisible to Hydra and any
   override naming them fails with "Key not in struct". Options: merge missing keys on
   load, version the copied tree, or document a manual refresh. Out of scope for this
   change, but it gates discoverability of every new config option, not just these two.
