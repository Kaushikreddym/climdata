# Copilot Instructions for climdata

## Build & Install

```bash
# Editable install (development)
pip install -e .

# With all optional ML/imputation extras
pip install -e ".[full]"

# PyTorch-based imputation (install separately)
pip install torch==2.3.1 --index-url https://download.pytorch.org/whl/cpu
pip install torch-cluster -f https://data.pyg.org/whl/torch-2.3.1+cpu.html
```

## Testing

```bash
# Run the full test suite (excludes dataset-integration tests)
python -m unittest discover -p 'test_[!d]*.py' tests/

# Run a single test file
python -m unittest tests/test_climdata.py

# Run a single test
python -m unittest tests.test_climdata.TestClimdata.test_version
```

> `tests/test_datasets_extraction.py` and `tests/test_point_workflow.py` require HPC data paths and PyTorch — they are skipped in CI (`CI` env var set).

## Linting

```bash
flake8 climdata/
```

Max line length is 88 (black-compatible). Config is in `[tool.flake8]` in `pyproject.toml`.

## Version Bumping

```bash
bump2version patch   # or minor / major
```

Updates `pyproject.toml` and `climdata/__init__.py` in sync (see `[tool.bumpversion]`).

---

## Architecture

### Central class: `ClimData` / `ClimateExtractor`

`climdata.ClimData` is an alias for `climdata.utils.wrapper_workflow.ClimateExtractor`. It is the entry point for all data extraction, index computation, and I/O.

```python
from climdata import ClimData

extractor = ClimData(overrides=["dataset=mswx", "lat=52.5", "lon=13.4", ...])
ds = extractor.extract()
ds_idx = extractor.calc_index(ds)
df = extractor.to_dataframe(ds_idx)
extractor.to_csv(df)

# or chain steps:
result = extractor.run_workflow(actions=["extract", "impute", "calc_index", "to_nc"])
```

`run_workflow` returns a `WorkflowResult` dataclass holding `dataset`, `dataframe`, `filename`, `index_ds`, `impute_ds`, etc.

### Configuration system (Hydra)

All configuration is managed by **Hydra**. The config is loaded from `climdata/conf/config.yaml` (copied to the CWD on first run by `climdata.utils.config._ensure_local_conf`).

Config sub-trees (all in `climdata/conf/mappings/`):
- `parameters.yaml` — per-dataset access params (folder IDs, API keys, variable mappings)
- `variables.yaml` — CF-standard variable definitions
- `indices.yaml` — extreme climate index definitions (dispatched to xclim)
- `actions.yaml` — workflow action definitions
- `imputation.yaml` — imputation method configurations

Key top-level config fields:
```yaml
dataset: MSWX          # which provider to use
lat: null              # point extraction
lon: null
aoi: null              # GeoJSON string for polygon/point extraction
shapefile: null
variables: [tasmin, tasmax, pr, ...]
data_dir: ./data
time_range: {start_date, end_date}
index: null            # e.g. "tn10p"
impute: null           # e.g. "BRITS"
experiment_id: historical   # for CMIP
source_id: MIROC6
output: {out_dir, filename_csv, filename_zarr, filename_nc, fmt}
```

Overrides are passed as Hydra-style strings: `"dataset=cmip"`, `"lat=52.5"`, `"variables=[tasmin,pr]"`.

### Dataset providers

Each provider is a class in `climdata/datasets/` accepting a `DictConfig` and exposing `extract()` / `fetch()` methods. Providers exported from `climdata/__init__.py`:

| Import name | Class | Notes |
|---|---|---|
| `MSWX` | `MSWXmirror` | Requires Google Drive service account JSON (`conf/service.json`) |
| `DWD` | `DWDmirror` | Uses `wetterdienst` library; returns nearest station data |
| `ERA5` | `ERA5Mirror` | Requires `~/.cdsapirc` with CDS API key |
| `CMIP` | `CMIPCloud` | CMIP6 from cloud storage |
| `W5E5` | `W5E5` | |
| `CMIPW5E5` | `CMIPW5E5` | Bias-corrected CMIP with W5E5 baseline |
| `NEXGDDP` | `NEXGDDP` | NASA NEX-GDDP |
| `HYRAS` | `HYRASmirror` | DWD HYRAS gridded dataset |
| `HOSTRADA` | `HOSTRADAmirror` | |
| `POWER` | `POWER` | NASA POWER |
| `AGRI_ISIMIP` | `AGRI_ISIMIP` | |

### Variable naming

Variables use **CF standard names** throughout: `tas`, `tasmin`, `tasmax`, `pr`, `rsds`, `rlds`, `sfcWind`, `hurs`, `huss`, `ps`. The `CF_TO_DWD_NAMES` dict in `wrapper_workflow.py` maps these to DWD internal names where needed.

### Extreme indices

Indices are declared in `climdata/conf/mappings/indices.yaml` and dispatched dynamically to `xclim` functions by `climdata.extremes.indices.extreme_index`. The `link` key in an index config allows intermediate variables to be computed via `function_call` or `operation` before calling the main index function.

### Decorators on `ClimateExtractor` methods

- `@update_ds(attr_name)` — after the wrapped method returns an `xr.Dataset`, stores it as `self.current_ds` (and optionally `self.<attr_name>`), then regenerates output filenames via `self._gen_fn_cfg()`.
- `@update_df(attr_name)` — same pattern for `pd.DataFrame` → `self.current_df`.

### Explore / catalog

`climdata.explore` provides `list_available_data`, `explore`, `find`, `inspect`, `list_esm_experiments`, `list_esm_models`, and `DatasetRegistry` for browsing available datasets and their metadata without extracting data.

### Vendored dependencies

`climdata/_vendor/` contains vendored copies of:
- `imputegap` — gap-filling / imputation algorithms (ML and statistical methods)
- `isimip3basd` — ISIMIP3 bias adjustment and statistical downscaling

### GUI

`climdata_gui/` is a separate PySide6 desktop application. It is packaged with PyInstaller (see `climdata_gui.spec`). The GUI bundles the `climdata` package via `collect_all('climdata')` in the spec file.
