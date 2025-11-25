# climdata


[![image](https://img.shields.io/pypi/v/climdata.svg)](https://pypi.python.org/pypi/climdata)
[![image](https://img.shields.io/conda/vn/conda-forge/climdata.svg)](https://anaconda.org/conda-forge/climdata)

`climdata` is a Python package designed to automate fetching, extraction, and processing of climate data from various sources, including MSWX, DWD HYRAS, ERA5-Land, and NASA-NEX-GDDP. It provides tools to retrieve data for specific locations and time ranges, facilitating climate analysis and research.

---

## Key features
- Fetch and load datasets: MSWX, CMIP (cloud via intake), DWD, HYRAS
- Spatial extraction: point, box (region via config bounds), or shapefile (GeoJSON/Feature)
- Temporal subsetting via config or programmatic call
- Multi-format export: NetCDF, Zarr, CSV (standardized long format: variable, value, units)
- Hydra configuration + easy CLI overrides
- Helper to normalize AOI (GeoJSON → point / bbox / polygon)
- Provenance-friendly workflow (designed to be used with CI/CD workflows)

## Install (development)
1. Clone repository
```bash
git clone <repo-url>
cd climdata
```
2. Create virtualenv and install deps
```bash
python -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install -e ".[dev]"   # or pip install -r requirements.txt
```

## Quick CLI (Hydra) usage
Hydra reads configs from `conf/`. Override any config value on the CLI.

Examples:
```bash
# Region extraction (saves NetCDF by default when region is used)
python examples/climdata_cli.py dataset=CMIP region=europe time_range.start_date=2010-01-01 time_range.end_date=2010-12-31

# Point extraction (saves CSV)
python examples/climdata_cli.py dataset=MSWX lat=52.5 lon=13.4 variables=['tas','pr'] time_range.start_date=2000-01-01

# HYRAS / DWD (point only)
python examples/climdata_cli.py dataset=HYRAS lat=52 lon=10
```

Notes:
- Use `dataset=<MSWX|CMIP|DWD|HYRAS>` in CLI.
- Override any config key: e.g. `time_range.start_date=2000-01-01`.
- DWD/HYRAS: region (box) extraction is not supported — script will raise an error if attempted.

## Programmatic usage
Use the wrapper to compose configs, preprocess AOI, extract, and save.

```python
from climdata.utils.wrapper import extract_data

# returns (cfg, filename, ds, index) when save_to_file=True
cfg, filename, ds, index = extract_data(cfg_name="config", overrides=["dataset=MSWX","lat=52.5","lon=13.4"])
```

Or use the dataset classes directly:
```python
import climdata, xarray as xr
cmip = climdata.CMIP(cfg)
cmip.fetch()
cmip.load()
cmip.extract(box=cfg.bounds[cfg.region])
cmip.save_netcdf("output.nc")
```

## Configs
- Config files live in `climdata/conf/`. There are dataset-specific config entry points e.g. `config_cmip`, `config_mswx`, etc.  
- Filename templates are configurable in `cfg.output`:
  - `cfg.output.filename_nc`
  - `cfg.output.filename_csv`
  - `cfg.output.filename_zarr`

The wrapper generates filenames via `get_output_filename(cfg, output_type, ...)` using `cfg.bounds`, `cfg.time_range`, etc.

## Output CSV format
CSV produced by `save_csv` is standardized to the long form with columns (where available):
- source_id, experiment_id, table_id, time, lat, lon, variable, value, units

This ensures a single `value` column and a `variable` column for stacked variables.

## Common issues & tips
- NetCDF write ValueError (datetime encoding): call `ds["time"].encoding.clear()` before `to_netcdf()` (wrapper handles this).
- PermissionError writing files: ensure output directory is writable or write to `/tmp/` (or adjust permissions).
- CMIP cloud access requires network access — use the Pangeo intake catalog URL already referenced in code.

## AOI handling
`preprocess_aoi(cfg)` accepts:
- GeoJSON strings / Feature / FeatureCollection
- Point → sets `cfg.lat`, `cfg.lon`
- Polygon or bbox → sets `cfg.bounds['custom']` and `cfg.region='custom'`

## HYRAS support
HYRAS class mirrors MSWX design:
- `fetch()` / `load()` / `extract(point=...)` / `save_csv()` / `save_netcdf()`
- HYRAS extraction currently supports point extraction; attempt to use a region will raise an error.

## Development & provenance
- CI: add GitHub Actions workflows to run tests and build/publish to PyPI.
- Keep config and runtime overrides in Hydra to enable reproducible runs.
- Include `CITATION.cff`, license, and a changelog for FAIR discoverability.

## Contributing
- Run tests: `pytest`
- Style: follow repository linting config
- Open PRs against `main` with tests and a short changelog entry

## License
Specify the license (e.g. MIT or Apache 2.0) in `LICENSE`.

---

For further examples, see `examples/` and the `docs/` folder (usage, installation, faq).// filepath: /beegfs/muduchuru/pkgs_fnl/climdata/README.md
# climdata

Lightweight toolkit to fetch, subset and export climate data (MSWX, CMIP, DWD, HYRAS).  
Provides a Hydra-driven CLI, programmatic wrapper, cloud-native CMIP access, local dataset handling, and standardized CSV/NetCDF/Zarr exports.

## Key features
- Fetch and load datasets: MSWX, CMIP (cloud via intake), DWD, HYRAS
- Spatial extraction: point, box (region via config bounds), or shapefile (GeoJSON/Feature)
- Temporal subsetting via config or programmatic call
- Multi-format export: NetCDF, Zarr, CSV (standardized long format: variable, value, units)
- Hydra configuration + easy CLI overrides
- Helper to normalize AOI (GeoJSON → point / bbox / polygon)
- Provenance-friendly workflow (designed to be used with CI/CD workflows)

## Install (development)
1. Clone repository
```bash
git clone <repo-url>
cd climdata
```
2. Create virtualenv and install deps
```bash
python -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install -e ".[dev]"   # or pip install -r requirements.txt
```

## Quick CLI (Hydra) usage
Hydra reads configs from `conf/`. Override any config value on the CLI.

Examples:
```bash
# Region extraction (saves NetCDF by default when region is used)
python examples/climdata_cli.py dataset=CMIP region=europe time_range.start_date=2010-01-01 time_range.end_date=2010-12-31

# Point extraction (saves CSV)
python examples/climdata_cli.py dataset=MSWX lat=52.5 lon=13.4 variables=['tas','pr'] time_range.start_date=2000-01-01

# HYRAS / DWD (point only)
python examples/climdata_cli.py dataset=HYRAS lat=52 lon=10
```

Notes:
- Use `dataset=<MSWX|CMIP|DWD|HYRAS>` in CLI.
- Override any config key: e.g. `time_range.start_date=2000-01-01`.
- DWD/HYRAS: region (box) extraction is not supported — script will raise an error if attempted.

## Programmatic usage
Use the wrapper to compose configs, preprocess AOI, extract, and save.

```python
from climdata.utils.wrapper import extract_data

# returns (cfg, filename, ds, index) when save_to_file=True
cfg, filename, ds, index = extract_data(cfg_name="config", overrides=["dataset=MSWX","lat=52.5","lon=13.4"])
```

Or use the dataset classes directly:
```python
import climdata, xarray as xr
cmip = climdata.CMIP(cfg)
cmip.fetch()
cmip.load()
cmip.extract(box=cfg.bounds[cfg.region])
cmip.save_netcdf("output.nc")
```

## Configs
- Config files live in `climdata/conf/`. There are dataset-specific config entry points e.g. `config_cmip`, `config_mswx`, etc.  
- Filename templates are configurable in `cfg.output`:
  - `cfg.output.filename_nc`
  - `cfg.output.filename_csv`
  - `cfg.output.filename_zarr`

The wrapper generates filenames via `get_output_filename(cfg, output_type, ...)` using `cfg.bounds`, `cfg.time_range`, etc.

## Output CSV format
CSV produced by `save_csv` is standardized to the long form with columns (where available):
- source_id, experiment_id, table_id, time, lat, lon, variable, value, units

This ensures a single `value` column and a `variable` column for stacked variables.

## Common issues & tips
- NetCDF write ValueError (datetime encoding): call `ds["time"].encoding.clear()` before `to_netcdf()` (wrapper handles this).
- PermissionError writing files: ensure output directory is writable or write to `/tmp/` (or adjust permissions).
- CMIP cloud access requires network access — use the Pangeo intake catalog URL already referenced in code.

## AOI handling
`preprocess_aoi(cfg)` accepts:
- GeoJSON strings / Feature / FeatureCollection
- Point → sets `cfg.lat`, `cfg.lon`
- Polygon or bbox → sets `cfg.bounds['custom']` and `cfg.region='custom'`

## HYRAS support
HYRAS class mirrors MSWX design:
- `fetch()` / `load()` / `extract(point=...)` / `save_csv()` / `save_netcdf()`
- HYRAS extraction currently supports point extraction; attempt to use a region will raise an error.

## Development & provenance
- CI: add GitHub Actions workflows to run tests and build/publish to PyPI.
- Keep config and runtime overrides in Hydra to enable reproducible runs.
- Include `CITATION.cff`, license, and a changelog for FAIR discoverability.

## Contributing
- Run tests: `pytest`
- Style: follow repository linting config
- Open PRs against `main` with tests and a short changelog entry

## License
Specify the license (e.g. MIT or Apache 2.0) in `LICENSE`.

---

For further examples, see `examples/` and the `docs/` folder (usage, installation, faq).