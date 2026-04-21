# climdata.explore — Dataset Discovery & Browsing

A professional, intuitive fluent interface for browsing climate datasets without needing to know underlying file paths.

## Quick Start

```python
import climdata as cd

# See all available datasets
cd.list_available_data()

# Explore a specific dataset
cd.explore(dataset="NEXGDDP")

# Search datasets by criteria
cd.find(variable="pr", frequency="daily")
cd.find(type_filter="ESM")

# Check variable metadata and BASD compatibility
cd.inspect("ERA5", variable="tp")

# ESM-specific discovery
cd.list_esm_experiments("CMIP")        # Available scenarios/experiments
cd.list_esm_models("NEXGDDP", experiment="historical")  # Available models
```

## Core Concepts

### 1. **High-Level Overview** (`list_available_data()`)
Prints a formatted table of all 9 climate datasets in climdata:

```
AVAILABLE CLIMATE DATASETS
────────────────────────────────────────────────────────────────────
Abbr.    | Long Name                          | Type
────────────────────────────────────────────────────────────────────
ERA5     | ECMWF Reanalysis v5                | Observation
MSWX     | Multi-Source Weather – GloH2O      | Observation
CMIP     | CMIP6 via Pangeo Cloud Catalog     | ESM (Raw)
NEXGDDP  | NASA Earth Exchange GDDP           | ESM (Bias-Corrected)
... (5 more datasets)
```

### 2. **Deep-Dive into a Dataset** (`explore()`)
```python
cd.explore(dataset="NEXGDDP")
```

Outputs:
- Type, coverage, resolution, frequency
- Time range & data source
- Available experiments (for ESM)
- Available variables with long names
- Usage notes & BASD tips

### 3. **Targeted Queries** (`find()`)
Search by one or more criteria:

```python
# Find all daily precipitation datasets
cd.find(variable="pr", frequency="daily")

# Find all ESM models
cd.find(type_filter="ESM")

# Find global datasets
cd.find(coverage="Global")
```

### 4. **Variable-Level Inspection** (`inspect()`)
Check units, BASD compatibility, and conversion hints:

```python
cd.inspect("ERA5", variable="tp")
```

Output:
```
────────────────────────────────────────────
  VARIABLE INSPECTION: tp  in  ERA5
────────────────────────────────────────────
  Long name  : Total Precipitation
  BASD unit  : kg m-2 s-1
  Conversion : ⚠  Auto-conversion available (m → kg m⁻² s⁻¹).
────────────────────────────────────────────
```

### 5. **ESM Models & Experiments Discovery** (`list_esm_experiments()`, `list_esm_models()`)
Explore available CMIP6 models and climate scenarios:

```python
# See all climate scenarios (experiments) in a dataset
cd.list_esm_experiments("CMIP")
cd.list_esm_experiments("NEXGDDP")

# List all ESM models available for an experiment
cd.list_esm_models("CMIP", experiment="historical")
cd.list_esm_models("NEXGDDP", experiment="ssp585")

# If not provided, defaults to first available experiment
cd.list_esm_models("NEXGDDP")
```

Output:
```
────────────────────────────────────────────
  EXPERIMENTS AVAILABLE IN: NEXGDDP
────────────────────────────────────────────
  1. historical
────────────────────────────────────────────

────────────────────────────────────────────
  ESM MODELS IN: NEXGDDP  (experiment: historical)
────────────────────────────────────────────
  1. ACCESS-CM2
  2. ACCESS-ESM1-5
  3. BCC-CSM2-MR
  ... (28 more models)
────────────────────────────────────────────
  Total: 31 models
```

## Object-Oriented Interface

Use `DatasetRegistry` class for programmatic access:

```python
from climdata.explore import DatasetRegistry

reg = DatasetRegistry()

# Method calls (same as functional API)
reg.list_available_data()
reg.explore("ERA5")
reg.find(variable="pr")

# Dict-like access
print(reg["ERA5"]["variables"])
print(reg["NEXGDDP"]["experiments"])

# Print summary
print(reg)
```

## Direct Dataset Class Access

For ESM datasets, you can also work directly with dataset classes for more control:

```python
from climdata.datasets.CMIPCloud import CMIPCloud
from climdata.datasets.NEXGDDP import NEXGDDP
from omegaconf import DictConfig

# Create minimal config
cfg = DictConfig({
    "experiment_id": "historical",
    "source_id": "dummy"  # Just to initialize
})

# Query CMIP6 catalog
cmip = CMIPCloud(cfg)
experiments = cmip.get_experiment_ids()  # ['historical', 'ssp126', 'ssp245', 'ssp370', 'ssp585']
models = cmip.get_source_ids("historical")  # [30+ models]
variables = cmip.get_variables(experiment_id="historical", source_id="GFDL-ESM4")

# Query NEX-GDDP
nexgddp = NEXGDDP(cfg)
models = nexgddp.get_source_ids("historical")  # [31 models]
```

## Architecture

### Directory Structure
```
climdata/explore/
├── __init__.py          # Public API exports
├── registry.py          # Dynamic registry loading from conf/
├── queries.py           # Query & display functions
├── catalog.py           # OOP interface (DatasetRegistry)
└── variables.py         # BASD conversion metadata
```

### How It Works

1. **Registry Loading** (`registry.py`)
   - Loads dataset definitions from `climdata/conf/mappings/parameters.yaml`
   - Adds static datasets (ERA5, W5E5, NEXGDDP, NASAPOWER)
   - Variables for each dataset come from config files

2. **Variable Metadata** (`variables.py`)
   - Only stores BASD conversion hints (units, notes)
   - Links CF variable names to their expected units in BASD
   - Flags variables needing pre-processing

3. **Query Functions** (`queries.py`)
   - `list_available_data()` — summary table
   - `explore(dataset)` — detailed dataset info
   - `find(**kwargs)` — filtering by criteria
   - `inspect(dataset, variable)` — unit & BASD checks

4. **OOP Wrapper** (`catalog.py`)
   - `DatasetRegistry` class wraps query functions
   - Supports dict-like access for power users
   - Pretty-prints summary via `__repr__`

## Dataset Inventory

| Abbr. | Name | Type | Coverage | Resolution |
|-------|------|------|----------|------------|
| ERA5 | ECMWF Reanalysis v5 | Observation | Global | 0.25° |
| MSWX | Multi-Source Weather | Observation | Global | 0.1° |
| W5E5 | WFDE5 + ERA5 (ISIMIP) | Observation | Global | 0.5° |
| DWD | German Weather Service | Station | Germany | point |
| HYRAS | DWD Gridded | Gridded | Germany | 1 km |
| CMIP | CMIP6 (Pangeo) | ESM (Raw) | Global | ~1° |
| CMIPW5E5 | CMIP6 (W5E5 format) | ESM (Bias-Corrected) | Global | 0.5° |
| NEXGDDP | NASA Earth Exchange | ESM (Downscaled) | Global | 0.25° |
| NASAPOWER | NASA POWER | Reanalysis-derived | Global | 0.5° |

## BASD Unit Hints

The `inspect()` function flags variables that need unit conversion before BASD:

| Variable | Units in Data | BASD Needs | Conversion Required? |
|----------|---------------|-----------|----|
| pr | kg m-2 s-1 | kg m-2 s-1 | ✓ No |
| tp | m | kg m-2 s-1 | ⚠ Yes |
| srad | MJ m-2 day-1 | W m-2 | ⚠ Yes (÷ 0.0864) |
| tas, tasmax, tasmin | K | K | ✓ No |
| hurs | % | % | ✓ No |

## Examples

See [examples/catalog_discovery_demo.py](../examples/catalog_discovery_demo.py) for a complete walkthrough.

## Data Source

- **Dataset metadata**: Loaded from `climdata/conf/mappings/parameters.yaml` (explore metadata + variables)
- **Variable BASD metadata**: Loaded from `climdata/conf/mappings/parameters.yaml` (basd_unit, basd_note per variable)
- **Experiments**: Extracted from `experiment_id` field in YAML config where available
- **Everything is configuration-driven** — No hardcoded dataset profiles or BASD metadata in Python code
