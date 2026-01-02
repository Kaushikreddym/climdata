# Welcome to climdata



[![image](https://img.shields.io/pypi/v/climdata.svg)](https://pypi.python.org/pypi/climdata)
<!-- [![image](https://img.shields.io/conda/vn/conda-forge/climdata.svg)](https://anaconda.org/conda-forge/climdata) -->

# ClimData — Quickstart & Overview
<p align="center">
  <img src="assets/climdata_logo.png" alt="ClimData Logo" class="page-logo" width="200">
</p>
ClimData provides a unified interface for extracting climate data from multiple providers (MSWX, CMIP, POWER, DWD, HYRAS), computing extreme indices, and converting results to tabular form. The ClimData (or ClimateExtractor) class is central: it manages configuration, extraction, index computation, and common I/O.

## Key features
- Provider-agnostic extraction (point / region / shapefile)
- Unit normalization via xclim
- Compute extreme indices using package indices
- Convert xarray Datasets → long-form pandas DataFrames
- Simple workflow runner for chained actions

## Installation

1) Create and activate a conda environment:
```bash
# create
conda create -n climdata python=3.11 -y

# activate
conda activate climdata
```

2) Install via pip (PyPI, if available) or from source:
```bash
# from PyPI
pip install climdata

# or from local source (editable)
git clone <repo-url>
cd climdata
pip install -e .
```

Install optional extras as needed (e.g., xclim, shapely, hydra, dask):
```bash
pip install xarray xclim shapely hydra-core dask "pandas>=1.5"
```

## Quick example
```python
from climdata import ClimData  # or from climdata.utils.wrapper_workflow import ClimateExtractor

overrides = [
    "dataset=mswx",
    "lat=52.5",
    "lon=13.4",
    "time_range.start_date=2014-01-01",
    "time_range.end_date=2014-12-31",
    "variables=[tasmin,tasmax,pr]",
    "data_dir=/path/to/data",
    "index=tn10p",
]

# initialize
extractor = ClimData(overrides=overrides)

# extract data (returns xarray.Dataset and updates internal state)
ds = extractor.extract()

# compute index (uses cfg.index)
ds_index = extractor.calc_index(ds)

# convert to long-form dataframe and save
df = extractor.to_dataframe(ds_index)
extractor.to_csv(df, filename="index.csv")
```

## Workflow runner
Use `run_workflow` for multi-step sequences:
```python
result = extractor.run_workflow(actions=["extract", "calc_index", "to_dataframe", "to_csv"])
```
`WorkflowResult` contains produced dataset(s), dataframe(s), and filenames.

## Documentation & API
- See API docs under `docs/api/` for detailed descriptions of ClimData/ClimateExtractor methods.
- Examples and notebooks are under `examples/`.

## Contributing
- Run tests and lint locally.
- Follow project coding and documentation conventions; submit PRs with tests.

## License
Refer to the repository LICENSE file for terms.

### ⚡️ Tip

- Make sure `yq` is installed:
  ```bash
  brew install yq   # macOS
  # OR
  pip install yq
  ```

- To see available variables for a specific dataset (for example `mswx`), run:
  ```bash
  python download_location.py --cfg job | yq '.mappings.mswx.variables | keys'
  ```

---

---

## ⚙️ **Key Features**

- **Supports multiple weather data providers**
- **Uses `xarray` for robust gridded data extraction**
- **Handles curvilinear and rectilinear grids**
- **Uses a Google Drive Service Account for secure downloads**
- **Easily reproducible runs using Hydra**

---
## 📡 Google Drive API Setup

This project uses the **Google Drive API** with a **Service Account** to securely download weather data files from a shared Google Drive folder.

Follow these steps to set it up correctly:

---

### ✅ 1. Create a Google Cloud Project

- Go to [Google Cloud Console](https://console.cloud.google.com/).
- Click **“Select Project”** → **“New Project”**.
- Enter a project name (e.g. `WeatherDataDownloader`).
- Click **“Create”**.

---

### ✅ 2. Enable the Google Drive API

- In the left sidebar, go to **APIs & Services → Library**.
- Search for **“Google Drive API”**.
- Click it, then click **“Enable”**.

---

### ✅ 3. Create a Service Account

- Go to **IAM & Admin → Service Accounts**.
- Click **“Create Service Account”**.
- Enter a name (e.g. `weather-downloader-sa`).
- Click **“Create and Continue”**. You can skip assigning roles for read-only Drive access.
- Click **“Done”** to finish.

---

### ✅ 4. Create and Download a JSON Key

- After creating the Service Account, click on its email address to open its details.
- Go to the **“Keys”** tab.
- Click **“Add Key” → “Create new key”** → choose **`JSON`** → click **“Create”**.
- A `.json` key file will download automatically. **Store it securely!**

### ✅ 5. Store the JSON Key Securely

- Place the downloaded `.json` key in the conf folder with the name service.json. 


## Setup Instructions from ERA5 api

### 1. CDS API Key Setup

1. Create a free account on the
[Copernicus Climate Data Store](https://cds.climate.copernicus.eu/user/register)
2. Once logged in, go to your [user profile](https://cds.climate.copernicus.eu/user)
3. Click on the "Show API key" button
4. Create the file `~/.cdsapirc` with the following content:

   ```bash
   url: https://cds.climate.copernicus.eu/api/v2
   key: <your-api-key-here>
   ```

5. Make sure the file has the correct permissions: `chmod 600 ~/.cdsapirc`

