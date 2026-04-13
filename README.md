# Welcome to climdata

[![DOI](https://zenodo.org/badge/19554926/climdata.svg)](https://zenodo.org/record/19554926)
[![image](https://img.shields.io/pypi/v/climdata.svg)](https://pypi.python.org/pypi/climdata)
[![image](https://img.shields.io/conda/vn/conda-forge/climdata.svg)](https://anaconda.org/conda-forge/climdata)

# ClimData — Quickstart & Overview

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

### Optional: Imputation Dependencies

If you need the **imputation functionality** (gap filling with ML models), install PyTorch and related packages:

```bash
# Install PyTorch CPU version (recommended for most users)
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu

# Install additional ML packages for imputation
pip install torch-cluster -f https://data.pyg.org/whl/torch-2.3.1+cpu.html
pip install pytorch-lightning torchmetrics lightning torchcde reformer-pytorch
pip install tensorflow darts sktime prophet tsfel tsfresh transformers timm
```

For GPU support, see [PyTorch installation guide](https://pytorch.org/get-started/locally/).

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

## Citation

If you use **climdata** in your research or projects, please cite it using the following formats:

### BibTeX
```bibtex
@software{muduchuru2024climdata,
  title={climdata: Automated Climate Data Extraction and Processing},
  author={Muduchuru, Kaushik},
  year={2024},
  version={0.5.0},
  url={https://github.com/Kaushikreddym/climdata},
  note={Available at https://Kaushikreddym.github.io/climdata}
}
```

### APA
Muduchuru, K. (2024). climdata: Automated climate data extraction and processing (v0.5.0). Retrieved from https://github.com/Kaushikreddym/climdata

### Citation.cff Format
Our repository includes a `CITATION.cff` file. GitHub will automatically show a "Cite this repository" button with ready-to-use citation formats.

### DOI (Zenodo)
[![DOI](https://zenodo.org/badge/19554926/climdata.svg)](https://zenodo.org/record/19554926)

**DOI**: https://doi.org/10.5281/zenodo.19554926  
**Zenodo Record**: https://zenodo.org/record/19554926

For archival setup details, see [Zenodo & DOI Guide](docs/zenodo_guide.md).

---

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
## ⚖️ Data Licensing & Access

### MSWX Dataset — Non-Commercial Use Only

**MSWX** (Multi-Source Weather) is released under the **CC BY-NC 4.0 license**. This means:

✅ **Allowed uses:**
- Academic research
- Non-profit scientific studies
- Personal projects
- Government or NGO applications (non-commercial)

❌ **Not allowed:**
- Commercial use or products
- For-profit services

**To access MSWX data:**
1. Visit [https://www.gloh2o.org/mswx/](https://www.gloh2o.org/mswx/)
2. Submit a data request for non-commercial use
3. Once approved, follow the Google Drive API setup below to configure climdata

> ⚠️ **Important:** By using MSWX via climdata, you agree to the CC BY-NC 4.0 license terms. Unauthorized commercial use is prohibited.

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