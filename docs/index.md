# Welcome to climdata



[![image](https://img.shields.io/pypi/v/climdata.svg)](https://pypi.python.org/pypi/climdata)
[![image](https://img.shields.io/conda/vn/conda-forge/climdata.svg)](https://anaconda.org/conda-forge/climdata)


**This project automates the fetching and extraction of weather data from multiple sources — such as MSWX, DWD HYRAS, ERA5-Land, NASA-NEX-GDDP, and more — for a given location and time range.**


-   Free software: MIT License
-   Documentation: https://Kaushikreddym.github.io/climdata
    

## 📦 Data Sources

This project utilizes climate and weather datasets from a variety of data sources:

- **DWD Station Data**  
  Retrieved using the [DWD API](https://wetterdienst.readthedocs.io/en/latest/index.html). Provides high-resolution observational data from Germany's national meteorological service.

- **MSWX (Multi-Source Weather)**  
  Accessed via [GloH2O's Google Drive](https://www.gloh2o.org/mswx/). Combines multiple satellite and reanalysis datasets for global gridded weather variables.

- **DWD HYRAS**  
  Downloaded from the [DWD Open Data FTP Server](https://opendata.dwd.de/). Offers gridded observational data for Central Europe, useful for hydrological applications.

- **ERA5, ERA5-Land**  
  Accessed through the [Google Earth Engine](https://developers.google.com/earth-engine/datasets/catalog). Provides reanalysis datasets from ECMWF with high temporal and spatial resolution.

- **NASA NEX-GDDP**  
  Also retrieved via Earth Engine. Downscaled CMIP5/CMIP6 climate projections developed by NASA for local-scale impact assessment.

- **CMIP6**  
  Obtained using [ESGPull](https://esgf.github.io/esgf-download/) from the ESGF data nodes. Includes multi-model climate simulations following various future scenarios.

It supports:
✅ Automatic file download (e.g., from Google Drive or online servers)  
✅ Flexible configuration via `config.yaml`  
✅ Time series extraction for a user-specified latitude/longitude  
✅ Batch processing for many locations from a CSV file


## 🚀 How to Run and Explore Configurations

### ✅ Run a download job with custom overrides

You can run the data download script and override any configuration value directly in the command line using [Hydra](https://hydra.cc/).

For example, to download **ERA5-Land** data for **January 1–4, 2020**, run:

```bash
python download_location.py dataset='era5-land' \
  time_range.start_date='2020-01-01' \
  time_range.end_date='2020-01-04' \
  location.lat=52.5200 \
  location.lon=13.4050
```

For downloading multiple locations from a csv file `locations.csv`, run:

```bash
python download_csv.py dataset='era5-land' \
  time_range.start_date='2020-01-01' \
  time_range.end_date='2020-01-04' \
```

an example `locations.csv` can be

```csv
lat,lon,city
52.5200,13.4050,berlin
48.1351,11.5820,munich
53.5511,9.9937,hamburg
```

**What this does:**

- `dataset='era5-land'` tells the script which dataset to use.
- `time_range.start_date` and `time_range.end_date` override the default dates in your YAML config.
- All other settings use your existing `config.yaml` in the `conf` folder.

---

### ✅ List all available datasets defined in your configuration

To see what datasets are available (without running the downloader), you can dump the **resolved configuration** and filter it using [`yq`](https://github.com/mikefarah/yq).

Run:

```bash
python download_location.py --cfg job | yq '.mappings | keys'
```

**What this does:**

- `--cfg job` tells Hydra to output the final resolved configuration and exit.
- `| yq '.mappings | keys'` filters the output to show only the dataset names defined under the `mappings` section.

---

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


## Setup Instructions fro ERA5 api

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

