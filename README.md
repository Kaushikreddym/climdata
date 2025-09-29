# climdata


[![image](https://img.shields.io/pypi/v/climdata.svg)](https://pypi.python.org/pypi/climdata)
[![image](https://img.shields.io/conda/vn/conda-forge/climdata.svg)](https://anaconda.org/conda-forge/climdata)

`climdata` is a Python package designed to automate fetching, extraction, and processing of climate data from various sources, including MSWX, DWD HYRAS, ERA5-Land, and NASA-NEX-GDDP. It provides tools to retrieve data for specific locations and time ranges, facilitating climate analysis and research.

---

## 📄 Table of Contents

- [Installation](#installation)
- [Usage](#usage)
  - [Command-Line Interface (CLI)](#command-line-interface-cli)
  - [Python API](#python-api)
- [Configuration](#configuration)
- [Datasets](#datasets)
  - [MSWX](#mswx)
  - [DWD HYRAS](#dwd-hyras)
  - [ERA5-Land](#era5-land)
  - [NASA-NEX-GDDP](#nasa-nex-gddp)
- [Contributing](#contributing)
- [License](#license)

---

## 🛠️ Installation

Clone the repository and install dependencies:

```bash
git clone https://github.com/Kaushikreddym/climdata.git
cd climdata
pip install -r requirements.txt
````

---

## 🧭 Usage

### Command-Line Interface (CLI)

The `climdata_cli.py` script allows fetching and processing climate data via CLI:

```bash
python examples/climdata_cli.py dataset=MSWX lat=52.507 lon=13.137 time_range.start_date=2000-01-01 time_range.end_date=2000-12-31
```

This command retrieves data from the MSWX dataset for the specified latitude, longitude, and date range.

### Python API

You can also use `climdata` programmatically:

```python
from climdata import MSWX

mswx = MSWX(cfg)
mswx.fetch()
mswx.load()
df = mswx.extract(point=(13.137, 52.507))
mswx.save_csv("output.csv", df)
```

---

## ⚙️ Configuration

Configuration is managed via Hydra. A sample configuration (`config_mswx.yaml`) might look like:

```yaml
defaults:
  - _self_
  - mappings/parameters
  - mappings/variables

dataset: MSWX
lat: null
lon: null
variables: ["tasmin", "tasmax", "pr"]

data_dir: ./data
region: europe

experiment_id: historical
source_id: MIROC6
table_id: day

bounds:
  europe:
    lat_min: 34.0
    lat_max: 71.0
    lon_min: -25.0
    lon_max: 45.0

time_range:
  start_date: "1989-01-01"
  end_date: "2020-12-31"

output:
  out_dir: "./climdata/data/"
  filename_csv: "{provider}_{parameter}_LAT_{lat}_LON_{lon}_{start}_{end}.csv"
  filename_zarr: "{provider}_{parameter}_LAT{lat_range}_LON{lon_range}_{start}_{end}.zarr"
  filename_nc: "{provider}_{parameter}_LAT{lat_range}_LON{lon_range}_{start}_{end}.nc"
  fmt: "standard"
```

CLI overrides are automatically injected into the Hydra config.

---

## 📚 Datasets

### MSWX

The `MSWX` class interacts with the MSWX dataset. Key methods:

* `fetch()`: Download required data files.
* `load()`: Load data into an xarray Dataset.
* `extract()`: Extract data for a point or region.
* `save_csv()`: Save extracted data to CSV.

### DWD HYRAS

The `HYRASmirror` class manages DWD HYRAS data:

* `fetch(variable)`: Download NetCDF files for the variable and time range (only if not already present locally).
* `load(variable)`: Load NetCDF files into a single xarray Dataset.
* `extract(point, box, shapefile, buffer_km)`: Extract data at a point, bounding box, or shapefile area.
* `save_csv(filename, df)`: Save extracted data to CSV.

Supports **point extraction**, **box extraction**, and **shapefile-based clipping**.

### ERA5-Land

The `ERA5Land` class provides access to ERA5-Land reanalysis data:

* `fetch()`: Download specified variables and time ranges.
* `load()`: Load into xarray Dataset.
* `extract()`: Extract for location or region.
* `save_csv()`: Save extracted data to CSV.

### NASA-NEX-GDDP

The `NEXGDDP` class handles NASA-NEX GDDP climate projections:

* `fetch()`: Download necessary data files.
* `load()`: Load into xarray Dataset.
* `extract()`: Extract for point or region.
* `save_csv()`: Save to CSV.

---

## 🤝 Contributing

Contributions are welcome! Fork the repo, make changes, and submit a pull request. Ensure code style and tests match existing patterns.

---

## 📄 License

`climdata` is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

```

---

I can also create a **`docs/` folder structure** in Markdown for each `datasets/*.py` file with detailed class and method documentation, similar to a Sphinx-style doc, if you want.  

Do you want me to do that next?
```

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
