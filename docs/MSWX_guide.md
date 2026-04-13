# MSWX Dataset Guide

## Overview

**MSWX** (Multi-Source Weather) is a global, daily, 0.1° meteorological dataset based on ERA5 reanalysis data, bias-corrected and validated against in-situ station observations. It provides consistent, high-quality weather data for climate and hydrological applications worldwide.

## Data Access & Licensing

### License: CC BY-NC 4.0 (Non-Commercial Use Only)

MSWX is released under the **Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)** license.

#### ✅ Permitted Uses
- **Academic research** at universities and research institutions
- **Non-profit scientific studies**
- **Personal projects** and non-commercial analysis
- **Government and NGO applications** for public benefit
- Educational purposes and teaching

#### ❌ Prohibited Uses
- **Commercial products or services** using MSWX data
- **For-profit companies** incorporating MSWX into commercial offerings
- **Commercial consulting services** based on MSWX data
- Any revenue-generating use without explicit permission

### How to Request Access

1. **Visit the MSWX website**: [https://www.gloh2o.org/mswx/](https://www.gloh2o.org/mswx/)
2. **Submit a data request** describing your intended use
3. **Confirm non-commercial use** in your application
4. **Receive approval** and download instructions from the MSWX team
5. **Configure climdata** with your Google Service Account credentials (see below)

> ⚠️ **Compliance**: By using MSWX via climdata, you agree to the CC BY-NC 4.0 license terms. Violations may result in legal action.

---

## Google Drive API Setup

MSWX data is hosted on Google Drive and accessed through the **Google Drive API** using a Service Account. This ensures secure, traceable access without sharing passwords.

### Step 1: Create a Google Cloud Project

1. Go to [Google Cloud Console](https://console.cloud.google.com/)
2. Click **"Select Project"** → **"New Project"**
3. Enter a project name (e.g., `climdata-mswx`)
4. Click **"Create"**

### Step 2: Enable Google Drive API

1. In the left sidebar, go to **APIs & Services** → **Library**
2. Search for **"Google Drive API"**
3. Click on it and click **"Enable"**

### Step 3: Create a Service Account

1. Go to **IAM & Admin** → **Service Accounts**
2. Click **"Create Service Account"**
3. Enter a name (e.g., `climdata-mswx-sa`)
4. Click **"Create and Continue"**
5. Click **"Done"** (you can skip role assignments for read-only Drive access)

### Step 4: Generate and Download JSON Key

1. Click on the Service Account email you just created
2. Go to the **"Keys"** tab
3. Click **"Add Key"** → **"Create new key"**
4. Select **"JSON"** and click **"Create"**
5. A JSON file will download automatically — **store it securely!**

### Step 5: Store the Key Securely

```bash
# Create a hidden directory for credentials
mkdir -p ~/.climdata_conf

# Move the downloaded JSON file
mv ~/Downloads/service-account-key.json ~/.climdata_conf/service.json

# Restrict file permissions (Linux/macOS)
chmod 600 ~/.climdata_conf/service.json

# Add to .gitignore to prevent accidental version control commits
echo "~/.climdata_conf/" >> .gitignore
```

### Step 6: Configure climdata

Update your climdata configuration to point to your service account JSON:

```python
from climdata import ClimData

overrides = [
    "dataset=mswx",
    "lat=52.5",
    "lon=13.4",
    "time_range.start_date=2010-01-01",
    "time_range.end_date=2020-12-31",
    "variables=[tasmin,tasmax,pr]",
    "data_dir=./data",
    "dsinfo.mswx.params.google_service_account=~/.climdata_conf/service.json",
]

extractor = ClimData(overrides=overrides)
ds = extractor.extract()
```

---

## Available Variables

| Variable | Description | Units | Long Name |
|----------|-------------|-------|-----------|
| `tasmin` | Daily minimum air temperature | °C | Minimum Temperature |
| `tasmax` | Daily maximum air temperature | °C | Maximum Temperature |
| `tas` | Daily mean air temperature | °C | Mean Temperature |
| `pr` | Daily precipitation | mm/day | Precipitation |
| `rsds` | Downward shortwave radiation | W/m² | Solar Radiation |
| `hurs` | Relative humidity | % | Relative Humidity |
| `sfcWind` | Wind speed | m/s | Surface Wind Speed |

---

## Spatial and Temporal Coverage

- **Spatial Resolution**: 0.1° × 0.1° (approximately 10 km at equator)
- **Temporal Resolution**: Daily
- **Global Coverage**: Worldwide
- **Time Period**: 1979–present (continuously updated)

---

## Data Quality & Methodology

MSWX combines ERA5 reanalysis with bias correction using in-situ observations:

1. **ERA5 as baseline**: High-resolution atmospheric reanalysis
2. **Station bias correction**: Adjusted using quality-controlled weather station data
3. **Validation**: Cross-validated against independent observations
4. **Uncertainty estimates**: Available for selected variables

For detailed methodology, see the MSWX publications on [https://www.gloh2o.org/mswx/](https://www.gloh2o.org/mswx/)

---

## Usage Examples

### Point Extraction

```python
from climdata import ClimData
from climdata.datasets.MSWX import MSWXmirror

overrides = [
    "dataset=mswx",
    "lat=52.5",  # Berlin
    "lon=13.4",
    "time_range.start_date=2020-01-01",
    "time_range.end_date=2020-12-31",
    "variables=[tasmin,tasmax,pr]",
    "data_dir=./data",
    "dsinfo.mswx.params.google_service_account=~/.climdata_conf/service.json",
]

extractor = ClimData(overrides=overrides)
mswx = MSWXmirror(extractor.cfg)
mswx.extract(point=(extractor.cfg.lon, extractor.cfg.lat))

ds = mswx.load(variable="tasmax")
print(ds)
```

### Bounding Box Extraction

```python
# Extract a region (Germany)
mswx.extract(box={
    "lon_min": 5.8,
    "lon_max": 15.0,
    "lat_min": 47.3,
    "lat_max": 55.1
})

ds = mswx.load(variable="pr")
print(ds)  # Returns spatiotemporal data for the region
```

### Workflow with Climate Indices

```python
overrides_wf = [
    "dataset=mswx",
    "lat=52.5",
    "lon=13.4",
    "time_range.start_date=2010-01-01",
    "time_range.end_date=2020-12-31",
    "variables=[tasmin,tasmax,pr]",
    "data_dir=./data",
    "dsinfo.mswx.params.google_service_account=~/.climdata_conf/service.json",
    "index=tn10p",      # Cold nights index
    "impute=BRITS",     # Optional: gap-filling
]

extractor = ClimData(overrides=overrides_wf)
result = extractor.run_workflow(actions=["extract", "impute", "calc_index", "to_nc"])
print(result)
```

---

## Common Issues

### Error: "No authentication provided" or "Google Drive API error"

**Solution**: Ensure:
1. Service account JSON key exists at the specified path
2. File path is correctly spelled in the config
3. Service account email has been shared with MSWX Google Drive folders

### Error: "Service Account not authorized to access folder"

**Solution**: 
1. Check your MSWX data request approval status
2. Confirm the service account email was shared with the correct folder
3. Contact MSWX support if access was granted but still failing

### Downloads are slow

**Solution**: 
1. Use `time_range` to limit data to only what you need
2. Consider requesting multiple points instead of a large bounding box
3. MSWX files are large; parallel downloads may be rate-limited by Google Drive

---

## Data Citation

When publishing results using MSWX data accessed via climdata, cite:

1. **MSWX Dataset**: Include the MSWX publication and DOI from [https://www.gloh2o.org/mswx/](https://www.gloh2o.org/mswx/)
2. **climdata Package**: Include version and GitHub repository
3. **License Compliance**: State "Data used under CC BY-NC 4.0 license"

Example:
> We used Multi-Source Weather (MSWX) data [citation] accessed via the climdata package [version X.X.X] under the CC BY-NC 4.0 license. Data were extracted for [location] during [time period].

---

## Related Resources

- **MSWX Homepage**: [https://www.gloh2o.org/mswx/](https://www.gloh2o.org/mswx/)
- **Climate Data Store (ERA5)**: [https://cds.climate.copernicus.eu/](https://cds.climate.copernicus.eu/)
- **climdata GitHub**: [https://github.com/Kaushikreddym/climdata](https://github.com/Kaushikreddym/climdata)

---

## Support

- **MSWX Issues**: Contact the MSWX team via [https://www.gloh2o.org/mswx/](https://www.gloh2o.org/mswx/)
- **climdata Issues**: Open an issue on [GitHub](https://github.com/Kaushikreddym/climdata/issues)
