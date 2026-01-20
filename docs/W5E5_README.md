# W5E5 Dataset Support

The climdata package now includes support for the W5E5 (WFDE5 over land merged with ERA5 over ocean) global meteorological forcing dataset from ISIMIP.

## About W5E5

**W5E5** is a high-quality global meteorological forcing dataset that combines:
- WFDE5 (WATCH Forcing Data ERA5) over land
- ERA5 reanalysis over ocean
- Available through the ISIMIP (Inter-Sectoral Impact Model Intercomparison Project) data repository

### Dataset Characteristics

- **Spatial Resolution**: 0.5° × 0.5° (approximately 50 km)
- **Temporal Resolution**: Daily
- **Temporal Coverage**: 1979 - present
- **Global Coverage**: Complete global coverage
- **Data Source**: ISIMIP3a (observational/historical climate)

### Available Variables

| Variable | Name | Description | Units |
|----------|------|-------------|-------|
| `tas` | Near-Surface Air Temperature | Mean daily temperature at 2m height | K |
| `tasmax` | Maximum Temperature | Daily maximum temperature at 2m height | K |
| `tasmin` | Minimum Temperature | Daily minimum temperature at 2m height | K |
| `pr` | Precipitation | Total daily precipitation | kg m⁻² s⁻¹ |
| `rsds` | Shortwave Radiation | Surface downwelling shortwave radiation | W m⁻² |
| `rlds` | Longwave Radiation | Surface downwelling longwave radiation | W m⁻² |
| `hurs` | Relative Humidity | Near-surface relative humidity | % |
| `sfcWind` | Wind Speed | Near-surface wind speed | m s⁻¹ |
| `ps` | Surface Pressure | Surface air pressure | Pa |
| `huss` | Specific Humidity | Near-surface specific humidity | 1 |

## Installation

To use the W5E5 dataset, you need to install the `isimip-client` package:

```bash
pip install isimip-client
```

Or if installing climdata from the updated requirements:

```bash
pip install -e .
```

## Basic Usage

### 1. Using the Configuration System

```python
from climdata.utils.config import load_config
from climdata.datasets.W5E5 import W5E5

# Load configuration
cfg = load_config(
    config_name='config',
    overrides=[
        'dataset=w5e5',
        'lat=52.52',  # Berlin, Germany
        'lon=13.405',
        'variables=[tas,tasmax,tasmin,pr,rsds]',
        'time_range.start_date=2010-01-01',
        'time_range.end_date=2010-12-31',
        'data_dir=./data'
    ]
)

# Initialize and fetch data
w5e5 = W5E5(cfg)
w5e5.fetch()  # Download from ISIMIP
w5e5.load()   # Load into xarray
w5e5.extract(point=(cfg.lon, cfg.lat))  # Extract for location

# Save results
w5e5.save_csv('output.csv')
w5e5.save_netcdf('output.nc')
```

### 2. Direct Instantiation

```python
from omegaconf import OmegaConf
from climdata.datasets.W5E5 import W5E5

# Create configuration manually
cfg = OmegaConf.create({
    'dataset': 'w5e5',
    'lat': 40.7128,  # New York
    'lon': -74.0060,
    'variables': ['tas', 'pr'],
    'time_range': {
        'start_date': '2015-01-01',
        'end_date': '2015-12-31'
    },
    'data_dir': './w5e5_data'
})

w5e5 = W5E5(cfg)
w5e5.fetch()
w5e5.load()
w5e5.extract(point=(cfg.lon, cfg.lat))
```

## Advanced Usage

### Spatial Extraction Options

#### Point Extraction
```python
# Extract for a single point
w5e5.extract(point=(lon, lat))

# Extract with buffer (average over surrounding area)
w5e5.extract(point=(lon, lat), buffer_km=50)
```

#### Bounding Box Extraction
```python
# Extract for a rectangular region
w5e5.extract(box={
    'lon_min': 10.0,
    'lon_max': 15.0,
    'lat_min': 50.0,
    'lat_max': 55.0
})
```

#### Shapefile Extraction
```python
import geopandas as gpd

# Extract for a shapefile region
gdf = gpd.read_file('region.shp')
w5e5.extract(shapefile=gdf)
```

### Working with the xarray Dataset

```python
# Access the loaded dataset
ds = w5e5.ds

# Convert to pandas DataFrame
df = ds.to_dataframe()

# Unit conversions
df['tas_celsius'] = df['tas'] - 273.15  # K to °C
df['pr_mm_day'] = df['pr'] * 86400      # kg/m²/s to mm/day

# Basic statistics
print(df.describe())

# Plotting
import matplotlib.pyplot as plt
df['tas_celsius'].plot(title='Temperature Time Series')
plt.show()
```

## Data Access Details

### ISIMIP Repository Structure

W5E5 data is organized in the ISIMIP repository as:
```
ISIMIP3a/InputData/climate/atmosphere/obsclim/global/daily/historical/w5e5v2.0/
```

The W5E5 class automatically:
1. Searches for datasets matching your criteria
2. Downloads relevant files covering your time range
3. Caches files locally to avoid re-downloading
4. Loads and merges multiple files as needed

### File Naming Convention

W5E5 files follow this pattern:
```
w5e5v2.0_obsclim_{variable}_global_daily_{start_year}_{end_year}.nc
```

For example:
```
w5e5v2.0_obsclim_tas_global_daily_2010_2019.nc
```

## Performance Tips

1. **Download once**: Downloaded files are cached in `data_dir/w5e5/`
2. **Subset early**: Extract for your region immediately after loading to reduce memory usage
3. **Use appropriate time ranges**: W5E5 files cover multi-year periods, so choose ranges that minimize file downloads
4. **Parallel processing**: The `load()` method uses dask for parallel file reading

## Comparison with Other Datasets

| Feature | W5E5 | ERA5 | MSWX | NASA POWER |
|---------|------|------|------|------------|
| Resolution | 0.5° | 0.25° | 0.1° | 0.5° |
| Start Year | 1979 | 1950 | 1979 | 1981 |
| Coverage | Global | Global | Global | Global |
| Variables | 10+ | 100+ | 7 | 6+ |
| Update Freq | Annual | Monthly | Annual | Near real-time |
| Quality | High | Very High | High | Good |

### When to Use W5E5

- ✅ Need bias-corrected ERA5 data
- ✅ Working with ISIMIP climate impact models
- ✅ Require consistent land-ocean dataset
- ✅ 0.5° resolution is sufficient
- ✅ Post-1979 period is adequate

### When to Use Alternatives

- Use **ERA5** if you need higher resolution (0.25°) or more variables
- Use **MSWX** if you need higher resolution (0.1°) and focus on precipitation
- Use **NASA POWER** if you need near real-time data or solar energy applications

## Troubleshooting

### Installation Issues

```bash
# If isimip-client installation fails
pip install --upgrade pip
pip install isimip-client

# Or use conda
conda install -c conda-forge isimip-client
```

### Download Issues

```python
# Check ISIMIP API status
from isimip_client.client import ISIMIPClient
client = ISIMIPClient()
response = client.datasets(query='w5e5')
print(f"Found {len(response.get('results', []))} datasets")
```

### Memory Issues

```python
# Process data in chunks for large regions
w5e5.load()
# Extract smaller region first
w5e5.extract(box={'lon_min': 10, 'lon_max': 11, 
                   'lat_min': 50, 'lat_max': 51})
```

## References

1. **W5E5 Dataset**: Lange, S. (2019). WFDE5 over land merged with ERA5 over the ocean (W5E5). V. 1.0. DOI: [10.5880/pik.2019.023](https://doi.org/10.5880/pik.2019.023)

2. **ISIMIP Project**: [https://www.isimip.org/](https://www.isimip.org/)

3. **ISIMIP Data Repository**: [https://data.isimip.org/](https://data.isimip.org/)

4. **ISIMIP Client**: [https://github.com/ISI-MIP/isimip-client](https://github.com/ISI-MIP/isimip-client)

## Example Notebooks

See the following example notebooks:
- [`docs/examples/w5e5_example.ipynb`](examples/w5e5_example.ipynb) - Complete usage example

## Support

For issues related to:
- **W5E5 implementation**: Open an issue in the climdata repository
- **ISIMIP data access**: Visit [ISIMIP support](https://www.isimip.org/contact/)
- **Data quality/methodology**: See W5E5 documentation at [PIK](https://doi.org/10.5880/pik.2019.023)
