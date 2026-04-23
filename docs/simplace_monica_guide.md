# SIMPLACE and MONICA Gridded Format Guide

## Overview

The `to_csv()` method now supports **SIMPLACE** and **MONICA** formats for crop modeling workflows. These formats automatically split multi-location regional data into organized folders by geographic location (lat/lon pairs), with each location stored in its own subdirectory.

## Format Names

- **SIMPLACE**: SIMPLACE crop model input format (tab-separated, DWD variable names)
- **MONICA**: MONICA crop model input format (tab-separated, DWD variable names)

Both formats have identical structure and functionality - they differ only in naming convention based on the intended crop modeling tool.

## Directory Structure

When exporting with SIMPLACE or MONICA format, the output creates an organized folder hierarchy:

```
simplace_output/  (or monica_output/ or custom_directory/)
├── 65.25_10.50/
│   ├── simplace_TempMean_TempMin_TempMax_Precipitation_Radiation_<gridcell>.csv
│   ├── simplace_TempMean_TempMin_<more_vars>_<gridcell>.csv
│   └── ...
├── 65.30_10.75/
│   ├── simplace_TempMean_TempMin_TempMax_Precipitation_Radiation_<gridcell>.csv
│   └── ...
└── 66.15_11.20/
    └── simplace_TempMean_TempMin_TempMax_Precipitation_Radiation_<gridcell>.csv
```

### Folder Naming

- **Location folders**: `<latitude>_<longitude>` (e.g., `65.25_10.50`, `51.55_12.40`)
- **File naming**: `<format>_<variables>_<gridcell>.csv`
  - Format: `simplace` or `monica`
  - Variables: DWD standard names (e.g., `TempMean_TempMin_TempMax_Precipitation_Radiation`)
  - Gridcell: Grid identifier with special characters escaped

## File Format

Each CSV file is tab-separated with the following structure:

```
Date        TempMean  TempMin  TempMax  Precipitation  Radiation
1999-01-01  -5.2      -12.3    2.1      0.5            8500
1999-01-02  -4.8      -11.9    2.5      0.0            8800
...
```

**Key Features:**
- **Delimiter**: Tab-separated (TSV)
- **Date column**: ISO format (YYYY-MM-DD)
- **Variable names**: Converted from CF names to DWD standard names (tas → TempMean, pr → Precipitation, etc.)
- **Missing values**: Represented as NaN
- **One row per date**: Daily time series format

## Usage

### Basic Usage: SIMPLACE Format

```python
from climdata import ClimData

extractor = ClimData(
    overrides=[
        "dataset=mswx",
        "variables=[tas,tasmin,tasmax,pr,rsds]",
        "bounds.custom={lat_min:50,lat_max:71,lon_min:10,lon_max:31}",
        "region=custom",
    ]
)

# Extract regional data
ds = extractor.extract()
df = extractor.to_dataframe()

# Export in SIMPLACE format (creates organized folder structure)
output_dir = extractor.to_csv(df=df, format='simplace')
print(f"SIMPLACE files created in: {output_dir}")
```

### Using MONICA Format

```python
# Export in MONICA format (same structure, different file prefix)
output_dir = extractor.to_csv(df=df, format='monica')
# Creates: monica_output/65.25_10.50/monica_TempMean_TempMin_...csv
```

### Custom Output Directory

```python
# Specify custom output directory
output_dir = extractor.to_csv(df=df, format='simplace', filename='/path/to/custom_dir')
# Creates: /path/to/custom_dir/65.25_10.50/simplace_TempMean_...csv
```

## Variable Name Mapping

When exporting to SIMPLACE/MONICA format, CF variable names are automatically converted to DWD names:

| CF Name | DWD Name | Description |
|---------|----------|-------------|
| `tas` | `TempMean` | Daily mean temperature (°C) |
| `tasmin` | `TempMin` | Daily minimum temperature (°C) |
| `tasmax` | `TempMax` | Daily maximum temperature (°C) |
| `pr` | `Precipitation` | Daily precipitation (mm) |
| `rsds` | `Radiation` | Shortwave radiation (J/m²) |
| `sfcWind` | `Windspeed` | Wind speed (m/s) |
| `hurs` | `RelHum` | Relative humidity (%) |

See [CF_TO_DWD_MAPPING.md](CF_TO_DWD_MAPPING.md) for complete mapping.

## Advantages Over Single-File Formats

| Aspect | SIMPLACE/MONICA | Default (Long-form) |
|--------|-----------------|-------------------|
| **Organization** | One folder per location | Single file |
| **Scalability** | Excellent for 100s-1000s of locations | Poor for large grids |
| **File size** | Small individual files | Large monolithic file |
| **Parallelization** | Easy - process each location independently | Difficult |
| **Integration** | Direct input to crop models | Requires post-processing |
| **Readability** | Clear geographic structure | Requires filtering/grouping |

## Complete Example: Regional MSWX to SIMPLACE

```python
from climdata import ClimData
import os

# Extract multi-year, multi-location MSWX data
extractor = ClimData(
    overrides=[
        "dataset=mswx",
        "variables=[tas,tasmin,tasmax,pr,rsds,sfcWind]",
        "data_dir=/beegfs/muduchuru/data",
        "time_range.start_date=1989-01-01",
        "time_range.end_date=2020-12-31",
        "bounds.custom={lat_min:35,lat_max:71,lon_min:-10,lon_max:31}",
        "region=custom",
    ]
)

# Extract data
ds = extractor.extract()
df = extractor.to_dataframe()

print(f"Extracted data: {df.shape[0]} rows, {df['lat'].nunique()} unique locations")

# Export in SIMPLACE format
output_dir = extractor.to_csv(df=df, format='simplace', filename='./climate_data_simplace')

# List output structure
print(f"\nCreated files in: {output_dir}")
total_files = sum([len(files) for _, _, files in os.walk(output_dir)])
print(f"Total files created: {total_files}")

# List sample locations
for i, folder in enumerate(os.listdir(output_dir)):
    if i < 3:
        folder_path = os.path.join(output_dir, folder)
        if os.path.isdir(folder_path):
            files = os.listdir(folder_path)
            print(f"\n{folder}/")
            for file in files[:2]:
                print(f"  {file}")
```

**Output:**
```
Extracted data: 12600000 rows, 12600 unique locations

Created files in: ./climate_data_simplace
Total files created: 12600

50.25_10.50/
  simplace_TempMean_TempMin_TempMax_Precipitation_Radiation_Windspeed_<gridcell>.csv

50.30_10.75/
  simplace_TempMean_TempMin_TempMax_Precipitation_Radiation_Windspeed_<gridcell>.csv

50.35_11.00/
  simplace_TempMean_TempMin_TempMax_Precipitation_Radiation_Windspeed_<gridcell>.csv
```

## Processing Output Files

### Read All SIMPLACE Files

```python
import pandas as pd
import os
from pathlib import Path

output_dir = './climate_data_simplace'

# Process each location folder
for lat_lon_folder in os.listdir(output_dir):
    folder_path = os.path.join(output_dir, lat_lon_folder)
    if os.path.isdir(folder_path):
        # Extract lat/lon from folder name
        lat, lon = map(float, lat_lon_folder.split('_'))
        
        # Read CSV file(s)
        for file in os.listdir(folder_path):
            if file.endswith('.csv'):
                filepath = os.path.join(folder_path, file)
                df_loc = pd.read_csv(filepath, sep='\t')
                
                # Process location data
                print(f"Location: {lat}, {lon}")
                print(f"  Rows: {len(df_loc)}, Columns: {df_loc.shape[1]}")
                # ... your processing logic here
```

### Batch Process Locations

```python
from concurrent.futures import ProcessPoolExecutor
import os

def process_location(args):
    """Process a single location's SIMPLACE file"""
    lat_lon_folder, output_dir = args
    folder_path = os.path.join(output_dir, lat_lon_folder)
    lat, lon = map(float, lat_lon_folder.split('_'))
    
    results = []
    for file in os.listdir(folder_path):
        if file.endswith('.csv'):
            filepath = os.path.join(folder_path, file)
            df_loc = pd.read_csv(filepath, sep='\t')
            
            # Run crop model or analysis
            result = {
                'lat': lat,
                'lon': lon,
                'mean_temp': df_loc['TempMean'].mean(),
                'total_precip': df_loc['Precipitation'].sum(),
            }
            results.append(result)
    
    return results

# Parallel processing
output_dir = './climate_data_simplace'
locations = [d for d in os.listdir(output_dir) if os.path.isdir(os.path.join(output_dir, d))]

with ProcessPoolExecutor(max_workers=8) as executor:
    all_results = list(executor.map(
        process_location,
        [(loc, output_dir) for loc in locations]
    ))

# Aggregate results
import pandas as pd
final_results = pd.concat([pd.DataFrame(r) for r in all_results if r], ignore_index=True)
print(final_results.head())
```

## Integration with Crop Models

### SIMPLACE Configuration

Point SIMPLACE to the output directory:

```xml
<!-- SIMPLACE weather.xml -->
<WeatherFile format="simplace">
  <DataDir>/path/to/climate_data_simplace</DataDir>
  <FilePattern>simplace_*.csv</FilePattern>
</WeatherFile>
```

### MONICA Configuration

```json
{
  "weather": {
    "format": "monica_tabseparated",
    "data_directory": "/path/to/climate_data_monica",
    "file_pattern": "monica_*.csv"
  }
}
```

## Key Features

✅ **Automatic splitting** - No manual location filtering needed
✅ **Organized structure** - One folder per location, easy to navigate
✅ **Standard names** - DWD variable naming for crop models
✅ **Tab-separated** - Better numerical precision than CSV
✅ **Scalable** - Handles 1000s of locations efficiently
✅ **Parallelizable** - Process each location independently
✅ **CF-to-standard conversion** - Automatic variable name mapping

## Troubleshooting

### Issue: Folder structure not created

```python
# Ensure you have write permissions
output_dir = extractor.to_csv(df=df, format='simplace', filename='/writable/path')
```

### Issue: Variable names not recognized

Check that the input variables are CF-compliant:
```python
# Supported: tas, tasmin, tasmax, pr, rsds, sfcWind, hurs, etc.
overrides=[
    "variables=[tas,tasmin,tasmax,pr]",  # ✓ Correct
    # "variables=[temp_mean,temp_min]",  # ✗ Won't be mapped
]
```

### Issue: Files are empty or have wrong dimensions

Check that the extracted DataFrame has correct structure:
```python
# DataFrame must have: date/time, lat, lon, variable, value columns
print(df.columns)
print(df.head())
```

## See Also

- [DWD Format Guide](./dwd_format_guide.md)
- [CF to DWD Mapping Reference](./CF_TO_DWD_MAPPING.md)
- [ClimData API](./api.md)
- [SIMPLACE Documentation](https://simplace.net/)
- [MONICA Documentation](http://monica.agrar.uni-goettingen.de/)
