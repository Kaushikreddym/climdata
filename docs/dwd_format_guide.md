# DWD Format Guide for `to_csv()`

## Overview

The DWD (Deutscher Wetterdienst) format is a standardized tabular format for single-station weather time series. When extracting climate data using the `ClimData` extractor, you can export results in DWD format using the `format='dwd'` parameter.

## DWD Format Characteristics

### File Structure
- **Format**: Tab-separated values (TSV)
- **Rows**: One row per date (time series)
- **Columns**: Date, followed by climate variables as columns
- **Location**: Single gridcell/station per file
- **Time Range**: Full extracted date range (e.g., 1951-01-01 to 2024-08-30)

### Example Columns
```
Date    Precipitation   TempMin  TempMean  TempMax  Radiation  SunshineDuration  Windspeed  RelHumCalc  Gridcell
1951-01-01  1.4  -7.5  -2.4  1.1  1313  0.0  2.7  69  C_180:R_1
1951-01-02  3.8  -1.3  0.1  1.2  1326  0.0  2.5  87  C_180:R_1
...
```

## Standard DWD Variables
| Variable | Unit | Description | CF Name |
|----------|------|-------------|---------|
| Precipitation | mm | Daily precipitation | `pr` |
| TempMin | °C | Daily minimum temperature | `tasmin` |
| TempMean | °C | Daily mean temperature | `tas` |
| TempMax | °C | Daily maximum temperature | `tasmax` |
| Radiation | J/m² | Solar radiation | `rsds` |
| SunshineDuration | hours | Hours of sunshine | - |
| Windspeed | m/s | Mean wind speed | `sfcWind` |
| RelHumCalc | % | Calculated relative humidity | `hurs` |
| Gridcell | - | Grid cell identifier (C_<col>:R_<row>) | - |

## Usage

### Basic Usage
```python
from climdata import ClimData

# Extract DWD data for a region
extractor = ClimData(
    overrides=[
        "dataset=dwd",
        "variables=[tas,tasmin,tasmax,pr,rsds]",
        "time_range.start_date=2020-01-01",
        "time_range.end_date=2020-12-31",
        "bounds.custom={lat_min:50,lat_max:52,lon_min:10,lon_max:13}",
        "region=custom",
    ]
)

# Extract and convert to DWD format
ds = extractor.extract()
df = extractor.to_dataframe()

# Save with DWD format
output_path = extractor.to_csv(format='dwd')
print(f"Saved to: {output_path}")
```

### Output Filename Format
```
dwd_<n_lat>x<n_lon>_<variable1>_<variable2>_..._<gridcell_id>.csv
```

**Example**: `dwd_1x1_TempMean_Precipitation_Radiation_C_180_R_1.csv`

Where:
- `n_lat`: Number of unique latitude points (1 for point extraction)
- `n_lon`: Number of unique longitude points (1 for point extraction)
- `variable1`, `variable2`, ...: Extracted variables in DWD names (e.g., TempMean, Precipitation, Radiation)
- `gridcell_id`: Grid cell identifier with `:` replaced by `_`

## Comparison: Default vs DWD Format

### Default Format (Long-form)
```
time     lat    lon    variable    value    ...
2020-01-01  50.5  11.2  tas  5.3
2020-01-01  50.5  11.2  pr   2.1
2020-01-02  50.5  11.2  tas  6.1
...
```
- **Use case**: Multi-location extraction, analysis across many stations
- **File size**: Larger (one row per variable per date per location)
- **Columns**: time, lat, lon, variable, value

### DWD Format (Tabular)
```
Date       Precipitation  TempMean  TempMin  TempMax  Radiation  Gridcell
2020-01-01  2.1           5.3       0.5      10.1     15000      C_180_R_1
2020-01-02  0.5           6.1       1.2      11.2     16500      C_180_R_1
```
- **Use case**: Single-station analysis, compatibility with DWD tools
- **File size**: Smaller (one row per date)
- **Columns**: Date, variable1 (DWD name), variable2 (DWD name), ...
- **Variable Names**: DWD standard names (TempMean, Precipitation, Radiation, etc.)

## Implementation Details

### Conversion Process
1. **Input**: Long-form DataFrame with columns `date`, `lat`, `lon`, `variable`, `value`
2. **Validation**: 
   - Checks for single location (multiple locations will raise error)
   - Verifies required columns exist
3. **Pivot**: Transforms wide format (Date x Variables)
4. **Metadata**: Extracts gridcell identifier and dimension counts
5. **Output**: Tab-separated CSV with DWD standard structure

### Error Handling
```python
# Will raise ValueError for multi-station data
try:
    extractor.to_csv(format='dwd')
except ValueError as e:
    print(f"Error: {e}")
    # "DWD format supports single station only, but found 2 locations"
```

## Examples

### Extract single-location DWD data
```python
from climdata import ClimData
import pandas as pd

# Point extraction
extractor = ClimData(
    overrides=[
        "dataset=dwd",
        "variables=[tas,tasmin,tasmax,pr]",
        "lat=52.5",
        "lon=13.4",
        "time_range.start_date=2010-01-01",
        "time_range.end_date=2020-12-31",
    ]
)

ds = extractor.extract()
df = extractor.to_dataframe()

# Export as DWD format (uses DWD variable names: TempMean, TempMin, TempMax, Precipitation)
output_file = extractor.to_csv(format='dwd')
print(f"Saved DWD time series to: {output_file}")
# Output filename: dwd_1x1_TempMean_TempMin_TempMax_Precipitation_C_<col>_R_<row>.csv

# Read back and inspect
df_dwd = pd.read_csv(output_file, sep='\t')
print(df_dwd.head())
print(f"Rows: {len(df_dwd)}, Columns: {len(df_dwd.columns)}")
```

### Export MSWX data as DWD format
```python
# Extract multi-location MSWX, then save single station as DWD
extractor = ClimData(
    overrides=[
        "dataset=mswx",
        "variables=[tas,pr]",
        "bounds.custom={lat_min:50,lat_max:52,lon_min:10,lon_max:13}",
        "region=custom",
    ]
)

ds = extractor.extract()
df = extractor.to_dataframe()

# Filter to single location
single_station = df[df['lon'] == df['lon'].iloc[0]][df['lat'] == df['lat'].iloc[0]]

# Save as DWD
output_file = extractor.to_csv(df=single_station, format='dwd')
```

## Integration with SIMPLACE

The DWD format is optimized for crop modeling workflows like SIMPLACE:

```python
from climdata import ClimData

# Prepare climate data for SIMPLACE in DWD format
extractor = ClimData(
    overrides=[
        "dataset=dwd",
        "variables=[tas,tasmin,tasmax,pr,sfcWind,rsds]",
        "bounds.custom={lat_min:47,lat_max:55,lon_min:5,lon_max:16}",
        "time_range.start_date=1989-01-01",
        "time_range.end_date=2020-12-31",
    ]
)

ds = extractor.extract()
df = extractor.to_dataframe()

# Export single station in DWD format for SIMPLACE
dwd_file = extractor.to_csv(format='dwd')

# File structure with DWD names: dwd_1x1_TempMean_TempMin_TempMax_Precipitation_Windspeed_Radiation_<gridcell>.csv
# Columns: Date, TempMean, TempMin, TempMax, Precipitation, Windspeed, Radiation, ...
# Ready for SIMPLACE model input
```

## Notes

- **Delimiter**: DWD format uses tab-separator (TSV) for better numerical precision
- **Missing Values**: Represented as NaN in output
- **Gridcell Format**: Format "C_<col>:R_<row>" is converted to "C_<col>_R_<row>" in filenames for compatibility
- **Single Location Only**: DWD format does not support multiple gridcells; point extraction recommended

## See Also

- [ClimData API Documentation](./api.md)
- [Usage Guide](./usage.md)
- [DWD Dataset Documentation](https://www.dwd.de/)
