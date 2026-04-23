# CF to DWD Variable Name Mapping Reference

## Overview

When exporting climate data to DWD format using `to_csv(format='dwd')`, all variable names are automatically converted from CF (Climate and Forecast) standard names to DWD (Deutscher Wetterdienst) standard names for consistency and compatibility.

## Complete Mapping Table

| CF Name | DWD Name | Description | Unit |
|---------|----------|-------------|------|
| `tas` | `TempMean` | Daily mean temperature | °C |
| `tasmin` | `TempMin` | Daily minimum temperature | °C |
| `tasmax` | `TempMax` | Daily maximum temperature | °C |
| `pr` | `Precipitation` | Daily precipitation | mm |
| `rsds` | `Radiation` | Shortwave downwelling radiation | J/m² |
| `rlds` | `LongwaveRadiation` | Longwave downwelling radiation | J/m² |
| `sfcWind` | `Windspeed` | Surface wind speed | m/s |
| `hurs` | `RelHum` | Relative humidity | % |
| `huss` | `SpecHum` | Specific humidity | kg/kg |
| `ps` | `SurfPressure` | Surface air pressure | Pa |

## Examples

### Example 1: Temperature Variables
**Input (CF names)**: `tas`, `tasmin`, `tasmax`
**Output (DWD names)**: `TempMean`, `TempMin`, `TempMax`

```python
from climdata import ClimData

extractor = ClimData(
    overrides=[
        "dataset=dwd",
        "variables=[tas,tasmin,tasmax]",  # CF names
        "lat=52.5",
        "lon=13.4",
    ]
)

ds = extractor.extract()
df = extractor.to_dataframe()

# Export with automatic name conversion
output_file = extractor.to_csv(format='dwd')
# Filename: dwd_1x1_TempMean_TempMin_TempMax_C_<col>_R_<row>.csv
# Columns in file: Date, TempMean, TempMin, TempMax
```

### Example 2: Multi-variable Extraction
**Input (CF)**: `tas`, `pr`, `rsds`, `sfcWind`
**Output (DWD)**: `TempMean`, `Precipitation`, `Radiation`, `Windspeed`

```python
extractor = ClimData(
    overrides=[
        "dataset=mswx",
        "variables=[tas,pr,rsds,sfcWind]",
        "bounds.custom={lat_min:50,lat_max:52,lon_min:10,lon_max:13}",
        "region=custom",
    ]
)

ds = extractor.extract()
df = extractor.to_dataframe()

# Filter to single station if needed
df_single = df[(df['lat'] == df['lat'].iloc[0]) & (df['lon'] == df['lon'].iloc[0])]

# Export with automatic name conversion
output_file = extractor.to_csv(df=df_single, format='dwd')
# Filename: dwd_1x1_TempMean_Precipitation_Radiation_Windspeed_C_<col>_R_<row>.csv
# Columns: Date, TempMean, Precipitation, Radiation, Windspeed, Gridcell
```

## Unmapped Variables

Variables not in the CF_TO_DWD_NAMES mapping dictionary are preserved as-is in the output. For example:
- `sfcWindmax` → `sfcWindmax` (no mapping, kept unchanged)
- Custom variables → kept unchanged

## Filename Examples

### With CF Variables (Converted to DWD)
```
# Input extraction
extractor = ClimData(overrides=[
    "dataset=dwd",
    "variables=[tas,tasmin,tasmax,pr,rsds]",  # CF names
    "lat=52.5",
    "lon=13.4",
])

# Output filename
dwd_1x1_TempMean_TempMin_TempMax_Precipitation_Radiation_C_180_R_1.csv
    ↑     ↑                        ↑                         ↑
    |     |                        |                         └─ Gridcell (converted)
    |     |                        └─ DWD variable names (converted)
    |     └─ 1 lat × 1 lon point
    └─ DWD format identifier
```

### Regional Extraction with DWD Names
```
dwd_3x5_TempMean_TempMin_TempMax_Precipitation_Radiation_C_100_R_50.csv
```

## Key Points

1. **Automatic**: No user action required - conversion happens automatically when `format='dwd'`
2. **Consistent**: All CF variables have corresponding DWD names
3. **Safe**: Unmapped variables pass through unchanged
4. **Reversible**: Reverse mapping available in code via `DWD_TO_CF_NAMES`
5. **Documented**: All mappings include descriptions and units

## Integration with SIMPLACE

SIMPLACE crop models typically expect DWD format with standard variable names:

```python
# Prepare climate input for SIMPLACE
extractor = ClimData(
    overrides=[
        "dataset=dwd",
        "variables=[tas,tasmin,tasmax,pr,sfcWind,rsds]",
        "lat=50.5",
        "lon=10.5",
        "time_range.start_date=1989-01-01",
        "time_range.end_date=2020-12-31",
    ]
)

ds = extractor.extract()
df = extractor.to_dataframe()

# Export in SIMPLACE-compatible DWD format
simplace_input = extractor.to_csv(format='dwd')
# File contains: Date, TempMean, TempMin, TempMax, Precipitation, Windspeed, Radiation
```

## See Also

- [DWD Format Guide](./docs/dwd_format_guide.md)
- [DWD Format Implementation](./DWD_FORMAT_IMPLEMENTATION.md)
- [CF Conventions](http://cfconventions.org/)
- [DWD Documentation](https://www.dwd.de/)
