# Datasets

Every provider is reached the same way: name it with `dataset=`, and change
nothing else. This page shows that pattern once, then lists what is actually
different per provider — the extra overrides, the credentials, and the quirks
worth knowing before you commit to a download.

## The pattern

```python
import climdata as cd

extractor = cd.ClimData(overrides=[
    "dataset=hyras",                      # <- the only line that changes
    "lat=52.5", "lon=13.4",
    "variables=[tasmin,tasmax,pr]",
    "time_range.start_date=2014-01-01",
    "time_range.end_date=2014-12-31",
    "data_dir=./data",
])

ds = extractor.extract()
df = extractor.to_dataframe()
```

Whatever the provider, `extract()` returns an `xarray.Dataset` whose variables
carry CF names and the units declared in `conf/mappings/variables.yaml` — climdata
converts them on the way out through xclim, so `tasmax` is degrees Celsius from
DWD stations and from a CMIP6 Zarr store alike.

!!! tip "Look before you download"
    `cd.explore("NEXGDDP")` prints a provider's coverage, resolution and
    variables without fetching anything, and `cd.find(variable="pr",
    coverage="Germany")` searches across providers. See
    [Discovering data](usage.md#discovering-what-is-available).

## What is available

| `dataset=` | Type | Coverage | Resolution | Period | Needs |
|---|---|---|---|---|---|
| `mswx` | Observation | Global | 0.1° daily | 1979– | Google service-account key |
| `era5` | Reanalysis | Global | 0.25° hourly/daily | 1940– | `~/.cdsapirc` key |
| `w5e5` | Observation | Global | 0.5° daily | 1979–2019 | `isimip-client` |
| `power` | Reanalysis-derived | Global, point only | daily | 1981– | nothing |
| `hyras` | Observation, gridded | Germany | 1 km daily | 1951– | nothing |
| `hostrada` | Observation, gridded | Germany | 1 km hourly | 1995– | nothing |
| `dwd` | Observation, stations | Germany | per station | varies | `wetterdienst` |
| `cmip` | ESM, raw | Global | model-dependent daily | 1850–2100 | `intake-esm` |
| `cmip_w5e5` | ESM, bias-corrected | Global | 0.5° daily | 1850–2100 | `isimip-client` |
| `nexgddp` | ESM, bias-corrected | Global | 0.25° daily | 1950–2100 | nothing |
| `agri_isimip` | Crop-model output | Global | 0.5° daily | varies | `isimip-client` |

### Choosing between the three CMIP6 options

They differ in how much processing has already been done for you:

- **`cmip`** is raw model output on each model's native grid, streamed from the
  Pangeo cloud catalogue. Use it when you want the model as published, and are
  prepared to bias-correct yourself — see [BCSD](bcsd_guide.md).
- **`nexgddp`** and **`cmip_w5e5`** are already bias-corrected and downscaled,
  to 0.25° and 0.5° respectively. Use these for impact work, where a raw model's
  bias would swamp the signal.
- Pair `cmip_w5e5` with `w5e5`, its own bias-adjustment reference, when you want
  a projection and the observations it was calibrated against on the same grid.

## Per-provider notes

### Observations, global

=== "MSWX"

    Distributed as one global NetCDF per day per variable on Google Drive, so a
    year of one variable is 365 files. Downloads are incremental — an
    interrupted run resumes rather than restarting.

    ```python
    overrides = [
        "dataset=mswx", "lat=52.5", "lon=13.4",
        "variables=[tasmax,pr]",
        "time_range.start_date=2014-01-01",
        "time_range.end_date=2014-12-31",
        "dsinfo.MSWX.params.google_service_account=/path/to/service.json",
    ]
    ```

    The key is only needed for data that is not already on disk. See the
    [MSWX guide](MSWX_guide.md) for obtaining one.

=== "ERA5"

    `climdata.ERA5` mirrors whole global months into Zarr stores rather than
    extracting a region, so it is driven directly rather than through
    `ClimData(...).extract()`:

    ```python
    import datetime
    from climdata import ERA5

    mirror = ERA5("./era5")
    paths = mirror.download(
        ["2m_temperature", ("temperature", 500)],
        (datetime.date(2020, 1, 1), datetime.date(2020, 3, 1)),
    )
    ds = xr.open_zarr(paths[0])
    ```

    Progress is journalled beside the stores, so an interrupted mirror resumes.
    Requires a Copernicus CDS account and a `~/.cdsapirc` key.

=== "W5E5"

    The ISIMIP reference dataset — WFDE5 over land merged with ERA5 over ocean.
    Its main use is as the observational baseline for `cmip_w5e5`.

    ```python
    overrides = [
        "dataset=w5e5", "lat=52.5", "lon=13.4",
        "variables=[tasmin,tasmax,pr]",
        "time_range.start_date=2004-01-01",
        "time_range.end_date=2004-12-31",
    ]
    ```

    Coverage stops at 2019.

=== "NASA POWER"

    The only genuinely point-based provider: the API returns one coordinate, so
    there is no grid to subset and box or shapefile requests do not apply. No
    credentials.

    ```python
    overrides = [
        "dataset=power", "lat=52.5", "lon=13.4",
        "variables=[tasmax,pr]",
        "time_range.start_date=2020-01-01",
        "time_range.end_date=2020-12-31",
    ]
    ```

    POWER encodes gaps as `-999`; screen for that sentinel before computing
    statistics.

### Observations, Germany

=== "HYRAS"

    1 km gridded observations, in two file layouts that climdata reconciles for
    you: yearly NetCDF for `pr`, `tas`, `tasmin`, `tasmax`, `hurs`; monthly
    `.tgz` archives of daily ASCII rasters for `evpot`, `soilMoist`, `soilTemp`.
    Archive members are read in place, never extracted to disk.

    ```python
    overrides = [
        "dataset=hyras",
        "box.lat_min=47", "box.lat_max=55",
        "box.lon_min=6", "box.lon_max=15",
        "variables=[tasmax,pr]",
        "time_range.start_date=2014-01-01",
        "time_range.end_date=2014-12-31",
    ]
    ```

    The grid is Gauss-Krüger zone 3 with two-dimensional `lat`/`lon`
    coordinates, so output has `x`/`y` dimensions. Pass it through
    [`reproject`](examples/grid/reprojection.ipynb) for a geographic grid.

=== "HOSTRADA"

    1 km *hourly* grids, one NetCDF per variable per month. The grid is Lambert
    Conformal Conic (EPSG:3034), again with 2-D `lat`/`lon`, so a box request is
    resolved to the enclosing `Y`/`X` rectangle — you get slightly more area
    than you asked for.

    ```python
    overrides = [
        "dataset=hostrada",
        "box.lat_min=47", "box.lat_max=55",
        "box.lon_min=6", "box.lon_max=15",
        "variables=[rsds]",
        "time_range.start_date=2020-06-01",
        "time_range.end_date=2020-06-30",
    ]
    ```

=== "DWD"

    Station measurements, not a grid. Every spatial selection picks a *set of
    stations*: a point request finds the stations near it, a box request the
    stations inside it.

    ```python
    overrides = [
        "dataset=dwd", "lat=52.5", "lon=13.4",
        "variables=[tasmax,pr]",
        "time_range.start_date=2014-01-01",
        "time_range.end_date=2014-12-31",
    ]
    ```

    The result is indexed by `station_id`, so it has no regular `lat`/`lon` axes
    and cannot be plotted as a map or passed to `reproject` without being
    gridded first. A small search radius in a sparsely instrumented area finds
    no stations and raises — widen it rather than assuming the data is absent.
    See the [DWD format guide](dwd_format_guide.md) for the tabular export.

### Projections

=== "CMIP6 (raw)"

    Streamed from the Pangeo cloud catalogue as Zarr, so nothing is downloaded
    until the data is computed.

    ```python
    overrides = [
        "dataset=cmip", "lat=52.5", "lon=13.4",
        "variables=[tasmax,pr]",
        "experiment_id=ssp585", "source_id=MIROC6", "table_id=day",
        "time_range.start_date=2050-01-01",
        "time_range.end_date=2050-12-31",
    ]
    ```

    Two things to know. Historical runs cover 1850–2014 and SSP scenarios
    2015–2100; asking an SSP for 1990 is rejected at construction rather than
    after a catalogue query. And not every model publishes every variable at
    every frequency, so filter the model list before choosing:

    ```python
    from climdata import CMIP
    CMIP.get_source_ids("ssp585", table_id="day", variables=["tasmax", "pr"])
    ```

=== "NEX-GDDP-CMIP6"

    Downscaled to 0.25°, served as one NetCDF per variable per year from NASA's
    THREDDS server. Model names are uppercase.

    ```python
    overrides = [
        "dataset=nexgddp", "lat=52.5", "lon=13.4",
        "variables=[tasmax,pr]",
        "experiment_id=ssp585", "source_id=GFDL-ESM4",
        "time_range.start_date=2050-01-01",
        "time_range.end_date=2050-12-31",
    ]
    ```

    There is no catalogue API, so the realisation and grid label are discovered
    by probing URLs. A global annual file is ~1440×600 cells; the region of
    interest is applied per file as each is opened, which is why a country-sized
    request is tractable.

=== "CMIP6 bias-adjusted to W5E5"

    ISIMIP3b output on the W5E5 0.5° grid. Model names are **lowercase** here,
    unlike the CMIP6 catalogue.

    ```python
    overrides = [
        "dataset=cmip_w5e5", "lat=52.5", "lon=13.4",
        "variables=[tasmax,pr]",
        "experiment_id=ssp370", "source_id=gfdl-esm4",
        "time_range.start_date=2080-01-01",
        "time_range.end_date=2080-12-31",
    ]
    ```

    See the [worked example](examples/datasets/cmip_w5e5_example.ipynb) for
    looping over scenarios.

=== "Crop model output"

    `agri_isimip` serves the *output* of crop models — yields, biomass, planting
    dates — rather than climate variables. A file is identified by five
    coordinates (impact model, climate forcing, period, crop, irrigation
    regime), and which combinations exist cannot be guessed, so this is the one
    provider with an explicit discovery step:

    ```python
    from climdata import AGRI_ISIMIP

    agri = AGRI_ISIMIP(cfg)
    agri.explore()                    # what exists for the configured queries
    print(sorted(agri.AVAILABLE_CROPS))
    ```

    The `AVAILABLE_*` properties raise until `explore()` or
    `query_available_metadata()` has run.

## Selecting a region

The same three modes work for every gridded provider.

=== "Point"

    ```python
    overrides = ["dataset=hyras", "lat=52.5", "lon=13.4", ...]
    ```

    Selects the nearest grid cell. For DWD, the nearest station.

=== "Bounding box"

    ```python
    overrides = [
        "dataset=hyras",
        "box.lat_min=47", "box.lat_max=55",
        "box.lon_min=6", "box.lon_max=15",
        ...,
    ]
    ```

=== "Shapefile"

    ```python
    overrides = ["dataset=hyras", "shapefile=/path/to/region.shp", ...]
    ```

    Clipped to the polygon for most providers; HOSTRADA uses the bounding box
    of the geometry instead.

=== "GeoJSON"

    A GeoJSON `Point` sets `lat`/`lon`; a `Polygon` sets the bounds. Useful when
    the region comes from a web map or another tool.

    ```python
    import json

    aoi = {"type": "Point", "coordinates": [13.4, 52.5]}   # lon, lat
    overrides = ["dataset=hyras", f"aoi='{json.dumps(aoi)}'", ...]
    ```

## Adding a provider

Providers are discovered from `conf/mappings/parameters.yaml`, not listed in
code, so a new entry appears in `list_available_data()`, `explore()` and the
capability registry with no code change. Extraction itself does need a branch in
`ClimateExtractor.extract()`; if it is missing, the error names the providers
that do have one.
