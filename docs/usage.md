# Usage

How to drive climdata once it is [installed](installation.md). For what each
provider needs, see [Datasets](datasets.md); for the full signatures, the
[API reference](api.md).

## The idea

A climdata run is described by a **configuration**, not by code. `ClimData`
composes that configuration from Hydra defaults plus a list of overrides, then
carries a dataset through each step you ask for:

```python
import climdata as cd

extractor = cd.ClimData(overrides=[
    "dataset=hyras",
    "lat=52.5", "lon=13.4",
    "variables=[tasmin,tasmax,pr]",
    "time_range.start_date=2014-01-01",
    "time_range.end_date=2014-12-31",
])

ds = extractor.extract()
```

Each step stores its result on the extractor, so the next one needs no
arguments. `calc_index()` operates on what `extract()` produced; `to_csv()`
writes what `to_dataframe()` made.

```python
extractor.calc_index()      # uses the dataset extract() produced
df = extractor.to_dataframe()
extractor.to_csv(df)
```

Every stage is also addressable by name, should you want an earlier one back:
`extractor.ds`, `.index_ds`, `.impute_ds`, `.reprojected_ds`.

## Discovering what is available

Before extracting anything, ask what exists. None of this touches the network
except `list_esm_models`.

```python
cd.list_available_data()                      # the full catalogue
cd.explore("NEXGDDP")                         # one provider in detail
cd.find(variable="pr", coverage="Germany")    # search across providers
cd.inspect("HYRAS", variable="pr")            # units, and BASD conversion notes
cd.list_esm_experiments("CMIP")               # scenarios, from config
cd.list_esm_models("CMIP", experiment="ssp585")   # models, from the live catalogue
```

For programmatic access, `cd.explore.get_registry()` returns the same metadata
as a dictionary.

Within a configured extractor:

```python
extractor.get_datasets()              # providers this config knows
extractor.get_variables("HYRAS")      # variables that provider offers
extractor.get_varinfo("tasmax")       # CF name, long name, units
extractor.get_actions()               # workflow actions
extractor.get_indices(["tasmin"])     # indices computable from these variables
```

## Selecting a region

Point, bounding box, shapefile or GeoJSON — see
[Selecting a region](datasets.md#selecting-a-region) for all four with
examples. In brief:

```python
"lat=52.5", "lon=13.4"                                    # point
"box.lat_min=47", "box.lat_max=55", ...                   # bounding box
"shapefile=/path/to/region.shp"                           # polygon
f"aoi='{json.dumps(geojson)}'"                            # GeoJSON
```

## Chaining steps

`run_workflow` runs a named sequence and returns everything it produced:

```python
result = extractor.run_workflow(
    actions=["extract", "calc_index", "to_csv", "to_nc"]
)

print(result.keys())          # which fields were actually populated
print(result.dataframe.head())
print("Saved to:", result.filename)
```

The available actions are:

| Action | Does | Needs |
|---|---|---|
| `extract` | Fetch from the configured provider | `dataset` set |
| `upload_netcdf` | Read an existing `.nc` instead | `file=` argument |
| `upload_csv` | Read an existing long-form `.csv` | `file=` argument |
| `impute` | Fill gaps | a dataset, `impute=` set |
| `calc_index` | Compute the configured extreme index | a dataset, `index=` set |
| `reproject` | Regrid onto a target CRS/resolution | a dataset, a target set |
| `to_csv` | Convert to long form and write CSV | a dataset |
| `to_nc` | Write NetCDF | a dataset |
| `to_fair` | Write an RO-Crate | a dataset |

Steps that are configured as no-ops are skipped rather than failing:
`calc_index` with no `index=` returns `None` and leaves the dataset alone, as
does `reproject` with no target grid. That makes it safe to put them in a
standard sequence.

An action that cannot run says why:

```python
extractor.run_workflow(actions=["calc_index"])
# ValueError: Action 'calc_index' requires a dataset, but no dataset is
# available. Upload or extract a dataset before computing an index.
```

## Extreme indices

Indices are declared in `conf/mappings/indices.yaml` and computed by xclim.
Select one with `index=`:

```python
extractor = cd.ClimData(overrides=[..., "index=tn10p"])
extractor.extract()
index_ds = extractor.calc_index()
```

Percentile-based indices such as `tn10p` are defined against a 30-year
reference period. climdata warns when the record is shorter — the number it
returns is then a percentile of too little data, not a climate statistic.

To see which indices your variables support:

```python
extractor.get_indices(["tasmin", "tasmax"])              # needs all of them
extractor.get_indices(["tasmin"], require_all=False)     # needs any of them
```

## Regridding

`reproject` puts a dataset onto a target CRS and cell size:

```python
extractor = cd.ClimData(overrides=[
    ..., "target_projection=EPSG:3035", "target_resolution=10 km",
])
extractor.extract()
grid_ds = extractor.reproject()
```

The units of `target_resolution` must match the axis units of
`target_projection`. A metric resolution has no fixed angular size — 10 km spans
0.1395° of longitude but 0.0899° of latitude at 50°N — so combining `"10 km"`
with a geographic CRS raises `ResolutionCRSMismatch` rather than silently
approximating:

```python
cd.ClimData(overrides=[..., "target_projection=EPSG:4326",
                            "target_resolution=10 km"]).reproject()
# ResolutionCRSMismatch: ... Choose one:
#   - angular resolution   target_resolution="0.1 deg"
#   - projected CRS        target_projection="EPSG:3035"
```

If you want the approximation anyway, `cd.to_angular("10 km", latitude=52.5)`
converts it explicitly — exact at that latitude and only there.

See the [worked example](examples/grid/reprojection.ipynb).

## Filling gaps

```python
extractor = cd.ClimData(overrides=[..., "impute=SoftImpute"])
extractor.extract()
filled = extractor.impute()
```

Each grid cell is treated as an independent time series, so a cell with no valid
values of its own cannot be recovered from its neighbours. `SoftImpute`,
`CDRec` and `XGBOOST` work out of the box; `BRITS` needs PyTorch. Calling
`impute()` on data with no gaps returns it unchanged.

## Output

=== "Long-form CSV"

    ```python
    df = extractor.to_dataframe()
    extractor.to_csv(df, filename="out.csv")
    ```

    Columns: `time`, `lat`, `lon`, `variable`, `value`, `units`, `source`.
    Tab-separated.

=== "NetCDF"

    ```python
    extractor.to_nc(filename="out.nc")
    ```

=== "Crop-model formats"

    One file per grid cell, foldered by column index, with CF names mapped to
    the names SIMPLACE and MONICA expect (`tasmax` → `TempMax`):

    ```python
    extractor.to_csv(df, filename="./sim", format="simplace")
    extractor.to_csv(df, filename="./mon", format="monica")
    ```

=== "FAIR RO-Crate"

    A self-describing directory with the data, a JSON-LD manifest and an HTML
    preview:

    ```python
    extractor.to_fair(title="HYRAS Berlin 2014", creator="Your Name")
    ```

Output filenames are generated from the configuration and the current dataset,
so they follow the data through the workflow — after `calc_index()` the filename
names the index rather than the input variables. The templates live under
`output:` in `conf/config.yaml`.

## Working from existing files

```python
extractor.upload_netcdf("path/to/file.nc")
extractor.upload_csv("path/to/long_form.csv")
```

The CSV must have `time`, `lat`/`latitude`, `lon`/`longitude`, `variable` and
`value` columns; `units` and `source` are used if present. Both set the current
dataset, so the rest of the workflow follows as usual.

## Plotting

```python
extractor.plot(variable="tasmax")               # nearest time step
extractor.plot(variable="tasmax", reduce="mean")  # averaged over time
```

Draws on Cartopy axes with coastlines and borders; needs `matplotlib` and
`cartopy`. Point and station extractions have no grid to draw and raise a
message saying so.

## Parallelism

Extraction is single-threaded unless you opt in:

```python
extractor = cd.ClimData(overrides=[..., "load.dask.enabled=true"])
```

That starts a local Dask cluster and leaves the dataset lazy; without it, data
is materialised into memory as it loads. Close it with `extractor.close_dask()`.
An existing client in the process is reused rather than duplicated.
