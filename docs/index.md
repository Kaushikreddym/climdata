# climdata

[![PyPI](https://img.shields.io/pypi/v/climdata.svg)](https://pypi.python.org/pypi/climdata)
[![DOI](https://zenodo.org/badge/19554926/climdata.svg)](https://zenodo.org/record/19554926)

<p align="center">
  <img src="assets/climdata_logo.png" alt="climdata logo" class="page-logo" width="200">
</p>

**One interface for eleven climate data providers.** Extract from reanalyses,
observational grids, station networks and CMIP6 projections through a single
configuration, then compute extreme indices, regrid, bias-correct, fill gaps and
export — without rewriting the code for each source.

```python
import climdata as cd

extractor = cd.ClimData(overrides=[
    "dataset=hyras",                      # swap for mswx, cmip, nexgddp, ...
    "lat=52.5", "lon=13.4",
    "variables=[tasmin,tasmax,pr]",
    "time_range.start_date=2014-01-01",
    "time_range.end_date=2014-12-31",
])

ds = extractor.extract()                  # xarray.Dataset, CF names, unit-normalised
df = extractor.to_dataframe()             # long-form pandas
```

Change `dataset=` and nothing else. Variable names and units come out the same
whichever provider answered.

## Start here

<div class="grid cards" markdown>

- **[Install](installation.md)** — pip, conda, and the credentials three of the
  providers need.
- **[Usage](usage.md)** — the workflow, region selection, indices, regridding,
  imputation, output formats.
- **[Datasets](datasets.md)** — what each of the eleven providers offers, and
  its quirks.
- **[API reference](api.md)** — every public function and class.

</div>

## What it does

**Extraction** — point, bounding box, shapefile or GeoJSON, against
[eleven providers](datasets.md). Variables are renamed to CF conventions and
converted to declared units on the way out, so a `tasmax` from a DWD station and
a `tasmax` from a CMIP6 Zarr store are directly comparable.

**Extreme indices** — ETCCDI and threshold indices via xclim, declared in
configuration rather than code. climdata warns when a percentile-based index is
asked of a record too short to support it.

**Regridding** — reprojection and resampling onto a target CRS and cell size,
with a compatibility gate: a metric resolution combined with a geographic CRS is
[rejected rather than silently approximated](usage.md#regridding).

**Bias correction and downscaling** — the ISIMIP3BASD method, wrapped for
xarray. See the [BCSD guide](bcsd_guide.md).

**Gap filling** — imputation per grid cell, with several backends.

**Reproducible output** — long-form CSV, NetCDF, crop-model formats for SIMPLACE
and MONICA, or a FAIR RO-Crate with a JSON-LD manifest.

## Finding your way around the data

```python
import climdata as cd

cd.list_available_data()                     # the catalogue
cd.explore("NEXGDDP")                        # one provider in detail
cd.find(variable="pr", coverage="Germany")   # search across providers
```

None of these download anything.

## How configuration works

climdata is driven by [Hydra](https://hydra.cc). Defaults live in
`climdata/conf/`, and a working copy is placed in your current directory the
first time you run something, so you can edit it. Overrides passed to
`ClimData(overrides=[...])` are applied on top:

```python
"dataset=cmip"                    # replace a value
"variables=[tasmax,pr]"           # a list
"time_range.start_date=2050-01-01"   # a nested value
```

Providers themselves are declared in `conf/mappings/parameters.yaml`, so adding
one to the catalogue takes no code change.

## Citing

If climdata contributes to published work, please cite it via its
[Zenodo record](https://doi.org/10.5281/zenodo.19554926). `CITATION.cff` in the
repository has the machine-readable form.

## License

MIT — see the repository `LICENSE`. Note that individual **datasets** carry
their own terms; MSWX in particular is CC BY-NC 4.0 (non-commercial use only).
See the [MSWX guide](MSWX_guide.md).
