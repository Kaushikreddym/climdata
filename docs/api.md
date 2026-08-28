# API Reference

Generated from the source. For task-oriented documentation see
[Usage](usage.md); for provider specifics see [Datasets](datasets.md).

## Workflow

`ClimData` is the main entry point — an alias for `ClimateExtractor`. It
composes the configuration, drives each provider, and carries a dataset between
steps.

::: climdata.utils.wrapper_workflow.ClimateExtractor

::: climdata.utils.wrapper_workflow.WorkflowResult

::: climdata.utils.config.load_config

---

## Discovery

Browsing the catalogue without downloading anything.

::: climdata.explore.list_available_data

::: climdata.explore.explore

::: climdata.explore.find

::: climdata.explore.inspect

::: climdata.explore.list_esm_experiments

::: climdata.explore.list_esm_models

::: climdata.explore.get_registry

::: climdata.explore.resolve_dataset_key

::: climdata.explore.DatasetRegistry

---

## Providers

Each class can be used directly when the workflow wrapper is more than you need.
All follow a `fetch` → `load` → `extract` lifecycle, with the per-provider
deviations noted in each docstring.

### Observations, global

::: climdata.datasets.MSWX.MSWXmirror

::: climdata.datasets.ERA5.ERA5Mirror

::: climdata.datasets.W5E5.W5E5

::: climdata.datasets.NASAPOWER.POWER

### Observations, Germany

::: climdata.datasets.HYRAS.HYRASmirror

::: climdata.datasets.HOSTRADA.HOSTRADAmirror

::: climdata.datasets.DWD.DWDmirror

### Projections

::: climdata.datasets.CMIPCloud.CMIPCloud

::: climdata.datasets.NEXGDDP.NEXGDDP

::: climdata.datasets.CMIP_W5E5.CMIPW5E5

::: climdata.datasets.AGRI_ISIMIP.AGRI_ISIMIP

---

## Extreme indices

Indices are declared in `conf/mappings/indices.yaml` and evaluated through
xclim.

::: climdata.extremes.indices.extreme_index

---

## Grid transformation

Reprojection and resampling, plus the resolution/CRS compatibility gate. See the
[worked example](examples/grid/reprojection.ipynb).

::: climdata.grid.reproject.reproject

::: climdata.grid.reproject.resampling_from_name

### Resolution parsing

::: climdata.grid.units.parse_resolution

::: climdata.grid.units.resolution_in_crs_units

::: climdata.grid.units.to_angular

::: climdata.grid.units.Resolution

::: climdata.grid.units.ResolutionCRSMismatch

### CRS handling

::: climdata.grid.crs.parse_crs

::: climdata.grid.crs.crs_axis_unit

::: climdata.grid.crs.infer_spatial_dims

::: climdata.grid.crs.infer_src_crs

::: climdata.grid.crs.normalize_longitude

---

## Bias correction and downscaling

The ISIMIP3BASD method, wrapped for xarray. See the [BCSD guide](bcsd_guide.md).

::: climdata.sdba.BCSD

::: climdata.sdba.BiasCorrection

::: climdata.sdba.StatisticalDownscaling

::: climdata.sdba.regrid_to_coarse

::: climdata.sdba.utils.compute_daily_climatology

::: climdata.sdba.utils.smooth_fft

---

## Imputation

::: climdata.impute.impute_xarray.Imputer

---

## Visualisation

::: climdata.viz.plot_map
