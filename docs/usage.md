# Usage

Quick examples to get started with the ClimData workflow utilities.

## Quickstart

Install into a conda env (recommended) and then pip:
```bash
conda create -n climdata python=3.11 -y
conda activate climdata
pip install climdata
# or from source:
# pip install -e .
```

## Minimal example
```python
from climdata import ClimData

overrides = [
    "dataset=mswx",
    "lat=52.5",
    "lon=13.4",
    "time_range.start_date=2014-01-01",
    "time_range.end_date=2014-12-31",
    "variables=[tasmin,tasmax,pr]",
]

extractor = ClimData(overrides=overrides)

# Extract dataset (updates extractor.current_ds)
ds = extractor.extract()

# Compute configured extreme index (updates extractor.index_ds)
index_ds = extractor.calc_index(ds)

# Convert to long-form DataFrame (updates extractor.current_df)
df = extractor.to_dataframe(index_ds)

# Save DataFrame to CSV
extractor.to_csv(df, filename="index.csv")
```

## Single-call workflow
Use the high-level runner to chain common steps:
```python
result = extractor.run_workflow(actions=["extract", "calc_index", "to_dataframe", "to_csv"])
# result contains produced dataset/dataframe and filenames
print(result.dataframe.head())
print("Saved to:", result.filename)
```

## Uploading existing files
- Load NetCDF: extractor.upload_netcdf("path/to/file.nc")
- Load long-form CSV: extractor.upload_csv("path/to/file.csv")

## Introspection helpers
- extractor.get_datasets()
- extractor.get_variables(dataset_name)
- extractor.get_varinfo(varname)
- extractor.get_actions()

## Notes
- See `docs/index.md` for installation details and full examples.
- For provider-specific options (MSWX, CMIP, POWER, DWD, HYRAS) consult the configuration files under `conf/` and the API docs.
