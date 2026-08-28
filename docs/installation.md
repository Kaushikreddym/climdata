# Installation

## Stable release

A conda environment is recommended, because several of climdata's geospatial
dependencies (`cartopy`, `rasterio`, `geopandas`) are considerably easier to
install that way:

```bash
conda create -n climdata python=3.11 -y
conda activate climdata

pip install climdata
```

## From source

```bash
git clone https://github.com/Kaushikreddym/climdata
cd climdata
pip install -e .
```

Or directly:

```bash
pip install git+https://github.com/Kaushikreddym/climdata
```

## Optional extras

Most providers work with the base install. These add capability:

| Extra | Command | Enables |
|---|---|---|
| Regridding (Dask) | `pip install climdata[grid]` | Parallel reprojection via `odc-geo` |
| Conservative regridding | `conda install -c conda-forge xesmf` | Area-weighted regridding, and BCSD |
| Plotting | `conda install -c conda-forge cartopy matplotlib` | `plot_map` and `extractor.plot` |
| CMIP6 cloud | `pip install intake intake-esm` | `dataset=cmip` |
| ISIMIP | `pip install isimip-client` | `w5e5`, `cmip_w5e5`, `agri_isimip` |
| DWD stations | `pip install wetterdienst` | `dataset=dwd` |
| ERA5 | `pip install cdsapi` | `climdata.ERA5` |
| Bias correction | `conda install -c conda-forge iris` | `climdata.sdba` |

### Imputation backends

Gap filling with the matrix-completion methods (`SoftImpute`, `CDRec`) works out
of the box. The neural methods need PyTorch, which climdata does not install by
default because it is large and platform-specific:

```bash
# CPU build, sufficient for most gap filling
pip install torch torchvision torchaudio \
    --index-url https://download.pytorch.org/whl/cpu

pip install torch-cluster -f https://data.pyg.org/whl/torch-2.3.1+cpu.html
pip install pytorch-lightning torchmetrics torchcde reformer-pytorch

conda install -c conda-forge pycatch22       # feature extraction
```

For GPU builds see the [PyTorch installation guide](https://pytorch.org/get-started/locally/).

`pip install climdata[full]` pulls the remaining ML backends (scikit-learn,
xgboost, lightgbm, optuna and friends) but deliberately not PyTorch or
`pycatch22` — install those as above.

## Credentials

Three providers need an account. The rest need nothing.

### MSWX — Google Drive service account

MSWX is distributed from a shared Google Drive folder, so downloads go through
the Drive API with a service-account key. The [MSWX guide](MSWX_guide.md) walks
through creating one, and covers the CC BY-NC 4.0 licence terms.

Once you have the key:

```bash
mkdir -p ~/.climdata
mv ~/Downloads/your-key.json ~/.climdata/service.json
chmod 600 ~/.climdata/service.json
```

```python
overrides = [
    "dataset=mswx", ...,
    "dsinfo.MSWX.params.google_service_account=~/.climdata/service.json",
]
```

The key is only required for data not already on disk.

### ERA5 — Copernicus CDS

1. Register at the [Copernicus Climate Data Store](https://cds.climate.copernicus.eu/user/register).
2. Open your [user profile](https://cds.climate.copernicus.eu/user) and reveal your API key.
3. Write `~/.cdsapirc`:

    ```
    url: https://cds.climate.copernicus.eu/api
    key: <your-api-key>
    ```

4. Restrict the permissions: `chmod 600 ~/.cdsapirc`

You must also accept each dataset's licence once, in the CDS web interface,
before the API will serve it.

### ISIMIP

No account is needed, only the client:

```bash
pip install isimip-client
```

## Checking the install

```python
import climdata as cd

print(cd.__version__)
cd.list_available_data()        # should print the catalogue
```

If that prints eleven providers, the configuration was found and you are ready.

## Development

```bash
pip install -e .
pip install -r requirements_dev.txt

pytest tests/                   # the offline suite
```

The test suite is offline by default. `tests/test_datasets_extraction.py` and
`tests/test_point_workflow.py` need real data paths and credentials and skip
themselves when `CI` is set.
