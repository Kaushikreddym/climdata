import climdata
import xarray as xr
import xclim
import pandas as pd
from climdata.utils.config import _ensure_local_conf
from omegaconf import DictConfig
import hydra
from hydra.core.global_hydra import GlobalHydra
import sys

# Example usage:
#
# Run the CLI with overrides:
#
#   python climdata_cli.py \
#       dataset=mswx \
#       lat=52.507 \
#       lon=13.137 \
#       time_range.start_date=2000-01-01 \
#       time_range.end_date=2000-12-31 \
#       dsinfo.mswx.params.google_service_account=/home/muduchuru/.climdata_conf/service.json \
#       data_dir=/beegfs/muduchuru/data/ \
#       variables=['tas']
#
# All Hydra overrides follow the format key=value.


## uncomment the below snippet for parallel processing 
# import dask
# from dask.distributed import Client

# # Configure Dask
# client = Client(
#     n_workers=20,        # or match number of physical cores
#     threads_per_worker=2,
#     memory_limit="10GB"  # per worker (8 * 10GB = 80GB total)
# )
# from multiprocessing import freeze_support

_ensure_local_conf()
@hydra.main(config_path="./conf", config_name="config", version_base=None)
def main(cfg: DictConfig) -> None:
    overrides = sys.argv[1:]

    # Extract data
    cfg, filename, ds = climdata.extract_data(overrides=overrides)

if __name__ == "__main__":
    main()
