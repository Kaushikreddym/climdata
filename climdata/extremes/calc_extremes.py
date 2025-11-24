import climdata
from climdata.extremes.indices import extreme_index 
import sys
from omegaconf import DictConfig
import hydra
from climdata.utils.config import _ensure_local_conf   
@hydra.main(config_path="conf", config_name="config", version_base="1.3")
def main(cfg: DictConfig):
    overrides = sys.argv[1:]

    # Extract data
    clim_cfg, clim_ds = climdata.extract_data(save_to_file=False)
        
    indices = extreme_index(clim_cfg, clim_ds)
    indices.calculate('tn10p')
    # indices.run()
if __name__ == "__main__":
    # from dask.distributed import Client, performance_report
    # client = Client(n_workers=20, threads_per_worker=2, memory_limit="80GB")  # ~80 logical cores
    import multiprocessing
    multiprocessing.freeze_support()
    main()