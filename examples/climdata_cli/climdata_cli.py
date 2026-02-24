from climdata import ClimData
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
#       data_dir=./data \
#       variables=[tas]
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

if __name__ == "__main__":
    # Get command line overrides
    overrides = sys.argv[1:]
    
    # Initialize ClimData with overrides
    extractor = ClimData(overrides=overrides)
    
    # Extract data (returns xarray.Dataset)
    ds = extractor.extract()
    
    # Optional: compute index if configured
    if extractor.cfg.get('index'):
        ds_index = extractor.calc_index(ds)
        df = extractor.to_dataframe(ds_index)
        extractor.to_csv(df)
    else:
        # Convert to dataframe and save
        df = extractor.to_dataframe(ds)
        extractor.to_csv(df)
