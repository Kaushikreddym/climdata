import numpy as np
import xarray as xr
from imputegap.recovery.manager import TimeSeries

def contaminate_ds_mcar(ds, *, variables=None, time_dim="time",
                        rate_dataset=0.01, rate_series=0.01,
                        seed=None, inplace=False):
    """
    Apply MCAR contamination to an xarray Dataset or DataArray.

    Parameters
    - ds: xarray.Dataset or xarray.DataArray
    - variables: list of variable names (for Dataset). If None, all data variables are used.
    - time_dim: name of the time dimension (default: "time")
    - rate_dataset, rate_series: forwarded to TimeSeries.Contamination.mcar
    - seed: optional int for reproducible contamination
    - inplace: if True modifies the input `ds` in-place and returns it; otherwise returns a modified copy

    Returns
    - contaminated xarray.Dataset or xarray.DataArray
    """
    is_da = isinstance(ds, xr.DataArray)
    ds_out = ds if inplace else ds.copy(deep=True)

    targets = [None] if is_da else (variables or list(ds_out.data_vars.keys()))

    # preserve RNG state if seed provided
    state = None
    if seed is not None:
        state = np.random.get_state()
        np.random.seed(seed)

    try:
        for var in targets:
            arr = ds_out if is_da else ds_out[var]
            # stack all non-time dims into a single 'series' dimension
            series_dims = [d for d in arr.dims if d != time_dim]
            if not series_dims:
                # single series -> shape (1, seq_len)
                arr2 = arr.values[np.newaxis, ...]
                coords = None
                stacked = None
            else:
                stacked = arr.stack(series=series_dims)
                arr2 = stacked.values  # shape (n_series, seq_len)
            # Apply MCAR contamination using imputegap TimeSeries helper
            arr_cont = TimeSeries.Contamination.mcar(arr2, rate_dataset=rate_dataset, rate_series=rate_series)
            # ensure existing NaNs are preserved (keep original NaNs)
            arr_cont = np.where(np.isnan(arr2), np.nan, arr_cont)
            # rebuild xarray object
            if series_dims:
                new_stack = xr.DataArray(arr_cont, coords=stacked.coords, dims=stacked.dims)
                new_arr = new_stack.unstack("series")
            else:
                # single series -> drop leading axis
                new_arr = xr.DataArray(arr_cont[0], coords=arr.coords, dims=arr.dims)
            # assign back
            if is_da:
                ds_out = new_arr
            else:
                ds_out[var] = new_arr
    finally:
        if state is not None:
            np.random.set_state(state)

    return ds_out