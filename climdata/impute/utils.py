import numpy as np
import xarray as xr
from .._vendor.imputegap.recovery.manager import TimeSeries


def contaminate_ds_mcar(
    ds,
    *,
    variables=None,
    time_dim="time",
    rate_dataset=0.01,
    rate_series=0.01,
    seed=None,
    inplace=False,
):
    """
    Apply MCAR contamination on flattened spatial series (lat/lon),
    contaminating each series independently as shape (1, seq_len).
    """

    is_da = isinstance(ds, xr.DataArray)
    ds_out = ds if inplace else ds.copy(deep=True)

    targets = [None] if is_da else (variables or list(ds_out.data_vars))

    rng = np.random.default_rng(seed)

    for var in targets:
        arr = ds_out if is_da else ds_out[var]

        # all non-time dims become series
        series_dims = [d for d in arr.dims if d != time_dim]
        if not series_dims:
            raise ValueError("At least one non-time dimension is required")

        stacked = arr.stack(series=series_dims)
        stacked = stacked.transpose("series", time_dim)

        data = stacked.values  # (n_series, seq_len)
        out = data.copy()
        print(data.shape)
        for i in range(data.shape[0]):
            print(i)
            ts = data[i][None, :]  # shape (1, seq_len)

            ts_cont = TimeSeries.Contamination.mcar(
                ts,
                rate_dataset=rate_dataset,
                rate_series=rate_series,
                seed=rng.integers(0, 2**32 - 1),
            )

            # preserve original NaNs
            ts_cont = np.where(np.isnan(ts), np.nan, ts_cont)

            out[i] = ts_cont[0]

        da_cont = xr.DataArray(
            out,
            coords=stacked.coords,
            dims=stacked.dims,
        ).unstack("series")

        if is_da:
            ds_out = da_cont
        else:
            ds_out[var] = da_cont

    return ds_out
