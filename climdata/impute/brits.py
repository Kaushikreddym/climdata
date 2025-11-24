import numpy as np
import xarray as xr
from imputegap.wrapper.AlgoPython.BRITS.runnerBRITS import brits_recovery


def impute_brits_xarray(ds, vars, spatial_dims=("lat", "lon"),
                        model="brits", epoch=10, hidden_layers=64, seq_length=32,
                        batch_size=64, verbose=True):
    """
    Apply ImputeGAP BRITS to an xarray Dataset using variables as channels.
    
    Parameters
    ----------
    ds : xr.Dataset
        Dataset with dims (time, lat, lon) or (time, y, x)
    vars : list of str
        Variables to use as features/channels.
    spatial_dims : tuple
        The spatial dims to flatten.
    Returns
    -------
    xr.Dataset
        New Dataset with imputed values.
    """

    assert all(v in ds.data_vars for v in vars), "some variables missing"

    # -----------------------------
    # 1) Extract and stack variables
    # -----------------------------
    da = ds[vars].to_array("channel")  # shape: (channel, time, lat, lon)

    # move channel to last: (time, lat, lon, channels)
    if len(spatial_dims) == 2:
        da = da.transpose("time", spatial_dims[0], spatial_dims[1], "channel")
    elif len(spatial_dims) == 1:
        da = da.transpose("time", spatial_dims[0], "channel")
    else:
        # No spatial dims (just time, channel)
        da = da.transpose("time", "channel")

    # flatten spatial dims -> pixels
    flat = da.stack(pixel=spatial_dims)  # (time, pixel, channel)

    arr = flat.values  # numpy array

    T, P, C = arr.shape
    if verbose:
        print(f"BRITS input shape: time={T}, pixels={P}, channels={C}")

    # BRITS expects (samples, seq_length, nbr_features)
    # samples = pixels
    data_for_brits = np.transpose(arr, (1, 0, 2))  # (pixel, time, channel)

    # collapse channel into feature dimension (BRITS only supports 2D)
    data_for_brits = data_for_brits.reshape(P, T * C)

    # -----------------------------
    # 2) Run BRITS
    # -----------------------------
    recov = brits_recovery(
        incomp_data=data_for_brits,
        model=model,
        epoch=epoch,
        batch_size=batch_size,
        nbr_features=1,
        hidden_layers=hidden_layers,
        seq_length=seq_length,
        verbose=verbose
    )

    # reshape back
    recov = recov.reshape(P, T, C)
    recov = np.transpose(recov, (1, 0, 2))  # (time, pixel, channel)

    # -----------------------------
    # 3) Reconstruct xarray dataset
    # -----------------------------
    pixel_coord = np.arange(ds.sizes[spatial_dims[0]])

    imputed = xr.DataArray(
        recov,
        dims=("time", "pixel", "channel"),
        coords={
            "time": ds.time,
            "pixel": pixel_coord,
            "channel": vars,
        }
    )

    # unstack spatial dims
    imputed = imputed.unstack("pixel")  # back to (time, lat, lon, channel)

    # split channels back into variables
    out = imputed.to_dataset("channel")

    # rename channels to original var names
    out = out.rename({v: vars[i] for i, v in enumerate(out.data_vars)})

    return out
