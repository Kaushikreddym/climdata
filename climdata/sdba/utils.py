"""Smoothing helpers shared by the bias-adjustment and downscaling routines.

The day-of-year climatology used throughout :mod:`climdata.sdba` is smoothed in
the frequency domain: keeping only the first few Fourier harmonics removes the
sampling noise that a short record leaves in a raw day-of-year mean, without the
phase shift a running window introduces.
"""

import numpy as np
import xarray as xr
from scipy.fft import fft, ifft

__all__ = ["compute_daily_climatology", "smooth_fft"]


def smooth_fft(x, n_harmonics: int = 3):
    """Smooth a periodic series by truncating its Fourier spectrum.

    All harmonics above ``n_harmonics`` are zeroed and the series transformed
    back, leaving the annual cycle and its first overtones. The input is treated
    as exactly one period, so it is meant for a full day-of-year axis (365 or 366
    values), not for an arbitrary time slice.

    Args:
        x (numpy.ndarray): One-dimensional periodic series, e.g. a day-of-year
            climatology.
        n_harmonics (int): Number of harmonics to retain. ``1`` keeps only the
            annual cycle; larger values follow the seasonal shape more closely.
            Defaults to ``3``.

    Returns:
        numpy.ndarray: Smoothed series, the same shape and dtype-kind as ``x``.

    Example:
        >>> import numpy as np
        >>> doy = np.arange(365)
        >>> raw = 10 * np.sin(2 * np.pi * doy / 365) + np.random.default_rng(0).normal(size=365)
        >>> smooth = smooth_fft(raw, n_harmonics=2)
        >>> smooth.std() < raw.std()
        np.True_
    """
    f = fft(x)
    f[n_harmonics + 1: -n_harmonics] = 0
    return np.real(ifft(f))


def compute_daily_climatology(data: xr.DataArray, n_harmonics: int = 3) -> xr.DataArray:
    """Build a smoothed day-of-year climatology.

    Groups ``data`` by day of year, averages over all years, then applies
    :func:`smooth_fft` along the ``dayofyear`` axis. The computation is lazy and
    Dask-friendly: chunked inputs stay chunked.

    Args:
        data (xr.DataArray): Input with a ``time`` dimension spanning at least one
            full year. Any additional dimensions (``lat``, ``lon``, …) are
            broadcast over.
        n_harmonics (int): Harmonics retained by :func:`smooth_fft`. Defaults to ``3``.

    Returns:
        xr.DataArray: The climatology, with ``time`` replaced by ``dayofyear``.

    Example:
        >>> clim = compute_daily_climatology(ds.tasmax)     # doctest: +SKIP
        >>> clim.dims                                       # doctest: +SKIP
        ('dayofyear', 'lat', 'lon')
    """
    doy_clim = data.groupby("time.dayofyear").mean("time")
    return xr.apply_ufunc(
        smooth_fft,
        doy_clim,
        kwargs={"n_harmonics": n_harmonics},
        input_core_dims=[["dayofyear"]],
        output_core_dims=[["dayofyear"]],
        vectorize=True,
        dask="parallelized",
        output_dtypes=[float],
        dask_gufunc_kwargs={"allow_rechunk": True},
    )
