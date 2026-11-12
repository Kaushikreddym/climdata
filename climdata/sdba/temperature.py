"""Keep bias-corrected daily temperatures physically ordered.

Bias-correcting ``tasmax`` and ``tasmin`` as two separate variables does not
keep ``tasmin <= tas <= tasmax``. MBCn corrects them jointly, but its last step
is a univariate quantile mapping per variable, and ISIMIP3BASD downscaling is
univariate by construction, so on the fine grid each of the three fields is
mapped on its own. In the CMIP6 -> MSWX run this left ``tasmin > tasmax`` on
2-3 % of cell-days and ``tas`` outside ``[tasmin, tasmax]`` on about 5 %.

The protocol in this module:

1. **Before correction**, replace ``tasmax`` with the diurnal temperature
   range ``dtr = tasmax - tasmin`` (:func:`add_dtr`). ``dtr`` has a physical
   lower bound of 0, so it is corrected as a lower-bounded variable (the
   ISIMIP3BASD ``tasrange`` settings: Rice distribution, mixed trend
   preservation, lower bound 0, threshold 0.01). The variables corrected are
   ``tas``, ``tasmin`` and ``dtr``. ``tasmax`` itself is never corrected.
2. **After correction and downscaling**, reconcile the three fields
   (:func:`enforce_temperature_order`):

   * ``dtr`` is raised to at least ``dtr_min``;
   * ``tas`` is kept as it is. It is the best-constrained of the three, and the
     one most downstream uses depend on;
   * ``tasmin`` is moved only where it falls outside the band that
     ``tas`` and ``dtr`` allow, ``tas - (1 - skew_min) * dtr <= tasmin <=
     tas - skew_min * dtr``. It is clipped to the nearer edge of that band;
   * ``tasmax`` is derived as ``tasmin + dtr``.

   The result always satisfies ``tasmin < tas < tasmax`` (strictly, for
   ``skew_min > 0`` and ``dtr_min > 0``), with ``tasmax - tasmin = dtr``.

   A cell-day that is already consistent keeps its ``tasmin`` exactly, so the
   correction only touches the days that need it.

The band is written in terms of the skewness ``s = (tas - tasmin) / dtr``: the
position of the daily mean within the daily range. ``s`` must lie in ``(0, 1)``
for ``tasmin < tas < tasmax``. Over Europe, MSWX 1979-2014 has ``s`` between
0.21 and 0.80 on 99.8 % of cell-days, so the default ``skew_min = 0.05`` only
binds on values the observations practically never show.
"""

from typing import Tuple

import numpy as np
import xarray as xr

__all__ = ["DTR_MIN", "SKEW_MIN", "add_dtr", "enforce_temperature_order"]

DTR_MIN = 0.1    # degC. MSWX's 0.1th percentile DTR over Europe is ~0.19 degC.
SKEW_MIN = 0.05  # (tas - tasmin) / dtr is kept in [SKEW_MIN, 1 - SKEW_MIN].

DTR_ATTRS = {
    "long_name": "Daily Near-Surface Air Temperature Range (tasmax - tasmin)",
    "units": "degC",
}


def add_dtr(ds: xr.Dataset) -> xr.Dataset:
    """Add ``dtr = tasmax - tasmin`` to a dataset holding both.

    Negative values (a raw model or an observation with ``tasmax < tasmin``)
    are set to 0, the variable's lower bound, so the bias adjustment's bound
    handling sees a valid series. ``tasmax`` is left in the dataset; callers
    choose which variables to correct.

    Args:
        ds (xr.Dataset): Dataset with ``tasmax`` and ``tasmin`` in the same units.

    Returns:
        xr.Dataset: A copy of ``ds`` with ``dtr`` added. Lazy inputs stay lazy.

    Raises:
        KeyError: If ``tasmax`` or ``tasmin`` is missing.
    """
    missing = {"tasmax", "tasmin"} - set(ds.data_vars)
    if missing:
        raise KeyError(f"add_dtr needs tasmax and tasmin; missing {sorted(missing)}")
    units = ds["tasmax"].attrs.get("units", DTR_ATTRS["units"])
    dtr = (ds["tasmax"] - ds["tasmin"]).clip(min=0)
    return ds.assign(dtr=dtr.assign_attrs({**DTR_ATTRS, "units": units}))


def enforce_temperature_order(
    tas,
    tasmin,
    dtr,
    dtr_min: float = DTR_MIN,
    skew_min: float = SKEW_MIN,
) -> Tuple[object, object, dict]:
    """Reconcile corrected ``tas``, ``tasmin`` and ``dtr``, and derive ``tasmax``.

    ``tas`` is not changed. ``dtr`` is raised to ``dtr_min``. ``tasmin`` is
    clipped into ``[tas - (1 - skew_min) * dtr, tas - skew_min * dtr]``, then
    ``tasmax = tasmin + dtr``. NaNs propagate.

    Works on NumPy arrays and on aligned xarray DataArrays.

    Args:
        tas: Corrected daily mean temperature.
        tasmin: Corrected daily minimum temperature.
        dtr: Corrected diurnal temperature range, same units.
        dtr_min (float): Smallest allowed range. Must be > 0 for a strict ordering.
        skew_min (float): Smallest allowed ``(tas - tasmin) / dtr`` and
            ``(tasmax - tas) / dtr``. Must be in ``[0, 0.5)``.

    Returns:
        tuple: ``(tasmin, tasmax, stats)``. ``stats`` counts the valid
        cell-days (``n_valid``), those whose ``dtr`` was raised (``n_dtr_raised``),
        and those whose ``tasmin`` was moved (``n_tasmin_adjusted``).

    Example:
        >>> import numpy as np
        >>> tmin, tmax, st = enforce_temperature_order(
        ...     np.array([10.0, 10.0]), np.array([5.0, 11.0]), np.array([8.0, 8.0]))
        >>> tmin.tolist(), tmax.tolist(), st["n_tasmin_adjusted"]
        ([5.0, 9.6], [13.0, 17.6], 1)
    """
    if not 0 <= skew_min < 0.5:
        raise ValueError(f"skew_min must be in [0, 0.5), got {skew_min}")
    if dtr_min < 0:
        raise ValueError(f"dtr_min must be >= 0, got {dtr_min}")

    dtr_new = np.maximum(dtr, dtr_min)
    lo = tas - (1.0 - skew_min) * dtr_new     # coldest tasmin allowed
    hi = tas - skew_min * dtr_new             # warmest tasmin allowed
    tasmin_new = np.minimum(np.maximum(tasmin, lo), hi)
    tasmax_new = tasmin_new + dtr_new

    valid = np.isfinite(tas) & np.isfinite(tasmin) & np.isfinite(dtr)
    stats = {
        "n_valid": int(np.sum(valid)),
        "n_dtr_raised": int(np.sum(valid & (dtr < dtr_min))),
        "n_tasmin_adjusted": int(np.sum(valid & (tasmin_new != tasmin))),
    }
    return tasmin_new, tasmax_new, stats
