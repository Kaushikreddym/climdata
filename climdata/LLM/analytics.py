"""
ClimData Analytics Layer
========================

Phase 3b. Consumes EXTRACTED ClimData datasets (xarray) and emits structured,
machine-readable climate summaries (JSON dicts). No narrative text.

Capability policy — reuse ClimData first, implement only what is missing:

  ETCCDI indices / extremes  REUSE   -> climdata.utils.ClimateExtractor.calc_index
                                        (config-driven, 45 indices from Phase 1)
  climatology                IMPLEMENT (sdba.utils.compute_daily_climatology exists
                                        but runs argparse at import → not importable;
                                        re-implemented cleanly here)
  trends                     IMPLEMENT (not provided by ClimData)
  evaluation metrics         IMPLEMENT (only imputation RMSE/MAE exists; no model
                                        bias / evaluation — implemented here)

See ``capability_report()`` for the live provided-vs-implemented map.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional

import numpy as np
import xarray as xr

LOG = logging.getLogger("climdata.analytics")


# ===========================================================================
# Capability probe (reuse-first decision, derived from Phase 1 inventory)
# ===========================================================================

def capability_report() -> Dict[str, Dict]:
    """Map each analysis to whether it reuses ClimData or is implemented here.

    A reuse-first audit: extremes and ETCCDI indices delegate to ClimData's
    config-driven index engine, while trends, climatology and evaluation metrics
    have no ClimData equivalent and are implemented in this module.

    Returns:
        dict[str, dict]: One entry per analysis, with ``source``
        (``"climdata"`` or ``"analytics"``), ``mode`` (``"reuse"`` or
        ``"implement"``), and either an ``entrypoint`` or a ``reason``.
    """
    return {
        "etccdi_indices": {"source": "climdata", "mode": "reuse",
                           "entrypoint": "ClimateExtractor.calc_index"},
        "extremes": {"source": "climdata", "mode": "reuse",
                     "entrypoint": "ClimateExtractor.calc_index"},
        "climatology": {"source": "analytics", "mode": "implement",
                        "reason": "sdba.utils.compute_daily_climatology returns a smoothed "
                                  "day-of-year cycle; this layer needs plain monthly means"},
        "trends": {"source": "analytics", "mode": "implement",
                   "reason": "no trend function in climdata"},
        "evaluation_metrics": {"source": "analytics", "mode": "implement",
                               "reason": "only imputation RMSE/MAE exists; no model-evaluation metrics"},
    }


# ===========================================================================
# Helpers
# ===========================================================================

_SPATIAL_DIMS = ("lat", "lon", "latitude", "longitude", "x", "y")

# Variables the fixed-degree hot-day threshold is meaningful for.
_TEMPERATURE_VARIABLES = {"tas", "tasmax", "tasmin", "soilTemp"}
_TEMPERATURE_UNITS = {"degc", "°c", "c", "celsius", "degree_celsius",
                      "k", "kelvin", "degk"}


def _as_dataarray(data: xr.Dataset | xr.DataArray, variable: Optional[str]) -> xr.DataArray:
    if isinstance(data, xr.DataArray):
        return data
    if variable is None:
        if len(data.data_vars) != 1:
            raise ValueError(f"Specify `variable`; dataset has {list(data.data_vars)}.")
        variable = next(iter(data.data_vars))
    return data[variable]


def _reduce_space(da: xr.DataArray) -> xr.DataArray:
    """
    Average away any spatial dims, then any other leftover non-time dim, to a
    pure time series.

    The second pass matters beyond named lat/lon dims: some providers attach
    their own extra dimension during extraction (e.g. CMIPCloud's
    ``expand_dims("source_id")``), and a stray length-1 dim like that still
    makes every downstream stat 2-D — silently breaking the boolean indexing
    in ``_linear_slope`` with a numpy IndexError. Mean-reducing it is also the
    sane default for the rare case it is NOT length-1 (e.g. several ensemble
    members): an ensemble mean rather than a crash.
    """
    spatial = [d for d in da.dims if d in _SPATIAL_DIMS]
    da = da.mean(spatial) if spatial else da
    leftover = [d for d in da.dims if d != "time"]
    return da.mean(leftover) if leftover else da


def _annual_mean(da: xr.DataArray) -> xr.DataArray:
    return _reduce_space(da).groupby("time.year").mean()


def _linear_slope(years: np.ndarray, values: np.ndarray) -> float:
    """Least-squares slope per year, ignoring NaNs."""
    m = np.isfinite(values)
    if m.sum() < 2:
        return float("nan")
    return float(np.polyfit(years[m], values[m], 1)[0])


def is_temperature(data, variable: Optional[str] = None) -> bool:
    """Test whether a fixed degree-Celsius threshold means anything here.

    Guards the hot-day count in :func:`hot_days_change`, which would otherwise
    silently produce a number for precipitation or wind speed. Decided from the
    registry variable name when one is given, otherwise from the array's
    declared units.

    Args:
        data (xr.Dataset | xr.DataArray): Data whose units are inspected when
            ``variable`` is ``None``.
        variable (str, optional): Registry variable name, e.g. ``"tasmax"``.

    Returns:
        bool: ``True`` for a temperature variable or temperature units.
        ``False`` when neither can be established — the safe answer, since it
        suppresses the threshold statistic rather than fabricating one.
    """
    if variable is not None:
        return variable in _TEMPERATURE_VARIABLES
    try:
        units = str(_as_dataarray(data, None).attrs.get("units", "")).strip().lower()
    except (ValueError, KeyError):
        return False
    return units.replace(" ", "") in _TEMPERATURE_UNITS


def variable_units(variable: Optional[str]) -> Optional[str]:
    """Return the canonical display units for a variable, for the narrator to quote.

    Temperature is always reported in ``degC`` and precipitation in ``mm/day``
    regardless of a dataset's native units, because ClimData standardises to
    those on extraction. Fixing them here is what stops the narrator inventing
    or converting units of its own. Anything else falls back to the units
    declared in ``conf/mappings/variables.yaml``.

    Args:
        variable (str, optional): Registry variable name.

    Returns:
        str | None: The units string, or ``None`` if the variable is unknown or
        the registry cannot be read. Lookup failures are logged, never raised.
    """
    if variable in _TEMPERATURE_VARIABLES:
        return "degC"
    if variable == "pr":
        return "mm/day"
    try:
        from .capabilities import discover_variables
        for cap in discover_variables():
            if cap.name == variable:
                return cap.metadata.get("units")
    except Exception:  # noqa: BLE001 — units are supplementary, never fatal
        LOG.warning("could not resolve units for variable '%s'", variable)
    return None


# ===========================================================================
# TRENDS  (implemented — not in ClimData)
# ===========================================================================

def annual_trend(data, variable: Optional[str] = None) -> float:
    """Fit a linear trend to the annual-mean series.

    Any spatial dimensions are averaged away first, so a gridded input yields a
    single domain-mean trend.

    Args:
        data (xr.Dataset | xr.DataArray): Data with a ``time`` dimension.
        variable (str, optional): Variable to analyse. Required when ``data`` is
            a Dataset with more than one.

    Returns:
        float: Slope in units per year, rounded to four decimals. NaN if fewer
        than two years have finite values.
    """
    s = _annual_mean(_as_dataarray(data, variable))
    return round(_linear_slope(s["year"].values, s.values), 4)


def seasonal_trend(data, season: str = "JJA", variable: Optional[str] = None) -> float:
    """Fit a linear trend to one season's annual-mean series.

    Args:
        data (xr.Dataset | xr.DataArray): Data with a ``time`` dimension.
        season (str): Three-letter season code — ``"DJF"``, ``"MAM"``, ``"JJA"``
            or ``"SON"``. Defaults to ``"JJA"``.
        variable (str, optional): Variable to analyse.

    Returns:
        float: Slope in units per year, rounded to four decimals. NaN when the
        record holds no days in that season — logged, rather than failing on an
        empty grouping.
    """
    da = _reduce_space(_as_dataarray(data, variable))
    da = da.sel(time=da["time"].dt.season == season)
    if da["time"].size == 0:
        LOG.warning("seasonal_trend: no %s days in the record.", season)
        return float("nan")
    s = da.groupby("time.year").mean()
    return round(_linear_slope(s["year"].values, s.values), 4)


# ===========================================================================
# CLIMATOLOGY  (implemented — clean re-implementation)
# ===========================================================================

def climatology(data, variable: Optional[str] = None) -> Dict[str, float]:
    """Compute the long-term mean and the mean of each calendar month.

    Distinct from :func:`climdata.sdba.compute_daily_climatology`, which builds a
    Fourier-smoothed day-of-year cycle for bias adjustment. This is the coarser,
    directly quotable summary the narrator needs.

    Args:
        data (xr.Dataset | xr.DataArray): Data with a ``time`` dimension.
        variable (str, optional): Variable to analyse.

    Returns:
        dict: ``overall_mean`` (float) and ``monthly_mean`` (month number 1-12
        to mean), all rounded to four decimals.
    """
    da = _reduce_space(_as_dataarray(data, variable))
    monthly = da.groupby("time.month").mean()
    return {
        "overall_mean": round(float(da.mean()), 4),
        "monthly_mean": {int(m): round(float(v), 4)
                         for m, v in zip(monthly["month"].values, monthly.values)},
    }


def hot_days_change(data, thresh: float = 25.0, variable: Optional[str] = None,
                    window: int = 10) -> Optional[int]:
    """Measure how the annual count of days above a threshold has changed.

    Compares the mean count over the last ``window`` years against the first
    ``window`` years. A record too short for two non-overlapping windows falls
    back to a linear slope projected across its span, so a short series gets an
    estimate rather than nothing.

    Two guards keep the answer honest. Years with no valid observations are
    dropped rather than counted as zero exceedances, which would otherwise read
    as a spurious cooling; and a record with fewer than two usable years returns
    ``None`` rather than a fabricated ``0``.

    Only meaningful for temperature — check :func:`is_temperature` first.

    Args:
        data (xr.Dataset | xr.DataArray): Data with a ``time`` dimension.
        thresh (float): Threshold in the data's units. Defaults to ``25.0``.
        variable (str, optional): Variable to analyse.
        window (int): Years averaged at each end. Defaults to ``10``.

    Returns:
        int | None: Change in days per year, or ``None`` when the record is too
        short or too gappy to support one.
    """
    da = _reduce_space(_as_dataarray(data, variable))
    observed = da.notnull().groupby("time.year").sum()
    counts = (da > thresh).groupby("time.year").sum().where(observed > 0)

    yrs = counts["year"].values
    vals = counts.values.astype(float)
    finite = np.isfinite(vals)
    if finite.sum() < 2:
        LOG.warning("hot_days_change: only %d usable year(s) in the record.",
                    int(finite.sum()))
        return None

    if finite.sum() >= 2 * window:
        late, early = np.nanmean(vals[-window:]), np.nanmean(vals[:window])
        if np.isfinite(late) and np.isfinite(early):
            return int(round(late - early))
    # short series → slope projected over the full span
    span = float(yrs[finite][-1] - yrs[finite][0])
    slope = _linear_slope(yrs, vals)
    if not np.isfinite(slope) or span == 0:
        return None
    return int(round(slope * span))


# ===========================================================================
# EVALUATION METRICS  (implemented — not in ClimData)
# ===========================================================================

def evaluation_metrics(model, reference, variable: Optional[str] = None) -> Dict[str, float]:
    """Score a model series against a reference over the period they share.

    The two series are inner-joined on time first, so only overlapping steps are
    compared, and non-finite pairs are dropped. With no overlap at all every
    metric is NaN rather than an exception, so one unusable pairing degrades the
    assessment instead of ending the pipeline.

    Args:
        model (xr.Dataset | xr.DataArray): The simulation.
        reference (xr.Dataset | xr.DataArray): The observations.
        variable (str, optional): Variable to analyse.

    Returns:
        dict[str, float]: ``model_mean``, ``obs_mean``, ``model_bias``,
        ``rmse``, ``mae`` and ``correlation``, rounded to four decimals.
    """
    m = _reduce_space(_as_dataarray(model, variable))
    o = _reduce_space(_as_dataarray(reference, variable))
    m, o = xr.align(m, o, join="inner")
    mv, ov = m.values.astype(float), o.values.astype(float)
    mask = np.isfinite(mv) & np.isfinite(ov)
    mv, ov = mv[mask], ov[mask]
    if mv.size == 0:
        return {"model_mean": float("nan"), "obs_mean": float("nan"),
                "model_bias": float("nan"), "rmse": float("nan"),
                "mae": float("nan"), "correlation": float("nan")}
    diff = mv - ov
    return {
        "model_mean": round(float(mv.mean()), 4),
        "obs_mean": round(float(ov.mean()), 4),
        "model_bias": round(float(diff.mean()), 4),
        "rmse": round(float(np.sqrt((diff ** 2).mean())), 4),
        "mae": round(float(np.abs(diff).mean()), 4),
        "correlation": round(float(np.corrcoef(mv, ov)[0, 1]), 4) if mv.size > 1 else float("nan"),
    }


# ===========================================================================
# ETCCDI / EXTREMES  (REUSE ClimData)
# ===========================================================================

def compute_index(data, index: str, dataset: str = "MSWX",
                  reduce: str = "mean") -> Optional[float]:
    """Compute one ETCCDI or extreme index through ClimData's index engine.

    Reuses :meth:`~climdata.utils.wrapper_workflow.ClimateExtractor.calc_index`
    rather than reimplementing the index, so the definitions stay in
    ``conf/mappings/indices.yaml``.

    Every failure path returns ``None`` and logs why — an unimportable engine, a
    rejected dataset, an index the data cannot support. That keeps the pipeline
    fail-soft while leaving a trace, because a ``None`` here means the caller
    asked for something it will not get.

    Args:
        data (xr.Dataset | xr.DataArray): Input data.
        index (str): Index name, e.g. ``"tx90p"``, ``"rx1day"``.
        dataset (str): Provider name for the index engine's configuration.
            Defaults to ``"MSWX"``.
        reduce (str): ``"mean"`` for the average over the index series, or
            ``"change"`` for last minus first. Defaults to ``"mean"``.

    Returns:
        float | None: The scalar summary, or ``None`` if it could not be computed.
    """
    try:
        from climdata.utils.wrapper_workflow import ClimateExtractor
    except Exception as e:  # noqa: BLE001 — optional dependency boundary
        LOG.warning("index '%s': ClimData index engine unavailable (%s).", index, e)
        return None
    try:
        extractor = ClimateExtractor(overrides=[f"dataset={dataset}", f"index={index}"])
        idx_ds = extractor.calc_index(ds=data if isinstance(data, xr.Dataset) else data.to_dataset())
        if idx_ds is None:
            LOG.warning("index '%s': ClimData returned no result.", index)
            return None
        series = _reduce_space(idx_ds[index]) if index in idx_ds else _reduce_space(
            idx_ds[next(iter(idx_ds.data_vars))])
        vals = series.values.astype(float)
        if reduce == "change" and vals.size >= 2:
            return round(float(vals[-1] - vals[0]), 4)
        return round(float(np.nanmean(vals)), 4)
    except Exception as e:  # noqa: BLE001 — index engine boundary
        LOG.warning("index '%s' could not be computed: %s", index, e)
        return None


# ===========================================================================
# Orchestrator → structured JSON summary
# ===========================================================================

def analyze(
    data: xr.Dataset | xr.DataArray,
    intents: List[str],
    variable: Optional[str] = None,
    reference: Optional[xr.Dataset | xr.DataArray] = None,
    indices: Optional[List[str]] = None,
    hot_day_thresh: float = 25.0,
) -> Dict:
    """
    Run the requested analyses and return ONE flat, machine-readable summary.

    Keys are only present for the analyses actually requested/possible. Any
    requested statistic that could NOT be produced is recorded under the
    ``notes`` key — a list of plain-language caveats — so a missing figure is
    never silently indistinguishable from one that was never asked for. Notes
    carry no digits, so they cannot pollute the narrator's grounding check.
    """
    out: Dict[str, float] = {}
    notes: List[str] = []

    if "historical_trend" in intents:
        out["annual_trend"] = annual_trend(data, variable)
        out["summer_trend"] = seasonal_trend(data, "JJA", variable)

    if "climatology" in intents:
        clim = climatology(data, variable)
        out["overall_mean"] = clim["overall_mean"]
        out["monthly_mean"] = clim["monthly_mean"]

    if ("extremes" in intents) or ("etccdi_indices" in intents):
        for idx in (indices or []):
            val = compute_index(data, idx, reduce="change")
            if val is None:
                notes.append(f"the index '{idx}' could not be computed from the "
                             f"extracted data and is not reported")
            else:
                out[f"{idx}_change"] = val
        # threshold extreme available without ClimData — temperature only
        if is_temperature(data, variable):
            change = hot_days_change(data, hot_day_thresh, variable)
            if change is None:
                notes.append("the hot-day count needs at least two years of "
                             "observations and is not reported")
            else:
                out["hot_days_change"] = change
        else:
            notes.append(f"the fixed-degree hot-day threshold does not apply to "
                         f"'{variable or 'this variable'}' and was not evaluated")

    if "model_evaluation" in intents:
        if reference is None:
            notes.append("model evaluation was requested but no reference "
                         "dataset was available")
        else:
            out.update(evaluation_metrics(data, reference, variable))

    if notes:
        out["notes"] = notes
    return out


# ===========================================================================
# Example — synthetic demo producing an example-style summary (offline)
# ===========================================================================

if __name__ == "__main__":
    import json

    # Synthetic daily tasmax over Berlin: warming + seasonal cycle + noise.
    rng = np.random.default_rng(0)
    time = xr.cftime_range("1980-01-01", "2020-12-31", freq="D", calendar="standard")
    t = np.arange(len(time))
    years = np.array([d.year for d in time])
    season = 10 * np.sin(2 * np.pi * t / 365.25)
    warming = 0.04 * (years - years[0])          # ~0.04 degC/yr
    tasmax = 15 + season + warming + rng.normal(0, 2, len(time))
    da = xr.DataArray(tasmax, coords={"time": time}, dims="time", name="tasmax")
    da.attrs["units"] = "degC"

    # A "reference" that is 0.2 degC cooler on average → model_bias ≈ +0.2
    ref = (da - 0.2).rename("tasmax")

    summary = analyze(
        da, intents=["historical_trend", "extremes", "model_evaluation"],
        variable="tasmax", reference=ref,
    )
    print(json.dumps(summary, indent=2))
    print("\ncapability_report:", json.dumps(capability_report(), indent=2))
