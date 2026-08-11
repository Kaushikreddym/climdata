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
    """What ClimData already provides vs. what this layer implements."""
    return {
        "etccdi_indices": {"source": "climdata", "mode": "reuse",
                           "entrypoint": "ClimateExtractor.calc_index"},
        "extremes": {"source": "climdata", "mode": "reuse",
                     "entrypoint": "ClimateExtractor.calc_index"},
        "climatology": {"source": "analytics", "mode": "implement",
                        "reason": "sdba.utils.compute_daily_climatology not importable (argparse side-effect)"},
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
    """
    Whether a fixed degree-Celsius threshold is meaningful for these data.

    Decided from the registry variable key when one is given, otherwise from
    the array's declared units.
    """
    if variable is not None:
        return variable in _TEMPERATURE_VARIABLES
    try:
        units = str(_as_dataarray(data, None).attrs.get("units", "")).strip().lower()
    except (ValueError, KeyError):
        return False
    return units.replace(" ", "") in _TEMPERATURE_UNITS


def variable_units(variable: Optional[str]) -> Optional[str]:
    """
    Canonical display units for a variable, for the narrator to quote.

    Temperature variables are always reported in degC and precipitation
    (``pr``) is always mm/day, regardless of a dataset's native units —
    ClimData standardises to these on extraction, and the narrator must never
    invent or convert units on its own. Other variables fall back to the
    registry's declared units (``climdata/conf/mappings/variables.yaml``).
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
    """Linear trend of the annual-mean series, in <units> per year."""
    s = _annual_mean(_as_dataarray(data, variable))
    return round(_linear_slope(s["year"].values, s.values), 4)


def seasonal_trend(data, season: str = "JJA", variable: Optional[str] = None) -> float:
    """
    Linear trend of a single season's annual-mean series, in <units> per year.

    Returns NaN when the record contains no data for the season, rather than
    failing on an empty grouping.
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
    """Long-term mean and monthly climatology (machine-readable)."""
    da = _reduce_space(_as_dataarray(data, variable))
    monthly = da.groupby("time.month").mean()
    return {
        "overall_mean": round(float(da.mean()), 4),
        "monthly_mean": {int(m): round(float(v), 4)
                         for m, v in zip(monthly["month"].values, monthly.values)},
    }


def hot_days_change(data, thresh: float = 25.0, variable: Optional[str] = None,
                    window: int = 10) -> Optional[int]:
    """
    Change in the annual count of days above `thresh`, comparing the mean of the
    last `window` years to the first `window` years. Threshold-based extreme;
    falls back to ClimData ETCCDI for named indices via compute_index().

    Years with no valid observations are excluded rather than counted as zero
    exceedances, and the statistic is None — not a fabricated 0 — when the
    record holds fewer than two usable years.
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
    """Model-vs-reference metrics on the temporally-overlapping period."""
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
    """
    Reuse ClimData's config-driven index engine for an ETCCDI/extreme index.

    Returns a scalar summary (mean / change) of the index series, or None if the
    ClimData machinery is unavailable or rejects the data (graceful
    degradation). Every failure is logged — a None here means the caller asked
    for an index it will not get, and that must be traceable.
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
