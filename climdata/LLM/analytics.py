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

from typing import Dict, List, Optional

import numpy as np
import xarray as xr


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


def _as_dataarray(data: xr.Dataset | xr.DataArray, variable: Optional[str]) -> xr.DataArray:
    if isinstance(data, xr.DataArray):
        return data
    if variable is None:
        if len(data.data_vars) != 1:
            raise ValueError(f"Specify `variable`; dataset has {list(data.data_vars)}.")
        variable = next(iter(data.data_vars))
    return data[variable]


def _reduce_space(da: xr.DataArray) -> xr.DataArray:
    """Average away any spatial dims → a pure time series."""
    dims = [d for d in da.dims if d in _SPATIAL_DIMS]
    return da.mean(dims) if dims else da


def _annual_mean(da: xr.DataArray) -> xr.DataArray:
    return _reduce_space(da).groupby("time.year").mean()


def _linear_slope(years: np.ndarray, values: np.ndarray) -> float:
    """Least-squares slope per year, ignoring NaNs."""
    m = np.isfinite(values)
    if m.sum() < 2:
        return float("nan")
    return float(np.polyfit(years[m], values[m], 1)[0])


# ===========================================================================
# TRENDS  (implemented — not in ClimData)
# ===========================================================================

def annual_trend(data, variable: Optional[str] = None) -> float:
    """Linear trend of the annual-mean series, in <units> per year."""
    s = _annual_mean(_as_dataarray(data, variable))
    return round(_linear_slope(s["year"].values, s.values), 4)


def seasonal_trend(data, season: str = "JJA", variable: Optional[str] = None) -> float:
    """Linear trend of a single season's annual-mean series, in <units> per year."""
    da = _reduce_space(_as_dataarray(data, variable))
    da = da.sel(time=da["time"].dt.season == season)
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
                    window: int = 10) -> int:
    """
    Change in the annual count of days above `thresh`, comparing the mean of the
    last `window` years to the first `window` years. Threshold-based extreme;
    falls back to ClimData ETCCDI for named indices via compute_index().
    """
    da = _reduce_space(_as_dataarray(data, variable))
    counts = (da > thresh).groupby("time.year").sum()
    vals = counts.values.astype(float)
    if len(vals) >= 2 * window:
        return int(round(float(vals[-window:].mean() - vals[:window].mean())))
    # short series → slope projected over the full span
    yrs = counts["year"].values
    return int(round(_linear_slope(yrs, vals) * (yrs[-1] - yrs[0])))


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
        return {"model_bias": float("nan"), "rmse": float("nan"),
                "mae": float("nan"), "correlation": float("nan")}
    diff = mv - ov
    return {
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
    ClimData machinery is unavailable in this environment (graceful degradation).
    """
    try:
        from climdata.utils.wrapper_workflow import ClimateExtractor
    except Exception:
        return None
    try:
        extractor = ClimateExtractor(overrides=[f"dataset={dataset}", f"index={index}"])
        idx_ds = extractor.calc_index(ds=data if isinstance(data, xr.Dataset) else data.to_dataset())
        if idx_ds is None:
            return None
        series = _reduce_space(idx_ds[index]) if index in idx_ds else _reduce_space(
            idx_ds[next(iter(idx_ds.data_vars))])
        vals = series.values.astype(float)
        if reduce == "change" and vals.size >= 2:
            return round(float(vals[-1] - vals[0]), 4)
        return round(float(np.nanmean(vals)), 4)
    except Exception:
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

    Keys are only present for the analyses actually requested/possible.
    """
    out: Dict[str, float] = {}

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
            if val is not None:
                out[f"{idx}_change"] = val
        # threshold extreme always available without ClimData
        out["hot_days_change"] = hot_days_change(data, hot_day_thresh, variable)

    if "model_evaluation" in intents and reference is not None:
        out.update(evaluation_metrics(data, reference, variable))

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
