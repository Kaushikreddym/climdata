"""Side-by-side comparison of two climate datasets for the Comparison tab.

Extracts both selections over the shared area of interest, reduces each to a
spatial-mean time series, and returns summary statistics, overlap metrics, and
a rendered figure (annual means + monthly climatology).
"""

from __future__ import annotations

import base64
import io
from datetime import date
from typing import Callable, Dict, Optional

from .config_builder import build_overrides
from .runner import run_pipeline

_SPATIAL_DIMS = ("lat", "latitude", "lon", "longitude", "x", "y")

# Readable on both the light and dark dashboard themes.
_INK = "#8b94a7"
_SERIES_COLORS = ("#4f91ff", "#f6a54a")


def _extract(dataset: str, variable: str, start: date, end: date,
             lat, lon, box, data_dir, experiment_id, source_id, log):
    overrides = build_overrides(
        dataset=dataset,
        lat=lat, lon=lon, box=box,
        start_date=start, end_date=end,
        variables=[variable],
        data_dir=data_dir,
        experiment_id=experiment_id,
        source_id=source_id,
    )
    log(f"  extract {dataset} [{variable}] {start} → {end}")
    result = run_pipeline(overrides, ["extract"])
    ds = getattr(result, "dataset", None)
    if ds is None or variable not in ds:
        raise ValueError(
            f"Extraction of '{variable}' from '{dataset}' returned no data "
            f"for {start} → {end}."
        )
    return ds[variable]


def _reduce_space(da):
    """Average away every non-time dimension to a pure time series."""
    spatial = [d for d in da.dims if d in _SPATIAL_DIMS]
    da = da.mean(spatial) if spatial else da
    leftover = [d for d in da.dims if d != "time"]
    return da.mean(leftover) if leftover else da


def _stats(da, variable: str) -> Dict[str, object]:
    """Summary statistics of one spatially-reduced series."""
    import numpy as np
    from climdata.LLM.analytics import annual_trend

    values = da.values.astype(float)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        raise ValueError(f"'{variable}' contains no finite values.")
    times = da["time"].values
    try:
        trend = annual_trend(da, variable)
    except Exception:      # noqa: BLE001 — a trend is supplementary, never fatal
        trend = float("nan")
    return {
        "mean":  round(float(finite.mean()), 4),
        "std":   round(float(finite.std()), 4),
        "min":   round(float(finite.min()), 4),
        "max":   round(float(finite.max()), 4),
        "trend_per_decade": (None if trend != trend else round(trend * 10, 4)),
        "n":     int(finite.size),
        "start": str(times[0])[:10],
        "end":   str(times[-1])[:10],
        "units": da.attrs.get("units", ""),
    }


def _render_figure(series_a, series_b, label_a: str, label_b: str,
                   variable: str, units: str) -> str:
    """Annual means + monthly climatology for both series, as a base64 PNG."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return ""

    fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.6), dpi=110)
    fig.patch.set_alpha(0)

    for ax in axes:
        ax.patch.set_alpha(0)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        for spine in ("left", "bottom"):
            ax.spines[spine].set_color(_INK)
        ax.tick_params(colors=_INK, labelsize=8)
        ax.grid(True, alpha=0.15, linewidth=0.6)

    ylabel = f"{variable} [{units}]" if units else variable

    # ── Annual means ──
    for (series, label, color) in ((series_a, label_a, _SERIES_COLORS[0]),
                                   (series_b, label_b, _SERIES_COLORS[1])):
        annual = series.groupby("time.year").mean()
        axes[0].plot(annual["year"].values, annual.values,
                     label=label, color=color, linewidth=1.6)
    axes[0].set_title("Annual mean", color=_INK, fontsize=9.5, loc="left")
    axes[0].set_ylabel(ylabel, color=_INK, fontsize=8.5)
    legend = axes[0].legend(fontsize=8, frameon=False)
    for text in legend.get_texts():
        text.set_color(_INK)

    # ── Monthly climatology ──
    for (series, label, color) in ((series_a, label_a, _SERIES_COLORS[0]),
                                   (series_b, label_b, _SERIES_COLORS[1])):
        monthly = series.groupby("time.month").mean()
        axes[1].plot(monthly["month"].values, monthly.values,
                     label=label, color=color, linewidth=1.6, marker="o",
                     markersize=3)
    axes[1].set_title("Monthly climatology", color=_INK, fontsize=9.5, loc="left")
    axes[1].set_ylabel(ylabel, color=_INK, fontsize=8.5)
    axes[1].set_xticks(range(1, 13))
    axes[1].set_xticklabels(["J", "F", "M", "A", "M", "J",
                             "J", "A", "S", "O", "N", "D"])

    fig.tight_layout(pad=1.2)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", transparent=True, bbox_inches="tight", dpi=110)
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode()


def run_comparison(
    *,
    a_dataset: str,
    a_variable: str,
    a_start: date,
    a_end: date,
    b_dataset: str,
    b_variable: str,
    b_start: date,
    b_end: date,
    lat: Optional[float] = None,
    lon: Optional[float] = None,
    box: Optional[dict] = None,
    data_dir: Optional[str] = None,
    a_experiment: Optional[str] = None,
    a_source: Optional[str] = None,
    b_experiment: Optional[str] = None,
    b_source: Optional[str] = None,
    log: Callable[[str], None] = print,
) -> Dict:
    """Compare two dataset/variable/period selections over the current AOI.

    Returns:
        ``{"a": {...}, "b": {...}, "overlap": {...} | None, "plot_b64": str}``
        where each side carries its label and summary statistics, and
        ``overlap`` holds bias / RMSE / MAE / correlation when the two periods
        intersect.
    """
    log(f"Comparing {a_dataset} [{a_variable}] vs {b_dataset} [{b_variable}]")

    log("Extracting dataset A…")
    da_a = _reduce_space(_extract(a_dataset, a_variable, a_start, a_end,
                                  lat, lon, box, data_dir,
                                  a_experiment, a_source, log))
    log("Extracting dataset B…")
    da_b = _reduce_space(_extract(b_dataset, b_variable, b_start, b_end,
                                  lat, lon, box, data_dir,
                                  b_experiment, b_source, log))

    label_a = f"{a_dataset} · {a_variable}"
    label_b = f"{b_dataset} · {b_variable}"

    stats_a = _stats(da_a, a_variable)
    stats_b = _stats(da_b, b_variable)

    # The two sides often differ only by period, so the plot legend carries the
    # years even when the dataset and variable are identical.
    legend_a = f"{label_a} ({stats_a['start'][:4]}–{stats_a['end'][:4]})"
    legend_b = f"{label_b} ({stats_b['start'][:4]}–{stats_b['end'][:4]})"

    # ── Overlap metrics: only meaningful where the records coincide ──
    overlap = None
    try:
        import xarray as xr
        aligned_a, aligned_b = xr.align(da_a, da_b, join="inner")
        if aligned_a["time"].size > 0:
            from climdata.LLM.analytics import evaluation_metrics
            metrics = evaluation_metrics(da_a, da_b, a_variable)
            overlap = {
                "n_days": int(aligned_a["time"].size),
                "start": str(aligned_a["time"].values[0])[:10],
                "end": str(aligned_a["time"].values[-1])[:10],
                # evaluation_metrics names its arguments model/reference:
                # here A plays the model and B the reference.
                "bias": metrics["model_bias"],
                "rmse": metrics["rmse"],
                "mae": metrics["mae"],
                "correlation": metrics["correlation"],
            }
            log(f"  overlapping period: {overlap['start']} → {overlap['end']} "
                f"({overlap['n_days']} steps)")
        else:
            log("  periods do not overlap — reporting each series separately.")
    except Exception as exc:      # noqa: BLE001 — metrics are supplementary
        log(f"  overlap metrics unavailable: {exc}")

    log("Rendering comparison figure…")
    units = stats_a["units"] or stats_b["units"]
    plot_b64 = _render_figure(da_a, da_b, legend_a, legend_b, a_variable, units)

    return {
        "a": {"label": label_a, "stats": stats_a},
        "b": {"label": label_b, "stats": stats_b},
        "overlap": overlap,
        "difference": round(stats_a["mean"] - stats_b["mean"], 4)
                      if a_variable == b_variable else None,
        "units": units,
        "plot_b64": plot_b64,
    }
