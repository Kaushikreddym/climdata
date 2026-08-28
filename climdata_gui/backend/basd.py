"""Bias adjustment & statistical downscaling for the BASD tab.

Three extractions feed every method (the standard ISIMIP framing):

    obs_hist   observations over the reference period
    sim_hist   model output over the *same* reference period
    sim_fut    model output over the target period — the series being corrected

``qdm`` / ``dqm`` / ``qm`` are quantile mappings from ``xsdba``; ``bcsd``
delegates to climdata's ISIMIP3BASD :class:`~climdata.sdba.BiasCorrection`.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from datetime import date
from typing import Callable, Dict, List, Optional

from .config_builder import build_overrides
from .runner import run_pipeline

# Datasets whose runs are split into a `historical` simulation and an SSP
# scenario, so the reference and target windows need different experiments.
ESM_DATASETS = ("cmip", "cmip_w5e5", "nexgddp")

# CMIP6 protocol boundaries between the historical run and the scenario runs.
HISTORICAL_END = date(2014, 12, 31)
SCENARIO_START = date(2015, 1, 1)

METHOD_LABELS = {
    "bcsd": "BCSD — Bias Correction & Spatial Disaggregation (ISIMIP3BASD)",
    "qdm":  "QDM — Quantile Delta Mapping",
    "dqm":  "DQM — Detrended Quantile Mapping",
    "qm":   "QM — Empirical Quantile Mapping",
}

# Multiplicative correction suits strictly-positive, skewed variables.
_MULTIPLICATIVE_VARS = {"pr", "evpot", "et0"}

_SPATIAL_DIMS = ("lat", "latitude", "lon", "longitude", "x", "y")


@dataclass
class BasdResult:
    """Outcome of one BASD run."""
    dataset: object = None                     # corrected xr.Dataset (target period)
    filename: Optional[str] = None
    method: str = ""
    variable: str = ""
    metrics: Dict[str, Dict[str, float]] = field(default_factory=dict)
    notes: List[str] = field(default_factory=list)


def is_esm_dataset(dataset: Optional[str]) -> bool:
    """True for datasets split into historical + scenario experiments."""
    return bool(dataset) and dataset.lower() in ESM_DATASETS


# ---------------------------------------------------------------------------
# Extraction helpers
# ---------------------------------------------------------------------------

def _extract(dataset: str, variable: str, start: date, end: date,
             lat, lon, box, data_dir, experiment_id, source_id, log):
    """Run climdata's extract action and return the xarray Dataset."""
    overrides = build_overrides(
        dataset=dataset,
        lat=lat, lon=lon, box=box,
        start_date=start, end_date=end,
        variables=[variable],
        data_dir=data_dir,
        experiment_id=experiment_id,
        source_id=source_id,
    )
    log(f"  extract {dataset} [{variable}] {start} → {end}"
        + (f"  ({experiment_id})" if experiment_id else ""))
    result = run_pipeline(overrides, ["extract"])
    ds = getattr(result, "dataset", None)
    if ds is None or variable not in ds:
        raise ValueError(
            f"Extraction of '{variable}' from '{dataset}' returned no data "
            f"for {start} → {end}."
        )
    return ds


def _clamp_to_experiment(start: date, end: date, experiment: str,
                         notes: List[str], label: str) -> tuple:
    """Clip a window to the period its CMIP6 experiment actually covers."""
    if experiment == "historical" and end > HISTORICAL_END:
        notes.append(f"{label}: end clipped to {HISTORICAL_END} "
                     f"(the historical run stops there).")
        end = HISTORICAL_END
    if experiment != "historical" and start < SCENARIO_START:
        notes.append(f"{label}: start clipped to {SCENARIO_START} "
                     f"(scenario runs begin there).")
        start = SCENARIO_START
    if start > end:
        raise ValueError(
            f"{label}: no overlap between the requested window and the "
            f"'{experiment}' experiment period."
        )
    return start, end


# ---------------------------------------------------------------------------
# Grid / calendar harmonisation
# ---------------------------------------------------------------------------

def _spatial_dims(da) -> List[str]:
    return [d for d in da.dims if d in _SPATIAL_DIMS]


def _match_grid(sim, ref, notes: List[str], log):
    """Put *sim* on *ref*'s grid when the two differ.

    Only geographic (lat/lon) grids are interpolated — a projected grid such as
    HYRAS's LAEA x/y cannot be matched to degrees by interpolation alone, and
    silently doing so would be wrong.
    """
    sim_dims, ref_dims = _spatial_dims(sim), _spatial_dims(ref)
    if not sim_dims and not ref_dims:
        return sim                                    # point vs point
    if sorted(sim_dims) == sorted(ref_dims) and all(
        sim.sizes[d] == ref.sizes[d] for d in sim_dims
    ):
        return sim                                    # already aligned

    geographic = {"lat", "latitude", "lon", "longitude"}
    if not (set(sim_dims) <= geographic and set(ref_dims) <= geographic):
        raise ValueError(
            "Reference and model grids differ and at least one is a projected "
            f"grid (dims {ref_dims} vs {sim_dims}). Select a single point on "
            "the map, or pick datasets that share a grid."
        )
    coords = {d: ref[d] for d in ref_dims if d in sim.coords}
    if not coords:
        raise ValueError(
            f"Cannot align model grid {sim_dims} to reference grid {ref_dims}."
        )
    notes.append("Model interpolated onto the reference grid (bilinear); this "
                 "is a simple interpolation, not conservative regridding.")
    log("  interpolating model onto the reference grid…")
    return sim.interp(coords, method="linear")


def _match_calendar(arrays: List, notes: List[str]):
    """Convert every array to one calendar so quantile mapping can group them."""
    calendars = {getattr(a["time"].dt, "calendar", "standard") for a in arrays}
    if len(calendars) <= 1:
        return arrays
    notes.append(f"Calendars differed ({', '.join(sorted(calendars))}); "
                 f"all series converted to 'noleap'.")
    return [a.convert_calendar("noleap", align_on="date") for a in arrays]


# ---------------------------------------------------------------------------
# Methods
# ---------------------------------------------------------------------------

def _single_time_chunk(da):
    """Collapse the time dimension into one dask chunk (no-op if in memory)."""
    return da.chunk({"time": -1}) if da.chunks else da


def _adjust_xsdba(method: str, ref, hist, sim, variable: str,
                  n_quantiles: int, notes: List[str], log):
    """Train a quantile mapping on (ref, hist) and adjust both hist and sim."""
    from xsdba import Grouper
    from xsdba.adjustment import (
        DetrendedQuantileMapping,
        EmpiricalQuantileMapping,
        QuantileDeltaMapping,
    )

    classes = {
        "qdm": QuantileDeltaMapping,
        "dqm": DetrendedQuantileMapping,
        "qm":  EmpiricalQuantileMapping,
    }
    kind = "*" if variable in _MULTIPLICATIVE_VARS else "+"
    group = Grouper("time.dayofyear", window=31)

    # xsdba refuses to group over a dimension split across dask chunks.
    ref, hist, sim = (_single_time_chunk(a) for a in (ref, hist, sim))

    log(f"  training {method.upper()} (kind='{kind}', "
        f"{n_quantiles} quantiles, day-of-year ±15)…")
    adjuster = classes[method].train(
        ref, hist, nquantiles=n_quantiles, group=group, kind=kind,
    )
    log("  adjusting the target period…")
    corrected_fut = adjuster.adjust(sim, interp="linear")
    log("  adjusting the reference period (for validation)…")
    corrected_hist = adjuster.adjust(hist, interp="linear")
    return corrected_fut, corrected_hist


def _adjust_bcsd(ref_ds, hist_ds, sim_ds, variable: str, log):
    """ISIMIP3BASD trend-preserving quantile mapping via climdata.sdba."""
    try:
        from climdata.sdba import BiasCorrection
    except ImportError as exc:
        raise ImportError(
            f"BCSD needs climdata's sdba module: {exc}"
        ) from exc
    log("  running ISIMIP3BASD bias correction…")
    try:
        corrected = BiasCorrection(variable=variable).correct(
            obs_hist=ref_ds, sim_hist=hist_ds, sim_fut=sim_ds,
        )
    except ImportError as exc:
        raise ImportError(
            "BCSD (ISIMIP3BASD) requires scitools-iris, which is not "
            f"installed: {exc}\n"
            "Install it (conda install -c conda-forge scitools-iris) or choose "
            "QDM / DQM / QM, which need only xsdba."
        ) from exc
    return corrected[variable] if variable in corrected else corrected


# ---------------------------------------------------------------------------
# Metrics & output
# ---------------------------------------------------------------------------

def _metrics(model, reference, variable: str) -> Dict[str, float]:
    from climdata.LLM.analytics import evaluation_metrics
    return evaluation_metrics(model, reference, variable)


def _write(ds, out_dir: Optional[str], stem: str, fmt: str, log) -> Optional[str]:
    """Persist the corrected dataset as NetCDF or CSV."""
    target_dir = out_dir or "."
    os.makedirs(target_dir, exist_ok=True)
    if fmt.upper() == "CSV":
        path = os.path.join(target_dir, f"{stem}.csv")
        ds.to_dataframe().reset_index().to_csv(path, index=False)
    else:
        path = os.path.join(target_dir, f"{stem}.nc")
        ds.to_netcdf(path)
    log(f"  wrote {path}")
    return path


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def run_basd(
    *,
    ref_dataset: str,
    ref_start: date,
    ref_end: date,
    target_dataset: str,
    target_start: date,
    target_end: date,
    method: str,
    variable: str,
    lat: Optional[float] = None,
    lon: Optional[float] = None,
    box: Optional[dict] = None,
    data_dir: Optional[str] = None,
    out_dir: Optional[str] = None,
    out_format: str = "NetCDF",
    experiment_id: Optional[str] = None,
    source_id: Optional[str] = None,
    n_quantiles: int = 50,
    log: Callable[[str], None] = print,
) -> BasdResult:
    """Bias-adjust *target_dataset* against *ref_dataset* and save the result.

    Args:
        ref_dataset: Observational reference provider (e.g. ``"mswx"``).
        ref_start/ref_end: Reference (calibration) period.
        target_dataset: Model provider to correct (e.g. ``"cmip"``).
        target_start/target_end: Period of model output to correct.
        method: One of ``bcsd``, ``qdm``, ``dqm``, ``qm``.
        variable: Single climdata variable key (e.g. ``"tasmax"``).
        experiment_id: Scenario for the *target* window; the reference window
            of an ESM always uses ``historical``.
        source_id: CMIP6 model name for ESM targets.

    Returns:
        BasdResult with the corrected dataset, output path, and — for the
        quantile-mapping methods — before/after validation metrics.
    """
    method = (method or "qdm").lower()
    if method not in METHOD_LABELS:
        raise ValueError(f"Unknown method '{method}'. "
                         f"Choose from {', '.join(METHOD_LABELS)}.")

    notes: List[str] = []
    log(f"BASD: {METHOD_LABELS[method]}")
    log(f"  variable={variable}  reference={ref_dataset}  target={target_dataset}")

    # ── Experiments: the reference window of an ESM is its historical run ──
    hist_experiment = fut_experiment = None
    hist_start, hist_end = ref_start, ref_end
    fut_start, fut_end = target_start, target_end
    if is_esm_dataset(target_dataset):
        hist_experiment = "historical"
        fut_experiment = experiment_id or "ssp585"
        hist_start, hist_end = _clamp_to_experiment(
            hist_start, hist_end, hist_experiment, notes, "Reference period")
        fut_start, fut_end = _clamp_to_experiment(
            fut_start, fut_end, fut_experiment, notes, "Target period")

    # ── 1. Extract the three series ───────────────────────────────────
    log("Extracting reference observations…")
    obs_ds = _extract(ref_dataset, variable, ref_start, ref_end,
                      lat, lon, box, data_dir, None, None, log)
    log("Extracting model output over the reference period…")
    hist_ds = _extract(target_dataset, variable, hist_start, hist_end,
                       lat, lon, box, data_dir, hist_experiment, source_id, log)
    log("Extracting model output over the target period…")
    fut_ds = _extract(target_dataset, variable, fut_start, fut_end,
                      lat, lon, box, data_dir, fut_experiment, source_id, log)

    # ── 2. Harmonise grid & calendar ──────────────────────────────────
    ref_da = obs_ds[variable]
    hist_da = _match_grid(hist_ds[variable], ref_da, notes, log)
    fut_da = _match_grid(fut_ds[variable], ref_da, notes, log)
    ref_da, hist_da, fut_da = _match_calendar([ref_da, hist_da, fut_da], notes)

    # ── 3. Correct ────────────────────────────────────────────────────
    metrics: Dict[str, Dict[str, float]] = {}
    if method == "bcsd":
        corrected = _adjust_bcsd(
            ref_da.to_dataset(name=variable),
            hist_da.to_dataset(name=variable),
            fut_da.to_dataset(name=variable),
            variable, log,
        )
    else:
        corrected, corrected_hist = _adjust_xsdba(
            method, ref_da, hist_da, fut_da, variable, n_quantiles, notes, log)
        log("Validating over the reference period…")
        metrics["before"] = _metrics(hist_da, ref_da, variable)
        metrics["after"] = _metrics(corrected_hist, ref_da, variable)

    corrected = corrected.rename(variable)
    corrected.attrs.setdefault("units", ref_da.attrs.get("units", ""))
    corrected.attrs["bias_adjustment"] = (
        f"{method.upper()} against {ref_dataset} "
        f"({ref_start} → {ref_end})"
    )
    out_ds = corrected.to_dataset(name=variable)

    # ── 4. Save ───────────────────────────────────────────────────────
    stem = (f"{variable}_{method}_{target_dataset}_vs_{ref_dataset}_"
            f"{fut_start:%Y%m%d}_{fut_end:%Y%m%d}")
    filename = _write(out_ds, out_dir or data_dir, stem, out_format, log)

    for note in notes:
        log(f"  note: {note}")
    return BasdResult(dataset=out_ds, filename=filename, method=method,
                      variable=variable, metrics=metrics, notes=notes)
