"""
ClimData Climate-Assistant Orchestrator
=======================================

Production orchestration layer wiring the full pipeline:

    Planner  →  Execution Layer  →  ClimData  →  Analytics  →  Narrator

Features
--------
* structured logging (per stage, timed)
* graceful error handling (fail-soft: a failed stage degrades output, never the
  process; clarification short-circuits cleanly)
* retry with exponential backoff for transient LLM / extraction failures
* validation at every boundary (planner, execution, narrator are self-validating)
* extensible tool registry (tool_registry.REGISTRY) for datasets / indices /
  diagnostics — future capabilities plug in without touching this file
* point analyses (lat/lon) and area analyses (lat/lon extent) via the same path
* multi-dataset role resolution (observation / model_historical / model_future)
  so "do CMIP models capture the signal?" and "future conditions" are supported

Data access is pluggable via the DataProvider protocol:
  * RealClimDataProvider  — runs ClimData(overrides=...).extract()
  * SyntheticProvider     — deterministic synthetic series, so the end-to-end
                            examples are reproducible offline and tests don't hit
                            the network. Swap in the real provider in production.
"""

from __future__ import annotations

import logging
import time
import traceback
from dataclasses import dataclass, field
from enum import Enum
from functools import lru_cache
from typing import Callable, Dict, List, Optional, Protocol

import numpy as np
import xarray as xr

from . import analytics
from . import planner as planner_mod
from . import execution as execution_mod
from . import narrator as narrator_mod
from .execution import build_execution_plan, ExecutionStatus
from .planner import ReadyPlan, ClarificationRequest, TimeRange
from .tool_registry import REGISTRY


# ===========================================================================
# Logging
# ===========================================================================

def configure_logging(level: int = logging.INFO) -> logging.Logger:
    logger = logging.getLogger("climdata.orchestrator")
    if not logger.handlers:
        h = logging.StreamHandler()
        h.setFormatter(logging.Formatter(
            "%(asctime)s | %(levelname)-7s | %(name)s | %(message)s",
            datefmt="%H:%M:%S"))
        logger.addHandler(h)
    logger.setLevel(level)
    logger.propagate = False
    return logger


LOG = configure_logging()

# Default LLM for both the planner and the narrator. Each stage can still be
# pointed at a different model via ClimateAssistant(planner_model=/narrator_model=);
# a single constant keeps them from silently drifting apart, which is what
# forced callers to monkeypatch plan()/narrate() to change the model.
DEFAULT_MODEL = "openai-gpt-oss-120b"


# ===========================================================================
# Retry
# ===========================================================================

def with_retry(fn: Callable, *args, retries: int = 3, backoff: float = 1.5,
               exceptions=(Exception,), logger: logging.Logger = LOG, what: str = "",
               **kwargs):
    """Call fn with exponential-backoff retry on transient failure."""
    delay = 1.0
    last: Optional[Exception] = None
    for attempt in range(1, retries + 1):
        try:
            return fn(*args, **kwargs)
        except exceptions as e:  # noqa: BLE001 — retry boundary
            last = e
            logger.warning("retry %d/%d for %s: %s", attempt, retries, what or fn.__name__, e)
            if attempt < retries:
                time.sleep(delay)
                delay *= backoff
    raise last


# ===========================================================================
# Data providers  (pluggable)
# ===========================================================================

class DataProvider(Protocol):
    def extract(self, overrides: List[str], role: str, spatial: Dict,
                time_range: Optional[Dict]) -> xr.Dataset: ...


class RealClimDataProvider:
    """
    Runs the actual ClimData extraction.

    `extra_overrides` is appended to every role's Hydra overrides — use it for
    settings the planner/execution layer never produces, such as dataset
    credentials (e.g. ``dsinfo.mswx.params.google_service_account=...``) or a
    non-default `data_dir`.
    """
    def __init__(self, extra_overrides: Optional[List[str]] = None):
        self.extra_overrides = list(extra_overrides or [])

    def extract(self, overrides, role, spatial, time_range) -> xr.Dataset:
        import climdata
        return climdata.ClimData(overrides=overrides + self.extra_overrides).extract()


class SyntheticProvider:
    """
    Deterministic synthetic daily data for reproducible demos/tests.

    Encodes a defensible toy signal per role so downstream diagnostics are
    meaningful: observations warm ~0.04/yr; the historical model is ~0.2 warm-
    biased but well correlated; the future scenario warms faster (~0.06/yr).

    Every variable named in the plan's ``variables=[...]`` override is emitted,
    each with its own offset and seed, so a multi-variable request is served
    with multi-variable data.
    """
    ROLE_SIGNAL = {
        "observation":      dict(base=15.0, warming=0.04, bias=0.0, noise=2.0, seed=1),
        "model_historical": dict(base=15.0, warming=0.038, bias=0.2, noise=2.1, seed=2),
        "model_future":     dict(base=16.5, warming=0.06, bias=0.0, noise=2.2, seed=3),
    }
    # per-variable offset from the role's base level (degC-like toy values)
    VARIABLE_OFFSET = {"tasmax": 0.0, "tas": -4.0, "tasmin": -8.0}

    @staticmethod
    def _requested_variables(overrides) -> List[str]:
        for ov in overrides or []:
            if ov.startswith("variables=["):
                names = [v.strip() for v in ov[len("variables=["):-1].split(",")]
                return [v for v in names if v] or ["tasmax"]
        return ["tasmax"]

    def extract(self, overrides, role, spatial, time_range) -> xr.Dataset:
        sig = self.ROLE_SIGNAL.get(role, self.ROLE_SIGNAL["observation"])
        start = (time_range or {}).get("start") or ("2015-01-01" if role == "model_future"
                                                    else "1980-01-01")
        end = (time_range or {}).get("end") or ("2050-12-31" if role == "model_future"
                                               else "2014-12-31")
        time = xr.cftime_range(start, end, freq="D", calendar="standard")
        t = np.arange(len(time))
        years = np.array([d.year for d in time])
        seasonal = 10 * np.sin(2 * np.pi * t / 365.25)
        warming = sig["warming"] * (years - years[0])

        ds = xr.Dataset()
        for n, name in enumerate(self._requested_variables(overrides)):
            rng = np.random.default_rng(sig["seed"] * 100 + n)
            series = (sig["base"] + sig["bias"] + self.VARIABLE_OFFSET.get(name, 0.0)
                      + seasonal + warming + rng.normal(0, sig["noise"], len(time)))
            da = xr.DataArray(series, coords={"time": time}, dims="time", name=name)
            da.attrs["units"] = "degC"
            ds[name] = da
        return ds


# ===========================================================================
# Pipeline result model
# ===========================================================================

class StageStatus(str, Enum):
    OK = "ok"
    SKIPPED = "skipped"
    FAILED = "failed"
    NEEDS_CLARIFICATION = "needs_clarification"


@dataclass
class StageResult:
    name: str
    status: StageStatus
    seconds: float = 0.0
    error: Optional[str] = None
    detail: Dict = field(default_factory=dict)


@dataclass
class PipelineResult:
    request: str
    status: str
    stages: List[StageResult] = field(default_factory=list)
    plan: Optional[Dict] = None
    execution: Optional[Dict] = None
    summary: Optional[Dict] = None
    narrative: Optional[str] = None
    clarification: Optional[Dict] = None
    provenance: Dict = field(default_factory=dict)


# Role → (dataset, experiment, time window) defaults. Extensible / overridable.
#
# The windows matter: a CMIP historical simulation ends in 2014 and a scenario
# run starts in 2015, so a request that pins a single period cannot be applied
# verbatim to every role. `_role_time_range` intersects the request with the
# role's window and falls back to the role window when they do not overlap.
DEFAULT_ROLE_DATASETS = {
    "observation":      {"dataset": "MSWX", "experiment": None},
    "model_historical": {"dataset": "CMIP", "experiment": "historical",
                         "time_range": {"start": None, "end": "2014-12-31"}},
    "model_future":     {"dataset": "CMIP", "experiment": "ssp585",
                         "time_range": {"start": "2015-01-01", "end": "2100-12-31"}},
}

# Which registry `type` values can serve which role. The registry already
# classifies every dataset ("Observation", "Observation (Gridded)", "ESM (Raw)",
# "ESM (Bias-Corrected)", "Simulated (Impact Model)"), so role capability is
# derived from that rather than duplicated as new config.
_ROLE_DATASET_KINDS = {
    "observation":      ("observation",),
    "model_historical": ("esm",),
    "model_future":     ("esm",),
}


@lru_cache(maxsize=1)
def _dataset_kinds() -> Dict[str, str]:
    """Map dataset key -> coarse kind ('observation' | 'esm' | 'other')."""
    from .capabilities import discover_datasets
    kinds: Dict[str, str] = {}
    for cap in discover_datasets():
        raw = str(cap.metadata.get("type") or "").strip().lower()
        if raw.startswith("observation"):
            kinds[cap.name] = "observation"
        elif raw.startswith("esm"):
            kinds[cap.name] = "esm"
        else:
            kinds[cap.name] = "other"
    return kinds


def _role_dataset(role: str, planner_dataset: Optional[str], cfg: Dict,
                  logger: logging.Logger = LOG) -> tuple[Optional[str], Optional[str]]:
    """
    Resolve the dataset for one role.

    Precedence: the planner's explicit pin (when it can actually serve the
    role) -> the role default -> None, which lets the execution layer
    auto-resolve by priority.

    The pin has to win: the role default used to be applied unconditionally, so
    "use ERA5 for observations" silently produced MSWX. But a pin cannot be
    applied blindly either — an observation product asked to serve
    `model_future` cannot run a scenario. In that case fall back to the role
    default and return a note, matching the fail-soft contract used everywhere
    else in this layer.
    """
    default = cfg.get("dataset")
    if not planner_dataset:
        return default, None

    allowed = _ROLE_DATASET_KINDS.get(role)
    kind = _dataset_kinds().get(planner_dataset)
    if allowed is None or kind is None or kind in allowed:
        return planner_dataset, None

    logger.warning("role '%s': dataset '%s' (%s) cannot serve this role; using %s",
                   role, planner_dataset, kind, default)
    return default, (f"the requested dataset '{planner_dataset}' cannot serve the "
                     f"'{role}' role, so '{default}' was used instead")


def _role_time_range(requested: Optional[Dict], role_window: Optional[Dict],
                     logger: logging.Logger = LOG, role: str = "") -> Optional[Dict]:
    """
    Intersect the requested period with the role's valid window.

    Returns the role window when the two do not overlap — a future scenario
    asked for over 1980-2014 must not be simulated over 1980-2014.
    """
    if not role_window:
        return requested
    if not requested:
        return {k: v for k, v in role_window.items() if v}

    start = max(x for x in (requested.get("start"), role_window.get("start")) if x) \
        if any((requested.get("start"), role_window.get("start"))) else None
    end = min(x for x in (requested.get("end"), role_window.get("end")) if x) \
        if any((requested.get("end"), role_window.get("end"))) else None

    if start and end and start > end:
        logger.warning("role '%s': requested period does not overlap its valid "
                       "window; using %s", role, role_window)
        return {k: v for k, v in role_window.items() if v}
    return {k: v for k, v in {"start": start, "end": end}.items() if v}


# ===========================================================================
# Orchestrator
# ===========================================================================

class ClimateAssistant:
    """
    The full pipeline, wired.

    `provider` is required and has no default: it used to fall back to
    SyntheticProvider, so omitting it in production yielded fabricated numbers
    presented as a real assessment. Pass RealClimDataProvider() to extract, or
    SyntheticProvider() explicitly for offline demos and tests.

    `planner_model` / `narrator_model` select the LLM for each stage. Both
    default to DEFAULT_MODEL; pass `model=` to set both at once.
    """

    def __init__(self, llm_client, provider: DataProvider,
                 registry=REGISTRY, logger: logging.Logger = LOG,
                 role_datasets: Dict = None, model: Optional[str] = None,
                 planner_model: Optional[str] = None,
                 narrator_model: Optional[str] = None):
        self.llm = llm_client
        self.provider = provider
        self.registry = registry
        self.log = logger
        self.role_datasets = role_datasets or DEFAULT_ROLE_DATASETS
        self.planner_model = planner_model or model or DEFAULT_MODEL
        self.narrator_model = narrator_model or model or DEFAULT_MODEL
        self._role_notes: List[str] = []

    # -- staged execution helper -------------------------------------------
    def _stage(self, name: str, fn: Callable, result: PipelineResult,
               critical: bool = True):
        t0 = time.time()
        try:
            out = fn()
            sec = round(time.time() - t0, 2)
            result.stages.append(StageResult(name, StageStatus.OK, sec))
            self.log.info("stage '%s' ok (%.2fs)", name, sec)
            return out
        except Exception as e:  # noqa: BLE001 — stage boundary
            sec = round(time.time() - t0, 2)
            self.log.error("stage '%s' FAILED: %s", name, e)
            self.log.debug(traceback.format_exc())
            result.stages.append(StageResult(name, StageStatus.FAILED, sec, error=str(e)))
            if critical:
                raise
            return None

    # -- role resolution ----------------------------------------------------
    def _required_roles(self, intents: List[str]) -> List[str]:
        roles = set()
        for intent in intents:
            for tool in self.registry.diagnostics_for_intent(intent):
                roles.update(tool.metadata.get("required_roles", []))
        roles.add("observation")  # always have a baseline
        # order matters for synthetic baseline-dependency
        order = ["observation", "model_historical", "model_future"]
        return [r for r in order if r in roles]

    # -- main entry point ---------------------------------------------------
    def run(self, nl_request: str) -> PipelineResult:
        """Run the full pipeline. Never raises — a fatal stage yields status='failed'."""
        result = PipelineResult(request=nl_request, status="running")
        try:
            return self._run(nl_request, result)
        except Exception as e:  # noqa: BLE001 — top-level safety net
            self.log.error("pipeline aborted: %s", e)
            result.status = "failed"
            result.provenance["error"] = str(e)
            return result

    def _run(self, nl_request: str, result: PipelineResult) -> PipelineResult:
        self.log.info("=" * 70)
        self.log.info("REQUEST: %s", nl_request)

        # 1) PLAN -----------------------------------------------------------
        # retries=1 here: plan() runs its own repair loop, feeding the
        # validation error back to the model. Re-invoking it wholesale only
        # duplicates that (3 x 3 = 9 LLM calls) and cannot fix a deterministic
        # ValidationError, so this layer retries transport failures only.
        plan_obj = self._stage(
            "plan",
            lambda: with_retry(planner_mod.plan, nl_request, client=self.llm,
                               model=self.planner_model,
                               retries=1, what="planner"),
            result,
        )
        if isinstance(plan_obj, ClarificationRequest):
            self.log.warning("planner needs clarification: %s", plan_obj.reason)
            result.stages[-1].status = StageStatus.NEEDS_CLARIFICATION
            result.status = "needs_clarification"
            result.clarification = plan_obj.model_dump()
            return result
        result.plan = plan_obj.model_dump(exclude_none=True)

        # 2) EXECUTION PLANS (one per dataset role) -------------------------
        roles = self._required_roles(list(plan_obj.analysis))
        self.log.info("roles required: %s", roles)
        exec_plans = self._stage(
            "execution", lambda: self._build_role_plans(plan_obj, roles), result)
        result.execution = {r: p.model_dump(exclude_none=True) for r, p in exec_plans.items()}

        # 3) CLIMDATA EXTRACTION (per role; model roles are non-critical) ---
        role_data = self._stage(
            "extract", lambda: self._extract_roles(plan_obj, exec_plans), result)

        # 4+5) ANALYTICS → structured summary -------------------------------
        summary = self._stage(
            "analytics", lambda: self._run_analytics(plan_obj, role_data), result)
        result.summary = summary

        # 6) NARRATOR -------------------------------------------------------
        # retries=1 for the same reason as the planner: narrate() owns the
        # grounding-repair loop.
        narrative = self._stage(
            "narrate",
            lambda: with_retry(narrator_mod.narrate, summary, client=self.llm,
                               model=self.narrator_model,
                               retries=1, what="narrator"),
            result, critical=False)
        if narrative is not None:
            result.narrative = narrative.text
            result.provenance["narrative_grounded"] = narrative.grounded
            if not narrative.grounded:
                self.log.warning("narrative ungrounded numbers: %s",
                                 narrative.ungrounded_numbers)

        result.status = "ok"
        result.provenance.update({
            "plan_hashes": {r: p.provenance.get("plan_hash") for r, p in exec_plans.items()},
            "provider": type(self.provider).__name__,
        })
        self.log.info("DONE: status=%s", result.status)
        return result

    # -- helpers ------------------------------------------------------------
    def _build_role_plans(self, plan_obj: ReadyPlan, roles: List[str]) -> Dict:
        plans: Dict[str, execution_mod.ExecutionPlan] = {}
        self._role_notes = []
        for role in roles:
            cfg = self.role_datasets.get(role, {})
            dataset, note = _role_dataset(role, plan_obj.dataset, cfg, self.log)
            if note:
                self._role_notes.append(note)

            # The future/scenario role gets its own explicitly-named period
            # when the user gave one (e.g. "2005-2014 as reference, 2040-2059
            # as future") — a single time_range cannot hold both at once.
            source_range = (plan_obj.future_time_range
                            if role == "model_future" and plan_obj.future_time_range
                            else plan_obj.time_range)
            requested = ({"start": source_range.start, "end": source_range.end}
                         if source_range else None)
            window = _role_time_range(requested, cfg.get("time_range"),
                                      self.log, role)

            # Rebuild through the schema rather than assigning onto a copy:
            # ReadyPlan has no validate_assignment, so a mutated dataset would
            # never be re-checked against the requested variables.
            req = ReadyPlan.model_validate({
                **plan_obj.model_dump(),
                "dataset": dataset,
                "time_range": window or None,
            })

            # Scenario/GCM travel INTO the builder so the plan hash covers them:
            # two runs differing only in experiment are different extractions.
            #
            # A user-named scenario applies to the FUTURE role only —
            # model_historical is by definition the `historical` simulation, so
            # letting "SSP2-4.5" reach it would run the historical baseline as a
            # scenario. The GCM, in contrast, applies to both model roles: the
            # comparison is only meaningful within one model.
            extra_ov, extra_cfg = [], {}
            experiment = cfg.get("experiment")
            if role == "model_future" and plan_obj.experiment_id:
                experiment = plan_obj.experiment_id
            if experiment and role != "observation":
                extra_ov.append(f"experiment_id={experiment}")
                extra_cfg["experiment_id"] = experiment
            if plan_obj.source_id and role != "observation":
                extra_ov.append(f"source_id={plan_obj.source_id}")
                extra_cfg["source_id"] = plan_obj.source_id

            ep = build_execution_plan(req, extra_overrides=extra_ov,
                                      extra_config=extra_cfg)
            if ep.status == ExecutionStatus.ERROR.value:
                self.log.warning("role '%s' execution plan has gaps: %s",
                                 role, [m.detail for m in ep.missing])
            plans[role] = ep
        return plans

    def _extract_roles(self, plan_obj: ReadyPlan, exec_plans: Dict) -> Dict:
        data: Dict[str, xr.Dataset] = {}
        for role, ep in exec_plans.items():
            critical = role == "observation"
            try:
                ds = with_retry(self.provider.extract, ep.overrides, role,
                                ep.spatial, ep.time_range, retries=2,
                                what=f"extract:{role}")
                data[role] = ds
                self.log.info("extracted role '%s' (%s)", role, ep.dataset)
            except Exception as e:  # noqa: BLE001
                if critical:
                    raise
                self.log.warning("non-critical role '%s' extraction failed: %s — skipping", role, e)
        return data

    @staticmethod
    def _is_unusable(value) -> bool:
        """True for a non-finite scalar — NaN/inf carry no information."""
        return isinstance(value, float) and not np.isfinite(value)

    @classmethod
    def _drop_unusable(cls, value):
        """
        Strip non-finite numbers, recursing into nested containers.

        Returns None when nothing usable survives, so the caller can turn the
        statistic into a note instead of reporting it.
        """
        if isinstance(value, dict):
            cleaned = {k: v for k, v in value.items() if not cls._is_unusable(v)}
            return cleaned or None
        return None if cls._is_unusable(value) else value

    @classmethod
    def _merge(cls, summary: Dict, produced: Dict, prefix: str = "") -> None:
        """
        Fold one diagnostic's output into the summary, accumulating notes.

        Non-finite statistics are dropped and recorded as notes rather than
        published. Analytics legitimately returns NaN for a degenerate record,
        but a bare `annual_trend = nan` in the narrator prompt is read as a
        statistic: `narrator._quotes` cannot match any number against NaN, so
        every figure the model writes gets flagged ungrounded and the repair
        loop burns all its retries before failing closed.
        """
        label = prefix.rstrip("_")
        for key, value in produced.items():
            if key == "notes":
                # notes are prose, so they take a readable label, not a key prefix
                summary.setdefault("notes", []).extend(
                    f"{label}: {v}" if label else v for v in value)
                continue
            usable = cls._drop_unusable(value)
            if usable is None and value is not None:
                note = f"{key} could not be computed from the extracted data"
                summary.setdefault("notes", []).append(
                    f"{label}: {note}" if label else note)
                continue
            summary[f"{prefix}{key}" if prefix else key] = usable

    def _run_analytics(self, plan_obj: ReadyPlan, role_data: Dict) -> Dict:
        """
        Run every requested diagnostic for every requested variable.

        With more than one variable each key is prefixed with the variable name,
        so a two-variable request yields two sets of statistics instead of
        silently reporting only the first.
        """
        variables = plan_obj.variables
        summary: Dict = {}
        for variable in variables:
            prefix = f"{variable}_" if len(variables) > 1 else ""
            produced_any = False
            for intent in plan_obj.analysis:
                for tool in self.registry.diagnostics_for_intent(intent):
                    needed = tool.metadata.get("required_roles", [])
                    if not all(r in role_data for r in needed):
                        self.log.warning("skip diagnostic '%s' (missing roles %s)",
                                         tool.name,
                                         [r for r in needed if r not in role_data])
                        continue
                    kw = {"variable": variable}
                    if tool.name in ("extremes", "etccdi_indices"):
                        kw["indices"] = plan_obj.indices or []
                    try:
                        self._merge(summary, tool.fn(role_data, **kw), prefix)
                        produced_any = True
                    except Exception as e:  # noqa: BLE001 — one diagnostic must not sink the rest
                        self.log.warning("diagnostic '%s' failed for '%s': %s",
                                         tool.name, variable, e, exc_info=True)
            # Units are supplementary context for whatever WAS computed — an
            # intent with no registered diagnostic must still yield an empty
            # summary, not a units key with nothing to attach it to.
            if produced_any:
                units = analytics.variable_units(variable)
                if units:
                    summary[f"{prefix}units"] = units

        # A dataset substitution changes what the numbers describe, so it is a
        # caveat the assessment has to carry.
        if summary and self._role_notes:
            summary.setdefault("notes", []).extend(self._role_notes)
        return summary


# ===========================================================================
# End-to-end examples
# ===========================================================================

SAIA_BASE_URL = "https://chat-ai.academiccloud.de/v1"


def build_client(api_key: Optional[str] = None):
    """
    Build the SAIA client.

    Prefers an explicitly passed `api_key`; otherwise falls back to the
    SAIA_API_KEY / OPENAI_API_KEY environment variables. Never hardcode a key
    at the call site — read it from a secrets manager or env var and pass it
    through.
    """
    import os
    from openai import OpenAI

    api_key = api_key or os.environ.get("SAIA_API_KEY") or os.environ.get("OPENAI_API_KEY")
    if not api_key:
        raise RuntimeError(
            "No API key found. Pass api_key=... or export SAIA_API_KEY "
            "(or OPENAI_API_KEY) before running the examples.")
    return OpenAI(api_key=api_key,
                  base_url=os.environ.get("SAIA_BASE_URL", SAIA_BASE_URL))


def print_result(result: PipelineResult):
    """Human-readable dump of a finished pipeline run."""
    import json
    print("\n" + "#" * 72)
    print("STATUS:", result.status)
    if result.clarification:
        print("CLARIFICATION:", json.dumps(result.clarification, indent=2)); return
    print("\nPLAN:", json.dumps(result.plan, indent=2))
    print("\nSUMMARY:", json.dumps(result.summary, indent=2))
    print("\nNARRATIVE:\n", result.narrative)
    print("\nPROVENANCE:", json.dumps(result.provenance, indent=2))
    print("STAGES:", [(s.name, s.status.value, s.seconds) for s in result.stages])


# Backwards-compatible aliases — these were the de-facto public API (notebooks
# call them) before they were given non-underscored names.
_client = build_client
_print = print_result


if __name__ == "__main__":
    assistant = ClimateAssistant(llm_client=build_client(),
                                 provider=SyntheticProvider())

    print("\n>>> EXAMPLE 1 — BERLIN (point)")
    berlin = (
        "Describe the historical trend of daily maximum temperature observations "
        "over Berlin (lat=52.5200, lon=13.4050). Check whether CMIP models capture "
        "the broad climate change signal. Also describe plausible future climate "
        "conditions for Berlin."
    )
    print_result(assistant.run(berlin))

    print("\n>>> EXAMPLE 2 — MÄRKISCH-ODERLAND (area)")
    mol = (
        "Describe historical and projected daily maximum temperature changes over "
        "the Märkisch-Oderland district using the extent lat_min=52.20, lat_max=53.55, "
        "lon_min=13.60, lon_max=15.25."
    )
    print_result(assistant.run(mol))
