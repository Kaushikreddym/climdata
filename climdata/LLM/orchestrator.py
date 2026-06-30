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
from typing import Callable, Dict, List, Optional, Protocol

import numpy as np
import xarray as xr

import planner as planner_mod
import execution as execution_mod
import narrator as narrator_mod
from execution import build_execution_plan, ExecutionStatus
from planner import ReadyPlan, ClarificationRequest
from tool_registry import REGISTRY


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
    """Runs the actual ClimData extraction."""
    def extract(self, overrides, role, spatial, time_range) -> xr.Dataset:
        import climdata
        return climdata.ClimData(overrides=overrides).extract()


class SyntheticProvider:
    """
    Deterministic synthetic daily tasmax for reproducible demos/tests.

    Encodes a defensible toy signal per role so downstream diagnostics are
    meaningful: observations warm ~0.04/yr; the historical model is ~0.2 warm-
    biased but well correlated; the future scenario warms faster (~0.06/yr).
    """
    ROLE_SIGNAL = {
        "observation":      dict(base=15.0, warming=0.04, bias=0.0, noise=2.0, seed=1),
        "model_historical": dict(base=15.0, warming=0.038, bias=0.2, noise=2.1, seed=2),
        "model_future":     dict(base=16.5, warming=0.06, bias=0.0, noise=2.2, seed=3),
    }

    def extract(self, overrides, role, spatial, time_range) -> xr.Dataset:
        sig = self.ROLE_SIGNAL.get(role, self.ROLE_SIGNAL["observation"])
        start = (time_range or {}).get("start") or ("2015-01-01" if role == "model_future"
                                                    else "1980-01-01")
        end = (time_range or {}).get("end") or ("2050-12-31" if role == "model_future"
                                               else "2014-12-31")
        time = xr.cftime_range(start, end, freq="D", calendar="standard")
        t = np.arange(len(time))
        years = np.array([d.year for d in time])
        rng = np.random.default_rng(sig["seed"])
        seasonal = 10 * np.sin(2 * np.pi * t / 365.25)
        warming = sig["warming"] * (years - years[0])
        series = (sig["base"] + sig["bias"] + seasonal + warming
                  + rng.normal(0, sig["noise"], len(time)))
        da = xr.DataArray(series, coords={"time": time}, dims="time", name="tasmax")
        da.attrs["units"] = "degC"
        return da.to_dataset()


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


# Role → (dataset, experiment) defaults. Extensible / overridable.
DEFAULT_ROLE_DATASETS = {
    "observation":      {"dataset": None,  "experiment": None},        # planner's or auto
    "model_historical": {"dataset": "CMIP", "experiment": "historical"},
    "model_future":     {"dataset": "CMIP", "experiment": "ssp585"},
}


# ===========================================================================
# Orchestrator
# ===========================================================================

class ClimateAssistant:
    def __init__(self, llm_client, provider: DataProvider = None,
                 registry=REGISTRY, logger: logging.Logger = LOG,
                 role_datasets: Dict = None):
        self.llm = llm_client
        self.provider = provider or SyntheticProvider()
        self.registry = registry
        self.log = logger
        self.role_datasets = role_datasets or DEFAULT_ROLE_DATASETS

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
        plan_obj = self._stage(
            "plan",
            lambda: with_retry(planner_mod.plan, nl_request, client=self.llm,
                               retries=3, what="planner"),
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
        narrative = self._stage(
            "narrate",
            lambda: with_retry(narrator_mod.narrate, summary, client=self.llm,
                               retries=2, what="narrator"),
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
        for role in roles:
            cfg = self.role_datasets.get(role, {})
            req = plan_obj.model_copy(deep=True)
            if cfg.get("dataset"):
                req.dataset = cfg["dataset"]
            ep = build_execution_plan(req)
            if cfg.get("experiment"):
                ep.overrides.append(f"experiment_id={cfg['experiment']}")
                ep.config["experiment_id"] = cfg["experiment"]
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

    def _run_analytics(self, plan_obj: ReadyPlan, role_data: Dict) -> Dict:
        variable = plan_obj.variables[0]
        summary: Dict = {}
        for intent in plan_obj.analysis:
            for tool in self.registry.diagnostics_for_intent(intent):
                needed = tool.metadata.get("required_roles", [])
                if not all(r in role_data for r in needed):
                    self.log.warning("skip diagnostic '%s' (missing roles %s)",
                                     tool.name, [r for r in needed if r not in role_data])
                    continue
                kw = {"variable": variable}
                if tool.name in ("extremes", "etccdi_indices"):
                    kw["indices"] = plan_obj.indices or []
                try:
                    summary.update(tool.fn(role_data, **kw))
                except Exception as e:  # noqa: BLE001 — one diagnostic must not sink the rest
                    self.log.warning("diagnostic '%s' failed: %s", tool.name, e)
        return summary


# ===========================================================================
# End-to-end examples
# ===========================================================================

def _client():
    from openai import OpenAI
    return OpenAI(api_key="a80c8e8293a9826a995643cc46b9b015",
                  base_url="https://chat-ai.academiccloud.de/v1")


def _print(result: PipelineResult):
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


if __name__ == "__main__":
    assistant = ClimateAssistant(llm_client=_client(), provider=SyntheticProvider())

    print("\n>>> EXAMPLE 1 — BERLIN (point)")
    berlin = (
        "Describe the historical trend of daily maximum temperature observations "
        "over Berlin (lat=52.5200, lon=13.4050). Check whether CMIP models capture "
        "the broad climate change signal. Also describe plausible future climate "
        "conditions for Berlin."
    )
    _print(assistant.run(berlin))

    print("\n>>> EXAMPLE 2 — MÄRKISCH-ODERLAND (area)")
    mol = (
        "Describe historical and projected daily maximum temperature changes over "
        "the Märkisch-Oderland district using the extent lat_min=52.20, lat_max=53.55, "
        "lon_min=13.60, lon_max=15.25."
    )
    _print(assistant.run(mol))
