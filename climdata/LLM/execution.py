"""
ClimData Execution Layer
========================

Phase 3a. Turns a validated planner request (``ReadyPlan`` from ``planner.py``)
into a **reproducible ClimData execution plan**: the exact Hydra overrides and
ordered workflow steps needed to run the extraction, plus provenance.

Responsibilities:
  * Convert planner output into valid ClimData overrides (Hydra dot-notation).
  * Validate parameters against the Phase 1 capability registry.
  * Handle missing capabilities gracefully (status="error" + structured report,
    never a crash / never an invented capability).
  * Generate reproducible extraction configurations (content hash + version +
    full resolved config snapshot).

This layer does NOT run ClimData, compute anything, or narrate. It produces a
plan that a runner (or the analytics layer) can execute deterministically.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from enum import Enum
from functools import lru_cache
from typing import Dict, List, Optional, Union

from pydantic import BaseModel, Field

from capabilities import discover_all
from planner import ReadyPlan, AnalysisIntent

# The planner's ReadyPlan IS the execution layer's input contract.
ClimateRequest = ReadyPlan


# ===========================================================================
# Registry helpers
# ===========================================================================

@lru_cache(maxsize=1)
def _caps() -> dict:
    reg = discover_all()
    datasets = {c.name: c for c in reg.datasets}
    return {
        "dataset_names": list(datasets),
        "dataset_vars": {k: set(v.metadata.get("variables", [])) for k, v in datasets.items()},
        "variables": {c.name for c in reg.variables},
        "indices": {c.name for c in reg.indices},
        "index_inputs": {c.name: set(c.metadata.get("input_variables", [])) for c in reg.indices},
        "version": reg.version,
    }


# Preference order when the planner did not pin a dataset. Global, well-supported
# observation/reanalysis sources first; impact/model sources last.
DEFAULT_DATASET_PRIORITY = [
    "ERA5", "MSWX", "W5E5", "HYRAS", "POWER", "DWD",
    "NEXGDDP", "CMIP_W5E5", "CMIP", "HOSTRADA", "AGRI_ISIMIP",
]

# Analysis intents ClimData can satisfy natively vs. those deferred to analytics.
_CLIMDATA_NATIVE_INTENTS = {AnalysisIntent.EXTREMES.value, AnalysisIntent.ETCCDI_INDICES.value}


# ===========================================================================
# 1. EXECUTION SCHEMA
# ===========================================================================

class ExecutionStatus(str, Enum):
    READY = "ready"
    ERROR = "error"


class MissingCapability(BaseModel):
    kind: str                    # "dataset" | "variable" | "index" | "intent"
    name: str
    detail: str


class WorkflowStep(BaseModel):
    """One ordered action in the ClimData pipeline."""
    action: str                  # extract | calc_index | analyze
    params: Dict = Field(default_factory=dict)
    provided_by: str             # "climdata" | "analytics"


class ExecutionPlan(BaseModel):
    status: ExecutionStatus
    dataset: Optional[str] = None
    variables: List[str] = Field(default_factory=list)
    spatial: Dict = Field(default_factory=dict)          # {mode: point|box, ...}
    time_range: Optional[Dict] = None
    overrides: List[str] = Field(default_factory=list)   # reproducible Hydra overrides
    steps: List[WorkflowStep] = Field(default_factory=list)
    config: Dict = Field(default_factory=dict)           # full resolved config snapshot
    missing: List[MissingCapability] = Field(default_factory=list)
    provenance: Dict = Field(default_factory=dict)

    model_config = {"use_enum_values": True}

    def climdata_call(self) -> str:
        """Reproducible one-liner to recreate the extraction."""
        return f"ClimData(overrides={self.overrides!r})"


# ===========================================================================
# 2 + 3. CONFIG BUILDER  +  VALIDATION LOGIC
# ===========================================================================

def _resolve_dataset(
    variables: List[str], requested: Optional[str], missing: List[MissingCapability]
) -> Optional[str]:
    """
    Pick a dataset that provides every requested variable.

    If the planner pinned one, validate it. Otherwise choose by priority. Records
    a MissingCapability (and returns None) instead of raising when nothing fits.
    """
    caps = _caps()
    need = set(variables)

    if requested is not None:
        if requested not in caps["dataset_names"]:
            missing.append(MissingCapability(
                kind="dataset", name=requested,
                detail=f"Unknown dataset. Known: {caps['dataset_names']}"))
            return None
        unsupported = need - caps["dataset_vars"][requested]
        if unsupported:
            missing.append(MissingCapability(
                kind="dataset", name=requested,
                detail=f"Does not provide {sorted(unsupported)}."))
            return None
        return requested

    # auto-resolve by priority
    for cand in DEFAULT_DATASET_PRIORITY:
        if cand in caps["dataset_vars"] and need.issubset(caps["dataset_vars"][cand]):
            return cand
    missing.append(MissingCapability(
        kind="dataset", name="<auto>",
        detail=f"No single dataset provides all of {sorted(need)}."))
    return None


def _spatial(req: ClimateRequest) -> Dict:
    if req.lat is not None and req.lon is not None:
        return {"mode": "point", "lat": req.lat, "lon": req.lon}
    return {
        "mode": "box",
        "lat_min": req.lat_min, "lat_max": req.lat_max,
        "lon_min": req.lon_min, "lon_max": req.lon_max,
    }


def _overrides(dataset: Optional[str], variables: List[str], spatial: Dict,
               time_range: Optional[Dict]) -> List[str]:
    ov: List[str] = []
    if dataset:
        ov.append(f"dataset={dataset}")
    ov.append(f"variables=[{','.join(variables)}]")
    if spatial["mode"] == "point":
        ov += [f"lat={spatial['lat']}", f"lon={spatial['lon']}"]
    else:
        ov += [
            f"box.lat_min={spatial['lat_min']}", f"box.lat_max={spatial['lat_max']}",
            f"box.lon_min={spatial['lon_min']}", f"box.lon_max={spatial['lon_max']}",
        ]
    if time_range:
        if time_range.get("start"):
            ov.append(f"time_range.start_date={time_range['start']}")
        if time_range.get("end"):
            ov.append(f"time_range.end_date={time_range['end']}")
    return ov


def _validate_indices(req: ClimateRequest, variables: List[str],
                      missing: List[MissingCapability]) -> List[str]:
    """Validate requested indices exist and their input variables are covered."""
    caps = _caps()
    valid: List[str] = []
    for idx in (req.indices or []):
        if idx not in caps["indices"]:
            missing.append(MissingCapability(
                kind="index", name=idx, detail="Not a registry index."))
            continue
        needed = caps["index_inputs"].get(idx, set())
        gap = needed - set(variables)
        if gap:
            missing.append(MissingCapability(
                kind="index", name=idx,
                detail=f"Needs variables {sorted(needed)}; missing {sorted(gap)}."))
            continue
        valid.append(idx)
    return valid


def build_execution_plan(request: Union[ClimateRequest, dict]) -> ExecutionPlan:
    """Convert a validated planner request into a reproducible ExecutionPlan."""
    req = ReadyPlan.model_validate(request) if isinstance(request, dict) else request
    caps = _caps()
    missing: List[MissingCapability] = []

    variables = req.variables
    dataset = _resolve_dataset(variables, req.dataset, missing)
    spatial = _spatial(req)
    time_range = (
        {"start": req.time_range.start, "end": req.time_range.end}
        if req.time_range else None
    )

    # indices (only meaningful for the extremes / ETCCDI intents)
    wants_index = any(a in _CLIMDATA_NATIVE_INTENTS for a in req.analysis)
    valid_indices = _validate_indices(req, variables, missing) if (req.indices or wants_index) else []
    if wants_index and not valid_indices and not req.indices:
        missing.append(MissingCapability(
            kind="index", name="<none>",
            detail="extremes/etccdi_indices intent requires at least one valid index."))

    overrides = _overrides(dataset, variables, spatial, time_range)

    # ----- ordered workflow steps -----
    steps: List[WorkflowStep] = [
        WorkflowStep(action="extract", provided_by="climdata",
                     params={"overrides": overrides})
    ]
    for idx in valid_indices:
        steps.append(WorkflowStep(action="calc_index", provided_by="climdata",
                                  params={"index": idx}))
    # intents ClimData doesn't compute natively → analytics layer
    for intent in req.analysis:
        if intent not in _CLIMDATA_NATIVE_INTENTS:
            steps.append(WorkflowStep(action="analyze", provided_by="analytics",
                                      params={"intent": intent}))

    # ----- resolved config snapshot (reproducible) -----
    config = {
        "dataset": dataset,
        "variables": variables,
        "spatial": spatial,
        "time_range": time_range,
        "indices": valid_indices,
        "analysis": list(req.analysis),
    }

    status = ExecutionStatus.ERROR if missing else ExecutionStatus.READY
    plan_hash = hashlib.sha256(
        json.dumps({"overrides": overrides, "config": config}, sort_keys=True,
                   default=str).encode()
    ).hexdigest()[:16]

    return ExecutionPlan(
        status=status,
        dataset=dataset,
        variables=variables,
        spatial=spatial,
        time_range=time_range,
        overrides=overrides,
        steps=steps,
        config=config,
        missing=missing,
        provenance={
            "climdata_version": caps["version"],
            "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "plan_hash": plan_hash,
        },
    )


# ===========================================================================
# 4. EXAMPLE — Berlin request
# ===========================================================================

if __name__ == "__main__":
    # A planner ReadyPlan for Berlin (coordinates supplied by the user, NOT geocoded).
    berlin = ReadyPlan.model_validate({
        "status": "ready",
        "lat": 52.52, "lon": 13.405,
        "variable": "tasmax",
        "dataset": "HYRAS",
        "time_range": {"start": "1980-01-01", "end": "2020-12-31"},
        "analysis": ["historical_trend"],
    })
    plan = build_execution_plan(berlin)
    print(plan.model_dump_json(indent=2, exclude_none=True))
    print("\nReproduce with:", plan.climdata_call())
