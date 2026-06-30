"""
ClimData Natural-Language Planner
=================================

Phase 2. Converts a natural-language climate question into a **validated,
structured request**, grounded in the Phase 1 capability registry
(``capabilities.py``).

Scope — this component ONLY plans and validates. It does not:
  * geocode / infer / estimate / invent coordinates
  * perform climate calculations or compute trends
  * generate Python code or ClimData overrides
  * narrate or analyse

Spatial contract
----------------
ClimData accepts ONLY machine coordinates — never a place name. A request is
spatially complete iff it carries EITHER:
  * a point:  lat AND lon
  * an extent: lat_min AND lat_max AND lon_min AND lon_max

If the user names a city / state / country / region without coordinates, the
planner MUST NOT guess them. It returns ``status = needs_clarification`` asking
for coordinates.

Output is one of two shapes, discriminated on ``status``:
  * ReadyPlan            (status="ready")
  * ClarificationRequest (status="needs_clarification")
"""

from __future__ import annotations

import json
from enum import Enum
from functools import lru_cache
from typing import List, Optional, Union

from pydantic import (
    BaseModel, Field, TypeAdapter, field_validator, model_validator,
)
from typing_extensions import Annotated, Literal

from capabilities import discover_all, CapabilityRegistry


# ===========================================================================
# Registry-derived allowed vocabularies (single source of truth = Phase 1)
# ===========================================================================

@lru_cache(maxsize=1)
def _registry() -> CapabilityRegistry:
    return discover_all()


@lru_cache(maxsize=1)
def vocab() -> dict:
    reg = _registry()
    datasets = {c.name: c for c in reg.datasets}
    return {
        "dataset_names": set(datasets),
        "dataset_vars": {k: set(v.metadata.get("variables", [])) for k, v in datasets.items()},
        "variables": {c.name for c in reg.variables},
        "indices": {c.name for c in reg.indices},
    }


# ===========================================================================
# Controlled vocabularies
# ===========================================================================

class Status(str, Enum):
    READY = "ready"
    NEEDS_CLARIFICATION = "needs_clarification"


class AnalysisIntent(str, Enum):
    """Recognised analysis intents (closed set)."""
    HISTORICAL_TREND = "historical_trend"
    CLIMATOLOGY = "climatology"
    FUTURE_PROJECTIONS = "future_projections"
    MODEL_EVALUATION = "model_evaluation"
    EXTREMES = "extremes"
    ETCCDI_INDICES = "etccdi_indices"
    DATASET_COMPARISON = "dataset_comparison"
    BIAS_CORRECTION = "bias_correction"


class ClarificationReason(str, Enum):
    MISSING_COORDINATES = "missing_coordinates"
    MISSING_VARIABLE = "missing_variable"
    UNSUPPORTED_VARIABLE = "unsupported_variable"
    UNSUPPORTED_DATASET = "unsupported_dataset"
    MISSING_ANALYSIS = "missing_analysis"


_MISSING_COORDS_MSG = (
    "ClimData requires either (1) lat/lon coordinates or (2) a lat/lon extent."
)


class TimeRange(BaseModel):
    start: Optional[str] = Field(None, description="ISO date, e.g. 1980-01-01")
    end: Optional[str] = Field(None, description="ISO date, e.g. 2020-12-31")
    model_config = {"extra": "forbid"}


# ===========================================================================
# 1. PLANNER SCHEMA  (the two structured output shapes)
# ===========================================================================

class ReadyPlan(BaseModel):
    """A spatially-complete, registry-valid request ready to hand downstream."""
    status: Literal[Status.READY] = Status.READY

    # --- point coordinates ---
    lat: Optional[float] = None
    lon: Optional[float] = None
    # --- OR bounding-box extent ---
    lat_min: Optional[float] = None
    lat_max: Optional[float] = None
    lon_min: Optional[float] = None
    lon_max: Optional[float] = None

    variable: Union[str, List[str]] = Field(..., description="Registry variable key(s).")
    dataset: Optional[str] = Field(None, description="Registry dataset key, if named.")
    time_range: Optional[TimeRange] = None
    analysis: List[AnalysisIntent] = Field(..., min_length=1)
    indices: Optional[List[str]] = Field(
        None, description="Registry index keys, when extremes/ETCCDI indices are requested."
    )

    model_config = {"use_enum_values": True, "extra": "forbid"}

    # ---- tolerate common LLM phrasings before validation -------------------
    @model_validator(mode="before")
    @classmethod
    def _flatten_spatial(cls, data):
        """Lift nested {extent:{...}} / {point:{...}} blocks to top-level fields."""
        if isinstance(data, dict):
            for key in ("extent", "bbox", "box", "point", "coordinates"):
                nested = data.get(key)
                if isinstance(nested, dict):
                    data = {**data, **nested}
                    data.pop(key, None)
        return data

    # ---- helpers -----------------------------------------------------------
    @property
    def variables(self) -> List[str]:
        return [self.variable] if isinstance(self.variable, str) else list(self.variable)

    # ---- 2. VALIDATION (registry-grounded; never invent) -------------------
    @field_validator("variable")
    @classmethod
    def _check_variable(cls, v):
        names = [v] if isinstance(v, str) else list(v)
        if not names:
            raise ValueError("At least one variable is required.")
        bad = [n for n in names if n not in vocab()["variables"]]
        if bad:
            raise ValueError(f"Unsupported variable(s) {bad}. Allowed: {sorted(vocab()['variables'])}")
        return v

    @field_validator("dataset")
    @classmethod
    def _check_dataset(cls, v):
        if v is not None and v not in vocab()["dataset_names"]:
            raise ValueError(f"Unknown dataset '{v}'. Allowed: {sorted(vocab()['dataset_names'])}")
        return v

    @field_validator("indices")
    @classmethod
    def _check_indices(cls, v):
        if v:
            bad = [i for i in v if i not in vocab()["indices"]]
            if bad:
                raise ValueError(f"Unsupported index/indices {bad}.")
        return v

    @model_validator(mode="after")
    def _check_spatial_and_consistency(self):
        point = [self.lat, self.lon]
        extent = [self.lat_min, self.lat_max, self.lon_min, self.lon_max]
        has_point = all(c is not None for c in point)
        any_point = any(c is not None for c in point)
        has_extent = all(c is not None for c in extent)
        any_extent = any(c is not None for c in extent)

        # exactly one fully-specified spatial mode
        if has_point and any_extent:
            raise ValueError("Provide a point OR an extent, not both.")
        if has_extent and any_point:
            raise ValueError("Provide a point OR an extent, not both.")
        if not has_point and not has_extent:
            raise ValueError(
                "ReadyPlan needs complete coordinates: lat+lon, or "
                "lat_min+lat_max+lon_min+lon_max. " + _MISSING_COORDS_MSG
            )
        if any_point and not has_point:
            raise ValueError("Point requires BOTH lat and lon.")
        if any_extent and not has_extent:
            raise ValueError("Extent requires all of lat_min, lat_max, lon_min, lon_max.")

        # range / ordering checks
        def _rng(name, val, lo, hi):
            if val is not None and not (lo <= val <= hi):
                raise ValueError(f"{name} {val} out of range [{lo}, {hi}].")

        for n, val in (("lat", self.lat), ("lat_min", self.lat_min), ("lat_max", self.lat_max)):
            _rng(n, val, -90.0, 90.0)
        for n, val in (("lon", self.lon), ("lon_min", self.lon_min), ("lon_max", self.lon_max)):
            _rng(n, val, -180.0, 180.0)
        if has_extent:
            if self.lat_min >= self.lat_max:
                raise ValueError("lat_min must be < lat_max.")
            if self.lon_min >= self.lon_max:
                raise ValueError("lon_min must be < lon_max.")

        # dataset must actually carry the requested variables
        if self.dataset:
            supported = vocab()["dataset_vars"].get(self.dataset, set())
            bad = [v for v in self.variables if v not in supported]
            if bad:
                raise ValueError(
                    f"Dataset '{self.dataset}' does not provide {bad}. Provides: {sorted(supported)}"
                )
        return self


class ClarificationRequest(BaseModel):
    """Returned when the request cannot be planned without more input."""
    status: Literal[Status.NEEDS_CLARIFICATION] = Status.NEEDS_CLARIFICATION
    reason: ClarificationReason
    message: str
    model_config = {"use_enum_values": True, "extra": "forbid"}


# Discriminated union — the planner always returns exactly one of these.
PlannerResponse = Annotated[
    Union[ReadyPlan, ClarificationRequest],
    Field(discriminator="status"),
]
_RESPONSE_ADAPTER: TypeAdapter = TypeAdapter(PlannerResponse)


# ===========================================================================
# 3. VALIDATION ENTRY POINT
# ===========================================================================

def validate_response(raw: Union[str, dict]) -> Union[ReadyPlan, ClarificationRequest]:
    """Parse + validate raw model output into a PlannerResponse (raises on violation)."""
    data = json.loads(raw) if isinstance(raw, str) else raw
    return _RESPONSE_ADAPTER.validate_python(data)


# ===========================================================================
# 4. PLANNER PROMPT + EXAMPLE IMPLEMENTATION
# ===========================================================================

def build_planner_prompt() -> str:
    reg = _registry()
    variable_lines = "\n".join(f"  - {c.name}: {c.description}" for c in reg.variables)
    datasets = ", ".join(c.name for c in reg.datasets)
    # ETCCDI / extreme index keys, grouped by family, for when indices are named
    fam: dict[str, list[str]] = {}
    for c in reg.indices:
        fam.setdefault(c.metadata.get("family", "other"), []).append(c.name)
    index_lines = "\n".join(f"  - {f}: {sorted(n)}" for f, n in sorted(fam.items()))
    intents = ", ".join(i.value for i in AnalysisIntent)

    return f"""You are the PLANNER for the ClimData climate toolkit. Convert the
user's request into ONE JSON object. Output JSON ONLY — no prose, no markdown.

You produce a structured request. You DO NOT compute anything, write code, or
generate ClimData overrides.

SPATIAL RULE (critical)
- ClimData accepts ONLY coordinates, never a place name.
- A request is spatially complete only if it has EITHER:
    point : lat AND lon
    extent: lat_min AND lat_max AND lon_min AND lon_max
- You MUST NOT geocode, infer, estimate, or invent coordinates.
- If the user gives a city/state/country/region name (e.g. "Berlin", "Germany")
  WITHOUT numeric coordinates, return:
    {{"status":"needs_clarification","reason":"missing_coordinates",
      "message":"{_MISSING_COORDS_MSG}"}}

WHEN COORDINATES ARE PRESENT, return (put ALL coordinate fields at the TOP
LEVEL — never nest them under "extent", "point", or "coordinates"):
  point :  {{"status":"ready","lat":<num>,"lon":<num>, ...}}
  extent:  {{"status":"ready","lat_min":<num>,"lat_max":<num>,
            "lon_min":<num>,"lon_max":<num>, ...}}
  full shape:
  {{"status":"ready",
    "lat":<num>,"lon":<num>,                # OR lat_min/lat_max/lon_min/lon_max
    "variable":"<key>" | ["<key>",...],
    "dataset":<KEY|null>,                    # only if the user named a source
    "time_range":{{"start":<ISO|null>,"end":<ISO|null>}}|null,
    "analysis":["<intent>",...],
    "indices":["<index_key>",...]|null}}     # only for extremes / ETCCDI

RULES
- Use ONLY the variables, datasets, indices and intents listed below. NEVER
  invent any of them.
- Map words to variable keys: "max temperature"->tasmax, "min temperature"->tasmin,
  "temperature"->tas, "rainfall"/"precipitation"->pr, "humidity"->hurs,
  "wind"->sfcWind, "solar"/"radiation"->rsds.
- Set "dataset": null unless the user names a source/coverage that pins it.
- "indices" is only filled when the intent is extremes or etccdi_indices and the
  user (or the index family) identifies specific indices.

ANALYSIS INTENTS
  {intents}

VARIABLES
{variable_lines}

DATASETS
  {datasets}

EXTREME / ETCCDI INDICES (by family)
{index_lines}

Return only the JSON object."""


def plan(
    nl_request: str,
    client=None,
    model: str = "qwen3-coder-30b-a3b-instruct",
    max_retries: int = 2,
) -> Union[ReadyPlan, ClarificationRequest]:
    """
    NL → LLM → validated PlannerResponse.

    `client` is an OpenAI-compatible client (e.g. the SAIA client). On a
    validation failure the error is fed back so the model can repair its JSON;
    the repair passes through the SAME validator, so nothing unsupported and no
    invented coordinate can slip through. Use validate_response() for offline
    testing without a client.
    """
    if client is None:
        raise ValueError("Provide an OpenAI-compatible `client` (e.g. the SAIA client).")

    messages = [
        {"role": "system", "content": build_planner_prompt()},
        {"role": "user", "content": nl_request},
    ]
    last_err: Optional[Exception] = None
    for _ in range(max_retries + 1):
        resp = client.chat.completions.create(
            model=model, messages=messages, temperature=0,
            response_format={"type": "json_object"},
        )
        content = resp.choices[0].message.content
        try:
            return validate_response(content)
        except Exception as e:
            last_err = e
            messages += [
                {"role": "assistant", "content": content},
                {"role": "user", "content":
                    f"That JSON failed validation:\n{e}\nReturn corrected JSON only."},
            ]
    raise ValueError(f"Planner could not produce a valid response: {last_err}")


if __name__ == "__main__":
    # Offline smoke test of the two canonical shapes from the spec.
    ready_point = validate_response({
        "status": "ready", "lat": 52.5, "lon": 13.4,
        "variable": "tasmax", "analysis": ["historical_trend"],
    })
    ready_extent = validate_response({
        "status": "ready", "lat_min": 52.0, "lat_max": 53.0,
        "lon_min": 13.0, "lon_max": 14.0,
        "variable": "tasmax", "analysis": ["historical_trend"],
    })
    clarify = validate_response({
        "status": "needs_clarification", "reason": "missing_coordinates",
        "message": _MISSING_COORDS_MSG,
    })
    for label, r in [("point", ready_point), ("extent", ready_extent), ("clarify", clarify)]:
        print(f"[OK] {label}: {r.model_dump_json(exclude_none=True)}")
