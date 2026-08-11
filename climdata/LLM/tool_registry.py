"""
Extensible Tool Registry
========================

The plugin backbone of the orchestration layer. Everything the orchestrator
dispatches — diagnostics, dataset-role resolvers, and the bridge to the Phase 1
capability registry (datasets + indices) — is looked up here by name.

New datasets and indices are already auto-discovered from ClimData's YAML config
(Phase 1). New *analytics* (diagnostics) plug in by decorating a function:

    from climdata.LLM.tool_registry import REGISTRY

    @REGISTRY.diagnostic("cold_days_change", intent="extremes",
                         required_roles=["observation"], provided_by="analytics")
    def cold_days_change(roles, variable=None):
        ...
        return {"cold_days_change": value}

A diagnostic takes ``roles`` (a dict of role-name -> extracted xarray data) plus
the target ``variable`` and returns a flat dict of summary keys. This uniform
contract is what lets the orchestrator stay generic.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional

from . import analytics
from .capabilities import discover_all


# ===========================================================================
# Tool model + registry
# ===========================================================================

@dataclass
class Tool:
    name: str
    kind: str                       # "diagnostic" | "provider" | ...
    fn: Optional[Callable]
    description: str = ""
    metadata: Dict = field(default_factory=dict)


class ToolRegistry:
    def __init__(self) -> None:
        self._tools: Dict[tuple, Tool] = {}

    # ---- generic registration --------------------------------------------
    def add(self, kind: str, name: str, fn: Optional[Callable] = None,
            description: str = "", **metadata) -> Tool:
        tool = Tool(name=name, kind=kind, fn=fn, description=description, metadata=metadata)
        self._tools[(kind, name)] = tool
        return tool

    def register(self, kind: str, name: str, description: str = "", **metadata):
        def deco(fn: Callable) -> Callable:
            self.add(kind, name, fn, description, **metadata)
            return fn
        return deco

    # ---- convenience: diagnostics ----------------------------------------
    def diagnostic(self, name: str, *, intent: str, required_roles: List[str],
                   provided_by: str = "analytics", description: str = ""):
        return self.register("diagnostic", name, description=description,
                             intent=intent, required_roles=required_roles,
                             provided_by=provided_by)

    # ---- lookup ----------------------------------------------------------
    def get(self, kind: str, name: str) -> Tool:
        if (kind, name) not in self._tools:
            raise KeyError(f"No {kind} tool named '{name}'.")
        return self._tools[(kind, name)]

    def has(self, kind: str, name: str) -> bool:
        return (kind, name) in self._tools

    def list(self, kind: Optional[str] = None) -> List[Tool]:
        return [t for (k, _), t in self._tools.items() if kind is None or k == kind]

    def diagnostics_for_intent(self, intent: str) -> List[Tool]:
        return [t for t in self.list("diagnostic") if t.metadata.get("intent") == intent]

    # ---- bridge to Phase 1 capability registry ---------------------------
    def datasets(self) -> List[str]:
        return [c.name for c in discover_all().datasets]

    def indices(self) -> List[str]:
        return [c.name for c in discover_all().indices]


REGISTRY = ToolRegistry()


# ===========================================================================
# Built-in diagnostics  (uniform contract: roles dict -> flat summary dict)
# ===========================================================================
# Each adapter pulls the role(s) it needs and reuses analytics.py functions.
# `provided_by="climdata"` marks diagnostics that reuse ClimData-native code.

@REGISTRY.diagnostic("historical_trend", intent="historical_trend",
                     required_roles=["observation"],
                     description="Annual and summer linear trends of the observed series.")
def _diag_historical_trend(roles, variable=None) -> Dict:
    da = roles["observation"]
    return {
        "annual_trend": analytics.annual_trend(da, variable),
        "summer_trend": analytics.seasonal_trend(da, "JJA", variable),
    }


@REGISTRY.diagnostic("climatology", intent="climatology",
                     required_roles=["observation"],
                     description="Long-term mean and monthly climatology.")
def _diag_climatology(roles, variable=None) -> Dict:
    clim = analytics.climatology(roles["observation"], variable)
    return {"overall_mean": clim["overall_mean"], "monthly_mean": clim["monthly_mean"]}


@REGISTRY.diagnostic("extremes", intent="extremes",
                     required_roles=["observation"], provided_by="climdata",
                     description="ETCCDI/threshold extremes (reuses ClimData calc_index).")
def _diag_extremes(roles, variable=None, indices=None) -> Dict:
    # Delegates to analytics.analyze so the threshold/index rules — including
    # which variables a fixed-degree threshold applies to, and the notes emitted
    # for indices that could not be computed — live in exactly one place.
    return analytics.analyze(roles["observation"], intents=["extremes"],
                             variable=variable, indices=indices)


# etccdi_indices shares the extremes adapter
REGISTRY.add("diagnostic", "etccdi_indices", _diag_extremes,
             intent="etccdi_indices", required_roles=["observation"],
             provided_by="climdata")


@REGISTRY.diagnostic("model_evaluation", intent="model_evaluation",
                     required_roles=["observation", "model_historical"],
                     description="Model-vs-observation bias / error / correlation.")
def _diag_model_evaluation(roles, variable=None) -> Dict:
    return analytics.evaluation_metrics(roles["model_historical"], roles["observation"], variable)


@REGISTRY.diagnostic("future_projections", intent="future_projections",
                     required_roles=["model_future"],
                     description="Projected trend in the future-scenario simulation.")
def _diag_future_projections(roles, variable=None) -> Dict:
    fut = roles["model_future"]
    out = {
        "future_annual_trend": analytics.annual_trend(fut, variable),
        "future_summer_trend": analytics.seasonal_trend(fut, "JJA", variable),
    }
    # projected change = future mean minus historical-model mean (no new math beyond means)
    if "model_historical" in roles:
        import numpy as np
        f = analytics._reduce_space(analytics._as_dataarray(fut, variable))
        h = analytics._reduce_space(analytics._as_dataarray(roles["model_historical"], variable))
        out["projected_mean_change"] = round(float(f.mean() - h.mean()), 4)
    return out
