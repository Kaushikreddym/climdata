"""Extensible tool registry — the plugin backbone of the orchestration layer.

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
    """One registered capability, addressed by ``(kind, name)``.

    Attributes:
        name (str): Unique name within its kind, e.g. ``"historical_trend"``.
        kind (str): Category — ``"diagnostic"``, ``"provider"``, and so on.
        fn (Callable | None): The implementation. ``None`` for entries that are
            metadata only.
        description (str): One-line summary, surfaced to the planning LLM.
        metadata (dict): Kind-specific extras. Diagnostics carry ``intent``,
            ``required_roles`` and ``provided_by``.
    """

    name: str
    kind: str                       # "diagnostic" | "provider" | ...
    fn: Optional[Callable]
    description: str = ""
    metadata: Dict = field(default_factory=dict)


class ToolRegistry:
    """A name-addressed table of the capabilities the orchestrator can dispatch.

    Keyed by ``(kind, name)``, so the same name may exist under several kinds.
    Registration is last-write-wins: adding a tool with an existing key replaces
    it, which is how a built-in diagnostic can be overridden.

    The module-level :data:`REGISTRY` is the instance everything shares; there is
    rarely a reason to construct another.

    Example:
        >>> from climdata.LLM.tool_registry import REGISTRY
        >>> @REGISTRY.diagnostic("cold_days", intent="extremes",       # doctest: +SKIP
        ...                      required_roles=["observation"])
        ... def cold_days(roles, variable=None):
        ...     return {"cold_days": 12}
    """

    def __init__(self) -> None:
        """Create an empty registry."""
        self._tools: Dict[tuple, Tool] = {}

    # ---- generic registration --------------------------------------------
    def add(self, kind: str, name: str, fn: Optional[Callable] = None,
            description: str = "", **metadata) -> Tool:
        """Register a tool directly, replacing any entry with the same key.

        The imperative counterpart to :meth:`register`, for wiring an existing
        function under a second name — as ``etccdi_indices`` reuses the
        ``extremes`` adapter.

        Args:
            kind (str): Category, e.g. ``"diagnostic"``.
            name (str): Name, unique within ``kind``.
            fn (Callable, optional): The implementation.
            description (str): One-line summary.
            **metadata: Kind-specific extras stored on the tool.

        Returns:
            Tool: The registered tool.
        """
        tool = Tool(name=name, kind=kind, fn=fn, description=description, metadata=metadata)
        self._tools[(kind, name)] = tool
        return tool

    def register(self, kind: str, name: str, description: str = "", **metadata):
        """Return a decorator that registers the function it wraps.

        Args:
            kind (str): Category, e.g. ``"diagnostic"``.
            name (str): Name, unique within ``kind``.
            description (str): One-line summary.
            **metadata: Kind-specific extras stored on the tool.

        Returns:
            Callable: A decorator returning its argument unchanged, so the
            function stays directly callable.
        """
        def deco(fn: Callable) -> Callable:
            self.add(kind, name, fn, description, **metadata)
            return fn
        return deco

    # ---- convenience: diagnostics ----------------------------------------
    def diagnostic(self, name: str, *, intent: str, required_roles: List[str],
                   provided_by: str = "analytics", description: str = ""):
        """Register a diagnostic, the kind the orchestrator dispatches on intent.

        A diagnostic takes ``(roles, variable=None)`` — where ``roles`` maps role
        names to extracted xarray data — and returns a flat dict of summary
        values. That uniform contract is what lets the orchestrator stay generic.

        Args:
            name (str): Diagnostic name.
            intent (str): Analysis intent it serves, matched against the planner's
                output. Keyword-only.
            required_roles (list[str]): Role names that must be present in
                ``roles``, e.g. ``["observation", "model_historical"]``. The
                orchestrator skips the diagnostic when any is missing.
                Keyword-only.
            provided_by (str): Which subsystem implements it — ``"analytics"`` or
                ``"climdata"``. Keyword-only. Defaults to ``"analytics"``.
            description (str): One-line summary. Keyword-only.

        Returns:
            Callable: A decorator for the diagnostic function.
        """
        return self.register("diagnostic", name, description=description,
                             intent=intent, required_roles=required_roles,
                             provided_by=provided_by)

    # ---- lookup ----------------------------------------------------------
    def get(self, kind: str, name: str) -> Tool:
        """Look up one tool.

        Args:
            kind (str): Category.
            name (str): Name within that category.

        Returns:
            Tool: The registered tool.

        Raises:
            KeyError: If no tool matches.
        """
        if (kind, name) not in self._tools:
            raise KeyError(f"No {kind} tool named '{name}'.")
        return self._tools[(kind, name)]

    def has(self, kind: str, name: str) -> bool:
        """Test whether a tool is registered.

        Args:
            kind (str): Category.
            name (str): Name within that category.

        Returns:
            bool: ``True`` if present.
        """
        return (kind, name) in self._tools

    def list(self, kind: Optional[str] = None) -> List[Tool]:
        """List registered tools, optionally of one kind.

        Args:
            kind (str, optional): Category to filter by. ``None`` returns
                everything.

        Returns:
            list[Tool]: Matching tools, in registration order.
        """
        return [t for (k, _), t in self._tools.items() if kind is None or k == kind]

    def diagnostics_for_intent(self, intent: str) -> List[Tool]:
        """List the diagnostics that serve one analysis intent.

        Args:
            intent (str): Intent name, e.g. ``"historical_trend"``.

        Returns:
            list[Tool]: Matching diagnostics, empty if the intent is unserved.
        """
        return [t for t in self.list("diagnostic") if t.metadata.get("intent") == intent]

    # ---- bridge to Phase 1 capability registry ---------------------------
    def datasets(self) -> List[str]:
        """List the dataset providers the capability layer discovered.

        Returns:
            list[str]: Provider names, from ``conf/mappings/parameters.yaml``.
        """
        return [c.name for c in discover_all().datasets]

    def indices(self) -> List[str]:
        """List the extreme indices the capability layer discovered.

        Returns:
            list[str]: Index names, from ``conf/mappings/indices.yaml``.
        """
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
    """Report annual and summer linear trends in the observed series.

    Args:
        roles (dict): Must contain ``"observation"``.
        variable (str, optional): Variable to analyse when the role holds a
            Dataset with more than one.

    Returns:
        dict: ``annual_trend`` and ``summer_trend``, each in units per year.
    """
    da = roles["observation"]
    return {
        "annual_trend": analytics.annual_trend(da, variable),
        "summer_trend": analytics.seasonal_trend(da, "JJA", variable),
    }


@REGISTRY.diagnostic("climatology", intent="climatology",
                     required_roles=["observation"],
                     description="Long-term mean and monthly climatology.")
def _diag_climatology(roles, variable=None) -> Dict:
    """Report the long-term mean and the monthly climatology.

    Args:
        roles (dict): Must contain ``"observation"``.
        variable (str, optional): Variable to analyse.

    Returns:
        dict: ``overall_mean`` and ``monthly_mean``.
    """
    clim = analytics.climatology(roles["observation"], variable)
    return {"overall_mean": clim["overall_mean"], "monthly_mean": clim["monthly_mean"]}


@REGISTRY.diagnostic("extremes", intent="extremes",
                     required_roles=["observation"], provided_by="climdata",
                     description="ETCCDI/threshold extremes (reuses ClimData calc_index).")
def _diag_extremes(roles, variable=None, indices=None) -> Dict:
    """Compute ETCCDI and threshold extremes for the observed series.

    Args:
        roles (dict): Must contain ``"observation"``.
        variable (str, optional): Variable to analyse.
        indices (list[str], optional): Specific indices to compute. Defaults to
            those appropriate for the variable.

    Returns:
        dict: One entry per computed index, plus notes for any that could not be
        computed.
    """
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
    """Score a historical simulation against observations.

    Args:
        roles (dict): Must contain ``"observation"`` and ``"model_historical"``.
        variable (str, optional): Variable to analyse.

    Returns:
        dict: Bias, error and correlation metrics.
    """
    return analytics.evaluation_metrics(roles["model_historical"], roles["observation"], variable)


@REGISTRY.diagnostic("future_projections", intent="future_projections",
                     required_roles=["model_future"],
                     description="Projected trend in the future-scenario simulation.")
def _diag_future_projections(roles, variable=None) -> Dict:
    """Report trends in a future-scenario simulation, and its change from historical.

    Args:
        roles (dict): Must contain ``"model_future"``. When ``"model_historical"``
            is also present, the projected mean change is added.
        variable (str, optional): Variable to analyse.

    Returns:
        dict: ``future_annual_trend`` and ``future_summer_trend``, plus
        ``projected_mean_change`` when a historical run is available.
    """
    fut = roles["model_future"]
    out = {
        "future_annual_trend": analytics.annual_trend(fut, variable),
        "future_summer_trend": analytics.seasonal_trend(fut, "JJA", variable),
    }
    # projected change = future mean minus historical-model mean (no new math beyond means)
    if "model_historical" in roles:
        f = analytics._reduce_space(analytics._as_dataarray(fut, variable))
        h = analytics._reduce_space(analytics._as_dataarray(roles["model_historical"], variable))
        out["projected_mean_change"] = round(float(f.mean() - h.mean()), 4)
    return out
