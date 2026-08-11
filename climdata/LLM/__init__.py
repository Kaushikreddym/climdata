"""
ClimData Climate-Assistant (LLM layer)
======================================

Turns a natural-language climate question into a grounded written assessment:

    Planner  →  Execution  →  ClimData  →  Analytics  →  Narrator

Each stage validates its own output, and the pipeline is fail-soft: a failed
stage degrades the assessment rather than the process.

Typical use::

    from climdata.LLM import ClimateAssistant, RealClimDataProvider, build_client

    assistant = ClimateAssistant(
        llm_client=build_client(),               # reads SAIA_API_KEY
        provider=RealClimDataProvider(),         # or SyntheticProvider() offline
    )
    result = assistant.run("Trend in daily max temperature at lat=52.52, lon=13.405?")
    print(result.narrative)

Each module also runs standalone as a demo, e.g.
``python -m climdata.LLM.planner``.
"""

from .analytics import analyze, capability_report, variable_units
from .capabilities import CapabilityRegistry, discover_all
from .execution import ExecutionPlan, ExecutionStatus, build_execution_plan
from .narrator import ClimateSummary, NarrativeAssessment, grounding_check, narrate
from .orchestrator import (
    DEFAULT_MODEL,
    ClimateAssistant,
    PipelineResult,
    RealClimDataProvider,
    SyntheticProvider,
    build_client,
    print_result,
)
from .planner import ClarificationRequest, ReadyPlan, plan, validate_response
from .tool_registry import REGISTRY, ToolRegistry

__all__ = [
    "ClimateAssistant",
    "PipelineResult",
    "RealClimDataProvider",
    "SyntheticProvider",
    "build_client",
    "print_result",
    "DEFAULT_MODEL",
    "plan",
    "ReadyPlan",
    "ClarificationRequest",
    "validate_response",
    "build_execution_plan",
    "ExecutionPlan",
    "ExecutionStatus",
    "narrate",
    "NarrativeAssessment",
    "ClimateSummary",
    "grounding_check",
    "analyze",
    "capability_report",
    "variable_units",
    "discover_all",
    "CapabilityRegistry",
    "REGISTRY",
    "ToolRegistry",
]
