"""Thin wrapper around ``climdata.ClimData`` for GUI consumption."""

from __future__ import annotations

from typing import Any, List, Optional, Sequence


def run_pipeline(overrides: Sequence[str], seq: Optional[Sequence[str]] = None) -> Any:
    """Run a climdata workflow with the given Hydra overrides.

    Args:
        overrides: Hydra override strings, e.g. ``["dataset=mswx", "lat=52.5"]``.
        seq: Ordered list of workflow actions. Defaults to
            ``["extract", "to_csv"]``.

    Returns:
        The ``WorkflowResult`` from ``ClimData.run_workflow``.

    Notes:
        Imports ``climdata`` lazily so the GUI module can be imported (and
        the scaffold inspected) even when the scientific stack is not yet
        installed.
    """
    from climdata import ClimData  # lazy import

    actions: List[str] = list(seq) if seq else ["extract", "to_csv"]
    extractor = ClimData(overrides=list(overrides))
    return extractor.run_workflow(actions=actions)
