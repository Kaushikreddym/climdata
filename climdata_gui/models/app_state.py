"""Shared GUI state container."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import date
from typing import Dict, List, Optional, Tuple


@dataclass
class AppState:
    """Holds user-selected inputs that flow into a climdata pipeline run."""

    dataset: str = "mswx"
    lat: Optional[float] = None
    lon: Optional[float] = None
    box: Optional[dict] = None          # {lat_min, lat_max, lon_min, lon_max}
    start_date: Optional[date] = None
    end_date: Optional[date] = None
    variables: List[str] = field(default_factory=lambda: ["tasmin", "tasmax", "pr"])
    actions: List[str] = field(default_factory=lambda: ["extract", "to_csv"])
    extra_overrides: List[str] = field(default_factory=list)
    data_dir: Optional[str] = None
    experiment_id: Optional[str] = None   # CMIP6 scenario, e.g. "ssp585"
    source_id: Optional[str] = None       # CMIP6 model, e.g. "MIROC6"

    # CMIP6 pickers outside the Download tab, keyed by slot name
    # ("basd-target", "cmp-a", "cmp-b"): {"experiment_id": ..., "source_id": ...}
    cmip_slots: Dict[str, Dict[str, Optional[str]]] = field(default_factory=dict)

    def set_cmip(self, slot: str, experiment_id: Optional[str] = None,
                 source_id: Optional[str] = None) -> None:
        """Record a slot's CMIP6 selection (``None`` leaves a value untouched)."""
        if slot == "download":
            if experiment_id is not None:
                self.experiment_id = experiment_id or None
            if source_id is not None:
                self.source_id = source_id or None
            return
        entry = self.cmip_slots.setdefault(
            slot, {"experiment_id": None, "source_id": None})
        if experiment_id is not None:
            entry["experiment_id"] = experiment_id or None
        if source_id is not None:
            entry["source_id"] = source_id or None

    def get_cmip(self, slot: str) -> Tuple[Optional[str], Optional[str]]:
        """Return ``(experiment_id, source_id)`` for a slot."""
        if slot == "download":
            return self.experiment_id, self.source_id
        entry = self.cmip_slots.get(slot, {})
        return entry.get("experiment_id"), entry.get("source_id")

    def clear_cmip(self, slot: str) -> None:
        """Forget a slot's selection (its dataset is no longer CMIP-driven)."""
        if slot == "download":
            self.experiment_id = None
            self.source_id = None
        else:
            self.cmip_slots.pop(slot, None)

    def is_ready(self) -> bool:
        """True once the minimum inputs for a run are populated."""
        has_spatial = (
            (self.lat is not None and self.lon is not None)
            or self.box is not None
        )
        return (
            self.dataset is not None
            and has_spatial
            and self.start_date is not None
            and self.end_date is not None
        )
