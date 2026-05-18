"""Shared GUI state container."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import date
from typing import List, Optional


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
