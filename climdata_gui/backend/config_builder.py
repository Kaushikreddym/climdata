"""Build Hydra override lists from the GUI state."""

from __future__ import annotations

from datetime import date
from typing import List, Optional, Sequence


def _fmt(d: date) -> str:
    return d.strftime("%Y-%m-%d")


def build_overrides(
    dataset: str,
    lat: Optional[float],
    lon: Optional[float],
    start_date: Optional[date],
    end_date: Optional[date],
    variables: Optional[Sequence[str]] = None,
) -> List[str]:
    """Construct a list of Hydra overrides for ``ClimData(overrides=...)``.

    Args:
        dataset: Provider name (e.g. ``"mswx"``, ``"cmip"``, ``"dwd"``).
        lat: Latitude of the point of interest.
        lon: Longitude of the point of interest.
        start_date: Start of the requested time range.
        end_date: End of the requested time range.
        variables: Optional list of variables; if omitted the defaults
            from ``conf/config.yaml`` are used.
    """
    overrides: List[str] = [f"dataset={dataset}"]

    if lat is not None:
        overrides.append(f"lat={lat}")
    if lon is not None:
        overrides.append(f"lon={lon}")

    if start_date is not None:
        overrides.append(f'time_range.start_date="{_fmt(start_date)}"')
    if end_date is not None:
        overrides.append(f'time_range.end_date="{_fmt(end_date)}"')

    if variables:
        vars_str = "[" + ",".join(variables) + "]"
        overrides.append(f"variables={vars_str}")

    return overrides
