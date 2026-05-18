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
    extra_overrides: Optional[Sequence[str]] = None,
    box: Optional[dict] = None,
    data_dir: Optional[str] = None,
) -> List[str]:
    """Construct a list of Hydra overrides for ``ClimData(overrides=...)``.

    Args:
        dataset: Provider name (e.g. ``"mswx"``, ``"cmip"``, ``"dwd"``).
        lat: Latitude of the point of interest (mutually exclusive with box).
        lon: Longitude of the point of interest (mutually exclusive with box).
        start_date: Start of the requested time range.
        end_date: End of the requested time range.
        variables: Optional list of variables; if omitted the defaults
            from ``conf/config.yaml`` are used.
        extra_overrides: Additional Hydra override strings appended last.
        box: Bounding box dict with keys lat_min/lat_max/lon_min/lon_max.
    """
    overrides: List[str] = [f"dataset={dataset}"]

    if box is not None:
        overrides.append(f"+box.lat_min={box['lat_min']}")
        overrides.append(f"+box.lat_max={box['lat_max']}")
        overrides.append(f"+box.lon_min={box['lon_min']}")
        overrides.append(f"+box.lon_max={box['lon_max']}")
    else:
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

    if data_dir:
        overrides.append(f"data_dir={data_dir}")

    if extra_overrides:
        overrides.extend(extra_overrides)

    return overrides
