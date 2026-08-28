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
    experiment_id: Optional[str] = None,
    source_id: Optional[str] = None,
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
        experiment_id: CMIP6 experiment/scenario (e.g. ``"ssp585"``); only
            meaningful for model-driven datasets such as ``"cmip"``.
        source_id: CMIP6 model / GCM name (e.g. ``"MIROC6"``).
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

    if experiment_id:
        overrides.append(f"experiment_id={experiment_id}")
    if source_id:
        overrides.append(f"source_id={source_id}")

    if variables:
        vars_str = "[" + ",".join(variables) + "]"
        overrides.append(f"variables={vars_str}")

    if data_dir:
        overrides.append(f"data_dir={data_dir}")
        overrides.append(f"output.out_dir={data_dir}")

    if extra_overrides:
        overrides.extend(extra_overrides)

    return overrides
