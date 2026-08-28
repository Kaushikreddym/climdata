"""Resolution parsing and the projection/resolution compatibility gate.

Unit parsing is delegated to :func:`xclim.core.units.str2pint` (pint), which already
understands ``"10 km"``, ``"0.1 deg"``, ``"30 arcsec"``, ``"0.5 mi"`` and ``"30 ft"``.
Two behaviours of that stack need handling here:

1. ``str2pint`` rejects unspaced input -- ``"10km"`` raises
   ``ValueError: Unit expression cannot have a scaling factor.`` Unspaced is exactly
   what gets typed on a Hydra command line, so :func:`_normalize` inserts the space.
2. In pint, degrees are *dimensionless*, so a bare number and an angle share a
   dimensionality. ``str2pint("0.1").to("degree")`` returns **5.7296** -- pint reads
   the bare value as radians. Bare numbers are therefore classified by *unit*, not by
   dimensionality, and are never passed through a conversion.

The compatibility gate itself has no upstream equivalent: a linear resolution has no
fixed size in a geographic CRS (10 km spans 0.1395 deg of longitude but 0.0899 deg of
latitude at 50N, and 0.2619 deg at 70N), so the combination is rejected rather than
silently approximated.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from typing import Any, Optional, Tuple, Union

from .crs import crs_axis_unit, parse_crs

__all__ = [
    "Resolution",
    "ResolutionCRSMismatch",
    "parse_resolution",
    "resolution_in_crs_units",
    "to_angular",
]

_NUM = r"^([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)\s*"
_DIMENSIONLESS = ("dimensionless", "1", "")
_DEG_PER_RAD = 180.0 / math.pi


class ResolutionCRSMismatch(ValueError):
    """``target_resolution`` units are incompatible with ``target_projection``.

    Subclasses :class:`ValueError` so existing ``except ValueError`` handlers --
    including the LLM execution layer -- keep working.
    """


@dataclass(frozen=True)
class Resolution:
    """A parsed grid resolution.

    Attributes:
        x: Magnitude along the x/longitude axis, in :attr:`unit`.
        y: Magnitude along the y/latitude axis, in :attr:`unit`.
        kind: ``"linear"``, ``"angular"`` or ``"unitless"``.
        unit: Canonical unit -- ``"m"`` for linear, ``"deg"`` for angular,
            ``""`` for unitless (adopts the target CRS axis unit).
    """

    x: float
    y: float
    kind: str
    unit: str

    def as_tuple(self) -> Tuple[float, float]:
        """Return ``(x, y)`` magnitudes."""
        return (self.x, self.y)

    def __str__(self) -> str:
        if self.x == self.y:
            return f"{self.x:g} {self.unit}".strip()
        return f"{self.x:g} x {self.y:g} {self.unit}".strip()


def _normalize(text: str) -> str:
    """Make a resolution string palatable to ``str2pint``.

    Inserts the space ``str2pint`` requires between magnitude and unit, and expands
    the degree/minute/second symbols. The ``'`` and ``"`` substitutions are anchored
    to end-of-string: unanchored, ``"'x'"`` corrupts to ``arcminx arcmin``.
    """
    text = text.strip().replace("°", " deg")
    text = re.sub(r"'$", " arcmin", text)
    text = re.sub(r'"$', " arcsec", text)
    return re.sub(_NUM, r"\1 ", text).strip()


def _classify(quantity) -> str:
    """Return ``"linear"``, ``"angular"`` or ``"unitless"`` for a pint quantity."""
    if str(quantity.dimensionality) == "[length]":
        return "linear"
    if str(quantity.units) in _DIMENSIONLESS:
        return "unitless"  # never convert: pint would read it as radians
    if str(quantity.dimensionality) in _DIMENSIONLESS:
        return "angular"
    raise ValueError(
        f"Resolution {quantity!r} has dimensionality "
        f"{quantity.dimensionality!s}, which is neither a length nor an angle."
    )


def _parse_scalar(value: Any) -> Tuple[float, str]:
    """Parse one component into ``(magnitude_in_canonical_unit, kind)``."""
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return float(value), "unitless"

    if hasattr(value, "magnitude") and hasattr(value, "units"):
        quantity = value
    elif isinstance(value, str):
        from xclim.core.units import str2pint

        text = _normalize(value)
        if not re.match(_NUM, text):
            # "km" parses as 1 kilometre and "" as dimensionless 1; neither is a
            # resolution, so require an explicit magnitude.
            raise ValueError(
                f"Resolution {value!r} has no magnitude. Write the number and the "
                "unit together, e.g. '10 km' or '0.1 deg'."
            )
        try:
            quantity = str2pint(text)
        except Exception as err:
            raise ValueError(
                f"Could not parse resolution {value!r}. Expected a magnitude with a "
                "length unit (e.g. '10 km', '1000 m', '30 ft') or an angular unit "
                "(e.g. '0.1 deg', '30 arcsec', '5 arcmin')."
            ) from err
    else:
        raise TypeError(
            f"Resolution must be a str, number, pint Quantity or 2-tuple; "
            f"got {type(value).__name__}."
        )

    kind = _classify(quantity)
    if kind == "linear":
        return float(quantity.to("meter").magnitude), kind
    if kind == "angular":
        return float(quantity.to("degree").magnitude), kind
    return float(quantity.magnitude), kind


def parse_resolution(value: Any) -> Resolution:
    """Parse a resolution into canonical units.

    Args:
        value: A string (``"10 km"``, ``"0.1 deg"``, ``"30 arcsec"``, ``"10km"``,
            ``"0.1°"``), a bare number, a pint ``Quantity``, a
            :class:`Resolution`, or a 2-tuple ``(x, y)`` of any of those for an
            anisotropic grid.

    Returns:
        Resolution: Linear values in metres, angular values in degrees, bare numbers
        left alone.

    Raises:
        ValueError: If the string cannot be parsed, or if the two components of a
            tuple are of different kinds.
        TypeError: If ``value`` is of an unsupported type.

    Example:
        >>> parse_resolution("10 km")
        Resolution(x=10000.0, y=10000.0, kind='linear', unit='m')
        >>> parse_resolution("0.1 deg").as_tuple()
        (0.1, 0.1)
    """
    if isinstance(value, Resolution):
        return value

    if isinstance(value, (tuple, list)):
        if len(value) != 2:
            raise ValueError(
                f"An anisotropic resolution needs exactly 2 components, got {len(value)}."
            )
        x_mag, x_kind = _parse_scalar(value[0])
        y_mag, y_kind = _parse_scalar(value[1])
        if x_kind != y_kind:
            raise ValueError(
                f"Mixed resolution kinds: {value[0]!r} is {x_kind} but {value[1]!r} "
                f"is {y_kind}. Both axes must use the same kind of unit."
            )
        kind = x_kind
    else:
        x_mag, kind = _parse_scalar(value)
        y_mag = x_mag

    unit = {"linear": "m", "angular": "deg", "unitless": ""}[kind]
    return Resolution(x=x_mag, y=y_mag, kind=kind, unit=unit)


def _mismatch_message(
    resolution: Resolution,
    crs,
    crs_kind: str,
    latitude: Optional[float],
    original: Any = None,
) -> str:
    """Build the ResolutionCRSMismatch text, with numbers for the actual domain.

    ``original`` is echoed verbatim so the message quotes what the caller typed
    ("10 km") rather than the canonicalised form ("10000 m").
    """
    crs_name = crs.to_string() if hasattr(crs, "to_string") else str(crs)
    shown = original if isinstance(original, str) else str(resolution)

    if resolution.kind == "linear" and crs_kind == "angular":
        lat = 50.0 if latitude is None else float(latitude)
        d_lon, d_lat = _geod_deltas(resolution.x, resolution.y, lat)
        far_lon, far_lat = _geod_deltas(resolution.x, resolution.y, 70.0)
        return (
            f'target_resolution="{shown}" is a linear resolution, but '
            f'target_projection="{crs_name}" is a geographic CRS with axis units '
            f"of degree.\n\n"
            f"A metric resolution has no fixed angular size: {shown} spans "
            f"{d_lon:.4f} deg of longitude but {d_lat:.4f} deg of latitude at "
            f"{lat:.1f}N, and {far_lon:.4f} deg x {far_lat:.4f} deg at 70N. There is "
            f"no single correct value for this grid.\n\n"
            f"Choose one:\n"
            f'  - angular resolution   target_resolution="0.1 deg"\n'
            f'  - projected CRS        target_projection="EPSG:3035"   '
            f"(ETRS89-LAEA, metres - Europe)\n"
            f'                         target_projection="EPSG:32632"  (UTM 32N, metres)'
        )

    return (
        f'target_resolution="{shown}" is an angular resolution, but '
        f'target_projection="{crs_name}" is a projected CRS with linear axis units.\n\n'
        f"A degree is not a fixed distance, so it cannot size a projected grid.\n\n"
        f"Choose one:\n"
        f'  - linear resolution     target_resolution="10 km"\n'
        f'  - geographic CRS        target_projection="EPSG:4326"'
    )


def _geod_deltas(x_metres: float, y_metres: float, latitude: float) -> Tuple[float, float]:
    """Return ``(delta_lon, delta_lat)`` in degrees for a metric step at ``latitude``."""
    from pyproj import Geod

    geod = Geod(ellps="WGS84")
    lon_east, _, _ = geod.fwd(0.0, latitude, 90.0, x_metres)
    _, lat_north, _ = geod.fwd(0.0, latitude, 0.0, y_metres)
    return abs(lon_east), abs(lat_north - latitude)


def resolution_in_crs_units(
    value: Any, crs: Any, *, latitude: Optional[float] = None
) -> Tuple[float, float]:
    """Express a resolution in the axis units of ``crs``, or refuse to.

    This is the single enforcement point for the projection/resolution compatibility
    rules, so every entry point -- direct call, ``ClimateExtractor.reproject``, a
    Hydra override, an LLM tool call -- inherits the check. It runs before any data is
    read, so a bad configuration fails immediately rather than after a download.

    Args:
        value: Anything :func:`parse_resolution` accepts.
        crs: Anything :func:`~climdata.grid.crs.parse_crs` accepts.
        latitude: Reference latitude used only to make the error message concrete.

    Returns:
        tuple[float, float]: ``(x, y)`` in the axis units of ``crs``.

    Raises:
        ResolutionCRSMismatch: For a linear resolution with a geographic CRS, or an
            angular resolution with a projected CRS.

    Example:
        >>> resolution_in_crs_units("10 km", "EPSG:3035")
        (10000.0, 10000.0)
        >>> resolution_in_crs_units("10 km", "EPSG:4326")
        Traceback (most recent call last):
            ...
        climdata.grid.units.ResolutionCRSMismatch: ...
    """
    resolution = parse_resolution(value)
    crs = parse_crs(crs)
    crs_kind, factor = crs_axis_unit(crs)

    if resolution.kind == "unitless":
        return resolution.as_tuple()  # adopt the CRS axis units verbatim

    if resolution.kind != crs_kind:
        raise ResolutionCRSMismatch(
            _mismatch_message(resolution, crs, crs_kind, latitude, original=value)
        )

    if crs_kind == "linear":
        # `factor` is metres per CRS unit (1.0 for metre, 0.3048006096 for US survey foot).
        return (resolution.x / factor, resolution.y / factor)

    # `factor` is radians per CRS unit; for a degree axis that is pi/180.
    deg_per_unit = factor * _DEG_PER_RAD
    return (resolution.x / deg_per_unit, resolution.y / deg_per_unit)


def to_angular(value: Any, latitude: float) -> Resolution:
    """Convert a linear resolution to degrees at an explicit latitude.

    The deliberate escape hatch from :class:`ResolutionCRSMismatch`, for a caller who
    wants the approximation and is prepared to own it. ``latitude`` has no default, so
    the approximation cannot be made by accident, and there is intentionally no config
    flag that reaches this -- a YAML toggle is how such an approximation creeps back
    into production runs unnoticed.

    The result is exact at ``latitude`` and only there.

    Args:
        value: A linear resolution, e.g. ``"10 km"``.
        latitude: Reference latitude in degrees. Required.

    Returns:
        Resolution: An angular resolution in degrees.

    Raises:
        ValueError: If ``value`` is not a linear resolution.

    Example:
        >>> res = to_angular("10 km", latitude=50.0)
        >>> round(res.x, 4), round(res.y, 4)
        (0.1395, 0.0899)
    """
    resolution = parse_resolution(value)
    if resolution.kind != "linear":
        raise ValueError(
            f"to_angular() converts a linear resolution to degrees, but {value!r} is "
            f"{resolution.kind}."
        )
    d_lon, d_lat = _geod_deltas(resolution.x, resolution.y, float(latitude))
    return Resolution(x=d_lon, y=d_lat, kind="angular", unit="deg")
