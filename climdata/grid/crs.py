"""CRS coercion and inference for :mod:`climdata.grid`.

Projection parsing delegates entirely to :func:`pyproj.CRS.from_user_input`, which
already accepts ``"EPSG:4326"``, integer EPSG codes, WKT1/WKT2, PROJ strings,
CF ``grid_mapping`` dicts, :class:`rasterio.crs.CRS` and :class:`pyproj.CRS`.

What is genuinely specific to ``climdata`` lives here too: inferring the source CRS
from the conventions this package already writes to disk (notably the ``crs_grid``
attribute written by the HYRAS reader), inferring spatial dimension names, and
normalising 0-360 longitudes that arrive from CMIP / NEX-GDDP.
"""

from __future__ import annotations

import math
import warnings
from typing import Any, Optional, Sequence, Tuple

__all__ = [
    "parse_crs",
    "crs_axis_unit",
    "infer_spatial_dims",
    "infer_src_crs",
    "normalize_longitude",
]

#: Attribute keys that may carry a CRS, most authoritative first. ``esri_pe_string``
#: matters in practice: real HYRAS files carry the full ESRI WKT there, on the data
#: variable, while the CF ``crs`` variable it points at is dropped during subsetting.
_CRS_ATTR_KEYS = (
    "epsg_code",
    "crs_wkt",
    "spatial_ref",
    "esri_pe_string",
    "proj4",
    "proj4_params",
    "crs_grid",
    "crs",
)

#: Candidate dimension names, most specific first.
_X_NAMES = ("x", "lon", "longitude", "rlon", "nav_lon", "X")
_Y_NAMES = ("y", "lat", "latitude", "rlat", "nav_lat", "Y")


def parse_crs(value: Any):
    """Coerce ``value`` to a :class:`pyproj.CRS`.

    Args:
        value: Anything :func:`pyproj.CRS.from_user_input` understands --
            ``"EPSG:4326"``, ``4326``, an OGC WKT string, a PROJ string,
            a CF ``grid_mapping`` dict, a ``rasterio.crs.CRS`` or a ``pyproj.CRS``.
            ``None`` passes through as ``None``.

    Returns:
        pyproj.CRS | None: The parsed CRS.

    Raises:
        ValueError: If ``value`` cannot be interpreted as a CRS.
    """
    if value is None:
        return None

    import pyproj

    if isinstance(value, pyproj.CRS):
        return value
    try:
        return pyproj.CRS.from_user_input(value)
    except Exception as err:  # pyproj raises CRSError, a subclass of RuntimeError
        raise ValueError(
            f"Could not interpret {value!r} as a coordinate reference system. "
            "Expected an EPSG code (e.g. 'EPSG:4326' or 4326), an OGC WKT string, "
            "a PROJ string (e.g. '+proj=laea +lat_0=52 ...'), or a pyproj/rasterio "
            "CRS object. Look codes up at https://epsg.io."
        ) from err


def _crs_from_attrs(attrs) -> Optional[Any]:
    """Return the first attribute of ``attrs`` that parses as a CRS, else ``None``."""
    for key in _CRS_ATTR_KEYS:
        value = attrs.get(key)
        if not value or not isinstance(value, str):
            continue
        try:
            return parse_crs(value)
        except ValueError:
            continue
    return None


def _horizontal(crs):
    """Reduce a compound/3D CRS to its horizontal part."""
    if getattr(crs, "is_compound", False) and crs.sub_crs_list:
        return crs.sub_crs_list[0]
    return crs


def crs_axis_unit(crs) -> Tuple[str, float]:
    """Return the axis unit kind and conversion factor of ``crs``.

    The kind is derived from :attr:`pyproj.CRS.is_geographic` -- never from the
    position of an axis. EPSG:4326 lists *latitude* first, so positional indexing
    silently mis-assigns anisotropic resolutions.

    Args:
        crs: Anything :func:`parse_crs` accepts.

    Returns:
        tuple[str, float]: ``("angular", radians_per_unit)`` for a geographic CRS
        or ``("linear", metres_per_unit)`` for a projected one. For EPSG:4326 the
        factor is ``pi/180``; for EPSG:3035 it is ``1.0``; for EPSG:2263
        (US survey foot) it is ``0.3048006096``.

    Raises:
        ValueError: If the two horizontal axes carry different units, or if the CRS
            is neither geographic nor projected.
    """
    crs = _horizontal(parse_crs(crs))
    axes = crs.axis_info[:2]
    if len(axes) < 2:
        raise ValueError(f"CRS {crs.name!r} does not expose two horizontal axes.")
    if axes[0].unit_name != axes[1].unit_name:
        raise ValueError(
            f"CRS {crs.name!r} has mixed horizontal axis units "
            f"({axes[0].unit_name!r} and {axes[1].unit_name!r}); a single "
            "resolution cannot describe such a grid."
        )

    factor = float(axes[0].unit_conversion_factor)
    if crs.is_geographic:
        return "angular", factor
    if crs.is_projected:
        return "linear", factor
    raise ValueError(
        f"CRS {crs.name!r} is neither geographic nor projected, so it has no "
        "well-defined grid resolution."
    )


def infer_spatial_dims(
    obj,
    x_dim: Optional[str] = None,
    y_dim: Optional[str] = None,
    dsinfo_dims: Optional[Sequence[str]] = None,
) -> Tuple[str, str]:
    """Infer the horizontal dimension names of ``obj``.

    Resolution order: explicit arguments, then whatever ``rioxarray`` already knows,
    then a table of conventional names, then ``dsinfo_dims``.

    Args:
        obj: An :class:`xarray.Dataset` or :class:`xarray.DataArray`.
        x_dim: Explicit x/longitude dimension name.
        y_dim: Explicit y/latitude dimension name.
        dsinfo_dims: Optional ``(x_dim, y_dim)`` fallback, typically read from
            ``cfg.dsinfo[DATASET].lon_dim`` / ``.lat_dim``.

    Returns:
        tuple[str, str]: ``(x_dim, y_dim)``.

    Raises:
        ValueError: If the dimensions cannot be determined.
    """
    if x_dim and y_dim:
        return x_dim, y_dim

    try:  # rioxarray may already have them set
        import rioxarray  # noqa: F401

        return obj.rio.x_dim, obj.rio.y_dim
    except Exception:
        pass

    dims = set(obj.dims) | set(obj.coords)
    found_x = x_dim or next((n for n in _X_NAMES if n in dims), None)
    found_y = y_dim or next((n for n in _Y_NAMES if n in dims), None)

    if (not found_x or not found_y) and dsinfo_dims:
        cand_x, cand_y = dsinfo_dims
        found_x = found_x or (cand_x if cand_x in dims else None)
        found_y = found_y or (cand_y if cand_y in dims else None)

    if not found_x or not found_y:
        raise ValueError(
            f"Could not infer spatial dimensions from {sorted(map(str, obj.dims))}. "
            "Pass x_dim= and y_dim= explicitly."
        )
    return found_x, found_y


def infer_src_crs(obj, default: Any = "EPSG:4326"):
    """Infer the CRS of ``obj``, following the conventions found in real files.

    Order: an existing ``rio.crs``; dataset-level attributes; a CF ``grid_mapping``
    variable, when it survived subsetting; CRS attributes carried on the data variables
    themselves (real HYRAS keeps its WKT in ``esri_pe_string`` there, and the ``crs``
    variable its ``grid_mapping`` points at is dropped); a ``spatial_ref`` coordinate;
    then a geographic heuristic on the coordinate ranges.

    Args:
        obj: An :class:`xarray.Dataset` or :class:`xarray.DataArray`.
        default: CRS assumed when the coordinates look geographic. Pass ``None``
            to make that case raise instead.

    Returns:
        pyproj.CRS: The inferred CRS.

    Raises:
        ValueError: If no CRS can be determined.
    """
    try:
        import rioxarray  # noqa: F401

        if obj.rio.crs is not None:
            return parse_crs(obj.rio.crs)
    except Exception:
        pass

    found = _crs_from_attrs(obj.attrs)
    if found is not None:
        return found

    # A CF grid_mapping pointer, when the variable it names still exists.
    variables = getattr(obj, "variables", {})
    candidates = list(getattr(obj, "data_vars", {}).values()) or [obj]
    for var in candidates:
        target = var.attrs.get("grid_mapping")
        if target and target in variables:
            found = _crs_from_attrs(variables[target].attrs)
            if found is not None:
                return found

    # CRS attributes carried directly on the data variables.
    for var in candidates:
        found = _crs_from_attrs(var.attrs)
        if found is not None:
            return found

    if "spatial_ref" in getattr(obj, "coords", {}):
        found = _crs_from_attrs(obj.coords["spatial_ref"].attrs)
        if found is not None:
            return found

    if default is not None:
        try:
            x_dim, y_dim = infer_spatial_dims(obj)
            x, y = obj[x_dim], obj[y_dim]
            if float(abs(y).max()) <= 90.0 and float(abs(x).max()) <= 360.0:
                warnings.warn(
                    f"No CRS found on the dataset; assuming {default} because the "
                    f"{y_dim!r}/{x_dim!r} coordinates fall in geographic ranges. "
                    "Pass src_crs= to be explicit.",
                    UserWarning,
                    stacklevel=2,
                )
                return parse_crs(default)
        except (ValueError, KeyError):
            pass

    raise ValueError(
        "Could not determine the source CRS of the dataset. Set it with "
        "ds.rio.write_crs(...) or pass src_crs= explicitly."
    )


def normalize_longitude(obj, x_dim: str):
    """Roll a 0-360 longitude axis onto -180-180 and re-sort.

    CMIP and NEX-GDDP data arrive on 0-360. ``rasterio`` will happily warp such an
    axis without raising and produce nonsense, so this runs unconditionally for
    geographic sources whose longitudes exceed 180.

    Args:
        obj: An :class:`xarray.Dataset` or :class:`xarray.DataArray`.
        x_dim: Name of the longitude dimension.

    Returns:
        The object with longitudes in -180..180, sorted ascending. Returned
        unchanged when no rolling is needed.
    """
    if x_dim not in obj.coords:
        return obj
    x = obj[x_dim]
    if float(x.max()) <= 180.0:
        return obj
    rolled = ((x + 180.0) % 360.0) - 180.0
    return obj.assign_coords({x_dim: rolled}).sortby(x_dim)
