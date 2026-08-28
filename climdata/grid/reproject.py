"""Reprojection and grid alignment, composed from rasterio/rioxarray.

The heavy lifting is all upstream: :func:`rasterio.warp.calculate_default_transform`
derives the destination transform and shape from a resolution,
:func:`rasterio.warp.aligned_target` snaps the origin so that independently-cropped
inputs land on identical cell centres, and ``DataArray.rio.reproject`` performs the
warp. What this module adds is validation (the compatibility gate), this package's
CRS/dimension inference, and metadata preservation across the warp.
"""

from __future__ import annotations

import warnings
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
import xarray as xr

from .crs import (
    crs_axis_unit,
    infer_spatial_dims,
    infer_src_crs,
    normalize_longitude,
    parse_crs,
)
from .units import resolution_in_crs_units

__all__ = ["reproject", "resampling_from_name"]


def _ensure_rioxarray():
    """Register the ``.rio`` accessor, with an actionable error if it is missing."""
    try:
        import rioxarray  # noqa: F401
    except ImportError as err:  # pragma: no cover - depends on the environment
        raise ImportError(
            "climdata.grid.reproject requires rioxarray. Install it with "
            "`pip install rioxarray` (or `conda install -c conda-forge rioxarray`)."
        ) from err

#: Coordinate names written on output, by target CRS kind.
_OUT_NAMES = {"angular": ("lon", "lat"), "linear": ("x", "y")}

#: Latitude span beyond which Web Mercator scale distortion is worth flagging.
_MERCATOR_SPAN_WARN = 20.0


def resampling_from_name(name: Any):
    """Look up a :class:`rasterio.enums.Resampling` member by name.

    Args:
        name: A method name such as ``"bilinear"`` or ``"average"``. A ``Resampling``
            member passes through unchanged.

    Returns:
        rasterio.enums.Resampling: The matching member.

    Raises:
        ValueError: If ``name`` is not a valid resampling method.
    """
    from rasterio.enums import Resampling

    if isinstance(name, Resampling):
        return name
    try:
        return Resampling[str(name).lower()]
    except KeyError as err:
        raise ValueError(
            f"Unknown resampling method {name!r}. Valid options: "
            f"{', '.join(Resampling.__members__)}."
        ) from err


def _prepare(obj, src_crs, x_dim, y_dim, dsinfo_dims, nodata):
    """Attach CRS/dims, wrap longitudes, and set nodata so gaps do not become zeros."""
    _ensure_rioxarray()
    x_dim, y_dim = infer_spatial_dims(obj, x_dim, y_dim, dsinfo_dims)
    src_crs = parse_crs(src_crs) if src_crs is not None else infer_src_crs(obj)

    if src_crs.is_geographic:
        obj = normalize_longitude(obj, x_dim)

    renamed = {}
    if x_dim != "x":
        renamed[x_dim] = "x"
    if y_dim != "y":
        renamed[y_dim] = "y"
    if renamed:
        obj = obj.rename(renamed)

    obj = obj.rio.set_spatial_dims(x_dim="x", y_dim="y", inplace=False)
    obj = obj.rio.write_crs(src_crs, inplace=False)

    # Without an explicit nodata, rasterio fills untouched cells with 0, which
    # silently corrupts masked domains such as HYRAS-over-Germany.
    if nodata is not None:
        obj = _write_nodata(obj, nodata)
    return obj, src_crs


def _write_nodata(obj, nodata):
    """Write ``nodata`` on every floating-point variable."""
    if isinstance(obj, xr.Dataset):
        for name, var in obj.data_vars.items():
            if np.issubdtype(var.dtype, np.floating):
                obj[name] = var.rio.write_nodata(nodata, inplace=False)
        return obj
    if np.issubdtype(obj.dtype, np.floating):
        return obj.rio.write_nodata(nodata, inplace=False)
    return obj


def _target_grid(obj, src_crs, dst_crs, resolution, bounds, align):
    """Return ``(transform, width, height)`` for the destination grid."""
    from affine import Affine
    from rasterio.warp import aligned_target, calculate_default_transform

    if bounds is not None:
        left, bottom, right, top = (float(v) for v in bounds)
        res_x, res_y = resolution
        width = max(1, int(round((right - left) / res_x)))
        height = max(1, int(round((top - bottom) / res_y)))
        transform = Affine.translation(left, top) * Affine.scale(res_x, -res_y)
    else:
        transform, width, height = calculate_default_transform(
            src_crs,
            dst_crs,
            obj.rio.width,
            obj.rio.height,
            *obj.rio.bounds(),
            resolution=resolution,
        )

    if align:
        transform, width, height = aligned_target(transform, width, height, resolution)
    return transform, width, height


def _warn_if_mercator(dst_crs, obj, src_crs):
    """Flag Web Mercator, whose metre axis is only true-scale at the equator."""
    if dst_crs.to_epsg() != 3857:
        return
    try:
        if src_crs.is_geographic:
            span = float(obj["y"].max()) - float(obj["y"].min())
        else:
            from pyproj import Transformer

            transformer = Transformer.from_crs(src_crs, "EPSG:4326", always_xy=True)
            left, bottom, right, top = obj.rio.bounds()
            _, lat_min = transformer.transform(left, bottom)
            _, lat_max = transformer.transform(right, top)
            span = abs(lat_max - lat_min)
    except Exception:
        return
    if span > _MERCATOR_SPAN_WARN:
        warnings.warn(
            "Reprojecting to EPSG:3857 (Web Mercator) across "
            f"{span:.0f} degrees of latitude. Its metre axis is only true-scale at "
            "the equator, so cell areas vary by more than a factor of two across "
            "this domain. Use an equal-area CRS (e.g. EPSG:3035) for area analysis.",
            UserWarning,
            stacklevel=3,
        )


def _restore_attrs(result, source):
    """Put back variable attributes that rioxarray drops during the warp."""
    result.attrs.update({k: v for k, v in source.attrs.items() if k not in result.attrs})
    if isinstance(result, xr.Dataset) and isinstance(source, xr.Dataset):
        for name, var in result.data_vars.items():
            if name in source.data_vars:
                keep = {
                    k: v
                    for k, v in source[name].attrs.items()
                    if k not in ("grid_mapping", "_FillValue")
                }
                var.attrs.update({k: v for k, v in keep.items() if k not in var.attrs})
    return result


def _rename_output(result, dst_crs):
    """Name output coordinates ``lon``/``lat`` for geographic CRSs, ``x``/``y`` otherwise."""
    kind, _ = crs_axis_unit(dst_crs)
    x_name, y_name = _OUT_NAMES[kind]
    renamed = {}
    if x_name != "x" and "x" in result.dims:
        renamed["x"] = x_name
    if y_name != "y" and "y" in result.dims:
        renamed["y"] = y_name
    if renamed:
        result = result.rename(renamed)
        result = result.rio.set_spatial_dims(x_dim=x_name, y_dim=y_name, inplace=False)
    return result


def reproject(
    obj: Union[xr.Dataset, xr.DataArray],
    target_projection: Any = None,
    target_resolution: Any = None,
    *,
    like: Optional[Union[xr.Dataset, xr.DataArray]] = None,
    method: Union[str, Mapping[str, str]] = "bilinear",
    bounds: Optional[Sequence[float]] = None,
    align: bool = True,
    src_crs: Any = None,
    x_dim: Optional[str] = None,
    y_dim: Optional[str] = None,
    nodata: Optional[float] = np.nan,
    engine: str = "rasterio",
    dsinfo_dims: Optional[Sequence[str]] = None,
) -> Union[xr.Dataset, xr.DataArray]:
    """Reproject and/or resample a dataset onto a new grid.

    Args:
        obj: The :class:`xarray.Dataset` or :class:`xarray.DataArray` to transform.
        target_projection: Destination CRS -- ``"EPSG:4326"``, an integer EPSG code,
            OGC WKT, a PROJ string, or a ``pyproj``/``rasterio`` CRS object. Defaults
            to the source CRS (resample without reprojecting).
        target_resolution: Destination cell size -- ``"10 km"``, ``"0.1 deg"``,
            ``"30 arcsec"``, a bare number (interpreted in the target CRS axis units),
            or a 2-tuple for an anisotropic grid. Its units must match the target CRS:
            a metric resolution with a geographic CRS raises
            :class:`~climdata.grid.units.ResolutionCRSMismatch`.
        like: A template whose grid the output should match exactly. Cannot be
            combined with ``target_projection``/``target_resolution``.
        method: Resampling method name, or a ``{variable: method}`` mapping. Use
            ``"bilinear"`` for state variables, ``"average"`` for fluxes,
            ``"nearest"`` for categorical fields.
        bounds: Optional ``(left, bottom, right, top)`` target extent, expressed in
            the units of the target CRS.
        align: Snap the grid origin to whole multiples of the resolution so that
            separately-processed domains stay co-registered. Default ``True``.
        src_crs: Source CRS, when it cannot be inferred from the data.
        x_dim: Source x/longitude dimension name.
        y_dim: Source y/latitude dimension name.
        nodata: Fill value written before warping. Keep the ``np.nan`` default for
            float data; passing ``None`` lets rasterio fill gaps with zeros.
        engine: ``"rasterio"`` (default) or ``"xesmf"`` for true area-weighted
            conservative remapping.
        dsinfo_dims: Optional ``(x_dim, y_dim)`` fallback, typically from
            ``cfg.dsinfo[DATASET]``.

    Returns:
        The reprojected object, of the same type as ``obj``, with coordinates named
        ``lon``/``lat`` for a geographic target and ``x``/``y`` for a projected one.

    Raises:
        ResolutionCRSMismatch: If the resolution units contradict the target CRS.
        ValueError: If ``like`` is combined with an explicit target, or if inputs
            cannot be interpreted.

    Example:
        >>> out = reproject(ds, "EPSG:3035", "10 km")            # doctest: +SKIP
        >>> out = reproject(ds, "EPSG:4326", "0.1 deg",          # doctest: +SKIP
        ...                 method={"pr": "average", "tas": "bilinear"})
    """
    if like is not None and (target_projection is not None or target_resolution is not None):
        raise ValueError(
            "`like` already fixes both the projection and the resolution; passing "
            "target_projection or target_resolution alongside it is contradictory. "
            "Drop one."
        )

    prepared, source_crs = _prepare(obj, src_crs, x_dim, y_dim, dsinfo_dims, nodata)

    if like is not None:
        template, _ = _prepare(like, None, None, None, None, nodata)
        result = prepared.rio.reproject_match(
            template, resampling=resampling_from_name(method if isinstance(method, str) else "bilinear")
        )
        return _rename_output(_restore_attrs(result, obj), infer_src_crs(template))

    dst_crs = parse_crs(target_projection) if target_projection is not None else source_crs

    # The gate runs on the *effective* target CRS, after defaulting: passing no
    # projection at all must fail exactly as naming the source CRS would.
    if target_resolution is not None:
        try:
            mid_lat = _mid_latitude(prepared, source_crs)
        except Exception:
            mid_lat = None
        resolution = resolution_in_crs_units(target_resolution, dst_crs, latitude=mid_lat)
    else:
        resolution = None

    _warn_if_mercator(dst_crs, prepared, source_crs)

    if resolution is None:
        result = _reproject_dispatch(prepared, dst_crs, None, None, None, method, nodata, engine)
    else:
        transform, width, height = _target_grid(
            prepared, source_crs, dst_crs, resolution, bounds, align
        )
        result = _reproject_dispatch(
            prepared, dst_crs, transform, width, height, method, nodata, engine
        )

    return _rename_output(_restore_attrs(result, obj), dst_crs)


def _mid_latitude(obj, src_crs) -> Optional[float]:
    """Mid-latitude of the domain in degrees, for concrete error messages."""
    left, bottom, right, top = obj.rio.bounds()
    if src_crs.is_geographic:
        return (bottom + top) / 2.0
    from pyproj import Transformer

    transformer = Transformer.from_crs(src_crs, "EPSG:4326", always_xy=True)
    _, lat_min = transformer.transform(left, bottom)
    _, lat_max = transformer.transform(right, top)
    return (lat_min + lat_max) / 2.0


def _reproject_dispatch(obj, dst_crs, transform, width, height, method, nodata, engine):
    """Warp with rasterio, or hand off to xESMF for conservative remapping."""
    if engine == "xesmf":
        return _reproject_xesmf(obj, dst_crs, transform, width, height, method)
    if engine != "rasterio":
        raise ValueError(f"Unknown engine {engine!r}. Use 'rasterio' or 'xesmf'.")

    kwargs: Dict[str, Any] = {}
    if transform is not None:
        kwargs["transform"] = transform
        kwargs["shape"] = (height, width)

    if isinstance(method, Mapping):
        if not isinstance(obj, xr.Dataset):
            raise ValueError("A per-variable `method` mapping requires a Dataset.")
        warped = {}
        for name, var in obj.data_vars.items():
            resampling = resampling_from_name(method.get(name, "bilinear"))
            warped[name] = var.rio.reproject(dst_crs, resampling=resampling, **kwargs)
        return xr.Dataset(warped, attrs=obj.attrs)

    return obj.rio.reproject(
        dst_crs, resampling=resampling_from_name(method), **kwargs
    )


def _reproject_xesmf(obj, dst_crs, transform, width, height, method):
    """Delegate to the existing xESMF regridder for area-weighted remapping."""
    dst_crs = parse_crs(dst_crs)
    if not dst_crs.is_geographic:
        raise ValueError(
            "engine='xesmf' regrids on a latitude/longitude mesh, so it needs a "
            "geographic target_projection (e.g. 'EPSG:4326'). Use engine='rasterio' "
            "for projected targets."
        )
    if transform is None:
        raise ValueError("engine='xesmf' requires an explicit target_resolution.")

    from climdata.sdba.bcsd import regrid_to_coarse

    lon = transform.c + transform.a * (np.arange(width) + 0.5)
    lat = transform.f + transform.e * (np.arange(height) + 0.5)
    template = xr.Dataset({"lat": ("lat", lat), "lon": ("lon", lon)})

    source = obj.rename({"x": "lon", "y": "lat"})
    if isinstance(source, xr.DataArray):
        source = source.to_dataset(name=source.name or "data")
    xesmf_method = method if isinstance(method, str) else "conservative"
    if xesmf_method in ("average", "sum"):
        xesmf_method = "conservative"
    return regrid_to_coarse(source, template, method=xesmf_method)
