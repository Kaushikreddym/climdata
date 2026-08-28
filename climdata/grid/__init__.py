"""Resolution and projection transformation for climate datasets.

Composes existing libraries rather than reimplementing them: unit parsing comes from
:mod:`xclim.core.units` (pint), CRS handling from :mod:`pyproj`, and warping plus grid
alignment from :mod:`rasterio` / :mod:`rioxarray`.

Example:
    >>> import climdata                                          # doctest: +SKIP
    >>> out = climdata.reproject(ds, "EPSG:3035", "10 km")       # doctest: +SKIP
"""

from .crs import (
    crs_axis_unit,
    infer_spatial_dims,
    infer_src_crs,
    normalize_longitude,
    parse_crs,
)
from .reproject import reproject, resampling_from_name
from .units import (
    Resolution,
    ResolutionCRSMismatch,
    parse_resolution,
    resolution_in_crs_units,
    to_angular,
)

__all__ = [
    "reproject",
    "resampling_from_name",
    "parse_crs",
    "crs_axis_unit",
    "infer_spatial_dims",
    "infer_src_crs",
    "normalize_longitude",
    "Resolution",
    "ResolutionCRSMismatch",
    "parse_resolution",
    "resolution_in_crs_units",
    "to_angular",
]
