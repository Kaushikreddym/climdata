"""Tests for :mod:`climdata.grid`.

The philosophy is to test *this package's* code, not its dependencies. pint's unit
table, pyproj's CRS database and rasterio's warp arithmetic are tested upstream;
duplicating that here would buy nothing. What is exercised instead is the
normalisation shim, the kind classification, the compatibility gate, the CRS/dimension
inference chain, and the end-to-end composition.
"""

import math

import numpy as np
import pytest
import xarray as xr

from climdata.grid import (
    ResolutionCRSMismatch,
    crs_axis_unit,
    parse_crs,
    parse_resolution,
    resolution_in_crs_units,
    to_angular,
)

rioxarray = pytest.importorskip("rioxarray")


# --------------------------------------------------------------------------
# Fixtures
# --------------------------------------------------------------------------
@pytest.fixture
def geo_ds():
    """A small 0.25 deg EPSG:4326 field over Germany."""
    lat = np.arange(47.0, 55.01, 0.25)
    lon = np.arange(5.0, 15.01, 0.25)
    rng = np.random.default_rng(42)
    ds = xr.Dataset(
        {
            "tas": (("lat", "lon"), rng.random((len(lat), len(lon))) + 280.0),
            "pr": (("lat", "lon"), rng.random((len(lat), len(lon)))),
        },
        coords={"lat": lat, "lon": lon},
    )
    ds.tas.attrs = {"units": "K", "long_name": "air temperature"}
    ds.pr.attrs = {"units": "mm d-1"}
    return ds.rio.write_crs("EPSG:4326")


@pytest.fixture
def hyras_v6_like():
    """The metadata shape of a real HYRAS v6-1 extraction.

    HYRAS is ETRS89-LAEA (EPSG:3035), not Gauss-Krueger. Subsetting drops the CF ``crs``
    variable but leaves ``grid_mapping`` pointing at it, so the only surviving CRS record
    is the ESRI WKT on the data variable.
    """
    esri = (
        'PROJCS["ETRS_1989_LAEA",GEOGCS["GCS_ETRS_1989",DATUM["D_ETRS_1989",'
        'SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],'
        'UNIT["Degree",0.0174532925199433]],PROJECTION["Lambert_Azimuthal_Equal_Area"],'
        'PARAMETER["False_Easting",4321000.0],PARAMETER["False_Northing",3210000.0],'
        'PARAMETER["Central_Meridian",10.0],PARAMETER["Latitude_Of_Origin",52.0],'
        'UNIT["Meter",1.0]]'
    )
    x = np.arange(4132500.0, 4156501.0, 1000.0)
    y = np.arange(2712500.0, 2736501.0, 1000.0)
    rng = np.random.default_rng(11)
    ds = xr.Dataset(
        {"tas": (("y", "x"), rng.random((len(y), len(x))) + 2.0)},
        coords={"x": x, "y": y},
    )
    ds.tas.attrs = {
        "units": "\u00b0C",
        "grid_mapping": "crs",          # dangling: the variable did not survive
        "esri_pe_string": esri,
    }
    return ds


@pytest.fixture
def hyras_like():
    """A dataset carrying its CRS in the ``crs_grid`` attribute (HYRAS reader path)."""
    x = np.arange(3280000.0, 3480001.0, 5000.0)
    y = np.arange(5240000.0, 5440001.0, 5000.0)
    rng = np.random.default_rng(7)
    ds = xr.Dataset(
        {"tasmin": (("y", "x"), rng.random((len(y), len(x))) + 270.0)},
        coords={"x": x, "y": y},
    )
    ds.attrs["crs_grid"] = "EPSG:31467"
    return ds


# --------------------------------------------------------------------------
# 1. Normalisation shim
# --------------------------------------------------------------------------
@pytest.mark.parametrize(
    "text, expected_x, expected_unit",
    [
        ("10 km", 10000.0, "m"),
        ("10km", 10000.0, "m"),          # str2pint alone rejects this
        ("1000 m", 1000.0, "m"),
        ("1e3 m", 1000.0, "m"),
        ("0.5 mi", 804.672, "m"),
        ("30 ft", 9.144, "m"),
        ("0.1 deg", 0.1, "deg"),
        ("0.1°", 0.1, "deg"),       # degree symbol
        ("0.25degree", 0.25, "deg"),
        ("30 arcsec", 1.0 / 120.0, "deg"),
        ("30arcsec", 1.0 / 120.0, "deg"),
        ("5 arcmin", 5.0 / 60.0, "deg"),
    ],
)
def test_parse_resolution_strings(text, expected_x, expected_unit):
    res = parse_resolution(text)
    assert res.x == pytest.approx(expected_x, rel=1e-9)
    assert res.y == pytest.approx(expected_x, rel=1e-9)
    assert res.unit == expected_unit


def test_symbol_substitution_is_anchored():
    """Unanchored ' / " replacement corrupts input: "'x'" -> "arcminx arcmin"."""
    with pytest.raises(ValueError):
        parse_resolution("'x'")


@pytest.mark.parametrize("bad", ["", "km", "10 bananas", "ten km"])
def test_parse_resolution_rejects_garbage(bad):
    with pytest.raises((ValueError, TypeError)):
        parse_resolution(bad)


def test_anisotropic_and_mixed_kinds():
    res = parse_resolution(("10 km", "5 km"))
    assert (res.x, res.y) == (10000.0, 5000.0)
    with pytest.raises(ValueError, match="Mixed resolution kinds"):
        parse_resolution(("10 km", "0.1 deg"))


# --------------------------------------------------------------------------
# 2. Kind classification -- the pint radians trap
# --------------------------------------------------------------------------
@pytest.mark.parametrize(
    "value, kind",
    [
        ("10 km", "linear"),
        ("30 ft", "linear"),
        ("0.1 deg", "angular"),
        ("30 arcsec", "angular"),
        ("0.1", "unitless"),
        (0.1, "unitless"),
    ],
)
def test_kind_classification(value, kind):
    assert parse_resolution(value).kind == kind


def test_bare_number_is_never_read_as_radians():
    """pint converts a bare 0.1 to 5.7296 degrees. It must stay 0.1."""
    assert parse_resolution("0.1").x == pytest.approx(0.1)
    assert resolution_in_crs_units("0.1", "EPSG:4326") == pytest.approx((0.1, 0.1))
    assert resolution_in_crs_units(0.1, "EPSG:4326")[0] != pytest.approx(5.7295779, abs=1e-3)


# --------------------------------------------------------------------------
# 3. The compatibility gate
# --------------------------------------------------------------------------
@pytest.mark.parametrize(
    "resolution, crs, expected",
    [
        ("0.1 deg", "EPSG:4326", (0.1, 0.1)),
        ("30 arcsec", "EPSG:4326", (1 / 120.0, 1 / 120.0)),
        ("10 km", "EPSG:3035", (10000.0, 10000.0)),
        ("10 km", "EPSG:32632", (10000.0, 10000.0)),
        ("0.1", "EPSG:4326", (0.1, 0.1)),
        ("0.1", "EPSG:3035", (0.1, 0.1)),
    ],
)
def test_gate_accepts_compatible(resolution, crs, expected):
    assert resolution_in_crs_units(resolution, crs) == pytest.approx(expected)


@pytest.mark.parametrize(
    "resolution, crs",
    [
        ("10 km", "EPSG:4326"),   # linear into geographic -- the headline case
        ("1000 m", "EPSG:4258"),
        ("0.1 deg", "EPSG:3035"),  # angular into projected -- the symmetric case
        ("30 arcsec", "EPSG:32632"),
    ],
)
def test_gate_rejects_incompatible(resolution, crs):
    with pytest.raises(ResolutionCRSMismatch):
        resolution_in_crs_units(resolution, crs)


def test_gate_runs_on_effective_crs(geo_ds):
    """No target_projection means the source CRS applies -- and must still fail."""
    from climdata.grid import reproject

    with pytest.raises(ResolutionCRSMismatch):
        reproject(geo_ds, target_resolution="10 km")


def test_error_message_is_actionable():
    with pytest.raises(ResolutionCRSMismatch) as excinfo:
        resolution_in_crs_units("10 km", "EPSG:4326", latitude=51.2)
    message = str(excinfo.value)
    assert '"10 km"' in message           # echoes what the caller typed
    assert "51.2N" in message             # the actual domain latitude
    assert "0.1 deg" in message           # both suggested fixes
    assert "EPSG:3035" in message


def test_non_metre_linear_crs():
    """EPSG:2263 is in US survey feet; conversion must be exact, not assumed metres."""
    x, y = resolution_in_crs_units("10 km", "EPSG:2263")
    assert x == pytest.approx(32808.333, abs=0.01)
    assert y == pytest.approx(32808.333, abs=0.01)


# --------------------------------------------------------------------------
# 4. CRS introspection
# --------------------------------------------------------------------------
@pytest.mark.parametrize("value", ["EPSG:4326", 4326, "epsg:4326"])
def test_parse_crs_round_trip(value):
    assert parse_crs(value).to_epsg() == 4326


def test_parse_crs_accepts_wkt_and_proj():
    import pyproj

    wkt = pyproj.CRS.from_epsg(3035).to_wkt()
    assert parse_crs(wkt).to_epsg() == 3035
    assert parse_crs("+proj=longlat +datum=WGS84 +no_defs").is_geographic


def test_parse_crs_rejects_garbage():
    with pytest.raises(ValueError, match="Could not interpret"):
        parse_crs("not a crs")


@pytest.mark.parametrize(
    "crs, kind, factor",
    [
        ("EPSG:4326", "angular", math.pi / 180.0),
        ("EPSG:3035", "linear", 1.0),
        ("EPSG:2263", "linear", 0.3048006096),
    ],
)
def test_crs_axis_unit(crs, kind, factor):
    got_kind, got_factor = crs_axis_unit(crs)
    assert got_kind == kind
    assert got_factor == pytest.approx(factor, rel=1e-9)


def test_axis_order_is_not_assumed(geo_ds):
    """EPSG:4326 lists latitude first; anisotropic resolution must not be swapped."""
    from climdata.grid import reproject

    out = reproject(geo_ds, "EPSG:4326", ("1.0 deg", "0.5 deg"), align=False)
    res_x, res_y = out.rio.resolution()
    assert abs(res_x) == pytest.approx(1.0)
    assert abs(res_y) == pytest.approx(0.5)


# --------------------------------------------------------------------------
# 5. The escape hatch
# --------------------------------------------------------------------------
def test_to_angular_matches_geod():
    from pyproj import Geod

    res = to_angular("10 km", latitude=50.0)
    geod = Geod(ellps="WGS84")
    lon_east, _, _ = geod.fwd(0.0, 50.0, 90.0, 10000.0)
    _, lat_north, _ = geod.fwd(0.0, 50.0, 0.0, 10000.0)
    assert res.x == pytest.approx(abs(lon_east), abs=1e-4)
    assert res.y == pytest.approx(abs(lat_north - 50.0), abs=1e-4)
    assert res.x == pytest.approx(0.1395, abs=1e-3)
    assert res.y == pytest.approx(0.0899, abs=1e-3)


def test_to_angular_requires_explicit_latitude():
    with pytest.raises(TypeError):
        to_angular("10 km")


def test_to_angular_output_passes_the_gate():
    res = to_angular("10 km", latitude=50.0)
    assert resolution_in_crs_units(res, "EPSG:4326")[0] == pytest.approx(0.1395, abs=1e-3)


def test_to_angular_rejects_angular_input():
    with pytest.raises(ValueError):
        to_angular("0.1 deg", latitude=50.0)


# --------------------------------------------------------------------------
# 6. Reprojection
# --------------------------------------------------------------------------
def test_known_answer_reprojection(geo_ds):
    """Values confirmed against the real stack while planning this module."""
    from climdata.grid import reproject

    out = reproject(geo_ds, "EPSG:3035", "10 km")
    assert out.rio.crs.to_epsg() == 3035
    assert (out.sizes["x"], out.sizes["y"]) == (79, 94)
    transform = out.rio.transform()
    assert transform.c == pytest.approx(3930000.0)
    assert transform.f == pytest.approx(3570000.0)
    assert out.rio.resolution() == pytest.approx((10000.0, -10000.0))


def test_alignment_invariant(geo_ds):
    """Differently-cropped inputs must land on identical cell centres."""
    from climdata.grid import reproject

    cropped = geo_ds.isel(lon=slice(3, None), lat=slice(2, None))
    full_out = reproject(geo_ds, "EPSG:3035", "10 km")
    crop_out = reproject(cropped, "EPSG:3035", "10 km")

    shared = set(np.round(full_out.x.values, 3)) & set(np.round(crop_out.x.values, 3))
    assert len(shared) == min(full_out.sizes["x"], crop_out.sizes["x"])


def test_align_false_preserves_extent(geo_ds):
    from climdata.grid import reproject

    out = reproject(geo_ds, "EPSG:4326", "0.25 deg", align=False)
    assert out.sizes["lon"] == geo_ds.sizes["lon"]
    assert out.sizes["lat"] == geo_ds.sizes["lat"]


def test_idempotent_on_native_grid(geo_ds):
    from climdata.grid import reproject

    out = reproject(geo_ds, "EPSG:4326", "0.25 deg", align=False).sortby("lat")
    got = out.tas.assign_coords(
        lat=np.round(out.lat.values, 6), lon=np.round(out.lon.values, 6)
    )
    want = geo_ds.tas.assign_coords(
        lat=np.round(geo_ds.lat.values, 6), lon=np.round(geo_ds.lon.values, 6)
    )
    assert np.nanmax(np.abs(got.values - want.values)) == pytest.approx(0.0, abs=1e-6)


def test_nodata_does_not_become_zero(geo_ds):
    """Without an explicit nodata, rasterio fills gaps with 0 and corrupts masks."""
    from climdata.grid import reproject

    masked = geo_ds.copy()
    masked["tas"] = masked.tas.where(masked.lon > 8.0)
    out = reproject(masked, "EPSG:3035", "20 km")
    assert not (out.tas.values == 0).any()
    assert np.isnan(out.tas.values).any()


def test_constant_field_survives_averaging(geo_ds):
    from climdata.grid import reproject

    const = xr.Dataset(
        {"t": (("lat", "lon"), np.full((geo_ds.sizes["lat"], geo_ds.sizes["lon"]), 5.0))},
        coords={"lat": geo_ds.lat, "lon": geo_ds.lon},
    ).rio.write_crs("EPSG:4326")
    out = reproject(const, "EPSG:4326", "1.0 deg", method="average")
    values = out.t.values
    assert np.allclose(values[~np.isnan(values)], 5.0)


def test_longitude_wrap(geo_ds):
    from climdata.grid import reproject

    shifted = geo_ds.assign_coords(lon=((geo_ds.lon + 360) % 360)).sortby("lon")
    shifted = shifted.rio.write_crs("EPSG:4326")
    out = reproject(shifted, "EPSG:4326", "0.5 deg")
    assert float(out.lon.min()) < 180.0
    assert float(out.lon.max()) < 180.0


def test_per_variable_method(geo_ds):
    from climdata.grid import reproject

    out = reproject(geo_ds, "EPSG:3035", "10 km", method={"pr": "average", "tas": "bilinear"})
    assert set(out.data_vars) == {"tas", "pr"}


def test_invalid_resampling_name(geo_ds):
    from climdata.grid import reproject

    with pytest.raises(ValueError, match="Unknown resampling method"):
        reproject(geo_ds, "EPSG:3035", "10 km", method="telepathy")


def test_like_matches_template(geo_ds):
    from climdata.grid import reproject

    template = reproject(geo_ds, "EPSG:3035", "10 km")
    out = reproject(geo_ds, like=template)
    assert out.sizes["x"] == template.sizes["x"]
    assert out.sizes["y"] == template.sizes["y"]


def test_like_with_explicit_target_is_contradictory(geo_ds):
    from climdata.grid import reproject

    template = reproject(geo_ds, "EPSG:3035", "10 km")
    with pytest.raises(ValueError, match="contradictory"):
        reproject(geo_ds, "EPSG:3035", "10 km", like=template)


def test_attributes_survive(geo_ds):
    from climdata.grid import reproject

    out = reproject(geo_ds, "EPSG:3035", "10 km")
    assert out.tas.attrs.get("units") == "K"
    assert out.pr.attrs.get("units") == "mm d-1"


def test_output_coord_names(geo_ds):
    from climdata.grid import reproject

    assert "lon" in reproject(geo_ds, "EPSG:4326", "0.5 deg").coords
    assert "x" in reproject(geo_ds, "EPSG:3035", "10 km").coords


# --------------------------------------------------------------------------
# 7. Integration: the full inference chain on HYRAS-shaped data
# --------------------------------------------------------------------------
def test_real_hyras_crs_from_esri_pe_string(hyras_v6_like):
    """Real HYRAS v6-1: EPSG:3035, recorded only as ESRI WKT on the data variable."""
    from climdata.grid import infer_src_crs

    assert "crs" not in hyras_v6_like.variables       # the grid_mapping target is gone
    assert infer_src_crs(hyras_v6_like).to_epsg() == 3035


def test_real_hyras_reprojects_to_geographic(hyras_v6_like):
    from climdata.grid import reproject

    out = reproject(hyras_v6_like, "EPSG:4326", "0.05 deg")
    assert out.rio.crs.to_epsg() == 4326
    assert {"lon", "lat"} <= set(out.coords)
    assert 7.0 < float(out.lon.mean()) < 9.0
    assert 46.5 < float(out.lat.mean()) < 48.5


def test_grid_mapping_variable_is_followed():
    """When the CF grid_mapping variable survives, its attributes are read."""
    from climdata.grid import infer_src_crs

    ds = xr.Dataset(
        {
            "tas": (("y", "x"), np.zeros((3, 3)), {"grid_mapping": "crs"}),
            "crs": ((), 0, {"epsg_code": "EPSG:3035"}),
        },
        coords={"x": [4.1e6, 4.2e6, 4.3e6], "y": [2.7e6, 2.8e6, 2.9e6]},
    )
    assert infer_src_crs(ds).to_epsg() == 3035


def test_hyras_like_inference_chain(hyras_like):
    """CRS comes only from ``attrs['crs_grid']``; dims are x/y, not lat/lon."""
    from climdata.grid import infer_src_crs, reproject

    assert infer_src_crs(hyras_like).to_epsg() == 31467

    out = reproject(hyras_like, "EPSG:4326", "0.05 deg")
    assert out.rio.crs.to_epsg() == 4326
    assert {"lon", "lat"} <= set(out.coords)
    assert 5.0 < float(out.lon.mean()) < 15.0
    assert 46.0 < float(out.lat.mean()) < 56.0


def test_missing_crs_raises_when_no_default():
    ds = xr.Dataset(
        {"v": (("y", "x"), np.zeros((3, 3)))},
        coords={"x": [1e6, 2e6, 3e6], "y": [1e6, 2e6, 3e6]},
    )
    from climdata.grid import infer_src_crs

    with pytest.raises(ValueError, match="Could not determine the source CRS"):
        infer_src_crs(ds, default=None)
