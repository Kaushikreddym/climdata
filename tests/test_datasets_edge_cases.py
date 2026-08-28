"""Offline edge cases for every dataset provider in :mod:`climdata.datasets`.

These never touch the network. What they exercise is the logic that runs
*before* any request goes out — time-range validation, filename parsing,
variable-name mapping, response-shape normalisation, raster decoding — plus the
error each provider raises when its optional dependency is absent.

That boundary matters: it is where a wrong configuration should be caught, and
catching it here costs nothing while catching it after a multi-gigabyte download
costs a great deal.

Live extraction is covered separately by ``test_datasets_extraction.py``.
"""

import importlib
from datetime import datetime

import numpy as np
import pytest
import xarray as xr
from omegaconf import OmegaConf

from climdata.datasets.HYRAS import find_nearest_xy, read_asc_file
from climdata.datasets.HOSTRADA import _first_hour_of_month, _last_hour_of_month


# ==========================================================================
# Every provider is importable, and names its missing dependency
# ==========================================================================
PROVIDERS = [
    ("climdata.datasets.AGRI_ISIMIP", "AGRI_ISIMIP"),
    ("climdata.datasets.CMIPCloud", "CMIPCloud"),
    ("climdata.datasets.CMIP_W5E5", "CMIPW5E5"),
    ("climdata.datasets.DWD", "DWDmirror"),
    ("climdata.datasets.ERA5", "ERA5Mirror"),
    ("climdata.datasets.HOSTRADA", "HOSTRADAmirror"),
    ("climdata.datasets.HYRAS", "HYRASmirror"),
    ("climdata.datasets.MSWX", "MSWXmirror"),
    ("climdata.datasets.NASAPOWER", "POWER"),
    ("climdata.datasets.NEXGDDP", "NEXGDDP"),
    ("climdata.datasets.W5E5", "W5E5"),
]


@pytest.mark.parametrize("module,cls", PROVIDERS)
def test_provider_imports_without_optional_dependencies(module, cls):
    """Importing climdata must not require every provider's extras."""
    assert hasattr(importlib.import_module(module), cls)


@pytest.mark.parametrize("module,cls", PROVIDERS)
def test_provider_class_is_documented(module, cls):
    """mkdocstrings renders these; an undocumented provider is an empty page."""
    doc = getattr(importlib.import_module(module), cls).__doc__
    assert doc and len(doc.strip()) > 40


# ==========================================================================
# Time-range validation (CMIPCloud, CMIPW5E5, NEXGDDP)
# ==========================================================================
def _cmip_cfg(experiment, start, end, **extra):
    return OmegaConf.create({
        "experiment_id": experiment,
        "source_id": "GFDL-ESM4",
        "table_id": "day",
        "variables": ["tas"],
        "time_range": {"start_date": start, "end_date": end},
        "data_dir": "/tmp",
        **extra,
    })


def _validator(module, cls, cfg):
    """Build the object far enough to call _validate_time_range, no network."""
    mod = importlib.import_module(module)
    obj = object.__new__(getattr(mod, cls))
    obj.cfg = cfg
    obj.experiment_id = cfg.experiment_id
    return obj


VALIDATING = [
    ("climdata.datasets.CMIPCloud", "CMIPCloud"),
    ("climdata.datasets.CMIP_W5E5", "CMIPW5E5"),
    ("climdata.datasets.NEXGDDP", "NEXGDDP"),
]


@pytest.mark.parametrize("module,cls", VALIDATING)
def test_historical_rejects_a_future_period(module, cls):
    """Historical runs stop at 2014; asking for 2050 would return nothing."""
    v = _validator(module, cls, _cmip_cfg("historical", "2050-01-01", "2050-12-31"))
    with pytest.raises(ValueError, match="[Tt]ime range"):
        v._validate_time_range()


@pytest.mark.parametrize("module,cls", VALIDATING)
def test_ssp_rejects_a_historical_period(module, cls):
    """SSP scenarios start in 2015."""
    v = _validator(module, cls, _cmip_cfg("ssp585", "1990-01-01", "1990-12-31"))
    with pytest.raises(ValueError, match="[Tt]ime range"):
        v._validate_time_range()


@pytest.mark.parametrize("module,cls", VALIDATING)
def test_a_valid_period_is_accepted(module, cls):
    v = _validator(module, cls, _cmip_cfg("historical", "1990-01-01", "1990-12-31"))
    v._validate_time_range()          # must not raise


@pytest.mark.parametrize("module,cls", VALIDATING)
def test_picontrol_skips_validation(module, cls):
    """Pre-industrial control runs are long and unbounded; no rule applies."""
    v = _validator(module, cls, _cmip_cfg("picontrol", "1000-01-01", "3000-12-31"))
    v._validate_time_range()


@pytest.mark.parametrize("module,cls", VALIDATING)
def test_an_unknown_experiment_skips_validation(module, cls):
    """Unrecognised experiments must not be blocked by a stale rule table."""
    v = _validator(module, cls, _cmip_cfg("some-new-mip", "1000-01-01", "3000-12-31"))
    v._validate_time_range()


@pytest.mark.parametrize("module,cls", VALIDATING)
def test_a_period_overrunning_one_end_warns_but_proceeds(module, cls, capsys):
    """Partial overlap is usable, so it warns rather than raising."""
    v = _validator(module, cls, _cmip_cfg("historical", "1990-01-01", "2020-12-31"))
    v._validate_time_range()
    assert "Warning" in capsys.readouterr().out


@pytest.mark.parametrize("module,cls", VALIDATING)
def test_boundary_year_is_inclusive(module, cls):
    """2014 is the last historical year and must be accepted, not rejected."""
    v = _validator(module, cls, _cmip_cfg("historical", "2014-01-01", "2014-12-31"))
    v._validate_time_range()


# ==========================================================================
# Filename year-range overlap (W5E5, CMIPW5E5)
# ==========================================================================
OVERLAP = [
    ("climdata.datasets.W5E5", "W5E5", "w5e5v2.0_obsclim_tas_global_daily_1979_1989.nc"),
    ("climdata.datasets.CMIP_W5E5", "CMIPW5E5",
     "gfdl-esm4_r1i1p1f1_ssp585_tas_global_daily_2015_2024.nc"),
]


@pytest.mark.parametrize("module,cls,filename", OVERLAP)
@pytest.mark.parametrize("start,end,expected", [
    ("1979-01-01", "1989-12-31", True),    # exactly the file's span
    ("1985-01-01", "1985-12-31", True),    # wholly inside
    ("1970-01-01", "1980-01-01", True),    # overlaps the start
    ("1988-01-01", "1995-01-01", True),    # overlaps the end
    ("1960-01-01", "1970-01-01", False),   # entirely before
    ("2050-01-01", "2060-01-01", False),   # entirely after
])
def test_date_range_overlap(module, cls, filename, start, end, expected):
    mod = importlib.import_module(module)
    obj = object.__new__(getattr(mod, cls))
    # Shift the CMIP_W5E5 window to that file's own span.
    if "2015_2024" in filename:
        shift = 36
        start = f"{int(start[:4]) + shift}{start[4:]}"
        end = f"{int(end[:4]) + shift}{end[4:]}"
    got = obj._is_file_in_date_range(
        filename, datetime.fromisoformat(start), datetime.fromisoformat(end)
    )
    assert got is expected


@pytest.mark.parametrize("module,cls,_f", OVERLAP)
def test_an_unparseable_filename_is_kept(module, cls, _f):
    """Better to download a file we cannot date than to silently drop it."""
    mod = importlib.import_module(module)
    obj = object.__new__(getattr(mod, cls))
    assert obj._is_file_in_date_range(
        "no_years_here.nc", datetime(2000, 1, 1), datetime(2001, 1, 1)
    ) is True


# ==========================================================================
# Variable-name mapping (W5E5, CMIPW5E5)
# ==========================================================================
@pytest.mark.parametrize("module,cls", [
    ("climdata.datasets.W5E5", "W5E5"),
    ("climdata.datasets.CMIP_W5E5", "CMIPW5E5"),
])
@pytest.mark.parametrize("cf_name,expected", [
    ("tas", "tas"),
    ("tasmax", "tasmax"),
    ("pr", "pr"),
    ("sfcWind", "sfcwind"),          # the one that is not identity
    ("not_a_variable", "not_a_variable"),   # unknown passes through
])
def test_variable_name_mapping(module, cls, cf_name, expected):
    mod = importlib.import_module(module)
    obj = object.__new__(getattr(mod, cls))
    assert obj._map_variable_name(cf_name) == expected


# ==========================================================================
# ISIMIP response-shape normalisation (W5E5, CMIPW5E5)
# ==========================================================================
@pytest.mark.parametrize("module,cls", [
    ("climdata.datasets.W5E5", "W5E5"),
    ("climdata.datasets.CMIP_W5E5", "CMIPW5E5"),
])
@pytest.mark.parametrize("response,expected_len", [
    (None, 0),
    ([], 0),
    ([{"name": "a"}], 1),
    ({"results": [{"name": "a"}, {"name": "b"}]}, 2),
    ({"results": None}, 0),
    ({"count": 0}, 0),
    ("unexpected string", 0),
])
def test_dataset_response_normalisation(module, cls, response, expected_len):
    """The client returns a bare list from some endpoints, an envelope from others."""
    mod = importlib.import_module(module)
    obj = object.__new__(getattr(mod, cls))
    assert len(obj._normalize_dataset_results(response)) == expected_len


@pytest.mark.parametrize("module,cls", [
    ("climdata.datasets.W5E5", "W5E5"),
    ("climdata.datasets.CMIP_W5E5", "CMIPW5E5"),
])
@pytest.mark.parametrize("record,expected_len", [
    ({"files": [{"name": "a.nc"}, {"name": "b.nc"}]}, 2),
    ({"file": {"name": "a.nc"}}, 1),
    ({}, 0),
    ({"files": None}, 0),
])
def test_file_list_extraction(module, cls, record, expected_len):
    mod = importlib.import_module(module)
    obj = object.__new__(getattr(mod, cls))
    assert len(obj._extract_files_list(record)) == expected_len


# ==========================================================================
# HYRAS ASCII raster decoding
# ==========================================================================
def _write_asc(path, values, nodata=-999.0, cellsize=1000.0):
    rows = "\n".join(" ".join(str(v) for v in row) for row in values)
    path.write_text(
        f"ncols {len(values[0])}\n"
        f"nrows {len(values)}\n"
        f"xllcorner 3400000.0\n"
        f"yllcorner 5200000.0\n"
        f"cellsize {cellsize}\n"
        f"NODATA_value {nodata}\n"
        f"{rows}\n"
    )
    return path


def test_read_asc_parses_the_header_geometry(tmp_path):
    da = read_asc_file(str(_write_asc(tmp_path / "a.asc", [[1, 2, 3], [4, 5, 6]])))

    assert da.dims == ("y", "x")
    assert da.shape == (2, 3)
    assert da.attrs["cellsize"] == 1000.0
    assert da.attrs["crs_grid"] == "EPSG:31467"
    assert da.attrs["xllcorner"] == 3400000.0


def test_read_asc_y_is_descending(tmp_path):
    """ASCII rasters store rows north-to-south; y must match that order."""
    da = read_asc_file(str(_write_asc(tmp_path / "a.asc", [[1, 2], [3, 4]])))
    assert da["y"].values[0] > da["y"].values[-1]
    # first row of the file is the northernmost
    assert da.isel(y=0).values.tolist() == [1.0, 2.0]


def test_read_asc_converts_nodata_to_nan(tmp_path):
    da = read_asc_file(str(_write_asc(tmp_path / "a.asc", [[1, -999.0], [3, 4]])))
    assert np.isnan(da.values[0, 1])
    assert not np.isnan(da.values[0, 0])


@pytest.mark.parametrize("varname,factor", [
    ("evpot", 0.1),        # published scaled by ten
    ("soilTemp", 0.1),     # published scaled by ten
    ("soilMoist", 1.0),    # not scaled
    (None, 1.0),           # unknown variable: leave alone
])
def test_read_asc_applies_only_the_documented_scaling(tmp_path, varname, factor):
    """Getting varname wrong yields values off by ten, so this must be exact."""
    da = read_asc_file(str(_write_asc(tmp_path / "a.asc", [[100, 200]])), varname=varname)
    assert da.values[0, 0] == pytest.approx(100 * factor)


def test_read_asc_records_units_when_given(tmp_path):
    da = read_asc_file(str(_write_asc(tmp_path / "a.asc", [[1]])), units="mm")
    assert da.attrs["units"] == "mm"


def test_read_asc_accepts_a_file_object(tmp_path):
    """Members are read straight out of a .tgz, never extracted to disk."""
    import io

    path = _write_asc(tmp_path / "a.asc", [[1, 2]])
    da = read_asc_file(io.BytesIO(path.read_bytes()))
    assert da.shape == (1, 2)


# ==========================================================================
# HYRAS nearest-cell lookup
# ==========================================================================
@pytest.fixture
def curvilinear_ds():
    """A 2-D lat/lon grid, as HYRAS NetCDF files carry."""
    y, x = np.arange(3), np.arange(4)
    lon2d, lat2d = np.meshgrid(np.linspace(6.0, 9.0, 4), np.linspace(50.0, 52.0, 3))
    return xr.Dataset(
        {"tas": (("y", "x"), np.zeros((3, 4)))},
        coords={"y": y, "x": x,
                "lat": (("y", "x"), lat2d), "lon": (("y", "x"), lon2d)},
    )


def test_find_nearest_xy_hits_an_exact_grid_point(curvilinear_ds):
    iy, ix = find_nearest_xy(curvilinear_ds, 51.0, 7.0)
    assert curvilinear_ds["lat"].values[iy, ix] == pytest.approx(51.0)
    assert curvilinear_ds["lon"].values[iy, ix] == pytest.approx(7.0)


def test_find_nearest_xy_snaps_an_off_grid_point(curvilinear_ds):
    iy, ix = find_nearest_xy(curvilinear_ds, 50.9, 7.1)
    assert (iy, ix) == (1, 1)


def test_find_nearest_xy_clamps_outside_the_domain(curvilinear_ds):
    """A point far outside still returns a valid corner index, never an error."""
    iy, ix = find_nearest_xy(curvilinear_ds, 0.0, 0.0)
    assert 0 <= iy < curvilinear_ds.sizes["y"]
    assert 0 <= ix < curvilinear_ds.sizes["x"]


# ==========================================================================
# HOSTRADA filename timestamps
# ==========================================================================
@pytest.mark.parametrize("year,month,expected", [
    (2020, 1, "2020013123"),
    (2020, 2, "2020022923"),     # leap year
    (2021, 2, "2021022823"),     # non-leap
    (2020, 4, "2020043023"),     # 30-day month
    (2020, 12, "2020123123"),
])
def test_last_hour_of_month(year, month, expected):
    assert _last_hour_of_month(year, month) == expected


@pytest.mark.parametrize("year,month,expected", [
    (2020, 1, "2020010100"),
    (2020, 12, "2020120100"),
])
def test_first_hour_of_month(year, month, expected):
    assert _first_hour_of_month(year, month) == expected


# ==========================================================================
# NASA POWER
# ==========================================================================
def test_power_load_without_fetch():
    """load() before fetch() must fail on the missing response, not silently."""
    from climdata.datasets.NASAPOWER import POWER

    power = POWER(OmegaConf.create({"lat": 52.5, "lon": 13.4}))
    assert power.raw is None
    with pytest.raises(TypeError):
        power.load()


def test_power_extract_without_bounds_is_a_noop():
    """The generic workflow calls extract() unconditionally."""
    from climdata.datasets.NASAPOWER import POWER

    power = POWER(OmegaConf.create({}))
    power.ds = xr.Dataset({"tas": ("time", [1.0, 2.0])},
                          coords={"time": [0, 1]})
    before = power.ds
    power.extract()
    assert power.ds is before


# ==========================================================================
# NEXGDDP static metadata
# ==========================================================================
def test_nexgddp_advertises_models_and_variables():
    from climdata.datasets.NEXGDDP import NEXGDDP

    assert "GFDL-ESM4" in NEXGDDP.AVAILABLE_MODELS
    assert "tasmax" in NEXGDDP.AVAILABLE_VARIABLES
    assert "historical" in NEXGDDP.AVAILABLE_EXPERIMENTS


def test_nexgddp_warns_rather_than_rejecting_unknown_names(capsys):
    """NASA adds models faster than the hardcoded list; rejecting would block them."""
    from climdata.datasets.NEXGDDP import NEXGDDP

    obj = object.__new__(NEXGDDP)
    obj.cfg = OmegaConf.create({"variables": ["not_a_variable"]})
    obj.source_id = "NOT-A-MODEL"
    obj.experiment_id = "not-an-experiment"
    obj._validate_inputs()                     # must not raise

    out = capsys.readouterr().out
    assert "may not be available" in out


# ==========================================================================
# Packaging
# ==========================================================================
def test_every_subpackage_is_shippable():
    """climdata.datasets and climdata.extremes once lacked __init__.py, so
    setuptools left them out of the wheel entirely and `pip install climdata`
    produced a package whose first import failed."""
    from setuptools import find_packages

    shipped = find_packages(include=["climdata*"],
                            exclude=["docs*", "climdata.sdba_usecase*"])
    for required in ("climdata", "climdata.datasets", "climdata.extremes",
                     "climdata.explore", "climdata.grid", "climdata.impute",
                     "climdata.sdba", "climdata.utils", "climdata.LLM"):
        assert required in shipped, f"{required} would not be installed"


def test_the_datasets_package_exports_every_provider():
    import climdata.datasets as pkg

    for _, cls in PROVIDERS:
        assert cls in pkg.__all__
        assert hasattr(pkg, cls)


def test_research_scratch_is_not_shipped():
    """sdba_usecase holds notebooks and a geopackage; keep it out of the wheel."""
    from setuptools import find_packages

    shipped = find_packages(include=["climdata*"],
                            exclude=["docs*", "climdata.sdba_usecase*"])
    assert not any(p.startswith("climdata.sdba_usecase") for p in shipped)
