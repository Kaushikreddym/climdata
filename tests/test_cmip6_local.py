"""Edge cases for reading a local CMIP6 mirror (``climdata.utils.cmip6_local``).

Real mirrors are not clean DRS trees, and every case here is one that was
actually hit on disk: macOS sidecar files that look like NetCDF but are not,
version directories with the ESGF host glued on and split across extra levels,
SSP extension files whose dates fall outside pandas' nanosecond range, and the
same variable published twice — once mirrored from a second host, once as a
genuinely disjoint time slice.

The theme is that the indexer must not lose data quietly. Dropping a file is
fine when it is not CMIP6; silently dropping half an experiment because a
newer version directory happens to cover different years is not.
"""

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from climdata.utils.cmip6_local import (
    _select_files,
    cmip6_coverage,
    index_cmip6,
    open_cmip6,
)


# ==========================================================================
# Fixtures
# ==========================================================================

def _write_nc(path, var, start, periods):
    """A minimal but genuinely openable DRS file at ``path``."""
    path.parent.mkdir(parents=True, exist_ok=True)
    time = pd.date_range(start, periods=periods, freq="D")
    ds = xr.Dataset(
        {var: (("time", "lat", "lon"),
               np.arange(periods * 2 * 3, dtype="float32").reshape(periods, 2, 3))},
        coords={"time": time, "lat": [10.0, 20.0], "lon": [0.0, 120.0, 240.0]},
    )
    ds.to_netcdf(path)


@pytest.fixture
def mirror(tmp_path):
    """A small mirror carrying every wart the real ones have."""
    root = tmp_path / "cmip6_data" / "CMIP6"

    # Normal, well-formed model.
    good = root / "CMIP/INST/GOOD-MODEL/historical/r1i1p1f1/day/tas/gn/v20200101"
    _write_nc(good / "tas_day_GOOD-MODEL_historical_r1i1p1f1_gn_20000101-20000110.nc",
              "tas", "2000-01-01", 10)
    _write_nc(good / "tas_day_GOOD-MODEL_historical_r1i1p1f1_gn_20000111-20000120.nc",
              "tas", "2000-01-11", 10)
    # AppleDouble sidecar: not NetCDF, must never be indexed.
    (good / "._tas_day_GOOD-MODEL_historical_r1i1p1f1_gn_20000121-20000130.nc"
     ).write_bytes(b"\x00\x05\x16\x07not netcdf")
    # Not DRS-named at all (a postproc/ style file).
    _write_nc(root / "postproc/tas/tas_GOOD-MODEL_2000-2001.nc", "tas", "2000-01-01", 5)

    # Version directory with the ESGF host appended and split across levels.
    messy = (root / "DAMIP/INST/MESSY-MODEL/hist-GHG/r1i1p1f1/day/pr/gn"
             / "v20190903esgf" / "ceda" / "ac" / "uk")
    _write_nc(messy / "pr_day_MESSY-MODEL_hist-GHG_r1i1p1f1_gn_20000101-20000110.nc",
              "pr", "2000-01-01", 10)
    # Same version, same slice, mirrored from a second host.
    messy2 = (root / "DAMIP/INST/MESSY-MODEL/hist-GHG/r1i1p1f1/day/pr/gn"
              / "v20190903esgf3" / "dkrz" / "de")
    _write_nc(messy2 / "pr_day_MESSY-MODEL_hist-GHG_r1i1p1f1_gn_20000101-20000110.nc",
              "pr", "2000-01-01", 10)

    return root.parent


# ==========================================================================
# Indexing
# ==========================================================================

def test_index_skips_sidecars_and_non_drs_files(mirror):
    """``._`` sidecars and postproc/ files are not CMIP6 and must not be rows."""
    df = index_cmip6(mirror, verbose=False)
    names = [p.rsplit("/", 1)[-1] for p in df["path"]]
    assert not any(name.startswith("._") for name in names)
    assert not any("postproc" in p for p in df["path"])
    assert set(df["source_id"]) == {"GOOD-MODEL", "MESSY-MODEL"}


def test_index_recovers_facets_from_a_garbled_version_dir(mirror):
    """The host glued onto the version must not shift the DRS facets."""
    df = index_cmip6(mirror, models="MESSY-MODEL", verbose=False)
    row = df.iloc[0]
    assert row["activity_id"] == "DAMIP"
    assert row["institution_id"] == "INST"
    assert row["experiment_id"] == "hist-GHG"
    assert row["grid_label"] == "gn"
    assert row["version"] == "v20190903"       # host suffix stripped


def test_index_filters_are_applied(mirror):
    df = index_cmip6(mirror, variables="pr", verbose=False)
    assert set(df["variable_id"]) == {"pr"}
    assert index_cmip6(mirror, experiments="ssp585", verbose=False).empty


def test_require_all_experiments_drops_incomplete_models(mirror):
    """GOOD-MODEL has only historical, MESSY-MODEL only hist-GHG, so asking for
    both must leave nothing rather than silently keeping half of each."""
    df = index_cmip6(mirror, experiments=["historical", "hist-GHG"],
                     require_all_experiments=True, verbose=False)
    assert df.empty


def test_require_all_experiments_needs_experiments(mirror):
    with pytest.raises(ValueError, match="experiments"):
        index_cmip6(mirror, require_all_experiments=True, verbose=False)


def test_missing_root_is_an_error(tmp_path):
    with pytest.raises(FileNotFoundError):
        index_cmip6(tmp_path / "nope", verbose=False)


def test_empty_result_keeps_the_schema(tmp_path):
    """An empty catalogue still has to be a usable DataFrame, or every caller
    needs its own emptiness special-case."""
    df = index_cmip6(tmp_path, verbose=False)
    assert df.empty
    for col in ("path", "source_id", "experiment_id", "variable_id", "years"):
        assert col in df.columns


# ==========================================================================
# Date stamps
# ==========================================================================

@pytest.mark.parametrize("stamp, start, end", [
    ("20000101-20001231", "2000-01-01", "2000-12-31"),   # daily
    ("200001-200012",     "2000-01-01", "2000-12-31"),   # monthly, end-of-month
    ("2000-2001",         "2000-01-01", "2001-12-31"),   # yearly, end-of-year
])
def test_time_range_ends_on_the_last_day_covered(stamp, start, end):
    from climdata.utils.cmip6_local import _parse_time_range
    got_start, got_end = _parse_time_range(stamp)
    assert got_start == pd.Timestamp(start)
    assert got_end == pd.Timestamp(end)


def test_dates_beyond_the_nanosecond_range_survive():
    """SSP extension files run to 2300, past pandas' default Timestamp bound —
    they must parse, not become NaT and silently lose their coverage."""
    from climdata.utils.cmip6_local import _parse_time_range
    start, end = _parse_time_range("22510101-23001231")
    assert start == pd.Timestamp("2251-01-01")
    assert end == pd.Timestamp("2300-12-31")


def test_unparseable_time_range_is_not_fatal():
    from climdata.utils.cmip6_local import _parse_time_range
    assert _parse_time_range("garbage") == (pd.NaT, pd.NaT)
    assert _parse_time_range(None) == (pd.NaT, pd.NaT)


# ==========================================================================
# Resolving duplicates before opening
# ==========================================================================

def _index_frame(rows):
    """Build a catalogue-shaped frame from (version, start, end, grid, member)."""
    return pd.DataFrame([
        {"source_id": "M", "experiment_id": "E", "variable_id": "pr",
         "member_id": member, "grid_label": grid, "version": version,
         "start_date": pd.Timestamp(start), "end_date": pd.Timestamp(end),
         "path": f"/{version}/{grid}/{member}/{start}.nc"}
        for version, start, end, grid, member in rows
    ])


def test_same_slice_mirrored_twice_is_opened_once():
    df = _index_frame([("v1", "2000-01-01", "2000-12-31", "gn", "r1i1p1f1"),
                       ("v2", "2000-01-01", "2000-12-31", "gn", "r1i1p1f1")])
    got = _select_files(df, "M", "E", "pr", None, None)
    assert len(got) == 1
    assert got["version"].iloc[0] == "v2"      # newest version of that slice wins


def test_disjoint_slices_from_different_versions_are_all_kept():
    """The regression this guards: taking only the newest version would drop
    2000-2100 and leave a store that starts in 2200."""
    df = _index_frame([("v1", "2000-01-01", "2100-12-31", "gn", "r1i1p1f1"),
                       ("v2", "2200-01-01", "2250-12-31", "gn", "r1i1p1f1")])
    got = _select_files(df, "M", "E", "pr", None, None)
    assert len(got) == 2
    assert set(got["version"]) == {"v1", "v2"}


def test_one_grid_and_one_member_are_chosen():
    df = _index_frame([("v1", "2000-01-01", "2000-12-31", "gn", "r1i1p1f1"),
                       ("v1", "2000-01-01", "2000-12-31", "gr2", "r1i1p1f1"),
                       ("v1", "2000-01-01", "2000-12-31", "gn", "r2i1p1f1")])
    got = _select_files(df, "M", "E", "pr", None, None)
    assert set(got["grid_label"]) == {"gn"}        # native preferred
    assert set(got["member_id"]) == {"r1i1p1f1"}   # r1i1p1f1 preferred


def test_explicit_member_and_grid_win():
    df = _index_frame([("v1", "2000-01-01", "2000-12-31", "gn", "r1i1p1f1"),
                       ("v1", "2000-01-01", "2000-12-31", "gr2", "r2i1p1f1")])
    got = _select_files(df, "M", "E", "pr", "r2i1p1f1", "gr2")
    assert got["grid_label"].tolist() == ["gr2"]
    assert got["member_id"].tolist() == ["r2i1p1f1"]


# ==========================================================================
# Opening
# ==========================================================================

def test_open_concatenates_and_normalises_the_grid(mirror):
    ds = open_cmip6("GOOD-MODEL", "historical", "tas", base_dir=mirror, verbose=False)
    assert ds.sizes["time"] == 20                       # both files, concatenated
    assert not ds.get_index("time").duplicated().any()
    assert float(ds.lon.max()) <= 180                   # 0-360 -> -180..180
    assert (ds.lon.diff("lon") > 0).all()               # and sorted
    assert ds.attrs["source_id"] == "GOOD-MODEL"


def test_open_dedupes_a_doubly_mirrored_file(mirror):
    """Both host copies hold the same 10 days; opening both would double them."""
    ds = open_cmip6("MESSY-MODEL", "hist-GHG", "pr", base_dir=mirror, verbose=False)
    assert ds.sizes["time"] == 10


def test_open_spans_the_datetime64_boundary(tmp_path):
    """An SSP whose extension runs past 2262 decodes its early files as numpy
    datetimes and its late ones as cftime; combining the two raises unless the
    whole call is decoded the same way.
    """
    d = (tmp_path / "CMIP6/ScenarioMIP/INST/EXT-MODEL/ssp126/r1i1p1f1/day/tas/gn/v20200101")
    _write_nc(d / "tas_day_EXT-MODEL_ssp126_r1i1p1f1_gn_20150101-20150110.nc",
              "tas", "2015-01-01", 10)
    d.mkdir(parents=True, exist_ok=True)
    time = xr.cftime_range("2290-01-01", periods=10, freq="D",
                           calendar="proleptic_gregorian")
    xr.Dataset(
        {"tas": (("time", "lat", "lon"),
                 np.zeros((10, 2, 3), dtype="float32"))},
        coords={"time": time, "lat": [10.0, 20.0], "lon": [0.0, 120.0, 240.0]},
    ).to_netcdf(d / "tas_day_EXT-MODEL_ssp126_r1i1p1f1_gn_22900101-22900110.nc")

    ds = open_cmip6("EXT-MODEL", "ssp126", "tas", base_dir=tmp_path, verbose=False)
    assert ds.sizes["time"] == 20
    assert str(ds.get_index("time")[-1])[:10] == "2290-01-10"


def test_open_handles_a_non_standard_calendar(tmp_path):
    """noleap/360_day models have no datetime64 representation at all, so the
    auto-detection must leave decoding to xarray rather than forcing
    datetime64 on the grounds that the dates are in range."""
    d = tmp_path / "CMIP6/CMIP/INST/NOLEAP-MODEL/historical/r1i1p1f1/day/tas/gn/v20200101"
    d.mkdir(parents=True)
    time = xr.cftime_range("2000-01-01", periods=10, freq="D", calendar="noleap")
    xr.Dataset(
        {"tas": (("time", "lat", "lon"), np.zeros((10, 2, 3), dtype="float32"))},
        coords={"time": time, "lat": [10.0, 20.0], "lon": [0.0, 120.0, 240.0]},
    ).to_netcdf(d / "tas_day_NOLEAP-MODEL_historical_r1i1p1f1_gn_20000101-20000110.nc")

    ds = open_cmip6("NOLEAP-MODEL", "historical", "tas", base_dir=tmp_path,
                    verbose=False)
    assert ds.sizes["time"] == 10
    assert ds.get_index("time")[0].calendar == "noleap"


def test_open_crops_in_space_and_time(mirror):
    ds = open_cmip6("GOOD-MODEL", "historical", "tas", base_dir=mirror,
                    bounds={"lat_min": 5, "lat_max": 15,
                            "lon_min": -180, "lon_max": 180},
                    time_range=("2000-01-05", "2000-01-15"), verbose=False)
    assert ds.sizes["lat"] == 1
    assert ds.sizes["time"] == 11


def test_open_warns_but_returns_what_exists(mirror):
    with pytest.warns(UserWarning, match="no files"):
        ds = open_cmip6("GOOD-MODEL", "historical", ["tas", "nosuchvar"],
                        base_dir=mirror, verbose=False)
    assert list(ds.data_vars) == ["tas"]


def test_open_strict_refuses_a_partial_result(mirror):
    with pytest.raises(FileNotFoundError, match="nosuchvar"):
        open_cmip6("GOOD-MODEL", "historical", ["tas", "nosuchvar"],
                   base_dir=mirror, strict=True, verbose=False)


def test_open_with_nothing_present_is_an_error(mirror):
    with pytest.raises(FileNotFoundError):
        open_cmip6("GOOD-MODEL", "ssp585", "tas", base_dir=mirror, verbose=False)


# ==========================================================================
# Coverage summary
# ==========================================================================

def test_coverage_reports_the_best_choice_not_the_sum():
    """Two mirrored copies of the same decade is ten years of data, not twenty."""
    df = _index_frame([("v1", "2000-01-01", "2009-12-31", "gn", "r1i1p1f1"),
                       ("v2", "2000-01-01", "2009-12-31", "gr2", "r1i1p1f1")])
    df["years"] = (df["end_date"] - df["start_date"]).dt.days / 365.25
    cov = cmip6_coverage(df)
    assert cov.loc[("M", "E"), "pr"] == pytest.approx(9.99, abs=0.02)


def test_coverage_of_an_empty_index_is_empty():
    assert cmip6_coverage(pd.DataFrame()).empty


def test_coverage_rejects_an_unknown_value(mirror):
    df = index_cmip6(mirror, verbose=False)
    with pytest.raises(ValueError, match="years"):
        cmip6_coverage(df, value="nope")
