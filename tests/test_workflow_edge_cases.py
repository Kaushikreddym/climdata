"""Edge cases for :class:`climdata.ClimData` (``ClimateExtractor``).

These run entirely offline. Nothing here extracts from a provider — that needs
network access and is covered by ``test_datasets_extraction.py``. What is
exercised instead is everything around extraction: configuration composition,
AOI parsing, the upload paths, filename generation, format conversion, and the
error messages a user actually hits when a request is malformed.

Several tests pin behaviour that was previously broken, and say so, so a
regression is recognisable rather than merely red.
"""

import json

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from climdata import ClimData
from climdata.grid import ResolutionCRSMismatch
from climdata.utils.wrapper_workflow import (
    WorkflowResult,
    _EXTRACTABLE,
    _WORKFLOW_ACTIONS,
)


BASE = [
    "dataset=mswx",
    "lat=52.5",
    "lon=13.4",
    "variables=[tasmax]",
    "time_range.start_date=2014-01-01",
    "time_range.end_date=2014-12-31",
]


# --------------------------------------------------------------------------
# Fixtures
# --------------------------------------------------------------------------
@pytest.fixture
def extractor(tmp_path):
    """A configured extractor writing into a temporary directory."""
    return ClimData(overrides=BASE + [f"output.out_dir={tmp_path}"])


@pytest.fixture
def gridded_ds():
    """A tiny 3-day, 2x2 gridded dataset with units and a source attribute."""
    time = pd.date_range("2014-01-01", periods=3, freq="D")
    lat = np.array([52.0, 52.5])
    lon = np.array([13.0, 13.5])
    rng = np.random.default_rng(0)
    ds = xr.Dataset(
        {"tasmax": (("time", "lat", "lon"), rng.normal(15, 3, (3, 2, 2)))},
        coords={"time": time, "lat": lat, "lon": lon},
    )
    ds.tasmax.attrs["units"] = "degC"
    ds.attrs["source"] = "synthetic"
    return ds


# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
def test_config_loads_and_applies_overrides(extractor):
    assert extractor.cfg.dataset == "mswx"
    assert extractor.cfg.lat == 52.5
    assert list(extractor.cfg.variables) == ["tasmax"]


def test_unknown_hydra_key_is_rejected():
    """A typo in an override must fail loudly, not be silently ignored."""
    with pytest.raises(Exception):
        ClimData(overrides=BASE + ["not_a_real_key=1"])


def test_repeated_construction_is_independent(tmp_path):
    """Hydra's GlobalHydra is process-wide; two extractors must not share state."""
    a = ClimData(overrides=BASE + [f"output.out_dir={tmp_path}"])
    b = ClimData(overrides=BASE + ["lat=48.0", f"output.out_dir={tmp_path}"])
    assert a.cfg.lat == 52.5
    assert b.cfg.lat == 48.0


# --------------------------------------------------------------------------
# AOI preprocessing
# --------------------------------------------------------------------------
def _aoi(extractor, geojson):
    extractor.cfg.aoi = json.dumps(geojson)
    return extractor.preprocess_aoi(extractor.cfg)


def test_aoi_point_sets_lat_lon(extractor):
    cfg = _aoi(extractor, {"type": "Point", "coordinates": [13.4, 52.5]})
    assert cfg.lat == 52.5 and cfg.lon == 13.4      # GeoJSON is (lon, lat)
    assert cfg.bounds is None


def test_aoi_polygon_sets_bounds_and_clears_point(extractor):
    poly = {
        "type": "Polygon",
        "coordinates": [[[9.0, 51.0], [11.0, 51.0], [11.0, 53.0], [9.0, 53.0], [9.0, 51.0]]],
    }
    cfg = _aoi(extractor, poly)
    assert cfg.region == "custom"
    assert cfg.lat is None and cfg.lon is None
    b = cfg.bounds["custom"]
    assert (b["lat_min"], b["lat_max"]) == (51.0, 53.0)
    assert (b["lon_min"], b["lon_max"]) == (9.0, 11.0)


def test_aoi_feature_collection_uses_first_feature(extractor):
    fc = {
        "type": "FeatureCollection",
        "features": [
            {"type": "Feature", "properties": {},
             "geometry": {"type": "Point", "coordinates": [13.4, 52.5]}},
            {"type": "Feature", "properties": {},
             "geometry": {"type": "Point", "coordinates": [0.0, 0.0]}},
        ],
    }
    cfg = _aoi(extractor, fc)
    assert (cfg.lat, cfg.lon) == (52.5, 13.4)


def test_aoi_bare_feature(extractor):
    feat = {"type": "Feature", "properties": {},
            "geometry": {"type": "Point", "coordinates": [13.4, 52.5]}}
    assert _aoi(extractor, feat).lat == 52.5


def test_aoi_invalid_json_raises(extractor):
    extractor.cfg.aoi = "{not valid json"
    with pytest.raises(ValueError, match="Invalid AOI JSON"):
        extractor.preprocess_aoi(extractor.cfg)


def test_aoi_unsupported_geometry_type_raises(extractor):
    line = {"type": "LineString", "coordinates": [[9.0, 51.0], [11.0, 53.0]]}
    with pytest.raises(ValueError, match="Unknown geometry type"):
        _aoi(extractor, line)


def test_aoi_without_type_key_raises(extractor):
    extractor.cfg.aoi = json.dumps({"coordinates": [13.4, 52.5]})
    with pytest.raises(ValueError, match="Unsupported AOI format"):
        extractor.preprocess_aoi(extractor.cfg)


def test_aoi_none_is_a_noop(extractor):
    extractor.cfg.aoi = None
    assert extractor.preprocess_aoi(extractor.cfg).lat == 52.5


# --------------------------------------------------------------------------
# Uploads
# --------------------------------------------------------------------------
def test_upload_netcdf_missing_file(extractor):
    with pytest.raises(FileNotFoundError):
        extractor.upload_netcdf("/nonexistent/file.nc")


def test_upload_netcdf_round_trip(extractor, gridded_ds, tmp_path):
    """Regression: _gen_fn is keyword-only, and used to be called positionally,
    so every upload raised TypeError before reaching the caller."""
    path = tmp_path / "in.nc"
    gridded_ds.to_netcdf(path)

    ds = extractor.upload_netcdf(str(path))

    assert "tasmax" in ds.data_vars
    assert extractor.current_ds is ds
    assert extractor.filename_nc.endswith(".nc")


def test_upload_csv_missing_file(extractor):
    with pytest.raises(FileNotFoundError):
        extractor.upload_csv("/nonexistent/file.csv")


def test_upload_csv_without_coordinate_columns(extractor, tmp_path):
    path = tmp_path / "bad.csv"
    pd.DataFrame({
        "time": pd.date_range("2014-01-01", periods=2),
        "variable": ["tasmax"] * 2,
        "value": [1.0, 2.0],
    }).to_csv(path, index=False)

    with pytest.raises(ValueError, match="lat.*lon"):
        extractor.upload_csv(str(path))


def test_upload_csv_round_trip(extractor, tmp_path):
    path = tmp_path / "in.csv"
    pd.DataFrame({
        "time": list(pd.date_range("2014-01-01", periods=2)) * 2,
        "lat": [52.0, 52.0, 52.5, 52.5],
        "lon": [13.0, 13.0, 13.5, 13.5],
        "variable": ["tasmax"] * 4,
        "value": [1.0, 2.0, 3.0, 4.0],
        "units": ["degC"] * 4,
        "source": ["synthetic"] * 4,
    }).to_csv(path, index=False)

    ds = extractor.upload_csv(str(path))

    assert "tasmax" in ds.data_vars
    assert ds.tasmax.attrs["units"] == "degC"
    assert ds.attrs["source"] == "synthetic"


def test_upload_csv_units_default_to_unknown(extractor, tmp_path):
    path = tmp_path / "nounits.csv"
    pd.DataFrame({
        "time": pd.date_range("2014-01-01", periods=2),
        "lat": [52.0, 52.0],
        "lon": [13.0, 13.0],
        "variable": ["tasmax"] * 2,
        "value": [1.0, 2.0],
        "units": [None, None],
    }).to_csv(path, index=False)

    ds = extractor.upload_csv(str(path))
    assert ds.tasmax.attrs["units"] == "unknown"


# --------------------------------------------------------------------------
# extract() dispatch
# --------------------------------------------------------------------------
def test_unhandled_provider_names_the_supported_ones(tmp_path):
    """Regression: an unhandled provider used to leave ds as None and fail with a
    bare 'NoneType is not subscriptable' that named nothing."""
    ex = ClimData(overrides=[
        "dataset=era5", "lat=52.5", "lon=13.4", "variables=[tas]",
        f"output.out_dir={tmp_path}",
    ])
    with pytest.raises(ValueError) as err:
        ex.extract()

    message = str(err.value)
    assert "era5" in message
    assert "MSWX" in message and "NEXGDDP" in message


def test_extractable_list_matches_the_dispatch_branches():
    """The error message must not drift from the code it describes."""
    import inspect as _inspect
    from climdata.utils.wrapper_workflow import ClimateExtractor

    source = _inspect.getsource(ClimateExtractor.extract)
    for name in _EXTRACTABLE:
        assert f'== "{name}"' in source, f"{name} is advertised but not dispatched"


# --------------------------------------------------------------------------
# to_dataframe
# --------------------------------------------------------------------------
def test_to_dataframe_without_a_dataset(extractor):
    with pytest.raises(ValueError, match="No dataset provided"):
        extractor.to_dataframe()


def test_to_dataframe_long_form_shape(extractor, gridded_ds):
    extractor.current_ds = gridded_ds
    df = extractor.to_dataframe(gridded_ds)

    assert {"time", "lat", "lon", "variable", "value", "units", "source"} <= set(df.columns)
    assert set(df["variable"]) == {"tasmax"}
    assert set(df["units"]) == {"degC"}
    assert len(df) == 3 * 2 * 2


def test_to_dataframe_on_a_dataset_with_no_variables(extractor):
    empty = xr.Dataset(coords={"time": pd.date_range("2014-01-01", periods=2)})
    extractor.current_ds = empty
    with pytest.raises(ValueError, match="No variables"):
        extractor.to_dataframe(empty)


def test_to_dataframe_sets_current_df(extractor, gridded_ds):
    extractor.current_ds = gridded_ds
    df = extractor.to_dataframe(gridded_ds)
    assert extractor.current_df is df
    assert extractor.df is df


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------
def test_to_csv_rejects_unknown_format(extractor, gridded_ds):
    extractor.current_ds = gridded_ds
    df = extractor.to_dataframe(gridded_ds)
    with pytest.raises(ValueError, match="Unsupported format"):
        extractor.to_csv(df, format="parquet")


def test_to_csv_default_is_tab_separated(extractor, gridded_ds, tmp_path):
    extractor.current_ds = gridded_ds
    df = extractor.to_dataframe(gridded_ds)
    out = tmp_path / "out.csv"

    extractor.to_csv(df, filename=str(out))

    assert out.exists()
    assert "\t" in out.read_text().splitlines()[0]


def test_to_csv_creates_missing_parent_directories(extractor, gridded_ds, tmp_path):
    extractor.current_ds = gridded_ds
    df = extractor.to_dataframe(gridded_ds)
    out = tmp_path / "deep" / "nested" / "out.csv"

    extractor.to_csv(df, filename=str(out))
    assert out.exists()


def test_to_csv_without_a_filename(extractor, gridded_ds):
    extractor.current_ds = gridded_ds
    df = extractor.to_dataframe(gridded_ds)
    extractor.filename_csv = None
    with pytest.raises(ValueError, match="No filename provided"):
        extractor.to_csv(df)


def test_to_nc_without_a_dataset(extractor):
    extractor.current_ds = None
    with pytest.raises(ValueError, match="No dataset available"):
        extractor.to_nc()


def test_to_nc_round_trips(extractor, gridded_ds, tmp_path):
    out = tmp_path / "out.nc"
    written = extractor.to_nc(gridded_ds, filename=str(out))

    assert written == str(out)
    reopened = xr.open_dataset(out)
    assert "tasmax" in reopened.data_vars
    reopened.close()


@pytest.mark.parametrize("fmt", ["simplace", "monica"])
def test_gridded_format_splits_by_location(extractor, gridded_ds, tmp_path, fmt):
    """One file per grid cell, foldered by column index."""
    extractor.current_ds = gridded_ds
    extractor.ds = gridded_ds
    df = extractor.to_dataframe(gridded_ds)

    base = extractor.to_csv(df, filename=str(tmp_path / fmt), format=fmt)

    written = sorted(p for p in (tmp_path / fmt).rglob("*.csv"))
    assert len(written) == 4                      # 2 lat x 2 lon
    header = written[0].read_text().splitlines()[0]
    assert "Date" in header
    assert "TempMax" in header                    # tasmax renamed to the DWD name
    assert base == str(tmp_path / fmt)


def test_gridded_format_requires_a_time_column(extractor):
    df = pd.DataFrame({"lat": [1.0], "lon": [2.0], "variable": ["tasmax"], "value": [1.0]})
    with pytest.raises(ValueError, match="'date' or 'time' column"):
        extractor.to_csv(df, format="simplace")


def test_gridded_format_requires_coordinate_columns(extractor):
    df = pd.DataFrame({"time": ["2014-01-01"], "variable": ["tasmax"], "value": [1.0]})
    with pytest.raises(ValueError, match="lat.*lon"):
        extractor.to_csv(df, format="simplace")


# --------------------------------------------------------------------------
# calc_index
# --------------------------------------------------------------------------
def test_calc_index_returns_none_when_unconfigured(extractor, gridded_ds):
    """index=None is the default, and must be a no-op rather than an error."""
    extractor.current_ds = gridded_ds
    assert extractor.cfg.index is None
    assert extractor.calc_index() is None


def test_calc_index_without_a_dataset(extractor):
    extractor.current_ds = None
    with pytest.raises(ValueError, match="No dataset provided"):
        extractor.calc_index()


def test_calc_index_warns_on_a_short_record(tmp_path, gridded_ds):
    """Percentile indices need ~30 years; three days must warn, not fail silently."""
    ex = ClimData(overrides=BASE + ["index=tn10p", f"output.out_dir={tmp_path}"])
    ex.current_ds = gridded_ds
    with pytest.warns(UserWarning, match="30 years"):
        with pytest.raises(Exception):
            ex.calc_index()


# --------------------------------------------------------------------------
# Config introspection
# --------------------------------------------------------------------------
def test_get_datasets_lists_providers(extractor):
    datasets = extractor.get_datasets()
    assert "MSWX" in datasets and "HYRAS" in datasets


def test_get_variables_for_a_known_dataset(extractor):
    assert "pr" in extractor.get_variables("HYRAS")


def test_get_variables_for_an_unknown_dataset(extractor):
    with pytest.raises(ValueError, match="No variable info"):
        extractor.get_variables("NOT_A_DATASET")


def test_get_varinfo_for_an_unknown_variable(extractor):
    with pytest.raises(ValueError, match="not found in varinfo"):
        extractor.get_varinfo("not_a_variable")


def test_get_varinfo_returns_units(extractor):
    assert "units" in extractor.get_varinfo("tasmax")


def test_get_indices_require_all_is_stricter_than_any(extractor):
    variables = ["tasmin"]
    strict = extractor.get_indices(variables, require_all=True)
    loose = extractor.get_indices(variables, require_all=False)
    assert set(strict) <= set(loose)


def test_get_indices_with_no_matching_variables(extractor):
    assert extractor.get_indices(["not_a_variable"], require_all=True) == {}


def test_get_actions_is_a_mapping(extractor):
    actions = extractor.get_actions()
    assert "extract" in actions


def test_get_impute_methods_is_a_mapping(extractor):
    assert isinstance(extractor.get_impute_methods(), dict)


# --------------------------------------------------------------------------
# reproject
# --------------------------------------------------------------------------
def test_reproject_without_a_target_is_a_noop(extractor, gridded_ds):
    """Neither target configured: skip, leaving the dataset untouched."""
    extractor.current_ds = gridded_ds
    assert extractor.reproject() is None
    assert extractor.current_ds is gridded_ds


def test_reproject_without_a_dataset(extractor):
    extractor.current_ds = None
    with pytest.raises(ValueError, match="No dataset provided"):
        extractor.reproject()


def test_reproject_rejects_metric_resolution_on_a_geographic_crs(tmp_path, gridded_ds):
    """The gate must fire before any data is touched: 10 km has no fixed angular size."""
    pytest.importorskip("rioxarray")
    ex = ClimData(overrides=BASE + [
        "target_projection=EPSG:4326", "target_resolution=10 km",
        f"output.out_dir={tmp_path}",
    ])
    ex.current_ds = gridded_ds
    with pytest.raises(ResolutionCRSMismatch):
        ex.reproject()


# --------------------------------------------------------------------------
# run_workflow
# --------------------------------------------------------------------------
def test_unknown_action_is_rejected(extractor):
    with pytest.raises(ValueError, match="Unknown action"):
        extractor.run_workflow(actions=["not_an_action"])


def test_upload_actions_require_a_file(extractor):
    with pytest.raises(ValueError, match="requires argument"):
        extractor.run_workflow(actions=["upload_netcdf"])
    with pytest.raises(ValueError, match="requires argument"):
        extractor.run_workflow(actions=["upload_csv"])


@pytest.mark.parametrize("action,bad", [
    ("upload_netcdf", "data.csv"),
    ("upload_csv", "data.nc"),
])
def test_upload_actions_validate_the_extension(extractor, action, bad):
    """Checked before opening, so a swapped argument fails clearly."""
    with pytest.raises(ValueError, match="Invalid file format"):
        extractor.run_workflow(actions=[action], file=bad)


def test_actions_requiring_a_dataset_fail_clearly(extractor):
    extractor.current_ds = None
    for action in ("calc_index", "to_csv", "to_nc", "to_fair", "impute"):
        with pytest.raises(ValueError, match="requires a dataset"):
            extractor.run_workflow(actions=[action])


def test_workflow_from_upload_to_files(extractor, gridded_ds, tmp_path):
    """The offline half of the pipeline, end to end."""
    src = tmp_path / "in.nc"
    gridded_ds.to_netcdf(src)

    result = extractor.run_workflow(
        actions=["upload_netcdf", "to_csv", "to_nc"], file=str(src)
    )

    assert isinstance(result, WorkflowResult)
    assert result.dataset is not None
    assert result.dataframe is not None
    assert "dataset" in result.keys() and "dataframe" in result.keys()


def test_workflow_result_keys_reports_only_populated_fields():
    result = WorkflowResult(cfg=None)
    assert result.keys() == []
    result.filename = "x.nc"
    assert result.keys() == ["filename"]


# --------------------------------------------------------------------------
# to_fair
# --------------------------------------------------------------------------
def test_to_fair_without_a_dataset(extractor):
    extractor.current_ds = None
    with pytest.raises(ValueError, match="No dataset available"):
        extractor.to_fair()


def test_to_fair_writes_a_complete_crate(extractor, gridded_ds, tmp_path):
    crate = tmp_path / "crate"
    path = extractor.to_fair(ds=gridded_ds, output_dir=str(crate), title="Test crate")

    assert (crate / "data.nc").exists()
    assert (crate / "ro-crate-metadata.json").exists()
    assert (crate / "ro-crate-preview.html").exists()
    assert path == str(crate)

    meta = json.loads((crate / "ro-crate-metadata.json").read_text())
    assert "@graph" in meta or "@context" in meta


def test_to_fair_can_zip(extractor, gridded_ds, tmp_path):
    path = extractor.to_fair(
        ds=gridded_ds, output_dir=str(tmp_path / "crate"), zip_crate=True
    )
    assert path.endswith(".zip")


# --------------------------------------------------------------------------
# Dask opt-in
# --------------------------------------------------------------------------
def test_dask_is_opt_in(extractor):
    """Disabled by default, so importing climdata never starts a cluster."""
    assert extractor._ensure_dask_client() is None


def test_close_dask_is_safe_without_a_cluster(extractor):
    extractor.close_dask()          # must not raise


# --------------------------------------------------------------------------
# Config / code agreement
# --------------------------------------------------------------------------
def test_every_configured_action_is_dispatchable(extractor):
    """get_actions() reads actions.yaml while run_workflow() dispatches in code.
    They drifted once — actions.yaml advertised 'reproject' and run_workflow
    raised 'Unknown action' for it."""
    configured = set(extractor.get_actions())
    assert configured <= set(_WORKFLOW_ACTIONS), (
        f"advertised but not dispatched: {configured - set(_WORKFLOW_ACTIONS)}"
    )


def test_every_dispatchable_action_is_reachable(extractor):
    """Each name in _WORKFLOW_ACTIONS must reach a real branch, not the else."""
    for action in _WORKFLOW_ACTIONS:
        extractor.current_ds = None
        try:
            extractor.run_workflow(actions=[action])
        except ValueError as err:
            assert "Unknown action" not in str(err), action
        except Exception:
            pass          # any other failure means the branch was reached


def test_the_unknown_action_message_lists_the_valid_ones(extractor):
    with pytest.raises(ValueError, match="Supported:") as err:
        extractor.run_workflow(actions=["nope"])
    assert "calc_index" in str(err.value)


def test_reproject_action_is_dispatched(extractor, gridded_ds):
    """No target grid configured, so it is a no-op — but it must not raise."""
    extractor.current_ds = gridded_ds
    result = extractor.run_workflow(actions=["reproject"])
    assert result.dataset is not None
