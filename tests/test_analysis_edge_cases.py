"""Edge cases for the analysis modules: extremes, imputation, sdba and viz.

Grouped together because they share a shape — each takes an xarray object and
returns a derived one, so the interesting cases are degenerate inputs: a record
too short, a series with no gaps, a field with no spatial dimensions.

The recurring theme is that these functions must degrade honestly. A statistic
that cannot be computed should be absent or ``None``, never a plausible-looking
zero.
"""

import numpy as np
import pandas as pd
import pytest
import xarray as xr
from omegaconf import OmegaConf

from climdata.extremes.indices import extreme_index
from climdata.sdba import compute_daily_climatology, smooth_fft


# ==========================================================================
# Fixtures
# ==========================================================================
@pytest.fixture
def daily_grid():
    """Two years of daily data on a 2x2 grid, with a seasonal cycle."""
    time = pd.date_range("2000-01-01", "2001-12-31", freq="D")
    doy = time.dayofyear.values
    base = 10 + 10 * np.sin(2 * np.pi * doy / 365.25)
    rng = np.random.default_rng(0)
    values = base[:, None, None] + rng.normal(0, 1, (len(time), 2, 2))
    ds = xr.Dataset(
        {"tasmin": (("time", "lat", "lon"), values)},
        coords={"time": time, "lat": [52.0, 52.5], "lon": [13.0, 13.5]},
    )
    ds.tasmin.attrs["units"] = "degC"
    return ds


# ==========================================================================
# extreme_index
# ==========================================================================
def _index_cfg(indices):
    return OmegaConf.create({"extinfo": {"indices": indices}})


def test_unknown_index_raises(daily_grid):
    idx = extreme_index(_index_cfg({}), daily_grid)
    with pytest.raises(KeyError):
        idx.calculate("no_such_index")


def test_link_without_a_recipe_is_rejected(daily_grid):
    """A link must say how to build the intermediate, or it is a config error."""
    cfg = _index_cfg({
        "bad": {
            "function": "numpy.mean",
            "args": {},
            "link": {"intermediate": {"inputs": ["tasmin"]}},   # no function_call/operation
        }
    })
    with pytest.raises(ValueError, match="function_call.*operation"):
        extreme_index(cfg, daily_grid).calculate("bad")


def test_unsupported_postprocess_is_rejected(daily_grid):
    cfg = _index_cfg({
        "bad": {
            "function": "numpy.mean",
            "args": {},
            "link": {
                "intermediate": {
                    "input": "tasmin",
                    "operation": "mean",
                    "postprocess": {"resample": {"time": "YS"}},
                }
            },
        }
    })
    with pytest.raises(NotImplementedError, match="not supported"):
        extreme_index(cfg, daily_grid).calculate("bad")


def test_a_link_referencing_a_missing_variable_raises(daily_grid):
    cfg = _index_cfg({
        "bad": {
            "function": "numpy.mean",
            "args": {},
            "link": {"intermediate": {"input": "not_a_variable", "operation": "mean"}},
        }
    })
    with pytest.raises(KeyError):
        extreme_index(cfg, daily_grid).calculate("bad")


def test_a_link_result_is_stored_for_later_reference(daily_grid):
    """Linked intermediates go back into climate_data so args can reference them."""
    cfg = _index_cfg({
        "demo": {
            "function": "xarray.ufuncs.fabs" if hasattr(xr, "ufuncs") else "numpy.fabs",
            "variables": ["intermediate"],
            "args": {},
            "link": {"intermediate": {"input": "tasmin", "operation": "mean"}},
        }
    })
    idx = extreme_index(cfg, daily_grid)
    idx.calculate("demo")
    assert "intermediate" in idx.climate_data


def test_the_result_is_renamed_to_the_index(daily_grid):
    cfg = _index_cfg({
        "my_index": {"function": "numpy.fabs", "variables": ["tasmin"], "args": {}}
    })
    result = extreme_index(cfg, daily_grid).calculate("my_index")
    assert isinstance(result, xr.DataArray)
    assert result.name == "my_index"


def test_a_real_xclim_index_computes(daily_grid):
    """End to end through the config-driven path, with a genuine indicator."""
    pytest.importorskip("xclim")
    cfg = _index_cfg({
        "frost_days": {
            "function": "xclim.indicators.atmos.frost_days",
            "variables": ["tasmin"],
            "args": {"freq": "YS"},
        }
    })
    result = extreme_index(cfg, daily_grid).calculate("frost_days").compute()
    assert result.name == "frost_days"
    assert result.sizes["time"] == 2          # two years in the fixture


# ==========================================================================
# Imputation
# ==========================================================================
imputegap = pytest.importorskip(
    "climdata._vendor.imputegap.recovery.manager", reason="imputegap unavailable"
)


@pytest.fixture
def gappy():
    """A small gridded dataset with a handful of holes."""
    time = pd.date_range("2000-01-01", periods=120, freq="D")
    rng = np.random.default_rng(1)
    values = rng.normal(10, 2, (120, 2, 2))
    values[10:15, 0, 0] = np.nan
    return xr.Dataset(
        {"tas": (("time", "lat", "lon"), values)},
        coords={"time": time, "lat": [52.0, 52.5], "lon": [13.0, 13.5]},
    )


@pytest.fixture
def complete(gappy):
    return gappy.fillna(10.0)


def test_missing_fraction_counts_gaps(gappy):
    from climdata.impute.impute_xarray import Imputer

    frac = Imputer(gappy, method="SoftImpute").missing_fraction()
    assert frac["tas"] == pytest.approx(5 / (120 * 4))
    assert frac["global"] == pytest.approx(frac["tas"])


def test_missing_fraction_is_zero_for_complete_data(complete):
    from climdata.impute.impute_xarray import Imputer

    assert Imputer(complete, method="SoftImpute").missing_fraction()["global"] == 0.0


def test_impute_is_a_noop_when_nothing_is_missing(complete):
    """Calling impute() unconditionally in a workflow must be safe and cheap."""
    from climdata.impute.impute_xarray import Imputer

    imputer = Imputer(complete, method="SoftImpute")
    out = imputer.impute()
    xr.testing.assert_allclose(out.tas, complete.tas)


def test_metrics_before_impute_raises(gappy):
    from climdata.impute.impute_xarray import Imputer

    with pytest.raises(RuntimeError, match="impute\\(\\) first"):
        Imputer(gappy, method="SoftImpute").metrics()


def test_unknown_method_is_rejected(gappy):
    from climdata.impute.impute_xarray import Imputer

    with pytest.raises(ValueError, match="Unknown method"):
        Imputer(gappy, method="NotAMethod").impute()


def test_summary_describes_the_setup_without_running(gappy):
    from climdata.impute.impute_xarray import Imputer

    summary = Imputer(gappy, method="SoftImpute", normalize=False).summary()
    assert summary["method"] == "SoftImpute"
    assert summary["normalize"] is False
    assert summary["variables"] == ["tas"]
    assert summary["missing_fraction"]["global"] > 0


def test_a_torch_method_reports_the_missing_backend(gappy):
    """The check runs at construction, before any reshaping work."""
    from climdata.impute.impute_xarray import Imputer

    try:
        import torch  # noqa: F401
    except ImportError:
        with pytest.raises(ImportError, match="requires PyTorch"):
            Imputer(gappy, method="BRITS")
    else:
        Imputer(gappy, method="BRITS")          # available: must construct fine


# ==========================================================================
# sdba smoothing
# ==========================================================================
def test_smooth_fft_preserves_length():
    x = np.random.default_rng(0).normal(size=365)
    assert smooth_fft(x).shape == x.shape


def test_smooth_fft_reduces_variance():
    """Truncating the spectrum can only remove variance, never add it."""
    rng = np.random.default_rng(0)
    doy = np.arange(365)
    x = 10 * np.sin(2 * np.pi * doy / 365) + rng.normal(0, 3, 365)
    assert smooth_fft(x, n_harmonics=2).std() < x.std()


def test_smooth_fft_preserves_the_mean():
    """Harmonic 0 is the mean and is always retained."""
    x = np.random.default_rng(0).normal(5, 2, 365)
    assert smooth_fft(x).mean() == pytest.approx(x.mean())


def test_more_harmonics_follow_the_signal_more_closely():
    rng = np.random.default_rng(0)
    doy = np.arange(365)
    x = 10 * np.sin(2 * np.pi * doy / 365) + 3 * np.sin(6 * np.pi * doy / 365)
    err = lambda n: np.abs(smooth_fft(x, n_harmonics=n) - x).mean()
    assert err(5) <= err(1)


def test_smooth_fft_is_exact_for_a_pure_harmonic():
    doy = np.arange(365)
    x = np.sin(2 * np.pi * doy / 365)
    assert np.allclose(smooth_fft(x, n_harmonics=3), x, atol=1e-8)


def test_compute_daily_climatology_replaces_time_with_dayofyear(daily_grid):
    clim = compute_daily_climatology(daily_grid.tasmin)
    assert "dayofyear" in clim.dims
    assert "time" not in clim.dims
    assert set(clim.dims) == {"dayofyear", "lat", "lon"}


def test_compute_daily_climatology_is_smoother_than_the_raw_cycle(daily_grid):
    raw = daily_grid.tasmin.groupby("time.dayofyear").mean("time")
    clim = compute_daily_climatology(daily_grid.tasmin)
    # day-to-day jitter, measured as the mean absolute first difference
    jitter = lambda da: float(np.abs(np.diff(da.isel(lat=0, lon=0).values)).mean())
    assert jitter(clim) < jitter(raw)


# ==========================================================================
# viz
# ==========================================================================
matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")


@pytest.fixture
def field():
    lat = np.linspace(47, 55, 9)
    lon = np.linspace(5, 15, 11)
    rng = np.random.default_rng(0)
    ds = xr.Dataset(
        {
            "tas": (("time", "lat", "lon"), rng.normal(10, 3, (4, 9, 11))),
            "pr": (("time", "lat", "lon"), rng.random((4, 9, 11))),
        },
        coords={"time": pd.date_range("2020-01-01", periods=4, freq="D"),
                "lat": lat, "lon": lon},
    )
    ds.tas.attrs["units"] = "degC"
    return ds


def test_multi_variable_dataset_requires_an_explicit_variable(field):
    from climdata.viz import _as_dataarray

    with pytest.raises(ValueError, match="pass `variable=`"):
        _as_dataarray(field)


def test_single_variable_dataset_needs_no_variable(field):
    from climdata.viz import _as_dataarray

    assert _as_dataarray(field[["tas"]]).name == "tas"


def test_a_dataarray_passes_through(field):
    from climdata.viz import _as_dataarray

    da = field.tas
    assert _as_dataarray(da) is da


def test_a_non_xarray_input_is_rejected():
    from climdata.viz import _as_dataarray

    with pytest.raises(TypeError, match="Expected an xarray"):
        _as_dataarray([1, 2, 3])


def test_missing_spatial_coordinates_are_reported():
    from climdata.viz import _spatial_dim_names

    da = xr.DataArray([1.0, 2.0], dims="time", coords={"time": [0, 1]})
    with pytest.raises(ValueError, match="longitude/latitude"):
        _spatial_dim_names(da)


@pytest.mark.parametrize("x_name,y_name", [
    ("lon", "lat"), ("longitude", "latitude"), ("x", "y"),
])
def test_spatial_coordinate_aliases_are_recognised(x_name, y_name):
    from climdata.viz import _spatial_dim_names

    da = xr.DataArray(
        np.zeros((2, 2)), dims=(y_name, x_name),
        coords={y_name: [0.0, 1.0], x_name: [0.0, 1.0]},
    )
    assert _spatial_dim_names(da) == (x_name, y_name)


def test_plot_map_draws_a_field(field):
    cartopy = pytest.importorskip("cartopy")          # noqa: F841
    from climdata.viz import plot_map

    ax = plot_map(field, variable="tas")
    assert ax.get_title().startswith("tas")


def test_plot_map_reduce_labels_the_title(field):
    pytest.importorskip("cartopy")
    from climdata.viz import plot_map

    ax = plot_map(field, variable="tas", reduce="mean")
    assert "mean over time" in ax.get_title()


def test_plot_map_rejects_point_data(field):
    """A point extraction has no grid to draw; say so rather than failing obscurely."""
    pytest.importorskip("cartopy")
    from climdata.viz import plot_map

    point = field.isel(lat=0, lon=0)
    with pytest.raises(ValueError, match="point/station extractions"):
        plot_map(point, variable="tas")
