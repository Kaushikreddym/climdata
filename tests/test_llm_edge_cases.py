"""Edge cases for the LLM layer: schema validation, analytics and grounding.

No network and no real model: the ``fake_llm`` fixture in ``conftest.py`` stands
in for the client, so the planner's repair loop and the narrator's grounding
check can be driven deterministically.

The load-bearing property of this layer is that it never invents anything. The
planner validates the model's output against the live capability registry, and
the narrator checks that every number in the prose traces back to a computed
statistic. Both are tested here against deliberately wrong model output.
"""

import numpy as np
import pytest
import xarray as xr

from climdata.LLM import analytics
from climdata.LLM.narrator import grounding_check
from climdata.LLM.planner import ReadyPlan, validate_response
from climdata.LLM.tool_registry import REGISTRY


# ==========================================================================
# ReadyPlan: spatial validation
# ==========================================================================
BASE = {"variable": "tasmax", "analysis": ["historical_trend"]}


def test_a_point_plan_validates():
    plan = ReadyPlan(lat=52.5, lon=13.4, **BASE)
    assert plan.lat == 52.5


def test_an_extent_plan_validates():
    plan = ReadyPlan(lat_min=47.0, lat_max=55.0, lon_min=5.0, lon_max=15.0, **BASE)
    assert plan.lat_min == 47.0


def test_no_coordinates_at_all_is_rejected():
    with pytest.raises(ValueError, match="needs complete coordinates"):
        ReadyPlan(**BASE)


def test_a_point_and_an_extent_together_are_rejected():
    """Ambiguous: the downstream extraction can only honour one."""
    with pytest.raises(ValueError, match="point OR an extent"):
        ReadyPlan(lat=52.5, lon=13.4, lat_min=47.0, lat_max=55.0,
                  lon_min=5.0, lon_max=15.0, **BASE)


def test_half_a_point_is_rejected():
    """lat without lon is neither a complete point nor an extent."""
    with pytest.raises(ValueError, match="needs complete coordinates"):
        ReadyPlan(lat=52.5, **BASE)


def test_half_an_extent_is_rejected():
    with pytest.raises(ValueError, match="needs complete coordinates"):
        ReadyPlan(lat_min=47.0, lat_max=55.0, **BASE)


@pytest.mark.parametrize("lat", [-91.0, 91.0, 1000.0])
def test_out_of_range_latitude_is_rejected(lat):
    with pytest.raises(ValueError, match="out of range"):
        ReadyPlan(lat=lat, lon=13.4, **BASE)


@pytest.mark.parametrize("lon", [-181.0, 181.0])
def test_out_of_range_longitude_is_rejected(lon):
    with pytest.raises(ValueError, match="out of range"):
        ReadyPlan(lat=52.5, lon=lon, **BASE)


@pytest.mark.parametrize("lat", [-90.0, 90.0])
def test_the_poles_are_accepted(lat):
    """Range checks are inclusive; +/-90 is a valid latitude."""
    assert ReadyPlan(lat=lat, lon=0.0, **BASE).lat == lat


def test_an_inverted_extent_is_rejected():
    with pytest.raises(ValueError, match="lat_min must be < lat_max"):
        ReadyPlan(lat_min=55.0, lat_max=47.0, lon_min=5.0, lon_max=15.0, **BASE)


def test_a_degenerate_extent_is_rejected():
    """Equal bounds enclose no cells; that is a mistake, not a point request."""
    with pytest.raises(ValueError, match="lon_min must be < lon_max"):
        ReadyPlan(lat_min=47.0, lat_max=55.0, lon_min=10.0, lon_max=10.0, **BASE)


# ==========================================================================
# ReadyPlan: registry grounding
# ==========================================================================
def test_an_invented_variable_is_rejected():
    """The whole point of validation: the model may not make names up."""
    with pytest.raises(ValueError):
        ReadyPlan(lat=52.5, lon=13.4, variable="unobtainium",
                  analysis=["historical_trend"])


def test_an_invented_dataset_is_rejected():
    with pytest.raises(ValueError):
        ReadyPlan(lat=52.5, lon=13.4, dataset="NOT_A_DATASET", **BASE)


def test_an_invented_index_is_rejected():
    with pytest.raises(ValueError):
        ReadyPlan(lat=52.5, lon=13.4, variable="tasmax",
                  analysis=["extremes"], indices=["not_an_index"])


def test_several_variables_are_accepted():
    plan = ReadyPlan(lat=52.5, lon=13.4, variable=["tasmax", "tasmin"],
                     analysis=["historical_trend"])
    assert plan.variables == ["tasmax", "tasmin"]


def test_variables_normalises_a_single_name_to_a_list():
    assert ReadyPlan(lat=52.5, lon=13.4, **BASE).variables == ["tasmax"]


# ==========================================================================
# validate_response
# ==========================================================================
def test_a_valid_ready_payload_parses(ready_point_payload):
    parsed = validate_response(ready_point_payload)
    assert isinstance(parsed, ReadyPlan)
    assert parsed.lat == 52.52


def test_a_payload_missing_coordinates_does_not_parse_as_ready():
    bad = {"status": "ready", "variable": "tasmax", "analysis": ["historical_trend"]}
    with pytest.raises(Exception):
        validate_response(bad)


# ==========================================================================
# Analytics: honest degradation
# ==========================================================================
def _series(years=35, warming=0.05, seed=0):
    time = xr.cftime_range("1980-01-01", periods=years * 365, freq="D",
                           calendar="noleap")
    t = np.arange(len(time))
    rng = np.random.default_rng(seed)
    values = 15 + warming * (t / 365.0) + rng.normal(0, 0.5, len(time))
    da = xr.DataArray(values, coords={"time": time}, dims="time", name="tasmax")
    da.attrs["units"] = "degC"
    return da


def test_annual_trend_recovers_a_known_slope():
    """A synthetic 0.05 degC/yr warming must come back as roughly 0.05."""
    assert analytics.annual_trend(_series(warming=0.05)) == pytest.approx(0.05, abs=0.01)


def test_annual_trend_on_a_flat_series_is_near_zero():
    assert analytics.annual_trend(_series(warming=0.0)) == pytest.approx(0.0, abs=0.01)


def test_seasonal_trend_with_no_data_for_the_season_returns_nan():
    """NaN, not an exception on an empty grouping."""
    jan = _series(years=2).sel(time=slice("1980-01-01", "1980-01-31"))
    assert np.isnan(analytics.seasonal_trend(jan, season="JJA"))


def test_climatology_covers_every_month():
    clim = analytics.climatology(_series(years=3))
    assert set(clim["monthly_mean"]) == set(range(1, 13))
    assert isinstance(clim["overall_mean"], float)


def test_hot_days_change_is_none_for_a_single_year():
    """One year cannot support a change; None beats a fabricated 0."""
    assert analytics.hot_days_change(_series(years=1)) is None


def test_hot_days_change_ignores_years_with_no_observations():
    """An all-NaN year must not read as zero exceedances, i.e. as cooling."""
    da = _series(years=25, warming=0.2)
    full = analytics.hot_days_change(da, thresh=15.0)

    gapped = da.copy()
    gapped.loc[{"time": slice("1990-01-01", "1990-12-31")}] = np.nan
    with_gap = analytics.hot_days_change(gapped, thresh=15.0)

    assert full is not None and with_gap is not None
    assert abs(full - with_gap) <= max(3, abs(full) * 0.5)


@pytest.mark.parametrize("variable,expected", [
    ("tas", True), ("tasmax", True), ("tasmin", True), ("soilTemp", True),
    ("pr", False), ("sfcWind", False), ("hurs", False),
])
def test_is_temperature_guards_the_degree_threshold(variable, expected):
    assert analytics.is_temperature(None, variable) is expected


def test_is_temperature_falls_back_to_units():
    da = xr.DataArray([1.0], dims="time", attrs={"units": "degC"})
    assert analytics.is_temperature(da) is True
    da.attrs["units"] = "mm/day"
    assert analytics.is_temperature(da) is False


def test_is_temperature_is_false_when_nothing_can_be_established():
    """The safe answer: suppress the statistic rather than invent one."""
    assert analytics.is_temperature(xr.DataArray([1.0], dims="time")) is False


@pytest.mark.parametrize("variable,units", [
    ("tas", "degC"), ("tasmax", "degC"), ("pr", "mm/day"),
])
def test_canonical_units_are_fixed_regardless_of_source(variable, units):
    assert analytics.variable_units(variable) == units


def test_units_for_an_unknown_variable_are_none():
    assert analytics.variable_units("not_a_variable") is None


def test_evaluation_metrics_on_identical_series_are_perfect():
    da = _series(years=5)
    m = analytics.evaluation_metrics(da, da)
    assert m["model_bias"] == pytest.approx(0.0, abs=1e-9)
    assert m["rmse"] == pytest.approx(0.0, abs=1e-9)
    assert m["correlation"] == pytest.approx(1.0, abs=1e-9)


def test_evaluation_metrics_with_no_overlap_are_all_nan():
    """No shared period: NaN throughout rather than an exception."""
    early = _series(years=2)
    late = _series(years=2)
    late = late.assign_coords(
        time=xr.cftime_range("2050-01-01", periods=late.sizes["time"],
                             freq="D", calendar="noleap")
    )
    m = analytics.evaluation_metrics(late, early)
    assert all(np.isnan(v) for v in m.values())


def test_evaluation_metrics_report_a_known_bias():
    obs = _series(years=5)
    model = obs + 2.0
    assert analytics.evaluation_metrics(model, obs)["model_bias"] == pytest.approx(2.0, abs=1e-6)


def test_capability_report_covers_every_analysis():
    report = analytics.capability_report()
    for key in ("etccdi_indices", "extremes", "climatology", "trends",
                "evaluation_metrics"):
        assert report[key]["mode"] in {"reuse", "implement"}


def test_analyze_only_returns_requested_intents():
    out = analytics.analyze(_series(years=5), intents=["climatology"],
                            variable="tasmax")
    assert "overall_mean" in out
    assert "annual_trend" not in out


def test_analyze_with_no_intents_returns_almost_nothing():
    out = analytics.analyze(_series(years=5), intents=[], variable="tasmax")
    assert "annual_trend" not in out and "overall_mean" not in out


# ==========================================================================
# Narrator grounding
# ==========================================================================
def test_prose_with_no_numbers_is_grounded():
    """grounding_check returns the ungrounded tokens; empty means clean."""
    assert grounding_check("Temperatures have risen over the record.", {}) == []


def test_a_quoted_statistic_is_grounded():
    assert grounding_check("The trend is 0.05 degC per year.",
                           {"annual_trend": 0.05}) == []


def test_a_faithful_rounding_is_grounded():
    """Quoting 0.05 for 0.0483 is reporting, not inventing."""
    assert grounding_check("The trend is 0.05 degC per year.",
                           {"annual_trend": 0.0483}) == []


def test_an_invented_number_is_caught():
    """The whole point: a figure with no matching statistic must be flagged."""
    bad = grounding_check("The trend is 0.42 degC per year.",
                          {"annual_trend": 0.05})
    assert "0.42" in bad


def test_a_derived_number_is_caught():
    """Per-decade rescaling is arithmetic the narrator must not do silently."""
    bad = grounding_check("Warming of 0.5 degC per decade.",
                          {"annual_trend": 0.05})
    assert bad


def test_nested_statistics_are_collected():
    """Monthly means live one level down and must still count as grounded."""
    assert grounding_check(
        "January averaged 2.5 degC.", {"monthly_mean": {1: 2.5}}
    ) == []


# ==========================================================================
# Tool registry
# ==========================================================================
def test_the_builtin_diagnostics_are_registered():
    names = {t.name for t in REGISTRY.list("diagnostic")}
    assert {"historical_trend", "climatology", "extremes",
            "model_evaluation", "future_projections"} <= names


def test_an_unknown_tool_lookup_raises():
    with pytest.raises(KeyError, match="No diagnostic tool named"):
        REGISTRY.get("diagnostic", "not_a_tool")


def test_has_does_not_raise():
    assert REGISTRY.has("diagnostic", "climatology")
    assert not REGISTRY.has("diagnostic", "not_a_tool")


def test_diagnostics_are_indexed_by_intent():
    found = REGISTRY.diagnostics_for_intent("historical_trend")
    assert [t.name for t in found] == ["historical_trend"]


def test_an_unserved_intent_yields_nothing():
    assert REGISTRY.diagnostics_for_intent("not_an_intent") == []


def test_every_diagnostic_declares_its_required_roles():
    """The orchestrator skips on missing roles; an undeclared one would crash it."""
    for tool in REGISTRY.list("diagnostic"):
        assert tool.metadata.get("required_roles")
        assert tool.metadata.get("intent")


def test_registering_a_diagnostic_makes_it_discoverable():
    @REGISTRY.diagnostic("test_only_diagnostic", intent="test_only_intent",
                         required_roles=["observation"])
    def _tool(roles, variable=None):
        return {"value": 1}

    assert REGISTRY.has("diagnostic", "test_only_diagnostic")
    assert REGISTRY.get("diagnostic", "test_only_diagnostic").fn(roles={}) == {"value": 1}


def test_the_capability_bridge_reports_datasets_and_indices():
    assert len(REGISTRY.datasets()) > 0
    assert len(REGISTRY.indices()) > 0


# ==========================================================================
# The fake client (fixtures moved here from the package)
# ==========================================================================
def test_fake_llm_replays_scripted_replies(fake_llm):
    client = fake_llm(plan_replies=[{"status": "ready"}, {"status": "second"}])
    system = [{"role": "system", "content": "PLANNER"}]
    first = client.chat.completions.create(model="m", messages=system)
    second = client.chat.completions.create(model="m", messages=system)
    assert "ready" in first.choices[0].message.content
    assert "second" in second.choices[0].message.content


def test_fake_llm_repeats_its_last_reply_once_exhausted(fake_llm):
    client = fake_llm(plan_replies=[{"status": "only"}])
    system = [{"role": "system", "content": "PLANNER"}]
    for _ in range(3):
        reply = client.chat.completions.create(model="m", messages=system)
    assert "only" in reply.choices[0].message.content


def test_fake_llm_records_which_model_each_call_used(fake_llm):
    client = fake_llm(plan_replies=[{"status": "ready"}])
    system = [{"role": "system", "content": "PLANNER"}]
    client.chat.completions.create(model="model-a", messages=system)
    client.chat.completions.create(model="model-b", messages=system)
    assert client.models_seen == ["model-a", "model-b"]


def test_fake_llm_can_simulate_a_transport_failure(fake_llm):
    client = fake_llm(plan_replies=[{}], raise_on={"plan"})
    with pytest.raises(RuntimeError, match="simulated plan"):
        client.chat.completions.create(
            model="m", messages=[{"role": "system", "content": "PLANNER"}]
        )


def test_synthetic_observation_fixture_has_a_warming_signal(obs):
    """The offline provider must produce data with a known, recoverable trend."""
    assert analytics.annual_trend(obs) > 0.01


def test_synthetic_grid_fixture_is_three_dimensional(obs_grid):
    assert set(obs_grid.dims) == {"time", "lat", "lon"}
