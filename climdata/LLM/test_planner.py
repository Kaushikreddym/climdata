"""
Unit tests for the ClimData planner (offline — no LLM, no network).

Run:  pytest test_planner.py -q     (or: python -m pytest test_planner.py)
"""

import pytest
from pydantic import ValidationError

from planner import (
    validate_response, ReadyPlan, ClarificationRequest,
    Status, ClarificationReason, _MISSING_COORDS_MSG,
)


# --------------------------------------------------------------------------
# READY — point
# --------------------------------------------------------------------------

def test_ready_point_ok():
    r = validate_response({
        "status": "ready", "lat": 52.5, "lon": 13.4,
        "variable": "tasmax", "analysis": ["historical_trend"],
    })
    assert isinstance(r, ReadyPlan)
    assert r.status == Status.READY.value
    assert r.lat == 52.5 and r.lon == 13.4
    assert r.variables == ["tasmax"]


def test_ready_extent_ok():
    r = validate_response({
        "status": "ready", "lat_min": 52.0, "lat_max": 53.0,
        "lon_min": 13.0, "lon_max": 14.0,
        "variable": "tasmax", "analysis": ["historical_trend"],
    })
    assert isinstance(r, ReadyPlan)
    assert r.lat_min == 52.0 and r.lon_max == 14.0


def test_multi_variable_list_ok():
    r = validate_response({
        "status": "ready", "lat": 48.1, "lon": 11.6,
        "variable": ["tasmin", "tasmax"], "analysis": ["extremes"],
        "indices": ["heat_wave_frequency"],
    })
    assert r.variables == ["tasmin", "tasmax"]


# --------------------------------------------------------------------------
# NEEDS_CLARIFICATION
# --------------------------------------------------------------------------

def test_clarification_parses_and_routes():
    r = validate_response({
        "status": "needs_clarification",
        "reason": "missing_coordinates",
        "message": _MISSING_COORDS_MSG,
    })
    assert isinstance(r, ClarificationRequest)
    assert r.reason == ClarificationReason.MISSING_COORDINATES.value


def test_place_name_without_coords_cannot_be_ready():
    # A ReadyPlan with a place name but no coordinates is invalid — the planner
    # must instead emit needs_clarification.
    with pytest.raises(ValidationError):
        validate_response({
            "status": "ready", "variable": "tasmax",
            "analysis": ["historical_trend"],
        })


# --------------------------------------------------------------------------
# SPATIAL validation
# --------------------------------------------------------------------------

def test_partial_point_rejected():
    with pytest.raises(ValidationError):
        validate_response({
            "status": "ready", "lat": 52.5,  # lon missing
            "variable": "tasmax", "analysis": ["historical_trend"],
        })


def test_partial_extent_rejected():
    with pytest.raises(ValidationError):
        validate_response({
            "status": "ready", "lat_min": 52.0, "lat_max": 53.0, "lon_min": 13.0,
            "variable": "tasmax", "analysis": ["historical_trend"],  # lon_max missing
        })


def test_point_and_extent_together_rejected():
    with pytest.raises(ValidationError):
        validate_response({
            "status": "ready", "lat": 52.5, "lon": 13.4,
            "lat_min": 52.0, "lat_max": 53.0, "lon_min": 13.0, "lon_max": 14.0,
            "variable": "tasmax", "analysis": ["historical_trend"],
        })


def test_lat_out_of_range_rejected():
    with pytest.raises(ValidationError):
        validate_response({
            "status": "ready", "lat": 99.0, "lon": 13.4,
            "variable": "tasmax", "analysis": ["historical_trend"],
        })


def test_extent_ordering_rejected():
    with pytest.raises(ValidationError):
        validate_response({
            "status": "ready", "lat_min": 53.0, "lat_max": 52.0,  # min >= max
            "lon_min": 13.0, "lon_max": 14.0,
            "variable": "tasmax", "analysis": ["historical_trend"],
        })


# --------------------------------------------------------------------------
# REGISTRY grounding — never invent
# --------------------------------------------------------------------------

def test_unsupported_variable_rejected():
    with pytest.raises(ValidationError):
        validate_response({
            "status": "ready", "lat": 52.5, "lon": 13.4,
            "variable": "snowfall", "analysis": ["historical_trend"],
        })


def test_unsupported_dataset_rejected():
    with pytest.raises(ValidationError):
        validate_response({
            "status": "ready", "lat": 52.5, "lon": 13.4, "dataset": "FAKESET",
            "variable": "tasmax", "analysis": ["historical_trend"],
        })


def test_unsupported_index_rejected():
    with pytest.raises(ValidationError):
        validate_response({
            "status": "ready", "lat": 52.5, "lon": 13.4,
            "variable": "pr", "analysis": ["etccdi_indices"],
            "indices": ["MegaDrought99"],
        })


def test_dataset_missing_variable_rejected():
    # POWER does not provide soil moisture.
    with pytest.raises(ValidationError):
        validate_response({
            "status": "ready", "lat": 52.5, "lon": 13.4, "dataset": "POWER",
            "variable": "soilMoist", "analysis": ["climatology"],
        })


def test_unknown_intent_rejected():
    with pytest.raises(ValidationError):
        validate_response({
            "status": "ready", "lat": 52.5, "lon": 13.4,
            "variable": "tasmax", "analysis": ["teleport"],
        })


def test_empty_analysis_rejected():
    with pytest.raises(ValidationError):
        validate_response({
            "status": "ready", "lat": 52.5, "lon": 13.4,
            "variable": "tasmax", "analysis": [],
        })


def test_extra_fields_rejected():
    with pytest.raises(ValidationError):
        validate_response({
            "status": "ready", "lat": 52.5, "lon": 13.4,
            "variable": "tasmax", "analysis": ["historical_trend"],
            "make_me_coffee": True,
        })


def test_dataset_named_ok():
    r = validate_response({
        "status": "ready", "lat": 52.5, "lon": 13.4, "dataset": "HYRAS",
        "variable": "tasmax", "analysis": ["historical_trend"],
    })
    assert r.dataset == "HYRAS"


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-q"]))
