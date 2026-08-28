"""Shared pytest fixtures for the LLM layer.

``climdata.LLM`` is a real subpackage, so the modules under test are imported
normally (``from climdata.LLM.planner import ...``) and no sys.path
manipulation is needed here.
"""

from __future__ import annotations

import json

import numpy as np
import pytest
import xarray as xr


# ===========================================================================
# Synthetic climate data
# ===========================================================================

def make_series(start="1980-01-01", end="2014-12-31", base=15.0, warming=0.04,
                noise=2.0, seed=0, name="tasmax", spatial=False):
    """Daily series with a linear warming signal, seasonal cycle and noise."""
    time = xr.cftime_range(start, end, freq="D", calendar="standard")
    t = np.arange(len(time))
    years = np.array([d.year for d in time])
    rng = np.random.default_rng(seed)
    values = (base + 10 * np.sin(2 * np.pi * t / 365.25)
              + warming * (years - years[0]) + rng.normal(0, noise, len(time)))
    if spatial:
        # broadcast onto a small lat/lon grid with a mild spatial gradient
        lat, lon = np.array([52.0, 52.5]), np.array([13.0, 13.5])
        grid = values[:, None, None] + np.array([[0.0, 0.5], [-0.5, 0.0]])[None]
        da = xr.DataArray(grid, coords={"time": time, "lat": lat, "lon": lon},
                          dims=("time", "lat", "lon"), name=name)
    else:
        da = xr.DataArray(values, coords={"time": time}, dims="time", name=name)
    da.attrs["units"] = "degC"
    return da


@pytest.fixture
def obs():
    return make_series(seed=1)


@pytest.fixture
def obs_grid():
    return make_series(seed=1, spatial=True)


# ===========================================================================
# Fake OpenAI-compatible LLM client
# ===========================================================================

class FakeLLM:
    """Minimal stand-in for the SAIA / OpenAI client.

    ``plan_replies`` and ``narrate_replies`` are consumed one per call, so a
    retry/repair loop can be driven deterministically. The last reply repeats
    once the list is exhausted.
    """

    def __init__(self, plan_replies=None, narrate_replies=None, raise_on=None):
        self.plan_replies = list(plan_replies or [])
        self.narrate_replies = list(narrate_replies or ["No numbers here."])
        self.raise_on = raise_on or set()
        self.plan_calls = 0
        self.narrate_calls = 0
        self.messages_seen = []
        # every call's kwargs, model included — without this the model id is
        # discarded and no test can assert which LLM a stage actually used.
        self.kwargs_seen = []
        self.chat = type("_Chat", (), {"completions": self})()

    @property
    def models_seen(self):
        return [kw.get("model") for kw in self.kwargs_seen]

    # the client surface used by planner/narrator: client.chat.completions.create
    def create(self, model=None, messages=None, **kwargs):
        self.messages_seen.append(messages)
        self.kwargs_seen.append({"model": model, "messages": messages, **kwargs})
        system = messages[0]["content"]
        kind = "plan" if "PLANNER" in system else "narrate"
        if kind in self.raise_on:
            raise RuntimeError(f"simulated {kind} transport failure")
        if kind == "plan":
            replies, self.plan_calls = self.plan_replies, self.plan_calls + 1
            n = self.plan_calls
        else:
            replies, self.narrate_calls = self.narrate_replies, self.narrate_calls + 1
            n = self.narrate_calls
        content = replies[min(n, len(replies)) - 1] if replies else "{}"
        if isinstance(content, (dict, list)):
            content = json.dumps(content)
        msg = type("_Msg", (), {"content": content})
        return type("_Resp", (), {"choices": [type("_C", (), {"message": msg})]})


@pytest.fixture
def fake_llm():
    return FakeLLM


# ===========================================================================
# Canonical planner payloads
# ===========================================================================

@pytest.fixture
def ready_point_payload():
    return {
        "status": "ready", "lat": 52.52, "lon": 13.405,
        "variable": "tasmax",
        "time_range": {"start": "1980-01-01", "end": "2014-12-31"},
        "analysis": ["historical_trend"],
    }
