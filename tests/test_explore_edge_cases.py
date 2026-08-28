"""Edge cases for :mod:`climdata.explore`.

The discovery layer is the first thing a new user touches, and its contract is
unusual: the query functions *print* and return ``None``, and they report a bad
name rather than raising. These tests pin that behaviour down, because a
lookup miss that starts raising would break exploratory notebook sessions.

Everything here is offline — the registry is built from the packaged YAML at
import time and never touches the network.
"""

import pytest

from climdata.explore import (
    DatasetRegistry,
    explore,
    find,
    get_registry,
    inspect,
    list_available_data,
    list_esm_experiments,
    list_esm_models,
    resolve_dataset_key,
)
from climdata.explore.registry import REGISTRY, _load_parameters_yaml


# --------------------------------------------------------------------------
# Registry construction
# --------------------------------------------------------------------------
def test_registry_is_not_empty():
    """A registry built from a readable config must find the providers."""
    assert len(REGISTRY) > 0
    # ERA5 is appended statically, so it is present even if the YAML is empty.
    assert "ERA5" in REGISTRY


@pytest.mark.parametrize("key", sorted(REGISTRY))
def test_every_entry_has_the_full_metadata_shape(key):
    """Consumers index these keys unconditionally; a missing one is a KeyError."""
    required = {
        "full_name", "type", "coverage", "resolution", "frequency",
        "time_range", "source", "notes", "variables", "experiments", "models",
    }
    assert required <= set(REGISTRY[key])


def test_registry_keys_are_upper_cased():
    """Lookup is case-folded through upper(), so stored keys must be upper."""
    assert all(k == k.upper() for k in REGISTRY)


def test_parameters_yaml_never_raises(monkeypatch):
    """A broken config degrades the catalogue to empty, never breaks the import."""
    def boom(*a, **k):
        raise OSError("corrupt config")

    from climdata.explore import registry as registry_mod

    monkeypatch.setattr(registry_mod.yaml, "safe_load", boom)
    assert _load_parameters_yaml() == {}


# --------------------------------------------------------------------------
# resolve_dataset_key
# --------------------------------------------------------------------------
@pytest.mark.parametrize("name", ["era5", "ERA5", "Era5", "eRa5"])
def test_resolve_is_case_insensitive(name):
    assert resolve_dataset_key(name) == "ERA5"


def test_resolve_returns_none_for_unknown():
    """None, not an exception — callers report the miss themselves."""
    assert resolve_dataset_key("not-a-dataset") is None


def test_resolve_does_not_match_on_prefix():
    """Substring matching would make 'CMIP' ambiguous with 'CMIP_W5E5'."""
    assert resolve_dataset_key("CMIP_") is None


def test_resolve_empty_string():
    assert resolve_dataset_key("") is None


# --------------------------------------------------------------------------
# Query functions: unknown input prints, never raises
# --------------------------------------------------------------------------
def test_explore_unknown_dataset_prints_hint(capsys):
    explore("nonexistent")
    out = capsys.readouterr().out
    assert "not found in the registry" in out
    assert "ERA5" in out          # the hint lists what *is* available


def test_inspect_unknown_dataset_prints_hint(capsys):
    inspect("nonexistent", variable="pr")
    assert "not found in the registry" in capsys.readouterr().out


def test_inspect_unknown_variable_lists_available(capsys):
    inspect("ERA5", variable="not_a_variable")
    out = capsys.readouterr().out
    assert "is not available in ERA5" in out
    assert "Available variables" in out


def test_list_esm_experiments_unknown_dataset(capsys):
    list_esm_experiments("nonexistent")
    assert "not found in the registry" in capsys.readouterr().out


def test_list_esm_experiments_on_non_esm_dataset(capsys):
    """ERA5 is an observation product; asking for its scenarios is a no-op."""
    list_esm_experiments("ERA5")
    assert "is not an ESM dataset" in capsys.readouterr().out


def test_list_esm_models_on_non_esm_dataset(capsys):
    """Must short-circuit before any network call."""
    list_esm_models("ERA5")
    assert "is not an ESM dataset" in capsys.readouterr().out


def test_list_available_data_prints_every_dataset(capsys):
    list_available_data()
    out = capsys.readouterr().out
    for key in REGISTRY:
        assert key in out
    assert f"Total: {len(REGISTRY)} datasets" in out


# --------------------------------------------------------------------------
# find()
# --------------------------------------------------------------------------
def test_find_with_no_criteria_lists_everything(capsys):
    find()
    assert f"{len(REGISTRY)} result(s)" in capsys.readouterr().out


def test_find_impossible_combination_reports_no_match(capsys):
    find(variable="pr", coverage="Mars")
    assert "No datasets matched" in capsys.readouterr().out


def test_find_unknown_variable_matches_nothing(capsys):
    find(variable="not_a_variable")
    assert "No datasets matched" in capsys.readouterr().out


def test_find_criteria_are_anded(capsys):
    """Each extra criterion can only narrow the result, never widen it."""
    find(type_filter="Observation")
    broad = capsys.readouterr().out
    find(type_filter="Observation", coverage="Germany")
    narrow = capsys.readouterr().out

    def n(text):
        return int(text.split("result(s)")[0].strip().split()[-1])

    assert n(narrow) <= n(broad)


def test_find_substring_matching_is_case_insensitive(capsys):
    """Only the echoed query differs between casings; the matches must not."""
    def matched(text):
        return [l for l in text.splitlines() if l.strip().startswith(tuple("0123456789"))]

    find(type_filter="observation")
    lower = matched(capsys.readouterr().out)
    find(type_filter="OBSERVATION")
    upper = matched(capsys.readouterr().out)
    assert lower and lower == upper


def test_find_variable_match_is_exact_not_substring(capsys):
    """'p' must not match 'pr' — that would make the filter useless."""
    find(variable="p")
    assert "No datasets matched" in capsys.readouterr().out


# --------------------------------------------------------------------------
# get_registry
# --------------------------------------------------------------------------
def test_get_registry_top_level_is_a_copy():
    reg = get_registry()
    reg["FAKE"] = {}
    assert "FAKE" not in REGISTRY


def test_get_registry_is_only_a_shallow_copy():
    """Documented behaviour: nested dicts are shared, so callers must not mutate."""
    reg = get_registry()
    assert reg["ERA5"] is REGISTRY["ERA5"]


# --------------------------------------------------------------------------
# DatasetRegistry
# --------------------------------------------------------------------------
def test_registry_getitem_is_case_insensitive():
    reg = DatasetRegistry()
    assert reg["era5"] == reg["ERA5"]


def test_registry_getitem_raises_on_unknown():
    """Dict-style access is the one place a miss raises rather than prints."""
    with pytest.raises(KeyError, match="not found in registry"):
        DatasetRegistry()["nonexistent"]


def test_registry_keys_match_module_registry():
    assert set(DatasetRegistry().keys()) == set(REGISTRY)


def test_registry_repr_is_a_table():
    text = repr(DatasetRegistry())
    assert "AVAILABLE CLIMATE DATASETS" in text
    for key in REGISTRY:
        assert key in text


def test_registry_methods_delegate(capsys):
    reg = DatasetRegistry()
    reg.explore("ERA5")
    assert "EXPLORING: ERA5" in capsys.readouterr().out
    reg.find(variable="pr")
    assert "MATCHING DATASETS" in capsys.readouterr().out
