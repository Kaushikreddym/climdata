"""
climdata.explore
================
Dataset discovery and exploration utilities.

Provides a fluent interface for exploring climdata's datasets without 
needing to know the underlying file paths. Dataset metadata is loaded 
from climdata/conf/mappings/ (Hydra configuration).

Public API (module-level functions):
    - list_available_data()           – print summary table
    - explore(dataset="ERA5")         – deep-dive into a dataset
    - find(variable="pr", ...)        – filter by criteria
    - inspect(dataset="ERA5", var)    – variable-level metadata + BASD hints
    - list_esm_experiments(dataset)   – ESM experiments/scenarios available
    - list_esm_models(dataset, ...)   – ESM models available

Alternative: Object-Oriented Interface
    - DatasetRegistry                 – class-based access
"""

from climdata.explore.queries import (
    list_available_data,
    explore,
    find,
    inspect,
    list_esm_experiments,
    list_esm_models,
)
from climdata.explore.catalog import DatasetRegistry
from climdata.explore.registry import get_registry, resolve_dataset_key

__all__ = [
    "list_available_data",
    "explore",
    "find",
    "inspect",
    "list_esm_experiments",
    "list_esm_models",
    "DatasetRegistry",
    "get_registry",
    "resolve_dataset_key",
]
