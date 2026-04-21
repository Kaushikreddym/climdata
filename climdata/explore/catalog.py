"""
Object-oriented interface for dataset browsing.
"""

from __future__ import annotations

from climdata.explore.registry import REGISTRY, resolve_dataset_key
from climdata.explore.queries import (
    list_available_data,
    explore,
    find,
    inspect,
)


class DatasetRegistry:
    """Object-oriented wrapper around the climdata dataset registry.

    Can be instantiated to get a printable catalogue object, or used
    programmatically via its methods.

    Examples
    --------
    >>> from climdata.explore import DatasetRegistry
    >>> reg = DatasetRegistry()
    >>> reg          # pretty-prints the summary table
    >>> reg.explore("ERA5")
    >>> reg.find(variable="pr")
    >>> reg.inspect("ERA5", variable="tp")
    """

    def __init__(self) -> None:
        self._registry = REGISTRY

    # ------------------------------------------------------------------
    def __repr__(self) -> str:
        col_abbr = 12
        col_name = 44
        col_type = 26
        sep = "-" * (col_abbr + col_name + col_type + 6)
        header = f"{'Abbr.':<{col_abbr}} | {'Long Name':<{col_name}} | {'Type':<{col_type}}"
        rows = [
            f"{k:<{col_abbr}} | {v['full_name']:<{col_name}} | {v['type']:<{col_type}}"
            for k, v in self._registry.items()
        ]
        return (
            "\nAVAILABLE CLIMATE DATASETS\n"
            + sep + "\n"
            + header + "\n"
            + sep + "\n"
            + "\n".join(rows) + "\n"
            + sep
        )

    # Delegate to module-level functions
    def list_available_data(self) -> None:
        list_available_data()

    def explore(self, dataset: str) -> None:
        explore(dataset)

    def find(self, **kwargs) -> None:
        find(**kwargs)

    def inspect(self, dataset: str, variable: str) -> None:
        inspect(dataset, variable)

    # Allow dict-like access for power users
    def __getitem__(self, key: str) -> dict:
        resolved = resolve_dataset_key(key)
        if resolved is None:
            raise KeyError(f"Dataset '{key}' not found in registry.")
        return self._registry[resolved]

    def keys(self):
        return self._registry.keys()
