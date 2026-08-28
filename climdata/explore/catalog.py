"""Object-oriented front end to the dataset registry."""

from __future__ import annotations

from climdata.explore.registry import REGISTRY, resolve_dataset_key
from climdata.explore.queries import (
    list_available_data,
    explore,
    find,
    inspect,
)


class DatasetRegistry:
    """A browsable handle on the climdata dataset catalogue.

    Wraps the module-level functions of :mod:`climdata.explore.queries` so the
    catalogue can be held in a variable, and adds dict-style lookup for callers
    who want the metadata rather than a printed report. Echoing the object at a
    REPL prints the summary table.

    Attributes:
        keys: Bound method returning the canonical dataset abbreviations.

    Example:
        >>> from climdata.explore import DatasetRegistry
        >>> reg = DatasetRegistry()
        >>> reg                                    # doctest: +SKIP
        AVAILABLE CLIMATE DATASETS ...
        >>> reg.explore("ERA5")                    # doctest: +SKIP
        >>> reg["era5"]["resolution"]              # dict-style, case-insensitive
        '0.25°'
    """

    def __init__(self) -> None:
        """Bind the process-wide registry snapshot."""
        self._registry = REGISTRY

    # ------------------------------------------------------------------
    def __repr__(self) -> str:
        """Render the catalogue as a summary table.

        Returns:
            str: Abbreviation, long name and type for every dataset.
        """
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
        """Print the summary table. See :func:`climdata.explore.list_available_data`."""
        list_available_data()

    def explore(self, dataset: str) -> None:
        """Print a dataset profile. See :func:`climdata.explore.explore`.

        Args:
            dataset (str): Dataset abbreviation, case-insensitive.
        """
        explore(dataset)

    def find(self, **kwargs) -> None:
        """Search the catalogue. See :func:`climdata.explore.find`.

        Args:
            **kwargs: Any of ``variable``, ``frequency``, ``type_filter`` or
                ``coverage``.
        """
        find(**kwargs)

    def inspect(self, dataset: str, variable: str) -> None:
        """Print variable metadata. See :func:`climdata.explore.inspect`.

        Args:
            dataset (str): Dataset abbreviation, case-insensitive.
            variable (str): CF variable name.
        """
        inspect(dataset, variable)

    # Allow dict-like access for power users
    def __getitem__(self, key: str) -> dict:
        """Return the raw metadata dict for one dataset.

        Args:
            key (str): Dataset abbreviation, case-insensitive.

        Returns:
            dict: The registry entry, as described in
            :func:`climdata.explore.get_registry`.

        Raises:
            KeyError: If no dataset matches ``key``.
        """
        resolved = resolve_dataset_key(key)
        if resolved is None:
            raise KeyError(f"Dataset '{key}' not found in registry.")
        return self._registry[resolved]

    def keys(self):
        """Return the canonical dataset abbreviations.

        Returns:
            KeysView[str]: The registry keys, e.g. ``"ERA5"``, ``"NEXGDDP"``.
        """
        return self._registry.keys()
