"""The dataset registry, built from the Hydra configuration at import time.

Providers are discovered from ``conf/mappings/parameters.yaml`` rather than
listed in code, so adding a dataset to the configuration is enough to make it
appear in :func:`climdata.explore.list_available_data` and friends.

The registry is a snapshot: :data:`REGISTRY` is built once when this module is
first imported. Editing ``parameters.yaml`` in a live session has no effect
until the interpreter is restarted.
"""

from __future__ import annotations
from typing import Optional
from pathlib import Path
import yaml


def _load_parameters_yaml() -> dict:
    """Read ``conf/mappings/parameters.yaml`` from the installed package.

    Reads the packaged copy next to this module, not the working copy that
    :func:`climdata.utils.config._ensure_local_conf` drops into the current
    directory — registry metadata should not change with the caller's cwd.

    Returns:
        dict: Parsed YAML, or an empty dict if the file is missing, empty or
        unreadable. Never raises: a broken config degrades the catalogue to
        empty rather than breaking ``import climdata``.
    """
    try:
        conf_dir = Path(__file__).parent.parent / "conf" / "mappings"
        params_file = conf_dir / "parameters.yaml"

        if not params_file.exists():
            return {}

        with open(params_file, "r") as f:
            return yaml.safe_load(f) or {}
    except Exception:
        return {}


def _load_registry() -> dict[str, dict]:
    """Build the dataset registry from configuration plus static entries.

    Every top-level mapping in ``parameters.yaml`` is treated as a dataset and
    registered under its upper-cased key (``dwd`` → ``DWD``, ``cmip_w5e5`` →
    ``CMIP_W5E5``), so newly configured providers appear without code changes.
    Descriptive fields come from each entry's optional ``explore`` block and fall
    back to ``"Unknown"``; variables come from its ``variables`` mapping.

    ERA5 is appended afterwards as a static entry: it is reached through the
    Copernicus CDS API rather than through ``parameters.yaml``, so it has no
    configuration block to discover.

    Returns:
        dict[str, dict]: Mapping of dataset abbreviation to a metadata dict with
        keys ``full_name``, ``type``, ``coverage``, ``resolution``,
        ``frequency``, ``time_range``, ``source``, ``notes``, ``variables``,
        ``experiments`` and ``models``.
    """
    params = _load_parameters_yaml()
    registry = {}

    # Extract every dataset from config (auto-discovered, upper-cased keys)
    for config_key, dataset_cfg in params.items():
        if not isinstance(dataset_cfg, dict):
            continue

        registry_key = config_key.upper()

        # Get explore metadata if available
        explore_meta = dataset_cfg.get("explore", {}) or {}

        # Extract variables from config
        variables = list((dataset_cfg.get("variables", {}) or {}).keys())

        # Extract experiments if available
        experiments = []
        if "experiment_id" in dataset_cfg:
            exp = dataset_cfg["experiment_id"]
            experiments = [exp] if isinstance(exp, str) else list(exp)

        # Build registry entry
        registry[registry_key] = {
            "full_name": explore_meta.get("full_name", config_key),
            "type": explore_meta.get("type", "Unknown"),
            "coverage": explore_meta.get("coverage", "Unknown"),
            "resolution": explore_meta.get("resolution", "Unknown"),
            "frequency": explore_meta.get("frequency", "Unknown"),
            "time_range": explore_meta.get("time_range", "Unknown"),
            "source": explore_meta.get("source", "Unknown"),
            "notes": explore_meta.get("notes", ""),
            "variables": variables,
            "experiments": experiments,
            "models": [],
        }

    # Add static datasets not in config files (ERA5)
    static_datasets = {
        "ERA5": {
            "full_name": "ECMWF Reanalysis v5",
            "type": "Observation",
            "coverage": "Global",
            "resolution": "0.25°",
            "frequency": "hourly / daily",
            "time_range": "1940-01-01 to present",
            "source": "Copernicus CDS (cdsapi)",
            "variables": ["tas", "tasmax", "tasmin", "pr", "tp", "hurs", "sfcWind", "rsds", "rlds"],
            "experiments": [],
            "models": [],
            "notes": "Requires a valid ~/.cdsapirc key.",
        },
    }

    registry.update(static_datasets)
    return registry


#: Dataset catalogue, built once at import time by :func:`_load_registry`.
REGISTRY = _load_registry()


def get_registry() -> dict[str, dict]:
    """Return the dataset registry for programmatic use.

    Returns:
        dict[str, dict]: A shallow copy of :data:`REGISTRY`, keyed by dataset
        abbreviation. The copy protects the top level only — the per-dataset
        metadata dicts are shared, so do not mutate them.

    Example:
        >>> from climdata.explore import get_registry
        >>> sorted(get_registry())[:3]                     # doctest: +SKIP
        ['AGRI_ISIMIP', 'CMIP', 'CMIP_W5E5']
    """
    return REGISTRY.copy()


def resolve_dataset_key(name: str) -> Optional[str]:
    """Resolve a dataset name to its canonical registry key, case-insensitively.

    Args:
        name (str): Dataset name in any casing, e.g. ``"era5"`` or ``"ERA5"``.

    Returns:
        str | None: The canonical key (``"ERA5"``), or ``None`` if no dataset
        matches. Callers report the miss themselves rather than catching an
        exception.

    Example:
        >>> from climdata.explore import resolve_dataset_key
        >>> resolve_dataset_key("nexgddp")
        'NEXGDDP'
        >>> resolve_dataset_key("not-a-dataset") is None
        True
    """
    name_up = name.upper()
    for k in REGISTRY:
        if k.upper() == name_up:
            return k
    return None
