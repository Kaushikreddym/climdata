"""Console-formatted queries over the climdata dataset registry.

Every function here prints a human-readable report and returns ``None`` — they
exist for interactive use at a REPL or in a notebook, to answer "what data can I
get, and what is in it?" before any download is attempted. For programmatic
access to the same metadata use :func:`climdata.explore.get_registry`, which
returns a plain dictionary.
"""

from __future__ import annotations
from typing import Optional

from climdata.explore.registry import REGISTRY, resolve_dataset_key


def list_available_data() -> None:
    """Print a summary table of every dataset in the registry.

    One row per provider, with its abbreviation, long name, type, spatial
    coverage and native resolution. The abbreviation in the first column is the
    key accepted by :func:`explore`, :func:`inspect` and the ``dataset=`` Hydra
    override.

    Returns:
        None: The table is written to stdout.

    Example:
        >>> import climdata as cd
        >>> cd.list_available_data()                       # doctest: +SKIP
    """
    col_abbr  = 12
    col_name  = 44
    col_type  = 26
    col_cov   = 10
    col_res   = 14
    sep = "-" * (col_abbr + col_name + col_type + col_cov + col_res + 12)

    header = (
        f"{'Abbr.':<{col_abbr}} | "
        f"{'Long Name':<{col_name}} | "
        f"{'Type':<{col_type}} | "
        f"{'Coverage':<{col_cov}} | "
        f"{'Resolution':<{col_res}}"
    )

    print("\nAVAILABLE CLIMATE DATASETS")
    print(sep)
    print(header)
    print(sep)
    for key, meta in REGISTRY.items():
        print(
            f"{key:<{col_abbr}} | "
            f"{meta['full_name']:<{col_name}} | "
            f"{meta['type']:<{col_type}} | "
            f"{meta['coverage']:<{col_cov}} | "
            f"{meta['resolution']:<{col_res}}"
        )
    print(sep)
    print(f"  Total: {len(REGISTRY)} datasets  |  Use cd.explore(dataset=...) for details.\n")


def explore(dataset: str) -> None:
    """Print a detailed profile for one dataset.

    Covers type, coverage, resolution, frequency, temporal extent and source,
    followed by the experiments, models and variables the registry knows about.
    An unknown name prints the list of valid ones rather than raising.

    Args:
        dataset (str): Dataset abbreviation as shown by :func:`list_available_data`
            (e.g. ``"ERA5"``, ``"NEXGDDP"``, ``"CMIP"``). Matched case-insensitively.

    Returns:
        None: The profile is written to stdout.

    Example:
        >>> import climdata as cd
        >>> cd.explore(dataset="NEXGDDP")                  # doctest: +SKIP
    """
    key = resolve_dataset_key(dataset)
    if key is None:
        _unknown_dataset_hint(dataset)
        return

    meta = REGISTRY[key]
    sep  = "=" * 68
    sep2 = "-" * 68

    print(f"\n{sep}")
    print(f"  EXPLORING: {key}  —  {meta['full_name']}")
    print(sep)
    print(f"  Type        : {meta['type']}")
    print(f"  Coverage    : {meta['coverage']}")
    print(f"  Resolution  : {meta['resolution']}")
    print(f"  Frequency   : {meta['frequency']}")
    print(f"  Time range  : {meta['time_range']}")
    print(f"  Data source : {meta['source']}")

    if meta.get("experiments"):
        exps = ", ".join(meta["experiments"])
        print(f"\n  Experiments ({len(meta['experiments'])}): [{exps}]")

    if meta.get("models"):
        print(f"\n  Models available ({len(meta['models'])}):")
        for m in meta["models"]:
            print(f"    • {m}")

    if meta.get("variables"):
        print(f"\n  Variables stored:")
        print(sep2)
        print(f"  {'Variable':<12} {'Long Name':<46} {'Units'}")
        print(sep2)
        for v in meta["variables"]:
            long  = v  # Variable short name as fallback
            units = "—"
            print(f"  {v:<12} {long:<46} {units}")
        print(sep2)

    if meta.get("notes"):
        print(f"\n  ℹ  Note: {meta['notes']}")

    print(f"\n  Tip: cd.inspect('{key}', variable='<var>') for unit/BASD details.\n")


def find(
    variable:    Optional[str] = None,
    frequency:   Optional[str] = None,
    type_filter: Optional[str] = None,
    coverage:    Optional[str] = None,
) -> None:
    """Search the registry for datasets matching the given criteria.

    Criteria combine with AND: a dataset is listed only if it satisfies every
    argument that was supplied. ``variable`` requires an exact CF-name match;
    the other three match as case-insensitive substrings. Passing no arguments
    lists everything.

    Args:
        variable (str, optional): CF variable name that the dataset must provide
            (e.g. ``"pr"``, ``"tas"``). Exact match.
        frequency (str, optional): Temporal resolution keyword, e.g. ``"daily"``,
            ``"hourly"``. Substring match.
        type_filter (str, optional): Dataset type, e.g. ``"Observation"``,
            ``"ESM"``. Substring match.
        coverage (str, optional): Coverage region, e.g. ``"Global"``,
            ``"Germany"``. Substring match.

    Returns:
        None: Matches are written to stdout.

    Example:
        >>> import climdata as cd
        >>> cd.find(variable="pr", frequency="daily")      # doctest: +SKIP
        >>> cd.find(type_filter="ESM")                     # doctest: +SKIP
    """
    matches = []
    for key, meta in REGISTRY.items():
        if variable and variable not in meta.get("variables", []):
            continue
        if frequency and frequency.lower() not in meta["frequency"].lower():
            continue
        if type_filter and type_filter.lower() not in meta["type"].lower():
            continue
        if coverage and coverage.lower() not in meta["coverage"].lower():
            continue
        matches.append((key, meta))

    if not matches:
        print("\n  No datasets matched the given criteria.\n")
        return

    # Build a compact query description for the header
    criteria = []
    if variable:    criteria.append(f"variable='{variable}'")
    if frequency:   criteria.append(f"frequency='{frequency}'")
    if type_filter: criteria.append(f"type='{type_filter}'")
    if coverage:    criteria.append(f"coverage='{coverage}'")
    query_str = ", ".join(criteria) if criteria else "all"

    print(f"\nMATCHING DATASETS  (query: {query_str})")
    print("-" * 68)
    for i, (key, meta) in enumerate(matches, 1):
        models_info = f"  ({len(meta['models'])} models)" if meta.get("models") else ""
        print(
            f"  {i:>2}. {key:<12} ({meta['type']:<26})  "
            f"{meta['coverage']}, {meta['resolution']}{models_info}"
        )
    print("-" * 68)
    print(f"  {len(matches)} result(s).  Use cd.explore(dataset=...) for full details.\n")


def inspect(dataset: str, variable: str) -> None:
    """Print variable-level metadata for one variable of one dataset.

    Reports the long name, the unit BASD expects for the variable, and any
    conversion note attached to it — the detail that matters before feeding data
    to :mod:`climdata.sdba`, since a unit mismatch there fails silently rather
    than loudly. Metadata comes from ``conf/mappings/parameters.yaml``.

    Unknown datasets, and variables the dataset does not carry, print the valid
    alternatives rather than raising.

    Args:
        dataset (str): Dataset abbreviation (e.g. ``"ERA5"``). Case-insensitive.
        variable (str): CF variable name (e.g. ``"tp"``, ``"pr"``).

    Returns:
        None: The report is written to stdout.

    Example:
        >>> import climdata as cd
        >>> cd.inspect("ERA5", variable="tp")              # doctest: +SKIP
    """
    key = resolve_dataset_key(dataset)
    if key is None:
        _unknown_dataset_hint(dataset)
        return

    meta = REGISTRY[key]

    if variable not in meta.get("variables", []):
        available = ", ".join(meta.get("variables", []))
        print(
            f"\n  ✗  Variable '{variable}' is not available in {key}.\n"
            f"     Available variables: [{available}]\n"
        )
        return

    # Try to load metadata from YAML config first
    vmeta = _get_variable_metadata(variable)

    if vmeta is None:
        print(f"\n  Variable '{variable}' found in {key} but has no detailed metadata yet.\n")
        return

    sep = "-" * 68
    print(f"\n{sep}")
    print(f"  VARIABLE INSPECTION: {variable}  in  {key} ({meta['full_name']})")
    print(sep)
    print(f"  Long name  : {vmeta['long_name']}")
    print(f"  BASD unit  : {vmeta['basd_unit']}")

    if vmeta.get("basd_note"):
        print(f"  Conversion : ⚠  {vmeta['basd_note']}")
    else:
        print("  Conversion : ✓  No special conversion needed for BASD.")

    print(sep + "\n")


# -----------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------

def list_esm_experiments(dataset: str) -> None:
    """Print the experiments (scenarios) available for an ESM dataset.

    Reads the experiment list straight from the registry, so it costs nothing and
    needs no credentials. Non-ESM datasets, and ESM datasets whose experiments
    are not declared in the configuration, print a hint pointing at the dataset
    class's own ``get_experiment_ids()``.

    Args:
        dataset (str): Dataset abbreviation (e.g. ``"CMIP"``, ``"NEXGDDP"``,
            ``"CMIP_W5E5"``). Case-insensitive.

    Returns:
        None: The list is written to stdout.

    Example:
        >>> import climdata as cd
        >>> cd.list_esm_experiments("CMIP")                # doctest: +SKIP
    """
    key = resolve_dataset_key(dataset)
    if key is None:
        _unknown_dataset_hint(dataset)
        return

    meta = REGISTRY[key]

    # Check if dataset is ESM-based
    if "ESM" not in meta.get("type", ""):
        print(f"\n  ℹ  {key} is not an ESM dataset (Type: {meta.get('type')}).\n")
        return

    experiments = meta.get("experiments", [])

    if not experiments:
        print(f"\n  ℹ  No experiment metadata found for {key}.\n"
              f"     Try accessing the dataset directly:\n"
              f"       from climdata.datasets import {key.upper()}\n"
              f"       ds = {key.upper()}(cfg)\n"
              f"       experiments = ds.get_experiment_ids()\n")
        return

    sep = "-" * 60
    print(f"\n{sep}")
    print(f"  EXPERIMENTS AVAILABLE IN: {key}")
    print(sep)
    for i, exp in enumerate(experiments, 1):
        print(f"  {i}. {exp}")
    print(sep + "\n")


def list_esm_models(dataset: str, experiment: Optional[str] = None) -> None:
    """Print the ESM models available for a dataset, optionally for one experiment.

    Unlike :func:`list_esm_experiments`, this queries the provider's live
    catalogue through the dataset class, so it needs network access and can be
    slow. Failures at any step — unimportable class, unreachable catalogue,
    unknown experiment — are reported and swallowed rather than raised, so an
    exploratory call never aborts a session.

    Args:
        dataset (str): Dataset abbreviation (e.g. ``"CMIP"``, ``"NEXGDDP"``,
            ``"CMIP_W5E5"``). Case-insensitive.
        experiment (str, optional): Experiment ID to filter by, e.g.
            ``"historical"`` or ``"ssp585"``. Defaults to the provider's first
            experiment.

    Returns:
        None: The list is written to stdout.

    Example:
        >>> import climdata as cd
        >>> cd.list_esm_models("NEXGDDP")                          # doctest: +SKIP
        >>> cd.list_esm_models("CMIP", experiment="ssp585")        # doctest: +SKIP
    """
    key = resolve_dataset_key(dataset)
    if key is None:
        _unknown_dataset_hint(dataset)
        return

    meta = REGISTRY[key]

    # Check if dataset is ESM-based
    if "ESM" not in meta.get("type", ""):
        print(f"\n  ℹ  {key} is not an ESM dataset (Type: {meta.get('type')}).\n")
        return

    # Import the dataset class dynamically
    try:
        if key == "CMIP":
            from climdata.datasets.CMIPCloud import CMIPCloud as DatasetClass
        elif key == "NEXGDDP":
            from climdata.datasets.NEXGDDP import NEXGDDP as DatasetClass
        elif key == "CMIP_W5E5":
            from climdata.datasets.CMIP_W5E5 import CMIPW5E5 as DatasetClass
        else:
            print(f"\n  ✗  {key} has no direct model access method yet.\n")
            return
    except ImportError as e:
        print(f"\n  ✗  Could not import {key} dataset class: {e}\n")
        return

    # Instantiate the dataset class to access methods
    try:
        from omegaconf import DictConfig
        dummy_cfg = DictConfig({"experiment_id": "historical", "source_id": "dummy"})
        ds_instance = DatasetClass(dummy_cfg)
    except Exception as e:
        print(f"\n  ✗  Could not initialize {key}: {e}\n")
        return

    # Get experiment IDs if not provided
    if experiment is None:
        try:
            all_experiments = ds_instance.get_experiment_ids()
            experiment = all_experiments[0] if all_experiments else "historical"
        except Exception as e:
            print(f"\n  ✗  Could not retrieve experiments: {e}\n")
            return

    # Get models for the specified experiment
    try:
        models = ds_instance.get_source_ids(experiment)
        sep = "-" * 60
        print(f"\n{sep}")
        print(f"  ESM MODELS IN: {key}  (experiment: {experiment})")
        print(sep)
        for i, model in enumerate(models, 1):
            print(f"  {i}. {model}")
        print(sep)
        print(f"  Total: {len(models)} models\n")
    except Exception as e:
        print(f"\n  ✗  Could not retrieve models for {experiment}: {e}\n")


def _unknown_dataset_hint(name: str) -> None:
    """Print the available dataset keys after a lookup miss.

    Args:
        name (str): The name the caller asked for.
    """
    available = ", ".join(REGISTRY.keys())
    print(
        f"\n  ✗  Dataset '{name}' not found in the registry.\n"
        f"     Available datasets: [{available}]\n"
        f"     Run cd.list_available_data() to see the full catalogue.\n"
    )


def _get_variable_metadata(variable: str) -> Optional[dict]:
    """Look up BASD metadata for a variable across every dataset in the config.

    The same CF variable is declared under several providers in
    ``parameters.yaml``, and only some carry a ``basd_note``. Entries that have
    one are preferred, since a note is what makes the difference between a
    correct and a silently wrong BASD run.

    Args:
        variable (str): CF variable name, e.g. ``"pr"``.

    Returns:
        dict | None: Mapping with ``long_name``, ``basd_unit`` and ``basd_note``,
        or ``None`` if no declaration carries a ``basd_unit`` (or the config
        cannot be read).
    """
    try:
        from climdata.explore.registry import _load_parameters_yaml
        params = _load_parameters_yaml()

        # Search all datasets for this variable's BASD metadata
        # Prefer entries with basd_note (most complete metadata)
        best_match = None

        for dataset_cfg in params.values():
            if isinstance(dataset_cfg, dict) and "variables" in dataset_cfg:
                var_cfg = dataset_cfg["variables"].get(variable, {})
                if var_cfg and "basd_unit" in var_cfg:
                    candidate = {
                        "long_name": var_cfg.get("long_name", variable),
                        "basd_unit": var_cfg.get("basd_unit", "Unknown"),
                        "basd_note": var_cfg.get("basd_note"),
                    }
                    # Prefer matches with basd_note if available
                    if var_cfg.get("basd_note") and not (best_match and best_match.get("basd_note")):
                        best_match = candidate
                    elif not best_match:
                        best_match = candidate

        return best_match
    except Exception:
        pass

    return None
