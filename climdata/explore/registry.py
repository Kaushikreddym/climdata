"""
Dynamic dataset registry loaded from Hydra configuration.
"""

from __future__ import annotations
from typing import Optional
from pathlib import Path
import yaml


def _load_parameters_yaml() -> dict:
    """Load parameters.yaml from climdata/conf/mappings/."""
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
    """
    Load dataset registry directly from YAML configuration files,
    then add static/hardcoded datasets.
    
    Returns a dictionary mapping dataset abbreviations to metadata dictionaries.
    """
    params = _load_parameters_yaml()
    registry = {}
    
    # Map config keys to uppercase registry keys
    key_mapping = {
        "dwd": "DWD",
        "mswx": "MSWX",
        "hyras": "HYRAS",
        "cmip": "CMIP",
        "power": "NASAPOWER",
        "w5e5": "W5E5",
        "cmip_w5e5": "CMIPW5E5",
        "nexgddp": "NEXGDDP",
    }
    
    # Extract datasets from config
    for config_key, registry_key in key_mapping.items():
        if config_key in params:
            dataset_cfg = params[config_key]
            
            # Get explore metadata if available
            explore_meta = dataset_cfg.get("explore", {})
            
            # Extract variables from config
            variables = list(dataset_cfg.get("variables", {}).keys())
            
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


# Load registry once at module import time
REGISTRY = _load_registry()


def get_registry() -> dict[str, dict]:
    """Get the loaded dataset registry."""
    return REGISTRY.copy()


def resolve_dataset_key(name: str) -> Optional[str]:
    """Case-insensitive lookup; returns canonical key or None."""
    name_up = name.upper()
    for k in REGISTRY:
        if k.upper() == name_up:
            return k
    return None
