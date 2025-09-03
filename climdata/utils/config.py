from hydra import initialize, compose
from omegaconf import OmegaConf
from typing import Optional, List
import os

def load_config(
    config_path: str = "../conf",
    config_name: str = "config",
    overrides: Optional[List[str]] = None,
    verbose: bool = False
):
    """
    Load a Hydra config file.

    Args:
        config_path (str): Path to the config directory.
        config_name (str): Name of the config YAML file (without `.yaml`).
        overrides (List[str], optional): List of override strings.
        verbose (bool): Whether to print the loaded config.

    Returns:
        OmegaConf.DictConfig: The loaded config object.
    """
    # config_path = os.path.abspath(config_path)

    with initialize(config_path=config_path, version_base=None):
        cfg = compose(config_name=config_name, overrides=overrides or [])
        if verbose:
            print(OmegaConf.to_yaml(cfg))
        return cfg
