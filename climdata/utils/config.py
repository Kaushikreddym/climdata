import os
import shutil
from pathlib import Path
from hydra import initialize, compose
from omegaconf import OmegaConf
import importlib.resources as resources

def _ensure_local_conf(package="climdata", local_dir="conf"):
    """
    Copy package conf/ to cwd if not exists.
    Returns the relative path "conf" for Hydra.
    """
    local_dir_path = Path(os.getcwd()) / local_dir
    if not local_dir_path.exists():
        # Get conf inside the installed package
        conf_src = resources.files(package).joinpath("conf")
        shutil.copytree(conf_src, local_dir_path)
    return local_dir  # relative for Hydra

def load_config(config_name="config", overrides=None, verbose=False):
    """
    Load Hydra config using ./conf in cwd.
    """
    config_path = _ensure_local_conf()
    print(config_path+config_name)
    with initialize(config_path=config_path, version_base=None):
        cfg = compose(config_name=config_name, overrides=overrides or [])
        if verbose:
            print(OmegaConf.to_yaml(cfg))
        return cfg