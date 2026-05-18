import os
import sys
import shutil
from pathlib import Path
from hydra import initialize_config_dir, compose
from hydra.core.global_hydra import GlobalHydra
from omegaconf import OmegaConf
import importlib.resources as resources

def _ensure_local_conf(package="climdata", local_dir="conf"):
    """
    Copy package conf/ to cwd if not exists.
    Returns the full path to the conf directory.
    
    Raises:
        FileNotFoundError: If conf directory cannot be found in package.
    """
    # PyInstaller bundle: conf is already extracted under _MEIPASS/climdata/conf.
    # os.getcwd() inside a Finder-launched .app is '/' (read-only), so never
    # attempt a copy there.
    if getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS"):
        bundled = Path(sys._MEIPASS) / package / "conf"
        if bundled.exists():
            return str(bundled)
        raise FileNotFoundError(
            f"Bundled conf not found at {bundled}. "
            "Ensure collect_all('climdata') is in the .spec file."
        )

    local_dir_path = Path(os.getcwd()) / local_dir
    if not local_dir_path.exists():
        conf_src = resources.files(package).joinpath("conf")
        shutil.copytree(conf_src, local_dir_path)
    return str(local_dir_path)

def load_config(config_name="config", overrides=None, verbose=False):
    """
    Load Hydra config from ./conf directory.
    
    Uses initialize_config_dir with absolute path to ensure Hydra finds the config.
    
    Args:
        config_name (str): Name of the config file (default: "config")
        overrides (list[str], optional): Hydra overrides to apply
        verbose (bool): Print the loaded config to stdout
        
    Returns:
        DictConfig: Hydra configuration object
        
    Raises:
        FileNotFoundError: If conf directory cannot be found
        RuntimeError: If Hydra initialization or config composition fails
    """
    # Ensure conf directory exists locally
    conf_dir_path = _ensure_local_conf()
    conf_dir_abs = os.path.abspath(conf_dir_path)
    
    # Reset Hydra if it was already initialized
    GlobalHydra.instance().clear()
    
    # Initialize and compose config using absolute path
    with initialize_config_dir(config_dir=conf_dir_abs, version_base=None):
        cfg = compose(config_name=config_name, overrides=overrides or [])
        if verbose:
            print(OmegaConf.to_yaml(cfg))
        return cfg