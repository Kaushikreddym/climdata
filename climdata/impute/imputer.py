import importlib
from pathlib import Path
from typing import List, Dict

import xarray as xr

class Imputer:
    """
    Imputation manager class.

    Usage:
        imputer = Imputer(cfg, ds)
        out = imputer.calculate("simple")  # runs the named imputation config
        imputer.run()  # runs all configured imputations

    The class expects imputation definitions under `cfg.mappings.imputes` (dict),
    where each impute entry looks like:
        name:
            function: "module.path.func"
            variables: ["tas", "tasmin"]
            kwargs: {...}
            postprocess: { "sel": {"lat": 52.0} }
    """
    def __init__(self, cfg, climate_data: xr.Dataset):
        self.cfg = cfg
        self.climate_data = climate_data
        self.successful = []
        self.failed = []
        self.tasks: List[str] = []

    def _get_imputes_map(self) -> Dict:
        """Return imputes mapping from config, support common key names."""
        if hasattr(self.cfg, "mappings") and isinstance(self.cfg.mappings, dict):
            return self.cfg.mappings.get("imputes", {})
        # fallback
        return getattr(self.cfg, "imputes", {}) or {}

    def calculate(self, name: str) -> xr.Dataset:
        """
        Perform the named imputation and return the imputed Dataset.

        Raises:
            KeyError if the named imputation is not found.
            Exception if the underlying imputation function fails.
        """
        imputes = self._get_imputes_map()
        if name not in imputes:
            raise KeyError(f"Imputation '{name}' not found in configuration")

        impute_cfg = imputes[name]
        func_path = impute_cfg["function"]
        module_name, func_name = func_path.rsplit(".", 1)
        func = getattr(importlib.import_module(module_name), func_name)

        variables = impute_cfg.get("variables", list(self.climate_data.data_vars))
        kwargs = dict(impute_cfg.get("kwargs", {}))

        # Call imputation function (common signature: ds, vars, **kwargs)
        result = func(self.climate_data, variables, **kwargs)

        # Optional postprocessing: simple 'sel' and 'rename' support
        post = impute_cfg.get("postprocess", {})
        if "sel" in post:
            result = result.sel(**post["sel"])
        if "rename" in post:
            result = result.rename(post["rename"])

        # If returned a DataArray, put into Dataset with the imputation name
        if isinstance(result, xr.DataArray):
            ds_out = result.to_dataset(name=name)
        else:
            ds_out = result

        # Ensure the output is a Dataset
        if not isinstance(ds_out, xr.Dataset):
            raise TypeError("Imputation function must return xr.Dataset or xr.DataArray")

        return ds_out

    def run(self):
        """Run all configured imputations, collect successes and failures."""
        imputes = self._get_imputes_map()
        for name in list(imputes.keys()):
            try:
                self.calculate(name)
                self.successful.append(name)
            except Exception as exc:
                self.failed.append(name)
                # Keep going, but log/print the error
                print(f"Imputer: error running '{name}': {exc}")