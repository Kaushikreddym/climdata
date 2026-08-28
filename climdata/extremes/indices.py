"""Configuration-driven computation of climate extreme indices.

Indices are declared in ``conf/mappings/indices.yaml`` rather than in code. Each
entry names the callable to invoke (usually an :mod:`xclim` indicator), the
dataset variables to feed it, and any literal keyword arguments. Optional
``link`` blocks build intermediate inputs — a day-of-year percentile, say —
before the index itself is evaluated.

Example:
    >>> from climdata.extremes.indices import extreme_index
    >>> indices = extreme_index(cfg, ds)          # doctest: +SKIP
    >>> tn10p = indices.calculate("tn10p")        # doctest: +SKIP
"""

import importlib

import xarray as xr


class extreme_index:
    """Evaluate a configured extreme index against a loaded dataset.

    Attributes:
        cfg (DictConfig): Hydra configuration; ``cfg.extinfo.indices`` holds the
            index definitions.
        climate_data (xr.Dataset): Dataset supplying the input variables. Linked
            intermediate variables are written back into it as they are built.

    Example:
        >>> indices = extreme_index(cfg, ds)              # doctest: +SKIP
        >>> da = indices.calculate("tx90p")               # doctest: +SKIP
    """

    def __init__(self, cfg, climate_data):
        """Bind a configuration and a dataset.

        Args:
            cfg (DictConfig): Hydra configuration containing ``extinfo.indices``.
            climate_data (xr.Dataset): Dataset holding the required input variables.
        """
        self.cfg = cfg
        self.climate_data = climate_data

    def calculate(self, index):
        """Compute a single index by name.

        Resolution proceeds in three steps:

        1. **Linked variables** — each entry under ``index_cfg.link`` produces an
           intermediate DataArray, either by calling an importable function
           (``function_call``) or by invoking a method on one input variable
           (``operation``). An optional ``postprocess`` block may then subset the
           result. Each is stored back into :attr:`climate_data` under its key so
           later links and the index itself can reference it.
        2. **Argument resolution** — any argument written as an OmegaConf-style
           reference (``"${...}"``) is replaced by the matching entry of
           :attr:`climate_data`.
        3. **Evaluation** — ``index_cfg.function`` is imported and called with the
           variables listed in ``index_cfg.variables`` (falling back to ``pr``).

        Args:
            index (str): Index name, matching a key of ``cfg.extinfo.indices``
                (e.g. ``"tn10p"``, ``"tx90p"``, ``"rx1day"``).

        Returns:
            xr.DataArray: The computed index, renamed to ``index``. The result is
            lazy if the inputs are lazy — call ``.compute()`` to materialise it.

        Raises:
            KeyError: If ``index`` is not defined in ``cfg.extinfo.indices``, or a
                referenced input variable is absent from :attr:`climate_data`.
            ValueError: If a ``link`` entry defines neither ``function_call`` nor
                ``operation``.
            NotImplementedError: If a ``postprocess`` block names an unsupported
                operation.

        Example:
            >>> indices = extreme_index(cfg, ds)          # doctest: +SKIP
            >>> indices.calculate("tn10p").compute()      # doctest: +SKIP
        """
        index_cfg = self.cfg.extinfo.indices[index]
        args = dict(index_cfg.args)

        # Handle linked intermediate variables
        if "link" in index_cfg:
            for var_name, link_cfg in index_cfg.link.items():
                # --- External function call ---
                if "function_call" in link_cfg:

                    inputs = [self.climate_data[name] for name in link_cfg["inputs"]]
                    module_path, func_name = link_cfg["function_call"].rsplit(".", 1)
                    func = getattr(importlib.import_module(module_path), func_name)
                    result = func(*inputs, **link_cfg.get("kwargs", {}))

                # --- Method call on single input ---
                elif "operation" in link_cfg:
                    input_var = self.climate_data[link_cfg["input"]]
                    method = getattr(input_var, link_cfg["operation"])
                    result = method(**link_cfg.get("kwargs", {}))

                else:
                    raise ValueError(f"Link for '{var_name}' must define either 'function_call' or 'operation'.")

                # Optional postprocessing
                if "postprocess" in link_cfg:
                    for op_name, op_args in link_cfg["postprocess"].items():
                        if op_name == "sel":
                            result = result.sel(**op_args)
                        else:
                            raise NotImplementedError(f"Postprocess operation '{op_name}' not supported.")

                self.climate_data[var_name] = result

        # Resolve references in args
        for key, val in args.items():
            if isinstance(val, str) and val.startswith("${"):
                ref = val.strip("${}").split(".")[-1]
                args[key] = self.climate_data[ref]

        # Load and call final function
        module_name, func_name = index_cfg.function.rsplit(".", 1)
        func = getattr(importlib.import_module(module_name), func_name)

        if hasattr(index_cfg, "variables"):
            inputs = [self.climate_data[v] for v in index_cfg.variables]
            result = func(*inputs, **args)
        else:
            result = func(self.climate_data["pr"], **args)

        if isinstance(result, xr.DataArray):
            result.name = index

        return result
