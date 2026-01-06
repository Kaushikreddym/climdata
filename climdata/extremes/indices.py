import os
import importlib
from copy import deepcopy
from pathlib import Path
import xarray as xr
from xclim.core.calendar import percentile_doy

class extreme_index:
    def __init__(self, cfg, climate_data):
        self.cfg = cfg
        self.climate_data = climate_data
        self.successful_indices = []
        self.failed_indices = []
        self.tasks = []

    def zarr_output_path(self, index):
        zarr_filename = self.cfg.output.filename.format(
            index=index,
            dataset=self.cfg.dataset,
            region=self.cfg.region,
            start=self.cfg.time_range.start_date,
            end=self.cfg.time_range.end_date,
            freq='YS',
        )
        return Path('data/10km/' + zarr_filename)

    def calculate(self, index):
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

        # ✅ Rename output DataArray
        if isinstance(result, xr.DataArray):
            result.name = index

        # result = result.chunk({'time':12,'lat':-1,'lon':-1})  # trigger chunking if needed
        
        # Save result
        # zarr_path = self.zarr_output_path(index)
        # os.makedirs(zarr_path.parent, exist_ok=True)
        # result.to_zarr(str(zarr_path), mode="w")
        return result

    def run(self):
        for index in list(self.cfg.mappings.indices.keys()):
            zarr_path = self.zarr_output_path(index)
            if zarr_path.exists():
                print(f"Skipping {index}: Zarr already exists at {zarr_path}")
                continue
            try:
                print(f"Processing index: {index}")
                self.calculate(index)
                self.successful_indices.append(index)
            except Exception as e:
                print(f"Error processing {index}: {e}")
                self.failed_indices.append(index)