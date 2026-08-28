"""
ClimData Capability Discovery
=============================

Automatically discovers the capabilities exposed by the ``climdata`` package
and returns them as a structured, validated registry.

This module performs **discovery only** — it does not run any workflow, build
prompts, or call an LLM. It is the inventory layer an agent would sit on top of.

Capability sources (all discovered at runtime, nothing hard-coded that the
package already declares):

    datasets            <- climdata/conf/mappings/parameters.yaml + explore.registry
    variables           <- climdata/conf/mappings/variables.yaml
    indices             <- climdata/conf/mappings/indices.yaml   (categorised)
    bias_correction     <- climdata.sdba  (class introspection)
    imputation          <- climdata/conf/mappings/imputation.yaml + Imputer
    workflows/actions   <- climdata/conf/mappings/actions.yaml
    extraction_modes    <- ClimateExtractor.extract (+ regions.yaml)
    output_formats      <- ClimateExtractor.to_csv / to_nc / to_dataframe

Usage
-----
    from climdata.LLM.capabilities import discover_all
    reg = discover_all()
    print(reg.summary())
    reg.to_json("climdata_capabilities.json")
"""

from __future__ import annotations

import inspect
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml
from pydantic import BaseModel, Field


# ===========================================================================
# 1. PYDANTIC SCHEMA  (capability metadata model)
# ===========================================================================

class CapabilityKind(str, Enum):
    """Top-level taxonomy of a ClimData capability."""
    DATASET = "dataset"
    VARIABLE = "variable"
    INDEX = "index"               # climate / ETCCDI / bioclim index
    BIAS_CORRECTION = "bias_correction"
    REGRIDDING = "regridding"     # projection / resolution transformation
    IMPUTATION = "imputation"
    WORKFLOW = "workflow"         # high-level action (extract, calc_index, ...)
    EXTRACTION_MODE = "extraction_mode"
    OUTPUT_FORMAT = "output_format"


class IndexFamily(str, Enum):
    """Sub-classification for climate indices."""
    ETCCDI_PRECIP = "etccdi_precipitation"
    ETCCDI_TEMP = "etccdi_temperature"
    DROUGHT = "drought_spi"
    BIOCLIM = "bioclimatic"
    AGROCLIMATIC = "agroclimatic"
    OTHER = "other"


class ParamSpec(BaseModel):
    """A single parameter accepted by a capability."""
    name: str
    type: str = "Any"
    required: bool = False
    default: Optional[Any] = None
    description: Optional[str] = None


class Capability(BaseModel):
    """
    Uniform metadata record for any ClimData capability.

    The same schema describes a dataset, a variable, an index, a bias-correction
    method, etc. — so an agent can reason over a single flat registry.
    """
    name: str = Field(..., description="Canonical identifier (key).")
    kind: CapabilityKind
    description: str = ""
    required_params: List[ParamSpec] = Field(default_factory=list)
    optional_params: List[ParamSpec] = Field(default_factory=list)
    example: Optional[str] = Field(None, description="Minimal runnable example.")
    # Free-form, kind-specific extras (units, coverage, distribution, ...)
    metadata: Dict[str, Any] = Field(default_factory=dict)


class CapabilityRegistry(BaseModel):
    """The complete discovered inventory, grouped by kind."""
    version: str = ""
    datasets: List[Capability] = Field(default_factory=list)
    variables: List[Capability] = Field(default_factory=list)
    indices: List[Capability] = Field(default_factory=list)
    bias_correction: List[Capability] = Field(default_factory=list)
    regridding: List[Capability] = Field(default_factory=list)
    imputation: List[Capability] = Field(default_factory=list)
    workflows: List[Capability] = Field(default_factory=list)
    extraction_modes: List[Capability] = Field(default_factory=list)
    output_formats: List[Capability] = Field(default_factory=list)

    # -- convenience ----------------------------------------------------
    def all(self) -> List[Capability]:
        """Flatten every capability group into one list.

        Returns:
            list[Capability]: All capabilities, grouped in declaration order.
        """
        return (
            self.datasets + self.variables + self.indices + self.bias_correction
            + self.regridding + self.imputation + self.workflows
            + self.extraction_modes + self.output_formats
        )

    def summary(self) -> str:
        """Render a counts-per-group table of the registry.

        Returns:
            str: A printable table with one row per capability group and a total.
        """
        groups = {
            "datasets": self.datasets, "variables": self.variables,
            "indices": self.indices, "bias_correction": self.bias_correction,
            "regridding": self.regridding, "imputation": self.imputation, "workflows": self.workflows,
            "extraction_modes": self.extraction_modes,
            "output_formats": self.output_formats,
        }
        lines = [f"ClimData capability registry (v{self.version})", "-" * 44]
        for g, items in groups.items():
            lines.append(f"{g:<18}: {len(items):>3}")
        lines.append("-" * 44)
        lines.append(f"{'TOTAL':<18}: {len(self.all()):>3}")
        return "\n".join(lines)

    def to_json(self, path: str | Path) -> Path:
        """Write the whole registry to a JSON file.

        This is what an LLM is handed as its inventory of what climdata can do.

        Args:
            path (str | Path): Destination path. The parent directory must exist.

        Returns:
            Path: The path written.
        """
        path = Path(path)
        path.write_text(self.model_dump_json(indent=2))
        return path


# ===========================================================================
# 2. DISCOVERY HELPERS
# ===========================================================================

def _conf_dir() -> Path:
    """Locate climdata/conf, regardless of where this script runs from."""
    import climdata
    return Path(climdata.__file__).parent / "conf"


def _load_yaml(rel: str) -> dict:
    f = _conf_dir() / rel
    if not f.exists():
        return {}
    with open(f) as fh:
        return yaml.safe_load(fh) or {}


def _params_from_signature(func) -> tuple[List[ParamSpec], List[ParamSpec]]:
    """Split a callable's signature into required / optional ParamSpecs."""
    required, optional = [], []
    for p in inspect.signature(func).parameters.values():
        if p.name in ("self", "cls", "args", "kwargs") or p.kind in (
            p.VAR_POSITIONAL, p.VAR_KEYWORD
        ):
            continue
        ann = p.annotation
        type_name = getattr(ann, "__name__", str(ann)) if ann is not p.empty else "Any"
        if p.default is p.empty:
            required.append(ParamSpec(name=p.name, type=type_name, required=True))
        else:
            optional.append(ParamSpec(
                name=p.name, type=type_name, required=False, default=p.default
            ))
    return required, optional


# ===========================================================================
# 3. PER-KIND DISCOVERY
# ===========================================================================

def discover_datasets() -> List[Capability]:
    """
    Datasets come from climdata.explore.registry, which auto-discovers every
    provider in parameters.yaml under its UPPER-CASED key and appends static
    code-registered datasets (ERA5).
    """
    try:
        from climdata.explore.registry import get_registry
        registry = get_registry()
    except Exception:
        registry = {}

    caps: List[Capability] = []
    for key, meta in registry.items():
        caps.append(Capability(
            name=key,
            kind=CapabilityKind.DATASET,
            description=meta.get("full_name", key),
            required_params=[
                ParamSpec(name="dataset", type="str", required=True,
                          default=key, description="Provider key for ClimData(overrides=...)"),
            ],
            optional_params=[
                ParamSpec(name="variables", type="list[str]", default=meta.get("variables", [])),
                ParamSpec(name="time_range", type="dict"),
                ParamSpec(name="experiment_id", type="str"),
            ],
            example=f'ClimData(overrides=["dataset={key}", "lat=52", "lon=13"]).extract()',
            metadata={
                "type": meta.get("type"),
                "coverage": meta.get("coverage"),
                "resolution": meta.get("resolution"),
                "frequency": meta.get("frequency"),
                "time_range": meta.get("time_range"),
                "source": meta.get("source"),
                "variables": meta.get("variables", []),
                "experiments": meta.get("experiments", []),
                "notes": meta.get("notes", ""),
            },
        ))
    return caps


def discover_variables() -> List[Capability]:
    """Canonical CF variables from variables.yaml."""
    var_map = _load_yaml("mappings/variables.yaml")
    caps: List[Capability] = []
    for name, meta in var_map.items():
        caps.append(Capability(
            name=name,
            kind=CapabilityKind.VARIABLE,
            description=meta.get("long_name", name),
            example=f'ClimData(overrides=["variables=[{name}]"])',
            metadata={
                "cf_name": meta.get("cf_name"),
                "units": meta.get("units"),
            },
        ))
    return caps


# ETCCDI / bioclim classification by xclim function or index name.
_ETCCDI_PRECIP = {
    "Rx1day", "Rx3day", "max_cdd", "max_cwd", "prcptot", "r75p", "f75p",
    "wetdays", "wetdays_prop", "daily_pr_intensity", "dry_days",
    "dry_spell_frequency", "dry_spell_max_length", "dry_spell_total_length",
    "wet_spell_frequency", "wet_spell_max_length", "wet_spell_total_length",
}
_ETCCDI_TEMP = {
    "tg10p", "tg90p", "tn10p", "ice_days", "heat_wave_index", "heat_wave_frequency",
    "heat_wave_max_length", "heat_wave_total_length", "hot_spell_frequency",
    "hot_spell_max_length", "hot_spell_total_length", "hot_spell_max_magnitude",
    "heating_degree_days", "maximum_consecutive_frost_days",
    "maximum_consecutive_frost_free_days", "maximum_consecutive_tx_days",
    "tg_days_above_10degC", "tg_days_below_10degC",
}
_DROUGHT = {"SPI1", "SPI3", "SPI6", "keetch_byram_drought_index"}
_BIOCLIM = {
    "isothermality", "temperature_seasonality", "precip_seasonality",
    "warmest_quarter", "coldest_quarter",
}
_AGROCLIM = {"huglin_index"}


def _index_family(name: str) -> IndexFamily:
    if name in _ETCCDI_PRECIP:
        return IndexFamily.ETCCDI_PRECIP
    if name in _ETCCDI_TEMP:
        return IndexFamily.ETCCDI_TEMP
    if name in _DROUGHT:
        return IndexFamily.DROUGHT
    if name in _BIOCLIM:
        return IndexFamily.BIOCLIM
    if name in _AGROCLIM:
        return IndexFamily.AGROCLIMATIC
    return IndexFamily.OTHER


def discover_indices() -> List[Capability]:
    """
    Climate indices from indices.yaml. ETCCDI vs bioclim vs drought is
    encoded in each Capability's metadata['family'] so they are queryable
    as one registry but separable on demand.
    """
    indices_cfg = _load_yaml("mappings/indices.yaml").get("indices", {})
    caps: List[Capability] = []
    for name, cfg in indices_cfg.items():
        args = cfg.get("args", {}) or {}
        req_vars = cfg.get("variables", [])
        caps.append(Capability(
            name=name,
            kind=CapabilityKind.INDEX,
            description=f"{cfg.get('function', '')} on {req_vars}",
            required_params=[
                ParamSpec(name="index", type="str", required=True, default=name),
                ParamSpec(name="variables", type="list[str]", required=True,
                          default=req_vars, description="Input variables this index needs"),
            ],
            optional_params=[
                ParamSpec(name=k, type=type(v).__name__, default=v)
                for k, v in args.items()
            ],
            example=f'ClimData(overrides=["dataset=MSWX", "index={name}"]).calc_index()',
            metadata={
                "family": _index_family(name).value,
                "function": cfg.get("function"),
                "input_variables": req_vars,
                "default_args": args,
                "has_linked_inputs": "link" in cfg,
            },
        ))
    return caps


def discover_bias_correction() -> List[Capability]:
    """Bias-correction / downscaling capabilities by introspecting climdata.sdba."""
    caps: List[Capability] = []
    try:
        from climdata import sdba
    except Exception:
        return caps

    class_examples = {
        "BiasCorrection": "BiasCorrection(variable='tas').correct(obs_hist, sim_hist, sim_fut)",
        "StatisticalDownscaling": "StatisticalDownscaling(variable='tas').downscale(obs_fine, sim_coarse)",
        "BCSD": "BCSD(variable='tas').run(obs_fine, sim_hist_coarse, sim_fut_coarse)",
    }
    for cls_name in ("BiasCorrection", "StatisticalDownscaling", "BCSD"):
        cls = getattr(sdba, cls_name, None)
        if cls is None:
            continue
        req, opt = _params_from_signature(cls.__init__)
        # Per-variable defaults (distribution / trend preservation) are valuable metadata
        default_configs = getattr(cls, "DEFAULT_CONFIGS", {})
        caps.append(Capability(
            name=cls_name,
            kind=CapabilityKind.BIAS_CORRECTION,
            description=(inspect.getdoc(cls) or "").split("\n\n")[0].strip(),
            required_params=req,
            optional_params=opt,
            example=class_examples.get(cls_name),
            metadata={
                "method": "ISIMIP3BASD (trend-preserving quantile mapping / MBCn)",
                "supported_variables": sorted(default_configs.keys()),
                "per_variable_defaults": default_configs,
            },
        ))

    # Functional capability: regridding
    if hasattr(sdba, "regrid_to_coarse"):
        req, opt = _params_from_signature(sdba.regrid_to_coarse)
        caps.append(Capability(
            name="regrid_to_coarse",
            kind=CapabilityKind.BIAS_CORRECTION,
            description="Regrid fine-resolution data to a coarse target grid (xESMF or CDO).",
            required_params=req,
            optional_params=opt,
            example="regrid_to_coarse(fine_data, coarse_template, method='conservative')",
            metadata={"tools": ["xesmf", "cdo"]},
        ))
    return caps


def discover_regridding() -> List[Capability]:
    """Projection and resolution transformation exposed by :mod:`climdata.grid`."""
    return [
        Capability(
            name="reproject",
            kind=CapabilityKind.REGRIDDING,
            description=(
                "Reproject and/or resample a dataset onto a target CRS and resolution "
                "(rasterio/rioxarray), with the grid origin snapped so that separately "
                "processed domains stay co-registered."
            ),
            required_params=[],
            optional_params=[
                ParamSpec(
                    name="target_projection", type="str", required=False,
                    description=(
                        "Destination CRS: 'EPSG:4326', an integer EPSG code, OGC WKT or "
                        "a PROJ string. Defaults to the source CRS."
                    ),
                ),
                ParamSpec(
                    name="target_resolution", type="str", required=False,
                    description=(
                        "Destination cell size. Units MUST match the target CRS: linear "
                        "('10 km', '1000 m') for a projected CRS, angular ('0.1 deg', "
                        "'30 arcsec') for a geographic one. Combining a metric "
                        "resolution with a geographic CRS such as EPSG:4326 raises "
                        "ResolutionCRSMismatch, because 10 km spans 0.1395 deg of "
                        "longitude but 0.0899 deg of latitude at 50N."
                    ),
                ),
                ParamSpec(
                    name="regrid.method", type="str", required=False,
                    description=(
                        "Resampling: nearest | bilinear | cubic | cubic_spline | lanczos "
                        "| average | mode | min | max | med | q1 | q3 | sum | rms. Use "
                        "'bilinear' for state variables, 'average' for fluxes, "
                        "'nearest' for categorical fields."
                    ),
                ),
                ParamSpec(
                    name="regrid.align", type="bool", required=False,
                    description="Snap the grid origin to whole multiples of the resolution.",
                ),
                ParamSpec(
                    name="regrid.engine", type="str", required=False,
                    description="'rasterio' (default) or 'xesmf' for area-weighted conservative remapping.",
                ),
            ],
            example=(
                'ClimData(overrides=["target_projection=EPSG:3035", '
                '"target_resolution=10 km"])'
            ),
            metadata={
                "tools": ["rasterio", "rioxarray", "pyproj", "xclim.core.units"],
                "valid_combinations": {
                    "projected_crs": "linear resolution (10 km, 1000 m, 30 ft)",
                    "geographic_crs": "angular resolution (0.1 deg, 30 arcsec, 5 arcmin)",
                },
                "raises": "ResolutionCRSMismatch on linear+geographic or angular+projected",
            },
        ),
    ]


def discover_imputation() -> List[Capability]:
    """
    Imputation methods declared in imputation.yaml, cross-checked against the
    methods actually wired into climdata.impute.Imputer.impute().
    """
    impute_cfg = _load_yaml("mappings/imputation.yaml")
    # Methods the Imputer.impute() branch actually dispatches on:
    wired = {"BRITS", "SoftImpute", "CDRec", "XGBOOST"}
    torch_methods = set()
    try:
        from climdata.impute.impute_xarray import Imputer
        torch_methods = getattr(Imputer, "_TORCH_METHODS", set())
    except Exception:
        pass

    caps: List[Capability] = []
    for name, meta in impute_cfg.items():
        opt = [ParamSpec(name="normalize", type="bool", default=True)]
        if name == "BRITS":
            opt.append(ParamSpec(name="epochs", type="int",
                                 default=meta.get("epochs", 300)))
        caps.append(Capability(
            name=name,
            kind=CapabilityKind.IMPUTATION,
            description=meta.get("name", name),
            required_params=[
                ParamSpec(name="impute", type="str", required=True, default=name),
            ],
            optional_params=opt,
            example=f'ClimData(overrides=["impute={name}"]).impute()',
            metadata={
                "wired_in_impute": name in wired,
                "requires_pytorch": name in torch_methods,
            },
        ))
    return caps


def discover_workflows() -> List[Capability]:
    """High-level actions from actions.yaml."""
    actions = _load_yaml("mappings/actions.yaml")
    caps: List[Capability] = []
    for name, meta in actions.items():
        caps.append(Capability(
            name=name,
            kind=CapabilityKind.WORKFLOW,
            description=meta.get("description", ""),
            example=f"extractor.{name}(...)",
            metadata={
                "action_type": meta.get("type"),     # primary | secondary
                "output": meta.get("output"),        # dataset | filename
            },
        ))
    return caps


def discover_extraction_modes() -> List[Capability]:
    """
    Spatial-subset modes accepted by ClimateExtractor.extract(), inferred from
    the config fields it reads (point / box / region / shapefile / aoi).
    """
    regions = list(_load_yaml("regions/regions.yaml").keys())
    modes = [
        Capability(
            name="point",
            kind=CapabilityKind.EXTRACTION_MODE,
            description="Nearest grid cell to a lat/lon point (DWD adds a 30 km buffer).",
            required_params=[
                ParamSpec(name="lat", type="float", required=True),
                ParamSpec(name="lon", type="float", required=True),
            ],
            example='ClimData(overrides=["lat=52", "lon=13"])',
        ),
        Capability(
            name="box",
            kind=CapabilityKind.EXTRACTION_MODE,
            description="Rectangular bounding box subset.",
            required_params=[
                ParamSpec(name="box.lat_min", type="float", required=True),
                ParamSpec(name="box.lat_max", type="float", required=True),
                ParamSpec(name="box.lon_min", type="float", required=True),
                ParamSpec(name="box.lon_max", type="float", required=True),
            ],
            example='ClimData(overrides=["box.lat_min=47", "box.lat_max=55", "box.lon_min=5", "box.lon_max=15"])',
        ),
        Capability(
            name="region",
            kind=CapabilityKind.EXTRACTION_MODE,
            description="Named region whose bounds come from regions.yaml.",
            required_params=[
                ParamSpec(name="region", type="str", required=True,
                          description=f"One of: {regions}"),
            ],
            example='ClimData(overrides=["region=germany"])',
            metadata={"available_regions": regions},
        ),
        Capability(
            name="shapefile",
            kind=CapabilityKind.EXTRACTION_MODE,
            description="Clip to the geometry of a shapefile / GeoJSON (rioxarray).",
            required_params=[ParamSpec(name="shapefile", type="str", required=True)],
            example='ClimData(overrides=["shapefile=/path/to/aoi.shp"])',
        ),
    ]
    return modes


def discover_output_formats() -> List[Capability]:
    """
    Output formats supported by the writer methods on ClimateExtractor.
    CSV sub-formats are read from the to_csv docstring/branches.
    """
    return [
        Capability(
            name="dataframe",
            kind=CapabilityKind.OUTPUT_FORMAT,
            description="Long-form pandas DataFrame (time, lat, lon, variable, value, units, source).",
            example="extractor.to_dataframe(ds)",
            metadata={"method": "to_dataframe"},
        ),
        Capability(
            name="csv:default",
            kind=CapabilityKind.OUTPUT_FORMAT,
            description="Tab-separated long-form CSV.",
            example='extractor.to_csv(df, format="default")',
            metadata={"method": "to_csv", "format": "default"},
        ),
        Capability(
            name="csv:simplace",
            kind=CapabilityKind.OUTPUT_FORMAT,
            description="SIMPLACE crop-model format (per-location, tab-separated).",
            example='extractor.to_csv(df, format="simplace")',
            metadata={"method": "to_csv", "format": "simplace"},
        ),
        Capability(
            name="csv:monica",
            kind=CapabilityKind.OUTPUT_FORMAT,
            description="MONICA crop-model format (per-location, tab-separated).",
            example='extractor.to_csv(df, format="monica")',
            metadata={"method": "to_csv", "format": "monica"},
        ),
        Capability(
            name="netcdf",
            kind=CapabilityKind.OUTPUT_FORMAT,
            description="CF-style NetCDF file.",
            example="extractor.to_nc(ds)",
            metadata={"method": "to_nc"},
        ),
        Capability(
            name="zarr",
            kind=CapabilityKind.OUTPUT_FORMAT,
            description="Zarr store (filename template configured in output.filename_zarr).",
            example="ds.to_zarr(extractor.filename_zarr)",
            metadata={"method": "to_zarr"},
        ),
    ]


# ===========================================================================
# 4. TOP-LEVEL ENTRY POINT
# ===========================================================================

def discover_all() -> CapabilityRegistry:
    """Discover every ClimData capability and return a validated registry."""
    try:
        import climdata
        version = getattr(climdata, "__version__", "")
    except Exception:
        version = ""

    return CapabilityRegistry(
        version=version,
        datasets=discover_datasets(),
        variables=discover_variables(),
        indices=discover_indices(),
        bias_correction=discover_bias_correction(),
        regridding=discover_regridding(),
        imputation=discover_imputation(),
        workflows=discover_workflows(),
        extraction_modes=discover_extraction_modes(),
        output_formats=discover_output_formats(),
    )


if __name__ == "__main__":
    reg = discover_all()
    print(reg.summary())
    out = reg.to_json(Path(__file__).parent / "climdata_capabilities.json")
    print(f"\nWrote full registry to: {out}")
