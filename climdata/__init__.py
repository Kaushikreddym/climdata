"""climdata — one interface for many climate data providers.

Extract climate data from reanalyses, observational grids, station networks and
CMIP6 projections through a single configuration-driven API, then compute
extreme indices, regrid, bias-correct, impute gaps and export the result.

The entry point for most work is :class:`ClimData`, which drives an entire
workflow from a Hydra configuration:

    >>> import climdata as cd
    >>> extractor = cd.ClimData(overrides=[        # doctest: +SKIP
    ...     "dataset=mswx", "lat=52.5", "lon=13.4",
    ...     "variables=[tasmin,tasmax,pr]",
    ...     "time_range.start_date=2014-01-01",
    ...     "time_range.end_date=2014-12-31",
    ... ])
    >>> ds = extractor.extract()                   # doctest: +SKIP
    >>> df = extractor.to_dataframe()              # doctest: +SKIP

To see what is available before extracting anything:

    >>> cd.list_available_data()                   # doctest: +SKIP
    >>> cd.explore("NEXGDDP")                      # doctest: +SKIP
    >>> cd.find(variable="pr", frequency="daily")  # doctest: +SKIP

The provider classes (:class:`MSWX`, :class:`CMIP`, :class:`HYRAS`, …) are also
exported for direct use when the workflow wrapper is more than you need.
"""

__author__ = """Kaushik Muduchuru"""
__email__ = "kaushik.reddy.m@gmail.com"
__version__ = "1.0.0"

from .utils.config import load_config
from .datasets.DWD import DWDmirror as DWD
from .datasets.MSWX import MSWXmirror as MSWX
from .datasets.ERA5 import ERA5Mirror as ERA5
from .datasets.CMIPCloud import CMIPCloud as CMIP
from .datasets.W5E5 import W5E5 as W5E5
from .datasets.CMIP_W5E5 import CMIPW5E5 as CMIPW5E5
from .datasets.NEXGDDP import NEXGDDP as NEXGDDP
from .datasets.HYRAS import HYRASmirror as HYRAS
from .datasets.HOSTRADA import HOSTRADAmirror as HOSTRADA
from .datasets.NASAPOWER import POWER as POWER
from .datasets.AGRI_ISIMIP import AGRI_ISIMIP as AGRI_ISIMIP
from .extremes.indices import extreme_index as extreme_index
from .utils.wrapper_workflow import ClimateExtractor as ClimData
from .viz import plot_map as plot_map
from .grid import (
    reproject,
    parse_crs,
    parse_resolution,
    to_angular,
    Resolution,
    ResolutionCRSMismatch,
)
from ._vendor import imputegap  # noqa: F401 — re-exported vendored dependency
from .explore import (
    list_available_data,
    explore,
    find,
    inspect,
    list_esm_experiments,
    list_esm_models,
    DatasetRegistry,
)

__all__ = [
    "__version__",
    # workflow
    "ClimData",
    "load_config",
    # providers
    "AGRI_ISIMIP",
    "CMIP",
    "CMIPW5E5",
    "DWD",
    "ERA5",
    "HOSTRADA",
    "HYRAS",
    "MSWX",
    "NEXGDDP",
    "POWER",
    "W5E5",
    # analysis
    "extreme_index",
    "plot_map",
    # grid
    "reproject",
    "parse_crs",
    "parse_resolution",
    "to_angular",
    "Resolution",
    "ResolutionCRSMismatch",
    # discovery
    "list_available_data",
    "explore",
    "find",
    "inspect",
    "list_esm_experiments",
    "list_esm_models",
    "DatasetRegistry",
]
