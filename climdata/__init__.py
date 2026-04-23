"""Top-level package for climdata."""

__author__ = """Kaushik Muduchuru"""
__email__ = "kaushik.reddy.m@gmail.com"
__version__ = "0.6.0"

from .utils.config import load_config
from .datasets.DWD import DWDmirror as DWD
from .datasets.MSWX import MSWXmirror as MSWX
from .datasets.ERA5 import ERA5Mirror as ERA5
# from .datasets.CMIPlocal import CMIPmirror as CMIPlocal
from .datasets.CMIPCloud import CMIPCloud as CMIP
from .datasets.W5E5 import W5E5 as W5E5
from .datasets.CMIP_W5E5 import CMIPW5E5 as CMIPW5E5
from .datasets.NEXGDDP import NEXGDDP as NEXGDDP
from .datasets.HYRAS import HYRASmirror as HYRAS
from .datasets.NASAPOWER import POWER as POWER
from .extremes.indices import extreme_index as extreme_index
from .utils.wrapper_workflow import ClimateExtractor as ClimData
from ._vendor import imputegap
from .explore import (
    list_available_data,
    explore,
    find,
    inspect,
    list_esm_experiments,
    list_esm_models,
    DatasetRegistry,
)
# from .impute.impute_xarray import Imputer as imputer_xarray

