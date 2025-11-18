"""Top-level package for climdata."""

__author__ = """Kaushik Muduchuru"""
__email__ = "kaushik.reddy.m@gmail.com"
__version__ = "0.2.2"

from .utils.utils_download import * # etc.
from .utils.config import load_config
from .utils.wrapper import extract_data
from .datasets.DWD import DWDmirror as DWD
from .datasets.MSWX import MSWXmirror as MSWX
from .datasets.ERA5 import ERA5Mirror as ERA5
from .datasets.CMIPlocal import CMIPmirror as CMIPlocal
from .datasets.CMIPCloud import CMIPCloud as CMIP
from .datasets.HYRAS import HYRASmirror as HYRAS

