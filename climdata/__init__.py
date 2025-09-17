"""Top-level package for climdata."""

__author__ = """Kaushik Muduchuru"""
__email__ = "kaushik.reddy.m@gmail.com"
__version__ = "0.0.2"

from .utils.utils_download import * # etc.
from .datasets.DWD import DWDmirror as DWD
