"""Path helpers for resources, compatible with PyInstaller frozen builds."""

from __future__ import annotations

import os
import sys
from pathlib import Path


def _base_dir() -> Path:
    """Return the base directory for bundled resources.

    When frozen by PyInstaller, ``sys._MEIPASS`` points at the temporary
    extraction directory. Otherwise we resolve relative to this file so that
    ``resource_path("resources/html/map.html")`` works from anywhere.
    """
    meipass = getattr(sys, "_MEIPASS", None)
    if meipass:
        return Path(meipass)
    return Path(__file__).resolve().parent.parent


def resource_path(path: str | os.PathLike) -> str:
    """Resolve a resource path relative to the GUI package root.

    Args:
        path: Relative path under the ``climdata_gui`` package
            (e.g. ``"resources/html/map.html"``).

    Returns:
        An absolute filesystem path as a string.
    """
    return str((_base_dir() / Path(path)).resolve())


def project_root() -> Path:
    """Return the repository root (one level above the GUI package)."""
    return _base_dir().parent
