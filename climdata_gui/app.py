"""QApplication factory — centralizes app-wide setup (styles, org info)."""

from __future__ import annotations

import sys
from typing import Optional

from PySide6.QtWidgets import QApplication

from utils.certs import configure_ca_bundle
from utils.paths import resource_path


def create_app(argv: Optional[list] = None) -> QApplication:
    """Return a singleton ``QApplication`` configured for climdata_gui."""
    # Idempotent safety net for entry points that bypass main.py.
    configure_ca_bundle()

    app = QApplication.instance()
    if app is None:
        app = QApplication(argv if argv is not None else sys.argv)

    app.setApplicationName("climdata")
    app.setOrganizationName("climdata")

    try:
        with open(resource_path("gui/styles/theme.qss"), "r", encoding="utf-8") as fh:
            app.setStyleSheet(fh.read())
    except OSError:
        # Styling is optional — fall back to defaults if the QSS is missing.
        pass

    return app
