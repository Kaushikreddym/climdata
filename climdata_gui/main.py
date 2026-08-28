"""Entry point for climdata_gui.

Run from the ``climdata_gui`` directory:

    python main.py
"""

from __future__ import annotations

import os
import sys


def _ensure_on_path() -> None:
    """Make sibling packages (gui/, backend/, ...) importable when launched
    via ``python main.py`` from inside ``climdata_gui``."""
    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.insert(0, here)


def main() -> int:
    _ensure_on_path()

    # Before anything can open an HTTPS connection: a frozen build's OpenSSL
    # looks for CA certificates where it was built, not where it now runs.
    from utils.certs import configure_ca_bundle
    configure_ca_bundle()

    from app import create_app
    from gui.main_window import MainWindow

    app = create_app(sys.argv)
    window = MainWindow()
    window.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
