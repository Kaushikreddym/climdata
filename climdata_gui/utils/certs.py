"""TLS trust-store setup for frozen builds.

A PyInstaller bundle carries its own OpenSSL, and that library still looks for
CA certificates at the path it was *compiled* with — a path that belongs to the
build machine's environment and generally does not exist on the machine running
the .app. Every HTTPS call then fails with::

    SSLCertVerificationError: unable to get local issuer certificate

The bundle already ships ``certifi``'s CA file, so the fix is to point OpenSSL
at it explicitly before anything opens a socket. ``SSL_CERT_FILE`` is read by
``ssl.create_default_context()`` (and therefore by aiohttp, fsspec, urllib) at
the moment the context is built, so setting it early in ``main()`` covers every
library the pipeline uses.
"""

from __future__ import annotations

import os
import ssl
import sys
from pathlib import Path
from typing import List, Optional

# Set by configure_ca_bundle(); surfaced in the GUI log for diagnosis.
_ca_bundle: Optional[str] = None
_ca_status: str = "not configured"


def _candidates() -> List[Path]:
    """Possible locations of a CA bundle, best first."""
    paths: List[Path] = []

    try:
        import certifi
        paths.append(Path(certifi.where()))
    except Exception:      # noqa: BLE001 — certifi is optional at runtime
        pass

    meipass = getattr(sys, "_MEIPASS", None)
    if meipass:
        paths.append(Path(meipass) / "certifi" / "cacert.pem")

    if getattr(sys, "frozen", False):
        # macOS .app layout: MacOS/<exe>, Resources/certifi/cacert.pem
        contents = Path(sys.executable).resolve().parent.parent
        paths.append(contents / "Resources" / "certifi" / "cacert.pem")
        paths.append(contents / "Frameworks" / "certifi" / "cacert.pem")

    return paths


def _default_store_is_usable() -> bool:
    """Whether OpenSSL's built-in CA location exists on this machine."""
    paths = ssl.get_default_verify_paths()
    if paths.cafile and os.path.isfile(paths.cafile):
        return True
    return bool(paths.capath and os.path.isdir(paths.capath)
                and os.listdir(paths.capath))


def configure_ca_bundle() -> Optional[str]:
    """Make sure HTTPS verification has a CA bundle it can actually read.

    Call once, as early as possible in ``main()`` — before any network library
    builds an SSL context.

    Returns:
        The CA bundle path now in force, or ``None`` when the platform default
        is already usable (or nothing could be found).
    """
    global _ca_bundle, _ca_status

    existing = os.environ.get("SSL_CERT_FILE")
    if existing and os.path.isfile(existing):
        _ca_bundle, _ca_status = existing, f"using SSL_CERT_FILE from the environment: {existing}"
        return existing

    frozen = getattr(sys, "frozen", False)
    # Running from source with a working system store: leave well alone.
    if not frozen and _default_store_is_usable():
        _ca_status = "using the interpreter's default CA store"
        return None

    for candidate in _candidates():
        try:
            if candidate.is_file():
                path = str(candidate)
                os.environ["SSL_CERT_FILE"] = path
                os.environ["REQUESTS_CA_BUNDLE"] = path
                # A stale SSL_CERT_DIR would otherwise still be consulted.
                if not os.path.isdir(os.environ.get("SSL_CERT_DIR", "")):
                    os.environ.pop("SSL_CERT_DIR", None)
                _ca_bundle, _ca_status = path, f"CA bundle: {path}"
                return path
        except OSError:
            continue

    _ca_status = ("no CA bundle found — HTTPS downloads may fail with "
                  "'certificate verify failed'")
    return None


def ca_bundle_status() -> str:
    """One-line description of the trust store in force, for the GUI log."""
    return _ca_status
