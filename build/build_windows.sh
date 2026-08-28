#!/usr/bin/env bash
# build_windows.sh — Build climdata_gui.exe using PyInstaller in a conda env
#
# MUST be run on Windows: PyInstaller does not cross-compile, so a Windows
# .exe can only be produced by a Windows Python interpreter. Use Git Bash,
# MSYS2 or the "Anaconda Prompt (bash)" shell.
#
# Usage (from anywhere in the repo):
#   bash build/build_windows.sh                  # uses the climdata_dev env
#   bash build/build_windows.sh --env sdba       # pick the conda env
#   CLIMDATA_ENV=sdba bash build/build_windows.sh
#   bash build/build_windows.sh --check          # pre-flight only, no build
#   bash build/build_windows.sh --zip            # also produce a .zip to ship
#
# Output: dist/climdata_gui/climdata_gui.exe

set -euo pipefail

# ── 0. Always operate from the repository root ────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

CONDA_ENV="${CLIMDATA_ENV:-climdata_dev}"
SPEC_FILE="climdata_gui_win.spec"
WORK_DIR="build_work"
DIST_DIR="dist"
APP_DIR="${DIST_DIR}/climdata_gui"
EXE_PATH="${APP_DIR}/climdata_gui.exe"

CHECK_ONLY=0
MAKE_ZIP=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --env)     CONDA_ENV="${2:?--env needs an environment name}"; shift 2 ;;
        --env=*)   CONDA_ENV="${1#*=}"; shift ;;
        --check)   CHECK_ONLY=1; shift ;;
        --zip)     MAKE_ZIP=1; shift ;;
        -h|--help) # print the header comment block, stopping at the first code line
                   sed -n '2,${/^#/!q;s/^# \{0,1\}//;p;}' "${BASH_SOURCE[0]}"; exit 0 ;;
        *)         echo "ERROR: unknown option '$1' (try --help)"; exit 2 ;;
    esac
done

# ── 1. Verify we are on Windows ───────────────────────────────────────────
# PyInstaller freezes for the host platform only; running this on Linux/macOS
# would silently produce a non-Windows binary.
case "$(uname -s)" in
    MINGW*|MSYS*|CYGWIN*|Windows_NT) IS_WINDOWS=1 ;;
    *)                               IS_WINDOWS=0 ;;
esac

if [[ ${IS_WINDOWS} -eq 0 ]]; then
    echo "ERROR: this script must run on Windows (Git Bash / MSYS2 / Cygwin)."
    echo "       Detected host: $(uname -s). PyInstaller cannot cross-compile;"
    echo "       build the .app with build/build_mac.sh on macOS instead."
    if [[ ${CHECK_ONLY} -eq 0 ]]; then
        exit 1
    fi
    echo "       (--check given: continuing with the platform-independent checks.)"
fi

# ── 2. Verify the conda env exists ────────────────────────────────────────
if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: 'conda' not on PATH. Open an Anaconda/Miniforge shell, or run"
    echo "       'conda init bash' once and restart Git Bash."
    exit 1
fi

if ! conda env list | grep -qE "^${CONDA_ENV}[[:space:]]"; then
    echo "ERROR: conda environment '${CONDA_ENV}' not found."
    echo "       Create it first:  conda env create -n ${CONDA_ENV} ..."
    exit 1
fi

echo "==> Using conda environment: ${CONDA_ENV}"

# ── 3. Verify the runtime dependencies are importable in that env ─────────
# A missing package here surfaces as a cryptic PyInstaller traceback later,
# so fail fast with a readable list instead.
echo "==> Checking GUI dependencies in ${CONDA_ENV}..."
conda run --no-capture-output -n "${CONDA_ENV}" python - <<'PY'
import importlib.util
import sys

required = [
    "PySide6", "hydra", "omegaconf", "matplotlib", "cartopy", "pyproj",
    "xclim", "climdata", "googleapiclient", "intake", "intake_esm",
    "wetterdienst", "cdsapi", "rasterio", "cf_xarray", "yamale",
]
missing = [m for m in required if importlib.util.find_spec(m) is None]
for m in required:
    print(f"    {m:<18} {'ok' if m not in missing else 'MISSING'}")
if missing:
    print("\nERROR: missing packages: " + ", ".join(missing), file=sys.stderr)
    sys.exit(1)
PY

# ── 4. Verify the spec and entry point are present ────────────────────────
for path in "${SPEC_FILE}" climdata_gui/main.py climdata_gui/resources climdata_gui/gui/styles climdata/conf; do
    if [[ ! -e "${path}" ]]; then
        echo "ERROR: expected '${path}' under ${REPO_ROOT} but it is missing."
        exit 1
    fi
done
echo "==> Spec, entry point and resources found."

if [[ ${CHECK_ONLY} -eq 1 ]]; then
    echo ""
    echo "✓ Pre-flight checks passed. Re-run without --check to build."
    exit 0
fi

# ── 5. Install / upgrade PyInstaller inside the env ───────────────────────
echo "==> Installing PyInstaller into ${CONDA_ENV}..."
conda run --no-capture-output -n "${CONDA_ENV}" \
    python -m pip install --quiet --upgrade pyinstaller

# ── 6. Sync conf into climdata_gui/ ───────────────────────────────────────
# The spec bundles climdata_gui/conf; keep it in step with the library conf.
echo "==> Syncing climdata/conf → climdata_gui/conf..."
rm -rf climdata_gui/conf
cp -r climdata/conf climdata_gui/conf

# ── 7. Build ──────────────────────────────────────────────────────────────
echo "==> Running PyInstaller (this takes several minutes)..."
rm -rf "${APP_DIR}"
conda run --no-capture-output -n "${CONDA_ENV}" pyinstaller "${SPEC_FILE}" \
    --noconfirm \
    --workpath "${WORK_DIR}" \
    --distpath "${DIST_DIR}"

# ── 8. Sanity-check output ────────────────────────────────────────────────
if [[ ! -f "${EXE_PATH}" ]]; then
    echo "ERROR: expected ${EXE_PATH} was not created. Check the PyInstaller output above."
    exit 1
fi

# QtWebEngine renders the Leaflet map; without its helper process the map
# panel comes up blank with no error, so check for it explicitly.
if ! find "${APP_DIR}" -name 'QtWebEngineProcess.exe' -print -quit | grep -q .; then
    echo "WARNING: QtWebEngineProcess.exe not found in the bundle."
    echo "         The map panel will not render. Verify PySide6 ships"
    echo "         QtWebEngine in the '${CONDA_ENV}' environment."
fi

SIZE=$(du -sh "${APP_DIR}" | cut -f1)
echo ""
echo "✓ Build succeeded!"
echo "  Output : ${APP_DIR}  (${SIZE})"
echo "  Exe    : ${EXE_PATH}"

# ── 9. Optional distributable zip ─────────────────────────────────────────
if [[ ${MAKE_ZIP} -eq 1 ]]; then
    ZIP_PATH="${DIST_DIR}/climdata_gui-windows-x64.zip"
    echo ""
    echo "==> Packaging ${ZIP_PATH}..."
    rm -f "${ZIP_PATH}"
    if command -v powershell.exe >/dev/null 2>&1; then
        # PowerShell needs native Windows paths, not the MSYS-style ones bash uses.
        WIN_SRC="$(cygpath -w "${APP_DIR}" 2>/dev/null || echo "${APP_DIR}")"
        WIN_ZIP="$(cygpath -w "${REPO_ROOT}/${ZIP_PATH}" 2>/dev/null || echo "${ZIP_PATH}")"
        powershell.exe -NoProfile -Command \
            "Compress-Archive -Path '${WIN_SRC}\\*' -DestinationPath '${WIN_ZIP}' -Force"
    elif command -v zip >/dev/null 2>&1; then
        ( cd "${DIST_DIR}" && zip -qr "$(basename "${ZIP_PATH}")" climdata_gui )
    else
        echo "WARNING: neither powershell.exe nor zip available; skipping the archive."
        ZIP_PATH=""
    fi
    [[ -n "${ZIP_PATH}" && -f "${ZIP_PATH}" ]] && echo "  Zip    : ${ZIP_PATH}"
fi

echo ""
echo "  To launch:"
echo "    ./${EXE_PATH}"
echo ""
echo "  NOTE: The exe is unsigned. SmartScreen will show"
echo "  'Windows protected your PC' on first run — choose"
echo "  More info → Run anyway."
