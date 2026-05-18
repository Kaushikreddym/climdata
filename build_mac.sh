#!/usr/bin/env bash
# build_mac.sh — Build climdata_gui.app using PyInstaller in the climdata_dev conda env
#
# Usage (from workspace root):
#   bash build_mac.sh
#
# Output: dist/climdata_gui.app

set -euo pipefail

CONDA_ENV="climdata_dev"
SPEC_FILE="climdata_gui.spec"

# ── 1. Verify the conda env exists ────────────────────────────────────────
if ! conda env list | grep -q "^${CONDA_ENV}[[:space:]]"; then
    echo "ERROR: conda environment '${CONDA_ENV}' not found."
    echo "       Create it first:  conda env create -n ${CONDA_ENV} ..."
    exit 1
fi

echo "==> Using conda environment: ${CONDA_ENV}"

# ── 2. Install / upgrade PyInstaller inside the env ───────────────────────
echo "==> Installing PyInstaller into ${CONDA_ENV}..."
conda run -n "${CONDA_ENV}" pip install --quiet --upgrade pyinstaller

# ── 3. Sync conf into climdata_gui/ ───────────────────────────────────────
echo "==> Syncing climdata/conf → climdata_gui/conf..."
rm -rf climdata_gui/conf
cp -r climdata/conf climdata_gui/conf

# ── 4. Build the .app bundle ──────────────────────────────────────────────
echo "==> Running PyInstaller..."
conda run -n "${CONDA_ENV}" pyinstaller "${SPEC_FILE}" --noconfirm

# ── 4. Fix QtWebEngineProcess search paths ────────────────────────────────
# conda-forge PySide6 places QtWebEngineProcess under
#   Contents/Frameworks/PySide6/Qt/lib/qt6/QtWebEngineProcess
# but Qt's WebEngineLibraryInfo looks in libexec/, bin/, and MacOS/.
# Create symlinks so all three searched paths resolve correctly.
APP_CONTENTS="dist/climdata_gui.app/Contents"
QT_BASE="${APP_CONTENTS}/Frameworks/PySide6/Qt"
HELPER_REL="../lib/qt6/QtWebEngineProcess"

echo "==> Fixing QtWebEngineProcess paths..."
mkdir -p "${QT_BASE}/libexec" "${QT_BASE}/bin"
ln -sf "${HELPER_REL}" "${QT_BASE}/libexec/QtWebEngineProcess"
ln -sf "${HELPER_REL}" "${QT_BASE}/bin/QtWebEngineProcess"
ln -sf "../Frameworks/PySide6/Qt/lib/qt6/QtWebEngineProcess" \
       "${APP_CONTENTS}/MacOS/QtWebEngineProcess"

# ── 5. Sanity-check output ────────────────────────────────────────────────
APP_PATH="dist/climdata_gui.app"
if [[ -d "${APP_PATH}" ]]; then
    SIZE=$(du -sh "${APP_PATH}" | cut -f1)
    echo ""
    echo "✓ Build succeeded!"
    echo "  Output : ${APP_PATH}  (${SIZE})"
    echo ""
    echo "  To launch:"
    echo "    open ${APP_PATH}"
    echo ""
    echo "  NOTE: The app is unsigned. On a new Mac, right-click → Open"
    echo "  to bypass Gatekeeper the first time."
else
    echo "ERROR: Expected ${APP_PATH} was not created. Check the PyInstaller output above."
    exit 1
fi
