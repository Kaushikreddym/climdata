# climdata_gui — Installation Guide

## Requirements

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda
- An existing conda environment (e.g. `climdata_dev`) with Python 3.12

---

## Install PySide6 and Qt6 WebEngine

**Do not use `pip install PySide6` alone.** PySide6 6.8+ pip wheels on macOS ship Qt
libraries as `.framework` bundles, but the Python binding `.so` files link against
`libQt6*.dylib` — causing an `ImportError` at runtime. Use conda-forge instead.

```bash
conda install -n climdata_dev -c conda-forge pyside6 qt6-webengine --yes
```

This installs:

| Package | Provides |
|---|---|
| `pyside6` | Python bindings: `QtWidgets`, `QtCore`, `QtWebChannel`, `QtWebEngineWidgets`, … |
| `qt6-webengine` | Native `libQt6WebEngine*.dylib` libraries (required by `map_widget.py`) |

> **Note:** `qt6-main` is pulled in automatically as a dependency of `pyside6`.

---

## Run the GUI

```bash
cd climdata_gui
conda run -n climdata_dev python main.py
```

Or activate the environment first:

```bash
conda activate climdata_dev
cd climdata_gui
python main.py
```

---

## Smoke test (no climdata stack needed)

```bash
conda run -n climdata_dev python test_gui/app.py
```

Opens a bare Qt window — confirms PySide6 and the display are working before
testing the full app.

---

## Troubleshooting

### `ImportError: Library not loaded: @rpath/libQt6WebEngineWidgets.6.dylib`

The pip-installed `PySide6-Addons` is present and conflicting with the conda-forge
version. Remove it and reinstall cleanly:

```bash
# Remove pip packages (no RECORD file — must delete manually)
SITE=$(conda run -n climdata_dev python -c "import site; print(site.getsitepackages()[0])")
rm -rf "$SITE/PySide6" "$SITE/shiboken6"

# Reinstall from conda-forge
conda install -n climdata_dev -c conda-forge pyside6 qt6-webengine --force-reinstall --yes
```

### Verify installed packages

```bash
conda list -n climdata_dev | grep -i "pyside\|qt6"
```

Expected output:
```
pyside6          6.x.x   py312...   conda-forge
qt6-main         6.x.x   ...        conda-forge
qt6-webengine    6.x.x   ...        conda-forge
```
