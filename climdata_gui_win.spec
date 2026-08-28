# -*- mode: python ; coding: utf-8 -*-
# PyInstaller spec for climdata_gui — Windows x86_64 one-folder distribution
#
# Run from the workspace root (Windows, inside the conda env):
#   pyinstaller climdata_gui_win.spec --noconfirm
# or via the wrapper:
#   bash build/build_windows.sh
#
# Output: dist/climdata_gui/climdata_gui.exe
#
# Mirrors climdata_gui.spec (macOS) but drops the .app BUNDLE, the arm64
# target_arch, and adds Windows-only DLL/plugin handling.

import os

from PyInstaller.utils.hooks import collect_all, collect_data_files, collect_submodules

block_cipher = None

# ── Collect packages that ship non-Python data files ──────────────────────
# collect_all() returns (datas, binaries, hiddenimports) for each package.
hydra_datas,     hydra_bins,     hydra_hidden     = collect_all('hydra')
omegaconf_datas, omegaconf_bins, omegaconf_hidden = collect_all('omegaconf')
mpl_datas,       mpl_bins,       mpl_hidden       = collect_all('matplotlib')
cartopy_datas,   cartopy_bins,   cartopy_hidden    = collect_all('cartopy')
pyproj_datas,    pyproj_bins,    pyproj_hidden     = collect_all('pyproj')
xclim_datas,     xclim_bins,     xclim_hidden      = collect_all('xclim')
climdata_datas,  climdata_bins,  climdata_hidden   = collect_all('climdata')

# ── Dataset-specific packages (namespace / plugin / data-file heavy) ───────
# google-api-python-client + google-auth are namespace packages; PyInstaller's
# static tracer cannot follow them without explicit collect_all.
google_datas,       google_bins,       google_hidden       = collect_all('googleapiclient')
# google.auth / google.oauth2 live in the same namespace; collect via submodules
googleauth_datas,   googleauth_bins,   googleauth_hidden   = [], [], collect_submodules('google.auth')
googleoauth_datas,  googleoauth_bins,  googleoauth_hidden  = [], [], collect_submodules('google.oauth2')

intake_datas,       intake_bins,       intake_hidden       = collect_all('intake')
intakeesm_datas,    intakeesm_bins,    intakeesm_hidden    = collect_all('intake_esm')
wetterdienst_datas, wetterdienst_bins, wetterdienst_hidden = collect_all('wetterdienst')
cdsapi_datas,       cdsapi_bins,       cdsapi_hidden       = collect_all('cdsapi')
rasterio_datas,     rasterio_bins,     rasterio_hidden     = collect_all('rasterio')
cfxarray_datas,     cfxarray_bins,     cfxarray_hidden     = collect_all('cf_xarray')
yamale_datas,       yamale_bins,       yamale_hidden       = collect_all('yamale')

# ── Windows-only: pull in the GDAL/PROJ data trees that conda ships in the
# environment root (Library/share) rather than inside the wheel. rasterio and
# pyproj resolve these at runtime via GDAL_DATA / PROJ_DATA; without them the
# frozen app fails on the first CRS lookup.
conda_datas = []
_env_root = os.environ.get('CONDA_PREFIX', '')
if _env_root:
    for _sub, _dest in (
        (os.path.join('Library', 'share', 'gdal'), 'Library/share/gdal'),
        (os.path.join('Library', 'share', 'proj'), 'Library/share/proj'),
    ):
        _src = os.path.join(_env_root, _sub)
        if os.path.isdir(_src):
            conda_datas.append((_src, _dest))

# Ship an icon only if one has been added to the repo.
_icon = os.path.join('climdata_gui', 'resources', 'icons', 'climdata.ico')
icon_path = _icon if os.path.isfile(_icon) else None

a = Analysis(
    # Entry point — climdata_gui/main.py
    ['climdata_gui/main.py'],

    # Add climdata_gui/ to sys.path so top-level imports like
    # "from app import create_app" and "from gui.main_window import …" work.
    # Mirrors what _ensure_on_path() does at runtime.
    pathex=['climdata_gui'],

    binaries=(
        hydra_bins + omegaconf_bins + mpl_bins +
        cartopy_bins + pyproj_bins + xclim_bins + climdata_bins +
        google_bins + googleauth_bins + googleoauth_bins +
        intake_bins + intakeesm_bins +
        wetterdienst_bins + cdsapi_bins +
        rasterio_bins + cfxarray_bins +
        yamale_bins
    ),

    datas=[
        # GUI resources (Leaflet HTML/JS/CSS, marker PNGs)
        # Destination matches paths.resource_path("resources/html/…") in _MEIPASS
        ('climdata_gui/resources',   'resources'),

        # Hydra config (dataset list, variable mappings, regions)
        ('climdata_gui/conf',        'conf'),

        # QSS theme loaded by app.py via resource_path("gui/styles/theme.qss")
        ('climdata_gui/gui/styles',  'gui/styles'),
    ] + collect_data_files('certifi')   # CA bundle — see utils/certs.py
      + hydra_datas + omegaconf_datas + mpl_datas
      + cartopy_datas + pyproj_datas + xclim_datas + climdata_datas
      + google_datas + googleauth_datas + googleoauth_datas
      + intake_datas + intakeesm_datas
      + wetterdienst_datas + cdsapi_datas
      + rasterio_datas + cfxarray_datas
      + yamale_datas
      + conda_datas,

    hiddenimports=[
        # ── PySide6 / QtWebEngine ──────────────────────────────────────────
        # These are loaded dynamically by map_widget.py; not detected by
        # PyInstaller's static import scanner.
        'PySide6.QtWebEngineWidgets',
        'PySide6.QtWebEngineCore',
        'PySide6.QtWebEngineQuick',
        'PySide6.QtWebChannel',
        'PySide6.QtNetwork',
        'PySide6.QtPositioning',
        'PySide6.QtPrintSupport',

        # ── GUI top-level modules (resolved via pathex=['climdata_gui']) ───
        'app',
        'gui',
        'gui.main_window',
        'gui.map_widget',
        'gui.controls',
        'gui.controls.dataset_selector',
        'gui.controls.date_range_picker',
        'gui.controls.yaml_editor',
        'backend',
        'backend.config_builder',
        'backend.runner',
        'backend.worker',
        'models',
        'models.app_state',
        'utils',
        'utils.paths',

        # ── Windows-only stdlib bits pulled in dynamically ────────────────
        # multiprocessing uses the spawn start method on Windows; dask and
        # the worker thread rely on it being importable in the frozen app.
        'multiprocessing.popen_spawn_win32',
        'multiprocessing.spawn',
        'encodings.idna',
        'win32timezone',

        # ── Scientific stack — lazy/plugin-loaded submodules ──────────────
        *collect_submodules('geopandas'),
        *collect_submodules('rioxarray'),
        *collect_submodules('xarray'),
        *collect_submodules('dask'),
        *collect_submodules('dask.dataframe'),
        *collect_submodules('fiona'),
        *collect_submodules('shapely'),
        *collect_submodules('pandas'),
        *collect_submodules('numpy'),
        *collect_submodules('scipy'),
        *collect_submodules('zarr'),
        *collect_submodules('wetterdienst'),
        *collect_submodules('pint'),
        *collect_submodules('cdsapi'),
        *collect_submodules('aiohttp'),
        *collect_submodules('rasterio'),
        *collect_submodules('cf_xarray'),
        *collect_submodules('intake'),
        *collect_submodules('intake_esm'),
        *collect_submodules('seaborn'),
        *collect_submodules('sklearn'),
        *collect_submodules('statsmodels'),
        *collect_submodules('xsdba'),          # BASD tab: QDM / DQM / QM
        'certifi',                             # TLS trust store (utils/certs.py)
        *collect_submodules('pyarrow'),

    ] + hydra_hidden + omegaconf_hidden + mpl_hidden
      + cartopy_hidden + pyproj_hidden + xclim_hidden + climdata_hidden
      + google_hidden + googleauth_hidden + googleoauth_hidden
      + intake_hidden + intakeesm_hidden
      + wetterdienst_hidden + cdsapi_hidden
      + rasterio_hidden + cfxarray_hidden
      + yamale_hidden,

    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],

    # Strip unused GUI toolkits to reduce bundle size
    excludes=['tkinter', 'PyQt5', 'PyQt6', 'wx', '_tkinter'],

    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='climdata_gui',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,          # `strip` needs binutils; not available on Windows
    upx=False,            # UPX corrupts Qt6 DLLs
    console=False,        # Suppress console window (windowed app)
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,     # Follows the interpreter (x86_64)
    codesign_identity=None,
    entitlements_file=None,
    icon=icon_path,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name='climdata_gui',
)
