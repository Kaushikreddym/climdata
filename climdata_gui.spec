# -*- mode: python ; coding: utf-8 -*-
# PyInstaller spec for climdata_gui — macOS arm64 .app bundle
#
# Run from the workspace root:
#   pyinstaller climdata_gui.spec --noconfirm
#
# All dependencies are bundled; no Python/conda installation required
# on the target machine.

from PyInstaller.utils.hooks import collect_all, collect_submodules, collect_data_files

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

a = Analysis(
    # Entry point — climdata_gui/main.py
    ['climdata_gui/main.py'],

    # Add climdata_gui/ to sys.path so top-level imports like
    # "from app import create_app" and "from gui.main_window import …" work.
    # Mirrors what _ensure_on_path() does at runtime.
    pathex=['climdata_gui'],

    binaries=(
        hydra_bins + omegaconf_bins + mpl_bins +
        cartopy_bins + pyproj_bins + xclim_bins + climdata_bins
    ),

    datas=[
        # GUI resources (Leaflet HTML/JS/CSS, marker PNGs)
        # Destination matches paths.resource_path("resources/html/…") in _MEIPASS
        ('climdata_gui/resources',   'resources'),

        # Hydra config (dataset list, variable mappings, regions)
        ('climdata_gui/conf',        'conf'),

        # QSS theme loaded by app.py via resource_path("gui/styles/theme.qss")
        ('climdata_gui/gui/styles',  'gui/styles'),
    ] + hydra_datas + omegaconf_datas + mpl_datas
      + cartopy_datas + pyproj_datas + xclim_datas + climdata_datas,

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

    ] + hydra_hidden + omegaconf_hidden + mpl_hidden
      + cartopy_hidden + pyproj_hidden + xclim_hidden + climdata_hidden,

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
    strip=False,
    upx=False,          # UPX is unreliable with macOS arm64 binaries
    console=False,      # Suppress terminal window (windowed app)
    disable_windowed_traceback=False,
    target_arch='arm64',
    codesign_identity=None,
    entitlements_file=None,
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

app = BUNDLE(
    coll,
    name='climdata_gui.app',
    icon=None,          # No .icns available — macOS uses a generic icon
    bundle_identifier='de.zalf.climdata',
    info_plist={
        'CFBundleDisplayName': 'climdata GUI',
        'CFBundleShortVersionString': '0.6.0',
        'NSHighResolutionCapable': True,
        'NSRequiresAquaSystemAppearance': False,   # Supports dark mode
        # Suppress "App downloaded from internet" quarantine warnings in dev
        'LSEnvironment': {
            'QTWEBENGINE_DISABLE_SANDBOX': '1',
        },
    },
)
