# Packaging scripts

Standalone desktop bundles of `climdata_gui`, built with PyInstaller. Each
script runs on the platform it targets — PyInstaller does not cross-compile,
so the Windows `.exe` must be built on Windows and the macOS `.app` on macOS.

| Script             | Host    | Spec                    | Output                             |
| ------------------ | ------- | ----------------------- | ---------------------------------- |
| `build_mac.sh`     | macOS   | `climdata_gui.spec`     | `dist/climdata_gui.app`            |
| `build_windows.sh` | Windows | `climdata_gui_win.spec` | `dist/climdata_gui/climdata_gui.exe` |

Both scripts `cd` to the repository root themselves, so they can be invoked
from anywhere. Both sync `climdata/conf` into `climdata_gui/conf` before
building, and both write PyInstaller's scratch tree to `build_work/` so it
does not collide with this directory.

## Windows

Run from **Git Bash**, MSYS2, or any bash shell where `conda` is on `PATH`:

```bash
bash build/build_windows.sh --env sdba        # build with the `sdba` env
bash build/build_windows.sh --check --env sdba  # pre-flight only, no build
bash build/build_windows.sh --env sdba --zip    # also emit a shippable .zip
```

The environment can also be set with `CLIMDATA_ENV`. The default is
`climdata_dev`.

`--check` verifies the conda environment, that every runtime dependency
imports, and that the spec, entry point and resource directories exist —
without running PyInstaller. Use it first; it takes seconds, whereas a full
build takes several minutes.

Notes specific to the Windows build:

- The spec bundles the conda environment's `Library/share/gdal` and
  `Library/share/proj` data trees. `rasterio` and `pyproj` look these up at
  runtime, and conda keeps them outside the wheels, so without them the
  frozen app fails on the first CRS lookup.
- After the build the script checks that `QtWebEngineProcess.exe` made it
  into the bundle. It renders the Leaflet map panel, and a missing helper
  shows up as a silently blank map rather than an error.
- The `.exe` is unsigned, so SmartScreen shows "Windows protected your PC"
  on first launch — *More info → Run anyway*.
- Drop a `climdata_gui/resources/icons/climdata.ico` into the repo and the
  spec picks it up automatically; without one Windows uses a generic icon.

## macOS

```bash
bash build/build_mac.sh                  # uses the `climdata_dev` env
CLIMDATA_ENV=sdba bash build/build_mac.sh
```

The script fixes up `QtWebEngineProcess` symlinks after the build:
conda-forge's PySide6 puts the helper under `Qt/lib/qt6/`, but Qt searches
`libexec/`, `bin/` and `Contents/MacOS/`.

The `.app` is unsigned — right-click → Open on first launch to get past
Gatekeeper.
