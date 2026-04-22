## 🏗️ TASK: Create GUI Subproject (`climdata_gui`)

Claude must create a new subdirectory inside the repository root:

```
climdata_gui/
```

This directory contains a standalone GUI application built with PySide6 that interfaces with the existing `climdata` package.

---

## 📁 Required Folder Structure

Create the following structure exactly:

```
climdata_gui/
│
├── main.py
├── app.py
│
├── gui/
│   ├── main_window.py
│   ├── map_widget.py
│   ├── controls/
│   │   ├── dataset_selector.py
│   │   └── date_range_picker.py
│   └── styles/
│       └── theme.qss
│
├── backend/
│   ├── runner.py
│   ├── worker.py
│   └── config_builder.py
│
├── models/
│   └── app_state.py
│
├── resources/
│   └── html/
│       └── map.html
│
├── utils/
│   └── paths.py
│
└── requirements.txt
```

---

## ⚙️ Implementation Requirements

### 1. `main.py`

* Entry point
* Initializes QApplication
* Launches MainWindow

---

### 2. GUI Requirements

#### `main_window.py`

Must include:

* Dataset dropdown
* Map widget (for AOI selection)
* Date range picker
* Run button
* Log/output panel

---

#### `map_widget.py`

* Use Qt WebEngine
* Load `resources/html/map.html`
* Use Qt WebChannel bridge
* Emit selected `(lat, lon)` to Python

---

### 3. Backend Layer

#### `config_builder.py`

* Build Hydra overrides from:

  * dataset
  * lat/lon
  * date range

---

#### `runner.py`

Wrap:

```python
from climdata import ClimData
```

Expose:

```python
run_pipeline(overrides, seq)
```

---

#### `worker.py`

* Use Qt threading (QObject + Signal)
* Run pipeline asynchronously
* Emit:

  * finished
  * error

---

### 4. Map HTML (`resources/html/map.html`)

* Use Leaflet
* Allow click-to-select point
* Send coordinates via Qt WebChannel

---

### 5. Utilities

#### `utils/paths.py`

Provide:

```python
resource_path(path)
```

Compatible with PyInstaller

---

### 6. Requirements File

Must include:

```
PySide6
PySide6-WebEngine
```

---

## 🔌 Integration Rules

* Import `climdata` as an installed package
* Do NOT duplicate backend logic
* GUI must call backend via `runner.py`

---

## 🧪 Minimum Working Requirement

After creation, the following must work:

```bash
cd climdata_gui
python main.py
```

This should:

* Open GUI window
* Allow selecting AOI from map
* Allow selecting date range
* Print constructed overrides in console/log

---

## 🚫 Constraints

* Do NOT modify existing `climdata` code
* Do NOT hardcode absolute paths
* Keep modules small and modular
* Ensure code is runnable without packaging

---

## 🎯 Goal of This Task

Bootstrap a **fully functional GUI scaffold** that:

* Is cleanly separated from backend
* Can be extended incrementally
* Is ready for packaging into `.exe`
