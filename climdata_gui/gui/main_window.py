"""Top-level GUI window wiring controls, map, and pipeline worker."""

from __future__ import annotations

from typing import Optional

from PySide6.QtCore import QThread, QTimer, Qt
from PySide6.QtWidgets import (
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QPlainTextEdit,
    QPushButton,
    QSplitter,
    QStatusBar,
    QVBoxLayout,
    QWidget,
)

from backend.config_builder import build_overrides
from backend.runner import get_datasets, get_overlay_meta, render_variable
from backend.worker import PipelineWorker
from models.app_state import AppState

from .controls import DatasetSelector, DateRangePicker
from .controls.yaml_editor import YamlEditorDialog
from .map_widget import MapWidget


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("climdata")
        self.resize(1100, 700)
        self.showMaximized()

        self._state = AppState()
        self._thread: Optional[QThread] = None
        self._worker: Optional[PipelineWorker] = None
        self._last_result = None

        self._build_ui()
        self._connect_signals()

    # ------------------------------------------------------------------ UI
    def _build_ui(self) -> None:
        central = QWidget(self)
        outer = QVBoxLayout(central)

        # Top controls row
        controls_row = QHBoxLayout()
        datasets = get_datasets() or None  # None → DatasetSelector uses DEFAULT_DATASETS
        self._dataset = DatasetSelector(datasets=datasets)
        self._dates = DateRangePicker()
        self._adv_btn = QPushButton("⚙")
        self._adv_btn.setFixedWidth(32)
        self._adv_btn.setToolTip("Advanced: edit config YAML")
        self._format_combo = QComboBox()
        self._format_combo.addItems(["CSV", "NetCDF"])
        self._format_combo.setToolTip("Output format")
        self._run_btn = QPushButton("Run")
        self._plot_btn = QPushButton("Plot")
        self._plot_btn.setToolTip("Render time-mean overlay on map")
        self._plot_btn.setEnabled(False)
        controls_row.addWidget(self._dataset, 1)
        controls_row.addWidget(self._adv_btn, 0)
        controls_row.addWidget(self._dates, 2)
        controls_row.addWidget(QLabel("Format:"), 0)
        controls_row.addWidget(self._format_combo, 0)
        controls_row.addWidget(self._run_btn, 0)
        controls_row.addWidget(self._plot_btn, 0)
        outer.addLayout(controls_row)

        # Data directory row
        dir_row = QHBoxLayout()
        dir_row.addWidget(QLabel("Data dir:"), 0)
        self._data_dir_edit = QLineEdit()
        self._data_dir_edit.setPlaceholderText("./data  (default)")
        self._data_dir_edit.setReadOnly(True)
        self._data_dir_btn = QPushButton("\U0001f4c2")
        self._data_dir_btn.setFixedWidth(36)
        self._data_dir_btn.setMinimumWidth(36)
        self._data_dir_btn.setToolTip("Choose data directory")
        dir_row.addWidget(self._data_dir_edit, 1)
        dir_row.addWidget(self._data_dir_btn, 0)
        outer.addLayout(dir_row)

        # Coordinate label
        self._coord_label = QLabel("No point selected")
        self._coord_label.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        outer.addWidget(self._coord_label)

        # Splitter: map on top, log panel below
        splitter = QSplitter(Qt.Vertical)
        self._map = MapWidget()
        splitter.addWidget(self._map)

        self._log = QPlainTextEdit()
        self._log.setReadOnly(True)
        self._log.setPlaceholderText("Pipeline log / output will appear here...")
        splitter.addWidget(self._log)
        splitter.setStretchFactor(0, 10)
        splitter.setStretchFactor(1, 1)
        outer.addWidget(splitter, 1)

        self.setCentralWidget(central)
        self.setStatusBar(QStatusBar(self))

        # Seed initial state from defaults.
        self._state.dataset = self._dataset.current_dataset()
        self._state.start_date = self._dates.start_date()
        self._state.end_date = self._dates.end_date()

    def _connect_signals(self) -> None:
        self._dataset.dataset_changed.connect(self._on_dataset_changed)
        self._dates.range_changed.connect(self._on_dates_changed)
        self._map.point_selected.connect(self._on_point_selected)
        self._map.box_selected.connect(self._on_box_selected)
        self._adv_btn.clicked.connect(self._on_advanced_clicked)
        self._format_combo.currentTextChanged.connect(self._on_format_changed)
        self._run_btn.clicked.connect(self._on_run_clicked)
        self._plot_btn.clicked.connect(self._on_plot_clicked)
        self._map.render_requested.connect(self._on_render_requested)
        self._data_dir_btn.clicked.connect(self._on_data_dir_clicked)

    # -------------------------------------------------------------- Logging
    def log(self, msg: str) -> None:
        self._log.appendPlainText(msg)

    # --------------------------------------------------------------- Slots
    def _on_format_changed(self, fmt: str) -> None:
        if fmt == "NetCDF":
            self._state.actions = ["extract", "to_nc"]
        else:
            self._state.actions = ["extract", "to_csv"]
        self.log(f"Output format set to: {fmt}")

    def _on_dataset_changed(self, dataset: str) -> None:
        self._state.dataset = dataset
        self.log(f"Dataset set to: {dataset}")

    def _on_dates_changed(self, start, end) -> None:
        self._state.start_date = start
        self._state.end_date = end
        self.log(f"Date range: {start.isoformat()} .. {end.isoformat()}")

    def _on_point_selected(self, lat: float, lon: float) -> None:
        self._state.lat = lat
        self._state.lon = lon
        self._state.box = None           # clear any previous box
        self._coord_label.setText(f"Point: lat={lat:.4f}, lon={lon:.4f}")
        self.log(f"AOI point: lat={lat:.4f}, lon={lon:.4f}")

    def _on_box_selected(self, lat_min: float, lat_max: float,
                         lon_min: float, lon_max: float) -> None:
        self._state.box = {
            "lat_min": lat_min, "lat_max": lat_max,
            "lon_min": lon_min, "lon_max": lon_max,
        }
        self._state.lat = None           # clear any previous point
        self._state.lon = None
        self._coord_label.setText(
            f"Box: ({lat_min:.3f}, {lon_min:.3f}) \u2192 ({lat_max:.3f}, {lon_max:.3f})"
        )
        self.log(f"AOI box: lat [{lat_min:.3f}, {lat_max:.3f}]  lon [{lon_min:.3f}, {lon_max:.3f}]")

    def _on_advanced_clicked(self) -> None:
        # Only pass simple scalar overrides — complex flattened values break Hydra grammar
        current_overrides = build_overrides(
            dataset=self._state.dataset,
            lat=self._state.lat,
            lon=self._state.lon,
            start_date=self._state.start_date,
            end_date=self._state.end_date,
            variables=self._state.variables,
        )
        dlg = YamlEditorDialog(self, overrides=current_overrides)
        if dlg.exec() == YamlEditorDialog.DialogCode.Accepted:
            self._state.extra_overrides = dlg.extra_overrides()
            fields = dlg.parsed_fields()

            if "dataset" in fields:
                self._dataset.set_dataset(fields["dataset"])
                # set_dataset triggers dataset_changed signal which updates state

            if "lat" in fields and "lon" in fields:
                self._state.lat = fields["lat"]
                self._state.lon = fields["lon"]
                self._coord_label.setText(
                    f"Selected point: lat={fields['lat']:.4f}, lon={fields['lon']:.4f}"
                )

            if "start_date" in fields or "end_date" in fields:
                from datetime import date as _date
                def _parse(s):
                    parts = str(s).split("-")
                    return _date(int(parts[0]), int(parts[1]), int(parts[2]))
                try:
                    start = _parse(fields["start_date"]) if "start_date" in fields else self._state.start_date
                    end = _parse(fields["end_date"]) if "end_date" in fields else self._state.end_date
                    self._dates.set_dates(start, end)
                    # set_dates triggers range_changed signal which updates state
                except Exception:
                    pass

            self.log(f"Advanced config saved ({len(self._state.extra_overrides)} overrides).")

    def _on_run_clicked(self) -> None:
        if not self._state.is_ready():
            self.log("Cannot run: select a dataset, a point on the map, and a date range.")
            return
        if self._thread is not None:
            self.log("Pipeline already running.")
            return

        overrides = build_overrides(
            dataset=self._state.dataset,
            lat=self._state.lat,
            lon=self._state.lon,
            start_date=self._state.start_date,
            end_date=self._state.end_date,
            variables=self._state.variables,
            extra_overrides=self._state.extra_overrides or None,
            box=self._state.box,
            data_dir=self._state.data_dir,
        )
        self.log("Constructed Hydra overrides:")
        for o in overrides:
            self.log(f"  {o}")

        self._start_worker(overrides)

    # -------------------------------------------------------------- Worker
    def _start_worker(self, overrides) -> None:
        self._thread = QThread(self)
        self._worker = PipelineWorker(overrides, self._state.actions)
        self._worker.moveToThread(self._thread)

        self._thread.started.connect(self._worker.run)
        self._worker.log.connect(self.log)
        self._worker.error.connect(self._on_worker_error)
        self._worker.finished.connect(self._on_worker_finished)
        self._worker.finished.connect(self._thread.quit)
        self._worker.error.connect(self._thread.quit)
        self._thread.finished.connect(self._cleanup_worker)

        self._run_btn.setEnabled(False)
        self._thread.start()

    def _on_worker_finished(self, result) -> None:
        self.log("Pipeline completed successfully.")
        if result is not None and hasattr(result, "filename") and result.filename:
            self.log(f"Output: {result.filename}")
        self._last_result = result
        self._plot_btn.setEnabled(result is not None)

    def _on_plot_clicked(self) -> None:
        if self._last_result is None:
            self.log("No result available — run the pipeline first.")
            return
        meta = get_overlay_meta(self._last_result)
        if meta is None:
            self.log("Cannot plot: dataset has no spatial grid dimensions.")
            return
        self.log(f"Initialising overlay for variables: {', '.join(meta['vars'])}")
        self._map.init_overlay_ui(meta)
        # _on_render_requested is called by JS immediately after initOverlayUI

    def _on_render_requested(self, var_name: str, colormap: str, vmin_str: str, vmax_str: str) -> None:
        """Render one variable on demand and push the result to the map."""
        if self._last_result is None or self._last_result.dataset is None:
            return
        vmin = None if vmin_str == 'auto' else float(vmin_str)
        vmax = None if vmax_str == 'auto' else float(vmax_str)
        self.log(f"Rendering '{var_name}' ({colormap}) vmin={vmin_str} vmax={vmax_str}...")
        payload = render_variable(self._last_result.dataset, var_name, colormap, vmin=vmin, vmax=vmax)
        if payload:
            # Defer runJavaScript out of the QWebChannel callback so Qt's
            # JS engine is not in a re-entrant state when we push the update.
            QTimer.singleShot(0, lambda p=payload: self._map.update_overlay(p))
            self.log(f"Overlay updated: {var_name} [{colormap}]")
        else:
            self.log(f"Render failed for '{var_name}' — no spatial dims.")

    def _on_data_dir_clicked(self) -> None:
        current = self._state.data_dir or ""
        chosen = QFileDialog.getExistingDirectory(
            self, "Select Data Directory", current or "."
        )
        if chosen:
            self._state.data_dir = chosen
            self._data_dir_edit.setText(chosen)
            self.log(f"Data directory: {chosen}")

    def _on_worker_error(self, tb: str) -> None:
        self.log("Pipeline failed:")
        self.log(tb)

    def _cleanup_worker(self) -> None:
        if self._worker is not None:
            self._worker.deleteLater()
            self._worker = None
        if self._thread is not None:
            self._thread.deleteLater()
            self._thread = None
        self._run_btn.setEnabled(True)
