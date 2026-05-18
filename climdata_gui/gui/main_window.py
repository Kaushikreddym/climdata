"""Top-level GUI window wiring controls, map, and pipeline worker."""

from __future__ import annotations

from typing import Optional

from PySide6.QtCore import QThread, Qt
from PySide6.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QPlainTextEdit,
    QPushButton,
    QSplitter,
    QStatusBar,
    QVBoxLayout,
    QWidget,
)

from backend.config_builder import build_overrides
from backend.worker import PipelineWorker
from models.app_state import AppState

from .controls import DatasetSelector, DateRangePicker
from .map_widget import MapWidget


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("climdata")
        self.resize(1100, 700)

        self._state = AppState()
        self._thread: Optional[QThread] = None
        self._worker: Optional[PipelineWorker] = None

        self._build_ui()
        self._connect_signals()

    # ------------------------------------------------------------------ UI
    def _build_ui(self) -> None:
        central = QWidget(self)
        outer = QVBoxLayout(central)

        # Top controls row
        controls_row = QHBoxLayout()
        self._dataset = DatasetSelector()
        self._dates = DateRangePicker()
        self._run_btn = QPushButton("Run")
        controls_row.addWidget(self._dataset, 1)
        controls_row.addWidget(self._dates, 2)
        controls_row.addWidget(self._run_btn, 0)
        outer.addLayout(controls_row)

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
        splitter.setStretchFactor(0, 3)
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
        self._run_btn.clicked.connect(self._on_run_clicked)

    # -------------------------------------------------------------- Logging
    def log(self, msg: str) -> None:
        self._log.appendPlainText(msg)

    # --------------------------------------------------------------- Slots
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
        self._coord_label.setText(f"Selected point: lat={lat:.4f}, lon={lon:.4f}")
        self.log(f"AOI point: lat={lat:.4f}, lon={lon:.4f}")

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
