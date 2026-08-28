"""Top-level window — thin QMainWindow wrapper around the web dashboard."""

from __future__ import annotations

from datetime import date
from typing import Optional

from PySide6.QtCore import QThread, QTimer
from PySide6.QtWidgets import QFileDialog, QMainWindow, QStatusBar

from backend.config_builder import build_overrides
from backend.runner import (
    CMIP_DATASETS,
    get_datasets,
    get_overlay_meta,
    is_cmip_dataset,
    render_variable,
)
from backend.worker import CmipOptionsWorker, PipelineWorker
from gui.controls.dataset_selector import DEFAULT_DATASETS
from models.app_state import AppState

from .controls.yaml_editor import YamlEditorDialog
from .map_widget import MapWidget


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("climdata")
        self.showMaximized()

        self._state = AppState()
        # Match the web toolbar's default values so state is in sync before first interaction
        self._state.start_date = date(1989, 1, 1)
        self._state.end_date   = date(2020, 12, 31)

        self._thread: Optional[QThread] = None
        self._worker: Optional[PipelineWorker] = None
        self._last_result = None

        # CMIP6 catalogue lookup (model / experiment pickers)
        self._opt_thread: Optional[QThread] = None
        self._opt_worker: Optional[CmipOptionsWorker] = None
        self._pending_experiment: Optional[str] = None

        self._map = MapWidget()
        self.setCentralWidget(self._map)
        self.setStatusBar(QStatusBar(self))

        self._connect_signals()

    def _connect_signals(self) -> None:
        self._map.page_ready.connect(self._on_page_ready)
        self._map.point_selected.connect(self._on_point_selected)
        self._map.box_selected.connect(self._on_box_selected)
        self._map.render_requested.connect(self._on_render_requested)
        self._map.dataset_changed.connect(self._on_dataset_changed)
        self._map.dates_changed.connect(self._on_dates_changed)
        self._map.format_changed.connect(self._on_format_changed)
        self._map.experiment_changed.connect(self._on_experiment_changed)
        self._map.model_changed.connect(self._on_model_changed)
        self._map.run_clicked.connect(self._on_run_clicked)
        self._map.plot_clicked.connect(self._on_plot_clicked)
        self._map.advanced_clicked.connect(self._on_advanced_clicked)
        self._map.data_dir_browse_requested.connect(self._on_data_dir_browse)

    # ── Logging ───────────────────────────────────────────────────────
    def log(self, msg: str) -> None:
        self._map.send_log(msg)

    # ── Page ready — push initial state to web frontend ───────────────
    def _on_page_ready(self) -> None:
        # Normalise to lowercase so syncToolbar's querySelector always matches
        # (get_datasets() returns uppercase keys from parameters.yaml)
        datasets = [d.lower() for d in (get_datasets() or DEFAULT_DATASETS)]
        if datasets and self._state.dataset not in datasets:
            self._state.dataset = datasets[0]
        self._map.dashboard_ready(
            datasets=datasets,
            dataset=self._state.dataset,
            start=self._state.start_date.isoformat(),
            end=self._state.end_date.isoformat(),
            data_dir=self._state.data_dir or "",
            cmip_datasets=list(CMIP_DATASETS),
        )
        self._map.send_plot_enabled(False)
        self.log("climdata dashboard ready.")

    # ── Toolbar slots (JS → Python) ───────────────────────────────────
    def _on_dataset_changed(self, dataset: str) -> None:
        previous = self._state.dataset
        self._state.dataset = dataset
        if is_cmip_dataset(dataset):
            # Only (re)query when switching in — avoids a catalogue hit on every
            # echo of the current selection from the web frontend.
            if not is_cmip_dataset(previous) or not self._state.source_id:
                self._load_cmip_options(self._state.experiment_id)
        else:
            self._state.experiment_id = None
            self._state.source_id = None

    def _on_experiment_changed(self, experiment_id: str) -> None:
        if not experiment_id or experiment_id == self._state.experiment_id:
            return
        self._state.experiment_id = experiment_id
        # The model list is experiment-dependent — refresh it.
        self._load_cmip_options(experiment_id)

    def _on_model_changed(self, source_id: str) -> None:
        if source_id:
            self._state.source_id = source_id

    # ── CMIP6 catalogue lookup ────────────────────────────────────────
    def _load_cmip_options(self, experiment: Optional[str] = None) -> None:
        """Fetch experiment/model lists for CMIP on a background thread."""
        if self._opt_thread is not None:
            # A lookup is in flight — remember the newest request and run it next.
            self._pending_experiment = experiment
            return
        self._pending_experiment = None
        self._map.set_cmip_options([], experiment or "", [], "",
                                   loading=True, note="Loading CMIP6 catalogue…")
        self.log("Querying CMIP6 catalogue for models and experiments…")

        self._opt_thread = QThread(self)
        self._opt_worker = CmipOptionsWorker(experiment)
        self._opt_worker.moveToThread(self._opt_thread)

        self._opt_thread.started.connect(self._opt_worker.run)
        self._opt_worker.log.connect(self.log)
        self._opt_worker.finished.connect(self._on_cmip_options_ready)
        self._opt_worker.error.connect(self._on_cmip_options_error)
        self._opt_worker.finished.connect(self._opt_thread.quit)
        self._opt_worker.error.connect(self._opt_thread.quit)
        self._opt_thread.finished.connect(self._cleanup_opt_worker)
        self._opt_thread.start()

    def _on_cmip_options_ready(self, opts: dict) -> None:
        experiments = opts["experiments"]
        models      = opts["models"]
        experiment  = opts["experiment"]

        model = self._state.source_id if self._state.source_id in models else None
        if model is None:
            model = models[0] if models else None

        self._state.experiment_id = experiment
        self._state.source_id = model

        note = (
            "Catalogue unavailable — showing common models/scenarios."
            if opts["fallback"] else
            f"{len(models)} models available for {experiment}."
        )
        self._map.set_cmip_options(experiments, experiment, models, model or "",
                                   loading=False, note=note)
        self.log(f"CMIP6: {note}")

    def _on_cmip_options_error(self, tb: str) -> None:
        self.log("✗ Could not load the CMIP6 catalogue:")
        self.log(tb)
        self._map.set_cmip_options([], self._state.experiment_id or "", [], "",
                                   loading=False,
                                   note="Catalogue lookup failed — see log.")

    def _cleanup_opt_worker(self) -> None:
        if self._opt_worker is not None:
            self._opt_worker.deleteLater()
            self._opt_worker = None
        if self._opt_thread is not None:
            self._opt_thread.deleteLater()
            self._opt_thread = None
        if self._pending_experiment is not None:
            pending, self._pending_experiment = self._pending_experiment, None
            QTimer.singleShot(0, lambda e=pending: self._load_cmip_options(e))

    def _on_dates_changed(self, start_str: str, end_str: str) -> None:
        try:
            self._state.start_date = date.fromisoformat(start_str)
            self._state.end_date   = date.fromisoformat(end_str)
        except ValueError:
            pass

    def _on_format_changed(self, fmt: str) -> None:
        self._state.actions = ["extract", "to_nc"] if fmt == "NetCDF" else ["extract", "to_csv"]

    # ── Map interaction slots ─────────────────────────────────────────
    def _on_point_selected(self, lat: float, lon: float) -> None:
        self._state.lat = lat
        self._state.lon = lon
        self._state.box = None

    def _on_box_selected(self, lat_min: float, lat_max: float,
                         lon_min: float, lon_max: float) -> None:
        self._state.box = {
            "lat_min": lat_min, "lat_max": lat_max,
            "lon_min": lon_min, "lon_max": lon_max,
        }
        self._state.lat = None
        self._state.lon = None

    # ── Pipeline ──────────────────────────────────────────────────────
    def _on_run_clicked(self) -> None:
        if not self._state.is_ready():
            self.log("Cannot run: select a dataset, a point/box on the map, and a date range.")
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
            experiment_id=self._state.experiment_id,
            source_id=self._state.source_id,
        )
        self.log("▶ Starting pipeline — overrides:")
        for o in overrides:
            self.log(f"  {o}")
        self._start_worker(overrides)

    def _start_worker(self, overrides: list) -> None:
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

        self._map.send_run_state(True)
        self._thread.start()

    def _on_worker_finished(self, result) -> None:
        self.log("✓ Pipeline completed.")
        if result is not None and hasattr(result, "filename") and result.filename:
            self.log(f"Output: {result.filename}")
        self._last_result = result
        self._map.send_plot_enabled(result is not None)

    def _on_worker_error(self, tb: str) -> None:
        self.log("✗ Pipeline failed:")
        self.log(tb)

    def _cleanup_worker(self) -> None:
        if self._worker is not None:
            self._worker.deleteLater()
            self._worker = None
        if self._thread is not None:
            self._thread.deleteLater()
            self._thread = None
        self._map.send_run_state(False)

    # ── Overlay ───────────────────────────────────────────────────────
    def _on_plot_clicked(self) -> None:
        if self._last_result is None:
            self.log("No result — run the pipeline first.")
            return
        meta = get_overlay_meta(self._last_result)
        if meta is None:
            self.log("Cannot plot: dataset has no spatial dimensions.")
            return
        self.log(f"Initialising overlay: {', '.join(meta['vars'])}")
        self._map.init_overlay_ui(meta)

    def _on_render_requested(self, var_name: str, colormap: str,
                              vmin_str: str, vmax_str: str) -> None:
        if self._last_result is None or self._last_result.dataset is None:
            return
        vmin = None if vmin_str == 'auto' else float(vmin_str)
        vmax = None if vmax_str == 'auto' else float(vmax_str)
        self.log(f"Rendering '{var_name}' ({colormap})…")
        payload = render_variable(self._last_result.dataset, var_name, colormap,
                                  vmin=vmin, vmax=vmax)
        if payload:
            QTimer.singleShot(0, lambda p=payload: self._map.update_overlay(p))
            self.log(f"Overlay updated: {var_name} [{colormap}]")
        else:
            self.log(f"Render failed for '{var_name}' — no spatial dims.")

    # ── Advanced config dialog ────────────────────────────────────────
    def _on_advanced_clicked(self) -> None:
        current_overrides = build_overrides(
            dataset=self._state.dataset,
            lat=self._state.lat,
            lon=self._state.lon,
            start_date=self._state.start_date,
            end_date=self._state.end_date,
            variables=self._state.variables,
            experiment_id=self._state.experiment_id,
            source_id=self._state.source_id,
        )
        dlg = YamlEditorDialog(self, overrides=current_overrides)
        if dlg.exec() == YamlEditorDialog.DialogCode.Accepted:
            self._state.extra_overrides = dlg.extra_overrides()
            fields = dlg.parsed_fields()

            if "dataset" in fields:
                self._state.dataset = str(fields["dataset"])
            if "experiment_id" in fields:
                self._state.experiment_id = str(fields["experiment_id"])
            if "source_id" in fields:
                self._state.source_id = str(fields["source_id"])
            if "lat" in fields and "lon" in fields:
                self._state.lat = float(fields["lat"])
                self._state.lon = float(fields["lon"])
            if "start_date" in fields or "end_date" in fields:
                def _parse(s):
                    p = str(s).split("-")
                    return date(int(p[0]), int(p[1]), int(p[2]))
                try:
                    if "start_date" in fields:
                        self._state.start_date = _parse(fields["start_date"])
                    if "end_date" in fields:
                        self._state.end_date = _parse(fields["end_date"])
                except Exception:
                    pass

            self._map.sync_toolbar(
                dataset=self._state.dataset,
                start=self._state.start_date.isoformat() if self._state.start_date else "",
                end=self._state.end_date.isoformat() if self._state.end_date else "",
                data_dir=self._state.data_dir or "",
                experiment=self._state.experiment_id or "",
                model=self._state.source_id or "",
            )
            if is_cmip_dataset(self._state.dataset):
                self._load_cmip_options(self._state.experiment_id)
            self.log(f"Advanced config saved ({len(self._state.extra_overrides)} extra overrides).")

    # ── Shutdown ──────────────────────────────────────────────────────
    def closeEvent(self, event) -> None:
        """Let an in-flight catalogue lookup unwind before the window dies."""
        if self._opt_thread is not None:
            self._pending_experiment = None
            self._opt_thread.quit()
            self._opt_thread.wait(3000)
        super().closeEvent(event)

    # ── Data directory ────────────────────────────────────────────────
    def _on_data_dir_browse(self) -> None:
        # Defer out of QWebChannel callback before opening a blocking dialog
        QTimer.singleShot(0, self._do_data_dir_browse)

    def _do_data_dir_browse(self) -> None:
        chosen = QFileDialog.getExistingDirectory(
            self, "Select Data Directory", self._state.data_dir or "."
        )
        if chosen:
            self._state.data_dir = chosen
            self._map.send_data_dir(chosen)
            self.log(f"Data directory: {chosen}")

