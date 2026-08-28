"""Top-level window — thin QMainWindow wrapper around the web dashboard."""

from __future__ import annotations

import json
import os
import sys
from datetime import date
from typing import Dict, List, Optional, Tuple

from PySide6.QtCore import QProcess, QThread, QTimer
from PySide6.QtWidgets import QFileDialog, QMainWindow, QMessageBox, QStatusBar

from backend.config_builder import build_overrides
from backend.runner import (
    CMIP_DATASETS,
    get_datasets,
    get_overlay_meta,
    is_cmip_dataset,
    render_variable,
)
from backend.worker import (
    BasdWorker,
    CmipOptionsWorker,
    ComparisonWorker,
    PipelineWorker,
)
from gui.controls.dataset_selector import DEFAULT_DATASETS
from models.app_state import AppState

from utils.certs import ca_bundle_status

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

        # CMIP6 catalogue lookup (model / experiment pickers, one slot at a time)
        self._opt_thread: Optional[QThread] = None
        self._opt_worker: Optional[CmipOptionsWorker] = None
        self._opt_slot: str = "download"
        self._opt_queue: List[Tuple[str, Optional[str]]] = []
        self._slot_datasets: Dict[str, str] = {}

        # BASD / Comparison jobs
        self._basd_thread: Optional[QThread] = None
        self._basd_worker: Optional[BasdWorker] = None
        self._cmp_thread: Optional[QThread] = None
        self._cmp_worker: Optional[ComparisonWorker] = None

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
        self._map.slot_dataset_changed.connect(self._on_slot_dataset_changed)
        self._map.experiment_changed.connect(self._on_experiment_changed)
        self._map.model_changed.connect(self._on_model_changed)
        self._map.basd_run_requested.connect(self._on_basd_run)
        self._map.comparison_run_requested.connect(self._on_comparison_run)
        self._map.restart_requested.connect(self._on_restart_requested)
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
        self.log(f"TLS {ca_bundle_status()}")

    # ── Toolbar slots (JS → Python) ───────────────────────────────────
    def _on_dataset_changed(self, dataset: str) -> None:
        previous = self._state.dataset
        self._state.dataset = dataset
        self._sync_slot_dataset("download", dataset, previous)

    def _on_slot_dataset_changed(self, slot: str, dataset: str) -> None:
        """A CMIP-capable dataset select outside the Download tab changed."""
        previous = self._slot_datasets.get(slot)
        self._slot_datasets[slot] = dataset
        self._sync_slot_dataset(slot, dataset, previous)

    def _sync_slot_dataset(self, slot: str, dataset: str,
                           previous: Optional[str]) -> None:
        if is_cmip_dataset(dataset):
            # Only (re)query when switching in — avoids a catalogue hit on every
            # echo of the current selection from the web frontend.
            _, model = self._state.get_cmip(slot)
            if not is_cmip_dataset(previous) or not model:
                experiment, _ = self._state.get_cmip(slot)
                self._load_cmip_options(slot, experiment)
        else:
            self._state.clear_cmip(slot)

    def _on_experiment_changed(self, slot: str, experiment_id: str) -> None:
        current, _ = self._state.get_cmip(slot)
        if not experiment_id or experiment_id == current:
            return
        self._state.set_cmip(slot, experiment_id=experiment_id)
        # The model list is experiment-dependent — refresh it.
        self._load_cmip_options(slot, experiment_id)

    def _on_model_changed(self, slot: str, source_id: str) -> None:
        if source_id:
            self._state.set_cmip(slot, source_id=source_id)

    # ── CMIP6 catalogue lookup ────────────────────────────────────────
    def _load_cmip_options(self, slot: str,
                           experiment: Optional[str] = None) -> None:
        """Fetch experiment/model lists for one slot on a background thread."""
        if self._opt_thread is not None:
            # A lookup is in flight — queue this one behind it.
            self._opt_queue = [q for q in self._opt_queue if q[0] != slot]
            self._opt_queue.append((slot, experiment))
            return
        self._opt_slot = slot
        self._map.set_cmip_options(slot, [], experiment or "", [], "",
                                   loading=True, note="Loading CMIP6 catalogue…")
        self.log("Querying CMIP6 catalogue for models and experiments…")

        self._opt_thread = QThread(self)
        # Only the Download tab knows its variables up front, so only its
        # picker can be narrowed to models that serve all of them. The BASD
        # and comparison tabs pick a single variable at run time; their lists
        # stay filtered by frequency alone rather than being over-restricted
        # to models carrying every Download variable.
        variables = self._state.variables if slot == "download" else None
        self._opt_worker = CmipOptionsWorker(experiment, variables)
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
        slot        = self._opt_slot
        experiments = opts["experiments"]
        models      = opts["models"]
        experiment  = opts["experiment"]

        _, current_model = self._state.get_cmip(slot)
        model = current_model if current_model in models else None
        if model is None:
            model = models[0] if models else None

        self._state.set_cmip(slot, experiment_id=experiment, source_id=model or "")

        if opts["fallback"]:
            note = ("Catalogue unreachable — showing common models/scenarios. "
                    + opts.get("reason", ""))
        elif slot == "download":
            vars_txt = ", ".join(self._state.variables)
            note = f"{len(models)} models provide daily {vars_txt} for {experiment}."
        else:
            note = f"{len(models)} models provide daily data for {experiment}."
        self._map.set_cmip_options(slot, experiments, experiment, models,
                                   model or "", loading=False, note=note)
        self.log(f"CMIP6 [{slot}]: {note}")

    def _on_cmip_options_error(self, tb: str) -> None:
        slot = self._opt_slot
        experiment, _ = self._state.get_cmip(slot)
        self.log("✗ Could not load the CMIP6 catalogue:")
        self.log(tb)
        self._map.set_cmip_options(slot, [], experiment or "", [], "",
                                   loading=False,
                                   note="Catalogue lookup failed — see log.")

    def _cleanup_opt_worker(self) -> None:
        if self._opt_worker is not None:
            self._opt_worker.deleteLater()
            self._opt_worker = None
        if self._opt_thread is not None:
            self._opt_thread.deleteLater()
            self._opt_thread = None
        if self._opt_queue:
            slot, experiment = self._opt_queue.pop(0)
            QTimer.singleShot(
                0, lambda s=slot, e=experiment: self._load_cmip_options(s, e))

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
                self._load_cmip_options("download", self._state.experiment_id)
            self.log(f"Advanced config saved ({len(self._state.extra_overrides)} extra overrides).")

    # ── BASD tab ──────────────────────────────────────────────────────
    def _on_basd_run(self, payload_json: str) -> None:
        if self._basd_thread is not None:
            self.log("BASD run already in progress.")
            return
        try:
            form = json.loads(payload_json)
        except ValueError as exc:
            self.log(f"✗ Could not read the BASD form: {exc}")
            return

        aoi = self._aoi_or_warning("BASD")
        if aoi is None:
            return
        try:
            params = dict(
                ref_dataset=form["ref_dataset"],
                ref_start=date.fromisoformat(form["ref_start"]),
                ref_end=date.fromisoformat(form["ref_end"]),
                target_dataset=form["target_dataset"],
                target_start=date.fromisoformat(form["target_start"]),
                target_end=date.fromisoformat(form["target_end"]),
                method=form["method"],
                variable=form["variable"],
                out_format=form.get("out_format", "NetCDF"),
                experiment_id=form.get("experiment_id") or None,
                source_id=form.get("source_id") or None,
                data_dir=self._state.data_dir,
                out_dir=self._state.data_dir,
                **aoi,
            )
        except (KeyError, ValueError) as exc:
            self.log(f"✗ Invalid BASD settings: {exc}")
            return

        self._basd_thread = QThread(self)
        self._basd_worker = BasdWorker(params)
        self._basd_worker.moveToThread(self._basd_thread)

        self._basd_thread.started.connect(self._basd_worker.run)
        self._basd_worker.log.connect(self.log)
        self._basd_worker.finished.connect(self._on_basd_finished)
        self._basd_worker.error.connect(self._on_basd_error)
        self._basd_worker.finished.connect(self._basd_thread.quit)
        self._basd_worker.error.connect(self._basd_thread.quit)
        self._basd_thread.finished.connect(self._cleanup_basd_worker)

        self._map.send_basd_run_state(True)
        self.log("▶ Starting bias adjustment…")
        self._basd_thread.start()

    def _on_basd_finished(self, result) -> None:
        self.log("✓ BASD completed.")
        if result.filename:
            self.log(f"Output: {result.filename}")
        self._last_result = result           # enables the Plot button's overlay
        self._map.send_plot_enabled(result.dataset is not None)
        self._map.set_basd_result({
            "method":   result.method,
            "variable": result.variable,
            "filename": result.filename,
            "metrics":  result.metrics,
            "notes":    result.notes,
        })

    def _on_basd_error(self, tb: str) -> None:
        self.log("✗ BASD failed:")
        self.log(tb)

    def _cleanup_basd_worker(self) -> None:
        if self._basd_worker is not None:
            self._basd_worker.deleteLater()
            self._basd_worker = None
        if self._basd_thread is not None:
            self._basd_thread.deleteLater()
            self._basd_thread = None
        self._map.send_basd_run_state(False)

    # ── Comparison tab ────────────────────────────────────────────────
    def _on_comparison_run(self, payload_json: str) -> None:
        if self._cmp_thread is not None:
            self.log("Comparison already in progress.")
            return
        try:
            form = json.loads(payload_json)
        except ValueError as exc:
            self.log(f"✗ Could not read the comparison form: {exc}")
            return

        aoi = self._aoi_or_warning("Comparison")
        if aoi is None:
            return
        try:
            params = dict(
                a_dataset=form["a_dataset"], a_variable=form["a_variable"],
                a_start=date.fromisoformat(form["a_start"]),
                a_end=date.fromisoformat(form["a_end"]),
                b_dataset=form["b_dataset"], b_variable=form["b_variable"],
                b_start=date.fromisoformat(form["b_start"]),
                b_end=date.fromisoformat(form["b_end"]),
                a_experiment=form.get("a_experiment") or None,
                a_source=form.get("a_source") or None,
                b_experiment=form.get("b_experiment") or None,
                b_source=form.get("b_source") or None,
                data_dir=self._state.data_dir,
                **aoi,
            )
        except (KeyError, ValueError) as exc:
            self.log(f"✗ Invalid comparison settings: {exc}")
            return

        self._cmp_thread = QThread(self)
        self._cmp_worker = ComparisonWorker(params)
        self._cmp_worker.moveToThread(self._cmp_thread)

        self._cmp_thread.started.connect(self._cmp_worker.run)
        self._cmp_worker.log.connect(self.log)
        self._cmp_worker.finished.connect(self._on_comparison_finished)
        self._cmp_worker.error.connect(self._on_comparison_error)
        self._cmp_worker.finished.connect(self._cmp_thread.quit)
        self._cmp_worker.error.connect(self._cmp_thread.quit)
        self._cmp_thread.finished.connect(self._cleanup_cmp_worker)

        self._map.send_comparison_run_state(True)
        self.log("▶ Starting comparison…")
        self._cmp_thread.start()

    def _on_comparison_finished(self, result: dict) -> None:
        self.log("✓ Comparison completed.")
        self._map.set_comparison_result(result)

    def _on_comparison_error(self, tb: str) -> None:
        self.log("✗ Comparison failed:")
        self.log(tb)

    def _cleanup_cmp_worker(self) -> None:
        if self._cmp_worker is not None:
            self._cmp_worker.deleteLater()
            self._cmp_worker = None
        if self._cmp_thread is not None:
            self._cmp_thread.deleteLater()
            self._cmp_thread = None
        self._map.send_comparison_run_state(False)

    # ── Shared helpers ────────────────────────────────────────────────
    def _aoi_or_warning(self, what: str) -> Optional[dict]:
        """Return the current AOI as kwargs, or None (with a log line) if unset."""
        if self._state.box is not None:
            return {"lat": None, "lon": None, "box": self._state.box}
        if self._state.lat is not None and self._state.lon is not None:
            return {"lat": self._state.lat, "lon": self._state.lon, "box": None}
        self.log(f"Cannot run {what}: select a point or bounding box on the map first.")
        return None

    # ── Restart ───────────────────────────────────────────────────────
    def _on_restart_requested(self) -> None:
        # Defer out of the QWebChannel callback before opening a modal dialog
        QTimer.singleShot(0, self._confirm_restart)

    def _busy_jobs(self) -> List[str]:
        """Names of the background jobs currently running."""
        running = []
        if self._thread is not None:
            running.append("download pipeline")
        if self._basd_thread is not None:
            running.append("BASD run")
        if self._cmp_thread is not None:
            running.append("comparison")
        if self._opt_thread is not None:
            running.append("CMIP6 catalogue lookup")
        return running

    def _confirm_restart(self) -> None:
        busy = self._busy_jobs()
        detail = (f"This will interrupt: {', '.join(busy)}.\n\n"
                  if busy else "")
        answer = QMessageBox.question(
            self,
            "Restart climdata",
            f"Restart the application?\n\n{detail}"
            "Loaded datasets, results and the log are discarded. "
            "Files already written to disk are kept.",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        if answer != QMessageBox.StandardButton.Yes:
            self.log("Restart cancelled.")
            return
        self._restart_application()

    def _restart_application(self) -> None:
        """Start a fresh instance and kill this one — jobs and all.

        A climdata extraction blocks inside network/IO calls that no Qt thread
        can be asked to abandon, so an in-process "interrupt" cannot be honest.
        Replacing the process is the only way to guarantee the running work
        actually stops — the same bargain a Jupyter kernel restart makes.
        """
        self._map.send_restarting(True)
        self.log("⟳ Restarting — interrupting all running jobs…")

        if getattr(sys, "frozen", False):     # PyInstaller bundle
            program, args = sys.executable, sys.argv[1:]
        else:
            program, args = sys.executable, list(sys.argv)

        started = QProcess.startDetached(program, args, os.getcwd())
        if not started:
            self.log("✗ Could not start a new instance — staying open.")
            self._map.send_restarting(False)
            return

        # Give the frontend a moment to paint the "Restarting…" state, then go.
        QTimer.singleShot(150, self._exit_now)

    def _exit_now(self) -> None:
        sys.stdout.flush()
        sys.stderr.flush()
        # os._exit skips interpreter teardown on purpose: a worker thread parked
        # in a socket read would otherwise keep this process alive.
        os._exit(0)

    # ── Shutdown ──────────────────────────────────────────────────────
    def closeEvent(self, event) -> None:
        """Let in-flight background jobs unwind before the window dies."""
        self._opt_queue.clear()
        for thread in (self._opt_thread, self._basd_thread, self._cmp_thread):
            if thread is not None:
                thread.quit()
                thread.wait(3000)
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

