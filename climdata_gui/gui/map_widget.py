"""Leaflet-backed map widget — thin QWebEngineView container for the web dashboard."""

from __future__ import annotations

import json

from PySide6.QtCore import QObject, QTimer, QUrl, Qt, Signal, Slot
from PySide6.QtWebChannel import QWebChannel
from PySide6.QtWebEngineCore import QWebEnginePage, QWebEngineSettings
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWidgets import QVBoxLayout, QWidget

from utils.paths import resource_path


class _DebugPage(QWebEnginePage):
    """Forwards JS console messages to Python stdout."""

    def javaScriptConsoleMessage(self, level, message, line_number, source_id):
        level_str = {
            QWebEnginePage.JavaScriptConsoleMessageLevel.InfoMessageLevel:    "INFO",
            QWebEnginePage.JavaScriptConsoleMessageLevel.WarningMessageLevel: "WARN",
            QWebEnginePage.JavaScriptConsoleMessageLevel.ErrorMessageLevel:   "ERROR",
        }.get(level, "LOG")
        print(f"[JS {level_str}] {message}  ({source_id}:{line_number})", flush=True)


class _Bridge(QObject):
    """QObject registered with QWebChannel — all signals/slots are exposed to JS."""

    # ── Python → JS signals (QWebChannel transmits these to JS) ──────
    log_message       = Signal(str)
    run_state_changed = Signal(bool)   # True = running
    plot_enabled      = Signal(bool)   # True = result ready
    data_dir_chosen   = Signal(str)    # path from QFileDialog

    # ── Relay signals (bubble up through MapWidget to MainWindow) ─────
    point_selected            = Signal(float, float)
    box_selected              = Signal(float, float, float, float)
    render_requested          = Signal(str, str, str, str)
    dataset_changed           = Signal(str)
    dates_changed             = Signal(str, str)   # ISO-date strings
    format_changed            = Signal(str)
    experiment_changed        = Signal(str)        # CMIP6 experiment_id
    model_changed             = Signal(str)        # CMIP6 source_id
    run_clicked               = Signal()
    plot_clicked              = Signal()
    advanced_clicked          = Signal()
    data_dir_browse_requested = Signal()
    page_ready                = Signal()

    # ── Slots: called by JS via QWebChannel ───────────────────────────
    @Slot(float, float)
    def on_point_selected(self, lat: float, lon: float) -> None:
        self.point_selected.emit(lat, lon)

    @Slot(float, float, float, float)
    def on_box_selected(self, lat_min: float, lat_max: float,
                        lon_min: float, lon_max: float) -> None:
        self.box_selected.emit(lat_min, lat_max, lon_min, lon_max)

    @Slot(str, str, str, str)
    def on_render_request(self, var_name: str, colormap: str,
                          vmin_str: str, vmax_str: str) -> None:
        self.render_requested.emit(var_name, colormap, vmin_str, vmax_str)

    @Slot(str)
    def on_dataset_changed(self, dataset: str) -> None:
        self.dataset_changed.emit(dataset)

    @Slot(str, str)
    def on_dates_changed(self, start: str, end: str) -> None:
        self.dates_changed.emit(start, end)

    @Slot(str)
    def on_format_changed(self, fmt: str) -> None:
        self.format_changed.emit(fmt)

    @Slot(str)
    def on_experiment_changed(self, experiment_id: str) -> None:
        self.experiment_changed.emit(experiment_id)

    @Slot(str)
    def on_model_changed(self, source_id: str) -> None:
        self.model_changed.emit(source_id)

    @Slot()
    def on_run_clicked(self) -> None:
        self.run_clicked.emit()

    @Slot()
    def on_plot_clicked(self) -> None:
        self.plot_clicked.emit()

    @Slot()
    def on_advanced_clicked(self) -> None:
        self.advanced_clicked.emit()

    @Slot()
    def on_data_dir_browse(self) -> None:
        self.data_dir_browse_requested.emit()

    @Slot()
    def on_page_ready(self) -> None:
        # Defer out of the QWebChannel callback to avoid re-entrancy
        QTimer.singleShot(0, self.page_ready.emit)


class MapWidget(QWidget):
    """Full-screen web dashboard container with bidirectional Qt ↔ JS bridge."""

    # Surface all bridge signals for MainWindow wiring
    point_selected            = Signal(float, float)
    box_selected              = Signal(float, float, float, float)
    render_requested          = Signal(str, str, str, str)
    dataset_changed           = Signal(str)
    dates_changed             = Signal(str, str)
    format_changed            = Signal(str)
    experiment_changed        = Signal(str)
    model_changed             = Signal(str)
    run_clicked               = Signal()
    plot_clicked              = Signal()
    advanced_clicked          = Signal()
    data_dir_browse_requested = Signal()
    page_ready                = Signal()

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self._view = QWebEngineView(self)
        self._page = _DebugPage(self._view)
        self._view.setPage(self._page)
        layout.addWidget(self._view)

        self._channel = QWebChannel(self._page)
        self._bridge = _Bridge()

        # Relay every bridge signal up to MapWidget
        for sig_name in (
            "point_selected", "box_selected", "render_requested",
            "dataset_changed", "dates_changed", "format_changed",
            "experiment_changed", "model_changed",
            "run_clicked", "plot_clicked", "advanced_clicked",
            "data_dir_browse_requested", "page_ready",
        ):
            getattr(self._bridge, sig_name).connect(
                getattr(self, sig_name).emit
            )

        self._channel.registerObject("bridge", self._bridge)
        self._page.setWebChannel(self._channel)

        settings = self._page.settings()
        settings.setAttribute(QWebEngineSettings.WebAttribute.LocalContentCanAccessRemoteUrls, True)
        settings.setAttribute(QWebEngineSettings.WebAttribute.LocalContentCanAccessFileUrls, True)
        settings.setAttribute(QWebEngineSettings.WebAttribute.TouchEventsApiEnabled, True)
        self._view.setAttribute(Qt.WidgetAttribute.WA_AcceptTouchEvents, True)

        html_path = resource_path("resources/html/index.html")
        self._view.load(QUrl.fromLocalFile(html_path))

    # ── Python → JS helpers ───────────────────────────────────────────

    def init_overlay_ui(self, config: dict) -> None:
        self._page.runJavaScript(f"initOverlayUI({json.dumps(config)});")

    def update_overlay(self, payload: dict) -> None:
        self._page.runJavaScript(f"updateOverlay({json.dumps(payload)});")

    def clear_overlay(self) -> None:
        self._page.runJavaScript("clearOverlays();")

    def send_log(self, msg: str) -> None:
        self._bridge.log_message.emit(msg)

    def send_run_state(self, running: bool) -> None:
        self._bridge.run_state_changed.emit(running)

    def send_plot_enabled(self, enabled: bool) -> None:
        self._bridge.plot_enabled.emit(enabled)

    def send_data_dir(self, path: str) -> None:
        self._bridge.data_dir_chosen.emit(path)

    def dashboard_ready(self, datasets: list, dataset: str,
                        start: str, end: str, data_dir: str,
                        cmip_datasets: list = None) -> None:
        """Push initial toolbar state to the web frontend after page_ready."""
        payload = {
            "datasets": datasets,
            "dataset":  dataset,
            "start":    start,
            "end":      end,
            "data_dir": data_dir,
            "cmip_datasets": list(cmip_datasets or []),
        }
        self._page.runJavaScript(f"onDashboardReady({json.dumps(payload)});")

    def sync_toolbar(self, dataset: str, start: str, end: str, data_dir: str,
                     experiment: str = None, model: str = None) -> None:
        """Sync toolbar controls to current Python state (e.g. after advanced dialog)."""
        payload = {"dataset": dataset, "start": start, "end": end, "data_dir": data_dir}
        if experiment is not None:
            payload["experiment"] = experiment
        if model is not None:
            payload["model"] = model
        self._page.runJavaScript(f"syncToolbar({json.dumps(payload)});")

    def set_cmip_options(self, experiments: list, experiment: str,
                         models: list, model: str,
                         loading: bool = False, note: str = "") -> None:
        """Fill (or put in a loading state) the CMIP6 experiment/model pickers."""
        payload = {
            "experiments": list(experiments),
            "experiment":  experiment or "",
            "models":      list(models),
            "model":       model or "",
            "loading":     loading,
            "note":        note,
        }
        self._page.runJavaScript(f"setCmipOptions({json.dumps(payload)});")

