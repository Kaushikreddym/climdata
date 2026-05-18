"""Leaflet-backed map widget that emits ``(lat, lon)`` on click."""

from __future__ import annotations

from PySide6.QtCore import QObject, QUrl, Qt, Signal, Slot
from PySide6.QtWebChannel import QWebChannel
from PySide6.QtWebEngineCore import QWebEnginePage, QWebEngineSettings
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWidgets import QVBoxLayout, QWidget

from utils.paths import resource_path


class _DebugPage(QWebEnginePage):
    """QWebEnginePage that forwards JS console messages to Python stdout."""

    def javaScriptConsoleMessage(self, level, message, line_number, source_id):
        level_str = {
            QWebEnginePage.JavaScriptConsoleMessageLevel.InfoMessageLevel:    "INFO",
            QWebEnginePage.JavaScriptConsoleMessageLevel.WarningMessageLevel: "WARN",
            QWebEnginePage.JavaScriptConsoleMessageLevel.ErrorMessageLevel:   "ERROR",
        }.get(level, "LOG")
        print(f"[JS {level_str}] {message}  ({source_id}:{line_number})", flush=True)


class _Bridge(QObject):
    """Object exposed to JavaScript via QWebChannel."""

    point_selected   = Signal(float, float)
    box_selected     = Signal(float, float, float, float)  # lat_min, lat_max, lon_min, lon_max
    render_requested = Signal(str, str, str, str)          # var_name, colormap, vmin_str, vmax_str

    @Slot(float, float)
    def on_point_selected(self, lat: float, lon: float) -> None:
        self.point_selected.emit(lat, lon)

    @Slot(float, float, float, float)
    def on_box_selected(self, lat_min: float, lat_max: float, lon_min: float, lon_max: float) -> None:
        self.box_selected.emit(lat_min, lat_max, lon_min, lon_max)

    @Slot(str, str, str, str)
    def on_render_request(self, var_name: str, colormap: str, vmin_str: str, vmax_str: str) -> None:
        self.render_requested.emit(var_name, colormap, vmin_str, vmax_str)


class MapWidget(QWidget):
    """Embedded Leaflet map with Qt <-> JS bridge."""

    point_selected   = Signal(float, float)
    box_selected     = Signal(float, float, float, float)  # lat_min, lat_max, lon_min, lon_max
    render_requested = Signal(str, str, str, str)          # var_name, colormap, vmin_str, vmax_str

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self._view = QWebEngineView(self)
        # Use a custom page that mirrors JS console output to Python stdout.
        self._page = _DebugPage(self._view)
        self._view.setPage(self._page)
        layout.addWidget(self._view)

        self._channel = QWebChannel(self._page)
        self._bridge = _Bridge()
        self._bridge.point_selected.connect(self.point_selected.emit)
        self._bridge.box_selected.connect(self.box_selected.emit)
        self._bridge.render_requested.connect(self.render_requested.emit)
        self._channel.registerObject("bridge", self._bridge)
        self._page.setWebChannel(self._channel)

        settings = self._page.settings()
        settings.setAttribute(QWebEngineSettings.WebAttribute.LocalContentCanAccessRemoteUrls, True)
        settings.setAttribute(QWebEngineSettings.WebAttribute.LocalContentCanAccessFileUrls, True)
        settings.setAttribute(QWebEngineSettings.WebAttribute.TouchEventsApiEnabled, True)

        self._view.setAttribute(Qt.WidgetAttribute.WA_AcceptTouchEvents, True)

        html_path = resource_path("resources/html/map.html")
        self._view.load(QUrl.fromLocalFile(html_path))

    def init_overlay_ui(self, config: dict) -> None:
        import json
        js = f"initOverlayUI({json.dumps(config)});"
        print(f"[PY] init_overlay_ui: vars={config.get('vars')} bounds=({config.get('lat_min')},{config.get('lon_min')})->({config.get('lat_max')},{config.get('lon_max')})", flush=True)
        self._page.runJavaScript(js)

    def update_overlay(self, payload: dict) -> None:
        import json
        b64_len = len(payload.get("b64", ""))
        print(f"[PY] update_overlay: b64_len={b64_len}", flush=True)
        self._page.runJavaScript(f"updateOverlay({json.dumps(payload)});")

    def clear_overlay(self) -> None:
        self._page.runJavaScript("clearOverlays();")
