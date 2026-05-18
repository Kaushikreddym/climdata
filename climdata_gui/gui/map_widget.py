"""Leaflet-backed map widget that emits ``(lat, lon)`` on click."""

from __future__ import annotations

from PySide6.QtCore import QObject, QUrl, Signal, Slot
from PySide6.QtWebChannel import QWebChannel
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWidgets import QVBoxLayout, QWidget

from utils.paths import resource_path


class _Bridge(QObject):
    """Object exposed to JavaScript via QWebChannel."""

    point_selected = Signal(float, float)

    @Slot(float, float)
    def on_point_selected(self, lat: float, lon: float) -> None:
        self.point_selected.emit(lat, lon)


class MapWidget(QWidget):
    """Embedded Leaflet map with Qt <-> JS bridge."""

    point_selected = Signal(float, float)

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self._view = QWebEngineView(self)
        layout.addWidget(self._view)

        self._channel = QWebChannel(self._view.page())
        self._bridge = _Bridge()
        self._bridge.point_selected.connect(self.point_selected.emit)
        self._channel.registerObject("bridge", self._bridge)
        self._view.page().setWebChannel(self._channel)

        html_path = resource_path("resources/html/map.html")
        self._view.load(QUrl.fromLocalFile(html_path))
