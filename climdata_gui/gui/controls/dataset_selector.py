"""Dropdown of available climdata providers."""

from __future__ import annotations

from typing import List

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QComboBox, QHBoxLayout, QLabel, QWidget

DEFAULT_DATASETS: List[str] = ["mswx", "dwd", "hyras", "cmip", "power", "nexgddp"]


class DatasetSelector(QWidget):
    """Labeled combo box for selecting the dataset provider."""

    dataset_changed = Signal(str)

    def __init__(self, datasets: List[str] = None, parent=None) -> None:
        super().__init__(parent)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        layout.addWidget(QLabel("Dataset:"))
        self._combo = QComboBox()
        self._combo.addItems(datasets or DEFAULT_DATASETS)
        self._combo.currentTextChanged.connect(self.dataset_changed.emit)
        layout.addWidget(self._combo, 1)

    def current_dataset(self) -> str:
        return self._combo.currentText()
