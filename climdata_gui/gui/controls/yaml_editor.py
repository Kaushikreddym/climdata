"""YAML config editor dialog with syntax highlighting and add-field support."""

from __future__ import annotations

from typing import Dict, List, Optional

import yaml
from PySide6.QtCore import QRegularExpression, Qt
from PySide6.QtGui import QColor, QFont, QSyntaxHighlighter, QTextCharFormat
from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QHBoxLayout,
    QInputDialog,
    QLabel,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QVBoxLayout,
)


# ---------------------------------------------------------------------------
# Syntax highlighter
# ---------------------------------------------------------------------------

class _YamlHighlighter(QSyntaxHighlighter):
    """Minimal YAML syntax highlighter."""

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._rules: list = []

        def rule(pattern: str, color: str, bold: bool = False) -> None:
            fmt = QTextCharFormat()
            fmt.setForeground(QColor(color))
            if bold:
                fmt.setFontWeight(700)
            self._rules.append((QRegularExpression(pattern), fmt))

        rule(r"#[^\n]*", "#6272a4")                          # comments  (slate-blue)
        rule(r"(?m)^\s*[\w\.\-]+(?=\s*:)", "#79b8ff", True)  # keys       (bright blue)
        rule(r"(?<=:\s)['\"]?[^#\n]+", "#9ecbff")             # values     (sky blue)
        rule(r"(?m)^\s*-\s", "#f97583")                       # list marks (rose)
        rule(r"(?m)^---", "#b392f0", True)                    # separator  (purple)

    def highlightBlock(self, text: str) -> None:
        for pattern, fmt in self._rules:
            it = pattern.globalMatch(text)
            while it.hasNext():
                m = it.next()
                self.setFormat(m.capturedStart(), m.capturedLength(), fmt)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_config_yaml(overrides: List[str] | None = None) -> str:
    """Return the effective climdata config as a YAML string."""
    try:
        from climdata import ClimData
        from omegaconf import OmegaConf
        cfg = ClimData(overrides=overrides or []).cfg
        return OmegaConf.to_yaml(cfg)
    except Exception as exc:
        return f"# Could not load config: {exc}\n"


def _flatten_to_overrides(d: dict, prefix: str = "") -> List[str]:
    """Recursively flatten a dict to Hydra override strings.
    Lists containing dicts or other complex types are skipped.
    """
    items: List[str] = []
    for k, v in d.items():
        key = f"{prefix}.{k}" if prefix else str(k)
        if isinstance(v, dict):
            items.extend(_flatten_to_overrides(v, key))
        elif isinstance(v, list):
            # Only emit lists of plain scalars — Hydra cannot parse lists of dicts
            if v and all(isinstance(x, (str, int, float, bool)) for x in v):
                vals = "[" + ",".join(str(x) for x in v) + "]"
                items.append(f"{key}={vals}")
        elif isinstance(v, str):
            items.append(f'{key}="{v}"')
        elif v is None:
            pass
        else:
            items.append(f"{key}={v}")
    return items


# ---------------------------------------------------------------------------
# Dialog
# ---------------------------------------------------------------------------

class YamlEditorDialog(QDialog):
    """Modal YAML config editor.

    Call :meth:`extra_overrides` after ``exec()`` returns ``Accepted`` to
    retrieve the Hydra override list derived from the edited YAML.
    """

    def __init__(self, parent=None, overrides: List[str] | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Advanced Config Editor")
        self.resize(720, 560)

        layout = QVBoxLayout(self)

        # Header
        header = QLabel(
            "Edit the config below. Changes are applied as Hydra overrides when you click OK."
        )
        header.setWordWrap(True)
        layout.addWidget(header)

        # Editor
        self._editor = QPlainTextEdit()
        font = QFont("Menlo, Courier New, monospace")
        font.setPointSize(11)
        self._editor.setFont(font)
        self._editor.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        self._highlighter = _YamlHighlighter(self._editor.document())
        layout.addWidget(self._editor, 1)

        # Toolbar row
        toolbar = QHBoxLayout()
        add_btn = QPushButton("＋ Add Field")
        add_btn.setToolTip("Append a new key: value pair at the end")
        add_btn.clicked.connect(self._add_field)
        toolbar.addWidget(add_btn)
        toolbar.addStretch()
        layout.addLayout(toolbar)

        # OK / Cancel
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Save | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self._on_ok)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

        # Load config with current GUI state
        self._original_yaml = _load_config_yaml(overrides)
        self._editor.setPlainText(self._original_yaml)
        self._overrides: List[str] = []
        self._parsed_fields: Dict = {}

    # ------------------------------------------------------------------
    def extra_overrides(self) -> List[str]:
        """Return the Hydra override list from the last accepted edit."""
        return self._overrides

    def parsed_fields(self) -> Dict:
        """Return a dict with any of dataset/lat/lon/start_date/end_date found in YAML."""
        return self._parsed_fields

    # ------------------------------------------------------------------
    def _add_field(self) -> None:
        key, ok = QInputDialog.getText(self, "Add Field", "Key (dot-notation, e.g. time_range.start_date):")
        if not ok or not key.strip():
            return
        value, ok = QInputDialog.getText(self, "Add Field", f"Value for  '{key.strip()}':")
        if not ok:
            return
        line = f"{key.strip()}: {value.strip()}\n"
        cursor = self._editor.textCursor()
        self._editor.moveCursor(cursor.MoveOperation.End)
        self._editor.insertPlainText(line)

    def _on_ok(self) -> None:
        text = self._editor.toPlainText()
        try:
            data = yaml.safe_load(text)
        except yaml.YAMLError as exc:
            QMessageBox.critical(self, "YAML Parse Error", str(exc))
            return
        if not isinstance(data, dict):
            QMessageBox.warning(self, "Invalid Config", "Top-level YAML must be a mapping.")
            return
        self._overrides = _flatten_to_overrides(data)
        # Extract well-known GUI fields so main window can sync controls
        fields: Dict = {}
        if "dataset" in data:
            fields["dataset"] = str(data["dataset"])
        if "lat" in data:
            try:
                fields["lat"] = float(data["lat"])
            except (TypeError, ValueError):
                pass
        if "lon" in data:
            try:
                fields["lon"] = float(data["lon"])
            except (TypeError, ValueError):
                pass
        if "experiment_id" in data and data["experiment_id"] is not None:
            fields["experiment_id"] = str(data["experiment_id"])
        if "source_id" in data and data["source_id"] is not None:
            fields["source_id"] = str(data["source_id"])
        tr = data.get("time_range", {})
        if isinstance(tr, dict):
            if "start_date" in tr:
                fields["start_date"] = str(tr["start_date"])
            if "end_date" in tr:
                fields["end_date"] = str(tr["end_date"])
        self._parsed_fields = fields
        self.accept()
