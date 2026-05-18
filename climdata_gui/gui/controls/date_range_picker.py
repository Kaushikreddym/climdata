"""Two ``QDateEdit`` widgets for selecting a start/end date range."""

from __future__ import annotations

from datetime import date

from PySide6.QtCore import QDate, Signal
from PySide6.QtWidgets import QDateEdit, QHBoxLayout, QLabel, QWidget


def _to_qdate(d: date) -> QDate:
    return QDate(d.year, d.month, d.day)


class DateRangePicker(QWidget):
    """Start/end date picker."""

    range_changed = Signal(object, object)  # (date, date)

    def __init__(
        self,
        start: date = date(1989, 1, 1),
        end: date = date(2020, 12, 31),
        parent=None,
    ) -> None:
        super().__init__(parent)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        layout.addWidget(QLabel("Start:"))
        self._start = QDateEdit(_to_qdate(start))
        self._start.setCalendarPopup(True)
        self._start.setDisplayFormat("yyyy-MM-dd")
        layout.addWidget(self._start)

        layout.addWidget(QLabel("End:"))
        self._end = QDateEdit(_to_qdate(end))
        self._end.setCalendarPopup(True)
        self._end.setDisplayFormat("yyyy-MM-dd")
        layout.addWidget(self._end)

        self._start.dateChanged.connect(self._emit)
        self._end.dateChanged.connect(self._emit)

    def _emit(self, *_):
        self.range_changed.emit(self.start_date(), self.end_date())

    def start_date(self) -> date:
        d = self._start.date()
        return date(d.year(), d.month(), d.day())

    def end_date(self) -> date:
        d = self._end.date()
        return date(d.year(), d.month(), d.day())

    def set_dates(self, start: date, end: date) -> None:
        self._start.setDate(QDate(start.year, start.month, start.day))
        self._end.setDate(QDate(end.year, end.month, end.day))
