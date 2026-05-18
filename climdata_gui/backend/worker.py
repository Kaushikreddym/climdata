"""QObject worker that runs the pipeline on a background thread."""

from __future__ import annotations

import traceback
from typing import List, Optional, Sequence

from PySide6.QtCore import QObject, Signal, Slot

from .runner import run_pipeline


class PipelineWorker(QObject):
    """Run ``run_pipeline`` off the GUI thread.

    Signals:
        started: emitted when the run begins.
        finished(object): emitted with the ``WorkflowResult`` on success.
        error(str): emitted with a traceback string on failure.
        log(str): emitted for human-readable progress messages.
    """

    started = Signal()
    finished = Signal(object)
    error = Signal(str)
    log = Signal(str)

    def __init__(
        self,
        overrides: Sequence[str],
        seq: Optional[Sequence[str]] = None,
        parent: Optional[QObject] = None,
    ) -> None:
        super().__init__(parent)
        self._overrides: List[str] = list(overrides)
        self._seq: Optional[List[str]] = list(seq) if seq else None

    @Slot()
    def run(self) -> None:
        self.started.emit()
        self.log.emit(f"Starting pipeline with overrides: {self._overrides}")
        try:
            result = run_pipeline(self._overrides, self._seq)
        except Exception:
            self.error.emit(traceback.format_exc())
            return
        self.log.emit("Pipeline finished.")
        self.finished.emit(result)
