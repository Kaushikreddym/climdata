"""QObject worker that runs the pipeline on a background thread."""

from __future__ import annotations

import traceback
from typing import Dict, List, Optional, Sequence

from PySide6.QtCore import QObject, Signal, Slot

from .basd import run_basd
from .comparison import run_comparison
from .runner import get_cmip_experiments, get_cmip_models, run_pipeline


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


class BasdWorker(QObject):
    """Run a bias-adjustment job off the GUI thread.

    Signals:
        finished(object): emitted with the :class:`~backend.basd.BasdResult`.
        error(str): emitted with a traceback string on failure.
        log(str): emitted for human-readable progress messages.
    """

    finished = Signal(object)
    error = Signal(str)
    log = Signal(str)

    def __init__(self, params: Dict, parent: Optional[QObject] = None) -> None:
        super().__init__(parent)
        self._params = dict(params)

    @Slot()
    def run(self) -> None:
        try:
            result = run_basd(log=self.log.emit, **self._params)
        except Exception:
            self.error.emit(traceback.format_exc())
            return
        self.finished.emit(result)


class ComparisonWorker(QObject):
    """Run a two-dataset comparison off the GUI thread.

    Signals:
        finished(object): emitted with the comparison result dict.
        error(str): emitted with a traceback string on failure.
        log(str): emitted for human-readable progress messages.
    """

    finished = Signal(object)
    error = Signal(str)
    log = Signal(str)

    def __init__(self, params: Dict, parent: Optional[QObject] = None) -> None:
        super().__init__(parent)
        self._params = dict(params)

    @Slot()
    def run(self) -> None:
        try:
            result = run_comparison(log=self.log.emit, **self._params)
        except Exception:
            self.error.emit(traceback.format_exc())
            return
        self.finished.emit(result)


class CmipOptionsWorker(QObject):
    """Fetch CMIP6 experiment / model lists off the GUI thread.

    The Pangeo catalogue is a network resource that takes several seconds on a
    cold cache, so the Download panel's model and experiment pickers are filled
    asynchronously.

    Signals:
        finished(object): emitted with
            ``{"experiments": [...], "experiment": str, "models": [...],
            "fallback": bool}``.
        error(str): emitted with a traceback string on failure.
        log(str): emitted for human-readable progress messages.
    """

    finished = Signal(object)
    error = Signal(str)
    log = Signal(str)

    def __init__(
        self,
        experiment: Optional[str] = None,
        parent: Optional[QObject] = None,
    ) -> None:
        super().__init__(parent)
        self._experiment = experiment

    @Slot()
    def run(self) -> None:
        try:
            experiments, exp_reason = get_cmip_experiments()
            experiment = self._experiment
            if experiment not in experiments:
                experiment = experiments[0] if experiments else "historical"
            models, mod_reason = get_cmip_models(experiment)
        except Exception:
            self.error.emit(traceback.format_exc())
            return
        reason = exp_reason or mod_reason
        if reason:
            self.log.emit(f"CMIP6 catalogue unreachable — {reason}")
        self.finished.emit({
            "experiments": experiments,
            "experiment":  experiment,
            "models":      models,
            "fallback":    bool(reason),
            "reason":      reason or "",
        })
