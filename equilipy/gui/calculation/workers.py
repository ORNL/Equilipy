"""Thread worker helpers for GUI calculations."""

from __future__ import annotations

from threading import Event

from PySide6.QtCore import QObject, Signal, Slot


class CalculationCanceled(RuntimeError):
    """Raised when the user requests cancellation through the progress dialog."""


class CalculationWorker(QObject):
    """Run one calculation callable away from the UI thread."""

    finished = Signal(object, object)
    progress = Signal(int, int, str)

    def __init__(self, work):
        super().__init__()
        self._work = work
        self._cancel_event = Event()

    @Slot()
    def cancel(self) -> None:
        """Request cancellation at the next calculation progress checkpoint."""
        self._cancel_event.set()

    @Slot()
    def run(self) -> None:
        """Execute the calculation and emit either result or exception."""
        try:
            self.finished.emit(self._work(self._emit_progress), None)
        except Exception as exc:  # pragma: no cover - exercised through GUI.
            self.finished.emit(None, exc)

    def _emit_progress(self, current: int, total: int, message: str = "") -> None:
        """Emit one progress update from the worker thread."""
        if self._cancel_event.is_set():
            raise CalculationCanceled("Calculation canceled.")
        self.progress.emit(int(current), int(total), str(message))
        if self._cancel_event.is_set():
            raise CalculationCanceled("Calculation canceled.")


class CalculationResultReceiver(QObject):
    """Main-thread receiver for calculation worker results."""

    finished = Signal(object, object)
    progress = Signal(int, int, str)

    @Slot(object, object)
    def receive(self, result, error) -> None:
        """Forward worker results through a main-thread signal."""
        self.finished.emit(result, error)

    @Slot(int, int, str)
    def receive_progress(self, current: int, total: int, message: str) -> None:
        """Forward worker progress through a main-thread signal."""
        self.progress.emit(current, total, message)
