"""Log display widgets for the Equilipy GUI."""

from __future__ import annotations

from PySide6.QtGui import QTextCursor
from PySide6.QtWidgets import QPlainTextEdit


class CalculationLog(QPlainTextEdit):
    """Append-style log that also remembers the latest status message."""

    def __init__(self):
        super().__init__()
        self.last_status = ""
        self.status_sink = None

    def setPlainText(self, text: str) -> None:
        """Append status text so the bottom log keeps calculation history."""
        self.last_status = text
        if not text:
            super().setPlainText(text)
            return
        if self.status_sink is not None:
            self.status_sink(text)
        if self.toPlainText().strip():
            self.appendPlainText(text)
        else:
            super().setPlainText(text)

    def append_stream_text(self, text: str) -> None:
        """Append raw stdout/stderr text without adding extra newlines."""
        if not text:
            return
        self.moveCursor(QTextCursor.MoveOperation.End)
        self.insertPlainText(text)
        self.moveCursor(QTextCursor.MoveOperation.End)
