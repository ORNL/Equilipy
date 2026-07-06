"""Bottom status bar widget for the Equilipy GUI."""

from __future__ import annotations

from PySide6.QtCore import QSize, Qt
from PySide6.QtSvgWidgets import QSvgWidget
from PySide6.QtWidgets import QFrame, QHBoxLayout, QLabel, QSizePolicy

from equilipy.gui.assets import ornl_logo_path
from equilipy.gui.models import function_error_icon, function_warning_icon


def status_level_from_message(message: str) -> str:
    """Infer a compact status severity from a user-facing message."""
    lowered = message.lower()
    if any(token in lowered for token in ("failed", "error", "invalid")):
        return "error"
    if any(
        token in lowered
        for token in ("warning", "issue", "diagnostic", "canceled", "dropped")
    ):
        return "warning"
    if any(
        token in lowered
        for token in ("running", "preparing", "finalizing", "calculating")
    ):
        return "running"
    return "ready"


class AppStatusBar(QFrame):
    """Application status bar with severity icon and ORNL branding."""

    def __init__(self) -> None:
        super().__init__()
        self.setObjectName("AppStatusBar")
        self.setFixedHeight(28)

        layout = QHBoxLayout(self)
        layout.setContentsMargins(12, 0, 12, 0)
        layout.setSpacing(8)

        self.status_symbol = QLabel()
        self.status_symbol.setObjectName("AppStatusSymbol")
        self.status_symbol.setFixedSize(9, 9)
        self.status_symbol.setAlignment(Qt.AlignmentFlag.AlignCenter)

        self.status_text = QLabel("Ready.")
        self.status_text.setObjectName("AppStatusText")
        self.status_text.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Fixed,
        )

        ornl_logo = QSvgWidget(str(ornl_logo_path()))
        ornl_logo.setObjectName("StatusOrnlLogo")
        ornl_logo.setFixedSize(18, 18)

        ornl_text = QLabel("Oak Ridge National Laboratory")
        ornl_text.setObjectName("StatusOrnlText")

        layout.addWidget(self.status_symbol)
        layout.addWidget(self.status_text, 1)
        layout.addWidget(ornl_logo)
        layout.addWidget(ornl_text)

    def set_message(self, message: str, level: str | None = None) -> None:
        """Set a user-facing status message and severity."""
        cleaned = " ".join(str(message or "Ready.").split())
        if not cleaned:
            cleaned = "Ready."
        status_level = level or status_level_from_message(cleaned)
        self.status_text.setText(cleaned)
        self.status_symbol.setProperty("status", status_level)
        self.status_symbol.clear()
        if status_level == "warning":
            self.status_symbol.setFixedSize(18, 18)
            self.status_symbol.setPixmap(
                function_warning_icon(16).pixmap(QSize(16, 16))
            )
        elif status_level == "error":
            self.status_symbol.setFixedSize(18, 18)
            self.status_symbol.setPixmap(function_error_icon(16).pixmap(QSize(16, 16)))
        else:
            self.status_symbol.setFixedSize(9, 9)
        self.status_symbol.style().unpolish(self.status_symbol)
        self.status_symbol.style().polish(self.status_symbol)
        self.status_symbol.update()
