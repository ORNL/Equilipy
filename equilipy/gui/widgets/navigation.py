"""Navigation and chrome widgets for the Equilipy GUI."""

from __future__ import annotations

from PySide6.QtCore import QModelIndex, QPersistentModelIndex, QRectF, Qt
from PySide6.QtGui import QColor, QPainter, QPen
from PySide6.QtWidgets import QPushButton, QTreeView, QTreeWidget


class LayoutToggleButton(QPushButton):
    """Chrome-free layout toggle with a VS Code-style panel icon."""

    def __init__(self, panel: str, tooltip: str, checked: bool = True):
        super().__init__("")
        self.panel = panel
        self.setObjectName("LayoutToggleButton")
        self.setToolTip(tooltip)
        self.setCheckable(True)
        self.setChecked(checked)
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setFixedSize(28, 28)

    def paintEvent(self, event) -> None:
        """Paint the panel-position icon with distinct checked state."""
        del event

        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)

        icon_rect = QRectF(
            (self.width() - 17) / 2,
            (self.height() - 15) / 2,
            17,
            15,
        )

        checked = self.isChecked()
        outline = QColor("#cfd6d8") if checked else QColor("#7a8588")
        fill = QColor("#8b9699") if checked else Qt.BrushStyle.NoBrush
        divider = QColor("#30383a") if checked else QColor("#687376")

        panel_rect = QRectF()
        if self.panel == "left":
            panel_rect = QRectF(
                icon_rect.left() + 1.5,
                icon_rect.top() + 1.5,
                4.7,
                12,
            )
        elif self.panel == "bottom":
            panel_rect = QRectF(
                icon_rect.left() + 1.5,
                icon_rect.bottom() - 5.8,
                14,
                4.3,
            )
        elif self.panel == "right":
            panel_rect = QRectF(
                icon_rect.right() - 6.2,
                icon_rect.top() + 1.5,
                4.7,
                12,
            )

        if checked and panel_rect.isValid():
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(fill)
            painter.drawRoundedRect(panel_rect, 1.4, 1.4)

        if panel_rect.isValid():
            painter.setPen(QPen(divider, 1.0))
            if self.panel in {"left", "right"}:
                x_pos = (
                    panel_rect.right()
                    if self.panel == "left"
                    else panel_rect.left()
                )
                painter.drawLine(
                    int(round(x_pos)),
                    int(round(icon_rect.top() + 2)),
                    int(round(x_pos)),
                    int(round(icon_rect.bottom() - 2)),
                )
            else:
                y_pos = panel_rect.top()
                painter.drawLine(
                    int(round(icon_rect.left() + 2)),
                    int(round(y_pos)),
                    int(round(icon_rect.right() - 2)),
                    int(round(y_pos)),
                )

        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.setPen(QPen(outline, 1.2))
        painter.drawRoundedRect(icon_rect, 3, 3)


class ClearableTreeView(QTreeView):
    """Tree view that clears row selection when empty space is clicked."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._drag_source_index: QPersistentModelIndex | None = None

    def mousePressEvent(self, event) -> None:
        """Clear visual selection on blank clicks without changing the editor."""
        index = self.indexAt(event.position().toPoint())
        if not index.isValid():
            self.selectionModel().clearSelection()
            return
        super().mousePressEvent(event)

    def startDrag(self, supported_actions) -> None:
        """Remember the source index so database drops can be copy-only."""
        current = self.currentIndex()
        self._drag_source_index = (
            QPersistentModelIndex(current) if current.isValid() else None
        )
        super().startDrag(Qt.DropAction.CopyAction)

    def dropEvent(self, event) -> None:
        """Route database tree drops through the main window copy handler."""
        source = (
            QModelIndex(self._drag_source_index)
            if self._drag_source_index is not None
            and self._drag_source_index.isValid()
            else QModelIndex()
        )
        target = self.indexAt(event.position().toPoint())
        handler = getattr(self.window(), "_copy_database_tree_payload", None)
        if source.isValid() and target.isValid() and callable(handler):
            if handler(source, target):
                event.setDropAction(Qt.DropAction.CopyAction)
                event.accept()
                return
        super().dropEvent(event)


class ClearableTreeWidget(QTreeWidget):
    """Tree widget that clears row selection when empty space is clicked."""

    def mousePressEvent(self, event) -> None:
        """Clear visual selection on blank clicks without changing the editor."""
        item = self.itemAt(event.position().toPoint())
        if item is None:
            self.clearSelection()
            return
        super().mousePressEvent(event)
