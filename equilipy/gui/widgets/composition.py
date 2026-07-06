"""Composition editor widget used by calculation setup forms."""

from __future__ import annotations

from typing import Any

from PySide6.QtCore import QEvent, Qt, Signal
from PySide6.QtGui import QKeyEvent
from PySide6.QtWidgets import (
    QFrame,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QTableWidgetItem,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

COMPOSITION_ROW_HEIGHT = 48
COMPOSITION_INPUT_MIN_HEIGHT = 38


class _CompositionIndex:
    """Small row index shim for composition row selection."""

    def __init__(self, row: int):
        self._row = row

    def row(self) -> int:
        return self._row


class CompositionEditor(QFrame):
    """Spreadsheet-free species/amount editor with table-like helper methods."""

    itemChanged = Signal(object)

    def __init__(self):
        super().__init__()
        self.setObjectName("CompositionEditor")
        self.setProperty("composition_table", True)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self._rows: list[dict[str, Any]] = []
        self._selected_row: int | None = None
        self._next_focus_widget: QWidget | None = None
        self._previous_focus_widget: QWidget | None = None

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        header = QFrame()
        header.setObjectName("CompositionEditorHeader")
        header_layout = QHBoxLayout(header)
        header_layout.setContentsMargins(12, 8, 12, 8)
        header_layout.setSpacing(8)
        species_label = QLabel("Species")
        species_label.setObjectName("CompositionEditorHeaderLabel")
        species_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        amount_label = QLabel("Amount")
        amount_label.setObjectName("CompositionEditorHeaderLabel")
        amount_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        header_divider = QFrame()
        header_divider.setObjectName("CompositionEditorHeaderDivider")
        header_divider.setFixedWidth(1)
        header_layout.addWidget(species_label, 1)
        header_layout.addWidget(header_divider)
        header_layout.addWidget(amount_label, 3)
        layout.addWidget(header)

        self._rows_container = QWidget()
        self._rows_layout = QVBoxLayout(self._rows_container)
        self._rows_layout.setContentsMargins(10, 10, 10, 10)
        self._rows_layout.setSpacing(8)
        layout.addWidget(self._rows_container)

    def rowCount(self) -> int:
        """Return the number of composition rows."""
        return len(self._rows)

    def columnCount(self) -> int:
        """Return the fixed Species/Amount column count."""
        return 2

    def insertRow(self, row: int) -> None:
        """Insert a new blank composition row."""
        row = max(0, min(row, len(self._rows)))
        frame = QFrame()
        frame.setObjectName("CompositionRow")
        frame.setProperty("selected", False)
        frame.setCursor(Qt.CursorShape.PointingHandCursor)
        frame.installEventFilter(self)
        row_layout = QHBoxLayout(frame)
        row_layout.setContentsMargins(0, 0, 0, 0)
        row_layout.setSpacing(8)
        frame.setMinimumHeight(COMPOSITION_ROW_HEIGHT)
        self._rows.insert(
            row,
            {
                "frame": frame,
                "items": [None, None],
                "widgets": [None, None],
            },
        )
        self._rebuild_rows()

    def removeRow(self, row: int) -> None:
        """Remove one composition row."""
        if row < 0 or row >= len(self._rows):
            return
        data = self._rows.pop(row)
        frame = data.get("frame")
        if isinstance(frame, QWidget):
            frame.hide()
            frame.deleteLater()
        if self._selected_row is not None:
            if self._selected_row == row:
                self._selected_row = None
            elif self._selected_row > row:
                self._selected_row -= 1
        self._rebuild_rows()

    def setRowCount(self, count: int) -> None:
        """Resize the composition row list."""
        count = max(0, count)
        while len(self._rows) > count:
            self.removeRow(len(self._rows) - 1)
        while len(self._rows) < count:
            self.insertRow(len(self._rows))

    def setItem(self, row: int, column: int, item: QTableWidgetItem) -> None:
        """Store the hidden item used by shared table helpers."""
        if not self._is_valid_cell(row, column):
            return
        self._rows[row]["items"][column] = item
        widget = self.cellWidget(row, column)
        if isinstance(widget, QLineEdit) and widget.text() != item.text():
            previous_blocked = widget.blockSignals(True)
            widget.setText(item.text())
            widget.blockSignals(previous_blocked)

    def item(self, row: int, column: int) -> QTableWidgetItem | None:
        """Return the hidden item for a cell."""
        if not self._is_valid_cell(row, column):
            return None
        item = self._rows[row]["items"][column]
        return item if isinstance(item, QTableWidgetItem) else None

    def setCellWidget(self, row: int, column: int, widget: QWidget) -> None:
        """Set the visible editor widget for a cell."""
        if not self._is_valid_cell(row, column):
            return
        old_widget = self._rows[row]["widgets"][column]
        if isinstance(old_widget, QWidget) and old_widget is not widget:
            old_widget.hide()
            old_widget.deleteLater()
        self._rows[row]["widgets"][column] = widget
        self._wire_cell_widget(row, column, widget)
        self._refresh_row_widget(row)
        self._refresh_row_indices()

    def cellWidget(self, row: int, column: int) -> QWidget | None:
        """Return the visible editor widget for a cell."""
        if not self._is_valid_cell(row, column):
            return None
        widget = self._rows[row]["widgets"][column]
        return widget if isinstance(widget, QWidget) else None

    def selectedIndexes(self) -> list[_CompositionIndex]:
        """Return selected row indexes in a table-compatible shape."""
        if self._selected_row is None:
            return []
        return [_CompositionIndex(self._selected_row)]

    def clearSelection(self) -> None:
        """Clear the selected composition row."""
        self._selected_row = None
        self._sync_selection_styles()

    def selectRow(self, row: int) -> None:
        """Select one composition row."""
        if row < 0 or row >= len(self._rows):
            self.clearSelection()
            return
        self._selected_row = row
        self._sync_selection_styles()

    def setCurrentCell(self, row: int, column: int) -> None:
        """Select a row and focus the requested editor."""
        if not self._is_valid_cell(row, column):
            return
        self.selectRow(row)
        widget = self.cellWidget(row, column)
        if isinstance(widget, QLineEdit):
            widget.setFocus()
        elif isinstance(widget, QWidget):
            edits = widget.findChildren(QLineEdit)
            if edits:
                edits[0].setFocus()

    def set_boundary_focus_widgets(
        self,
        *,
        next_widget: QWidget | None = None,
        previous_widget: QWidget | None = None,
    ) -> None:
        """Set explicit focus handoff targets outside the composition editor."""
        self._next_focus_widget = next_widget
        self._previous_focus_widget = previous_widget

    def rowHeight(self, row: int) -> int:
        """Return the effective row height."""
        if row < 0 or row >= len(self._rows):
            return 0
        frame = self._rows[row].get("frame")
        if isinstance(frame, QWidget):
            return max(frame.sizeHint().height(), COMPOSITION_ROW_HEIGHT)
        return COMPOSITION_ROW_HEIGHT

    def setRowHeight(self, row: int, height: int) -> None:
        """Set a minimum row height."""
        if row < 0 or row >= len(self._rows):
            return
        frame = self._rows[row].get("frame")
        if isinstance(frame, QWidget):
            frame.setMinimumHeight(height)

    def resizeRowsToContents(self) -> None:
        """Keep table helper compatibility; row widgets size themselves."""
        return

    def fit_to_rows(self, min_rows: int = 1) -> None:
        """Resize the editor to fit all row widgets."""
        row_count = max(min_rows, len(self._rows))
        rows_height = (
            sum(self.rowHeight(row) for row in range(len(self._rows)))
            + max(0, row_count - len(self._rows)) * COMPOSITION_ROW_HEIGHT
        )
        spacing = self._rows_layout.spacing() * max(0, row_count - 1)
        margins = self._rows_layout.contentsMargins()
        header_height = 36
        total = (
            header_height
            + rows_height
            + spacing
            + margins.top()
            + margins.bottom()
            + 2
        )
        self.setFixedHeight(total)

    def eventFilter(self, watched, event) -> bool:
        """Select the row when the user clicks or focuses inside it."""
        if (
            event.type() == QEvent.Type.KeyPress
            and isinstance(event, QKeyEvent)
            and event.key() in {Qt.Key.Key_Tab, Qt.Key.Key_Backtab}
        ):
            forward = event.key() == Qt.Key.Key_Tab and not (
                event.modifiers() & Qt.KeyboardModifier.ShiftModifier
            )
            if self._move_composition_focus(watched, forward):
                return True

        if event.type() in {
            QEvent.Type.FocusIn,
            QEvent.Type.MouseButtonPress,
        }:
            row = watched.property("composition_row")
            if isinstance(row, int):
                self.selectRow(row)
        return super().eventFilter(watched, event)

    def _is_valid_cell(self, row: int, column: int) -> bool:
        return 0 <= row < len(self._rows) and 0 <= column < 2

    def _rebuild_rows(self) -> None:
        while self._rows_layout.count():
            self._rows_layout.takeAt(0)
        for row, _data in enumerate(self._rows):
            self._refresh_row_widget(row)
            frame = self._rows[row]["frame"]
            self._rows_layout.addWidget(frame)
        self._rows_layout.addStretch(1)
        self._refresh_row_indices()
        self._sync_selection_styles()

    def _refresh_row_widget(self, row: int) -> None:
        if row < 0 or row >= len(self._rows):
            return
        data = self._rows[row]
        frame = data["frame"]
        row_layout = frame.layout()
        if row_layout is None:
            return
        active_widgets = {
            widget for widget in data["widgets"] if isinstance(widget, QWidget)
        }
        while row_layout.count():
            item = row_layout.takeAt(0)
            widget = item.widget()
            if widget is not None and widget not in active_widgets:
                widget.hide()
                widget.deleteLater()
        for column, stretch in ((0, 1), (1, 3)):
            widget = data["widgets"][column]
            if isinstance(widget, QWidget):
                row_layout.addWidget(widget, stretch)
            else:
                filler = QLabel("")
                filler.setMinimumHeight(COMPOSITION_INPUT_MIN_HEIGHT)
                row_layout.addWidget(filler, stretch)

    def _wire_cell_widget(self, row: int, column: int, widget: QWidget) -> None:
        del column
        self._mark_widget_row(widget, row)
        if isinstance(widget, QLineEdit):
            widget.textChanged.connect(lambda _value, r=row: self._emit_row_changed(r))
        else:
            for child in widget.findChildren(QWidget):
                self._mark_widget_row(child, row)
            for editor in widget.findChildren(QLineEdit):
                editor.textChanged.connect(
                    lambda _value, r=row: self._emit_row_changed(r)
                )
            for button in widget.findChildren(QToolButton):
                button.toggled.connect(
                    lambda _checked, r=row: self._emit_row_changed(r)
                )

    def _mark_widget_row(self, widget: QWidget, row: int) -> None:
        widget.setProperty("composition_row", row)
        widget.installEventFilter(self)

    def _refresh_row_indices(self) -> None:
        for row, data in enumerate(self._rows):
            frame = data["frame"]
            if isinstance(frame, QWidget):
                frame.setProperty("composition_row", row)
            for widget in data["widgets"]:
                if isinstance(widget, QWidget):
                    self._mark_widget_row(widget, row)
                    for child in widget.findChildren(QWidget):
                        self._mark_widget_row(child, row)

    def _emit_row_changed(self, row: int) -> None:
        if self.signalsBlocked() or row < 0 or row >= len(self._rows):
            return
        item = self.item(row, 0) or self.item(row, 1)
        self.itemChanged.emit(item)

    def _sync_selection_styles(self) -> None:
        for row, data in enumerate(self._rows):
            frame = data.get("frame")
            if not isinstance(frame, QWidget):
                continue
            frame.setProperty("selected", row == self._selected_row)
            frame.style().unpolish(frame)
            frame.style().polish(frame)
            frame.update()

    def _composition_focus_widgets(self) -> list[QLineEdit]:
        widgets: list[QLineEdit] = []
        for row in range(len(self._rows)):
            species_widget = self.cellWidget(row, 0)
            if isinstance(species_widget, QLineEdit):
                widgets.append(species_widget)
            amount_widget = self.cellWidget(row, 1)
            if isinstance(amount_widget, QLineEdit):
                widgets.append(amount_widget)
            elif isinstance(amount_widget, QWidget):
                widgets.extend(
                    editor
                    for editor in amount_widget.findChildren(QLineEdit)
                    if editor.isVisible() and editor.isEnabled()
                )
        return widgets

    def _move_composition_focus(self, watched: QWidget, forward: bool) -> bool:
        if not isinstance(watched, QLineEdit):
            return False
        widgets = self._composition_focus_widgets()
        if watched not in widgets:
            return False
        index = widgets.index(watched)
        next_index = index + (1 if forward else -1)
        if 0 <= next_index < len(widgets):
            widgets[next_index].setFocus(Qt.FocusReason.TabFocusReason)
            widgets[next_index].selectAll()
            return True
        boundary_widget = (
            self._next_focus_widget if forward else self._previous_focus_widget
        )
        return self._focus_boundary_widget(boundary_widget)

    @staticmethod
    def _focus_boundary_widget(widget: QWidget | None) -> bool:
        """Focus a configured widget or its first visible line edit."""
        if widget is None or not widget.isEnabled() or not widget.isVisible():
            return False
        if isinstance(widget, QLineEdit):
            widget.setFocus(Qt.FocusReason.TabFocusReason)
            widget.selectAll()
            return True
        editors = [
            editor
            for editor in widget.findChildren(QLineEdit)
            if editor.isEnabled() and editor.isVisible()
        ]
        if not editors:
            return False
        editors[0].setFocus(Qt.FocusReason.TabFocusReason)
        editors[0].selectAll()
        return True
