"""Reusable Qt form and table helpers for the GUI."""

# ruff: noqa: D103

from __future__ import annotations

from typing import Any

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractButton,
    QAbstractItemView,
    QFormLayout,
    QFrame,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPlainTextEdit,
    QSizePolicy,
    QTableWidget,
    QTableWidgetItem,
    QToolButton,
    QTreeWidget,
    QVBoxLayout,
    QWidget,
)

from .composition import (
    COMPOSITION_ROW_HEIGHT,
    CompositionEditor,
)


def apply_button_cursor(widget: QWidget) -> None:
    """Use a pointing cursor for clickable controls under one widget."""
    if isinstance(widget, QAbstractButton):
        widget.setCursor(Qt.CursorShape.PointingHandCursor)
    for button in widget.findChildren(QAbstractButton):
        button.setCursor(Qt.CursorShape.PointingHandCursor)


def form_layout() -> QFormLayout:
    form = QFormLayout()
    form.setLabelAlignment(Qt.AlignmentFlag.AlignRight)
    form.setFormAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
    form.setHorizontalSpacing(16)
    form.setVerticalSpacing(10)
    return form


def form_section(title: str) -> QFrame:
    section = QFrame()
    section.setObjectName("FormSection")
    section.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum)
    layout = QVBoxLayout(section)
    layout.setContentsMargins(16, 14, 16, 16)
    layout.setSpacing(12)

    heading = QLabel(title)
    heading.setObjectName("SectionTitle")
    layout.addWidget(heading)
    return section


def untitled_form_section() -> QFrame:
    """Return a form section without a visible heading."""
    section = QFrame()
    section.setObjectName("FormSection")
    section.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum)
    layout = QVBoxLayout(section)
    layout.setContentsMargins(16, 14, 16, 16)
    layout.setSpacing(12)
    return section


def collapsible_form_section(
    title: str,
    actions: list[QWidget] | None = None,
    *,
    collapsed: bool = False,
) -> tuple[QFrame, QVBoxLayout]:
    section = QFrame()
    section.setObjectName("FormSection")
    section.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum)
    outer_layout = QVBoxLayout(section)
    outer_layout.setContentsMargins(16, 14, 16, 16)
    outer_layout.setSpacing(12)

    header = QHBoxLayout()
    header.setSpacing(8)
    toggle = QToolButton()
    toggle.setObjectName("SectionCollapseButton")
    heading = QLabel(title)
    heading.setObjectName("SectionTitle")
    header.addWidget(toggle)
    header.addWidget(heading)
    header.addStretch(1)
    for action in actions or []:
        header.addWidget(action)
    outer_layout.addLayout(header)

    body = QWidget()
    body.setObjectName("CollapsibleSectionBody")
    body_layout = QVBoxLayout(body)
    body_layout.setContentsMargins(0, 0, 0, 0)
    body_layout.setSpacing(12)
    outer_layout.addWidget(body)

    def set_collapsed(value: bool) -> None:
        section.setProperty("collapsed", value)
        body.setVisible(not value)
        toggle.setArrowType(
            Qt.ArrowType.RightArrow if value else Qt.ArrowType.DownArrow
        )
        toggle.setToolTip("Expand section" if value else "Collapse section")

    toggle.clicked.connect(
        lambda _checked=False: set_collapsed(
            not bool(section.property("collapsed"))
        )
    )
    set_collapsed(collapsed)
    return section, body_layout


def form_section_with_actions(title: str, actions: list[QWidget]) -> QFrame:
    section = QFrame()
    section.setObjectName("FormSection")
    section.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum)
    layout = QVBoxLayout(section)
    layout.setContentsMargins(16, 14, 16, 16)
    layout.setSpacing(12)

    header = QHBoxLayout()
    header.setSpacing(8)
    heading = QLabel(title)
    heading.setObjectName("SectionTitle")
    header.addWidget(heading)
    header.addStretch(1)
    for action in actions:
        header.addWidget(action)
    layout.addLayout(header)
    return section


def text_section(title: str, text: str) -> QFrame:
    section = form_section(title)
    section.layout().addWidget(read_only_text(text, 160))
    return section


def table_section(
    title: str,
    headers: tuple[str, ...],
    rows: list[tuple[Any, ...]],
) -> QFrame:
    section = form_section(title)
    table = QTableWidget(len(rows), len(headers))
    table.setObjectName("FormTable")
    table.setHorizontalHeaderLabels(headers)
    table.verticalHeader().setVisible(False)
    table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
    table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
    table.setAlternatingRowColors(False)
    table.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)

    for row_index, row in enumerate(rows):
        for column_index, value in enumerate(row):
            item = QTableWidgetItem(str(value))
            item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            table.setItem(row_index, column_index, item)

    table.resizeColumnsToContents()
    fit_table_columns_to_headers(table)
    table.horizontalHeader().setStretchLastSection(True)
    disable_inner_scrollbars(table)
    fit_table_to_rows(table)
    section.layout().addWidget(table)
    return section


def fit_table_columns_to_headers(table: QTableWidget) -> None:
    """Ensure compact columns are at least wide enough for their header text."""
    metrics = table.horizontalHeader().fontMetrics()
    for column in range(table.columnCount()):
        header = table.horizontalHeaderItem(column)
        if header is None:
            continue
        header_width = metrics.horizontalAdvance(header.text()) + 28
        table.setColumnWidth(column, max(table.columnWidth(column), header_width))


def inline_table(
    label: str,
    headers: tuple[str, ...],
    rows: list[tuple[Any, ...]],
) -> QWidget:
    wrapper = QWidget()
    wrapper_layout = QVBoxLayout(wrapper)
    wrapper_layout.setContentsMargins(0, 0, 0, 0)
    wrapper_layout.setSpacing(4)

    heading = QLabel(label)
    heading.setObjectName("DatabaseListLabel")
    wrapper_layout.addWidget(heading)

    table = QTableWidget(len(rows), len(headers))
    table.setObjectName("FormTable")
    table.setHorizontalHeaderLabels(headers)
    table.verticalHeader().setVisible(False)
    table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
    table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
    table.setAlternatingRowColors(False)
    table.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)

    for row_index, row in enumerate(rows):
        for column_index, value in enumerate(row):
            item = QTableWidgetItem(str(value))
            item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            table.setItem(row_index, column_index, item)

    table.resizeColumnsToContents()
    table.horizontalHeader().setStretchLastSection(True)
    disable_inner_scrollbars(table)
    fit_table_to_rows(table)
    wrapper_layout.addWidget(table)
    return wrapper


def read_only_line(value: Any) -> QLineEdit:
    field = QLineEdit(str(value))
    field.setReadOnly(True)
    return field


def read_only_text(value: Any, height: int) -> QPlainTextEdit:
    field = QPlainTextEdit(str(value))
    field.setReadOnly(True)
    field.setMinimumHeight(height)
    return field


def disable_inner_scrollbars(widget: QWidget) -> None:
    widget.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
    widget.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)


def latest_status_text(status: QPlainTextEdit) -> str:
    value = getattr(status, "last_status", "")
    if isinstance(value, str) and value.strip():
        return value.strip()
    return status.toPlainText().strip()


def fit_table_to_rows(
    table: QTableWidget | CompositionEditor,
    min_rows: int = 1,
) -> None:
    if isinstance(table, CompositionEditor):
        table.fit_to_rows(min_rows=min_rows)
        return
    table.resizeRowsToContents()
    if bool(table.property("composition_table")):
        table.verticalHeader().setDefaultSectionSize(COMPOSITION_ROW_HEIGHT)
        table.verticalHeader().setMinimumSectionSize(COMPOSITION_ROW_HEIGHT)
        for row in range(table.rowCount()):
            table.setRowHeight(
                row,
                max(table.rowHeight(row), COMPOSITION_ROW_HEIGHT),
            )
    row_count = max(min_rows, table.rowCount())
    row_height = table.verticalHeader().defaultSectionSize()
    rows_height = sum(table.rowHeight(row) for row in range(table.rowCount()))
    rows_height += max(0, row_count - table.rowCount()) * row_height
    header_height = table.horizontalHeader().height()
    frame = table.frameWidth() * 2
    margin = 14
    table.setFixedHeight(header_height + rows_height + frame + margin)


def fit_tree_to_items(tree: QTreeWidget, min_rows: int = 4) -> None:
    row_count = max(min_rows, tree_visible_item_count(tree))
    row_height = tree.sizeHintForRow(0)
    if row_height <= 0:
        row_height = 30
    frame = tree.frameWidth() * 2
    header_height = 0 if tree.isHeaderHidden() else tree.header().height()
    tree.setFixedHeight(header_height + row_count * row_height + frame + 8)


def tree_visible_item_count(tree: QTreeWidget) -> int:
    count = 0
    for row in range(tree.topLevelItemCount()):
        parent = tree.topLevelItem(row)
        count += 1
        if parent.isExpanded():
            count += parent.childCount()
    return count
