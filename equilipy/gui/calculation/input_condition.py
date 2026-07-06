"""Composition, condition, and phase-selection GUI helpers."""
from __future__ import annotations

# ruff: noqa: F401,F403,F405,F811,F821,E402,I001,E501

import base64
import copy
import csv
import json
import os
import pickle
import re
import shutil
from dataclasses import asdict
from itertools import product
from pathlib import Path
from typing import Any, Callable, Iterable, Iterator

import numpy as np
from PySide6.QtCore import QModelIndex, QSize, Qt
from PySide6.QtGui import QAction, QBrush, QColor, QFontDatabase, QIcon
from PySide6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QCheckBox,
    QComboBox,
    QDialog,
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMenu,
    QMessageBox,
    QRadioButton,
    QScrollArea,
    QSizePolicy,
    QStackedWidget,
    QStyle,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QToolButton,
    QTreeWidget,
    QTreeWidgetItem,
    QPlainTextEdit,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from equilipy.composition import expand_condition_species as _expand_condition_species
from equilipy.composition import (
    formula_or_database_species_stoichiometry as _shared_formula_or_database_species_stoichiometry,
)
from equilipy.database_ir import (
    DatabaseIR,
    Diagnostic,
    FunctionDefinition,
    GibbsRange,
    Parameter,
    Phase,
    SourceRef,
    Species,
)
from equilipy.database_ir.tdb_syntax import parameter_command_prefix_match
from equilipy.phase_selection import (
    default_scheil_phase_selection as _metadata_scheil_phase_selection,
)
from equilipy.results import ResultTable
from equilipy.results.equilib import EquilibResult
from equilipy.results.scheil import ScheilResult
from equilipy.utils import G2HSCp, HSCp2G, NeumanKoppHSCp

from ..assets import calculation_icon_path, math_font_path
from .state import CalculationDatabase, CalculationModule, CalculationSession
from ..models import (
    CompoundPhasePayload,
    CompoundRecord,
    DatabaseTreeModel,
    ThermoRangePayload,
    function_error_icon,
    function_warning_icon,
    thermo_range_warning_key,
)
from ..units import AMOUNT_UNIT_OPTIONS as _AMOUNT_UNIT_OPTIONS
from ..widgets import forms as _forms
from ..widgets.composition import (
    COMPOSITION_INPUT_MIN_HEIGHT,
    COMPOSITION_ROW_HEIGHT,
    CompositionEditor,
)

_apply_button_cursor = _forms.apply_button_cursor
_collapsible_form_section = _forms.collapsible_form_section
_disable_inner_scrollbars = _forms.disable_inner_scrollbars
_fit_table_columns_to_headers = _forms.fit_table_columns_to_headers
_fit_table_to_rows = _forms.fit_table_to_rows
_fit_tree_to_items = _forms.fit_tree_to_items
_form_layout = _forms.form_layout
_form_section = _forms.form_section
_form_section_with_actions = _forms.form_section_with_actions
_inline_table = _forms.inline_table
_latest_status_text = _forms.latest_status_text
_read_only_line = _forms.read_only_line
_read_only_text = _forms.read_only_text
_table_section = _forms.table_section
_text_section = _forms.text_section
_untitled_form_section = _forms.untitled_form_section

_RESULT_OVERVIEW_ITEM_ID = "__result_overview__"
_THERMO_COEFFICIENT_BOX_WIDTH = 170
_THERMO_POWER_BOX_WIDTH = 80
_THERMO_TEMPERATURE_BOX_WIDTH = round(_THERMO_COEFFICIENT_BOX_WIDTH * 0.5)
_GIBBS_COEFFICIENT_SLOTS = (0, 1, 2, 3, 4, 5, 6, 8, 10, 12)
_FUNCTION_TOKEN_RE = re.compile(r"(?<![A-Za-z0-9_])([A-Za-z_][A-Za-z0-9_]*)")
_PORTABLE_TDB_FUNCTION_NAME_WIDTH = 8
_MAX_TDB_FUNCTION_NAME_WIDTH = 25
_GIBBS_CONTINUITY_TOLERANCE = 1e-5
_MATH_FONT_FAMILY: str | None = None


def _pure_compound_names(data: Any) -> list[str]:
    if not isinstance(data, dict):
        return []
    species_names = data.get("cSpeciesNameCS")
    if species_names is None:
        return []
    try:
        n_pure_species = int(data.get("nPureSpeciesCS", 0) or 0)
    except (TypeError, ValueError):
        return []
    if n_pure_species <= 0:
        return []

    start_index = max(0, len(species_names) - n_pure_species)
    return [
        str(species).strip()
        for species in species_names[start_index:]
        if str(species).strip()
        and not _is_reference_ser_phase(str(species).strip())
    ]


def _database_phase_names(data: Any) -> list[str]:
    if not isinstance(data, dict):
        return []

    phase_sources = data.get("cPhaseNames")
    if phase_sources is None:
        phase_sources = list(data.get("cSolnPhaseNameCS") or [])
        phase_sources.extend(_pure_compound_names(data))

    phase_names: list[str] = []
    for phase in phase_sources:
        name = str(phase).strip()
        if name and name not in phase_names and not _is_reference_ser_phase(name):
            phase_names.append(name)
    return phase_names


def _is_reference_ser_phase(phase_name: str) -> bool:
    return bool(re.search(r"_SER\s*\(s\)\s*$", phase_name.strip(), re.IGNORECASE))


def _solution_phase_overview_rows(data: Any) -> list[tuple[Any, ...]]:
    if not isinstance(data, dict):
        return []
    phase_names = data.get("cSolnPhaseNameCS")
    if phase_names is None:
        return []

    phase_types = data.get("cSolnPhaseTypeCS")
    species_counts = data.get("nSolnPhaseCS")
    rows: list[tuple[Any, ...]] = []
    species_start = 0
    for index, phase_name in enumerate(phase_names):
        name = str(phase_name).strip()
        if not name:
            continue

        model = ""
        if phase_types is not None and index < len(phase_types):
            model = str(phase_types[index]).strip()

        species_count = 0
        if species_counts is not None and index < len(species_counts):
            species_count = int(species_counts[index])
        species = _solution_phase_species(data, species_start, species_count)
        species_start += species_count

        rows.append(
            (
                len(rows) + 1,
                name,
                model,
                _solution_phase_structure(data, index),
                _format_species_list(species),
            )
        )
    return rows


def _solution_phase_species(
    data: dict[str, Any],
    start_index: int,
    count: int,
) -> list[str]:
    species_names = data.get("cSpeciesNameCS")
    if species_names is None or count <= 0:
        return []
    return [
        str(species).strip()
        for species in species_names[start_index : start_index + count]
        if str(species).strip()
    ]


def _solution_phase_structure(
    data: dict[str, Any],
    phase_index: int,
) -> str:
    sublattice_indices = data.get("iPhaseSublatticeCS")
    if sublattice_indices is None or phase_index >= len(sublattice_indices):
        return ""
    sublattice_index = int(sublattice_indices[phase_index])
    if sublattice_index <= 0:
        return ""

    sublattice_row = sublattice_index - 1
    sublattice_counts = data.get("nSublatticePhaseCS")
    constituent_counts = data.get("nConstituentSublatticeCS")
    constituent_names = data.get("cConstituentNameSUBCS")
    site_ratios = data.get("dStoichSublatticeCS")
    if (
        sublattice_counts is None
        or constituent_counts is None
        or constituent_names is None
        or site_ratios is None
        or sublattice_row >= len(sublattice_counts)
    ):
        return ""

    parts: list[str] = []
    for sublattice in range(int(sublattice_counts[sublattice_row])):
        constituent_count = int(constituent_counts[sublattice_row, sublattice])
        names = [
            _char_array_text(
                constituent_names[sublattice_row, sublattice, constituent]
            )
            for constituent in range(constituent_count)
        ]
        names = [name for name in names if name]
        if not names:
            continue
        ratio = _format_site_ratio(float(site_ratios[sublattice_row, sublattice]))
        parts.append(f"({', '.join(names)}){ratio}")
    return "".join(parts)


def _format_species_list(species: list[str]) -> str:
    if not species:
        return ""
    return f"{', '.join(species)} ({len(species)})"


def _char_array_text(value: Any) -> str:
    if hasattr(value, "tobytes"):
        return value.tobytes().decode(errors="ignore").strip()
    return str(value).strip()


def _format_site_ratio(value: float) -> str:
    if value.is_integer():
        return str(int(value))
    return f"{value:g}"


def _is_composition_table(widget: Any) -> bool:
    return isinstance(widget, (CompositionEditor, QTableWidget)) and bool(
        widget.property("composition_table")
    )


def _combo_box(values: list[str]) -> QComboBox:
    combo = QComboBox()
    combo.addItems(values)
    return combo


def _check_box(text: str, checked: bool) -> QCheckBox:
    check_box = QCheckBox(text)
    check_box.setChecked(checked)
    return check_box


def _composition_table() -> CompositionEditor:
    table = CompositionEditor()
    _add_composition_row(table)
    _fit_table_to_rows(table, min_rows=3)
    return table


def _add_composition_row(table: QTableWidget | CompositionEditor) -> None:
    row = table.rowCount()
    table.insertRow(row)
    _set_composition_cell(table, row, 0, "", "e.g. Mg2Si")
    if bool(table.property("allow_amount_ranges")):
        _set_composition_amount_range_cell(table, row, "")
    else:
        _set_composition_cell(table, row, 1, "", _composition_amount_placeholder(table))
    _fit_table_to_rows(table, min_rows=3)


def _composition_amount_placeholder(table: QTableWidget | CompositionEditor) -> str:
    if bool(table.property("allow_amount_ranges")):
        return "[initial] [Final] [Step]"
    if _fixed_total_for_table(table) is not None:
        return "Amount"
    return "e.g. 25"


def _amount_unit_fixed_total(amount_unit: str) -> float | None:
    unit = amount_unit.strip().lower()
    if unit in {"mol%", "at%", "wt%", "wt.%"}:
        return 100.0
    if unit in {
        "mole fraction",
        "atom fraction",
        "mass fraction",
        "weight fraction",
    }:
        return 1.0
    return None


def _fixed_total_for_table(table: QTableWidget | CompositionEditor) -> float | None:
    return _amount_unit_fixed_total(str(table.property("amount_unit") or "moles"))


def _balance_total_text(total: float | None) -> str:
    if total is None:
        return ""
    return _format_amount_value(total)


def _balance_placeholder(total: float | None) -> str:
    total_text = _balance_total_text(total)
    return f"Balance species: {total_text} total amount"


def _balance_expression_from_amount_text(text: str) -> str:
    cleaned = text.strip()
    if cleaned.startswith("?"):
        return cleaned if cleaned != "?" else "?100"
    try:
        total = float(cleaned)
    except ValueError:
        return "?100"
    if total <= 0:
        return "?100"
    return f"?{_format_amount_value(total)}"


def _amount_widget_text_property(widget: QWidget, name: str) -> str:
    return str(widget.property(name) or "").strip()


def _set_line_edit_text_silent(editor: QLineEdit, text: str) -> None:
    previous = editor.blockSignals(True)
    editor.setText(text)
    editor.blockSignals(previous)


def _other_open_total_balance_row(
    table: QTableWidget | CompositionEditor,
    active_row: int,
) -> int | None:
    for row in range(table.rowCount()):
        if row == active_row:
            continue
        if _table_cell_text(table, row, 1).strip().startswith("?"):
            return row
    return None


def _warn_duplicate_open_total_balance(
    table: QTableWidget | CompositionEditor,
) -> None:
    QMessageBox.warning(
        table.window(),
        "Balance Species Already Set",
        "Only one ?Total balance species is allowed for batch calculations.",
    )


def _set_open_total_range_editor_balance_state(
    editors: list[QLineEdit],
    is_balance: bool,
) -> None:
    if len(editors) < 3:
        return
    if is_balance:
        for editor in editors[1:3]:
            _set_line_edit_text_silent(editor, "")
            editor.setEnabled(False)
            editor.setPlaceholderText("Disabled for ?Total")
        return

    editors[1].setEnabled(True)
    editors[1].setPlaceholderText("Final")
    editors[2].setEnabled(True)
    editors[2].setPlaceholderText("Step")


def _set_amount_row_balance_expression(
    table: QTableWidget | CompositionEditor,
    row: int,
    total_text: str | None = None,
) -> None:
    total = _fixed_total_for_table(table)
    text = f"?{_balance_total_text(total)}" if total is not None else total_text
    if text is None:
        text = _balance_expression_from_amount_text(_table_cell_text(table, row, 1))
    if (
        text.strip().startswith("?")
        and total is None
        and bool(table.property("allow_amount_ranges"))
        and _other_open_total_balance_row(table, row) is not None
    ):
        _warn_duplicate_open_total_balance(table)
        return
    _set_table_cell_text(table, row, 1, text)


def _apply_amount_balance_from_editor(
    editor: QLineEdit,
    table: QTableWidget | CompositionEditor,
    row: int,
) -> None:
    total = _fixed_total_for_table(table)
    total_text = (
        f"?{_balance_total_text(total)}"
        if total is not None
        else _balance_expression_from_amount_text(editor.text())
    )
    _set_amount_row_balance_expression(table, row, total_text)


def _composition_editor_row(
    table: QTableWidget | CompositionEditor,
    editor: QLineEdit,
    fallback_row: int,
) -> int:
    for row in range(table.rowCount()):
        for column in range(table.columnCount()):
            widget = table.cellWidget(row, column)
            if widget is editor:
                return row
            if isinstance(widget, QWidget) and editor in widget.findChildren(QLineEdit):
                return row
    return fallback_row


def _amount_balance_editor_for_row(
    table: QTableWidget | CompositionEditor,
    row: int,
) -> QLineEdit | None:
    widget = table.cellWidget(row, 1)
    if isinstance(widget, QLineEdit):
        return widget
    if isinstance(widget, QWidget):
        editors = widget.findChildren(QLineEdit)
        if editors:
            return editors[0]
    return None


def _apply_amount_balance_from_row(
    table: QTableWidget | CompositionEditor,
    row: int,
) -> None:
    editor = _amount_balance_editor_for_row(table, row)
    if editor is None:
        _set_amount_row_balance_expression(table, row)
        return
    _apply_amount_balance_from_editor(editor, table, row)


def _duplicate_composition_row(
    table: QTableWidget | CompositionEditor,
    row: int,
) -> None:
    if row < 0 or row >= table.rowCount():
        return
    insert_row = row + 1
    species = _table_cell_text(table, row, 0)
    amount = _table_cell_text(table, row, 1)
    table.insertRow(insert_row)
    _set_composition_cell(table, insert_row, 0, species, "e.g. Mg2Si")
    if bool(table.property("allow_amount_ranges")):
        _set_composition_amount_range_cell(table, insert_row, amount)
    else:
        _set_composition_amount_single_cell(table, insert_row, amount)
    _refresh_composition_amount_editors(table)
    table.clearSelection()
    table.selectRow(insert_row)
    table.setCurrentCell(insert_row, 0)
    _fit_table_to_rows(table, min_rows=3)


def _delete_composition_row(
    table: QTableWidget | CompositionEditor,
    row: int,
) -> None:
    if row < 0 or row >= table.rowCount():
        return
    next_row = min(row, max(0, table.rowCount() - 2))
    table.removeRow(row)
    if table.rowCount() == 0:
        _add_composition_row(table)
        next_row = 0
    _ensure_required_balance_button(table)
    table.clearSelection()
    table.selectRow(next_row)
    table.setCurrentCell(next_row, 0)
    _fit_table_to_rows(table, min_rows=3)


def _line_edit_undo_available(editor: QLineEdit) -> bool:
    checker = getattr(editor, "isUndoAvailable", None)
    return bool(checker()) if callable(checker) else True


def _line_edit_redo_available(editor: QLineEdit) -> bool:
    checker = getattr(editor, "isRedoAvailable", None)
    return bool(checker()) if callable(checker) else True


def _composition_input_context_menu(
    editor: QLineEdit,
    table: QTableWidget | CompositionEditor,
    row: int,
) -> QMenu:
    menu = QMenu(editor)
    undo_action = QAction("Undo", menu)
    undo_action.setEnabled(_line_edit_undo_available(editor))
    undo_action.triggered.connect(editor.undo)
    menu.addAction(undo_action)

    redo_action = QAction("Redo", menu)
    redo_action.setEnabled(_line_edit_redo_available(editor))
    redo_action.triggered.connect(editor.redo)
    menu.addAction(redo_action)
    menu.addSeparator()

    duplicate_action = QAction("Duplicate", menu)
    duplicate_action.triggered.connect(
        lambda _checked=False: _duplicate_composition_row(
            table,
            _composition_editor_row(table, editor, row),
        )
    )
    menu.addAction(duplicate_action)

    delete_action = QAction("Delete", menu)
    delete_action.triggered.connect(
        lambda _checked=False: _delete_composition_row(
            table,
            _composition_editor_row(table, editor, row),
        )
    )
    menu.addAction(delete_action)
    menu.addSeparator()

    balance_action = QAction("Balance", menu)
    balance_action.triggered.connect(
        lambda _checked=False: _apply_amount_balance_from_row(
            table,
            _composition_editor_row(table, editor, row),
        )
    )
    menu.addAction(balance_action)
    return menu


def _install_composition_input_context_menu(
    editor: QLineEdit,
    table: QTableWidget | CompositionEditor,
    row: int,
) -> None:
    editor.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)

    def open_menu(position) -> None:
        menu = _composition_input_context_menu(
            editor,
            table,
            _composition_editor_row(table, editor, row),
        )
        menu.exec(editor.mapToGlobal(position))

    editor.customContextMenuRequested.connect(open_menu)


def _install_amount_balance_context_menu(
    editor: QLineEdit,
    table: QTableWidget | CompositionEditor,
    row: int,
) -> None:
    _install_composition_input_context_menu(editor, table, row)


def _refresh_composition_amount_placeholders(
    table: QTableWidget | CompositionEditor,
) -> None:
    placeholder = _composition_amount_placeholder(table)
    for row in range(table.rowCount()):
        widget = table.cellWidget(row, 1)
        if isinstance(widget, QLineEdit):
            widget.setPlaceholderText(placeholder)
        elif isinstance(widget, QWidget) and bool(
            widget.property("amount_single_editor")
        ):
            editor = _amount_editor_line_edit(widget)
            button = _amount_range_balance_button(widget)
            if editor is not None and (button is None or not button.isChecked()):
                editor.setPlaceholderText(placeholder)


def _refresh_composition_amount_editors(
    table: QTableWidget | CompositionEditor,
) -> None:
    use_range_editor = bool(table.property("allow_amount_ranges"))
    needs_balance_button = _fixed_total_for_table(table) is not None
    for row in range(table.rowCount()):
        text = _table_cell_text(table, row, 1)
        widget = table.cellWidget(row, 1)
        if use_range_editor:
            has_matching_editor = (
                isinstance(widget, QWidget)
                and bool(widget.property("amount_range_editor"))
                and (_amount_range_balance_button(widget) is not None)
                == needs_balance_button
            )
            if not has_matching_editor:
                _set_composition_amount_range_cell(table, row, text)
            else:
                _sync_amount_range_editor_state(table, row)
        else:
            has_matching_editor = (
                isinstance(widget, QWidget)
                and bool(widget.property("amount_single_editor"))
                and (_amount_range_balance_button(widget) is not None)
                == needs_balance_button
            )
            if not has_matching_editor:
                _set_composition_amount_single_cell(table, row, text)
            else:
                _sync_amount_single_editor_state(table, row)
    _refresh_composition_amount_placeholders(table)
    _ensure_required_balance_button(table)
    _fit_table_to_rows(table, min_rows=3)


def _set_composition_cell(
    table: QTableWidget | CompositionEditor,
    row: int,
    column: int,
    text: str,
    placeholder: str,
) -> None:
    if column == 1:
        _set_composition_amount_single_cell(table, row, text)
        return

    item = QTableWidgetItem(text)
    table.setItem(row, column, item)
    editor = QLineEdit(text)
    editor.setObjectName("CompositionCellInput")
    editor.setPlaceholderText(placeholder)
    editor.setMinimumHeight(COMPOSITION_INPUT_MIN_HEIGHT)

    def sync_item(value: str, target: QTableWidgetItem = item) -> None:
        if target.text() != value:
            target.setText(value)

    editor.textChanged.connect(sync_item)
    _install_composition_input_context_menu(editor, table, row)
    if column == 1:
        editor.editingFinished.connect(
            lambda r=row: _commit_balance_amount_expression(
                table,
                r,
            )
        )
    table.setCellWidget(row, column, editor)
    table.setRowHeight(row, COMPOSITION_ROW_HEIGHT)


def _set_composition_amount_single_cell(
    table: QTableWidget | CompositionEditor,
    row: int,
    text: str,
) -> None:
    fixed_total = _fixed_total_for_table(table)
    use_balance_button = fixed_total is not None
    previous_widget = table.cellWidget(row, 1)
    previous_normal = (
        _amount_widget_text_property(previous_widget, "normal_amount_text")
        if isinstance(previous_widget, QWidget)
        else ""
    )
    previous_balance_total = (
        _amount_widget_text_property(previous_widget, "balance_total_text")
        if isinstance(previous_widget, QWidget)
        else ""
    )

    item = table.item(row, 1)
    if item is None:
        item = QTableWidgetItem(text)
        table.setItem(row, 1, item)
    else:
        item.setText(text)

    container = QWidget()
    container.setObjectName("CompositionAmountInputContainer")
    container.setProperty("amount_single_editor", True)
    layout = QHBoxLayout(container)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.setSpacing(6)

    is_balance = text.strip().startswith("?")
    if not use_balance_button:
        editor = QLineEdit(text)
        editor.setObjectName("CompositionCellInput")
        editor.setPlaceholderText(_composition_amount_placeholder(table))
        editor.setMinimumHeight(COMPOSITION_INPUT_MIN_HEIGHT)
        layout.addWidget(editor, 1)

        def sync_item(value: str, target: QTableWidgetItem = item) -> None:
            if target.text() != value:
                target.setText(value)

        editor.textChanged.connect(sync_item)
        editor.editingFinished.connect(
            lambda r=row: _commit_balance_amount_expression(table, r)
        )
        _install_amount_balance_context_menu(editor, table, row)

        table.setCellWidget(row, 1, container)
        table.setRowHeight(row, COMPOSITION_ROW_HEIGHT)
        return

    balance_button = QToolButton(container)
    balance_button.setObjectName("BalanceAmountButton")
    balance_button.setText("Bal.")
    balance_button.setCheckable(True)
    balance_button.setToolTip("Balance species")
    layout.addWidget(balance_button)

    initial_balance_total = (
        text.strip()[1:].strip() if is_balance else previous_balance_total
    )
    initial_normal_text = previous_normal if is_balance else text
    editor_text = initial_balance_total if is_balance else initial_normal_text
    container.setProperty("normal_amount_text", initial_normal_text)
    container.setProperty("balance_total_text", initial_balance_total)
    editor = QLineEdit(editor_text)
    editor.setObjectName("CompositionCellInput")
    editor.setPlaceholderText(_composition_amount_placeholder(table))
    editor.setMinimumHeight(COMPOSITION_INPUT_MIN_HEIGHT)
    layout.addWidget(editor, 1)
    _install_amount_balance_context_menu(editor, table, row)

    def sync_item(target: QTableWidgetItem = item) -> None:
        total = _fixed_total_for_table(table)
        if balance_button.isChecked():
            if total is not None:
                combined = f"?{_balance_total_text(total)}"
            else:
                entered_total = _amount_widget_text_property(
                    container,
                    "balance_total_text",
                )
                combined = f"?{entered_total}" if entered_total else "?"
        else:
            combined = _amount_widget_text_property(container, "normal_amount_text")
        if target.text() != combined:
            target.setText(combined)

    def sync_editor_text(value: str) -> None:
        if balance_button.isChecked():
            container.setProperty("balance_total_text", value.strip())
        else:
            container.setProperty("normal_amount_text", value.strip())
        sync_item()

    def sync_balance_mode(checked: bool) -> None:
        total = _fixed_total_for_table(table)
        if checked:
            if not bool(container.property("initializing_amount_editor")):
                container.setProperty("normal_amount_text", editor.text().strip())
            _clear_other_balance_amount_buttons(table, row)
        elif (
            total is not None
            and bool(container.property("warn_on_balance_uncheck"))
            and not bool(table.property("clearing_balance_buttons"))
        ):
            if _checked_balance_button_count(table) == 0:
                QMessageBox.warning(
                    table.window(),
                    "Balance Species Required",
                    (
                        "A balance species is required for mol%, wt%, mole "
                        "fraction, and weight fraction amount units."
                    ),
                )
                previous = balance_button.blockSignals(True)
                balance_button.setChecked(True)
                balance_button.blockSignals(previous)
                checked = True

        if checked and total is not None:
            _set_line_edit_text_silent(editor, "")
            editor.setEnabled(False)
            editor.setPlaceholderText(_balance_placeholder(total))
        elif checked:
            balance_total = _amount_widget_text_property(
                container,
                "balance_total_text",
            )
            if not balance_total:
                balance_total = "100"
                container.setProperty("balance_total_text", balance_total)
            editor.setEnabled(True)
            _set_line_edit_text_silent(editor, balance_total)
            editor.setPlaceholderText("Total")
        else:
            normal_text = _amount_widget_text_property(container, "normal_amount_text")
            if not normal_text and not bool(table.property("clearing_balance_buttons")):
                normal_text = _amount_widget_text_property(
                    container,
                    "balance_total_text",
                )
                container.setProperty("normal_amount_text", normal_text)
            editor.setEnabled(True)
            _set_line_edit_text_silent(editor, normal_text)
            editor.setPlaceholderText(_composition_amount_placeholder(table))
        sync_item()

    editor.textChanged.connect(sync_editor_text)
    editor.editingFinished.connect(
        lambda r=row: _commit_balance_amount_expression(table, r)
    )
    balance_button.toggled.connect(sync_balance_mode)
    container.setProperty("initializing_amount_editor", True)
    balance_button.setChecked(is_balance)
    sync_balance_mode(is_balance)
    container.setProperty("initializing_amount_editor", False)
    container.setProperty("warn_on_balance_uncheck", True)

    table.setCellWidget(row, 1, container)
    table.setRowHeight(row, COMPOSITION_ROW_HEIGHT)


def _set_composition_amount_range_cell(
    table: QTableWidget | CompositionEditor,
    row: int,
    text: str,
) -> None:
    fixed_total = _fixed_total_for_table(table)
    use_balance_button = fixed_total is not None
    previous_widget = table.cellWidget(row, 1)
    previous_normal = (
        _amount_widget_text_property(previous_widget, "normal_amount_text")
        if isinstance(previous_widget, QWidget)
        else ""
    )
    previous_balance_total = (
        _amount_widget_text_property(previous_widget, "balance_total_text")
        if isinstance(previous_widget, QWidget)
        else ""
    )

    item = table.item(row, 1)
    if item is None:
        item = QTableWidgetItem(text)
        table.setItem(row, 1, item)
    else:
        item.setText(text)

    container = QWidget()
    container.setObjectName("CompositionRangeInputContainer")
    container.setProperty("amount_range_editor", True)
    layout = QHBoxLayout(container)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.setSpacing(6)

    is_balance = text.strip().startswith("?")
    balance_button: QToolButton | None = None
    if use_balance_button:
        balance_button = QToolButton(container)
        balance_button.setObjectName("BalanceAmountButton")
        balance_button.setText("Bal.")
        balance_button.setCheckable(True)
        balance_button.setToolTip("Balance species")
        layout.addWidget(balance_button)

    initial_balance_total = (
        text.strip()[1:].strip() if is_balance else previous_balance_total
    )
    initial_normal_text = previous_normal if is_balance else text
    container.setProperty("normal_amount_text", initial_normal_text)
    container.setProperty("balance_total_text", initial_balance_total)
    if is_balance and not use_balance_button:
        parts = [text.strip(), "", ""]
    else:
        parts = _amount_range_editor_parts(
            initial_balance_total if is_balance else initial_normal_text
        )
    placeholders = ("Initial", "Final", "Step")
    editors: list[QLineEdit] = []
    for index, placeholder in enumerate(placeholders):
        editor = QLineEdit(parts[index])
        editor.setObjectName("CompositionCellInput")
        editor.setPlaceholderText(placeholder)
        editor.setMinimumHeight(COMPOSITION_INPUT_MIN_HEIGHT)
        editors.append(editor)
        layout.addWidget(editor)

    if not use_balance_button:
        def sync_plain_range(
            _value: str = "",
            target: QTableWidgetItem = item,
        ) -> None:
            values = [editor.text().strip() for editor in editors]
            is_question_total = values[0].startswith("?")
            if is_question_total and _other_open_total_balance_row(
                table,
                row,
            ) is not None:
                _warn_duplicate_open_total_balance(table)
                values[0] = values[0][1:].strip()
                _set_line_edit_text_silent(editors[0], values[0])
                is_question_total = False
            if is_question_total:
                if not target.text().strip().startswith("?"):
                    container.setProperty("normal_amount_text", target.text().strip())
                container.setProperty(
                    "balance_total_text",
                    values[0][1:].strip(),
                )
                _set_open_total_range_editor_balance_state(editors, True)
                text_value = values[0]
            else:
                _set_open_total_range_editor_balance_state(editors, False)
                values = [editor.text().strip() for editor in editors]
                text_value = " ".join(value for value in values if value)
                container.setProperty("normal_amount_text", text_value)
            if target.text() != text_value:
                target.setText(text_value)

        for editor in editors:
            editor.textChanged.connect(sync_plain_range)
            editor.editingFinished.connect(
                lambda r=row: _commit_balance_amount_expression(table, r)
            )
            _install_amount_balance_context_menu(editor, table, row)

        table.setCellWidget(row, 1, container)
        table.setRowHeight(row, COMPOSITION_ROW_HEIGHT)
        _set_open_total_range_editor_balance_state(editors, is_balance)
        return

    def sync_item(target: QTableWidgetItem = item) -> None:
        total = _fixed_total_for_table(table)
        if balance_button is not None and balance_button.isChecked():
            if total is not None:
                combined = f"?{_balance_total_text(total)}"
            else:
                entered_total = _amount_widget_text_property(
                    container,
                    "balance_total_text",
                )
                combined = f"?{entered_total}" if entered_total else "?"
        else:
            combined = _amount_widget_text_property(container, "normal_amount_text")
        if target.text() != combined:
            target.setText(combined)

    def sync_editor_text() -> None:
        if balance_button.isChecked():
            container.setProperty("balance_total_text", editors[0].text().strip())
        else:
            values = [editor.text().strip() for editor in editors]
            container.setProperty(
                "normal_amount_text",
                " ".join(value for value in values if value),
            )
        sync_item()

    def sync_balance_mode(checked: bool) -> None:
        total = _fixed_total_for_table(table)
        if checked:
            if not bool(container.property("initializing_amount_editor")):
                values = [editor.text().strip() for editor in editors]
                container.setProperty(
                    "normal_amount_text",
                    " ".join(value for value in values if value),
                )
            _clear_other_balance_amount_buttons(table, row)
        elif (
            total is not None
            and bool(container.property("warn_on_balance_uncheck"))
            and not bool(table.property("clearing_balance_buttons"))
        ):
            if _checked_balance_button_count(table) == 0:
                QMessageBox.warning(
                    table.window(),
                    "Balance Species Required",
                    (
                        "A balance species is required for mol%, wt%, mole "
                        "fraction, and weight fraction amount units."
                    ),
                )
                previous = balance_button.blockSignals(True)
                balance_button.setChecked(True)
                balance_button.blockSignals(previous)
                checked = True
        if checked and total is not None:
            _set_line_edit_text_silent(editors[0], "")
            editors[0].setEnabled(False)
            editors[0].setPlaceholderText(_balance_placeholder(total))
        elif checked:
            balance_total = _amount_widget_text_property(
                container,
                "balance_total_text",
            )
            if not balance_total:
                balance_total = "100"
                container.setProperty("balance_total_text", balance_total)
            editors[0].setEnabled(True)
            _set_line_edit_text_silent(editors[0], balance_total)
            editors[0].setPlaceholderText("Total")
        else:
            normal_text = _amount_widget_text_property(container, "normal_amount_text")
            if not normal_text and not bool(table.property("clearing_balance_buttons")):
                normal_text = _amount_widget_text_property(
                    container,
                    "balance_total_text",
                )
                container.setProperty("normal_amount_text", normal_text)
            normal_parts = _amount_range_editor_parts(
                normal_text
            )
            for editor, part in zip(editors, normal_parts, strict=False):
                _set_line_edit_text_silent(editor, part)
            editors[0].setEnabled(True)
            editors[0].setPlaceholderText("Initial")
        editors[1].setVisible(not checked)
        editors[2].setVisible(not checked)
        sync_item()

    for editor in editors:
        editor.textChanged.connect(lambda _value, sync=sync_editor_text: sync())
    balance_button.toggled.connect(sync_balance_mode)
    container.setProperty("initializing_amount_editor", True)
    balance_button.setChecked(is_balance)
    sync_balance_mode(is_balance)
    container.setProperty("initializing_amount_editor", False)
    container.setProperty("warn_on_balance_uncheck", True)

    table.setCellWidget(row, 1, container)
    table.setRowHeight(row, COMPOSITION_ROW_HEIGHT)


def _amount_range_editor_parts(text: str) -> list[str]:
    cleaned = text.strip().strip("[]")
    if not cleaned:
        return ["", "", ""]
    if cleaned.startswith("?"):
        return [cleaned[1:].strip(), "", ""]
    parts = [part for part in re.split(r"[\s,]+", cleaned) if part]
    if len(parts) == 1:
        return [parts[0], parts[0], ""]
    return (parts + ["", "", ""])[:3]


def _table_cell_text(
    table: QTableWidget | CompositionEditor,
    row: int,
    column: int,
) -> str:
    widget = table.cellWidget(row, column)
    if isinstance(widget, QLineEdit):
        return widget.text().strip()
    if isinstance(widget, QWidget) and bool(widget.property("amount_single_editor")):
        balance_button = _amount_range_balance_button(widget)
        editor = _amount_editor_line_edit(widget)
        if balance_button is not None and balance_button.isChecked():
            total = _fixed_total_for_table(table)
            if total is not None:
                return f"?{_balance_total_text(total)}"
            total_text = editor.text().strip() if editor is not None else ""
            return f"?{total_text}" if total_text else "?"
        return editor.text().strip() if editor is not None else ""
    if isinstance(widget, QWidget) and bool(widget.property("amount_range_editor")):
        balance_button = _amount_range_balance_button(widget)
        editors = widget.findChildren(QLineEdit)
        if balance_button is not None and balance_button.isChecked():
            total = _fixed_total_for_table(table)
            if total is not None:
                return f"?{_balance_total_text(total)}"
            total_text = editors[0].text().strip() if editors else ""
            return f"?{total_text}" if total_text else "?"
        values = [editor.text().strip() for editor in editors]
        return " ".join(value for value in values if value)
    item = table.item(row, column)
    return item.text().strip() if item is not None else ""


def _set_table_cell_text(
    table: QTableWidget | CompositionEditor,
    row: int,
    column: int,
    text: str,
) -> None:
    widget = table.cellWidget(row, column)
    if isinstance(widget, QLineEdit):
        widget.setText(text)
        return
    if isinstance(widget, QWidget) and bool(widget.property("amount_single_editor")):
        balance_button = _amount_range_balance_button(widget)
        editor = _amount_editor_line_edit(widget)
        if balance_button is None:
            if editor is not None:
                editor.setText(text)
            _set_table_item_text_only(table, row, column, text)
            return
        is_balance = text.strip().startswith("?")
        if editor is not None:
            if is_balance:
                widget.setProperty("balance_total_text", text.strip()[1:].strip())
                editor.setText(
                    _amount_widget_text_property(widget, "balance_total_text")
                )
            else:
                widget.setProperty("normal_amount_text", text.strip())
                editor.setText(text)
        if balance_button is not None:
            balance_button.setChecked(is_balance)
        return
    if isinstance(widget, QWidget) and bool(widget.property("amount_range_editor")):
        balance_button = _amount_range_balance_button(widget)
        editors = widget.findChildren(QLineEdit)
        is_balance = text.strip().startswith("?")
        parts = [text.strip(), "", ""] if balance_button is None and is_balance else (
            _amount_range_editor_parts(text)
        )
        if is_balance:
            widget.setProperty("balance_total_text", text.strip()[1:].strip())
        else:
            widget.setProperty("normal_amount_text", text.strip())
        for editor, part in zip(editors, parts, strict=False):
            editor.setText(part)
        if balance_button is None:
            _set_open_total_range_editor_balance_state(editors, is_balance)
        if balance_button is not None:
            balance_button.setChecked(is_balance)
        return
    item = table.item(row, column)
    if item is not None:
        item.setText(text)


def _amount_range_balance_button(widget: QWidget) -> QToolButton | None:
    for button in widget.findChildren(QToolButton):
        if button.objectName() == "BalanceAmountButton":
            return button
    return None


def _amount_editor_line_edit(widget: QWidget) -> QLineEdit | None:
    editors = widget.findChildren(QLineEdit)
    return editors[0] if editors else None


def _sync_amount_single_editor_state(
    table: QTableWidget | CompositionEditor,
    row: int,
) -> None:
    widget = table.cellWidget(row, 1)
    if not isinstance(widget, QWidget) or not bool(
        widget.property("amount_single_editor")
    ):
        return
    button = _amount_range_balance_button(widget)
    editor = _amount_editor_line_edit(widget)
    if button is None or editor is None:
        return
    total = _fixed_total_for_table(table)
    if button.isChecked() and total is not None:
        if editor.isEnabled() and not _amount_widget_text_property(
            widget,
            "normal_amount_text",
        ):
            widget.setProperty("normal_amount_text", editor.text().strip())
        _set_line_edit_text_silent(editor, "")
        editor.setEnabled(False)
        editor.setPlaceholderText(_balance_placeholder(total))
        _set_table_item_text_only(table, row, 1, f"?{_balance_total_text(total)}")
    else:
        editor.setEnabled(True)
        editor.setPlaceholderText(
            "Total" if button.isChecked() else _composition_amount_placeholder(table)
        )


def _sync_amount_range_editor_state(
    table: QTableWidget | CompositionEditor,
    row: int,
) -> None:
    widget = table.cellWidget(row, 1)
    if not isinstance(widget, QWidget) or not bool(
        widget.property("amount_range_editor")
    ):
        return
    button = _amount_range_balance_button(widget)
    editors = widget.findChildren(QLineEdit)
    if button is None or not editors:
        return
    total = _fixed_total_for_table(table)
    if button.isChecked() and total is not None:
        if editors[1].isVisible() and not _amount_widget_text_property(
            widget,
            "normal_amount_text",
        ):
            values = [editor.text().strip() for editor in editors]
            widget.setProperty(
                "normal_amount_text",
                " ".join(value for value in values if value),
            )
        _set_line_edit_text_silent(editors[0], "")
        editors[0].setEnabled(False)
        editors[0].setPlaceholderText(_balance_placeholder(total))
        editors[1].setVisible(False)
        editors[2].setVisible(False)
        _set_table_item_text_only(table, row, 1, f"?{_balance_total_text(total)}")
    else:
        editors[0].setEnabled(True)
        editors[0].setPlaceholderText(
            "Total" if button.isChecked() else "Initial"
        )
        editors[1].setVisible(not button.isChecked())
        editors[2].setVisible(not button.isChecked())


def _set_table_item_text_only(
    table: QTableWidget | CompositionEditor,
    row: int,
    column: int,
    text: str,
) -> None:
    item = table.item(row, column)
    if item is None:
        table.setItem(row, column, QTableWidgetItem(text))
    elif item.text() != text:
        item.setText(text)


def _checked_balance_button_count(table: QTableWidget | CompositionEditor) -> int:
    count = 0
    for row in range(table.rowCount()):
        widget = table.cellWidget(row, 1)
        if not isinstance(widget, QWidget):
            continue
        button = _amount_range_balance_button(widget)
        if button is not None and button.isChecked():
            count += 1
    return count


def _ensure_required_balance_button(table: QTableWidget | CompositionEditor) -> None:
    total = _fixed_total_for_table(table)
    if total is None or _checked_balance_button_count(table) > 0:
        return
    for row in range(table.rowCount()):
        widget = table.cellWidget(row, 1)
        if not isinstance(widget, QWidget):
            continue
        button = _amount_range_balance_button(widget)
        if button is None:
            continue
        previous = button.blockSignals(True)
        button.setChecked(True)
        button.blockSignals(previous)
        if bool(widget.property("amount_single_editor")):
            _sync_amount_single_editor_state(table, row)
        else:
            _sync_amount_range_editor_state(table, row)
        return


def _clear_other_balance_amount_buttons(
    table: QTableWidget | CompositionEditor,
    active_row: int,
) -> None:
    previous_clearing = bool(table.property("clearing_balance_buttons"))
    table.setProperty("clearing_balance_buttons", True)
    try:
        for row in range(table.rowCount()):
            if row == active_row:
                continue
            widget = table.cellWidget(row, 1)
            if not isinstance(widget, QWidget):
                continue
            button = _amount_range_balance_button(widget)
            if button is not None and button.isChecked():
                button.setChecked(False)
    finally:
        table.setProperty("clearing_balance_buttons", previous_clearing)


def _commit_balance_amount_expression(
    table: QTableWidget | CompositionEditor,
    row: int,
) -> None:
    if bool(table.property("allow_amount_ranges")):
        return
    if row < 0 or row >= table.rowCount():
        return
    amount_text = _table_cell_text(table, row, 1)
    if not amount_text.startswith("?"):
        return
    amount_values, amount_errors = _resolved_composition_amounts(table)
    if amount_errors or row not in amount_values:
        return
    _set_table_cell_text(table, row, 1, _format_amount_value(amount_values[row]))


def _set_table_cell_error(
    table: QTableWidget | CompositionEditor,
    row: int,
    column: int,
    message: str,
    *,
    warning: bool = False,
) -> None:
    item = table.item(row, column)
    if item is not None:
        _set_table_item_error(item, message, warning=warning)
    widget = table.cellWidget(row, column)
    if isinstance(widget, QLineEdit):
        if message:
            background = "#8a5a18" if warning else "#6f2429"
            foreground = "#fff4de" if warning else "#ffffff"
            widget.setStyleSheet(
                "QLineEdit { "
                f"background: {background}; "
                f"color: {foreground}; "
                "border: 1px solid #566063; "
                "border-radius: 6px; "
                "padding: 6px 10px; "
                f"min-height: {COMPOSITION_INPUT_MIN_HEIGHT}px; "
                "}"
            )
            widget.setToolTip(message)
        else:
            widget.setStyleSheet("")
            widget.setToolTip("")
        return
    if isinstance(widget, QWidget) and bool(widget.property("amount_range_editor")):
        for editor in widget.findChildren(QLineEdit):
            if message:
                background = "#8a5a18" if warning else "#6f2429"
                foreground = "#fff4de" if warning else "#ffffff"
                editor.setStyleSheet(
                    "QLineEdit { "
                    f"background: {background}; "
                    f"color: {foreground}; "
                    "border: 1px solid #566063; "
                    "border-radius: 6px; "
                    "padding: 6px 10px; "
                    f"min-height: {COMPOSITION_INPUT_MIN_HEIGHT}px; "
                    "}"
                )
                editor.setToolTip(message)
            else:
                editor.setStyleSheet("")
                editor.setToolTip("")


def _remove_selected_table_rows(table: QTableWidget | CompositionEditor) -> None:
    rows = sorted({index.row() for index in table.selectedIndexes()}, reverse=True)
    if not rows:
        return
    next_row = min(rows)
    for row in rows:
        table.removeRow(row)
    if table.rowCount() == 0:
        _add_composition_row(table)
        next_row = 0
    else:
        next_row = min(next_row, table.rowCount() - 1)
    table.clearSelection()
    table.selectRow(next_row)
    table.setCurrentCell(next_row, 0)
    _fit_table_to_rows(table, min_rows=3)


def _validate_composition_items(
    table: QTableWidget,
    database: Any | None,
) -> list[str]:
    diagnostics: list[str] = []
    previous_blocked = table.blockSignals(True)
    amount_values, amount_errors = _resolved_composition_amounts(
        table,
        update_table=False,
        allow_ranges=bool(table.property("allow_amount_ranges")),
    )
    for row in range(table.rowCount()):
        species = _table_cell_text(table, row, 0)
        amount_text = _table_cell_text(table, row, 1)
        if not species and not amount_text:
            _set_table_cell_error(table, row, 0, "")
            _set_table_cell_error(table, row, 1, "")
            continue

        species_error = ""
        amount_error = amount_errors.get(row, "")
        if not species:
            species_error = "Species name cannot be empty."
            species_warning = False
        else:
            _composition, warnings = _formula_or_database_species_stoichiometry(
                species,
                database,
            )
            species_warning = all(
                warning.startswith("Element(s) not found in selected database")
                for warning in warnings
            )
            if warnings:
                species_error = "; ".join(warnings)

        parsed_amount = amount_values.get(row)
        if parsed_amount is not None and parsed_amount < 0:
            amount_error = "Amount must be non-negative."

        _set_table_cell_error(
            table,
            row,
            0,
            species_error,
            warning=species_warning,
        )
        _set_table_cell_error(table, row, 1, amount_error)
        if species_error:
            diagnostics.append(f"row {row + 1}: {species_error}")
        if amount_error:
            diagnostics.append(f"row {row + 1}: {amount_error}")
    table.blockSignals(previous_blocked)
    return diagnostics


def _set_table_item_error(
    item: QTableWidgetItem,
    message: str,
    *,
    warning: bool = False,
) -> None:
    if message:
        item.setBackground(QBrush(QColor("#8a5a18" if warning else "#6f2429")))
        item.setForeground(QBrush(QColor("#fff4de" if warning else "#ffffff")))
        item.setToolTip(message)
        return
    item.setData(Qt.ItemDataRole.BackgroundRole, None)
    item.setData(Qt.ItemDataRole.ForegroundRole, None)
    item.setToolTip("")


def _resolved_composition_amounts(
    table: QTableWidget | CompositionEditor,
    *,
    update_table: bool = False,
    allow_ranges: bool = False,
    amount_unit: str | None = None,
) -> tuple[dict[int, float], dict[int, str]]:
    values: dict[int, float] = {}
    errors: dict[int, str] = {}
    balance_rows: list[tuple[int, float]] = []
    known_total = 0.0
    fixed_total = (
        _amount_unit_fixed_total(amount_unit)
        if amount_unit is not None
        else _fixed_total_for_table(table)
    )
    first_used_row: int | None = None

    for row in range(table.rowCount()):
        species_text = _table_cell_text(table, row, 0)
        amount_text = _table_cell_text(table, row, 1)
        if not species_text and not amount_text:
            continue
        if first_used_row is None:
            first_used_row = row
        if amount_text.startswith("?"):
            target_text = amount_text[1:].strip()
            if fixed_total is not None:
                target_total = fixed_total
            elif not target_text:
                errors[row] = "Enter a total amount after '?'."
                continue
            else:
                try:
                    target_total = float(target_text)
                except ValueError:
                    errors[row] = "Amount total after '?' must be numeric."
                    continue
                if target_total <= 0:
                    errors[row] = "Total amount after '?' must be positive."
                    continue
            balance_rows.append((row, target_total))
            continue

        try:
            if allow_ranges:
                amount_values = _amount_condition_values(amount_text)
                amount_value = amount_values[0]
            else:
                amount_value = float(amount_text)
        except ValueError as exc:
            errors[row] = str(exc) if allow_ranges else "Amount must be numeric."
            continue
        values[row] = amount_value
        known_total += amount_value

    if len(balance_rows) > 1:
        for row, _target_total in balance_rows:
            errors[row] = "Only one balance species can be used."
    elif fixed_total is not None and not balance_rows and first_used_row is not None:
        errors[first_used_row] = (
            "Select one balance species for mol%, wt%, mole fraction, "
            "or weight fraction."
        )
    elif len(balance_rows) == 1 and not errors:
        row, target_total = balance_rows[0]
        residual = target_total - known_total
        if residual <= 0:
            errors[row] = "The total amount is smaller than the sum."
        else:
            values[row] = residual
            if update_table and not allow_ranges:
                _set_table_cell_text(table, row, 1, _format_amount_value(residual))

    return values, errors


def _format_amount_value(value: float) -> str:
    return f"{value:.12g}"


def _composition_elements(
    table: QTableWidget,
    database: Any | None = None,
    *,
    amount_unit: str | None = None,
    active_only: bool = False,
    allow_ranges: bool | None = None,
) -> list[str]:
    if active_only:
        if allow_ranges is None:
            allow_ranges = bool(table.property("allow_amount_ranges"))
        if amount_unit is None:
            amount_unit = str(table.property("amount_unit") or "moles")
        if allow_ranges:
            return _composition_elements_from_amount_ranges(table, database)
        composition = _expanded_element_composition(table, amount_unit, database)
        return [
            element
            for element, amount in composition.items()
            if not np.isclose(float(amount), 0.0)
        ]

    elements: list[str] = []
    element_map = _database_element_map(database)
    for row in range(table.rowCount()):
        name = _table_cell_text(table, row, 0)
        if not name:
            continue
        composition, warnings = _formula_or_database_species_stoichiometry(
            name,
            database,
        )
        if warnings:
            continue
        for element in composition:
            db_element = element_map.get(element, element)
            if db_element not in elements:
                elements.append(db_element)
    return elements


def _composition_elements_from_amount_ranges(
    table: QTableWidget,
    database: Any | None,
) -> list[str]:
    elements: list[str] = []
    element_map = _database_element_map(database)
    for row in range(table.rowCount()):
        species = _table_cell_text(table, row, 0)
        amount_text = _table_cell_text(table, row, 1)
        if not species and not amount_text:
            continue
        if not species:
            raise ValueError(f"row {row + 1}: species name cannot be empty")
        if not amount_text:
            raise ValueError(f"row {row + 1}: amount cannot be empty")
        if not _amount_text_may_be_nonzero(amount_text):
            continue

        composition, warnings = _formula_or_database_species_stoichiometry(
            species,
            database,
        )
        if warnings:
            raise ValueError("; ".join(warnings))
        for element in composition:
            db_element = element_map.get(element, element)
            if db_element not in elements:
                elements.append(db_element)
    return elements


def _amount_text_may_be_nonzero(text: str) -> bool:
    cleaned = text.strip()
    if cleaned.startswith("?"):
        return True
    values = _amount_condition_values(cleaned)
    return any(not np.isclose(value, 0.0) for value in values)


def _active_condition_elements(condition: dict[str, Any]) -> list[str]:
    elements: list[str] = []
    for key, values in list(condition.items())[2:]:
        array = np.asarray(values, dtype=float)
        if np.any(~np.isclose(array, 0.0)):
            elements.append(key)
    return elements


def _module_condition_elements(
    module: CalculationModule | None,
    composition_table: QTableWidget,
    database: Any | None,
) -> list[str]:
    if (
        module is not None
        and module.kind in {"equilibrium", "solidification"}
        and module.calculation_type == "batch"
    ):
        if module.batch_condition:
            expanded = _expand_condition_species(
                database,
                module.batch_condition,
                module.amount_unit,
            )
            return _active_condition_elements(expanded)
        return _composition_elements(
            composition_table,
            database,
            amount_unit=module.amount_unit,
            active_only=True,
            allow_ranges=True,
        )
    return _composition_elements(
        composition_table,
        database,
        amount_unit=module.amount_unit if module is not None else None,
        active_only=True,
    )


def _calculation_condition(
    table: QTableWidget,
    temperature: float,
    pressure: float,
    amount_unit: str = "moles",
    database: Any | None = None,
) -> dict[str, float]:
    condition: dict[str, float] = {"T": temperature, "P": pressure}
    composition = {
        element: amount
        for element, amount in _expanded_element_composition(
            table,
            amount_unit,
            database,
        ).items()
        if not np.isclose(float(amount), 0.0)
    }
    if not composition:
        raise ValueError("composition is empty")
    condition.update(composition)
    return condition


def _manual_equilibrium_condition(
    database: Any,
    table: QTableWidget,
    temperature_text: str,
    pressure_text: str,
    units: list[str],
    *,
    phases: list[str] | None,
    include_transitions: bool,
) -> dict[str, list[float]]:
    temperatures = _temperature_condition_values(temperature_text)
    try:
        pressure = float(pressure_text)
    except ValueError as exc:
        raise ValueError("pressure must be numeric") from exc

    base_condition = _calculation_condition(
        table,
        temperatures[0],
        pressure,
        units[2],
        database,
    )
    if include_transitions:
        if len(temperatures) < 2:
            raise ValueError(
                "transition search requires a temperature range: initial final step"
            )
        import equilipy as eq

        t_start = temperatures[0]
        t_stop = temperatures[-1]
        transition_temperatures = eq.find_transitions(
            database,
            dict(base_condition),
            max(t_start, t_stop),
            min(t_start, t_stop),
            unit=units,
            phases=phases,
        )
        temperatures = _merge_temperature_values(
            temperatures,
            transition_temperatures,
        )
    return _condition_for_temperatures(base_condition, temperatures)


def _temperature_condition_values(text: str) -> list[float]:
    cleaned = text.strip().strip("[]")
    if not cleaned:
        raise ValueError("temperature cannot be empty")
    parts = [part for part in re.split(r"[\s,]+", cleaned) if part]
    try:
        values = [float(part) for part in parts]
    except ValueError as exc:
        raise ValueError("temperature values must be numeric") from exc

    if len(values) == 1:
        return values
    if len(values) != 3:
        raise ValueError("temperature range must be: initial final step")

    start, stop, step = values
    if step == 0:
        raise ValueError("temperature step cannot be zero")
    if start == stop:
        return [start]
    if (stop - start) * step < 0:
        raise ValueError("temperature step must move from initial toward final")

    temperatures = [start]
    current = start + step
    if step > 0:
        while current < stop:
            temperatures.append(current)
            current += step
    else:
        while current > stop:
            temperatures.append(current)
            current += step
    if not np.isclose(temperatures[-1], stop):
        temperatures.append(stop)
    return temperatures


def _temperature_text_is_range(text: str) -> bool:
    cleaned = text.strip().strip("[]")
    if not cleaned:
        return False
    parts = [part for part in re.split(r"[\s,]+", cleaned) if part]
    if len(parts) != 3:
        return False
    try:
        start, stop, step = (float(part) for part in parts)
    except ValueError:
        return False
    return step != 0 and start != stop and (stop - start) * step > 0


def _amount_condition_values(text: str) -> list[float]:
    cleaned = text.strip().strip("[]")
    if not cleaned:
        raise ValueError("Amount cannot be empty.")
    if cleaned.startswith("?"):
        raise ValueError("'?' balance amount cannot also be a range.")
    parts = [part for part in re.split(r"[\s,]+", cleaned) if part]
    try:
        values = [float(part) for part in parts]
    except ValueError as exc:
        raise ValueError("Amount values must be numeric.") from exc

    if len(values) == 1:
        if values[0] < 0:
            raise ValueError("Amount must be non-negative.")
        return values
    if len(values) != 3:
        raise ValueError("Amount range must be: initial final step.")

    start, stop, step = values
    if start < 0 or stop < 0:
        raise ValueError("Amount range values must be non-negative.")
    if step == 0:
        return [start]
    if start == stop:
        return [start]
    if (stop - start) * step < 0:
        raise ValueError("Amount step must move from initial toward final.")

    amounts = [start]
    current = start + step
    if step > 0:
        while current < stop:
            amounts.append(current)
            current += step
    else:
        while current > stop:
            amounts.append(current)
            current += step
    if not np.isclose(amounts[-1], stop):
        amounts.append(stop)
    return amounts


def _manual_batch_point_count(
    table: QTableWidget | CompositionEditor,
    temperature_text: str,
    pressure_text: str,
) -> int:
    temperatures = _temperature_condition_values(temperature_text)
    try:
        float(pressure_text)
    except ValueError as exc:
        raise ValueError("pressure must be numeric") from exc
    amount_grid_count = 1
    used_rows = 0
    balance_count = 0
    fixed_total = _fixed_total_for_table(table)
    for row in range(table.rowCount()):
        species = _table_cell_text(table, row, 0)
        amount_text = _table_cell_text(table, row, 1)
        if not species and not amount_text:
            continue
        if not species:
            raise ValueError(f"row {row + 1} species is empty")
        if amount_text.startswith("?"):
            balance_count += 1
            if fixed_total is None:
                target_text = amount_text[1:].strip()
                if not target_text:
                    raise ValueError(f"row {row + 1} balance total is empty")
                try:
                    target_total = float(target_text)
                except ValueError as exc:
                    raise ValueError(
                        f"row {row + 1} balance total must be numeric"
                    ) from exc
                if target_total <= 0:
                    raise ValueError(f"row {row + 1} balance total must be positive")
            used_rows += 1
            continue
        amount_grid_count *= len(_amount_condition_values(amount_text))
        used_rows += 1
    if used_rows == 0:
        raise ValueError("enter at least one species")
    if balance_count > 1:
        raise ValueError("only one balance species can be used")
    if fixed_total is not None and balance_count == 0:
        raise ValueError(
            "select one balance species for mol%, wt%, mole fraction, "
            "or weight fraction"
        )
    return len(temperatures) * amount_grid_count


def _manual_batch_condition(
    database: Any,
    table: QTableWidget | CompositionEditor,
    temperature_text: str,
    pressure_text: str,
    units: list[str],
) -> dict[str, list[float]]:
    temperatures = _temperature_condition_values(temperature_text)
    try:
        pressure = float(pressure_text)
    except ValueError as exc:
        raise ValueError("pressure must be numeric") from exc

    species_order: list[str] = []
    ranged_rows: list[tuple[str, list[float]]] = []
    balance_row: tuple[str, float] | None = None
    fixed_total = _amount_unit_fixed_total(units[2])
    for row in range(table.rowCount()):
        species = _table_cell_text(table, row, 0)
        amount_text = _table_cell_text(table, row, 1)
        if not species and not amount_text:
            continue
        if not species:
            raise ValueError(f"row {row + 1}: species name cannot be empty")
        if species not in species_order:
            species_order.append(species)
        if amount_text.startswith("?"):
            if balance_row is not None:
                raise ValueError("Only one balance species can be used.")
            if fixed_total is not None:
                target_total = fixed_total
            else:
                target_text = amount_text[1:].strip()
                if not target_text:
                    raise ValueError("Enter a total amount after '?'.")
                try:
                    target_total = float(target_text)
                except ValueError as exc:
                    raise ValueError("Amount total after '?' must be numeric.") from exc
                if target_total <= 0:
                    raise ValueError("Total amount after '?' must be positive.")
            balance_row = (species, target_total)
            continue
        ranged_rows.append((species, _amount_condition_values(amount_text)))

    if not species_order:
        raise ValueError("composition is empty")
    if fixed_total is not None and balance_row is None:
        raise ValueError(
            "Select one balance species for mol%, wt%, mole fraction, "
            "or weight fraction."
        )

    amount_combinations = product(
        *(values for _species, values in ranged_rows),
    )
    raw_condition: dict[str, list[float]] = {"T": [], "P": []}
    for species in species_order:
        raw_condition[species] = []

    for temperature_value in temperatures:
        for amounts in amount_combinations:
            row_amounts: dict[str, float] = {}
            known_total = 0.0
            for (species, _values), amount in zip(
                ranged_rows,
                amounts,
                strict=False,
            ):
                row_amounts[species] = row_amounts.get(species, 0.0) + amount
                known_total += amount
            if balance_row is not None:
                balance_species, target_total = balance_row
                residual = target_total - known_total
                if residual <= 0:
                    raise ValueError("The total amount is smaller than the sum.")
                row_amounts[balance_species] = (
                    row_amounts.get(balance_species, 0.0) + residual
                )
            raw_condition["T"].append(temperature_value)
            raw_condition["P"].append(pressure)
            for species in species_order:
                raw_condition[species].append(row_amounts.get(species, 0.0))
        amount_combinations = product(
            *(values for _species, values in ranged_rows),
        )

    expanded = _expand_condition_species(database, raw_condition, units[2])
    return {
        key: [float(value) for value in np.asarray(values, dtype=float).tolist()]
        for key, values in expanded.items()
    }


def _merge_temperature_values(
    base_temperatures: list[float],
    extra_temperatures: Any,
) -> list[float]:
    lower = min(base_temperatures[0], base_temperatures[-1])
    upper = max(base_temperatures[0], base_temperatures[-1])
    merged = list(base_temperatures)
    for value in np.asarray(extra_temperatures, dtype=float).flatten():
        temperature = float(value)
        if lower <= temperature <= upper:
            merged.append(temperature)

    reverse = base_temperatures[-1] < base_temperatures[0]
    merged = sorted(merged, reverse=reverse)
    unique: list[float] = []
    for value in merged:
        if not any(np.isclose(value, existing) for existing in unique):
            unique.append(value)
    return unique


def _condition_for_temperatures(
    base_condition: dict[str, float],
    temperatures: list[float],
) -> dict[str, list[float]]:
    condition: dict[str, list[float]] = {"T": temperatures}
    row_count = len(temperatures)
    for key, value in base_condition.items():
        if key == "T":
            continue
        condition[key] = [float(value)] * row_count
    return condition


def _read_condition_csv(path: str) -> dict[str, list[float]]:
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        try:
            headers = [header.strip() for header in next(reader)]
        except StopIteration as exc:
            raise ValueError("Conditions CSV is empty.") from exc

        if len(headers) < 3 or headers[:2] != ["T", "P"]:
            raise ValueError("Conditions CSV header must be T, P, Species...")
        if any(not header for header in headers):
            raise ValueError("Conditions CSV contains an empty header.")
        if len(set(headers)) != len(headers):
            raise ValueError("Conditions CSV contains duplicate headers.")

        condition = {header: [] for header in headers}
        for line_number, row in enumerate(reader, start=2):
            if not row or all(not value.strip() for value in row):
                continue
            if len(row) != len(headers):
                raise ValueError(
                    f"Line {line_number}: expected {len(headers)} values, "
                    f"found {len(row)}."
                )
            for header, value in zip(headers, row, strict=False):
                text = value.strip()
                if not text:
                    raise ValueError(f"Line {line_number}: missing value for {header}.")
                try:
                    condition[header].append(float(text))
                except ValueError as exc:
                    raise ValueError(
                        f"Line {line_number}: {header} must be numeric."
                    ) from exc

    if not condition["T"]:
        raise ValueError("Conditions CSV has no data rows.")
    return condition


def _populate_batch_condition_table(
    table: QTableWidget,
    condition: dict[str, list[float]],
) -> None:
    headers = list(condition)
    row_count = _condition_row_count(condition)
    table.clear()
    table.setColumnCount(len(headers))
    table.setHorizontalHeaderLabels(headers)
    table.setRowCount(row_count)
    for row in range(row_count):
        for column, header in enumerate(headers):
            value = condition[header][row]
            item = QTableWidgetItem(_format_amount_value(value))
            item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            table.setItem(row, column, item)
    table.horizontalHeader().setStretchLastSection(True)
    _disable_inner_scrollbars(table)
    _fit_table_to_rows(table, min_rows=1)
    table.setVisible(True)


def _condition_row_count(condition: dict[str, list[float]]) -> int:
    if not condition:
        return 0
    lengths = {len(values) for values in condition.values()}
    if len(lengths) != 1:
        raise ValueError("condition columns must have the same length")
    return lengths.pop()


def _expanded_element_composition(
    table: QTableWidget,
    amount_unit: str,
    database: Any | None,
) -> dict[str, float]:
    composition: dict[str, float] = {}
    element_map = _database_element_map(database)
    atomic_masses = _database_atomic_mass_map(database)
    amount_values, amount_errors = _resolved_composition_amounts(
        table,
        amount_unit=amount_unit,
    )
    if amount_errors:
        first_row = min(amount_errors)
        raise ValueError(amount_errors[first_row])
    for row in range(table.rowCount()):
        species = _table_cell_text(table, row, 0)
        if not species:
            continue
        amount_value = amount_values.get(row, 0.0)
        if amount_value < 0:
            raise ValueError("composition amount must be non-negative")
        stoichiometry, warnings = _formula_or_database_species_stoichiometry(
            species,
            database,
        )
        if warnings:
            raise ValueError("; ".join(warnings))
        if _is_mass_amount_unit(amount_unit):
            molar_mass = _molar_mass(stoichiometry, atomic_masses)
            if molar_mass <= 0:
                raise ValueError(f"molar mass is unavailable for {species}")
            contributions = {
                element: (
                    amount_value
                    * coefficient
                    * atomic_masses[element]
                    / molar_mass
                )
                for element, coefficient in stoichiometry.items()
            }
        else:
            contributions = {
                element: amount_value * coefficient
                for element, coefficient in stoichiometry.items()
            }
        for element, element_amount in contributions.items():
            db_element = element_map.get(element, element)
            composition[db_element] = composition.get(db_element, 0.0) + element_amount
    return composition


def _formula_or_database_species_stoichiometry(
    formula: str,
    database: Any | None,
) -> tuple[dict[str, float], list[str]]:
    try:
        return _shared_formula_or_database_species_stoichiometry(formula, database), []
    except ValueError as exc:
        return {}, [str(exc)]


def _database_element_map(database: Any | None) -> dict[str, str]:
    if not isinstance(database, dict):
        return {}
    elements = database.get("cElementNameCS")
    if elements is None:
        return {}
    element_map: dict[str, str] = {}
    for element in elements:
        raw = str(element).strip()
        if not raw:
            continue
        element_map[_canonical_element_symbol(raw)] = raw
    return element_map


def _database_atomic_mass_map(database: Any | None) -> dict[str, float]:
    if not isinstance(database, dict):
        return {}
    element_map = _database_element_map(database)
    masses = database.get("dAtomicMass")
    if masses is None:
        return {}
    atomic_masses: dict[str, float] = {}
    for element, mass in zip(element_map, masses, strict=False):
        atomic_masses[element] = float(mass)
    return atomic_masses


def _canonical_element_symbol(symbol: str) -> str:
    stripped = symbol.strip()
    if stripped.lower() == "e-":
        return "e-"
    if not stripped:
        return ""
    return stripped[:1].upper() + stripped[1:].lower()


def _is_mass_amount_unit(amount_unit: str) -> bool:
    unit = amount_unit.strip().lower()
    return unit in {"g", "kg", "lbs", "wt%", "wt.%"} or unit.startswith(
        ("gram", "kilo", "pound", "mass fraction", "weight fraction")
    )


def _batch_cpu_count() -> int:
    return max(1, (os.cpu_count() or 1) - 2)


def _max_cpu_count() -> int:
    return max(1, os.cpu_count() or 1)


def _normalized_batch_cpu_count(value: Any) -> int:
    try:
        count = int(value)
    except (TypeError, ValueError):
        return 0
    return max(0, min(count, _max_cpu_count()))


def _populate_batch_cpu_combo(combo: QComboBox, selected: int = 0) -> None:
    combo.clear()
    auto_count = _batch_cpu_count()
    combo.addItem(f"Auto ({auto_count})", 0)
    for count in range(1, _max_cpu_count() + 1):
        combo.addItem(str(count), count)
    selected = _normalized_batch_cpu_count(selected)
    index = combo.findData(selected)
    combo.setCurrentIndex(index if index >= 0 else 0)


def _batch_cpu_combo_value(combo: QComboBox) -> int:
    return _normalized_batch_cpu_count(combo.currentData())


def _condition_row_count(condition: dict[str, Any]) -> int:
    row_count = 1
    for value in condition.values():
        if isinstance(value, np.ndarray):
            if value.ndim > 0:
                row_count = max(row_count, int(value.shape[0]))
            continue
        if isinstance(value, (list, tuple)):
            row_count = max(row_count, len(value))
    return max(1, row_count)


def _effective_batch_cpu_count(
    condition: dict[str, Any],
    configured_count: int = 0,
) -> int:
    requested_count = _normalized_batch_cpu_count(configured_count)
    if requested_count <= 0:
        requested_count = _batch_cpu_count()
    return max(1, min(requested_count, _condition_row_count(condition)))


def _run_equilibrium_batch_with_serial_retry(
    eq_module: Any,
    database: Any,
    condition: dict[str, Any],
    units: list[str],
    phases: list[str] | None,
    n_cpu: int,
    status: QPlainTextEdit | None,
    progress_callback=None,
) -> Any:
    try:
        return eq_module.equilib_batch(
            database,
            condition,
            unit=units,
            phases=phases,
            n_cpu=n_cpu,
            progress_callback=progress_callback,
        )
    except Exception:
        if n_cpu <= 1:
            raise
        if status is not None:
            status.setPlainText(
                "Parallel equilibrium calculation failed; retrying with one CPU."
            )
        return eq_module.equilib_batch(
            database,
            condition,
            unit=units,
            phases=phases,
            n_cpu=1,
            progress_callback=progress_callback,
        )


def _run_parallel_batch_with_serial_retry(
    run_batch,
    n_cpu: int,
    status: QPlainTextEdit | None,
    progress_callback,
    *,
    calculation_name: str,
) -> Any:
    """Retry a process-pool batch serially after a worker-startup failure."""
    try:
        return run_batch(n_cpu)
    except Exception as exc:
        if n_cpu <= 1 or not _is_parallel_worker_failure(exc):
            raise
        message = (
            f"Parallel {calculation_name} calculation failed; "
            "retrying with one CPU."
        )
        if status is not None:
            status.setPlainText(message)
        if progress_callback is not None:
            progress_callback(0, 0, message)
        return run_batch(1)


def _is_parallel_worker_failure(exc: Exception) -> bool:
    """Return whether an exception looks like a process-pool worker failure."""
    exc_type = type(exc).__name__.lower()
    message = str(exc).lower()
    return (
        "terminatedworkererror" in exc_type
        or "brokenprocesspool" in exc_type
        or "worker process managed by the executor" in message
        or "process died unexpectedly" in message
    )


def _molar_mass(
    stoichiometry: dict[str, float],
    atomic_masses: dict[str, float],
) -> float:
    total = 0.0
    for element, coefficient in stoichiometry.items():
        if element not in atomic_masses:
            raise ValueError(f"atomic mass is unavailable for {element}")
        total += coefficient * atomic_masses[element]
    return total


def _populate_phase_tree(
    tree: QTreeWidget,
    phases: list[str],
    selected_phases: list[str] | None = None,
    *,
    fit_to_contents: bool = True,
) -> None:
    tree.clear()
    tree.setColumnCount(3)
    _disable_inner_scrollbars(tree)
    solution_parent = QTreeWidgetItem(["Solutions", "", ""])
    compound_parent = QTreeWidgetItem(["Compounds", "", ""])
    solution_parent.setFirstColumnSpanned(True)
    compound_parent.setFirstColumnSpanned(True)
    solution_phases: list[str] = []
    compound_phases: list[str] = []
    for phase in phases:
        if _is_reference_ser_phase(phase):
            continue
        phase_list = compound_phases if _phase_is_compound(phase) else solution_phases
        phase_list.append(phase)
    selected = set(selected_phases) if selected_phases is not None else None
    _add_phase_columns(solution_parent, solution_phases, selected)
    _add_phase_columns(compound_parent, compound_phases, selected)
    tree.addTopLevelItem(solution_parent)
    tree.addTopLevelItem(compound_parent)
    tree.expandAll()
    if fit_to_contents:
        _fit_tree_to_items(tree)


def _populate_nucleation_undercooling_tree(
    tree: QTreeWidget,
    phases: list[str],
    undercooling: dict[str, float],
    liquid_phase_name: str,
    temperature_unit: str = "K",
) -> None:
    previous_blocked = tree.blockSignals(True)
    tree.clear()
    tree.setColumnCount(2)
    unit_label = temperature_unit.strip() or "K"
    tree.setHeaderLabels(["Phase", f"Undercooling [{unit_label}]"])
    _disable_inner_scrollbars(tree)

    solution_parent = QTreeWidgetItem(["Solutions", ""])
    compound_parent = QTreeWidgetItem(["Compounds", ""])
    solution_parent.setFirstColumnSpanned(True)
    compound_parent.setFirstColumnSpanned(True)

    for phase in phases:
        phase_name = str(phase).strip()
        if (
            not phase_name
            or phase_name == liquid_phase_name
            or _is_reference_ser_phase(phase_name)
        ):
            continue
        value = float(undercooling.get(phase_name, 0.0) or 0.0)
        item = QTreeWidgetItem([phase_name, _format_amount_value(value)])
        item.setFlags(item.flags() | Qt.ItemFlag.ItemIsEditable)
        parent = compound_parent if _phase_is_compound(phase_name) else solution_parent
        parent.addChild(item)

    tree.addTopLevelItem(solution_parent)
    tree.addTopLevelItem(compound_parent)
    tree.expandAll()
    _fit_tree_to_items(tree)
    tree.blockSignals(previous_blocked)


def _store_nucleation_undercooling_tree(
    module: CalculationModule,
    tree: QTreeWidget,
) -> None:
    try:
        module.nucleation_undercooling.update(
            _nucleation_undercooling_tree_values(tree, strict=False)
        )
    except ValueError:
        return


def _nucleation_undercooling_tree_values(
    tree: QTreeWidget,
    *,
    strict: bool,
) -> dict[str, float]:
    values: dict[str, float] = {}
    for parent_row in range(tree.topLevelItemCount()):
        parent = tree.topLevelItem(parent_row)
        for row in range(parent.childCount()):
            child = parent.child(row)
            phase = child.text(0).strip()
            if not phase:
                continue
            text = child.text(1).strip()
            if not text:
                value = 0.0
            else:
                try:
                    value = float(text)
                except ValueError as exc:
                    if strict:
                        raise ValueError(f"{phase} must be numeric.") from exc
                    continue
            if value < 0:
                if strict:
                    raise ValueError(f"{phase} must be non-negative.")
                continue
            values[phase] = value
    return values


def _nucleation_undercooling_for_phases(
    liquid_phase_name: str,
    phases: list[str] | None,
    undercooling: dict[str, float],
) -> dict[str, float]:
    if not phases:
        return {
            phase: float(value)
            for phase, value in undercooling.items()
            if phase != liquid_phase_name
        }
    return {
        phase: float(undercooling.get(phase, 0.0) or 0.0)
        for phase in phases
        if phase != liquid_phase_name
    }


def _missing_nucleation_undercooling_phases(
    liquid_phase_name: str,
    phases: list[str] | None,
    undercooling: dict[str, float],
) -> list[str]:
    """Return selected solid phases with no user-provided undercooling entry."""
    if not phases:
        return []
    return [
        phase
        for phase in phases
        if phase != liquid_phase_name and phase not in undercooling
    ]


def _phase_is_compound(phase: str) -> bool:
    return "(s" in phase or "_s" in phase


def _default_scheil_phase_selection(
    phases: list[str],
    database: dict[str, Any] | None = None,
) -> list[str]:
    """Return GUI default Scheil phases using order/disorder metadata."""
    visible_phases = [
        phase
        for phase in phases
        if not _is_reference_ser_phase(phase)
    ]
    if database is None:
        return visible_phases
    selected, _excluded = _metadata_scheil_phase_selection(database, visible_phases)
    return selected


def _scheil_default_excluded_ordered_phases(
    phases: list[str],
    database: dict[str, Any] | None,
) -> list[str]:
    """Return ordered phases hidden from the GUI's default Scheil selection."""
    if database is None:
        return []
    visible_phases = [
        phase
        for phase in phases
        if not _is_reference_ser_phase(phase)
    ]
    _selected, excluded = _metadata_scheil_phase_selection(database, visible_phases)
    return excluded


def _add_phase_columns(
    parent: QTreeWidgetItem,
    phases: list[str],
    selected_phases: set[str] | None,
) -> None:
    for row_start in range(0, len(phases), 3):
        row_phases = phases[row_start : row_start + 3]
        item = QTreeWidgetItem(row_phases + [""] * (3 - len(row_phases)))
        item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
        for column, phase in enumerate(row_phases):
            item.setText(column, phase)
            check_state = (
                Qt.CheckState.Checked
                if selected_phases is None or phase in selected_phases
                else Qt.CheckState.Unchecked
            )
            item.setCheckState(column, check_state)
        parent.addChild(item)


def _set_phase_tree_selection(
    tree: QTreeWidget,
    selected_phases: list[str],
) -> None:
    selected = set(selected_phases)
    for parent_row in range(tree.topLevelItemCount()):
        parent = tree.topLevelItem(parent_row)
        for row in range(parent.childCount()):
            child = parent.child(row)
            for column in range(tree.columnCount()):
                phase = child.text(column).strip()
                if not phase:
                    continue
                child.setCheckState(
                    column,
                    Qt.CheckState.Checked
                    if phase in selected
                    else Qt.CheckState.Unchecked,
                )


def _selected_phase_names(tree: QTreeWidget) -> list[str]:
    phases = []
    for parent_row in range(tree.topLevelItemCount()):
        parent = tree.topLevelItem(parent_row)
        for row in range(parent.childCount()):
            child = parent.child(row)
            for column in range(tree.columnCount()):
                phase = child.text(column).strip()
                if phase and child.checkState(column) == Qt.CheckState.Checked:
                    phases.append(phase)
    return phases


def _all_phase_names(tree: QTreeWidget) -> list[str]:
    phases = []
    for parent_row in range(tree.topLevelItemCount()):
        parent = tree.topLevelItem(parent_row)
        for row in range(parent.childCount()):
            child = parent.child(row)
            for column in range(tree.columnCount()):
                phase = child.text(column).strip()
                if phase:
                    phases.append(phase)
    return phases


def _resolve_liquid_phase_name(
    requested: str,
    phase_names: list[str],
) -> str:
    requested = requested.strip()
    for phase in phase_names:
        if phase == requested:
            return phase
    for phase in phase_names:
        if phase.upper() == requested.upper():
            return phase

    preferred = _preferred_liquid_phase(phase_names)
    if preferred and (not requested or requested.upper() == "LIQUID"):
        return preferred
    return requested or preferred


def _preferred_liquid_phase(phase_names: list[str]) -> str:
    for phase in phase_names:
        if phase.strip().upper() == "LIQUID":
            return phase
    for phase in phase_names:
        normalized = phase.strip().replace("_", "").replace("-", "").upper()
        if normalized.startswith("LIQ") or "LIQUID" in normalized:
            return phase
    return ""


def _toggle_phase_group(
    tree: QTreeWidget,
    group_name: str,
) -> None:
    parent = _phase_group_parent(tree, group_name)
    if parent is None:
        return
    checked = _phase_group_all_checked(parent, tree.columnCount())
    new_state = Qt.CheckState.Unchecked if checked else Qt.CheckState.Checked
    for row in range(parent.childCount()):
        child = parent.child(row)
        for column in range(tree.columnCount()):
            if child.text(column).strip():
                child.setCheckState(column, new_state)


def _phase_group_parent(
    tree: QTreeWidget,
    group_name: str,
) -> QTreeWidgetItem | None:
    for row in range(tree.topLevelItemCount()):
        parent = tree.topLevelItem(row)
        if parent.text(0) == group_name:
            return parent
    return None


def _phase_group_all_checked(
    parent: QTreeWidgetItem,
    column_count: int,
) -> bool:
    for row in range(parent.childCount()):
        child = parent.child(row)
        for column in range(column_count):
            if (
                child.text(column).strip()
                and child.checkState(column) != Qt.CheckState.Checked
            ):
                return False
    return True

__all__ = [name for name in globals() if name.startswith("_") and not name.startswith("__")]
