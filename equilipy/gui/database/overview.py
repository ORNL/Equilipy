"""DatabaseIR overview and thermodynamic record view helpers."""
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
from PySide6.QtGui import QBrush, QColor, QFontDatabase, QIcon
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
from equilipy.results import ResultTable
from equilipy.results.equilib import EquilibResult
from equilipy.results.scheil import ScheilResult
from equilipy.utils import G2HSCp, HSCp2G, NeumanKoppHSCp

from ..assets import calculation_icon_path, math_font_path
from ..calculation.state import CalculationDatabase, CalculationModule, CalculationSession
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
_MAGNETIC_PARAMETER_TYPES = {"BM", "BMAGN", "TC"}


from ..calculation.input_condition import _format_site_ratio
from ..calculation.session_serialization import _safe_file_stem
from ..calculation.results import _is_solution_phase
from .cp_editor import *

def _database_ir_element_names(database: DatabaseIR) -> list[str]:
    """Return sorted visible element symbols for a DatabaseIR overview."""
    return sorted(
        {
            element.symbol
            for element in database.elements
            if str(element.symbol).strip()
            and str(element.symbol).strip() not in {"/-", "VA"}
        }
    )


def _database_ir_has_content(database: DatabaseIR) -> bool:
    """Return whether a DatabaseIR contains displayable loaded content."""
    return any(
        (
            database.elements,
            database.species,
            database.functions,
            database.phases,
            database.parameters,
            database.diagnostics,
        )
    )


def _database_ir_phase_overview_rows(
    phases: list[Phase],
) -> list[tuple[Any, ...]]:
    """Return compact phase rows for the DatabaseIR overview."""
    rows: list[tuple[Any, ...]] = []
    for phase in phases:
        rows.append(
            (
                len(rows) + 1,
                phase.name,
                _phase_model_display_name(phase),
                _database_ir_phase_structure(phase),
                _format_species_list(_database_ir_phase_species(phase)),
            )
        )
    return rows


def _phase_model_display_name(phase: Phase) -> str:
    """Return the user-facing solution/compound model name."""
    model = str(phase.model or "").strip()
    model_upper = model.upper()
    if model_upper in {
        "CEF",
        "IDMX",
        "QKTO",
        "RKMP",
        "RKMPM",
        "SOLUTION",
        "SUBG",
        "SUBI",
        "SUBL",
        "SUBLM",
        "SUBM",
        "SUBQ",
    }:
        if len(phase.constituents) == 1:
            return "Bragg-Williams"
        return "CEF"
    return model


def _database_ir_phase_structure(phase: Phase) -> str:
    """Return a compact sublattice structure for a DatabaseIR phase."""
    parts: list[str] = []
    for constituent_set in phase.constituents:
        species = [name for name in constituent_set.species if str(name).strip()]
        if not species:
            continue
        ratio = _format_site_ratio(float(constituent_set.site_ratio))
        parts.append(f"({', '.join(species)}){ratio}")
    return "".join(parts)


def _database_ir_phase_species(phase: Phase) -> list[str]:
    """Return unique species listed across all phase sublattices."""
    species_names: list[str] = []
    for constituent_set in phase.constituents:
        for species_name in constituent_set.species:
            name = str(species_name).strip()
            if name and name not in species_names:
                species_names.append(name)
    return species_names


def _phase_formula(phase: Phase) -> str:
    """Return a compact stoichiometric sublattice formula."""
    parts: list[str] = []
    for constituent_set in phase.constituents:
        species = [
            _display_formula_species(species_name)
            for species_name in constituent_set.species
            if str(species_name).strip()
        ]
        if not species:
            continue
        species_text = ", ".join(species)
        ratio = _format_site_ratio(float(constituent_set.site_ratio))
        parts.append(f"({species_text}){ratio}")
    return "".join(parts)


def _display_formula_species(species_name: str) -> str:
    """Return a conventional display spelling for simple element symbols."""
    name = str(species_name).strip()
    if re.fullmatch(r"[A-Za-z]{1,2}", name):
        return name[0].upper() + name[1:].lower()
    return name


def _thermo_property_key(parameter_type: str) -> str:
    """Map a parameter type onto the GUI thermodynamic expression tabs."""
    normalized = str(parameter_type).strip().upper().replace(" ", "")
    if normalized in {"CP", "C_P", "HEATCAPACITY", "HEAT_CAPACITY"}:
        return "Cp"
    if normalized in {"H", "ENTHALPY"}:
        return "H"
    if normalized in {"S", "ENTROPY"}:
        return "S"
    return "G"


def _parameter_label(parameter: Parameter) -> str:
    """Return a compact parameter label."""
    target = ",".join(parameter.target)
    return f"{parameter.parameter_type}({target}; {parameter.order})"


def _function_cp_range_rows(function: FunctionDefinition) -> list[tuple[Any, ...]]:
    """Return parent-editor Cp range rows for one function."""
    ranges = list(function.temperature_ranges[:6])
    if not ranges and function.expression.strip():
        ranges = [("", "")]
    return [
        (f"Cp {_range_suffix(index, upper)}", lower, upper)
        for index, (lower, upper) in enumerate(ranges, start=1)
    ]


def _add_pertinent_expression_editor(
    section: QFrame,
    expression: str,
    object_name: str,
    validate_callback: Callable[[QLineEdit], None] | None,
) -> None:
    """Add the editable pertinent-function expression row."""
    header = QHBoxLayout()
    header.setContentsMargins(0, 0, 0, 0)
    header.setSpacing(8)
    pertinent_label = QLabel("Pertinent Functions")
    pertinent_label.setObjectName("DatabaseListLabel")
    header.addWidget(pertinent_label)
    validate_button = QPushButton("Validate")
    validate_button.setObjectName("SmallActionButton")
    header.addWidget(validate_button)
    header.addStretch(1)
    section.layout().addLayout(header)

    editor = QLineEdit(expression)
    editor.setObjectName(object_name)
    editor.setAlignment(Qt.AlignmentFlag.AlignLeft)
    editor.setPlaceholderText("-300+2*T+3*#Function")
    section.layout().addWidget(editor)
    if validate_callback is None:
        validate_button.setEnabled(False)
    else:
        validate_button.clicked.connect(lambda _checked=False: validate_callback(editor))


def _function_pertinent_expression(
    database: DatabaseIR,
    function: FunctionDefinition,
) -> str:
    """Return the editable expression that references pertinent Functions."""
    source_expression = str(function.tdb_src.source_expression or "").strip()
    if source_expression:
        return _expression_without_temperature_ranges(source_expression)
    if _expression_references_any_function(
        function.expression,
        database,
        exclude_name=function.name,
    ):
        return _expression_without_temperature_ranges(function.expression)
    return ""


def _compound_phase_pertinent_expression(
    database: DatabaseIR,
    phase: Phase,
) -> str:
    """Return the editable pertinent-function expression for one compound phase."""
    for parameter in _compound_phase_g_parameters(database, phase):
        expression = str(parameter.expression or "").strip()
        if _expression_references_any_function(expression, database):
            return _expression_without_temperature_ranges(expression)
        function_expression = _parameter_as_function_expression(parameter)
        if _expression_references_any_function(function_expression, database):
            return _expression_without_temperature_ranges(function_expression)
    return ""


def _expression_without_temperature_ranges(expression: str) -> str:
    """Return a single Gibbs expression without TDB piecewise range markers."""
    segments = _function_expression_segments(expression)
    if not segments:
        return str(expression or "").strip()
    segment_expressions = [
        segment_expression.strip()
        for _lower, _upper, segment_expression in segments
        if segment_expression.strip()
    ]
    if not segment_expressions:
        return ""
    first_expression = segment_expressions[0]
    if all(segment == first_expression for segment in segment_expressions):
        return first_expression
    return first_expression


def _expression_references_any_function(
    expression: str,
    database: DatabaseIR,
    *,
    exclude_name: str = "",
) -> bool:
    """Return whether an expression references a known Function name."""
    if not expression:
        return False
    excluded = exclude_name.upper()
    functions_by_name = {
        function.name.upper(): function
        for function in database.functions
        if function.name.upper() != excluded
    }
    return bool(_referenced_function_names(expression, functions_by_name))


def _function_overview_section(
    database: DatabaseIR,
    function: FunctionDefinition,
    validate_pertinent_expression: Callable[[QLineEdit], None] | None = None,
    *,
    title: str = "Function Overview",
) -> QFrame:
    """Return function Gibbs ranges plus editable pertinent-function expression."""
    _h298, _s298, _cp_rows, gibbs_rows = _function_thermo_data(function, database)
    name_warnings = _function_name_warning_messages(function)
    name_errors = _function_name_error_messages(function)
    continuity_warnings = _function_gibbs_continuity_warnings(gibbs_rows)
    action_icons = [
        icon
        for icon in (
            _function_name_warning_icon(name_warnings),
            _function_name_error_icon(name_errors),
            _function_continuity_warning_icon(continuity_warnings),
        )
        if icon is not None
    ]
    section = (
        _form_section_with_actions(title, action_icons)
        if action_icons
        else _form_section(title)
    )
    expression_label = QLabel("Gibbs Energy Expressions")
    expression_label.setObjectName("DatabaseListLabel")
    section.layout().addWidget(expression_label)

    expression_rows = _function_gibbs_expression_rows(function, database)
    expression_table = _function_expression_table(expression_rows)
    section.layout().addWidget(expression_table)

    _add_pertinent_expression_editor(
        section,
        _function_pertinent_expression(database, function),
        "FunctionPertinentExpressionInput",
        validate_pertinent_expression,
    )
    return section


def _function_gibbs_expression_rows(
    function: FunctionDefinition,
    database: DatabaseIR,
) -> list[tuple[Any, ...]]:
    """Return low/high temperature and Gibbs expression rows for a Function."""
    _h298, _s298, _cp_rows, gibbs_rows = _function_thermo_data(function, database)
    rows: list[tuple[Any, ...]] = []
    for row in gibbs_rows:
        rows.append(
            (
                _temperature_range_bound_label(float(row[0])),
                _temperature_range_bound_label(float(row[1])),
                _gibbs_row_expression(row),
            )
        )
    if rows:
        return rows
    return [("", "", function.expression)]


def _function_expression_table(rows: list[tuple[Any, ...]]) -> QTableWidget:
    """Return a compact read-only expression table for a Function overview."""
    headers = ("Low T [K]", "High T [K]", "Gibbs energy expression")
    table = QTableWidget(max(1, len(rows)), len(headers))
    table.setObjectName("FormTable")
    table.setHorizontalHeaderLabels(headers)
    table.verticalHeader().setVisible(False)
    table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
    table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
    table.setAlternatingRowColors(False)
    table.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
    for row_index, row in enumerate(rows or [("", "", "")]):
        for column_index, value in enumerate(row):
            item = QTableWidgetItem(str(value))
            item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            table.setItem(row_index, column_index, item)
    table.resizeColumnsToContents()
    _fit_table_columns_to_headers(table)
    table.horizontalHeader().setStretchLastSection(True)
    _disable_inner_scrollbars(table)
    _fit_table_to_rows(table, min_rows=1)
    return table


def _compound_cp_range_rows(parameters: list[Parameter]) -> list[tuple[Any, ...]]:
    """Return parent-editor Cp range rows for one compound phase."""
    cp_parameters = [
        parameter
        for parameter in parameters
        if _thermo_property_key(parameter.parameter_type) == "Cp"
    ]
    range_parameters = cp_parameters or [
        parameter
        for parameter in parameters
        if _thermo_property_key(parameter.parameter_type) == "G"
    ]
    rows: list[tuple[Any, ...]] = []
    for index, parameter in enumerate(range_parameters[:6], start=1):
        lower, upper = _parameter_temperature_range(parameter)
        property_name = _thermo_property_key(parameter.parameter_type)
        rows.append((f"{property_name} {_range_suffix(index, upper)}", lower, upper))
    return rows


def _compound_phase_display_label(phase: Phase, index: int) -> str:
    """Return the default compound phase label used in the sidebar."""
    label = str(getattr(phase, "metadata", {}).get("phase_label", "")).strip()
    if label:
        return label
    if _is_explicit_liquid_compound_phase(phase):
        return "L"
    return f"S{min(max(index, 1), 10)}"


def _compound_phase_header_widget(phase_label: str) -> QWidget:
    """Return an inline editable phase label for the editor subtitle line."""
    wrapper = QWidget()
    layout = QHBoxLayout(wrapper)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.setSpacing(4)
    label = QLabel("Phase:")
    label.setObjectName("EditorSubtitle")
    editor = QLineEdit(phase_label)
    editor.setObjectName("CompoundPhaseNameInlineInput")
    editor.setFrame(False)
    editor.setFixedWidth(70)
    editor.setCursorPosition(0)
    editor.setToolTip("Edit compound phase label")
    layout.addWidget(label)
    layout.addWidget(editor)
    layout.addStretch(1)
    return wrapper


def _is_explicit_liquid_compound_phase(phase: Phase) -> bool:
    """Return whether a compound phase was explicitly tagged as liquid."""
    command = str(phase.source.command).upper()
    if re.search(r"\bPHASE\s+\S+:L\b", command):
        return True
    return str(phase.name).upper() == "LIQUID"


def _compound_phase_overview_section(
    database: DatabaseIR,
    phase: Phase,
    phase_label: str,
) -> QFrame:
    """Return phase-state overview with reference values and expressions."""
    del phase_label
    h298, s298, _cp_rows, _gibbs_rows = _compound_phase_thermo_data(
        database,
        phase,
    )
    transition = _compound_phase_transition(database, phase)
    section = _form_section("Phase Overview")
    h_editor = _editable_line("" if h298 is None else _thermo_number(h298))
    h_editor.setObjectName("CompoundPhaseH298Input")
    s_editor = _editable_line("" if s298 is None else _thermo_number(s298))
    s_editor.setObjectName("CompoundPhaseS298Input")
    h_editor.setFixedWidth(_THERMO_COEFFICIENT_BOX_WIDTH)
    s_editor.setFixedWidth(_THERMO_COEFFICIENT_BOX_WIDTH)
    h_editor.editingFinished.connect(
        lambda: _commit_compound_phase_reference_state(
            database,
            phase,
            h_editor,
            s_editor,
        )
    )
    s_editor.editingFinished.connect(
        lambda: _commit_compound_phase_reference_state(
            database,
            phase,
            h_editor,
            s_editor,
        )
    )

    mode_row = QHBoxLayout()
    mode_row.setSpacing(8)
    heat_298_button = QPushButton("Heat at 298")
    heat_298_button.setObjectName("ThermoModeButton")
    heat_298_button.setCheckable(True)
    heat_transition_button = QPushButton("Heat at transition")
    heat_transition_button.setObjectName("ThermoModeButton")
    heat_transition_button.setCheckable(True)
    for button in (heat_298_button, heat_transition_button):
        button.setMinimumWidth(150)
    mode_row.addStretch(1)
    mode_row.addWidget(heat_298_button)
    mode_row.addWidget(heat_transition_button)
    mode_row.addStretch(1)
    section.layout().addLayout(mode_row)

    stack = QStackedWidget()

    reference_page = QWidget()
    reference_row = QHBoxLayout(reference_page)
    reference_row.setContentsMargins(0, 0, 0, 0)
    reference_row.setSpacing(36)
    reference_row.addWidget(
        _thermo_reference_editor_column(
            "<b>&Delta;H<sub>298.15 K</sub></b>",
            "[J/mol]",
            h_editor,
        )
    )
    reference_row.addWidget(
        _thermo_reference_editor_column(
            "<b>S<sub>298.15 K</sub></b>",
            "[J/(mol K)]",
            s_editor,
        )
    )
    reference_row.addStretch(1)
    stack.addWidget(reference_page)

    transition_page = QWidget()
    transition_row = QHBoxLayout(transition_page)
    transition_row.setContentsMargins(0, 0, 0, 0)
    transition_row.setSpacing(36)
    saved_transition_delta_h = phase.metadata.get("transition_delta_h")
    saved_transition_temperature = phase.metadata.get("transition_temperature")
    delta_h_text = (
        _thermo_number(transition["delta_h"])
        if transition is not None
        else str(saved_transition_delta_h or "")
    )
    transition_t_text = (
        _compact_number(transition["temperature"])
        if transition is not None
        else str(saved_transition_temperature or "")
    )
    transition_h_editor = _editable_line(delta_h_text)
    transition_h_editor.setObjectName("CompoundPhaseTransitionHInput")
    transition_t_editor = _editable_line(transition_t_text)
    transition_t_editor.setObjectName("CompoundPhaseTransitionTInput")
    transition_h_editor.setFixedWidth(_THERMO_COEFFICIENT_BOX_WIDTH)
    transition_t_editor.setFixedWidth(_THERMO_TEMPERATURE_BOX_WIDTH)
    if transition is None:
        transition_h_editor.setPlaceholderText("no transition")
        transition_t_editor.setPlaceholderText("no transition")
    transition_h_editor.editingFinished.connect(
        lambda: _commit_compound_phase_transition_state(
            database,
            phase,
            transition,
            transition_h_editor,
            transition_t_editor,
            h_editor,
            s_editor,
        )
    )
    transition_t_editor.editingFinished.connect(
        lambda: _commit_compound_phase_transition_state(
            database,
            phase,
            transition,
            transition_h_editor,
            transition_t_editor,
            h_editor,
            s_editor,
        )
    )
    transition_row.addWidget(
        _thermo_reference_editor_column(
            "<b>&Delta;H<sup>trans</sup></b>",
            "[J/mol]",
            transition_h_editor,
        )
    )
    transition_row.addWidget(
        _thermo_reference_editor_column(
            "<b>T<sup>trans</sup></b>",
            "[K]",
            transition_t_editor,
        )
    )
    transition_row.addStretch(1)
    stack.addWidget(transition_page)

    def set_heat_mode(mode: str) -> None:
        if mode == "298":
            updated_h298, updated_s298, _cp_rows, _gibbs_rows = (
                _compound_phase_thermo_data(database, phase)
            )
            if updated_h298 is not None:
                h_editor.setText(_thermo_number(updated_h298))
            if updated_s298 is not None:
                s_editor.setText(_thermo_number(updated_s298))
        phase.metadata["heat_mode"] = mode
        heat_298_button.setChecked(mode == "298")
        heat_transition_button.setChecked(mode == "transition")
        stack.setCurrentIndex(1 if mode == "transition" else 0)

    heat_298_button.clicked.connect(lambda _checked=False: set_heat_mode("298"))
    heat_transition_button.clicked.connect(
        lambda _checked=False: set_heat_mode("transition")
    )
    initial_mode = phase.metadata.get("heat_mode", "298")
    if initial_mode not in {"298", "transition"}:
        initial_mode = "298"
    set_heat_mode(initial_mode)
    section.layout().addWidget(stack)
    return section


def _compound_phase_parameters_section(
    database: DatabaseIR,
    phase: Phase,
    validate_pertinent_expression: Callable[[QLineEdit], None] | None = None,
) -> QFrame:
    """Return the compound phase parameter expression section."""
    section = _form_section("Parameters")
    expression_table = _function_expression_table(
        _compound_phase_gibbs_expression_rows(database, phase),
    )
    section.layout().addWidget(expression_table)

    _add_pertinent_expression_editor(
        section,
        _compound_phase_pertinent_expression(database, phase),
        "CompoundPhasePertinentExpressionInput",
        validate_pertinent_expression,
    )
    return section


def _phase_magnetic_parameter_rows(
    database: DatabaseIR,
    phase: Phase,
) -> list[tuple[Any, ...]]:
    """Return magnetic TC/BMAGN parameter rows attached to one phase."""
    phase_name = str(phase.name).strip().upper()
    phase_id = str(phase.metadata.get("phase_id", "")).strip()
    rows: list[tuple[Any, ...]] = []
    for parameter in database.parameters:
        if str(parameter.phase_name).strip().upper() != phase_name:
            continue
        if str(parameter.parameter_type).strip().upper() not in _MAGNETIC_PARAMETER_TYPES:
            continue
        parameter_phase_id = str(parameter.metadata.get("phase_id", "")).strip()
        if phase_id and parameter_phase_id and parameter_phase_id != phase_id:
            continue
        rows.append(
            (
                _canonical_magnetic_parameter_type(parameter.parameter_type),
                _magnetic_parameter_target_text(parameter),
                parameter.order,
                _parameter_as_function_expression(parameter),
            )
        )
    return rows


def _phase_magnetic_properties_section(
    database: DatabaseIR,
    phase: Phase,
) -> QFrame | None:
    """Return a magnetic-property table for a phase, when records exist."""
    rows = _phase_magnetic_parameter_rows(database, phase)
    if not rows:
        return None
    return _table_section(
        "Magnetic Properties",
        ("Parameter", "Target", "Order", "Expression"),
        rows,
    )


def _canonical_magnetic_parameter_type(parameter_type: str) -> str:
    """Return the GUI display name for a magnetic parameter type."""
    normalized = str(parameter_type).strip().upper()
    if normalized == "BM":
        return "BMAGN"
    return normalized


def _magnetic_parameter_target_text(parameter: Parameter) -> str:
    """Return a compact target label for magnetic parameters."""
    target = list(parameter.target)
    if len(target) == 1 and ":" in target[0]:
        return target[0]
    return ":".join(target) or "-"


def _compound_record_parameter_rows(
    database: DatabaseIR,
    record: CompoundRecord,
) -> list[tuple[Any, ...]]:
    """Return phase-scoped thermodynamic parameter rows for one compound."""
    rows: list[tuple[Any, ...]] = []
    for index, phase in enumerate(record.phases, start=1):
        phase_label = _compound_phase_display_label(phase, index)
        parameter_rows = _compound_phase_parameter_rows(database, phase)
        for row_index, row in enumerate(parameter_rows):
            parameter, low, high, expression = row
            rows.append(
                (
                    phase_label if row_index == 0 else "",
                    parameter,
                    low,
                    high,
                    expression,
                )
            )
    return rows or [("", "", "", "", "")]


def _compound_phase_parameter_rows(
    database: DatabaseIR,
    phase: Phase,
) -> list[tuple[Any, ...]]:
    """Return Gibbs-energy parameter rows attached to one compound phase."""
    rows: list[tuple[Any, ...]] = []
    for parameter in _compound_phase_g_parameters(database, phase):
        parameter_label = _parameter_label(parameter)
        expression = _parameter_as_function_expression(parameter)
        segments = _function_expression_segments(expression)
        if segments:
            for segment_index, (lower, upper, segment_expression) in enumerate(
                segments,
            ):
                rows.append(
                    (
                        parameter_label if segment_index == 0 else "",
                        _temperature_range_bound_label(lower),
                        _temperature_range_bound_label(upper),
                        segment_expression,
                    )
                )
            continue
        lower, upper = _parameter_temperature_range(parameter)
        rows.append(
            (
                parameter_label,
                _temperature_range_bound_label(float(lower))
                if lower != ""
                else "",
                _temperature_range_bound_label(float(upper))
                if upper != ""
                else "",
                parameter.expression,
            )
        )
    return rows or [("", "", "", "")]


def _compound_phase_thermo_data(
    database: DatabaseIR,
    phase: Phase,
) -> tuple[float | None, float | None, list[list[float]], list[list[float]]]:
    """Return best-effort H/S/Cp/G rows for one compound phase."""
    for parameter in _compound_phase_g_parameters(database, phase):
        expression = _parameter_as_function_expression(parameter)
        if not expression:
            continue
        function = FunctionDefinition(
            name=f"__{phase.name}_G__",
            expression=expression,
        )
        h298, s298, cp_rows, gibbs_rows = _function_thermo_data(
            function,
            database,
        )
        if h298 is not None or s298 is not None or cp_rows or gibbs_rows:
            return h298, s298, cp_rows, gibbs_rows
    return None, None, [], []


def _compound_phase_g_parameters(database: DatabaseIR, phase: Phase) -> list[Parameter]:
    """Return Gibbs-energy parameters attached to one compound phase."""
    phase_id = _compound_phase_state_id(phase)
    if phase_id:
        matched = [
            parameter
            for parameter in database.parameters
            if parameter.phase_name == phase.name
            and parameter.metadata.get("phase_id") == phase_id
            and _thermo_property_key(parameter.parameter_type) == "G"
        ]
        if matched:
            return matched
        # Existing databases may have gained a GUI phase id before their original
        # parameter was tagged. Use a single unclaimed same-name parameter.
        unclaimed = [
            parameter
            for parameter in database.parameters
            if parameter.phase_name == phase.name
            and not parameter.metadata.get("phase_id")
            and _thermo_property_key(parameter.parameter_type) == "G"
        ]
        same_name_phases = [
            candidate for candidate in database.phases if candidate.name == phase.name
        ]
        if len(unclaimed) == 1 and len(same_name_phases) == 1:
            unclaimed[0].metadata["phase_id"] = phase_id
            return unclaimed
        return []
    return [
        parameter
        for parameter in database.parameters
        if parameter.phase_name == phase.name
        and not parameter.metadata.get("phase_id")
        and _thermo_property_key(parameter.parameter_type) == "G"
    ]


def _compound_phase_state_id(phase: Phase) -> str:
    """Return the persistent GUI phase-state id for a compound phase."""
    return str(getattr(phase, "metadata", {}).get("phase_id", "")).strip()


def _new_compound_phase_state_id(database: DatabaseIR, compound_name: str) -> str:
    """Return a phase-state id unique within one DatabaseIR."""
    prefix = _safe_file_stem(compound_name) or "compound"
    existing = {
        str(getattr(phase, "metadata", {}).get("phase_id", "")).strip()
        for phase in database.phases
    }
    existing.update(
        str(getattr(parameter, "metadata", {}).get("phase_id", "")).strip()
        for parameter in database.parameters
    )
    for index in range(1, 1000):
        candidate = f"{prefix}_phase_{index}"
        if candidate not in existing:
            return candidate
    return f"{prefix}_phase_{len(existing) + 1}"


def _compound_phase_gibbs_expression_rows(
    database: DatabaseIR,
    phase: Phase,
) -> list[tuple[Any, ...]]:
    """Return Gibbs expression rows for a compound phase overview."""
    rows: list[tuple[Any, ...]] = []
    for parameter in _compound_phase_g_parameters(database, phase):
        expression = _parameter_as_function_expression(parameter)
        segments = _function_expression_segments(expression)
        if segments:
            rows.extend(
                (
                    _temperature_range_bound_label(lower),
                    _temperature_range_bound_label(upper),
                    segment_expression,
                )
                for lower, upper, segment_expression in segments
            )
            continue
        lower, upper = _parameter_temperature_range(parameter)
        rows.append(
            (
                _temperature_range_bound_label(float(lower))
                if lower != ""
                else "",
                _temperature_range_bound_label(float(upper))
                if upper != ""
                else "",
                parameter.expression,
            )
        )
    return rows or [("", "", "")]


def _parameter_as_function_expression(parameter: Parameter) -> str:
    """Return a parameter expression in FUNCTION-like piecewise form."""
    expression = str(parameter.expression).strip()
    if _function_expression_segments(expression):
        return expression

    command_expression = _parameter_command_expression(parameter.source.command)
    if _function_expression_segments(command_expression):
        return command_expression

    lower, upper = _parameter_temperature_range(parameter)
    lower_value = lower if lower != "" else 298.15
    upper_value = upper if upper != "" else 6000
    if expression:
        return f"{lower_value} {expression}; {upper_value} N"
    return ""


def _new_compound_g_parameter(phase: Phase) -> Parameter:
    """Return a default Gibbs parameter for a stoichiometric compound phase."""
    target = _compound_parameter_target_species(phase)
    target_text = ":".join(target)
    expression = "298.15 0; 6000 N"
    source = SourceRef(
        command=f"PARAMETER G({phase.name},{target_text};0) {expression} !",
    )
    return Parameter(
        phase_name=phase.name,
        parameter_type="G",
        target=target,
        order=0,
        expression=expression,
        source=source,
        metadata=(
            {"phase_id": _compound_phase_state_id(phase)}
            if _compound_phase_state_id(phase)
            else {}
        ),
    )


def _compound_parameter_target_species(phase: Phase) -> list[str]:
    """Return flattened constituent species for a compound parameter target."""
    target: list[str] = []
    for constituent_set in phase.constituents:
        for species_name in constituent_set.species:
            name = str(species_name).strip()
            if name and name not in target:
                target.append(name)
    return target or [phase.name]


def _set_parameter_gibbs_rows(
    parameter: Parameter,
    gibbs_rows: list[list[float]],
) -> None:
    """Store serialized Gibbs rows on a compound parameter."""
    expression = _gibbs_rows_to_tdb_expression(gibbs_rows)
    parameter.expression = expression
    parameter.metadata.pop("source_expression", None)
    command = str(parameter.source.command).strip()
    if not command:
        target = (
            ":".join(parameter.target)
            if parameter.target
            else parameter.phase_name
        )
        parameter.source.command = (
            f"PARAMETER G({parameter.phase_name},{target};0) {expression} !"
        )
        return
    match = _parameter_command_prefix_match(command)
    if match is None:
        parameter.source.command = f"{command} {expression} !"
        return
    parameter.source.command = f"{command[: match.end()]}{expression} !"


def _parameter_command_prefix_match(command: str) -> re.Match[str] | None:
    """Return the PARAMETER command prefix before the temperature expression."""
    return parameter_command_prefix_match(command)


def _parameter_command_expression(command: str) -> str:
    """Return expression text from a TDB PARAMETER command."""
    command_text = str(command or "").strip()
    match = _parameter_command_prefix_match(command_text)
    if match is None:
        return ""
    expression = command_text[match.end() :].strip()
    return re.sub(r"\s*!\s*$", "", expression).strip()


def _compound_parameter_for_range_payload(
    database: DatabaseIR,
    payload: ThermoRangePayload,
) -> Parameter | None:
    """Return the compound parameter backing a tree Cp payload."""
    if isinstance(payload.source, SourceRef):
        for parameter in database.parameters:
            if parameter.source is payload.source:
                return parameter
        source_command = str(payload.source.command).strip()
        if source_command:
            for parameter in database.parameters:
                if str(parameter.source.command).strip() == source_command:
                    return parameter

    compound_name = payload.owner_name.split(":", 1)[0]
    phase_label = (
        payload.owner_name.split(":", 1)[1]
        if ":" in payload.owner_name
        else ""
    )
    if phase_label:
        matching_phase = _compound_phase_by_label(
            database,
            compound_name,
            phase_label,
        )
        if matching_phase is not None:
            candidates = _compound_phase_g_parameters(database, matching_phase)
            for parameter in candidates:
                if _parameter_contains_payload_range(parameter, payload):
                    return parameter
            return candidates[0] if candidates else None

    candidates = [
        parameter
        for parameter in database.parameters
        if parameter.phase_name == compound_name
        and not parameter.metadata.get("phase_id")
        and _thermo_property_key(parameter.parameter_type) == "G"
    ]
    for parameter in candidates:
        if _parameter_contains_payload_range(parameter, payload):
            return parameter
    return candidates[0] if candidates else None


def _compound_phase_by_label(
    database: DatabaseIR,
    compound_name: str,
    phase_label: str,
) -> Phase | None:
    """Return a compound phase matching a sidebar phase label."""
    for index, phase in enumerate(
        [candidate for candidate in database.phases if candidate.name == compound_name],
        start=1,
    ):
        if _compound_phase_display_label(phase, index) == phase_label:
            return phase
    return None


def _parameter_contains_payload_range(
    parameter: Parameter,
    payload: ThermoRangePayload,
) -> bool:
    """Return whether a parameter expression contains a displayed Cp range."""
    try:
        payload_lower = float(payload.lower)
        payload_upper = float(payload.upper)
    except (TypeError, ValueError):
        return False
    for lower, upper, _segment in _function_expression_segments(
        _parameter_as_function_expression(parameter),
    ):
        if np.isclose(lower, payload_lower, rtol=1e-10, atol=1e-8) and np.isclose(
            upper,
            payload_upper,
            rtol=1e-10,
            atol=1e-8,
        ):
            return True
    return False


def _compound_pertinent_function_rows(
    database: DatabaseIR,
    phase: Phase,
) -> list[tuple[Any, ...]]:
    """Return referenced functions for a compound phase."""
    expressions = [
        _parameter_as_function_expression(parameter)
        for parameter in database.parameters
        if parameter.phase_name == phase.name
    ]
    rows: list[tuple[Any, ...]] = []
    for function in database.functions:
        factor = sum(
            _expression_multiplier_for_name(expression, function.name) or 0.0
            for expression in expressions
        )
        if factor == 0.0:
            continue
        rows.append(
            (
                function.name,
                _compact_number(factor),
                _temperature_ranges_label(function.temperature_ranges),
            )
        )
    return rows


def _range_suffix(index: int, upper: Any) -> str:
    """Return upper-bound range suffix when available, otherwise row index."""
    if upper != "":
        return _compact_number(upper)
    return str(index)


def _compact_number(value: Any) -> str:
    """Format a numeric value compactly for GUI labels."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if number.is_integer():
        return str(int(number))
    return f"{number:g}"


def _thermo_number(value: Any) -> str:
    """Format a thermodynamic coefficient with enough visible precision."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not np.isfinite(number):
        return str(value)
    if number == 0.0:
        return "0"
    return f"{number:.16g}"


def _parameter_temperature_range(
    parameter: Parameter,
) -> tuple[float | str, float | str]:
    """Best-effort extraction of TDB parameter lower/upper temperatures."""
    command = parameter.source.command
    lower_match = re.search(r"\)\s*([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)", command)
    upper_matches = re.findall(
        r";\s*([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)\s+[A-Z]\b",
        command,
        flags=re.IGNORECASE,
    )
    lower: float | str = ""
    upper: float | str = ""
    if lower_match:
        lower = float(lower_match.group(1))
    if upper_matches:
        upper = float(upper_matches[-1])
    return lower, upper


def _function_building_block_rows(
    database: DatabaseIR,
    function: FunctionDefinition,
) -> list[tuple[Any, ...]]:
    """Return function dependencies with multipliers and Cp ranges."""
    rows: list[tuple[Any, ...]] = []
    for referenced_function, factor in _function_dependency_terms(database, function):
        rows.append(
            (
                referenced_function.name,
                _compact_number(factor),
                _temperature_ranges_label(
                    referenced_function.temperature_ranges,
                ),
            )
        )
    return rows


def _temperature_ranges_label(ranges: list[tuple[float, float]]) -> str:
    """Format temperature ranges as compact GUI text."""
    return "; ".join(
        f"{_temperature_range_bound_label(lower)}-{_temperature_range_bound_label(upper)}"
        for lower, upper in ranges
    )


def _temperature_range_bound_label(value: float) -> str:
    """Format TDB range bounds without losing the 298.15 reference marker."""
    if abs(float(value) - 298.0) < 1e-9:
        return "298.15"
    return _compact_number(value)


def _expression_multiplier_for_name(expression: str, name: str) -> float | None:
    """Return the summed simple multiplier for a named function reference."""
    if not expression or not name:
        return None
    pattern = (
        rf"(?P<sign>[+\-]?)\s*"
        rf"(?:(?P<number>(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+\-]?\d+)?)\s*\*\s*)?"
        rf"(?<![A-Za-z0-9_]){re.escape(name)}(?![A-Za-z0-9_])"
    )
    factors: list[float] = []
    for match in re.finditer(pattern, expression, flags=re.IGNORECASE):
        sign = -1.0 if match.group("sign") == "-" else 1.0
        number_text = match.group("number")
        factor = float(number_text) if number_text is not None else 1.0
        factors.append(sign * factor)
    if not factors:
        return None
    return sum(factors)


def _compound_building_block_rows(
    database: DatabaseIR,
    parameters: list[Parameter],
) -> list[tuple[Any, ...]]:
    """Return parameters and referenced functions for one compound."""
    rows: list[tuple[Any, ...]] = []
    expressions: list[str] = []
    for parameter in parameters:
        rows.append(
            (
                _parameter_label(parameter),
                f"{_thermo_property_key(parameter.parameter_type)} parameter",
                parameter.expression,
            )
        )
        expressions.append(parameter.expression)

    joined_expression = "\n".join(expressions)
    for function in database.functions:
        if _expression_references_name(joined_expression, function.name):
            rows.append((function.name, "Referenced function", function.expression))
    return rows


def _expression_references_name(expression: str, name: str) -> bool:
    """Return whether an expression references a named thermodynamic function."""
    if not expression or not name:
        return False
    pattern = rf"(?<![A-Za-z0-9_]){re.escape(name)}(?![A-Za-z0-9_])"
    return re.search(pattern, expression, re.IGNORECASE) is not None


def _load_math_font_family() -> str | None:
    """Load the bundled math font and return its family name."""
    global _MATH_FONT_FAMILY
    if _MATH_FONT_FAMILY is not None:
        return _MATH_FONT_FAMILY
    font_path = math_font_path()
    if not font_path.exists():
        _MATH_FONT_FAMILY = ""
        return None
    font_id = QFontDatabase.addApplicationFont(str(font_path))
    if font_id < 0:
        _MATH_FONT_FAMILY = ""
        return None
    families = QFontDatabase.applicationFontFamilies(font_id)
    _MATH_FONT_FAMILY = families[0] if families else ""
    return _MATH_FONT_FAMILY or None


def _unique_messages(messages: Iterable[str]) -> list[str]:
    """Return messages with duplicates removed while preserving order."""
    seen: set[str] = set()
    unique: list[str] = []
    for message in messages:
        if message in seen:
            continue
        seen.add(message)
        unique.append(message)
    return unique


def _unique_diagnostics(diagnostics: Iterable[Diagnostic]) -> list[Diagnostic]:
    """Return diagnostics with duplicates removed while preserving order."""
    seen: set[tuple[str, str, str, int, str]] = set()
    unique: list[Diagnostic] = []
    for diagnostic in diagnostics:
        source = diagnostic.source
        key = (
            diagnostic.severity,
            diagnostic.message,
            source.file,
            source.line,
            source.command,
        )
        if key in seen:
            continue
        seen.add(key)
        unique.append(diagnostic)
    return unique


def _database_validation_object_count(database: DatabaseIR) -> int:
    """Return the number of objects shown in validation progress."""
    function_count = len(database.functions)
    compound_count = sum(
        1 for phase in database.phases if not _is_solution_phase(phase)
    )
    solution_count = sum(1 for phase in database.phases if _is_solution_phase(phase))
    return function_count + compound_count + solution_count


def _validation_progress_message(current: int, total: int, label: str) -> str:
    """Return a validation progress label with two-decimal percentage."""
    if total <= 0:
        percent = 100.0
        count_text = "0 / 0"
    else:
        bounded = min(max(current, 0), total)
        percent = 100.0 * bounded / total
        count_text = f"{bounded} / {total}"
    return f"{label}: {count_text} objects checked ({percent:.2f}%)"


def _iter_database_tree_payloads(model: DatabaseTreeModel) -> Iterator[Any]:
    """Yield every payload in a DatabaseTreeModel."""
    def visit(parent: QModelIndex) -> Iterator[Any]:
        for row in range(model.rowCount(parent)):
            index = model.index(row, 0, parent)
            yield model.payload(index)
            yield from visit(index)

    yield from visit(QModelIndex())


def _thermo_expression_tab_section(
    expressions_by_property: dict[str, list[tuple[str, str]]],
) -> QFrame:
    """Build the shared Cp/H/S/G expression section for functions and compounds."""
    section = _form_section("Thermodynamic Expressions")
    tabs = QTabWidget()
    tabs.setObjectName("ThermoExpressionTabs")
    for property_name in ("Cp", "H", "S", "G"):
        rows = expressions_by_property.get(property_name, [])
        table = QTableWidget(max(1, len(rows)), 2)
        table.setObjectName("FormTable")
        table.setHorizontalHeaderLabels(("Name", "Expression"))
        table.verticalHeader().setVisible(False)
        table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        table.setAlternatingRowColors(False)
        table.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        for row_index, row in enumerate(rows or [("", "")]):
            for column_index, value in enumerate(row):
                item = QTableWidgetItem(str(value))
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                table.setItem(row_index, column_index, item)
        table.resizeColumnsToContents()
        table.horizontalHeader().setStretchLastSection(True)
        _disable_inner_scrollbars(table)
        _fit_table_to_rows(table, min_rows=1)
        tabs.addTab(table, property_name)
    section.layout().addWidget(tabs)
    return section


def _object_label(value: Any) -> str:
    if isinstance(value, CompoundRecord):
        return value.name
    if isinstance(value, CompoundPhasePayload):
        return value.phase_label
    for attribute in ("name", "symbol", "phase_name", "severity"):
        if hasattr(value, attribute):
            return str(getattr(value, attribute))
    return type(value).__name__


def _object_subtitle(value: Any) -> str:
    if isinstance(value, CompoundRecord):
        return f"{len(value.phases)} phase(s)"
    if isinstance(value, CompoundPhasePayload):
        return value.phase.name
    if isinstance(value, Parameter):
        return f"{value.parameter_type}, order {value.order}"
    if isinstance(value, Phase):
        return _phase_model_display_name(value)
    if isinstance(value, FunctionDefinition):
        return value.expression
    if isinstance(value, Diagnostic):
        return value.message
    if hasattr(value, "__dataclass_fields__"):
        data = asdict(value)
        return ", ".join(f"{key}={item}" for key, item in list(data.items())[:2])
    return ""

__all__ = [name for name in globals() if name.startswith("_") and not name.startswith("__")]
