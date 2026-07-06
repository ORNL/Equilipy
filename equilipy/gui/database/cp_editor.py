"""Thermodynamic expression editor and conversion helpers."""
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


def _function_reference_state_section(
    database: DatabaseIR,
    function: FunctionDefinition,
) -> QFrame:
    """Return the compact Function reference-state editor section."""
    h298, s298, _cp_rows, gibbs_rows = _function_thermo_data(function, database)
    continuity_warnings = _function_gibbs_continuity_warnings(gibbs_rows)
    warning_icon = _function_continuity_warning_icon(continuity_warnings)
    section = (
        _form_section_with_actions("Properties", [warning_icon])
        if warning_icon is not None
        else _form_section("Properties")
    )
    grid = QGridLayout()
    grid.setHorizontalSpacing(14)
    grid.setVerticalSpacing(6)
    grid.setColumnStretch(0, 3)
    grid.setColumnStretch(1, 3)
    grid.setColumnStretch(2, 2)
    h_editor = _editable_line("" if h298 is None else _thermo_number(h298))
    s_editor = _editable_line("" if s298 is None else _thermo_number(s298))
    h_editor.editingFinished.connect(
        lambda: _commit_function_reference_state(
            database,
            function,
            h_editor,
            s_editor,
        )
    )
    s_editor.editingFinished.connect(
        lambda: _commit_function_reference_state(
            database,
            function,
            h_editor,
            s_editor,
        )
    )

    rows = [
        (
            _latex_label("&Delta;H<sub>298.15 K</sub> [J/mol]"),
            _latex_label("S<sub>298.15 K</sub> [J/(mol K)]"),
            QLabel("Refs"),
        ),
        (
            h_editor,
            s_editor,
            _read_only_line(_reference_tokens_from_source(function.source)),
        ),
    ]
    for row_index, row in enumerate(rows):
        for column_index, widget in enumerate(row):
            grid.addWidget(widget, row_index, column_index)
    section.layout().addLayout(grid)
    return section


def _function_continuity_warning_icon(warnings: list[str]) -> QLabel | None:
    """Return a warning icon for discontinuous piecewise Gibbs records."""
    if not warnings:
        return None
    label = QLabel()
    label.setObjectName("FunctionContinuityWarning")
    label.setPixmap(function_warning_icon().pixmap(QSize(22, 22)))
    label.setToolTip("Gibbs energy discontinuity:\n" + "\n".join(warnings))
    return label


def _function_name_warning_icon(warnings: list[str]) -> QLabel | None:
    """Return a warning icon for non-portable TDB function names."""
    if not warnings:
        return None
    label = QLabel()
    label.setObjectName("FunctionNameWarning")
    label.setPixmap(function_warning_icon().pixmap(QSize(22, 22)))
    label.setToolTip("\n".join(warnings))
    return label


def _function_name_warning_messages(function: FunctionDefinition) -> list[str]:
    """Return portability warnings for long TDB function names."""
    if len(function.name) <= _PORTABLE_TDB_FUNCTION_NAME_WIDTH:
        return []
    if len(function.name) > _MAX_TDB_FUNCTION_NAME_WIDTH:
        return []
    return [
        "Name is longer than "
        f"{_PORTABLE_TDB_FUNCTION_NAME_WIDTH} characters; some TDB tools may "
        "require a shorter Function name"
    ]


def _function_name_error_icon(errors: list[str]) -> QLabel | None:
    """Return an error icon for unsupported TDB function-name violations."""
    if not errors:
        return None
    label = QLabel()
    label.setObjectName("FunctionNameError")
    label.setPixmap(function_error_icon().pixmap(QSize(22, 22)))
    label.setToolTip("\n".join(errors))
    return label


def _function_name_error_messages(function: FunctionDefinition) -> list[str]:
    """Return hard function-name errors for GUI display."""
    if len(function.name) <= _MAX_TDB_FUNCTION_NAME_WIDTH:
        return []
    return [f"Name must be maximum {_MAX_TDB_FUNCTION_NAME_WIDTH} characters"]


def _function_gibbs_continuity_warnings(
    gibbs_rows: list[list[float]],
    *,
    tolerance: float = _GIBBS_CONTINUITY_TOLERANCE,
) -> list[str]:
    """Return warnings for adjacent Gibbs ranges with mismatched boundary G."""
    warnings: list[str] = []
    if len(gibbs_rows) <= 1:
        return warnings
    for previous, current in zip(gibbs_rows[:-1], gibbs_rows[1:], strict=True):
        transition_temperature = float(previous[1])
        if not np.isclose(
            transition_temperature,
            float(current[0]),
            rtol=1e-10,
            atol=1e-8,
        ):
            continue
        previous_gibbs = _gibbs_value_at_temperature(previous, transition_temperature)
        current_gibbs = _gibbs_value_at_temperature(current, transition_temperature)
        delta = abs(previous_gibbs - current_gibbs)
        if delta > tolerance:
            warnings.append(
                "Mismatching Gibbs energy at "
                f"{_compact_number(transition_temperature)} K: "
                f"Delta G = {_compact_number(delta)} J/mol"
            )
    return warnings


def _selected_gibbs_continuity_warnings(
    gibbs_rows: list[list[float]],
    range_index: int,
    *,
    tolerance: float = _GIBBS_CONTINUITY_TOLERANCE,
) -> list[str]:
    """Return warnings for the selected higher-temperature Gibbs range."""
    warnings: list[str] = []
    boundary = _selected_lower_boundary_gibbs_delta(
        gibbs_rows,
        range_index,
        tolerance=tolerance,
    )
    if boundary is not None:
        transition_temperature, delta = boundary
        warnings.append(
            "Mismatching Gibbs energy at "
            f"{_compact_number(transition_temperature)} K: "
            f"Delta G = {_compact_number(abs(delta))} J/mol"
        )
    return warnings


def _selected_lower_boundary_gibbs_delta(
    gibbs_rows: list[list[float]],
    range_index: int,
    *,
    tolerance: float = _GIBBS_CONTINUITY_TOLERANCE,
) -> tuple[float, float] | None:
    """Return transition T and correction delta for a selected higher G range."""
    if range_index <= 0 or range_index >= len(gibbs_rows):
        return None
    previous = gibbs_rows[range_index - 1]
    current = gibbs_rows[range_index]
    transition_temperature = float(current[0])
    if not np.isclose(
        float(previous[1]),
        transition_temperature,
        rtol=1e-10,
        atol=1e-8,
    ):
        return None
    previous_gibbs = _gibbs_value_at_temperature(previous, transition_temperature)
    current_gibbs = _gibbs_value_at_temperature(current, transition_temperature)
    delta = previous_gibbs - current_gibbs
    if abs(delta) <= tolerance:
        return None
    return transition_temperature, delta


def _selected_upper_boundary_gibbs_delta(
    gibbs_rows: list[list[float]],
    range_index: int,
    *,
    tolerance: float = _GIBBS_CONTINUITY_TOLERANCE,
) -> tuple[float, float] | None:
    """Return transition T and mismatch for a selected lower G range."""
    if range_index < 0 or range_index >= len(gibbs_rows) - 1:
        return None
    current = gibbs_rows[range_index]
    next_row = gibbs_rows[range_index + 1]
    transition_temperature = float(current[1])
    if not np.isclose(
        transition_temperature,
        float(next_row[0]),
        rtol=1e-10,
        atol=1e-8,
    ):
        return None
    current_gibbs = _gibbs_value_at_temperature(current, transition_temperature)
    next_gibbs = _gibbs_value_at_temperature(next_row, transition_temperature)
    delta = current_gibbs - next_gibbs
    if abs(delta) <= tolerance:
        return None
    return transition_temperature, delta


def _gibbs_coefficient_basis_value(
    row: list[float],
    coefficient_slot: int,
    temperature: float,
) -> float:
    """Return the G basis multiplier for one editable coefficient slot."""
    if temperature <= 0:
        return float("nan")
    coefficients = [float(value) for value in row[2:]]
    while len(coefficients) < 14:
        coefficients.append(0.0)
    if coefficient_slot == 0:
        return 1.0
    if coefficient_slot == 1:
        return float(temperature)
    if coefficient_slot == 2:
        return float(temperature * np.log(temperature))
    if coefficient_slot == 3:
        return float(temperature**2)
    if coefficient_slot == 4:
        return float(temperature**3)
    if coefficient_slot == 5:
        return float(1.0 / temperature)
    if coefficient_slot in {6, 8, 10, 12}:
        power = coefficients[coefficient_slot + 1]
        return float(temperature**power)
    return float("nan")


def _balanced_gibbs_rows_at_lower_boundary(
    gibbs_rows: list[list[float]],
    row_index: int,
    coefficient_slot: int,
    *,
    target_tolerance: float = 1e-10,
) -> list[list[float]] | None:
    """Return G rows with one higher-range coefficient balanced at its lower bound."""
    candidate_rows = [_complete_gibbs_row(row) for row in gibbs_rows]
    best_rows: list[list[float]] | None = None
    best_error = float("inf")
    for _iteration in range(8):
        boundary = _selected_lower_boundary_gibbs_delta(
            candidate_rows,
            row_index,
            tolerance=0.0,
        )
        if boundary is None:
            return candidate_rows
        transition_temperature, delta = boundary
        error = abs(delta)
        if error < best_error:
            best_error = error
            best_rows = [[float(value) for value in row] for row in candidate_rows]
        if error <= target_tolerance:
            return candidate_rows
        row = _complete_gibbs_row(candidate_rows[row_index])
        basis_value = _gibbs_coefficient_basis_value(
            row,
            coefficient_slot,
            transition_temperature,
        )
        if not np.isfinite(basis_value) or abs(basis_value) < 1e-14:
            return best_rows
        row[2 + coefficient_slot] += delta / basis_value
        candidate_rows[row_index] = row
        parsed_rows = _function_gibbs_rows(
            _gibbs_rows_to_tdb_expression(candidate_rows)
        )
        if not parsed_rows:
            return best_rows
        candidate_rows = [_complete_gibbs_row(row) for row in parsed_rows]
    return best_rows


def _balance_gibbs_rows_a_coefficients_sequentially(
    gibbs_rows: list[list[float]],
    *,
    tolerance: float = _GIBBS_CONTINUITY_TOLERANCE,
) -> tuple[list[list[float]], int]:
    """Balance discontinuous G ranges by shifting each higher range's ``a`` term."""
    candidate_rows, _selected_index = _normalized_gibbs_rows(
        [_complete_gibbs_row(row) for row in gibbs_rows]
    )
    if not candidate_rows:
        return [], 0
    balanced_count = 0
    for row_index in range(1, len(candidate_rows)):
        if (
            _selected_lower_boundary_gibbs_delta(
                candidate_rows,
                row_index,
                tolerance=tolerance,
            )
            is None
        ):
            continue
        balanced_rows = _balanced_gibbs_rows_at_lower_boundary(
            candidate_rows,
            row_index,
            coefficient_slot=0,
        )
        if balanced_rows is None:
            continue
        candidate_rows, _selected_index = _normalized_gibbs_rows(balanced_rows)
        balanced_count += 1
    return candidate_rows, balanced_count


def _complete_gibbs_row(row: list[float]) -> list[float]:
    """Return a Gibbs coefficient row padded to the full G2HSCp width."""
    values = [float(value) for value in row]
    if len(values) < 16:
        values.extend([0.0] * (16 - len(values)))
    return values[:16]


def _gibbs_value_at_temperature(row: list[float], temperature: float) -> float:
    """Evaluate one G2HSCp Gibbs coefficient row at temperature."""
    if temperature <= 0:
        return float("nan")
    coefficients = [float(value) for value in row[2:]]
    while len(coefficients) < 14:
        coefficients.append(0.0)
    value = (
        coefficients[0]
        + coefficients[1] * temperature
        + coefficients[2] * temperature * np.log(temperature)
        + coefficients[3] * temperature**2
        + coefficients[4] * temperature**3
        + coefficients[5] / temperature
    )
    for index in range(4):
        coefficient = coefficients[6 + 2 * index]
        power = coefficients[7 + 2 * index]
        if coefficient == 0.0:
            continue
        if _is_log_t_power(power):
            value += coefficient * np.log(temperature)
            continue
        value += coefficient * temperature**power
    return float(value)


def _enthalpy_value_at_temperature(row: list[float], temperature: float) -> float:
    """Evaluate H(T) for one G2HSCp Gibbs coefficient row."""
    if temperature <= 0:
        return float("nan")
    coefficients = [float(value) for value in row[2:]]
    while len(coefficients) < 14:
        coefficients.append(0.0)
    value = (
        coefficients[0]
        - coefficients[2] * temperature
        - coefficients[3] * temperature**2
        - 2.0 * coefficients[4] * temperature**3
        + 2.0 * coefficients[5] / temperature
    )
    for index in range(4):
        coefficient = coefficients[6 + 2 * index]
        power = coefficients[7 + 2 * index]
        if coefficient == 0.0:
            continue
        if _is_log_t_power(power):
            value += coefficient * (np.log(temperature) - 1.0)
            continue
        value += coefficient * (1.0 - power) * temperature**power
    return float(value)


def _entropy_value_at_temperature(row: list[float], temperature: float) -> float:
    """Evaluate S(T) for one G2HSCp Gibbs coefficient row."""
    if temperature <= 0:
        return float("nan")
    coefficients = [float(value) for value in row[2:]]
    while len(coefficients) < 14:
        coefficients.append(0.0)
    value = (
        -(coefficients[1] + coefficients[2])
        - coefficients[2] * np.log(temperature)
        - 2.0 * coefficients[3] * temperature
        - 3.0 * coefficients[4] * temperature**2
        + coefficients[5] / temperature**2
    )
    for index in range(4):
        coefficient = coefficients[6 + 2 * index]
        power = coefficients[7 + 2 * index]
        if coefficient == 0.0:
            continue
        if _is_log_t_power(power):
            value -= coefficient / temperature
            continue
        value -= coefficient * power * temperature ** (power - 1.0)
    return float(value)


def _commit_function_reference_state(
    database: DatabaseIR,
    function: FunctionDefinition,
    h_editor: QLineEdit,
    s_editor: QLineEdit,
) -> None:
    """Update the Function expression after editing H298/S298 in overview."""
    h298, s298, cp_rows, _gibbs_rows = _function_thermo_data(function, database)
    if not cp_rows:
        return
    updated_h298 = _editor_float(h_editor, h298 or 0.0)
    updated_s298 = _editor_float(s_editor, s298 or 0.0)
    try:
        gibbs_rows = HSCp2G(updated_h298, updated_s298, cp_rows).tolist()
    except (AssertionError, ValueError, TypeError):
        return
    _set_function_gibbs_rows(function, gibbs_rows)
    h_editor.setText(_compact_number(updated_h298))
    s_editor.setText(_compact_number(updated_s298))


def _commit_compound_phase_reference_state(
    database: DatabaseIR,
    phase: Phase,
    h_editor: QLineEdit,
    s_editor: QLineEdit,
) -> None:
    """Update one compound phase after editing H298/S298."""
    h298, s298, cp_rows, _gibbs_rows = _compound_phase_thermo_data(database, phase)
    if not cp_rows:
        return
    updated_h298 = _editor_float(h_editor, h298 or 0.0)
    updated_s298 = _editor_float(s_editor, s298 or 0.0)
    try:
        gibbs_rows = HSCp2G(updated_h298, updated_s298, cp_rows).tolist()
    except (AssertionError, ValueError, TypeError):
        return
    parameter = _first_compound_phase_g_parameter(database, phase)
    if parameter is None:
        return
    _set_parameter_gibbs_rows(parameter, gibbs_rows)
    phase.metadata["heat_mode"] = "298"
    h_editor.setText(_compact_number(updated_h298))
    s_editor.setText(_compact_number(updated_s298))


def _commit_compound_phase_transition_state(
    database: DatabaseIR,
    phase: Phase,
    transition: dict[str, Any] | None,
    delta_h_editor: QLineEdit,
    transition_temperature_editor: QLineEdit,
    h298_editor: QLineEdit | None = None,
    s298_editor: QLineEdit | None = None,
) -> None:
    """Update one compound phase from transition heat and temperature."""
    h298, s298, cp_rows, gibbs_rows = _compound_phase_thermo_data(database, phase)
    if h298 is None or s298 is None or not cp_rows or not gibbs_rows:
        return
    fallback_temperature = (
        float(transition["temperature"])
        if transition is not None
        else _safe_float(phase.metadata.get("transition_temperature"), 0.0)
    )
    fallback_delta_h = (
        float(transition["delta_h"])
        if transition is not None
        else _safe_float(phase.metadata.get("transition_delta_h"), 0.0)
    )
    transition_temperature = _editor_float(
        transition_temperature_editor,
        fallback_temperature,
    )
    if transition_temperature <= 0:
        return
    target_delta_h = _editor_float(
        delta_h_editor,
        fallback_delta_h,
    )
    phase.metadata["heat_mode"] = "transition"
    phase.metadata["transition_temperature"] = _compact_number(
        transition_temperature
    )
    phase.metadata["transition_delta_h"] = _compact_number(target_delta_h)

    from_phase = _compound_phase_transition_reference_phase(
        database,
        phase,
        transition,
        transition_temperature,
    )
    if from_phase is None:
        delta_h_editor.setText(_thermo_number(target_delta_h))
        transition_temperature_editor.setText(_compact_number(transition_temperature))
        return
    from_enthalpy = _phase_property_value_at_temperature(
        database,
        from_phase,
        transition_temperature,
        _enthalpy_value_at_temperature,
    )
    from_entropy = _phase_property_value_at_temperature(
        database,
        from_phase,
        transition_temperature,
        _entropy_value_at_temperature,
    )
    current_enthalpy = _rows_property_value_at_temperature(
        gibbs_rows,
        transition_temperature,
        _enthalpy_value_at_temperature,
    )
    current_entropy = _rows_property_value_at_temperature(
        gibbs_rows,
        transition_temperature,
        _entropy_value_at_temperature,
    )
    if any(
        value is None
        for value in (from_enthalpy, from_entropy, current_enthalpy, current_entropy)
    ):
        return
    target_enthalpy = float(from_enthalpy) + target_delta_h
    target_entropy = float(from_entropy) + target_delta_h / transition_temperature
    updated_h298 = float(h298) + target_enthalpy - float(current_enthalpy)
    updated_s298 = float(s298) + target_entropy - float(current_entropy)
    try:
        updated_gibbs_rows = HSCp2G(updated_h298, updated_s298, cp_rows).tolist()
    except (AssertionError, ValueError, TypeError):
        return
    parameter = _first_compound_phase_g_parameter(database, phase)
    if parameter is None:
        return
    _set_parameter_gibbs_rows(parameter, updated_gibbs_rows)
    delta_h_editor.setText(_thermo_number(target_delta_h))
    transition_temperature_editor.setText(_compact_number(transition_temperature))
    if h298_editor is not None:
        h298_editor.setText(_thermo_number(updated_h298))
    if s298_editor is not None:
        s298_editor.setText(_thermo_number(updated_s298))


def _compound_phase_transition_reference_phase(
    database: DatabaseIR,
    phase: Phase,
    transition: dict[str, Any] | None,
    transition_temperature: float,
) -> Phase | None:
    """Return the reference phase for a user-entered transition state."""
    if transition is not None:
        from_phase = transition.get("from_phase")
        if isinstance(from_phase, Phase):
            return from_phase

    records = [
        record
        for record in _compound_phase_stability_records(database, phase.name)
        if record["phase"] is not phase
    ]
    if not records:
        return None
    epsilon = max(1e-5, abs(transition_temperature) * 1e-8)
    reference_record = _stable_record_at_temperature(
        records,
        transition_temperature - epsilon,
    )
    if reference_record is None:
        reference_record = _stable_record_at_temperature(
            records,
            transition_temperature,
        )
    if reference_record is None:
        return None
    reference_phase = reference_record.get("phase")
    return reference_phase if isinstance(reference_phase, Phase) else None


def _first_compound_phase_g_parameter(
    database: DatabaseIR,
    phase: Phase,
) -> Parameter | None:
    """Return the primary Gibbs parameter for a compound phase."""
    parameters = _compound_phase_g_parameters(database, phase)
    return parameters[0] if parameters else None


def _compound_phase_transition(
    database: DatabaseIR,
    phase: Phase,
) -> dict[str, Any] | None:
    """Return the stable lower-envelope transition into a compound phase."""
    records = _compound_phase_stability_records(database, phase.name)
    if len(records) < 2:
        return None
    phase_record = next(
        (record for record in records if record["phase"] is phase),
        None,
    )
    if phase_record is None:
        return None

    candidates: set[float] = set()
    for left_index, left in enumerate(records):
        for right in records[left_index + 1 :]:
            candidates.update(_phase_pair_crossing_temperatures(left, right))
    for temperature in sorted(candidates):
        lower_record, higher_record = _stable_records_around_temperature(
            records,
            temperature,
        )
        if lower_record is None or higher_record is None:
            continue
        if lower_record["phase"] is higher_record["phase"]:
            continue
        if higher_record["phase"] is not phase:
            continue
        from_enthalpy = _rows_property_value_at_temperature(
            lower_record["rows"],
            temperature,
            _enthalpy_value_at_temperature,
        )
        to_enthalpy = _rows_property_value_at_temperature(
            higher_record["rows"],
            temperature,
            _enthalpy_value_at_temperature,
        )
        if from_enthalpy is None or to_enthalpy is None:
            continue
        return {
            "temperature": temperature,
            "delta_h": float(to_enthalpy) - float(from_enthalpy),
            "from_phase": lower_record["phase"],
            "from_label": lower_record["label"],
            "to_phase": higher_record["phase"],
            "to_label": higher_record["label"],
        }
    return None


def _compound_phase_stability_records(
    database: DatabaseIR,
    compound_name: str,
) -> list[dict[str, Any]]:
    """Return compound phase Gibbs rows for transition analysis."""
    phases = [
        phase
        for phase in database.phases
        if phase.name == compound_name and not _is_solution_phase(phase)
    ]
    records: list[dict[str, Any]] = []
    for index, phase in enumerate(phases, start=1):
        _h298, _s298, _cp_rows, gibbs_rows = _compound_phase_thermo_data(
            database,
            phase,
        )
        if not gibbs_rows:
            continue
        records.append(
            {
                "phase": phase,
                "label": _compound_phase_display_label(phase, index),
                "rows": [[float(value) for value in row] for row in gibbs_rows],
            }
        )
    return records


def _phase_pair_crossing_temperatures(
    left: dict[str, Any],
    right: dict[str, Any],
) -> set[float]:
    """Return pairwise Gibbs crossing temperatures for two phase records."""
    boundaries = sorted(
        {
            point
            for record in (left, right)
            for row in record["rows"]
            for point in (float(row[0]), float(row[1]))
        }
    )
    roots: set[float] = set()
    for lower, upper in zip(boundaries[:-1], boundaries[1:], strict=True):
        if upper <= lower:
            continue
        midpoint = 0.5 * (lower + upper)
        if _rows_gibbs_value_at_temperature(left["rows"], midpoint) is None:
            continue
        if _rows_gibbs_value_at_temperature(right["rows"], midpoint) is None:
            continue
        root = _bisect_phase_pair_crossing(left, right, lower, upper)
        if root is not None:
            roots.add(root)
    return roots


def _bisect_phase_pair_crossing(
    left: dict[str, Any],
    right: dict[str, Any],
    lower: float,
    upper: float,
) -> float | None:
    """Return one Gibbs crossing in a bracket, if present."""
    span = upper - lower
    edge = max(1e-8, abs(span) * 1e-9)
    lo = lower + edge
    hi = upper - edge
    if hi <= lo:
        return None

    def diff(temperature: float) -> float | None:
        left_value = _rows_gibbs_value_at_temperature(left["rows"], temperature)
        right_value = _rows_gibbs_value_at_temperature(right["rows"], temperature)
        if left_value is None or right_value is None:
            return None
        return float(left_value) - float(right_value)

    f_lo = diff(lo)
    f_hi = diff(hi)
    if f_lo is None or f_hi is None:
        return None
    if abs(f_lo) <= _GIBBS_CONTINUITY_TOLERANCE:
        return lo
    if abs(f_hi) <= _GIBBS_CONTINUITY_TOLERANCE:
        return hi
    if f_lo * f_hi > 0:
        return None
    for _iteration in range(80):
        mid = 0.5 * (lo + hi)
        f_mid = diff(mid)
        if f_mid is None:
            return None
        if abs(f_mid) <= _GIBBS_CONTINUITY_TOLERANCE:
            return mid
        if f_lo * f_mid <= 0:
            hi = mid
            f_hi = f_mid
        else:
            lo = mid
            f_lo = f_mid
    return 0.5 * (lo + hi)


def _stable_records_around_temperature(
    records: list[dict[str, Any]],
    temperature: float,
) -> tuple[dict[str, Any] | None, dict[str, Any] | None]:
    """Return stable phase records just below and above a temperature."""
    epsilon = max(1e-5, abs(temperature) * 1e-8)
    lower = _stable_record_at_temperature(records, temperature - epsilon)
    higher = _stable_record_at_temperature(records, temperature + epsilon)
    return lower, higher


def _stable_record_at_temperature(
    records: list[dict[str, Any]],
    temperature: float,
) -> dict[str, Any] | None:
    """Return the lowest-Gibbs record at a temperature."""
    candidates: list[tuple[float, dict[str, Any]]] = []
    for record in records:
        value = _rows_gibbs_value_at_temperature(record["rows"], temperature)
        if value is not None and np.isfinite(value):
            candidates.append((float(value), record))
    if not candidates:
        return None
    return min(candidates, key=lambda item: item[0])[1]


def _phase_property_value_at_temperature(
    database: DatabaseIR,
    phase: Phase,
    temperature: float,
    evaluator: Callable[[list[float], float], float],
) -> float | None:
    """Evaluate a phase property at temperature."""
    _h298, _s298, _cp_rows, gibbs_rows = _compound_phase_thermo_data(
        database,
        phase,
    )
    return _rows_property_value_at_temperature(gibbs_rows, temperature, evaluator)


def _rows_gibbs_value_at_temperature(
    rows: list[list[float]],
    temperature: float,
) -> float | None:
    """Evaluate Gibbs rows at temperature."""
    return _rows_property_value_at_temperature(
        rows,
        temperature,
        _gibbs_value_at_temperature,
    )


def _rows_property_value_at_temperature(
    rows: list[list[float]],
    temperature: float,
    evaluator: Callable[[list[float], float], float],
) -> float | None:
    """Evaluate a property over piecewise Gibbs rows."""
    row = _row_for_temperature(rows, temperature)
    if row is None:
        return None
    value = evaluator(row, temperature)
    return value if np.isfinite(value) else None


def _row_for_temperature(
    rows: list[list[float]],
    temperature: float,
) -> list[float] | None:
    """Return the piecewise row covering a temperature."""
    for row in rows:
        lower = float(row[0])
        upper = float(row[1])
        if lower - 1e-8 <= temperature <= upper + 1e-8:
            return row
    return None


class _FunctionThermoEditor(QWidget):
    """Interactive Cp/G/H/S editor for one Function record."""

    _modes = ("G", "Cp", "H", "S")

    def __init__(
        self,
        function: FunctionDefinition,
        *,
        database: DatabaseIR | None = None,
        range_index: int = 0,
        range_changed_callback: Any | None = None,
        continuity_changed_callback: Any | None = None,
        gibbs_changed_callback: Any | None = None,
    ) -> None:
        super().__init__()
        self.function = function
        self.database = database
        self.range_index = max(0, range_index)
        self.range_changed_callback = range_changed_callback
        self.continuity_changed_callback = continuity_changed_callback
        self.gibbs_changed_callback = gibbs_changed_callback
        self.current_mode = "G"
        self.h298, self.s298, self.cp_rows, self.gibbs_rows = (
            _function_thermo_data(function, database)
        )
        self._editor_refs: dict[str, Any] = {}
        self._buttons: dict[str, QPushButton] = {}
        self._dirty_modes: set[str] = set()
        self.range_widget = self._build_range_widget()

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(12)

        selector = QFrame()
        selector.setObjectName("ModeSelector")
        selector_layout = QHBoxLayout(selector)
        selector_layout.setContentsMargins(0, 0, 0, 0)
        selector_layout.setSpacing(8)
        selector_layout.addStretch(1)
        for mode in self._modes:
            button = QPushButton(mode)
            button.setObjectName("ThermoModeButton")
            button.setCheckable(True)
            button.setMinimumWidth(64)
            button.clicked.connect(
                lambda _checked=False, selected=mode: self._select_mode(selected)
            )
            self._buttons[mode] = button
            selector_layout.addWidget(button)
        selector_layout.addStretch(1)

        self.stack = QStackedWidget()
        self.stack.setObjectName("FunctionThermoStack")
        layout.addWidget(selector)
        layout.addWidget(self.stack)
        self._refresh_pages()
        self._update_continuity_warning(self.gibbs_rows)

    def _build_range_widget(self) -> QWidget:
        """Return the header-owned temperature range editor."""
        widget = QWidget()
        layout = QHBoxLayout(widget)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)
        lower, upper = self._current_temperature_range()
        self.range_lower_editor = _editable_line(_compact_number(lower))
        self.range_upper_editor = _editable_line(_compact_number(upper))
        self.range_lower_editor.setObjectName("ThermoRangeLowerInput")
        self.range_upper_editor.setObjectName("ThermoRangeUpperInput")
        self.range_lower_editor.setFixedWidth(_THERMO_TEMPERATURE_BOX_WIDTH)
        self.range_upper_editor.setFixedWidth(_THERMO_TEMPERATURE_BOX_WIDTH)
        self.range_lower_editor.editingFinished.connect(self._commit_temperature_range)
        self.range_upper_editor.editingFinished.connect(self._commit_temperature_range)
        self.continuity_warning_label = QLabel()
        self.continuity_warning_label.setObjectName("FunctionContinuityWarning")
        self.continuity_warning_label.setPixmap(
            function_warning_icon().pixmap(QSize(22, 22))
        )
        self.continuity_warning_label.hide()
        layout.addWidget(QLabel("Range from"))
        layout.addWidget(self.range_lower_editor)
        layout.addWidget(QLabel("[K] to"))
        layout.addWidget(self.range_upper_editor)
        layout.addWidget(QLabel("[K]"))
        layout.addWidget(self.continuity_warning_label)
        layout.addStretch(1)
        return widget

    def _current_temperature_range(self) -> tuple[float, float]:
        rows = self.cp_rows or self.gibbs_rows
        if rows:
            index = min(self.range_index, len(rows) - 1)
            return float(rows[index][0]), float(rows[index][1])
        return 298.15, 6000.0

    def _commit_temperature_range(self) -> None:
        current_lower, current_upper = self._current_temperature_range()
        lower = _editor_float(self.range_lower_editor, current_lower)
        upper = _editor_float(self.range_upper_editor, current_upper)
        index = min(self.range_index, max(0, len(self.gibbs_rows) - 1))
        for rows in (self.cp_rows, self.gibbs_rows):
            if rows:
                row_index = min(self.range_index, len(rows) - 1)
                rows[row_index][0] = lower
                rows[row_index][1] = upper
        if self.gibbs_rows:
            self._set_gibbs_rows(
                self.gibbs_rows,
                selected_index=index,
            )
        lower, upper = self._current_temperature_range()
        self.range_lower_editor.setText(_compact_number(lower))
        self.range_upper_editor.setText(_compact_number(upper))
        if self.range_changed_callback is not None:
            self.range_changed_callback(self.function, self.range_index)

    def _select_mode(self, mode: str) -> None:
        if mode == self.current_mode:
            return
        self._capture_current_mode()
        self.current_mode = mode
        self._refresh_pages()

    def _refresh_pages(self) -> None:
        while self.stack.count():
            widget = self.stack.widget(0)
            self.stack.removeWidget(widget)
            widget.deleteLater()
        self._editor_refs.clear()
        selected_cp_rows = _selected_thermo_rows(self.cp_rows, self.range_index)
        selected_gibbs_rows = _selected_thermo_rows(self.gibbs_rows, self.range_index)

        pages = {
            "G": _function_coefficient_tab(
                selected_gibbs_rows,
                _gibbs_term_rows,
                _gibbs_expression_text(),
                editor_refs=self._editor_refs.setdefault("G", {}),
                show_temperature_rows=False,
                balance_callbacks=self._gibbs_balance_callbacks(),
            ),
            "Cp": _function_cp_tab(
                self.h298,
                self.s298,
                selected_cp_rows,
                editor_refs=self._editor_refs.setdefault("Cp", {}),
            ),
            "H": _function_coefficient_tab(
                selected_gibbs_rows,
                _enthalpy_term_rows,
                _enthalpy_expression_text(),
                editor_refs=self._editor_refs.setdefault("H", {}),
                show_temperature_rows=False,
            ),
            "S": _function_coefficient_tab(
                selected_gibbs_rows,
                _entropy_term_rows,
                _entropy_expression_text(),
                editor_refs=self._editor_refs.setdefault("S", {}),
                show_temperature_rows=False,
            ),
        }
        for mode in self._modes:
            page = pages[mode]
            for editor in page.findChildren(QLineEdit):
                editor.textEdited.connect(
                    lambda _text="", edited_mode=mode: (
                        self._mark_mode_dirty_and_update_warning(edited_mode)
                    )
                )
            self.stack.addWidget(page)
            self._buttons[mode].setChecked(mode == self.current_mode)
        self.stack.setCurrentIndex(self._modes.index(self.current_mode))
        self._update_continuity_warning(self.gibbs_rows)

    def _gibbs_balance_callbacks(self) -> dict[int, Any] | None:
        """Return Bal. callbacks for G coefficient rows."""
        return {
            term_index: (
                lambda _checked=False, slot=coefficient_slot: (
                    self._balance_gibbs_coefficient(slot)
                )
            )
            for term_index, coefficient_slot in enumerate(_GIBBS_COEFFICIENT_SLOTS)
        }

    def _balance_gibbs_coefficient(self, coefficient_slot: int) -> None:
        """Adjust one higher-range G coefficient to remove boundary mismatch."""
        candidate_rows = self._candidate_gibbs_rows_from_current_editor()
        row_index = min(self.range_index, len(candidate_rows) - 1)
        optimized_rows = self._optimized_gibbs_balance_rows(
            candidate_rows,
            row_index,
            coefficient_slot,
        )
        if optimized_rows is None:
            return
        self._set_gibbs_rows(optimized_rows, selected_index=self.range_index)
        self._dirty_modes.discard("G")
        self.current_mode = "G"
        self._refresh_pages()

    def _optimized_gibbs_balance_rows(
        self,
        gibbs_rows: list[list[float]],
        row_index: int,
        coefficient_slot: int,
    ) -> list[list[float]] | None:
        """Balance a serialized G expression to the warning tolerance."""
        return _balanced_gibbs_rows_at_lower_boundary(
            gibbs_rows,
            row_index,
            coefficient_slot,
        )

    def _mark_mode_dirty_and_update_warning(self, mode: str) -> None:
        """Mark an editor tab dirty and refresh the continuity warning preview."""
        self._dirty_modes.add(mode)
        self._update_continuity_warning(
            self._candidate_gibbs_rows_from_current_editor()
        )

    def _candidate_gibbs_rows_from_current_editor(self) -> list[list[float]]:
        """Return prospective Gibbs rows from the active tab without committing."""
        if self.current_mode == "Cp":
            captured = _capture_cp_editor(
                self._editor_refs.get("Cp", {}),
                _selected_thermo_rows(self.cp_rows, self.range_index),
            )
            if captured is None:
                return self.gibbs_rows
            h298, s298, cp_rows = captured
            try:
                full_cp_rows = [list(row) for row in self.cp_rows]
                if full_cp_rows and cp_rows:
                    full_cp_rows[self.range_index] = cp_rows[0]
                else:
                    full_cp_rows = cp_rows
                return HSCp2G(h298, s298, full_cp_rows).tolist()
            except (AssertionError, ValueError, TypeError):
                return self.gibbs_rows

        captured_rows = _capture_expression_rows(
            self._editor_refs.get(self.current_mode, {}),
            self.current_mode,
            _selected_thermo_rows(self.gibbs_rows, self.range_index),
        )
        if captured_rows is None:
            return self.gibbs_rows
        if self.current_mode == "G":
            return self._replace_selected_gibbs_row(captured_rows[0])
        if self.current_mode == "H":
            selected_row = _gibbs_rows_from_enthalpy_rows(
                captured_rows,
                _selected_thermo_rows(self.gibbs_rows, self.range_index),
            )[0]
            return self._replace_selected_gibbs_row(selected_row)
        if self.current_mode == "S":
            selected_row = _gibbs_rows_from_entropy_rows(
                captured_rows,
                _selected_thermo_rows(self.gibbs_rows, self.range_index),
            )[0]
            return self._replace_selected_gibbs_row(selected_row)
        return self.gibbs_rows

    def _capture_current_mode(self) -> None:
        if self.current_mode not in self._dirty_modes:
            return
        if not self.gibbs_rows and self.current_mode != "Cp":
            return

        if self.current_mode == "Cp":
            captured = _capture_cp_editor(
                self._editor_refs.get("Cp", {}),
                _selected_thermo_rows(self.cp_rows, self.range_index),
            )
            if captured is None:
                return
            h298, s298, cp_rows = captured
            try:
                full_cp_rows = [list(row) for row in self.cp_rows]
                if full_cp_rows and cp_rows:
                    full_cp_rows[self.range_index] = cp_rows[0]
                else:
                    full_cp_rows = cp_rows
                gibbs_rows = HSCp2G(h298, s298, full_cp_rows).tolist()
                self._set_gibbs_rows(gibbs_rows, selected_index=self.range_index)
                self._dirty_modes.discard("Cp")
            except (AssertionError, ValueError, TypeError):
                return
            return

        captured_rows = _capture_expression_rows(
            self._editor_refs.get(self.current_mode, {}),
            self.current_mode,
            _selected_thermo_rows(self.gibbs_rows, self.range_index),
        )
        if captured_rows is None:
            return
        if self.current_mode == "G":
            gibbs_rows = self._replace_selected_gibbs_row(captured_rows[0])
            self._set_gibbs_rows(gibbs_rows, selected_index=self.range_index)
        elif self.current_mode == "H":
            selected_row = _gibbs_rows_from_enthalpy_rows(
                captured_rows,
                _selected_thermo_rows(self.gibbs_rows, self.range_index),
            )[0]
            self._set_gibbs_rows(
                self._replace_selected_gibbs_row(selected_row),
                selected_index=self.range_index,
            )
        elif self.current_mode == "S":
            selected_row = _gibbs_rows_from_entropy_rows(
                captured_rows,
                _selected_thermo_rows(self.gibbs_rows, self.range_index),
            )[0]
            self._set_gibbs_rows(
                self._replace_selected_gibbs_row(selected_row),
                selected_index=self.range_index,
            )
        self._dirty_modes.discard(self.current_mode)

    def _replace_selected_gibbs_row(self, row: list[float]) -> list[list[float]]:
        """Return full Gibbs rows with the selected segment replaced."""
        rows = [list(value) for value in self.gibbs_rows]
        if rows:
            rows[self.range_index] = list(row)
            return rows
        return [list(row)]

    def _set_gibbs_rows(
        self,
        gibbs_rows: list[list[float]],
        *,
        selected_index: int | None = None,
    ) -> None:
        selected_row = self.range_index if selected_index is None else selected_index
        gibbs_rows, normalized_index = _normalized_gibbs_rows(
            gibbs_rows,
            selected_index=selected_row,
        )
        try:
            h298, s298, cp = G2HSCp(gibbs_rows)
        except (AssertionError, ValueError, TypeError):
            return
        self.gibbs_rows = [[float(value) for value in row] for row in gibbs_rows]
        self.range_index = normalized_index
        self.h298 = float(h298)
        self.s298 = float(s298)
        self.cp_rows = cp.tolist()
        _set_function_gibbs_rows(self.function, self.gibbs_rows)
        if self.gibbs_changed_callback is not None:
            self.gibbs_changed_callback(self.function, self.gibbs_rows)
        self._update_continuity_warning(self.gibbs_rows)

    def _update_continuity_warning(self, gibbs_rows: list[list[float]]) -> None:
        """Update the active Function editor continuity warning indicator."""
        warnings = _function_gibbs_continuity_warnings(gibbs_rows)
        selected_warnings = _selected_gibbs_continuity_warnings(
            gibbs_rows,
            self.range_index,
        )
        self.continuity_warning_label.setVisible(bool(selected_warnings))
        self.continuity_warning_label.setToolTip(
            "Gibbs energy discontinuity:\n" + "\n".join(selected_warnings)
            if selected_warnings
            else ""
        )
        if self.continuity_changed_callback is not None:
            self.continuity_changed_callback(
                self.function,
                warnings,
                selected_warnings,
                gibbs_rows,
            )
        self._update_balance_buttons(gibbs_rows)

    def _update_balance_buttons(self, gibbs_rows: list[list[float]]) -> None:
        """Show Bal. buttons only for discontinuous higher-temperature G rows."""
        higher_side = _selected_lower_boundary_gibbs_delta(
            gibbs_rows,
            self.range_index,
        )
        lower_side = _selected_upper_boundary_gibbs_delta(
            gibbs_rows,
            self.range_index,
        )
        visible = self.current_mode == "G" and higher_side is not None
        enabled = self.current_mode == "G" and higher_side is not None
        if enabled:
            tooltip = "Balance Gibbs continuity at transition to 1e-10 J/mol"
        elif lower_side is not None:
            tooltip = (
                "Select the higher-temperature Cp range to balance this transition"
            )
        else:
            tooltip = "No active lower-boundary discontinuity to balance"
        for row_ref in self._editor_refs.get("G", {}).get("rows", []):
            for term_ref in row_ref.get("terms", []):
                button = term_ref.get("balance_button")
                if isinstance(button, QToolButton):
                    button.setVisible(visible)
                    button.setEnabled(enabled)
                    button.setToolTip(tooltip)


def _function_thermo_expression_section(editor: _FunctionThermoEditor) -> QFrame:
    """Return the interactive H/S/Cp/G function expression editor."""
    section = _untitled_form_section()
    section.layout().addWidget(editor)
    return section


def _selected_thermo_rows(
    rows: list[list[float]],
    range_index: int,
) -> list[list[float]]:
    """Return a one-row view of the selected temperature segment."""
    if not rows:
        return []
    index = min(max(0, range_index), len(rows) - 1)
    return [list(rows[index])]


def _function_cp_tab(
    h298: float | None,
    s298: float | None,
    cp_rows: list[list[float]],
    *,
    editor_refs: dict[str, Any] | None = None,
) -> QWidget:
    """Build the Cp tab using fixed and optional powers."""
    tab = QWidget()
    layout = QVBoxLayout(tab)
    layout.setContentsMargins(10, 10, 10, 10)
    layout.setSpacing(12)

    if not cp_rows:
        layout.addWidget(
            QLabel(
                "Cp conversion is available for direct TDB Gibbs polynomials. "
                "Referenced functions are listed in the Function Overview."
            )
        )
        layout.addStretch(1)
        return tab
    if editor_refs is not None:
        editor_refs.clear()
        editor_refs["rows"] = []

    h_editor = _editable_line("" if h298 is None else _thermo_number(h298))
    s_editor = _editable_line("" if s298 is None else _thermo_number(s298))
    h_editor.setObjectName("FunctionCpH298Input")
    s_editor.setObjectName("FunctionCpS298Input")
    h_editor.setFixedWidth(_THERMO_COEFFICIENT_BOX_WIDTH)
    s_editor.setFixedWidth(_THERMO_COEFFICIENT_BOX_WIDTH)
    if editor_refs is not None:
        editor_refs["h298"] = h_editor
        editor_refs["s298"] = s_editor

    reference_row = QHBoxLayout()
    reference_row.setSpacing(36)
    reference_row.addStretch(1)
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
    layout.addLayout(reference_row)
    layout.addWidget(_expression_basis_label(_cp_expression_text()))

    for index, cp_row in enumerate(cp_rows, start=1):
        if len(cp_rows) > 1:
            range_label = QLabel(f"Range {index}")
            range_label.setObjectName("DatabaseListLabel")
            layout.addWidget(range_label)
        row_refs: dict[str, Any] = {}
        if editor_refs is not None:
            row_refs["lower"] = None
            row_refs["upper"] = None
        if editor_refs is not None:
            editor_refs["rows"].append(row_refs)
        layout.addWidget(
            _coefficient_range_editor(
                cp_row[0],
                cp_row[1],
                _cp_term_rows(cp_row),
                show_temperature_row=index > 1,
                editor_refs=row_refs if editor_refs is not None else None,
                compact_fields=True,
            )
        )
    layout.addStretch(1)
    return tab


def _thermo_reference_editor_column(
    label_text: str,
    unit_text: str,
    editor: QLineEdit,
) -> QWidget:
    """Return a two-row reference-state editor column."""
    wrapper = QWidget()
    layout = QVBoxLayout(wrapper)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.setSpacing(6)

    header = QHBoxLayout()
    header.setContentsMargins(0, 0, 0, 0)
    header.setSpacing(6)
    label = _latex_label(label_text)
    unit = _latex_label(unit_text)
    header.addWidget(label)
    header.addWidget(unit)
    header.addStretch(1)

    layout.addLayout(header)
    layout.addWidget(editor, alignment=Qt.AlignmentFlag.AlignLeft)
    return wrapper


def _function_coefficient_tab(
    data_rows: list[list[float]],
    term_builder: Any,
    expression_text: str,
    *,
    editor_refs: dict[str, Any] | None = None,
    show_temperature_rows: bool = True,
    balance_callbacks: dict[int, Any] | None = None,
) -> QWidget:
    """Build a generic coefficient/power expression tab."""
    tab = QWidget()
    layout = QVBoxLayout(tab)
    layout.setContentsMargins(10, 10, 10, 10)
    layout.setSpacing(12)
    if not data_rows:
        layout.addWidget(
            QLabel(
                "Expression editing is available for direct TDB Gibbs polynomials. "
                "Referenced functions are listed in the Function Overview."
            )
        )
        layout.addStretch(1)
        return tab
    if editor_refs is not None:
        editor_refs.clear()
        editor_refs["rows"] = []
    for index, data_row in enumerate(data_rows, start=1):
        if len(data_rows) > 1:
            range_label = QLabel(f"Range {index}")
            range_label.setObjectName("DatabaseListLabel")
            layout.addWidget(range_label)
        row_refs: dict[str, Any] = {}
        if show_temperature_rows:
            range_layout, lower_editor, upper_editor = _temperature_range_row(
                data_row[0],
                data_row[1],
            )
            row_refs["lower"] = lower_editor
            row_refs["upper"] = upper_editor
            layout.addLayout(range_layout)
        else:
            row_refs["lower"] = None
            row_refs["upper"] = None
        if editor_refs is not None:
            editor_refs["rows"].append(row_refs)
        layout.addWidget(_expression_basis_label(expression_text))
        layout.addWidget(
            _coefficient_range_editor(
                data_row[0],
                data_row[1],
                term_builder(data_row),
                show_temperature_row=False,
                editor_refs=row_refs if editor_refs is not None else None,
                compact_fields=True,
                balance_callbacks=balance_callbacks,
            )
        )
    layout.addStretch(1)
    return tab


def _capture_cp_editor(
    editor_refs: dict[str, Any],
    fallback_rows: list[list[float]],
) -> tuple[float, float, list[list[float]]] | None:
    """Read H298, S298, and Cp rows from the current Cp editor."""
    row_refs = editor_refs.get("rows", [])
    if not row_refs:
        return None
    h298 = _editor_float(editor_refs.get("h298"), 0.0)
    s298 = _editor_float(editor_refs.get("s298"), 0.0)
    rows: list[list[float]] = []
    for index, row_ref in enumerate(row_refs):
        fallback = (
            fallback_rows[index]
            if index < len(fallback_rows)
            else [298.15, 6000]
        )
        lower = _editor_float(row_ref.get("lower"), fallback[0])
        upper = _editor_float(row_ref.get("upper"), fallback[1])
        terms = row_ref.get("terms", [])
        coefficients = _capture_fixed_and_custom_terms(terms, fixed_count=4)
        rows.append([lower, upper, *coefficients])
    return h298, s298, rows


def _capture_expression_rows(
    editor_refs: dict[str, Any],
    mode: str,
    fallback_rows: list[list[float]],
) -> list[list[float]] | None:
    """Read coefficient rows from a G/H/S editor."""
    row_refs = editor_refs.get("rows", [])
    if not row_refs:
        return None
    fixed_count = 6 if mode == "G" else 5
    rows: list[list[float]] = []
    for index, row_ref in enumerate(row_refs):
        fallback = (
            fallback_rows[index]
            if index < len(fallback_rows)
            else [298.15, 6000]
        )
        lower = _editor_float(row_ref.get("lower"), fallback[0])
        upper = _editor_float(row_ref.get("upper"), fallback[1])
        coefficients = _capture_fixed_and_custom_terms(
            row_ref.get("terms", []),
            fixed_count=fixed_count,
        )
        rows.append([lower, upper, *coefficients])
    return rows


def _capture_fixed_and_custom_terms(
    terms: list[dict[str, Any]],
    *,
    fixed_count: int,
) -> list[float]:
    """Return fixed coefficients followed by four custom coefficient/power pairs."""
    fixed = [
        _editor_float(term.get("coefficient"), 0.0)
        for term in terms[:fixed_count]
    ]
    custom: list[float] = []
    for term in terms[fixed_count : fixed_count + 4]:
        coefficient = _editor_float(term.get("coefficient"), 0.0)
        exponent = term.get("fixed_exponent")
        if exponent is None:
            exponent = _editor_float(term.get("exponent"), 0.0)
        else:
            exponent = _safe_float(exponent, 0.0)
        custom.extend([coefficient, exponent])
    while len(custom) < 8:
        custom.extend([0.0, 0.0])
    return [*fixed, *custom]


def _gibbs_rows_from_enthalpy_rows(
    enthalpy_rows: list[list[float]],
    previous_gibbs_rows: list[list[float]],
) -> list[list[float]]:
    """Update Gibbs rows from editable H coefficients where the inverse is defined."""
    rows: list[list[float]] = []
    for index, enthalpy in enumerate(enthalpy_rows):
        previous = (
            previous_gibbs_rows[index]
            if index < len(previous_gibbs_rows)
            else [enthalpy[0], enthalpy[1], *([0.0] * 14)]
        )
        gibbs = [float(value) for value in previous]
        gibbs[0] = enthalpy[0]
        gibbs[1] = enthalpy[1]
        if len(enthalpy) >= 7:
            gibbs[2] = enthalpy[2]
            gibbs[4] = -enthalpy[3]
            gibbs[5] = -enthalpy[4]
            gibbs[6] = -0.5 * enthalpy[5]
            gibbs[7] = 0.5 * enthalpy[6]
        _copy_enthalpy_custom_terms_to_gibbs(enthalpy, gibbs)
        rows.append(gibbs)
    return rows


def _gibbs_rows_from_entropy_rows(
    entropy_rows: list[list[float]],
    previous_gibbs_rows: list[list[float]],
) -> list[list[float]]:
    """Update Gibbs rows from editable S coefficients where the inverse is defined."""
    rows: list[list[float]] = []
    for index, entropy in enumerate(entropy_rows):
        previous = (
            previous_gibbs_rows[index]
            if index < len(previous_gibbs_rows)
            else [entropy[0], entropy[1], *([0.0] * 14)]
        )
        gibbs = [float(value) for value in previous]
        gibbs[0] = entropy[0]
        gibbs[1] = entropy[1]
        if len(entropy) >= 7:
            gibbs[4] = -entropy[3]
            gibbs[3] = -entropy[2] - gibbs[4]
            gibbs[5] = -0.5 * entropy[4]
            gibbs[6] = -entropy[5] / 3.0
            gibbs[7] = entropy[6]
        _copy_entropy_custom_terms_to_gibbs(entropy, gibbs)
        rows.append(gibbs)
    return rows


def _copy_enthalpy_custom_terms_to_gibbs(
    enthalpy: list[float],
    gibbs: list[float],
) -> None:
    for index in range(4):
        source_index = 7 + 2 * index
        target_index = 8 + 2 * index
        if source_index + 1 >= len(enthalpy):
            break
        coefficient = enthalpy[source_index]
        power = enthalpy[source_index + 1]
        gibbs[target_index + 1] = power
        denominator = 1.0 - power
        gibbs[target_index] = (
            coefficient / denominator if abs(denominator) > 1e-14 else 0.0
        )


def _copy_entropy_custom_terms_to_gibbs(
    entropy: list[float],
    gibbs: list[float],
) -> None:
    for index in range(4):
        source_index = 7 + 2 * index
        target_index = 8 + 2 * index
        if source_index + 1 >= len(entropy):
            break
        coefficient = entropy[source_index]
        gibbs_power = entropy[source_index + 1] + 1.0
        gibbs[target_index + 1] = gibbs_power
        gibbs[target_index] = (
            -coefficient / gibbs_power if abs(gibbs_power) > 1e-14 else 0.0
        )


def _gibbs_rows_to_tdb_expression(gibbs_rows: list[list[float]]) -> str:
    """Serialize G2HSCp Gibbs coefficient rows back to a direct TDB expression."""
    if not gibbs_rows:
        return ""
    expression = (
        f"{_compact_number(gibbs_rows[0][0])} "
        f"{_gibbs_row_expression(gibbs_rows[0])}; "
    )
    for index, row in enumerate(gibbs_rows):
        is_last = index == len(gibbs_rows) - 1
        expression += f"{_compact_number(row[1])} {'N' if is_last else 'Y'}"
        if not is_last:
            expression += f" {_gibbs_row_expression(gibbs_rows[index + 1])}; "
    return expression


def _set_function_gibbs_rows(
    function: FunctionDefinition,
    gibbs_rows: list[list[float]],
) -> None:
    """Store sorted Gibbs rows and synchronized temperature ranges on a Function."""
    normalized_rows, _selected_index = _normalized_gibbs_rows(gibbs_rows)
    function.expression = _gibbs_rows_to_tdb_expression(normalized_rows)
    function.temperature_ranges = [
        (float(row[0]), float(row[1]))
        for row in normalized_rows
    ]
    function.gibbs_ranges = [
        GibbsRange(
            T_min=float(row[0]),
            T_max=float(row[1]),
            Gibbs=_gibbs_row_expression(row),
            status="N" if index == len(normalized_rows) - 1 else "Y",
        )
        for index, row in enumerate(normalized_rows)
    ]


def _normalized_gibbs_rows(
    gibbs_rows: list[list[float]],
    *,
    selected_index: int = 0,
) -> tuple[list[list[float]], int]:
    """Sort ranges by upper temperature and make adjacent ranges continuous."""
    indexed_rows = [
        (index, [float(value) for value in row])
        for index, row in enumerate(gibbs_rows)
    ]
    if not indexed_rows:
        return [], 0
    initial_lower = min(row[0] for _index, row in indexed_rows)
    indexed_rows.sort(key=lambda item: item[1][1])
    normalized: list[list[float]] = []
    selected_new_index = 0
    previous_upper: float | None = initial_lower
    for new_index, (old_index, row) in enumerate(indexed_rows):
        row[0] = previous_upper
        if row[1] < row[0]:
            row[1] = row[0]
        previous_upper = row[1]
        if old_index == selected_index:
            selected_new_index = new_index
        normalized.append(row)
    return normalized, selected_new_index


def _normalize_function_temperature_ranges(function: FunctionDefinition) -> None:
    """Normalize a Function's direct Gibbs expression and Cp range metadata."""
    gibbs_rows = _direct_function_gibbs_rows(function)
    if gibbs_rows:
        _set_function_gibbs_rows(function, gibbs_rows)
        return
    ranges = sorted(
        [
            (float(lower), float(upper))
            for lower, upper in function.temperature_ranges
        ],
        key=lambda item: item[1],
    )
    normalized_ranges: list[tuple[float, float]] = []
    previous_upper: float | None = min(
        (lower for lower, _upper in ranges),
        default=None,
    )
    for lower, upper in ranges:
        if previous_upper is not None:
            lower = previous_upper
        if upper < lower:
            upper = lower
        normalized_ranges.append((lower, upper))
        previous_upper = upper
    function.temperature_ranges = normalized_ranges


def _function_by_name(
    database: DatabaseIR,
    name: str,
) -> FunctionDefinition | None:
    for function in database.functions:
        if function.name == name:
            return function
    return None


def _gibbs_row_expression(row: list[float]) -> str:
    """Return a compact Gibbs polynomial expression for one coefficient row."""
    coefficients = row[2:]
    terms: list[str] = []
    _append_polynomial_term(terms, coefficients[0], "")
    _append_polynomial_term(terms, coefficients[1], "*T")
    _append_polynomial_term(terms, coefficients[2], "*T*LN(T)")
    _append_polynomial_term(terms, coefficients[3], "*T**2")
    _append_polynomial_term(terms, coefficients[4], "*T**3")
    _append_polynomial_term(terms, coefficients[5], "*T**(-1)")
    for index in range(4):
        coefficient = coefficients[6 + 2 * index]
        power = coefficients[7 + 2 * index]
        if coefficient == 0.0:
            continue
        if _is_log_t_power(power):
            _append_polynomial_term(terms, coefficient, "*LN(T)")
            continue
        _append_polynomial_term(
            terms,
            coefficient,
            f"*T**({_compact_number(power)})",
        )
    return "".join(terms).lstrip("+") or "0"


def _append_polynomial_term(terms: list[str], coefficient: float, suffix: str) -> None:
    if coefficient == 0.0:
        return
    terms.append(f"{_signed_tdb_number(coefficient)}{suffix}")


def _signed_tdb_number(value: float) -> str:
    """Format a signed thermodynamic coefficient without losing precision."""
    return f"{float(value):+.16g}"


def _editor_float(editor: Any, fallback: float) -> float:
    if isinstance(editor, QLineEdit):
        return _safe_float(editor.text(), fallback)
    return _safe_float(editor, fallback)


def _safe_float(value: Any, fallback: float) -> float:
    text = str(value).strip()
    if text == "":
        return float(fallback)
    try:
        return float(text)
    except (TypeError, ValueError):
        return float(fallback)


def _coefficient_range_editor(
    lower: float,
    upper: float,
    term_rows: list[tuple[str, str, str, str, bool]],
    *,
    show_temperature_row: bool = True,
    editor_refs: dict[str, Any] | None = None,
    compact_fields: bool = False,
    balance_callbacks: dict[int, Any] | None = None,
) -> QWidget:
    """Return one editable coefficient/power range editor."""
    wrapper = QWidget()
    if compact_fields:
        outer_layout = QHBoxLayout(wrapper)
        outer_layout.setContentsMargins(0, 0, 0, 0)
        outer_layout.setSpacing(0)
        grid_holder = QWidget()
        grid = QGridLayout(grid_holder)
        outer_layout.addStretch(1)
        outer_layout.addWidget(grid_holder)
        outer_layout.addStretch(1)
    else:
        grid = QGridLayout(wrapper)
    grid.setContentsMargins(0, 0, 0, 0)
    grid.setHorizontalSpacing(8)
    grid.setVerticalSpacing(6)
    balance_column_offset = 1 if balance_callbacks else 0
    grid.setColumnStretch(1 + balance_column_offset, 0 if compact_fields else 2)
    grid.setColumnStretch(3 + balance_column_offset, 0)
    if compact_fields:
        grid.setColumnStretch(4 + balance_column_offset, 1)

    row_offset = 1 if show_temperature_row else 0
    if show_temperature_row:
        lower_editor = _editable_line(_compact_number(lower))
        upper_editor = _editable_line(_compact_number(upper))
        if compact_fields:
            lower_editor.setFixedWidth(_THERMO_TEMPERATURE_BOX_WIDTH)
            upper_editor.setFixedWidth(_THERMO_TEMPERATURE_BOX_WIDTH)
        grid.addWidget(QLabel("from"), 0, 0)
        grid.addWidget(lower_editor, 0, 1)
        grid.addWidget(QLabel("[K] to"), 0, 2)
        grid.addWidget(upper_editor, 0, 3)
        grid.addWidget(QLabel("[K]"), 0, 4)
        if editor_refs is not None:
            editor_refs["lower"] = lower_editor
            editor_refs["upper"] = upper_editor
    if editor_refs is not None:
        editor_refs["terms"] = []

    for offset, (
        name,
        coefficient,
        term,
        exponent,
        fixed_exponent,
    ) in enumerate(term_rows):
        row = offset + row_offset
        grid.addWidget(QLabel(f"{name}: "), row, 0)
        balance_button: QToolButton | None = None
        if balance_callbacks:
            balance_button = QToolButton()
            balance_button.setObjectName("GibbsBalanceCoefficientButton")
            balance_button.setText("Bal.")
            balance_button.setToolTip(
                "Balance Gibbs continuity at transition to 1e-10 J/mol"
            )
            balance_button.setCursor(Qt.CursorShape.PointingHandCursor)
            balance_button.setVisible(False)
            balance_button.setEnabled(False)
            callback = balance_callbacks.get(offset)
            if callback is not None:
                balance_button.clicked.connect(callback)
            grid.addWidget(balance_button, row, 1)
        coefficient_editor = _editable_line(coefficient)
        if compact_fields:
            coefficient_editor.setFixedWidth(_THERMO_COEFFICIENT_BOX_WIDTH)
        grid.addWidget(coefficient_editor, row, 1 + balance_column_offset)
        grid.addWidget(QLabel(term), row, 2 + balance_column_offset)
        exponent_widget = (
            _fixed_exponent_label(exponent)
            if fixed_exponent
            else _editable_line(exponent)
        )
        if compact_fields and isinstance(exponent_widget, QLineEdit):
            exponent_widget.setFixedWidth(_THERMO_POWER_BOX_WIDTH)
        grid.addWidget(exponent_widget, row, 3 + balance_column_offset)
        if editor_refs is not None:
            editor_refs["terms"].append(
                {
                    "name": name,
                    "coefficient": coefficient_editor,
                    "exponent": None if fixed_exponent else exponent_widget,
                    "fixed_exponent": exponent if fixed_exponent else None,
                    "balance_button": balance_button,
                }
            )
    return wrapper


def _temperature_range_row(
    lower: float,
    upper: float,
) -> tuple[QHBoxLayout, QLineEdit, QLineEdit]:
    """Return a compact temperature range row."""
    row = QHBoxLayout()
    row.setSpacing(8)
    row.addStretch(1)
    row.addWidget(QLabel("Range from"))
    lower_editor = _editable_line(_compact_number(lower))
    upper_editor = _editable_line(_compact_number(upper))
    lower_editor.setFixedWidth(_THERMO_TEMPERATURE_BOX_WIDTH)
    upper_editor.setFixedWidth(_THERMO_TEMPERATURE_BOX_WIDTH)
    row.addWidget(lower_editor)
    row.addWidget(QLabel("[K] to"))
    row.addWidget(upper_editor)
    row.addWidget(QLabel("[K]"))
    row.addStretch(1)
    return row, lower_editor, upper_editor


def _fixed_exponent_label(exponent: str) -> QLabel:
    """Return a static label for non-editable basis exponents."""
    label = QLabel(exponent)
    label.setObjectName("FixedExponentLabel")
    label.setAlignment(
        Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft
    )
    return label


def _expression_basis_label(text: str) -> QLabel:
    """Return a compact expression-basis label."""
    label = _math_label(text)
    label.setAlignment(Qt.AlignmentFlag.AlignCenter)
    label.setWordWrap(True)
    font = label.font()
    font.setPointSize(font.pointSize() + 2)
    label.setFont(font)
    label.setObjectName("ExpressionBasisLabel")
    return label


def _cp_expression_text() -> str:
    """Return the Cp expression basis used by G2HSCp output."""
    return (
        "C<sub>p</sub>(T) = "
        "a + bT + cT<sup>-2</sup> "
        "+ dT<sup>2</sup> "
        "+ c<sub>1</sub>T<sup>p1</sup> "
        "+ c<sub>2</sub>T<sup>p2</sup> "
        "+ c<sub>3</sub>T<sup>p3</sup> "
        "+ c<sub>4</sub>T<sup>p4</sup>"
    )


def _gibbs_expression_text() -> str:
    """Return the Gibbs expression basis used by G2HSCp input."""
    return (
        "G(T) = "
        "a + bT + cT ln(T) "
        "+ dT<sup>2</sup> + eT<sup>3</sup> "
        "+ fT<sup>-1</sup> "
        "+ c<sub>1</sub>T<sup>p1</sup> "
        "+ c<sub>2</sub>T<sup>p2</sup> "
        "+ c<sub>3</sub>T<sup>p3</sup> "
        "+ c<sub>4</sub>T<sup>p4</sup>"
    )


def _enthalpy_expression_text() -> str:
    """Return the H expression basis derived from G."""
    return (
        "H(T) = "
        "a + bT + cT<sup>2</sup> "
        "+ dT<sup>3</sup> + eT<sup>-1</sup> "
        "+ c<sub>1</sub>T<sup>p1</sup> "
        "+ c<sub>2</sub>T<sup>p2</sup> "
        "+ c<sub>3</sub>T<sup>p3</sup> "
        "+ c<sub>4</sub>T<sup>p4</sup>"
    )


def _entropy_expression_text() -> str:
    """Return the S expression basis derived from G."""
    return (
        "S(T) = "
        "a + b ln(T) + cT "
        "+ dT<sup>2</sup> + eT<sup>-2</sup> "
        "+ c<sub>1</sub>T<sup>p1</sup> "
        "+ c<sub>2</sub>T<sup>p2</sup> "
        "+ c<sub>3</sub>T<sup>p3</sup> "
        "+ c<sub>4</sub>T<sup>p4</sup>"
    )


def _cp_term_rows(cp_row: list[float]) -> list[tuple[str, str, str, str, bool]]:
    """Return Cp coefficient rows."""
    rows = [
        (name, _thermo_number(coefficient), "T^", f"{power:.2f}", True)
        for name, coefficient, power in zip(
            ["a", "b", "c", "d"],
            cp_row[2:6],
            [0.0, 1.0, -2.0, 2.0],
            strict=True,
        )
    ]
    rows.extend(_custom_power_rows(cp_row[6:14]))
    return rows


def _gibbs_term_rows(gibbs_row: list[float]) -> list[tuple[str, str, str, str, bool]]:
    """Return G coefficient rows from G2HSCp Gibbs input basis."""
    coefficients = gibbs_row[2:]
    rows = [
        ("a", _thermo_number(coefficients[0]), "T^", "0.00", True),
        ("b", _thermo_number(coefficients[1]), "T^", "1.00", True),
        ("c", _thermo_number(coefficients[2]), "T*ln(T)", "", True),
        ("d", _thermo_number(coefficients[3]), "T^", "2.00", True),
        ("e", _thermo_number(coefficients[4]), "T^", "3.00", True),
        ("f", _thermo_number(coefficients[5]), "T^", "-1.00", True),
    ]
    rows.extend(_custom_power_rows(coefficients[6:14]))
    return rows


def _enthalpy_term_rows(
    gibbs_row: list[float],
) -> list[tuple[str, str, str, str, bool]]:
    """Return H coefficient rows derived from G = H - TS."""
    coefficients = gibbs_row[2:]
    rows = [
        ("a", _thermo_number(coefficients[0]), "T^", "0.00", True),
        ("b", _thermo_number(-coefficients[2]), "T^", "1.00", True),
        ("c", _thermo_number(-coefficients[3]), "T^", "2.00", True),
        ("d", _thermo_number(-2.0 * coefficients[4]), "T^", "3.00", True),
        ("e", _thermo_number(2.0 * coefficients[5]), "T^", "-1.00", True),
    ]
    rows.extend(
        _derived_custom_power_rows(
            coefficients[6:14],
            coefficient_factor=lambda power: 1.0 - power,
            power_shift=0.0,
        )
    )
    return rows


def _entropy_term_rows(gibbs_row: list[float]) -> list[tuple[str, str, str, str, bool]]:
    """Return S coefficient rows derived from S = -dG/dT."""
    coefficients = gibbs_row[2:]
    rows = [
        (
            "a",
            _thermo_number(-(coefficients[1] + coefficients[2])),
            "T^",
            "0.00",
            True,
        ),
        ("b", _thermo_number(-coefficients[2]), "ln(T)", "", True),
        ("c", _thermo_number(-2.0 * coefficients[3]), "T^", "1.00", True),
        ("d", _thermo_number(-3.0 * coefficients[4]), "T^", "2.00", True),
        ("e", _thermo_number(coefficients[5]), "T^", "-2.00", True),
    ]
    rows.extend(
        _derived_custom_power_rows(
            coefficients[6:14],
            coefficient_factor=lambda power: -power,
            power_shift=-1.0,
        )
    )
    return rows


def _custom_power_rows(values: list[float]) -> list[tuple[str, str, str, str, bool]]:
    """Return four editable custom coefficient/power rows."""
    rows: list[tuple[str, str, str, str, bool]] = []
    for index in range(4):
        coefficient = values[2 * index]
        power = values[2 * index + 1]
        name = f"c{index + 1}"
        if coefficient == 0.0:
            rows.append((name, "", "T^", "", False))
        else:
            rows.append(
                (
                    name,
                    _thermo_number(coefficient),
                    "T^",
                    _compact_number(power),
                    False,
                )
            )
    return rows


def _derived_custom_power_rows(
    values: list[float],
    *,
    coefficient_factor: Any,
    power_shift: float,
) -> list[tuple[str, str, str, str, bool]]:
    """Return custom rows after derivative/integral power transform."""
    rows: list[tuple[str, str, str, str, bool]] = []
    for index in range(4):
        coefficient = values[2 * index]
        power = values[2 * index + 1]
        name = f"c{index + 1}"
        if coefficient == 0.0:
            rows.append((name, "", "T^", "", False))
        else:
            rows.append(
                (
                    name,
                    _thermo_number(coefficient * coefficient_factor(power)),
                    "T^",
                    _compact_number(power + power_shift),
                    False,
                )
            )
    return rows


def _function_thermo_data(
    function: FunctionDefinition,
    database: DatabaseIR | None = None,
    seen: set[str] | None = None,
) -> tuple[float | None, float | None, list[list[float]], list[list[float]]]:
    """Convert direct or pertinent-function Gibbs records through G2HSCp."""
    gibbs_rows = _direct_function_gibbs_rows(function)
    if not gibbs_rows:
        return _referenced_function_thermo_data(function, database, seen)
    try:
        h298, s298, cp = G2HSCp(gibbs_rows)
    except (AssertionError, ValueError, TypeError):
        return _referenced_function_thermo_data(function, database, seen)
    return float(h298), float(s298), cp.tolist(), gibbs_rows


def _direct_function_gibbs_rows(function: FunctionDefinition) -> list[list[float]]:
    """Return direct numeric Gibbs rows, preferring parsed DatabaseIR ranges."""
    return _function_gibbs_rows_from_ranges(function) or _function_gibbs_rows(
        function.expression,
    )


def _function_gibbs_rows_from_ranges(
    function: FunctionDefinition,
) -> list[list[float]]:
    """Convert expanded DatabaseIR Gibbs ranges into G2HSCp row format."""
    rows: list[list[float]] = []
    for gibbs_range in function.gibbs_ranges:
        coefficients = _gibbs_coefficients_from_segment(gibbs_range.Gibbs)
        if coefficients is None:
            return []
        rows.append(
            [
                float(gibbs_range.T_min),
                float(gibbs_range.T_max),
                *coefficients,
            ],
        )
    return rows


def _referenced_function_thermo_data(
    function: FunctionDefinition,
    database: DatabaseIR | None,
    seen: set[str] | None,
) -> tuple[float | None, float | None, list[list[float]], list[list[float]]]:
    """Combine referenced functions by Neumann-Kopp over a split range grid."""
    if database is None:
        return None, None, [], []
    seen = set() if seen is None else set(seen)
    function_key = function.name.upper()
    if function_key in seen:
        return None, None, [], []
    seen.add(function_key)

    dependencies = _function_dependency_terms(database, function)
    dataset: dict[str, dict[str, Any]] = {}
    names: list[str] = []
    multipliers: list[float] = []
    local_gibbs_rows = _function_local_gibbs_rows(function.expression, dependencies)
    if local_gibbs_rows:
        try:
            h298, s298, cp_rows = G2HSCp(local_gibbs_rows)
        except (AssertionError, ValueError, TypeError):
            return None, None, [], []
        local_name = f"__{function.name}_local__"
        names.append(local_name)
        multipliers.append(1.0)
        dataset[local_name] = {
            "H298": float(h298),
            "S298": float(s298),
            "Cp": np.asarray(cp_rows, dtype=float),
        }

    for referenced_function, multiplier in dependencies:
        h298, s298, cp_rows, _gibbs_rows = _function_thermo_data(
            referenced_function,
            database,
            seen,
        )
        if h298 is None or s298 is None or not cp_rows:
            return None, None, [], []
        names.append(referenced_function.name)
        multipliers.append(multiplier)
        dataset[referenced_function.name] = {
            "H298": float(h298),
            "S298": float(s298),
            "Cp": np.asarray(cp_rows, dtype=float),
        }
    if not names:
        return None, None, [], []

    try:
        aligned_dataset = _align_neumann_kopp_temperature_ranges(dataset, names)
        combined = NeumanKoppHSCp(aligned_dataset, names, multipliers)
        cp_rows_array = np.asarray(combined["Cp"], dtype=float)
        gibbs_rows_array = HSCp2G(
            float(combined["H298"]),
            float(combined["S298"]),
            cp_rows_array,
        )
        cp_rows_array, gibbs_rows_array = _merge_identical_thermo_rows(
            cp_rows_array,
            gibbs_rows_array,
        )
    except (AssertionError, ValueError, TypeError, IndexError):
        return None, None, [], []
    return (
        float(combined["H298"]),
        float(combined["S298"]),
        cp_rows_array.tolist(),
        gibbs_rows_array.tolist(),
    )


def _merge_identical_thermo_rows(
    cp_rows: np.ndarray,
    gibbs_rows: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Merge adjacent ranges when Cp and Gibbs coefficient rows are unchanged."""
    cp_rows = np.asarray(cp_rows, dtype=float)
    gibbs_rows = np.asarray(gibbs_rows, dtype=float)
    if cp_rows.ndim != 2 or gibbs_rows.ndim != 2:
        return cp_rows, gibbs_rows
    if len(cp_rows) <= 1 or len(cp_rows) != len(gibbs_rows):
        return cp_rows, gibbs_rows

    merged_cp: list[np.ndarray] = [cp_rows[0].copy()]
    merged_gibbs: list[np.ndarray] = [gibbs_rows[0].copy()]
    for cp_row, gibbs_row in zip(cp_rows[1:], gibbs_rows[1:], strict=True):
        cp_previous = merged_cp[-1]
        gibbs_previous = merged_gibbs[-1]
        same_boundary = np.isclose(
            cp_previous[1],
            cp_row[0],
            rtol=1e-10,
            atol=1e-8,
        ) and np.isclose(
            gibbs_previous[1],
            gibbs_row[0],
            rtol=1e-10,
            atol=1e-8,
        )
        same_cp = np.allclose(
            cp_previous[2:],
            cp_row[2:],
            rtol=1e-10,
            atol=1e-8,
        )
        same_gibbs = np.allclose(
            gibbs_previous[2:],
            gibbs_row[2:],
            rtol=1e-10,
            atol=1e-8,
        )
        if same_boundary and same_cp and same_gibbs:
            cp_previous[1] = cp_row[1]
            gibbs_previous[1] = gibbs_row[1]
        else:
            merged_cp.append(cp_row.copy())
            merged_gibbs.append(gibbs_row.copy())
    return np.vstack(merged_cp), np.vstack(merged_gibbs)


def _function_dependency_terms(
    database: DatabaseIR,
    function: FunctionDefinition,
) -> list[tuple[FunctionDefinition, float]]:
    """Return referenced functions and simple multipliers, excluding self."""
    functions_by_name = {
        referenced_function.name.upper(): referenced_function
        for referenced_function in database.functions
    }
    dependencies: list[tuple[FunctionDefinition, float]] = []
    for name in _referenced_function_names(function.expression, functions_by_name):
        referenced_function = functions_by_name[name]
        if referenced_function.name.upper() == function.name.upper():
            continue
        factor = _expression_multiplier_for_name(
            function.expression,
            referenced_function.name,
        )
        if factor is not None and factor != 0.0:
            dependencies.append((referenced_function, factor))
    return dependencies


def _referenced_function_names(
    expression: str,
    functions_by_name: dict[str, FunctionDefinition],
) -> list[str]:
    """Return function names referenced as expression tokens."""
    if not expression or not functions_by_name:
        return []
    names: list[str] = []
    seen: set[str] = set()
    for match in _FUNCTION_TOKEN_RE.finditer(expression):
        name = match.group(1).upper()
        if name in seen or name not in functions_by_name:
            continue
        seen.add(name)
        names.append(name)
    return names


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


def _expression_references_name(expression: str, name: str) -> bool:
    """Return whether an expression references a named thermodynamic function."""
    if not expression or not name:
        return False
    pattern = rf"(?<![A-Za-z0-9_]){re.escape(name)}(?![A-Za-z0-9_])"
    return re.search(pattern, expression, re.IGNORECASE) is not None


def _function_local_gibbs_rows(
    expression: str,
    dependencies: list[tuple[FunctionDefinition, float]],
) -> list[list[float]]:
    """Return direct Gibbs rows after removing referenced-function terms."""
    segments = _function_expression_segments(expression)
    if not segments:
        return []
    names = [function.name for function, _factor in dependencies]
    rows: list[list[float]] = []
    for lower, upper, segment_expression in segments:
        local_expression = _remove_referenced_function_terms(segment_expression, names)
        coefficients = _gibbs_coefficients_from_segment(local_expression)
        if coefficients is None:
            return []
        rows.append([lower, upper, *coefficients])
    return rows


def _remove_referenced_function_terms(expression: str, names: list[str]) -> str:
    """Remove additive terms that contain referenced function names."""
    local_terms: list[str] = []
    for term in _signed_expression_terms(expression):
        if any(_expression_references_name(term, name) for name in names):
            continue
        local_terms.append(term)
    return "".join(local_terms) or "0"


def _align_neumann_kopp_temperature_ranges(
    dataset: dict[str, dict[str, Any]],
    names: list[str],
) -> dict[str, dict[str, Any]]:
    """Split all referenced Cp arrays onto the union of their temperature bounds."""
    bounds: list[float] = []
    for name in names:
        cp_rows = np.asarray(dataset[name]["Cp"], dtype=float)
        bounds.extend(float(value) for value in cp_rows[:, :2].flatten())
    unique_bounds = sorted({round(value, 6) for value in bounds})
    if len(unique_bounds) < 2:
        return dataset

    aligned: dict[str, dict[str, Any]] = {}
    for name in names:
        original_cp = np.asarray(dataset[name]["Cp"], dtype=float)
        aligned_rows: list[np.ndarray] = []
        for lower, upper in zip(unique_bounds[:-1], unique_bounds[1:], strict=True):
            source_row = _cp_row_for_temperature_range(original_cp, lower, upper)
            row = source_row.copy()
            row[0] = lower
            row[1] = upper
            aligned_rows.append(row)
        aligned[name] = {
            "H298": dataset[name]["H298"],
            "S298": dataset[name]["S298"],
            "Cp": np.vstack(aligned_rows),
        }
    return aligned


def _cp_row_for_temperature_range(
    cp_rows: np.ndarray,
    lower: float,
    upper: float,
) -> np.ndarray:
    """Return the source Cp row that should represent one split interval."""
    midpoint = 0.5 * (lower + upper)
    tolerance = 1e-8
    for row in cp_rows:
        if row[0] - tolerance <= midpoint <= row[1] + tolerance:
            return row
    if midpoint < cp_rows[0, 0]:
        return cp_rows[0]
    return cp_rows[-1]


def _function_gibbs_rows(expression: str) -> list[list[float]]:
    """Parse direct TDB Gibbs polynomial segments for G2HSCp."""
    segments = _function_expression_segments(expression)
    rows: list[list[float]] = []
    for lower, upper, segment_expression in segments:
        coefficients = _gibbs_coefficients_from_segment(segment_expression)
        if coefficients is None:
            return []
        rows.append([lower, upper, *coefficients])
    return rows


def _function_expression_segments(expression: str) -> list[tuple[float, float, str]]:
    """Split a TDB function expression into lower/upper/expression segments."""
    lower_match = re.match(
        r"\s*([-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?)",
        expression,
    )
    if lower_match is None:
        return []
    lower = float(lower_match.group(1).replace("D", "E").replace("d", "E"))
    remainder = expression[lower_match.end() :]
    segments: list[tuple[float, float, str]] = []
    while ";" in remainder:
        segment_text, remainder = remainder.split(";", 1)
        upper_match = re.match(
            r"\s*([-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?)",
            remainder,
        )
        if upper_match is None:
            break
        upper = float(upper_match.group(1).replace("D", "E").replace("d", "E"))
        segments.append((lower, upper, segment_text.strip()))
        lower = upper
        remainder = remainder[upper_match.end() :].lstrip()
        if remainder[:1].upper() in {"Y", "N"}:
            remainder = remainder[1:].lstrip()
    return segments


def _gibbs_coefficients_from_segment(expression: str) -> list[float] | None:
    """Return G2HSCp coefficients for one direct Gibbs expression segment."""
    coefficients = [0.0] * 14
    custom_pairs: list[tuple[float, float]] = []
    for term in _signed_expression_terms(expression):
        parsed = _parse_gibbs_term(term)
        if parsed is None:
            return None
        power, coefficient = parsed
        if power == "TLOGT":
            coefficients[2] += coefficient
        elif abs(power) < 1e-12:
            coefficients[0] += coefficient
        elif abs(power - 1.0) < 1e-12:
            coefficients[1] += coefficient
        elif abs(power - 2.0) < 1e-12:
            coefficients[3] += coefficient
        elif abs(power - 3.0) < 1e-12:
            coefficients[4] += coefficient
        elif abs(power + 1.0) < 1e-12:
            coefficients[5] += coefficient
        else:
            for index, (_existing_coefficient, existing_power) in enumerate(
                custom_pairs,
            ):
                if abs(power - existing_power) < 1e-12:
                    custom_pairs[index] = (
                        custom_pairs[index][0] + coefficient,
                        existing_power,
                    )
                    break
            else:
                custom_pairs.append((coefficient, power))
    if len(custom_pairs) > 4:
        return None
    for index, (coefficient, power) in enumerate(custom_pairs):
        coefficients[6 + 2 * index] = coefficient
        coefficients[7 + 2 * index] = power
    return coefficients


def _signed_expression_terms(expression: str) -> list[str]:
    """Split an expression into additive terms without splitting exponents."""
    compact = expression.replace(" ", "").replace("D", "E").replace("d", "E")
    if not compact:
        return []
    terms: list[str] = []
    start = 0
    depth = 0
    for index, character in enumerate(compact):
        if character == "(":
            depth += 1
        elif character == ")":
            depth = max(0, depth - 1)
        elif (
            character in "+-"
            and index > start
            and depth == 0
            and compact[index - 1] not in "Ee"
            and compact[max(start, index - 2) : index] != "**"
        ):
            terms.append(compact[start:index])
            start = index
    terms.append(compact[start:])
    return [term for term in terms if term]


def _parse_gibbs_term(term: str) -> tuple[float | str, float] | None:
    """Parse one simple Gibbs polynomial term."""
    term = term.rstrip("#")
    term = re.sub(r"(?i)\bLN\(T\)", "LOG(T)", term)
    term = re.sub(r"(?i)\bLOG\(T\)", "LOG(T)", term)
    if term.endswith("LOG(T)") and not term.endswith("T*LOG(T)"):
        coefficient_text = term[: -len("LOG(T)")].rstrip("*")
        try:
            coefficient = _coefficient_from_prefix(coefficient_text)
        except ValueError:
            return None
        return 99.0, coefficient
    if "T" not in term:
        constant = _constant_with_r_value(term)
        return (0.0, constant) if constant is not None else None
    if term.endswith("T*LOG(T)"):
        coefficient_text = term[: -len("T*LOG(T)")].rstrip("*")
        try:
            coefficient = _coefficient_from_prefix(coefficient_text)
        except ValueError:
            return None
        return "TLOGT", coefficient
    if "*T" in term:
        coefficient_text, suffix = term.split("*T", 1)
        try:
            coefficient = _coefficient_from_prefix(coefficient_text)
        except ValueError:
            return None
    elif term.endswith("T") or term.startswith(("+T", "-T", "T")):
        coefficient = -1.0 if term.startswith("-") else 1.0
        suffix = term[2:] if term[:2] in {"+T", "-T"} else term[1:]
    else:
        return None
    if not suffix:
        return 1.0, coefficient
    if suffix.startswith("**"):
        power_text = suffix[2:].strip("()")
        return (float(power_text), coefficient) if _is_float_text(power_text) else None
    return None


def _coefficient_from_prefix(prefix: str) -> float:
    """Return numeric coefficient from an optional multiplication prefix."""
    if prefix in {"", "+"}:
        return 1.0
    if prefix == "-":
        return -1.0
    coefficient = _constant_with_r_value(prefix.rstrip("*"))
    if coefficient is None:
        raise ValueError(f"non-numeric coefficient prefix: {prefix}")
    return coefficient


def _constant_with_r_value(value: str) -> float | None:
    """Return a numeric value for constants that may include TDB's R symbol."""
    text = value.strip().rstrip("*")
    if not text:
        return 1.0
    sign = 1.0
    if text[0] == "+":
        text = text[1:]
    elif text[0] == "-":
        sign = -1.0
        text = text[1:]
    if not text:
        return sign
    product = sign
    for factor in text.split("*"):
        factor = factor.strip()
        if not factor:
            continue
        if factor.upper() == "R":
            product *= 8.31451
        elif _is_float_text(factor):
            product *= float(factor)
        else:
            return None
    return product


def _is_float_text(value: str) -> bool:
    try:
        float(value)
    except ValueError:
        return False
    return True


def _is_log_t_power(power: float) -> bool:
    """Return whether a custom Gibbs term represents standalone ln(T)."""
    return abs(float(power) - 99.0) < 1e-8


def _editable_line(value: Any) -> QLineEdit:
    """Return an editable line edit for database fields."""
    return QLineEdit(str(value))


def _math_label(text: str) -> QLabel:
    """Return a QLabel using the bundled math font."""
    label = QLabel(text)
    label.setTextFormat(Qt.TextFormat.RichText)
    math_family = _load_math_font_family()
    if math_family:
        font = label.font()
        font.setFamily(math_family)
        label.setFont(font)
    return label


def _latex_label(text: str) -> QLabel:
    """Return a QLabel that renders a small HTML/LaTeX-like label."""
    return _math_label(text)


def _reference_tokens_from_source(source: SourceRef) -> str:
    """Extract compact REF tokens from a TDB source command."""
    tokens: list[str] = []
    for match in re.finditer(
        r"\bREF(?:[:#\s]*([A-Za-z0-9_.-]+)|([A-Za-z0-9_.-]*))",
        source.command,
        flags=re.IGNORECASE,
    ):
        token = match.group(1) or match.group(2) or ""
        token = token.strip()
        if token and token.upper() != "REF" and token not in tokens:
            tokens.append(token)
    return " ".join(tokens[:2])


def _source_section(source: SourceRef) -> QFrame:
    rows = [
        ("File", source.file),
        ("Line", source.line),
        ("Column", source.column),
    ]
    return _table_section("Source Mapping", ("Field", "Value"), rows)

__all__ = [name for name in globals() if name.startswith("_") and not name.startswith("__")]
