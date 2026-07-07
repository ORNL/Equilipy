"""Calculation result table and display helpers."""
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
from PySide6.QtGui import (
    QAction,
    QBrush,
    QColor,
    QFont,
    QFontDatabase,
    QFontMetrics,
    QIcon,
    QPainter,
    QPen,
)
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
try:
    from PySide6.QtCharts import QChart, QChartView, QLineSeries, QValueAxis
except ImportError:  # pragma: no cover - depends on optional QtCharts packaging
    QChart = None
    QChartView = None
    QLineSeries = None
    QValueAxis = None

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
_RESULT_COLUMN_GROUP_PREFIX = "@group:"
_THERMO_COEFFICIENT_BOX_WIDTH = 170
_THERMO_POWER_BOX_WIDTH = 80
_THERMO_TEMPERATURE_BOX_WIDTH = round(_THERMO_COEFFICIENT_BOX_WIDTH * 0.5)
_GIBBS_COEFFICIENT_SLOTS = (0, 1, 2, 3, 4, 5, 6, 8, 10, 12)
_FUNCTION_TOKEN_RE = re.compile(r"(?<![A-Za-z0-9_])([A-Za-z_][A-Za-z0-9_]*)")
_PORTABLE_TDB_FUNCTION_NAME_WIDTH = 8
_MAX_TDB_FUNCTION_NAME_WIDTH = 25
_GIBBS_CONTINUITY_TOLERANCE = 1e-5
_RESULT_TABLE_PREVIEW_ROW_LIMIT = 500
_RESULT_TABLE_INLINE_ROW_LIMIT = 12
_MATH_FONT_FAMILY: str | None = None


def _populate_module_results(
    module: CalculationModule,
    table: QTableWidget,
    result: Any,
) -> None:
    if (
        module.kind == "equilibrium"
        and not module.result_columns
        and _populate_default_equilibrium_result_columns(module, table, result)
    ):
        _refresh_module_result_view(module)
        return

    result_table = _result_table_for_module(module, result)
    selected_columns = (
        list(module.result_columns)
        if module.result_columns
        else (
            _default_result_column_keys(module, result_table, result)
            if result_table is not None
            else []
        )
    )
    if selected_columns and _populate_selected_result_columns(
        table,
        result,
        selected_columns,
        module,
    ):
        _refresh_module_result_view(module)
        return

    if module.result_columns:
        module.result_columns = []

    if module.kind == "solidification":
        _populate_scheil_results(table, result)
    elif module.kind == "equilibrium":
        _populate_equilibrium_results(table, result)
    _refresh_module_result_view(module)


def _populate_selected_result_columns(
    table: QTableWidget,
    result: Any,
    selected_columns: list[str],
    module: CalculationModule,
) -> bool:
    result_table = _result_table_for_module(module, result)
    if result_table is None:
        return False
    return _populate_result_table_columns(table, result_table, selected_columns, module)


def _populate_result_table_columns(
    table: QTableWidget,
    result_table: ResultTable,
    selected_columns: list[str],
    module: CalculationModule,
) -> bool:
    available = set(result_table.available_columns())
    expanded_columns = _expand_result_column_defaults(
        selected_columns,
        result_table,
        module,
    )
    columns = [column for column in expanded_columns if column in available]
    if not columns:
        return False

    try:
        selected_table = result_table.select(columns)
    except KeyError:
        return False

    selected_column_metadata = selected_table.columns
    headers = tuple(
        _result_column_display_header(module.kind, column.key, column.label)
        for column in selected_column_metadata
    )
    total_row_count = _result_table_row_count(selected_table)
    preview_row_count = min(
        total_row_count,
        _result_table_display_row_limit(table, total_row_count),
    )
    selected_data = selected_table.to_dict()
    rows = [
        tuple(
            _format_result_column_value(
                module.kind,
                column.key,
                _sequence_value(selected_data.get(column.key), row_index),
            )
            for column in selected_column_metadata
        )
        for row_index in range(preview_row_count)
    ]
    _populate_result_rows(table, headers, rows, total_row_count=total_row_count)
    return True


def _populate_default_equilibrium_result_columns(
    module: CalculationModule,
    table: QTableWidget,
    result: Any,
) -> bool:
    result_table = _default_equilibrium_result_table(module, result)
    if result_table is None:
        return False
    return _populate_result_table_columns(
        table,
        result_table,
        result_table.available_columns(),
        module,
    )


def _default_equilibrium_result_table(
    module: CalculationModule,
    result: Any,
) -> ResultTable | None:
    points = list(getattr(result, "points", []) or [])
    if not points:
        return None

    temperature_unit = _temperature_unit_label(module.temperature_unit)
    amount_is_mass = _is_mass_amount_unit(module.amount_unit)
    amount_attribute = "w_i" if amount_is_mass else "n_i"
    amount_prefix = "w" if amount_is_mass else "n"
    amount_unit = "g" if amount_is_mass else "sp-mol"
    phase_amount_suffix = (
        "_amount_w [g]" if amount_is_mass else "_amount_n [sp-mol]"
    )
    phase_amount_attribute = "amount_w" if amount_is_mass else "amount_n"

    data: dict[str, list[Any]] = {
        f"T [{temperature_unit}]": [
            _temperature_from_kelvin_value(
                getattr(point, "T", np.nan),
                temperature_unit,
            )
            for point in points
        ],
    }

    for component in _ordered_point_amount_keys(points, amount_attribute):
        data[f"{amount_prefix}_{component} [{amount_unit}]"] = [
            _point_amount_value(point, amount_attribute, component) for point in points
        ]

    for phase_name in _stable_phase_names_from_points(points):
        data[f"{phase_name}{phase_amount_suffix}"] = [
            _point_phase_amount_value(point, phase_name, phase_amount_attribute)
            for point in points
        ]

    return ResultTable.from_dict(data)


def _ordered_point_amount_keys(points: list[Any], attribute_name: str) -> list[str]:
    names: list[str] = []
    for point in points:
        mapping = getattr(point, attribute_name, {})
        if not hasattr(mapping, "keys"):
            continue
        for key in mapping.keys():
            text = str(key).strip()
            if text and text not in names:
                names.append(text)
    return names


def _point_amount_value(point: Any, attribute_name: str, key: str) -> Any:
    mapping = getattr(point, attribute_name, {})
    if hasattr(mapping, "get"):
        return mapping.get(key, 0.0)
    return 0.0


def _stable_phase_names_from_points(points: list[Any]) -> list[str]:
    names: list[str] = []
    for point in points:
        stable = getattr(point, "stable_phases", None)
        for name in list(getattr(stable, "names", []) or []):
            text = str(name).strip()
            if text and text.lower() != "nan" and text not in names:
                names.append(text)
    return names


def _point_phase_amount_value(
    point: Any,
    phase_name: str,
    attribute_name: str,
) -> Any:
    phase_map = getattr(point, "phase_map", {})
    phase = phase_map.get(phase_name) if hasattr(phase_map, "get") else None
    if phase is None:
        phases = getattr(point, "phases", None)
        phase = phases.get(phase_name) if hasattr(phases, "get") else None
    if phase is None:
        return 0.0
    return getattr(phase, attribute_name, 0.0)


def _result_table_for_object(result: Any) -> ResultTable | None:
    if result is None:
        return None
    try:
        if isinstance(result, dict) and result:
            return ResultTable.from_dict(result)
        if hasattr(result, "to_table"):
            table = result.to_table()
            return table if isinstance(table, ResultTable) else None
        if hasattr(result, "to_dict"):
            result_data = result.to_dict()
            if isinstance(result_data, dict) and result_data:
                return ResultTable.from_dict(result_data)
    except Exception:
        return None
    return None


def _result_table_for_module(
    module: CalculationModule,
    result: Any,
) -> ResultTable | None:
    base_table = _result_table_for_object(result)
    if base_table is None:
        return None
    base_table = _filter_compound_pseudo_result_columns(base_table, result)
    return _augment_result_table_units(
        base_table,
        temperature_unit=module.temperature_unit,
        pressure_unit=module.pressure_unit,
    )


_COMPOUND_PSEUDO_RESULT_MARKERS = (
    "_endmembers_",
    "_elements_",
    "_partial_gibbs_",
    "_standard_gibbs_energy_",
    "_activity_",
    "_partial_enthalpy_",
    "_partial_entropy_",
    "_partial_heat_capacity_",
)


def _filter_compound_pseudo_result_columns(
    table: ResultTable,
    result: Any,
) -> ResultTable:
    data = table.to_dict()
    filtered = {
        key: value
        for key, value in data.items()
        if not _is_compound_pseudo_result_column(key, result)
    }
    if len(filtered) == len(data):
        return table
    return ResultTable.from_dict(filtered)


def _is_compound_pseudo_result_column(key: str, result: Any) -> bool:
    base_key = key.rsplit(" [", 1)[0]
    compound_phase_names = _compound_phase_names_from_result(result)
    for marker in _COMPOUND_PSEUDO_RESULT_MARKERS:
        if marker not in base_key:
            continue
        phase_name, species_name = base_key.split(marker, 1)
        if phase_name and phase_name in compound_phase_names:
            return True
        if phase_name and species_name and phase_name == species_name:
            return True
    return False


def _compound_phase_names_from_result(result: Any) -> set[str]:
    context = getattr(result, "context", None)
    names = getattr(context, "compound_phase_names", None)
    if names is None:
        return set()
    return {str(name) for name in names}


def _augment_result_table_units(
    table: ResultTable,
    *,
    temperature_unit: str,
    pressure_unit: str,
) -> ResultTable:
    data = table.to_dict()
    augmented = dict(data)

    kelvin_values = _result_temperature_values_in_kelvin(data)
    if kelvin_values is not None:
        for unit in _temperature_output_units(temperature_unit):
            header = f"T [{unit}]"
            augmented[header] = _map_result_column_values(
                kelvin_values,
                lambda value, target_unit=unit: _temperature_from_kelvin_value(
                    value,
                    target_unit,
                ),
            )

    atm_values = _result_pressure_values_in_atm(data)
    if atm_values is not None:
        for unit in _pressure_output_units(pressure_unit):
            header = f"P [{unit}]"
            augmented[header] = _map_result_column_values(
                atm_values,
                lambda value, target_unit=unit: _pressure_from_atm_value(
                    value,
                    target_unit,
                ),
            )

    return ResultTable.from_dict(augmented)


def _temperature_output_units(default_unit: str) -> list[str]:
    return _ordered_with_default(
        ["K", "C", "F", "R"],
        _temperature_unit_label(default_unit),
    )


def _pressure_output_units(default_unit: str) -> list[str]:
    return _ordered_with_default(
        ["atm", "psi", "bar", "Pa", "kPa"],
        _pressure_unit_label(default_unit),
    )


def _ordered_with_default(values: list[str], default: str) -> list[str]:
    ordered = [value for value in values if value == default]
    ordered.extend(value for value in values if value != default)
    return ordered


def _result_temperature_values_in_kelvin(
    data: dict[str, Any],
) -> Any | None:
    for unit in ("K", "C", "F", "R"):
        for prefix in ("T", "Temperature"):
            key = f"{prefix} [{unit}]"
            if key in data:
                return _map_result_column_values(
                    data[key],
                    lambda value, source_unit=unit: _temperature_to_kelvin_value(
                        value,
                        source_unit,
                    ),
                )
    return None


def _result_pressure_values_in_atm(data: dict[str, Any]) -> Any | None:
    for unit in ("atm", "psi", "bar", "Pa", "kPa"):
        key = f"P [{unit}]"
        if key in data:
            return _map_result_column_values(
                data[key],
                lambda value, source_unit=unit: _pressure_to_atm_value(
                    value,
                    source_unit,
                ),
            )
    return None


def _map_result_column_values(value: Any, converter: Any) -> Any:
    if isinstance(value, np.ndarray):
        return [converter(item) for item in value.tolist()]
    if isinstance(value, (list, tuple)):
        return [converter(item) for item in value]
    return converter(value)


def _temperature_to_kelvin_value(value: Any, unit: str) -> Any:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return value
    normalized = _temperature_unit_label(unit)
    if normalized == "C":
        return numeric + 273.15
    if normalized == "F":
        return (numeric - 32.0) * 5.0 / 9.0 + 273.15
    if normalized == "R":
        return numeric * 5.0 / 9.0
    return numeric


def _temperature_from_kelvin_value(value: Any, unit: str) -> Any:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return value
    normalized = _temperature_unit_label(unit)
    if normalized == "C":
        return numeric - 273.15
    if normalized == "F":
        return (numeric - 273.15) * 9.0 / 5.0 + 32.0
    if normalized == "R":
        return numeric * 9.0 / 5.0
    return numeric


def _pressure_to_atm_value(value: Any, unit: str) -> Any:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return value
    normalized = _pressure_unit_label(unit)
    if normalized == "psi":
        return numeric / 14.6959487755
    if normalized == "bar":
        return numeric / 1.01325
    if normalized == "Pa":
        return numeric / 101325.0
    if normalized == "kPa":
        return numeric / 101.325
    return numeric


def _pressure_from_atm_value(value: Any, unit: str) -> Any:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return value
    normalized = _pressure_unit_label(unit)
    if normalized == "psi":
        return numeric * 14.6959487755
    if normalized == "bar":
        return numeric * 1.01325
    if normalized == "Pa":
        return numeric * 101325.0
    if normalized == "kPa":
        return numeric * 101.325
    return numeric


def _pressure_unit_label(unit: str) -> str:
    normalized = unit.strip().lower()
    if normalized == "psi":
        return "psi"
    if normalized == "bar":
        return "bar"
    if normalized == "pa":
        return "Pa"
    if normalized == "kpa":
        return "kPa"
    return "atm"


def _result_column_display_header(
    module_kind: str,
    column: str,
    label: str | None = None,
) -> str:
    if label and label != column:
        return label
    if module_kind == "solidification":
        return _scheil_display_header(column)
    return label or column


def _result_column_group_label(group: str) -> str:
    """Return the display label used for a result-column property set."""
    return group.replace("_", " ").title()


def _result_column_basis(column: str) -> str:
    if (
        column.startswith("n_")
        or "_amount_n" in column
        or "_endmembers_x_" in column
        or "_elements_x_" in column
        or column in {"fl", "fs"}
    ):
        return "mole"
    if (
        column.startswith("w_")
        or "_amount_w" in column
        or "_endmembers_w_" in column
        or "_elements_w_" in column
        or column in {"fl_w", "fs_w"}
    ):
        return "mass"
    return ""


def _result_output_unit_filter_options(
    columns: list[str],
) -> dict[str, list[tuple[str, str]]]:
    options: dict[str, list[tuple[str, str]]] = {
        "Temperature": [],
        "Pressure": [],
        "Amount": [],
    }
    if any(_result_column_basis(column) == "mole" for column in columns):
        options["Amount"].append(("Mole", "basis:mole"))
    if any(_result_column_basis(column) == "mass" for column in columns):
        options["Amount"].append(("Mass", "basis:mass"))
    for unit in ("K", "C", "F", "R"):
        if f"T [{unit}]" in columns:
            options["Temperature"].append((unit, f"temperature:{unit}"))
    for unit in ("atm", "psi", "bar", "Pa", "kPa"):
        if f"P [{unit}]" in columns:
            options["Pressure"].append((unit, f"pressure:{unit}"))
    return {section: values for section, values in options.items() if values}


def _result_unit_filters_from_columns(
    columns: set[str],
    module: CalculationModule,
    available_columns: list[str],
) -> set[str]:
    filters: set[str] = set()
    for column in columns:
        basis = _result_column_basis(column)
        if basis:
            filters.add(f"basis:{basis}")
        temperature_unit = _result_temperature_unit_column(column)
        if temperature_unit:
            filters.add(f"temperature:{temperature_unit}")
        pressure_unit = _result_pressure_unit_column(column)
        if pressure_unit:
            filters.add(f"pressure:{pressure_unit}")

    if not any(item.startswith("basis:") for item in filters):
        default_basis = "mass" if _is_mass_amount_unit(module.amount_unit) else "mole"
        if any(
            _result_column_basis(column) == default_basis
            for column in available_columns
        ):
            filters.add(f"basis:{default_basis}")
    if not any(item.startswith("temperature:") for item in filters):
        default_temperature = _temperature_unit_label(module.temperature_unit)
        if f"T [{default_temperature}]" in available_columns:
            filters.add(f"temperature:{default_temperature}")
    if not any(item.startswith("pressure:") for item in filters):
        default_pressure = _pressure_unit_label(module.pressure_unit)
        if f"P [{default_pressure}]" in available_columns:
            filters.add(f"pressure:{default_pressure}")
    return filters


def _result_column_visible_for_unit_filters(
    column: str,
    filters: set[str],
) -> bool:
    basis = _result_column_basis(column)
    if basis:
        return f"basis:{basis}" in filters
    temperature_unit = _result_temperature_unit_column(column)
    if temperature_unit:
        return f"temperature:{temperature_unit}" in filters
    pressure_unit = _result_pressure_unit_column(column)
    if pressure_unit:
        return f"pressure:{pressure_unit}" in filters
    return True


def _result_column_group_token(group: str) -> str:
    """Return the saved-default token for a result column group."""
    return f"{_RESULT_COLUMN_GROUP_PREFIX}{group}"


def _result_column_default_group(value: str) -> str:
    """Return the group encoded in a saved-default token."""
    if value.startswith(_RESULT_COLUMN_GROUP_PREFIX):
        return value[len(_RESULT_COLUMN_GROUP_PREFIX) :]
    return ""


def _expand_result_column_defaults(
    defaults: list[str],
    table: ResultTable,
    module: CalculationModule,
) -> list[str]:
    """Expand saved group defaults into concrete columns for this result."""
    available_columns = table.available_columns()
    literal_defaults = {
        column for column in defaults if not _result_column_default_group(column)
    }
    filters = _result_unit_filters_from_columns(
        literal_defaults,
        module,
        available_columns,
    )
    expanded: list[str] = []
    seen: set[str] = set()

    def add_column(column: str) -> None:
        if column not in seen and column in available_columns:
            expanded.append(column)
            seen.add(column)

    for item in defaults:
        group = _result_column_default_group(item)
        if not group:
            add_column(item)
            continue
        for column in table.columns:
            if column.group != group:
                continue
            if _result_column_is_stability(column.key):
                continue
            if not _result_column_visible_for_unit_filters(column.key, filters):
                continue
            add_column(column.key)
    return expanded


def _generalized_result_column_defaults(
    selected_columns: list[str],
    table: ResultTable,
    unit_filters: set[str],
) -> list[str]:
    """Store full selected groups as reusable defaults plus literal columns."""
    selected = set(selected_columns)
    visible_columns_by_group: dict[str, list[str]] = {}
    for column in table.columns:
        if _result_column_is_stability(column.key):
            continue
        if not _result_column_visible_for_unit_filters(column.key, unit_filters):
            continue
        visible_columns_by_group.setdefault(column.group, []).append(column.key)

    token_groups = {
        group
        for group, keys in visible_columns_by_group.items()
        if len(keys) > 1 and set(keys).issubset(selected)
    }
    emitted_groups: set[str] = set()
    defaults: list[str] = []

    for column in table.columns:
        if column.key not in selected:
            continue
        if not _result_column_visible_for_unit_filters(column.key, unit_filters):
            continue
        group = column.group
        if group in token_groups and group not in emitted_groups:
            defaults.append(_result_column_group_token(group))
            emitted_groups.add(group)
        defaults.append(column.key)
    return defaults


def _result_temperature_unit_column(column: str) -> str:
    for unit in ("K", "C", "F", "R"):
        if column == f"T [{unit}]":
            return unit
    return ""


def _result_pressure_unit_column(column: str) -> str:
    for unit in ("atm", "psi", "bar", "Pa", "kPa"):
        if column == f"P [{unit}]":
            return unit
    return ""


def _default_result_column_keys(
    module: CalculationModule,
    table: ResultTable,
    result: Any,
) -> list[str]:
    if module.kind == "equilibrium":
        return _default_equilibrium_result_column_keys(module, table, result)
    if module.kind == "solidification":
        return _default_solidification_result_column_keys(module, table, result)

    return table.available_columns()


def _default_solidification_result_column_keys(
    module: CalculationModule,
    table: ResultTable,
    result: Any,
) -> list[str]:
    available = set(table.available_columns())
    default_temperature = f"T [{_temperature_unit_label(module.temperature_unit)}]"
    default_basis = "mass" if _is_mass_amount_unit(module.amount_unit) else "mole"
    default_fraction = "fs_w" if default_basis == "mass" else "fs"
    stable_phases = set(_stable_phase_names_from_result(result))

    selected: list[str] = []

    for column in table.columns:
        key = column.key
        if key == default_temperature:
            selected.append(key)
            continue
        if key in {"label", "Label"}:
            selected.append(key)
            continue
        if key == default_fraction:
            selected.append(key)
            continue
        if _result_column_basis(key) != default_basis:
            continue
        if column.group == "composition":
            selected.append(key)
            continue
        if column.group != "phase_amount":
            continue
        phase_name = _result_column_phase_name(key)
        if not stable_phases or phase_name in stable_phases:
            selected.append(key)

    if default_fraction not in selected:
        for fallback in ("fs", "fs_w"):
            if fallback in available:
                selected.append(fallback)
                break

    if selected:
        return selected
    if default_temperature in available:
        return [default_temperature]
    return []


def _default_equilibrium_result_column_keys(
    module: CalculationModule,
    table: ResultTable,
    result: Any,
) -> list[str]:
    default_temperature = f"T [{_temperature_unit_label(module.temperature_unit)}]"
    default_basis = "mass" if _is_mass_amount_unit(module.amount_unit) else "mole"
    stable_phases = set(_stable_phase_names_from_result(result))
    selected: list[str] = []

    for column in table.columns:
        if column.key == default_temperature:
            selected.append(column.key)
            continue
        if _result_column_basis(column.key) != default_basis:
            continue
        if column.group == "composition":
            selected.append(column.key)
            continue
        if column.group != "phase_amount":
            continue
        phase_name = _result_column_phase_name(column.key)
        if not stable_phases or phase_name in stable_phases:
            selected.append(column.key)

    return selected or table.available_columns()


def _stable_result_column_keys(table: ResultTable, result: Any) -> list[str]:
    available = table.available_columns()
    stable_phases = set(_stable_phase_names_from_result(result))
    if not stable_phases:
        return available
    selected = []
    for column in available:
        if _result_column_is_stability(column):
            continue
        phase_name = _result_column_phase_name(column)
        if phase_name is None or phase_name in stable_phases:
            selected.append(column)
    return selected


def _result_column_phase_name(column: str) -> str | None:
    if column.startswith("stable_phase_"):
        return None
    markers = (
        "_amount_n",
        "_amount_w",
        "_stability",
        "_endmembers_",
        "_elements_",
        "_partial_gibbs_",
        "_standard_gibbs_energy_",
        "_activity_",
        "_partial_enthalpy_",
        "_partial_entropy_",
        "_partial_heat_capacity_",
    )
    for marker in markers:
        if marker in column:
            return column.split(marker, 1)[0]
    return None


def _result_column_is_stability(column: str) -> bool:
    return column.endswith("_stability")


def _stable_phase_names_from_result(result: Any) -> list[str]:
    names: list[str] = []

    def add_name(value: Any) -> None:
        text = str(value).strip()
        if text and text.lower() != "nan" and text not in names:
            names.append(text)

    for point in list(getattr(result, "points", []) or []):
        stable = getattr(point, "stable_phases", None)
        for name in list(getattr(stable, "names", []) or []):
            add_name(name)
        phases = getattr(point, "phases", None)
        if hasattr(phases, "items"):
            for name, phase in phases.items():
                try:
                    stability = float(getattr(phase, "stability", 0.0))
                except (TypeError, ValueError):
                    stability = 0.0
                if np.isclose(stability, 1.0):
                    add_name(name)

    try:
        stable = getattr(result, "stable_phases", None)
    except Exception:
        stable = None
    for name in list(getattr(stable, "names", []) or []):
        add_name(name)

    nested = getattr(result, "equilib_result", None)
    if nested is not None and nested is not result:
        for name in _stable_phase_names_from_result(nested):
            add_name(name)

    return names


def _format_result_column_value(module_kind: str, column: str, value: Any) -> str:
    if module_kind == "solidification":
        return _format_scheil_result_value(column, value)
    if column.startswith("T ["):
        return _format_float_value(value, 2)
    if "_endmembers_" in column or "_elements_" in column:
        return _format_float_value(value, 5)
    return _format_result_value(value)


def _populate_equilibrium_results(table: QTableWidget, result: Any) -> None:
    points = list(getattr(result, "points", []) or [])
    if points:
        if len(points) > 1:
            _populate_batch_equilibrium_point_results(table, points)
            return
        stable_phases = getattr(points[0], "stable_phases", None)
        if hasattr(stable_phases, "values"):
            rows = [
                (phase.name, phase.amount_n, phase.amount_w)
                for phase in stable_phases.values()
                if str(phase.name).strip()
                and str(phase.name).strip().lower() != "nan"
            ]
        else:
            rows = [
                (name, "", "")
                for name in list(getattr(stable_phases, "names", []) or [])
                if str(name).strip() and str(name).strip().lower() != "nan"
            ]
        _populate_result_rows(table, ("Phase", "n(phase)", "w(phase)"), rows)
        return

    _populate_result_rows(table, ("Phase", "n(phase)", "w(phase)"), [])


def _populate_batch_equilibrium_point_results(
    table: QTableWidget,
    points: list[Any],
) -> None:
    rows = []
    preview_points = points[: _result_table_display_row_limit(table, len(points))]
    for index, point in enumerate(preview_points):
        names = [
            str(name).strip()
            for name in getattr(point.stable_phases, "names", [])
            if str(name).strip() and str(name).strip().lower() != "nan"
        ]
        rows.append(
            (
                index + 1,
                _format_result_value(getattr(point, "T", "")),
                _format_result_value(getattr(point, "G", "")),
                " + ".join(names),
            )
        )
    _populate_result_rows(
        table,
        ("#", "T [K]", "G [J]", "Stable phases"),
        rows,
        total_row_count=len(points),
    )


def _as_result_sequence(value: Any) -> list[Any]:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, list):
        return value
    if isinstance(value, tuple):
        return list(value)
    if value is None:
        return []
    return [value]


def _calculation_error_message(result: Any) -> str | None:
    if isinstance(result, dict):
        error = result.get("Error")
        if isinstance(error, list):
            messages = [str(item) for item in error if str(item).strip()]
            return "; ".join(messages) if messages else None
        if error:
            return str(error)
        return None

    gibbs_energy = getattr(result, "G", None)
    if _contains_nan(gibbs_energy):
        return "Calculation failed: equilibrium solver returned NaN Gibbs energy."

    points = list(getattr(result, "points", []) or [])
    if points:
        if any(not getattr(point.stable_phases, "names", []) for point in points):
            return "Calculation failed: no stable phases were found."
        return None

    return "Calculation failed: no stable phases were found."


def _solidification_model_options() -> list[tuple[str, str]]:
    return [
        ("scheil", "Scheil-Gulliver"),
        ("nucleoscheil", "NucleoScheil"),
        ("equilibrium_cooling", "Equilib cooling"),
    ]


def _solidification_model_value(value: str) -> str:
    values = {key for key, _label in _solidification_model_options()}
    text = str(value or "").strip().lower().replace("-", "_").replace(" ", "_")
    aliases = {
        "equilibrium": "equilibrium_cooling",
        "equilib": "equilibrium_cooling",
        "equilib_cooling": "equilibrium_cooling",
        "scheil_gulliver": "scheil",
        "scheil_gulliver_cooling": "scheil",
        "nucleo_scheil": "nucleoscheil",
        "nucleoscheil_cooling": "nucleoscheil",
    }
    normalized = aliases.get(text, text)
    return normalized if normalized in values else "scheil"


def _default_solidification_delta_t(module: CalculationModule) -> str:
    if _solidification_model_value(module.solidification_model) == "nucleoscheil":
        return "0.1"
    return "5"


def _batch_condition_rows(condition: dict[str, Any]) -> list[dict[str, Any]]:
    if not condition:
        return []
    headers = list(condition)
    values_by_header = {
        header: _as_result_sequence(condition[header]) for header in headers
    }
    row_count = max((len(values) for values in values_by_header.values()), default=0)
    rows: list[dict[str, Any]] = []
    for row_index in range(row_count):
        rows.append(
            {
                header: _sequence_value(values, row_index)
                for header, values in values_by_header.items()
            }
        )
    return rows


def _equilib_cooling_batch(
    eq_module: Any,
    liquid_phase_name: str,
    database: Any,
    condition: dict[str, Any],
    *,
    delta_T: float,
    unit: list[str],
    phases: list[str] | None,
    start_from_liquidus: bool,
    progress_callback=None,
) -> dict[str, Any]:
    rows = _batch_condition_rows(condition)
    results: list[dict[str, Any]] = []
    if progress_callback is not None:
        progress_callback(0, len(rows), "Equilibrium Cooling Batch")
    for task_id, row_condition in enumerate(rows, start=1):
        result = eq_module.equilib_cooling(
            liquid_phase_name,
            database,
            row_condition,
            delta_T=delta_T,
            unit=unit,
            phases=phases,
            start_from_liquidus=start_from_liquidus,
        )
        result_data = result.to_dict()
        row_count = max(
            (len(_as_result_sequence(value)) for value in result_data.values()),
            default=0,
        )
        result_data["task_id"] = [task_id] * row_count
        results.append(result_data)
        if progress_callback is not None:
            progress_callback(task_id, len(rows), "Equilibrium Cooling Batch")
    return _combine_result_dicts(results)


def _combine_result_dicts(results: list[dict[str, Any]]) -> dict[str, Any]:
    combined: dict[str, list[Any]] = {}
    for result in results:
        for key, value in result.items():
            combined.setdefault(key, []).extend(_as_result_sequence(value))
    return combined


def _contains_nan(value: Any) -> bool:
    try:
        values = np.asarray(value, dtype=float)
    except (TypeError, ValueError):
        return False
    if values.size == 0:
        return False
    return bool(np.any(np.isnan(values)))


def _populate_scheil_results(table: QTableWidget, result: Any) -> None:
    result_data = result.to_dict() if hasattr(result, "to_dict") else {}
    if not result_data:
        _populate_result_rows(table, ("Temperature", "fs", "Label"), [])
        return

    input_unit = getattr(result, "input_unit", ["K", "atm", "moles"])
    temperature_header = _scheil_temperature_header(result_data, input_unit)
    fs_header = _scheil_fs_header(result_data, input_unit)
    liquid_phase_name = _scheil_liquid_phase_name(result, result_data)
    liquid_endmember_headers = _scheil_liquid_endmember_headers(
        result_data,
        liquid_phase_name,
        input_unit,
    )
    label_header = "label" if "label" in result_data else "Label"
    data_headers = [
        temperature_header,
        fs_header,
        label_header,
        *liquid_endmember_headers,
    ]
    display_headers = [_scheil_display_header(header) for header in data_headers]
    row_count = max(len(result_data.get(header, [])) for header in data_headers)
    preview_row_count = min(row_count, _result_table_display_row_limit(table, row_count))
    rows = []
    for row_index in range(preview_row_count):
        rows.append(
            tuple(
                _format_scheil_result_value(
                    header,
                    _sequence_value(result_data.get(header), row_index)
                )
                for header in data_headers
            )
        )
    _populate_result_rows(
        table,
        tuple(display_headers),
        rows,
        total_row_count=row_count,
    )


def _scheil_temperature_header(
    result_data: dict[str, Any],
    input_unit: list[str],
) -> str:
    unit = input_unit[0] if input_unit else "K"
    header = f"T [{_temperature_unit_label(unit)}]"
    if header in result_data:
        return header
    return "T [K]" if "T [K]" in result_data else header


def _scheil_fs_header(result_data: dict[str, Any], input_unit: list[str]) -> str:
    amount_unit = input_unit[2] if len(input_unit) > 2 else "moles"
    preferred = "fs_w" if _is_mass_amount_unit(amount_unit) else "fs"
    if preferred in result_data:
        return preferred
    fallbacks = ("fs_w", "fs") if _is_mass_amount_unit(amount_unit) else ("fs", "fs_w")
    for fallback in fallbacks:
        if fallback in result_data:
            return fallback
    return preferred


def _scheil_liquid_phase_name(result: Any, result_data: dict[str, Any]) -> str:
    liquid_phase_name = str(getattr(result, "liquid_phase_name", "") or "").strip()
    if liquid_phase_name:
        return liquid_phase_name
    phase_names = {
        key.split("_endmembers_", 1)[0]
        for key in result_data
        if "_endmembers_" in key
    }
    for phase_name in sorted(phase_names):
        if "LIQ" in phase_name.upper():
            return phase_name
    return "LIQUID"


def _scheil_liquid_endmember_headers(
    result_data: dict[str, Any],
    liquid_phase_name: str,
    input_unit: list[str],
) -> list[str]:
    amount_unit = input_unit[2] if len(input_unit) > 2 else "moles"
    preferred = "w" if _is_mass_amount_unit(amount_unit) else "x"
    prefix = f"{liquid_phase_name}_endmembers_{preferred}_"
    headers = sorted(key for key in result_data if key.startswith(prefix))
    if headers:
        return headers

    fallback = "x" if preferred == "w" else "w"
    fallback_prefix = f"{liquid_phase_name}_endmembers_{fallback}_"
    headers = sorted(key for key in result_data if key.startswith(fallback_prefix))
    if headers:
        return headers

    return []


def _scheil_display_header(header: str) -> str:
    if header == "label":
        return "Label"
    if header == "fs":
        return "fs"
    if header == "fs_w":
        return "fs_w"

    marker = "_endmembers_"
    if marker not in header:
        return header

    phase_name, suffix = header.split(marker, 1)
    fraction_type, element = suffix.split("_", 1)
    if fraction_type == "w":
        return f"w({element}@{phase_name})"
    return f"x({element}@{phase_name})"


def _temperature_unit_label(unit: str) -> str:
    normalized_unit = unit.strip().upper()
    if normalized_unit in {"C", "CELSIUS"}:
        return "C"
    if normalized_unit in {"F", "FAHRENHEIT"}:
        return "F"
    if normalized_unit in {"R", "RANKINE"}:
        return "R"
    return "K"


def _sequence_value(values: Any, index: int) -> Any:
    if values is None:
        return ""
    if isinstance(values, np.ndarray):
        return values[index].item() if index < len(values) else ""
    if isinstance(values, (list, tuple)):
        return values[index] if index < len(values) else ""
    return values if index == 0 else ""


def _format_scheil_result_value(header: str, value: Any) -> str:
    if header.startswith("T ["):
        return _format_float_value(value, 2)
    if (
        header in {"fs", "fs_w"}
        or "_endmembers_" in header
    ):
        return _format_float_value(value, 5)
    return _format_result_value(value)


def _format_float_value(value: Any, decimals: int) -> str:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return str(value)
    if np.isnan(numeric):
        return "nan"
    text = f"{numeric:.{decimals}f}".rstrip("0").rstrip(".")
    if text in {"", "-0"}:
        return "0"
    return text


def _format_result_value(value: Any) -> str:
    if isinstance(value, float):
        return f"{value:.8g}"
    return str(value)


def _populate_result_rows(
    table: QTableWidget,
    headers: tuple[str, ...],
    rows: list[tuple[Any, ...]],
    *,
    total_row_count: int | None = None,
) -> None:
    total_rows = len(rows) if total_row_count is None else int(total_row_count)
    display_rows = rows[: _result_table_display_row_limit(table, total_rows)]
    if total_rows > len(display_rows):
        display_rows.append(_result_table_preview_notice_row(headers, total_rows))
    table.setUpdatesEnabled(False)
    try:
        table.clear()
        table.setColumnCount(len(headers))
        table.setHorizontalHeaderLabels(headers)
        table.setRowCount(len(display_rows))
        table.setVisible(bool(display_rows))
        for row_index, row in enumerate(display_rows):
            for column_index, value in enumerate(row):
                table.setItem(row_index, column_index, QTableWidgetItem(str(value)))
        table.horizontalHeader().setStretchLastSection(True)
        table.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        table.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        if len(display_rows) <= _RESULT_TABLE_INLINE_ROW_LIMIT:
            _fit_table_to_rows(table, min_rows=3)
        else:
            _fit_large_result_table_to_rows(table, min_rows=3)
    finally:
        table.setUpdatesEnabled(True)


def _result_table_preview_notice_row(
    headers: tuple[str, ...],
    total_row_count: int,
) -> tuple[str, ...]:
    """Return a table row explaining that only a preview is shown."""
    limit = _result_table_preview_row_limit()
    message = (
        f"Showing first {limit} of {total_row_count} rows. "
        "Export results for the full table."
    )
    return (message, *("" for _ in headers[1:]))


def _result_table_preview_row_limit() -> int:
    """Return the maximum number of rows to render into a Qt preview table."""
    try:
        value = int(os.environ.get("EQUILIPY_RESULT_TABLE_PREVIEW_ROWS", ""))
    except ValueError:
        value = 0
    if value <= 0:
        value = _RESULT_TABLE_PREVIEW_ROW_LIMIT
    return max(1, value)


def _result_table_display_row_limit(
    table: QTableWidget,
    total_row_count: int,
) -> int:
    """Return the row count to render for this result table."""
    if bool(table.property("show_all_results")):
        return max(0, int(total_row_count))
    return _result_table_preview_row_limit()


def _fit_large_result_table_to_rows(table: QTableWidget, min_rows: int = 3) -> None:
    """Give large result tables a bounded viewport with internal scrolling."""
    row_count = max(
        min_rows,
        min(table.rowCount(), _RESULT_TABLE_INLINE_ROW_LIMIT),
    )
    row_height = table.verticalHeader().defaultSectionSize()
    header_height = table.horizontalHeader().height()
    frame = table.frameWidth() * 2
    margin = 14
    table.setFixedHeight(header_height + row_count * row_height + frame + margin)


def _result_table_row_count(table: ResultTable) -> int:
    """Return the row count implied by a column-oriented result table."""
    data = table.to_dict()
    lengths = [
        _result_column_length(data.get(column))
        for column in table.available_columns()
        if _result_column_length(data.get(column)) is not None
    ]
    return max(lengths) if lengths else 1


def _result_column_length(value: Any) -> int | None:
    """Return sequence length for result columns, or None for scalars."""
    if isinstance(value, np.ndarray):
        return len(value)
    if isinstance(value, (list, tuple)):
        return len(value)
    return None


_RESULT_FIGURE_AXIS_TITLE_SIZE = 16
_RESULT_FIGURE_AXIS_LABEL_SIZE = 12
_RESULT_FIGURE_CHART_MIN_HEIGHT = 448
_RESULT_FIGURE_CHART_MAX_HEIGHT = 16777215
_EQUILIBRIUM_FIGURE_CHART_HEIGHT = 448
_RESULT_FIGURE_LEGEND_MAX_COLUMNS = 5
_RESULT_FIGURE_PALETTES: dict[str, list[str]] = {
    "Ocean Sunset": [
        "#001219",
        "#005f73",
        "#0a9396",
        "#94d2bd",
        "#e9d8a6",
        "#ee9b00",
        "#ca6702",
        "#bb3e03",
        "#ae2012",
        "#9b2226",
    ],
    "Electric Rainbow Burst": [
        "#4d86a5",
        "#00bcd4",
        "#2dd4bf",
        "#8be000",
        "#d000ff",
        "#ff8c42",
        "#ffd166",
        "#ef476f",
        "#06d6a0",
        "#118ab2",
    ],
    "Toasty Earth Tones": [
        "#3d2c2e",
        "#7f4f24",
        "#936639",
        "#a68a64",
        "#b6ad90",
        "#c2c5aa",
        "#656d4a",
        "#414833",
    ],
    "Rainbow Fiesta Fun": [
        "#ff006e",
        "#fb5607",
        "#ffbe0b",
        "#8338ec",
        "#3a86ff",
        "#06d6a0",
        "#ef476f",
        "#118ab2",
    ],
    "Whimsical Garden Party": [
        "#386641",
        "#6a994e",
        "#a7c957",
        "#f2e8cf",
        "#bc4749",
        "#7b2cbf",
        "#2a9d8f",
        "#e76f51",
    ],
    "Earthy Tones Adventure": [
        "#264653",
        "#2a9d8f",
        "#8ab17d",
        "#e9c46a",
        "#f4a261",
        "#e76f51",
        "#6d597a",
        "#355070",
    ],
    "Magical Seaside": [
        "#003049",
        "#1d3557",
        "#457b9d",
        "#a8dadc",
        "#f1faee",
        "#ffb703",
        "#fb8500",
        "#c1121f",
    ],
    "Vibrant Tones": [
        "#0077b6",
        "#00b4d8",
        "#90be6d",
        "#f9c74f",
        "#f8961e",
        "#f3722c",
        "#f94144",
        "#577590",
    ],
}
_RESULT_FIGURE_PALETTES["All"] = [
    color
    for palette_name, palette_colors in _RESULT_FIGURE_PALETTES.items()
    if palette_name != "All"
    for color in palette_colors
]
_RESULT_FIGURE_DEFAULT_PALETTE = "Ocean Sunset"


def _new_scheil_figure_chart_view() -> QWidget:
    """Return the shared result figure chart view."""
    if QChart is None or QChartView is None:
        fallback = QLabel("QtCharts is not available.")
        fallback.setObjectName("ResultFigureChart")
        fallback.setMinimumHeight(_RESULT_FIGURE_CHART_MIN_HEIGHT)
        fallback.setMaximumHeight(_RESULT_FIGURE_CHART_MAX_HEIGHT)
        fallback.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )
        return fallback

    chart = QChart()
    chart.setTitle("")
    chart.legend().setVisible(False)
    chart.setBackgroundBrush(QBrush(QColor("#ffffff")))
    chart.setPlotAreaBackgroundBrush(QBrush(QColor("#ffffff")))
    chart.setPlotAreaBackgroundVisible(True)
    chart.setMargins(chart.margins())

    chart_view = QChartView(chart)
    chart_view.setObjectName("ResultFigureChart")
    chart_view.setRenderHint(QPainter.RenderHint.Antialiasing)
    chart_view.setMinimumHeight(_RESULT_FIGURE_CHART_MIN_HEIGHT)
    chart_view.setMaximumHeight(_RESULT_FIGURE_CHART_MAX_HEIGHT)
    chart_view.setSizePolicy(
        QSizePolicy.Policy.Expanding,
        QSizePolicy.Policy.Expanding,
    )
    return chart_view


def _new_result_figure_legend_widget() -> QWidget:
    """Return the wrapping legend used below result figures."""
    return _ResultFigureLegend()


def _refresh_module_result_view(module: CalculationModule) -> None:
    """Sync the table/figure result controls with the current module result."""
    runtime = getattr(module, "runtime", {}) or {}
    result = runtime.get("result") or module.result
    stack = runtime.get("result_view_stack")
    table_button = runtime.get("result_table_button")
    figure_button = runtime.get("result_figure_button")
    y_axis = runtime.get("result_figure_y_axis")
    y_columns_button = runtime.get("result_figure_y_columns_button")
    qdot_field = runtime.get("result_figure_qdot_field")
    qdot_label = runtime.get("result_figure_qdot_label")
    qdot = runtime.get("result_figure_qdot")

    figure_kind = _result_figure_kind(module, result)
    figure_available = bool(figure_kind)
    if isinstance(figure_button, QPushButton):
        figure_button.setEnabled(figure_available)
    if isinstance(stack, QStackedWidget) and not figure_available and stack.currentIndex() == 1:
        stack.setCurrentIndex(0)

    figure_is_active = (
        isinstance(stack, QStackedWidget)
        and stack.currentIndex() == 1
        and figure_available
    )
    if isinstance(table_button, QPushButton):
        table_button.setChecked(not figure_is_active)
    if isinstance(figure_button, QPushButton):
        figure_button.setChecked(figure_is_active)

    if isinstance(y_axis, QComboBox):
        y_axis.setVisible(figure_kind != "equilibrium")
    if isinstance(y_columns_button, QPushButton):
        y_columns_button.setVisible(figure_kind == "equilibrium")
        y_columns_button.setText("select")
    if isinstance(qdot_field, QWidget):
        qdot_field.setVisible(False)
    if isinstance(qdot_label, QLabel):
        qdot_label.setVisible(False)
    if isinstance(qdot, QLineEdit):
        qdot.setVisible(False)

    if not figure_available:
        _clear_result_figure(module)
        return

    _sync_result_figure_controls(module, result, figure_kind)
    if figure_is_active:
        _populate_result_figure(module)


def _result_figure_kind(module: CalculationModule, result: Any) -> str:
    if _solidification_figure_available(module, result):
        return "solidification"
    if _equilibrium_figure_available(module, result):
        return "equilibrium"
    return ""


def _solidification_figure_available(
    module: CalculationModule,
    result: Any,
) -> bool:
    if module.kind != "solidification" or result is None:
        return False
    table = _result_table_for_object(result)
    if table is None or _result_table_row_count(table) <= 1:
        return False
    data = table.to_dict()
    return bool(_solidification_figure_x_axis_options(module, result)) and bool(
        _scheil_temperature_celsius_values(data)
    )


def _equilibrium_figure_available(
    module: CalculationModule,
    result: Any,
) -> bool:
    if module.kind != "equilibrium" or result is None:
        return False
    table = _result_table_for_module(module, result)
    if table is None or _result_table_row_count(table) <= 1:
        return False
    if not _equilibrium_figure_x_axis_pairs(module, result, table):
        return False
    return bool(_equilibrium_figure_y_column_options(module, table))


def _sync_result_figure_controls(
    module: CalculationModule,
    result: Any,
    figure_kind: str,
) -> None:
    runtime = getattr(module, "runtime", {}) or {}
    x_axis = runtime.get("result_figure_x_axis")
    y_axis = runtime.get("result_figure_y_axis")
    palette = runtime.get("result_figure_palette")
    qdot_field = runtime.get("result_figure_qdot_field")
    if isinstance(palette, QComboBox):
        _set_combo_items(
            palette,
            list(_RESULT_FIGURE_PALETTES),
            current=_RESULT_FIGURE_DEFAULT_PALETTE,
        )
    _set_result_figure_chart_size(module, figure_kind)

    if figure_kind == "solidification":
        if isinstance(x_axis, QComboBox):
            _set_combo_items(
                x_axis,
                _solidification_figure_x_axis_options(module, result),
            )
        if isinstance(y_axis, QComboBox):
            _set_combo_items(y_axis, ["Temperature \N{DEGREE SIGN}C"])
            y_axis.setVisible(True)
        y_button = runtime.get("result_figure_y_columns_button")
        if isinstance(y_button, QPushButton):
            y_button.setVisible(False)
        if isinstance(x_axis, QComboBox):
            show_qdot = x_axis.currentText() == "t [s]"
            if isinstance(qdot_field, QWidget):
                qdot_field.setVisible(show_qdot)
            qdot_label = runtime.get("result_figure_qdot_label")
            qdot = runtime.get("result_figure_qdot")
            if isinstance(qdot_label, QLabel):
                qdot_label.setVisible(show_qdot)
            if isinstance(qdot, QLineEdit):
                qdot.setVisible(show_qdot)
        return

    if figure_kind == "equilibrium":
        table = _result_table_for_module(module, result)
        if isinstance(x_axis, QComboBox) and table is not None:
            _set_combo_items(
                x_axis,
                [
                    label
                    for label, _key in _equilibrium_figure_x_axis_pairs(
                        module,
                        result,
                        table,
                    )
                ],
            )
        if isinstance(y_axis, QComboBox):
            y_axis.setVisible(False)
        if isinstance(qdot_field, QWidget):
            qdot_field.setVisible(False)
        if table is not None:
            _sync_equilibrium_figure_y_columns(module, table)


def _set_combo_items(
    combo: QComboBox,
    items: list[str],
    *,
    current: str | None = None,
) -> None:
    old_text = combo.currentText()
    target_text = old_text if old_text in items else current if current in items else ""
    if not target_text and items:
        target_text = items[0]
    existing = [combo.itemText(index) for index in range(combo.count())]
    if existing == items and combo.currentText() == target_text:
        return
    combo.blockSignals(True)
    try:
        combo.clear()
        combo.addItems(items)
        if target_text:
            combo.setCurrentText(target_text)
    finally:
        combo.blockSignals(False)


def _set_result_figure_chart_size(
    module: CalculationModule,
    figure_kind: str,
) -> None:
    chart_view = getattr(module, "runtime", {}).get("result_figure_chart")
    if not isinstance(chart_view, QWidget):
        return
    if figure_kind == "equilibrium":
        chart_view.setMinimumHeight(_EQUILIBRIUM_FIGURE_CHART_HEIGHT)
        chart_view.setMaximumHeight(_EQUILIBRIUM_FIGURE_CHART_HEIGHT)
        chart_view.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Fixed,
        )
        return
    chart_view.setMinimumHeight(_RESULT_FIGURE_CHART_MIN_HEIGHT)
    chart_view.setMaximumHeight(_RESULT_FIGURE_CHART_MAX_HEIGHT)
    chart_view.setSizePolicy(
        QSizePolicy.Policy.Expanding,
        QSizePolicy.Policy.Expanding,
    )


def _populate_result_figure(module: CalculationModule) -> None:
    runtime = getattr(module, "runtime", {}) or {}
    result = runtime.get("result") or module.result
    figure_kind = _result_figure_kind(module, result)
    if figure_kind == "solidification":
        _populate_solidification_result_figure(module, result)
        return
    if figure_kind == "equilibrium":
        _populate_equilibrium_result_figure(module, result)
        return
    _clear_result_figure(module)


def _clear_result_figure(module: CalculationModule) -> None:
    chart_view = getattr(module, "runtime", {}).get("result_figure_chart")
    chart = _result_figure_chart(chart_view)
    if chart is not None:
        chart.removeAllSeries()
        for axis in list(chart.axes()):
            chart.removeAxis(axis)
        chart.setTitle("")
    legend = getattr(module, "runtime", {}).get("result_figure_legend")
    if isinstance(legend, _ResultFigureLegend):
        legend.set_entries([])


def _populate_solidification_result_figure(
    module: CalculationModule,
    result: Any,
) -> None:
    runtime = getattr(module, "runtime", {}) or {}
    x_axis = runtime.get("result_figure_x_axis")
    qdot = runtime.get("result_figure_qdot")
    x_label = x_axis.currentText() if isinstance(x_axis, QComboBox) else ""
    qdot_text = qdot.text() if isinstance(qdot, QLineEdit) else ""
    segments = _scheil_figure_segments(result, x_label, qdot_text)
    if not segments:
        _show_result_figure_message(module, "No solidification figure data.")
        return

    labels = _unique_preserve_order([label for label, _xs, _ys in segments])
    colors = _result_figure_colors(module, len(labels))
    color_by_label = dict(zip(labels, colors, strict=False))
    series_specs = [
        (label, xs, ys, color_by_label.get(label, colors[0]))
        for label, xs, ys in segments
    ]
    legend_entries = [
        (label, color_by_label.get(label, colors[0]))
        for label in labels
    ]
    x_values = [value for _label, xs, _ys, _color in series_specs for value in xs]
    x_bounds: tuple[float, float] | None = None
    if x_label in {"fs [mol/mol]", "fs_w [g/g]"}:
        x_bounds = (0.0, 1.0)
    elif x_label == "t [s]":
        x_bounds = (0.0, max([0.0, *x_values]))
        if x_bounds[1] <= x_bounds[0]:
            x_bounds = (0.0, 1.0)

    _render_result_figure_series(
        module,
        series_specs,
        legend_entries=legend_entries,
        x_title=x_label or "fs",
        y_title="Temperature [\N{DEGREE SIGN}C]",
        x_bounds=x_bounds,
    )


def _populate_equilibrium_result_figure(
    module: CalculationModule,
    result: Any,
) -> None:
    table = _result_table_for_module(module, result)
    if table is None:
        _show_result_figure_message(module, "No equilibrium figure data.")
        return
    runtime = getattr(module, "runtime", {}) or {}
    _sync_equilibrium_figure_y_columns(module, table)
    selected_y = list(runtime.get("result_figure_y_columns") or [])
    if not selected_y:
        _show_result_figure_message(module, "Select at least one y column.")
        return

    x_axis = runtime.get("result_figure_x_axis")
    x_label = x_axis.currentText() if isinstance(x_axis, QComboBox) else ""
    x_key = _equilibrium_figure_x_axis_key(module, result, table, x_label)
    row_count = _result_table_row_count(table)
    data = table.to_dict()
    x_values = _equilibrium_figure_x_values(data, x_key, row_count)
    if not x_values:
        _show_result_figure_message(module, "No x-axis data.")
        return

    colors = _result_figure_colors(module, len(selected_y))
    metadata_by_key = {column.key: column for column in table.columns}
    series_specs: list[tuple[str, list[float], list[float], str]] = []
    for index, key in enumerate(selected_y):
        y_values = _numeric_column_values(data.get(key), row_count)
        label = _result_column_display_header(
            module.kind,
            key,
            metadata_by_key.get(key).label if key in metadata_by_key else key,
        )
        series_specs.append((label, x_values, y_values, colors[index]))

    y_title = _equilibrium_figure_y_axis_title(table, selected_y, series_specs)
    _render_result_figure_series(
        module,
        series_specs,
        legend_entries=[(label, color) for label, _xs, _ys, color in series_specs],
        x_title=x_label or x_key,
        y_title=y_title,
    )


def _equilibrium_figure_y_axis_title(
    table: ResultTable,
    selected_columns: list[str],
    series_specs: list[tuple[str, list[float], list[float], str]],
) -> str:
    """Return a group-aware y-axis title for equilibrium result figures."""
    metadata_by_key = {column.key: column for column in table.columns}
    selected_metadata = [
        metadata_by_key[key] for key in selected_columns if key in metadata_by_key
    ]
    selected_groups = {column.group for column in selected_metadata}
    if len(selected_groups) == 1:
        group_label = _result_column_group_label(selected_metadata[0].group)
        unit_label = _equilibrium_figure_y_axis_basis_label(selected_columns)
        if unit_label:
            return f"{group_label} [{unit_label}]"
        return group_label

    if len(series_specs) == 1:
        return series_specs[0][0]
    return "Value"


def _equilibrium_figure_y_axis_basis_label(selected_columns: list[str]) -> str:
    """Return a basis label only when every selected y column shares one basis."""
    basis_by_column = [_result_column_basis(column) for column in selected_columns]
    basis_values = {basis for basis in basis_by_column if basis}
    if len(basis_values) != 1:
        return ""
    if any(not basis for basis in basis_by_column):
        return ""
    return next(iter(basis_values))


def _render_result_figure_series(
    module: CalculationModule,
    series_specs: list[tuple[str, list[float], list[float], str]],
    *,
    legend_entries: list[tuple[str, str]],
    x_title: str,
    y_title: str,
    x_bounds: tuple[float, float] | None = None,
) -> None:
    chart = _result_figure_chart(
        getattr(module, "runtime", {}).get("result_figure_chart")
    )
    if chart is None or QLineSeries is None or QValueAxis is None:
        _show_result_figure_message(module, "QtCharts is not available.")
        return

    chart.removeAllSeries()
    for axis in list(chart.axes()):
        chart.removeAxis(axis)
    chart.setTitle("")
    chart.legend().setVisible(False)

    plotted_x: list[float] = []
    plotted_y: list[float] = []
    for label, xs, ys, color in series_specs:
        pairs = _finite_xy_pairs(xs, ys)
        if not pairs:
            continue
        series = QLineSeries()
        series.setName(label)
        pen = QPen(QColor(color))
        pen.setWidthF(2.4)
        series.setPen(pen)
        for x_value, y_value in pairs:
            series.append(x_value, y_value)
            plotted_x.append(x_value)
            plotted_y.append(y_value)
        chart.addSeries(series)

    if not chart.series():
        _show_result_figure_message(module, "No plottable figure data.")
        return

    x_axis = QValueAxis()
    y_axis = QValueAxis()
    x_range = _result_figure_axis_range(plotted_x, forced=x_bounds)
    y_range = _result_figure_axis_range(plotted_y)
    _configure_result_figure_axis(x_axis, x_title, value_range=x_range)
    _configure_result_figure_axis(y_axis, y_title, value_range=y_range)
    x_axis.setRange(*x_range)
    y_axis.setRange(*y_range)

    chart.addAxis(x_axis, Qt.AlignmentFlag.AlignBottom)
    chart.addAxis(y_axis, Qt.AlignmentFlag.AlignLeft)
    for series in chart.series():
        series.attachAxis(x_axis)
        series.attachAxis(y_axis)

    legend = getattr(module, "runtime", {}).get("result_figure_legend")
    if isinstance(legend, _ResultFigureLegend):
        legend.set_entries(legend_entries)
    _show_result_figure_message(module, "")


def _result_figure_chart(chart_view: Any) -> Any | None:
    try:
        return chart_view.chart()
    except Exception:
        return None


def _configure_result_figure_axis(
    axis: Any,
    title: str,
    value_range: tuple[float, float] | None = None,
) -> None:
    axis.setTitleText(title)
    title_font = QFont()
    title_font.setPointSizeF(float(_RESULT_FIGURE_AXIS_TITLE_SIZE))
    title_font.setBold(True)
    label_font = QFont()
    label_font.setPointSizeF(float(_RESULT_FIGURE_AXIS_LABEL_SIZE))
    axis.setTitleFont(title_font)
    axis.setLabelsFont(label_font)
    if value_range is None:
        axis.setLabelFormat("%.3g")
    else:
        axis.setLabelFormat(_result_figure_axis_label_format(*value_range))


def _result_figure_axis_label_format(low: float, high: float) -> str:
    """Return a plain printf label format sized to the axis span.

    Wide axes (temperatures, energies) read as integers; narrow axes keep
    just enough decimals to distinguish adjacent ticks.  Scientific notation
    is never used.
    """
    span = abs(high - low)
    if not np.isfinite(span) or span <= 0.0:
        span = max(abs(low), abs(high), 1.0)
    if span >= 20.0:
        return "%.0f"
    decimals = 2
    if span < 1.0:
        decimals = min(6, int(np.ceil(-np.log10(span))) + 2)
    return f"%.{decimals}f"


def _result_figure_axis_range(
    values: list[float],
    *,
    forced: tuple[float, float] | None = None,
) -> tuple[float, float]:
    if forced is not None:
        low, high = forced
        if high <= low:
            high = low + 1.0
        return low, high
    finite_values = [value for value in values if np.isfinite(value)]
    if not finite_values:
        return 0.0, 1.0
    low = min(finite_values)
    high = max(finite_values)
    if np.isclose(low, high):
        pad = abs(low) * 0.05 or 1.0
        return low - pad, high + pad
    pad = (high - low) * 0.04
    return low - pad, high + pad


def _show_result_figure_message(module: CalculationModule, message: str) -> None:
    status = getattr(module, "runtime", {}).get("result_figure_status")
    if isinstance(status, QLabel):
        status.setText(message)
        status.setVisible(bool(message))


def _solidification_figure_x_axis_options(
    module: CalculationModule,
    result: Any,
) -> list[str]:
    table = _result_table_for_object(result)
    if table is None:
        return []
    available = set(table.available_columns())
    preferred = ["fs_w", "fs"]
    labels = {"fs_w": "fs_w [g/g]", "fs": "fs [mol/mol]"}
    options = [
        labels[column]
        for column in preferred
        if column in available
    ]
    if "Q [J]" in available:
        options.append("t [s]")
    return options


def _scheil_figure_segments(
    result: Any,
    x_axis: str,
    qdot_text: str,
) -> list[tuple[str, list[float], list[float]]]:
    table = _result_table_for_object(result)
    if table is None:
        return []
    data = table.to_dict()
    x_values = _scheil_figure_x_values(data, x_axis, qdot_text)
    y_values = _scheil_temperature_celsius_values(data)
    labels = _string_column_values(data.get("label"), min(len(x_values), len(y_values)))
    if not labels:
        labels = ["Result" for _ in range(min(len(x_values), len(y_values)))]
    count = min(len(x_values), len(y_values), len(labels))
    if count <= 0:
        return []

    segments: list[tuple[str, list[float], list[float]]] = []
    current_label = labels[0]
    current_x = [x_values[0]]
    current_y = [y_values[0]]
    for index in range(1, count):
        label = labels[index]
        if label != current_label:
            segments.append((current_label, current_x, current_y))
            current_label = label
            current_x = [x_values[index - 1], x_values[index]]
            current_y = [y_values[index - 1], y_values[index]]
            continue
        current_x.append(x_values[index])
        current_y.append(y_values[index])
    segments.append((current_label, current_x, current_y))
    return segments


def _scheil_figure_x_values(
    data: dict[str, Any],
    x_axis: str,
    qdot_text: str,
) -> list[float]:
    if x_axis == "t [s]":
        q_values = _numeric_column_values(data.get("Q [J]"), _column_value_count(data))
        qdot = _positive_float(qdot_text, default=1.0)
        if not q_values:
            return []
        q0 = q_values[0]
        return [abs(value - q0) / qdot for value in q_values]
    if x_axis == "fs [mol/mol]":
        return _numeric_column_values(data.get("fs"), _column_value_count(data))
    return _numeric_column_values(data.get("fs_w"), _column_value_count(data))


def _scheil_temperature_celsius_values(data: dict[str, Any]) -> list[float]:
    for unit in ("C", "K", "F", "R"):
        key = f"T [{unit}]"
        if key not in data:
            continue
        values = _numeric_column_values(data.get(key), _column_value_count(data))
        if unit == "C":
            return values
        return [
            float(_temperature_from_kelvin_value(
                _temperature_to_kelvin_value(value, unit),
                "C",
            ))
            for value in values
        ]
    return []


def _equilibrium_figure_x_axis_pairs(
    module: CalculationModule,
    result: Any,
    table: ResultTable,
) -> list[tuple[str, str]]:
    available = set(table.available_columns())
    pairs: list[tuple[str, str]] = []
    preferred_temperature = f"T [{_temperature_unit_label(module.temperature_unit)}]"
    if preferred_temperature in available:
        pairs.append((preferred_temperature, preferred_temperature))
    else:
        for column in table.columns:
            if column.key.startswith("T ["):
                pairs.append((column.key, column.key))
                break

    for column in table.columns:
        if column.key.startswith("w_") and column.key.endswith(" [g]"):
            label = column.label or column.key
            if label != column.key and not label.endswith(" [g]"):
                label = f"{label} [g]"
            pairs.append((label, column.key))

    if _result_table_row_count(table) > 1:
        pairs.append(("task_id", "task_id"))
    return pairs


def _equilibrium_figure_x_axis_key(
    module: CalculationModule,
    result: Any,
    table: ResultTable,
    label: str,
) -> str:
    for option_label, key in _equilibrium_figure_x_axis_pairs(module, result, table):
        if option_label == label:
            return key
    pairs = _equilibrium_figure_x_axis_pairs(module, result, table)
    return pairs[0][1] if pairs else ""


def _equilibrium_figure_x_values(
    data: dict[str, Any],
    key: str,
    row_count: int,
) -> list[float]:
    if key == "task_id":
        return [float(index + 1) for index in range(row_count)]
    return _numeric_column_values(data.get(key), row_count)


def _equilibrium_figure_y_column_options(
    module: CalculationModule,
    table: ResultTable,
) -> list[str]:
    data = table.to_dict()
    row_count = _result_table_row_count(table)
    candidates: list[str] = []
    for column in table.columns:
        key = column.key
        if column.group == "stable_phase_summary" or key.startswith("stable_phase_"):
            continue
        if column.group == "condition" or key == "task_id":
            continue
        if _result_column_is_stability(key):
            continue
        if _result_column_has_numeric_values(data.get(key), row_count):
            candidates.append(key)
    return candidates


def _sync_equilibrium_figure_y_columns(
    module: CalculationModule,
    table: ResultTable,
) -> None:
    runtime = getattr(module, "runtime", {}) or {}
    candidates = _equilibrium_figure_y_column_options(module, table)
    candidate_set = set(candidates)
    selected = [
        key
        for key in list(runtime.get("result_figure_y_columns") or [])
        if key in candidate_set
    ]
    if not selected:
        preferred = ["G [J]", "H [J]", "S [J/K]", "Cp [J/K]"]
        selected = [key for key in preferred if key in candidate_set]
    if not selected and candidates:
        selected = [candidates[0]]
    runtime["result_figure_y_columns"] = selected
    button = runtime.get("result_figure_y_columns_button")
    if isinstance(button, QPushButton):
        button.setText("select")
        button.setToolTip(
            f"{len(selected)} selected"
            if selected
            else "Select equilibrium figure columns"
        )


def _open_result_figure_y_columns_dialog(module: CalculationModule) -> None:
    """Open a Result Columns style selector for equilibrium figure lines."""
    result = getattr(module, "runtime", {}).get("result") or module.result
    table = _result_table_for_module(module, result)
    if table is None:
        QMessageBox.information(
            None,
            "Result Columns",
            "Run or load a calculation result before selecting columns.",
        )
        return

    candidates = _equilibrium_figure_y_column_options(module, table)
    if not candidates:
        QMessageBox.information(
            None,
            "Result Columns",
            "No plottable equilibrium columns are available.",
        )
        return

    _sync_equilibrium_figure_y_columns(module, table)
    runtime = getattr(module, "runtime", {}) or {}
    parent = _result_figure_dialog_parent(module)
    dialog = QDialog(parent)
    dialog.setWindowTitle(f"{module.name} Result Columns")
    dialog.setModal(True)
    dialog.resize(560, 560)

    layout = QVBoxLayout(dialog)
    layout.setContentsMargins(18, 16, 18, 16)
    layout.setSpacing(12)

    candidate_set = set(candidates)
    current_columns = set(runtime.get("result_figure_y_columns") or [])
    default_columns = set(runtime.get("result_figure_y_columns_default") or [])
    if not default_columns:
        default_columns = set(current_columns)
    unit_filters = _result_unit_filters_from_columns(
        current_columns or set(candidates),
        module,
        candidates,
    )
    checkboxes: list[tuple[str, QCheckBox]] = []
    sync_callbacks: list[Any] = []
    syncing_section_checks = False

    scroll = QScrollArea(dialog)
    scroll.setWidgetResizable(True)
    scroll.setFrameShape(QFrame.Shape.NoFrame)
    checkbox_host = QWidget()
    checkbox_layout = QVBoxLayout(checkbox_host)
    checkbox_layout.setContentsMargins(0, 0, 0, 0)
    checkbox_layout.setSpacing(6)
    scroll.setWidget(checkbox_host)
    layout.addWidget(scroll)

    def set_checkboxes(
        target_checkboxes: list[QCheckBox],
        checked: bool,
    ) -> None:
        for checkbox in target_checkboxes:
            checkbox.setChecked(checked)

    def set_section_state(
        checkbox: QCheckBox,
        target_checkboxes: list[QCheckBox],
    ) -> None:
        nonlocal syncing_section_checks
        checked_count = sum(1 for target in target_checkboxes if target.isChecked())
        syncing_section_checks = True
        checkbox.setTristate(True)
        if checked_count == 0:
            checkbox.setCheckState(Qt.CheckState.Unchecked)
        elif checked_count == len(target_checkboxes):
            checkbox.setCheckState(Qt.CheckState.Checked)
        else:
            checkbox.setCheckState(Qt.CheckState.PartiallyChecked)
        syncing_section_checks = False

    def sync_all_headers() -> None:
        for sync in sync_callbacks:
            sync()

    def clear_checkbox_layout() -> None:
        while checkbox_layout.count():
            item = checkbox_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()

    def visible_column_keys() -> set[str]:
        return {
            key
            for key in candidate_set
            if _result_column_visible_for_unit_filters(key, unit_filters)
        }

    def rebuild_column_groups() -> None:
        nonlocal checkboxes, sync_callbacks, current_columns
        if checkboxes:
            current_columns = {
                key for key, checkbox in checkboxes if checkbox.isChecked()
            }
        visible_keys = visible_column_keys()
        current_columns &= visible_keys
        checkboxes = []
        sync_callbacks = []
        clear_checkbox_layout()

        columns_by_group: dict[str, list[Any]] = {}
        for column in table.columns:
            if column.key not in visible_keys:
                continue
            if column.group == "stable_phase_summary":
                continue
            columns_by_group.setdefault(column.group, []).append(column)

        for group, columns in columns_by_group.items():
            group_frame = QFrame()
            group_frame.setObjectName("ResultColumnGroup")
            group_layout = QVBoxLayout(group_frame)
            group_layout.setContentsMargins(0, 0, 0, 0)
            group_layout.setSpacing(6)

            group_header = QHBoxLayout()
            group_header.setSpacing(6)
            group_toggle = QToolButton()
            group_toggle.setObjectName("ResultColumnGroupToggle")
            group_toggle.setText(group.replace("_", " ").title())
            group_toggle.setToolButtonStyle(
                Qt.ToolButtonStyle.ToolButtonTextBesideIcon
            )
            group_toggle.setArrowType(Qt.ArrowType.DownArrow)
            group_toggle.setCheckable(True)
            group_toggle.setChecked(True)
            group_header.addWidget(group_toggle)

            section_box = QCheckBox("")
            section_box.setObjectName("ResultSectionCheck")
            section_box.setToolTip("Select all visible columns in this section")
            group_header.addWidget(section_box)
            group_header.addStretch(1)

            content = QWidget()
            content_layout = QVBoxLayout(content)
            content_layout.setContentsMargins(20, 0, 0, 0)
            content_layout.setSpacing(4)

            group_checkboxes: list[QCheckBox] = []
            for column in columns:
                checkbox = QCheckBox(column.label or column.key)
                checkbox.setObjectName("ResultColumnCheck")
                checkbox.setToolTip(column.key)
                checkbox.setChecked(column.key in current_columns)
                content_layout.addWidget(checkbox)
                checkboxes.append((column.key, checkbox))
                group_checkboxes.append(checkbox)

            def sync_header_checks(
                all_checkbox: QCheckBox = section_box,
                all_targets: list[QCheckBox] = group_checkboxes,
            ) -> None:
                set_section_state(all_checkbox, all_targets)

            def on_section_state_changed(
                state: int,
                targets: list[QCheckBox] = group_checkboxes,
                sync: Any = sync_header_checks,
            ) -> None:
                if syncing_section_checks:
                    return
                set_checkboxes(targets, state != 0)
                sync()

            section_box.stateChanged.connect(on_section_state_changed)
            for checkbox in group_checkboxes:
                checkbox.stateChanged.connect(
                    lambda _state, sync=sync_header_checks: sync()
                )
            sync_callbacks.append(sync_header_checks)
            group_toggle.toggled.connect(
                lambda checked, widget=content, toggle=group_toggle: (
                    widget.setVisible(checked),
                    toggle.setArrowType(
                        Qt.ArrowType.DownArrow
                        if checked
                        else Qt.ArrowType.RightArrow
                    ),
                )
            )
            sync_header_checks()
            group_layout.addLayout(group_header)
            group_layout.addWidget(content)
            checkbox_layout.addWidget(group_frame)

        checkbox_layout.addStretch(1)
        for _key, checkbox in checkboxes:
            checkbox.stateChanged.connect(lambda _state: sync_all_headers())
        sync_all_headers()

    def set_selected_columns(keys: set[str]) -> None:
        for key, checkbox in checkboxes:
            checkbox.setChecked(key in keys)
        sync_all_headers()

    bulk_actions = QHBoxLayout()
    output_unit_button = QPushButton("Output Unit")
    output_unit_button.setObjectName("SmallActionButton")
    select_button = QToolButton()
    select_button.setText("Select")
    select_button.setObjectName("SmallMenuButton")
    select_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
    select_menu = QMenu(select_button)
    select_all_action = QAction("Select All", select_menu)
    select_stable_action = QAction("Select Stable", select_menu)
    select_menu.addAction(select_all_action)
    select_menu.addAction(select_stable_action)
    select_button.setMenu(select_menu)
    clear_all = QPushButton("Clear")
    clear_all.setObjectName("SmallActionButton")
    save_default = QPushButton("Save as Default")
    save_default.setObjectName("SmallActionButton")
    reset_default = QPushButton("Default")
    reset_default.setObjectName("SmallActionButton")

    select_all_action.triggered.connect(
        lambda: set_selected_columns({key for key, _checkbox in checkboxes})
    )
    select_stable_action.triggered.connect(
        lambda: set_selected_columns(
            set(_stable_result_column_keys(table, result)) & visible_column_keys()
        )
    )
    clear_all.clicked.connect(lambda: set_selected_columns(set()))
    output_unit_button.clicked.connect(
        lambda: _open_result_output_unit_dialog(
            dialog,
            candidates,
            unit_filters,
            lambda filters: (
                unit_filters.clear(),
                unit_filters.update(filters),
                rebuild_column_groups(),
            ),
        )
    )
    bulk_actions.addWidget(output_unit_button)
    bulk_actions.addWidget(select_button)
    bulk_actions.addWidget(clear_all)
    bulk_actions.addStretch(1)
    bulk_actions.addWidget(save_default)
    bulk_actions.addWidget(reset_default)
    layout.addLayout(bulk_actions)
    rebuild_column_groups()

    actions = QHBoxLayout()
    actions.addStretch(1)
    apply = QPushButton("Apply")
    apply.setObjectName("PrimaryButton")
    cancel = QPushButton("Cancel")
    cancel.setObjectName("SmallActionButton")
    actions.addWidget(apply)
    actions.addWidget(cancel)
    layout.addLayout(actions)

    def apply_selection() -> None:
        selected = [key for key, checkbox in checkboxes if checkbox.isChecked()]
        if not selected:
            QMessageBox.warning(
                dialog,
                "Result Columns",
                "Select at least one result column.",
            )
            return
        runtime["result_figure_y_columns"] = selected
        _sync_equilibrium_figure_y_columns(module, table)
        _populate_result_figure(module)
        dialog.accept()

    def save_selection_as_default() -> None:
        selected = [key for key, checkbox in checkboxes if checkbox.isChecked()]
        if not selected:
            QMessageBox.warning(
                dialog,
                "Result Columns",
                "Select at least one result column before saving a default.",
            )
            return
        runtime["result_figure_y_columns_default"] = selected

    apply.clicked.connect(apply_selection)
    save_default.clicked.connect(save_selection_as_default)
    reset_default.clicked.connect(lambda: set_selected_columns(default_columns))
    cancel.clicked.connect(dialog.reject)
    dialog.exec()


def _open_result_output_unit_dialog(
    parent: QWidget,
    available_columns: list[str],
    current_filters: set[str],
    apply_callback: Callable[[set[str]], None],
) -> None:
    options = _result_output_unit_filter_options(available_columns)
    if not options:
        QMessageBox.information(
            parent,
            "Output Unit",
            "No alternate output units are available for this result.",
        )
        return

    dialog = QDialog(parent)
    dialog.setWindowTitle("Output Unit")
    dialog.setModal(True)
    dialog.resize(420, 360)

    layout = QVBoxLayout(dialog)
    layout.setContentsMargins(18, 16, 18, 16)
    layout.setSpacing(14)
    boxes: list[tuple[str, QCheckBox]] = []

    for section, section_options in options.items():
        if not section_options:
            continue
        section_label = QLabel(section)
        section_label.setObjectName("SmallSectionLabel")
        layout.addWidget(section_label)
        for label, filter_key in section_options:
            checkbox = QCheckBox(label)
            checkbox.setObjectName("ResultUnitCheck")
            checkbox.setChecked(filter_key in current_filters)
            layout.addWidget(checkbox)
            boxes.append((filter_key, checkbox))

    layout.addStretch(1)
    actions = QHBoxLayout()
    actions.addStretch(1)
    apply = QPushButton("Apply")
    apply.setObjectName("PrimaryButton")
    cancel = QPushButton("Cancel")
    cancel.setObjectName("SmallActionButton")
    actions.addWidget(apply)
    actions.addWidget(cancel)
    layout.addLayout(actions)

    def apply_units() -> None:
        selected = {key for key, checkbox in boxes if checkbox.isChecked()}
        apply_callback(selected)
        dialog.accept()

    apply.clicked.connect(apply_units)
    cancel.clicked.connect(dialog.reject)
    dialog.exec()


def _result_figure_dialog_parent(module: CalculationModule) -> QWidget | None:
    for key in (
        "result_figure_y_columns_button",
        "result_figure_chart",
        "results_table",
    ):
        widget = getattr(module, "runtime", {}).get(key)
        if isinstance(widget, QWidget):
            return widget.window()
    return None


def _result_figure_colors(module: CalculationModule, count: int) -> list[str]:
    runtime = getattr(module, "runtime", {}) or {}
    palette_widget = runtime.get("result_figure_palette")
    palette_name = (
        palette_widget.currentText()
        if isinstance(palette_widget, QComboBox) and palette_widget.currentText()
        else _RESULT_FIGURE_DEFAULT_PALETTE
    )
    palette = _RESULT_FIGURE_PALETTES.get(
        palette_name,
        _RESULT_FIGURE_PALETTES[_RESULT_FIGURE_DEFAULT_PALETTE],
    )
    if count <= 0:
        return []
    return [palette[index % len(palette)] for index in range(count)]


def _numeric_column_values(value: Any, row_count: int) -> list[float]:
    if value is None:
        return []
    if isinstance(value, np.ndarray):
        values = value.tolist()
    elif isinstance(value, (list, tuple)):
        values = list(value)
    else:
        values = [value for _ in range(max(1, row_count))]
    numeric: list[float] = []
    for item in values[: max(0, row_count)]:
        try:
            number = float(item)
        except (TypeError, ValueError):
            continue
        if np.isfinite(number):
            numeric.append(number)
    return numeric


def _string_column_values(value: Any, row_count: int) -> list[str]:
    if value is None:
        return []
    if isinstance(value, np.ndarray):
        values = value.tolist()
    elif isinstance(value, (list, tuple)):
        values = list(value)
    else:
        values = [value for _ in range(max(1, row_count))]
    return [str(item) for item in values[: max(0, row_count)]]


def _result_column_has_numeric_values(value: Any, row_count: int) -> bool:
    return bool(_numeric_column_values(value, row_count))


def _column_value_count(data: dict[str, Any]) -> int:
    lengths = [
        _result_column_length(value)
        for value in data.values()
        if _result_column_length(value) is not None
    ]
    return max(lengths) if lengths else 1


def _finite_xy_pairs(xs: list[float], ys: list[float]) -> list[tuple[float, float]]:
    pairs: list[tuple[float, float]] = []
    for x_value, y_value in zip(xs, ys, strict=False):
        if np.isfinite(x_value) and np.isfinite(y_value):
            pairs.append((x_value, y_value))
    return pairs


def _positive_float(text: str, *, default: float) -> float:
    try:
        value = float(text)
    except (TypeError, ValueError):
        return default
    return value if value > 0.0 else default


def _unique_preserve_order(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    unique: list[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        unique.append(value)
    return unique


def _wrapped_legend_label(label: str, *, max_line_length: int) -> str:
    """Wrap long legend labels at phase separators."""
    text = str(label)
    if len(text) <= max_line_length:
        return text
    parts = text.split("+")
    if len(parts) <= 1:
        return text
    lines: list[str] = []
    current = ""
    for part in parts:
        candidate = part if not current else f"{current}+{part}"
        if len(candidate) <= max_line_length or not current:
            current = candidate
            continue
        lines.append(current)
        current = part
    if current:
        lines.append(current)
    return "\n".join(lines)


class _ResultFigureLegend(QWidget):
    """Compact wrapping legend for result figures."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._entries: list[tuple[str, str]] = []
        self.setObjectName("ResultFigureLegend")
        self.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Preferred,
        )
        self.setVisible(False)

    def set_entries(self, entries: list[tuple[str, str]]) -> None:
        self._entries = list(entries)
        self.setVisible(bool(self._entries))
        self.updateGeometry()
        self.update()

    def hasHeightForWidth(self) -> bool:
        return True

    def heightForWidth(self, width: int) -> int:
        if not self._entries:
            return 0
        columns = self._legend_column_count(width)
        rows = self._legend_rows(width, columns)
        return sum(row_height for _row_entries, row_height in rows) + 8

    def sizeHint(self) -> QSize:
        return QSize(420, self.heightForWidth(max(420, self.width())))

    def _legend_column_count(self, width: int) -> int:
        if not self._entries:
            return 0
        if width >= 650:
            return min(_RESULT_FIGURE_LEGEND_MAX_COLUMNS, len(self._entries))
        if width >= 520:
            return min(4, len(self._entries))
        if width >= 380:
            return min(3, len(self._entries))
        if width >= 240:
            return min(2, len(self._entries))
        return 1

    def _legend_row_count(self, width: int) -> int:
        columns = self._legend_column_count(width)
        return len(self._legend_rows(width, columns))

    def _legend_label_wrap_length(self, columns: int) -> int:
        if columns >= 5:
            return 18
        if columns >= 4:
            return 22
        if columns >= 3:
            return 26
        if columns >= 2:
            return 34
        return 60

    def _legend_rows(
        self,
        width: int,
        columns: int,
    ) -> list[tuple[list[tuple[str, str]], int]]:
        if columns <= 0:
            return []
        metrics = QFontMetrics(self.font())
        wrap_length = self._legend_label_wrap_length(columns)
        rows: list[tuple[list[tuple[str, str]], int]] = []
        for start in range(0, len(self._entries), columns):
            row_entries = self._entries[start : start + columns]
            max_lines = max(
                1,
                *[
                    _wrapped_legend_label(label, max_line_length=wrap_length).count("\n")
                    + 1
                    for label, _color in row_entries
                ],
            )
            row_height = max(22, metrics.lineSpacing() * max_lines + 8)
            rows.append((row_entries, row_height))
        return rows

    def paintEvent(self, event: Any) -> None:
        super().paintEvent(event)
        if not self._entries:
            return
        painter = QPainter(self)
        try:
            painter.setRenderHint(QPainter.RenderHint.Antialiasing)
            width = max(1, self.width())
            columns = self._legend_column_count(width)
            column_width = max(1, width // max(1, columns))
            wrap_length = self._legend_label_wrap_length(columns)
            y = 4
            for row_entries, row_height in self._legend_rows(width, columns):
                for column_index, (label, color) in enumerate(row_entries):
                    x = column_index * column_width + 8
                    square_size = 8
                    square_y = y + max(0, (row_height - square_size) // 2)
                    painter.setPen(QPen(QColor("#1a1f20")))
                    painter.setBrush(QBrush(QColor(color)))
                    painter.drawRect(x, square_y, square_size, square_size)
                    painter.setPen(QColor("#eef2f3"))
                    wrapped = _wrapped_legend_label(
                        label,
                        max_line_length=wrap_length,
                    )
                    text_rect = self.rect().adjusted(
                        x + square_size + 6,
                        y,
                        -(width - (column_index + 1) * column_width + 4),
                        0,
                    )
                    text_rect.setHeight(row_height)
                    painter.drawText(
                        text_rect,
                        int(
                            Qt.AlignmentFlag.AlignLeft
                            | Qt.AlignmentFlag.AlignVCenter
                        ),
                        wrapped,
                    )
                y += row_height
        finally:
            painter.end()


def _is_solution_phase(phase: Phase) -> bool:
    """Return whether a DatabaseIR phase should use the solution editor path."""
    return phase.model.upper() in {
        "CEF",
        "IDMX",
        "QKTO",
        "SUBL",
        "RKMP",
        "RKMPM",
        "SUBLM",
        "SUBG",
        "SUBQ",
        "SUBI",
        "SUBM",
        "SOLUTION",
    }

__all__ = [name for name in globals() if name.startswith("_") and not name.startswith("__")]
