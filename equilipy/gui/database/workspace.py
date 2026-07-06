"""Database workspace behavior: tree actions, validation, import, and export."""

from __future__ import annotations

# ruff: noqa: F401,F403,F405,I001,E501

import copy
import json
import os
import re
import sys
from itertools import product
from pathlib import Path
from typing import Any

import numpy as np
from PySide6.QtCore import (
    QItemSelectionModel,
    QModelIndex,
    QPersistentModelIndex,
    Qt,
)
from PySide6.QtGui import QAction, QBrush, QColor
from PySide6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QHeaderView,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMenu,
    QMessageBox,
    QProgressDialog,
    QPushButton,
    QRadioButton,
    QStyle,
    QTableView,
    QTableWidget,
    QTableWidgetItem,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from equilipy.database_ir import (
    ConstituentSet,
    DatabaseIR,
    Diagnostic,
    Element,
    FunctionDefinition,
    GibbsRange,
    Parameter,
    Phase,
    SourceRef,
    Species,
    TdbFunctionTerm,
    load_database_ir,
    remove_redundant_disordered_phase_aliases,
    split_database,
    write_eqdb,
    write_tdb,
)
from equilipy.database_ir.tdb_writer import _canonical_tdb_gibbs_expression
from equilipy.gui.widgets.periodic_table import PeriodicTableSelector
from equilipy.utils import HSCp2G

from equilipy.gui.compat_exports import *
from equilipy.gui.database.dragdrop import DatabaseDragDropCopyMixin
from equilipy.gui.models import (
    CompoundPhasePayload,
    CompoundRecord,
    DatabaseTreeModel,
    DiagnosticsModel,
    SolutionGroupPayload,
    SolutionParameterPayload,
    SolutionSpeciesPayload,
    SolutionSublatticePayload,
    ThermoRangePayload,
    thermo_range_warning_key,
)
from equilipy.gui.widgets import forms as _forms

_collapsible_form_section = _forms.collapsible_form_section
_fit_table_columns_to_headers = _forms.fit_table_columns_to_headers
_fit_table_to_rows = _forms.fit_table_to_rows
_fit_tree_to_items = _forms.fit_tree_to_items
_form_section = _forms.form_section
_form_section_with_actions = _forms.form_section_with_actions
_inline_table = _forms.inline_table
_read_only_line = _forms.read_only_line
_read_only_text = _forms.read_only_text
_table_section = _forms.table_section
_text_section = _forms.text_section
_untitled_form_section = _forms.untitled_form_section


def _main_window_attr(name: str, fallback: Any) -> Any:
    """Return compatibility-patched attributes from ``main_window`` when present."""
    module = sys.modules.get("equilipy.gui.main_window")
    return getattr(module, name, fallback) if module is not None else fallback


def _runtime_function_thermo_data(
    function: FunctionDefinition,
    database: DatabaseIR | None = None,
    seen: set[str] | None = None,
) -> tuple[float, float, list[list[float]], list[list[float]]]:
    thermo_data = _main_window_attr("_function_thermo_data", _function_thermo_data)
    return thermo_data(function, database, seen)


def _runtime_progress_dialog(*args: Any, **kwargs: Any) -> Any:
    progress_dialog = _main_window_attr("QProgressDialog", QProgressDialog)
    return progress_dialog(*args, **kwargs)


def _strip_pertinent_function_markers(expression: str) -> str:
    """Normalize GUI-only `#Function` markers before parser validation."""
    text = " ".join(str(expression or "").strip().split())
    return re.sub(r"#\s*([A-Za-z_][A-Za-z0-9_]*)", r"\1", text)


def _wrapped_single_range_expression(expression: str) -> str:
    """Wrap one expression line in a 298.15-6000 K TDB-like range."""
    return f"298.15 {expression}; 6000 N"


def _missing_pertinent_function_names(
    expression: str,
    database: DatabaseIR,
) -> list[str]:
    """Return referenced names that do not exist in the database Functions."""
    available = {function.name.upper() for function in database.functions}
    ignored = {"T", "LN", "LOG", "EXP", "R"}
    missing: list[str] = []
    seen: set[str] = set()
    for match in re.finditer(
        r"(?<![A-Za-z0-9_])([A-Za-z_][A-Za-z0-9_]*)",
        expression,
    ):
        name = match.group(1)
        normalized = name.upper()
        if normalized in ignored or normalized in seen:
            continue
        seen.add(normalized)
        if normalized not in available:
            missing.append(name)
    return missing


def _restore_function_name_case(expression: str, database: DatabaseIR) -> str:
    """Restore DatabaseIR Function name casing after strict math normalization."""
    restored = str(expression or "")
    for function in database.functions:
        restored = re.sub(
            rf"(?<![A-Za-z0-9_]){re.escape(function.name.upper())}(?![A-Za-z0-9_])",
            function.name,
            restored,
            flags=re.IGNORECASE,
        )
    return restored


def _solution_parameters_for_phase(
    database: DatabaseIR,
    phase: Phase,
) -> list[Parameter]:
    """Return parameters attached to one solution phase, preserving order."""
    phase_name = phase.name.upper()
    return [
        parameter
        for parameter in database.parameters
        if parameter.phase_name.upper() == phase_name
    ]


def _solution_parameter_is_interaction(parameter: Parameter) -> bool:
    """Return whether a solution parameter is a mixing or non-endmember term."""
    if str(parameter.parameter_type).strip().upper() != "G":
        return True
    return any("," in section for section in _solution_parameter_sections(parameter))


def _solution_parameter_sections(parameter: Parameter) -> list[str]:
    """Return target sections split by sublattice."""
    target = list(parameter.target)
    if len(target) == 1 and ":" in target[0]:
        return target[0].split(":")
    return target


def _solution_parameter_target_text(parameter: Parameter) -> str:
    """Return a compact target string for solution tables."""
    return ":".join(_solution_parameter_sections(parameter)) or "-"


def _solution_parameter_edit_text(parameter: Parameter) -> str:
    """Return the editable expression body for one solution parameter."""
    source_expression = str(parameter.metadata.get("source_expression", "")).strip()
    if source_expression:
        return source_expression
    expression = str(parameter.expression or "").strip()
    return _single_tdb_range_body(expression) or expression


def _single_tdb_range_body(expression: str) -> str:
    """Return the body from a single ``Tlow expr; Thigh status`` range."""
    text = " ".join(str(expression or "").split())
    if text.count(";") != 1:
        return ""
    number = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
    match = re.fullmatch(
        rf"{number}\s+(?P<body>.+?)\s*;\s*{number}\s+[A-Za-z]",
        text,
        flags=re.IGNORECASE,
    )
    if match is None:
        return ""
    return match.group("body").strip()


def _solution_parameter_command_with_expression(
    parameter: Parameter,
    expression: str,
) -> str:
    """Return a solution PAR command with an updated expression range."""
    command = str(parameter.source.command or "").strip()
    match = _parameter_command_prefix_match(command)
    if match is not None:
        return f"{command[: match.end()]}{expression} !"
    parameter_type = str(parameter.parameter_type or "G").strip() or "G"
    target = _solution_parameter_target_text(parameter)
    return (
        f"PAR {parameter_type}({parameter.phase_name},{target};"
        f"{parameter.order}) {expression} !"
    )


def _solution_sublattice_label(
    index: int,
    species_count: int,
    phase: Phase | None = None,
) -> str:
    """Return FactSage-like sublattice labels."""
    if phase is not None:
        return f"{_solution_sublattice_name(phase, index)} ({species_count})"
    if 1 <= index <= 26:
        return f"{chr(ord('A') + index - 1)} ({species_count})"
    return f"Sublattice {index} ({species_count})"


def _solution_sublattice_rows(phase: Phase) -> list[tuple[str, str, str]]:
    """Return rows describing solution sublattices."""
    rows: list[tuple[str, str, str]] = []
    for index, constituent_set in enumerate(phase.constituents, start=1):
        rows.append(
            (
                _solution_sublattice_label(index, len(constituent_set.species), phase),
                _compact_gui_number(constituent_set.site_ratio),
                ", ".join(constituent_set.species),
            )
        )
    return rows


def _solution_expected_endmembers(phase: Phase) -> list[tuple[str, ...]]:
    """Return cartesian-product endmembers implied by solution constituents."""
    species_sets = [
        tuple(species for species in constituent_set.species if str(species).strip())
        for constituent_set in phase.constituents
    ]
    if not species_sets or any(not species for species in species_sets):
        return []
    return [tuple(str(part) for part in target) for target in product(*species_sets)]


def _solution_species_is_vacancy(species_name: str) -> bool:
    """Return whether a solution constituent is a vacancy-like placeholder."""
    return str(species_name).strip().upper() in {"VA", "VACUUM", "/-", "ELECTRON_GAS"}


def _solution_species_composition(
    database: DatabaseIR,
    species_name: str,
) -> dict[str, float]:
    """Return elemental composition for one solution species."""
    canonical_name = _canonical_solution_species_name(species_name)
    if _solution_species_is_vacancy(canonical_name):
        return {}
    species = _solution_species_record(database, canonical_name)
    if species is not None and species.composition:
        return dict(species.composition)
    return _species_composition_guess(canonical_name)


def _solution_endmember_stoichiometry_key(
    database: DatabaseIR,
    phase: Phase,
    target: tuple[str, ...],
) -> tuple[tuple[str, float], ...]:
    """Return a canonical stoichiometry key for matching endmember targets."""
    composition: dict[str, float] = {}
    for index, species_name in enumerate(target):
        if _solution_species_is_vacancy(species_name):
            continue
        site_ratio = (
            phase.constituents[index].site_ratio
            if 0 <= index < len(phase.constituents)
            else 1.0
        )
        for symbol, amount in _solution_species_composition(database, species_name).items():
            canonical_symbol = _canonical_solution_species_name(symbol)
            composition[canonical_symbol] = (
                composition.get(canonical_symbol, 0.0)
                + float(amount) * float(site_ratio)
            )
    return tuple(
        (symbol, round(amount, 12))
        for symbol, amount in sorted(composition.items())
        if abs(amount) > 1e-12
    )


def _solution_endmember_match_key(
    database: DatabaseIR,
    phase: Phase,
    target: tuple[str, ...],
) -> tuple[tuple[str, ...], tuple[tuple[str, float], ...]]:
    """Return a stable key that preserves ordered sublattice identity."""
    canonical_target = tuple(
        _canonical_solution_species_name(part)
        for part in target
        if str(part).strip()
    )
    non_vacancy_target = tuple(
        part for part in canonical_target if not _solution_species_is_vacancy(part)
    )
    return (
        non_vacancy_target,
        _solution_endmember_stoichiometry_key(database, phase, canonical_target),
    )


def _solution_endmember_parameter(
    parameter: Parameter,
) -> bool:
    """Return whether a solution parameter is a concrete endmember."""
    return (
        str(parameter.parameter_type).strip().upper() == "G"
        and bool(parameter.target)
        and not _solution_parameter_is_interaction(parameter)
    )


def _solution_endmember_source_range(parameter: Parameter) -> str:
    """Return the TDB range used when writing an endmember PAR command."""
    expression = str(parameter.expression or "").strip()
    if _single_tdb_range_body(expression):
        return expression
    body = _solution_parameter_edit_text(parameter).strip() or "ZERO#"
    return f"298.15 {body}; 6000 N"


def _solution_endmember_owner_name(phase: Phase, parameter: Parameter) -> str:
    """Return the tree owner name for a solution endmember parameter."""
    return f"{phase.name}:{_solution_parameter_target_text(parameter)}"


def _solution_endmember_parameter_for_range_payload(
    database: DatabaseIR,
    payload: ThermoRangePayload,
) -> Parameter | None:
    """Return the endmember parameter represented by one Cp range payload."""
    if ":" not in payload.owner_name:
        return None
    phase_name, target_text = payload.owner_name.split(":", 1)
    return next(
        (
            parameter
            for parameter in database.parameters
            if parameter.phase_name.upper() == phase_name.upper()
            and _solution_endmember_parameter(parameter)
            and _solution_parameter_target_text(parameter).upper()
            == target_text.upper()
        ),
        None,
    )


def _solution_metadata_value(phase: Phase, *keys: str) -> str:
    """Return the first non-empty phase metadata value from candidate keys."""
    for key in keys:
        value = str(phase.metadata.get(key, "")).strip()
        if value:
            return value
    return ""


_SOLUTION_MODEL_CHOICES = ("Bragg-Williams", "CEF")
_MAX_SOLUTION_SUBLATTICES = 5
_CHARGED_SPECIES_RE = re.compile(r"^\s*([^\[\],]+?)\s*(?:\[\s*([^\]]+)\s*\])?\s*$")


def _solution_model_editor_value(phase: Phase) -> str:
    """Return the editable two-choice solution model name."""
    if str(phase.model or "").strip().upper() == "CEF":
        return "CEF"
    display_name = _phase_model_display_name(phase)
    if display_name in _SOLUTION_MODEL_CHOICES:
        return display_name
    return "CEF"


def _default_sublattice_name(index: int) -> str:
    """Return the default A-E style sublattice name."""
    if 1 <= index <= 26:
        return chr(ord("A") + index - 1)
    return f"S{index}"


def _solution_sublattice_name(phase: Phase, index: int) -> str:
    """Return the stored editable sublattice name for a 1-based index."""
    labels = phase.metadata.get("sublattice_labels")
    if isinstance(labels, list) and 0 <= index - 1 < len(labels):
        label = str(labels[index - 1]).strip()
        if label:
            return label
    return _default_sublattice_name(index)


def _set_solution_sublattice_name(phase: Phase, index: int, name: str) -> None:
    """Store a user-edited sublattice name in phase metadata."""
    labels = phase.metadata.get("sublattice_labels")
    if not isinstance(labels, list):
        labels = []
    while len(labels) < len(phase.constituents):
        labels.append(_default_sublattice_name(len(labels) + 1))
    labels[index - 1] = str(name or _default_sublattice_name(index)).strip()
    phase.metadata["sublattice_labels"] = labels[: len(phase.constituents)]


def _trim_solution_sublattice_labels(phase: Phase) -> None:
    """Keep stored sublattice labels aligned with the current sublattice count."""
    labels = phase.metadata.get("sublattice_labels")
    if isinstance(labels, list):
        phase.metadata["sublattice_labels"] = labels[: len(phase.constituents)]


def _canonical_solution_species_name(name: str) -> str:
    """Normalize common species casing without changing custom names."""
    text = str(name or "").strip()
    if text.upper() == "VA":
        return "Va"
    if text in {"/-", "ELECTRON_GAS"}:
        return text
    if re.fullmatch(r"[A-Za-z]{1,2}", text):
        return text[0].upper() + text[1:].lower()
    return text


def _parse_solution_charge(text: str) -> float | None:
    """Parse charge notations like 3+, +3, 2-, or -2."""
    token = str(text or "").strip().replace(" ", "")
    if not token:
        return None
    sign = 1.0
    if token[-1:] in {"+", "-"}:
        sign = 1.0 if token[-1] == "+" else -1.0
        magnitude = token[:-1]
    elif token[:1] in {"+", "-"}:
        sign = 1.0 if token[0] == "+" else -1.0
        magnitude = token[1:]
    else:
        magnitude = token
    if not magnitude:
        return None
    try:
        return sign * float(magnitude)
    except ValueError:
        return None


def _parse_solution_constituent_token(token: str) -> tuple[str, float | None] | None:
    """Parse one comma-separated constituent token with optional charge."""
    match = _CHARGED_SPECIES_RE.match(str(token or ""))
    if match is None:
        return None
    name = _canonical_solution_species_name(match.group(1))
    if not name:
        return None
    charge = _parse_solution_charge(match.group(2) or "")
    return name, charge


def _species_composition_guess(species_name: str) -> dict[str, float]:
    """Guess a simple elemental composition for new GUI-created species."""
    name = _canonical_solution_species_name(species_name)
    if name in {"Va", "/-"}:
        return {}
    if re.fullmatch(r"[A-Z][a-z]?", name):
        return {name: 1.0}
    return {name: 1.0}


def _solution_charge_text(charge: float) -> str:
    """Return a compact bracket charge suffix."""
    if abs(charge) < 1e-12:
        return ""
    magnitude = _compact_gui_number(abs(charge))
    sign = "+" if charge > 0 else "-"
    return f" [{magnitude}{sign}]"


def _solution_species_record(
    database: DatabaseIR,
    species_name: str,
) -> Species | None:
    """Return a species record by case-insensitive name."""
    target = species_name.upper()
    return next(
        (species for species in database.species if species.name.upper() == target),
        None,
    )


def _upsert_solution_species(
    database: DatabaseIR,
    species_name: str,
    charge: float | None,
) -> None:
    """Create or update a species record from a solution constituent editor."""
    species = _solution_species_record(database, species_name)
    if species is None:
        database.species.append(
            Species(
                species_name,
                _species_composition_guess(species_name),
                0.0 if charge is None else charge,
            )
        )
        return
    if charge is not None:
        species.charge = charge


def _solution_constituent_display(
    database: DatabaseIR,
    species_name: str,
) -> str:
    """Return a constituent display name including charge when known."""
    canonical_name = _canonical_solution_species_name(species_name)
    species = _solution_species_record(database, canonical_name)
    if species is None:
        return canonical_name
    return f"{canonical_name}{_solution_charge_text(species.charge)}"


def _solution_constituents_display(
    database: DatabaseIR,
    constituent_set: ConstituentSet,
) -> str:
    """Return comma-separated constituents for the editable solution table."""
    return ", ".join(
        _solution_constituent_display(database, species_name)
        for species_name in constituent_set.species
    )


def _solution_species_formula(database: DatabaseIR, species_name: str) -> str:
    """Return a readable formula for one species."""
    species = _solution_species_record(database, species_name)
    if species is None or not species.composition:
        return species_name
    parts: list[str] = []
    for symbol, amount in species.composition.items():
        suffix = "" if amount == 1 else _compact_gui_number(amount)
        parts.append(f"{symbol}{suffix}")
    if species.charge:
        charge = _compact_gui_number(abs(species.charge))
        sign = "+" if species.charge > 0 else "-"
        parts.append(f"[{charge}{sign}]")
    return "".join(parts)


def _compact_gui_number(value: Any) -> str:
    """Render numeric values compactly for read-only GUI fields."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if number.is_integer():
        return str(int(number))
    return f"{number:g}"


def _compound_phase_primary_g_parameter(
    database: DatabaseIR,
    phase: Phase,
) -> Parameter | None:
    """Return the editable primary Gibbs parameter for a compound phase."""
    parameters = _compound_phase_g_parameters(database, phase)
    return parameters[0] if parameters else None


def _parameter_command_with_expression(
    parameter: Parameter,
    expression: str,
) -> str:
    """Return a PARAMETER command with an updated expression body."""
    command = str(parameter.source.command or "").strip()
    match = _parameter_command_prefix_match(command)
    if match is not None:
        return f"{command[: match.end()]}{expression} !"
    target = ":".join(parameter.target) if parameter.target else parameter.phase_name
    return f"PARAMETER G({parameter.phase_name},{target};0) {expression} !"


class _FunctionDropLineEdit(QLineEdit):
    """Line edit that accepts Function drags as plain expression text."""

    def __init__(self, value: str = ""):
        super().__init__(value)
        self.setAcceptDrops(True)

    def dragEnterEvent(self, event) -> None:
        if event.mimeData().hasText():
            event.acceptProposedAction()
            return
        super().dragEnterEvent(event)

    def dragMoveEvent(self, event) -> None:
        if event.mimeData().hasText():
            event.acceptProposedAction()
            return
        super().dragMoveEvent(event)

    def dropEvent(self, event) -> None:
        text = event.mimeData().text().strip()
        if text:
            self.insert(text)
            event.acceptProposedAction()
            return
        super().dropEvent(event)


_EXPORT_PSEUDO_ELEMENTS = {"VA", "/-"}


def _ask_export_element_subset(
    parent,
    database: DatabaseIR,
) -> list[str] | None:
    """Element-subset dialog for database export.

    Lists the real elements alphabetically with one checkbox each, all
    checked by default.  ``VA`` and the electron pseudo-element are always
    kept and are not listed.  Returns the checked symbols in alphabetical
    order, or ``None`` when the dialog is cancelled.  Checking every element
    means "export unchanged".
    """
    symbols = sorted(
        (
            element.symbol
            for element in database.elements
            if element.symbol.upper() not in _EXPORT_PSEUDO_ELEMENTS
        ),
        key=str.upper,
    )
    if not symbols:
        return []

    dialog = QDialog(parent)
    dialog.setWindowTitle("Export Elements")
    layout = QVBoxLayout(dialog)
    layout.addWidget(
        QLabel(
            "Select the elements to export.  Unchecked elements are removed\n"
            "by an exact subsystem split (VA is kept automatically).",
            dialog,
        )
    )

    # Periodic-table selector: database elements are clickable cells with an
    # accent edge box when selected; elements outside the database render 70%
    # transparent and disabled.
    selector = PeriodicTableSelector(symbols, dialog)
    layout.addWidget(selector)
    checkboxes = selector.buttons

    select_row = QHBoxLayout()
    select_all_button = QPushButton("Select All", dialog)
    deselect_all_button = QPushButton("Deselect All", dialog)
    select_row.addWidget(select_all_button)
    select_row.addWidget(deselect_all_button)
    select_row.addStretch(1)
    layout.addLayout(select_row)

    buttons = QDialogButtonBox(
        QDialogButtonBox.StandardButton.Ok
        | QDialogButtonBox.StandardButton.Cancel,
        dialog,
    )
    buttons.accepted.connect(dialog.accept)
    buttons.rejected.connect(dialog.reject)
    layout.addWidget(buttons)
    ok_button = buttons.button(QDialogButtonBox.StandardButton.Ok)

    def _refresh_ok_state() -> None:
        ok_button.setEnabled(
            any(checkbox.isChecked() for checkbox in checkboxes.values())
        )

    for checkbox in checkboxes.values():
        checkbox.toggled.connect(_refresh_ok_state)

    def _set_all(checked: bool) -> None:
        for checkbox in checkboxes.values():
            checkbox.setChecked(checked)

    select_all_button.clicked.connect(lambda: _set_all(True))
    deselect_all_button.clicked.connect(lambda: _set_all(False))

    # Exposed for tests.
    dialog._element_checkboxes = checkboxes
    dialog._export_ok_button = ok_button

    if dialog.exec() != QDialog.DialogCode.Accepted:
        return None
    return [symbol for symbol in symbols if checkboxes[symbol].isChecked()]


class DatabaseWorkspaceLogicMixin(DatabaseDragDropCopyMixin):
    """Database tree, rendering, validation, and file I/O controller methods."""

    def _show_database_summary(self, database: DatabaseIR | None = None) -> None:
        database = database or self.database
        self.current_database = database
        if not _database_ir_has_content(database) and not database.name.strip():
            self._set_editor_heading("", "")
            self.editor_badge.setText("")
            self._clear_form()
            self._set_source_preview(None)
            return

        element_names = _database_ir_element_names(database)
        subtitle = (
            f"Available elements: {', '.join(element_names)}"
            if element_names
            else ""
        )
        self._set_editor_heading("Database Overview", subtitle)
        self._clear_form()

        solution_rows = _database_ir_phase_overview_rows(
            [phase for phase in database.phases if _is_solution_phase(phase)]
        )
        compound_rows = _database_ir_phase_overview_rows(
            [phase for phase in database.phases if not _is_solution_phase(phase)]
        )
        self.editor_form_layout.insertWidget(
            0,
            _table_section(
                f"Solution phases ({len(solution_rows)})",
                ("#", "Phase", "Solution model", "Structure", "Species (count)"),
                solution_rows or [("", "", "", "", "")],
            ),
        )
        self.editor_form_layout.insertWidget(
            1,
            _table_section(
                f"Stoichiometric compounds ({len(compound_rows)})",
                ("#", "Phase", "Model", "Structure", "Species (count)"),
                compound_rows or [("", "", "", "", "")],
            ),
        )
        self._set_source_preview(None)

    def _show_selected(self, current, _previous) -> None:
        self._commit_active_database_thermo_editor()
        if current.isValid():
            self.last_tree_index = QPersistentModelIndex(current)
        payload = self.tree_model.payload(current)
        display_name = self.tree_model.label(current)
        self.current_database = self.tree_model.database_for(current)
        self._update_database_selection_status(payload, display_name)
        self._render_payload(payload, display_name)

    def _commit_active_database_thermo_editor(self) -> None:
        """Persist dirty Cp/G/H/S editor values before replacing the editor."""
        editor = getattr(self, "_active_database_thermo_editor", None)
        if editor is None:
            return
        try:
            editor._capture_current_mode()
        except RuntimeError:
            return

    def _update_database_selection_status(
        self,
        payload: Any,
        display_name: str,
    ) -> None:
        """Reflect the selected database object in the bottom status bar."""
        selected_function = self._function_for_payload(payload)
        if selected_function is not None:
            self._set_function_status_bar_message(selected_function)
            return
        if isinstance(payload, ThermoRangePayload) and payload.owner_kind == "Compound":
            self._set_compound_range_status_bar_message(payload, display_name)
            return
        if isinstance(payload, CompoundPhasePayload):
            warnings = self._compound_phase_status_warnings(payload)
            if warnings:
                self._set_status_bar_message(
                    f"Compound {payload.compound_name} phase "
                    f"{payload.phase_label}: {warnings[0]}",
                    "warning",
                )
                return
        if isinstance(payload, CompoundRecord):
            warnings = self._compound_record_status_warnings(payload)
            if warnings:
                self._set_status_bar_message(
                    f"Compound {payload.name}: {warnings[0]}",
                    "warning",
                )
                return

        if isinstance(payload, DatabaseIR):
            name = payload.name or "Database"
            self._set_status_bar_message(f"Database {name} selected.", "ready")
            return
        if isinstance(payload, list):
            self._set_status_bar_message(
                f"{display_name or 'Database group'} selected.",
                "ready",
            )
            return
        if display_name:
            self._set_status_bar_message(f"{display_name} selected.", "ready")
            return
        self._set_status_bar_message("Ready.", "ready")

    def _function_for_payload(self, payload: Any) -> FunctionDefinition | None:
        """Return the function represented by a selected tree payload, if any."""
        if isinstance(payload, FunctionDefinition):
            return payload
        if isinstance(payload, ThermoRangePayload) and payload.owner_kind == "Function":
            return _function_by_name(self.current_database, payload.owner_name)
        return None

    def _current_selection_function(self) -> FunctionDefinition | None:
        """Return the currently selected function or function range owner."""
        if self.tree is None or self.tree_model is None:
            return None
        current = self.tree.currentIndex()
        if not current.isValid():
            return None
        return self._function_for_payload(self.tree_model.payload(current))

    def _set_function_status_bar_message(self, function: FunctionDefinition) -> None:
        """Show strict-name errors or Gibbs warnings for a selected function."""
        errors = self.function_name_errors.get(
            id(function),
            _function_name_error_messages(function),
        )
        if errors:
            self._set_status_bar_message(
                f"Function {function.name}: {errors[0]}",
                "error",
            )
            return

        name_warnings = self.function_name_warnings.get(
            id(function),
            _function_name_warning_messages(function),
        )
        if name_warnings:
            self._set_status_bar_message(
                f"Function {function.name}: {name_warnings[0]}",
                "warning",
            )
            return

        warnings = self.function_continuity_warnings.get(id(function), [])
        if warnings:
            self._set_status_bar_message(
                f"Function {function.name}: {warnings[0]}",
                "warning",
            )
            return

        self._set_status_bar_message(f"Function {function.name}: OK.", "ready")

    def _set_compound_range_status_bar_message(
        self,
        payload: ThermoRangePayload,
        display_name: str,
    ) -> None:
        """Show compound range Gibbs warnings in the status bar."""
        key = thermo_range_warning_key(payload)
        warnings = self.compound_range_warnings.get(key or "", [])
        label = display_name or f"{payload.property_name} {payload.index}"
        if warnings:
            self._set_status_bar_message(f"{label}: {warnings[0]}", "warning")
            return
        self._set_status_bar_message(f"{label} selected.", "ready")

    def _compound_phase_status_warnings(
        self,
        payload: CompoundPhasePayload,
    ) -> list[str]:
        """Return warnings for a compound phase and its Cp ranges."""
        prefix = f"Compound|{payload.compound_name}:{payload.phase_label}|"
        return _unique_messages(
            warning
            for key, warnings in self.compound_range_warnings.items()
            if key.startswith(prefix)
            for warning in warnings
        )

    def _compound_record_status_warnings(
        self,
        payload: CompoundRecord,
    ) -> list[str]:
        """Return warnings for all phase states under a compound."""
        prefix = f"Compound|{payload.name}:"
        return _unique_messages(
            warning
            for key, warnings in self.compound_range_warnings.items()
            if key.startswith(prefix)
            for warning in warnings
        )

    def _show_database_tree_menu(self, position) -> None:
        index = self.tree.indexAt(position)
        if not index.isValid():
            return
        self.tree.setCurrentIndex(index)
        payload = self.tree_model.payload(index)
        label = self.tree_model.label(index)
        database = self.tree_model.database_for(index)
        menu = QMenu(self)
        menu_handlers: dict[str, Any] = {}

        if isinstance(payload, DatabaseIR):
            edit_metadata = QAction("Edit Metadata", self)
            edit_metadata.setData("edit_database_metadata")
            menu.addAction(edit_metadata)
            menu_handlers["edit_database_metadata"] = lambda target_database=payload: (
                self._edit_database_metadata(target_database)
            )
            menu.addSeparator()

            resolve_all = QAction("Resolve All", self)
            resolve_all.triggered.connect(
                lambda _checked=False, target_database=payload: (
                    self._resolve_database_gibbs_mismatches(target_database)
                )
            )
            menu.addAction(resolve_all)
            menu.addSeparator()

            save_and_close = QAction("Save and Close", self)
            save_and_close.triggered.connect(
                lambda _checked=False, target_database=payload: (
                    self._save_and_close_database(target_database)
                )
            )
            menu.addAction(save_and_close)

            close_database = QAction(
                self._theme_icon(
                    "window-close",
                    QStyle.StandardPixmap.SP_DialogCloseButton,
                ),
                "Close",
                self,
            )
            close_database.triggered.connect(
                lambda _checked=False, target_database=payload: (
                    self._close_database(target_database)
                )
            )
            menu.addAction(close_database)
        elif isinstance(payload, list) and label == "Functions":
            add_function = QAction("Add New Function", self)
            add_function.triggered.connect(
                lambda _checked=False, target_database=database: (
                    self._add_database_function(target_database)
                )
            )
            menu.addAction(add_function)
        elif isinstance(payload, list) and label == "Compounds":
            add_compound = QAction("Add New Compound", self)
            add_compound.triggered.connect(
                lambda _checked=False, target_database=database: (
                    self._add_database_compound(target_database)
                )
            )
            menu.addAction(add_compound)
        elif isinstance(payload, list) and label == "Solutions":
            add_solution = QAction("Add New Solution", self)
            add_solution.triggered.connect(
                lambda _checked=False, target_database=database: (
                    self._add_database_solution(target_database)
                )
            )
            menu.addAction(add_solution)
        elif isinstance(payload, FunctionDefinition):
            add_cp = QAction("Add New Cp", self)
            add_cp.setEnabled(len(payload.temperature_ranges) < 6)
            add_cp.triggered.connect(
                lambda _checked=False, function=payload, target_database=database: (
                    self._add_function_cp_range(target_database, function)
                )
            )
            menu.addAction(add_cp)
            menu.addSeparator()

            duplicate_function = QAction("Duplicate Function", self)
            duplicate_function.triggered.connect(
                lambda _checked=False, function=payload, target_database=database: (
                    self._duplicate_database_function(target_database, function)
                )
            )
            menu.addAction(duplicate_function)

            rename_function = QAction("Rename Function", self)
            rename_function.triggered.connect(
                lambda _checked=False, function=payload, target_database=database: (
                    self._begin_rename_database_function(target_database, function)
                )
            )
            menu.addAction(rename_function)

            validate_function = QAction("Validate Function", self)
            validate_function.triggered.connect(
                lambda _checked=False, function=payload, target_database=database: (
                    self._validate_database_function(target_database, function)
                )
            )
            menu.addAction(validate_function)

            copy_expression = QAction("Copy Expression", self)
            copy_expression.triggered.connect(
                lambda _checked=False, function=payload: (
                    self._copy_database_function_expression(function)
                )
            )
            menu.addAction(copy_expression)
            menu.addSeparator()

            delete_function = QAction(
                self._theme_icon("edit-delete", QStyle.StandardPixmap.SP_TrashIcon),
                "Delete Function",
                self,
            )

            def delete_selected_function(
                _checked: bool = False,
                target_database: DatabaseIR = database,
                function: FunctionDefinition = payload,
            ) -> None:
                self._delete_database_function(target_database, function)

            delete_function.triggered.connect(delete_selected_function)
            menu.addAction(delete_function)
        elif (
            isinstance(payload, ThermoRangePayload)
            and payload.owner_kind == "Function"
            and payload.property_name == "Cp"
        ):
            delete_cp = QAction(
                self._theme_icon("edit-delete", QStyle.StandardPixmap.SP_TrashIcon),
                "Delete",
                self,
            )

            def delete_selected_cp(
                _checked: bool = False,
                target_database: DatabaseIR = database,
                range_payload: ThermoRangePayload = payload,
            ) -> None:
                self._delete_function_cp_range(target_database, range_payload)

            delete_cp.triggered.connect(delete_selected_cp)
            menu.addAction(delete_cp)
        elif isinstance(payload, Phase) and _is_solution_phase(payload):
            duplicate_solution = QAction("Duplicate Solution", self)
            duplicate_solution.triggered.connect(
                lambda _checked=False,
                target_database=database,
                phase=payload: self._duplicate_database_solution(
                    target_database,
                    phase,
                )
            )
            menu.addAction(duplicate_solution)

            validate_solution = QAction("Validate Solution", self)
            validate_solution.triggered.connect(
                lambda _checked=False,
                target_database=database,
                phase=payload: self._validate_database_solution_phase(
                    target_database,
                    phase,
                )
            )
            menu.addAction(validate_solution)
            menu.addSeparator()

            delete_solution = QAction(
                self._theme_icon("edit-delete", QStyle.StandardPixmap.SP_TrashIcon),
                "Delete Solution",
                self,
            )
            delete_solution.triggered.connect(
                lambda _checked=False,
                target_database=database,
                phase=payload: self._delete_database_solution(
                    target_database,
                    phase,
                )
            )
            menu.addAction(delete_solution)
        elif isinstance(payload, SolutionParameterPayload):
            if payload.kind == "Endmember":
                add_cp = QAction("Add New Cp", self)
                function = self._solution_endmember_function_for_parameter(
                    payload.phase,
                    payload.parameter,
                )
                _h298, _s298, _cp_rows, gibbs_rows = _runtime_function_thermo_data(
                    function,
                    database,
                )
                add_cp.setEnabled(len(gibbs_rows or _cp_rows) < 6)
                add_cp.triggered.connect(
                    lambda _checked=False,
                    target_database=database,
                    parameter_payload=payload: self._add_solution_endmember_cp_range(
                        target_database,
                        parameter_payload,
                    )
                )
                menu.addAction(add_cp)
                menu.addSeparator()

            validate_parameter = QAction("Validate Parameter", self)
            validate_parameter.triggered.connect(
                lambda _checked=False,
                target_database=database,
                parameter_payload=payload: self._validate_solution_parameter_payload(
                    target_database,
                    parameter_payload,
                )
            )
            menu.addAction(validate_parameter)

            copy_expression = QAction("Copy Expression", self)
            copy_expression.triggered.connect(
                lambda _checked=False,
                parameter_payload=payload: self._copy_solution_parameter_expression(
                    parameter_payload,
                )
            )
            menu.addAction(copy_expression)
            menu.addSeparator()

            delete_parameter = QAction(
                self._theme_icon("edit-delete", QStyle.StandardPixmap.SP_TrashIcon),
                "Delete Parameter",
                self,
            )
            delete_parameter.triggered.connect(
                lambda _checked=False,
                target_database=database,
                parameter_payload=payload: self._delete_solution_parameter(
                    target_database,
                    parameter_payload,
                )
            )
            menu.addAction(delete_parameter)
        elif (
            isinstance(payload, ThermoRangePayload)
            and payload.owner_kind == "SolutionEndmember"
            and payload.property_name == "Cp"
        ):
            delete_cp = QAction(
                self._theme_icon("edit-delete", QStyle.StandardPixmap.SP_TrashIcon),
                "Delete",
                self,
            )
            delete_cp.triggered.connect(
                lambda _checked=False,
                target_database=database,
                range_payload=payload: self._delete_solution_endmember_cp_range(
                    target_database,
                    range_payload,
                )
            )
            menu.addAction(delete_cp)
        elif isinstance(payload, CompoundRecord):
            add_phase = QAction("Add New Phase", self)
            add_phase.setEnabled(len(payload.phases) < 10)
            add_phase.triggered.connect(
                lambda _checked=False, target_database=database, record=payload: (
                    self._add_compound_phase(target_database, record)
                )
            )
            menu.addAction(add_phase)
            menu.addSeparator()
            delete_compound = QAction(
                self._theme_icon("edit-delete", QStyle.StandardPixmap.SP_TrashIcon),
                "Delete Compound",
                self,
            )
            delete_compound.triggered.connect(
                lambda _checked=False,
                target_database=database,
                record=payload: self._delete_compound_record(
                    target_database,
                    record,
                )
            )
            menu.addAction(delete_compound)
        elif isinstance(payload, CompoundPhasePayload):
            add_cp = QAction("Add New Cp", self)
            _h298, _s298, _cp_rows, gibbs_rows = _compound_phase_thermo_data(
                database,
                payload.phase,
            )
            add_cp.setEnabled(len(gibbs_rows or _cp_rows) < 6)
            add_cp.triggered.connect(
                lambda _checked=False,
                target_database=database,
                phase_payload=payload: self._add_compound_cp_range(
                    target_database,
                    phase_payload,
                )
            )
            menu.addAction(add_cp)
            menu.addSeparator()
            delete_phase = QAction(
                self._theme_icon("edit-delete", QStyle.StandardPixmap.SP_TrashIcon),
                "Delete Phase",
                self,
            )
            delete_phase.triggered.connect(
                lambda _checked=False,
                target_database=database,
                phase_payload=payload: self._delete_compound_phase(
                    target_database,
                    phase_payload,
                )
            )
            menu.addAction(delete_phase)
        elif (
            isinstance(payload, ThermoRangePayload)
            and payload.owner_kind == "Compound"
            and payload.property_name == "Cp"
        ):
            delete_cp = QAction(
                self._theme_icon("edit-delete", QStyle.StandardPixmap.SP_TrashIcon),
                "Delete",
                self,
            )
            delete_cp.triggered.connect(
                lambda _checked=False,
                target_database=database,
                range_payload=payload: self._delete_compound_cp_range(
                    target_database,
                    range_payload,
                )
            )
            menu.addAction(delete_cp)

        if not menu.actions():
            return
        selected_action = menu.exec(self.tree.viewport().mapToGlobal(position))
        action_key = str(selected_action.data()) if selected_action is not None else ""
        handler = menu_handlers.get(action_key)
        if handler is not None:
            handler()

    def _edit_database_metadata(self, database: DatabaseIR) -> None:
        values = self._database_metadata_dialog_values(database)
        if values is None:
            return
        self._apply_database_metadata(database, values)
        self._rebuild_database_tree(database, select_payload=database)
        self._set_status_bar_message(
            f"{database.name}: metadata updated.",
            level="ready",
        )

    def _database_metadata_dialog_values(
        self,
        database: DatabaseIR,
    ) -> dict[str, str] | None:
        dialog = QDialog(self)
        dialog.setWindowTitle("Database Metadata")
        dialog.setModal(True)
        dialog.setMinimumWidth(520)
        layout = QVBoxLayout(dialog)
        form = QFormLayout()

        title_edit = QLineEdit(
            str(database.metadata.get("title") or database.name or "").strip()
        )
        authors_edit = QLineEdit(
            str(
                database.metadata.get("authors")
                or database.metadata.get("author")
                or ""
            ).strip()
        )
        date_edit = QLineEdit(
            str(
                database.metadata.get("date_modified")
                or database.metadata.get("date")
                or ""
            ).strip()
        )
        description_edit = QTextEdit()
        description_edit.setPlainText(str(database.metadata.get("description", "")))
        description_edit.setMinimumHeight(120)

        form.addRow("Title", title_edit)
        form.addRow("Authors", authors_edit)
        form.addRow("Date modified", date_edit)
        form.addRow("Description", description_edit)
        layout.addLayout(form)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok
            | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        layout.addWidget(buttons)

        if dialog.exec() != QDialog.DialogCode.Accepted:
            return None
        return {
            "title": title_edit.text().strip(),
            "authors": authors_edit.text().strip(),
            "date_modified": date_edit.text().strip(),
            "description": description_edit.toPlainText().strip(),
        }

    def _apply_database_metadata(
        self,
        database: DatabaseIR,
        values: dict[str, str],
    ) -> None:
        for key in ("title", "authors", "date_modified", "description"):
            value = str(values.get(key, "")).strip()
            if value:
                database.metadata[key] = value
            else:
                database.metadata.pop(key, None)

    def _resolve_database_gibbs_mismatches(self, database: DatabaseIR) -> None:
        """Apply TDB structural cleanup and balance Gibbs ranges low to high T."""
        self._commit_active_database_thermo_editor()
        structural_corrections = remove_redundant_disordered_phase_aliases(database)
        resolved_boundaries = 0

        for function in database.functions:
            _h298, _s298, _cp_rows, gibbs_rows = _runtime_function_thermo_data(
                function,
                database,
            )
            balanced_rows, count = _balance_gibbs_rows_a_coefficients_sequentially(
                gibbs_rows,
            )
            if count:
                _set_function_gibbs_rows(function, balanced_rows)
                resolved_boundaries += count

        phase_counts: dict[str, int] = {}
        for phase in database.phases:
            if _is_solution_phase(phase):
                continue
            compound_name = str(phase.name).strip() or "Compound"
            phase_counts[compound_name] = phase_counts.get(compound_name, 0) + 1
            phase_label = _compound_phase_display_label(
                phase,
                phase_counts[compound_name],
            )
            owner_name = f"{compound_name}:{phase_label}"
            for parameter in _compound_phase_g_parameters(database, phase):
                function = FunctionDefinition(
                    name=owner_name.replace(":", "_"),
                    expression=_parameter_as_function_expression(parameter),
                    source=parameter.source,
                )
                _h298, _s298, _cp_rows, gibbs_rows = _runtime_function_thermo_data(
                    function,
                    database,
                )
                balanced_rows, count = (
                    _balance_gibbs_rows_a_coefficients_sequentially(gibbs_rows)
                )
                if count:
                    _set_parameter_gibbs_rows(parameter, balanced_rows)
                    resolved_boundaries += count

        self._run_global_database_validation(database)
        self._refresh_database_diagnostics()
        self._rebuild_database_tree(database, select_payload=database)
        if resolved_boundaries or structural_corrections:
            pieces: list[str] = []
            if structural_corrections:
                pieces.append(f"{structural_corrections} structural correction(s)")
            if resolved_boundaries:
                pieces.append(f"{resolved_boundaries} Gibbs boundary mismatch(es)")
            self._set_status_bar_message(
                f"{database.name}: resolved " + " and ".join(pieces) + ".",
                level="ready",
            )
            return
        self._set_status_bar_message(
            f"{database.name}: no Gibbs boundary mismatches to resolve.",
            level="ready",
        )

    def _save_and_close_database(self, database: DatabaseIR) -> None:
        """Save a database as native .eqdb, then close it from the sidebar."""
        default_name = f"{_safe_file_stem(database.name)}.eqdb"
        path, _selected_filter = QFileDialog.getSaveFileName(
            self,
            "Save and Close Database",
            default_name,
            "Equilipy DB (*.eqdb);;JSON files (*.json);;All files (*)",
        )
        if not path:
            return
        output_path = Path(path)
        if output_path.suffix.lower() not in {".eqdb", ".json"}:
            output_path = output_path.with_suffix(".eqdb")
        if self._write_database_eqdb(
            output_path,
            title="Save and Close Database",
            database=database,
        ):
            self._close_database(database)

    def _close_database(self, database: DatabaseIR) -> None:
        """Remove one loaded database from the database workspace."""
        try:
            closed_index = next(
                index
                for index, candidate in enumerate(self.databases)
                if candidate is database
            )
        except StopIteration:
            return

        self._clear_validation_notices_for_database(database)
        self.databases.pop(closed_index)
        if not self.databases:
            blank = DatabaseIR(name="")
            self.databases = [blank]
            self.database = blank
            self.current_database = blank
            self._rebuild_database_tree(blank)
            self._show_database_summary(blank)
            self._refresh_database_diagnostics()
            return

        if self.current_database is database:
            replacement_index = min(closed_index, len(self.databases) - 1)
            self.current_database = self.databases[replacement_index]
        self.database = self.databases[0]
        self._rebuild_database_tree(self.current_database)
        self._refresh_database_diagnostics()

    def _add_database_function(self, database: DatabaseIR) -> None:
        existing = [function.name for function in database.functions]
        name = _unique_numbered_name("NEW_FUNCTION", existing, fallback="NEW_FUNCTION")
        function = FunctionDefinition(
            name=name,
            expression="298.15 0; 6000 N",
            temperature_ranges=[(298.15, 6000.0)],
        )
        database.functions.append(function)
        self._rebuild_database_tree(database, select_payload=function)

    def _duplicate_database_function(
        self,
        database: DatabaseIR,
        function: FunctionDefinition,
    ) -> FunctionDefinition | None:
        """Duplicate one Function beside the original without changing order."""
        self._commit_active_database_thermo_editor()
        if function not in database.functions:
            return None
        copied = copy.deepcopy(function)
        copied.name = self._unique_function_duplicate_name(
            function.name,
            [candidate.name for candidate in database.functions],
        )
        copied.dirty = True
        insert_at = database.functions.index(function) + 1
        database.functions.insert(insert_at, copied)
        self._rebuild_database_tree(database, select_payload=copied)
        self._set_status_bar_message(f"Duplicated function {copied.name}.", "ready")
        return copied

    def _unique_function_duplicate_name(
        self,
        name: str,
        existing_names: list[str],
    ) -> str:
        """Return a duplicate Function name that stays within TDB's 8-char limit."""
        base = re.sub(r"[^A-Za-z0-9_]", "", str(name or "").strip()) or "FUNC"
        existing = set(existing_names)
        index = 2
        while True:
            suffix = str(index)
            candidate = f"{base[: max(1, 8 - len(suffix))]}{suffix}"
            if candidate not in existing:
                return candidate
            index += 1

    def _begin_rename_database_function(
        self,
        database: DatabaseIR,
        function: FunctionDefinition,
    ) -> None:
        """Select one Function and put its title editor into rename mode."""
        self._rebuild_database_tree(database, select_payload=function)
        self._render_function(function, function.name)
        self.editor_title.setFocus(Qt.FocusReason.OtherFocusReason)
        self.editor_title.selectAll()
        self._set_status_bar_message(f"Editing function {function.name}.", "ready")

    def _validate_database_function(
        self,
        database: DatabaseIR,
        function: FunctionDefinition,
    ) -> list[Diagnostic]:
        """Validate one Function and preserve per-database validation state."""
        self._commit_active_database_thermo_editor()
        diagnostics: list[Diagnostic] = []
        name_warnings = _function_name_warning_messages(function)
        self._set_function_name_warnings(function, name_warnings)
        diagnostics.extend(
            Diagnostic(
                "warning",
                f"Function {function.name}: {message}",
                function.source,
            )
            for message in name_warnings
        )
        name_errors = _function_name_error_messages(function)
        self._set_function_name_errors(function, name_errors)
        diagnostics.extend(
            Diagnostic(
                "error",
                f"Function {function.name}: {message}",
                function.source,
            )
            for message in name_errors
        )

        try:
            _h298, _s298, _cp_rows, gibbs_rows = _runtime_function_thermo_data(
                function,
                database,
            )
        except Exception as exc:
            self._set_function_continuity_warnings(function, [], gibbs_rows=[])
            diagnostics.append(
                Diagnostic(
                    "error",
                    f"Function {function.name}: {exc}",
                    function.source,
                )
            )
            self._set_status_bar_message(
                f"Function {function.name}: validation failed.",
                "error",
            )
            return diagnostics

        if not gibbs_rows and (function.expression.strip() or function.gibbs_ranges):
            diagnostics.append(
                Diagnostic(
                    "warning",
                    (
                        f"Function {function.name}: could not be converted "
                        "into H/S/Cp/G rows."
                    ),
                    function.source,
                )
            )
        warnings = _function_gibbs_continuity_warnings(gibbs_rows)
        self._set_function_continuity_warnings(
            function,
            warnings,
            gibbs_rows=gibbs_rows,
        )
        diagnostics.extend(
            Diagnostic(
                "warning",
                f"Function {function.name}: {message}",
                function.source,
            )
            for message in warnings
        )

        if name_errors:
            self._set_status_bar_message(
                f"Function {function.name}: {name_errors[0]}",
                "error",
            )
        elif diagnostics:
            self._set_status_bar_message(
                f"Function {function.name}: {diagnostics[0].message}",
                diagnostics[0].severity,
            )
        else:
            self._set_status_bar_message(
                f"Function {function.name}: validated.",
                "ready",
            )
        return diagnostics

    def _copy_database_function_expression(
        self,
        function: FunctionDefinition,
    ) -> str:
        """Copy the editable expression text for one Function to the clipboard."""
        expression = str(function.tdb_src.source_expression or function.expression)
        QApplication.clipboard().setText(expression)
        self._set_status_bar_message(
            f"Function {function.name}: expression copied.",
            "ready",
        )
        return expression

    def _delete_database_function(
        self,
        database: DatabaseIR,
        function: FunctionDefinition,
    ) -> None:
        """Remove a function from a database and keep the tree selection valid."""
        if function in database.functions:
            database.functions.remove(function)
        else:
            database.functions = [
                candidate
                for candidate in database.functions
                if candidate.name != function.name
            ]
        if self._editable_editor_title_record is function:
            self._editable_editor_title_record = None
        self._rebuild_database_tree(database, select_payload=database.functions)

    def _add_database_solution(self, database: DatabaseIR) -> None:
        """Add a minimal solution phase that can be edited further."""
        existing = [phase.name for phase in database.phases]
        name = _unique_numbered_name("NEW_SOLUTION", existing, fallback="NEW_SOLUTION")
        phase = Phase(
            name=name,
            model="CEF",
            constituents=[ConstituentSet(1, ["Va"], 1.0)],
        )
        database.phases.append(phase)
        self._rebuild_database_tree(database, select_payload=phase)
        self._set_status_bar_message(f"Added solution {phase.name}.", "ready")

    def _duplicate_database_solution(
        self,
        database: DatabaseIR,
        phase: Phase,
    ) -> Phase | None:
        """Duplicate one solution phase and its attached parameters."""
        self._commit_active_database_thermo_editor()
        if phase not in database.phases or not _is_solution_phase(phase):
            return None
        copied = copy.deepcopy(phase)
        copied.name = _unique_numbered_name(
            f"{phase.name}_COPY",
            [candidate.name for candidate in database.phases],
            fallback=f"{phase.name}_COPY",
        )
        insert_at = database.phases.index(phase) + 1
        database.phases.insert(insert_at, copied)
        copied_parameters = []
        for parameter in database.parameters:
            if parameter.phase_name.upper() != phase.name.upper():
                continue
            copied_parameter = copy.deepcopy(parameter)
            copied_parameter.phase_name = copied.name
            copied_parameters.append(copied_parameter)
        database.parameters.extend(copied_parameters)
        self._rebuild_database_tree(database, select_payload=copied)
        self._set_status_bar_message(f"Duplicated solution {copied.name}.", "ready")
        return copied

    def _delete_database_solution(
        self,
        database: DatabaseIR,
        phase: Phase,
    ) -> None:
        """Delete one solution phase and its attached parameters."""
        if phase in database.phases:
            database.phases.remove(phase)
        else:
            database.phases = [
                candidate
                for candidate in database.phases
                if candidate.name.upper() != phase.name.upper()
            ]
        database.parameters = [
            parameter
            for parameter in database.parameters
            if parameter.phase_name.upper() != phase.name.upper()
        ]
        self._rebuild_database_tree(database, select_payload=database.phases)
        self._set_status_bar_message(f"Deleted solution {phase.name}.", "ready")

    def _validate_database_solution_phase(
        self,
        database: DatabaseIR,
        phase: Phase,
    ) -> list[Diagnostic]:
        """Validate required endmember parameters for one solution phase."""
        diagnostics: list[Diagnostic] = []
        parameters = _solution_parameters_for_phase(database, phase)
        defined_endmembers = {
            _solution_parameter_target_text(parameter).upper()
            for parameter in parameters
            if not _solution_parameter_is_interaction(parameter)
        }
        for target in _solution_expected_endmembers(phase):
            normalized_target = ":".join(target).upper()
            if normalized_target in defined_endmembers:
                continue
            diagnostics.append(
                Diagnostic(
                    "warning",
                    (
                        f"Solution {phase.name}: missing endmember "
                        f"{':'.join(target)}."
                    ),
                    phase.source,
                )
            )
        if diagnostics:
            first = diagnostics[0]
            self._set_status_bar_message(first.message, first.severity)
        else:
            self._set_status_bar_message(
                f"Solution {phase.name}: validated.",
                "ready",
            )
        return diagnostics

    def _add_function_cp_range(
        self,
        database: DatabaseIR,
        function: FunctionDefinition,
    ) -> None:
        gibbs_rows = _direct_function_gibbs_rows(function)
        if not gibbs_rows:
            lower = (
                function.temperature_ranges[-1][1]
                if function.temperature_ranges
                else 298.15
            )
            gibbs_rows = [[float(lower), float(lower) + 100.0, *([0.0] * 14)]]
        if len(gibbs_rows) >= 6:
            return
        last = [float(value) for value in gibbs_rows[-1]]
        lower = float(last[1])
        last[0] = lower
        last[1] = lower + 100.0
        gibbs_rows.append(last)
        _set_function_gibbs_rows(function, gibbs_rows)
        _normalize_function_temperature_ranges(function)
        self._select_function_cp_range(
            database,
            function,
            len(function.temperature_ranges) - 1,
        )

    def _delete_function_cp_range(
        self,
        database: DatabaseIR,
        payload: ThermoRangePayload,
    ) -> None:
        function = _function_by_name(database, payload.owner_name)
        if function is None:
            return
        gibbs_rows = _direct_function_gibbs_rows(function)
        remove_index = max(0, payload.index - 1)
        if remove_index < len(gibbs_rows):
            gibbs_rows.pop(remove_index)
        if gibbs_rows:
            _set_function_gibbs_rows(function, gibbs_rows)
            _normalize_function_temperature_ranges(function)
            self._select_function_cp_range(
                database,
                function,
                min(remove_index, len(function.temperature_ranges) - 1),
            )
            return
        else:
            function.expression = ""
            function.temperature_ranges.clear()
            function.gibbs_ranges.clear()
        self._rebuild_database_tree(database, select_payload=function)

    def _add_database_compound(self, database: DatabaseIR) -> None:
        """Append a new stoichiometric compound scaffold to a database."""
        existing = [str(phase.name).strip() for phase in database.phases]
        name = _unique_numbered_name("NEW_COMPOUND", existing, fallback="NEW_COMPOUND")
        source = SourceRef(
            command=f"PHASE {name} % 1 1 ! CONST {name} :{name}: !",
        )
        phase = Phase(
            name=name,
            model="COMP",
            constituents=[ConstituentSet(1, [name], 1.0)],
            source=source,
            metadata={
                "phase_id": _new_compound_phase_state_id(database, name),
                "phase_label": "S1",
            },
        )
        database.species.append(Species(name, {}, source=source))
        database.phases.append(phase)
        database.parameters.append(_new_compound_g_parameter(phase))
        self._rebuild_database_tree(database, show_summary=False)
        self._select_compound_phase(phase)

    def _add_compound_phase(
        self,
        database: DatabaseIR,
        record: CompoundRecord,
    ) -> None:
        """Append one phase state under an existing compound record."""
        if len(record.phases) >= 10:
            return
        template = record.phases[0] if record.phases else None
        phase = Phase(
            name=record.name,
            model=(template.model if template is not None else "COMP") or "COMP",
            constituents=copy.deepcopy(
                template.constituents
                if template is not None
                else [ConstituentSet(1, [record.name], 1.0)]
            ),
            source=SourceRef(command=f"PHASE {record.name} % 1 1 !"),
            metadata={
                "phase_id": _new_compound_phase_state_id(database, record.name),
                "phase_label": f"S{len(record.phases) + 1}",
            },
        )
        database.phases.append(phase)
        database.parameters.append(
            self._new_compound_g_parameter_from_base(
                database,
                phase,
                template,
            )
        )
        self._rebuild_database_tree(database, show_summary=False)
        self._select_compound_phase(phase)

    def _new_compound_g_parameter_from_base(
        self,
        database: DatabaseIR,
        phase: Phase,
        base_phase: Phase | None,
    ) -> Parameter:
        """Return a new compound Gibbs parameter, inheriting base Cp ranges if present."""
        if base_phase is None:
            return _new_compound_g_parameter(phase)
        base_parameters = _compound_phase_g_parameters(database, base_phase)
        if not base_parameters:
            return _new_compound_g_parameter(phase)

        base_parameter = base_parameters[0]
        expression = _parameter_as_function_expression(base_parameter)
        if not expression:
            return _new_compound_g_parameter(phase)
        function = FunctionDefinition(
            name=f"{phase.name}_base_cp",
            expression=expression,
            source=base_parameter.source,
        )
        _h298, _s298, cp_rows, _gibbs_rows = _runtime_function_thermo_data(
            function,
            database,
        )
        try:
            gibbs_rows = (
                HSCp2G(0.0, 0.0, np.asarray(cp_rows, dtype=float)).tolist()
                if cp_rows
                else []
            )
        except (AssertionError, ValueError, TypeError):
            gibbs_rows = []
        target = _compound_parameter_target_species(phase)
        target_text = ":".join(target)
        source = SourceRef(
            command=f"PARAMETER G({phase.name},{target_text};0) {expression} !",
        )
        parameter = Parameter(
            phase_name=phase.name,
            parameter_type="G",
            target=target,
            order=base_parameter.order,
            expression=expression,
            source=source,
            metadata=(
                {"phase_id": _compound_phase_state_id(phase)}
                if _compound_phase_state_id(phase)
                else {}
            ),
        )
        if gibbs_rows:
            _set_parameter_gibbs_rows(parameter, gibbs_rows)
        return parameter

    def _delete_compound_phase(
        self,
        database: DatabaseIR,
        payload: CompoundPhasePayload,
    ) -> None:
        """Remove one compound phase state from the database."""
        phase = payload.phase
        try:
            removed_index = database.phases.index(phase)
        except ValueError:
            removed_index = None
        if phase in database.phases:
            database.phases.remove(phase)
        phase_id = _compound_phase_state_id(phase)
        if phase_id:
            database.parameters = [
                parameter
                for parameter in database.parameters
                if parameter.metadata.get("phase_id") != phase_id
            ]
            self._preserve_compound_order_after_phase_delete(
                database,
                phase.name,
                removed_index,
            )
            self._rebuild_database_tree(database, show_summary=False)
            self._show_database_summary(database)
            return
        remaining_same_name = [
            candidate for candidate in database.phases if candidate.name == phase.name
        ]
        if not remaining_same_name:
            database.parameters = [
                parameter
                for parameter in database.parameters
                if parameter.phase_name != phase.name
            ]
        self._preserve_compound_order_after_phase_delete(
            database,
            phase.name,
            removed_index,
        )
        self._rebuild_database_tree(database, show_summary=False)
        self._show_database_summary(database)

    def _preserve_compound_order_after_phase_delete(
        self,
        database: DatabaseIR,
        compound_name: str,
        removed_index: int | None,
    ) -> None:
        """Keep a compound at its old tree position when deleting one phase state."""
        if removed_index is None:
            return
        same_name_phases = [
            phase
            for phase in database.phases
            if phase.name == compound_name and not _is_solution_phase(phase)
        ]
        if not same_name_phases:
            return
        first_remaining_index = min(database.phases.index(phase) for phase in same_name_phases)
        if first_remaining_index <= removed_index:
            return
        database.phases = [
            phase
            for phase in database.phases
            if not (phase.name == compound_name and not _is_solution_phase(phase))
        ]
        insert_index = min(removed_index, len(database.phases))
        for offset, phase in enumerate(same_name_phases):
            database.phases.insert(insert_index + offset, phase)

    def _delete_compound_record(
        self,
        database: DatabaseIR,
        record: CompoundRecord,
    ) -> None:
        """Remove one compound and all attached phase-state parameters."""
        self._remove_compound_record(database, record)
        self._rebuild_database_tree(database, show_summary=False)
        self._show_database_summary(database)

    def _commit_compound_phase_label(
        self,
        payload: CompoundPhasePayload,
        editor: QLineEdit,
    ) -> None:
        """Persist an edited compound phase-state label."""
        label = str(editor.text()).strip() or payload.phase_label
        payload.phase.metadata["phase_label"] = label
        editor.setText(label)
        self._rebuild_database_tree(self.current_database, show_summary=False)
        self._select_compound_phase(payload.phase)

    def _add_compound_cp_range(
        self,
        database: DatabaseIR,
        payload: CompoundPhasePayload,
    ) -> None:
        """Append a new Cp/G range to one compound phase parameter."""
        parameter = _compound_phase_g_parameters(database, payload.phase)
        if parameter:
            target_parameter = parameter[0]
        else:
            target_parameter = _new_compound_g_parameter(payload.phase)
            database.parameters.append(target_parameter)
        function = FunctionDefinition(
            name=f"{payload.compound_name}_{payload.phase_label}",
            expression=_parameter_as_function_expression(target_parameter),
            source=target_parameter.source,
        )
        gibbs_rows = _direct_function_gibbs_rows(function)
        if not gibbs_rows:
            _h298, _s298, _cp_rows, gibbs_rows = _runtime_function_thermo_data(
                function,
                database,
            )
        if len(gibbs_rows) >= 6:
            return
        if not gibbs_rows:
            gibbs_rows = [[298.15, 6000.0, *([0.0] * 14)]]
        last = [float(value) for value in gibbs_rows[-1]]
        lower = float(last[1])
        last[0] = lower
        last[1] = lower + 100.0
        gibbs_rows.append(last)
        _set_parameter_gibbs_rows(target_parameter, gibbs_rows)
        self._rebuild_database_tree(database, show_summary=False)
        self._select_compound_cp_range(
            f"{payload.compound_name}:{payload.phase_label}",
            len(gibbs_rows),
        )

    def _delete_compound_cp_range(
        self,
        database: DatabaseIR,
        payload: ThermoRangePayload,
    ) -> None:
        """Delete one compound Cp/G range from its backing parameter."""
        parameter = _compound_parameter_for_range_payload(database, payload)
        if parameter is None:
            return
        function = FunctionDefinition(
            name=payload.owner_name.replace(":", "_"),
            expression=_parameter_as_function_expression(parameter),
            source=parameter.source,
        )
        gibbs_rows = _direct_function_gibbs_rows(function)
        if not gibbs_rows:
            _h298, _s298, _cp_rows, gibbs_rows = _runtime_function_thermo_data(
                function,
                database,
            )
        remove_index = max(0, payload.index - 1)
        if remove_index < len(gibbs_rows):
            gibbs_rows.pop(remove_index)
        if gibbs_rows:
            _set_parameter_gibbs_rows(parameter, gibbs_rows)
        else:
            database.parameters.remove(parameter)
        self._rebuild_database_tree(database, show_summary=False)
        self._select_compound_cp_range(
            payload.owner_name,
            min(remove_index + 1, max(1, len(gibbs_rows))),
        )

    def _rebuild_database_tree(
        self,
        database: DatabaseIR | None = None,
        *,
        select_payload: Any | None = None,
        show_summary: bool = True,
    ) -> None:
        previous_expanded_paths = self._database_tree_expanded_paths()
        had_visible_database_roots = self.tree_model.rowCount(QModelIndex()) > 0
        self.database = self.databases[0]
        self.current_database = database or self.current_database
        self.tree_model = self._new_database_tree_model()
        self.tree.setModel(self.tree_model)
        self.tree.setHeaderHidden(True)
        if had_visible_database_roots:
            self._restore_database_tree_expansion(previous_expanded_paths)
        else:
            self._expand_database_roots()
        self.tree.selectionModel().currentChanged.connect(self._show_selected)
        if select_payload is not None:
            index = self._find_database_tree_index(select_payload)
            if index.isValid():
                self.tree.setCurrentIndex(index)
                self.tree.scrollTo(index)
                return
        if database is not None and show_summary:
            self._show_database_summary(database)

    def _database_tree_expanded_paths(self) -> set[tuple[str, ...]]:
        """Return stable label paths for currently expanded database tree nodes."""
        if self.tree is None or self.tree_model is None:
            return set()
        expanded: set[tuple[str, ...]] = set()

        def visit(parent: QModelIndex, parent_path: tuple[str, ...] = ()) -> None:
            for row in range(self.tree_model.rowCount(parent)):
                index = self.tree_model.index(row, 0, parent)
                path = (*parent_path, self.tree_model.label(index))
                if self.tree.isExpanded(index):
                    expanded.add(path)
                visit(index, path)

        visit(QModelIndex())
        return expanded

    def _restore_database_tree_expansion(
        self,
        expanded_paths: set[tuple[str, ...]],
    ) -> None:
        """Restore expanded database tree branches after replacing the model."""

        def visit(parent: QModelIndex, parent_path: tuple[str, ...] = ()) -> None:
            for row in range(self.tree_model.rowCount(parent)):
                index = self.tree_model.index(row, 0, parent)
                path = (*parent_path, self.tree_model.label(index))
                self.tree.setExpanded(index, path in expanded_paths)
                visit(index, path)

        visit(QModelIndex())

    def _expand_database_roots(self) -> None:
        """Default to expanded database roots with folded object groups."""
        for row in range(self.tree_model.rowCount(QModelIndex())):
            index = self.tree_model.index(row, 0, QModelIndex())
            self.tree.setExpanded(index, True)

    def _new_database_tree_model(self) -> DatabaseTreeModel:
        """Return a database tree model with current function warning metadata."""
        return DatabaseTreeModel(
            self.databases,
            function_warnings=self._tree_function_warnings(),
            function_errors=self.function_name_errors,
            range_warnings=self._tree_range_warnings(),
        )

    def _tree_function_warnings(self) -> dict[int, list[str]]:
        """Return Function warnings shown by the database tree."""
        warnings_by_id = {
            key: list(warnings)
            for key, warnings in self.function_name_warnings.items()
        }
        for key, warnings in self.function_continuity_warnings.items():
            warnings_by_id.setdefault(key, []).extend(warnings)
        return warnings_by_id

    def _tree_range_warnings(self) -> dict[str, list[str]]:
        """Return all Cp-range warnings shown by the database tree."""
        return {
            **self.function_range_warnings,
            **self.compound_range_warnings,
        }

    def _refresh_tree_model_range_warnings(self) -> None:
        """Synchronize live range-warning dictionaries into the tree model."""
        self.tree_model.range_warnings = self._tree_range_warnings()

    def _function_name_errors(self) -> dict[int, list[str]]:
        """Return strict TDB function-name errors keyed by Function object identity."""
        errors_by_id: dict[int, list[str]] = {}
        for database in self.databases:
            for function in database.functions:
                errors = _function_name_error_messages(function)
                if errors:
                    errors_by_id[id(function)] = errors
        return errors_by_id

    def _function_name_warnings(self) -> dict[int, list[str]]:
        """Return TDB function-name warnings keyed by Function object identity."""
        warnings_by_id: dict[int, list[str]] = {}
        for database in self.databases:
            for function in database.functions:
                warnings = _function_name_warning_messages(function)
                if warnings:
                    warnings_by_id[id(function)] = warnings
        return warnings_by_id

    def _set_function_name_warnings(
        self,
        function: FunctionDefinition,
        warnings: list[str],
    ) -> None:
        """Update live function-name warning state for one Function row."""
        if warnings:
            self.function_name_warnings[id(function)] = list(warnings)
        else:
            self.function_name_warnings.pop(id(function), None)
        index = self._find_database_tree_index(function)
        if index.isValid():
            self.tree_model.dataChanged.emit(index, index)
        if self._current_selection_function() is function:
            self._set_function_status_bar_message(function)

    def _set_function_name_errors(
        self,
        function: FunctionDefinition,
        errors: list[str],
    ) -> None:
        """Update live function-name error state for one Function row in the tree."""
        if errors:
            self.function_name_errors[id(function)] = list(errors)
        else:
            self.function_name_errors.pop(id(function), None)
        index = self._find_database_tree_index(function)
        if index.isValid():
            self.tree_model.dataChanged.emit(index, index)
        if self._current_selection_function() is function:
            self._set_function_status_bar_message(function)

    def _function_continuity_warnings(self) -> dict[int, list[str]]:
        """Return Gibbs continuity warnings keyed by Function object identity."""
        warnings_by_id: dict[int, list[str]] = {}
        for database in self.databases:
            for function in database.functions:
                gibbs_rows = _direct_function_gibbs_rows(function)
                warnings = _function_gibbs_continuity_warnings(gibbs_rows)
                if warnings:
                    warnings_by_id[id(function)] = warnings
        return warnings_by_id

    def _function_range_continuity_warnings(self) -> dict[str, list[str]]:
        """Return Gibbs warnings keyed by function Cp-range tree payload."""
        warnings_by_key: dict[str, list[str]] = {}
        for database in self.databases:
            for function in database.functions:
                warnings_by_key.update(
                    self._function_range_warnings_for_function(function)
                )
        return warnings_by_key

    def _function_range_warnings_for_function(
        self,
        function: FunctionDefinition,
        gibbs_rows: list[list[float]] | None = None,
    ) -> dict[str, list[str]]:
        """Return selected-boundary warnings for each Cp range of one function."""
        warnings_by_key: dict[str, list[str]] = {}
        gibbs_rows = (
            _direct_function_gibbs_rows(function)
            if gibbs_rows is None
            else gibbs_rows
        )
        ranges = list(function.temperature_ranges[:6])
        if not ranges and function.expression.strip():
            ranges = [("", "")]
        for index, (lower, upper) in enumerate(ranges, start=1):
            payload = ThermoRangePayload(
                owner_name=function.name,
                owner_kind="Function",
                property_name="Cp",
                index=index,
                expression=function.expression,
                lower=lower,
                upper=upper,
                source=function.source,
            )
            key = thermo_range_warning_key(payload)
            if not key:
                continue
            warnings = _selected_gibbs_continuity_warnings(gibbs_rows, index - 1)
            if warnings:
                warnings_by_key[key] = warnings
        return warnings_by_key

    def _compound_range_direct_continuity_warnings(self) -> dict[str, list[str]]:
        """Return compound warnings that do not require H/S/Cp/G conversion."""
        warnings_by_key: dict[str, list[str]] = {}
        for database in self.databases:
            phase_counts: dict[str, int] = {}
            for phase in database.phases:
                if _is_solution_phase(phase):
                    continue
                compound_name = str(phase.name).strip() or "Compound"
                phase_counts[compound_name] = phase_counts.get(compound_name, 0) + 1
                phase_label = _compound_phase_display_label(
                    phase,
                    phase_counts[compound_name],
                )
                owner_name = f"{compound_name}:{phase_label}"
                for parameter in _compound_phase_g_parameters(database, phase):
                    function = FunctionDefinition(
                        name=owner_name.replace(":", "_"),
                        expression=_parameter_as_function_expression(parameter),
                        source=parameter.source,
                    )
                    gibbs_rows = _direct_function_gibbs_rows(function)
                    if not gibbs_rows:
                        continue
                    for index, row in enumerate(gibbs_rows, start=1):
                        payload = ThermoRangePayload(
                            owner_name=owner_name,
                            owner_kind="Compound",
                            property_name="Cp",
                            index=index,
                            expression=function.expression,
                            lower=float(row[0]),
                            upper=float(row[1]),
                            source=parameter.source,
                        )
                        key = thermo_range_warning_key(payload)
                        if not key:
                            continue
                        warnings = _selected_gibbs_continuity_warnings(
                            gibbs_rows,
                            index - 1,
                        )
                        if warnings:
                            warnings_by_key[key] = warnings
        return warnings_by_key

    def _set_function_continuity_warnings(
        self,
        function: FunctionDefinition,
        warnings: list[str],
        _selected_warnings: list[str] | None = None,
        gibbs_rows: list[list[float]] | None = None,
    ) -> None:
        """Update live warning state for one Function row in the tree."""
        if warnings:
            self.function_continuity_warnings[id(function)] = list(warnings)
        else:
            self.function_continuity_warnings.pop(id(function), None)
        self._set_function_range_continuity_warnings(
            function,
            gibbs_rows=gibbs_rows,
        )
        self._refresh_tree_model_range_warnings()
        index = self._find_database_tree_index(function)
        if index.isValid():
            self._emit_tree_index_and_ancestors_changed(index)
            for row in range(self.tree_model.rowCount(index)):
                child_index = self.tree_model.index(row, 0, index)
                self.tree_model.dataChanged.emit(child_index, child_index)
        if self._current_selection_function() is function:
            self._set_function_status_bar_message(function)

    def _set_function_range_continuity_warnings(
        self,
        function: FunctionDefinition,
        *,
        gibbs_rows: list[list[float]] | None = None,
    ) -> None:
        """Refresh range-warning metadata for one function."""
        prefix = f"Function|{function.name}|"
        self.function_range_warnings = {
            key: warnings
            for key, warnings in self.function_range_warnings.items()
            if not key.startswith(prefix)
        }
        self.function_range_warnings.update(
            self._function_range_warnings_for_function(
                function,
                gibbs_rows=gibbs_rows,
            )
        )

    def _set_compound_range_continuity_warnings(
        self,
        payload: ThermoRangePayload,
        function: FunctionDefinition,
        *,
        gibbs_rows: list[list[float]] | None = None,
    ) -> None:
        """Update live warning state for one compound Cp range in the tree."""
        if gibbs_rows is None:
            function = self._compound_function_for_range_payload(payload)
            _h298, _s298, _cp_rows, gibbs_rows = _runtime_function_thermo_data(
                function,
                self.current_database,
            )
            self._set_compound_range_warnings_from_gibbs_rows(
                payload,
                function,
                gibbs_rows,
            )
        else:
            self._set_compound_range_warnings_from_gibbs_rows(
                payload,
                function,
                gibbs_rows,
            )
        self._refresh_tree_model_range_warnings()
        key = thermo_range_warning_key(payload)
        if not key:
            return
        index = self._find_thermo_range_tree_index(key)
        if index.isValid():
            self._emit_tree_index_and_ancestors_changed(index)
            parent = index.parent()
            if parent.isValid():
                self._emit_tree_subtree_changed(parent)
            if self.tree.currentIndex() == index:
                self._set_compound_range_status_bar_message(
                    self.tree_model.payload(index),
                    self.tree_model.label(index),
                )

    def _set_compound_range_warnings_from_gibbs_rows(
        self,
        payload: ThermoRangePayload,
        function: FunctionDefinition,
        gibbs_rows: list[list[float]],
    ) -> None:
        """Refresh one compound phase's tree warnings from live editor rows."""
        prefix = f"Compound|{payload.owner_name}|"
        self.compound_range_warnings = {
            key: warnings
            for key, warnings in self.compound_range_warnings.items()
            if not key.startswith(prefix)
        }
        source = payload.source if isinstance(payload.source, SourceRef) else None
        for index, row in enumerate(gibbs_rows, start=1):
            range_payload = ThermoRangePayload(
                owner_name=payload.owner_name,
                owner_kind="Compound",
                property_name="Cp",
                index=index,
                expression=function.expression,
                lower=float(row[0]),
                upper=float(row[1]),
                source=source,
            )
            key = thermo_range_warning_key(range_payload)
            if not key:
                continue
            warnings = _selected_gibbs_continuity_warnings(
                gibbs_rows,
                index - 1,
            )
            if warnings:
                self.compound_range_warnings[key] = warnings

    def _refresh_compound_phase_warning_cache(
        self,
        database: DatabaseIR,
        phase: Phase,
        owner_name: str,
        *,
        emit: bool = True,
    ) -> None:
        """Compute compound warnings only for the phase being inspected."""
        parameters = _compound_phase_g_parameters(database, phase)
        if not parameters:
            return
        parameter = parameters[0]
        function = FunctionDefinition(
            name=owner_name.replace(":", "_"),
            expression=_parameter_as_function_expression(parameter),
            source=parameter.source,
        )
        _h298, _s298, _cp_rows, gibbs_rows = _runtime_function_thermo_data(
            function,
            database,
        )
        payload = ThermoRangePayload(
            owner_name=owner_name,
            owner_kind="Compound",
            property_name="Cp",
            index=1,
            expression=function.expression,
            source=parameter.source,
        )
        self._set_compound_range_warnings_from_gibbs_rows(
            payload,
            function,
            gibbs_rows,
        )
        self._refresh_tree_model_range_warnings()
        if not emit:
            return
        phase_index = self._find_compound_phase_tree_index(phase)
        if phase_index.isValid():
            self._emit_tree_index_and_ancestors_changed(phase_index)
            self._emit_tree_subtree_changed(phase_index)
        current = self.tree.currentIndex()
        if current.isValid():
            current_payload = self.tree_model.payload(current)
            if (
                isinstance(current_payload, CompoundPhasePayload)
                and current_payload.phase is phase
            ):
                warnings = self._compound_phase_status_warnings(current_payload)
                if warnings:
                    self._set_status_bar_message(
                        f"Compound {current_payload.compound_name} phase "
                        f"{current_payload.phase_label}: {warnings[0]}",
                        "warning",
                    )

    def _emit_tree_index_and_ancestors_changed(self, index: QModelIndex) -> None:
        """Tell the view to redraw one tree node and its visible parents."""
        while index.isValid():
            self.tree_model.dataChanged.emit(index, index)
            index = index.parent()

    def _emit_tree_subtree_changed(self, index: QModelIndex) -> None:
        """Tell the view to redraw one tree node and all visible descendants."""
        if not index.isValid():
            return
        self.tree_model.dataChanged.emit(index, index)
        for row in range(self.tree_model.rowCount(index)):
            child = self.tree_model.index(row, 0, index)
            self._emit_tree_subtree_changed(child)

    def _find_database_tree_index(self, payload: Any) -> QModelIndex:
        def visit(parent: QModelIndex) -> QModelIndex:
            for row in range(self.tree_model.rowCount(parent)):
                index = self.tree_model.index(row, 0, parent)
                if self.tree_model.payload(index) is payload:
                    return index
                found = visit(index)
                if found.isValid():
                    return found
            return QModelIndex()

        return visit(QModelIndex())

    def _find_thermo_range_tree_index(self, key: str) -> QModelIndex:
        """Return the tree index for one thermo range warning key."""
        def visit(parent: QModelIndex) -> QModelIndex:
            for row in range(self.tree_model.rowCount(parent)):
                index = self.tree_model.index(row, 0, parent)
                if thermo_range_warning_key(self.tree_model.payload(index)) == key:
                    return index
                found = visit(index)
                if found.isValid():
                    return found
            return QModelIndex()

        return visit(QModelIndex())

    def _find_function_cp_tree_index(
        self,
        function_name: str,
        range_index: int,
    ) -> QModelIndex:
        target_index = range_index + 1

        def visit(parent: QModelIndex) -> QModelIndex:
            for row in range(self.tree_model.rowCount(parent)):
                index = self.tree_model.index(row, 0, parent)
                payload = self.tree_model.payload(index)
                if (
                    isinstance(payload, ThermoRangePayload)
                    and payload.owner_kind == "Function"
                    and payload.property_name == "Cp"
                    and payload.owner_name == function_name
                    and payload.index == target_index
                ):
                    return index
                found = visit(index)
                if found.isValid():
                    return found
            return QModelIndex()

        return visit(QModelIndex())

    def _find_compound_phase_tree_index(self, phase: Phase) -> QModelIndex:
        """Return the tree index wrapping one compound phase object."""
        def visit(parent: QModelIndex) -> QModelIndex:
            for row in range(self.tree_model.rowCount(parent)):
                index = self.tree_model.index(row, 0, parent)
                payload = self.tree_model.payload(index)
                if (
                    isinstance(payload, CompoundPhasePayload)
                    and payload.phase is phase
                ):
                    return index
                found = visit(index)
                if found.isValid():
                    return found
            return QModelIndex()

        return visit(QModelIndex())

    def _find_compound_cp_tree_index(
        self,
        owner_name: str,
        range_index: int,
    ) -> QModelIndex:
        """Return one compound Cp-range index by owner and 1-based range index."""
        def visit(parent: QModelIndex) -> QModelIndex:
            for row in range(self.tree_model.rowCount(parent)):
                index = self.tree_model.index(row, 0, parent)
                payload = self.tree_model.payload(index)
                if (
                    isinstance(payload, ThermoRangePayload)
                    and payload.owner_kind == "Compound"
                    and payload.property_name == "Cp"
                    and payload.owner_name == owner_name
                    and payload.index == range_index
                ):
                    return index
                found = visit(index)
                if found.isValid():
                    return found
            return QModelIndex()

        return visit(QModelIndex())

    def _find_solution_endmember_cp_tree_index(
        self,
        owner_name: str,
        range_index: int,
    ) -> QModelIndex:
        """Return one solution endmember Cp-range index by owner and range."""
        def visit(parent: QModelIndex) -> QModelIndex:
            for row in range(self.tree_model.rowCount(parent)):
                index = self.tree_model.index(row, 0, parent)
                payload = self.tree_model.payload(index)
                if (
                    isinstance(payload, ThermoRangePayload)
                    and payload.owner_kind == "SolutionEndmember"
                    and payload.property_name == "Cp"
                    and payload.owner_name == owner_name
                    and payload.index == range_index
                ):
                    return index
                found = visit(index)
                if found.isValid():
                    return found
            return QModelIndex()

        return visit(QModelIndex())

    def _find_solution_parameter_tree_index(self, parameter: Parameter) -> QModelIndex:
        """Return the tree index wrapping one solution parameter object."""
        def visit(parent: QModelIndex) -> QModelIndex:
            for row in range(self.tree_model.rowCount(parent)):
                index = self.tree_model.index(row, 0, parent)
                payload = self.tree_model.payload(index)
                if (
                    isinstance(payload, SolutionParameterPayload)
                    and payload.parameter is parameter
                ):
                    return index
                found = visit(index)
                if found.isValid():
                    return found
            return QModelIndex()

        return visit(QModelIndex())

    def _select_function_cp_range(
        self,
        database: DatabaseIR,
        function: FunctionDefinition,
        range_index: int,
    ) -> None:
        self._rebuild_database_tree(database, show_summary=False)
        index = self._find_function_cp_tree_index(function.name, range_index)
        if index.isValid():
            self.tree.setCurrentIndex(index)
            self.tree.scrollTo(index)
            return
        function_index = self._find_database_tree_index(function)
        if function_index.isValid():
            self.tree.setCurrentIndex(function_index)
            self.tree.scrollTo(function_index)
            return
        self._show_database_summary(database)

    def _select_compound_phase(self, phase: Phase) -> None:
        """Select a compound phase wrapper after rebuilding the tree."""
        index = self._find_compound_phase_tree_index(phase)
        if index.isValid():
            self.tree.setCurrentIndex(index)
            self.tree.scrollTo(index)

    def _select_compound_cp_range(self, owner_name: str, range_index: int) -> None:
        """Select a compound Cp range after rebuilding the tree."""
        index = self._find_compound_cp_tree_index(owner_name, range_index)
        if index.isValid():
            self.tree.setCurrentIndex(index)
            self.tree.scrollTo(index)

    def _select_solution_endmember_cp_range(
        self,
        owner_name: str,
        range_index: int,
    ) -> None:
        """Select a solution endmember Cp range after rebuilding the tree."""
        index = self._find_solution_endmember_cp_tree_index(owner_name, range_index)
        if index.isValid():
            self.tree.setCurrentIndex(index)
            self.tree.scrollTo(index)

    def _select_solution_parameter(self, parameter: Parameter) -> None:
        """Select a solution parameter after rebuilding the tree."""
        index = self._find_solution_parameter_tree_index(parameter)
        if index.isValid():
            self.tree.setCurrentIndex(index)
            self.tree.scrollTo(index)



    def _restore_tree_selection(self) -> None:
        if self.last_tree_index is None or not self.last_tree_index.isValid():
            return
        index = QModelIndex(self.last_tree_index)
        self.tree.selectionModel().select(
            index,
            QItemSelectionModel.SelectionFlag.ClearAndSelect
            | QItemSelectionModel.SelectionFlag.Rows,
        )
        self.tree.setCurrentIndex(index)



    def _render_payload(self, payload: Any, display_name: str = "") -> None:
        if isinstance(payload, DatabaseIR):
            self._show_database_summary(payload)
        elif isinstance(payload, list):
            self._render_list(payload, display_name)
        elif isinstance(payload, Element):
            self._render_element(payload, display_name)
        elif isinstance(payload, Species):
            self._render_species(payload, display_name)
        elif isinstance(payload, ThermoRangePayload):
            self._render_thermo_range(payload, display_name)
        elif isinstance(payload, CompoundRecord):
            self._render_compound_record(payload, display_name)
        elif isinstance(payload, CompoundPhasePayload):
            self._render_compound_phase_payload(payload, display_name)
        elif isinstance(payload, FunctionDefinition):
            self._render_function(payload, display_name)
        elif isinstance(payload, SolutionGroupPayload):
            self._render_solution_group(payload, display_name)
        elif isinstance(payload, SolutionSublatticePayload):
            self._render_solution_sublattice(payload, display_name)
        elif isinstance(payload, SolutionSpeciesPayload):
            self._render_solution_species(payload, display_name)
        elif isinstance(payload, SolutionParameterPayload):
            self._render_solution_parameter(payload, display_name)
        elif isinstance(payload, Phase):
            self._render_phase(payload, display_name)
        elif isinstance(payload, Parameter):
            self._render_parameter(payload, display_name)
        elif isinstance(payload, Diagnostic):
            self._render_diagnostic(payload, display_name)
        else:
            self._set_editor_heading(
                display_name or "Unsupported object",
                type(payload).__name__,
            )
            self._clear_form()
            self.editor_form_layout.insertWidget(
                0,
                _text_section("Raw Value", str(payload)),
            )
            self._set_source_preview(None)

    def _render_list(self, payload: list[Any], display_name: str) -> None:
        item_type = type(payload[0]).__name__ if payload else "Object"
        self._set_editor_heading(
            display_name or f"{item_type} list",
            f"{len(payload)} item(s)",
        )
        self._clear_form()
        rows = [
            (index + 1, _object_label(item), _object_subtitle(item))
            for index, item in enumerate(payload)
        ]
        section_title = "Items"
        headers = ("#", "Name", "Details")
        if display_name == "Functions":
            section_title = "Function Expressions"
            headers = ("#", "Name", "Expression")
        elif display_name == "Compounds":
            section_title = "Stoichiometric compounds"
            headers = ("#", "Name", "Phase count")
        elif display_name == "Solutions":
            section_title = "Solution phases"
            headers = ("#", "Name", "Solution model")
        self.editor_form_layout.insertWidget(
            0,
            _table_section(section_title, headers, rows),
        )
        self._set_source_preview(None)

    def _render_element(self, element: Element, display_name: str) -> None:
        self._set_editor_heading(
            display_name or element.symbol,
            element.reference_state,
        )
        self._clear_form()
        section = _form_section("Element Properties")
        form = _form_layout()
        form.addRow("Symbol:", _read_only_line(element.symbol))
        form.addRow("Atomic mass:", _read_only_line(element.atomic_mass))
        form.addRow("Reference state:", _read_only_line(element.reference_state))
        section.layout().addLayout(form)
        self.editor_form_layout.insertWidget(0, section)
        self._set_source_preview(element.source)

    def _render_species(self, species: Species, display_name: str) -> None:
        self._set_editor_heading(
            display_name or species.name,
            "Elemental composition",
        )
        self._clear_form()
        section = _form_section("Species Properties")
        form = _form_layout()
        form.addRow("Name:", _read_only_line(species.name))
        form.addRow("Charge:", _read_only_line(species.charge))
        section.layout().addLayout(form)
        self.editor_form_layout.insertWidget(0, section)
        self.editor_form_layout.insertWidget(
            1,
            _table_section(
                "Composition",
                ("Element", "Amount"),
                list(species.composition.items()) or [("", "")],
            ),
        )
        self._set_source_preview(species.source)

    def _render_function(self, function: FunctionDefinition, display_name: str) -> None:
        self._set_editor_heading(
            display_name or function.name,
            "Stoichiometric Gibbs function",
        )
        self._set_function_title_editor(function)
        self._clear_form()

        self.editor_form_layout.insertWidget(
            0,
            _function_overview_section(
                self.current_database,
                function,
                lambda editor, selected_function=function: (
                    self._validate_function_pertinent_expression(
                        selected_function,
                        editor,
                    )
                ),
            ),
        )
        self._set_source_preview(function.source)

    def _set_function_title_editor(self, function: FunctionDefinition) -> None:
        """Allow the database editor header title to edit a Function name."""
        self._editable_editor_title_record = function
        self.editor_title.setReadOnly(False)
        self.editor_title.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.editor_title.setToolTip("Edit function name")

    def _set_solution_title_editor(self, phase: Phase) -> None:
        """Allow the database editor header title to edit a solution phase name."""
        self._editable_editor_title_record = phase
        self.editor_title.setReadOnly(False)
        self.editor_title.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.editor_title.setToolTip("Edit solution phase name")

    def _commit_function_name(self) -> None:
        """Update the selected editable database object name from the header."""
        record = self._editable_editor_title_record
        if record is None:
            return
        if isinstance(record, Phase):
            self._commit_solution_phase_name(record)
            return
        if not isinstance(record, FunctionDefinition):
            return

        updated_name = self.editor_title.text().strip()
        if not updated_name:
            self.editor_title.setText(record.name)
            return

        old_name = record.name
        record.name = updated_name
        self.editor_title.setText(updated_name)
        self._set_function_name_warnings(
            record,
            _function_name_warning_messages(record),
        )
        self._set_function_name_errors(
            record,
            _function_name_error_messages(record),
        )
        current_index = self.tree.currentIndex()
        is_current_function = (
            current_index.isValid()
            and self.tree_model.payload(current_index) is record
        )
        if is_current_function:
            node = current_index.internalPointer()
            if hasattr(node, "label"):
                node.label = updated_name
            self.tree_model.dataChanged.emit(current_index, current_index)
        for table in self.editor_form.findChildren(QTableWidget, "FormTable"):
            if table.columnCount() < 2:
                continue
            header = table.horizontalHeaderItem(0)
            factor_header = table.horizontalHeaderItem(1)
            if (
                header is None
                or factor_header is None
                or header.text() != "Name"
                or factor_header.text() != "Factor"
            ):
                continue
            name_item = table.item(0, 0)
            if name_item is not None:
                if name_item.text() == old_name:
                    name_item.setText(updated_name)
        if is_current_function:
            self._render_function(record, updated_name)

    def _commit_solution_phase_name(self, phase: Phase) -> None:
        """Update a solution phase name and keep attached parameters in sync."""
        updated_name = self.editor_title.text().strip()
        if not updated_name:
            self.editor_title.setText(phase.name)
            return

        database = self.current_database
        duplicate = next(
            (
                candidate
                for candidate in database.phases
                if candidate is not phase
                and str(candidate.name).strip().upper() == updated_name.upper()
            ),
            None,
        )
        if duplicate is not None:
            QMessageBox.warning(
                self,
                "Rename Solution Phase",
                f"A phase named {updated_name} already exists.",
            )
            self.editor_title.setText(phase.name)
            return

        old_name = phase.name
        if updated_name == old_name:
            return
        phase.name = updated_name
        for parameter in database.parameters:
            if str(parameter.phase_name).strip().upper() == old_name.upper():
                parameter.phase_name = updated_name
        self.editor_title.setText(updated_name)

        current_index = self.tree.currentIndex()
        is_current_phase = (
            current_index.isValid()
            and self.tree_model.payload(current_index) is phase
        )
        if is_current_phase:
            node = current_index.internalPointer()
            if hasattr(node, "label"):
                node.label = updated_name
            self.tree_model.dataChanged.emit(current_index, current_index)
            self._render_phase(phase, updated_name)

    def _validate_function_pertinent_expression(
        self,
        function: FunctionDefinition,
        editor: QLineEdit,
    ) -> None:
        """Validate and save a Function pertinent-expression edit."""
        normalized = self._normalized_pertinent_expression(editor.text())
        if normalized is None:
            return
        wrapped_expression = _wrapped_single_range_expression(normalized)
        candidate = FunctionDefinition(
            name=function.name,
            expression=wrapped_expression,
        )
        _h298, _s298, _cp_rows, gibbs_rows = _runtime_function_thermo_data(
            candidate,
            self.current_database,
        )
        if not gibbs_rows:
            QMessageBox.warning(
                self,
                "Invalid Pertinent Function Expression",
                "The expression could not be converted to Gibbs energy ranges.",
            )
            self._set_status_bar_message(
                f"Function {function.name}: invalid pertinent expression.",
                "error",
            )
            return
        _set_function_gibbs_rows(function, gibbs_rows)
        function.tdb_src.source_expression = normalized
        function.tdb_src.pertinent_functions = [
            TdbFunctionTerm(referenced_function.name, factor)
            for referenced_function, factor in _function_dependency_terms(
                self.current_database,
                FunctionDefinition(function.name, normalized),
            )
        ]
        function.dirty = True
        editor.setText(normalized)
        self._set_status_bar_message(
            f"Function {function.name}: pertinent expression validated.",
            "ready",
        )
        self._rebuild_database_tree(self.current_database, select_payload=function)

    def _validate_compound_phase_pertinent_expression(
        self,
        payload: CompoundPhasePayload,
        editor: QLineEdit,
    ) -> None:
        """Validate and save a compound-phase pertinent-expression edit."""
        normalized = self._normalized_pertinent_expression(editor.text())
        if normalized is None:
            return
        parameter = _compound_phase_primary_g_parameter(
            self.current_database,
            payload.phase,
        )
        if parameter is None:
            parameter = _new_compound_g_parameter(payload.phase)
            self.current_database.parameters.append(parameter)
        wrapped_expression = _wrapped_single_range_expression(normalized)
        candidate = FunctionDefinition(
            name=f"{payload.compound_name}_{payload.phase_label}",
            expression=wrapped_expression,
        )
        _h298, _s298, _cp_rows, gibbs_rows = _runtime_function_thermo_data(
            candidate,
            self.current_database,
        )
        if not gibbs_rows:
            QMessageBox.warning(
                self,
                "Invalid Pertinent Function Expression",
                "The expression could not be converted to Gibbs energy ranges.",
            )
            self._set_status_bar_message(
                (
                    f"Compound {payload.compound_name} phase "
                    f"{payload.phase_label}: invalid pertinent expression."
                ),
                "error",
            )
            return
        parameter.expression = wrapped_expression
        parameter.source.command = _parameter_command_with_expression(
            parameter,
            wrapped_expression,
        )
        parameter.metadata["phase_id"] = _compound_phase_state_id(payload.phase)
        parameter.metadata["source_expression"] = wrapped_expression
        editor.setText(normalized)
        self._set_status_bar_message(
            (
                f"Compound {payload.compound_name} phase "
                f"{payload.phase_label}: pertinent expression validated."
            ),
            "ready",
        )
        self._rebuild_database_tree(self.current_database, show_summary=False)
        self._select_compound_phase(payload.phase)

    def _validate_solution_parameter_expression(
        self,
        payload: SolutionParameterPayload,
        editor: QLineEdit,
    ) -> None:
        """Validate and save a solution endmember/interaction expression edit."""
        normalized = self._normalized_pertinent_expression(editor.text())
        if normalized is None:
            return
        parameter = payload.parameter
        wrapped_expression = _wrapped_single_range_expression(normalized)
        parameter.expression = wrapped_expression
        parameter.metadata["source_expression"] = normalized
        parameter.source.command = _solution_parameter_command_with_expression(
            parameter,
            wrapped_expression,
        )
        editor.setText(normalized)
        self._set_status_bar_message(
            (
                f"Solution {payload.phase.name} {payload.kind.lower()} "
                f"{_solution_parameter_target_text(parameter)}: expression validated."
            ),
            "ready",
        )
        self._set_source_preview(parameter.source)
        if payload.kind == "Endmember":
            self._rebuild_database_tree(self.current_database, show_summary=False)
            self._select_solution_parameter(parameter)

    def _validate_solution_parameter_payload(
        self,
        database: DatabaseIR,
        payload: SolutionParameterPayload,
    ) -> None:
        """Validate and canonicalize a solution parameter from the tree menu."""
        previous_database = self.current_database
        self.current_database = database
        try:
            normalized = self._normalized_pertinent_expression(
                _solution_parameter_edit_text(payload.parameter),
            )
        finally:
            self.current_database = previous_database
        if normalized is None:
            return
        parameter = payload.parameter
        wrapped_expression = _wrapped_single_range_expression(normalized)
        parameter.expression = wrapped_expression
        parameter.metadata["source_expression"] = normalized
        parameter.source.command = _solution_parameter_command_with_expression(
            parameter,
            wrapped_expression,
        )
        self._set_status_bar_message(
            (
                f"Solution {payload.phase.name} {payload.kind.lower()} "
                f"{_solution_parameter_target_text(parameter)}: parameter validated."
            ),
            "ready",
        )
        self._rebuild_database_tree(database, show_summary=False)
        self._select_solution_parameter(parameter)

    def _copy_solution_parameter_expression(
        self,
        payload: SolutionParameterPayload,
    ) -> None:
        """Copy one solution parameter expression for reuse in editors."""
        QApplication.clipboard().setText(_solution_parameter_edit_text(payload.parameter))
        self._set_status_bar_message(
            (
                f"Copied {payload.phase.name} {payload.kind.lower()} "
                f"{_solution_parameter_target_text(payload.parameter)} expression."
            ),
            "ready",
        )

    def _delete_solution_parameter(
        self,
        database: DatabaseIR,
        payload: SolutionParameterPayload,
    ) -> None:
        """Delete one solution endmember or interaction parameter."""
        if payload.parameter not in database.parameters:
            return
        database.parameters.remove(payload.parameter)
        self._set_status_bar_message(
            (
                f"Deleted {payload.phase.name} {payload.kind.lower()} "
                f"{_solution_parameter_target_text(payload.parameter)}."
            ),
            "ready",
        )
        self._rebuild_database_tree(database, select_payload=payload.phase)

    def _solution_endmember_function_for_parameter(
        self,
        phase: Phase,
        parameter: Parameter,
    ) -> FunctionDefinition:
        """Return a transient Function wrapper for a solution endmember."""
        expression = _solution_endmember_source_range(parameter)
        ranges = [
            (float(lower), float(upper))
            for lower, upper, _segment in _function_expression_segments(expression)
        ]
        function = FunctionDefinition(
            name=_solution_parameter_target_text(parameter).replace(":", "_"),
            expression=expression,
            temperature_ranges=ranges,
            source=parameter.source,
        )
        function.tdb_src.source_expression = str(
            parameter.metadata.get("source_expression", "")
        )
        return function

    def _solution_endmember_function_for_range_payload(
        self,
        payload: ThermoRangePayload,
        *,
        database: DatabaseIR | None = None,
    ) -> FunctionDefinition:
        """Return a transient Function wrapper for one endmember Cp range."""
        active_database = database or self.current_database
        parameter = _solution_endmember_parameter_for_range_payload(
            active_database,
            payload,
        )
        if parameter is None:
            ranges = []
            if payload.lower != "" and payload.upper != "":
                ranges = [(float(payload.lower), float(payload.upper))]
            source = (
                payload.source
                if isinstance(payload.source, SourceRef)
                else SourceRef()
            )
            return FunctionDefinition(
                name=payload.owner_name.replace(":", "_"),
                expression=str(payload.expression or ""),
                temperature_ranges=ranges,
                source=source,
            )
        phase = next(
            (
                candidate
                for candidate in active_database.phases
                if candidate.name.upper() == parameter.phase_name.upper()
            ),
            Phase(parameter.phase_name),
        )
        return self._solution_endmember_function_for_parameter(phase, parameter)

    def _add_solution_endmember_cp_range(
        self,
        database: DatabaseIR,
        payload: SolutionParameterPayload,
    ) -> None:
        """Append a Cp/G range to one solution endmember parameter."""
        function = self._solution_endmember_function_for_parameter(
            payload.phase,
            payload.parameter,
        )
        gibbs_rows = _direct_function_gibbs_rows(function)
        if not gibbs_rows:
            _h298, _s298, _cp_rows, gibbs_rows = _runtime_function_thermo_data(
                function,
                database,
            )
        if len(gibbs_rows) >= 6:
            return
        if not gibbs_rows:
            gibbs_rows = [[298.15, 6000.0, *([0.0] * 14)]]
        last = [float(value) for value in gibbs_rows[-1]]
        lower = float(last[1])
        last[0] = lower
        last[1] = lower + 100.0
        gibbs_rows.append(last)
        _set_parameter_gibbs_rows(payload.parameter, gibbs_rows)
        payload.parameter.source.command = _solution_parameter_command_with_expression(
            payload.parameter,
            payload.parameter.expression,
        )
        self._rebuild_database_tree(database, show_summary=False)
        self._select_solution_endmember_cp_range(
            _solution_endmember_owner_name(payload.phase, payload.parameter),
            len(gibbs_rows),
        )

    def _delete_solution_endmember_cp_range(
        self,
        database: DatabaseIR,
        payload: ThermoRangePayload,
    ) -> None:
        """Delete one Cp/G range from a solution endmember parameter."""
        parameter = _solution_endmember_parameter_for_range_payload(database, payload)
        if parameter is None:
            return
        function = self._solution_endmember_function_for_range_payload(
            payload,
            database=database,
        )
        gibbs_rows = _direct_function_gibbs_rows(function)
        if not gibbs_rows:
            _h298, _s298, _cp_rows, gibbs_rows = _runtime_function_thermo_data(
                function,
                database,
            )
        remove_index = max(0, payload.index - 1)
        if remove_index < len(gibbs_rows):
            gibbs_rows.pop(remove_index)
        if gibbs_rows:
            _set_parameter_gibbs_rows(parameter, gibbs_rows)
            parameter.source.command = _solution_parameter_command_with_expression(
                parameter,
                parameter.expression,
            )
        else:
            database.parameters.remove(parameter)
        self._rebuild_database_tree(database, show_summary=False)
        if gibbs_rows:
            self._select_solution_endmember_cp_range(
                payload.owner_name,
                min(remove_index + 1, len(gibbs_rows)),
            )

    def _normalized_pertinent_expression(self, expression: str) -> str | None:
        """Return validated normalized expression text, or show an error."""
        expression_text = _strip_pertinent_function_markers(expression)
        if not expression_text:
            QMessageBox.warning(
                self,
                "Invalid Pertinent Function Expression",
                "Enter a Gibbs energy expression before validating.",
            )
            return None
        if ";" in expression_text:
            QMessageBox.warning(
                self,
                "Invalid Pertinent Function Expression",
                "Enter one expression without Cp temperature ranges.",
            )
            return None
        missing = _missing_pertinent_function_names(
            expression_text,
            self.current_database,
        )
        if missing:
            QMessageBox.warning(
                self,
                "Missing Function",
                (
                    "The expression references Function name(s) that are not "
                    f"available: {', '.join(missing)}"
                ),
            )
            self._set_status_bar_message(
                f"Missing Function: {', '.join(missing)}",
                "error",
            )
            return None
        try:
            normalized = _canonical_tdb_gibbs_expression(expression_text)
        except ValueError as exc:
            QMessageBox.warning(
                self,
                "Invalid Pertinent Function Expression",
                str(exc),
            )
            self._set_status_bar_message(
                "Invalid pertinent expression.",
                "error",
            )
            return None
        return _restore_function_name_case(normalized, self.current_database)

    def _render_solution_group(
        self,
        payload: SolutionGroupPayload,
        display_name: str,
    ) -> None:
        """Render one group node under a solution phase."""
        phase = payload.phase
        group_name = payload.group_name
        self._set_editor_heading(display_name or group_name, phase.name)
        self._clear_form()

        if group_name == "Sublattices":
            rows = _solution_sublattice_rows(phase)
            self.editor_form_layout.insertWidget(
                0,
                _table_section(
                    "Sublattices",
                    ("Label", "Site Ratio", "Species"),
                    rows or [("", "", "")],
                ),
            )
            self._set_source_preview(phase.source)
            return

        parameters = _solution_parameters_for_phase(self.current_database, phase)
        if group_name == "Endmembers":
            rows = [
                (
                    index,
                    _solution_parameter_target_text(parameter),
                    parameter.expression,
                )
                for index, parameter in enumerate(parameters)
                if not _solution_parameter_is_interaction(parameter)
            ]
            self.editor_form_layout.insertWidget(
                0,
                _table_section(
                    "Endmembers",
                    ("#", "Target", "Gibbs Energy Function"),
                    rows or [("", "", "")],
                ),
            )
            self._set_source_preview(phase.source)
            return

        if group_name == "Interactions":
            rows = [
                (
                    index,
                    parameter.parameter_type,
                    _solution_parameter_target_text(parameter),
                    parameter.order,
                    parameter.expression,
                )
                for index, parameter in enumerate(parameters)
                if _solution_parameter_is_interaction(parameter)
            ]
            self.editor_form_layout.insertWidget(
                0,
                _table_section(
                    "Interactions",
                    ("#", "Type", "Target", "Order", "Expression"),
                    rows or [("", "", "", "", "")],
                ),
            )
            self._set_source_preview(phase.source)
            return

        self.editor_form_layout.insertWidget(
            0,
            _text_section("Solution Group", group_name),
        )
        self._set_source_preview(phase.source)

    def _render_solution_sublattice(
        self,
        payload: SolutionSublatticePayload,
        display_name: str,
    ) -> None:
        """Render a solution sublattice and its species list."""
        self._set_editor_heading(
            display_name or payload.label,
            f"{payload.phase.name} sublattice",
        )
        self._clear_form()

        identity = _form_section("Sublattice")
        form = _form_layout()
        form.addRow("Label:", _read_only_line(payload.label))
        form.addRow("Index:", _read_only_line(payload.sublattice))
        form.addRow("Site Ratio:", _read_only_line(payload.site_ratio))
        form.addRow("Species Count:", _read_only_line(len(payload.species)))
        identity.layout().addLayout(form)
        self.editor_form_layout.insertWidget(0, identity)

        rows = [
            (species_name, _solution_species_formula(self.current_database, species_name))
            for species_name in payload.species
        ]
        self.editor_form_layout.insertWidget(
            1,
            _table_section(
                "Species",
                ("Name", "Formula"),
                rows or [("", "")],
            ),
        )
        self._set_source_preview(payload.phase.source)

    def _render_solution_species(
        self,
        payload: SolutionSpeciesPayload,
        display_name: str,
    ) -> None:
        """Render one species inside a solution sublattice."""
        species = _solution_species_record(self.current_database, payload.species_name)
        self._set_editor_heading(
            display_name or payload.species_name,
            f"{payload.phase.name} {payload.sublattice_label}",
        )
        self._clear_form()

        section = _form_section("Species")
        form = _form_layout()
        form.addRow("Species Name:", _read_only_line(payload.species_name))
        form.addRow("Sublattice:", _read_only_line(payload.sublattice_label))
        form.addRow("Site Ratio:", _read_only_line(payload.site_ratio))
        form.addRow(
            "Valence:",
            _read_only_line("" if species is None else _compact_gui_number(species.charge)),
        )
        form.addRow("Chem. Group:", _read_only_line(""))
        form.addRow(
            "Formula:",
            _read_only_line(_solution_species_formula(self.current_database, payload.species_name)),
        )
        section.layout().addLayout(form)
        self.editor_form_layout.insertWidget(0, section)
        self._set_source_preview(species.source if species is not None else payload.phase.source)

    def _render_solution_parameter(
        self,
        payload: SolutionParameterPayload,
        display_name: str,
    ) -> None:
        """Render a solution endmember or interaction parameter."""
        parameter = payload.parameter
        target = _solution_parameter_target_text(parameter)
        if payload.kind == "Endmember":
            self._render_solution_endmember_parameter(payload, display_name)
            return
        self._set_editor_heading(
            display_name or target,
            f"{payload.phase.name} {payload.kind.lower()}",
        )
        self._clear_form()

        section = _form_section(payload.kind)
        form = _form_layout()
        form.addRow("Phase:", _read_only_line(payload.phase.name))
        form.addRow("Type:", _read_only_line(parameter.parameter_type))
        form.addRow("Target:", _read_only_line(target))
        form.addRow("Order:", _read_only_line(parameter.order))
        form.addRow("Status:", _read_only_line(parameter.metadata.get("status", "")))
        section.layout().addLayout(form)
        self.editor_form_layout.insertWidget(0, section)

        expression_editor = _FunctionDropLineEdit(
            _solution_parameter_edit_text(parameter),
        )
        expression_editor.setObjectName("SolutionParameterExpressionInput")
        expression_editor.setToolTip("Drop a Function here to insert FUNCTION#.")
        validate_button = QPushButton("Validate")
        validate_button.setObjectName("SolutionParameterValidateButton")
        validate_button.clicked.connect(
            lambda _checked=False,
            selected_payload=payload,
            editor=expression_editor: self._validate_solution_parameter_expression(
                selected_payload,
                editor,
            )
        )
        expression_section = _form_section_with_actions(
            "Expression",
            [validate_button],
        )
        expression_section.layout().addWidget(expression_editor)
        self.editor_form_layout.insertWidget(1, expression_section)
        self._set_source_preview(parameter.source)

    def _render_solution_endmember_parameter(
        self,
        payload: SolutionParameterPayload,
        display_name: str,
    ) -> None:
        """Render a solution endmember as a stoichiometric Gibbs/Cp record."""
        parameter = payload.parameter
        target = _solution_parameter_target_text(parameter)
        function = self._solution_endmember_function_for_parameter(
            payload.phase,
            parameter,
        )
        self._set_editor_heading(
            display_name or target,
            f"{payload.phase.name} endmember",
        )
        self._clear_form()

        self.editor_form_layout.insertWidget(
            0,
            _function_overview_section(
                self.current_database,
                function,
                lambda editor, selected_payload=payload: (
                    self._validate_solution_parameter_expression(
                        selected_payload,
                        editor,
                    )
                ),
                title="Endmember Overview",
            ),
        )
        self._set_source_preview(parameter.source)

    def _solution_properties_section(self, phase: Phase) -> QFrame:
        """Return the editable top-row solution properties section."""
        add_button = QPushButton("Add Sublattice")
        add_button.setObjectName("SmallActionButton")
        add_button.setEnabled(
            _solution_model_editor_value(phase) == "CEF"
            and len(phase.constituents) < _MAX_SOLUTION_SUBLATTICES
        )
        add_button.clicked.connect(
            lambda _checked=False, selected_phase=phase: (
                self._add_solution_sublattice(selected_phase)
            )
        )
        remove_button = QPushButton("Remove Sublattice")
        remove_button.setObjectName("SmallActionButton")
        remove_button.setEnabled(
            _solution_model_editor_value(phase) == "CEF"
            and len(phase.constituents) > 1
        )
        remove_button.clicked.connect(
            lambda _checked=False, selected_phase=phase: (
                self._remove_solution_sublattice(selected_phase)
            )
        )
        populate_button = QPushButton("Populate Endmembers")
        populate_button.setObjectName("PrimaryButton")
        populate_button.clicked.connect(
            lambda _checked=False, selected_phase=phase: (
                self._populate_solution_endmembers(selected_phase)
            )
        )

        section = _form_section_with_actions(
            "Solution Properties",
            [add_button, remove_button, populate_button],
        )
        row = QHBoxLayout()
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(12)

        model_label = QLabel("Solution model:")
        model_combo = QComboBox()
        model_combo.setObjectName("SolutionModelCombo")
        model_combo.addItems(_SOLUTION_MODEL_CHOICES)
        model_combo.setCurrentText(_solution_model_editor_value(phase))

        lattice_label = QLabel("Lattices:")
        lattice_input = QLineEdit(str(len(phase.constituents)))
        lattice_input.setObjectName("SolutionLatticesInput")
        lattice_input.setReadOnly(True)
        lattice_input.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        lattice_input.setFixedWidth(72)

        model_combo.currentTextChanged.connect(
            lambda value, selected_phase=phase: (
                self._update_solution_model(selected_phase, value)
            )
        )

        row.addWidget(model_label)
        row.addWidget(model_combo, 1)
        row.addSpacing(24)
        row.addWidget(lattice_label)
        row.addWidget(lattice_input)
        row.addStretch(2)
        section.layout().addLayout(row)
        section.layout().addWidget(self._solution_sublattice_table(phase))
        return section

    def _update_solution_model(
        self,
        phase: Phase,
        model_name: str,
    ) -> None:
        """Update the solution model from the overview dropdown."""
        if model_name == "CEF":
            phase.model = "CEF"
        else:
            phase.model = "SOLUTION"
            if len(phase.constituents) != 1:
                self._set_solution_lattice_count(phase, 1)
        self._refresh_solution_phase_editor(phase)

    def _solution_sublattice_table(self, phase: Phase) -> QTableWidget:
        """Return the editable sublattice table for a solution phase."""
        table = QTableWidget(0, 3)
        table.setObjectName("SolutionSublatticeTable")
        table.setHorizontalHeaderLabels(
            ["Sublattice", "Stoichiometric coefficients", "Constituents"]
        )
        table.verticalHeader().setVisible(False)
        table.verticalHeader().setDefaultSectionSize(62)
        table.verticalHeader().setMinimumSectionSize(62)
        header = table.horizontalHeader()
        header.setStretchLastSection(True)
        header.setDefaultAlignment(Qt.AlignmentFlag.AlignCenter)
        header.setMinimumSectionSize(130)
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        table.setShowGrid(False)
        table.setAlternatingRowColors(False)
        table.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        table.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)

        for row, constituent_set in enumerate(phase.constituents):
            table.insertRow(row)
            table.setRowHeight(row, 62)

            name_editor = QLineEdit(_solution_sublattice_name(phase, row + 1))
            name_editor.setObjectName("SolutionSublatticeNameInput")
            name_editor.setPlaceholderText("A")
            name_editor.setMinimumHeight(38)
            name_editor.editingFinished.connect(
                lambda r=row, editor=name_editor, selected_phase=phase: (
                    self._commit_solution_sublattice_name(
                        selected_phase,
                        r,
                        editor.text(),
                    )
                )
            )
            table.setCellWidget(row, 0, name_editor)

            coefficient_editor = QLineEdit(
                _compact_gui_number(constituent_set.site_ratio)
            )
            coefficient_editor.setObjectName("SolutionSublatticeCoefficientInput")
            coefficient_editor.setPlaceholderText("e.g. 1")
            coefficient_editor.setMinimumHeight(38)
            coefficient_editor.editingFinished.connect(
                lambda r=row, editor=coefficient_editor, selected_phase=phase: (
                    self._commit_solution_sublattice_coefficient(
                        selected_phase,
                        r,
                        editor.text(),
                    )
                )
            )
            table.setCellWidget(row, 1, coefficient_editor)

            constituents_editor = QLineEdit(
                _solution_constituents_display(
                    self.current_database,
                    constituent_set,
                )
            )
            constituents_editor.setObjectName("SolutionSublatticeConstituentsInput")
            constituents_editor.setPlaceholderText("e.g. Al [3+], O [-2]")
            constituents_editor.setMinimumHeight(38)
            constituents_editor.setToolTip(
                "Separate constituents with commas. Charges can be written as Al [3+] or Al [+3]."
            )
            constituents_editor.editingFinished.connect(
                lambda r=row, editor=constituents_editor, selected_phase=phase: (
                    self._commit_solution_sublattice_constituents(
                        selected_phase,
                        r,
                        editor.text(),
                    )
                )
            )
            table.setCellWidget(row, 2, constituents_editor)

        _fit_table_columns_to_headers(table)
        _fit_table_to_rows(table, min_rows=max(4, len(phase.constituents)))
        return table

    def _commit_solution_sublattice_name(
        self,
        phase: Phase,
        row: int,
        name: str,
    ) -> None:
        """Commit an edited sublattice name and refresh sidebar labels."""
        if not 0 <= row < len(phase.constituents):
            return
        _set_solution_sublattice_name(phase, row + 1, name)
        self._refresh_solution_phase_editor(phase)

    def _commit_solution_sublattice_coefficient(
        self,
        phase: Phase,
        row: int,
        coefficient: str,
    ) -> None:
        """Commit an edited sublattice stoichiometric coefficient."""
        if not 0 <= row < len(phase.constituents):
            return
        try:
            phase.constituents[row].site_ratio = float(str(coefficient).strip())
        except ValueError:
            self._refresh_solution_phase_editor(phase)
            return
        self._refresh_solution_phase_editor(phase)

    def _commit_solution_sublattice_constituents(
        self,
        phase: Phase,
        row: int,
        constituents: str,
    ) -> None:
        """Commit comma-separated constituents, including optional charges."""
        if not 0 <= row < len(phase.constituents):
            return
        species_names: list[str] = []
        for token in str(constituents or "").split(","):
            parsed = _parse_solution_constituent_token(token)
            if parsed is None:
                continue
            species_name, charge = parsed
            _upsert_solution_species(self.current_database, species_name, charge)
            species_names.append(species_name)
        phase.constituents[row].species = species_names
        self._refresh_solution_phase_editor(phase)

    def _set_solution_lattice_count(self, phase: Phase, lattice_count: int) -> None:
        """Resize the solution sublattice list."""
        target_count = max(1, min(_MAX_SOLUTION_SUBLATTICES, int(lattice_count)))
        current_count = len(phase.constituents)
        if target_count == current_count:
            return
        if target_count > current_count:
            for sublattice in range(current_count + 1, target_count + 1):
                phase.constituents.append(ConstituentSet(sublattice, [], 1.0))
        else:
            del phase.constituents[target_count:]
        for index, constituent_set in enumerate(phase.constituents, start=1):
            constituent_set.sublattice = index
        _trim_solution_sublattice_labels(phase)

    def _add_solution_sublattice(self, phase: Phase) -> None:
        """Add one CEF sublattice row, up to the supported GUI limit."""
        if _solution_model_editor_value(phase) != "CEF":
            return
        if len(phase.constituents) >= _MAX_SOLUTION_SUBLATTICES:
            return
        self._set_solution_lattice_count(phase, len(phase.constituents) + 1)
        self._refresh_solution_phase_editor(phase)

    def _remove_solution_sublattice(self, phase: Phase) -> None:
        """Remove the selected sublattice row, or the last row if none selected."""
        if _solution_model_editor_value(phase) != "CEF" or len(phase.constituents) <= 1:
            return
        table = self.editor_form.findChild(QTableWidget, "SolutionSublatticeTable")
        remove_index = (
            table.currentRow()
            if table is not None and table.currentRow() >= 0
            else len(phase.constituents) - 1
        )
        remove_index = max(0, min(remove_index, len(phase.constituents) - 1))
        del phase.constituents[remove_index]
        labels = phase.metadata.get("sublattice_labels")
        if isinstance(labels, list) and remove_index < len(labels):
            del labels[remove_index]
        for index, constituent_set in enumerate(phase.constituents, start=1):
            constituent_set.sublattice = index
        _trim_solution_sublattice_labels(phase)
        self._refresh_solution_phase_editor(phase)

    def _populate_solution_endmembers(self, phase: Phase) -> None:
        """Synchronize endmember parameters implied by the sublattice table."""
        expected_endmembers = _solution_expected_endmembers(phase)
        if not expected_endmembers:
            QMessageBox.warning(
                self,
                "Populate Endmembers",
                "Enter at least one constituent in each sublattice first.",
            )
            return

        phase_parameters = _solution_parameters_for_phase(
            self.current_database,
            phase,
        )
        existing_endmembers = [
            parameter
            for parameter in phase_parameters
            if _solution_endmember_parameter(parameter)
        ]
        existing_by_key = {
            _solution_endmember_match_key(
                self.current_database,
                phase,
                tuple(_solution_parameter_sections(parameter)),
            ): parameter
            for parameter in existing_endmembers
        }
        expected_keys = {
            _solution_endmember_match_key(
                self.current_database,
                phase,
                endmember,
            )
            for endmember in expected_endmembers
        }
        existing_keys = set(existing_by_key)
        preserve_existing = (
            bool(existing_endmembers)
            and existing_keys == expected_keys
            and len(existing_endmembers) == len(expected_endmembers)
        )
        if existing_endmembers and not preserve_existing:
            QMessageBox.warning(
                self,
                "Populate Endmembers",
                (
                    "Existing endmembers do not match the current sublattice "
                    "stoichiometries. They will be replaced with freshly "
                    "generated endmembers."
                ),
            )

        new_endmembers: list[Parameter] = []
        for endmember in expected_endmembers:
            target_text = ":".join(endmember)
            key = _solution_endmember_match_key(
                self.current_database,
                phase,
                endmember,
            )
            if preserve_existing and key in existing_by_key:
                parameter = copy.deepcopy(existing_by_key[key])
                parameter.phase_name = phase.name
                parameter.parameter_type = "G"
                parameter.target = [target_text]
                parameter.order = 0
                expression = _solution_endmember_source_range(parameter)
            else:
                expression = "298.15 ZERO#; 6000 N"
                parameter = Parameter(
                    phase.name,
                    "G",
                    [target_text],
                    0,
                    expression,
                )
                parameter.metadata["source_expression"] = "ZERO#"
            parameter.source.command = (
                f"PAR G({phase.name},{target_text};0) {expression} !"
            )
            new_endmembers.append(parameter)

        old_parameters = list(self.current_database.parameters)
        phase_endmember_indexes = [
            index
            for index, parameter in enumerate(old_parameters)
            if parameter.phase_name.upper() == phase.name.upper()
            and _solution_endmember_parameter(parameter)
        ]
        phase_parameter_indexes = [
            index
            for index, parameter in enumerate(old_parameters)
            if parameter.phase_name.upper() == phase.name.upper()
        ]
        insert_index = (
            phase_endmember_indexes[0]
            if phase_endmember_indexes
            else phase_parameter_indexes[0]
            if phase_parameter_indexes
            else len(old_parameters)
        )
        self.current_database.parameters = [
            parameter
            for parameter in old_parameters
            if not (
                parameter.phase_name.upper() == phase.name.upper()
                and _solution_endmember_parameter(parameter)
            )
        ]
        for offset, parameter in enumerate(new_endmembers):
            self.current_database.parameters.insert(insert_index + offset, parameter)

        self._set_status_bar_message(
            f"Populated {len(new_endmembers)} endmember parameter(s).",
            "ready",
        )
        self._refresh_solution_phase_editor(phase)

    def _refresh_solution_phase_editor(self, phase: Phase) -> None:
        """Refresh the solution tree labels and current editor after an edit."""
        current_index = self.tree.currentIndex()
        if current_index.isValid() and self.tree_model.payload(current_index) is phase:
            self._rebuild_database_tree(
                self.current_database,
                select_payload=phase,
                show_summary=False,
            )
            return
        self._render_phase(phase, phase.name)

    def _render_phase(self, phase: Phase, display_name: str) -> None:
        if not _is_solution_phase(phase):
            self._render_compound_phase_payload(
                CompoundPhasePayload(
                    compound_name=display_name or phase.name,
                    phase_label="S1",
                    phase=phase,
                    index=1,
                ),
                "S1",
            )
            return

        self._set_editor_heading(
            display_name or phase.name,
            "Solution Phase",
        )
        self._set_solution_title_editor(phase)
        self._clear_form()
        self.editor_form_layout.addWidget(self._solution_properties_section(phase))
        magnetic_section = _phase_magnetic_properties_section(
            self.current_database,
            phase,
        )
        if magnetic_section is not None:
            self.editor_form_layout.addWidget(magnetic_section)
        self._set_source_preview(phase.source)

    def _render_compound_record(
        self,
        record: CompoundRecord,
        display_name: str,
    ) -> None:
        """Render one compound parent record with its phase-state list."""
        for index, phase in enumerate(record.phases[:10], start=1):
            self._refresh_compound_phase_warning_cache(
                self.current_database,
                phase,
                f"{record.name}:{_compound_phase_display_label(phase, index)}",
                emit=False,
            )
        self._refresh_tree_model_range_warnings()
        self._set_editor_heading(
            display_name or record.name,
            "Stoichiometric Compound",
        )
        self._clear_form()
        rows = [
            (
                _compound_phase_display_label(phase, index),
                "Berman&Brown",
                _phase_formula(phase),
            )
            for index, phase in enumerate(record.phases, start=1)
        ]
        self.editor_form_layout.insertWidget(
            0,
            _table_section(
                "Phases",
                ("Phase", "Model", "Formula"),
                rows or [("", "", "")],
            ),
        )
        self.editor_form_layout.insertWidget(
            1,
            _table_section(
                "Parameters",
                ("Phase", "Parameter", "Low T [K]", "High T [K]", "Expression"),
                _compound_record_parameter_rows(self.current_database, record),
            ),
        )
        source = record.phases[0].source if record.phases else SourceRef()
        self._set_source_preview(source)

    def _render_compound_phase_payload(
        self,
        payload: CompoundPhasePayload,
        display_name: str,
    ) -> None:
        """Render the editable overview for one compound phase state."""
        phase = payload.phase
        phase_label = display_name or payload.phase_label
        self._refresh_compound_phase_warning_cache(
            self.current_database,
            phase,
            f"{payload.compound_name}:{payload.phase_label}",
        )
        self._set_editor_heading(payload.compound_name, "")
        header_widget = _compound_phase_header_widget(phase_label)
        phase_label_editor = header_widget.findChild(
            QLineEdit,
            "CompoundPhaseNameInlineInput",
        )
        if phase_label_editor is not None:
            phase_label_editor.editingFinished.connect(
                lambda phase_payload=payload, editor=phase_label_editor: (
                    self._commit_compound_phase_label(phase_payload, editor)
                )
            )
        self._set_editor_subtitle_widget(header_widget)
        self._clear_form()
        self.editor_form_layout.insertWidget(
            0,
            _compound_phase_overview_section(
                self.current_database,
                phase,
                payload.phase_label,
            ),
        )
        insert_index = 1
        magnetic_section = _phase_magnetic_properties_section(
            self.current_database,
            phase,
        )
        if magnetic_section is not None:
            self.editor_form_layout.insertWidget(insert_index, magnetic_section)
            insert_index += 1
        self.editor_form_layout.insertWidget(
            insert_index,
            _compound_phase_parameters_section(
                self.current_database,
                phase,
                lambda editor, selected_payload=payload: (
                    self._validate_compound_phase_pertinent_expression(
                        selected_payload,
                        editor,
                    )
                ),
            ),
        )
        self._set_source_preview(phase.source)

    def _render_thermo_record(
        self,
        *,
        record_type: str,
        name: str,
        formula: str,
        model: str,
        cp_range_rows: list[tuple[Any, ...]],
        building_block_rows: list[tuple[Any, ...]],
        source: SourceRef,
    ) -> None:
        self._set_editor_heading(name, f"{record_type} thermodynamic record")
        self._clear_form()

        identity = _form_section("Identity")
        identity_form = _form_layout()
        identity_form.addRow("Name:", _read_only_line(name))
        identity_form.addRow("Formula:", _read_only_line(formula))
        identity_form.addRow("Type:", _read_only_line(record_type))
        identity_form.addRow("Model:", _read_only_line(model))
        identity.layout().addLayout(identity_form)
        self.editor_form_layout.insertWidget(0, identity)

        reference = _form_section("Reference State")
        reference_form = _form_layout()
        reference_form.addRow("Delta H298:", _read_only_line(""))
        reference_form.addRow("S298:", _read_only_line(""))
        reference_form.addRow("Density:", _read_only_line(""))
        reference_form.addRow("References:", _read_only_line(""))
        reference.layout().addLayout(reference_form)
        self.editor_form_layout.insertWidget(1, reference)

        self.editor_form_layout.insertWidget(
            2,
            _table_section(
                "Cp Ranges",
                ("Range", "Lower", "Upper"),
                cp_range_rows or [("", "", "")],
            ),
        )
        self.editor_form_layout.insertWidget(
            3,
            _table_section(
                "Building Blocks",
                ("Name", "Role", "Expression"),
                building_block_rows or [("", "", "")],
            ),
        )
        self._set_source_preview(source)

    def _render_thermo_range(
        self,
        payload: ThermoRangePayload,
        display_name: str,
    ) -> None:
        if (
            payload.owner_kind == "Function"
            and payload.property_name == "Cp"
        ):
            self._render_function_thermo_range(payload, display_name)
            return
        if (
            payload.owner_kind == "Compound"
            and payload.property_name == "Cp"
        ):
            self._render_compound_thermo_range(payload, display_name)
            return
        if (
            payload.owner_kind == "SolutionEndmember"
            and payload.property_name == "Cp"
        ):
            self._render_solution_endmember_thermo_range(payload, display_name)
            return

        title = display_name or f"{payload.property_name} {payload.index}"
        self._set_editor_heading(
            title,
            f"{payload.owner_kind} range for {payload.owner_name}",
        )
        self._clear_form()

        range_section = _form_section("Cp Range")
        range_form = _form_layout()
        range_form.addRow("Owner:", _read_only_line(payload.owner_name))
        range_form.addRow("Type:", _read_only_line(payload.owner_kind))
        range_form.addRow("Property:", _read_only_line(payload.property_name))
        range_form.addRow("Lower T:", _read_only_line(payload.lower))
        range_form.addRow("Upper T:", _read_only_line(payload.upper))
        range_section.layout().addLayout(range_form)
        self.editor_form_layout.insertWidget(0, range_section)

        self.editor_form_layout.insertWidget(
            1,
            _thermo_expression_tab_section(
                {payload.property_name: [(title, payload.expression)]}
            ),
        )
        if isinstance(payload.source, SourceRef):
            self._set_source_preview(payload.source)
        else:
            self._set_source_preview(None)

    def _render_compound_thermo_range(
        self,
        payload: ThermoRangePayload,
        display_name: str,
    ) -> None:
        """Render a compound Gibbs parameter through the Cp/G/H/S editor."""
        function = self._compound_function_for_range_payload(payload)
        title = display_name or f"Cp {payload.index}"
        range_index = self._compound_range_index_for_payload(function, payload)
        editor = _FunctionThermoEditor(
            function,
            database=self.current_database,
            range_index=range_index,
            continuity_changed_callback=(
                lambda edited_function,
                _warnings,
                _selected_warnings=None,
                _gibbs_rows=None,
                range_payload=payload: (
                    self._set_compound_range_continuity_warnings(
                        range_payload,
                        edited_function,
                        gibbs_rows=_gibbs_rows,
                    )
                )
            ),
            gibbs_changed_callback=(
                lambda _function, gibbs_rows, range_payload=payload: (
                    self._compound_thermo_range_changed(
                        range_payload,
                        gibbs_rows,
                    )
                )
            ),
        )
        self._set_editor_heading(title, "")
        self._set_editor_subtitle_widget(editor.range_widget)
        self._clear_form()
        self._active_database_thermo_editor = editor
        self.editor_form_layout.insertWidget(
            0,
            _function_thermo_expression_section(editor),
        )
        if isinstance(payload.source, SourceRef):
            self._set_source_preview(payload.source)
        else:
            self._set_source_preview(None)

    def _render_solution_endmember_thermo_range(
        self,
        payload: ThermoRangePayload,
        display_name: str,
    ) -> None:
        """Render a solution endmember through the Cp/G/H/S editor."""
        function = self._solution_endmember_function_for_range_payload(payload)
        title = display_name or f"Cp {payload.index}"
        range_index = self._solution_endmember_range_index_for_payload(
            function,
            payload,
        )
        editor = _FunctionThermoEditor(
            function,
            database=self.current_database,
            range_index=range_index,
            gibbs_changed_callback=(
                lambda _function, gibbs_rows, range_payload=payload: (
                    self._solution_endmember_thermo_range_changed(
                        range_payload,
                        gibbs_rows,
                    )
                )
            ),
        )
        self._set_editor_heading(title, "")
        self._set_editor_subtitle_widget(editor.range_widget)
        self._clear_form()
        self._active_database_thermo_editor = editor
        self.editor_form_layout.insertWidget(
            0,
            _function_thermo_expression_section(editor),
        )
        if isinstance(payload.source, SourceRef):
            self._set_source_preview(payload.source)
        else:
            self._set_source_preview(None)

    def _render_function_thermo_range(
        self,
        payload: ThermoRangePayload,
        display_name: str,
    ) -> None:
        """Render the tabbed H/S/Cp/G editor for a Function Cp child."""
        function = self._function_for_range_payload(payload)
        title = display_name or f"{payload.property_name} {payload.index}"
        editor = _FunctionThermoEditor(
            function,
            database=self.current_database,
            range_index=payload.index - 1,
            range_changed_callback=self._function_cp_range_changed,
            continuity_changed_callback=self._set_function_continuity_warnings,
        )
        self._set_editor_heading(title, "")
        self._set_editor_subtitle_widget(editor.range_widget)
        self._clear_form()
        self._active_database_thermo_editor = editor
        self.editor_form_layout.insertWidget(
            0,
            _function_thermo_expression_section(editor),
        )
        if isinstance(payload.source, SourceRef):
            self._set_source_preview(payload.source)
        else:
            self._set_source_preview(None)

    def _compound_function_for_range_payload(
        self,
        payload: ThermoRangePayload,
        *,
        database: DatabaseIR | None = None,
    ) -> FunctionDefinition:
        """Return a transient Function wrapper for one compound parameter."""
        active_database = database or self.current_database
        parameter = _compound_parameter_for_range_payload(active_database, payload)
        expression = (
            _parameter_as_function_expression(parameter)
            if parameter is not None
            else str(payload.expression).strip()
        )
        if not _function_expression_segments(expression):
            lower_value = payload.lower if payload.lower != "" else 298.15
            upper_value = payload.upper if payload.upper != "" else 6000
            expression = f"{lower_value} {expression}; {upper_value} N"
        ranges = [
            (float(lower), float(upper))
            for lower, upper, _segment in _function_expression_segments(expression)
        ]
        source = parameter.source if parameter is not None else payload.source
        if not isinstance(source, SourceRef):
            source = SourceRef()
        return FunctionDefinition(
            name=payload.owner_name.replace(":", "_"),
            expression=expression,
            temperature_ranges=ranges,
            source=source,
        )

    def _compound_range_index_for_payload(
        self,
        function: FunctionDefinition,
        payload: ThermoRangePayload,
        *,
        database: DatabaseIR | None = None,
    ) -> int:
        """Return the thermo editor segment represented by a compound tree child."""
        fallback = max(0, payload.index - 1)
        try:
            payload_upper = float(payload.upper)
        except (TypeError, ValueError):
            return fallback
        _h298, _s298, _cp_rows, gibbs_rows = _runtime_function_thermo_data(
            function,
            database or self.current_database,
        )
        for index, row in enumerate(gibbs_rows):
            if np.isclose(float(row[1]), payload_upper, rtol=1e-10, atol=1e-8):
                return index
        return fallback

    def _solution_endmember_range_index_for_payload(
        self,
        function: FunctionDefinition,
        payload: ThermoRangePayload,
        *,
        database: DatabaseIR | None = None,
    ) -> int:
        """Return the thermo editor segment represented by an endmember child."""
        fallback = max(0, payload.index - 1)
        try:
            payload_upper = float(payload.upper)
        except (TypeError, ValueError):
            return fallback
        _h298, _s298, _cp_rows, gibbs_rows = _runtime_function_thermo_data(
            function,
            database or self.current_database,
        )
        for index, row in enumerate(gibbs_rows):
            if np.isclose(float(row[1]), payload_upper, rtol=1e-10, atol=1e-8):
                return index
        return fallback

    def _compound_thermo_range_changed(
        self,
        payload: ThermoRangePayload,
        gibbs_rows: list[list[float]],
    ) -> None:
        """Persist edits from a compound Cp editor back to its parameter."""
        parameter = _compound_parameter_for_range_payload(
            self.current_database,
            payload,
        )
        if parameter is None:
            return
        _set_parameter_gibbs_rows(parameter, gibbs_rows)

    def _solution_endmember_thermo_range_changed(
        self,
        payload: ThermoRangePayload,
        gibbs_rows: list[list[float]],
    ) -> None:
        """Persist Cp editor edits back to a solution endmember parameter."""
        parameter = _solution_endmember_parameter_for_range_payload(
            self.current_database,
            payload,
        )
        if parameter is None:
            return
        _set_parameter_gibbs_rows(parameter, gibbs_rows)
        parameter.source.command = _solution_parameter_command_with_expression(
            parameter,
            parameter.expression,
        )

    def _function_for_range_payload(
        self,
        payload: ThermoRangePayload,
    ) -> FunctionDefinition:
        """Return the Function object represented by a range payload."""
        for function in self.current_database.functions:
            if function.name == payload.owner_name:
                return function
        ranges = []
        if payload.lower != "" and payload.upper != "":
            ranges = [(float(payload.lower), float(payload.upper))]
        source = (
            payload.source
            if isinstance(payload.source, SourceRef)
            else SourceRef()
        )
        return FunctionDefinition(
            name=payload.owner_name,
            expression=payload.expression,
            temperature_ranges=ranges,
            source=source,
        )

    def _function_cp_range_changed(
        self,
        function: FunctionDefinition,
        range_index: int,
    ) -> None:
        ranges = function.temperature_ranges
        if range_index < len(ranges):
            self.editor_title.setText(f"Cp {_compact_number(ranges[range_index][1])}")
        self._select_function_cp_range(
            self.current_database,
            function,
            min(range_index, max(0, len(ranges) - 1)),
        )

    def _render_parameter(self, parameter: Parameter, display_name: str) -> None:
        self._set_editor_heading(
            display_name or f"{parameter.parameter_type}: {parameter.phase_name}",
            ", ".join(parameter.target),
        )
        self._clear_form()
        section = _form_section("Parameter Properties")
        form = _form_layout()
        form.addRow("Phase:", _read_only_line(parameter.phase_name))
        form.addRow("Type:", _combo_box([parameter.parameter_type]))
        form.addRow("Target:", _read_only_line(", ".join(parameter.target)))
        form.addRow("Order:", _read_only_line(parameter.order))
        form.addRow("Expression:", _read_only_text(parameter.expression, 100))
        section.layout().addLayout(form)
        self.editor_form_layout.insertWidget(0, section)
        self._set_source_preview(parameter.source)

    def _render_diagnostic(self, diagnostic: Diagnostic, display_name: str) -> None:
        self._set_editor_heading(
            display_name or diagnostic.severity,
            diagnostic.message,
        )
        self._clear_form()
        section = _form_section("Diagnostic")
        form = _form_layout()
        form.addRow("Severity:", _read_only_line(diagnostic.severity))
        form.addRow("Message:", _read_only_text(diagnostic.message, 90))
        section.layout().addLayout(form)
        self.editor_form_layout.insertWidget(0, section)
        self._set_source_preview(diagnostic.source)

    def _clear_form(self) -> None:
        self._active_database_thermo_editor = None
        while self.editor_form_layout.count():
            item = self.editor_form_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.setParent(None)
                widget.deleteLater()

    def _set_editor_heading(self, title: str, subtitle: str = "") -> None:
        self._editable_editor_title_record = None
        self.editor_badge.setText(self.current_database.name)
        self.editor_title.setText(title)
        self.editor_title.setReadOnly(True)
        self.editor_title.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.editor_title.setToolTip("")
        self.editor_subtitle.setText(subtitle)
        self.editor_subtitle.show()
        self._clear_editor_subtitle_widget()

    def _set_editor_subtitle_widget(self, widget: QWidget) -> None:
        self.editor_subtitle.hide()
        if self.editor_subtitle_control_layout is None:
            return
        self._clear_editor_subtitle_widget()
        self.editor_subtitle_control_layout.addWidget(widget)
        self.editor_subtitle_control.show()

    def _clear_editor_subtitle_widget(self) -> None:
        if self.editor_subtitle_control_layout is None:
            return
        while self.editor_subtitle_control_layout.count():
            item = self.editor_subtitle_control_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.setParent(None)
                widget.deleteLater()
        self.editor_subtitle_control.hide()

    def _toggle_left_sidebar(self) -> None:
        current_sidebar = (
            self.calculation_sidebar
            if self.workspace_stack.currentIndex() == 0
            else self.database_sidebar
        )
        if current_sidebar is not None:
            current_sidebar.setVisible(not current_sidebar.isVisible())

    def _toggle_bottom_panel(self) -> None:
        if self.workspace_stack.currentIndex() == 0 and self.calculation_log_tabs:
            self.calculation_log_tabs.setVisible(
                not self.calculation_log_tabs.isVisible()
            )
        elif self.workspace_stack.currentIndex() == 1 and self.database_inspector_tabs:
            self.database_inspector_tabs.setVisible(
                not self.database_inspector_tabs.isVisible()
            )

    def _toggle_right_chat(self) -> None:
        current_chat = (
            self.calculation_chat_sidebar
            if self.workspace_stack.currentIndex() == 0
            else self.database_chat_sidebar
        )
        if current_chat is not None:
            current_chat.setVisible(not current_chat.isVisible())

    def _set_source_preview(self, source: SourceRef | None) -> None:
        if source and source.file:
            self.source_preview.setPlainText(
                f"{source.file}:{source.line}:{source.column}\n\n{source.command}"
            )
        else:
            self.source_preview.setPlainText(
                "No source mapping for this object yet.\n\n"
                "The parser will populate this panel with command text and line "
                "numbers once TDB reading is connected."
            )

    def _show_validation_summary(self) -> None:
        database = self.current_database
        total = _database_validation_object_count(database)
        progress = _runtime_progress_dialog(
            _validation_progress_message(0, total, "Preparing validation"),
            None,
            0,
            max(total, 1),
            self,
        )
        progress.setWindowTitle("Validate Database")
        progress.setWindowModality(Qt.WindowModality.ApplicationModal)
        progress.setMinimumDuration(0)
        progress.setAutoClose(False)
        progress.setAutoReset(False)
        progress.setValue(0)
        progress.show()
        self._set_status_bar_message(
            _validation_progress_message(0, total, "Validating database"),
        )
        QApplication.processEvents()

        def update_progress(current: int, total_items: int, label: str) -> None:
            message = _validation_progress_message(current, total_items, label)
            progress.setMaximum(max(total_items, 1))
            progress.setValue(min(current, max(total_items, 1)))
            progress.setLabelText(message)
            self._set_status_bar_message(message)
            QApplication.processEvents()

        self._run_global_database_validation(database, update_progress)
        progress.setValue(max(total, 1))
        progress.close()
        diagnostics = self._database_diagnostics(database)
        self._refresh_database_diagnostics()
        error_count = sum(
            1
            for diagnostic in diagnostics
            if diagnostic.severity.lower() == "error"
        )
        warning_count = sum(
            1
            for diagnostic in diagnostics
            if diagnostic.severity.lower() == "warning"
        )
        if diagnostics:
            message = (
                f"{database.name} validation complete: "
                f"{error_count} error(s), {warning_count} warning(s), "
                f"{len(diagnostics)} diagnostic(s) total."
            )
        else:
            message = f"{database.name} validation complete: no diagnostics found."
        QMessageBox.information(self, "DatabaseIR Validation", message)
        level = "error" if error_count else "warning" if warning_count else "ready"
        self._set_status_bar_message(
            message,
            level=level,
        )

    def _database_diagnostics(
        self,
        database: DatabaseIR | None = None,
    ) -> list[Diagnostic]:
        """Return cached full validation diagnostics or lightweight messages."""
        target = database or self.current_database
        diagnostics = self.database_validation_diagnostics.get(id(target))
        if diagnostics is not None:
            return list(diagnostics)
        return target.validation_messages()

    def _run_global_database_validation(
        self,
        database: DatabaseIR | None = None,
        progress_callback: Callable[[int, int, str], None] | None = None,
    ) -> None:
        """Run the explicit Validate-button checks for one selected database."""
        target = database or self.current_database
        self._clear_validation_notices_for_database(target)
        diagnostics = self._validate_database_thermodynamics(
            target,
            progress_callback=progress_callback,
        )
        self.database_validation_diagnostics[id(target)] = diagnostics
        self._refresh_tree_notice_state()

    def _clear_validation_notices_for_database(self, database: DatabaseIR) -> None:
        """Clear cached validation notices only for one database."""
        for function in database.functions:
            self.function_name_warnings.pop(id(function), None)
            self.function_name_errors.pop(id(function), None)
            self.function_continuity_warnings.pop(id(function), None)
        function_prefixes = [
            f"Function|{function.name}|"
            for function in database.functions
        ]
        compound_prefixes: list[str] = []
        phase_counts: dict[str, int] = {}
        for phase in database.phases:
            if _is_solution_phase(phase):
                continue
            compound_name = str(phase.name).strip() or "Compound"
            phase_counts[compound_name] = phase_counts.get(compound_name, 0) + 1
            phase_label = _compound_phase_display_label(
                phase,
                phase_counts[compound_name],
            )
            compound_prefixes.append(f"Compound|{compound_name}:{phase_label}|")
        self.function_range_warnings = {
            key: warnings
            for key, warnings in self.function_range_warnings.items()
            if not any(key.startswith(prefix) for prefix in function_prefixes)
        }
        self.compound_range_warnings = {
            key: warnings
            for key, warnings in self.compound_range_warnings.items()
            if not any(key.startswith(prefix) for prefix in compound_prefixes)
        }
        self.database_validation_diagnostics.pop(id(database), None)

    def _validate_database_thermodynamics(
        self,
        database: DatabaseIR,
        progress_callback: Callable[[int, int, str], None] | None = None,
    ) -> list[Diagnostic]:
        """Return full parser, Function, and compound thermo diagnostics."""
        total = _database_validation_object_count(database)
        checked = 0

        def advance(label: str) -> None:
            nonlocal checked
            checked += 1
            if progress_callback is not None:
                progress_callback(checked, total, label)

        if progress_callback is not None:
            progress_callback(0, total, "Starting validation")
        diagnostics = list(database.validation_messages())
        diagnostics.extend(
            self._validate_function_thermodynamics(database, progress_callback=advance)
        )
        diagnostics.extend(
            self._validate_compound_thermodynamics(database, progress_callback=advance)
        )
        diagnostics.extend(
            self._validate_solution_thermodynamics(database, progress_callback=advance)
        )
        return _unique_diagnostics(diagnostics)

    def _validate_function_thermodynamics(
        self,
        database: DatabaseIR,
        progress_callback: Callable[[str], None] | None = None,
    ) -> list[Diagnostic]:
        """Run strict Function checks, including recursive thermo conversion."""
        diagnostics: list[Diagnostic] = []
        for function in database.functions:
            name_warnings = _function_name_warning_messages(function)
            if name_warnings:
                self.function_name_warnings[id(function)] = list(name_warnings)
                diagnostics.extend(
                    Diagnostic(
                        "warning",
                        f"Function {function.name}: {message}",
                        function.source,
                    )
                    for message in name_warnings
                )
            name_errors = _function_name_error_messages(function)
            if name_errors:
                self.function_name_errors[id(function)] = list(name_errors)
                diagnostics.extend(
                    Diagnostic(
                        "error",
                        f"Function {function.name}: {message}",
                        function.source,
                    )
                    for message in name_errors
                )

            _h298, _s298, _cp_rows, gibbs_rows = _runtime_function_thermo_data(
                function,
                database,
            )
            if not gibbs_rows:
                if function.expression.strip() or function.gibbs_ranges:
                    diagnostics.append(
                        Diagnostic(
                            "warning",
                            (
                                f"Function {function.name}: could not be "
                                "converted into H/S/Cp/G rows."
                            ),
                            function.source,
                        )
                    )
                self._set_function_range_continuity_warnings(
                    function,
                    gibbs_rows=[],
                )
                if progress_callback is not None:
                    progress_callback(f"Function {function.name}")
                continue

            warnings = _function_gibbs_continuity_warnings(gibbs_rows)
            if warnings:
                self.function_continuity_warnings[id(function)] = list(warnings)
                diagnostics.extend(
                    Diagnostic(
                        "warning",
                        f"Function {function.name}: {message}",
                        function.source,
                    )
                    for message in warnings
                )
            self._set_function_range_continuity_warnings(
                function,
                gibbs_rows=gibbs_rows,
            )
            if progress_callback is not None:
                progress_callback(f"Function {function.name}")
        return diagnostics

    def _validate_compound_thermodynamics(
        self,
        database: DatabaseIR,
        progress_callback: Callable[[str], None] | None = None,
    ) -> list[Diagnostic]:
        """Run full compound Gibbs checks with Function dependencies resolved."""
        diagnostics: list[Diagnostic] = []
        phase_counts: dict[str, int] = {}
        for phase in database.phases:
            if _is_solution_phase(phase):
                continue
            compound_name = str(phase.name).strip() or "Compound"
            phase_counts[compound_name] = phase_counts.get(compound_name, 0) + 1
            phase_label = _compound_phase_display_label(
                phase,
                phase_counts[compound_name],
            )
            owner_name = f"{compound_name}:{phase_label}"
            for parameter in _compound_phase_g_parameters(database, phase):
                function = FunctionDefinition(
                    name=owner_name.replace(":", "_"),
                    expression=_parameter_as_function_expression(parameter),
                    source=parameter.source,
                )
                _h298, _s298, _cp_rows, gibbs_rows = _runtime_function_thermo_data(
                    function,
                    database,
                )
                payload = ThermoRangePayload(
                    owner_name=owner_name,
                    owner_kind="Compound",
                    property_name="Cp",
                    index=1,
                    expression=function.expression,
                    source=parameter.source,
                )
                if not gibbs_rows:
                    if function.expression.strip():
                        diagnostics.append(
                            Diagnostic(
                                "warning",
                                (
                                    f"Compound {compound_name} phase "
                                    f"{phase_label}: could not be converted "
                                    "into H/S/Cp/G rows."
                                ),
                                parameter.source,
                            )
                        )
                    continue
                warnings = _function_gibbs_continuity_warnings(gibbs_rows)
                if warnings:
                    diagnostics.extend(
                        Diagnostic(
                            "warning",
                            (
                                f"Compound {compound_name} phase "
                                f"{phase_label}: {message}"
                            ),
                            parameter.source,
                        )
                        for message in warnings
                    )
                self._set_compound_range_warnings_from_gibbs_rows(
                    payload,
                    function,
                    gibbs_rows,
                )
            if progress_callback is not None:
                progress_callback(f"Compound phase {compound_name}:{phase_label}")
        return diagnostics

    def _validate_solution_thermodynamics(
        self,
        database: DatabaseIR,
        progress_callback: Callable[[str], None] | None = None,
    ) -> list[Diagnostic]:
        """Count solution phases in full validation without eager solution checks."""
        for phase in database.phases:
            if not _is_solution_phase(phase):
                continue
            if progress_callback is not None:
                progress_callback(f"Solution phase {phase.name}")
        return []

    def _refresh_tree_notice_state(self) -> None:
        """Push current error/warning maps into the visible database tree."""
        self.tree_model.function_warnings = self._tree_function_warnings()
        self.tree_model.function_errors = self.function_name_errors
        self.tree_model.range_warnings = self._tree_range_warnings()
        for row in range(self.tree_model.rowCount(QModelIndex())):
            index = self.tree_model.index(row, 0, QModelIndex())
            self._emit_tree_subtree_changed(index)
        self.tree.viewport().update()

    def _open_database_ir_file(self) -> None:
        paths, _selected_filter = QFileDialog.getOpenFileNames(
            self,
            "Open DatabaseIR or TDB",
            "",
            "Equilipy DB or TDB (*.eqdb *.json *.tdb);;Equilipy DB (*.eqdb);;"
            "TDB files (*.tdb);;JSON files (*.json);;All files (*)",
        )
        if not paths:
            return

        databases: list[DatabaseIR] = []
        failed: list[str] = []
        skipped: list[str] = []
        loaded_paths = _loaded_database_file_identities(self.databases)
        selected_paths: set[str] = set()
        for path in paths:
            identity = _database_file_identity_from_path(path)
            if identity in loaded_paths or identity in selected_paths:
                skipped.append(Path(path).name)
                continue
            selected_paths.add(identity)
            try:
                database = load_database_ir(path)
                _remember_database_loaded_path(database, path)
                databases.append(database)
            except Exception as exc:  # pragma: no cover - GUI error path.
                failed.append(f"{Path(path).name}: {exc}")
        if failed:
            QMessageBox.critical(
                self,
                "Open Database Failed",
                "\n".join(failed),
            )
            self._set_status_bar_message(
                f"Failed to open {len(failed)} database file(s).",
                level="error",
            )
        if not databases:
            if skipped:
                self._set_status_bar_message(
                    f"Skipped {len(skipped)} already loaded database file(s).",
                    level="ready",
                )
            return
        self._add_databases(databases)
        if skipped:
            self._set_status_bar_message(
                (
                    f"Loaded {len(databases)} database file(s); "
                    f"skipped {len(skipped)} duplicate file(s)."
                ),
                level="ready",
            )

    def _add_databases(self, databases: list[DatabaseIR]) -> None:
        existing_paths = _loaded_database_file_identities(self.databases)
        unique_databases: list[DatabaseIR] = []
        for database in databases:
            identity = _database_file_identity(database)
            if identity and identity in existing_paths:
                continue
            if identity:
                existing_paths.add(identity)
            unique_databases.append(database)
        if not unique_databases:
            self._set_status_bar_message(
                "Skipped already loaded database file(s).",
                level="ready",
            )
            return
        if (
            len(self.databases) == 1
            and not _database_ir_has_content(self.databases[0])
            and not self.databases[0].name.strip()
        ):
            self.databases = list(unique_databases)
        else:
            self.databases.extend(unique_databases)
        self.database = self.databases[0]
        self.current_database = unique_databases[-1]
        self._rebuild_database_tree(self.current_database)
        if self.database_diagnostics_table is not None:
            self.database_diagnostics_table.setModel(
                DiagnosticsModel(self._database_diagnostics())
            )
            self.database_diagnostics_table.resizeColumnsToContents()
        self._set_status_bar_message(
            f"Loaded {len(unique_databases)} database file(s).",
            level="ready",
        )
        self._show_pressure_function_drop_warnings(unique_databases)

    def _show_pressure_function_drop_warnings(
        self,
        databases: list[DatabaseIR],
    ) -> None:
        """Warn when TDB pressure functions were intentionally skipped."""
        messages: list[str] = []
        for database in databases:
            dropped = [
                name.strip()
                for name in database.metadata.get(
                    "dropped_pressure_functions",
                    "",
                ).split(",")
                if name.strip()
            ]
            if not dropped:
                continue
            shown = ", ".join(dropped[:8])
            if len(dropped) > 8:
                shown = f"{shown}, ..."
            messages.append(
                f"{database.name}: {shown}"
            )
        if not messages:
            return
        self._set_status_bar_message(
            "Pressure-dependent TDB functions were dropped during import.",
            level="warning",
        )
        QMessageBox.warning(
            self,
            "Pressure Functions Dropped",
            (
                "The following pressure-dependent TDB function(s) were skipped "
                "during import because Equilipy already applies the gas pressure "
                "correction internally:\n\n"
                + "\n".join(messages)
            ),
        )

    def _replace_database(self, database: DatabaseIR) -> None:
        """Replace loaded databases with one database."""
        self.databases = [database]
        self.database = database
        self.current_database = database
        self._rebuild_database_tree(database)
        self._refresh_database_diagnostics()

    def _refresh_database_diagnostics(self) -> None:
        if self.database_diagnostics_table is not None:
            self.database_diagnostics_table.setModel(
                DiagnosticsModel(self._database_diagnostics())
            )
            self.database_diagnostics_table.resizeColumnsToContents()

    def _placeholder(self, label: str) -> None:
        QMessageBox.information(
            self,
            "Not Connected Yet",
            f"{label} will be connected after parser/writer integration.",
        )

    def _sidebar_action(self, label: str) -> None:
        if label == "Open":
            self._open_database_ir_file()
            return
        if label == "Validate":
            self._show_validation_summary()
            return
        if label == "Export":
            self._export_database_file()
            return
        self._placeholder(label)

    def _export_database_file(self) -> None:
        database = self.current_database
        if not _database_ir_has_content(database):
            QMessageBox.warning(
                self,
                "Export Database",
                "Open a database before exporting.",
            )
            return
        export_database = self._database_for_export()
        if export_database is None:
            return
        default_name = f"{_safe_file_stem(export_database.name)}.eqdb"
        path, selected_filter = QFileDialog.getSaveFileName(
            self,
            "Export Database",
            default_name,
            (
                "Equilipy DB (*.eqdb);;TDB files (*.tdb);;"
                "JSON files (*.json);;All files (*)"
            ),
        )
        if not path:
            return
        output_path = Path(path)
        filter_text = selected_filter.lower()
        suffix = output_path.suffix.lower()
        if suffix == ".tdb" or ("tdb" in filter_text and suffix == ""):
            if output_path.suffix.lower() != ".tdb":
                output_path = output_path.with_suffix(".tdb")
            export_options = self._tdb_export_options_dialog()
            if export_options is None:
                return
            self._write_database_tdb(
                output_path, database=export_database, **export_options
            )
            return
        if suffix not in {".eqdb", ".json"}:
            output_path = output_path.with_suffix(".eqdb")
        self._write_database_eqdb(
            output_path, title="Export Equilipy DB", database=export_database
        )

    def _export_database_tdb(self) -> None:
        database = self.current_database
        if not _database_ir_has_content(database):
            QMessageBox.warning(
                self,
                "Export TDB",
                "Open a database before exporting.",
            )
            return
        export_database = self._database_for_export()
        if export_database is None:
            return
        default_name = f"{_safe_file_stem(export_database.name)}.tdb"
        path, _selected_filter = QFileDialog.getSaveFileName(
            self,
            "Export TDB",
            default_name,
            "TDB files (*.tdb);;All files (*)",
        )
        if not path:
            return
        output_path = Path(path)
        if output_path.suffix.lower() != ".tdb":
            output_path = output_path.with_suffix(".tdb")
        export_options = self._tdb_export_options_dialog()
        if export_options is None:
            return
        self._write_database_tdb(
            output_path, database=export_database, **export_options
        )

    def _database_for_export(self) -> DatabaseIR | None:
        """Apply the element-subset dialog to the current database.

        Returns the database unchanged when every element stays checked, a
        split database for a proper subset, or None when the user cancels.
        """
        database = self.current_database
        selected = _ask_export_element_subset(self, database)
        if selected is None:
            return None
        real_element_count = sum(
            1
            for element in database.elements
            if element.symbol.upper() not in _EXPORT_PSEUDO_ELEMENTS
        )
        if not selected or len(selected) == real_element_count:
            return database
        try:
            return split_database(database, selected)
        except Exception as exc:  # pragma: no cover - GUI error path.
            QMessageBox.critical(self, "Export Element Split Failed", str(exc))
            return None

    def _tdb_export_options_dialog(self) -> dict[str, str] | None:
        """Ask which TDB export style to use."""
        dialog = QDialog(self)
        dialog.setWindowTitle("Export TDB")
        layout = QVBoxLayout(dialog)
        layout.addWidget(
            QLabel(
                "Choose the TDB export style.",
                dialog,
            )
        )
        normal_button = QRadioButton(
            "Normal TDB",
            dialog,
        )
        normal_button.setChecked(True)
        factsage_button = QRadioButton(
            "FactSage-style TDB",
            dialog,
        )
        layout.addWidget(normal_button)
        layout.addWidget(factsage_button)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok
            | QDialogButtonBox.StandardButton.Cancel,
            dialog,
        )
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        layout.addWidget(buttons)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return None
        return {
            "export_style": "factsage" if factsage_button.isChecked() else "native"
        }

    def _write_database_tdb(
        self,
        output_path: Path,
        *,
        include_endmember_compounds: bool = False,
        export_style: str = "native",
        database: DatabaseIR | None = None,
    ) -> None:
        database = self.current_database if database is None else database
        try:
            write_tdb(
                database,
                output_path,
                include_endmember_compounds=include_endmember_compounds,
                export_style=export_style,
            )
        except Exception as exc:  # pragma: no cover - GUI error path.
            QMessageBox.critical(self, "Export TDB Failed", str(exc))
            return
        message = f"Exported TDB to:\n{output_path}"
        subset = database.metadata.get("split_elements", "")
        if subset:
            message += f"\nElement subset: {subset}"
        QMessageBox.information(
            self,
            "Export TDB",
            message,
        )

    def _save_database_eqdb(self) -> None:
        database = self.current_database
        if not _database_ir_has_content(database):
            QMessageBox.warning(
                self,
                "Save Equilipy DB",
                "Open a database before saving.",
            )
            return
        default_name = f"{_safe_file_stem(database.name)}.eqdb"
        path, _selected_filter = QFileDialog.getSaveFileName(
            self,
            "Save Equilipy DB",
            default_name,
            "Equilipy DB (*.eqdb);;JSON files (*.json);;All files (*)",
        )
        if not path:
            return
        output_path = Path(path)
        if output_path.suffix.lower() not in {".eqdb", ".json"}:
            output_path = output_path.with_suffix(".eqdb")
        self._write_database_eqdb(output_path, title="Save Equilipy DB")

    def _write_database_eqdb(
        self,
        output_path: Path,
        *,
        title: str,
        database: DatabaseIR | None = None,
    ) -> bool:
        database = self.current_database if database is None else database
        try:
            write_eqdb(database, output_path)
        except Exception as exc:  # pragma: no cover - GUI error path.
            QMessageBox.critical(self, f"{title} Failed", str(exc))
            return False
        QMessageBox.information(
            self,
            title,
            f"Saved DatabaseIR to:\n{output_path}",
        )
        return True
