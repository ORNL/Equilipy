"""Calculation payload, script, and session serialization helpers."""
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
from pprint import pformat
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


from .input_condition import *
from .results import *

def _safe_file_stem(value: str) -> str:
    stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip()).strip("._")
    return stem or "equilipy"


def _eq_source_path(path: str | Path) -> str:
    return str(Path(path).expanduser().resolve())


def _database_file_identity_from_path(path: str | Path) -> str:
    """Return a stable absolute identity for a database file path."""
    try:
        return str(Path(path).expanduser().resolve())
    except OSError:
        return str(Path(path).expanduser().absolute())


def _remember_database_loaded_path(database: DatabaseIR, path: str | Path) -> None:
    """Store the runtime file identity used to prevent duplicate GUI loads."""
    database.metadata["_loaded_path"] = _database_file_identity_from_path(path)


def _database_file_identity(database: DatabaseIR) -> str:
    """Return the runtime file identity for a loaded database, if known."""
    for key in ("_loaded_path", "source_file"):
        value = str(database.metadata.get(key, "")).strip()
        if value:
            return _database_file_identity_from_path(value)
    return ""


def _loaded_database_file_identities(databases: list[DatabaseIR]) -> set[str]:
    """Return known loaded file identities for duplicate-load prevention."""
    return {
        identity
        for identity in (_database_file_identity(database) for database in databases)
        if identity
    }


def _unique_numbered_name(
    value: str,
    existing_names: list[str],
    *,
    fallback: str,
) -> str:
    candidate = str(value or "").strip() or fallback
    existing = set(existing_names)
    if candidate not in existing:
        return candidate

    match = re.match(r"^(.*?#)(\d+)$", candidate)
    if match is not None:
        prefix = match.group(1)
        index = int(match.group(2)) + 1
    else:
        prefix = f"{candidate}#"
        index = 2

    while f"{prefix}{index}" in existing:
        index += 1
    return f"{prefix}{index}"


def _eq_path_for_save(path: str | Path) -> Path:
    """Return the normalized `.eq` path used for a save target."""
    eq_path = Path(path).expanduser()
    if eq_path.suffix.lower() != ".eq":
        eq_path = eq_path.with_suffix(".eq")
    return eq_path


def _resolve_saved_database_path(
    path: str | Path,
    base_dir: Path | None = None,
) -> Path:
    """Resolve an absolute or `.eq`-relative saved database path."""
    database_path = Path(path).expanduser()
    if not database_path.is_absolute() and base_dir is not None:
        database_path = base_dir / database_path
    return database_path.resolve()


def _unique_database_copy_path(
    database_dir: Path,
    source_path: Path,
    used_names: set[str],
    *,
    avoid_existing: bool = True,
) -> Path:
    """Return a collision-free database copy path under `database_dir`."""
    stem = _safe_file_stem(source_path.stem or "database")
    suffix = source_path.suffix or ".dat"
    candidate = f"{stem}{suffix}"
    index = 2
    while candidate in used_names or (
        avoid_existing and (database_dir / candidate).exists()
    ):
        existing_path = database_dir / candidate
        try:
            if existing_path.resolve() == source_path.resolve():
                break
        except FileNotFoundError:
            pass
        candidate = f"{stem}_{index}{suffix}"
        index += 1
    used_names.add(candidate)
    return database_dir / candidate


def _copy_database_for_directory(
    database_path: str | Path,
    target_directory: str | Path,
) -> Path:
    """Copy one database into `target_directory/databases` and return the copy."""
    source_path = Path(database_path).expanduser().resolve()
    if not source_path.exists():
        raise FileNotFoundError(f"Database file does not exist: {source_path}")
    database_dir = Path(target_directory).expanduser().resolve() / "databases"
    database_dir.mkdir(parents=True, exist_ok=True)
    target_path = _unique_database_copy_path(database_dir, source_path, set())
    if source_path != target_path.resolve():
        shutil.copy2(source_path, target_path)
    return target_path.resolve()


def _self_contained_eq_payload(
    payload: dict[str, Any],
    save_path: str | Path,
) -> dict[str, Any]:
    """Copy saved databases beside an `.eq` file and rewrite paths as relative."""
    prepared = copy.deepcopy(payload)
    eq_path = _eq_path_for_save(save_path)
    base_dir = eq_path.parent.resolve()
    database_dir = base_dir / "databases"
    database_entries = (prepared.get("session") or {}).get("databases") or []
    if not database_entries:
        return prepared

    database_dir.mkdir(parents=True, exist_ok=True)
    used_names: set[str] = set()
    copied_paths: dict[str, str] = {}

    for entry in database_entries:
        if not isinstance(entry, dict):
            continue
        original_text = str(entry.get("path") or "").strip()
        if not original_text:
            continue
        source_path = _resolve_saved_database_path(original_text, base_dir)
        if not source_path.exists():
            raise FileNotFoundError(f"Database file does not exist: {source_path}")
        target_path = _unique_database_copy_path(
            database_dir,
            source_path,
            used_names,
            avoid_existing=False,
        )
        if source_path != target_path.resolve():
            shutil.copy2(source_path, target_path)
        relative_path = target_path.relative_to(base_dir).as_posix()
        copied_paths[original_text] = relative_path
        copied_paths[str(source_path)] = relative_path
        copied_paths[str(source_path.resolve())] = relative_path
        entry["original_path"] = str(source_path.resolve())
        entry["path"] = relative_path

    selected_database_path = ""
    for entry in database_entries:
        if isinstance(entry, dict) and entry.get("selected"):
            selected_database_path = str(entry.get("path") or "")
            break
    if not selected_database_path and isinstance(database_entries[0], dict):
        selected_database_path = str(database_entries[0].get("path") or "")

    for module in prepared.get("modules") or []:
        if not isinstance(module, dict):
            continue
        script = module.get("script")
        if not isinstance(script, dict):
            continue
        database_path = str(script.get("database_path") or "")
        script["database_path"] = copied_paths.get(
            database_path,
            selected_database_path,
        )

    return prepared


def _write_eq_payload(path: str, payload: dict[str, Any]) -> Path:
    eq_path = _eq_path_for_save(path)
    eq_path.write_text(
        json.dumps(_json_safe(payload), indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return eq_path


def _save_eq_and_script(
    path: str | Path,
    payload: dict[str, Any],
    script_text: str,
    *,
    script_path: str | Path | None = None,
) -> tuple[Path, Path]:
    eq_path = _write_eq_payload(path, payload)
    output_script_path = (
        Path(script_path).expanduser()
        if script_path is not None
        else eq_path.with_suffix(".py")
    )
    output_script_path.write_text(script_text, encoding="utf-8")
    return eq_path, output_script_path


def _write_result_table_csv(table: QTableWidget, path: str | Path) -> Path:
    output_path = Path(path).expanduser()
    headers = [
        table.horizontalHeaderItem(column).text()
        if table.horizontalHeaderItem(column) is not None
        else ""
        for column in range(table.columnCount())
    ]
    rows = [
        [
            _table_cell_text(table, row, column)
            for column in range(table.columnCount())
        ]
        for row in range(table.rowCount())
    ]
    with output_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(headers)
        writer.writerows(rows)
    return output_path


def _read_eq_payload(path: str) -> dict[str, Any]:
    payload = json.loads(Path(path).expanduser().read_text(encoding="utf-8"))
    if payload.get("format") != "equilipy.eq":
        raise ValueError("The selected file is not an Equilipy .eq file.")
    if int(payload.get("version", 0)) > 1:
        raise ValueError("This .eq file was created by a newer Equilipy GUI.")
    return payload


def _eq_payloads_in_directory(directory: Path) -> list[tuple[Path, dict[str, Any]]]:
    payloads: list[tuple[Path, dict[str, Any]]] = []
    for path in sorted(directory.glob("*.eq")):
        try:
            payloads.append((path, _read_eq_payload(str(path))))
        except Exception:
            continue
    return payloads


def _payload_contains_pickled_results(payload: dict[str, Any]) -> bool:
    for module in payload.get("modules") or []:
        result = module.get("result") if isinstance(module, dict) else None
        if isinstance(result, dict) and result.get("pickle"):
            return True
    return False


def _module_state_result_score(state: dict[str, Any]) -> int:
    """Rank saved module state by how much restorable result data it contains."""
    table_payload = state.get("results_table") or {}
    result_payload = state.get("result") or {}
    score = 0
    if isinstance(table_payload, dict):
        score += 1000 * len(table_payload.get("rows") or [])
        score += 10 * len(table_payload.get("headers") or [])
    if isinstance(result_payload, dict):
        bundle = result_payload.get("bundle")
        if isinstance(bundle, dict):
            score += 100000
        to_dict = result_payload.get("to_dict")
        if isinstance(to_dict, dict):
            score += len(to_dict)
        if result_payload.get("pickle"):
            score += 1
    return score


def _session_state_payload(session: CalculationSession) -> dict[str, Any]:
    return {
        "id": session.id,
        "name": session.name,
        "temperature_unit": session.temperature_unit,
        "pressure_unit": session.pressure_unit,
        "amount_unit": session.amount_unit,
        "default_phase_names": list(session.default_phase_names),
        "result_column_defaults": {
            str(kind): list(columns)
            for kind, columns in session.result_column_defaults.items()
        },
        "databases": [
            {
                "name": database.name,
                "path": database.path,
                "selected": database.selected,
            }
            for database in session.databases.values()
        ],
    }


def _module_state_payload(
    session: CalculationSession,
    module: CalculationModule,
) -> dict[str, Any]:
    runtime = module.runtime
    composition_table = runtime.get("composition_table")
    phase_tree = runtime.get("phase_tree")
    undercooling_tree = runtime.get("nucleation_undercooling_tree")
    results_table = runtime.get("results_table")

    phase_names = list(module.phase_names)
    selected_phase_names: list[str] = []
    if isinstance(phase_tree, QTreeWidget):
        if not phase_names:
            phase_names = _all_phase_names(phase_tree)
        selected_phase_names = _selected_phase_names(phase_tree)
    nucleation_undercooling = dict(module.nucleation_undercooling)
    if isinstance(undercooling_tree, QTreeWidget):
        try:
            nucleation_undercooling.update(
                _nucleation_undercooling_tree_values(
                    undercooling_tree,
                    strict=False,
                )
            )
        except ValueError:
            pass

    state = {
        "id": module.id,
        "kind": module.kind,
        "name": module.name,
        "temperature_unit": module.temperature_unit,
        "pressure_unit": module.pressure_unit,
        "amount_unit": module.amount_unit,
        "calculation_type": module.calculation_type,
        "solidification_model": module.solidification_model,
        "nucleation_undercooling": nucleation_undercooling,
        "batch_condition": {
            key: list(values) for key, values in module.batch_condition.items()
        },
        "batch_condition_path": module.batch_condition_path,
        "batch_cpu_count": module.batch_cpu_count,
        "transition_search": module.transition_search,
        "start_from_liquidus": module.start_from_liquidus,
        "result_columns": list(module.result_columns),
        "phase_names": phase_names,
        "selected_phase_names": selected_phase_names,
        "composition": (
            _composition_table_rows(composition_table)
            if _is_composition_table(composition_table)
            else []
        ),
        "temperature": _line_edit_text(runtime.get("temperature"), "1600"),
        "pressure": _line_edit_text(runtime.get("pressure"), "1"),
        "delta_t": _line_edit_text(
            runtime.get("delta_t"),
            _default_solidification_delta_t(module),
        ),
        "liquid_phase": _line_edit_text(runtime.get("liquid_phase"), "LIQUID"),
        "result": _result_payload(runtime.get("result") or module.result),
        "results_table": (
            _table_payload(results_table)
            if isinstance(results_table, QTableWidget)
            else dict(module.results_table_payload)
        ),
    }
    state["script"] = _script_payload(session, module, state)
    return state


def _line_edit_text(widget: Any, default: str = "") -> str:
    if isinstance(widget, QLineEdit):
        return widget.text()
    return default


def _composition_table_rows(
    table: QTableWidget | CompositionEditor,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for row in range(table.rowCount()):
        species = _table_cell_text(table, row, 0)
        amount = _table_cell_text(table, row, 1)
        if species or amount:
            rows.append({"species": species, "amount": amount})
    return rows


def _set_composition_table_rows(
    table: QTableWidget | CompositionEditor,
    rows: list[dict[str, Any]],
) -> None:
    previous_blocked = table.blockSignals(True)
    table.setRowCount(0)
    for row_data in rows:
        row = table.rowCount()
        table.insertRow(row)
        _set_composition_cell(
            table,
            row,
            0,
            str(row_data.get("species") or ""),
            "e.g. Mg2Si",
        )
        _set_composition_cell(
            table,
            row,
            1,
            str(row_data.get("amount") or ""),
            _composition_amount_placeholder(table),
        )
    if table.rowCount() == 0:
        _add_composition_row(table)
    table.blockSignals(previous_blocked)
    _refresh_composition_amount_editors(table)
    _fit_table_to_rows(table, min_rows=3)


def _table_payload(table: QTableWidget) -> dict[str, Any]:
    headers = [
        table.horizontalHeaderItem(column).text()
        if table.horizontalHeaderItem(column) is not None
        else ""
        for column in range(table.columnCount())
    ]
    rows = []
    for row in range(table.rowCount()):
        rows.append(
            [
                _table_cell_text(table, row, column)
                for column in range(table.columnCount())
            ]
        )
    return {"headers": headers, "rows": rows, "visible": bool(rows)}


def _module_has_displayable_result(module: CalculationModule) -> bool:
    if module.result is not None or module.runtime.get("result") is not None:
        return True
    payload = module.results_table_payload
    return bool(payload.get("headers") and payload.get("rows"))


def _restore_table_payload(
    table: QTableWidget,
    payload: dict[str, Any],
) -> None:
    headers = tuple(str(header) for header in payload.get("headers") or [])
    rows = [tuple(row) for row in payload.get("rows") or []]
    if headers:
        _populate_result_rows(table, headers, rows)
        table.setVisible(bool(rows))
    else:
        table.clear()
        table.setColumnCount(0)
        table.setRowCount(0)
        table.setVisible(False)


def _result_payload(result: Any) -> dict[str, Any]:
    if result is None:
        return {}
    payload: dict[str, Any] = {
        "class": f"{type(result).__module__}.{type(result).__qualname__}",
    }
    if hasattr(result, "to_bundle"):
        try:
            payload["bundle"] = _json_safe(result.to_bundle())
            return payload
        except Exception:
            payload.pop("bundle", None)
    if hasattr(result, "to_dict"):
        try:
            payload["to_dict"] = _json_safe(result.to_dict())
        except Exception:
            payload["to_dict"] = {}
    if payload.get("bundle"):
        return payload
    encoded = _encode_result_object(result)
    if encoded:
        payload["pickle"] = encoded
    return payload


def _encode_result_object(result: Any) -> str:
    try:
        return base64.b64encode(pickle.dumps(result)).decode("ascii")
    except Exception:
        return ""


def _decode_result_object(payload: dict[str, Any]) -> Any | None:
    encoded = payload.get("pickle") if isinstance(payload, dict) else ""
    if not encoded:
        return None
    return pickle.loads(base64.b64decode(encoded.encode("ascii")))


def _decode_result_bundle(payload: dict[str, Any]) -> Any | None:
    bundle = payload.get("bundle") if isinstance(payload, dict) else None
    if not isinstance(bundle, dict):
        return None
    kind = bundle.get("kind")
    if kind in {"equilib", "equilibrium"}:
        return EquilibResult.from_bundle(bundle)
    if kind == "scheil":
        return ScheilResult.from_bundle(bundle)
    return None


def _apply_module_state_to_runtime(
    module: CalculationModule,
    runtime: dict[str, Any],
    state: dict[str, Any],
    *,
    allow_pickle_results: bool = False,
) -> None:
    composition_table = runtime.get("composition_table")
    if _is_composition_table(composition_table):
        _set_composition_table_rows(composition_table, state.get("composition") or [])

    for key, value in (
        ("temperature", state.get("temperature") or "1600"),
        ("pressure", state.get("pressure") or "1"),
        ("delta_t", state.get("delta_t") or _default_solidification_delta_t(module)),
        ("liquid_phase", state.get("liquid_phase") or "LIQUID"),
    ):
        widget = runtime.get(key)
        if isinstance(widget, QLineEdit):
            widget.setText(str(value))

    transition_search = runtime.get("transition_search")
    if isinstance(transition_search, QCheckBox):
        transition_search.setChecked(module.transition_search)

    from_liquidus = runtime.get("from_liquidus")
    if isinstance(from_liquidus, QCheckBox):
        from_liquidus.setChecked(module.start_from_liquidus)

    phase_tree = runtime.get("phase_tree")
    selected_phase_names = list(state.get("selected_phase_names") or [])
    if isinstance(phase_tree, QTreeWidget) and module.phase_names:
        _populate_phase_tree(
            phase_tree,
            module.phase_names,
            selected_phase_names or None,
        )

    _restore_module_result(
        module,
        runtime,
        state,
        allow_pickle_results=allow_pickle_results,
    )


def _restore_module_result(
    module: CalculationModule,
    runtime: dict[str, Any],
    state: dict[str, Any],
    *,
    allow_pickle_results: bool = False,
) -> None:
    table = runtime.get("results_table")
    if not isinstance(table, QTableWidget):
        return

    result = None
    result_payload = state.get("result") or {}
    module.result_payload = (
        dict(result_payload) if isinstance(result_payload, dict) else {}
    )
    module.results_table_payload = dict(state.get("results_table") or {})
    try:
        result = _decode_result_bundle(result_payload)
    except Exception:
        result = None
    if allow_pickle_results:
        if result is None:
            try:
                result = _decode_result_object(result_payload)
            except Exception:
                result = None
    if result is not None:
        module.result = result
        runtime["result"] = result
    else:
        module.result = None
        runtime.pop("result", None)

    table_payload = state.get("results_table") or {}
    if module.result_columns and result is not None:
        _populate_module_results(module, table, result)
        module.results_table_payload = _table_payload(table)
        return

    if table_payload.get("headers"):
        _restore_table_payload(table, table_payload)
        module.results_table_payload = _table_payload(table)
        return

    if isinstance(result_payload, dict):
        result_data = result_payload.get("to_dict")
        if isinstance(result_data, dict) and result_data:
            result_table = _augment_result_table_units(
                ResultTable.from_dict(result_data),
                temperature_unit=module.temperature_unit,
                pressure_unit=module.pressure_unit,
            )
            selected_columns = (
                _expand_result_column_defaults(
                    list(module.result_columns),
                    result_table,
                    module,
                )
                if module.result_columns
                else _default_result_column_keys(module, result_table, result)
            )
            if _populate_result_table_columns(
                table,
                result_table,
                selected_columns,
                module,
            ):
                module.results_table_payload = _table_payload(table)
                return

    if module.kind == "solidification" and isinstance(result_payload, dict):
        result_data = result_payload.get("to_dict")
        if isinstance(result_data, dict) and result_data:
            snapshot = _SavedScheilResultSnapshot(
                result_data,
                [
                    module.temperature_unit,
                    module.pressure_unit,
                    module.amount_unit,
                ],
                str(state.get("liquid_phase") or "LIQUID"),
            )
            _populate_module_results(module, table, snapshot)
            module.results_table_payload = _table_payload(table)
            return

    if result is not None:
        try:
            _populate_module_results(module, table, result)
            module.results_table_payload = _table_payload(table)
            return
        except Exception:
            pass

    _restore_table_payload(table, table_payload)
    module.results_table_payload = _table_payload(table)


def _restore_cached_module_result(
    module: CalculationModule,
    table: QTableWidget,
) -> None:
    if module.result is not None:
        module.runtime["result"] = module.result
    if module.result is not None and module.result_columns:
        try:
            _populate_module_results(module, table, module.result)
        except Exception:
            return
        module.results_table_payload = _table_payload(table)
        _refresh_module_result_view(module)
        return
    if module.results_table_payload.get("headers"):
        _restore_table_payload(table, module.results_table_payload)
        _refresh_module_result_view(module)
        return
    if module.result is None:
        _refresh_module_result_view(module)
        return
    try:
        _populate_module_results(module, table, module.result)
    except Exception:
        return
    module.results_table_payload = _table_payload(table)
    _refresh_module_result_view(module)


class _SavedScheilResultSnapshot:
    """Tiny adapter for saved Scheil result dictionaries."""

    def __init__(
        self,
        data: dict[str, Any],
        input_unit: list[str],
        liquid_phase_name: str,
    ):
        self._data = data
        self.input_unit = input_unit
        self.liquid_phase_name = liquid_phase_name

    def to_dict(self) -> dict[str, Any]:
        """Return the saved flattened result data without touching globals."""
        return self._data


def _script_payload(
    session: CalculationSession,
    module: CalculationModule,
    state: dict[str, Any],
) -> dict[str, Any]:
    database = next(
        (db for db in session.databases.values() if db.selected),
        next(iter(session.databases.values()), None),
    )
    selected_phases = (
        state.get("selected_phase_names") or state.get("phase_names") or []
    )
    payload: dict[str, Any] = {
        "database_path": database.path if database else "",
        "units": [
            module.temperature_unit,
            module.pressure_unit,
            module.amount_unit,
        ],
        "phases": selected_phases,
    }

    try:
        if module.kind == "equilibrium":
            payload.update(_equilibrium_script_payload(database, module, state))
        elif module.kind == "solidification":
            payload.update(_solidification_script_payload(database, state))
        else:
            payload["error"] = "Nucleation is to be implemented."
    except Exception as exc:
        payload["error"] = str(exc)
    return payload


def _equilibrium_script_payload(
    database: CalculationDatabase | None,
    module: CalculationModule,
    state: dict[str, Any],
) -> dict[str, Any]:
    database_object = database.database if database else None
    if module.calculation_type == "batch" and module.batch_condition:
        return {
            "mode": "batch",
            "condition": {
                key: list(values) for key, values in module.batch_condition.items()
            },
            "n_cpu": _effective_batch_cpu_count(
                module.batch_condition,
                module.batch_cpu_count,
            ),
            "include_transitions": False,
        }
    if module.calculation_type == "batch":
        condition = _manual_batch_condition(
            database_object,
            _composition_table_from_state(state),
            str(state.get("temperature") or ""),
            str(state.get("pressure") or ""),
            [module.temperature_unit, module.pressure_unit, module.amount_unit],
        )
        return {
            "mode": "batch",
            "manual_batch": True,
            "manual_batch_temperature": str(state.get("temperature") or ""),
            "manual_batch_pressure": str(state.get("pressure") or ""),
            "manual_batch_composition": list(state.get("composition") or []),
            "condition": condition,
            "n_cpu": _effective_batch_cpu_count(
                condition,
                module.batch_cpu_count,
            ),
            "include_transitions": False,
        }

    temperatures = _temperature_condition_values(str(state.get("temperature") or ""))
    pressure = float(state.get("pressure") or 1)
    base_condition = _condition_from_composition_rows(
        state.get("composition") or [],
        temperatures[0],
        pressure,
        module.amount_unit,
        database_object,
    )
    return {
        "mode": "single",
        "base_condition": base_condition,
        "temperature_values": temperatures,
        "include_transitions": bool(
            module.transition_search
            and module.kind == "equilibrium"
            and _temperature_text_is_range(str(state.get("temperature") or ""))
        ),
    }


def _solidification_script_payload(
    database: CalculationDatabase | None,
    state: dict[str, Any],
) -> dict[str, Any]:
    database_object = database.database if database else None
    is_batch = str(state.get("calculation_type") or "single") == "batch"
    if is_batch and bool(state.get("batch_condition")):
        condition: dict[str, Any] = {
            key: list(values)
            for key, values in (state.get("batch_condition") or {}).items()
        }
    elif is_batch:
        condition = _manual_batch_condition(
            database_object,
            _composition_table_from_state(state),
            str(state.get("temperature") or ""),
            str(state.get("pressure") or ""),
            [
                str(state.get("temperature_unit") or "K"),
                str(state.get("pressure_unit") or "atm"),
                str(state.get("amount_unit") or "moles"),
            ],
        )
        manual_batch_payload = {
            "manual_batch": True,
            "manual_batch_temperature": str(state.get("temperature") or ""),
            "manual_batch_pressure": str(state.get("pressure") or ""),
            "manual_batch_composition": list(state.get("composition") or []),
        }
    else:
        condition = _condition_from_composition_rows(
            state.get("composition") or [],
            float(state.get("temperature") or 1600),
            float(state.get("pressure") or 1),
            str(state.get("amount_unit") or "moles"),
            database_object,
        )
        manual_batch_payload = {}
    return {
        "mode": "batch" if is_batch else "single",
        "solidification_model": _solidification_model_value(
            str(state.get("solidification_model") or "scheil")
        ),
        "nucleation_undercooling": {
            str(phase): float(value)
            for phase, value in (
                state.get("nucleation_undercooling") or {}
            ).items()
        },
        "condition": condition,
        "n_cpu": (
            _effective_batch_cpu_count(
                condition,
                _normalized_batch_cpu_count(state.get("batch_cpu_count", 0)),
            )
            if is_batch
            else 1
        ),
        "delta_t": float(state.get("delta_t") or 5),
        "liquid_phase": str(state.get("liquid_phase") or "LIQUID"),
        "start_from_liquidus": bool(state.get("start_from_liquidus", True)),
        **manual_batch_payload,
    }


def _composition_table_from_state(state: dict[str, Any]) -> QTableWidget:
    table = _composition_table()
    table.setProperty("allow_amount_ranges", True)
    _set_composition_table_rows(table, state.get("composition") or [])
    return table


def _condition_from_composition_rows(
    rows: list[dict[str, Any]],
    temperature: float,
    pressure: float,
    amount_unit: str,
    database: Any | None,
) -> dict[str, Any]:
    table = _composition_table()
    _set_composition_table_rows(table, rows)
    return _calculation_condition(table, temperature, pressure, amount_unit, database)


def _script_header_lines(title: str, *, uses_manual_batch: bool = False) -> list[str]:
    lines = [
        "#!/usr/bin/env python3",
        f'"""{title}"""',
        "",
        "import os",
    ]
    if uses_manual_batch:
        lines.append("from itertools import product")
    lines.extend(
        [
            "import equilipy as eq",
            "import polars as pl",
            "",
            "",
            "fpath = os.path.dirname(os.path.abspath(__file__))",
            "",
        ]
    )
    return lines


def _script_assignment_lines(
    name: str,
    value: Any,
    *,
    indent: str = "",
) -> list[str]:
    literal = pformat(_json_safe(value), width=88, sort_dicts=False)
    return [f"{indent}{line}" for line in f"{name} = {literal}".splitlines()]


def _script_database_lines(
    database_path: str,
    *,
    indent: str = "",
) -> list[str]:
    return [
        f"{indent}# Step 1: Parse database",
        *(
            _script_assignment_lines(
                "database_path",
                database_path,
                indent=indent,
            )
        ),
        f"{indent}if not database_path:",
        f"{indent}    raise RuntimeError('No database path was saved with this module.')",
        f"{indent}if not os.path.isabs(database_path):",
        f"{indent}    database_path = os.path.join(fpath, database_path)",
        f"{indent}if database_path.lower().endswith('.tdb'):",
        f"{indent}    DB = eq.read_tdb(database_path, strict=False, auto_correct=False)",
        f"{indent}else:",
        f"{indent}    DB = eq.read_dat(database_path)",
        "",
    ]


def _script_result_lines(
    result_filename: str,
    *,
    indent: str = "",
    return_result: bool = False,
) -> list[str]:
    lines = [
        "",
        f"{indent}# Step 4: Save and print result",
        *(
            _script_assignment_lines(
                "result_filename",
                result_filename,
                indent=indent,
            )
        ),
        f"{indent}result_path = os.path.join(fpath, result_filename)",
        f"{indent}if hasattr(res, 'to_bundle'):",
        f"{indent}    eq.save_result(res, result_path)",
        f"{indent}res_table = res.to_dict() if hasattr(res, 'to_dict') else res",
        f"{indent}print(pl.DataFrame(res_table) if isinstance(res_table, dict) else res_table)",
    ]
    if return_result:
        lines.append(f"{indent}return res")
    return lines


def _script_range_condition_lines(
    script: dict[str, Any],
    units: list[Any],
    phases: list[Any] | None,
    *,
    indent: str = "",
) -> list[str]:
    base_condition = script.get("base_condition") or {}
    temperatures = script.get("temperature_values") or []
    include_transitions = bool(script.get("include_transitions"))
    lines = [
        *(
            _script_assignment_lines(
                "NPTBase",
                base_condition,
                indent=indent,
            )
        ),
        *(
            _script_assignment_lines(
                "TemperatureValues",
                temperatures,
                indent=indent,
            )
        ),
    ]
    if include_transitions:
        lines.extend(
            [
                *(
                    _script_assignment_lines(
                        "IncludeTransitions",
                        include_transitions,
                        indent=indent,
                    )
                ),
                f"{indent}if IncludeTransitions:",
                f"{indent}    TransitionTemperatures = eq.find_transitions(",
                f"{indent}        DB,",
                f"{indent}        dict(NPTBase),",
                f"{indent}        max(TemperatureValues[0], TemperatureValues[-1]),",
                f"{indent}        min(TemperatureValues[0], TemperatureValues[-1]),",
                f"{indent}        unit=UnitIn,",
                f"{indent}        phases=PhaseSelection,",
                f"{indent}    )",
                f"{indent}    lower_temperature = min(TemperatureValues[0], TemperatureValues[-1])",
                f"{indent}    upper_temperature = max(TemperatureValues[0], TemperatureValues[-1])",
                f"{indent}    reverse_temperature = TemperatureValues[-1] < TemperatureValues[0]",
                f"{indent}    merged_temperatures = list(TemperatureValues)",
                f"{indent}    for transition_temperature in TransitionTemperatures:",
                f"{indent}        if lower_temperature <= transition_temperature <= upper_temperature:",
                f"{indent}            merged_temperatures.append(float(transition_temperature))",
                f"{indent}    TemperatureValues = []",
                f"{indent}    for temperature in sorted(merged_temperatures, reverse=reverse_temperature):",
                f"{indent}        if not any(abs(temperature - existing) <= 1e-8 for existing in TemperatureValues):",
                f"{indent}            TemperatureValues.append(temperature)",
            ]
        )
    lines.extend(
        [
            f"{indent}NPT = {{'T': TemperatureValues}}",
            f"{indent}n_points = len(TemperatureValues)",
            f"{indent}for species, amount in NPTBase.items():",
            f"{indent}    if species == 'T':",
            f"{indent}        continue",
            f"{indent}    NPT[species] = [float(amount)] * n_points",
        ]
    )
    return lines


def _script_uses_manual_batch(module_state: dict[str, Any]) -> bool:
    script = module_state.get("script") or {}
    return bool(script.get("manual_batch"))


def _script_manual_batch_condition_lines(
    script: dict[str, Any],
    *,
    indent: str = "",
) -> list[str]:
    rows = [
        [str(row.get("species") or ""), str(row.get("amount") or "")]
        for row in script.get("manual_batch_composition") or []
        if isinstance(row, dict)
    ]
    lines = [
        *(
            _script_assignment_lines(
                "TemperatureText",
                str(script.get("manual_batch_temperature") or ""),
                indent=indent,
            )
        ),
        *(
            _script_assignment_lines(
                "PressureText",
                str(script.get("manual_batch_pressure") or ""),
                indent=indent,
            )
        ),
        *(_script_assignment_lines("CompositionRows", rows, indent=indent)),
        "",
        f"{indent}def _batch_range_values(text):",
        f"{indent}    values = [float(item) for item in str(text).strip().strip('[]').replace(',', ' ').split()]",
        f"{indent}    if len(values) == 1:",
        f"{indent}        return values",
        f"{indent}    if len(values) != 3:",
        f"{indent}        raise ValueError('Batch ranges must be: initial final step.')",
        f"{indent}    start, stop, step = values",
        f"{indent}    if step == 0 or start == stop:",
        f"{indent}        return [start]",
        f"{indent}    if (stop - start) * step < 0:",
        f"{indent}        raise ValueError('Batch step must move from initial toward final.')",
        f"{indent}    result = [start]",
        f"{indent}    current = start + step",
        f"{indent}    if step > 0:",
        f"{indent}        while current < stop:",
        f"{indent}            result.append(current)",
        f"{indent}            current += step",
        f"{indent}    else:",
        f"{indent}        while current > stop:",
        f"{indent}            result.append(current)",
        f"{indent}            current += step",
        f"{indent}    if abs(result[-1] - stop) > 1e-8:",
        f"{indent}        result.append(stop)",
        f"{indent}    return result",
        "",
        f"{indent}def _fixed_total(amount_unit):",
        f"{indent}    unit = str(amount_unit).strip().lower()",
        f"{indent}    if unit in {{'mol%', 'at%', 'wt%', 'wt.%'}}:",
        f"{indent}        return 100.0",
        f"{indent}    if unit in {{'mole fraction', 'atom fraction', 'mass fraction', 'weight fraction'}}:",
        f"{indent}        return 1.0",
        f"{indent}    return None",
        "",
        f"{indent}def _manual_batch_condition(temperature_text, pressure_text, rows, amount_unit):",
        f"{indent}    temperatures = _batch_range_values(temperature_text)",
        f"{indent}    pressure = float(pressure_text)",
        f"{indent}    fixed_total = _fixed_total(amount_unit)",
        f"{indent}    species_order = []",
        f"{indent}    ranged_rows = []",
        f"{indent}    balance_row = None",
        f"{indent}    for species, amount_text in rows:",
        f"{indent}        species = str(species).strip()",
        f"{indent}        amount_text = str(amount_text).strip()",
        f"{indent}        if not species and not amount_text:",
        f"{indent}            continue",
        f"{indent}        if not species:",
        f"{indent}            raise ValueError('Species name cannot be empty.')",
        f"{indent}        if species not in species_order:",
        f"{indent}            species_order.append(species)",
        f"{indent}        if amount_text.startswith('?'):",
        f"{indent}            if balance_row is not None:",
        f"{indent}                raise ValueError('Only one balance species can be used.')",
        f"{indent}            if fixed_total is None:",
        f"{indent}                total_text = amount_text[1:].strip()",
        f"{indent}                if not total_text:",
        f"{indent}                    raise ValueError(\"Enter a total amount after '?'.\")",
        f"{indent}                target_total = float(total_text)",
        f"{indent}            else:",
        f"{indent}                target_total = fixed_total",
        f"{indent}            balance_row = (species, target_total)",
        f"{indent}        else:",
        f"{indent}            ranged_rows.append((species, _batch_range_values(amount_text)))",
        f"{indent}    if fixed_total is not None and balance_row is None:",
        f"{indent}        raise ValueError('Select one balance species for fraction or percent units.')",
        f"{indent}    condition = {{'T': [], 'P': []}}",
        f"{indent}    for species in species_order:",
        f"{indent}        condition[species] = []",
        f"{indent}    for temperature in temperatures:",
        f"{indent}        for amounts in product(*(values for _species, values in ranged_rows)):",
        f"{indent}            row_amounts = {{}}",
        f"{indent}            known_total = 0.0",
        f"{indent}            for (species, _values), amount in zip(ranged_rows, amounts, strict=False):",
        f"{indent}                row_amounts[species] = row_amounts.get(species, 0.0) + amount",
        f"{indent}                known_total += amount",
        f"{indent}            if balance_row is not None:",
        f"{indent}                species, target_total = balance_row",
        f"{indent}                residual = target_total - known_total",
        f"{indent}                if residual <= 0:",
        f"{indent}                    raise ValueError('The total amount is smaller than the sum.')",
        f"{indent}                row_amounts[species] = row_amounts.get(species, 0.0) + residual",
        f"{indent}            condition['T'].append(temperature)",
        f"{indent}            condition['P'].append(pressure)",
        f"{indent}            for species in species_order:",
        f"{indent}                condition[species].append(row_amounts.get(species, 0.0))",
        f"{indent}    return condition",
        "",
        f"{indent}NPT = _manual_batch_condition(",
        f"{indent}    TemperatureText,",
        f"{indent}    PressureText,",
        f"{indent}    CompositionRows,",
        f"{indent}    UnitIn[2],",
        f"{indent})",
    ]
    return lines


def _script_equilibrium_lines(
    module_state: dict[str, Any],
    script: dict[str, Any],
    units: list[Any],
    phases: list[Any] | None,
    *,
    indent: str = "",
) -> list[str]:
    lines = [
        f"{indent}# Step 2: Set input data",
        *(_script_assignment_lines("UnitIn", units, indent=indent)),
        *(_script_assignment_lines("PhaseSelection", phases, indent=indent)),
    ]
    mode = str(script.get("mode") or "single")
    temperatures = script.get("temperature_values") or []
    include_transitions = bool(script.get("include_transitions"))

    if mode == "batch":
        if script.get("manual_batch"):
            lines.extend(_script_manual_batch_condition_lines(script, indent=indent))
        else:
            lines.extend(
                _script_assignment_lines(
                    "NPT",
                    script.get("condition") or {},
                    indent=indent,
                )
            )
        function_name = "equilib_batch"
    elif len(temperatures) == 1 and not include_transitions:
        lines.extend(
            _script_assignment_lines(
                "NPT",
                script.get("base_condition") or {},
                indent=indent,
            )
        )
        function_name = "equilib_single"
    else:
        lines.extend(
            _script_range_condition_lines(
                script,
                units,
                phases,
                indent=indent,
            )
        )
        function_name = "equilib_batch"

    lines.extend(["", f"{indent}# Step 3: Calculate equilibrium"])
    if function_name == "equilib_single":
        lines.extend(
            [
                f"{indent}res = eq.equilib_single(",
                f"{indent}    DB,",
                f"{indent}    NPT,",
                f"{indent}    unit=UnitIn,",
                f"{indent}    phases=PhaseSelection,",
                f"{indent})",
            ]
        )
    else:
        lines.extend(
            [
                f"{indent}res = eq.equilib_batch(",
                f"{indent}    DB,",
                f"{indent}    NPT,",
                f"{indent}    unit=UnitIn,",
                f"{indent}    phases=PhaseSelection,",
                f"{indent}    n_cpu={int(script.get('n_cpu', 1))!r},",
                f"{indent})",
            ]
        )
    return lines


def _script_solidification_lines(
    module_state: dict[str, Any],
    script: dict[str, Any],
    units: list[Any],
    phases: list[Any] | None,
    *,
    indent: str = "",
) -> list[str]:
    solidification_model = _solidification_model_value(
        str(script.get("solidification_model") or "scheil")
    )
    mode = str(script.get("mode") or "single")
    if solidification_model == "equilibrium_cooling":
        function_name = "equilib_cooling"
    elif solidification_model == "nucleoscheil":
        function_name = "nucleoscheil_batch" if mode == "batch" else "nucleoscheil_cooling"
    else:
        function_name = "scheil_batch" if mode == "batch" else "scheil_cooling"

    lines = [
        f"{indent}# Step 2: Set input data",
        *(_script_assignment_lines("UnitIn", units, indent=indent)),
        *(_script_assignment_lines("PhaseSelection", phases, indent=indent)),
    ]
    if mode == "batch" and script.get("manual_batch"):
        lines.extend(_script_manual_batch_condition_lines(script, indent=indent))
    else:
        lines.extend(
            _script_assignment_lines("NPT", script.get("condition") or {}, indent=indent)
        )
    lines.extend(
        [
            *(
                _script_assignment_lines(
                    "TargetPhase",
                    str(script.get("liquid_phase") or "LIQUID"),
                    indent=indent,
                )
            ),
            "",
            f"{indent}# Step 3: Calculate solidification",
        ]
    )

    if solidification_model == "nucleoscheil":
        lines.extend(
            [
                *(
                    _script_assignment_lines(
                        "CriticalUndercooling",
                        script.get("nucleation_undercooling") or None,
                        indent=indent,
                    )
                ),
                f"{indent}res = eq.{function_name}(",
                f"{indent}    TargetPhase,",
                f"{indent}    DB,",
                f"{indent}    NPT,",
                f"{indent}    CriticalUndercooling,",
                f"{indent}    delta_T={float(script.get('delta_t', 5))!r},",
                f"{indent}    unit=UnitIn,",
                f"{indent}    phases=PhaseSelection,",
            ]
        )
        if mode == "batch":
            lines.append(f"{indent}    n_cpu={int(script.get('n_cpu', 1))!r},")
        lines.append(f"{indent})")
        return lines

    lines.extend(
        [
            f"{indent}res = eq.{function_name}(",
            f"{indent}    TargetPhase,",
            f"{indent}    DB,",
            f"{indent}    NPT,",
            f"{indent}    delta_T={float(script.get('delta_t', 5))!r},",
            f"{indent}    unit=UnitIn,",
            f"{indent}    phases=PhaseSelection,",
            f"{indent}    start_from_liquidus={bool(script.get('start_from_liquidus', True))!r},",
        ]
    )
    if mode == "batch" and solidification_model != "equilibrium_cooling":
        lines.append(f"{indent}    n_cpu={int(script.get('n_cpu', 1))!r},")
    lines.append(f"{indent})")
    return lines


def _module_script_step_lines(
    module_state: dict[str, Any],
    session_state: dict[str, Any],
    *,
    indent: str = "",
    return_result: bool = False,
) -> list[str]:
    del session_state
    script = module_state.get("script") or {}
    database_path = script.get("database_path") or ""
    units = script.get("units") or [
        module_state.get("temperature_unit", "K"),
        module_state.get("pressure_unit", "atm"),
        module_state.get("amount_unit", "moles"),
    ]
    phases = script.get("phases") or None
    module_kind = module_state.get("kind") or "equilibrium"
    module_name = str(module_state.get("name") or module_kind)
    result_filename = f"{_safe_file_stem(module_name)}.eqres"
    error = script.get("error")

    lines = [f"{indent}# Module: {module_name}", ""]
    lines.extend(_script_database_lines(str(database_path), indent=indent))

    if error:
        lines.extend(
            [
                f"{indent}# This module could not be exported as an executable calculation.",
                f"{indent}raise RuntimeError({error!r})",
            ]
        )
        return lines

    if module_kind == "solidification":
        lines.extend(
            _script_solidification_lines(
                module_state,
                script,
                units,
                phases,
                indent=indent,
            )
        )
    else:
        lines.extend(
            _script_equilibrium_lines(
                module_state,
                script,
                units,
                phases,
                indent=indent,
            )
        )

    lines.extend(
        _script_result_lines(
            result_filename,
            indent=indent,
            return_result=return_result,
        )
    )
    return lines


def _module_script_text(
    module_state: dict[str, Any],
    session_state: dict[str, Any],
) -> str:
    lines = _script_header_lines(
        "Run an Equilipy GUI-saved calculation module.",
        uses_manual_batch=_script_uses_manual_batch(module_state),
    )
    lines.extend(_module_script_step_lines(module_state, session_state))
    lines.append("")
    return "\n".join(lines)


def _session_script_text(payload: dict[str, Any]) -> str:
    session_state = payload.get("session") or {}
    modules = payload.get("modules") or []
    if not modules:
        return "\n".join(
            [
                "#!/usr/bin/env python3",
                '"""Run an Equilipy GUI-saved calculation session."""',
                "",
                "def run():",
                "    return []",
                "",
                "",
                "if __name__ == '__main__':",
                "    run()",
                "",
            ]
        )

    lines = _script_header_lines(
        "Run an Equilipy GUI-saved calculation session.",
        uses_manual_batch=any(_script_uses_manual_batch(module) for module in modules),
    )
    for index, module_state in enumerate(modules, start=1):
        function_name = f"run_module_{index}"
        lines.append(f"def {function_name}():")
        lines.extend(
            _module_script_step_lines(
                module_state,
                session_state,
                indent="    ",
                return_result=True,
            )
        )
        lines.append("")
    lines.extend(
        [
            "def run():",
            "    results = []",
        ]
    )
    for index, _module_state in enumerate(modules, start=1):
        lines.append(f"    results.append(run_module_{index}())")
    lines.extend(
        [
            "    return results",
            "",
            "",
            "if __name__ == '__main__':",
            "    run()",
            "",
        ]
    )
    return "\n".join(lines)


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return _json_safe(value.tolist())
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    if hasattr(value, "model_dump"):
        return _json_safe(value.model_dump())
    if hasattr(value, "__dataclass_fields__"):
        return _json_safe(asdict(value))
    return str(value)

__all__ = [name for name in globals() if name.startswith("_") and not name.startswith("__")]
