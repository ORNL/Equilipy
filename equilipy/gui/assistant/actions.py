"""Whitelisted GUI actions requested by assistant providers."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from itertools import product
from pathlib import Path
from typing import Any

from PySide6.QtWidgets import (
    QCheckBox,
    QLineEdit,
    QMessageBox,
    QPlainTextEdit,
    QTreeWidget,
    QWidget,
)

from equilipy.composition import expand_condition_species as _expand_condition_species
from equilipy.gui.calculation.input_condition import (
    _all_phase_names,
    _database_phase_names,
    _is_composition_table,
    _populate_phase_tree,
    _selected_phase_names,
    _set_phase_tree_selection,
)
from equilipy.gui.calculation.results import _solidification_model_value
from equilipy.gui.calculation.session_serialization import _set_composition_table_rows
from equilipy.gui.calculation.state import (
    CalculationDatabase,
    CalculationModule,
    CalculationSession,
)
from equilipy.gui.units import AMOUNT_UNIT_OPTIONS

from .action_contract import (
    ACTION_BLOCK_END,
    ACTION_BLOCK_START,
    strip_json_code_fence,
)


@dataclass(frozen=True)
class AssistantGuiAction:
    """One parsed assistant-requested GUI action."""

    name: str
    args: dict[str, Any] = field(default_factory=dict)


class AssistantActionError(RuntimeError):
    """Raised when a requested assistant GUI action is unsupported or unsafe."""


def parse_action_requests(text: str) -> list[AssistantGuiAction]:
    """Parse action requests embedded in assistant text."""
    actions: list[AssistantGuiAction] = []
    for payload_text in _iter_action_payloads(text):
        payload_text = strip_json_code_fence(payload_text)
        if not payload_text.strip():
            continue
        try:
            payload = json.loads(payload_text)
        except json.JSONDecodeError as exc:
            raise AssistantActionError(f"Invalid assistant action JSON: {exc}") from exc
        raw_actions = payload.get("actions") if isinstance(payload, dict) else payload
        if isinstance(raw_actions, dict):
            raw_actions = [raw_actions]
        if not isinstance(raw_actions, list):
            raise AssistantActionError("Assistant action payload must contain actions.")
        for raw in raw_actions:
            if not isinstance(raw, dict):
                raise AssistantActionError("Assistant action entries must be objects.")
            name = str(raw.get("name") or "").strip()
            if not name:
                raise AssistantActionError("Assistant action entry is missing a name.")
            raw_args = raw.get("args") or {}
            if not isinstance(raw_args, dict):
                raise AssistantActionError(f"Arguments for {name} must be an object.")
            actions.append(AssistantGuiAction(name=name, args=dict(raw_args)))
    return actions


def strip_action_markup(text: str) -> str:
    """Remove assistant action blocks from user-visible transcript text."""
    cleaned = str(text)
    while True:
        start = cleaned.find(ACTION_BLOCK_START)
        if start < 0:
            return cleaned
        end = cleaned.find(ACTION_BLOCK_END, start)
        if end < 0:
            return cleaned[:start].rstrip()
        cleaned = cleaned[:start] + cleaned[end + len(ACTION_BLOCK_END) :]


def execute_assistant_actions(
    owner: Any,
    actions: list[AssistantGuiAction],
) -> list[str]:
    """Execute safe GUI actions and return user-visible status lines."""
    context: dict[str, Any] = {}
    summaries: list[str] = []
    for action in actions:
        if action.name == "calculation.create_session":
            session = _create_calculation_session(owner, action.args)
            context["session"] = session
            summaries.append(f"Created calculation session {session.name}.")
            continue
        if action.name in {
            "calculation.attach_current_database",
            "calculation.load_database",
        }:
            session = _target_session(owner, context, action.args)
            loaded = _load_database_from_project_directory(owner, session, action.args)
            summaries.append(f"Loaded database {loaded.name} into {session.name}.")
            continue
        if action.name == "calculation.add_module":
            session = _target_session(owner, context, action.args)
            module_kind = _normalize_module_kind(action.args.get("kind"))
            module = owner._add_calculation_module(session, module_kind)
            owner._show_calculation_module(session, module)
            context["module"] = module
            summaries.append(f"Added {_module_receipt_label(module_kind)} module.")
            continue
        if action.name == "calculation.set_condition":
            session, module = _target_module(owner, context, action.args)
            _set_module_condition(owner, session, module, action.args)
            summaries.append(f"Updated condition for {module.name}.")
            continue
        if action.name == "calculation.set_units":
            session, module = _target_module(owner, context, action.args)
            _set_module_units(owner, session, module, action.args)
            summaries.append(f"Updated units for {module.name}.")
            continue
        if action.name == "calculation.set_type":
            _session, module = _target_module(owner, context, action.args)
            _set_module_type(owner, module, action.args)
            summaries.append(f"Updated type for {module.name}.")
            continue
        if action.name in {
            "calculation.set_batch_condition",
            "calculation.set_batch_grid",
            "solidification.set_batch_grid",
        }:
            session, module = _target_module(owner, context, action.args)
            rows = _set_module_batch_condition(owner, session, module, action.args)
            summaries.append(
                f"Updated batch condition for {module.name}: {rows} row(s)."
            )
            continue
        if action.name == "calculation.set_nucleation_undercooling":
            _session, module = _target_module(owner, context, action.args)
            count = _set_module_nucleation_undercooling(owner, module, action.args)
            summaries.append(
                f"Updated nucleation undercooling for {module.name}: {count} phase(s)."
            )
            continue
        if action.name == "calculation.set_result_columns":
            _session, module = _target_module(owner, context, action.args)
            count = _set_module_result_columns(module, action.args)
            summaries.append(
                f"Updated result columns for {module.name}: {count} column(s)."
            )
            continue
        if action.name == "calculation.select_phases":
            session, module = _target_module(owner, context, action.args)
            selected = _set_module_phase_selection(owner, session, module, action.args)
            summaries.append(f"Selected {len(selected)} phase(s) for {module.name}.")
            continue
        if action.name == "calculation.calculate":
            session, module = _target_module(owner, context, action.args)
            _show_module(owner, session, module)
            if not _assistant_calculation_confirmed(action.args):
                summaries.append(
                    f"Prepared {module.name}. Review it, then press Calculate to run."
                )
                _notify_calculation_not_started(owner, module.name)
                continue
            owner._calculate_current_module()
            summaries.append(f"Started calculation for {module.name}.")
            continue
        if action.name == "calculation.calculate_all":
            if not _assistant_calculation_confirmed(action.args):
                summaries.append(
                    "Prepared calculation modules. "
                    "Review them, then press Calculate All to run."
                )
                _notify_calculation_not_started(owner, "all modules")
                continue
            owner._calculate_all_modules()
            summaries.append("Started all calculations.")
            continue
        raise AssistantActionError(f"Unsupported assistant action: {action.name}")
    return summaries


def _iter_action_payloads(text: str):
    cursor = 0
    source = str(text or "")
    while True:
        start = source.find(ACTION_BLOCK_START, cursor)
        if start < 0:
            return
        payload_start = start + len(ACTION_BLOCK_START)
        end = source.find(ACTION_BLOCK_END, payload_start)
        if end < 0:
            raise AssistantActionError("Assistant action block is missing its end tag.")
        yield source[payload_start:end].strip()
        cursor = end + len(ACTION_BLOCK_END)


def _create_calculation_session(
    owner: Any,
    args: dict[str, Any],
) -> CalculationSession:
    _switch_to_calculation_workspace(owner)
    name = str(args.get("name") or "").strip()
    if name:
        session = owner._create_calculation_session(name)
        owner._add_session_tree_item(session)
        owner._show_calculation_overview(session)
        return session
    return owner._add_calculation_session(False)


def _target_session(
    owner: Any,
    context: dict[str, Any],
    args: dict[str, Any] | None = None,
) -> CalculationSession:
    requested = _requested_session_text(args or {})
    if requested:
        sessions = getattr(owner, "calculation_sessions", {})
        for candidate in reversed(list(sessions.values())):
            if isinstance(candidate, CalculationSession) and _session_text_matches(
                candidate,
                requested,
            ):
                previous = context.get("session")
                context["session"] = candidate
                if previous is not candidate:
                    context.pop("module", None)
                return candidate
        raise AssistantActionError(f"Could not find calculation session {requested}.")

    session = context.get("session")
    if isinstance(session, CalculationSession):
        return session
    selected = owner._selected_session_for_new_module()
    if selected is not None:
        context["session"] = selected
        return selected
    return _create_calculation_session(owner, {})


def _target_module(
    owner: Any,
    context: dict[str, Any],
    args: dict[str, Any],
) -> tuple[CalculationSession, CalculationModule]:
    requested_session = _requested_session_text(args)
    session = (
        _target_session(owner, context, args)
        if requested_session
        else context.get("session")
    )
    module = context.get("module")
    if (
        isinstance(session, CalculationSession)
        and isinstance(module, CalculationModule)
        and module.id in session.modules
        and _module_matches(module, args)
    ):
        return session, module

    current_context = getattr(owner, "_current_calculation_context", None)
    if not requested_session and callable(current_context):
        current_session, current_module = current_context()
        if (
            isinstance(current_session, CalculationSession)
            and isinstance(current_module, CalculationModule)
            and _module_matches(current_module, args)
        ):
            context["session"] = current_session
            context["module"] = current_module
            return current_session, current_module

    if not isinstance(session, CalculationSession):
        session = _target_session(owner, context, args)
    modules = list(session.modules.values())
    if not modules:
        raise AssistantActionError("Create or select a calculation module first.")

    requested = _requested_module_text(args)
    if requested:
        for candidate in reversed(modules):
            if _module_text_matches(candidate, requested):
                context["session"] = session
                context["module"] = candidate
                return session, candidate
        raise AssistantActionError(f"Could not find calculation module {requested}.")

    module = modules[-1]
    context["session"] = session
    context["module"] = module
    return session, module


def _requested_module_text(args: dict[str, Any]) -> str:
    for key in ("module", "module_name", "name", "id", "kind"):
        value = str(args.get(key) or "").strip()
        if value:
            return value
    return ""


def _requested_session_text(args: dict[str, Any]) -> str:
    for key in ("session", "session_name", "session_id"):
        value = str(args.get(key) or "").strip()
        if value:
            return value
    return ""


def _session_text_matches(session: CalculationSession, requested: str) -> bool:
    normalized = requested.strip().lower().replace("-", "_").replace(" ", "_")
    if session.id.lower() == normalized:
        return True
    return (
        session.name.strip().lower().replace("-", "_").replace(" ", "_")
        == normalized
    )


def _module_matches(module: CalculationModule, args: dict[str, Any]) -> bool:
    requested = _requested_module_text(args)
    return not requested or _module_text_matches(module, requested)


def _module_text_matches(module: CalculationModule, requested: str) -> bool:
    normalized = requested.strip().lower().replace("-", "_").replace(" ", "_")
    if module.id.lower() == normalized:
        return True
    if module.name.strip().lower().replace("-", "_").replace(" ", "_") == normalized:
        return True
    try:
        return module.kind == _normalize_module_kind(requested)
    except AssistantActionError:
        return module.kind.lower() == normalized


def _switch_to_calculation_workspace(owner: Any) -> None:
    handler = getattr(owner, "_handle_workspace_mode_clicked", None)
    if callable(handler):
        handler(0)
        return
    workspace_stack = getattr(owner, "workspace_stack", None)
    if workspace_stack is not None:
        workspace_stack.setCurrentIndex(0)


def _show_module(
    owner: Any,
    session: CalculationSession,
    module: CalculationModule,
) -> None:
    show = getattr(owner, "_show_calculation_module", None)
    if callable(show):
        show(session, module)


def _assistant_calculation_confirmed(args: dict[str, Any]) -> bool:
    """Return whether an assistant action may launch a calculation."""
    for key in (
        "confirmed",
        "user_confirmed",
        "explicit",
        "run",
        "start",
        "launch",
        "calculate",
    ):
        if _truthy_action_arg(args.get(key)):
            return True
    return False


def _truthy_action_arg(value: Any) -> bool:
    """Interpret common structured-response truth values."""
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return value != 0
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "y", "run", "start"}
    return False


def _notify_calculation_not_started(owner: Any, target: str) -> None:
    """Leave a non-blocking GUI trace when assistant setup stops before run."""
    message = (
        f"Assistant prepared {target}; calculation was not launched without "
        "explicit confirmation."
    )
    logger = getattr(owner, "_append_calculation_log_message", None)
    if callable(logger):
        logger(message)
    status = getattr(owner, "_set_status_bar_message", None)
    if callable(status):
        status(message, level="warning")


def _set_module_condition(
    owner: Any,
    session: CalculationSession,
    module: CalculationModule,
    args: dict[str, Any],
) -> None:
    runtime = owner._ensure_calculation_module_runtime(session, module)
    _set_runtime_line(runtime, "temperature", _first_arg(args, "T", "temperature"))
    _set_runtime_line(runtime, "pressure", _first_arg(args, "P", "pressure"))
    _set_runtime_line(runtime, "delta_t", _first_arg(args, "delta_t", "dT"))
    _set_runtime_line(
        runtime,
        "liquid_phase",
        _first_arg(args, "liquid_phase", "liquid"),
    )

    composition = _composition_rows_from_args(args)
    table = runtime.get("composition_table")
    if composition is not None:
        if not _is_composition_table(table):
            raise AssistantActionError("The target module has no composition table.")
        _set_composition_table_rows(table, composition)

    transition_search = _first_arg(args, "transition_search", "transitions")
    widget = runtime.get("transition_search")
    if transition_search is not None and isinstance(widget, QCheckBox):
        module.transition_search = bool(transition_search)
        widget.setChecked(module.transition_search)

    from_liquidus = _first_arg(args, "from_liquidus", "start_from_liquidus")
    widget = runtime.get("from_liquidus")
    if from_liquidus is not None and isinstance(widget, QCheckBox):
        module.start_from_liquidus = bool(from_liquidus)
        widget.setChecked(module.start_from_liquidus)

    status = runtime.get("status")
    validator = getattr(owner, "_validate_composition_table", None)
    if callable(validator) and _is_composition_table(table):
        validator(
            session,
            table,
            status if isinstance(status, QPlainTextEdit) else None,
            quiet=True,
        )
    updater = getattr(owner, "_update_manual_batch_point_count", None)
    if callable(updater):
        updater(module)
    refresher = getattr(owner, "_refresh_nucleation_undercooling_runtime", None)
    if callable(refresher):
        refresher(module)


def _set_runtime_line(
    runtime: dict[str, Any],
    key: str,
    value: Any,
) -> None:
    if value is None:
        return
    widget = runtime.get(key)
    if isinstance(widget, QLineEdit):
        widget.setText(_action_value_text(value))


def _first_arg(args: dict[str, Any], *keys: str) -> Any:
    for key in keys:
        if key in args:
            return args[key]
    return None


def _action_value_text(value: Any) -> str:
    if isinstance(value, (list, tuple)):
        return " ".join(_action_value_text(item) for item in value)
    return str(value)


def _composition_rows_from_args(args: dict[str, Any]) -> list[dict[str, str]] | None:
    if "composition" in args:
        raw = args["composition"]
    elif "species" in args:
        raw = args["species"]
    else:
        return None

    if isinstance(raw, dict):
        return [
            {"species": str(species), "amount": _action_value_text(amount)}
            for species, amount in raw.items()
        ]
    if isinstance(raw, list):
        rows: list[dict[str, str]] = []
        for entry in raw:
            if isinstance(entry, dict):
                species = (
                    entry.get("species") or entry.get("name") or entry.get("element")
                )
                amount = entry.get("amount")
                if amount is None:
                    amount = entry.get("value")
            elif isinstance(entry, (list, tuple)) and len(entry) >= 2:
                species, amount = entry[0], entry[1]
            else:
                raise AssistantActionError(
                    "Composition entries must name species and amount."
                )
            if not str(species or "").strip():
                raise AssistantActionError(
                    "Composition entries must include species names."
                )
            rows.append({"species": str(species), "amount": _action_value_text(amount)})
        return rows
    raise AssistantActionError("Composition must be an object or a list of entries.")


def _set_module_units(
    owner: Any,
    session: CalculationSession,
    module: CalculationModule,
    args: dict[str, Any],
) -> None:
    temperature = _canonical_choice(
        _first_arg(args, "temperature", "temperature_unit", "T"),
        {"K", "C", "F", "R"},
        module.temperature_unit,
        "temperature unit",
    )
    pressure = _canonical_choice(
        _first_arg(args, "pressure", "pressure_unit", "P"),
        {"atm", "psi", "bar", "Pa", "kPa"},
        module.pressure_unit,
        "pressure unit",
    )
    amount = _canonical_choice(
        _first_arg(args, "amount", "amount_unit", "composition"),
        set(AMOUNT_UNIT_OPTIONS),
        module.amount_unit,
        "amount unit",
    )
    if str(args.get("scope") or "").strip().lower() == "session":
        setter = getattr(owner, "_set_session_default_units", None)
        if callable(setter):
            setter(session, temperature, pressure, amount)
        else:
            session.temperature_unit = temperature
            session.pressure_unit = pressure
            session.amount_unit = amount
        return

    module.temperature_unit = temperature
    module.pressure_unit = pressure
    module.amount_unit = amount
    if module.runtime:
        owner._apply_calculation_type_to_runtime(module)
    header = getattr(owner, "_set_header_unit_controls", None)
    if callable(header):
        header(module)


def _canonical_choice(
    value: Any,
    allowed: set[str],
    fallback: str,
    label: str,
) -> str:
    if value is None:
        return fallback
    text = str(value).strip()
    for option in allowed:
        if option.lower() == text.lower():
            return option
    raise AssistantActionError(f"Unsupported {label}: {value}")


def _set_module_type(
    owner: Any,
    module: CalculationModule,
    args: dict[str, Any],
) -> None:
    mode = str(
        _first_arg(args, "mode", "type", "calculation_type") or module.calculation_type
    )
    normalized = mode.strip().lower().replace("-", "_").replace(" ", "_")
    aliases = {
        "single": "single",
        "single_condition": "single",
        "batch": "batch",
        "batch_condition": "batch",
    }
    if normalized not in aliases:
        raise AssistantActionError(f"Unsupported calculation type: {mode}")
    previous_signature = (
        module.calculation_type,
        module.solidification_model,
        module.batch_condition_path,
    )
    module.calculation_type = aliases[normalized]
    if module.calculation_type == "single":
        module.batch_condition = {}
        module.batch_condition_path = ""
    if module.kind == "solidification":
        module.solidification_model = _solidification_model_value(
            str(
                _first_arg(args, "model", "solidification_model")
                or module.solidification_model
            )
        )
    cpu_count = _first_arg(args, "cpu_count", "cpus")
    if cpu_count is not None:
        try:
            module.batch_cpu_count = max(0, int(cpu_count))
        except (TypeError, ValueError) as exc:
            raise AssistantActionError(
                f"CPU count must be an integer: {cpu_count}"
            ) from exc
    current_signature = (
        module.calculation_type,
        module.solidification_model,
        module.batch_condition_path,
    )
    if current_signature != previous_signature:
        module.phase_names = []
        phase_tree = module.runtime.get("phase_tree")
        if phase_tree is not None:
            phase_tree.clear()
    if module.runtime:
        owner._apply_calculation_type_to_runtime(module)


def _set_module_batch_condition(
    owner: Any,
    session: CalculationSession,
    module: CalculationModule,
    args: dict[str, Any],
) -> int:
    if module.kind not in {"equilibrium", "solidification"}:
        raise AssistantActionError(f"{module.name} does not support batch conditions.")
    active_database = owner._active_session_database(session)
    if active_database is None:
        raise AssistantActionError("Load a database before setting batch conditions.")

    condition = _batch_condition_from_args(
        args,
        module,
        active_database.database,
    )
    row_count = _condition_row_count(condition)
    module.calculation_type = "batch"
    module.batch_condition = condition
    module.batch_condition_path = str(args.get("label") or "Assistant generated")
    if "cpu_count" in args or "cpus" in args:
        cpu_count = _first_arg(args, "cpu_count", "cpus")
        try:
            module.batch_cpu_count = max(0, int(cpu_count))
        except (TypeError, ValueError) as exc:
            raise AssistantActionError(
                f"CPU count must be an integer: {cpu_count}"
            ) from exc

    runtime = owner._ensure_calculation_module_runtime(session, module)
    _set_runtime_line(runtime, "temperature", _first_arg(args, "T", "temperature"))
    _set_runtime_line(runtime, "pressure", _first_arg(args, "P", "pressure"))
    _set_runtime_line(runtime, "delta_t", _first_arg(args, "delta_t", "dT"))
    _set_runtime_line(
        runtime,
        "liquid_phase",
        _first_arg(args, "liquid_phase", "liquid"),
    )
    from_liquidus = _first_arg(args, "from_liquidus", "start_from_liquidus")
    widget = runtime.get("from_liquidus")
    if from_liquidus is not None and isinstance(widget, QCheckBox):
        module.start_from_liquidus = bool(from_liquidus)
        widget.setChecked(module.start_from_liquidus)

    table = runtime.get("composition_table")
    composition = _composition_rows_from_args(args)
    if composition is not None and _is_composition_table(table):
        _set_composition_table_rows(table, composition)
    if module.runtime:
        owner._apply_calculation_type_to_runtime(module)
    updater = getattr(owner, "_update_manual_batch_point_count", None)
    if callable(updater):
        updater(module)
    return row_count


def _batch_condition_from_args(
    args: dict[str, Any],
    module: CalculationModule,
    database: Any,
) -> dict[str, list[float]]:
    raw = _first_arg(args, "condition", "conditions", "batch_condition")
    if raw is not None:
        return _normalize_batch_condition(raw, database, module.amount_unit)
    return _batch_grid_condition(args, module, database)


def _normalize_batch_condition(
    raw: Any,
    database: Any,
    amount_unit: str,
) -> dict[str, list[float]]:
    if not isinstance(raw, dict):
        raise AssistantActionError("Batch condition must be an object.")
    condition: dict[str, list[float]] = {}
    for key, values in raw.items():
        name = str(key).strip()
        if not name:
            raise AssistantActionError("Batch condition contains an empty column name.")
        if isinstance(values, (list, tuple)):
            column = [_float_arg(value, f"condition column {name}") for value in values]
        else:
            column = [_float_arg(values, f"condition column {name}")]
        condition[name] = column
    _condition_row_count(condition)
    return _expand_batch_species(database, condition, amount_unit)


def _batch_grid_condition(
    args: dict[str, Any],
    module: CalculationModule,
    database: Any,
) -> dict[str, list[float]]:
    axes = _first_arg(args, "axes", "grid", "ranges")
    if not isinstance(axes, dict) or not axes:
        raise AssistantActionError(
            "Batch grid needs args.axes, for example "
            '{"Cu":[0,100,10],"Mg":[0,100,10]}.'
        )
    balance = str(_first_arg(args, "balance", "balance_species") or "").strip()
    if not balance:
        raise AssistantActionError("Batch grid needs a balance species.")
    fixed = _first_arg(args, "fixed", "constants")
    if fixed is None:
        fixed = {}
    if not isinstance(fixed, dict):
        raise AssistantActionError("Batch grid fixed values must be an object.")

    total = _first_arg(args, "total", "target_total")
    target_total = (
        _float_arg(total, "batch grid total")
        if total is not None
        else _default_batch_grid_total(module.amount_unit)
    )
    if target_total <= 0:
        raise AssistantActionError("Batch grid total must be positive.")

    temperatures = _condition_values(
        _first_arg(args, "T", "temperature"),
        "temperature",
        fallback=_runtime_line_float(module, "temperature", 1600.0),
    )
    pressures = _condition_values(
        _first_arg(args, "P", "pressure"),
        "pressure",
        fallback=_runtime_line_float(module, "pressure", 1.0),
    )
    if len(pressures) not in {1, len(temperatures)}:
        raise AssistantActionError("Pressure must be one value or match temperatures.")

    axis_names = [str(name).strip() for name in axes]
    if any(not name for name in axis_names):
        raise AssistantActionError("Batch grid contains an empty axis name.")
    axis_values = [
        _condition_values(values, f"axis {name}") for name, values in axes.items()
    ]
    fixed_values = {
        str(name).strip(): _float_arg(value, f"fixed value {name}")
        for name, value in fixed.items()
        if str(name).strip()
    }
    duplicate_axes = sorted(set(axis_names) & set(fixed_values))
    if duplicate_axes:
        raise AssistantActionError(
            "Batch grid species cannot be both axis and fixed: "
            + ", ".join(duplicate_axes)
        )
    if balance in axis_names or balance in fixed_values:
        raise AssistantActionError(
            "Balance species cannot also be an axis or fixed species."
        )

    species_order = [balance, *axis_names, *fixed_values]
    condition: dict[str, list[float]] = {"T": [], "P": []}
    for species in species_order:
        condition[species] = []

    for temperature_index, temperature in enumerate(temperatures):
        pressure = pressures[temperature_index] if len(pressures) > 1 else pressures[0]
        for combination in product(*axis_values):
            known_total = sum(combination) + sum(fixed_values.values())
            if known_total > target_total and not _nearly_equal(
                known_total, target_total
            ):
                continue
            balance_value = target_total - known_total
            if balance_value < 0 and _nearly_equal(balance_value, 0.0):
                balance_value = 0.0
            if balance_value < 0:
                continue
            condition["T"].append(float(temperature))
            condition["P"].append(float(pressure))
            condition[balance].append(float(balance_value))
            for species, amount in zip(axis_names, combination, strict=False):
                condition[species].append(float(amount))
            for species, amount in fixed_values.items():
                condition[species].append(float(amount))

    if not condition["T"]:
        raise AssistantActionError("Batch grid produced no valid rows.")
    return _expand_batch_species(database, condition, module.amount_unit)


def _condition_values(
    value: Any, label: str, fallback: float | None = None
) -> list[float]:
    if value is None:
        if fallback is None:
            raise AssistantActionError(f"{label} is required.")
        return [float(fallback)]
    if isinstance(value, dict):
        value = [
            value.get("start"),
            value.get("stop", value.get("end")),
            value.get("step"),
        ]
    if isinstance(value, str):
        stripped = value.strip().strip("[]")
        if not stripped:
            raise AssistantActionError(f"{label} cannot be empty.")
        parts = [part for part in stripped.replace(",", " ").split() if part]
        return _range_or_values(parts, label)
    if isinstance(value, (list, tuple)):
        return _range_or_values(list(value), label)
    return [_float_arg(value, label)]


def _range_or_values(values: list[Any], label: str) -> list[float]:
    if len(values) == 3:
        start = _float_arg(values[0], label)
        stop = _float_arg(values[1], label)
        step = _float_arg(values[2], label)
        return _inclusive_range(start, stop, step, label)
    if not values:
        raise AssistantActionError(f"{label} cannot be empty.")
    return [_float_arg(value, label) for value in values]


def _inclusive_range(start: float, stop: float, step: float, label: str) -> list[float]:
    if step == 0:
        return [float(start)]
    if start == stop:
        return [float(start)]
    if (stop - start) * step < 0:
        raise AssistantActionError(f"{label} step must move from start toward stop.")
    values = [float(start)]
    current = start + step
    if step > 0:
        while current < stop and not _nearly_equal(current, stop):
            values.append(float(current))
            current += step
    else:
        while current > stop and not _nearly_equal(current, stop):
            values.append(float(current))
            current += step
    if not _nearly_equal(values[-1], stop):
        values.append(float(stop))
    return values


def _expand_batch_species(
    database: Any,
    condition: dict[str, list[float]],
    amount_unit: str,
) -> dict[str, list[float]]:
    try:
        expanded = _expand_condition_species(database, condition, amount_unit)
    except Exception as exc:
        raise AssistantActionError(
            f"Could not expand batch condition species: {exc}"
        ) from exc
    return {
        str(key): _float_column(values)
        for key, values in expanded.items()
    }


def _condition_row_count(condition: dict[str, list[float]]) -> int:
    if not condition:
        raise AssistantActionError("Batch condition is empty.")
    lengths = {len(values) for values in condition.values()}
    if len(lengths) != 1:
        raise AssistantActionError("Batch condition columns must have the same length.")
    return lengths.pop()


def _float_column(values: Any) -> list[float]:
    if isinstance(values, (list, tuple)):
        return [float(value) for value in values]
    try:
        return [float(value) for value in values.tolist()]
    except AttributeError:
        return [float(values)]


def _default_batch_grid_total(amount_unit: str) -> float:
    normalized = str(amount_unit or "").strip().lower()
    if normalized in {"mol%", "wt%", "at%", "mass%"} or "%" in normalized:
        return 100.0
    if normalized in {"mole fraction", "weight fraction", "fraction"}:
        return 1.0
    return 1.0


def _runtime_line_float(
    module: CalculationModule,
    key: str,
    fallback: float,
) -> float:
    widget = module.runtime.get(key)
    if isinstance(widget, QLineEdit):
        text = widget.text().strip()
        try:
            return float(text.split()[0])
        except (IndexError, ValueError):
            return fallback
    return fallback


def _float_arg(value: Any, label: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise AssistantActionError(f"{label} must be numeric: {value}") from exc


def _nearly_equal(left: float, right: float) -> bool:
    return abs(left - right) <= 1e-10


def _set_module_nucleation_undercooling(
    owner: Any,
    module: CalculationModule,
    args: dict[str, Any],
) -> int:
    if module.kind != "solidification":
        raise AssistantActionError(
            "Nucleation undercooling is only available for solidification."
        )
    raw = _first_arg(args, "undercooling", "critical_undercooling", "values", "phases")
    if not isinstance(raw, dict) or not raw:
        raise AssistantActionError(
            "Nucleation undercooling must be a phase-to-value object."
        )
    values = {
        str(phase).strip(): _float_arg(value, f"undercooling for {phase}")
        for phase, value in raw.items()
        if str(phase).strip()
    }
    if not values:
        raise AssistantActionError("Nucleation undercooling has no phase names.")
    module.nucleation_undercooling.update(values)
    refresher = getattr(owner, "_refresh_nucleation_undercooling_runtime", None)
    if callable(refresher):
        refresher(module)
    tree = module.runtime.get("nucleation_undercooling_tree")
    if isinstance(tree, QTreeWidget):
        for parent_row in range(tree.topLevelItemCount()):
            parent = tree.topLevelItem(parent_row)
            for row in range(parent.childCount()):
                child = parent.child(row)
                phase = child.text(0).strip()
                if phase in values:
                    child.setText(1, str(values[phase]))
    return len(values)


def _set_module_result_columns(
    module: CalculationModule,
    args: dict[str, Any],
) -> int:
    raw = _first_arg(args, "columns", "result_columns", "selected")
    if raw is None:
        module.result_columns = []
        return 0
    if isinstance(raw, str):
        columns = [
            item.strip() for item in raw.replace(";", ",").split(",") if item.strip()
        ]
    elif isinstance(raw, list):
        columns = [str(item).strip() for item in raw if str(item).strip()]
    else:
        raise AssistantActionError("Result columns must be a string or list.")
    module.result_columns = columns
    return len(columns)


def _set_module_phase_selection(
    owner: Any,
    session: CalculationSession,
    module: CalculationModule,
    args: dict[str, Any],
) -> list[str]:
    runtime = owner._ensure_calculation_module_runtime(session, module)
    phase_tree = runtime.get("phase_tree")
    requested = _phase_list_from_args(args)
    available = _available_phase_names(owner, session, module, phase_tree)
    if not available:
        raise AssistantActionError("No phases are available for this module yet.")

    if not requested or _select_all_requested(args):
        selected = list(available)
    else:
        selected = _resolve_requested_phases(requested, available)

    module.phase_names = list(available)
    if phase_tree is not None:
        if not _all_phase_names(phase_tree):
            _populate_phase_tree(phase_tree, available, selected)
        else:
            _set_phase_tree_selection(phase_tree, selected)
        selected = _selected_phase_names(phase_tree)
    return selected


def _phase_list_from_args(args: dict[str, Any]) -> list[str]:
    raw = args.get("phases")
    if raw is None:
        raw = args.get("phase_names")
    if raw is None:
        raw = args.get("selected")
    if raw is None:
        return []
    if isinstance(raw, str):
        return [
            item.strip() for item in raw.replace(";", ",").split(",") if item.strip()
        ]
    if isinstance(raw, list):
        return [str(item).strip() for item in raw if str(item).strip()]
    raise AssistantActionError("Phase selection must be a string or list.")


def _select_all_requested(args: dict[str, Any]) -> bool:
    value = args.get("mode") or args.get("selection")
    return str(value or "").strip().lower() in {"all", "select_all", "*"}


def _available_phase_names(
    owner: Any,
    session: CalculationSession,
    module: CalculationModule,
    phase_tree: Any,
) -> list[str]:
    existing = _all_phase_names(phase_tree) if phase_tree is not None else []
    if existing:
        return existing
    if module.phase_names:
        return list(module.phase_names)
    active_database = owner._active_session_database(session)
    if active_database is None:
        raise AssistantActionError("Load a database before selecting phases.")
    phase_names = _database_phase_names(active_database.database)
    module.phase_names = list(phase_names)
    return phase_names


def _resolve_requested_phases(
    requested: list[str],
    available: list[str],
) -> list[str]:
    by_lower = {phase.lower(): phase for phase in available}
    selected: list[str] = []
    missing: list[str] = []
    for phase in requested:
        match = by_lower.get(phase.lower())
        if match is None:
            missing.append(phase)
        elif match not in selected:
            selected.append(match)
    if missing:
        choices = ", ".join(available[:8])
        if len(available) > 8:
            choices += f", ... {len(available) - 8} more"
        raise AssistantActionError(
            f"Unknown phase(s): {', '.join(missing)}. Available: {choices}"
        )
    return selected


def _load_database_from_project_directory(
    owner: Any,
    session: CalculationSession,
    args: dict[str, Any],
) -> CalculationDatabase:
    path = _database_path_from_project_directory(owner, args)
    database_name = path.name
    duplicate = _find_duplicate_database(session, path, database_name)
    if duplicate is not None:
        for loaded in session.databases.values():
            loaded.selected = loaded.id == duplicate.id
        _refresh_session_after_database_change(owner, session)
        return duplicate

    runtime_database = _runtime_database_from_path(path)
    session.database_counter += 1
    database_id = f"{session.id}:database:{session.database_counter}"
    for loaded in session.databases.values():
        loaded.selected = False
    loaded = CalculationDatabase(
        id=database_id,
        name=database_name,
        path=str(path),
        database=runtime_database,
        selected=True,
    )
    session.databases[database_id] = loaded
    _refresh_session_after_database_change(owner, session)
    return loaded


def _database_path_from_project_directory(
    owner: Any,
    args: dict[str, Any],
) -> Path:
    active_directory = getattr(owner, "active_calculation_directory", None)
    if active_directory is None:
        _prompt_for_active_calculation_directory(owner)
        active_directory = getattr(owner, "active_calculation_directory", None)
        if active_directory is None:
            raise AssistantActionError("Active directory was not set.")
    root = Path(active_directory).expanduser().resolve()
    preferred_dir = root / "database"
    preferred_dir.mkdir(parents=True, exist_ok=True)
    database_dirs = [preferred_dir]
    compatibility_dir = root / "databases"
    if compatibility_dir.exists() and compatibility_dir != preferred_dir:
        database_dirs.append(compatibility_dir)

    requested = str(args.get("file") or args.get("name") or "").strip()
    files = _project_database_files(database_dirs)
    if requested:
        requested_path = _resolve_requested_database_file(database_dirs, requested)
        if requested_path is None:
            raise AssistantActionError(
                f"Could not find {requested} in {preferred_dir}."
            )
        return requested_path
    if not files:
        raise AssistantActionError(
            f"Created database folder {preferred_dir}. Put one .dat or .tdb "
            "database there, then ask me to load it again."
        )
    if len(files) > 1:
        names = ", ".join(path.name for path in files[:8])
        if len(files) > 8:
            names += f", ... {len(files) - 8} more"
        raise AssistantActionError(
            "Multiple databases are available. Ask me to load one by file name: "
            f"{names}"
        )
    return files[0]


def _prompt_for_active_calculation_directory(owner: Any) -> None:
    """Ask the user whether to open the active-directory chooser."""
    parent = owner if isinstance(owner, QWidget) else None
    result = QMessageBox.warning(
        parent,
        "Set active directory",
        (
            "Equilipy needs an active calculation directory before it can load "
            "a database. Set the active directory now?"
        ),
        QMessageBox.StandardButton.Ok | QMessageBox.StandardButton.Cancel,
        QMessageBox.StandardButton.Ok,
    )
    if result != QMessageBox.StandardButton.Ok:
        return
    setter = getattr(owner, "_set_active_calculation_directory", None)
    if callable(setter):
        setter()


def _project_database_files(database_dirs: list[Path]) -> list[Path]:
    files: list[Path] = []
    seen: set[Path] = set()
    for database_dir in database_dirs:
        if not database_dir.exists():
            continue
        for path in sorted(database_dir.iterdir(), key=lambda item: item.name.lower()):
            if path.suffix.lower() not in {".dat", ".tdb"} or not path.is_file():
                continue
            resolved = path.resolve()
            if resolved not in seen:
                seen.add(resolved)
                files.append(resolved)
    return files


def _resolve_requested_database_file(
    database_dirs: list[Path],
    requested: str,
) -> Path | None:
    requested_path = Path(requested).expanduser()
    if requested_path.is_absolute() and requested_path.is_file():
        if requested_path.suffix.lower() in {".dat", ".tdb"}:
            return requested_path.resolve()
        return None
    for database_dir in database_dirs:
        for candidate in (
            database_dir / requested,
            database_dir / f"{requested}.dat",
            database_dir / f"{requested}.tdb",
        ):
            if candidate.is_file() and candidate.suffix.lower() in {".dat", ".tdb"}:
                return candidate.resolve()
    normalized = requested.lower()
    for path in _project_database_files(database_dirs):
        if path.name.lower() == normalized or path.stem.lower() == normalized:
            return path
    return None


def _find_duplicate_database(
    session: CalculationSession,
    path: Path | None,
    name: str,
) -> CalculationDatabase | None:
    for loaded in session.databases.values():
        if path is not None:
            try:
                if Path(loaded.path).expanduser().resolve() == path:
                    return loaded
            except OSError:
                pass
        if loaded.name == name:
            return loaded
    return None


def _runtime_database_from_path(path: Path) -> Any:
    suffix = path.suffix.lower()
    try:
        import equilipy as eq

        if suffix == ".dat":
            return eq.read_dat(str(path))
        if suffix == ".tdb":
            return eq.read_tdb(str(path))
    except Exception as exc:
        raise AssistantActionError(f"Could not load runtime database: {exc}") from exc
    raise AssistantActionError(f"Unsupported database file type: {path.name}")


def _refresh_session_after_database_change(
    owner: Any, session: CalculationSession
) -> None:
    owner._show_calculation_overview(session)
    owner._update_current_calculation_heading(session)
    revalidate = getattr(owner, "_revalidate_session_composition_tables", None)
    if callable(revalidate):
        revalidate(session)


def _normalize_module_kind(value: object) -> str:
    normalized = str(value or "").strip().lower().replace("-", "_").replace(" ", "_")
    aliases = {
        "equilib": "equilibrium",
        "equilibrium": "equilibrium",
        "single_equilibrium": "equilibrium",
        "solidification": "solidification",
        "solidify": "solidification",
        "scheil": "solidification",
        "nucleoscheil": "solidification",
        "nucleo_scheil": "solidification",
    }
    try:
        return aliases[normalized]
    except KeyError as exc:
        raise AssistantActionError(
            f"Unsupported calculation module kind: {value}"
        ) from exc


def _module_receipt_label(module_kind: str) -> str:
    labels = {
        "equilibrium": "Equilib",
        "solidification": "Solidification",
    }
    return labels.get(module_kind, module_kind.title())
