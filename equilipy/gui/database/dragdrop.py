"""Copy-only database tree drag/drop behavior.

This module owns object-level copy rules for the database sidebar. Keep future
merge behavior close to these routines so drag/drop and explicit merge workflows
share the same semantics.
"""

from __future__ import annotations

import copy
import re
from typing import Any

from PySide6.QtCore import QModelIndex
from PySide6.QtWidgets import QMessageBox

from equilipy.database_ir import DatabaseIR, FunctionDefinition, Parameter, Phase
from equilipy.gui.database.cp_editor import (
    _direct_function_gibbs_rows,
    _function_by_name,
    _function_thermo_data,
    _normalize_function_temperature_ranges,
    _normalized_gibbs_rows,
    _referenced_function_names,
    _set_function_gibbs_rows,
)
from equilipy.gui.database.overview import (
    _compound_parameter_for_range_payload,
    _compound_parameter_target_species,
    _compound_phase_by_label,
    _compound_phase_display_label,
    _compound_phase_g_parameters,
    _compound_phase_state_id,
    _new_compound_g_parameter,
    _new_compound_phase_state_id,
    _parameter_as_function_expression,
    _set_parameter_gibbs_rows,
)
from equilipy.gui.models import (
    CompoundPhasePayload,
    CompoundRecord,
    SolutionParameterPayload,
    ThermoRangePayload,
    _is_solution_phase,
)


class DatabaseDragDropCopyMixin:
    """Copy-only drag/drop operations for DatabaseIR tree payloads."""

    def _copy_database_tree_payload(
        self,
        source_index: QModelIndex,
        target_index: QModelIndex,
    ) -> bool:
        """Copy database tree payloads by drag-and-drop without moving them."""
        if not source_index.isValid() or not target_index.isValid():
            return False
        commit_editor = getattr(self, "_commit_active_database_thermo_editor", None)
        if callable(commit_editor):
            commit_editor()
        source_payload = self.tree_model.payload(source_index)
        source_database = self.tree_model.database_for(source_index)
        target_database = self.tree_model.database_for(target_index)

        if isinstance(source_payload, FunctionDefinition):
            target_payload = self.tree_model.payload(target_index)
            if isinstance(target_payload, CompoundPhasePayload):
                return self._replace_compound_phase_expression_from_function(
                    source_database,
                    source_payload,
                    target_database,
                    target_payload,
                )
            if isinstance(target_payload, SolutionParameterPayload):
                return self._replace_solution_parameter_expression_from_function(
                    source_database,
                    source_payload,
                    target_database,
                    target_payload,
                )
            return self._copy_function_to_database(
                source_database,
                source_payload,
                target_database,
            )
        if isinstance(source_payload, CompoundRecord):
            return self._copy_compound_to_database(
                source_database,
                source_payload,
                target_database,
            )
        if isinstance(source_payload, Phase) and _is_solution_phase(source_payload):
            return self._copy_solution_to_database(
                source_database,
                source_payload,
                target_database,
            )
        if isinstance(source_payload, CompoundPhasePayload):
            target_payload = self.tree_model.payload(target_index)
            if isinstance(target_payload, SolutionParameterPayload):
                return self._replace_solution_parameter_expression_from_compound_phase(
                    source_database,
                    source_payload,
                    target_database,
                    target_payload,
                )
            target_function = self._target_function_for_drop(
                target_database,
                target_index,
            )
            if target_function is not None:
                return self._replace_function_expression_from_compound_phase(
                    source_database,
                    source_payload,
                    target_database,
                    target_function,
                )
            if self._is_function_collection_drop_target(target_index):
                return self._copy_compound_phase_to_function_collection(
                    source_database,
                    source_payload,
                    target_database,
                )
            target_record = self._target_compound_record_for_drop(target_index)
            if target_record is None:
                return False
            return self._copy_compound_phase_to_record(
                source_database,
                source_payload,
                target_database,
                target_record,
            )
        if isinstance(source_payload, SolutionParameterPayload):
            if source_payload.kind != "Endmember":
                return False
            target_function = self._target_function_for_drop(
                target_database,
                target_index,
            )
            if target_function is not None:
                return self._replace_function_expression_from_solution_parameter(
                    source_database,
                    source_payload,
                    target_database,
                    target_function,
                )
            if self._is_function_collection_drop_target(target_index):
                return self._copy_solution_parameter_to_function_collection(
                    source_database,
                    source_payload,
                    target_database,
                )
            return False
        if isinstance(source_payload, ThermoRangePayload):
            return self._copy_thermo_range_to_target(
                source_database,
                source_payload,
                target_database,
                target_index,
            )
        return False

    def _copy_function_to_database(
        self,
        source_database: DatabaseIR,
        function: FunctionDefinition,
        target_database: DatabaseIR,
    ) -> bool:
        """Copy one function and all Cp ranges to a target database."""
        existing = _function_by_name(target_database, function.name)
        if existing is not None and not self._confirm_database_tree_replace(
            "function",
            function.name,
            target_database,
        ):
            return False
        copied = copy.deepcopy(function)
        if existing is not None:
            target_database.functions[
                target_database.functions.index(existing)
            ] = copied
        else:
            target_database.functions.append(copied)
        self._copy_referenced_functions_for_expressions(
            source_database,
            target_database,
            self._function_dependency_expressions(function),
            exclude_names={function.name},
        )
        self._rebuild_database_tree(target_database, select_payload=copied)
        return True

    def _copy_compound_to_database(
        self,
        source_database: DatabaseIR,
        record: CompoundRecord,
        target_database: DatabaseIR,
    ) -> bool:
        """Copy a compound record, including phases and Gibbs/Cp parameters."""
        existing = self._compound_record_by_name(target_database, record.name)
        if existing is not None and not self._confirm_database_tree_replace(
            "compound",
            record.name,
            target_database,
        ):
            return False
        if (
            source_database is target_database
            and existing is not None
            and self._compound_records_share_phase_objects(record, existing)
        ):
            self._rebuild_database_tree(target_database, show_summary=False)
            selected = self._compound_record_by_name(target_database, record.name)
            if selected is not None:
                index = self._find_database_tree_index(selected)
                if index.isValid():
                    self.tree.setCurrentIndex(index)
                    self.tree.scrollTo(index)
            return True

        replace_index = (
            self._compound_record_first_phase_index(target_database, existing)
            if existing is not None
            else None
        )
        phases = [copy.deepcopy(phase) for phase in record.phases]
        phase_id_map: dict[str, str] = {}
        for phase in phases:
            old_id = _compound_phase_state_id(phase)
            new_id = _new_compound_phase_state_id(target_database, phase.name)
            phase.metadata["phase_id"] = new_id
            if old_id:
                phase_id_map[old_id] = new_id

        source_parameters = [
            copy.deepcopy(parameter)
            for parameter in source_database.parameters
            if parameter.phase_name == record.name
        ]
        for parameter in source_parameters:
            old_id = str(parameter.metadata.get("phase_id", "")).strip()
            if old_id in phase_id_map:
                parameter.metadata["phase_id"] = phase_id_map[old_id]

        if existing is not None:
            self._remove_compound_record(target_database, existing)
        self._copy_referenced_functions_for_parameters(
            source_database,
            target_database,
            source_parameters,
        )
        self._copy_species_for_phases(source_database, target_database, phases)
        if replace_index is None:
            target_database.phases.extend(phases)
        else:
            insertion_index = min(replace_index, len(target_database.phases))
            target_database.phases[insertion_index:insertion_index] = phases
        target_database.parameters.extend(source_parameters)
        self._rebuild_database_tree(target_database, show_summary=False)
        selected = self._compound_record_by_name(target_database, record.name)
        if selected is not None:
            index = self._find_database_tree_index(selected)
            if index.isValid():
                self.tree.setCurrentIndex(index)
                self.tree.scrollTo(index)
        return True

    def _copy_solution_to_database(
        self,
        source_database: DatabaseIR,
        phase: Phase,
        target_database: DatabaseIR,
    ) -> bool:
        """Copy one solution phase and its parameters to a target database."""
        existing = next(
            (
                candidate
                for candidate in target_database.phases
                if _is_solution_phase(candidate)
                and candidate.name.upper() == phase.name.upper()
            ),
            None,
        )
        if existing is not None and not self._confirm_database_tree_replace(
            "solution",
            phase.name,
            target_database,
        ):
            return False

        copied_phase = copy.deepcopy(phase)
        copied_parameters = [
            copy.deepcopy(parameter)
            for parameter in source_database.parameters
            if parameter.phase_name == phase.name
        ]
        if existing is not None:
            target_database.phases = [
                candidate
                for candidate in target_database.phases
                if candidate is not existing
            ]
            target_database.parameters = [
                parameter
                for parameter in target_database.parameters
                if parameter.phase_name.upper() != phase.name.upper()
            ]
        self._copy_referenced_functions_for_parameters(
            source_database,
            target_database,
            copied_parameters,
        )
        self._copy_species_for_phases(source_database, target_database, [phase])
        target_database.phases.append(copied_phase)
        target_database.parameters.extend(copied_parameters)
        self._rebuild_database_tree(target_database, select_payload=copied_phase)
        return True

    def _replace_compound_phase_expression_from_function(
        self,
        source_database: DatabaseIR,
        function: FunctionDefinition,
        target_database: DatabaseIR,
        target_payload: CompoundPhasePayload,
    ) -> bool:
        """Replace a compound phase Gibbs expression with a Function expression."""
        expression = self._function_expression_for_compound_phase_drop(function)
        if not expression:
            return False
        target_phase = target_payload.phase
        parameters = _compound_phase_g_parameters(target_database, target_phase)
        target_parameter = parameters[0] if parameters else None
        if target_parameter is None:
            target_parameter = _new_compound_g_parameter(target_phase)
            target_database.parameters.append(target_parameter)
        self._set_compound_phase_parameter_expression(
            target_phase,
            target_parameter,
            expression,
        )
        self._copy_referenced_functions_for_expressions(
            source_database,
            target_database,
            [expression],
        )
        self._rebuild_database_tree(target_database, show_summary=False)
        self._select_compound_phase(target_phase)
        return True

    def _function_expression_for_compound_phase_drop(
        self,
        function: FunctionDefinition,
    ) -> str:
        """Return a Function expression suitable for a compound G parameter."""
        expression = str(function.expression or "").strip()
        if expression:
            return expression
        source_expression = str(function.tdb_src.source_expression or "").strip()
        if source_expression:
            return f"298.15 {source_expression}; 6000 N"
        return ""

    def _set_compound_phase_parameter_expression(
        self,
        phase: Phase,
        parameter: Parameter,
        expression: str,
    ) -> None:
        """Store a full Gibbs expression on a compound phase parameter only."""
        target = _compound_parameter_target_species(phase)
        parameter.phase_name = phase.name
        parameter.parameter_type = "G"
        parameter.target = target
        parameter.order = 0
        parameter.expression = expression
        parameter.metadata.pop("source_expression", None)
        phase_id = _compound_phase_state_id(phase)
        if phase_id:
            parameter.metadata["phase_id"] = phase_id
        target_text = ":".join(target)
        parameter.source.command = (
            f"PARAMETER G({phase.name},{target_text};0) {expression} !"
        )

    def _replace_solution_parameter_expression_from_function(
        self,
        source_database: DatabaseIR,
        function: FunctionDefinition,
        target_database: DatabaseIR,
        target_payload: SolutionParameterPayload,
    ) -> bool:
        """Replace a solution parameter expression with a Function reference."""
        parameter = target_payload.parameter
        reference = f"{function.name}#"
        expression = f"298.15 {reference}; 6000 N"
        parameter.expression = expression
        parameter.metadata["source_expression"] = reference
        parameter.source.command = self._solution_parameter_command_with_expression(
            parameter,
            expression,
        )
        self._copy_referenced_functions_for_expressions(
            source_database,
            target_database,
            [reference],
        )
        self._rebuild_database_tree(target_database, show_summary=False)
        return True

    def _replace_solution_parameter_expression_from_compound_phase(
        self,
        source_database: DatabaseIR,
        source_payload: CompoundPhasePayload,
        target_database: DatabaseIR,
        target_payload: SolutionParameterPayload,
    ) -> bool:
        """Replace a solution endmember expression from a compound phase."""
        if target_payload.kind != "Endmember":
            return False
        expression = self._compound_phase_expression_for_drop(
            source_database,
            source_payload,
        )
        if not expression:
            return False
        parameter = target_payload.parameter
        parameter.expression = expression
        parameter.metadata.pop("source_expression", None)
        parameter.source.command = self._solution_parameter_command_with_expression(
            parameter,
            expression,
        )
        self._copy_referenced_functions_for_expressions(
            source_database,
            target_database,
            [expression],
        )
        self._rebuild_database_tree(target_database, show_summary=False)
        return True

    def _replace_function_expression_from_compound_phase(
        self,
        source_database: DatabaseIR,
        source_payload: CompoundPhasePayload,
        target_database: DatabaseIR,
        target_function: FunctionDefinition,
    ) -> bool:
        """Replace a Function expression from a compound phase."""
        expression = self._compound_phase_expression_for_drop(
            source_database,
            source_payload,
        )
        if not expression:
            return False
        self._replace_function_expression(
            target_function,
            expression,
        )
        self._copy_referenced_functions_for_expressions(
            source_database,
            target_database,
            [expression],
            exclude_names={target_function.name},
        )
        self._rebuild_database_tree(target_database, select_payload=target_function)
        return True

    def _copy_compound_phase_to_function_collection(
        self,
        source_database: DatabaseIR,
        source_payload: CompoundPhasePayload,
        target_database: DatabaseIR,
    ) -> bool:
        """Create or replace a Function from a compound phase expression."""
        expression = self._compound_phase_expression_for_drop(
            source_database,
            source_payload,
        )
        if not expression:
            return False
        function_name = self._derived_function_name(
            source_payload.compound_name,
            source_payload.phase_label,
        )
        return self._copy_expression_to_function_collection(
            source_database,
            target_database,
            function_name,
            expression,
        )

    def _replace_function_expression_from_solution_parameter(
        self,
        source_database: DatabaseIR,
        source_payload: SolutionParameterPayload,
        target_database: DatabaseIR,
        target_function: FunctionDefinition,
    ) -> bool:
        """Replace a Function expression from a solution endmember."""
        expression = self._solution_parameter_expression_for_drop(source_payload)
        if not expression:
            return False
        self._replace_function_expression(target_function, expression)
        self._copy_referenced_functions_for_expressions(
            source_database,
            target_database,
            [expression],
            exclude_names={target_function.name},
        )
        self._rebuild_database_tree(target_database, select_payload=target_function)
        return True

    def _copy_solution_parameter_to_function_collection(
        self,
        source_database: DatabaseIR,
        source_payload: SolutionParameterPayload,
        target_database: DatabaseIR,
    ) -> bool:
        """Create or replace a Function from a solution endmember expression."""
        expression = self._solution_parameter_expression_for_drop(source_payload)
        if not expression:
            return False
        function_name = self._derived_function_name(
            source_payload.phase.name,
            *source_payload.parameter.target,
        )
        return self._copy_expression_to_function_collection(
            source_database,
            target_database,
            function_name,
            expression,
        )

    def _compound_phase_expression_for_drop(
        self,
        database: DatabaseIR,
        payload: CompoundPhasePayload,
    ) -> str:
        """Return a compound phase Gibbs expression suitable for expression drops."""
        parameters = _compound_phase_g_parameters(database, payload.phase)
        if not parameters:
            return ""
        return _parameter_as_function_expression(parameters[0]).strip()

    def _solution_parameter_expression_for_drop(
        self,
        payload: SolutionParameterPayload,
    ) -> str:
        """Return a solution endmember expression suitable for Function drops."""
        if payload.kind != "Endmember":
            return ""
        return _parameter_as_function_expression(payload.parameter).strip()

    def _replace_function_expression(
        self,
        function: FunctionDefinition,
        expression: str,
    ) -> None:
        """Replace a Function with a full-range expression."""
        function.expression = expression
        function.gibbs_ranges.clear()
        function.temperature_ranges = []
        function.source.command = f"FUNCTION {function.name} {expression} !"
        function.dirty = True
        _normalize_function_temperature_ranges(function)

    def _copy_expression_to_function_collection(
        self,
        source_database: DatabaseIR,
        target_database: DatabaseIR,
        function_name: str,
        expression: str,
    ) -> bool:
        """Create or replace a Function from a copied expression."""
        existing = _function_by_name(target_database, function_name)
        if existing is not None and not self._confirm_database_tree_replace(
            "function",
            function_name,
            target_database,
        ):
            return False
        if existing is None:
            function = FunctionDefinition(function_name, expression)
            target_database.functions.append(function)
        else:
            function = existing
            self._replace_function_expression(function, expression)
        if existing is None:
            self._replace_function_expression(function, expression)
        self._copy_referenced_functions_for_expressions(
            source_database,
            target_database,
            [expression],
            exclude_names={function.name},
        )
        self._rebuild_database_tree(target_database, select_payload=function)
        return True

    def _derived_function_name(self, *parts: str) -> str:
        """Return a stable Function name for expression-copy drops."""
        text = "_".join(str(part) for part in parts if str(part).strip())
        name = re.sub(r"[^A-Za-z0-9_]+", "_", text).strip("_").upper()
        if not name:
            return "COPIED_FUNCTION"
        if name[0].isdigit():
            return f"F_{name}"
        return name

    def _is_function_collection_drop_target(self, target_index: QModelIndex) -> bool:
        """Return whether a target should receive a copied Function object."""
        if not target_index.isValid():
            return False
        payload = self.tree_model.payload(target_index)
        if isinstance(payload, DatabaseIR):
            return True
        return self.tree_model.label(target_index) == "Functions"

    def _solution_parameter_command_with_expression(
        self,
        parameter: Parameter,
        expression: str,
    ) -> str:
        """Return a PAR command for a solution parameter expression."""
        target = (
            ":".join(parameter.target) if parameter.target else parameter.phase_name
        )
        parameter_type = str(parameter.parameter_type or "G").strip() or "G"
        return (
            f"PAR {parameter_type}({parameter.phase_name},{target};"
            f"{parameter.order}) {expression} !"
        )

    def _copy_compound_phase_to_record(
        self,
        source_database: DatabaseIR,
        payload: CompoundPhasePayload,
        target_database: DatabaseIR,
        target_record: CompoundRecord,
    ) -> bool:
        """Copy one compound phase state into a target compound record."""
        existing_phase = self._compound_phase_by_record_label(
            target_record,
            payload.phase_label,
        )
        replacing = existing_phase is not None
        if replacing and not self._confirm_database_tree_replace(
            "phase",
            f"{target_record.name}:{payload.phase_label}",
            target_database,
        ):
            return False
        if not replacing and len(target_record.phases) >= 10:
            return False

        copied_phase = copy.deepcopy(payload.phase)
        copied_phase.name = target_record.name
        copied_phase.metadata["phase_label"] = payload.phase_label
        copied_phase.metadata["phase_id"] = _new_compound_phase_state_id(
            target_database,
            target_record.name,
        )
        source_parameters = _compound_phase_g_parameters(
            source_database,
            payload.phase,
        )
        copied_parameters: list[Parameter] = []
        for parameter in source_parameters:
            copied = copy.deepcopy(parameter)
            copied.phase_name = target_record.name
            copied.target = _compound_parameter_target_species(copied_phase)
            copied.metadata["phase_id"] = _compound_phase_state_id(copied_phase)
            target_text = ":".join(copied.target)
            expression = _parameter_as_function_expression(parameter)
            copied.expression = expression
            copied.source.command = (
                f"PARAMETER G({target_record.name},{target_text};0) "
                f"{expression} !"
            )
            copied_parameters.append(copied)

        replace_index = None
        if replacing and existing_phase is not None:
            try:
                replace_index = target_database.phases.index(existing_phase)
            except ValueError:
                replace_index = None
            self._remove_compound_phase(target_database, existing_phase)
        self._copy_referenced_functions_for_parameters(
            source_database,
            target_database,
            copied_parameters,
        )
        self._copy_species_for_phases(source_database, target_database, [payload.phase])
        if replace_index is None:
            target_database.phases.append(copied_phase)
        else:
            target_database.phases.insert(
                min(replace_index, len(target_database.phases)),
                copied_phase,
            )
        target_database.parameters.extend(copied_parameters)
        if not copied_parameters:
            target_database.parameters.append(_new_compound_g_parameter(copied_phase))
        self._rebuild_database_tree(target_database, show_summary=False)
        self._select_compound_phase(copied_phase)
        return True

    def _copy_thermo_range_to_target(
        self,
        source_database: DatabaseIR,
        payload: ThermoRangePayload,
        target_database: DatabaseIR,
        target_index: QModelIndex,
    ) -> bool:
        """Copy one Cp/G range into a target function or compound phase."""
        source_row = self._gibbs_row_for_range_payload(source_database, payload)
        if source_row is None:
            return False
        target_payload = self.tree_model.payload(target_index)
        target_function = self._target_function_for_drop(
            target_database,
            target_index,
        )
        if target_function is not None:
            gibbs_rows = _direct_function_gibbs_rows(target_function)
            replace_index = (
                target_payload.index - 1
                if isinstance(target_payload, ThermoRangePayload)
                and target_payload.owner_kind == "Function"
                and target_payload.property_name == "Cp"
                else None
            )
            if replace_index is not None:
                if replace_index < 0 or replace_index >= len(gibbs_rows):
                    return False
                if not self._confirm_database_tree_replace(
                    "Cp range",
                    f"{target_payload.owner_name}:{target_payload.index}",
                    target_database,
                ):
                    return False
                gibbs_rows[replace_index] = source_row
                selected_index = replace_index
            else:
                if len(gibbs_rows) >= 6:
                    return False
                gibbs_rows.append(source_row)
                selected_index = len(gibbs_rows) - 1
            normalized_rows, selected_index = _normalized_gibbs_rows(
                gibbs_rows,
                selected_index=selected_index,
            )
            _set_function_gibbs_rows(target_function, normalized_rows)
            _normalize_function_temperature_ranges(target_function)
            self._select_function_cp_range(
                target_database,
                target_function,
                selected_index,
            )
            return True

        if (
            isinstance(target_payload, ThermoRangePayload)
            and target_payload.owner_kind == "Compound"
            and target_payload.property_name == "Cp"
        ):
            replace_index = target_payload.index - 1
        else:
            replace_index = None
        target_phase = self._target_compound_phase_for_drop(
            target_database,
            target_index,
        )
        if target_phase is None:
            return False
        parameters = _compound_phase_g_parameters(target_database, target_phase)
        target_parameter = parameters[0] if parameters else None
        if target_parameter is None:
            target_parameter = _new_compound_g_parameter(target_phase)
            target_database.parameters.append(target_parameter)
        function = FunctionDefinition(
            name=f"{target_phase.name}_copy_target",
            expression=_parameter_as_function_expression(target_parameter),
            source=target_parameter.source,
        )
        gibbs_rows = _direct_function_gibbs_rows(function)
        if not gibbs_rows:
            _h298, _s298, _cp_rows, gibbs_rows = _function_thermo_data(
                function,
                target_database,
            )
        if replace_index is not None:
            if replace_index < 0 or replace_index >= len(gibbs_rows):
                return False
            if not self._confirm_database_tree_replace(
                "Cp range",
                f"{target_payload.owner_name}:{target_payload.index}",
                target_database,
            ):
                return False
            gibbs_rows[replace_index] = source_row
            selected_index = replace_index
        else:
            if len(gibbs_rows) >= 6:
                return False
            gibbs_rows.append(source_row)
            selected_index = len(gibbs_rows) - 1
        normalized_rows, selected_index = _normalized_gibbs_rows(
            gibbs_rows,
            selected_index=selected_index,
        )
        _set_parameter_gibbs_rows(target_parameter, normalized_rows)
        self._rebuild_database_tree(target_database, show_summary=False)
        self._select_compound_cp_range(
            self._compound_phase_owner_name(target_database, target_phase),
            selected_index + 1,
        )
        return True

    def _gibbs_row_for_range_payload(
        self,
        database: DatabaseIR,
        payload: ThermoRangePayload,
    ) -> list[float] | None:
        """Return a copied Gibbs row for one function or compound Cp payload."""
        if payload.owner_kind == "Function":
            function = _function_by_name(database, payload.owner_name)
            if function is None:
                return None
            gibbs_rows = _direct_function_gibbs_rows(function)
            if not gibbs_rows:
                _h298, _s298, _cp_rows, gibbs_rows = _function_thermo_data(
                    function,
                    database,
                )
        else:
            parameter = _compound_parameter_for_range_payload(database, payload)
            if parameter is None:
                return None
            function = FunctionDefinition(
                name=f"{payload.owner_name}_copy_source",
                expression=_parameter_as_function_expression(parameter),
                source=parameter.source,
            )
            gibbs_rows = _direct_function_gibbs_rows(function)
            if not gibbs_rows:
                _h298, _s298, _cp_rows, gibbs_rows = _function_thermo_data(
                    function,
                    database,
                )
        index = payload.index - 1
        if index < 0 or index >= len(gibbs_rows):
            return None
        return [float(value) for value in gibbs_rows[index]]

    def _target_function_for_drop(
        self,
        database: DatabaseIR,
        target_index: QModelIndex,
    ) -> FunctionDefinition | None:
        """Return the target function for a function Cp drop."""
        payload = self.tree_model.payload(target_index)
        if isinstance(payload, FunctionDefinition):
            return payload
        if isinstance(payload, ThermoRangePayload) and payload.owner_kind == "Function":
            return _function_by_name(database, payload.owner_name)
        return self._ancestor_payload(target_index, FunctionDefinition)

    def _target_compound_phase_for_drop(
        self,
        database: DatabaseIR,
        target_index: QModelIndex,
    ) -> Phase | None:
        """Return the target compound phase for a compound Cp drop."""
        payload = self.tree_model.payload(target_index)
        if isinstance(payload, CompoundPhasePayload):
            return payload.phase
        if isinstance(payload, ThermoRangePayload) and payload.owner_kind == "Compound":
            compound_name, _sep, phase_label = payload.owner_name.partition(":")
            return _compound_phase_by_label(database, compound_name, phase_label)
        phase_payload = self._ancestor_payload(target_index, CompoundPhasePayload)
        return phase_payload.phase if phase_payload is not None else None

    def _target_compound_record_for_drop(
        self,
        target_index: QModelIndex,
    ) -> CompoundRecord | None:
        """Return the target compound record for a compound phase drop."""
        payload = self.tree_model.payload(target_index)
        if isinstance(payload, CompoundRecord):
            return payload
        return self._ancestor_payload(target_index, CompoundRecord)

    def _ancestor_payload(self, index: QModelIndex, payload_type: type) -> Any:
        """Return the nearest ancestor payload matching a type."""
        current = index
        while current.isValid():
            payload = self.tree_model.payload(current)
            if isinstance(payload, payload_type):
                return payload
            current = current.parent()
        return None

    def _compound_record_by_name(
        self,
        database: DatabaseIR,
        name: str,
    ) -> CompoundRecord | None:
        """Return a compound record from current database phases by name."""
        phases = [
            phase
            for phase in database.phases
            if not _is_solution_phase(phase) and phase.name.upper() == name.upper()
        ]
        return CompoundRecord(name, tuple(phases)) if phases else None

    def _compound_records_share_phase_objects(
        self,
        left: CompoundRecord,
        right: CompoundRecord,
    ) -> bool:
        """Return whether two compound records represent the same live phases."""
        return {id(phase) for phase in left.phases} == {
            id(phase) for phase in right.phases
        }

    def _compound_record_first_phase_index(
        self,
        database: DatabaseIR,
        record: CompoundRecord | None,
    ) -> int | None:
        """Return the first phase-list index occupied by a compound record."""
        if record is None:
            return None
        phase_ids = {id(phase) for phase in record.phases}
        for index, phase in enumerate(database.phases):
            if id(phase) in phase_ids:
                return index
        return None

    def _compound_phase_by_record_label(
        self,
        record: CompoundRecord,
        phase_label: str,
    ) -> Phase | None:
        """Return a phase in a record by visible phase label."""
        for index, phase in enumerate(record.phases, start=1):
            if _compound_phase_display_label(phase, index) == phase_label:
                return phase
        return None

    def _remove_compound_record(
        self,
        database: DatabaseIR,
        record: CompoundRecord,
    ) -> None:
        """Remove all phases and parameters for a compound record."""
        phase_ids = {
            _compound_phase_state_id(phase)
            for phase in record.phases
            if _compound_phase_state_id(phase)
        }
        database.phases = [
            phase
            for phase in database.phases
            if not (
                not _is_solution_phase(phase)
                and phase.name.upper() == record.name.upper()
            )
        ]
        database.parameters = [
            parameter
            for parameter in database.parameters
            if not (
                parameter.phase_name.upper() == record.name.upper()
                or (
                    parameter.metadata.get("phase_id")
                    and parameter.metadata.get("phase_id") in phase_ids
                )
            )
        ]

    def _remove_compound_phase(self, database: DatabaseIR, phase: Phase) -> None:
        """Remove one compound phase and its phase-scoped parameters."""
        phase_id = _compound_phase_state_id(phase)
        database.phases = [
            candidate for candidate in database.phases if candidate is not phase
        ]
        if phase_id:
            database.parameters = [
                parameter
                for parameter in database.parameters
                if parameter.metadata.get("phase_id") != phase_id
            ]
            return
        database.parameters = [
            parameter
            for parameter in database.parameters
            if parameter.phase_name.upper() != phase.name.upper()
        ]

    def _copy_species_for_phases(
        self,
        source_database: DatabaseIR,
        target_database: DatabaseIR,
        phases: list[Phase],
    ) -> None:
        """Copy species referenced by copied phases when missing in target."""
        names = {
            str(species_name).strip()
            for phase in phases
            for constituent_set in phase.constituents
            for species_name in constituent_set.species
            if str(species_name).strip()
        }
        existing = {species.name.upper() for species in target_database.species}
        for species in source_database.species:
            if species.name in names and species.name.upper() not in existing:
                target_database.species.append(copy.deepcopy(species))
                existing.add(species.name.upper())

    def _copy_referenced_functions_for_parameters(
        self,
        source_database: DatabaseIR,
        target_database: DatabaseIR,
        parameters: list[Parameter],
    ) -> None:
        """Copy functions referenced by copied parameter expressions."""
        expressions = [
            _parameter_as_function_expression(parameter)
            for parameter in parameters
        ]
        self._copy_referenced_functions_for_expressions(
            source_database,
            target_database,
            expressions,
        )

    def _copy_referenced_functions_for_expressions(
        self,
        source_database: DatabaseIR,
        target_database: DatabaseIR,
        expressions: list[str],
        *,
        exclude_names: set[str] | None = None,
    ) -> None:
        """Copy recursively referenced functions without replacing target ones."""
        functions_by_name = {
            function.name.upper(): function
            for function in source_database.functions
        }
        if not functions_by_name:
            return
        excluded = {name.upper() for name in (exclude_names or set())}
        visited: set[str] = set()
        target_function_names = {
            function.name.upper()
            for function in target_database.functions
        }

        def visit(expression: str) -> None:
            for function_name in _referenced_function_names(
                expression,
                functions_by_name,
            ):
                if function_name in visited or function_name in excluded:
                    continue
                visited.add(function_name)
                function = functions_by_name[function_name]
                if function_name not in target_function_names:
                    target_database.functions.append(copy.deepcopy(function))
                    target_function_names.add(function_name)
                for dependency_expression in self._function_dependency_expressions(
                    function,
                ):
                    visit(dependency_expression)

        for expression in expressions:
            visit(expression)

    def _function_dependency_expressions(
        self,
        function: FunctionDefinition,
    ) -> list[str]:
        """Return expressions that may contain function dependencies."""
        expressions = [str(function.expression or "")]
        source_expression = str(function.tdb_src.source_expression or "")
        if source_expression and source_expression not in expressions:
            expressions.append(source_expression)
        return expressions

    def _compound_phase_owner_name(
        self,
        database: DatabaseIR,
        phase: Phase,
    ) -> str:
        """Return tree owner name for one compound phase."""
        phases = [
            candidate
            for candidate in database.phases
            if candidate.name == phase.name and not _is_solution_phase(candidate)
        ]
        for index, candidate in enumerate(phases, start=1):
            if candidate is phase:
                return f"{phase.name}:{_compound_phase_display_label(phase, index)}"
        return f"{phase.name}:{phase.metadata.get('phase_label', 'S1')}"

    def _confirm_database_tree_replace(
        self,
        kind: str,
        name: str,
        target_database: DatabaseIR,
    ) -> bool:
        """Ask whether a drag/drop copy should replace an existing object."""
        result = QMessageBox.question(
            self,
            f"Replace {kind.title()}?",
            (
                f"{kind.title()} '{name}' already exists in "
                f"{target_database.name}. Replace it?"
            ),
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        return result == QMessageBox.StandardButton.Yes
