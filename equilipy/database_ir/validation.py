"""Headless validation helpers for DatabaseIR objects."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

from .model import DatabaseIR, Diagnostic, FunctionDefinition, Phase

_SOLUTION_MODELS = {
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


@dataclass(slots=True)
class ValidationReport:
    """Result from TDB validation and optional automatic correction."""

    diagnostics: list[Diagnostic]
    corrected_count: int = 0

    @property
    def errors(self) -> list[Diagnostic]:
        """Return error diagnostics."""
        return [
            diagnostic
            for diagnostic in self.diagnostics
            if diagnostic.severity.lower() == "error"
        ]

    @property
    def warnings(self) -> list[Diagnostic]:
        """Return warning diagnostics."""
        return [
            diagnostic
            for diagnostic in self.diagnostics
            if diagnostic.severity.lower() == "warning"
        ]

    @property
    def ok(self) -> bool:
        """Return whether validation has no errors."""
        return not self.errors


def validate_tdb(
    database: DatabaseIR,
    *,
    auto_correct: bool = False,
    progress_callback: Callable[[int, int, str], None] | None = None,
) -> ValidationReport:
    """Validate a DatabaseIR with TDB-oriented thermodynamic checks.

    When ``auto_correct`` is true, discontinuous Gibbs piecewise ranges are
    balanced in place by shifting the higher-temperature range's ``a`` term.
    Functions are corrected before compound parameters so compound expressions
    can benefit from corrected function dependencies.
    """
    corrected_count = 0
    if auto_correct:
        corrected_count += _auto_correct_function_gibbs(database)
        corrected_count += _auto_correct_compound_gibbs(database)

    diagnostics = _validate_database_thermodynamics(database, progress_callback)
    return ValidationReport(
        diagnostics=_unique_diagnostics(diagnostics),
        corrected_count=corrected_count,
    )


def _validate_database_thermodynamics(
    database: DatabaseIR,
    progress_callback: Callable[[int, int, str], None] | None = None,
) -> list[Diagnostic]:
    """Return parser, Function, compound, and solution validation diagnostics."""
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
    diagnostics.extend(_validate_function_thermodynamics(database, advance))
    diagnostics.extend(_validate_compound_thermodynamics(database, advance))
    diagnostics.extend(_validate_solution_thermodynamics(database, advance))
    return diagnostics


def _validate_function_thermodynamics(
    database: DatabaseIR,
    progress_callback: Callable[[str], None],
) -> list[Diagnostic]:
    """Run strict Function checks, including thermo conversion."""
    helpers = _validation_helpers()
    diagnostics: list[Diagnostic] = []
    for function in database.functions:
        for message in helpers["function_name_warning_messages"](function):
            diagnostics.append(
                Diagnostic(
                    "warning",
                    f"Function {function.name}: {message}",
                    function.source,
                )
            )
        for message in helpers["function_name_error_messages"](function):
            diagnostics.append(
                Diagnostic(
                    "error",
                    f"Function {function.name}: {message}",
                    function.source,
                )
            )

        _h298, _s298, _cp_rows, gibbs_rows = helpers["function_thermo_data"](
            function,
            database,
        )
        if not gibbs_rows:
            if function.expression.strip() or function.gibbs_ranges:
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
            progress_callback(f"Function {function.name}")
            continue

        diagnostics.extend(
            Diagnostic(
                "warning",
                f"Function {function.name}: {message}",
                function.source,
            )
            for message in helpers["function_gibbs_continuity_warnings"](gibbs_rows)
        )
        progress_callback(f"Function {function.name}")
    return diagnostics


def _validate_compound_thermodynamics(
    database: DatabaseIR,
    progress_callback: Callable[[str], None],
) -> list[Diagnostic]:
    """Run full compound Gibbs checks with Function dependencies resolved."""
    helpers = _validation_helpers()
    diagnostics: list[Diagnostic] = []
    phase_counts: dict[str, int] = {}
    for phase in database.phases:
        if _is_solution_phase(phase):
            continue
        compound_name = str(phase.name).strip() or "Compound"
        phase_counts[compound_name] = phase_counts.get(compound_name, 0) + 1
        phase_label = helpers["compound_phase_display_label"](
            phase,
            phase_counts[compound_name],
        )
        owner_name = f"{compound_name}:{phase_label}"
        for parameter in helpers["compound_phase_g_parameters"](database, phase):
            function = FunctionDefinition(
                name=owner_name.replace(":", "_"),
                expression=helpers["parameter_as_function_expression"](parameter),
                source=parameter.source,
            )
            _h298, _s298, _cp_rows, gibbs_rows = helpers["function_thermo_data"](
                function,
                database,
            )
            if not gibbs_rows:
                if function.expression.strip():
                    diagnostics.append(
                        Diagnostic(
                            "warning",
                            (
                                f"Compound {compound_name} phase {phase_label}: "
                                "could not be converted into H/S/Cp/G rows."
                            ),
                            parameter.source,
                        )
                    )
                continue
            diagnostics.extend(
                Diagnostic(
                    "warning",
                    f"Compound {compound_name} phase {phase_label}: {message}",
                    parameter.source,
                )
                for message in helpers["function_gibbs_continuity_warnings"](
                    gibbs_rows
                )
            )
        progress_callback(f"Compound phase {compound_name}:{phase_label}")
    return diagnostics


def _validate_solution_thermodynamics(
    database: DatabaseIR,
    progress_callback: Callable[[str], None],
) -> list[Diagnostic]:
    """Count solution phases in full validation without eager solution checks."""
    for phase in database.phases:
        if _is_solution_phase(phase):
            progress_callback(f"Solution phase {phase.name}")
    return []


def _auto_correct_function_gibbs(database: DatabaseIR) -> int:
    """Balance Function Gibbs discontinuities in place."""
    helpers = _validation_helpers()
    corrected_count = 0
    for function in database.functions:
        _h298, _s298, _cp_rows, gibbs_rows = helpers["function_thermo_data"](
            function,
            database,
        )
        balanced_rows, count = helpers[
            "balance_gibbs_rows_a_coefficients_sequentially"
        ](gibbs_rows)
        if count:
            helpers["set_function_gibbs_rows"](function, balanced_rows)
            corrected_count += count
    return corrected_count


def _auto_correct_compound_gibbs(database: DatabaseIR) -> int:
    """Balance compound Gibbs parameter discontinuities in place."""
    helpers = _validation_helpers()
    corrected_count = 0
    phase_counts: dict[str, int] = {}
    for phase in database.phases:
        if _is_solution_phase(phase):
            continue
        compound_name = str(phase.name).strip() or "Compound"
        phase_counts[compound_name] = phase_counts.get(compound_name, 0) + 1
        phase_label = helpers["compound_phase_display_label"](
            phase,
            phase_counts[compound_name],
        )
        owner_name = f"{compound_name}:{phase_label}"
        for parameter in helpers["compound_phase_g_parameters"](database, phase):
            function = FunctionDefinition(
                name=owner_name.replace(":", "_"),
                expression=helpers["parameter_as_function_expression"](parameter),
                source=parameter.source,
            )
            _h298, _s298, _cp_rows, gibbs_rows = helpers["function_thermo_data"](
                function,
                database,
            )
            balanced_rows, count = helpers[
                "balance_gibbs_rows_a_coefficients_sequentially"
            ](gibbs_rows)
            if count:
                helpers["set_parameter_gibbs_rows"](parameter, balanced_rows)
                corrected_count += count
    return corrected_count


def _database_validation_object_count(database: DatabaseIR) -> int:
    """Return the number of objects represented in validation progress."""
    compound_count = sum(
        1 for phase in database.phases if not _is_solution_phase(phase)
    )
    solution_count = sum(1 for phase in database.phases if _is_solution_phase(phase))
    return max(1, len(database.functions) + compound_count + solution_count)


def _is_solution_phase(phase: Phase) -> bool:
    """Return whether a phase should use solution-phase validation."""
    return phase.model.upper() in _SOLUTION_MODELS


def _unique_diagnostics(diagnostics: list[Diagnostic]) -> list[Diagnostic]:
    """Return diagnostics without duplicate severity/message/source triples."""
    unique: list[Diagnostic] = []
    seen: set[tuple[str, str, str, int]] = set()
    for diagnostic in diagnostics:
        key = (
            diagnostic.severity,
            diagnostic.message,
            diagnostic.source.file,
            diagnostic.source.line,
        )
        if key in seen:
            continue
        seen.add(key)
        unique.append(diagnostic)
    return unique


def _validation_helpers() -> dict[str, object]:
    """Load the existing thermo conversion helpers lazily.

    The logic still lives in the GUI modules today; importing lazily keeps the
    public DatabaseIR module importable without immediately loading Qt widgets.
    """
    import equilipy.gui.database.cp_editor as cp_editor
    import equilipy.gui.database.overview as overview

    # Some thermo helpers historically relied on main_window wildcard exports.
    # Provide the missing formatting helper explicitly for headless API usage.
    if not hasattr(cp_editor, "_compact_number"):
        cp_editor._compact_number = overview._compact_number

    return {
        "balance_gibbs_rows_a_coefficients_sequentially": (
            cp_editor._balance_gibbs_rows_a_coefficients_sequentially
        ),
        "compound_phase_display_label": overview._compound_phase_display_label,
        "compound_phase_g_parameters": overview._compound_phase_g_parameters,
        "function_gibbs_continuity_warnings": (
            cp_editor._function_gibbs_continuity_warnings
        ),
        "function_name_warning_messages": cp_editor._function_name_warning_messages,
        "function_name_error_messages": cp_editor._function_name_error_messages,
        "function_thermo_data": cp_editor._function_thermo_data,
        "parameter_as_function_expression": overview._parameter_as_function_expression,
        "set_function_gibbs_rows": cp_editor._set_function_gibbs_rows,
        "set_parameter_gibbs_rows": overview._set_parameter_gibbs_rows,
    }
