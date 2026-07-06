"""Element-subset operations on the format-neutral database IR.

`split_database` projects a multicomponent database onto an element subset.
The projection is thermodynamically exact for CEF-family models: any
parameter whose constituent array references a removed constituent multiplies
a site fraction that is identically zero inside the subsystem, so dropping it
cannot change any property evaluated there.  The differential test in
tests/test_devop/test19_tdb_split.py enforces this against full-database
equilibrium results.

Merge (`merge_database`) is designed but not implemented yet; see
_ongoing_development/tdb_SplitMerge.md for the design contract.
"""

from __future__ import annotations

import copy
import re
from collections.abc import Sequence

from .model import DatabaseIR, Diagnostic, FunctionDefinition, Parameter, Phase

_IDENTIFIER_RE = re.compile(r"(?<![A-Za-z0-9_])([A-Za-z_][A-Za-z0-9_]*)")

# Pseudo-elements that always survive a split when present in the parent.
_ALWAYS_KEPT_ELEMENTS = ("VA", "/-")

# tdb_command.parsed keys whose value names a phase; a command is dropped
# when any of these names a phase that did not survive the split.
_COMMAND_PHASE_KEYS = ("target_phase", "ordered_phase", "disordered_phase")


def split_database(database: DatabaseIR, elements: Sequence[str]) -> DatabaseIR:
    """Project a database IR onto an element subset.

    Input
    -----
    database : Editable ``DatabaseIR``, e.g. from
               ``eq.read_tdb(path, editable=True)``.
    elements : Element symbols to keep, e.g. ``["Al", "Fe", "Si"]``.
               ``VA`` and the electron pseudo-element are always kept.

    Output
    ------
    A new ``DatabaseIR`` containing only the phases, species, parameters,
    functions, commands, and references representable in the subsystem.
    Info diagnostics record what was dropped; the input is not modified.
    """
    if not isinstance(database, DatabaseIR):
        raise TypeError(
            "split_database expects a DatabaseIR; parse the file with "
            "eq.read_tdb(path, editable=True) first."
        )

    available = {element.symbol.upper(): element for element in database.elements}
    requested = [str(symbol) for symbol in elements]
    unknown = [symbol for symbol in requested if symbol.upper() not in available]
    if unknown:
        raise ValueError(
            f"Elements not in database: {unknown}. "
            f"Available: {sorted(available)}"
        )

    kept_elements = {symbol.upper() for symbol in requested}
    for pseudo in _ALWAYS_KEPT_ELEMENTS:
        if pseudo in available:
            kept_elements.add(pseudo)

    # Species survive when their full composition lies inside the subset.
    kept_species = {
        species.name.upper()
        for species in database.species
        if species.composition
        and all(symbol.upper() in kept_elements for symbol in species.composition)
    }
    allowed_constituents = kept_elements | kept_species

    # Order/disorder pairs from active DIS_PART commands: the ordered phase
    # may survive only if its disordered partner survives.
    dispart_partner = {}
    for command in database.tdb_commands:
        if not command.active:
            continue
        if command.parsed.get("action", "").upper() != "DIS_PART":
            continue
        ordered = command.parsed.get("ordered_phase", "").upper()
        disordered = command.parsed.get("disordered_phase", "").upper()
        if ordered and disordered:
            dispart_partner[ordered] = disordered

    kept_phases: list[Phase] = []
    dropped_phase_names: set[str] = set()
    for phase in database.phases:
        projected = _project_phase(phase, allowed_constituents)
        if projected is None:
            dropped_phase_names.add(phase.name.upper())
        else:
            kept_phases.append(projected)

    # Enforce DIS_PART atomicity after constituent projection.
    kept_phase_names = {phase.name.upper() for phase in kept_phases}
    dropped_ordered = {
        ordered
        for ordered, disordered in dispart_partner.items()
        if ordered in kept_phase_names and disordered not in kept_phase_names
    }
    if dropped_ordered:
        kept_phases = [
            phase for phase in kept_phases if phase.name.upper() not in dropped_ordered
        ]
        dropped_phase_names |= dropped_ordered
        kept_phase_names -= dropped_ordered

    kept_parameters: list[Parameter] = []
    for parameter in database.parameters:
        if parameter.phase_name.upper() not in kept_phase_names:
            continue
        if _parameter_targets_allowed(parameter, allowed_constituents):
            kept_parameters.append(copy.deepcopy(parameter))

    kept_functions = _function_closure(database.functions, kept_parameters)

    kept_commands = [
        copy.deepcopy(command)
        for command in database.tdb_commands
        if _command_survives(command, kept_phase_names)
    ]

    cited_references = {
        parameter.metadata.get("reference", "").upper()
        for parameter in kept_parameters
    }
    cited_references.discard("")
    kept_references = [
        copy.deepcopy(reference)
        for reference in database.references
        if reference.key.upper() in cited_references
    ]

    result = DatabaseIR(
        schema_version=database.schema_version,
        name=f"{database.name} (split: {'-'.join(requested)})",
        metadata=dict(database.metadata),
        elements=[
            copy.deepcopy(element)
            for element in database.elements
            if element.symbol.upper() in kept_elements
        ],
        species=[
            copy.deepcopy(species)
            for species in database.species
            if species.name.upper() in kept_species
        ],
        functions=kept_functions,
        phases=kept_phases,
        parameters=kept_parameters,
        tdb_commands=kept_commands,
        references=kept_references,
    )
    result.metadata["split_from"] = database.name
    result.metadata["split_elements"] = ",".join(requested)

    n_dropped_parameters = len(database.parameters) - len(kept_parameters)
    result.diagnostics.append(
        Diagnostic(
            severity="info",
            message=(
                f"split kept {len(kept_phases)}/{len(database.phases)} phases, "
                f"dropped {n_dropped_parameters} parameters outside "
                f"{'-'.join(requested)}"
            ),
        )
    )
    for ordered in sorted(dropped_ordered):
        result.diagnostics.append(
            Diagnostic(
                severity="warning",
                message=(
                    f"dropped ordered phase {ordered} because its DIS_PART "
                    "partner is outside the element subset"
                ),
            )
        )
    return result


def _project_phase(phase: Phase, allowed_constituents: set[str]) -> Phase | None:
    """Filter phase constituents; return None when a sublattice empties."""
    projected = copy.deepcopy(phase)
    if not projected.constituents:
        return projected
    for constituent_set in projected.constituents:
        constituent_set.species = [
            name
            for name in constituent_set.species
            if name.upper() in allowed_constituents
        ]
        if not constituent_set.species:
            return None
    return projected


def _parameter_targets_allowed(
    parameter: Parameter,
    allowed_constituents: set[str],
) -> bool:
    """Check that every constituent in the parameter target survives."""
    for target in parameter.target:
        for section in str(target).split(":"):
            for token in section.split(","):
                name = token.strip()
                if not name or name == "*":
                    continue
                # Interaction wildcards and default markers travel with the
                # constituent token in some dialects.
                name = name.rstrip("%").strip()
                if name and name.upper() not in allowed_constituents:
                    return False
    return True


def _function_closure(
    functions: list[FunctionDefinition],
    parameters: list[Parameter],
) -> list[FunctionDefinition]:
    """Keep functions transitively referenced by parameters or functions."""
    by_name = {function.name.upper(): function for function in functions}

    def references(text: str) -> set[str]:
        return {
            match.group(1).upper()
            for match in _IDENTIFIER_RE.finditer(text)
            if match.group(1).upper() in by_name
        }

    pending: set[str] = set()
    for parameter in parameters:
        pending |= references(parameter.expression)

    kept: set[str] = set()
    while pending:
        name = pending.pop()
        if name in kept:
            continue
        kept.add(name)
        function = by_name[name]
        referenced = references(function.expression)
        referenced |= references(function.tdb_src.source_expression)
        for gibbs_range in function.gibbs_ranges:
            referenced |= references(gibbs_range.Gibbs)
        for term in function.tdb_src.pertinent_functions:
            if term.name.upper() in by_name:
                referenced.add(term.name.upper())
        pending |= referenced - kept

    return [
        copy.deepcopy(function)
        for function in functions
        if function.name.upper() in kept
    ]


def _command_survives(command, kept_phase_names: set[str]) -> bool:
    """Drop commands that reference a phase removed by the split."""
    for key in _COMMAND_PHASE_KEYS:
        name = command.parsed.get(key, "").strip()
        if not name or name in ("@", "*"):
            continue
        if name.upper() not in kept_phase_names:
            return False
    return True
