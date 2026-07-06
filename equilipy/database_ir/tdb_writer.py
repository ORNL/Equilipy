"""TDB writer for the DatabaseIR development path.

The writer is intentionally narrow for now: it emits a strict, readable TDB
header plus the DatabaseIR sections currently supported by the parser/editor.
Later database sections should be added here rather than creating a second
writer path.
"""

from __future__ import annotations

import re
from dataclasses import replace
from datetime import date
from fractions import Fraction
from functools import reduce
from itertools import product
from math import gcd
from pathlib import Path
from typing import Literal

from .model import (
    ConstituentSet,
    DatabaseIR,
    Element,
    FunctionDefinition,
    GibbsRange,
    Parameter,
    Phase,
    TdbCommand,
    TdbReference,
)
from .tdb_canonical import (
    DISORDERED_PHASE_CANONICAL_NAMES,
    canonical_disordered_phase_name,
)

_SPECIAL_ELEMENT_SYMBOLS = {"/-", "*", "VA"}
_FLOAT_STATUS_RE = re.compile(
    r"^\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EDed][+-]?\d+)?)\s+([A-Za-z])"
    r"(?:\s+(.*))?\s*$"
)
_BARE_DECIMAL_EXPONENT_RE = re.compile(
    r"(?<![A-Za-z0-9_])([+\-]?\d+)\.(?=[EeDd][+\-]?\d+)"
)
_MAX_TDB_LINE_LENGTH = 78
_MAX_TDB_REFERENCE_LINE_LENGTH = 80
_CONTINUATION_PREFIX = "    "
_FUNCTION_NAME_WIDTH = 25
_COMMENT_BODY_WIDTH = _MAX_TDB_LINE_LENGTH - 2
_METADATA_LABEL_WIDTH = 20
_ENDMEMBER_COMPOUND_SUFFIX = "_ENDMBR"
_ENDMEMBER_COMPOUND_GENERATOR = "tdb_endmember_compound_export"
OrderDisorderAliasMode = Literal["canonical", "preserve"]
TdbExportStyle = Literal["native", "factsage"]
_ORDER_DISORDER_HELPER_NAMES = {
    canonical: alias
    for alias, canonical in DISORDERED_PHASE_CANONICAL_NAMES.items()
}
_ORDER_DISORDER_ALIAS_MODES = {"canonical", "preserve"}
_TDB_EXPORT_STYLES = {"native", "factsage"}
_ORDER_DISORDER_HELPER_OFFSET = 0.1
_FACTSAGE_SITE_RATIO_TOLERANCE = 5e-3
_FACTSAGE_SIMPLE_SITE_RATIO_MAX_DENOMINATOR = 32
_FACTSAGE_EXACT_SITE_RATIO_MAX_DENOMINATOR = 400
_FACTSAGE_SIMPLE_SITE_RATIO_MAX_SCALE = 64.0
_FACTSAGE_EXACT_SITE_RATIO_MAX_SCALE = 400.0
_FACTSAGE_ENERGY_PARAMETER_TYPES = {"G", "L"}


def write_tdb(
    database: DatabaseIR,
    path: str | Path,
    *,
    include_endmember_compounds: bool = False,
    order_disorder_alias_mode: OrderDisorderAliasMode = "preserve",
    export_style: TdbExportStyle = "native",
) -> Path:
    """Write a DatabaseIR as a TDB file with ELEMENT and FUNCTION sections."""
    output_path = Path(path)
    output_path.write_text(
        dumps_tdb(
            database,
            include_endmember_compounds=include_endmember_compounds,
            order_disorder_alias_mode=order_disorder_alias_mode,
            export_style=export_style,
        ),
        encoding="utf-8",
    )
    return output_path


def dumps_tdb(
    database: DatabaseIR,
    *,
    include_endmember_compounds: bool = False,
    order_disorder_alias_mode: OrderDisorderAliasMode = "preserve",
    export_style: TdbExportStyle = "native",
) -> str:
    """Return a strict TDB text representation of the supported IR sections."""
    style = _normalize_tdb_export_style(export_style)
    if style == "factsage":
        include_endmember_compounds = True
        order_disorder_alias_mode = "canonical"
    mode = _normalize_order_disorder_alias_mode(order_disorder_alias_mode)
    phases, source_parameters, tdb_commands = _order_disorder_export_inputs(
        database.phases,
        database.parameters,
        database.tdb_commands,
        mode,
    )
    source_parameters = _without_generated_compatibility_zero_parameters(
        source_parameters
    )
    alias_phases, alias_parameters = _without_endmember_compound_exports(
        phases,
        source_parameters,
    )
    phase_alias_map = _duplicate_phase_alias_map_for_export(
        alias_phases,
        alias_parameters,
        tdb_commands,
        extra_protected_phases=_order_disorder_alias_mode_protected_phases(
            alias_phases,
            mode,
        ),
    )
    phases, source_parameters = _endmember_compound_export_inputs(
        phases,
        source_parameters,
        include_endmember_compounds=include_endmember_compounds,
        phase_alias_map=phase_alias_map,
    )
    parameters = _parameter_exports(source_parameters, phases, phase_alias_map)
    if style == "factsage":
        phases, parameters, tdb_commands, phase_alias_map = _factsage_export_inputs(
            phases,
            parameters,
            tdb_commands,
            phase_alias_map,
        )
    setup_code_map = _tdb_setup_type_code_export_map(tdb_commands)
    phase_names = {_tdb_identifier(phase.name) for phase in phases}
    functions = _functions_with_zero_for_parameter_export(
        database.functions,
        parameters,
    )
    magnetic_codes = _magnetic_type_codes(database.tdb_commands, setup_code_map)
    phases_with_magnetic_parameters = _phases_with_magnetic_parameters(parameters)
    export_database = replace(database, phases=phases, parameters=parameters)
    lines: list[str] = _database_header(database)
    lines.append("")
    lines.extend(_element_section(database.elements))
    lines.append("")
    lines.extend(_function_section(functions))
    lines.append("")
    lines.extend(
        _tdb_setup_section(
            tdb_commands,
            setup_code_map,
            database.elements,
            phase_names,
            phase_alias_map,
            factsage_style=style == "factsage",
        )
    )
    lines.append("")
    lines.extend(
        _phase_section(
            phases,
            setup_code_map,
            phase_alias_map,
            magnetic_codes,
            phases_with_magnetic_parameters,
        )
    )
    lines.append("")
    lines.extend(_parameter_section(export_database, parameters))
    reference_lines = _reference_section(database.references)
    if reference_lines:
        lines.append("")
        lines.extend(reference_lines)
    return "\n".join(lines).rstrip() + "\n"


def _normalize_order_disorder_alias_mode(mode: str) -> str:
    """Return the writer dialect for known order/disorder helper aliases."""
    normalized = str(mode or "preserve").strip().lower().replace("_", "-")
    aliases = {
        "canonical": "canonical",
        "preserve": "preserve",
        "source": "preserve",
    }
    mapped = aliases.get(normalized)
    if mapped is None:
        valid = ", ".join(sorted(_ORDER_DISORDER_ALIAS_MODES))
        raise ValueError(
            f"Unknown order/disorder alias export mode {mode!r}; expected {valid}."
        )
    return mapped


def _normalize_tdb_export_style(style: str) -> str:
    """Return the named TDB export style."""
    normalized = (
        str(style or "native")
        .strip()
        .lower()
        .replace("-", "_")
        .replace(" ", "_")
    )
    aliases = {
        "default": "native",
        "normal": "native",
        "native": "native",
        "equilipy": "native",
        "fact_sage": "factsage",
        "fact_sage_style": "factsage",
        "factsage": "factsage",
        "factsage_style": "factsage",
    }
    mapped = aliases.get(normalized)
    if mapped is None:
        valid = ", ".join(sorted(_TDB_EXPORT_STYLES))
        raise ValueError(f"Unknown TDB export style {style!r}; expected {valid}.")
    return mapped


def _order_disorder_export_inputs(
    phases: list[Phase],
    parameters: list[Parameter],
    commands: list[TdbCommand],
    mode: str,
) -> tuple[list[Phase], list[Parameter], list[TdbCommand]]:
    """Return phase/parameter/setup inputs for one order/disorder dialect."""
    if mode == "canonical":
        return _canonical_order_disorder_alias_inputs(phases, parameters, commands)
    return _preserve_order_disorder_alias_inputs(phases, parameters, commands)


def _preserve_order_disorder_alias_inputs(
    phases: list[Phase],
    parameters: list[Parameter],
    commands: list[TdbCommand],
) -> tuple[list[Phase], list[Parameter], list[TdbCommand]]:
    """Preserve or synthesize known DIS_PART helper aliases on export."""
    export_phases = list(phases)
    source_parameters = list(parameters)
    export_parameters = list(source_parameters)
    phases_by_name = {_tdb_identifier(phase.name): phase for phase in export_phases}
    disordered_export_names: dict[str, str] = {}

    for command in commands:
        ordered_phase, disordered_phase = _dispart_phase_names(command)
        if not ordered_phase or not disordered_phase:
            continue
        canonical = canonical_disordered_phase_name(disordered_phase)
        helper = _ORDER_DISORDER_HELPER_NAMES.get(canonical)
        if not helper:
            disordered_export_names[ordered_phase] = disordered_phase
            continue
        source_phase = phases_by_name.get(canonical) or phases_by_name.get(
            disordered_phase
        )
        if helper not in phases_by_name:
            if source_phase is None:
                disordered_export_names[ordered_phase] = disordered_phase
                continue
            helper_phase = _clone_phase_for_order_disorder_alias(source_phase, helper)
            _insert_order_disorder_helper_phase(
                export_phases,
                helper_phase,
                ordered_phase,
            )
            phases_by_name[helper] = helper_phase
            export_parameters.extend(
                _with_order_disorder_helper_endmember_offsets(
                    _clone_parameters_for_order_disorder_alias(
                        source_parameters,
                        source_phase.name,
                        helper,
                    ),
                    helper_phase,
                    canonical,
                    helper,
                )
            )
        else:
            if (
                source_phase is not None
                and _tdb_identifier(source_phase.name) != helper
            ):
                export_parameters = _repair_degenerate_order_disorder_helper_offsets(
                    export_parameters,
                    source_phase.name,
                    phases_by_name[helper],
                    canonical,
                    helper,
                )
            _move_order_disorder_helper_before_ordered(
                export_phases,
                helper,
                ordered_phase,
            )
        disordered_export_names[ordered_phase] = helper

    export_commands = [
        _with_dispart_disordered_phase(
            command,
            disordered_export_names.get(
                _tdb_identifier(command.parsed.get("ordered_phase", "")),
                "",
            ),
        )
        for command in commands
    ]
    return export_phases, export_parameters, export_commands


def _insert_order_disorder_helper_phase(
    phases: list[Phase],
    helper_phase: Phase,
    ordered_phase_name: str,
) -> None:
    """Insert a synthesized DIS_PART helper before its ordered phase."""
    ordered_index = _phase_index(phases, ordered_phase_name)
    if ordered_index is None:
        phases.append(helper_phase)
    else:
        phases.insert(ordered_index, helper_phase)


def _move_order_disorder_helper_before_ordered(
    phases: list[Phase],
    helper_phase_name: str,
    ordered_phase_name: str,
) -> None:
    """Keep existing DIS_PART helper declarations before ordered phases."""
    helper_index = _phase_index(phases, helper_phase_name)
    ordered_index = _phase_index(phases, ordered_phase_name)
    if helper_index is None or ordered_index is None or helper_index < ordered_index:
        return
    helper_phase = phases.pop(helper_index)
    ordered_index = _phase_index(phases, ordered_phase_name)
    if ordered_index is None:
        phases.append(helper_phase)
    else:
        phases.insert(ordered_index, helper_phase)


def _phase_index(phases: list[Phase], phase_name: str) -> int | None:
    """Return the index of a phase by TDB-normalized name."""
    normalized = _tdb_identifier(phase_name)
    for index, phase in enumerate(phases):
        if _tdb_identifier(phase.name) == normalized:
            return index
    return None


def _canonical_order_disorder_alias_inputs(
    phases: list[Phase],
    parameters: list[Parameter],
    commands: list[TdbCommand],
) -> tuple[list[Phase], list[Parameter], list[TdbCommand]]:
    """Collapse known helper aliases to canonical disordered phase names."""
    phase_names = {_tdb_identifier(phase.name) for phase in phases}
    aliases = {
        alias: canonical
        for alias, canonical in DISORDERED_PHASE_CANONICAL_NAMES.items()
        if alias in phase_names
    }
    if not aliases:
        return (
            list(phases),
            list(parameters),
            [
                _with_dispart_disordered_phase(
                    command,
                    canonical_disordered_phase_name(
                        command.parsed.get("disordered_phase", "")
                    ),
                )
                for command in commands
            ],
        )

    export_phases: list[Phase] = []
    emitted_phases: set[str] = set()
    for phase in phases:
        phase_name = _tdb_identifier(phase.name)
        canonical = aliases.get(phase_name)
        if canonical and canonical in phase_names:
            continue
        export_phase = (
            _clone_phase_for_order_disorder_alias(phase, canonical)
            if canonical
            else phase
        )
        export_name = _tdb_identifier(export_phase.name)
        if export_name in emitted_phases:
            continue
        emitted_phases.add(export_name)
        export_phases.append(export_phase)

    canonical_keys = {
        _parameter_duplicate_key(parameter)
        for parameter in parameters
        if _tdb_identifier(parameter.phase_name) not in aliases
    }
    export_parameters: list[Parameter] = []
    emitted_parameters: set[tuple[str, str, tuple[str, ...], int]] = set()
    for parameter in parameters:
        phase_name = _tdb_identifier(parameter.phase_name)
        canonical = aliases.get(phase_name)
        export_parameter = (
            replace(
                parameter,
                phase_name=canonical,
                metadata={
                    **parameter.metadata,
                    "order_disorder_alias_collapsed_from": phase_name,
                },
            )
            if canonical
            else parameter
        )
        key = _parameter_duplicate_key(export_parameter)
        if canonical and key in canonical_keys:
            continue
        if key in emitted_parameters:
            continue
        emitted_parameters.add(key)
        export_parameters.append(export_parameter)

    export_commands = [
        _with_dispart_disordered_phase(
            command,
            canonical_disordered_phase_name(command.parsed.get("disordered_phase", "")),
        )
        for command in commands
    ]
    return export_phases, export_parameters, export_commands


def _dispart_phase_names(command: TdbCommand) -> tuple[str, str]:
    """Return ordered/disordered phase names for a DIS_PART TYPE_DEF command."""
    if _tdb_setup_keyword(command.command) != "TYPE_DEF":
        return "", ""
    if command.parsed.get("action", "").upper() != "DIS_PART":
        return "", ""
    return (
        _tdb_identifier(command.parsed.get("ordered_phase", "")),
        _tdb_identifier(command.parsed.get("disordered_phase", "")),
    )


def _clone_phase_for_order_disorder_alias(phase: Phase, name: str) -> Phase:
    """Return a phase copy with a different order/disorder helper name."""
    metadata = dict(phase.metadata)
    metadata["order_disorder_alias_export"] = "true"
    metadata["order_disorder_alias_source"] = _tdb_identifier(phase.name)
    return replace(phase, name=_tdb_identifier(name), metadata=metadata)


def _clone_parameters_for_order_disorder_alias(
    parameters: list[Parameter],
    source_phase_name: str,
    alias_phase_name: str,
) -> list[Parameter]:
    """Return parameter copies for a generated disordered helper phase."""
    source = _tdb_identifier(source_phase_name)
    alias = _tdb_identifier(alias_phase_name)
    cloned: list[Parameter] = []
    for parameter in parameters:
        if _tdb_identifier(parameter.phase_name) != source:
            continue
        metadata = dict(parameter.metadata)
        metadata["order_disorder_alias_export"] = "true"
        metadata["order_disorder_alias_source"] = source
        cloned.append(replace(parameter, phase_name=alias, metadata=metadata))
    return cloned


def _repair_degenerate_order_disorder_helper_offsets(
    parameters: list[Parameter],
    source_phase_name: str,
    helper_phase: Phase,
    canonical_phase_name: str,
    helper_phase_name: str,
) -> list[Parameter]:
    """Add the BCC helper offset when an existing helper duplicates BCC_A2."""
    if not _uses_bcc_order_disorder_helper_offset(
        canonical_phase_name, helper_phase_name
    ):
        return parameters
    source = _tdb_identifier(source_phase_name)
    helper = _tdb_identifier(helper_phase_name)
    source_parameters = {
        _phase_independent_parameter_key(parameter): parameter
        for parameter in parameters
        if _tdb_identifier(parameter.phase_name) == source
    }
    repaired: list[Parameter] = []
    for parameter in parameters:
        if (
            _tdb_identifier(parameter.phase_name) != helper
            or not _is_order_disorder_helper_offset_candidate(parameter, helper_phase)
        ):
            repaired.append(parameter)
            continue
        source_parameter = source_parameters.get(
            _phase_independent_parameter_key(parameter)
        )
        if source_parameter is None:
            repaired.append(parameter)
            continue
        helper_expression = (
            _parameter_source_expression_for_export(parameter) or parameter.expression
        )
        source_expression = (
            _parameter_source_expression_for_export(source_parameter)
            or source_parameter.expression
        )
        if _tdb_expressions_equivalent(helper_expression, source_expression):
            repaired.append(_with_order_disorder_helper_offset(parameter))
        else:
            repaired.append(parameter)
    return repaired


def _with_order_disorder_helper_endmember_offsets(
    parameters: list[Parameter],
    helper_phase: Phase,
    canonical_phase_name: str,
    helper_phase_name: str,
) -> list[Parameter]:
    """Add the BCC helper offset to synthesized helper endmember records."""
    if not _uses_bcc_order_disorder_helper_offset(
        canonical_phase_name, helper_phase_name
    ):
        return parameters
    return [
        _with_order_disorder_helper_offset(parameter)
        if _is_order_disorder_helper_offset_candidate(parameter, helper_phase)
        else parameter
        for parameter in parameters
    ]


def _uses_bcc_order_disorder_helper_offset(
    canonical_phase_name: str,
    helper_phase_name: str,
) -> bool:
    """Return whether this helper pair uses the BCC/B2 anti-degeneracy offset."""
    return (
        _tdb_identifier(canonical_phase_name) == "BCC_A2"
        and _tdb_identifier(helper_phase_name) == "A2_BCC"
    )


def _is_order_disorder_helper_offset_candidate(
    parameter: Parameter,
    helper_phase: Phase,
) -> bool:
    """Return whether a helper G endmember should receive the BCC offset."""
    if (
        _tdb_parameter_type(parameter) != "G"
        or int(parameter.order) != 0
        or len(parameter.target) != 1
        or _is_interaction_parameter(parameter)
        or _parameter_target_has_wildcard(parameter)
    ):
        return False
    target = _tdb_parameter_target_text(parameter.target[0])
    sections = target.split(":")
    if len(sections) != len(helper_phase.constituents):
        return False
    tokens = [_tdb_identifier(section.strip()) for section in sections]
    if not tokens or tokens[0] in {"", "VA", "*"}:
        return False
    return all(token == "VA" for token in tokens[1:])


def _phase_independent_parameter_key(
    parameter: Parameter,
) -> tuple[str, tuple[str, ...], int]:
    """Return a parameter identity without the phase name."""
    return (
        _tdb_parameter_type(parameter),
        tuple(
            _tdb_identifier(_tdb_parameter_target_text(target))
            for target in parameter.target
        ),
        int(parameter.order),
    )


def _with_order_disorder_helper_offset(parameter: Parameter) -> Parameter:
    """Return a parameter with the BCC helper anti-degeneracy offset added."""
    expression = _order_disorder_helper_offset_expression(
        _parameter_source_expression_for_export(parameter) or parameter.expression
    )
    metadata = dict(parameter.metadata)
    metadata.pop("source_expression", None)
    metadata["order_disorder_helper_offset"] = _tdb_number(
        _ORDER_DISORDER_HELPER_OFFSET
    )
    return replace(parameter, expression=expression, metadata=metadata)


def _order_disorder_helper_offset_expression(expression: str) -> str:
    """Add the BCC helper offset to a simple or piecewise TDB expression."""
    expression = _normalize_tdb_expression(expression or "0")
    ranges = _parse_function_ranges(expression)
    if not ranges:
        return _additive_gibbs_expression(
            _tdb_number(_ORDER_DISORDER_HELPER_OFFSET),
            expression,
        )
    offset_ranges = [
        GibbsRange(
            float(gibbs_range.T_min),
            float(gibbs_range.T_max),
            _additive_gibbs_expression(
                _tdb_number(_ORDER_DISORDER_HELPER_OFFSET),
                gibbs_range.Gibbs,
            ),
            gibbs_range.status,
        )
        for gibbs_range in ranges
    ]
    return _canonical_tdb_range_expression(offset_ranges, zero_as_function=True)


def _with_dispart_disordered_phase(
    command: TdbCommand,
    disordered_phase: str,
) -> TdbCommand:
    """Return a TYPE_DEF command with an updated DIS_PART target phase."""
    if not disordered_phase or _tdb_setup_keyword(command.command) != "TYPE_DEF":
        return command
    if command.parsed.get("action", "").upper() != "DIS_PART":
        return command
    args = list(command.args)
    upper_args = [_tdb_identifier(arg) for arg in args]
    for index, token in enumerate(upper_args):
        if token in {"DIS_PART", "DISORDER_PART"} and index + 1 < len(args):
            args[index] = "DIS_PART"
            args[index + 1] = _tdb_identifier(disordered_phase)
            break
    parsed = dict(command.parsed)
    parsed["action"] = "DIS_PART"
    parsed["disordered_phase"] = _tdb_identifier(disordered_phase)
    return replace(command, args=args, parsed=parsed)


def _order_disorder_alias_mode_protected_phases(
    phases: list[Phase],
    mode: str,
) -> set[str]:
    """Return helper/canonical aliases that duplicate suppression must keep."""
    if mode == "canonical":
        return set()
    phase_names = {_tdb_identifier(phase.name) for phase in phases}
    protected: set[str] = set()
    for canonical, helper in _ORDER_DISORDER_HELPER_NAMES.items():
        if canonical in phase_names and helper in phase_names:
            protected.update({canonical, helper})
    return protected


def _without_generated_compatibility_zero_parameters(
    parameters: list[Parameter],
) -> list[Parameter]:
    """Drop legacy parser-created non-vacancy zero filler from TDB export.

    Missing non-vacancy endmembers are runtime implicit zero-compatibility
    records when Equilipy reads a TDB.  Older DatabaseIR payloads may still
    contain parser-created filler, so keep source records and writer-created
    all-vacancy safety records, but do not serialize legacy ordinary filler.
    """
    return [
        parameter
        for parameter in parameters
        if not (
            str(parameter.metadata.get("generated_by", "")) == "tdb_parser"
            and parameter.metadata.get("reason")
            in {"missing_endmember", "missing_ordered_dispart_endmember"}
        )
    ]


def _factsage_export_inputs(
    phases: list[Phase],
    parameters: list[Parameter],
    commands: list[TdbCommand],
    phase_alias_map: dict[str, str],
) -> tuple[list[Phase], list[Parameter], list[TdbCommand], dict[str, str]]:
    """Return export inputs adjusted for FactSage-style TDB import."""
    phase_name_map = _factsage_phase_name_map(phases, phase_alias_map)
    export_phases: list[Phase] = []
    phase_site_scales: dict[str, float] = {}

    for phase in phases:
        original_name = _tdb_identifier(phase.name)
        export_name = phase_name_map.get(original_name, original_name)
        constituents, site_scale = _factsage_integerized_constituents(
            phase.constituents
        )
        metadata = dict(phase.metadata)
        if export_name != original_name:
            metadata["factsage_phase_name_source"] = original_name
        if not _float_close(site_scale, 1.0):
            metadata["factsage_site_ratio_scale"] = _tdb_number(site_scale)
        export_phases.append(
            replace(
                phase,
                name=export_name,
                constituents=constituents,
                metadata=metadata,
            )
        )
        phase_site_scales[export_name] = site_scale

    export_parameters = [
        _factsage_parameter_export(parameter, phase_name_map, phase_site_scales)
        for parameter in parameters
    ]
    export_commands = [
        _factsage_tdb_command_export(command, phase_name_map) for command in commands
    ]
    export_alias_map = {
        phase_name_map.get(_tdb_identifier(alias), _factsage_phase_identifier(alias)):
        phase_name_map.get(
            _tdb_identifier(canonical),
            _factsage_phase_identifier(canonical),
        )
        for alias, canonical in phase_alias_map.items()
    }
    return export_phases, export_parameters, export_commands, export_alias_map


def _factsage_phase_name_map(
    phases: list[Phase],
    phase_alias_map: dict[str, str],
) -> dict[str, str]:
    """Return unique phase-name aliases legal for FactSage ConvTDB."""
    mapping: dict[str, str] = {}
    used: set[str] = set()
    for phase in phases:
        original = _tdb_identifier(phase.name)
        sanitized = _factsage_unique_phase_identifier(
            _factsage_phase_identifier(original),
            used,
        )
        mapping[original] = sanitized
        used.add(sanitized)
    for name in set(phase_alias_map) | set(phase_alias_map.values()):
        original = _tdb_identifier(name)
        mapping.setdefault(original, _factsage_phase_identifier(original))
    return mapping


def _factsage_phase_identifier(name: str) -> str:
    """Return a phase identifier with characters accepted by FactSage."""
    identifier = _tdb_identifier(name)
    sanitized = re.sub(r"[^A-Z0-9_]", "", identifier)
    sanitized = re.sub(r"_+", "_", sanitized).strip("_")
    if not sanitized:
        sanitized = "PHASE"
    if sanitized[:1].isdigit():
        sanitized = f"PHASE_{sanitized}"
    return sanitized


def _factsage_unique_phase_identifier(name: str, used: set[str]) -> str:
    """Return a unique sanitized phase identifier."""
    candidate = name
    index = 2
    while candidate in used:
        candidate = f"{name}_{index}"
        index += 1
    return candidate


def _factsage_integerized_constituents(
    constituents: list[ConstituentSet],
) -> tuple[list[ConstituentSet], float]:
    """Return site ratios suitable for FactSage plus the energy scale factor."""
    if not constituents:
        return list(constituents), 1.0
    ratios = [float(constituent.site_ratio) for constituent in constituents]
    if any(ratio <= 0.0 for ratio in ratios):
        return list(constituents), 1.0

    integer_ratios = _factsage_reduced_integer_ratios(ratios)
    if integer_ratios is None:
        integer_ratios = _factsage_simple_fraction_ratios(ratios)
    if integer_ratios is None:
        return list(constituents), 1.0

    scale = _factsage_site_ratio_scale(ratios, integer_ratios)
    if _float_close(scale, 1.0):
        return list(constituents), 1.0
    export_constituents = [
        replace(constituent, site_ratio=float(site_ratio))
        for constituent, site_ratio in zip(constituents, integer_ratios, strict=True)
    ]
    return export_constituents, scale


def _factsage_reduced_integer_ratios(ratios: list[float]) -> list[int] | None:
    """Return reduced integer ratios when all source ratios are integers."""
    integers: list[int] = []
    for ratio in ratios:
        rounded = round(ratio)
        if not _float_close(ratio, float(rounded)):
            return None
        integers.append(int(rounded))
    divisor = _positive_gcd(integers)
    if divisor <= 1:
        return None
    return [value // divisor for value in integers]


def _factsage_simple_fraction_ratios(ratios: list[float]) -> list[int] | None:
    """Return small integer ratios for simple decimal site fractions."""
    for max_denominator, max_scale in (
        (
            _FACTSAGE_SIMPLE_SITE_RATIO_MAX_DENOMINATOR,
            _FACTSAGE_SIMPLE_SITE_RATIO_MAX_SCALE,
        ),
        (
            _FACTSAGE_EXACT_SITE_RATIO_MAX_DENOMINATOR,
            _FACTSAGE_EXACT_SITE_RATIO_MAX_SCALE,
        ),
    ):
        integer_ratios = _factsage_fraction_ratios(
            ratios,
            max_denominator=max_denominator,
            max_scale=max_scale,
        )
        if integer_ratios is not None:
            return integer_ratios
    return None


def _factsage_fraction_ratios(
    ratios: list[float],
    *,
    max_denominator: int,
    max_scale: float,
) -> list[int] | None:
    """Return integer ratios from bounded rational approximations."""
    fractions: list[Fraction] = []
    for ratio in ratios:
        fraction = Fraction(ratio).limit_denominator(max_denominator)
        if abs(float(fraction) - ratio) > _FACTSAGE_SITE_RATIO_TOLERANCE:
            return None
        fractions.append(fraction)

    denominator_lcm = reduce(_integer_lcm, (item.denominator for item in fractions), 1)
    integers = [
        item.numerator * (denominator_lcm // item.denominator)
        for item in fractions
    ]
    divisor = _positive_gcd(integers)
    if divisor > 1:
        integers = [value // divisor for value in integers]
    scale = _factsage_site_ratio_scale(ratios, integers)
    if scale > max_scale + 1e-8:
        return None
    pairs = zip(integers, ratios, strict=False)
    if all(_float_close(float(value), ratio) for value, ratio in pairs):
        return None
    return integers


def _factsage_site_ratio_scale(
    source_ratios: list[float],
    export_ratios: list[int],
) -> float:
    """Return the common multiplier from source to exported site ratios."""
    return sum(float(value) for value in export_ratios) / sum(source_ratios)


def _positive_gcd(values: list[int]) -> int:
    """Return the greatest common divisor of positive integer values."""
    positives = [abs(value) for value in values if value != 0]
    if not positives:
        return 1
    return reduce(gcd, positives)


def _integer_lcm(left: int, right: int) -> int:
    """Return least common multiple for positive integers."""
    if left == 0 or right == 0:
        return 0
    return abs(left * right) // gcd(left, right)


def _factsage_parameter_export(
    parameter: Parameter,
    phase_name_map: dict[str, str],
    phase_site_scales: dict[str, float],
) -> Parameter:
    """Return one parameter renamed and scaled for FactSage-style export."""
    original_phase = _tdb_identifier(parameter.phase_name)
    export_phase = phase_name_map.get(
        original_phase,
        _factsage_phase_identifier(original_phase),
    )
    metadata = dict(parameter.metadata)
    if export_phase != original_phase:
        metadata["factsage_phase_name_source"] = original_phase
    export_parameter = replace(parameter, phase_name=export_phase, metadata=metadata)

    site_scale = phase_site_scales.get(export_phase, 1.0)
    if _float_close(site_scale, 1.0):
        return export_parameter
    if _tdb_parameter_type(export_parameter) not in _FACTSAGE_ENERGY_PARAMETER_TYPES:
        return export_parameter

    metadata = dict(export_parameter.metadata)
    metadata.pop("source_expression", None)
    metadata["factsage_site_ratio_scale"] = _tdb_number(site_scale)
    expression = _factsage_scaled_parameter_expression(
        _parameter_source_expression_for_export(export_parameter)
        or export_parameter.expression,
        site_scale,
    )
    return replace(export_parameter, expression=expression, metadata=metadata)


def _factsage_scaled_parameter_expression(expression: str, scale: float) -> str:
    """Return a parameter expression scaled with transformed site ratios."""
    normalized = _normalize_tdb_expression(expression)
    ranges = _parse_function_ranges(normalized)
    if not ranges:
        return _scaled_gibbs_segment(normalized or "0", scale)
    return _canonical_tdb_range_expression(
        [
            replace(
                gibbs_range,
                Gibbs=_scaled_gibbs_segment(gibbs_range.Gibbs, scale),
            )
            for gibbs_range in ranges
        ],
        zero_as_function=False,
    )


def _scaled_gibbs_segment(expression: str, scale: float) -> str:
    """Scale one Gibbs-expression segment with existing canonical rules."""
    if _float_close(scale, 1.0):
        return _canonical_tdb_gibbs_expression(expression or "0")
    if _float_close(scale, 0.0):
        return "0"
    return _canonical_tdb_gibbs_expression(
        f"{_tdb_number(scale)}*({expression or '0'})"
    )


def _factsage_tdb_command_export(
    command: TdbCommand,
    phase_name_map: dict[str, str],
) -> TdbCommand:
    """Return a setup command with renamed phase references."""
    args = [
        phase_name_map.get(_tdb_identifier(arg), arg)
        for arg in command.args
    ]
    parsed = {
        key: (
            phase_name_map.get(
                _tdb_identifier(value),
                _factsage_phase_identifier(value),
            )
            if key in {"ordered_phase", "disordered_phase"}
            else value
        )
        for key, value in command.parsed.items()
    }
    return replace(command, args=args, parsed=parsed)


def _float_close(left: float, right: float, *, tolerance: float = 1e-8) -> bool:
    return abs(float(left) - float(right)) <= tolerance


def _database_header(database: DatabaseIR) -> list[str]:
    """Return an AlCuMgSi-style database metadata header."""
    title = _database_metadata_value(
        database,
        "title",
        "database_title",
    ) or "Database generated by Equilipy."
    program = _database_metadata_value(database, "program") or "Equilipy"
    authors = _database_metadata_value(database, "authors", "author")
    date_modified = _database_metadata_value(
        database,
        "date_modified",
        "date",
        "modified_date",
    ) or date.today().strftime("%m/%d/%Y")
    description = _database_metadata_value(database, "description")

    lines = [
        _comment_rule("="),
        "$",
        _comment_title(title),
        "$ ",
        _comment_rule("="),
        _comment_rule("-"),
    ]
    lines.append(
        _comment_metadata_line(
            "Program",
            program,
            label_width=_METADATA_LABEL_WIDTH,
        )
    )
    if authors:
        lines.append(
            _comment_metadata_line(
                "Author",
                authors,
                label_width=_METADATA_LABEL_WIDTH,
            )
        )
    lines.append(
        _comment_metadata_line(
            "Date modified",
            date_modified,
            label_width=_METADATA_LABEL_WIDTH,
        )
    )
    if description:
        description_lines = _wrap_comment_text(
            description,
            first_prefix="Description:".ljust(_METADATA_LABEL_WIDTH),
        )
        lines.extend(description_lines)
    lines.append("$")
    lines.append(_comment_rule("-"))
    return lines


def _database_metadata_value(database: DatabaseIR, *keys: str) -> str:
    for key in keys:
        value = str(database.metadata.get(key, "")).strip()
        if value:
            return value
    return ""


def _comment_title(text: str) -> str:
    normalized = " ".join(str(text).split())
    if len(normalized) >= _COMMENT_BODY_WIDTH:
        return _comment_line(normalized)
    left_pad = max(0, (_COMMENT_BODY_WIDTH - len(normalized)) // 2)
    return "$ " + " " * left_pad + normalized


def _comment_metadata_line(label: str, value: str, *, label_width: int) -> str:
    prefix = f"{label}:".ljust(label_width)
    return _comment_preformatted_line(f"{prefix}{value}")


def _wrap_comment_text(text: str, *, first_prefix: str = "") -> list[str]:
    """Wrap comment metadata without emitting non-comment continuation lines."""
    words = " ".join(str(text).split()).split()
    if not words:
        return []
    lines: list[str] = []
    prefix = first_prefix
    current = prefix
    for word in words:
        if current.endswith(" "):
            candidate = f"{current}{word}"
        else:
            candidate = f"{current} {word}".rstrip() if current else word
        if len(candidate) <= _COMMENT_BODY_WIDTH:
            current = candidate
            continue
        if current:
            lines.append("$ " + current)
        prefix = ""
        current = word
    if current:
        lines.append("$ " + current)
    return lines


def _comment_rule(character: str) -> str:
    return "$ " + character * _COMMENT_BODY_WIDTH


def _section_start(title: str, *headers: str) -> list[str]:
    lines = [
        _comment_arrow(">"),
        _comment_line(f"Define {title}:"),
        _comment_separator(),
    ]
    for header in headers:
        lines.append(_comment_preformatted_line(header))
    if headers:
        lines.append(_comment_separator())
    return lines


def _section_end() -> list[str]:
    return [_comment_arrow("<")]


def _element_section(elements: list[Element]) -> list[str]:
    lines = [
        _comment_arrow(">"),
        _comment_line("Define ELEMENTS:"),
        _comment_separator(),
        (
            "$ELEMENT <Symbol> <Phase SER.x24>          "
            "<molar_mass>    <H298>      <S298>!"
        ),
        _comment_separator(),
    ]
    for element in _ordered_elements(_elements_with_system_defaults(elements)):
        lines.append(format_tdb_element(element))
    lines.extend(_section_end())
    return lines


def _elements_with_system_defaults(elements: list[Element]) -> list[Element]:
    """Return elements plus pseudo-elements required by the system defaults."""
    by_symbol = {_tdb_identifier(element.symbol): element for element in elements}
    output = list(elements)
    if "/-" not in by_symbol:
        output.append(
            Element(
                symbol="/-",
                atomic_mass=0.0,
                reference_state="ELECTRON_GAS",
                h298=0.0,
                s298=0.0,
                tdb_symbol="/-",
            )
        )
    if "VA" not in by_symbol:
        output.append(
            Element(
                symbol="VA",
                atomic_mass=0.0,
                reference_state="VACUUM",
                h298=0.0,
                s298=0.0,
                tdb_symbol="VA",
            )
        )
    return output


def _ordered_elements(elements: list[Element]) -> list[Element]:
    special_order = {"/-": 0, "VA": 1, "*": 2}

    def sort_key(element: Element) -> tuple[int, int, str]:
        symbol = _tdb_identifier(element.symbol)
        if symbol in special_order:
            return (0, special_order[symbol], symbol)
        return (1, 0, symbol)

    return sorted(elements, key=sort_key)


def format_tdb_element(element: Element) -> str:
    """Format one ELEMENT command with stable alignment."""
    symbol = _tdb_element_symbol(element.symbol)
    reference_state = _tdb_identifier(element.reference_state)
    return (
        f" ELEMENT {symbol:>7}   {reference_state:<24}"
        f"{_tdb_scientific_number(element.atomic_mass):>10}  "
        f"{_tdb_scientific_number(element.h298):>10}  "
        f"{_tdb_scientific_number(element.s298):>10} !"
    )


def _function_section(functions: list[FunctionDefinition]) -> list[str]:
    lines = _section_start(
        "FUNCTIONS",
        "FUNCTION Name           LowTemp(K) Equation;        HighTemp(K) Y/N",
    )
    for function in _ordered_functions(functions):
        lines.extend(format_tdb_function(function))
    lines.extend(_section_end())
    return lines


def _functions_with_zero_for_parameter_export(
    functions: list[FunctionDefinition],
    parameters: list[Parameter],
) -> list[FunctionDefinition]:
    """Ensure ``ZERO#`` is defined when zero PAR segments are exported."""
    if not _parameters_need_zero_function(parameters):
        return functions
    existing = {_tdb_identifier(function.name) for function in functions}
    if "ZERO" in existing:
        return functions
    return [
        *functions,
        FunctionDefinition(
            name="ZERO",
            expression="",
            gibbs_ranges=[GibbsRange(298.15, 6000.0, "0", "N")],
        ),
    ]


def _parameters_need_zero_function(parameters: list[Parameter]) -> bool:
    """Return whether thermodynamic parameters contain literal zero segments."""
    return any(
        _parameter_contains_zero_gibbs_segment(parameter)
        for parameter in parameters
    )


def _parameter_contains_zero_gibbs_segment(parameter: Parameter) -> bool:
    """Return whether one G/L parameter has a segment that canonicalizes to zero."""
    parameter_type = _tdb_parameter_type(parameter)
    if parameter_type not in {"G", "L"}:
        return False
    expression = _normalize_tdb_expression(parameter.expression)
    ranges = _parse_function_ranges(expression)
    if not ranges:
        ranges = [GibbsRange(298.15, 6000.0, expression or "0", "N")]
    for gibbs_range in ranges:
        if _canonical_tdb_gibbs_expression(gibbs_range.Gibbs) == "0":
            return True
    return False


def format_tdb_function(function: FunctionDefinition) -> list[str]:
    """Format one FUNCTION command using continuation lines between ranges."""
    name = _tdb_function_name(function.name)
    range_lines = _canonical_function_range_lines(function)

    lines: list[str] = []
    for index, range_line in enumerate(range_lines):
        prefix = _function_prefix(name) if index == 0 else _CONTINUATION_PREFIX
        suffix = " !" if index == len(range_lines) - 1 else ""
        lines.extend(_wrap_tdb_command(prefix, f"{range_line}{suffix}"))
    return lines


def _tdb_setup_section(
    commands: list[TdbCommand],
    code_map: dict[str, str],
    elements: list[Element],
    phase_names: set[str],
    phase_alias_map: dict[str, str],
    *,
    factsage_style: bool = False,
) -> list[str]:
    """Return preserved system default/type-definition commands."""
    type_definition_keyword = "TYPE_DEFINITION" if factsage_style else "TYPE_DEF"
    lines = _section_start(
        "SYSTEM DEFAULTS, PHASE TYPE DEFINITIONS, AND ORDRER-DISORDER",
        "Command",
    )
    lines.extend(_wrap_tdb_command(" ", f"{type_definition_keyword} % SEQ * !"))
    lines.extend(_wrap_tdb_command(" ", "DEFINE_SYSTEM_DEFAULT SPECIE 2 !"))
    default_elements = _default_system_element_command(elements)
    lines.extend(_wrap_tdb_command(" ", f"{default_elements} !"))
    for args in _canonical_type_definition_args(
        commands,
        code_map,
        phase_names,
        phase_alias_map,
    ):
        if factsage_style:
            args = [
                "DISORDER_PART" if arg.upper() == "DIS_PART" else arg
                for arg in args
            ]
        text = " ".join([type_definition_keyword, *args])
        lines.extend(_wrap_tdb_command(" ", f"{text} !"))
    lines.extend(_section_end())
    return lines


def _default_system_element_command(elements: list[Element]) -> str:
    """Return the default system element command for exported pseudo-elements."""
    element_symbols = {
        _tdb_identifier(element.symbol)
        for element in _elements_with_system_defaults(elements)
    }
    default_symbols = []
    if "VA" in element_symbols:
        default_symbols.append("VA")
    if "/-" in element_symbols:
        default_symbols.append("/-")
    return " ".join(["DEFAULT_COMMAND", "DEF_SYS_ELEMENT", *default_symbols])


def format_tdb_setup_command(
    command: TdbCommand,
    code_map: dict[str, str] | None = None,
    phase_names: set[str] | None = None,
    phase_alias_map: dict[str, str] | None = None,
) -> list[str]:
    """Format one preserved setup/type command."""
    keyword = _tdb_setup_keyword(command.command)
    args = _tdb_setup_args(
        command,
        code_map or {},
        phase_names or set(),
        phase_alias_map or {},
    )
    text = " ".join([keyword, *args]).strip()
    prefix = " " if command.active else "$"
    return _wrap_tdb_command(prefix, f"{text} !")


def _tdb_setup_keyword(command: str) -> str:
    """Return the canonical TDB setup command keyword."""
    upper = _tdb_identifier(command).replace("_", "-")
    if upper == "DEFINE-SYSTEM-DEFAULT":
        return "DEFINE_SYSTEM_DEFAULT"
    if upper in {"DEFAULT-COM", "DEFAULT-COMMAND"}:
        return "DEFAULT_COMMAND"
    if upper in {"TYPE-DEF", "TYPE-DEFINITION"}:
        return "TYPE_DEF"
    return upper


def _tdb_setup_args(
    command: TdbCommand,
    code_map: dict[str, str],
    phase_names: set[str],
    phase_alias_map: dict[str, str] | None = None,
) -> list[str]:
    """Return setup/type-definition arguments in the selected export dialect."""
    args = list(command.args)
    if _tdb_setup_keyword(command.command) != "TYPE_DEF" or not args:
        return args

    original_code = args[0].upper()
    if original_code in code_map:
        args[0] = code_map[original_code]

    upper_args = [arg.upper() for arg in args]
    for index, token in enumerate(upper_args):
        if token == "DIS_PART":
            args[index] = "DIS_PART"
        elif token == "DISORDER_PART":
            args[index] = "DIS_PART"

    if "DIS_PART" in [arg.upper() for arg in args]:
        disorder_index = [arg.upper() for arg in args].index("DIS_PART")
        if disorder_index + 1 < len(args):
            args[disorder_index + 1] = _disordered_phase_name_for_export(
                args[disorder_index + 1],
                phase_names,
                phase_alias_map or {},
            )
    return args


def _tdb_setup_type_code_export_map(commands: list[TdbCommand]) -> dict[str, str]:
    """Map AMEND_PHASE_DESCRIPTION type codes to compact alphabetical codes."""
    code_map: dict[str, str] = {}
    group_codes: dict[tuple[str, ...], str] = {}
    available_codes = iter("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    for command in commands:
        if _tdb_setup_keyword(command.command) != "TYPE_DEF":
            continue
        group_key = _type_definition_group_key(command)
        if not group_key:
            continue
        code = command.parsed.get("code") or (command.args[0] if command.args else "")
        if not code:
            continue
        if group_key not in group_codes:
            group_codes[group_key] = next(available_codes)
        code_map[code.upper()] = group_codes[group_key]
    return code_map


def _magnetic_type_codes(
    commands: list[TdbCommand],
    code_map: dict[str, str],
) -> set[str]:
    """Return exported TYPE_DEF codes that activate magnetic models."""
    codes: set[str] = set()
    for command in commands:
        if _tdb_setup_keyword(command.command) != "TYPE_DEF":
            continue
        if command.parsed.get("action", "").upper() != "MAGNETIC":
            continue
        code = command.parsed.get("code") or (command.args[0] if command.args else "")
        if code:
            codes.add(code_map.get(code.upper(), code.upper()))
    return codes


def _phases_with_magnetic_parameters(parameters: list[Parameter]) -> set[str]:
    """Return phase names with explicit BM/BMAGN or TC parameters."""
    magnetic_parameter_types = {"BM", "BMAGN", "TC"}
    return {
        _tdb_identifier(parameter.phase_name)
        for parameter in parameters
        if _tdb_identifier(parameter.parameter_type) in magnetic_parameter_types
    }


def _canonical_type_definition_args(
    commands: list[TdbCommand],
    code_map: dict[str, str],
    phase_names: set[str],
    phase_alias_map: dict[str, str],
    *,
    include_actions: set[str] | None = None,
    exclude_actions: set[str] | None = None,
) -> list[list[str]]:
    """Return canonical alphabetical TYPE_DEF argument lists."""
    definitions: dict[str, list[str]] = {}
    for command in commands:
        if _tdb_setup_keyword(command.command) != "TYPE_DEF":
            continue
        group_key = _type_definition_group_key(command)
        if not group_key:
            continue
        action = group_key[0]
        if include_actions is not None and action not in include_actions:
            continue
        if exclude_actions is not None and action in exclude_actions:
            continue
        original_code = (
            command.parsed.get("code") or (command.args[0] if command.args else "")
        ).upper()
        code = code_map.get(original_code)
        if not code or code in definitions:
            continue
        definitions[code] = _canonical_type_definition_for_group(
            code,
            group_key,
            phase_names,
            phase_alias_map,
        )
    return [definitions[code] for code in sorted(definitions)]


def _type_definition_group_key(command: TdbCommand) -> tuple[str, ...]:
    """Return the semantic setup group for a TYPE_DEF command."""
    action = command.parsed.get("action", "").upper()
    if action == "MAGNETIC":
        model = command.parsed.get("magnetic_model", "")
        factor = command.parsed.get("magnetic_factor", "")
        if model and factor:
            return ("MAGNETIC", model, factor)
    if action == "DIS_PART":
        ordered = _tdb_identifier(command.parsed.get("ordered_phase", ""))
        disordered = _tdb_identifier(command.parsed.get("disordered_phase", ""))
        if ordered and disordered:
            return ("DIS_PART", ordered, disordered)
    if action in {"SEQ", ""}:
        return ()
    args = tuple(_tdb_identifier(arg) for arg in command.args[1:])
    return ("RAW", *args) if args else ()


def _canonical_type_definition_for_group(
    code: str,
    group_key: tuple[str, ...],
    phase_names: set[str],
    phase_alias_map: dict[str, str],
) -> list[str]:
    """Return one canonical TYPE_DEF payload."""
    if group_key[0] == "MAGNETIC":
        return [
            code,
            "GES",
            "AMEND_PHASE_DESCRIPTION",
            "@",
            "MAGNETIC",
            group_key[1],
            group_key[2],
        ]
    if group_key[0] == "DIS_PART":
        return [
            code,
            "GES",
            "AMEND_PHASE_DESCRIPTION",
            group_key[1],
            "DIS_PART",
            _disordered_phase_name_for_export(
                group_key[2],
                phase_names,
                phase_alias_map,
            ),
        ]
    return [code, *group_key[1:]]


def _disordered_phase_name_for_export(
    name: str,
    phase_names: set[str],
    phase_alias_map: dict[str, str] | None = None,
) -> str:
    """Return the exported disordered phase name from a TYPE_DEF record."""
    normalized = _tdb_identifier(name)
    mapped = (phase_alias_map or {}).get(normalized)
    if mapped:
        return mapped
    if normalized in phase_names:
        return normalized
    canonical = canonical_disordered_phase_name(normalized)
    if canonical and canonical in phase_names:
        return canonical
    return normalized


def _duplicate_phase_alias_map_for_export(
    phases: list[Phase],
    parameters: list[Parameter],
    commands: list[TdbCommand],
    *,
    extra_protected_phases: set[str] | None = None,
) -> dict[str, str]:
    """Return thermodynamic duplicate phases to suppress on export."""
    phase_names = {_tdb_identifier(phase.name) for phase in phases}
    preferred_phases = _preferred_disordered_phase_names_for_export(
        commands,
        phase_names,
    )
    protected_phases = _order_disorder_protected_phase_names(commands, phase_names)
    protected_phases.update(extra_protected_phases or set())
    aliases = _duplicate_phase_aliases_from_signatures(
        phases,
        parameters,
        preferred_phases,
        protected_phases,
    )
    return aliases


def _preferred_disordered_phase_names_for_export(
    commands: list[TdbCommand],
    phase_names: set[str],
) -> set[str]:
    """Return canonical disordered phase names referenced by TYPE_DEF records."""
    names: set[str] = set()
    for command in commands:
        if _tdb_setup_keyword(command.command) != "TYPE_DEF":
            continue
        if command.parsed.get("action", "").upper() != "DIS_PART":
            continue
        disordered_phase = command.parsed.get("disordered_phase", "")
        names.add(_disordered_phase_name_for_export(disordered_phase, phase_names))
    return names


def _duplicate_phase_aliases_from_signatures(
    phases: list[Phase],
    parameters: list[Parameter],
    preferred_phases: set[str],
    protected_phases: set[str],
) -> dict[str, str]:
    """Return duplicate phases with identical constituents and parameters."""
    parameters_by_phase = _parameters_by_phase(parameters)
    groups: dict[tuple, list[str]] = {}
    phase_order: list[str] = []
    for phase in phases:
        phase_name = _tdb_identifier(phase.name)
        if phase_name in phase_order:
            continue
        parameter_signatures = _phase_parameter_signatures(
            parameters_by_phase.get(phase_name, []),
        )
        if not parameter_signatures:
            continue
        signature = (
            _phase_constituent_signature(phase),
            tuple(parameter_signatures),
        )
        groups.setdefault(signature, []).append(phase_name)
        phase_order.append(phase_name)

    aliases: dict[str, str] = {}
    order_index = {name: index for index, name in enumerate(phase_order)}
    for names in groups.values():
        if len(names) < 2:
            continue
        canonical = _preferred_duplicate_phase_name(
            names,
            preferred_phases,
            order_index,
        )
        for name in names:
            if name == canonical:
                continue
            if name not in protected_phases:
                aliases[name] = canonical
    return aliases


def _order_disorder_protected_phase_names(
    commands: list[TdbCommand],
    phase_names: set[str],
) -> set[str]:
    """Return order-disorder phases that must not be duplicate-suppressed."""
    protected: set[str] = set()
    for command in commands:
        if _tdb_setup_keyword(command.command) != "TYPE_DEF":
            continue
        if command.parsed.get("action", "").upper() != "DIS_PART":
            continue
        ordered_phase = _tdb_identifier(command.parsed.get("ordered_phase", ""))
        disordered_phase = _disordered_phase_name_for_export(
            command.parsed.get("disordered_phase", ""),
            phase_names,
        )
        if ordered_phase in phase_names:
            protected.add(ordered_phase)
        if disordered_phase in phase_names:
            protected.add(disordered_phase)
    return protected


def _parameters_by_phase(parameters: list[Parameter]) -> dict[str, list[Parameter]]:
    """Return parameters keyed by normalized phase name."""
    by_phase: dict[str, list[Parameter]] = {}
    for parameter in parameters:
        by_phase.setdefault(_tdb_identifier(parameter.phase_name), []).append(parameter)
    return by_phase


def _phase_constituent_signature(phase: Phase) -> tuple:
    """Return a normalized constituent layout signature for duplicate checks."""
    return tuple(
        (
            _tdb_number(constituent.site_ratio),
            tuple(sorted(_tdb_constituent_name(name) for name in constituent.species)),
        )
        for constituent in phase.constituents
    )


def _phase_parameter_signatures(parameters: list[Parameter]) -> list[tuple]:
    """Return normalized parameter signatures for one phase."""
    signatures: list[tuple] = []
    for parameter in parameters:
        parameter_type = _tdb_parameter_type(parameter)
        signatures.append(
            (
                parameter_type,
                tuple(_tdb_identifier(target) for target in parameter.target),
                int(parameter.order),
                _canonical_tdb_parameter_expression(parameter, parameter_type),
            )
        )
    return sorted(signatures)


def _preferred_duplicate_phase_name(
    names: list[str],
    preferred_phases: set[str],
    order_index: dict[str, int],
) -> str:
    """Return the canonical phase name to keep from a duplicate group."""
    known_canonical = [
        canonical
        for alias, canonical in DISORDERED_PHASE_CANONICAL_NAMES.items()
        if alias in names and canonical in names
    ]
    if known_canonical:
        return min(known_canonical, key=lambda name: order_index[name])
    preferred = [name for name in names if name in preferred_phases]
    if preferred:
        return min(preferred, key=lambda name: order_index[name])
    return min(names, key=lambda name: order_index[name])


def _phase_section(
    phases: list[Phase],
    code_map: dict[str, str],
    phase_alias_map: dict[str, str],
    magnetic_codes: set[str],
    phases_with_magnetic_parameters: set[str],
) -> list[str]:
    """Return PHASE/CONST commands for phase and constituent definitions."""
    if not phases:
        return []

    lines = _section_start(
        "PHASES",
        "PHASE <name> <type-codes> <sublattice-count> <site-ratios>",
        "CONST <name> :<sublattice constituents>:",
    )
    wrote_phase = False
    for phase in phases:
        if _tdb_identifier(phase.name) in phase_alias_map:
            continue
        if wrote_phase:
            lines.append("")
        lines.extend(
            format_tdb_phase(
                phase,
                code_map,
                magnetic_codes=magnetic_codes,
                phase_has_magnetic_parameters=(
                    _tdb_identifier(phase.name) in phases_with_magnetic_parameters
                ),
            )
        )
        wrote_phase = True
    lines.extend(_section_end())
    return lines


def format_tdb_phase(
    phase: Phase,
    code_map: dict[str, str] | None = None,
    *,
    magnetic_codes: set[str] | None = None,
    phase_has_magnetic_parameters: bool = True,
) -> list[str]:
    """Format one phase as PHASE followed by CONST."""
    phase_name = _tdb_phase_name(phase)
    type_codes = _phase_type_codes(
        phase,
        code_map or {},
        magnetic_codes or set(),
        phase_has_magnetic_parameters,
    )
    constituents = list(phase.constituents)
    sublattice_count = len(constituents) if constituents else 1
    site_ratios = (
        [_tdb_site_ratio_number(constituent.site_ratio) for constituent in constituents]
        if constituents
        else ["1"]
    )
    phase_command = (
        f"PHASE {phase_name} {type_codes} {sublattice_count} "
        f"{' '.join(site_ratios)} !"
    )
    constituent_command = (
        f"CONST {phase_name} {_constituent_tdb_body(constituents, phase_name)} !"
    )
    lines = _wrap_tdb_command(" ", phase_command)
    lines.extend(_wrap_tdb_command(" ", constituent_command))
    return lines


def _tdb_phase_name(phase: Phase) -> str:
    """Return the exported phase name, preserving source suffixes like ``:L``."""
    source_name = _source_phase_name(phase.source.command)
    if source_name:
        bare_source_name = source_name.split(":", 1)[0]
        if bare_source_name.upper() == str(phase.name).strip().upper():
            return _tdb_identifier(source_name)
    return _tdb_identifier(phase.name)


def _source_phase_name(command: str) -> str:
    """Return the original PHASE name token when available."""
    match = re.search(r"(?i)(?:^|\s)PHASE\s+(\S+)", str(command or ""))
    return match.group(1).strip() if match is not None else ""


def _phase_type_codes(
    phase: Phase,
    code_map: dict[str, str],
    magnetic_codes: set[str],
    phase_has_magnetic_parameters: bool,
) -> str:
    """Return TDB type-code token from source, falling back to ``%``."""
    match = re.search(
        r"(?i)(?:^|\s)PHASE\s+\S+\s+(\S+)",
        str(phase.source.command or ""),
    )
    token = match.group(1).strip() if match is not None else ""
    if not token.startswith("%"):
        return "%"
    return _map_phase_type_codes(
        _tdb_identifier(token),
        code_map,
        magnetic_codes,
        phase_has_magnetic_parameters,
    )


def _map_phase_type_codes(
    token: str,
    code_map: dict[str, str],
    magnetic_codes: set[str],
    phase_has_magnetic_parameters: bool,
) -> str:
    """Apply setup type-code remapping and skip ineffective magnetic codes."""
    if not token.startswith("%"):
        return token
    remapped_codes = [code_map.get(character, character) for character in token[1:]]
    exported_codes = [
        code
        for code in remapped_codes
        if phase_has_magnetic_parameters or code not in magnetic_codes
    ]
    return "%" + "".join(exported_codes)


def _constituent_tdb_body(constituents: list, phase_name: str) -> str:
    """Return the colon-delimited CONST payload."""
    if not constituents:
        return f": {_tdb_identifier(phase_name)} :"
    sections: list[str] = []
    for constituent in constituents:
        species = [
            _tdb_constituent_name(name)
            for name in constituent.species
            if str(name).strip()
        ]
        sections.append(" ".join(species) if species else "*")
    return ": " + " : ".join(sections) + " :"


def _tdb_constituent_name(name: str) -> str:
    text = str(name).strip()
    if not text:
        return text
    upper = text.upper()
    if upper in _SPECIAL_ELEMENT_SYMBOLS:
        return upper
    if upper.isalpha() and len(upper) <= 2:
        return _conventional_element_symbol(upper)
    return upper


def _parameter_section(
    database: DatabaseIR,
    parameters: list[Parameter],
) -> list[str]:
    """Return grouped PAR commands with recoverable comments."""
    if not parameters:
        return []

    groups = _parameter_groups(database, parameters)
    lines = _section_start(
        "PARAMETERS",
        "PAR <type>(<phase>,<constituents>;<order>) <T_low> expression; <T_high> N",
    )
    for arity in sorted(groups):
        label = _arity_label(arity)
        lines.append(_comment_arrow(">"))
        lines.append(_comment_line(f"{label} parameters"))
        for elements in sorted(groups[arity]):
            parameters = groups[arity][elements]
            element_text = ",".join(elements) if elements else "Unknown"
            lines.append(_comment_separator())
            readable_elements = "-".join(elements) if elements else "Unknown"
            lines.append(_comment_line(readable_elements))
            lines.append(
                _comment_line(
                    f"Equilipy:PARAMETER_GROUP ARITY={arity} ELEMENTS={element_text}"
                )
            )
            lines.append(_comment_line("Comments:"))
            for parameter in parameters:
                lines.extend(format_tdb_parameter(parameter, database))
        lines.append(_comment_arrow("<"))
    lines.extend(_section_end())
    return lines


def _reference_section(references: list[TdbReference]) -> list[str]:
    """Return TDB reference table commands."""
    entries = _deduplicated_references(references)
    if not entries:
        return []

    lines = [
        _comment_arrow(">"),
        _comment_line("Define REFERENCES:"),
        _comment_separator(),
    ]
    list_entries = [
        reference
        for reference in entries
        if _reference_section_name(reference) != "ADD_REFERENCES"
    ]
    add_entries = [
        reference
        for reference in entries
        if _reference_section_name(reference) == "ADD_REFERENCES"
    ]
    if list_entries:
        lines.extend(_reference_command_block("LIST-OF-REFERENCE", list_entries))
    if add_entries:
        if list_entries:
            lines.append("")
        command = "ADD_REFERENCES" if list_entries else "LIST-OF-REFERENCE"
        lines.extend(_reference_command_block(command, add_entries))
    lines.extend(_section_end())
    return lines


def _reference_command_block(
    command: str,
    references: list[TdbReference],
) -> list[str]:
    """Return one LIST-OF-REFERENCE or ADD_REFERENCES command block."""
    lines = [f" {command}"]
    if command == "LIST-OF-REFERENCE":
        lines.append(" NUMBER SOURCE")
    for reference in references:
        lines.extend(format_tdb_reference(reference))
    lines.append("!")
    return lines


def _reference_section_name(reference: TdbReference) -> str:
    """Return the canonical source section for one reference entry."""
    section = str(reference.section).upper().replace("-", "_")
    if section.startswith("ADD_"):
        return "ADD_REFERENCES"
    return "LIST-OF-REFERENCE"


def _deduplicated_references(
    references: list[TdbReference],
) -> list[TdbReference]:
    """Return source-order references with later duplicate keys winning."""
    by_key: dict[str, TdbReference] = {}
    order: list[str] = []
    for reference in references:
        key = _tdb_reference_key(reference.key)
        if not key:
            continue
        lookup_key = key.upper()
        if lookup_key not in by_key:
            order.append(lookup_key)
        by_key[lookup_key] = replace(reference, key=key)
    return [by_key[key] for key in order]


def format_tdb_reference(reference: TdbReference) -> list[str]:
    """Format one TDB bibliographic reference entry."""
    key = _tdb_reference_key(reference.key)
    text = _escape_tdb_reference_text(reference.text)
    if not key:
        return []
    if not text:
        return [f"{_reference_prefix(key)}'"]

    first_prefix = _reference_prefix(key)
    continuation_prefix = "          "
    one_line = f"{first_prefix}{text}'"
    if len(one_line) <= _MAX_TDB_REFERENCE_LINE_LENGTH:
        return [one_line]

    parts = _wrap_reference_words(
        text,
        _MAX_TDB_REFERENCE_LINE_LENGTH - len(first_prefix),
        _MAX_TDB_REFERENCE_LINE_LENGTH - len(continuation_prefix) - 1,
    )
    if not parts:
        return [f"{first_prefix}'"]
    if len(parts) == 1:
        return [f"{first_prefix}{parts[0]}'"]

    lines = [f"{first_prefix}{parts[0]}"]
    lines.extend(f"{continuation_prefix}{part}" for part in parts[1:-1])
    lines.append(f"{continuation_prefix}{parts[-1]}'")
    return lines


def _reference_prefix(key: str) -> str:
    """Return Hallstedt-style reference table prefix up to the opening quote."""
    if len(key) <= 6:
        return f"  {key:<7}'"
    return f"  {key} '"


def _tdb_reference_key(value: str) -> str:
    """Return a portable TDB reference key while preserving source casing."""
    return re.sub(r"[^A-Za-z0-9_.+\-/]", "", str(value).strip())


def _escape_tdb_reference_text(value: str) -> str:
    """Escape one quoted TDB reference string."""
    return " ".join(str(value).split()).replace("'", "''")


def _wrap_reference_words(
    text: str,
    first_width: int,
    continuation_width: int,
) -> list[str]:
    """Wrap quoted reference text without splitting words when possible."""
    words = text.split()
    if not words:
        return []

    parts: list[str] = []
    width = max(8, first_width)
    current = ""
    for word in words:
        candidate = f"{current} {word}".strip()
        if current and len(candidate) > width:
            parts.append(current)
            width = max(8, continuation_width)
            current = word
            while len(current) > width:
                parts.append(current[:width])
                current = current[width:]
            continue
        while not current and len(candidate) > width:
            parts.append(candidate[:width])
            candidate = candidate[width:]
            width = max(8, continuation_width)
        current = candidate
    if current:
        parts.append(current)
    return parts


def _endmember_compound_export_inputs(
    phases: list[Phase],
    parameters: list[Parameter],
    *,
    include_endmember_compounds: bool,
    phase_alias_map: dict[str, str],
) -> tuple[list[Phase], list[Parameter]]:
    """Return export inputs with optional generated endmember compounds."""
    export_phases, export_parameters = (
        (list(phases), list(parameters))
        if include_endmember_compounds
        else _without_endmember_compound_exports(phases, parameters)
    )
    if not include_endmember_compounds:
        return export_phases, export_parameters

    generated_phases, generated_parameters = _generated_endmember_compounds(
        export_phases,
        export_parameters,
        phase_alias_map,
    )
    if not generated_phases:
        return export_phases, export_parameters
    return (
        [*export_phases, *generated_phases],
        [*export_parameters, *generated_parameters],
    )


def _without_endmember_compound_exports(
    phases: list[Phase],
    parameters: list[Parameter],
) -> tuple[list[Phase], list[Parameter]]:
    """Remove exported endmember-compound copies from the normal TDB path."""
    endmember_phase_names = {
        _tdb_identifier(phase.name)
        for phase in phases
        if _is_endmember_compound_phase(phase)
    }
    if not endmember_phase_names:
        return list(phases), list(parameters)
    return (
        [
            phase
            for phase in phases
            if _tdb_identifier(phase.name) not in endmember_phase_names
        ],
        [
            parameter
            for parameter in parameters
            if _tdb_identifier(parameter.phase_name) not in endmember_phase_names
        ],
    )


def _is_endmember_compound_phase(phase: Phase) -> bool:
    """Return whether a phase is a generated endmember-compound export copy."""
    metadata_value = str(phase.metadata.get("endmember_compound", "")).lower()
    return metadata_value in {"1", "true", "yes"} or _is_endmember_compound_name(
        phase.name
    )


def _is_endmember_compound_name(name: str) -> bool:
    """Return whether a phase name uses Equilipy's endmember-compound suffix."""
    identifier = _tdb_identifier(name)
    return identifier.endswith(_ENDMEMBER_COMPOUND_SUFFIX) or identifier.endswith(
        "_ENDMEMBER"
    )


def _generated_endmember_compounds(
    phases: list[Phase],
    parameters: list[Parameter],
    phase_alias_map: dict[str, str],
) -> tuple[list[Phase], list[Parameter]]:
    """Build stoichiometric compound phases from concrete solution endmembers."""
    existing_names = {_tdb_identifier(phase.name) for phase in phases}
    source_phase_names = {
        _tdb_identifier(phase.name)
        for phase in phases
        if not _is_endmember_compound_phase(phase)
        and _tdb_identifier(phase.name) not in phase_alias_map
        and str(phase.model).upper() != "COMPOUND"
    }
    phases_by_name = {
        _tdb_identifier(phase.name): phase
        for phase in phases
        if _tdb_identifier(phase.name) in source_phase_names
    }
    generated_phases: list[Phase] = []
    generated_parameters: list[Parameter] = []
    emitted_targets: set[tuple[str, str]] = set()

    for parameter in parameters:
        phase_name = _tdb_identifier(parameter.phase_name)
        phase = phases_by_name.get(phase_name)
        if phase is None:
            continue
        if (
            _tdb_parameter_type(parameter) != "G"
            or int(parameter.order) != 0
            or not parameter.target
        ):
            continue
        for target in _endmember_compound_source_targets(parameter, phase):
            key = (phase_name, _tdb_identifier(target))
            if key in emitted_targets:
                continue
            emitted_targets.add(key)
            stoichiometry = _endmember_compound_stoichiometry(phase, target)
            if not stoichiometry:
                continue

            compound_name = _unique_endmember_compound_name(
                phase.name,
                target,
                existing_names,
            )
            if not compound_name:
                continue
            compound_phase = Phase(
                name=compound_name,
                model="COMPOUND",
                constituents=[
                    ConstituentSet(index, [species], ratio)
                    for index, (species, ratio) in enumerate(stoichiometry, start=1)
                ],
                metadata={
                    "generated_by": _ENDMEMBER_COMPOUND_GENERATOR,
                    "endmember_compound": "true",
                    "source_phase": _tdb_identifier(phase.name),
                    "source_endmember": target,
                },
            )
            compound_target = ":".join(species for species, _ratio in stoichiometry)
            parameter_metadata = dict(parameter.metadata)
            parameter_metadata.update(
                {
                    "generated_by": _ENDMEMBER_COMPOUND_GENERATOR,
                    "endmember_compound": "true",
                    "source_phase": _tdb_identifier(phase.name),
                    "source_endmember": target,
                }
            )
            generated_phases.append(compound_phase)
            generated_parameters.append(
                replace(
                    parameter,
                    phase_name=compound_name,
                    target=[compound_target],
                    order=0,
                    metadata=parameter_metadata,
                )
            )
            existing_names.add(_tdb_identifier(compound_name))
    return generated_phases, generated_parameters


def _endmember_compound_source_targets(
    parameter: Parameter,
    phase: Phase,
) -> list[str]:
    """Return concrete source endmember targets for compound-copy generation."""
    if _parameter_target_has_wildcard(parameter):
        targets = _expanded_wildcard_targets(parameter, phase)
    else:
        targets = [_tdb_parameter_target_text(parameter.target[0])]
    return [target for target in targets if "," not in target]


def _endmember_compound_stoichiometry(
    phase: Phase,
    target: str,
) -> list[tuple[str, float]]:
    """Return normalized species/site-ratio pairs for one concrete endmember."""
    sections = str(target).split(":")
    constituents = sorted(phase.constituents, key=lambda item: item.sublattice)
    if len(sections) != len(constituents):
        return []
    ratios_by_species: dict[str, float] = {}
    order: list[str] = []
    for section, constituent in zip(sections, constituents, strict=True):
        species_names = [name.strip() for name in section.split(",") if name.strip()]
        if len(species_names) != 1:
            return []
        species = _tdb_constituent_name(species_names[0])
        if _tdb_identifier(species) in {"VA", "/-", "*"}:
            continue
        if species not in ratios_by_species:
            ratios_by_species[species] = 0.0
            order.append(species)
        ratios_by_species[species] += float(constituent.site_ratio)
    return [
        (species, ratios_by_species[species])
        for species in order
        if ratios_by_species[species] != 0.0
    ]


def _unique_endmember_compound_name(
    phase_name: str,
    target: str,
    existing_names: set[str],
) -> str:
    """Return a stable, recoverable generated phase name for one endmember."""
    target_part = "_".join(
        _tdb_identifier(token) or "EMPTY"
        for token in re.split(r"[:,]+", str(target))
        if token.strip()
    )
    base = _tdb_identifier(f"{phase_name}_{target_part}{_ENDMEMBER_COMPOUND_SUFFIX}")
    if base in existing_names:
        return ""
    name = base
    index = 2
    while name in existing_names:
        name = f"{base}_{index}"
        index += 1
    return name


def _parameter_exports(
    parameters: list[Parameter],
    phases: list[Phase],
    phase_alias_map: dict[str, str],
) -> list[Parameter]:
    """Return final parameters to write, including export safety records."""
    parameters = _parameters_with_export_phase_aliases(parameters, phase_alias_map)
    parameters = _deduplicate_parameter_exports(parameters)
    parameters = _parameters_compatible_with_phase_constituents(
        parameters,
        phases,
        phase_alias_map,
    )
    parameters = _parameters_without_symmetry_duplicate_zero_endmembers(
        parameters,
        phases,
        phase_alias_map,
    )
    return _parameters_with_missing_all_vacancy_endmembers(
        parameters,
        phases,
        phase_alias_map,
    )


def _parameters_with_expanded_wildcard_targets(
    parameters: list[Parameter],
    phases: list[Phase],
    phase_alias_map: dict[str, str],
) -> list[Parameter]:
    """Replace TDB wildcard target sections with explicit exported targets."""
    phases_by_name = {
        _tdb_identifier(phase.name): phase
        for phase in phases
        if _tdb_identifier(phase.name) not in phase_alias_map
    }
    expanded: list[Parameter] = []
    for parameter in parameters:
        if not _parameter_target_has_wildcard(parameter):
            expanded.append(parameter)
            continue
        phase = phases_by_name.get(_tdb_identifier(parameter.phase_name))
        if phase is None or not phase.constituents:
            raise ValueError(
                "Cannot export wildcard PAR target without phase constituents: "
                f"{parameter.phase_name!r}"
            )
        for target in _expanded_wildcard_targets(parameter, phase):
            metadata = dict(parameter.metadata)
            metadata["expanded_from_wildcard"] = _tdb_parameter_target(parameter)
            expanded.append(replace(parameter, target=[target], metadata=metadata))
    return _deduplicate_parameter_exports(expanded)


def _parameters_compatible_with_phase_constituents(
    parameters: list[Parameter],
    phases: list[Phase],
    phase_alias_map: dict[str, str],
) -> list[Parameter]:
    """Drop parameter records whose target species are not in the exported phase."""
    phases_by_name = {
        _tdb_identifier(phase.name): phase
        for phase in phases
        if _tdb_identifier(phase.name) not in phase_alias_map
    }
    return [
        parameter
        for parameter in parameters
        if _parameter_compatible_with_phase_constituents(
            parameter,
            phases_by_name.get(_tdb_identifier(parameter.phase_name)),
        )
    ]


def _parameter_compatible_with_phase_constituents(
    parameter: Parameter,
    phase: Phase | None,
) -> bool:
    """Return whether one parameter target is legal for the exported phase."""
    if phase is None or not phase.constituents or len(parameter.target) != 1:
        return True
    target = _tdb_parameter_target_text(parameter.target[0])
    sections = target.split(":")
    allowed = _phase_sublattice_species(phase)
    if len(sections) != len(allowed):
        return True

    for section, allowed_species in zip(sections, allowed, strict=True):
        allowed_ids = {_tdb_identifier(species) for species in allowed_species}
        tokens = [
            _tdb_constituent_name(token)
            for token in section.split(",")
            if token.strip()
        ]
        if "*" in tokens:
            continue
        for token in tokens:
            if _tdb_identifier(token) not in allowed_ids:
                return False
    return True


def _parameters_without_symmetry_duplicate_zero_endmembers(
    parameters: list[Parameter],
    phases: list[Phase],
    phase_alias_map: dict[str, str],
) -> list[Parameter]:
    """Drop zero G endmember permutations equivalent under ``:F``."""
    phases_by_name = {
        _tdb_identifier(phase.name): phase
        for phase in phases
        if _tdb_identifier(phase.name) not in phase_alias_map
    }
    authoritative_keys = {
        key
        for parameter in parameters
        if not _is_zero_endmember_parameter(parameter)
        for key in [
            _symmetry_reduced_endmember_parameter_key(
                parameter,
                phases_by_name.get(_tdb_identifier(parameter.phase_name)),
            )
        ]
        if key is not None
    }
    emitted_generated_keys: set[
        tuple[str, str, tuple[str, ...], int]
    ] = set()
    output: list[Parameter] = []
    for parameter in parameters:
        if not _is_zero_endmember_parameter(parameter):
            output.append(parameter)
            continue
        key = _symmetry_reduced_endmember_parameter_key(
            parameter,
            phases_by_name.get(_tdb_identifier(parameter.phase_name)),
        )
        if key is None:
            output.append(parameter)
            continue
        if key in authoritative_keys or key in emitted_generated_keys:
            continue
        emitted_generated_keys.add(key)
        output.append(parameter)
    return output


def _symmetry_reduced_endmember_parameter_key(
    parameter: Parameter,
    phase: Phase | None,
) -> tuple[str, str, tuple[str, ...], int] | None:
    """Return a parameter identity reduced by equivalent sublattice groups."""
    if phase is None or len(parameter.target) != 1:
        return None
    if (
        _tdb_parameter_type(parameter) != "G"
        or int(parameter.order) != 0
        or _is_interaction_parameter(parameter)
        or _parameter_target_has_wildcard(parameter)
    ):
        return None
    groups = _symmetry_equivalent_sublattice_groups(phase)
    if not groups:
        return None
    target = _tdb_parameter_target_text(parameter.target[0])
    if len(target.split(":")) != len(phase.constituents):
        return None
    return (
        _tdb_identifier(parameter.phase_name),
        "G",
        _symmetry_reduced_target_key(target, groups),
        0,
    )


def _is_zero_endmember_parameter(parameter: Parameter) -> bool:
    """Return whether a parameter is a zero G endmember record."""
    if (
        _tdb_parameter_type(parameter) != "G"
        or int(parameter.order) != 0
        or len(parameter.target) != 1
        or _is_interaction_parameter(parameter)
        or _parameter_target_has_wildcard(parameter)
    ):
        return False
    return _parameter_expression_is_zero(parameter)


def _parameter_expression_is_zero(parameter: Parameter) -> bool:
    """Return whether every Gibbs range for a parameter is zero-valued."""
    expression = _normalize_tdb_expression(
        _parameter_source_expression_for_export(parameter) or parameter.expression
    )
    ranges = _parse_function_ranges(expression)
    if not ranges:
        return _gibbs_segment_is_zero(expression or "0")
    return all(_gibbs_segment_is_zero(gibbs_range.Gibbs) for gibbs_range in ranges)


def _gibbs_segment_is_zero(expression: str) -> bool:
    try:
        canonical = _canonical_tdb_gibbs_expression(expression)
    except ValueError:
        return False
    return canonical in {"0", "ZERO#"}


def _expanded_wildcard_targets(parameter: Parameter, phase: Phase) -> list[str]:
    """Return concrete targets for one parameter target containing ``*``."""
    if len(parameter.target) != 1:
        raise ValueError(
            "Cannot export wildcard PAR target split across multiple targets: "
            f"{parameter.phase_name!r}"
        )
    allowed = _phase_sublattice_species(phase)
    sections = _tdb_parameter_target_text(parameter.target[0]).split(":")
    if len(sections) != len(allowed):
        raise ValueError(
            "Cannot export wildcard PAR target with mismatched sublattices: "
            f"{parameter.phase_name!r} {parameter.target[0]!r}"
        )

    options: list[list[str]] = []
    for section, sublattice_species in zip(sections, allowed, strict=True):
        tokens = [
            _tdb_constituent_name(token)
            for token in section.split(",")
            if token.strip()
        ]
        if "*" not in tokens:
            options.append([",".join(tokens)])
            continue
        replacements: list[str] = []
        for species in sublattice_species:
            replaced = [species if token == "*" else token for token in tokens]
            if _has_duplicate_identifiers(replaced):
                continue
            replacements.append(",".join(replaced))
        options.append(replacements)
    targets = _unique_preserve_order(
        [":".join(combination) for combination in product(*options)]
    )
    return _unique_symmetry_reduced_targets(targets, phase)


def _has_duplicate_identifiers(values: list[str]) -> bool:
    identifiers = [_tdb_identifier(value) for value in values]
    return len(set(identifiers)) != len(identifiers)


def _unique_symmetry_reduced_targets(targets: list[str], phase: Phase) -> list[str]:
    """Return one wildcard-expanded target per phase symmetry class."""
    groups = _symmetry_equivalent_sublattice_groups(phase)
    if not groups:
        return targets

    output: list[str] = []
    seen: set[tuple[str, ...]] = set()
    for target in targets:
        key = _symmetry_reduced_target_key(target, groups)
        if key in seen:
            continue
        seen.add(key)
        output.append(target)
    return output


def _symmetry_equivalent_sublattice_groups(phase: Phase) -> list[list[int]]:
    """Return zero-based sublattice groups equivalent under the ``:F`` option."""
    if not _phase_has_f_option(phase):
        return []
    groups_by_key: dict[tuple[float, tuple[str, ...]], list[int]] = {}
    for index, constituent in enumerate(phase.constituents):
        species_key = tuple(_tdb_identifier(species) for species in constituent.species)
        if not species_key:
            continue
        key = (round(float(constituent.site_ratio), 12), species_key)
        groups_by_key.setdefault(key, []).append(index)
    return [indices for indices in groups_by_key.values() if len(indices) > 1]


def _phase_has_f_option(phase: Phase) -> bool:
    source_name = _source_phase_name(phase.source.command)
    if ":" not in source_name:
        return False
    suffix = source_name.split(":", 1)[1]
    return "F" in suffix.upper()


def _symmetry_reduced_target_key(
    target: str,
    groups: list[list[int]],
) -> tuple[str, ...]:
    sections = target.split(":")
    keyed_sections = [_canonical_target_section_key(section) for section in sections]
    for group in groups:
        if any(index >= len(keyed_sections) for index in group):
            continue
        reduced_group = sorted(keyed_sections[index] for index in group)
        for index, value in zip(group, reduced_group, strict=True):
            keyed_sections[index] = value
    return tuple(keyed_sections)


def _canonical_target_section_key(section: str) -> str:
    tokens = [
        _tdb_identifier(token)
        for token in section.split(",")
        if token.strip()
    ]
    return ",".join(sorted(tokens))


def _unique_preserve_order(values: list[str]) -> list[str]:
    seen: set[str] = set()
    unique: list[str] = []
    for value in values:
        key = _tdb_identifier(value)
        if key in seen:
            continue
        seen.add(key)
        unique.append(value)
    return unique


def _deduplicate_parameter_exports(parameters: list[Parameter]) -> list[Parameter]:
    """Merge duplicate exported PAR identities created by wildcard expansion."""
    output: list[Parameter] = []
    by_key: dict[tuple[str, str, tuple[str, ...], int], int] = {}
    for parameter in parameters:
        key = _parameter_output_key(parameter)
        existing_index = by_key.get(key)
        if existing_index is None:
            by_key[key] = len(output)
            output.append(parameter)
            continue
        output[existing_index] = _merge_duplicate_parameter_export(
            output[existing_index],
            parameter,
        )
    return output


def _merge_duplicate_parameter_export(left: Parameter, right: Parameter) -> Parameter:
    """Return one PAR record for additive duplicate identities."""
    left_rank = _parameter_specificity_rank(left)
    right_rank = _parameter_specificity_rank(right)
    primary, secondary = (left, right) if left_rank <= right_rank else (right, left)
    if _tdb_identifier(left.parameter_type) != _tdb_identifier(right.parameter_type):
        return primary
    parameter_type = _tdb_parameter_type(primary)
    left_expression = _canonical_tdb_parameter_expression(
        left,
        _tdb_parameter_type(left),
    )
    right_expression = _canonical_tdb_parameter_expression(
        right,
        _tdb_parameter_type(right),
    )
    if left_expression == right_expression:
        return primary
    expression = _additive_parameter_expression(
        left_expression,
        right_expression,
        zero_as_function=parameter_type in {"G", "L"},
    )
    metadata = dict(primary.metadata)
    metadata["merged_duplicate_parameter"] = "true"
    metadata.pop("source_expression", None)
    if left.metadata.get("expanded_from_wildcard") or right.metadata.get(
        "expanded_from_wildcard"
    ):
        metadata["merged_expanded_wildcard"] = "true"
    return replace(primary, expression=expression, metadata=metadata)


def _additive_parameter_expression(
    left: str,
    right: str,
    *,
    zero_as_function: bool,
) -> str:
    """Return a TDB range expression that adds two same-identity parameters."""
    left_ranges = _parameter_expression_ranges(left)
    right_ranges = _parameter_expression_ranges(right)
    if len(left_ranges) != len(right_ranges):
        raise ValueError(
            "Cannot merge duplicate PAR records with different temperature ranges."
        )
    merged: list[GibbsRange] = []
    for left_range, right_range in zip(left_ranges, right_ranges, strict=True):
        if (
            abs(float(left_range.T_min) - float(right_range.T_min)) > 1e-8
            or abs(float(left_range.T_max) - float(right_range.T_max)) > 1e-8
            or left_range.status.upper()[:1] != right_range.status.upper()[:1]
        ):
            raise ValueError(
                "Cannot merge duplicate PAR records with different temperature ranges."
            )
        merged.append(
            GibbsRange(
                float(left_range.T_min),
                float(left_range.T_max),
                _additive_gibbs_expression(
                    left_range.Gibbs,
                    right_range.Gibbs,
                ),
                left_range.status,
            )
        )
    return _canonical_tdb_range_expression(
        merged,
        zero_as_function=zero_as_function,
    )


def _parameter_expression_ranges(expression: str) -> list[GibbsRange]:
    ranges = _parse_function_ranges(expression)
    if ranges:
        return ranges
    return [GibbsRange(298.15, 6000.0, expression or "0", "N")]


def _additive_gibbs_expression(left: str, right: str) -> str:
    """Return two Gibbs segment bodies joined as a readable sum."""
    pieces: list[str] = []
    for expression in (left, right):
        for term in _signed_expression_terms(expression):
            normalized = _normalize_reference_signed_term(term)
            if pieces and not normalized.startswith(("+", "-")):
                normalized = f"+{normalized}"
            pieces.append(normalized)
    result = "".join(pieces)
    return result[1:] if result.startswith("+") else result or "0"


def _parameter_output_key(
    parameter: Parameter,
) -> tuple[str, str, tuple[str, ...], int]:
    return (
        _tdb_identifier(parameter.phase_name),
        _tdb_parameter_type(parameter),
        tuple(
            _tdb_identifier(_tdb_parameter_target_text(target))
            for target in parameter.target
        ),
        int(parameter.order),
    )


def _parameter_specificity_rank(parameter: Parameter) -> tuple[int, int, int]:
    """Return smaller rank for more authoritative export records."""
    return (
        1 if parameter.metadata.get("expanded_from_wildcard") else 0,
        1 if parameter.metadata.get("generated_by") else 0,
        {"BMAGN": 0, "BMAG": 1, "BM": 2}.get(
            _tdb_identifier(parameter.parameter_type),
            0,
        ),
    )


def _phase_sublattice_species(phase: Phase) -> list[list[str]]:
    """Return normalized per-sublattice species in export order."""
    species_by_sublattice: list[list[str]] = []
    for constituent in sorted(phase.constituents, key=lambda item: item.sublattice):
        species = [
            _tdb_constituent_name(name)
            for name in constituent.species
            if str(name).strip() and str(name).strip() != "*"
        ]
        if not species:
            return []
        species_by_sublattice.append(species)
    return species_by_sublattice


def _parameters_with_missing_all_vacancy_endmembers(
    parameters: list[Parameter],
    phases: list[Phase],
    phase_alias_map: dict[str, str],
) -> list[Parameter]:
    """Add ZERO# records only for missing all-vacancy endmembers."""
    phases_by_name = {
        _tdb_identifier(phase.name): phase
        for phase in phases
        if _tdb_identifier(phase.name) not in phase_alias_map
    }
    existing_endmembers = {
        (_tdb_identifier(parameter.phase_name), _canonical_parameter_target(parameter))
        for parameter in parameters
        if _tdb_parameter_type(parameter) == "G"
        and not _is_interaction_parameter(parameter)
        and not _parameter_target_has_wildcard(parameter)
        and parameter.target
    }
    added: list[Parameter] = []
    for phase_name, phase in phases_by_name.items():
        if _phase_uses_symmetry_reduced_source(phase):
            continue
        target = _phase_all_vacancy_endmember_target(phase)
        if not target:
            continue
        key = (phase_name, target)
        if key in existing_endmembers:
            continue
        added.append(
            Parameter(
                phase_name=phase.name,
                parameter_type="G",
                target=[target],
                order=0,
                expression="298.15 0; 6000 N",
                metadata={
                    "generated_by": "tdb_writer",
                    "reason": "missing_all_vacancy_endmember",
                },
            )
        )
        existing_endmembers.add(key)
    return [*parameters, *added]


def _phase_all_vacancy_endmember_target(phase: Phase) -> str:
    """Return the all-vacancy endmember target when every sublattice allows it."""
    target: list[str] = []
    for constituent in sorted(phase.constituents, key=lambda item: item.sublattice):
        vacancy = next(
            (
                _tdb_constituent_name(name)
                for name in constituent.species
                if _tdb_identifier(name) == "VA"
            ),
            "",
        )
        if not vacancy:
            return ""
        target.append(vacancy)
    return ":".join(target)


def _parameter_target_has_wildcard(parameter: Parameter) -> bool:
    """Return whether a parameter target contains a TDB wildcard token."""
    return any(
        token.strip() == "*"
        for target in parameter.target
        for token in re.split(r"[:,]", str(target))
    )


def _canonical_parameter_target(parameter: Parameter) -> str:
    """Return the canonical single-target representation for endmember keys."""
    return ":".join(
        _tdb_constituent_name(token)
        for token in str(parameter.target[0]).split(":")
        if token.strip()
    )


def _phase_concrete_endmember_targets(phase: Phase) -> list[str]:
    """Return concrete endmember targets implied by one phase constitution."""
    if not phase.constituents:
        return []
    sublattice_species: list[list[str]] = []
    for constituent in sorted(phase.constituents, key=lambda item: item.sublattice):
        species = [
            _tdb_constituent_name(name)
            for name in constituent.species
            if str(name).strip() and str(name).strip() != "*"
        ]
        if not species:
            return []
        sublattice_species.append(species)

    targets: list[str] = []
    for combination in product(*sublattice_species):
        targets.append(":".join(combination))
    return targets


def _phase_uses_symmetry_reduced_source(phase: Phase) -> bool:
    """Return whether the source PHASE name carries a symmetry/model suffix."""
    return ":" in _tdb_phase_name(phase)


def _parameter_groups(
    database: DatabaseIR,
    parameters: list[Parameter],
) -> dict[int, dict[tuple[str, ...], list[Parameter]]]:
    species_elements = {
        species.name.upper(): tuple(sorted(species.composition))
        for species in database.species
        if species.composition
    }
    groups: dict[int, dict[tuple[str, ...], list[Parameter]]] = {}
    for parameter in parameters:
        elements = _parameter_elements(parameter, species_elements)
        arity = max(1, len(elements))
        groups.setdefault(arity, {}).setdefault(elements, []).append(parameter)
    return groups


def _parameters_with_export_phase_aliases(
    parameters: list[Parameter],
    phase_alias_map: dict[str, str],
) -> list[Parameter]:
    """Return parameters with duplicate disordered-alias phases suppressed."""
    if not phase_alias_map:
        return list(parameters)

    canonical_keys = {
        _parameter_duplicate_key(parameter)
        for parameter in parameters
        if _tdb_identifier(parameter.phase_name) not in phase_alias_map
    }
    exported: list[Parameter] = []
    emitted_keys: set[tuple[str, str, tuple[str, ...], int]] = set()
    for parameter in parameters:
        phase_name = _tdb_identifier(parameter.phase_name)
        canonical_phase = phase_alias_map.get(phase_name)
        export_parameter = (
            replace(parameter, phase_name=canonical_phase)
            if canonical_phase is not None
            else parameter
        )
        key = _parameter_duplicate_key(export_parameter)
        if canonical_phase is not None and key in canonical_keys:
            continue
        if key in emitted_keys:
            continue
        exported.append(export_parameter)
        emitted_keys.add(key)
    return exported


def _parameter_duplicate_key(
    parameter: Parameter,
) -> tuple[str, str, tuple[str, ...], int]:
    """Return the TDB identity of a parameter independent of expression text."""
    return (
        _tdb_identifier(parameter.phase_name),
        _tdb_identifier(parameter.parameter_type),
        tuple(_tdb_identifier(target) for target in parameter.target),
        int(parameter.order),
    )


def _parameter_elements(
    parameter: Parameter,
    species_elements: dict[str, tuple[str, ...]],
) -> tuple[str, ...]:
    names: set[str] = set()
    saw_vacancy = False
    for target in parameter.target:
        for token in re.split(r"[:,;\s]+", target):
            token = token.strip()
            if not token or token == "*":
                continue
            upper = token.upper()
            if upper == "VA":
                saw_vacancy = True
                continue
            if upper in species_elements:
                names.update(species_elements[upper])
            else:
                names.update(_element_symbols_from_name(token))
    if not names and saw_vacancy:
        return ("VA",)
    return tuple(sorted(names))


def _element_symbols_from_name(name: str) -> list[str]:
    """Best-effort element symbols from a target/endmember name."""
    text = name.replace("_", "")
    symbols = re.findall(r"[A-Z][a-z]?|[A-Z]{1,2}", text)
    if not symbols and text.isalpha() and len(text) <= 2:
        symbols = [text]
    return [_conventional_element_symbol(symbol) for symbol in symbols]


def format_tdb_parameter(
    parameter: Parameter,
    database: DatabaseIR | None = None,
) -> list[str]:
    """Format one PAR command with stable wrapping."""
    target = _tdb_parameter_target(parameter)
    parameter_type = _tdb_parameter_type(parameter)
    phase_name = _tdb_identifier(parameter.phase_name)
    order = int(parameter.order)
    expression = _canonical_tdb_parameter_expression(
        parameter,
        parameter_type,
        database,
    )
    reference = _tdb_parameter_reference(parameter)
    command = (
        f"PAR {parameter_type}({phase_name},{target};{order}) "
        f"{expression}{reference} !"
    )
    return _wrap_tdb_command(" ", command)


def _tdb_parameter_reference(parameter: Parameter) -> str:
    """Return the TDB reference suffix for a parameter command."""
    reference = _tdb_reference_key(str(parameter.metadata.get("reference", "")))
    return f" {reference}" if reference else ""


def _tdb_parameter_target(parameter: Parameter) -> str:
    """Return a PAR target with conventional element/species casing."""
    return ",".join(
        _tdb_parameter_target_text(target)
        for target in parameter.target
        if str(target).strip()
    )


def _tdb_parameter_target_text(target: str) -> str:
    """Normalize one colon/comma-delimited PAR target component."""
    sections: list[str] = []
    for section in str(target).split(":"):
        tokens = [
            _tdb_constituent_name(token)
            for token in section.split(",")
            if token.strip()
        ]
        sections.append(",".join(tokens) if tokens else "")
    return ":".join(sections)


def _tdb_parameter_type(parameter: Parameter) -> str:
    """Return the exported parameter type.

    Preserve the parsed TDB parameter family.  Some vendor TDBs use comma or
    wildcard targets in ``G`` records for ordered-phase terms, and rewriting
    them as ``L`` changes strict endmember coverage on roundtrip.  For
    programmatically-created IR records without source text, keep the legacy
    convenience rule that mixed ``G`` targets export as ``L``.  Magnetic moment
    aliases export as ``BMAGN``; other parameter families, such as ``TC``, keep
    their original type.
    """
    parameter_type = _tdb_identifier(parameter.parameter_type)
    if parameter_type in {"BM", "BMAG"}:
        return "BMAGN"
    source_type = _source_tdb_parameter_type(parameter)
    if source_type in {"G", "L"}:
        return parameter_type
    if parameter_type == "G":
        return "L" if _is_interaction_parameter(parameter) else "G"
    return parameter_type


def _source_tdb_parameter_type(parameter: Parameter) -> str:
    """Return the parameter family declared by a source TDB PAR command."""
    match = re.match(
        r"^\s*(?:PARAMETER|PARAM|PAR)\s+([A-Za-z]+)\s*\(",
        parameter.source.command,
        flags=re.IGNORECASE,
    )
    if not match:
        return ""
    return _tdb_identifier(match.group(1))


def _is_interaction_parameter(parameter: Parameter) -> bool:
    """Return whether a parameter target describes a mixed sublattice."""
    return len(parameter.target) > 1 or any(
        "," in str(target) for target in parameter.target
    )


def _canonical_tdb_parameter_expression(
    parameter: Parameter,
    parameter_type: str,
    database: DatabaseIR | None = None,
) -> str:
    """Return a strict TDB expression body for one PAR command."""
    source_expression = _parameter_source_expression_for_export(parameter)
    expression = _normalize_tdb_expression(
        source_expression or parameter.expression
    )
    reconstructed = _compound_reference_expression_for_export(
        parameter,
        parameter_type,
        expression,
        database,
    )
    if reconstructed:
        expression = reconstructed

    ranges = _parse_function_ranges(expression)
    zero_as_function = parameter_type in {"G", "L"}
    if ranges:
        return _canonical_tdb_range_expression(
            ranges,
            zero_as_function=zero_as_function,
        )

    segment_expression = expression or "0"
    segment_expression = _canonical_tdb_gibbs_expression(segment_expression)
    if zero_as_function:
        if segment_expression == "0":
            segment_expression = "ZERO#"
        return f"298.15 {segment_expression}; 6000 N"
    return segment_expression


def _parameter_source_expression_for_export(parameter: Parameter) -> str:
    """Return original composite PAR expression when it is still authoritative."""
    return str(parameter.metadata.get("source_expression", "")).strip()


def _compound_reference_expression_for_export(
    parameter: Parameter,
    parameter_type: str,
    expression: str,
    database: DatabaseIR | None,
) -> str:
    """Recover compact compound ``formation + GHSER#`` expressions when proven."""
    if database is None or parameter_type != "G":
        return ""
    phase = _phase_by_export_name(database).get(_tdb_identifier(parameter.phase_name))
    if phase is None or str(phase.model).upper() != "COMPOUND":
        return ""
    reference_terms = _compound_reference_terms(parameter, phase, database)
    if not reference_terms:
        return ""
    ranges = _parse_function_ranges(expression)
    if not ranges:
        return ""
    if any(_expression_has_symbolic_references(row.Gibbs) for row in ranges):
        return ""

    reference_ranges = {
        name: _function_range_records(function)
        for name, _coefficient, function in reference_terms
    }
    output_ranges: list[GibbsRange] = []
    for gibbs_range in ranges:
        coefficients = _gibbs_coefficients_from_segment(gibbs_range.Gibbs)
        if coefficients is None:
            return ""
        bounds = _compound_reference_segment_bounds(
            gibbs_range,
            reference_ranges,
        )
        if len(bounds) < 2:
            return ""
        for lower, upper in zip(bounds, bounds[1:], strict=False):
            residual = list(coefficients)
            for name, coefficient, _function in reference_terms:
                reference_range = _active_function_range(
                    reference_ranges[name],
                    lower,
                    upper,
                )
                if reference_range is None:
                    return ""
                reference_coefficients = _gibbs_coefficients_from_segment(
                    reference_range.Gibbs,
                )
                if reference_coefficients is None:
                    return ""
                residual = _subtract_gibbs_coefficients(
                    residual,
                    reference_coefficients,
                    coefficient,
                )
                if residual is None:
                    return ""
            body_terms = _signed_expression_terms(
                _gibbs_row_expression(_zero_small_coefficients(residual))
            )
            body_terms.extend(
                _format_reference_term(name, coefficient)
                for name, coefficient, _function in reference_terms
            )
            body = _join_signed_terms(body_terms)
            if not body:
                return ""
            output_ranges.append(GibbsRange(lower, upper, body, "Y"))
    if not output_ranges:
        return ""
    for index, output_range in enumerate(output_ranges):
        status = ranges[-1].status if index == len(output_ranges) - 1 else "Y"
        output_ranges[index] = replace(output_range, status=status)
    return _canonical_tdb_range_expression(output_ranges)


def _phase_by_export_name(database: DatabaseIR) -> dict[str, Phase]:
    """Return phases keyed by exported TDB identifier."""
    return {_tdb_identifier(phase.name): phase for phase in database.phases}


def _expression_has_symbolic_references(expression: str) -> bool:
    """Return whether an expression already contains non-math identifiers."""
    pattern = r"(?<![A-Za-z0-9_])([A-Za-z_][A-Za-z0-9_]*#?)"
    for match in re.finditer(pattern, expression):
        name = match.group(1).rstrip("#").upper()
        if name not in {"T", "LN", "LOG", "EXP", "R"}:
            return True
    return False


def _compound_reference_terms(
    parameter: Parameter,
    phase: Phase,
    database: DatabaseIR,
) -> list[tuple[str, float, FunctionDefinition]]:
    """Return stoichiometric reference functions for a compound endmember."""
    if not parameter.target:
        return []
    sections = str(parameter.target[0]).split(":")
    constituents = sorted(phase.constituents, key=lambda item: item.sublattice)
    if len(sections) != len(constituents):
        return []
    functions_by_name = {
        _tdb_identifier(function.name): function for function in database.functions
    }
    terms: list[tuple[str, float, FunctionDefinition]] = []
    for section, constituent in zip(sections, constituents, strict=False):
        species_names = [name.strip() for name in section.split(",") if name.strip()]
        if len(species_names) != 1:
            return []
        species_name = species_names[0]
        if species_name.upper() in {"VA", "/-", "*"}:
            continue
        element_symbol = _conventional_element_symbol(species_name)
        if not element_symbol.isalpha() or len(element_symbol) > 2:
            return []
        reference_name = f"GHSER{element_symbol.upper()}"
        reference_function = functions_by_name.get(reference_name)
        if reference_function is None:
            return []
        terms.append(
            (reference_name, float(constituent.site_ratio), reference_function)
        )
    return terms


def _active_function_range(
    ranges: list[GibbsRange],
    lower: float,
    upper: float,
) -> GibbsRange | None:
    """Return the reference function range covering or extrapolating a segment."""
    midpoint = (float(lower) + float(upper)) / 2.0
    for gibbs_range in ranges:
        if (
            float(gibbs_range.T_min) - 1e-8
            <= midpoint
            <= float(gibbs_range.T_max) + 1e-8
        ):
            return gibbs_range
    if not ranges:
        return None
    if midpoint < float(ranges[0].T_min):
        return ranges[0]
    return ranges[-1]


def _compound_reference_segment_bounds(
    gibbs_range: GibbsRange,
    reference_ranges: dict[str, list[GibbsRange]],
) -> list[float]:
    """Return source bounds split at active reference-function boundaries."""
    lower = float(gibbs_range.T_min)
    upper = float(gibbs_range.T_max)
    bounds = [lower, upper]
    for ranges in reference_ranges.values():
        for reference_range in ranges:
            for bound in (float(reference_range.T_min), float(reference_range.T_max)):
                if lower + 1e-8 < bound < upper - 1e-8:
                    bounds.append(bound)
    return _merged_temperature_bounds(bounds)


def _merged_temperature_bounds(bounds: list[float]) -> list[float]:
    merged: list[float] = []
    for bound in sorted(float(value) for value in bounds):
        if not merged or abs(bound - merged[-1]) > 1e-8:
            merged.append(bound)
    return merged


def _subtract_gibbs_coefficients(
    left: list[float],
    right: list[float],
    scale: float,
) -> list[float] | None:
    """Subtract Gibbs coefficient rows while matching custom powers by basis."""
    terms = _coefficient_terms(left)
    terms.extend(
        (power, -scale * coefficient)
        for power, coefficient in _coefficient_terms(right)
    )
    return _gibbs_coefficients_from_parsed_terms(terms)


def _coefficient_terms(
    coefficients: list[float],
) -> list[tuple[float | str, float]]:
    """Return parsed Gibbs terms from a compact coefficient row."""
    terms: list[tuple[float | str, float]] = [
        (0.0, coefficients[0]),
        (1.0, coefficients[1]),
        ("TLOGT", coefficients[2]),
        (2.0, coefficients[3]),
        (3.0, coefficients[4]),
        (-1.0, coefficients[5]),
    ]
    for index in range(4):
        coefficient = coefficients[6 + 2 * index]
        power = coefficients[7 + 2 * index]
        if coefficient != 0.0:
            terms.append((power, coefficient))
    return terms


def _zero_small_coefficients(coefficients: list[float]) -> list[float]:
    """Return coefficients without dropping small nonzero physical terms."""
    return [_clean_residual_coefficient(value) for value in coefficients]


def _clean_residual_coefficient(value: float) -> float:
    """Return a readable residual coefficient after cancellation arithmetic."""
    if value == 0.0:
        return 0.0
    return value


def _join_signed_terms(terms: list[str]) -> str:
    """Join additive terms and remove a leading plus sign."""
    result = "".join(
        _normalize_reference_signed_term(term)
        for term in terms
        if str(term).strip() and str(term).strip() != "0"
    )
    return result[1:] if result.startswith("+") else result


def _canonical_tdb_range_expression(
    ranges: list[GibbsRange],
    *,
    zero_as_function: bool = False,
) -> str:
    """Return canonical piecewise expression text for FUNCTION/PAR records."""
    lines: list[str] = []
    for index, gibbs_range in enumerate(ranges):
        expression = _canonical_tdb_gibbs_expression(gibbs_range.Gibbs)
        if zero_as_function and expression == "0":
            expression = "ZERO#"
        prefix = f"{_tdb_number(gibbs_range.T_min)} " if index == 0 else ""
        lines.append(
            f"{prefix}{expression}; "
            f"{_tdb_number(gibbs_range.T_max)} {gibbs_range.status.upper()[:1] or 'N'}"
        )
    return " ".join(lines)


def _arity_label(arity: int) -> str:
    if arity == 1:
        return "Unary"
    if arity == 2:
        return "Binary"
    if arity == 3:
        return "Ternary"
    if arity == 4:
        return "Quaternary"
    return f"{arity}-component"


def _comment_arrow(character: str) -> str:
    return "$ " + character * _COMMENT_BODY_WIDTH


def _comment_separator() -> str:
    return "$ " + "-" * _COMMENT_BODY_WIDTH


def _comment_line(text: str) -> str:
    """Return one TDB comment line capped at the TDB line width."""
    normalized = " ".join(str(text).split())
    return "$ " + normalized[:_COMMENT_BODY_WIDTH]


def _comment_preformatted_line(text: str) -> str:
    """Return one TDB comment line without collapsing table spacing."""
    return "$ " + str(text).rstrip()[:_COMMENT_BODY_WIDTH]


def _canonical_function_range_lines(function: FunctionDefinition) -> list[str]:
    """Return canonical TDB range lines from parsed Gibbs ranges."""
    ranges = _function_range_records(function)
    if not ranges:
        ranges = [GibbsRange(298.15, 6000.0, "0", "N")]

    lines: list[str] = []
    for index, gibbs_range in enumerate(ranges):
        expression = _canonical_tdb_gibbs_expression(gibbs_range.Gibbs)
        prefix = f"{_tdb_number(gibbs_range.T_min)} " if index == 0 else ""
        lines.append(
            f"{prefix}{expression}; "
            f"{_tdb_number(gibbs_range.T_max)} {gibbs_range.status.upper()[:1] or 'N'}"
        )
    return lines


def _function_range_records(function: FunctionDefinition) -> list[GibbsRange]:
    """Return function Gibbs ranges, falling back to parsing source expression."""
    source_expression = _function_source_expression_for_export(function)
    if source_expression:
        source_ranges = _parse_function_ranges(
            _normalize_tdb_expression(source_expression)
        )
        if source_ranges:
            return source_ranges
    if function.gibbs_ranges:
        return function.gibbs_ranges
    expression = _normalize_tdb_expression(function.expression)
    if not expression:
        return []
    return _parse_function_ranges(expression)


def _function_source_expression_for_export(function: FunctionDefinition) -> str:
    """Return original composite TDB expression when it still matches the function."""
    source_expression = function.tdb_src.source_expression.strip()
    if not source_expression:
        return ""
    current_expression = function.expression.strip()
    if current_expression and not _tdb_expressions_equivalent(
        current_expression,
        source_expression,
    ):
        return ""
    return source_expression


def _tdb_expressions_equivalent(left: str, right: str) -> bool:
    """Return whether two TDB range expressions canonicalize to the same text."""
    left_ranges = _parse_function_ranges(_normalize_tdb_expression(left))
    right_ranges = _parse_function_ranges(_normalize_tdb_expression(right))
    if left_ranges or right_ranges:
        if not left_ranges or not right_ranges:
            return False
        left_expression = _canonical_tdb_range_expression(left_ranges)
        right_expression = _canonical_tdb_range_expression(right_ranges)
        return left_expression == right_expression
    return _canonical_tdb_gibbs_expression(left) == _canonical_tdb_gibbs_expression(
        right
    )


def _parse_function_ranges(expression: str) -> list[GibbsRange]:
    """Parse FUNCTION piecewise ranges for writer fallback paths."""
    lower_match = re.match(
        r"\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EDed][+-]?\d+)?)",
        expression,
    )
    if lower_match is None:
        return []
    lower = float(lower_match.group(1).replace("D", "E").replace("d", "E"))
    remainder = expression[lower_match.end() :].lstrip()
    ranges: list[GibbsRange] = []
    while ";" in remainder:
        segment_expression, remainder = remainder.split(";", 1)
        remainder = remainder.lstrip()
        upper_match = re.match(
            r"\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EDed][+-]?\d+)?)",
            remainder,
        )
        if upper_match is None:
            break
        upper = float(upper_match.group(1).replace("D", "E").replace("d", "E"))
        tail = remainder[upper_match.end() :].lstrip()
        status = tail[:1].upper() if tail[:1].upper() in {"Y", "N"} else "N"
        if tail[:1].upper() in {"Y", "N"}:
            tail = tail[1:].lstrip()
        ranges.append(
            GibbsRange(
                T_min=lower,
                T_max=upper,
                Gibbs=segment_expression.strip(),
                status=status,
            )
        )
        lower = upper
        remainder = tail
    return ranges


def _function_prefix(name: str) -> str:
    """Return ``FUNCTION`` prefix with a fixed-width name field and 4 spaces."""
    return f" FUNCTION {name:<{_FUNCTION_NAME_WIDTH}}    "


def _ordered_functions(
    functions: list[FunctionDefinition],
) -> list[FunctionDefinition]:
    """Return root functions before functions that reference them."""
    by_name = {_tdb_identifier(function.name): function for function in functions}
    original_index = {
        _tdb_identifier(function.name): index
        for index, function in enumerate(functions)
    }
    dependencies = {
        name: _function_dependencies(function.expression, set(by_name)) - {name}
        for name, function in by_name.items()
    }

    ordered: list[FunctionDefinition] = []
    emitted: set[str] = set()
    while len(emitted) < len(by_name):
        ready = [
            name
            for name, required in dependencies.items()
            if name not in emitted and required <= emitted
        ]
        if not ready:
            ready = [
                name
                for name in by_name
                if name not in emitted
            ]
        ready.sort(key=lambda name: original_index[name])
        for name in ready:
            ordered.append(by_name[name])
            emitted.add(name)
    return ordered


def _function_dependencies(expression: str, function_names: set[str]) -> set[str]:
    expression_upper = expression.upper()
    dependencies: set[str] = set()
    for name in function_names:
        pattern = rf"(?<![A-Z0-9_]){re.escape(name)}(?![A-Z0-9_])"
        if re.search(pattern, expression_upper):
            dependencies.add(name)
    return dependencies


def _split_function_ranges(expression: str) -> list[str]:
    parts = [part.strip() for part in expression.split(";")]
    if len(parts) < 2:
        return []

    range_lines: list[str] = []
    current_expression = parts[0]
    for index, status_part in enumerate(parts[1:], start=1):
        if not status_part:
            continue
        match = _FLOAT_STATUS_RE.match(status_part)
        if match is None:
            return []
        upper, status, tail = match.groups()
        status = status.upper()
        tail = (tail or "").strip()
        is_last = index == len(parts) - 1
        suffix = f" {tail}" if is_last and tail else ""
        range_lines.append(
            f"{current_expression}; {_tdb_number_text(upper)} {status}{suffix}"
        )
        current_expression = "" if is_last else tail

    return [line for line in range_lines if line.strip()]


def _normalize_tdb_expression(expression: str) -> str:
    text = " ".join(str(expression).replace("\n", " ").split())
    if text.endswith("!"):
        text = text[:-1].rstrip()
    return text


def _normalize_bare_decimal_exponents(expression: str) -> str:
    """Return expression with ``1.E-4`` style literals as ``1.0E-4``."""
    return _BARE_DECIMAL_EXPONENT_RE.sub(r"\1.0", str(expression))


def _normalize_numeric_exponents(expression: str) -> str:
    """Convert TDB ``D`` exponents without touching identifiers.

    A plain string replacement would corrupt valid function names such as
    ``LDF0MNNI`` into ``LEF0MNNI``.  Only exponent markers adjacent to numeric
    text are normalized.
    """
    return re.sub(r"(?<=[0-9.])[Dd](?=[+\-]?\d)", "E", str(expression))


def _canonical_tdb_gibbs_expression(expression: str) -> str:
    """Canonicalize one Gibbs expression segment for strict TDB export."""
    expression = _normalize_bare_decimal_exponents(
        _simplify_symbolic_functions(expression)
    )
    numeric_parsed_terms: list[tuple[float | str, float]] = []
    referenced_terms: dict[str, float] = {}
    referenced_order: list[str] = []
    for term, scale in _expanded_scaled_terms(expression):
        if _contains_pressure_variable(term):
            raise ValueError(f"Unsupported TDB Gibbs expression: {expression}")
        reference = _reference_term(term)
        if reference is None:
            parsed_terms = _parsed_gibbs_terms(term, scale)
            if parsed_terms is None:
                raise ValueError(f"Unsupported TDB Gibbs expression: {expression}")
            numeric_parsed_terms.extend(parsed_terms)
            continue
        name, coefficient = reference
        if name not in referenced_terms:
            referenced_terms[name] = 0.0
            referenced_order.append(name)
        referenced_terms[name] += coefficient * scale

    terms: list[str] = []
    if numeric_parsed_terms:
        coefficients = _gibbs_coefficients_from_parsed_terms(numeric_parsed_terms)
        if coefficients is None:
            raise ValueError(f"Unsupported TDB Gibbs expression: {expression}")
        numeric_expression = _gibbs_row_expression(coefficients)
        if numeric_expression != "0":
            terms.extend(_signed_expression_terms(numeric_expression))
    for name in referenced_order:
        coefficient = referenced_terms[name]
        if coefficient == 0.0:
            continue
        terms.append(_format_reference_term(name, coefficient))

    if not terms:
        return "0"
    result = "".join(_normalize_reference_signed_term(term) for term in terms)
    if not result:
        return "0"
    return result[1:] if result.startswith("+") else result


def _expanded_scaled_terms(
    expression: str,
    scale: float = 1.0,
) -> list[tuple[str, float]]:
    """Return additive terms with simple scaled parentheses distributed."""
    terms: list[tuple[str, float]] = []
    for term in _signed_expression_terms(expression):
        grouped = _scaled_parenthesized_term(term)
        if grouped is not None:
            group_scale, inner = grouped
            terms.extend(_expanded_scaled_terms(inner, scale * group_scale))
            continue
        terms.append((term, scale))
    return terms


def _contains_pressure_variable(expression: str) -> bool:
    """Return whether a Gibbs expression term references pressure."""
    for match in re.finditer(r"(?<![A-Za-z0-9_])([A-Za-z_][A-Za-z0-9_]*)", expression):
        if match.group(1).upper() == "P":
            return True
    return False


def _reference_term(term: str) -> tuple[str, float] | None:
    """Return ``(function_name, coefficient)`` for simple referenced functions."""
    compact = _normalize_numeric_exponents(
        _normalize_bare_decimal_exponents(term.replace(" ", ""))
    )
    compact = re.sub(r"(?i)\bLOG\(", "LN(", compact)
    identifiers = [
        identifier.upper()
        for identifier in re.findall(
            r"(?<![A-Za-z0-9_])([A-Za-z_][A-Za-z0-9_]*#?)",
            compact,
        )
    ]
    references = [
        identifier
        for identifier in identifiers
        if identifier.rstrip("#") not in {"T", "LN", "LOG", "EXP", "R"}
    ]
    if not references:
        return None
    if len(references) != 1:
        raise ValueError(f"Unsupported TDB Gibbs expression: {term}")

    reference = references[0]
    reference_pattern = re.escape(reference)
    if reference.endswith("#"):
        reference_pattern = re.escape(reference[:-1]) + "#"
    match = re.fullmatch(
        rf"(?P<prefix>[+\-]?(?:(?:R|(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+\-]?\d+)?)\*)*)"
        rf"(?P<name>{reference_pattern})",
        compact,
        flags=re.IGNORECASE,
    )
    if match is None:
        raise ValueError(f"Unsupported TDB Gibbs expression: {term}")
    coefficient = _constant_with_r_value(match.group("prefix").rstrip("*"))
    if coefficient is None:
        raise ValueError(f"Unsupported TDB Gibbs expression: {term}")
    return reference.rstrip("#"), coefficient


def _format_reference_term(name: str, coefficient: float) -> str:
    """Return a signed TDB function-reference term."""
    reference_name = f"{_tdb_identifier(name).rstrip('#')}#"
    if coefficient == 1.0:
        return f"+{reference_name}"
    if coefficient == -1.0:
        return f"-{reference_name}"
    return f"{float(coefficient):+.16g}*{reference_name}"


def _normalize_reference_signed_term(term: str) -> str:
    """Normalize spaces and power syntax in a signed expression term."""
    normalized = re.sub(r"\s+", "", term)
    normalized = re.sub(r"(?i)\bLOG\(", "LN(", normalized)
    normalized = re.sub(r"(?i)\bLN\(", "LN(", normalized)
    normalized = re.sub(
        r"T\*\*(?!\()([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+\-]?\d+)?)",
        r"T**(\1)",
        normalized,
    )
    return normalized


def _simplify_symbolic_functions(expression: str) -> str:
    """Simplify supported symbolic function calls before strict export."""
    text = re.sub(r"(?i)\bLOG\(", "LN(", expression)
    text = re.sub(
        r"(?i)\bEXP\(\s*\(?\s*"
        r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+\-]?\d+)?)"
        r"\s*\*\s*LN\(T\)\s*\)?\s*\)",
        lambda match: f"T**({_tdb_number(float(match.group(1)))})",
        text,
    )
    return re.sub(r"(?i)\bEXP\(\s*\(?\s*LN\(T\)\s*\)?\s*\)", "T", text)


def _gibbs_coefficients_from_segment(expression: str) -> list[float] | None:
    """Return compact Gibbs coefficients for one supported TDB expression."""
    parsed_terms = _parsed_gibbs_terms(expression)
    if parsed_terms is None:
        return None
    return _gibbs_coefficients_from_parsed_terms(parsed_terms)


def _gibbs_coefficients_from_parsed_terms(
    parsed_terms: list[tuple[float | str, float]],
) -> list[float] | None:
    """Return compact Gibbs coefficients from already parsed basis terms."""
    coefficients = [0.0] * 14
    custom_pairs: list[tuple[float, float]] = []
    for parsed in parsed_terms:
        power, coefficient = parsed
        if power == "TLOGT":
            coefficients[2] += coefficient
        elif power == 0.0:
            coefficients[0] += coefficient
        elif power == 1.0:
            coefficients[1] += coefficient
        elif power == 2.0:
            coefficients[3] += coefficient
        elif power == 3.0:
            coefficients[4] += coefficient
        elif power == -1.0:
            coefficients[5] += coefficient
        else:
            for index, (_existing_coefficient, existing_power) in enumerate(
                custom_pairs,
            ):
                if power == existing_power:
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


def _parsed_gibbs_terms(
    expression: str,
    scale: float = 1.0,
) -> list[tuple[float | str, float]] | None:
    """Parse supported Gibbs terms, distributing simple scaled parentheses."""
    parsed_terms: list[tuple[float | str, float]] = []
    for term in _signed_expression_terms(expression):
        grouped = _scaled_parenthesized_term(term)
        if grouped is not None:
            group_scale, inner = grouped
            parsed_terms.extend(_parsed_gibbs_terms(inner, scale * group_scale))
            continue
        parsed = _parse_gibbs_term(term)
        if parsed is None:
            return None
        power, coefficient = parsed
        parsed_terms.append((power, coefficient * scale))
    return parsed_terms


def _scaled_parenthesized_term(term: str) -> tuple[float, str] | None:
    text = term.replace(" ", "")
    for index in range(len(text) - 1):
        if text[index : index + 2] != "*(":
            continue
        scale = _constant_with_r_value(text[:index])
        if scale is None:
            continue
        if _matching_parenthesis_index(text, index + 1) == len(text) - 1:
            return scale, text[index + 2 : -1]
    return None


def _matching_parenthesis_index(text: str, open_index: int) -> int | None:
    depth = 0
    for index in range(open_index, len(text)):
        if text[index] == "(":
            depth += 1
        elif text[index] == ")":
            depth -= 1
            if depth == 0:
                return index
    return None


def _signed_expression_terms(expression: str) -> list[str]:
    """Split expression into additive terms without splitting exponent signs."""
    compact = _normalize_numeric_exponents(
        _normalize_bare_decimal_exponents(expression.replace(" ", ""))
    )
    compact = re.sub(r"(?i)\bLOG\(T\)", "LN(T)", compact)
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
            and not _is_exponent_sign(compact, index)
            and compact[max(start, index - 2) : index] != "**"
        ):
            terms.append(compact[start:index])
            start = index
    terms.append(compact[start:])
    return [term for term in terms if term]


def _is_exponent_sign(text: str, index: int) -> bool:
    """Return whether ``text[index]`` is the sign in scientific notation."""
    if index <= 1 or text[index - 1] not in {"E", "e"}:
        return False
    return text[index - 2].isdigit() or text[index - 2] == "."


def _parse_gibbs_term(term: str) -> tuple[float | str, float] | None:
    """Parse one supported Gibbs basis term."""
    term = re.sub(r"(?i)\bLOG\(T\)", "LN(T)", term)
    term = re.sub(r"(?i)\bLN\(T\)", "LN(T)", term)
    if term.endswith("LN(T)") and not term.endswith("T*LN(T)"):
        coefficient_text = term[: -len("LN(T)")].rstrip("*")
        coefficient = _coefficient_from_prefix(coefficient_text)
        return (99.0, coefficient) if coefficient is not None else None
    if "T" not in term:
        constant = _constant_with_r_value(term)
        return (0.0, constant) if constant is not None else None
    if term.endswith("T*LN(T)"):
        coefficient = _coefficient_from_prefix(term[: -len("T*LN(T)")].rstrip("*"))
        return ("TLOGT", coefficient) if coefficient is not None else None
    if "*T" in term:
        coefficient_text, suffix = term.split("*T", 1)
        coefficient = _coefficient_from_prefix(coefficient_text)
        if coefficient is None:
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


def _coefficient_from_prefix(prefix: str) -> float | None:
    """Return numeric coefficient from an optional multiplication prefix."""
    if prefix in {"", "+"}:
        return 1.0
    if prefix == "-":
        return -1.0
    return _constant_with_r_value(prefix.rstrip("*"))


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


def _gibbs_row_expression(coefficients: list[float]) -> str:
    """Return canonical TDB text from Gibbs coefficients."""
    terms: list[str] = []
    _append_gibbs_term(terms, coefficients[0], "")
    _append_gibbs_term(terms, coefficients[1], "*T")
    _append_gibbs_term(terms, coefficients[2], "*T*LN(T)")
    _append_gibbs_term(terms, coefficients[3], "*T**(2)")
    _append_gibbs_term(terms, coefficients[4], "*T**(3)")
    _append_gibbs_term(terms, coefficients[5], "*T**(-1)")
    for index in range(4):
        coefficient = coefficients[6 + 2 * index]
        power = coefficients[7 + 2 * index]
        if coefficient == 0.0:
            continue
        if power == 99.0:
            _append_gibbs_term(terms, coefficient, "*LN(T)")
        else:
            _append_gibbs_term(terms, coefficient, f"*T**({_tdb_number(power)})")
    return "".join(terms).lstrip("+") or "0"


def _append_gibbs_term(terms: list[str], coefficient: float, suffix: str) -> None:
    if coefficient == 0.0:
        return
    terms.append(f"{float(coefficient):+.16g}{suffix}")


def _wrap_tdb_command(prefix: str, text: str) -> list[str]:
    """Wrap command text to Thermo-Calc's 78-character TDB line limit."""
    remaining = text.strip()
    line_prefix = prefix
    lines: list[str] = []

    while remaining:
        available = _MAX_TDB_LINE_LENGTH - len(line_prefix)
        if available <= 8:
            available = _MAX_TDB_LINE_LENGTH - len(_CONTINUATION_PREFIX)
            line_prefix = _CONTINUATION_PREFIX
        if len(remaining) <= available:
            lines.append(f"{line_prefix}{remaining}")
            break

        split_at = _tdb_wrap_index(remaining, available)
        lines.append(f"{line_prefix}{remaining[:split_at].rstrip()}")
        remaining = remaining[split_at:].lstrip()
        line_prefix = _CONTINUATION_PREFIX
    return lines


def _tdb_wrap_index(text: str, limit: int) -> int:
    """Choose a stable split point without breaking exponent syntax."""
    search_start = max(1, limit - 28)
    for index in range(limit, search_start - 1, -1):
        char = text[index - 1]
        if char.isspace():
            return index

    for index in range(limit, search_start - 1, -1):
        if _is_expression_break(text, index):
            return index

    return limit


def _is_expression_break(text: str, index: int) -> bool:
    """Return whether the line can break before ``text[index]``."""
    if index <= 0 or index >= len(text):
        return False
    char = text[index]
    previous = text[index - 1]
    if char in "+-":
        return previous not in "Ee("
    return False


def _conventional_element_symbol(symbol: str) -> str:
    text = str(symbol).strip()
    upper = text.upper()
    if upper in _SPECIAL_ELEMENT_SYMBOLS:
        return upper
    if upper.isalpha() and len(upper) <= 2:
        return upper[:1] + upper[1:].lower()
    return text


def _tdb_element_symbol(symbol: str) -> str:
    upper = str(symbol).strip().upper()
    if upper in {"/-", "*"}:
        return upper
    return _conventional_element_symbol(symbol)


def _tdb_identifier(value: str) -> str:
    return str(value).strip().upper()


def _tdb_function_name(value: str) -> str:
    name = _tdb_identifier(value)
    if len(name) > _FUNCTION_NAME_WIDTH:
        raise ValueError(
            f"TDB FUNCTION names must be {_FUNCTION_NAME_WIDTH} characters or fewer: "
            f"{name!r}"
        )
    return name


def _tdb_number(value: float) -> str:
    return f"{float(value):.16g}"


def _tdb_scientific_number(value: float) -> str:
    return f"{float(value):.4E}"


def _tdb_site_ratio_number(value: float) -> str:
    """Return a readable site-ratio number while avoiding binary float noise."""
    return f"{float(value):.12g}"


def _tdb_number_text(value: str) -> str:
    try:
        return _tdb_number(float(value.replace("D", "E").replace("d", "E")))
    except ValueError:
        return value.strip()
