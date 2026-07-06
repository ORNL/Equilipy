"""Initial TDB-to-DatabaseIR parser.

This parser is intentionally conservative: it keeps the original command text
as source metadata and extracts the pieces needed by the DatabaseIR browser.
It is not yet a full Thermo-Calc/Factsage compatibility layer.
"""

from __future__ import annotations

import re
from itertools import permutations, product
from pathlib import Path

from equilipy.exceptions import DatabaseParsingError

from . import tdb_syntax as _tdb_syntax
from .model import (
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
    TdbCommand,
    TdbFunctionSource,
    TdbFunctionTerm,
    TdbReference,
)
from .tdb_canonical import DISORDERED_PHASE_CANONICAL_NAMES

_COMMAND_RE = _tdb_syntax.COMMAND_RE
_CONSTITUENT_KEYWORDS = _tdb_syntax.CONSTITUENT_KEYWORDS
_DEFAULT_COMMAND_KEYWORDS = _tdb_syntax.DEFAULT_COMMAND_KEYWORDS
_ELEMENT_KEYWORDS = _tdb_syntax.ELEMENT_KEYWORDS
_FUNCTION_KEYWORDS = _tdb_syntax.FUNCTION_KEYWORDS
_PARAMETER_KEYWORDS = _tdb_syntax.PARAMETER_KEYWORDS
_PHASE_KEYWORDS = _tdb_syntax.PHASE_KEYWORDS
_PRESERVED_COMMENT_COMMANDS = _tdb_syntax.PRESERVED_COMMENT_COMMANDS
_REFERENCE_KEYWORDS = _tdb_syntax.REFERENCE_KEYWORDS
_SPECIES_KEYWORDS = _tdb_syntax.SPECIES_KEYWORDS
_SYSTEM_DEFAULT_KEYWORDS = _tdb_syntax.SYSTEM_DEFAULT_KEYWORDS
_TYPE_DEFINITION_KEYWORDS = _tdb_syntax.TYPE_DEFINITION_KEYWORDS
_FLOAT_RE = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EDed][+-]?\d+)?")
_FUNCTION_TOKEN_RE = re.compile(r"(?<![A-Za-z0-9_])([A-Za-z_][A-Za-z0-9_]*)")
_PARAMETER_REFERENCE_RE = re.compile(
    r";\s*(?:,,|[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EDed][+-]?\d+)?)"
    r"\s+[YNyn](?:\s+(?P<reference>[A-Za-z0-9_.+\-/]+))?\s*$"
)
_REFERENCE_ENTRY_RE = re.compile(
    r"(?P<key>[A-Za-z0-9_.+\-/]+)\s+'(?P<text>(?:''|[^'])*)'",
    re.DOTALL,
)
_PSEUDO_SPECIES = {"VA", "*"}
_COMMON_SOLUTION_PHASES = {
    "LIQUID",
    "FCC_A1",
    "A1_FCC",
    "FCC_4SL",
    "BCC_A2",
    "A2_BCC",
    "BCC_B2",
    "HCP_A3",
    "CBCC_A12",
    "CUB_A13",
    "DIAMOND_A4",
}
_DEFAULT_CALPHAD_T_MAX = 6000.0
_PORTABLE_TDB_FUNCTION_NAME_WIDTH = 8
_MAX_TDB_FUNCTION_NAME_WIDTH = 25
_EXACT_ENDMEMBER_DIAGNOSTIC_LIMIT = 4096
_ENDMEMBER_DIAGNOSTIC_PREVIEW_LIMIT = 8
_ENDMEMBER_COMPOUND_SUFFIX = "_ENDMBR"
_ENDMEMBER_COMPOUND_GENERATOR = "tdb_endmember_compound_export"


def load_tdb(path: str | Path) -> DatabaseIR:
    """Parse a TDB file into the current DatabaseIR browsing model."""
    source_path = Path(path)
    source_text = source_path.read_text(encoding="utf-8", errors="replace")
    commands = list(_iter_commands(source_path))
    database = DatabaseIR(
        name=source_path.stem,
        metadata={
            "source_format": "TDB",
            "source_file": str(source_path),
            "parser": "equilipy.database_ir.tdb.load_tdb",
            **_parse_tdb_header_metadata(source_text),
        },
    )
    phase_by_name: dict[str, Phase] = {}
    phase_site_ratios: dict[str, list[float]] = {}

    for command in commands:
        keyword, body = _split_keyword(command.text)
        if not command.active and not _is_preserved_setup_keyword(keyword):
            continue
        if keyword in _ELEMENT_KEYWORDS:
            element = _parse_element(body, command.source)
            if element is not None:
                database.elements.append(element)
        elif keyword in _SPECIES_KEYWORDS:
            species = _parse_species(body, command.source)
            if species is not None:
                database.species.append(species)
        elif keyword in _FUNCTION_KEYWORDS:
            function = _parse_function(body, command.source)
            if function is not None:
                database.functions.append(function)
                _record_long_function_name(database, function)
        elif keyword in _PHASE_KEYWORDS:
            phase, site_ratios = _parse_phase(body, command.source)
            if phase is not None:
                database.phases.append(phase)
                phase_by_name[phase.name] = phase
                phase_site_ratios[phase.name] = site_ratios
        elif keyword in _CONSTITUENT_KEYWORDS:
            phase_name, constituents = _parse_constituents(
                body,
                command.source,
                phase_site_ratios,
            )
            if phase_name and constituents:
                if phase_name not in phase_by_name:
                    phase = Phase(
                        name=phase_name,
                        model="COMPOUND",
                        source=command.source,
                        metadata=_endmember_compound_phase_metadata(phase_name),
                    )
                    database.phases.append(phase)
                    phase_by_name[phase.name] = phase
                phase_by_name[phase_name].constituents = constituents
        elif keyword in _PARAMETER_KEYWORDS:
            parameter = _parse_parameter(body, command.source)
            if parameter is not None:
                database.parameters.append(parameter)
        elif keyword in _REFERENCE_KEYWORDS:
            database.references.extend(
                _parse_references(keyword, body, command.source)
            )
        elif (
            keyword in _SYSTEM_DEFAULT_KEYWORDS
            or keyword in _DEFAULT_COMMAND_KEYWORDS
            or keyword in _TYPE_DEFINITION_KEYWORDS
        ):
            database.tdb_commands.append(_parse_tdb_command(keyword, body, command))

    _canonicalize_element_references(database)
    _infer_missing_species(database)
    _classify_phase_models(database)
    _record_missing_endmember_diagnostics(database)
    _drop_pressure_dependent_functions(database)
    _sync_function_tdb_sources(database)
    _sort_ir(database)
    return database


def read_tdb(
    path: str | Path,
    *,
    strict: bool = True,
    auto_correct: bool = False,
    remove_redundant_phases: bool = False,
) -> DatabaseIR:
    """Read a TDB file into a DatabaseIR object.

    ``load_tdb`` is intentionally tolerant for the GUI: it imports what it can
    and stores parser/validation problems as diagnostics. ``read_tdb`` is the
    calculation-facing API and defaults to strict mode, raising on parser
    errors while leaving missing concrete endmember parameters as runtime
    implicit zero Gibbs compatibility records.
    """
    database = load_tdb(path)
    if remove_redundant_phases:
        remove_redundant_disordered_phase_aliases(database)
    if auto_correct:
        database.validate_tdb(auto_correct=True)
    if strict:
        _raise_for_database_errors(database)
    return database


def remove_redundant_disordered_phase_aliases(database: DatabaseIR) -> int:
    """Collapse known disordered phase aliases to their canonical phase names.

    Thermo-Calc style databases often carry helper names such as ``A2_BCC`` or
    ``A1_FCC`` for the disordered part of an ordered phase.  In Equilipy these
    are not distinct phases: they are aliases for ``BCC_A2`` and ``FCC_A1``.
    Prefer already-defined canonical phase records and only migrate alias-only
    records when the canonical record is absent.
    """
    phase_names = {phase.name.upper(): phase.name for phase in database.phases}
    aliases: dict[str, str] = {}
    for alias, canonical in DISORDERED_PHASE_CANONICAL_NAMES.items():
        alias_name = phase_names.get(alias)
        if alias_name is None:
            continue
        aliases[alias_name] = phase_names.get(canonical, canonical)
    if not aliases:
        return 0

    alias_names = set(aliases)
    canonical_phase_names = {name.upper() for name in aliases.values()}
    collapsed_phases: list[Phase] = []
    emitted_phase_names: set[str] = set()
    for phase in database.phases:
        canonical_name = aliases.get(phase.name)
        if canonical_name is not None and canonical_name.upper() in phase_names:
            continue
        if canonical_name is not None:
            phase.name = canonical_name
        phase_key = phase.name.upper()
        if phase_key in emitted_phase_names:
            continue
        collapsed_phases.append(phase)
        emitted_phase_names.add(phase_key)
    database.phases = collapsed_phases

    canonical_keys = {
        _parameter_identity(parameter)
        for parameter in database.parameters
        if parameter.phase_name.upper() in canonical_phase_names
    }
    pruned_parameters: list[Parameter] = []
    emitted_keys: set[tuple[str, str, tuple[str, ...], int]] = set()
    for parameter in database.parameters:
        alias_target = aliases.get(parameter.phase_name)
        if alias_target is not None:
            parameter.phase_name = alias_target
        key = _parameter_identity(parameter)
        if alias_target is not None and key in canonical_keys:
            continue
        if key in emitted_keys:
            continue
        pruned_parameters.append(parameter)
        emitted_keys.add(key)
    database.parameters = pruned_parameters
    for command in database.tdb_commands:
        parsed = command.parsed
        if parsed.get("action", "").upper() != "DIS_PART":
            continue
        disordered = parsed.get("disordered_phase", "")
        canonical = DISORDERED_PHASE_CANONICAL_NAMES.get(disordered.upper())
        if not canonical:
            continue
        parsed["disordered_phase"] = canonical
        command.args = [
            canonical if arg.upper() == disordered.upper() else arg
            for arg in command.args
        ]
    _drop_missing_endmember_diagnostics(database, alias_names)
    return len(aliases)


def _parameter_identity(
    parameter: Parameter,
) -> tuple[str, str, tuple[str, ...], int]:
    return (
        parameter.phase_name.upper(),
        parameter.parameter_type.upper(),
        tuple(target.upper() for target in parameter.target),
        int(parameter.order),
    )


def _drop_missing_endmember_diagnostics(
    database: DatabaseIR,
    phase_names: set[str],
) -> None:
    upper_phase_names = {name.upper() for name in phase_names}
    filtered: list[Diagnostic] = []
    for diagnostic in database.diagnostics:
        message = diagnostic.message
        if " is missing G endmember parameters " not in message:
            filtered.append(diagnostic)
            continue
        if not upper_phase_names:
            continue
        match = re.search(r"Phase '([^']+)' is missing G endmember parameters", message)
        if match is not None and match.group(1).upper() in upper_phase_names:
            continue
        filtered.append(diagnostic)
    database.diagnostics = filtered


def _raise_for_database_errors(database: DatabaseIR) -> None:
    """Raise a parsing error when the imported TDB has hard diagnostics."""
    errors = [
        diagnostic
        for diagnostic in database.validation_messages()
        if diagnostic.severity.lower() == "error"
    ]
    if not errors:
        return
    details = "\n".join(
        _format_diagnostic_for_exception(diagnostic) for diagnostic in errors[:8]
    )
    if len(errors) > 8:
        details += f"\n... {len(errors) - 8} additional error(s)."
    raise DatabaseParsingError(
        f"TDB parsing failed with {len(errors)} error(s):\n{details}"
    )


def _format_diagnostic_for_exception(diagnostic: Diagnostic) -> str:
    """Return a compact source-aware diagnostic message."""
    source = diagnostic.source
    if source.file and source.line:
        return f"{source.file}:{source.line}: {diagnostic.message}"
    return diagnostic.message


def _is_preserved_setup_keyword(keyword: str) -> bool:
    """Return whether an inactive/commented command is safe to keep."""
    return (
        keyword in _SYSTEM_DEFAULT_KEYWORDS
        or keyword in _DEFAULT_COMMAND_KEYWORDS
        or keyword in _TYPE_DEFINITION_KEYWORDS
    )


def _parse_tdb_header_metadata(text: str) -> dict[str, str]:
    """Extract user-facing title/author/date/description comments."""
    comment_lines: list[str] = []
    for raw_line in text.splitlines():
        stripped = raw_line.lstrip()
        if stripped.startswith("$"):
            comment_lines.append(stripped[1:].rstrip())
            continue
        if raw_line.strip():
            break

    metadata: dict[str, str] = {}
    title = _parse_header_title(comment_lines)
    if title:
        metadata["title"] = title

    description_parts: list[str] = []
    collecting_description = False
    for line in comment_lines:
        stripped = line.strip()
        lower = stripped.lower()
        if lower.startswith("author:"):
            metadata["authors"] = stripped.split(":", 1)[1].strip()
            collecting_description = False
            continue
        if lower.startswith("program:"):
            metadata["program"] = stripped.split(":", 1)[1].strip()
            collecting_description = False
            continue
        if lower.startswith(("date:", "date modified:", "modified:")):
            metadata["date_modified"] = stripped.split(":", 1)[1].strip()
            collecting_description = False
            continue
        if lower.startswith("description:"):
            first = stripped.split(":", 1)[1].strip()
            if first:
                description_parts.append(first)
            collecting_description = True
            continue
        if collecting_description:
            if (
                not stripped
                or set(stripped) <= {"-", "="}
                or lower.startswith(("alloy elements", "attention:"))
            ):
                collecting_description = False
                continue
            description_parts.append(stripped)

    if description_parts:
        metadata["description"] = " ".join(description_parts)
    return {key: value for key, value in metadata.items() if value}


def _parse_header_title(comment_lines: list[str]) -> str:
    separator_indices = [
        index
        for index, line in enumerate(comment_lines)
        if line.strip() and set(line.strip()) == {"="}
    ]
    if len(separator_indices) < 2:
        return ""
    start, end = separator_indices[0], separator_indices[1]
    for line in comment_lines[start + 1 : end]:
        stripped = line.strip()
        if stripped:
            return stripped
    return ""


class _Command:
    def __init__(self, text: str, source: SourceRef, *, active: bool = True):
        self.text = text
        self.source = source
        self.active = active


def _iter_commands(path: Path) -> list[_Command]:
    """Yield TDB commands separated by ``!`` while skipping comment lines."""
    commands: list[_Command] = []
    buffer: list[str] = []
    start_line = 0

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            stripped = raw_line.lstrip()
            if stripped.startswith("$"):
                comment_text = stripped[1:].strip()
                if not _looks_like_preserved_comment_command(comment_text):
                    continue
                commands.extend(_comment_commands(comment_text, path, line_number))
                continue
            line = _strip_inline_comment(raw_line.rstrip("\n"))
            if not line.strip():
                continue
            if not buffer:
                start_line = line_number
            buffer.append(line)

            while "!" in buffer[-1]:
                joined = "\n".join(buffer)
                before, after = joined.split("!", 1)
                text = _clean_command_text(before)
                if text:
                    commands.append(
                        _Command(
                            text,
                            SourceRef(
                                file=str(path),
                                line=start_line,
                                column=1,
                                command=f"{text} !",
                            ),
                            active=True,
                        )
                    )
                if after.strip():
                    buffer = [after]
                    start_line = line_number
                else:
                    buffer = []
                    start_line = 0
                    break

    if buffer:
        text = _clean_command_text("\n".join(buffer))
        if text:
            commands.append(
                _Command(
                    text,
                    SourceRef(
                        file=str(path),
                        line=start_line,
                        column=1,
                        command=text,
                    ),
                    active=True,
                )
            )
    return commands


def _strip_inline_comment(line: str) -> str:
    """Strip ordinary TDB ``$`` comments from active command lines."""
    return line.split("$", 1)[0].rstrip()


def _looks_like_preserved_comment_command(text: str) -> bool:
    """Return whether a TDB comment line contains a setup command to preserve."""
    if not text:
        return False
    keyword, _body = _split_keyword(text.split("!", 1)[0])
    return keyword in _PRESERVED_COMMENT_COMMANDS


def _comment_commands(text: str, path: Path, line_number: int) -> list[_Command]:
    """Return inactive command records embedded in a TDB comment line."""
    commands: list[_Command] = []
    for before in text.split("!"):
        cleaned = _clean_command_text(before)
        if not cleaned:
            continue
        keyword, _body = _split_keyword(cleaned)
        if keyword not in _PRESERVED_COMMENT_COMMANDS:
            continue
        commands.append(
            _Command(
                cleaned,
                SourceRef(
                    file=str(path),
                    line=line_number,
                    column=1,
                    command=f"${cleaned} !",
                ),
                active=False,
            )
        )
    return commands


def _clean_command_text(text: str) -> str:
    """Normalize whitespace without erasing expression content."""
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    return " ".join(lines)


def _split_keyword(command: str) -> tuple[str, str]:
    match = _COMMAND_RE.match(command)
    if match is None:
        return "", command
    return match.group(1).upper(), match.group(2).strip()


def _parse_tdb_command(keyword: str, body: str, command: _Command) -> TdbCommand:
    """Preserve a TDB setup/type command with lightweight parsed metadata."""
    args = body.split()
    return TdbCommand(
        kind=_tdb_command_kind(keyword),
        command=_canonical_tdb_setup_command(keyword),
        args=args,
        active=command.active,
        raw=command.text,
        parsed=_parse_tdb_command_metadata(keyword, args),
        source=command.source,
    )


def _tdb_command_kind(keyword: str) -> str:
    if keyword in _SYSTEM_DEFAULT_KEYWORDS:
        return "system_default"
    if keyword in _DEFAULT_COMMAND_KEYWORDS:
        return "default_command"
    if keyword in _TYPE_DEFINITION_KEYWORDS:
        return "type_definition"
    return "generic"


def _canonical_tdb_setup_command(keyword: str) -> str:
    if keyword in _SYSTEM_DEFAULT_KEYWORDS:
        return "DEFINE_SYSTEM_DEFAULT"
    if keyword in _DEFAULT_COMMAND_KEYWORDS:
        return "DEFAULT_COMMAND"
    if keyword in _TYPE_DEFINITION_KEYWORDS:
        return "TYPE_DEF"
    return keyword


def _parse_tdb_command_metadata(keyword: str, args: list[str]) -> dict[str, str]:
    """Return recognized metadata for a preserved setup/type command."""
    if keyword in _SYSTEM_DEFAULT_KEYWORDS:
        return {
            "subject": args[0].upper() if args else "",
            "value": " ".join(args[1:]),
        }
    if keyword in _DEFAULT_COMMAND_KEYWORDS:
        return {
            "action": args[0].upper() if args else "",
            "arguments": " ".join(args[1:]),
        }
    if keyword in _TYPE_DEFINITION_KEYWORDS:
        return _parse_type_definition_metadata(args)
    return {}


def _parse_type_definition_metadata(args: list[str]) -> dict[str, str]:
    """Recognize common TYPE_DEFINITION payloads without executing them."""
    metadata: dict[str, str] = {
        "code": args[0] if args else "",
        "module": args[1].upper() if len(args) > 1 else "",
    }
    upper_args = [arg.upper() for arg in args]
    if "MAGNETIC" in upper_args:
        index = upper_args.index("MAGNETIC")
        metadata.update(
            {
                "action": "MAGNETIC",
                "target_phase": args[3] if len(args) > 3 else "",
                "magnetic_model": args[index + 1] if len(args) > index + 1 else "",
                "magnetic_factor": args[index + 2] if len(args) > index + 2 else "",
            }
        )
    elif "DIS_PART" in upper_args or "DISORDER_PART" in upper_args:
        keyword = "DIS_PART" if "DIS_PART" in upper_args else "DISORDER_PART"
        index = upper_args.index(keyword)
        metadata.update(
            {
                "action": "DIS_PART",
                "ordered_phase": args[3] if len(args) > 3 else "",
                "disordered_phase": args[index + 1] if len(args) > index + 1 else "",
                "dialect_keyword": keyword,
            }
        )
    elif "SEQ" in upper_args:
        index = upper_args.index("SEQ")
        metadata.update(
            {
                "action": "SEQ",
                "sequence": args[index + 1] if len(args) > index + 1 else "",
            }
        )
    else:
        metadata["action"] = upper_args[2] if len(upper_args) > 2 else ""
    return metadata


def _parse_references(
    keyword: str,
    body: str,
    source: SourceRef,
) -> list[TdbReference]:
    """Parse a TDB reference table command into individual entries."""
    section = _canonical_reference_section(keyword)
    references: list[TdbReference] = []
    for match in _REFERENCE_ENTRY_RE.finditer(body):
        key = match.group("key").strip()
        text = " ".join(match.group("text").replace("''", "'").split())
        if not key:
            continue
        references.append(
            TdbReference(
                key=key,
                text=text,
                section=section,
                source=source,
            )
        )
    return references


def _canonical_reference_section(keyword: str) -> str:
    """Return the canonical reference table command spelling."""
    normalized = keyword.upper().replace("_", "-")
    if normalized.startswith("ADD-"):
        return "ADD_REFERENCES"
    return "LIST-OF-REFERENCE"


def _parse_element(body: str, source: SourceRef) -> Element | None:
    parts = body.split()
    if len(parts) < 2:
        return None
    atomic_mass = _to_float(parts[2], 0.0) if len(parts) > 2 else 0.0
    h298 = _to_float(parts[3], 0.0) if len(parts) > 3 else 0.0
    s298 = _to_float(parts[4], 0.0) if len(parts) > 4 else 0.0
    return Element(
        symbol=_normalize_element_symbol(parts[0]),
        reference_state=_normalize_name(parts[1]),
        atomic_mass=atomic_mass,
        source=source,
        h298=h298,
        s298=s298,
        tdb_symbol=_normalize_name(parts[0]),
    )


def _parse_species(body: str, source: SourceRef) -> Species | None:
    parts = body.split(None, 1)
    if not parts:
        return None
    name = _normalize_name(parts[0])
    composition = _formula_composition(name)
    return Species(name=name, composition=composition, source=source)


def _parse_function(body: str, source: SourceRef) -> FunctionDefinition | None:
    parts = body.split(None, 1)
    if len(parts) < 2:
        return None
    name = _normalize_name(parts[0])
    expression = parts[1].strip()
    gibbs_ranges = _parse_gibbs_ranges(expression)
    return FunctionDefinition(
        name=name,
        expression=expression,
        temperature_ranges=[
            (gibbs_range.t_min, gibbs_range.t_max) for gibbs_range in gibbs_ranges
        ],
        source=source,
        gibbs_ranges=gibbs_ranges,
    )


def _function_depends_on_pressure(function: FunctionDefinition) -> bool:
    """Return whether a parsed TDB FUNCTION explicitly references pressure."""
    expressions = [function.expression]
    expressions.extend(gibbs_range.Gibbs for gibbs_range in function.gibbs_ranges)
    for expression in expressions:
        for match in _FUNCTION_TOKEN_RE.finditer(expression):
            if match.group(1).upper() == "P":
                return True
    return False


def _drop_pressure_dependent_functions(database: DatabaseIR) -> None:
    """Remove functions that directly or indirectly reference pressure."""
    by_name = {function.name.upper(): function for function in database.functions}
    pressure_names = {
        function.name.upper()
        for function in database.functions
        if _function_depends_on_pressure(function)
    }
    changed = True
    while changed:
        changed = False
        for function in database.functions:
            name = function.name.upper()
            if name in pressure_names:
                continue
            referenced = set(_referenced_function_names(function.expression, by_name))
            if referenced.intersection(pressure_names):
                pressure_names.add(name)
                changed = True

    if not pressure_names:
        return

    kept: list[FunctionDefinition] = []
    for function in database.functions:
        if function.name.upper() in pressure_names:
            _record_dropped_pressure_function(database, function)
        else:
            kept.append(function)
    database.functions = kept


def _record_dropped_pressure_function(
    database: DatabaseIR,
    function: FunctionDefinition,
) -> None:
    """Attach a diagnostic for a skipped pressure-dependent TDB function."""
    dropped = _metadata_name_list(database.metadata.get("dropped_pressure_functions"))
    if function.name not in dropped:
        dropped.append(function.name)
    database.metadata["dropped_pressure_functions"] = ", ".join(dropped)
    database.diagnostics.append(
        Diagnostic(
            "warning",
            (
                f"Dropped pressure-dependent function '{function.name}'. "
                "Equilipy applies the gas pressure correction internally, so this "
                "TDB pressure function was not imported."
            ),
            function.source,
        )
    )


def _record_long_function_name(
    database: DatabaseIR,
    function: FunctionDefinition,
) -> None:
    """Attach diagnostics for FUNCTION names outside portable TDB conventions."""
    name_length = len(function.name)
    if name_length <= _PORTABLE_TDB_FUNCTION_NAME_WIDTH:
        return
    if name_length > _MAX_TDB_FUNCTION_NAME_WIDTH:
        database.diagnostics.append(
            Diagnostic(
                severity="error",
                message=(
                    f"TDB FUNCTION name '{function.name}' has {name_length} "
                    "characters. Equilipy supports FUNCTION names up to "
                    f"{_MAX_TDB_FUNCTION_NAME_WIDTH} characters."
                ),
                source=function.source,
            )
        )
        return
    database.diagnostics.append(
        Diagnostic(
            severity="warning",
            message=(
                f"TDB FUNCTION name '{function.name}' has {name_length} "
                "characters. Some TDB tools expect FUNCTION names of at most "
                f"{_PORTABLE_TDB_FUNCTION_NAME_WIDTH} characters."
            ),
            source=function.source,
        )
    )


def _metadata_name_list(value: str | None) -> list[str]:
    """Return comma-separated metadata names as a clean list."""
    if not value:
        return []
    return [name.strip() for name in value.split(",") if name.strip()]


def _endmember_compound_phase_metadata(phase_name: str) -> dict[str, str]:
    """Return marker metadata for generated endmember-compound phase copies."""
    if not _is_endmember_compound_name(phase_name):
        return {}
    return {
        "generated_by": _ENDMEMBER_COMPOUND_GENERATOR,
        "endmember_compound": "true",
    }


def _is_endmember_compound_name(phase_name: str) -> bool:
    """Return whether a phase name carries Equilipy's generated suffix."""
    identifier = _normalize_name(phase_name).upper()
    return identifier.endswith(_ENDMEMBER_COMPOUND_SUFFIX) or identifier.endswith(
        "_ENDMEMBER"
    )


def _parse_phase(body: str, source: SourceRef) -> tuple[Phase | None, list[float]]:
    parts = body.split()
    if len(parts) < 3:
        return None, []

    phase_name = _phase_name(parts[0])
    try:
        sublattice_count = int(float(parts[2]))
    except ValueError:
        sublattice_count = 0
    site_tokens = parts[3 : 3 + sublattice_count]
    site_ratios = [_to_float(token, 1.0) for token in site_tokens]
    if sublattice_count and len(site_ratios) < sublattice_count:
        site_ratios.extend([1.0] * (sublattice_count - len(site_ratios)))

    phase = Phase(
        name=phase_name,
        model="COMPOUND",
        constituents=[
            ConstituentSet(index, [], site_ratios[index - 1])
            for index in range(1, len(site_ratios) + 1)
        ],
        source=source,
        metadata=_endmember_compound_phase_metadata(phase_name),
    )
    return phase, site_ratios


def _parse_constituents(
    body: str,
    source: SourceRef,
    phase_site_ratios: dict[str, list[float]],
) -> tuple[str, list[ConstituentSet]]:
    parts = body.split(None, 1)
    if not parts:
        return "", []
    phase_name = _phase_name(parts[0])
    constituent_part = parts[1] if len(parts) > 1 else ""
    constituent_part = constituent_part.strip()
    if constituent_part.startswith(":"):
        constituent_part = constituent_part[1:]
    site_ratios = phase_site_ratios.get(phase_name, [])
    constituents: list[ConstituentSet] = []

    sections = [section for section in constituent_part.split(":") if section.strip()]
    for index, section in enumerate(sections, start=1):
        species = [
            _normalize_name(name)
            for name in re.split(r"[,\s]+", section)
            if _normalize_name(name)
        ]
        site_ratio = site_ratios[index - 1] if index <= len(site_ratios) else 1.0
        constituents.append(ConstituentSet(index, species, site_ratio))

    return phase_name, constituents


def _parse_parameter(body: str, source: SourceRef) -> Parameter | None:
    match = re.match(
        r"^\s*([A-Za-z0-9_]+)\s*\(\s*([^,()]+)\s*,\s*(.*?)\s*\)\s*(.*)$",
        body,
        re.DOTALL,
    )
    if match is None:
        return None

    parameter_type = match.group(1).upper()
    phase_name = _phase_name(match.group(2))
    target_and_order = match.group(3).strip()
    expression_text, reference = _split_parameter_reference(match.group(4).strip())
    expression = _normalize_parameter_expression(expression_text)
    target, order = _split_target_order(target_and_order)
    metadata = {"source_expression": expression} if expression else {}
    if reference:
        metadata["reference"] = reference
    return Parameter(
        phase_name=phase_name,
        parameter_type=parameter_type,
        target=[target] if target else [],
        order=order,
        expression=expression,
        source=source,
        metadata=metadata,
    )


def _split_parameter_reference(expression: str) -> tuple[str, str]:
    """Return a parameter expression without its trailing TDB reference key."""
    match = _PARAMETER_REFERENCE_RE.search(expression)
    if match is None:
        return expression, ""
    reference = (match.group("reference") or "").strip()
    if not reference:
        return expression, ""
    return expression[: match.start("reference")].rstrip(), reference


def _normalize_parameter_expression(expression: str) -> str:
    """Return a strict numeric temperature range expression for TDB parameters."""
    ranges = _parse_gibbs_ranges(expression)
    if not ranges:
        return expression
    pieces: list[str] = []
    for index, gibbs_range in enumerate(ranges):
        prefix = f"{_tdb_number_text(gibbs_range.t_min)} " if index == 0 else ""
        pieces.append(
            f"{prefix}{gibbs_range.Gibbs.strip()}; "
            f"{_tdb_number_text(gibbs_range.t_max)} {gibbs_range.status.upper()}"
        )
    return " ".join(pieces)


def _tdb_number_text(value: float) -> str:
    """Return compact TDB numeric text."""
    return f"{float(value):.15g}"


def _split_target_order(value: str) -> tuple[str, int]:
    if ";" not in value:
        return value.strip(), 0
    target, order_text = value.rsplit(";", 1)
    return target.strip(), int(_to_float(order_text.strip(), 0.0))


def _parse_gibbs_ranges(expression: str) -> list[GibbsRange]:
    """Parse piecewise TDB Gibbs ranges from a FUNCTION expression."""
    text = expression.strip()
    if text.startswith(",,"):
        lower = 298.15
        remainder = text[2:].lstrip(" ,")
    else:
        lower_match = _FLOAT_RE.match(text)
        if lower_match is None:
            return []
        lower = _to_float(lower_match.group(0), 0.0)
        remainder = text[lower_match.end() :].lstrip()
    ranges: list[GibbsRange] = []

    while ";" in remainder:
        segment_expression, remainder = remainder.split(";", 1)
        remainder = remainder.lstrip()
        upper_match = _FLOAT_RE.match(remainder)
        if upper_match is None:
            upper = _DEFAULT_CALPHAD_T_MAX
            tail = remainder.lstrip(" ,")
        else:
            upper = _to_float(upper_match.group(0), lower)
            tail = remainder[upper_match.end() :].lstrip()
        status = tail[:1].upper() if tail else "N"
        if status not in {"Y", "N"}:
            status = "N"
        else:
            tail = tail[1:].lstrip()
        segment_expression = _strip_repeated_range_lower(
            segment_expression.strip(),
            lower,
        )
        ranges.append(
            GibbsRange(
                T_min=lower,
                T_max=upper,
                Gibbs=segment_expression,
                status=status,
            )
        )
        lower = upper
        remainder = tail

    return ranges


def _strip_repeated_range_lower(segment_expression: str, lower: float) -> str:
    """Remove a duplicated lower-bound token from a range body.

    Older development exports accidentally wrote each range as
    ``...; 700 Y 700 expr; 933 Y ...``.  That second ``700`` is not part of a
    valid Gibbs expression, and leaving it in the IR later makes the writer try
    to canonicalize ``700 expr`` as a single segment.  Strip only the exact
    whitespace-separated lower-bound duplicate so legitimate coefficients such
    as ``700*T`` remain untouched.
    """
    match = _FLOAT_RE.match(segment_expression)
    if match is None:
        return segment_expression
    tail = segment_expression[match.end() :]
    if not tail or not tail[0].isspace():
        return segment_expression
    token_value = _to_float(match.group(0), float("nan"))
    if abs(token_value - float(lower)) > 1e-8:
        return segment_expression
    stripped = tail.lstrip()
    return stripped or segment_expression


def _temperature_ranges(expression: str) -> list[tuple[float, float]]:
    return [
        (gibbs_range.t_min, gibbs_range.t_max)
        for gibbs_range in _parse_gibbs_ranges(expression)
    ]


def _sync_function_tdb_sources(database: DatabaseIR) -> None:
    """Attach TDB provenance and split composite functions onto dependency ranges."""
    functions_by_name = {
        function.name.upper(): function for function in database.functions
    }
    for function in database.functions:
        dependency_terms = _function_dependency_terms(function, functions_by_name)
        if not dependency_terms:
            function.tdb_src = TdbFunctionSource()
            continue
        function.tdb_src = TdbFunctionSource(
            source_expression=function.expression,
            pertinent_functions=[
                TdbFunctionTerm(name=name, factor=factor)
                for name, factor in dependency_terms
            ],
        )
        function.gibbs_ranges = _split_composite_gibbs_ranges(
            function,
            dependency_terms,
            functions_by_name,
        )
        function.temperature_ranges = [
            (gibbs_range.t_min, gibbs_range.t_max)
            for gibbs_range in function.gibbs_ranges
        ]


def _function_dependency_terms(
    function: FunctionDefinition,
    functions_by_name: dict[str, FunctionDefinition],
) -> list[tuple[str, float]]:
    """Return referenced functions and multipliers for one function."""
    terms: list[tuple[str, float]] = []
    referenced_names = _referenced_function_names(
        function.expression,
        functions_by_name,
    )
    for name in referenced_names:
        if name == function.name.upper():
            continue
        factor = _expression_multiplier_for_name(function.expression, name)
        if factor is not None and factor != 0.0:
            terms.append((name, factor))
    return terms


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


def _split_composite_gibbs_ranges(
    function: FunctionDefinition,
    dependency_terms: list[tuple[str, float]],
    functions_by_name: dict[str, FunctionDefinition],
) -> list[GibbsRange]:
    """Split a composite function over every referenced-function range boundary."""
    expansion_cache: dict[tuple[str, float, float], str] = {}
    source_ranges = function.gibbs_ranges or _parse_gibbs_ranges(function.expression)
    dependency_bounds = _dependency_temperature_bounds(
        dependency_terms,
        functions_by_name,
    )
    if not source_ranges and dependency_bounds:
        source_ranges = [
            GibbsRange(
                min(dependency_bounds),
                max(dependency_bounds),
                function.expression,
                "N",
            )
        ]
    if not dependency_bounds:
        return source_ranges

    split_ranges: list[GibbsRange] = []
    for source_range in source_ranges:
        bounds = _merged_bounds(
            [
                source_range.t_min,
                *[
                    bound
                    for bound in dependency_bounds
                    if source_range.t_min < bound < source_range.t_max
                ],
                source_range.t_max,
            ]
        )
        for lower, upper in zip(bounds, bounds[1:], strict=False):
            expression = _expand_composite_expression(
                source_range.expression,
                lower,
                upper,
                functions_by_name,
                seen={function.name.upper()},
                cache=expansion_cache,
            )
            expression = _simplify_gibbs_expression(expression)
            split_ranges.append(
                GibbsRange(
                    T_min=lower,
                    T_max=upper,
                    Gibbs=expression,
                    status=source_range.status,
                )
            )
    return split_ranges


def _expand_composite_expression(
    expression: str,
    lower: float,
    upper: float,
    functions_by_name: dict[str, FunctionDefinition],
    seen: set[str],
    cache: dict[tuple[str, float, float], str] | None = None,
) -> str:
    """Replace referenced function names with the active Gibbs expression."""
    cache = {} if cache is None else cache
    expanded = expression
    referenced_names = _referenced_function_names(expanded, functions_by_name)
    for name in sorted(referenced_names, key=len, reverse=True):
        if name in seen:
            continue
        dependency = functions_by_name[name]
        dependency_range = _active_gibbs_range(dependency, lower, upper)
        if dependency_range is None:
            continue
        cache_key = (name, round(float(lower), 8), round(float(upper), 8))
        if cache_key in cache:
            dependency_expression = cache[cache_key]
        else:
            dependency_expression = _expand_composite_expression(
                dependency_range.expression,
                lower,
                upper,
                functions_by_name,
                seen={*seen, name},
                cache=cache,
            )
            cache[cache_key] = dependency_expression
        pattern = rf"(?<![A-Za-z0-9_]){re.escape(name)}#?(?![A-Za-z0-9_])"
        expanded = re.sub(
            pattern,
            f"({dependency_expression})",
            expanded,
            flags=re.IGNORECASE,
        )
    return expanded


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


def _simplify_gibbs_expression(expression: str) -> str:
    """Combine safely identifiable constant and linear-T Gibbs terms."""
    expression = _simplify_symbolic_functions(expression)
    expression = _expand_scaled_parenthesized_terms(expression)
    flattened = _strip_additive_group_parentheses(expression)
    terms = _split_additive_terms(flattened)
    constant = 0.0
    linear_t = 0.0
    remainder: list[str] = []

    for term in terms:
        coefficient = _constant_term_value(term)
        if coefficient is not None:
            constant += coefficient
            continue
        coefficient = _linear_t_term_value(term)
        if coefficient is not None:
            linear_t += coefficient
            continue
        remainder.append(term)

    combined_terms: list[str] = []
    if constant != 0.0 or not remainder and linear_t == 0.0:
        combined_terms.append(_signed_number_term(constant))
    if linear_t != 0.0:
        combined_terms.append(_signed_number_term(linear_t, suffix="*T"))
    combined_terms.extend(_normalize_signed_term(term) for term in remainder)
    result = "".join(combined_terms)
    return result[1:] if result.startswith("+") else result


def _simplify_symbolic_functions(expression: str) -> str:
    """Simplify supported symbolic function calls in expanded TDB expressions."""
    text = re.sub(r"(?i)\bLOG\(", "LN(", expression)
    text = re.sub(
        r"(?i)\bEXP\(\s*\(?\s*"
        r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+\-]?\d+)?)"
        r"\s*\*\s*LN\(T\)\s*\)?\s*\)",
        lambda match: f"T**({_format_gibbs_number(float(match.group(1)))})",
        text,
    )
    text = re.sub(r"(?i)\bEXP\(\s*\(?\s*LN\(T\)\s*\)?\s*\)", "T", text)
    return text


def _expand_scaled_parenthesized_terms(expression: str) -> str:
    """Distribute simple numeric/R factors over parenthesized additive terms."""
    expanded: list[str] = []
    for term in _split_additive_terms(expression):
        expanded.extend(_expand_scaled_parenthesized_term(term))
    result = "".join(expanded)
    return result[1:] if result.startswith("+") else result


def _expand_scaled_parenthesized_term(term: str) -> list[str]:
    grouped = _scaled_parenthesized_term(term)
    if grouped is None:
        return [term]
    scale, inner = grouped
    expanded: list[str] = []
    for inner_term in _split_additive_terms(inner):
        expanded.append(_scale_signed_term(inner_term, scale))
    return expanded


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


def _scale_signed_term(term: str, scale: float) -> str:
    body = term[1:] if term[:1] in "+-" else term
    sign = -1.0 if term.startswith("-") else 1.0
    scaled = scale * sign
    prefix = _signed_number_term(scaled)
    if not body:
        return prefix
    return f"{prefix}*{body}"


def _strip_additive_group_parentheses(expression: str) -> str:
    """Remove grouping parentheses after + or - while preserving LN(T), T**(-1)."""
    remove_indices: set[int] = set()
    stack: list[tuple[int, bool]] = []
    for index, char in enumerate(expression):
        if char == "(":
            previous = _previous_nonspace(expression, index)
            removable = previous in {"", "+", "-"}
            stack.append((index, removable))
        elif char == ")" and stack:
            open_index, removable = stack.pop()
            if removable:
                remove_indices.add(open_index)
                remove_indices.add(index)
    return "".join(
        char for index, char in enumerate(expression) if index not in remove_indices
    )


def _previous_nonspace(text: str, index: int) -> str:
    for position in range(index - 1, -1, -1):
        if not text[position].isspace():
            return text[position]
    return ""


def _split_additive_terms(expression: str) -> list[str]:
    """Split expression into signed top-level additive terms."""
    text = re.sub(r"\+\s*-", "-", expression.strip())
    text = re.sub(r"-\s*-", "+", text)
    text = re.sub(r"\+\s*\+", "+", text)
    text = re.sub(r"-\s*\+", "-", text)
    if not text:
        return []
    starts: list[int] = [0]
    depth = 0
    for index, char in enumerate(text):
        if char == "(":
            depth += 1
        elif char == ")":
            depth = max(0, depth - 1)
        elif (
            char in "+-"
            and depth == 0
            and index > 0
            and text[index - 1] not in {"E", "e"}
            and text[max(0, index - 2) : index] != "**"
        ):
            starts.append(index)
    starts.append(len(text))

    terms: list[str] = []
    for start, end in zip(starts, starts[1:], strict=False):
        term = text[start:end].strip()
        if not term:
            continue
        if term[0] not in "+-":
            term = f"+{term}"
        terms.append(term)
    return terms


def _constant_term_value(term: str) -> float | None:
    body = term[1:].strip()
    value = _constant_with_r_value(body)
    if value is None:
        return None
    sign = -1.0 if term.startswith("-") else 1.0
    return sign * value


def _linear_t_term_value(term: str) -> float | None:
    sign = -1.0 if term.startswith("-") else 1.0
    body = term[1:].replace(" ", "")
    if body == "T":
        return sign
    if body.endswith("*T"):
        coefficient = _constant_with_r_value(body[:-2])
        if coefficient is not None:
            return sign * coefficient
    if body.startswith("T*"):
        coefficient = _constant_with_r_value(body[2:])
        if coefficient is not None:
            return sign * coefficient
    return None


def _constant_with_r_value(value: str) -> float | None:
    """Return numeric value for a simple constant expression, including R."""
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
        elif re.fullmatch(r"(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+\-]?\d+)?", factor):
            product *= float(factor)
        else:
            return None
    return product


def _signed_number_term(value: float, suffix: str = "") -> str:
    sign = "+" if value >= 0 else "-"
    return f"{sign}{_format_gibbs_number(abs(value))}{suffix}"


def _format_gibbs_number(value: float) -> str:
    return f"{value:.15g}"


def _normalize_signed_term(term: str) -> str:
    return re.sub(r"\s+", "", term)


def _active_gibbs_range(
    function: FunctionDefinition,
    lower: float,
    upper: float,
) -> GibbsRange | None:
    """Return the Gibbs range active over an interval, extrapolating at edges."""
    ranges = function.gibbs_ranges or _parse_gibbs_ranges(function.expression)
    if not ranges:
        return None
    midpoint = (lower + upper) / 2.0
    for gibbs_range in ranges:
        if gibbs_range.t_min - 1e-8 <= midpoint <= gibbs_range.t_max + 1e-8:
            return gibbs_range
    if midpoint < ranges[0].t_min:
        return ranges[0]
    return ranges[-1]


def _dependency_temperature_bounds(
    dependency_terms: list[tuple[str, float]],
    functions_by_name: dict[str, FunctionDefinition],
) -> list[float]:
    """Return all range boundaries from referenced functions."""
    bounds: list[float] = []
    for name, _factor in dependency_terms:
        referenced = functions_by_name.get(name)
        if referenced is None:
            continue
        if referenced.gibbs_ranges:
            for gibbs_range in referenced.gibbs_ranges:
                bounds.extend([gibbs_range.t_min, gibbs_range.t_max])
        else:
            for lower, upper in referenced.temperature_ranges:
                bounds.extend([float(lower), float(upper)])
    return _merged_bounds(bounds)


def _merged_bounds(bounds: list[float]) -> list[float]:
    """Return sorted temperature bounds with near-duplicates merged."""
    merged: list[float] = []
    for bound in sorted(float(value) for value in bounds):
        if not merged or abs(bound - merged[-1]) > 1e-8:
            merged.append(bound)
    return merged


def _canonicalize_element_references(database: DatabaseIR) -> None:
    """Use conventional element symbols throughout element-level IR fields."""
    element_symbols = {
        element.symbol.upper(): element.symbol for element in database.elements
    }
    if not element_symbols:
        return

    for species in database.species:
        species.name = _canonical_element_token(species.name, element_symbols)
        species.composition = {
            _canonical_element_token(symbol, element_symbols): amount
            for symbol, amount in species.composition.items()
        }

    for phase in database.phases:
        for constituent_set in phase.constituents:
            constituent_set.species = [
                _canonical_element_token(name, element_symbols)
                for name in constituent_set.species
            ]

    for parameter in database.parameters:
        parameter.target = [
            _canonicalize_parameter_target(target, element_symbols)
            for target in parameter.target
        ]


def _canonicalize_parameter_target(
    target: str,
    element_symbols: dict[str, str],
) -> str:
    """Normalize element symbols inside a parameter target string."""
    if ":" not in target:
        return _canonical_element_token(target, element_symbols)
    sections = target.split(":")
    return ":".join(
        ",".join(
            _canonical_element_token(token.strip(), element_symbols)
            for token in section.split(",")
        )
        for section in sections
    )


def _canonical_element_token(
    value: str,
    element_symbols: dict[str, str],
) -> str:
    """Return the conventional spelling for pure element tokens."""
    token = str(value).strip()
    return element_symbols.get(token.upper(), token)


def _infer_missing_species(database: DatabaseIR) -> None:
    species_by_name = {species.name: species for species in database.species}
    element_names = {element.symbol for element in database.elements}

    for element in database.elements:
        if element.symbol not in species_by_name:
            composition = {} if element.symbol == "VA" else {element.symbol: 1.0}
            species_by_name[element.symbol] = Species(
                name=element.symbol,
                composition=composition,
                source=element.source,
            )

    for phase in database.phases:
        for constituent_set in phase.constituents:
            for name in constituent_set.species:
                if name in _PSEUDO_SPECIES or name in species_by_name:
                    continue
                composition = (
                    {name: 1.0} if name in element_names else _formula_composition(name)
                )
                species_by_name[name] = Species(name=name, composition=composition)

    database.species = list(species_by_name.values())


def _classify_phase_models(database: DatabaseIR) -> None:
    parameters_by_phase: dict[str, list[Parameter]] = {}
    for parameter in database.parameters:
        parameters_by_phase.setdefault(parameter.phase_name, []).append(parameter)

    for phase in database.phases:
        parameters = parameters_by_phase.get(phase.name, [])
        if _looks_like_cef_solution(phase):
            phase.model = "CEF"
        elif _looks_like_solution(phase, parameters):
            phase.model = "SOLUTION"
        else:
            phase.model = "COMPOUND"


def _record_missing_endmember_diagnostics(database: DatabaseIR) -> None:
    """Record warnings for concrete endmembers without explicit G records."""
    phase_by_name = {phase.name: phase for phase in database.phases}
    ordered_dispart_phases = _ordered_dispart_phase_names(database)
    patterns_by_phase: dict[str, list[list[list[str]]]] = {}
    for parameter in database.parameters:
        if parameter.parameter_type.upper() != "G" or int(parameter.order) != 0:
            continue
        phase = phase_by_name.get(parameter.phase_name)
        if phase is None or not parameter.target:
            continue
        pattern = _parameter_endmember_pattern(
            parameter.target[0],
            len(phase.constituents),
        )
        if pattern:
            patterns_by_phase.setdefault(phase.name, []).append(pattern)

    for phase in database.phases:
        species_by_sublattice = _phase_concrete_endmember_species(phase)
        if not species_by_sublattice:
            continue
        total_endmembers = _endmember_space_size(species_by_sublattice)
        patterns = patterns_by_phase.get(phase.name, [])
        if total_endmembers <= _EXACT_ENDMEMBER_DIAGNOSTIC_LIMIT:
            missing = _exact_missing_endmembers(
                species_by_sublattice,
                patterns,
                phase,
            )
        else:
            missing = _missing_endmember_examples(
                species_by_sublattice,
                patterns,
                _ENDMEMBER_DIAGNOSTIC_PREVIEW_LIMIT,
                phase,
            )
        if not missing:
            continue
        is_ordered_dispart = phase.name.upper() in ordered_dispart_phases
        preview = ", ".join(
            _format_endmember_target(sections)
            for sections in missing[:_ENDMEMBER_DIAGNOSTIC_PREVIEW_LIMIT]
        )
        if len(missing) > _ENDMEMBER_DIAGNOSTIC_PREVIEW_LIMIT:
            preview += ", ..."
        phase_detail = (
            " is a DIS_PART ordered phase and"
            if is_ordered_dispart
            else ""
        )
        if total_endmembers <= _EXACT_ENDMEMBER_DIAGNOSTIC_LIMIT:
            message = (
                f"Phase '{phase.name}'{phase_detail} is missing "
                f"{len(missing)} G endmember parameter(s); runtime calculations "
                "will use implicit zero Gibbs compatibility record(s): "
                f"{preview}."
            )
        else:
            message = (
                f"Phase '{phase.name}'{phase_detail} is missing G endmember "
                "parameter coverage; runtime calculations will use implicit "
                f"zero Gibbs compatibility records. Examples: {preview}."
            )
        database.diagnostics.append(
            Diagnostic(
                "warning",
                message,
                phase.source,
            )
        )


def _ordered_dispart_phase_names(database: DatabaseIR) -> set[str]:
    """Return phase names that are ordered parts in DIS_PART definitions."""
    ordered: set[str] = set()
    phase_names = {phase.name.upper() for phase in database.phases}
    for command in database.tdb_commands:
        if not command.active:
            continue
        parsed = command.parsed
        if parsed.get("action", "").upper() != "DIS_PART":
            continue
        phase_name = _phase_name(parsed.get("ordered_phase", ""))
        if phase_name.upper() in phase_names:
            ordered.add(phase_name.upper())
    return ordered


def _phase_concrete_endmember_species(phase: Phase) -> list[list[str]]:
    """Return concrete species choices for each phase sublattice."""
    if not phase.constituents:
        return []
    sublattice_species: list[list[str]] = []
    for constituent_set in sorted(phase.constituents, key=lambda item: item.sublattice):
        species = [
            str(name).strip()
            for name in constituent_set.species
            if str(name).strip() and str(name).strip() != "*"
        ]
        if not species:
            return []
        sublattice_species.append(species)
    return sublattice_species


def _phase_all_vacancy_endmember_sections(
    species_by_sublattice: list[list[str]],
) -> list[str]:
    """Return the all-vacancy endmember when the constitution allows it."""
    sections: list[str] = []
    for species_choices in species_by_sublattice:
        vacancy = next(
            (
                species
                for species in species_choices
                if _normalize_name(species) == "VA"
            ),
            None,
        )
        if vacancy is None:
            return []
        sections.append(vacancy)
    return sections


def _endmember_space_size(species_by_sublattice: list[list[str]]) -> int:
    total = 1
    for species in species_by_sublattice:
        total *= len(species)
    return total


def _exact_missing_endmembers(
    species_by_sublattice: list[list[str]],
    patterns: list[list[list[str]]],
    phase: Phase | None = None,
) -> list[list[str]]:
    """Return all missing endmembers for a small phase constitution."""
    symmetry_groups = _phase_equivalent_sublattice_groups(phase)
    return [
        list(sections)
        for sections in product(*species_by_sublattice)
        if not any(
            _endmember_pattern_covers_symmetric(
                pattern,
                list(sections),
                symmetry_groups,
            )
            for pattern in patterns
        )
    ]


def _missing_endmember_examples(
    species_by_sublattice: list[list[str]],
    patterns: list[list[list[str]]],
    limit: int,
    phase: Phase | None = None,
) -> list[list[str]]:
    """Return representative missing endmembers without full enumeration."""
    examples: list[list[str]] = []
    symmetry_groups = _phase_equivalent_sublattice_groups(phase)

    def search(
        depth: int,
        prefix: list[str],
        active_patterns: list[list[list[str]]],
    ) -> None:
        if len(examples) >= limit:
            return
        if any(
            _pattern_covers_remaining_space(pattern, species_by_sublattice, depth)
            for pattern in active_patterns
        ):
            return
        if depth >= len(species_by_sublattice):
            if not active_patterns:
                examples.append(list(prefix))
            return

        for species in species_by_sublattice[depth]:
            next_active = [
                pattern
                for pattern in active_patterns
                if _pattern_allows_species_symmetric(
                    pattern,
                    depth,
                    species,
                    prefix,
                    symmetry_groups,
                )
            ]
            next_prefix = [*prefix, species]
            if not next_active:
                examples.append(
                    [
                        *next_prefix,
                        *(choices[0] for choices in species_by_sublattice[depth + 1 :]),
                    ]
                )
                if len(examples) >= limit:
                    return
                continue
            search(depth + 1, next_prefix, next_active)
            if len(examples) >= limit:
                return

    search(0, [], patterns)
    return examples


def _endmember_pattern_covers_symmetric(
    pattern: list[list[str]],
    sections: list[str],
    symmetry_groups: list[list[int]],
) -> bool:
    """Return whether a pattern covers sections, allowing ``:F`` permutations."""
    if _endmember_pattern_covers(pattern, sections):
        return True
    for permuted in _symmetry_equivalent_sections(sections, symmetry_groups):
        if _endmember_pattern_covers(pattern, permuted):
            return True
    return False


def _pattern_allows_species_symmetric(
    pattern: list[list[str]],
    depth: int,
    species: str,
    prefix: list[str],
    symmetry_groups: list[list[int]],
) -> bool:
    """Return whether a partial search branch can still match a pattern."""
    if _pattern_allows_species(pattern[depth], species):
        return True
    if not symmetry_groups:
        return False
    sections = [*prefix, species]
    for group in symmetry_groups:
        if depth not in group:
            continue
        known_group_values = [
            sections[index]
            for index in group
            if index < len(sections)
        ]
        if any(
            _pattern_allows_species(pattern[index], value)
            for index in group
            for value in known_group_values
        ):
            return True
    return False


def _symmetry_equivalent_sections(
    sections: list[str],
    symmetry_groups: list[list[int]],
) -> list[list[str]]:
    """Return representative section permutations for ``:F``-equivalent sites."""
    variants = [list(sections)]
    for group in symmetry_groups:
        if any(index >= len(sections) for index in group):
            continue
        next_variants: list[list[str]] = []
        for variant in variants:
            seen_permutations: set[tuple[str, ...]] = set()
            for permutation in permutations(
                [variant[index] for index in group],
                len(group),
            ):
                if permutation in seen_permutations:
                    continue
                seen_permutations.add(permutation)
                candidate = list(variant)
                for index, value in zip(group, permutation, strict=True):
                    candidate[index] = value
                next_variants.append(candidate)
        variants = next_variants
    unique: list[list[str]] = []
    seen: set[tuple[str, ...]] = set()
    for variant in variants:
        key = tuple(variant)
        if key in seen:
            continue
        seen.add(key)
        unique.append(variant)
    return unique


def _phase_equivalent_sublattice_groups(phase: Phase | None) -> list[list[int]]:
    """Return zero-based sublattice groups equivalent under a ``PHASE ...:F``."""
    if phase is None or not _phase_has_f_option(phase):
        return []
    groups_by_key: dict[tuple[float, tuple[str, ...]], list[int]] = {}
    for index, constituent in enumerate(
        sorted(phase.constituents, key=lambda item: item.sublattice)
    ):
        species_key = tuple(_normalize_name(species) for species in constituent.species)
        key = (float(constituent.site_ratio), species_key)
        groups_by_key.setdefault(key, []).append(index)
    return [indices for indices in groups_by_key.values() if len(indices) > 1]


def _phase_has_f_option(phase: Phase) -> bool:
    source_name = _source_phase_name(phase.source.command)
    if ":" not in source_name:
        return False
    return "F" in source_name.split(":", 1)[1].upper()


def _source_phase_name(command: str) -> str:
    """Return the original PHASE name token when available."""
    match = re.search(r"(?i)(?:^|\s)PHASE\s+(\S+)", str(command or ""))
    return match.group(1).strip() if match is not None else ""


def _pattern_covers_remaining_space(
    pattern: list[list[str]],
    species_by_sublattice: list[list[str]],
    depth: int,
) -> bool:
    """Return whether one active pattern covers all remaining combinations."""
    for index in range(depth, len(species_by_sublattice)):
        allowed = pattern[index]
        if "*" in allowed:
            continue
        species = {_normalize_name(name) for name in species_by_sublattice[index]}
        if not species.issubset(set(allowed)):
            return False
    return True


def _parameter_endmember_pattern(
    target: str,
    sublattice_count: int,
) -> list[list[str]] | None:
    """Return normalized allowed species per target sublattice.

    Only concrete endmember targets and single-site wildcards can cover
    endmembers.  CEF ``G`` records with comma lists, such as
    ``G(FCC_4SL,MN,NI:MN,NI:*:*:VA;0)``, are interaction/order parameters and
    must not satisfy missing-endmember checks.
    """
    sections = str(target).split(":")
    if len(sections) != sublattice_count:
        return None
    pattern: list[list[str]] = []
    for section in sections:
        if "," in section:
            return None
        tokens = [
            _normalize_name(token) for token in section.split(",") if token.strip()
        ]
        if not tokens:
            return None
        pattern.append(tokens)
    return pattern


def _endmember_pattern_covers(
    pattern: list[list[str]],
    endmember_sections: list[str],
) -> bool:
    """Return whether a parameter target pattern covers one endmember."""
    if len(pattern) != len(endmember_sections):
        return False
    for allowed, species in zip(pattern, endmember_sections, strict=True):
        if not _pattern_allows_species(allowed, species):
            return False
    return True


def _pattern_allows_species(allowed: list[str], species: str) -> bool:
    return "*" in allowed or _normalize_name(species) in allowed


def _format_endmember_target(sections: list[str]) -> str:
    return ":".join(str(section).strip() for section in sections)


def _looks_like_solution(phase: Phase, parameters: list[Parameter]) -> bool:
    if phase.name in _COMMON_SOLUTION_PHASES:
        return True
    if _looks_like_cef_solution(phase):
        return True
    if _has_mixing_constituents(phase):
        return True
    if any("," in target for parameter in parameters for target in parameter.target):
        return True
    return False


def _has_mixing_constituents(phase: Phase) -> bool:
    return any(
        len([species for species in constituent_set.species if str(species).strip()])
        > 1
        for constituent_set in phase.constituents
    )


def _looks_like_cef_solution(phase: Phase) -> bool:
    """Return whether a multi-sublattice phase allows sublattice mixing."""
    if len(phase.constituents) <= 1:
        return False
    return any(
        len([species for species in constituent_set.species if str(species).strip()])
        > 1
        for constituent_set in phase.constituents
    )


def _sort_ir(database: DatabaseIR) -> None:
    database.elements.sort(key=lambda item: item.symbol)
    database.species.sort(key=lambda item: item.name)
    database.functions.sort(key=lambda item: item.name)
    database.phases.sort(key=lambda item: item.name)


def _phase_name(raw_name: str) -> str:
    return _normalize_name(raw_name.split(":", 1)[0])


def _normalize_name(name: str) -> str:
    return name.strip().strip("'").replace("%", "").upper()


def _normalize_element_symbol(symbol: str) -> str:
    text = str(symbol).strip().strip("'").replace("%", "")
    upper = text.upper()
    if upper in {"VA", "/-", "*"}:
        return upper
    if upper.isalpha() and len(upper) <= 2:
        return upper[:1] + upper[1:].lower()
    return text


def _to_float(value: str, default: float) -> float:
    try:
        return float(value.replace("D", "E").replace("d", "e"))
    except ValueError:
        return default


def _formula_composition(name: str) -> dict[str, float]:
    """Best-effort composition for formula-like species names."""
    cleaned = re.sub(r"[^A-Za-z0-9]", "", name)
    matches = re.findall(r"([A-Z][a-z]?)([0-9.]*)", cleaned.title())
    composition: dict[str, float] = {}
    for element, amount_text in matches:
        amount = float(amount_text) if amount_text else 1.0
        composition[element.upper()] = composition.get(element.upper(), 0.0) + amount
    return composition
