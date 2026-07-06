"""Qt model/view wrappers for DatabaseIR."""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from PySide6.QtCore import (
    QAbstractItemModel,
    QAbstractTableModel,
    QByteArray,
    QMimeData,
    QModelIndex,
    Qt,
)
from PySide6.QtGui import QIcon

from equilipy.database_ir import (
    DatabaseIR,
    Diagnostic,
    FunctionDefinition,
    Parameter,
    Phase,
)
from equilipy.gui.assets import status_icon_path

_INVALID_INDEX = QModelIndex()
_DATABASE_TREE_DRAG_MIME = "application/x-equilipy-database-tree"
_MAX_CP_RANGES = 6
_FUNCTION_TOKEN_RE = re.compile(r"(?<![A-Za-z0-9_])([A-Za-z_][A-Za-z0-9_]*)")
_SUBSCRIPT_TRANSLATION = str.maketrans(
    {
        "0": "₀",
        "1": "₁",
        "2": "₂",
        "3": "₃",
        "4": "₄",
        "5": "₅",
        "6": "₆",
        "7": "₇",
        "8": "₈",
        "9": "₉",
        "+": "₊",
        "-": "₋",
    }
)
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
class _TreeBuildContext:
    """Per-database caches used while building the database tree."""

    functions_by_name: dict[str, FunctionDefinition]
    transition_cache: dict[str, set[float]]


@dataclass(frozen=True, slots=True)
class ThermoRangePayload:
    """GUI payload for one Cp temperature range under a function or compound."""

    owner_name: str
    owner_kind: str
    property_name: str
    index: int
    expression: str
    lower: float | str = ""
    upper: float | str = ""
    source: Any = None


@dataclass(frozen=True, slots=True)
class CompoundRecord:
    """GUI grouping for one compound with up to ten phase states."""

    name: str
    phases: tuple[Phase, ...]


@dataclass(frozen=True, slots=True)
class CompoundPhasePayload:
    """GUI payload for one compound phase state such as S1 or L."""

    compound_name: str
    phase_label: str
    phase: Phase
    index: int


@dataclass(frozen=True, slots=True)
class SolutionGroupPayload:
    """GUI grouping payload under one solution phase."""

    phase: Phase
    group_name: str


@dataclass(frozen=True, slots=True)
class SolutionSublatticePayload:
    """GUI payload for one solution sublattice."""

    phase: Phase
    sublattice: int
    label: str
    site_ratio: float
    species: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class SolutionSpeciesPayload:
    """GUI payload for one species within a solution sublattice."""

    phase: Phase
    sublattice: int
    sublattice_label: str
    species_name: str
    site_ratio: float


@dataclass(frozen=True, slots=True)
class SolutionParameterPayload:
    """GUI payload for one solution endmember or interaction parameter."""

    phase: Phase
    parameter: Parameter
    kind: str
    index: int


class TreeNode:
    """Simple tree node used by DatabaseTreeModel."""

    def __init__(
        self,
        label: str,
        payload: Any = None,
        parent: "TreeNode | None" = None,
    ):
        self.label = label
        self.payload = payload
        self.parent = parent
        self.children: list[TreeNode] = []

    def append(self, child: "TreeNode") -> None:
        """Append a child node."""
        child.parent = self
        self.children.append(child)

    def row(self) -> int:
        """Return this node's row within its parent."""
        if self.parent is None:
            return 0
        return self.parent.children.index(self)


class DatabaseTreeModel(QAbstractItemModel):
    """Tree model for browsing loaded databases by workflow category."""

    def __init__(
        self,
        database: DatabaseIR | list[DatabaseIR],
        *,
        function_warnings: dict[int, list[str]] | None = None,
        function_errors: dict[int, list[str]] | None = None,
        range_warnings: dict[str, list[str]] | None = None,
    ):
        super().__init__()
        self.databases = database if isinstance(database, list) else [database]
        self.database = self.databases[0]
        self.function_warnings = {} if function_warnings is None else function_warnings
        self.function_errors = {} if function_errors is None else function_errors
        self.range_warnings = {} if range_warnings is None else range_warnings
        self.root = self._build_tree(self.databases)

    def _build_tree(self, databases: list[DatabaseIR]) -> TreeNode:
        root = TreeNode("Databases")
        for database in databases:
            if not _database_has_content(database) and not database.name.strip():
                continue
            context = _tree_build_context(database)
            database_node = TreeNode(_database_tree_label(database), database)
            root.append(database_node)

            compounds, solutions = _split_phase_groups(database.phases)
            compound_records = _compound_records(compounds)
            groups = [
                ("Functions", database.functions, "Function"),
                ("Compounds", compound_records, "Compound"),
                ("Solutions", solutions, "Solution"),
            ]

            for group_name, items, prefix in groups:
                group = TreeNode(group_name, items)
                database_node.append(group)
                for index, item in enumerate(items, start=1):
                    object_node = TreeNode(_tree_item_label(item, prefix, index), item)
                    group.append(object_node)
                    if isinstance(item, CompoundRecord):
                        for phase_index, phase in enumerate(item.phases[:10], start=1):
                            phase_payload = CompoundPhasePayload(
                                compound_name=item.name,
                                phase_label=_compound_phase_label(phase, phase_index),
                                phase=phase,
                                index=phase_index,
                            )
                            phase_node = TreeNode(
                                phase_payload.phase_label,
                                phase_payload,
                            )
                            object_node.append(phase_node)
                            for range_payload in _cp_range_payloads(
                                database,
                                phase_payload,
                                context,
                            ):
                                phase_node.append(
                                    TreeNode(
                                        _cp_range_label(range_payload),
                                        range_payload,
                                    )
                                )
                        continue
                    if isinstance(item, Phase) and _is_solution_phase(item):
                        for child in _solution_phase_children(database, item, context):
                            object_node.append(child)
                        continue
                    for range_payload in _cp_range_payloads(database, item, context):
                        object_node.append(
                            TreeNode(_cp_range_label(range_payload), range_payload)
                        )
        return root

    def columnCount(self, parent: QModelIndex = _INVALID_INDEX) -> int:
        """Return column count."""
        return 1

    def rowCount(self, parent: QModelIndex = _INVALID_INDEX) -> int:
        """Return child count for parent."""
        node = self._node_from_index(parent)
        return len(node.children)

    def index(
        self, row: int, column: int, parent: QModelIndex = _INVALID_INDEX
    ) -> QModelIndex:
        """Return index for a child node."""
        if not self.hasIndex(row, column, parent):
            return QModelIndex()
        parent_node = self._node_from_index(parent)
        child = parent_node.children[row]
        return self.createIndex(row, column, child)

    def parent(self, index: QModelIndex) -> QModelIndex:
        """Return parent index for a node."""
        if not index.isValid():
            return QModelIndex()
        node = self._node_from_index(index)
        parent = node.parent
        if parent is None or parent is self.root:
            return QModelIndex()
        return self.createIndex(parent.row(), 0, parent)

    def data(self, index: QModelIndex, role: int = Qt.ItemDataRole.DisplayRole) -> Any:
        """Return display data."""
        if not index.isValid():
            return None
        node = self._node_from_index(index)
        if role == Qt.ItemDataRole.DisplayRole:
            return node.label
        if role == Qt.ItemDataRole.DecorationRole and self._has_error(node):
            return function_error_icon()
        if role == Qt.ItemDataRole.DecorationRole and self._has_warning(node):
            return function_warning_icon()
        if role == Qt.ItemDataRole.ToolTipRole and self._has_notice(node):
            parts: list[str] = []
            errors = self.function_errors.get(id(node.payload), [])
            warnings = self._warning_messages(node)
            if errors:
                parts.append("\n".join(errors))
            if warnings:
                parts.append("\n".join(warnings))
            return "\n\n".join(parts)
        return None

    def headerData(
        self,
        section: int,
        orientation: Qt.Orientation,
        role: int = Qt.ItemDataRole.DisplayRole,
    ) -> Any:
        """Return tree header."""
        if (
            orientation == Qt.Orientation.Horizontal
            and role == Qt.ItemDataRole.DisplayRole
        ):
            return "Database"
        return None

    def flags(self, index: QModelIndex) -> Qt.ItemFlag:
        """Return tree item flags, including copy-only drag/drop support."""
        flags = super().flags(index)
        if not index.isValid():
            return flags
        payload = self.payload(index)
        if _is_draggable_payload(payload):
            flags |= Qt.ItemFlag.ItemIsDragEnabled
        flags |= Qt.ItemFlag.ItemIsDropEnabled
        return flags

    def supportedDragActions(self) -> Qt.DropAction:
        """Drag copies database objects; it never moves them."""
        return Qt.DropAction.CopyAction

    def supportedDropActions(self) -> Qt.DropAction:
        """Drop copies database objects; it never moves them."""
        return Qt.DropAction.CopyAction

    def mimeTypes(self) -> list[str]:
        """Return MIME types for database tree drags."""
        return [_DATABASE_TREE_DRAG_MIME, "text/plain"]

    def mimeData(self, indexes: list[QModelIndex]) -> QMimeData:
        """Attach internal drag data, plus Function text for expression editors."""
        mime_data = QMimeData()
        for index in indexes:
            if not index.isValid() or index.column() != 0:
                continue
            payload = self.payload(index)
            if not _is_draggable_payload(payload):
                continue
            mime_data.setData(_DATABASE_TREE_DRAG_MIME, QByteArray(b"copy"))
            if isinstance(payload, FunctionDefinition):
                mime_data.setText(f"{payload.name}#")
            break
        return mime_data

    def payload(self, index: QModelIndex) -> Any:
        """Return payload for the selected node."""
        if index.isValid():
            return self._node_from_index(index).payload
        return self.database

    def label(self, index: QModelIndex) -> str:
        """Return display label for the selected node."""
        if index.isValid():
            return self._node_from_index(index).label
        return self.database.name

    def database_for(self, index: QModelIndex) -> DatabaseIR:
        """Return the loaded database that owns the selected node."""
        node = self._node_from_index(index)
        while node.parent is not None:
            if isinstance(node.payload, DatabaseIR):
                return node.payload
            node = node.parent
        return self.database

    def _node_from_index(self, index: QModelIndex) -> TreeNode:
        if index.isValid():
            return index.internalPointer()
        return self.root

    def _has_warning(self, node: TreeNode) -> bool:
        return bool(self._warning_messages(node))

    def _local_warning_messages(self, node: TreeNode) -> list[str]:
        """Return warnings attached directly to one tree node."""
        warnings = list(self.function_warnings.get(id(node.payload), []))
        key = thermo_range_warning_key(node.payload)
        if key:
            warnings.extend(self.range_warnings.get(key, []))
        return warnings

    def _warning_messages(self, node: TreeNode) -> list[str]:
        """Return warnings for one node, including relevant child ranges."""
        warnings = self._local_warning_messages(node)
        if self._bubbles_range_warnings(node):
            for child in node.children:
                warnings.extend(self._warning_messages(child))
        return _unique_messages(warnings)

    def _bubbles_range_warnings(self, node: TreeNode) -> bool:
        """Return whether child range warnings should decorate this node."""
        return isinstance(
            node.payload,
            (FunctionDefinition, CompoundRecord, CompoundPhasePayload),
        )

    def _has_error(self, node: TreeNode) -> bool:
        return (
            isinstance(node.payload, FunctionDefinition)
            and bool(self.function_errors.get(id(node.payload)))
        )

    def _has_notice(self, node: TreeNode) -> bool:
        return self._has_error(node) or self._has_warning(node)


def function_error_icon(size: int = 22) -> QIcon:
    """Return a red triangle error icon for database tree validation errors."""
    return _function_notice_icon(size, "error")


def function_warning_icon(size: int = 22) -> QIcon:
    """Return a yellow triangle warning icon for database tree validation warnings."""
    return _function_notice_icon(size, "warning")


def _function_notice_icon(_size: int, name: str) -> QIcon:
    """Return a vector triangle validation icon with consistent shape."""
    return QIcon(str(status_icon_path(name)))


def _is_draggable_payload(payload: Any) -> bool:
    """Return whether one database tree payload can be copied by dragging."""
    return isinstance(
        payload,
        (
            FunctionDefinition,
            CompoundRecord,
            CompoundPhasePayload,
            SolutionParameterPayload,
            ThermoRangePayload,
        ),
    ) or (isinstance(payload, Phase) and _is_solution_phase(payload))


def thermo_range_warning_key(payload: Any) -> str | None:
    """Return a stable warning key for a tree thermo range payload."""
    if not isinstance(payload, ThermoRangePayload):
        return None
    source_command = getattr(payload.source, "command", "")
    return "|".join(
        (
            payload.owner_kind,
            payload.owner_name,
            payload.property_name,
            str(payload.lower),
            str(payload.upper),
            str(source_command),
        )
    )


def _unique_messages(messages: list[str]) -> list[str]:
    """Return messages with duplicates removed while preserving order."""
    seen: set[str] = set()
    unique: list[str] = []
    for message in messages:
        if message in seen:
            continue
        seen.add(message)
        unique.append(message)
    return unique


def _split_phase_groups(phases: list[Phase]) -> tuple[list[Phase], list[Phase]]:
    """Split phases into temporary compound and solution GUI groups."""
    solutions = [phase for phase in phases if _is_solution_phase(phase)]
    compounds = [phase for phase in phases if phase not in solutions]

    if not compounds and not solutions and phases:
        compounds = [phases[0]]
        solutions = [phase for phase in phases[1:] if phase.model.upper()]

    if not solutions and len(phases) > 1:
        solutions = phases[1:]

    return compounds, solutions


def _database_has_content(database: DatabaseIR) -> bool:
    """Return whether a DatabaseIR should be displayed in the tree."""
    return any(
        (
            database.elements,
            database.species,
            database.functions,
            database.phases,
            database.parameters,
            database.diagnostics,
        )
    )


def _database_tree_label(database: DatabaseIR) -> str:
    """Return a database root label that includes the file extension when known."""
    loaded_path = str(
        database.metadata.get("_loaded_path")
        or database.metadata.get("source_file")
        or ""
    ).strip()
    path = Path(loaded_path) if loaded_path else None
    suffix = path.suffix if path is not None else ""
    name = str(database.name).strip()
    if not name and path is not None:
        return path.name
    if suffix and not name.lower().endswith(suffix.lower()):
        return f"{name}{suffix}"
    return name


def _tree_build_context(database: DatabaseIR) -> _TreeBuildContext:
    """Return lookup caches shared across one tree-model build."""
    return _TreeBuildContext(
        functions_by_name={
            function.name.upper(): function
            for function in database.functions
        },
        transition_cache={},
    )


def _is_solution_phase(phase: Phase) -> bool:
    """Return whether a phase should use the solution-phase editor."""
    return phase.model.upper() in _SOLUTION_MODELS


def _tree_item_label(item: Any, prefix: str, index: int) -> str:
    """Return a useful tree label for one database object."""
    if isinstance(item, CompoundRecord):
        return item.name
    if isinstance(item, CompoundPhasePayload):
        return item.phase_label
    for attribute in ("name", "symbol", "phase_name", "severity"):
        value = getattr(item, attribute, "")
        if str(value).strip():
            return str(value).strip()
    return f"{prefix}#{index}"


def _cp_range_payloads(
    database: DatabaseIR,
    item: Any,
    context: _TreeBuildContext,
) -> list[ThermoRangePayload]:
    """Return up to six Cp range child payloads for a function or compound."""
    if isinstance(item, FunctionDefinition):
        return _function_cp_range_payloads(item, context)
    if isinstance(item, CompoundPhasePayload):
        return _compound_cp_range_payloads(
            database,
            item.phase,
            context,
            owner_name=f"{item.compound_name}:{item.phase_label}",
        )
    if isinstance(item, Phase):
        return _compound_cp_range_payloads(database, item, context)
    return []


def _solution_phase_children(
    database: DatabaseIR,
    phase: Phase,
    context: _TreeBuildContext,
) -> list[TreeNode]:
    """Return structured tree children for one solution phase."""
    parameters = _solution_phase_parameters(database, phase)
    endmembers = [
        parameter
        for parameter in parameters
        if not _is_solution_interaction_parameter(parameter)
    ]
    interactions = [
        parameter
        for parameter in parameters
        if _is_solution_interaction_parameter(parameter)
    ]

    endmember_group = TreeNode(
        f"Endmembers ({len(endmembers)})",
        SolutionGroupPayload(phase, "Endmembers"),
    )
    for index, parameter in enumerate(endmembers):
        endmember_node = TreeNode(
            _solution_endmember_label(phase, parameter, index),
            SolutionParameterPayload(phase, parameter, "Endmember", index),
        )
        for range_payload in _solution_endmember_cp_range_payloads(
            phase,
            parameter,
            context,
        ):
            endmember_node.append(
                TreeNode(_cp_range_label(range_payload), range_payload)
            )
        endmember_group.append(endmember_node)

    interaction_group = TreeNode(
        f"Interactions ({len(interactions)})",
        SolutionGroupPayload(phase, "Interactions"),
    )
    for index, parameter in enumerate(interactions):
        interaction_group.append(
            TreeNode(
                f"({index}) {_parameter_target_text(parameter)}",
                SolutionParameterPayload(phase, parameter, "Interaction", index),
            )
        )

    return [endmember_group, interaction_group]


def _solution_endmember_cp_range_payloads(
    phase: Phase,
    parameter: Parameter,
    context: _TreeBuildContext,
) -> list[ThermoRangePayload]:
    """Return Cp child ranges for one solution endmember parameter."""
    payloads: list[ThermoRangePayload] = []
    for lower, upper in _compound_parameter_cp_ranges(parameter, context):
        payloads.append(
            ThermoRangePayload(
                owner_name=f"{phase.name}:{_parameter_target_text(parameter)}",
                owner_kind="SolutionEndmember",
                property_name="Cp",
                index=len(payloads) + 1,
                expression=parameter.expression,
                lower=lower,
                upper=upper,
                source=parameter.source,
            )
        )
        if len(payloads) >= _MAX_CP_RANGES:
            break
    return payloads


def _solution_sublattice_label(phase: Phase, index: int, species_count: int) -> str:
    """Return FactSage-like sublattice labels, A-Z then numbered fallback."""
    labels = phase.metadata.get("sublattice_labels")
    if isinstance(labels, list) and 0 <= index - 1 < len(labels):
        name = str(labels[index - 1]).strip()
        if name:
            return f"{name} ({species_count})"
    if 1 <= index <= 26:
        return f"{chr(ord('A') + index - 1)} ({species_count})"
    return f"S{index} ({species_count})"


def _solution_phase_parameters(
    database: DatabaseIR,
    phase: Phase,
) -> list[Parameter]:
    """Return parameters attached to one solution phase, preserving file order."""
    phase_name = phase.name.upper()
    return [
        parameter
        for parameter in database.parameters
        if parameter.phase_name.upper() == phase_name
    ]


def _is_solution_interaction_parameter(parameter: Parameter) -> bool:
    """Return whether a solution parameter describes mixing/interactions."""
    if str(parameter.parameter_type).strip().upper() != "G":
        return True
    return any("," in section for section in _parameter_target_sections(parameter))


def _parameter_target_sections(parameter: Parameter) -> list[str]:
    """Return target sections split by sublattice, not by mixed species."""
    target = list(parameter.target)
    if len(target) == 1 and ":" in target[0]:
        return target[0].split(":")
    return target


def _parameter_target_text(parameter: Parameter) -> str:
    """Return a compact target string for tree labels and tables."""
    return ":".join(_parameter_target_sections(parameter)) or "-"


def _solution_endmember_label(
    phase: Phase,
    parameter: Parameter,
    index: int,
) -> str:
    """Return an endmember tree label with sublattice stoichiometry."""
    sections = _parameter_target_sections(parameter)
    if not sections:
        return f"({index}) -"
    constituents = sorted(phase.constituents, key=lambda item: item.sublattice)
    pieces: list[str] = []
    for section_index, section in enumerate(sections):
        species = _solution_species_display_name(phase, section)
        if section_index < len(constituents):
            coefficient = _subscript_text(
                _compact_number(constituents[section_index].site_ratio)
            )
            pieces.append(f"{species}{coefficient}")
        else:
            pieces.append(species)
    return f"({index}) {''.join(pieces)}"


def _solution_species_display_name(phase: Phase, name: str) -> str:
    """Return a target species name using the phase's displayed casing."""
    normalized_name = str(name).strip()
    if normalized_name.upper() == "VA":
        return "Va"
    for constituent in phase.constituents:
        for species_name in constituent.species:
            if str(species_name).upper() == normalized_name.upper():
                return str(species_name)
    if len(normalized_name) <= 2 and normalized_name.isalpha():
        return normalized_name.title()
    return normalized_name


def _subscript_text(text: str) -> str:
    """Return text with digits rendered as unicode subscripts."""
    return str(text).translate(_SUBSCRIPT_TRANSLATION)


def _function_cp_range_payloads(
    function: FunctionDefinition,
    context: _TreeBuildContext,
) -> list[ThermoRangePayload]:
    """Return Cp child ranges for one function definition."""
    ranges = list(function.temperature_ranges[:_MAX_CP_RANGES])
    if not ranges and function.expression.strip():
        ranges = _referenced_function_cp_ranges(function.expression, context)
    if not ranges and function.expression.strip():
        ranges = [("", "")]
    return [
        ThermoRangePayload(
            owner_name=function.name,
            owner_kind="Function",
            property_name="Cp",
            index=index,
            expression=function.expression,
            lower=lower,
            upper=upper,
            source=function.source,
        )
        for index, (lower, upper) in enumerate(ranges, start=1)
    ]


def _referenced_function_cp_ranges(
    expression: str,
    context: _TreeBuildContext,
) -> list[tuple[float, float]]:
    """Return cheap dependency-derived ranges for composite Functions."""
    span = _referenced_function_temperature_span(expression, context, set())
    if span is None:
        return []
    lower, upper, transitions = span
    points = [lower, upper]
    points.extend(point for point in transitions if lower < point < upper)
    sorted_points = sorted(set(points))
    return [
        (sorted_points[index], sorted_points[index + 1])
        for index in range(min(len(sorted_points) - 1, _MAX_CP_RANGES))
    ]


def _referenced_function_temperature_span(
    expression: str,
    context: _TreeBuildContext,
    seen: set[str],
) -> tuple[float, float, set[float]] | None:
    """Return min lower, max upper, and transition points for dependencies."""
    lower_bounds: list[float] = []
    upper_bounds: list[float] = []
    transitions: set[float] = set()
    for function in _referenced_functions_in_expression(expression, context):
        normalized_name = function.name.upper()
        if normalized_name in seen:
            continue
        next_seen = {*seen, normalized_name}
        ranges = function.temperature_ranges or _piecewise_temperature_ranges(
            function.expression,
        )
        if ranges:
            lower_bounds.append(float(ranges[0][0]))
            upper_bounds.append(float(ranges[-1][1]))
            transitions.update(float(lower) for lower, _upper in ranges[1:])
        nested = _referenced_function_temperature_span(
            function.expression,
            context,
            next_seen,
        )
        if nested is not None:
            nested_lower, nested_upper, nested_transitions = nested
            lower_bounds.append(nested_lower)
            upper_bounds.append(nested_upper)
            transitions.update(nested_transitions)
    if not lower_bounds or not upper_bounds:
        return None
    return min(lower_bounds), max(upper_bounds), transitions


def _compound_cp_range_payloads(
    database: DatabaseIR,
    phase: Phase,
    context: _TreeBuildContext,
    owner_name: str | None = None,
) -> list[ThermoRangePayload]:
    """Return Cp child ranges for one compound phase."""
    if _is_solution_phase(phase):
        return []
    payloads: list[ThermoRangePayload] = []
    phase_id = str(getattr(phase, "metadata", {}).get("phase_id", "")).strip()
    g_parameters = [
        parameter
        for parameter in database.parameters
        if parameter.phase_name == phase.name
        and (
            parameter.metadata.get("phase_id") == phase_id
            if phase_id
            else not parameter.metadata.get("phase_id")
        )
        and _thermo_property_key(parameter.parameter_type) == "G"
    ]
    range_parameters = g_parameters or [
        parameter
        for parameter in database.parameters
        if parameter.phase_name == phase.name
        and (
            parameter.metadata.get("phase_id") == phase_id
            if phase_id
            else not parameter.metadata.get("phase_id")
        )
        and _thermo_property_key(parameter.parameter_type) == "Cp"
    ]
    seen_ranges: set[tuple[str, str, str]] = set()
    for parameter in range_parameters:
        property_name = "Cp"
        for lower, upper in _compound_parameter_cp_ranges(
            parameter,
            context,
        ):
            range_key = (property_name, str(lower), str(upper))
            if range_key in seen_ranges:
                continue
            seen_ranges.add(range_key)
            payloads.append(
                ThermoRangePayload(
                    owner_name=owner_name or phase.name,
                    owner_kind="Compound",
                    property_name=property_name,
                    index=len(payloads) + 1,
                    expression=parameter.expression,
                    lower=lower,
                    upper=upper,
                    source=parameter.source,
                )
            )
            if len(payloads) >= _MAX_CP_RANGES:
                break
        if len(payloads) >= _MAX_CP_RANGES:
            break
    return payloads


def _compound_records(phases: list[Phase]) -> list[CompoundRecord]:
    """Return compound GUI records grouped by compound name."""
    records: dict[str, list[Phase]] = {}
    order: list[str] = []
    for phase in phases:
        name = _compound_record_name(phase)
        if name not in records:
            records[name] = []
            order.append(name)
        if len(records[name]) < 10:
            records[name].append(phase)
    return [
        CompoundRecord(name=name, phases=tuple(records[name]))
        for name in order
    ]


def _compound_record_name(phase: Phase) -> str:
    """Return the visible compound name for a phase-state record."""
    return str(phase.name).strip() or "Compound"


def _compound_phase_label(phase: Phase, index: int) -> str:
    """Return default phase-state labels, using L for explicit liquids."""
    label = str(getattr(phase, "metadata", {}).get("phase_label", "")).strip()
    if label:
        return label
    if _is_explicit_liquid_phase(phase):
        return "L"
    return f"S{min(max(index, 1), 10)}"


def _is_explicit_liquid_phase(phase: Phase) -> bool:
    """Return whether a compound phase was explicitly tagged as liquid."""
    command = str(phase.source.command).upper()
    if re.search(r"\bPHASE\s+\S+:L\b", command):
        return True
    return str(phase.name).upper() == "LIQUID"


def _cp_range_label(payload: ThermoRangePayload) -> str:
    """Return a tree label for a Cp range, using the upper bound when known."""
    if payload.upper != "":
        return f"{payload.property_name} {_compact_number(payload.upper)}"
    return f"{payload.property_name} {payload.index}"


def _compact_number(value: Any) -> str:
    """Format a range boundary compactly for tree display."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if number.is_integer():
        return str(int(number))
    return f"{number:g}"


def _parameter_temperature_range(
    parameter: Parameter,
) -> tuple[float | str, float | str]:
    """Best-effort extraction of TDB parameter lower/upper temperatures."""
    ranges = _parameter_temperature_ranges(parameter)
    if ranges:
        return ranges[0][0], ranges[-1][1]
    return "", ""


def _compound_parameter_cp_ranges(
    parameter: Parameter,
    context: _TreeBuildContext,
) -> list[tuple[float | str, float | str]]:
    """Return displayed Cp ranges for a compound parameter."""
    parameter_ranges = _parameter_temperature_ranges(parameter)
    if not parameter_ranges:
        return [("", "")]

    split_points = _referenced_function_transition_points(
        parameter.expression,
        context,
        set(),
    )
    ranges: list[tuple[float | str, float | str]] = []
    for lower, upper in parameter_ranges:
        points = [lower, upper]
        points.extend(point for point in split_points if lower < point < upper)
        sorted_points = sorted(set(points))
        ranges.extend(
            (sorted_points[index], sorted_points[index + 1])
            for index in range(len(sorted_points) - 1)
        )
    return ranges or parameter_ranges


def _parameter_temperature_ranges(
    parameter: Parameter,
) -> list[tuple[float, float]]:
    """Best-effort extraction of all TDB parameter temperature ranges."""
    return _piecewise_temperature_ranges(parameter.expression) or (
        _piecewise_temperature_ranges(parameter.source.command)
    )


def _piecewise_temperature_ranges(text: str) -> list[tuple[float, float]]:
    """Return piecewise temperature ranges encoded in TDB-like text."""
    command = str(text or "")
    lower_match = re.search(r"\)\s*([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)", command)
    if lower_match is None:
        lower_match = re.search(
            r"^\s*([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)\b",
            command,
        )
    upper_matches = re.findall(
        r";\s*([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)\s+[A-Z]\b",
        command,
        flags=re.IGNORECASE,
    )
    if lower_match is None or not upper_matches:
        return []
    ranges: list[tuple[float, float]] = []
    lower = float(lower_match.group(1))
    for upper_text in upper_matches:
        upper = float(upper_text)
        if upper > lower:
            ranges.append((lower, upper))
        lower = upper
    return ranges


def _referenced_function_transition_points(
    expression: str,
    context: _TreeBuildContext,
    seen: set[str],
) -> set[float]:
    """Return transition temperatures introduced by referenced Functions."""
    transitions: set[float] = set()
    for function in _referenced_functions_in_expression(expression, context):
        normalized_name = function.name.upper()
        if normalized_name in seen:
            continue
        transitions.update(
            _function_transition_points(function, context, {*seen, normalized_name})
        )
    return transitions


def _referenced_functions_in_expression(
    expression: str,
    context: _TreeBuildContext,
) -> list[FunctionDefinition]:
    """Return function references found as exact expression tokens."""
    functions: list[FunctionDefinition] = []
    seen: set[str] = set()
    for match in _FUNCTION_TOKEN_RE.finditer(str(expression or "")):
        name = match.group(1).upper()
        if name in seen:
            continue
        function = context.functions_by_name.get(name)
        if function is None:
            continue
        seen.add(name)
        functions.append(function)
    return functions


def _function_transition_points(
    function: FunctionDefinition,
    context: _TreeBuildContext,
    seen: set[str],
) -> set[float]:
    """Return all Cp transition points for a referenced function tree."""
    normalized_name = function.name.upper()
    cached = context.transition_cache.get(normalized_name)
    if cached is not None:
        return set(cached)

    ranges = function.temperature_ranges or _piecewise_temperature_ranges(
        function.expression,
    )
    transitions = {float(lower) for lower, _upper in ranges[1:]}
    for referenced_function in _referenced_functions_in_expression(
        function.expression,
        context,
    ):
        referenced_name = referenced_function.name.upper()
        if referenced_name in seen:
            continue
        transitions.update(
            _function_transition_points(
                referenced_function,
                context,
                {*seen, referenced_name},
            )
        )

    context.transition_cache[normalized_name] = set(transitions)
    return transitions


def _thermo_property_key(parameter_type: str) -> str:
    """Map a parameter type onto the GUI thermodynamic expression tabs."""
    normalized = str(parameter_type).strip().upper().replace(" ", "")
    if normalized in {"CP", "C_P", "HEATCAPACITY", "HEAT_CAPACITY"}:
        return "Cp"
    if normalized in {"H", "ENTHALPY"}:
        return "H"
    if normalized in {"S", "ENTROPY"}:
        return "S"
    return "G"


class DiagnosticsModel(QAbstractTableModel):
    """Table model for diagnostics."""

    headers = ("Severity", "Message", "Source")

    def __init__(self, diagnostics: list[Diagnostic]):
        super().__init__()
        self.diagnostics = diagnostics

    def rowCount(self, parent: QModelIndex = _INVALID_INDEX) -> int:
        """Return diagnostic count."""
        return len(self.diagnostics)

    def columnCount(self, parent: QModelIndex = _INVALID_INDEX) -> int:
        """Return column count."""
        return len(self.headers)

    def data(self, index: QModelIndex, role: int = Qt.ItemDataRole.DisplayRole) -> Any:
        """Return diagnostic table data."""
        if not index.isValid() or role != Qt.ItemDataRole.DisplayRole:
            return None
        diagnostic = self.diagnostics[index.row()]
        if index.column() == 0:
            return diagnostic.severity
        if index.column() == 1:
            return diagnostic.message
        source = diagnostic.source
        if source.file:
            return f"{source.file}:{source.line}:{source.column}"
        return ""

    def headerData(
        self,
        section: int,
        orientation: Qt.Orientation,
        role: int = Qt.ItemDataRole.DisplayRole,
    ) -> Any:
        """Return diagnostic table headers."""
        if (
            orientation == Qt.Orientation.Horizontal
            and role == Qt.ItemDataRole.DisplayRole
        ):
            return self.headers[section]
        return None


def object_to_text(value: Any) -> str:
    """Render a DatabaseIR object for the first detail panel."""
    if value is None:
        return ""
    if isinstance(value, list):
        return "\n".join(object_to_text(item) for item in value)
    if hasattr(value, "__dataclass_fields__"):
        data = asdict(value)
    else:
        data = value

    if isinstance(data, dict):
        lines = []
        for key, item in data.items():
            if isinstance(item, (dict, list)):
                lines.append(f"{key}: {item}")
            else:
                lines.append(f"{key}: {item}")
        return "\n".join(lines)
    return str(data)
