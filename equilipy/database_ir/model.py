"""Dataclasses for a format-neutral thermodynamic database IR."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any


@dataclass(slots=True)
class SourceRef:
    """Location of an object in an original source file."""

    file: str = ""
    line: int = 0
    column: int = 0
    command: str = ""

    @classmethod
    def from_dict(cls, data: dict[str, Any] | None) -> "SourceRef":
        """Create a source reference from serialized data."""
        if data is None:
            return cls()
        return cls(
            file=str(data.get("file", "")),
            line=int(data.get("line", 0)),
            column=int(data.get("column", 0)),
            command=str(data.get("command", "")),
        )


@dataclass(slots=True)
class Diagnostic:
    """Validation or parsing diagnostic attached to the database."""

    severity: str
    message: str
    source: SourceRef = field(default_factory=SourceRef)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "Diagnostic":
        """Create a diagnostic from serialized data."""
        return cls(
            severity=str(data.get("severity", "info")),
            message=str(data.get("message", "")),
            source=SourceRef.from_dict(data.get("source")),
        )


@dataclass(slots=True)
class Element:
    """Chemical element or pseudo-element in a database."""

    symbol: str
    atomic_mass: float = 0.0
    reference_state: str = ""
    source: SourceRef = field(default_factory=SourceRef)
    h298: float = 0.0
    s298: float = 0.0
    tdb_symbol: str = ""

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "Element":
        """Create an element from serialized data."""
        return cls(
            symbol=str(data["symbol"]),
            atomic_mass=float(data.get("atomic_mass", 0.0)),
            reference_state=str(data.get("reference_state", "")),
            source=SourceRef.from_dict(data.get("source")),
            h298=float(data.get("h298", 0.0)),
            s298=float(data.get("s298", 0.0)),
            tdb_symbol=str(data.get("tdb_symbol", data.get("symbol", ""))),
        )


@dataclass(slots=True)
class Species:
    """Species or endmember with elemental composition."""

    name: str
    composition: dict[str, float] = field(default_factory=dict)
    charge: float = 0.0
    source: SourceRef = field(default_factory=SourceRef)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "Species":
        """Create a species from serialized data."""
        composition = {
            str(symbol): float(amount)
            for symbol, amount in data.get("composition", {}).items()
        }
        return cls(
            name=str(data["name"]),
            composition=composition,
            charge=float(data.get("charge", 0.0)),
            source=SourceRef.from_dict(data.get("source")),
        )


@dataclass(slots=True)
class GibbsRange:
    """One piecewise Gibbs energy expression over a temperature interval."""

    T_min: float
    T_max: float
    Gibbs: str
    status: str = "N"

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "GibbsRange":
        """Create a Gibbs range from serialized data."""
        return cls(
            T_min=float(data.get("T_min", data.get("t_min", 0.0))),
            T_max=float(data.get("T_max", data.get("t_max", 0.0))),
            Gibbs=str(data.get("Gibbs", data.get("expression", ""))),
            status=str(data.get("status", "N")).upper()[:1] or "N",
        )

    @property
    def t_min(self) -> float:
        """Backward-compatible lower-temperature alias."""
        return self.T_min

    @t_min.setter
    def t_min(self, value: float) -> None:
        self.T_min = value

    @property
    def t_max(self) -> float:
        """Backward-compatible upper-temperature alias."""
        return self.T_max

    @t_max.setter
    def t_max(self, value: float) -> None:
        self.T_max = value

    @property
    def expression(self) -> str:
        """Backward-compatible Gibbs-expression alias."""
        return self.Gibbs

    @expression.setter
    def expression(self, value: str) -> None:
        self.Gibbs = value


@dataclass(slots=True)
class TdbFunctionTerm:
    """A referenced TDB function term and its multiplier."""

    name: str
    factor: float = 1.0

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "TdbFunctionTerm":
        """Create a referenced-function term from serialized data."""
        return cls(
            name=str(data.get("name", "")),
            factor=float(data.get("factor", 1.0)),
        )


@dataclass(slots=True)
class TdbFunctionSource:
    """TDB-specific source/provenance for a function-like record."""

    source_expression: str = ""
    pertinent_functions: list[TdbFunctionTerm] = field(default_factory=list)

    @classmethod
    def from_dict(cls, data: dict[str, Any] | None) -> "TdbFunctionSource":
        """Create TDB source metadata from serialized data."""
        if data is None:
            return cls()
        return cls(
            source_expression=str(data.get("source_expression", "")),
            pertinent_functions=[
                TdbFunctionTerm.from_dict(item)
                for item in data.get("pertinent_functions", [])
            ],
        )


@dataclass(slots=True)
class TdbCommand:
    """TDB setup/type command preserved for later semantic expansion."""

    kind: str
    command: str
    args: list[str] = field(default_factory=list)
    active: bool = True
    raw: str = ""
    parsed: dict[str, str] = field(default_factory=dict)
    source: SourceRef = field(default_factory=SourceRef)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "TdbCommand":
        """Create a preserved TDB command from serialized data."""
        return cls(
            kind=str(data.get("kind", "generic")),
            command=str(data.get("command", "")),
            args=[str(item) for item in data.get("args", [])],
            active=bool(data.get("active", True)),
            raw=str(data.get("raw", "")),
            parsed={str(k): str(v) for k, v in data.get("parsed", {}).items()},
            source=SourceRef.from_dict(data.get("source")),
        )


@dataclass(slots=True)
class TdbReference:
    """Bibliographic reference entry from a TDB reference table."""

    key: str
    text: str
    section: str = "LIST-OF-REFERENCE"
    source: SourceRef = field(default_factory=SourceRef)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "TdbReference":
        """Create a TDB reference from serialized data."""
        return cls(
            key=str(data.get("key", "")),
            text=str(data.get("text", "")),
            section=str(data.get("section", "LIST-OF-REFERENCE")),
            source=SourceRef.from_dict(data.get("source")),
        )


@dataclass(slots=True)
class FunctionDefinition:
    """Named thermodynamic function expression."""

    name: str
    expression: str
    temperature_ranges: list[tuple[float, float]] = field(default_factory=list)
    source: SourceRef = field(default_factory=SourceRef)
    gibbs_ranges: list[GibbsRange] = field(default_factory=list)
    tdb_src: TdbFunctionSource = field(default_factory=TdbFunctionSource)
    dirty: bool = False

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "FunctionDefinition":
        """Create a function definition from serialized data."""
        gibbs_ranges = [
            GibbsRange.from_dict(item)
            for item in data.get("gibbs_ranges", [])
        ]
        ranges = [
            (float(item[0]), float(item[1]))
            for item in data.get("temperature_ranges", [])
        ] or [
            (gibbs_range.t_min, gibbs_range.t_max)
            for gibbs_range in gibbs_ranges
        ]
        tdb_src = TdbFunctionSource.from_dict(data.get("tdb_src"))
        expression = str(data.get("expression", "")) or tdb_src.source_expression
        return cls(
            name=str(data["name"]),
            expression=expression,
            temperature_ranges=ranges,
            source=SourceRef.from_dict(data.get("source")),
            gibbs_ranges=gibbs_ranges,
            tdb_src=tdb_src,
            dirty=bool(data.get("dirty", False)),
        )


@dataclass(slots=True)
class ConstituentSet:
    """Constituent list on one phase sublattice."""

    sublattice: int
    species: list[str] = field(default_factory=list)
    site_ratio: float = 1.0

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "ConstituentSet":
        """Create a constituent set from serialized data."""
        return cls(
            sublattice=int(data.get("sublattice", 1)),
            species=[str(name) for name in data.get("species", [])],
            site_ratio=float(data.get("site_ratio", 1.0)),
        )


@dataclass(slots=True)
class Phase:
    """Solution or pure phase definition."""

    name: str
    model: str = ""
    constituents: list[ConstituentSet] = field(default_factory=list)
    source: SourceRef = field(default_factory=SourceRef)
    metadata: dict[str, str] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "Phase":
        """Create a phase from serialized data."""
        return cls(
            name=str(data["name"]),
            model=str(data.get("model", "")),
            constituents=[
                ConstituentSet.from_dict(item)
                for item in data.get("constituents", [])
            ],
            source=SourceRef.from_dict(data.get("source")),
            metadata={str(k): str(v) for k, v in data.get("metadata", {}).items()},
        )


@dataclass(slots=True)
class Parameter:
    """Thermodynamic parameter attached to a phase."""

    phase_name: str
    parameter_type: str
    target: list[str] = field(default_factory=list)
    order: int = 0
    expression: str = ""
    source: SourceRef = field(default_factory=SourceRef)
    metadata: dict[str, str] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "Parameter":
        """Create a parameter from serialized data."""
        return cls(
            phase_name=str(data["phase_name"]),
            parameter_type=str(data.get("parameter_type", "G")),
            target=[str(name) for name in data.get("target", [])],
            order=int(data.get("order", 0)),
            expression=str(data.get("expression", "")),
            source=SourceRef.from_dict(data.get("source")),
            metadata={str(k): str(v) for k, v in data.get("metadata", {}).items()},
        )


@dataclass(slots=True)
class DatabaseIR:
    """Format-neutral model used by parsers, GUI, writers, and adapters."""

    schema_version: str = "database_ir.v2"
    name: str = "Untitled database"
    metadata: dict[str, str] = field(default_factory=dict)
    elements: list[Element] = field(default_factory=list)
    species: list[Species] = field(default_factory=list)
    functions: list[FunctionDefinition] = field(default_factory=list)
    phases: list[Phase] = field(default_factory=list)
    parameters: list[Parameter] = field(default_factory=list)
    tdb_commands: list[TdbCommand] = field(default_factory=list)
    references: list[TdbReference] = field(default_factory=list)
    diagnostics: list[Diagnostic] = field(default_factory=list)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "DatabaseIR":
        """Create a database IR from serialized data."""
        return cls(
            schema_version=str(data.get("schema_version", "database_ir.v1")),
            name=str(data.get("name", "Untitled database")),
            metadata={str(k): str(v) for k, v in data.get("metadata", {}).items()},
            elements=[Element.from_dict(item) for item in data.get("elements", [])],
            species=[Species.from_dict(item) for item in data.get("species", [])],
            functions=[
                FunctionDefinition.from_dict(item)
                for item in data.get("functions", [])
            ],
            phases=[Phase.from_dict(item) for item in data.get("phases", [])],
            parameters=[
                Parameter.from_dict(item) for item in data.get("parameters", [])
            ],
            tdb_commands=[
                TdbCommand.from_dict(item)
                for item in data.get("tdb_commands", [])
            ],
            references=[
                TdbReference.from_dict(item)
                for item in data.get("references", [])
            ],
            diagnostics=[
                Diagnostic.from_dict(item) for item in data.get("diagnostics", [])
            ],
        )

    def to_dict(self) -> dict[str, Any]:
        """Serialize the database IR to plain Python containers."""
        return asdict(self)

    def summary_counts(self) -> dict[str, int]:
        """Return object counts for high-level UI summaries."""
        return {
            "elements": len(self.elements),
            "species": len(self.species),
            "functions": len(self.functions),
            "phases": len(self.phases),
            "parameters": len(self.parameters),
            "diagnostics": len(self.diagnostics),
        }

    def validation_messages(self) -> list[Diagnostic]:
        """Return lightweight validation diagnostics for the initial GUI."""
        messages = list(self.diagnostics)
        phase_names = {phase.name for phase in self.phases}
        species_names = {species.name for species in self.species}

        for parameter in self.parameters:
            if parameter.phase_name not in phase_names:
                messages.append(
                    Diagnostic(
                        "error",
                        f"Parameter references unknown phase '{parameter.phase_name}'.",
                        parameter.source,
                    )
                )

        for phase in self.phases:
            for constituent_set in phase.constituents:
                for name in constituent_set.species:
                    if name != "VA" and name not in species_names:
                        messages.append(
                            Diagnostic(
                                "warning",
                                (
                                    f"Phase '{phase.name}' references unknown "
                                    f"species '{name}'."
                                ),
                                phase.source,
                            )
                        )

        return messages

    def validate_tdb(self, *, auto_correct: bool = False):
        """Validate TDB-oriented thermodynamic records on this database."""
        from .validation import validate_tdb

        return validate_tdb(self, auto_correct=auto_correct)

    def write_tdb(
        self,
        path: str | Path,
        *,
        include_endmember_compounds: bool = False,
        order_disorder_alias_mode: str = "preserve",
        export_style: str = "native",
    ) -> Path:
        """Write this database as a TDB file."""
        from .tdb_writer import write_tdb

        return write_tdb(
            self,
            path,
            include_endmember_compounds=include_endmember_compounds,
            order_disorder_alias_mode=order_disorder_alias_mode,
            export_style=export_style,
        )


def load_database_ir(path: str | Path) -> DatabaseIR:
    """Load a database IR JSON or TDB file."""
    source_path = Path(path)
    if source_path.suffix.lower() == ".tdb":
        from .tdb import read_tdb

        return read_tdb(
            source_path,
            strict=False,
            remove_redundant_phases=False,
        )

    with source_path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    return DatabaseIR.from_dict(data)
