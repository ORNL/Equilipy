"""UI-side calculation session state."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass
class CalculationDatabase:
    """A database loaded into one calculation session."""

    id: str
    name: str
    path: str
    database: Any
    selected: bool = True


@dataclass
class CalculationModule:
    """UI-side state for one calculation module instance."""

    id: str
    kind: str
    name: str
    phase_names: list[str] = field(default_factory=list)
    temperature_unit: str = "K"
    pressure_unit: str = "atm"
    amount_unit: str = "moles"
    calculation_type: str = "single"
    solidification_model: str = "scheil"
    nucleation_undercooling: dict[str, float] = field(default_factory=dict)
    batch_condition: dict[str, list[float]] = field(default_factory=dict)
    batch_condition_path: str = ""
    batch_cpu_count: int = 0
    transition_search: bool = True
    start_from_liquidus: bool = True
    result: Any | None = None
    result_payload: dict[str, Any] = field(default_factory=dict)
    results_table_payload: dict[str, Any] = field(default_factory=dict)
    result_columns: list[str] = field(default_factory=list)
    source_path: str = ""
    runtime: dict[str, Any] = field(default_factory=dict)


@dataclass
class CalculationSession:
    """UI-side state for one calculation session."""

    id: str
    name: str
    temperature_unit: str = "K"
    pressure_unit: str = "atm"
    amount_unit: str = "moles"
    default_phase_names: list[str] = field(default_factory=list)
    result_column_defaults: dict[str, list[str]] = field(default_factory=dict)
    database_counter: int = 0
    databases: dict[str, CalculationDatabase] = field(default_factory=dict)
    module_counts: dict[str, int] = field(default_factory=dict)
    modules: dict[str, CalculationModule] = field(default_factory=dict)
    source_path: str = ""
