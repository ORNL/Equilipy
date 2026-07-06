"""Column-oriented result table helpers."""

from __future__ import annotations

import csv
import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


@dataclass(frozen=True)
class ResultColumn:
    """Metadata for one result table column."""

    key: str
    label: str | None = None
    group: str = "other"
    unit: str | None = None
    default: bool = True

    def __post_init__(self) -> None:
        """Fill a missing display label with the column key."""
        if self.label is None:
            object.__setattr__(self, "label", self.key)


class ResultTable:
    """Column-oriented table with optional result-column metadata."""

    def __init__(
        self,
        data: Mapping[str, Any] | None = None,
        columns: Sequence[ResultColumn] | None = None,
    ) -> None:
        self._data = dict(data or {})
        if columns is None:
            columns = columns_from_dict(self._data)
        self._columns = _merge_columns_with_data(columns, self._data)

    @classmethod
    def from_dict(
        cls,
        data: Mapping[str, Any],
        columns: Sequence[ResultColumn] | None = None,
    ) -> "ResultTable":
        """Build a result table from an existing flattened result dictionary."""
        return cls(data, columns=columns)

    @property
    def columns(self) -> list[ResultColumn]:
        """Return column metadata in display/export order."""
        return list(self._columns)

    @property
    def column_keys(self) -> list[str]:
        """Return column keys in display/export order."""
        return [column.key for column in self._columns]

    def available_columns(self) -> list[str]:
        """Return every selectable column key."""
        return self.column_keys

    def default_columns(self) -> list[ResultColumn]:
        """Return columns marked for default display/export."""
        return [column for column in self._columns if column.default]

    def default_column_keys(self) -> list[str]:
        """Return column keys marked for default display/export."""
        return [column.key for column in self.default_columns()]

    def select(self, keys: Iterable[str]) -> "ResultTable":
        """Return a new table with only the requested columns, in that order."""
        selected_keys = list(keys)
        missing = [key for key in selected_keys if key not in self._data]
        if missing:
            raise KeyError(f"Unknown result column(s): {', '.join(missing)}")

        metadata_by_key = {column.key: column for column in self._columns}
        selected_columns = [
            metadata_by_key.get(key, ResultColumn(key=key)) for key in selected_keys
        ]
        selected_data = {key: self._data[key] for key in selected_keys}
        return ResultTable(selected_data, selected_columns)

    def select_groups(self, groups: str | Iterable[str]) -> "ResultTable":
        """Return a new table with columns whose group is selected."""
        if isinstance(groups, str):
            group_names = {groups}
        else:
            group_names = set(groups)
        keys = [column.key for column in self._columns if column.group in group_names]
        return self.select(keys)

    def to_dict(self) -> dict[str, Any]:
        """Return column-oriented data in this table's current column order."""
        return {column.key: self._data[column.key] for column in self._columns}

    def to_rows(self) -> list[dict[str, Any]]:
        """Return row-oriented data, broadcasting scalar columns if needed."""
        normalized: dict[str, tuple[bool, list[Any] | Any]] = {}
        sequence_lengths: list[tuple[str, int]] = []

        for key in self.column_keys:
            value = self._data[key]
            if _is_column_sequence(value):
                values = _as_list(value)
                normalized[key] = (True, values)
                sequence_lengths.append((key, len(values)))
            else:
                normalized[key] = (False, _scalar_value(value))

        if not normalized:
            return []

        if not sequence_lengths:
            row_count = 1
        else:
            row_count = sequence_lengths[0][1]
            inconsistent = [
                f"{key}={length}"
                for key, length in sequence_lengths
                if length != row_count
            ]
            if inconsistent:
                detail = ", ".join(inconsistent)
                raise ValueError(
                    "Result table columns have inconsistent lengths: "
                    f"{detail}"
                )

        rows: list[dict[str, Any]] = []
        for index in range(row_count):
            row: dict[str, Any] = {}
            for key, (is_sequence, value) in normalized.items():
                row[key] = value[index] if is_sequence else value
            rows.append(row)
        return rows

    def to_csv(self, path: str | Path) -> Path:
        """Write the selected table columns to a CSV file."""
        output_path = Path(path)
        rows = self.to_rows()
        with output_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=self.column_keys)
            writer.writeheader()
            writer.writerows(rows)
        return output_path

    def to_polars(self) -> Any:
        """Return a Polars DataFrame for the selected columns."""
        try:
            import polars as pl
        except ImportError as exc:
            raise ImportError(
                "ResultTable.to_polars() requires polars to be installed."
            ) from exc
        return pl.DataFrame(self.to_rows())

    def to_pandas(self) -> Any:
        """Return a pandas DataFrame for the selected columns."""
        try:
            import pandas as pd
        except ImportError as exc:
            raise ImportError(
                "ResultTable.to_pandas() requires pandas to be installed."
            ) from exc
        return pd.DataFrame(self.to_rows(), columns=self.column_keys)


def columns_from_dict(data: Mapping[str, Any]) -> list[ResultColumn]:
    """Infer basic result-column metadata from flattened result keys."""
    return [
        ResultColumn(
            key=key,
            label=_infer_label(key),
            group=_infer_group(key),
            unit=_infer_unit(key),
        )
        for key in data
    ]


def _merge_columns_with_data(
    columns: Sequence[ResultColumn],
    data: Mapping[str, Any],
) -> list[ResultColumn]:
    columns_by_key = {column.key: column for column in columns}
    ordered: list[ResultColumn] = []
    seen: set[str] = set()

    for column in columns:
        if column.key in data and column.key not in seen:
            ordered.append(column)
            seen.add(column.key)

    for key in data:
        if key not in seen:
            ordered.append(columns_by_key.get(key, ResultColumn(key=key)))

    return ordered


def _infer_group(key: str) -> str:
    if key == "task_id":
        return "condition"
    if key.startswith(("T [", "P [")):
        return "condition"
    if key in {"G [J]", "H [J]", "S [J/K]", "Cp [J/K]", "Q [J]"}:
        return "thermodynamic"
    if key in {"label", "Label"}:
        return "state"
    if key in {"fl", "fs", "fl_w", "fs_w"} or key.startswith(("fl ", "fs ")):
        return "solidification_fraction"
    if key.startswith("stable_phase_"):
        return "stable_phase_summary"
    if "_amount_n" in key or "_amount_w" in key:
        return "phase_amount"
    if key.endswith("_stability"):
        return "phase_stability"
    if "_endmembers_" in key:
        return "phase_endmembers"
    if "_elements_" in key:
        return "phase_elements"
    if (
        "_partial_gibbs_" in key
        or "_standard_gibbs_energy_" in key
        or "_activity_" in key
        or "_partial_enthalpy_" in key
        or "_partial_entropy_" in key
        or "_partial_heat_capacity_" in key
    ):
        return "phase_endmembers"
    if (
        key.startswith(("n_", "w_"))
        and (key.endswith(" [sp-mol]") or key.endswith(" [g]"))
    ):
        return "composition"
    return "other"


def _infer_label(key: str) -> str:
    """Infer a compact human-facing label from an internal result key."""
    if key == "task_id":
        return "Task ID"
    if key in {"label", "Label"}:
        return "Label"

    if key.startswith("n_") and key.endswith(" [sp-mol]"):
        return f"n({_strip_prefix_suffix(key, 'n_', ' [sp-mol]')})"
    if key.startswith("w_") and key.endswith(" [g]"):
        return f"w({_strip_prefix_suffix(key, 'w_', ' [g]')})"

    if key == "stable_phase_amount_n [sp-mol]":
        return "n(stable phases)"
    if key == "stable_phase_amount_w [g]":
        return "w(stable phases)"
    if key == "stable_phase_amount_n_basis":
        return "basis(n(stable phases))"

    if key.endswith("_amount_n_basis"):
        return f"basis(n({key.rsplit('_amount_n_basis', 1)[0]}))"
    if "_amount_n" in key:
        return f"n({key.split('_amount_n', 1)[0]})"
    if "_amount_w" in key:
        return f"w({key.split('_amount_w', 1)[0]})"

    endmember_label = _phase_fraction_label(
        key,
        marker="_endmembers_",
        mole_tag="x",
        mass_tag="w",
        mole_prefix="x",
        mass_prefix="w",
    )
    if endmember_label is not None:
        return endmember_label

    element_label = _phase_fraction_label(
        key,
        marker="_elements_",
        mole_tag="x",
        mass_tag="w",
        mole_prefix="x_el",
        mass_prefix="w_el",
    )
    if element_label is not None:
        return element_label

    property_label = _phase_property_label(key)
    if property_label is not None:
        return property_label

    return key


def _strip_prefix_suffix(key: str, prefix: str, suffix: str) -> str:
    return key[len(prefix) : -len(suffix)]


def _phase_fraction_label(
    key: str,
    *,
    marker: str,
    mole_tag: str,
    mass_tag: str,
    mole_prefix: str,
    mass_prefix: str,
) -> str | None:
    if marker not in key:
        return None
    phase_name, suffix = key.split(marker, 1)

    if "_" not in suffix:
        return None
    fraction_tag, species_name = suffix.split("_", 1)
    if fraction_tag == mole_tag:
        return f"{mole_prefix}({species_name}@{phase_name})"
    if fraction_tag == mass_tag:
        return f"{mass_prefix}({species_name}@{phase_name})"
    return None


def _phase_property_label(key: str) -> str | None:
    property_markers = (
        ("_partial_gibbs_", "g"),
        ("_standard_gibbs_energy_", "Gref"),
        ("_activity_", "a"),
        ("_partial_enthalpy_", "h"),
        ("_partial_entropy_", "s"),
        ("_partial_heat_capacity_", "cp"),
    )
    base_key = key.rsplit(" [", 1)[0]
    for marker, symbol in property_markers:
        if marker not in base_key:
            continue
        phase_name, species_name = base_key.split(marker, 1)
        if phase_name and species_name:
            return f"{symbol}({species_name}@{phase_name})"
    return None


def _infer_unit(key: str) -> str | None:
    match = re.search(r"\[([^\[\]]+)\]$", key)
    return match.group(1) if match else None


def _is_column_sequence(value: Any) -> bool:
    if isinstance(value, np.ndarray):
        return value.ndim > 0
    return isinstance(value, (list, tuple))


def _as_list(value: Any) -> list[Any]:
    if isinstance(value, np.ndarray):
        return value.tolist()
    return list(value)


def _scalar_value(value: Any) -> Any:
    if isinstance(value, np.ndarray) and value.ndim == 0:
        return value.item()
    if isinstance(value, np.generic):
        return value.item()
    return value
