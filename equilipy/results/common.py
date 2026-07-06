"""Shared result-model helper functions."""

from __future__ import annotations

from typing import Any, Dict, List, Mapping, Optional

import numpy as np

from equilipy.exceptions import PostProcessError

PHASE_AMOUNT_ZERO_TOL = 1e-13


def normalize_phase_amount(value: Any) -> float:
    """Return phase amounts with Fortran-scale numerical noise set to zero."""
    numeric_value = float(value)
    if np.isfinite(numeric_value) and abs(numeric_value) <= PHASE_AMOUNT_ZERO_TOL:
        return 0.0
    return numeric_value


def normalize_phase_amounts(values: Any) -> np.ndarray:
    """Return phase amount arrays with tiny numerical noise set to zero."""
    amount_array = np.asarray(values, dtype=float).copy()
    finite_values = np.isfinite(amount_array)
    amount_array[finite_values & (np.abs(amount_array) <= PHASE_AMOUNT_ZERO_TOL)] = 0.0
    return amount_array


def clean_result_quantity_key(key: str) -> str:
    """Return a result mapping key without its trailing unit label."""
    return str(key).split(" [", 1)[0].strip()


def clean_quantity_mapping(values: Mapping[str, Any]) -> Dict[str, float]:
    """Return amount mappings keyed by clean component names."""
    cleaned: Dict[str, float] = {}
    for key, value in values.items():
        try:
            cleaned[clean_result_quantity_key(key)] = float(value)
        except (TypeError, ValueError):
            cleaned[clean_result_quantity_key(key)] = value
    return cleaned


def as_list_for_segments(value: Any) -> List[Any]:
    """Return a scalar or array-like result series as a plain list."""
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, list):
        return value
    if isinstance(value, tuple):
        return list(value)
    if value is None:
        return []
    return [value]


def component_names_from_amount_maps(
    context: Any,
    *amount_maps: Any,
) -> List[str]:
    """Return ordered component names from context and stored amount maps."""
    names = list(context.component_names) if context is not None else []
    for amount_map in amount_maps:
        if isinstance(amount_map, list):
            for item in amount_map:
                names.extend(component_names_from_amount_maps(None, item))
            continue
        if not isinstance(amount_map, dict):
            continue
        for key in amount_map:
            component = component_name_from_amount_key(str(key))
            if component:
                names.append(component)
    return unique_preserve_order(names)


def component_name_from_amount_key(key: str) -> str:
    """Extract `Al` from keys such as `Al [sp-mol]` or `Al [g]`."""
    return key.split("[", 1)[0].strip()


def temperature_from_kelvin(temperature: float, unit: str) -> float:
    """Convert a temperature in Kelvin to the requested display unit."""
    normalized_unit = unit.strip().upper()
    if normalized_unit in {"C", "CELSIUS"}:
        return temperature - 273.15
    if normalized_unit in {"F", "FAHRENHEIT"}:
        return temperature * 9.0 / 5.0 - 459.67
    if normalized_unit in {"R", "RANKINE"}:
        return temperature * 9.0 / 5.0
    return temperature


def temperature_unit_label(unit: str) -> str:
    """Return the compact temperature unit label used in result headers."""
    normalized_unit = unit.strip().upper()
    if normalized_unit in {"C", "CELSIUS"}:
        return "C"
    if normalized_unit in {"F", "FAHRENHEIT"}:
        return "F"
    if normalized_unit in {"R", "RANKINE"}:
        return "R"
    return "K"


def temperature_axis_unit(
    unit: Optional[str],
    input_unit: Optional[List[str]],
) -> Optional[str]:
    """Resolve a requested plotting temperature unit."""
    if unit is None:
        return None
    if unit.strip().lower() == "input":
        return input_unit[0] if input_unit else "K"
    return unit


def apply_temperature_axis_unit(
    axis_name: str,
    values: List[Any],
    unit: Optional[str],
    input_unit: Optional[List[str]],
    option_name: str,
) -> List[Any]:
    """Convert a plotting axis from K to a requested display temperature unit."""
    resolved_unit = temperature_axis_unit(unit, input_unit)
    if resolved_unit is None:
        return values
    if axis_name != "T":
        raise PostProcessError(
            f"{option_name} can only be used when that axis is 'T'. "
            f"Received axis {axis_name!r}."
        )
    unit_label = temperature_unit_label(resolved_unit)
    return [
        None
        if value is None
        else temperature_from_kelvin(float(value), unit_label)
        for value in values
    ]


def unique_preserve_order(values: List[str]) -> List[str]:
    """Return values without duplicates while preserving first-seen order."""
    unique = []
    for value in values:
        if value not in unique:
            unique.append(value)
    return unique


def phase_is_liquid(
    phase_name: str,
    liquid_phase_name: str | None = None,
) -> bool:
    """Return whether a phase is the active liquid phase."""
    normalized_phase_name = phase_name.strip()
    if liquid_phase_name:
        return normalized_phase_name == liquid_phase_name.strip()
    return "LIQ" in normalized_phase_name.upper()


def pad_missing_cumulative_phase_amounts(
    last_cumulative_phases_n: Dict[str, float],
    last_cumulative_phases_w: Dict[str, float],
    new_cumulative_phases_n: Dict[str, float],
    new_cumulative_phases_w: Dict[str, float],
    liquid_phase_name: str | None = None,
) -> None:
    """Carry or zero cumulative phase amounts absent from the current step."""
    for phase_name, last_amount_n in last_cumulative_phases_n.items():
        if phase_name in new_cumulative_phases_n:
            continue

        last_amount_w = last_cumulative_phases_w.get(phase_name, 0.0)
        if phase_is_liquid(phase_name, liquid_phase_name):
            new_cumulative_phases_n[phase_name] = 0.0
            new_cumulative_phases_w[phase_name] = 0.0
        else:
            new_cumulative_phases_n[phase_name] = last_amount_n
            new_cumulative_phases_w[phase_name] = last_amount_w


def combine_dicts(dict_list: List[Dict[str, float]]) -> Dict[str, np.ndarray]:
    """
    Combine row dictionaries into a columnar dictionary.

    This is the high-performance, dependency-free alternative to
    `pd.DataFrame(list_of_dicts).fillna(0.0)`.
    """
    if not dict_list:
        return {}

    all_keys = set().union(*dict_list)
    col_dict: Dict[str, list] = {}

    for i, row_dict in enumerate(dict_list):
        n = 1
        for key, value in row_dict.items():
            if key.startswith("T"):
                try:
                    if hasattr(value, "__iter__") and not isinstance(value, str):
                        n = len(value)
                except TypeError:
                    pass
                break

        for key in all_keys:
            values = col_dict.setdefault(key, [])
            if key in row_dict:
                value = row_dict[key]
                if not (hasattr(value, "__iter__") and not isinstance(value, str)):
                    value = [value] * n
                values.extend(np.array(value))
            else:
                values.extend(np.zeros(n))

        col_dict.setdefault("task_id", []).extend(np.ones(n, dtype=int) * (i + 1))

    return {key: np.array(values) for key, values in col_dict.items()}
