"""Composition normalization helpers for calculation inputs."""

from __future__ import annotations

import re
from collections.abc import Sequence
from typing import Any

import numpy as np

from . import variables as var
from .variables import cPeriodicTable

_FORMULA_NUMBER = re.compile(r"(?:\d+(?:\.\d*)?|\.\d+)")


def expand_condition_species(
    database: dict[str, Any],
    condition: dict[str, Any],
    amount_unit: str = "moles",
) -> dict[str, Any]:
    """
    Expand formula/species condition columns into database element columns.

    The returned condition always keeps ``T`` and ``P`` first, followed by element
    columns. Element columns can come from direct element names, chemical formulae
    such as ``Al13Fe4``, or pure-compound names found in the database.
    """
    if "T" not in condition or "P" not in condition:
        raise ValueError("condition must include T and P columns")

    _clear_pseudo_component_metadata()
    columns = {key: _as_column(value) for key, value in condition.items()}
    row_count = _condition_row_count(columns)
    normalized: dict[str, np.ndarray] = {
        "T": columns["T"].astype(float),
        "P": columns["P"].astype(float),
    }
    element_map = database_element_map(database)
    atomic_masses = database_atomic_mass_map(database)
    component_names: list[str] = []
    component_stoichiometries: list[dict[str, float]] = []
    has_pseudo_component = False
    for species, values in columns.items():
        if species in {"T", "P"}:
            continue
        stoichiometry = formula_or_database_species_stoichiometry(species, database)
        component_names.append(str(species))
        component_stoichiometries.append(stoichiometry)
        has_pseudo_component = has_pseudo_component or not _is_direct_element_component(
            species,
            stoichiometry,
            element_map,
        )
        values = values.astype(float)
        if is_mass_amount_unit(amount_unit):
            molar_mass = _molar_mass(stoichiometry, atomic_masses)
            if molar_mass <= 0:
                raise ValueError(f"molar mass is unavailable for {species}")
            contributions = {
                element: values * coefficient * atomic_masses[element] / molar_mass
                for element, coefficient in stoichiometry.items()
            }
        else:
            contributions = {
                element: values * coefficient
                for element, coefficient in stoichiometry.items()
            }
        for element, contribution in contributions.items():
            db_element = element_map.get(element, element)
            if db_element not in normalized:
                normalized[db_element] = np.zeros(row_count, dtype=float)
            normalized[db_element] = normalized[db_element] + contribution

    if has_pseudo_component:
        _set_pseudo_component_metadata(
            component_names,
            component_stoichiometries,
            database,
        )

    if row_count == 1 and all(not _is_sequence(value) for value in condition.values()):
        return {key: value[0].item() for key, value in normalized.items()}
    return normalized


def formula_or_database_species_stoichiometry(
    species_name: str,
    database: dict[str, Any] | None,
) -> dict[str, float]:
    """Return elemental stoichiometry for a formula or database pure species."""
    element_map = database_element_map(database)
    try:
        stoichiometry = parse_formula_stoichiometry(species_name)
    except ValueError:
        formula_stoichiometry = _parse_formula_with_database_elements(
            species_name,
            element_map,
        )
        if formula_stoichiometry:
            return formula_stoichiometry
        species_stoichiometry = database_species_stoichiometry(database, species_name)
        if species_stoichiometry:
            return species_stoichiometry
        raise

    if element_map:
        missing = sorted(
            element for element in stoichiometry if element not in element_map
        )
        if missing:
            formula_stoichiometry = _parse_formula_with_database_elements(
                species_name,
                element_map,
            )
            if formula_stoichiometry:
                return formula_stoichiometry
            species_stoichiometry = database_species_stoichiometry(
                database,
                species_name,
            )
            if species_stoichiometry:
                return species_stoichiometry
            raise ValueError(
                "Element(s) not found in selected database: " + ", ".join(missing)
            )
    return stoichiometry


def parse_formula_stoichiometry(formula: str) -> dict[str, float]:
    """Parse a chemical formula into elemental stoichiometry."""
    text = formula.strip().replace(" ", "")
    if not text:
        raise ValueError("Species name cannot be empty.")
    position = 0

    def parse_number() -> float:
        nonlocal position
        match = _FORMULA_NUMBER.match(text, position)
        if not match:
            return 1.0
        position = match.end()
        return float(match.group(0))

    def parse_group(expect_close: bool = False) -> dict[str, float]:
        nonlocal position
        composition: dict[str, float] = {}
        while position < len(text):
            char = text[position]
            if char == ")":
                break
            if char == "(":
                position += 1
                inner = parse_group(expect_close=True)
                if position >= len(text) or text[position] != ")":
                    raise ValueError("Invalid formula: unmatched '('.")
                position += 1
                multiplier = parse_number()
                for element, coefficient in inner.items():
                    composition[element] = (
                        composition.get(element, 0.0) + coefficient * multiplier
                    )
                continue
            if not char.isalpha() or not char.isupper():
                raise ValueError(f"Invalid formula near '{char}'.")

            symbol = char
            position += 1
            if position < len(text) and text[position].islower():
                symbol += text[position]
                position += 1
            if symbol not in cPeriodicTable:
                raise ValueError(
                    f"Invalid element symbol '{symbol}'. "
                    "Use standard capitalization, for example Al."
                )
            coefficient = parse_number()
            composition[symbol] = composition.get(symbol, 0.0) + coefficient

        if expect_close and (position >= len(text) or text[position] != ")"):
            raise ValueError("Invalid formula: unmatched '('.")
        if not composition:
            raise ValueError("Invalid formula.")
        return composition

    composition = parse_group()
    if position != len(text):
        raise ValueError("Invalid formula: unmatched ')'.")
    return composition


def database_species_stoichiometry(
    database: dict[str, Any] | None,
    species_name: str,
) -> dict[str, float] | None:
    """Return database pure-species stoichiometry when available."""
    if not isinstance(database, dict):
        return None
    species_names = database.get("cSpeciesNameCS")
    stoichiometry = database.get("dStoichSpeciesCS")
    n_pure_species = int(database.get("nPureSpeciesCS", 0) or 0)
    if species_names is None or stoichiometry is None or n_pure_species <= 0:
        return None

    element_names = list(database_element_map(database))
    start_index = max(0, len(species_names) - n_pure_species)
    target = species_match_key(species_name)
    for index in range(start_index, len(species_names)):
        raw_name = str(species_names[index]).strip()
        candidates = {
            raw_name,
            raw_name.strip("'"),
            raw_name.split("_", 1)[0].strip("'"),
            raw_name.split("(", 1)[0].strip("'"),
        }
        if target not in {species_match_key(candidate) for candidate in candidates}:
            continue
        row = stoichiometry[index]
        composition: dict[str, float] = {}
        for element, coefficient in zip(element_names, row, strict=False):
            value = float(coefficient)
            if value:
                composition[element] = value
        return composition or None
    return None


def species_match_key(value: str) -> str:
    """Return a normalized species key for lookup comparisons."""
    return value.strip().strip("'").replace(" ", "")


def database_element_map(database: dict[str, Any] | None) -> dict[str, str]:
    """Return canonical element symbols mapped to database element names."""
    if not isinstance(database, dict):
        return {}
    elements = database.get("cElementNameCS")
    if elements is None:
        return {}
    element_map: dict[str, str] = {}
    for element in elements:
        raw = str(element).strip()
        if raw:
            element_map[canonical_element_symbol(raw)] = raw
    return element_map


def database_atomic_mass_map(database: dict[str, Any] | None) -> dict[str, float]:
    """Return canonical element symbols mapped to database atomic masses."""
    if not isinstance(database, dict):
        return {}
    element_map = database_element_map(database)
    masses = database.get("dAtomicMass")
    if masses is None:
        return {}
    atomic_masses: dict[str, float] = {}
    for element, mass in zip(element_map, masses, strict=False):
        atomic_masses[element] = float(mass)
    return atomic_masses


def canonical_element_symbol(symbol: str) -> str:
    """Return standard element capitalization."""
    stripped = symbol.strip()
    if stripped.lower() == "e-":
        return "e-"
    if not stripped:
        return ""
    return stripped[:1].upper() + stripped[1:].lower()


def is_mass_amount_unit(amount_unit: str) -> bool:
    """Return whether an amount unit is mass based."""
    unit = amount_unit.strip().lower()
    return unit in {"g", "kg", "lbs", "wt%", "wt.%"} or unit.startswith(
        ("gram", "kilo", "pound", "mass fraction", "weight fraction")
    )


def _clear_pseudo_component_metadata() -> None:
    var.iPseudoComponentDependentElementSys = []
    var.iPseudoComponentDependentElementDBIndex = []
    var.cPseudoComponentNameSys = []
    var.dPseudoComponentStoichSys = np.empty((0, 0), dtype=float)


def _set_pseudo_component_metadata(
    component_names: list[str],
    component_stoichiometries: list[dict[str, float]],
    database: dict[str, Any] | None,
) -> None:
    """Store input component stoichiometry for result-basis projections."""
    element_names = [
        str(element).strip()
        for element in (database or {}).get("cElementNameCS", [])
    ]
    if not element_names:
        return

    element_index = {element: index for index, element in enumerate(element_names)}
    element_map = database_element_map(database)
    stoich = np.zeros((len(component_stoichiometries), len(element_names)), dtype=float)
    for component_index, component_stoich in enumerate(component_stoichiometries):
        for element, coefficient in component_stoich.items():
            database_element = element_map.get(element, element)
            if database_element not in element_index:
                continue
            element_column = element_index[database_element]
            stoich[component_index, element_column] = float(coefficient)

    var.cPseudoComponentNameSys = list(component_names)
    var.dPseudoComponentStoichSys = stoich


def _is_direct_element_component(
    species: str,
    stoichiometry: dict[str, float],
    element_map: dict[str, str],
) -> bool:
    """Return whether an input column is a plain element amount."""
    if len(stoichiometry) != 1:
        return False
    element, coefficient = next(iter(stoichiometry.items()))
    if abs(float(coefficient) - 1.0) > 1e-12:
        return False
    return element_map.get(element, element).strip().lower() == species.strip().lower()


def _parse_condition_formula(
    species: str,
    element_map: dict[str, str],
) -> dict[str, float] | None:
    try:
        stoichiometry = parse_formula_stoichiometry(species)
    except ValueError:
        stoichiometry = _parse_formula_with_database_elements(species, element_map)
    else:
        if element_map and any(element not in element_map for element in stoichiometry):
            database_stoichiometry = _parse_formula_with_database_elements(
                species,
                element_map,
            )
            if database_stoichiometry:
                stoichiometry = database_stoichiometry
    return stoichiometry


def _parse_formula_with_database_elements(
    formula: str,
    element_map: dict[str, str],
) -> dict[str, float] | None:
    text = formula.strip().replace(" ", "")
    if (
        not text
        or not all(char.isalnum() or char == "." for char in text)
        or not element_map
    ):
        return None

    symbols_by_upper = {element.upper(): element for element in element_map}
    symbol_keys = sorted(symbols_by_upper, key=len, reverse=True)
    text_upper = text.upper()
    position = 0
    composition: dict[str, float] = {}

    while position < len(text):
        matched_key = next(
            (
                symbol_key
                for symbol_key in symbol_keys
                if text_upper.startswith(symbol_key, position)
            ),
            None,
        )
        if matched_key is None:
            return None
        position += len(matched_key)

        match = _FORMULA_NUMBER.match(text, position)
        if match:
            coefficient = float(match.group(0))
            position = match.end()
        else:
            coefficient = 1.0

        element = symbols_by_upper[matched_key]
        composition[element] = composition.get(element, 0.0) + coefficient

    return composition or None


def _molar_mass(
    stoichiometry: dict[str, float],
    atomic_masses: dict[str, float],
) -> float:
    total = 0.0
    for element, coefficient in stoichiometry.items():
        if element not in atomic_masses:
            raise ValueError(f"atomic mass is unavailable for {element}")
        total += coefficient * atomic_masses[element]
    return total


def _as_column(value: Any) -> np.ndarray:
    if not _is_sequence(value):
        return np.asarray([value])
    if isinstance(value, np.ndarray):
        return value
    if hasattr(value, "to_numpy"):
        return np.asarray(value.to_numpy())
    if hasattr(value, "to_list"):
        return np.asarray(value.to_list())
    if hasattr(value, "to_pylist"):
        return np.asarray(value.to_pylist())
    return np.asarray(value)


def _condition_row_count(columns: dict[str, np.ndarray]) -> int:
    row_counts = {key: len(value) for key, value in columns.items()}
    unique_counts = set(row_counts.values())
    if len(unique_counts) != 1:
        raise ValueError("condition columns must have the same length")
    return unique_counts.pop()


def _is_sequence(value: Any) -> bool:
    if isinstance(value, str):
        return False
    if isinstance(value, np.ndarray):
        return value.ndim > 0
    if np.isscalar(value):
        return False
    if isinstance(value, Sequence):
        return True
    if (
        hasattr(value, "to_numpy")
        or hasattr(value, "to_list")
        or hasattr(value, "to_pylist")
    ):
        return True
    try:
        len(value)
    except TypeError:
        return False
    return True
