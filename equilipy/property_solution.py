"""Solution-phase thermodynamic property evaluation."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

import equilipy.equilifort as fort
import equilipy.variables as var

from .composition import (
    expand_condition_species,
    formula_or_database_species_stoichiometry,
)
from .database_ir.tdb_canonical import DISORDERED_PHASE_CANONICAL_NAMES
from .equilib_single import _preprocess_single
from .exceptions import EquilibError, InputConditionError
from .minimize import _raise_fortran_status, minimize
from .results import capture_result_context
from .results.equilib import (
    _phase_weighted_property,
    endmembers2elements,
)
from .results.property import SlnPropertyPoint, SlnPropertyResult, SpeciesProperty
from .utils import _dict2np


class _ExplicitCompositionUnavailable(Exception):
    """Signal that elemental composition cannot define unique endmember fractions."""


_COMPOSITION_TOL = 1e-7


def property_solution(
    database: dict,
    phase_name: str,
    condition: dict,
    unit: list | None = None,
    *,
    endmember_fractions: Mapping[str, float]
    | Sequence[Mapping[str, float]]
    | None = None,
    include_heat_capacity: bool = True,
) -> SlnPropertyResult:
    """Evaluate one solution phase at the requested NPT condition.

    Callers use the same NPT style as equilibrium calculations.  If the
    non-T/P keys are phase endmember names such as ``Al:Fe:VA``, this function
    evaluates the phase at those fixed endmember fractions.  If the supplied
    element/component columns map one-to-one onto every phase endmember, the
    fixed-endmember path is also used.  Otherwise, the selected phase is
    minimized at the supplied bulk composition.
    """
    if unit is None:
        unit = ["K", "atm", "moles"]

    condition = dict(condition)
    if endmember_fractions is None:
        condition, endmember_fractions = _resolve_ntp_endmember_input(
            database,
            phase_name,
            condition,
            unit,
        )
    elif not _has_element_amounts(condition):
        condition = _merge_endmember_derived_composition(
            database,
            phase_name,
            condition,
            endmember_fractions,
            unit,
        )

    expanded_condition = expand_condition_species(database, condition, unit[2])
    rows = _condition_rows(expanded_condition)
    fraction_rows = _endmember_fraction_rows(endmember_fractions, len(rows))
    context = capture_result_context(
        database,
        expanded_condition,
        unit,
        [phase_name],
    )

    points: list[SlnPropertyPoint] = []
    for row_condition, row_fractions in zip(rows, fraction_rows, strict=True):
        try:
            points.append(
                _property_solution_point(
                    database,
                    phase_name,
                    row_condition,
                    unit,
                    endmember_fractions=row_fractions,
                    include_heat_capacity=include_heat_capacity,
                )
            )
        finally:
            fort.resetthermo()

    return SlnPropertyResult(
        points[0] if len(points) == 1 else points,
        context=context,
    )


def _condition_rows(condition: dict) -> list[dict[str, float]]:
    headers, values = _dict2np(condition)
    return [
        dict(zip(headers, row, strict=False))
        for row in np.asarray(values, dtype=float)
    ]


def _resolve_ntp_endmember_input(
    database: dict,
    phase_name: str,
    condition: dict,
    unit: list,
) -> tuple[dict, list[dict[str, float]] | None]:
    """Resolve direct endmember input encoded as normal NPT species columns."""
    species_keys = _condition_species_keys(condition)
    if not species_keys:
        return condition, None

    element_names, endmember_names, stoich = _phase_stoich_from_database(
        database,
        phase_name,
    )
    endmember_lookup = _endmember_name_lookup(endmember_names)
    if all(":" in key for key in species_keys) and all(
        _matched_endmember_index(key, endmember_lookup) is not None
        for key in species_keys
    ):
        rows = _endmember_amount_fraction_rows(
            database,
            condition,
            species_keys,
            endmember_names,
            stoich,
            unit,
        )
        base_condition = _tp_condition(condition)
        merged = _merge_endmember_derived_composition(
            database,
            phase_name,
            base_condition,
            rows,
            unit,
        )
        return merged, rows

    if len(species_keys) != len(endmember_names):
        return condition, None

    component_to_endmember = _map_components_to_endmembers(
        database,
        species_keys,
        element_names,
        endmember_names,
        stoich,
    )
    if component_to_endmember is None:
        return condition, None

    rows = _component_amount_fraction_rows(
        database,
        condition,
        species_keys,
        component_to_endmember,
        endmember_names,
        unit,
    )
    return condition, rows


def _condition_species_keys(condition: Mapping[str, Any]) -> list[str]:
    return [
        str(key)
        for key in condition
        if str(key).strip() not in {"T", "P"}
    ]


def _tp_condition(condition: Mapping[str, Any]) -> dict:
    return {key: condition[key] for key in condition if str(key).strip() in {"T", "P"}}


def _has_element_amounts(condition: Mapping[str, Any]) -> bool:
    return any(str(key).strip() not in {"T", "P"} for key in condition)


def _merge_endmember_derived_composition(
    database: dict,
    phase_name: str,
    condition: dict,
    endmember_payload: Mapping[str, float] | Sequence[Mapping[str, float]],
    unit: list,
) -> dict:
    """Derive elemental composition from inline endmember fractions."""
    merged = dict(condition)
    count = _condition_row_count(condition)
    fraction_rows = _endmember_fraction_rows(endmember_payload, count)
    element_names, endmember_names, stoich = _phase_stoich_from_database(
        database,
        phase_name,
    )
    atomic_masses = _database_atomic_mass_map(database)

    element_values = {element: [] for element in element_names}
    for fraction_row in fraction_rows:
        fractions = _normalize_endmember_fractions(fraction_row, endmember_names)
        values = fractions @ stoich
        if _is_mass_unit(unit[2]):
            values = np.array(
                [
                    value * atomic_masses[element]
                    for element, value in zip(element_names, values, strict=False)
                ],
                dtype=float,
            )
        for element, value in zip(element_names, values, strict=False):
            element_values[element].append(float(value))

    for element, values in element_values.items():
        merged[element] = values[0] if count == 1 else values
    return merged


def _condition_row_count(condition: Mapping[str, Any]) -> int:
    lengths: list[int] = []
    for value in condition.values():
        if isinstance(value, (str, bytes)):
            continue
        try:
            array = np.asarray(value)
        except (TypeError, ValueError):
            continue
        if array.ndim > 0 and array.size > 1:
            lengths.append(int(array.size))
    if not lengths:
        return 1
    if len(set(lengths)) != 1:
        raise InputConditionError(
            "Vectorized T/P entries must have the same length for endmember input."
        )
    return lengths[0]


def _phase_stoich_from_database(
    database: dict,
    phase_name: str,
) -> tuple[list[str], list[str], np.ndarray]:
    phase_names = [str(name).strip() for name in database["cSolnPhaseNameCS"]]
    try:
        phase_index = phase_names.index(phase_name)
    except ValueError as exc:
        raise InputConditionError(
            f"Unknown solution phase {phase_name!r} for inline endmember input."
        ) from exc

    species_boundaries = np.append(
        [0],
        np.asarray(database["nSpeciesPhaseCS"], dtype=int),
    )
    first = int(species_boundaries[phase_index])
    last = int(species_boundaries[phase_index + 1])
    endmember_names = [
        str(name).strip() for name in database["cEndmemberNameCS"][first:last]
    ]
    element_names = [str(name).strip() for name in database["cElementNameCS"]]
    stoich = np.asarray(database["dStoichSpeciesCS"], dtype=float)[first:last, :]
    particles = np.asarray(database["iParticlesPerMoleCS"], dtype=float)[first:last]
    return element_names, endmember_names, stoich / particles[:, None]


def _database_atomic_mass_map(database: dict) -> dict[str, float]:
    return {
        str(symbol).strip(): float(mass)
        for symbol, mass in zip(
            database["cElementNameCS"],
            database["dAtomicMass"],
            strict=False,
        )
    }


def _is_mass_unit(unit_name: str) -> bool:
    return str(unit_name).strip().lower() in {
        "grams",
        "kilograms",
        "pounds",
        "g",
        "kg",
        "lbs",
        "mass fraction",
        "weight fraction",
        "wt%",
        "wt.%",
    }


def _endmember_name_lookup(endmember_names: Sequence[str]) -> dict[str, int]:
    lookup: dict[str, int] = {}
    for index, name in enumerate(endmember_names):
        text = str(name).strip()
        lookup[text] = index
        lookup[text.upper()] = index
    return lookup


def _matched_endmember_index(
    name: str,
    lookup: Mapping[str, int],
) -> int | None:
    text = str(name).strip()
    return lookup.get(text, lookup.get(text.upper()))


def _endmember_amount_fraction_rows(
    database: dict,
    condition: Mapping[str, Any],
    species_keys: Sequence[str],
    endmember_names: Sequence[str],
    stoich: np.ndarray,
    unit: Sequence[str],
) -> list[dict[str, float]]:
    lookup = _endmember_name_lookup(endmember_names)
    mass_map = _database_atomic_mass_map(database)
    amount_rows = _species_amount_rows(condition, species_keys)
    fraction_rows: list[dict[str, float]] = []
    for amount_row in amount_rows:
        amounts = np.zeros(len(endmember_names), dtype=float)
        for species, value in amount_row.items():
            index = _matched_endmember_index(species, lookup)
            if index is None:
                raise InputConditionError(f"Unknown endmember species {species!r}.")
            amounts[index] = _amount_to_moles(
                value,
                unit[2],
                stoich[index],
                mass_map,
                species,
            )
        fraction_rows.append(
            _fraction_mapping_from_amounts(amounts, endmember_names)
        )
    return fraction_rows


def _component_amount_fraction_rows(
    database: dict,
    condition: Mapping[str, Any],
    species_keys: Sequence[str],
    component_to_endmember: Mapping[str, int],
    endmember_names: Sequence[str],
    unit: Sequence[str],
) -> list[dict[str, float]]:
    amount_rows = _species_amount_rows(condition, species_keys)
    fraction_rows: list[dict[str, float]] = []
    for amount_row in amount_rows:
        amounts = np.zeros(len(endmember_names), dtype=float)
        for species, value in amount_row.items():
            stoich = _component_stoich_vector(database, species)
            amounts[component_to_endmember[species]] = _amount_to_moles(
                value,
                unit[2],
                stoich,
                _database_atomic_mass_map(database),
                species,
            )
        fraction_rows.append(
            _fraction_mapping_from_amounts(amounts, endmember_names)
        )
    return fraction_rows


def _species_amount_rows(
    condition: Mapping[str, Any],
    species_keys: Sequence[str],
) -> list[dict[str, float]]:
    subset = {key: condition[key] for key in species_keys}
    headers, values = _dict2np(subset)
    return [
        {
            str(header).strip(): float(value)
            for header, value in zip(headers, row, strict=False)
        }
        for row in np.asarray(values, dtype=float)
    ]


def _fraction_mapping_from_amounts(
    amounts: np.ndarray,
    endmember_names: Sequence[str],
) -> dict[str, float]:
    if np.any(amounts < 0.0):
        raise InputConditionError("Endmember/component amounts must be non-negative.")
    total = float(np.sum(amounts))
    if total <= 0.0:
        raise InputConditionError("Endmember/component amounts must sum positive.")
    return dict(
        zip(
            [str(name).strip() for name in endmember_names],
            (amounts / total).tolist(),
            strict=False,
        )
    )


def _amount_to_moles(
    value: float,
    amount_unit: str,
    stoich: np.ndarray,
    atomic_masses: Mapping[str, float],
    species: str,
) -> float:
    amount = float(value)
    if amount < 0.0:
        raise InputConditionError(f"Negative amount supplied for {species!r}.")
    if not _is_mass_unit(amount_unit):
        return amount
    molar_mass = _stoich_molar_mass(stoich, atomic_masses)
    if molar_mass <= 0.0:
        if amount == 0.0:
            return 0.0
        raise InputConditionError(
            f"Cannot convert mass amount for zero-mass species {species!r}."
        )
    return amount / molar_mass


def _stoich_molar_mass(
    stoich: np.ndarray,
    atomic_masses: Mapping[str, float],
) -> float:
    total = 0.0
    for element, coefficient in zip(atomic_masses, stoich, strict=False):
        total += float(coefficient) * float(atomic_masses[element])
    return total


def _map_components_to_endmembers(
    database: dict,
    species_keys: Sequence[str],
    element_names: Sequence[str],
    endmember_names: Sequence[str],
    endmember_stoich: np.ndarray,
) -> dict[str, int] | None:
    used: set[int] = set()
    mapping: dict[str, int] = {}
    for species in species_keys:
        try:
            component_stoich = _component_stoich_vector(database, species)
        except ValueError:
            return None
        match = _matching_endmember_index(
            component_stoich,
            endmember_stoich,
            used,
        )
        if match is None:
            return None
        used.add(match)
        mapping[str(species).strip()] = match
    if len(mapping) != len(endmember_names):
        return None
    return mapping


def _component_stoich_vector(database: dict, species: str) -> np.ndarray:
    stoich_map = formula_or_database_species_stoichiometry(species, database)
    element_names = [str(name).strip() for name in database["cElementNameCS"]]
    return np.asarray(
        [float(stoich_map.get(element, 0.0)) for element in element_names],
        dtype=float,
    )


def _matching_endmember_index(
    component_stoich: np.ndarray,
    endmember_stoich: np.ndarray,
    used: set[int],
) -> int | None:
    component_shape = _normalized_stoich_shape(component_stoich)
    if component_shape is None:
        return None
    matches: list[int] = []
    for index, row in enumerate(endmember_stoich):
        if index in used:
            continue
        endmember_shape = _normalized_stoich_shape(row)
        if endmember_shape is None:
            continue
        if np.allclose(component_shape, endmember_shape, atol=1e-10, rtol=1e-10):
            matches.append(index)
    if len(matches) != 1:
        return None
    return matches[0]


def _normalized_stoich_shape(stoich: np.ndarray) -> np.ndarray | None:
    clean = np.clip(np.asarray(stoich, dtype=float), 0.0, None)
    total = float(np.sum(clean))
    if total <= 0.0:
        return None
    return clean / total


def _endmember_fraction_rows(
    endmember_fractions: Mapping[str, float]
    | Sequence[Mapping[str, float]]
    | None,
    count: int,
) -> list[Mapping[str, float] | None]:
    if endmember_fractions is None:
        return [None] * count
    if isinstance(endmember_fractions, Mapping):
        return [endmember_fractions] * count
    if len(endmember_fractions) != count:
        raise InputConditionError(
            "endmember_fractions must be one mapping or one mapping per condition."
        )
    return list(endmember_fractions)


def _property_solution_point(
    database: dict,
    phase_name: str,
    condition: dict,
    unit: list,
    *,
    endmember_fractions: Mapping[str, float] | None,
    include_heat_capacity: bool,
) -> SlnPropertyPoint:
    _preprocess_single(
        database,
        condition,
        unit,
        phases=[phase_name],
        include_heat_capacity=include_heat_capacity,
    )
    _run_solution_property_initialization()

    system_phase_index = _selected_solution_phase_index(phase_name)
    first, last = _solution_species_bounds(system_phase_index)
    endmember_names = _system_endmember_names(first, last)

    if endmember_fractions is None:
        try:
            fractions = _infer_explicit_endmember_fractions(first, last)
        except _ExplicitCompositionUnavailable:
            _raise_if_target_outside_endmember_range(
                _effective_stoich(first, last),
                _target_element_fraction(),
            )
            (
                system_phase_index,
                first,
                last,
                endmember_names,
                fractions,
            ) = _minimize_phase_local_model_view(
                database,
                phase_name,
                condition,
                unit,
                include_heat_capacity,
            )
    else:
        fractions = _normalize_endmember_fractions(endmember_fractions, endmember_names)

    _evaluate_solution_phase(system_phase_index, first, last, fractions)
    return _build_property_point(system_phase_index, first, last, endmember_names)


def _run_standard_property_initialization() -> None:
    for stage in (
        "checkthermoinput",
        "initthermo",
        "checksystem",
        "compthermodata",
        "checkthermodata",
    ):
        try:
            getattr(fort, stage)()
        except Exception as exc:
            info = fort.modulethermoio.infothermo
            raise EquilibError(
                f"Equilifort {stage} failed: infothermo={info}, Error: {exc}"
        ) from exc
        _raise_fortran_status(stage)


def _run_solution_property_initialization() -> None:
    for stage in (
        "checkthermoinput",
        "initthermo",
        "checksystem",
        "compthermodata",
    ):
        try:
            getattr(fort, stage)()
        except Exception as exc:
            info = fort.modulethermoio.infothermo
            raise EquilibError(
                f"Equilifort {stage} failed: infothermo={info}, Error: {exc}"
            ) from exc
        _raise_fortran_status(stage)

    try:
        fort.initgemsolver()
    except Exception as exc:
        info = fort.modulethermoio.infothermo
        raise EquilibError(
            f"Equilifort initgemsolver failed: infothermo={info}, Error: {exc}"
        ) from exc
    _raise_fortran_status("initgemsolver")


def _selected_solution_phase_index(phase_name: str) -> int:
    selected_names = [str(name).strip() for name in var.cPhaseNameSys]
    if phase_name in selected_names:
        index = selected_names.index(phase_name)
    else:
        prefix = f"{phase_name}#"
        matches = [
            i for i, name in enumerate(selected_names) if name.startswith(prefix)
        ]
        if not matches:
            selected_text = ", ".join(selected_names)
            raise InputConditionError(
                f"Selected phase {phase_name} is not a solution phase in this "
                f"system. Selected phases: {selected_text}"
            )
        index = matches[0]

    phase_type, _ = var.PhaseNameSys[selected_names[index]]
    if phase_type != "soln":
        raise InputConditionError(
            f"property_solution only supports solution phases; "
            f"{phase_name} is {phase_type}."
        )
    return index


def _solution_species_bounds(system_phase_index: int) -> tuple[int, int]:
    boundaries = np.asarray(fort.modulethermo.nspeciesphase, dtype=int)
    return int(boundaries[system_phase_index]), int(boundaries[system_phase_index + 1])


def _system_endmember_names(first: int, last: int) -> list[str]:
    database_indices = list(var.iSys2DBSpecies[first:last])
    return [str(var.cEndmemberNameCS[index]).strip() for index in database_indices]


def _infer_explicit_endmember_fractions(first: int, last: int) -> np.ndarray:
    stoich = _effective_stoich(first, last)
    element_amounts = np.asarray(fort.modulethermo.dmoleselement, dtype=float).copy()
    element_amounts = element_amounts[: int(fort.modulethermo.nelements)]
    if np.sum(element_amounts) <= 0.0:
        raise InputConditionError("At least one component amount must be non-zero.")
    target = element_amounts / np.sum(element_amounts)

    try:
        return _solve_endmember_fractions_from_elements(stoich, target)
    except _ExplicitCompositionUnavailable as exc:
        raise exc


def _effective_stoich(first: int, last: int) -> np.ndarray:
    stoich = np.asarray(fort.modulethermo.dstoichspecies, dtype=float).copy()
    particles = np.asarray(fort.modulethermo.iparticlespermole, dtype=float).copy()
    return stoich[first:last, :] / particles[first:last, None]


def _solve_endmember_fractions_from_elements(
    stoich: np.ndarray,
    target: np.ndarray,
) -> np.ndarray:
    n_endmembers, n_elements = stoich.shape
    if n_endmembers != n_elements:
        raise _ExplicitCompositionUnavailable()
    if n_endmembers == 1:
        return np.array([1.0], dtype=float)

    totals = np.sum(stoich, axis=1)
    rows = [
        stoich[:, j] - target[j] * totals
        for j in range(n_elements - 1)
    ]
    rows.append(np.ones(n_endmembers))
    rhs = np.zeros(n_endmembers)
    rhs[-1] = 1.0

    try:
        fractions = np.linalg.solve(np.vstack(rows), rhs)
    except np.linalg.LinAlgError as exc:
        raise _ExplicitCompositionUnavailable() from exc

    fractions[np.isclose(fractions, 0.0, atol=1e-12)] = 0.0
    if np.any(fractions < -1e-10):
        raise _ExplicitCompositionUnavailable()
    fractions = np.clip(fractions, 0.0, None)
    total = float(np.sum(fractions))
    if total <= 0.0:
        raise _ExplicitCompositionUnavailable()
    fractions = fractions / total

    reconstructed = fractions @ stoich
    reconstructed_total = float(np.sum(reconstructed))
    if reconstructed_total <= 0.0:
        raise _ExplicitCompositionUnavailable()
    reconstructed = reconstructed / reconstructed_total
    if not np.allclose(reconstructed, target, atol=1e-8, rtol=1e-8):
        raise _ExplicitCompositionUnavailable()
    return fractions


def _minimize_phase_local_model_view(
    database: dict,
    phase_name: str,
    condition: dict,
    unit: list,
    include_heat_capacity: bool,
) -> tuple[int, int, int, list[str], np.ndarray]:
    """Minimize with the standard stack, then report the requested model view."""
    fort.resetthermo()
    _preprocess_single(
        database,
        condition,
        unit,
        phases=[phase_name],
        include_heat_capacity=include_heat_capacity,
    )
    minimize()

    system_phase_index = _selected_solution_phase_index(phase_name)
    first, last = _solution_species_bounds(system_phase_index)
    endmember_names = _system_endmember_names(first, last)
    stoich = _effective_stoich(first, last)
    target = _target_element_fraction()
    _raise_if_target_outside_endmember_range(stoich, target)

    active_phase_indices = _active_solution_phase_indices()
    if system_phase_index in active_phase_indices:
        fractions = np.asarray(fort.modulethermo.dmolfraction, dtype=float).copy()
        fractions = fractions[first:last]
    else:
        companion_indices = [
            active_phase_index
            for active_phase_index in active_phase_indices
            if _is_order_disorder_companion(system_phase_index, active_phase_index)
        ]
        if len(companion_indices) != 1:
            active_names = [
                str(var.cPhaseNameSys[index]).strip()
                for index in active_phase_indices
            ]
            active_text = ", ".join(active_names) if active_names else "none"
            raise EquilibError(
                f"Selected phase {phase_name} cannot represent the requested "
                "composition on the standard minimization path; active "
                f"phase(s): {active_text}."
            )
        fractions = _requested_phase_random_state_from_companion(
            companion_indices[0],
            endmember_names,
        )

    total = float(np.sum(fractions))
    if total <= 0.0:
        raise EquilibError("Standard minimization produced zero phase fraction.")
    fractions = fractions / total
    composition = _composition_from_fractions(fractions, stoich)
    if composition is None:
        raise EquilibError("Standard minimization result has zero composition.")
    error = float(np.max(np.abs(composition - target)))
    if error > _COMPOSITION_TOL:
        raise EquilibError(
            "Standard minimization could not satisfy the selected-phase "
            f"requested composition; maximum atom-fraction residual is {error:.3e}."
        )
    return system_phase_index, first, last, endmember_names, fractions


def _active_solution_phase_indices() -> list[int]:
    """Return active solution phase indices from the converged assemblage."""
    mt = fort.modulethermo
    nelements = int(mt.nelements)
    assemblage = np.asarray(mt.iassemblage, dtype=int).copy()[:nelements]
    amounts = np.asarray(mt.dmolesphase, dtype=float).copy()[:nelements]
    return [
        -int(entry) - 1
        for entry, amount in zip(assemblage, amounts, strict=False)
        if int(entry) < 0 and float(amount) > 1e-12
    ]


def _is_order_disorder_companion(
    requested_phase_index: int,
    active_phase_index: int | None,
) -> bool:
    """Return True when the active phase is the requested phase's companion."""
    if active_phase_index is None:
        return False
    requested_db_id = int(var.iSys2DBSoln[requested_phase_index])
    active_db_id = int(var.iSys2DBSoln[active_phase_index])
    companion_ids = _requested_disordered_companion_db_ids(requested_db_id)
    if active_db_id in companion_ids:
        return True

    active_name = str(var.cPhaseNameSys[active_phase_index]).strip().upper()
    companion_names = {
        str(var.cSolnPhaseNameCS[index - 1]).strip().upper()
        for index in companion_ids
        if index > 0
    }
    companion_names.update(
        DISORDERED_PHASE_CANONICAL_NAMES.get(name, name).upper()
        for name in list(companion_names)
    )
    return active_name in companion_names


def _requested_disordered_companion_db_ids(requested_db_id: int) -> set[int]:
    """Return helper and canonical disordered database ids for an ordered phase."""
    disordered = np.asarray(var.iDisorderedPhaseCS, dtype=int)
    if requested_db_id < 1 or requested_db_id > disordered.size:
        return set()
    helper_id = int(disordered[requested_db_id - 1])
    if helper_id <= 0:
        return set()
    ids = {helper_id}
    helper_name = str(var.cSolnPhaseNameCS[helper_id - 1]).strip().upper()
    canonical_name = DISORDERED_PHASE_CANONICAL_NAMES.get(helper_name)
    if canonical_name is not None:
        for index, name in enumerate(var.cSolnPhaseNameCS, start=1):
            if str(name).strip().upper() == canonical_name.upper():
                ids.add(index)
    return ids


def _requested_phase_random_state_from_companion(
    active_phase_index: int,
    requested_endmember_names: list[str],
) -> np.ndarray:
    """Project a converged disordered branch onto the requested ordered model."""
    first, last = _solution_species_bounds(active_phase_index)
    fractions = np.asarray(fort.modulethermo.dmolfraction, dtype=float).copy()
    companion_x = fractions[first:last]
    companion_total = float(np.sum(companion_x))
    if companion_total <= 0.0:
        raise EquilibError("Companion phase has zero endmember fraction.")
    companion_x = companion_x / companion_total
    companion_stoich = _effective_stoich(first, last)
    element_fraction = _composition_from_fractions(companion_x, companion_stoich)
    if element_fraction is None:
        raise EquilibError("Companion phase has zero elemental composition.")
    return _ordered_random_endmember_fractions(
        requested_endmember_names,
        _selected_element_names(),
        element_fraction,
    )


def _ordered_random_endmember_fractions(
    endmember_names: Sequence[str],
    element_names: Sequence[str],
    element_fraction: np.ndarray,
) -> np.ndarray:
    """Build ordered endmember fractions on the disordered/random manifold."""
    parts = _endmember_sublattice_parts(endmember_names)
    if parts is None:
        raise EquilibError(
            "Requested companion mapping requires sublattice endmembers."
        )
    constituents_by_site = _sublattice_constituents(parts)
    element_map = {
        str(name).strip().upper(): float(value)
        for name, value in zip(element_names, element_fraction, strict=False)
    }

    fractions = np.zeros(len(parts), dtype=float)
    for index, endmember in enumerate(parts):
        product = 1.0
        for site, constituent in enumerate(endmember):
            name = str(constituent).strip().upper()
            if name == "VA":
                if len(constituents_by_site[site]) > 1:
                    product = 0.0
                    break
                continue
            if name not in element_map:
                product = 0.0
                break
            product *= max(element_map[name], 0.0)
        fractions[index] = product

    total = float(np.sum(fractions))
    if total <= 0.0:
        raise EquilibError(
            "Companion phase could not be projected onto the requested ordered model."
        )
    return fractions / total


def _raise_if_target_outside_endmember_range(
    stoich: np.ndarray,
    target: np.ndarray,
) -> None:
    row_totals = np.sum(stoich, axis=1)
    valid_rows = row_totals > 0.0
    if not np.any(valid_rows):
        raise InputConditionError(
            "Selected phase has no endmember with positive elemental content."
        )

    endmember_compositions = stoich[valid_rows] / row_totals[valid_rows, None]
    lower = np.min(endmember_compositions, axis=0) - _COMPOSITION_TOL
    upper = np.max(endmember_compositions, axis=0) + _COMPOSITION_TOL
    if np.any((target < lower) | (target > upper)):
        raise InputConditionError(
            "Requested element composition is outside the composition range "
            "spanned by the selected phase endmembers."
        )


def _target_element_fraction() -> np.ndarray:
    element_amounts = np.asarray(fort.modulethermo.dmoleselement, dtype=float).copy()
    element_amounts = element_amounts[: int(fort.modulethermo.nelements)]
    total = float(np.sum(element_amounts))
    if total <= 0.0:
        raise InputConditionError("At least one component amount must be non-zero.")
    return element_amounts / total


def _endmember_sublattice_parts(
    endmember_names: Sequence[str],
) -> list[tuple[str, ...]] | None:
    parts = [tuple(name.split(":")) for name in endmember_names]
    if not parts or len(parts[0]) <= 1:
        return None
    width = len(parts[0])
    if any(len(row) != width for row in parts):
        return None
    return parts


def _sublattice_constituents(
    parts: Sequence[tuple[str, ...]],
) -> list[list[str]]:
    constituents: list[list[str]] = []
    for sublattice in range(len(parts[0])):
        names = sorted({row[sublattice] for row in parts})
        constituents.append(names)
    return constituents


def _selected_element_names() -> list[str]:
    return [
        str(name).strip()
        for name in np.asarray(var.cElementNameCS)[var.iElementDBIndex]
    ]


def _composition_from_fractions(
    fractions: np.ndarray,
    stoich: np.ndarray,
) -> np.ndarray | None:
    amounts = np.asarray(fractions, dtype=float) @ stoich
    total = float(np.sum(amounts))
    if total <= 0.0:
        return None
    return amounts / total


def _normalize_endmember_fractions(
    endmember_fractions: Mapping[str, float],
    endmember_names: list[str],
) -> np.ndarray:
    exact_names = {name: i for i, name in enumerate(endmember_names)}
    upper_names = {name.upper(): i for i, name in enumerate(endmember_names)}
    fractions = np.zeros(len(endmember_names), dtype=float)

    for raw_name, raw_value in endmember_fractions.items():
        name = str(raw_name).strip()
        value = float(raw_value)
        index = exact_names.get(name)
        if index is None:
            index = upper_names.get(name.upper())
        if index is None:
            if abs(value) <= 1e-15:
                continue
            allowed = ", ".join(endmember_names)
            raise InputConditionError(
                f"Unknown endmember {raw_name!r} for selected phase. "
                f"Available endmembers: {allowed}"
            )
        fractions[index] = value

    if np.any(fractions < 0.0):
        raise InputConditionError("Endmember fractions must be non-negative.")
    total = float(np.sum(fractions))
    if total <= 0.0:
        raise InputConditionError("Endmember fractions must sum to a positive value.")
    return fractions / total


def _evaluate_solution_phase(
    system_phase_index: int,
    first: int,
    last: int,
    fractions: np.ndarray,
) -> None:
    fort.modulethermo.dmolfraction[first:last] = fractions
    fort.modulethermo.dmolesspecies[first:last] = fractions
    fort.modulethermo.dgibbssolnphase[system_phase_index] = 0.0
    try:
        fort.compexcessgibbsenergy(system_phase_index + 1)
    except Exception as exc:
        info = fort.modulethermoio.infothermo
        raise EquilibError(
            "Equilifort compexcessgibbsenergy failed: "
            f"infothermo={info}, Error: {exc}"
        ) from exc
    _raise_fortran_status("compexcessgibbsenergy")


def _build_property_point(
    system_phase_index: int,
    first: int,
    last: int,
    endmember_names: list[str],
) -> SlnPropertyPoint:
    ideal_constant = float(fort.modulethermo.didealconstant)
    temperature = float(fort.modulethermoio.dtemperature)
    energy_factor = ideal_constant * temperature
    species_slice = slice(first, last)

    fractions = np.asarray(fort.modulethermo.dmolfraction, dtype=float).copy()
    xi = fractions[species_slice]
    xi = xi / np.sum(xi)
    atom_denominator = _phase_occupied_atom_count(first, last, xi)

    partial_gibbs_values = (
        np.asarray(fort.modulethermo.dchemicalpotential, dtype=float).copy()
        * energy_factor
    )
    standard_gibbs_energies = (
        np.asarray(fort.modulethermo.dstdgibbsenergy, dtype=float).copy()
        * energy_factor
    )
    partial_enthalpies = (
        np.asarray(fort.modulethermo.dpartialenthalpy, dtype=float).copy()
        * energy_factor
    )
    partial_entropies = (
        np.asarray(fort.modulethermo.dpartialentropy, dtype=float).copy()
        * ideal_constant
    )
    partial_heat_capacities = (
        np.asarray(fort.modulethermo.dpartialheatcapacity, dtype=float).copy()
        * ideal_constant
    )
    activities = np.exp(
        np.clip(
            partial_gibbs_values[species_slice] / energy_factor
            - standard_gibbs_energies[species_slice] / energy_factor,
            -745.0,
            709.0,
        )
    )

    phase_name = str(var.cPhaseNameSys[system_phase_index]).strip()
    element_names = [
        str(name).strip()
        for name in np.asarray(var.cElementNameCS)[var.iElementDBIndex]
    ]
    endmembers_x = dict(zip(endmember_names, xi.tolist(), strict=False))
    elements_x = endmembers2elements(phase_name, endmembers_x, element_names)
    elements_w = endmembers2elements(
        phase_name,
        endmembers_x,
        element_names,
        unit_out=["K", "atm", "g"],
    )
    endmembers_w = _endmember_mass_fractions(first, last, xi)

    phase_G = _phase_weighted_property(
        xi,
        partial_gibbs_values[species_slice],
    )
    phase_H = _phase_weighted_property(
        xi,
        partial_enthalpies[species_slice],
    )
    phase_S = _phase_weighted_property(
        xi,
        partial_entropies[species_slice],
    )
    phase_Cp = _phase_weighted_property(
        xi,
        partial_heat_capacities[species_slice],
    )
    phase_G /= atom_denominator
    phase_H /= atom_denominator
    phase_S /= atom_denominator
    phase_Cp /= atom_denominator
    reported_partial_gibbs = partial_gibbs_values[species_slice] / atom_denominator
    reported_partial_enthalpies = partial_enthalpies[species_slice] / atom_denominator
    reported_partial_entropies = partial_entropies[species_slice] / atom_denominator
    reported_partial_heat_capacities = (
        partial_heat_capacities[species_slice] / atom_denominator
    )
    return SlnPropertyPoint(
        T=float(fort.modulethermoio.dtemperature),
        P=float(fort.modulethermoio.dpressure),
        phase_name=phase_name,
        G=phase_G,
        H=phase_H,
        S=phase_S,
        Cp=phase_Cp,
        n_i=_system_component_amounts(fort.modulethermo.dmoleselement),
        w_i=_system_component_mass_amounts(),
        endmembers_x=endmembers_x,
        endmembers_w=endmembers_w,
        elements_x=elements_x,
        elements_w=elements_w,
        partial_gibbs=SpeciesProperty(
            endmember_names,
            reported_partial_gibbs,
        ),
        partial_enthalpy=SpeciesProperty(
            endmember_names,
            reported_partial_enthalpies,
        ),
        partial_entropy=SpeciesProperty(
            endmember_names,
            reported_partial_entropies,
        ),
        activity=SpeciesProperty(endmember_names, activities),
        partial_heat_capacity=SpeciesProperty(
            endmember_names,
            reported_partial_heat_capacities,
        ),
    )


def _phase_occupied_atom_count(first: int, last: int, xi: np.ndarray) -> float:
    nelements = int(fort.modulethermo.nelements)
    stoich = np.asarray(fort.modulethermo.dstoichspecies, dtype=float).copy()
    occupied = xi @ stoich[first:last, :nelements]
    total = float(np.sum(occupied))
    if total <= 0.0 or not np.isfinite(total):
        return 1.0
    return total


def _endmember_mass_fractions(
    first: int,
    last: int,
    xi: np.ndarray,
) -> dict[str, float]:
    stoich = np.asarray(fort.modulethermo.dstoichspecies, dtype=float).copy()
    atomic_masses = np.asarray(fort.modulethermo.datomicmass, dtype=float).copy()
    endmember_masses = stoich[first:last, :] @ atomic_masses
    weighted = xi * endmember_masses
    total = float(np.sum(weighted))
    names = _system_endmember_names(first, last)
    if total <= 0.0:
        return dict(zip(names, np.zeros(len(names)), strict=False))
    return dict(zip(names, (weighted / total).tolist(), strict=False))


def _system_component_amounts(values: Any) -> dict[str, float]:
    component_names = [str(name).strip() for name in var.cComponentNameSys]
    component_values = np.asarray(values, dtype=float).copy()[var.iElementSysIndex]
    return dict(zip(component_names, component_values.tolist(), strict=False))


def _system_component_mass_amounts() -> dict[str, float]:
    component_names = [str(name).strip() for name in var.cComponentNameSys]
    mole_values = np.asarray(fort.modulethermo.dmoleselement, dtype=float).copy()
    mole_values = mole_values[var.iElementSysIndex]
    mass_values = [
        float(moles) * float(var.cPeriodicTable[name][1])
        for name, moles in zip(component_names, mole_values, strict=False)
    ]
    return dict(zip(component_names, mass_values, strict=False))
