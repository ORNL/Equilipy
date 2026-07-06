"""Stoichiometric species and endmember thermodynamic property evaluation."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np

import equilipy.equilifort as fort
import equilipy.variables as var

from .composition import species_match_key
from .equilib_single import _preprocess_single
from .exceptions import InputConditionError
from .property_solution import _run_standard_property_initialization
from .results import capture_result_context
from .results.equilib import (
    EquilibPoint,
    EquilibResult,
    PhaseResult,
    StablePhaseSummaryDict,
)


@dataclass(frozen=True)
class _SpeciesRecord:
    index: int
    name: str
    element_names: tuple[str, ...]
    stoichiometry: np.ndarray
    atomic_masses: np.ndarray

    @property
    def element_amounts(self) -> dict[str, float]:
        return {
            element: float(amount)
            for element, amount in zip(
                self.element_names,
                self.stoichiometry,
                strict=False,
            )
            if abs(float(amount)) > 0.0
        }

    @property
    def element_mass_amounts(self) -> dict[str, float]:
        return {
            element: float(amount * atomic_mass)
            for element, amount, atomic_mass in zip(
                self.element_names,
                self.stoichiometry,
                self.atomic_masses,
                strict=False,
            )
            if abs(float(amount)) > 0.0
        }

    @property
    def molar_mass(self) -> float:
        return float(np.dot(self.stoichiometry, self.atomic_masses))


def property_compound(
    database: dict,
    compound_name: str,
    temperature: float | Sequence[float],
    *,
    pressure: float | Sequence[float] = 1.0,
    unit: list | None = None,
    include_heat_capacity: bool = True,
) -> EquilibResult:
    """Evaluate one stoichiometric compound species or solution endmember.

    The composition is fixed to one mole of the requested database species.
    Result ``n_i`` therefore reports the species stoichiometry, for example
    one mole of ``Fe2Si`` gives ``{"Fe": 2.0, "Si": 1.0}``.
    """
    if unit is None:
        unit = ["K", "atm", "moles"]
    if len(unit) < 3 or unit[2] not in {"moles", "mol"}:
        raise InputConditionError(
            "property_compound fixes composition as one mole of species; "
            "the amount unit must be 'moles'."
        )

    species = _resolve_species_record(database, compound_name)
    rows = _temperature_pressure_rows(temperature, pressure)
    context_condition = _compound_condition(rows, species)
    result = EquilibResult(
        context=capture_result_context(
            database,
            context_condition,
            unit,
            [species.name],
        )
    )

    points: list[EquilibPoint] = []
    for row_temperature, row_pressure in rows:
        condition = {
            "T": row_temperature,
            "P": row_pressure,
            **species.element_amounts,
        }
        try:
            points.append(
                _property_compound_point(
                    database,
                    species,
                    condition,
                    unit,
                    include_heat_capacity=include_heat_capacity,
                )
            )
        finally:
            fort.resetthermo()

    result.data = points[0] if len(points) == 1 else points
    return result


def _temperature_pressure_rows(
    temperature: float | Sequence[float],
    pressure: float | Sequence[float],
) -> list[tuple[float, float]]:
    temperatures = _as_float_array(temperature)
    pressures = _as_float_array(pressure)
    if len(pressures) == 1 and len(temperatures) > 1:
        pressures = np.full(len(temperatures), pressures[0], dtype=float)
    if len(temperatures) != len(pressures):
        raise InputConditionError(
            "temperature and pressure must have the same length, or pressure "
            "must be scalar."
        )
    return [
        (float(row_temperature), float(row_pressure))
        for row_temperature, row_pressure in zip(
            temperatures,
            pressures,
            strict=True,
        )
    ]


def _as_float_array(value: float | Sequence[float]) -> np.ndarray:
    if isinstance(value, (str, bytes)):
        return np.asarray([float(value)], dtype=float)
    try:
        return np.asarray(list(value), dtype=float)
    except TypeError:
        return np.asarray([float(value)], dtype=float)


def _compound_condition(
    rows: list[tuple[float, float]],
    species: _SpeciesRecord,
) -> dict[str, Any]:
    temperatures = np.asarray([row[0] for row in rows], dtype=float)
    pressures = np.asarray([row[1] for row in rows], dtype=float)
    condition: dict[str, Any] = {"T": temperatures, "P": pressures}
    for element, amount in species.element_amounts.items():
        condition[element] = np.full(len(rows), amount, dtype=float)
    if len(rows) == 1:
        return {key: np.asarray(value).item() for key, value in condition.items()}
    return condition


def _property_compound_point(
    database: dict,
    species: _SpeciesRecord,
    condition: dict[str, float],
    unit: list,
    *,
    include_heat_capacity: bool,
) -> EquilibPoint:
    _preprocess_single(
        database,
        condition,
        unit,
        phases=None,
        include_heat_capacity=include_heat_capacity,
    )
    _run_standard_property_initialization()
    system_index = _system_species_index(species)
    return _build_compound_point(species, system_index)


def _system_species_index(species: _SpeciesRecord) -> int:
    try:
        return list(var.iSys2DBSpecies).index(species.index)
    except ValueError as exc:
        raise InputConditionError(
            f"Species {species.name} is not available for its element set."
        ) from exc


def _build_compound_point(
    species: _SpeciesRecord,
    system_index: int,
) -> EquilibPoint:
    ideal_constant = float(fort.modulethermo.didealconstant)
    temperature = float(fort.modulethermoio.dtemperature)
    energy_factor = ideal_constant * temperature

    G = float(fort.modulethermo.dstdgibbsenergy[system_index] * energy_factor)
    H = float(fort.modulethermo.dstdenthalpy[system_index] * energy_factor)
    S = float(fort.modulethermo.dstdentropy[system_index] * ideal_constant)
    Cp = float(fort.modulethermo.dstdheatcapacity[system_index] * ideal_constant)

    elements_x = _normalized_mapping(species.element_amounts)
    elements_w = _normalized_mapping(species.element_mass_amounts)
    species_property = {species.name: G}
    partial_h = {species.name: H}
    partial_s = {species.name: S}
    partial_cp = {species.name: Cp}
    phase = PhaseResult(
        id=system_index + 1,
        name=species.name,
        amount_n=1.0,
        amount_w=species.molar_mass,
        stability=1.0,
        endmembers_x={species.name: 1.0},
        endmembers_w={species.name: 1.0},
        elements_x=elements_x,
        elements_w=elements_w,
        G=G,
        H=H,
        S=S,
        Cp=Cp,
        partial_gibbs=species_property,
        standard_gibbs_energy=species_property,
        activity={species.name: 1.0},
        partial_enthalpy=partial_h,
        partial_entropy=partial_s,
        partial_heat_capacity=partial_cp,
    )

    return EquilibPoint(
        T=float(fort.modulethermoio.dtemperature),
        P=float(fort.modulethermoio.dpressure),
        G=G,
        H=H,
        S=S,
        Cp=Cp,
        n_i=species.element_amounts,
        w_i=species.element_mass_amounts,
        stable_phase_summary=StablePhaseSummaryDict(
            {
                "id": np.array([phase.id]),
                "name": np.array([phase.name]),
                "amount_n": np.array([phase.amount_n]),
                "amount_w": np.array([phase.amount_w]),
            }
        ),
        phase_map={phase.name: phase},
    )


def _normalized_mapping(values: dict[str, float]) -> dict[str, float]:
    total = float(sum(values.values()))
    if total <= 0.0:
        return {key: 0.0 for key in values}
    return {key: float(value) / total for key, value in values.items()}


def _resolve_species_record(database: dict, compound_name: str) -> _SpeciesRecord:
    matches = _matching_species_indices(database, compound_name)
    if not matches:
        raise InputConditionError(
            f"Unknown stoichiometric species or endmember: {compound_name}"
        )

    equivalent = _equivalent_species_indices(database, matches)
    if len(equivalent) > 1:
        names = ", ".join(_species_name(database, index) for index in equivalent)
        raise InputConditionError(
            f"Name {compound_name!r} matches multiple non-equivalent species: "
            f"{names}"
        )

    index = equivalent[0]
    element_names = tuple(str(name).strip() for name in database["cElementNameCS"])
    stoichiometry = np.asarray(database["dStoichSpeciesCS"][index], dtype=float)
    atomic_masses = np.asarray(database["dAtomicMass"], dtype=float)
    return _SpeciesRecord(
        index=index,
        name=_species_name(database, index),
        element_names=element_names,
        stoichiometry=stoichiometry,
        atomic_masses=atomic_masses,
    )


def _matching_species_indices(database: dict, compound_name: str) -> list[int]:
    target = _match_key(compound_name)
    names = [str(name).strip() for name in database["cSpeciesNameCS"]]

    exact = [index for index, name in enumerate(names) if _match_key(name) == target]
    if exact:
        return exact

    return [
        index
        for index, name in enumerate(names)
        if target in _species_lookup_keys(name)
    ]


def _species_lookup_keys(raw_name: str) -> set[str]:
    cleaned = raw_name.strip().strip("'")
    first_token = cleaned.split("_", 1)[0].strip("'")
    no_phase_suffix = cleaned.split("(", 1)[0].strip("'")
    first_no_suffix = first_token.split("(", 1)[0].strip("'")
    return {
        _match_key(candidate)
        for candidate in {
            cleaned,
            first_token,
            no_phase_suffix,
            first_no_suffix,
        }
        if candidate
    }


def _match_key(value: str) -> str:
    return species_match_key(value).upper()


def _equivalent_species_indices(database: dict, matches: list[int]) -> list[int]:
    reference = matches[0]
    equivalent = [reference]
    different = []
    for index in matches[1:]:
        if _species_rows_equal(database, reference, index):
            continue
        different.append(index)
    return equivalent + different


def _species_rows_equal(database: dict, left: int, right: int) -> bool:
    if int(database["iParticlesPerMoleCS"][left]) != int(
        database["iParticlesPerMoleCS"][right]
    ):
        return False
    if not np.allclose(
        database["dStoichSpeciesCS"][left],
        database["dStoichSpeciesCS"][right],
    ):
        return False
    if not np.array_equal(
        _gibbs_columns(database, left),
        _gibbs_columns(database, right),
    ):
        return False
    if "dGibbsMagneticCS" in database and not np.allclose(
        database["dGibbsMagneticCS"][left],
        database["dGibbsMagneticCS"][right],
    ):
        return False
    return True


def _gibbs_columns(database: dict, species_index: int) -> np.ndarray:
    counts = np.asarray(database["nGibbsEqSpecies"], dtype=int)
    start = int(np.sum(counts[:species_index]))
    stop = start + int(counts[species_index])
    return np.asarray(database["dGibbsCoeffSpeciesTemp"])[:, start:stop]


def _species_name(database: dict, index: int) -> str:
    return str(database["cSpeciesNameCS"][index]).strip()
