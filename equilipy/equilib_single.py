"""Single-condition equilibrium calculation helpers."""

from __future__ import annotations

import numpy as np

import equilipy.equilifort as fort
import equilipy.variables as var

from .composition import expand_condition_species
from .exceptions import EquilibError, InputConditionError
from .input_condition import input_condition
from .list_phases import list_phases
from .load_database import load_fortran_database
from .minimize import minimize
from .phase_selection import phase_selection
from .results import capture_result_context
from .results.equilib import EquilibResult
from .utils import _dict2np


def _preprocess_single(
    database: dict,
    condition,
    unit: list = None,
    phases=None,
    include_heat_capacity=True,
):
    if unit is None:
        unit = ["K", "atm", "moles"]

    condition = expand_condition_species(database, condition, unit[2])
    load_fortran_database(database)
    # Get info from input condition
    NPTheader, NPTvals = _dict2np(condition)

    NPTvals = np.squeeze(NPTvals)

    # Check zero input
    composition_condition = NPTvals
    comp = np.array(NPTvals[2:])
    phase_elements = list(NPTheader[2:])
    phase_selector = _validated_phase_selection
    if np.any(np.isclose(comp, 0.0)):
        elements = np.array(NPTheader[2:])
        active_mask = ~np.isclose(comp, 0.0)
        phase_elements = list(elements[active_mask])
        if not phase_elements:
            raise InputConditionError(
                "At least one component amount must be non-zero."
            )
        composition_condition = list(comp[active_mask])
        composition_condition = np.array(list(NPTvals[:2]) + composition_condition)
        phase_selector = _active_phase_selection

    all_phases = list_phases(database, phase_elements)
    if phases is not None:
        phases = phase_selector(phases, all_phases)
        phase_selection(phases)

    var.dConditionSys = composition_condition
    input_condition(unit, composition_condition, include_heat_capacity)

    return condition


def _validated_phase_selection(
    requested_phases,
    available_phases,
) -> list[str]:
    if not isinstance(requested_phases, list):
        raise InputConditionError("Make sure if the data type of phases is list.")
    unknown = [
        phase
        for phase in requested_phases
        if not _requested_phase_is_available(str(phase), available_phases)
    ]
    if unknown:
        unknown_text = ", ".join(str(phase) for phase in unknown)
        raise InputConditionError(
            f"Unknown phase(s) for selected system: {unknown_text}"
        )
    return [
        str(phase)
        for phase in requested_phases
        if _requested_phase_is_available(str(phase), available_phases)
    ]


def _active_phase_selection(
    requested_phases,
    available_phases,
) -> list[str]:
    """Return requested phases that remain available for active non-zero elements."""
    if not isinstance(requested_phases, list):
        raise InputConditionError("Make sure if the data type of phases is list.")
    selected = [
        str(phase)
        for phase in requested_phases
        if _requested_phase_is_available(str(phase), available_phases)
    ]
    if not selected:
        raise InputConditionError(
            "No selected phases are available for the active non-zero composition."
        )
    return selected


def _requested_phase_is_available(phase: str, available_phases) -> bool:
    """Return True when a requested phase or explicit ``#1`` alias is available."""
    if phase in available_phases:
        return True
    if "#" not in phase:
        return False
    base_phase, composition_set = phase.rsplit("#", 1)
    return composition_set == "1" and base_phase in available_phases


def _equilib_single(
    database: dict,
    condition,
    unit: list = None,
    phases=None,
    include_heat_capacity=True,
):
    if unit is None:
        unit = ["K", "atm", "moles"]

    _preprocess_single(
        database,
        condition,
        unit,
        phases=phases,
        include_heat_capacity=include_heat_capacity,
    )

    try:
        minimize()
    except Exception:
        fort.resetthermo()
        raise

    return None


def equilib_single(
    database: dict,
    condition,
    unit: list = None,
    phases=None,
    include_heat_capacity=True,
):
    """Calculate equilibrium for one condition.

    Description
    ===========
    This function conducts equilibrium calculations for one condition.


    Revisions
    =========

     Date            Programmer      Description of change
     ----            ----------      ---------------------
     12/29/2023      S.Y. KWON       Original code


    Variables
    =========

    Input
    datafile : A string of 'Directory/databasename.dat'
    units    : A list of units for temperature, pressure, mass e.g.['K','atm','moles'].
                   Temperature units, 'K'/'C'/'F'/'R'.
                   Pressure units, 'atm'/'psi'/'bar'/'Pa'/'kPa'.
                   Mass units, 'mass fraction'/'kilograms'/'grams'/'pounds',
                   'mole fraction'/'atom fraction'/'atoms'/'moles'.
    Components: list of strings containing element name[element1, element2 ...]
    condition : [T, P, element1 amount, element1 amount ...]

    Output
    Results dataclass
    """
    if unit is None:
        unit = ["K", "atm", "moles"]

    context_condition = expand_condition_species(database, condition, unit[2])
    res = EquilibResult(
        context=capture_result_context(
            database,
            context_condition,
            unit,
            phases,
        )
    )
    # As a default, synchronize input and output units.
    _preprocess_single(
        database,
        condition,
        unit,
        phases,
        include_heat_capacity,
    )

    try:
        minimize()

        res.append_output()
    except EquilibError as error:
        if not res.append_postprocess_warning(error):
            res.append_error(error)

    fort.resetthermo()

    return res
