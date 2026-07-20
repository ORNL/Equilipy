"""Nucleoscheil solidification routines."""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Union

import numpy as np
from tqdm import tqdm

import equilipy.equilifort as fort
import equilipy.variables as var

from ._parallel import starmap_joblib
from .composition import expand_condition_species
from .equilib_single import _equilib_single, _preprocess_single, equilib_single
from .exceptions import (
    EquilibError,
    InputConditionError,
    TransitionError,
    liquidus_search_failure_warning,
)
from .find_transition import _phase_row_groups, find_liquidus_transition
from .list_phases import list_phases
from .minimize import _raise_fortran_status, prepare_minimization
from .phase_selection import default_scheil_phase_selection
from .results import capture_result_context
from .results.common import combine_dicts
from .results.equilib import (
    _PHASE_CONSTITUTION_ATTRIBUTES,
    _mass_weighted_phase_constitution,
    _phase_constitution,
    _phase_with_constitution,
)
from .results.scheil import ScheilResult
from .scheil_cooling import (
    _emit_cooling_warning,
    _warn_if_truncated_with_liquid,
)
from .utils import _dict2np
from .variables import cPeriodicTable

# ============================================================================
# Nucleoscheil functions
_DF_ZERO_TOLERANCE = 1e-6
_GROWTH_MANIFOLD_TOLERANCE = 1e-4


class _CandidateDFStatus(str, Enum):
    """Classification of one phase in a liquid-state driving-force screen."""

    FAVORABLE = "favorable"
    UNFAVORABLE = "unfavorable"
    UNKNOWN = "unknown"


@dataclass(frozen=True, slots=True)
class CandidateDFFact:
    """Driving-force and undercooling facts for one candidate phase."""

    phase: str
    status: _CandidateDFStatus
    value: float | None
    saturation_temperature: float | None = None
    undercooling: float | None = None
    source_state: tuple = ()


@dataclass(frozen=True, slots=True)
class NucleationEventOutcome:
    """Accepted or rejected liquid-plus-newcomer nucleation transaction."""

    newcomer: str
    event_temperature: float
    accepted_liquid_state: dict[str, float] | None
    deposited_increment: float
    status: str


@dataclass(frozen=True, slots=True)
class SaturationEpisode:
    temperature: float
    liquid_update_count: int
    source_state: tuple


@dataclass(frozen=True, slots=True)
class DepositionTransaction:
    """Material ledger entry for one accepted liquid-to-solid transaction."""

    kind: str
    temperature: float
    deposited_elements_w: dict[str, float]
    liquid_before_elements_w: dict[str, float]
    liquid_after_elements_w: dict[str, float]


@dataclass(slots=True)
class NucleationState:
    """Mutable state for one NucleoScheil cooling path."""

    eligible_universe: tuple[str, ...]
    memory: dict[str, NucleationEventOutcome] = field(default_factory=dict)
    deposited_phases: dict[str, float] = field(default_factory=dict)
    saturation_episodes: dict[str, SaturationEpisode] = field(default_factory=dict)
    event_ledger: list[NucleationEventOutcome] = field(default_factory=list)
    transaction_ledger: list[DepositionTransaction] = field(default_factory=list)
    df_screens: dict[tuple, dict[str, CandidateDFFact]] = field(default_factory=dict)
    liquid_update_count: int = 0


def _condition_element_masses(
    condition: Dict[str, Union[float, str]],
) -> dict[str, float]:
    """Return elemental masses for an internal mole-basis liquid condition."""
    return {
        element: float(amount) * cPeriodicTable[element.strip()][1]
        for element, amount in condition.items()
        if element not in {"T", "P"}
    }


def _phase_element_masses(phase) -> dict[str, float]:
    """Return elemental masses carried by one equilibrium phase result."""
    mass = float(phase.amount_w)
    return {
        element: mass * float(fraction)
        for element, fraction in phase.elements.w_i.items()
    }


def _record_deposition_transaction(
    state: NucleationState,
    kind: str,
    condition_before: Dict[str, Union[float, str]],
    condition_after: Dict[str, Union[float, str]] | None,
    phases: list,
) -> None:
    """Record accepted deposits before any table-row accumulation."""
    deposited: dict[str, float] = {}
    for phase in phases:
        for element, mass in _phase_element_masses(phase).items():
            deposited[element] = deposited.get(element, 0.0) + mass
    before = _condition_element_masses(condition_before)
    after = (
        _condition_element_masses(condition_after)
        if condition_after is not None
        else {element: 0.0 for element in before}
    )
    state.transaction_ledger.append(
        DepositionTransaction(
            kind=kind,
            temperature=float(condition_before["T"]),
            deposited_elements_w=deposited,
            liquid_before_elements_w=before,
            liquid_after_elements_w=after,
        )
    )


def _default_critical_undercooling(
    liquid_phase_name: str,
    database: Dict,
    condition: Dict[str, Union[float, str]],
    unit: Optional[List[str]] = None,
    phases: Optional[List[str]] = None,
) -> Dict[str, float]:
    """Return zero nucleation undercooling for all eligible solid phases."""
    hidden_helpers = {
        str(name).strip() for name in database.get("cOrderDisorderHelperPhaseNames", [])
    }
    if phases is not None:
        phase_names = list(phases)
    else:
        amount_unit = unit[2] if unit and len(unit) > 2 else "moles"
        try:
            expanded = expand_condition_species(database, condition, amount_unit)
            elements = list(expanded.keys())[2:]
        except ValueError:
            elements = [
                str(key)
                for key in condition
                if str(key).strip() and str(key) not in {"T", "P"}
            ]
        phase_names = list(list_phases(database, elements))
        phase_names, _excluded = default_scheil_phase_selection(
            database,
            phase_names,
        )
    return {
        str(phase): 0.0
        for phase in phase_names
        if (
            str(phase).strip()
            and str(phase).strip() != liquid_phase_name
            and str(phase).strip() not in hidden_helpers
            and not _is_reference_ser_phase(str(phase).strip())
        )
    }


def _liquid_state_key(
    condition: Dict[str, Union[float, str]],
    universe: tuple[str, ...],
) -> tuple:
    """Return an exact, order-independent cache key for a liquid state."""
    coordinates = tuple(
        sorted((str(name), float(value)) for name, value in condition.items())
    )
    return coordinates + (("eligible_universe", universe),)


def _screen_candidate_driving_forces(
    liquid_phase_name: str,
    database: Dict,
    condition: Dict[str, Union[float, str]],
    state: NucleationState,
    unit: Optional[List[str]] = None,
    warnings: Optional[List[str]] = None,
    universe: Optional[tuple[str, ...]] = None,
) -> dict[str, CandidateDFFact]:
    """Recompute one batched driving-force screen for the current liquid."""
    if unit is None:
        unit = ["K", "atm", "moles"]
    candidates = state.eligible_universe if universe is None else tuple(universe)
    key = (
        ("liquid_update_count", state.liquid_update_count),
    ) + _liquid_state_key(condition, state.eligible_universe)[:-1]
    stored = state.df_screens.get(key, {})
    if all(phase in stored for phase in candidates):
        return {phase: stored[phase] for phase in candidates}
    source_state = key + (("eligible_universe", candidates),)
    try:
        values = _batched_liquid_driving_force_values(
            liquid_phase_name,
            database,
            condition,
            candidates,
            unit,
        )
    except Exception as exc:
        reason = getattr(getattr(exc, "reason", None), "value", "unknown")
        if warnings is not None:
            warnings.append(
                "NucleoScheil DF screen UNKNOWN at "
                f"T={float(condition['T']):g}: reason={reason}; {exc}"
            )
        facts = {
            phase: CandidateDFFact(
                phase=phase,
                status=_CandidateDFStatus.UNKNOWN,
                value=None,
                source_state=source_state,
            )
            for phase in candidates
        }
        stored.update(facts)
        state.df_screens[key] = stored
        return facts

    facts: dict[str, CandidateDFFact] = {}
    for phase in candidates:
        raw_value = values.get(phase)
        if raw_value is None or not np.isfinite(raw_value):
            value = None
            status = _CandidateDFStatus.UNKNOWN
        else:
            value = float(raw_value)
            status = (
                _CandidateDFStatus.FAVORABLE
                if value <= _DF_ZERO_TOLERANCE
                else _CandidateDFStatus.UNFAVORABLE
            )
        facts[phase] = CandidateDFFact(
            phase=phase,
            status=status,
            value=value,
            source_state=source_state,
        )

    unknown = [phase for phase, fact in facts.items() if fact.status == "unknown"]
    if unknown and warnings is not None:
        warnings.append(
            "NucleoScheil DF screen UNKNOWN for phase(s) "
            + ", ".join(unknown)
            + f" at T={float(condition['T']):g}; no nucleation or growth inferred."
        )
    stored.update(facts)
    state.df_screens[key] = stored
    return facts


def _active_element_potential_map() -> dict[str, float]:
    """Capture the current element-potential plane by active element name."""
    count = int(fort.modulethermo.nelements)
    indices = np.asarray(var.iElementDBIndex, dtype=int)[:count]
    names = [str(var.cElementNameCS[int(index)]).strip() for index in indices]
    values = np.asarray(fort.modulethermo.delementpotential, dtype=float)[:count]
    return dict(zip(names, values, strict=True))


def _batched_liquid_driving_force_values(
    liquid_phase_name: str,
    database: Dict,
    condition: Dict[str, Union[float, str]],
    candidates: tuple[str, ...],
    unit: List[str],
) -> dict[str, float | None]:
    """Evaluate every candidate against one metastable-liquid potential plane."""
    try:
        _equilib_single(
            database,
            condition.copy(),
            unit=unit,
            phases=[liquid_phase_name],
            include_heat_capacity=False,
        )
        parent_potentials = _active_element_potential_map()
    finally:
        fort.resetthermo()

    try:
        _preprocess_single(
            database,
            condition.copy(),
            unit,
            phases=[liquid_phase_name, *candidates],
            include_heat_capacity=False,
        )
        prepare_minimization()
        fort.initgemsolver()
        _raise_fortran_status("initgemsolver")

        active_names = [
            str(var.cElementNameCS[int(index)]).strip()
            for index in np.asarray(var.iElementDBIndex, dtype=int)[
                : int(fort.modulethermo.nelements)
            ]
        ]
        missing = [name for name in active_names if name not in parent_potentials]
        if missing:
            raise InputConditionError(
                "Candidate DF screen requires element potential(s) absent from "
                f"the liquid parent: {', '.join(missing)}"
            )
        fort.modulethermo.delementpotential[:] = [
            parent_potentials[name] for name in active_names
        ]
        fort.compdrivingforceall()
        _raise_fortran_status("compdrivingforceall")
        return _candidate_values_from_fortran(candidates)
    finally:
        fort.resetthermo()


def _candidate_values_from_fortran(
    candidates: tuple[str, ...],
) -> dict[str, float | None]:
    """Read signed candidate minima from the completed Fortran DF sweep."""
    phase_names = [str(name).strip() for name in var.cPhaseNameSys]
    solution_count = int(fort.modulethermo.nsolnphasessys)
    solution_values = np.asarray(
        fort.modulegemsolver.ddrivingforcesoln,
        dtype=float,
    )
    solution_status = np.asarray(
        fort.modulegemsolver.isubmincandidatestatussoln,
        dtype=int,
    )
    unknown_statuses = {
        int(getattr(fort.modulegemsolver, "submin_candidate_unknown", 0)),
        int(getattr(fort.modulegemsolver, "submin_candidate_max_iter", 3)),
    }
    row_groups = _phase_row_groups()
    phase_potentials = np.asarray(fort.modulegemsolver.dphasepotential, dtype=float)
    values: dict[str, float | None] = {}
    for candidate in candidates:
        matching = [
            index
            for index, phase_name in enumerate(phase_names)
            if phase_name == candidate
            or phase_name.split("#", 1)[0] == candidate.split("#", 1)[0]
        ]
        candidate_values: list[float] = []
        unknown = False
        for index in matching:
            if index < solution_count:
                if index >= solution_values.size:
                    unknown = True
                    continue
                if (
                    index < solution_status.size
                    and int(solution_status[index]) in unknown_statuses
                ):
                    unknown = True
                    continue
                candidate_values.append(float(solution_values[index]))
                continue
            rows = row_groups.get(phase_names[index], ())
            candidate_values.extend(
                float(phase_potentials[row])
                for row in rows
                if 0 <= row < phase_potentials.size
            )
        if candidate_values:
            values[candidate] = min(candidate_values)
        elif unknown or matching:
            values[candidate] = None
        else:
            values[candidate] = None
    return values


def _is_reference_ser_phase(phase_name: str) -> bool:
    return phase_name.strip().endswith("_SER(s)")


def _temperature_delta_to_kelvin(value: float, unit: str) -> float:
    unit_key = unit.strip().lower()
    if unit_key in {"f", "fahrenheit"}:
        return value * 5.0 / 9.0
    if unit_key in {"r", "rankine"}:
        return value * 5.0 / 9.0
    return value


def _transition_search_floor(unit: Optional[List[str]]) -> float:
    """Return a conservative lower bound for transition searches."""
    temperature_unit = unit[0] if unit else "K"
    unit_key = temperature_unit.strip().lower()
    if unit_key in {"c", "celsius"}:
        return 25.0
    if unit_key in {"f", "fahrenheit"}:
        return 77.0
    if unit_key in {"r", "rankine"}:
        return 536.67
    return 298.15


def _critical_undercooling_to_kelvin(
    critical_undercooling: Dict[str, float],
    unit: Optional[List[str]],
) -> Dict[str, float]:
    temperature_unit = unit[0] if unit else "K"
    return {
        phase: _temperature_delta_to_kelvin(float(value), temperature_unit)
        for phase, value in critical_undercooling.items()
    }


def _init_nuclei(
    liquid_phase_name: str,
    database: Dict,
    condition: Dict[str, Union[float, str]],
    critical_undercooling: Dict[str, float],
    unit: Optional[List[str]] = None,
    warnings: Optional[List[str]] = None,
    state: Optional[NucleationState] = None,
) -> (float, str):
    """Locate the first criterion-satisfied phase from a favorable DF screen."""
    if unit is None:
        unit = ["K", "atm", "moles"]
    if state is None:
        state = NucleationState(tuple(critical_undercooling))
    facts = _screen_candidate_driving_forces(
        liquid_phase_name,
        database,
        condition,
        state,
        unit=unit,
        warnings=warnings,
    )
    events = _criterion_satisfied_candidates(
        liquid_phase_name,
        database,
        condition,
        critical_undercooling,
        state,
        facts,
        unit=unit,
        warnings=warnings,
    )
    if not events:
        raise EquilibError(
            "No phase in the current liquid state satisfies its nucleation criterion."
        )
    return events[0]


def _target_driving_force(
    liquid_phase_name: str,
    database: Dict,
    condition: Dict[str, Union[float, str]],
    phase: str,
    state: NucleationState,
    unit: List[str],
    warnings: Optional[List[str]],
) -> CandidateDFFact:
    """Evaluate one favorable target without a pairwise transition search."""
    return _screen_candidate_driving_forces(
        liquid_phase_name,
        database,
        condition,
        state,
        unit=unit,
        warnings=warnings,
        universe=(phase,),
    )[phase]


def _locate_saturation_temperature(
    liquid_phase_name: str,
    database: Dict,
    condition: Dict[str, Union[float, str]],
    phase: str,
    state: NucleationState,
    unit: List[str],
    warnings: Optional[List[str]],
    tolerance: float = 1e-4,
) -> float | None:
    """Locate ``DF=0`` upward using target-only driving-force evaluations."""
    cold = float(condition["T"])
    cold_fact = _target_driving_force(
        liquid_phase_name,
        database,
        condition,
        phase,
        state,
        unit,
        warnings,
    )
    if cold_fact.status != _CandidateDFStatus.FAVORABLE:
        return None

    hot_limit = max(cold + 100.0, 2.0 * max(cold, 1.0))
    step = max(10.0, 0.02 * max(abs(cold), 1.0))
    hot = cold
    hot_fact = cold_fact
    while hot < hot_limit and hot_fact.status == _CandidateDFStatus.FAVORABLE:
        hot = min(hot + step, hot_limit)
        hot_condition = condition.copy()
        hot_condition["T"] = hot
        hot_fact = _target_driving_force(
            liquid_phase_name,
            database,
            hot_condition,
            phase,
            state,
            unit,
            warnings,
        )
        step *= 1.6
    if hot_fact.status == _CandidateDFStatus.UNKNOWN:
        return None
    if hot_fact.status != _CandidateDFStatus.UNFAVORABLE:
        if warnings is not None:
            warnings.append(
                "NucleoScheil saturation location UNKNOWN for "
                f"{phase}: no DF>0 structural upper bound below T={hot_limit:g}."
            )
        return None

    for _ in range(80):
        if hot - cold <= tolerance:
            break
        middle = 0.5 * (hot + cold)
        middle_condition = condition.copy()
        middle_condition["T"] = middle
        middle_fact = _target_driving_force(
            liquid_phase_name,
            database,
            middle_condition,
            phase,
            state,
            unit,
            warnings,
        )
        if middle_fact.status == _CandidateDFStatus.UNKNOWN:
            return None
        if middle_fact.status == _CandidateDFStatus.FAVORABLE:
            cold = middle
        else:
            hot = middle
    # The favorable side is the accepted ``DF=0`` numerical coordinate.
    return cold


def _criterion_satisfied_candidates(
    liquid_phase_name: str,
    database: Dict,
    condition: Dict[str, Union[float, str]],
    critical_undercooling: Dict[str, float],
    state: NucleationState,
    facts: dict[str, CandidateDFFact],
    unit: List[str],
    warnings: Optional[List[str]] = None,
    excluded_phases: Optional[set[str]] = None,
) -> list[tuple[float, str]]:
    """Return reached nucleation coordinates in deterministic event order."""
    current_temperature = float(condition["T"])
    events: list[tuple[float, str]] = []
    excluded = excluded_phases or set()
    for phase in state.eligible_universe:
        if phase in excluded:
            continue
        fact = facts.get(phase)
        if fact is None or fact.status == _CandidateDFStatus.UNKNOWN:
            continue
        if fact.status == _CandidateDFStatus.UNFAVORABLE:
            state.saturation_episodes.pop(phase, None)
            continue
        episode = state.saturation_episodes.get(phase)
        if episode is None or episode.liquid_update_count != state.liquid_update_count:
            saturation = _locate_saturation_temperature(
                liquid_phase_name,
                database,
                condition,
                phase,
                state,
                unit,
                warnings,
            )
            if saturation is None:
                continue
            episode = SaturationEpisode(
                float(saturation), state.liquid_update_count, fact.source_state
            )
            state.saturation_episodes[phase] = episode
        saturation = episode.temperature
        facts[phase] = CandidateDFFact(
            phase=fact.phase,
            status=fact.status,
            value=fact.value,
            saturation_temperature=float(saturation),
            undercooling=float(saturation - current_temperature),
            source_state=fact.source_state,
        )
        event_temperature = saturation - critical_undercooling[phase]
        if critical_undercooling[phase] == 0.0:
            # Put a zero-undercooling event on the solid-bearing side of the
            # transition.  The transaction is still solved and reported at
            # this one native coordinate.
            event_temperature -= 0.01
        if current_temperature <= event_temperature + 1e-8:
            # A displaced remembered phase re-enters on the first cooling
            # coordinate that crosses its renewed kinetic condition.  New
            # phases retain the exactly located first-nucleation coordinate.
            coordinate = (
                current_temperature if phase in state.memory else event_temperature
            )
            events.append((float(coordinate), phase))
    order = {phase: index for index, phase in enumerate(state.eligible_universe)}
    return sorted(events, key=lambda item: (-item[0], order[item[1]]))


def _update_liquid_composition(
    liquid_phase_name: str,
    latest_scheil_result: ScheilResult,
    condition: Dict[str, Union[float, str]],
) -> Dict[str, float]:

    res = condition.copy()
    liquid_phase = latest_scheil_result.equilib_result.points[-1].phase(
        liquid_phase_name
    )
    liquid_composition = liquid_phase.elements.x_i
    liquid_amount = latest_scheil_result.phase_amount_n(liquid_phase_name)[-1]
    if liquid_amount <= 0:
        liquid_amount = 0.0
    for element, fraction in liquid_composition.items():
        res[element] = np.round(float(fraction) * liquid_amount, 10)
    return res


def _update_liquid_composition_from_equilib(
    liquid_phase_name: str,
    equilib_result,
    condition: Dict[str, Union[float, str]],
) -> Dict[str, float]:
    """Return the next liquid condition from one equilibrium result."""
    res = condition.copy()
    elements = list(condition.keys())[2:]
    liquid_phase = equilib_result.phase(liquid_phase_name)
    liquid_amount = max(float(liquid_phase.amount_n), 0.0)
    for element in elements:
        fraction = liquid_phase.elements.x_i.get(element, 0.0)
        # A precipitate cannot create an element in the remaining liquid.
        # Clamp solver-scale positive residuals before they are carried into
        # the next transaction; the deposited solid is balanced to this same
        # outgoing state below.
        res[element] = min(
            float(fraction) * liquid_amount,
            float(condition[element]),
        )
    return res


def _stable_solid_candidates(
    equilib_result,
    candidates: List[str],
    liquid_phase_name: str,
) -> List[str]:
    """Return candidate solid phases that are stable in one equilibrium result."""
    stable_names = set(equilib_result.stable_phases.names)
    stable_solids = []
    for phase in candidates:
        if phase == liquid_phase_name or phase not in stable_names:
            continue
        stable_solids.append(phase)
    return stable_solids


def _append_unique_phases(target: List[str], phases: List[str]) -> None:
    """Append phases to a list while preserving the first-seen order."""
    for phase in phases:
        if phase not in target:
            target.append(phase)


def _set_result_points(result, points: List) -> None:
    """Replace a result object's internal point storage."""
    if not points:
        result.data = None
    elif len(points) == 1:
        result.data = points[0]
    else:
        result.data = points


def _point_with_zero_amount_phase(point, phase_template):
    """Add one zero-amount phase to an equilibrium point for an event anchor."""
    phase = phase_template.model_copy(
        update={"amount_n": 0.0, "amount_w": 0.0, "stability": 0.0}
    )
    phase_map = dict(point.phase_map)
    phase_map[phase.name] = phase
    names = [
        name
        for name in (str(raw_name).strip() for raw_name in point.stable_phase_summary["name"])
        if name
    ]
    if phase.name not in names:
        names.append(phase.name)
    summary = {
        "id": np.asarray([phase_map[name].id for name in names]),
        "name": np.asarray(names),
        "amount_n": np.asarray([phase_map[name].amount_n for name in names]),
        "amount_w": np.asarray([phase_map[name].amount_w for name in names]),
    }
    return point.model_copy(
        update={"phase_map": phase_map, "stable_phase_summary": summary}
    )


def _insert_nucleation_anchor(
    result: ScheilResult,
    start_index: int,
    newcomer: str,
    liquid_anchor_point,
) -> bool:
    """Insert or replace the pre-deposition row for one first phase event."""
    deposit_point = result.equilib_result.points[start_index]
    phase_template = deposit_point.phase(newcomer)
    reuse_previous = (
        start_index > 0
        and np.isclose(result.points[start_index - 1].T, liquid_anchor_point.T)
    )
    if reuse_previous:
        anchor_equilib = _point_with_zero_amount_phase(
            result.equilib_result.points[start_index - 1],
            phase_template,
        )
        previous = result.points[start_index - 1]
        cumulative_n = dict(previous.cumulative_phases_n)
        cumulative_w = dict(previous.cumulative_phases_w)
        cumulative_n[newcomer] = 0.0
        cumulative_w[newcomer] = 0.0
        label = "+".join(
            phase
            for phase in (newcomer, result.liquid_phase_name)
            if phase
        )
        anchor_scheil = previous.model_copy(
            update={
                "label": label,
                "cumulative_phases_n": cumulative_n,
                "cumulative_phases_w": cumulative_w,
            }
        )
        scheil_points = result.points
        equilib_points = result.equilib_result.points
        scheil_points[start_index - 1] = anchor_scheil
        equilib_points[start_index - 1] = anchor_equilib
    else:
        previous = result.points[start_index - 1] if start_index > 0 else None
        if previous is not None:
            previous_equilib = result.equilib_result.points[start_index - 1]
            anchor_equilib = _point_with_zero_amount_phase(
                previous_equilib.model_copy(
                    update={"T": liquid_anchor_point.T, "P": liquid_anchor_point.P}
                ),
                phase_template,
            )
            cumulative_n = dict(previous.cumulative_phases_n)
            cumulative_w = dict(previous.cumulative_phases_w)
            cumulative_n[newcomer] = 0.0
            cumulative_w[newcomer] = 0.0
            label = "+".join(
                phase
                for phase in (newcomer, result.liquid_phase_name)
                if phase
            )
            anchor_scheil = previous.model_copy(
                update={
                    "T": liquid_anchor_point.T,
                    "P": liquid_anchor_point.P,
                    "label": label,
                    "cumulative_phases_n": cumulative_n,
                    "cumulative_phases_w": cumulative_w,
                }
            )
        else:
            anchor_equilib = _point_with_zero_amount_phase(
                liquid_anchor_point,
                phase_template,
            )
            anchor_scheil = result._scheil_point_from_equilib_step(
                anchor_equilib,
                previous_point=None,
            )
        scheil_points = result.points
        equilib_points = result.equilib_result.points
        scheil_points.insert(start_index, anchor_scheil)
        equilib_points.insert(start_index, anchor_equilib)
    _set_result_points(result, scheil_points)
    _set_result_points(result.equilib_result, equilib_points)
    return reuse_previous


def _contains_nucleation_anchor(points: List) -> bool:
    """Return whether adjacent same-temperature rows are anchor then deposit."""
    for anchor, deposit in zip(points, points[1:]):
        for phase, amount in anchor.cumulative_phases_n.items():
            if amount == 0.0 and deposit.cumulative_phases_n.get(phase, 0.0) > 0.0:
                return True
    return False


def _combined_cascade_label(
    solid_phases: List[str],
    liquid_phase_name: str,
    include_liquid: bool,
) -> str:
    label_phases = list(solid_phases)
    if include_liquid and liquid_phase_name not in label_phases:
        label_phases.append(liquid_phase_name)
    return "+".join(label_phases)


def _collapse_same_temperature_cascade(
    result: ScheilResult,
    start_index: int,
    solid_phases: List[str],
    liquid_phase_name: str,
    end_index: Optional[int] = None,
) -> None:
    """Collapse same-temperature sub-steps into one visible Scheil row."""
    scheil_points = result.points
    equilib_points = result.equilib_result.points
    if start_index >= len(scheil_points) or start_index >= len(equilib_points):
        return
    if end_index is None:
        end_index = len(scheil_points)

    event_scheil_points = scheil_points[start_index:end_index]
    event_equilib_points = equilib_points[start_index:end_index]
    if not event_scheil_points or not event_equilib_points:
        return
    last_scheil_point = event_scheil_points[-1]
    last_equilib_point = event_equilib_points[-1]
    include_liquid = last_scheil_point.fl > 1e-12 or last_scheil_point.fl_w > 1e-12
    label = _combined_cascade_label(solid_phases, liquid_phase_name, include_liquid)

    merged_phase_map = dict(last_equilib_point.phase_map)
    for point in event_equilib_points:
        for phase_name, phase in point.phase_map.items():
            if phase_name == liquid_phase_name:
                continue
            if phase_name in solid_phases:
                merged_phase_map[phase_name] = phase

    stable_names = [
        phase_name
        for phase_name in solid_phases + ([liquid_phase_name] if include_liquid else [])
        if phase_name in merged_phase_map
    ]
    stable_summary = {
        "id": np.asarray([merged_phase_map[name].id for name in stable_names]),
        "name": np.asarray(stable_names),
        "amount_n": np.asarray(
            [merged_phase_map[name].amount_n for name in stable_names]
        ),
        "amount_w": np.asarray(
            [merged_phase_map[name].amount_w for name in stable_names]
        ),
    }

    merged_equilib_point = last_equilib_point.model_copy(
        update={
            "phase_map": merged_phase_map,
            "stable_phase_summary": stable_summary,
        }
    )
    merged_scheil_point = last_scheil_point.model_copy(update={"label": label})

    _set_result_points(
        result,
        scheil_points[:start_index] + [merged_scheil_point] + scheil_points[end_index:],
    )
    _set_result_points(
        result.equilib_result,
        equilib_points[:start_index]
        + [merged_equilib_point]
        + equilib_points[end_index:],
    )


def _collapse_same_temperature_output(
    result: ScheilResult,
    liquid_phase_name: str,
) -> None:
    """Merge adjacent same-temperature cascade rows for table readability."""
    index = 0
    while index < len(result.points):
        points = result.points
        start_index = index
        temperature = points[start_index].T
        end_index = start_index + 1
        while end_index < len(points) and np.isclose(points[end_index].T, temperature):
            end_index += 1

        if end_index - start_index > 1:
            same_temperature_points = points[start_index:end_index]
            if _contains_nucleation_anchor(same_temperature_points):
                index = end_index
                continue
            solid_phases: List[str] = []
            # Prefer the final transaction's phase order. An earlier connector
            # row at the same temperature may contain only the incumbent and
            # must not reorder the newcomer-led terminal cascade.
            for point in reversed(same_temperature_points):
                phases = [
                    phase
                    for phase in point.label.split("+")
                    if phase and phase != liquid_phase_name
                ]
                if not phases:
                    continue
                _append_unique_phases(solid_phases, phases)
            if len(solid_phases) > 1:
                _collapse_same_temperature_cascade(
                    result,
                    start_index,
                    solid_phases,
                    liquid_phase_name,
                    end_index=end_index,
                )
                index = start_index + 1
            else:
                index = end_index
        else:
            index += 1


def _get_nuclei_candiates(
    liquid_phase_name: str,
    database: Dict,
    condition: Dict[str, Union[float, str]],
    critical_undercooling: Dict[str, float],
    unit: Optional[List[str]] = None,
    state: Optional[NucleationState] = None,
    warnings: Optional[List[str]] = None,
) -> Dict[str, float]:
    """Compatibility wrapper over the cached DF/undercooling candidate screen."""
    if unit is None:
        unit = ["K", "atm", "moles"]
    if state is None:
        state = NucleationState(tuple(critical_undercooling))
    facts = _screen_candidate_driving_forces(
        liquid_phase_name,
        database,
        condition,
        state,
        unit=unit,
        warnings=warnings,
    )
    events = _criterion_satisfied_candidates(
        liquid_phase_name,
        database,
        condition,
        critical_undercooling,
        state,
        facts,
        unit,
        warnings,
    )
    return {phase: temperature for temperature, phase in events}


def _get_same_temperature_nuclei_candidates(
    liquid_phase_name: str,
    database: Dict,
    condition: Dict[str, Union[float, str]],
    critical_undercooling: Dict[str, float],
    unit: Optional[List[str]] = None,
    excluded_phases: Optional[List[str]] = None,
    state: Optional[NucleationState] = None,
    facts: Optional[dict[str, CandidateDFFact]] = None,
    warnings: Optional[List[str]] = None,
) -> Dict[str, float]:
    """Return criterion-satisfied newcomers from the exact-state DF cache."""
    if unit is None:
        unit = ["K", "atm", "moles"]
    if state is None:
        state = NucleationState(tuple(critical_undercooling))
    excluded = set(excluded_phases or ())
    if facts is None:
        facts = _screen_candidate_driving_forces(
            liquid_phase_name,
            database,
            condition,
            state,
            unit=unit,
            warnings=warnings,
        )
    events = _criterion_satisfied_candidates(
        liquid_phase_name,
        database,
        condition,
        critical_undercooling,
        state,
        facts,
        unit,
        warnings,
        excluded_phases=excluded,
    )
    return {
        phase: temperature for temperature, phase in events if phase not in excluded
    }


def _derive_growth_set(
    state: NucleationState,
    facts: dict[str, CandidateDFFact],
) -> set[str]:
    """Return the sole, most-recent grower when it remains on its manifold.

    A successful nucleation transaction displaces every resident grower.  The
    event ledger is already persistent nucleation history, so its last accepted
    newcomer identifies the only phase allowed to march away from the event.
    Older remembered phases can become the last newcomer again only by meeting
    their freshly recomputed nucleation condition.
    """
    accepted = (
        event.newcomer
        for event in reversed(state.event_ledger)
        if event.status == "accepted"
    )
    phase = next(accepted, next(reversed(state.memory), None) if state.memory else None)
    fact = facts.get(phase) if phase is not None else None
    if (
        fact is None
        or fact.value is None
        or abs(fact.value) > _GROWTH_MANIFOLD_TOLERANCE
    ):
        return set()
    return {phase}


def _append_equilibrium_transaction(
    result: ScheilResult,
    liquid_phase_name: str,
    newcomer: str,
    database: Dict,
    condition: Dict[str, Union[float, str]],
    unit: List[str],
) -> tuple[dict[str, float] | None, float, bool]:
    """Apply the pairwise liquid-plus-solid lever on the current liquid."""
    scheil_count = len(result.points)
    equilib_count = len(result.equilib_result.points)
    try:
        _equilib_single(
            database,
            condition,
            unit=unit,
            phases=[liquid_phase_name, newcomer],
        )
        result.append_output()
    except Exception:
        return None, 0.0, False

    point = result.equilib_result.points[-1]
    stable_names = set(point.stable_phases.names)
    if newcomer not in stable_names:
        _set_result_points(result, result.points[:scheil_count])
        _set_result_points(
            result.equilib_result, result.equilib_result.points[:equilib_count]
        )
        return None, 0.0, False

    deposited = max(float(point.phase(newcomer).amount_n), 0.0)
    if liquid_phase_name not in stable_names:
        _balance_deposition_transaction(
            result,
            point,
            liquid_phase_name,
            [newcomer],
            condition,
            None,
        )
        return None, deposited, True
    liquid_state = _update_liquid_composition_from_equilib(
        liquid_phase_name,
        point,
        condition,
    )
    _balance_deposition_transaction(
        result,
        point,
        liquid_phase_name,
        [newcomer],
        condition,
        liquid_state,
    )
    return liquid_state, deposited, True


def _balance_deposition_transaction(
    result: ScheilResult,
    point,
    liquid_phase_name: str,
    solid_phases: list[str],
    condition_before: Dict[str, Union[float, str]],
    condition_after: Dict[str, Union[float, str]] | None,
):
    """Make reported deposits equal the liquid consumed, element by element.

    Equilibrium amounts carry solver-scale residuals. Banking both those raw
    solid amounts and the independently reconstructed outgoing liquid repeats
    that residual at every Scheil transaction. Preserve the thermodynamic
    liquid path and phase split, but put the residual into the largest
    deposited solid so the transaction is exactly material balanced.
    """
    before = _condition_element_masses(condition_before)
    after = (
        _condition_element_masses(condition_after)
        if condition_after is not None
        else {element: 0.0 for element in before}
    )
    consumed = {
        element: before[element] - after.get(element, 0.0)
        for element in before
    }
    present = [phase for phase in solid_phases if phase in point.phase_map]
    if not present:
        return point
    closure_phase = max(present, key=lambda phase: float(point.phase(phase).amount_w))
    phase_map = dict(point.phase_map)
    accounted = {element: 0.0 for element in consumed}
    for phase in present:
        if phase == closure_phase:
            continue
        for element, mass in _phase_element_masses(point.phase(phase)).items():
            if element in accounted:
                accounted[element] += mass
    closure_elements = {
        element: consumed[element] - accounted[element] for element in consumed
    }
    closure_mass = sum(closure_elements.values())
    if closure_mass <= 0.0 or any(
        mass < -1e-12 for mass in closure_elements.values()
    ):
        return point
    closure_elements = {element: max(mass, 0.0) for element, mass in closure_elements.items()}
    closure_mass = sum(closure_elements.values())
    closure = point.phase(closure_phase).model_copy(
        update={
            "amount_w": closure_mass,
            "elements_w": {
                element: mass / closure_mass if closure_mass > 0.0 else 0.0
                for element, mass in closure_elements.items()
            },
        }
    )
    phase_map[closure_phase] = closure
    summary = dict(point.stable_phase_summary)
    names = [
        name
        for name in (str(raw_name).strip() for raw_name in summary["name"])
        if name
    ]
    summary["id"] = np.asarray([phase_map[name].id for name in names])
    summary["name"] = np.asarray(names)
    summary["amount_n"] = np.asarray(
        [phase_map[name].amount_n for name in names], dtype=float
    )
    summary["amount_w"] = np.asarray(
        [phase_map[name].amount_w for name in names], dtype=float
    )
    point = point.model_copy(update={"phase_map": phase_map, "stable_phase_summary": summary})

    equilib_points = result.equilib_result.points
    _set_result_points(result.equilib_result, equilib_points[:-1] + [point])
    previous = result.points[-2] if len(result.points) > 1 else None
    row = result._scheil_point_from_equilib_step(point, previous_point=previous)
    _set_result_points(result, result.points[:-1] + [row])
    return point


def _process_newcomer(
    liquid_phase_name: str,
    result: ScheilResult,
    newcomer: str,
    database: Dict,
    condition: Dict[str, Union[float, str]],
    state: NucleationState,
    unit: List[str],
) -> NucleationEventOutcome:
    """Apply one liquid-plus-solid deposition transactionally."""
    first_deposit = newcomer not in state.memory
    event_temperature = float(condition["T"])
    liquid_state, deposited, accepted = _append_equilibrium_transaction(
        result,
        liquid_phase_name,
        newcomer,
        database,
        condition.copy(),
        unit,
    )
    outcome = NucleationEventOutcome(
        newcomer=newcomer,
        event_temperature=event_temperature,
        accepted_liquid_state=liquid_state,
        deposited_increment=deposited,
        status="accepted" if accepted else "solve_failed",
    )
    state.event_ledger.append(outcome)
    if not accepted:
        return outcome
    if first_deposit:
        state.memory[newcomer] = outcome
    state.deposited_phases[newcomer] = (
        state.deposited_phases.get(newcomer, 0.0) + deposited
    )
    if first_deposit:
        state.saturation_episodes.pop(newcomer, None)
    if hasattr(result, "equilib_result") and result.equilib_result.points:
        point = result.equilib_result.points[-1]
        _record_deposition_transaction(
            state,
            "event",
            condition,
            liquid_state,
            [point.phase(newcomer)],
        )
    if liquid_state is not None and any(
        liquid_state.get(coordinate) != value
        for coordinate, value in condition.items()
    ):
        state.liquid_update_count += 1
    return outcome


def _constant_temperature_candidates(
    state: NucleationState,
    previous_facts: dict[str, CandidateDFFact],
    facts: dict[str, CandidateDFFact],
    new_candidates: Dict[str, float],
    previous_phase: str,
    component_count: int,
    firing_count: int,
) -> list[str]:
    """Return phases whose current-liquid nucleation condition is satisfied."""
    candidates = list(new_candidates)
    alternatives = [phase for phase in candidates if phase != previous_phase]
    alternatives.sort(
        key=lambda phase: (
            facts[phase].value is None,
            facts[phase].value if facts[phase].value is not None else np.inf,
            state.eligible_universe.index(phase),
        )
    )
    return alternatives


def _collapse_constant_temperature_rows(
    result: ScheilResult,
    start_index: int,
    deposited_phases: List[str],
    liquid_phase_name: str,
    deposited_masses: dict[str, float],
    deposited_element_masses: dict[str, dict[str, float]],
) -> None:
    """Expose one accumulated row for a constant-temperature transaction loop."""
    if len(result.points) - start_index <= 1:
        return
    _collapse_same_temperature_cascade(
        result,
        start_index,
        deposited_phases,
        liquid_phase_name,
    )
    point = result.equilib_result.points[-1]
    phase_map = dict(point.phase_map)
    for phase, mass in deposited_masses.items():
        if mass <= 0.0 or phase not in phase_map:
            continue
        elements_w = {
            element: element_mass / mass
            for element, element_mass in deposited_element_masses[phase].items()
        }
        phase_map[phase] = phase_map[phase].model_copy(
            update={"amount_w": mass, "elements_w": elements_w}
        )
    summary = dict(point.stable_phase_summary)
    names = [
        name
        for name in (str(raw_name).strip() for raw_name in summary["name"])
        if name
    ]
    summary["id"] = np.asarray([phase_map[name].id for name in names])
    summary["name"] = np.asarray(names)
    summary["amount_n"] = np.asarray(
        [phase_map[name].amount_n for name in names], dtype=float
    )
    summary["amount_w"] = np.asarray(
        [phase_map[name].amount_w for name in names], dtype=float
    )
    point = point.model_copy(
        update={"phase_map": phase_map, "stable_phase_summary": summary}
    )
    _set_result_points(
        result.equilib_result,
        result.equilib_result.points[:-1] + [point],
    )


def _advance_constant_temperature(
    liquid_phase_name: str,
    latest_scheil_result: ScheilResult,
    nuclei_current: List[str],
    database: Dict,
    condition: Dict[str, Union[float, str]],
    critical_undercooling: Dict[str, float],
    unit: Optional[List[str]] = None,
    state: Optional[NucleationState] = None,
    recheck_conditions: bool = True,
    allow_bare_newcomers: bool = False,
) -> tuple[ScheilResult, bool]:
    """Alternate pairwise deposits at one temperature until the loop settles."""
    if unit is None:
        unit = ["K", "atm", "moles"]
    res = latest_scheil_result
    if state is None:
        state = NucleationState(tuple(critical_undercooling))
    liquid_condition = condition.copy()
    queue = list(nuclei_current[:1])
    bare_newcomer_commotion = allow_bare_newcomers or len(nuclei_current) > 1
    start_index = len(res.points)
    deposited_phases: List[str] = []
    deposited_masses: dict[str, float] = {}
    deposited_element_masses: dict[str, dict[str, float]] = {}
    max_transactions = 10000
    transaction_count = 0
    previous_facts = _screen_candidate_driving_forces(
        liquid_phase_name,
        database,
        liquid_condition,
        state,
        unit=unit,
        warnings=res.warnings,
    )
    while queue:
        transaction_count += 1
        if transaction_count > max_transactions:
            raise EquilibError(
                "NucleoScheil constant-temperature loop exceeded "
                f"{max_transactions} pairwise deposits at "
                f"T={float(liquid_condition['T']):g}."
            )
        newcomer = queue[0]
        outcome = _process_newcomer(
            liquid_phase_name,
            res,
            newcomer,
            database,
            liquid_condition,
            state,
            unit,
        )
        if outcome.status != "accepted":
            res.warnings.append(
                f"NucleoScheil newcomer solve failed for {newcomer} at "
                f"T={float(liquid_condition['T']):g}; the last accepted liquid "
                "state was preserved."
            )
            condition.clear()
            condition.update(liquid_condition)
            return res, False
        _append_unique_phases(deposited_phases, [newcomer])
        point = res.equilib_result.points[-1]
        deposited_phase = point.phase(newcomer)
        deposited_mass = float(deposited_phase.amount_w)
        deposited_masses[newcomer] = (
            deposited_masses.get(newcomer, 0.0) + deposited_mass
        )
        phase_element_masses = deposited_element_masses.setdefault(newcomer, {})
        for element, fraction in deposited_phase.elements.w_i.items():
            phase_element_masses[element] = (
                phase_element_masses.get(element, 0.0)
                + deposited_mass * float(fraction)
            )
        _collapse_constant_temperature_rows(
            res,
            start_index,
            deposited_phases,
            liquid_phase_name,
            deposited_masses,
            deposited_element_masses,
        )
        if outcome.accepted_liquid_state is None:
            _collapse_constant_temperature_rows(
                res,
                start_index,
                deposited_phases,
                liquid_phase_name,
                deposited_masses,
                deposited_element_masses,
            )
            return res, True
        liquid_unchanged = all(
            outcome.accepted_liquid_state.get(coordinate) == value
            for coordinate, value in liquid_condition.items()
        )
        liquid_condition = outcome.accepted_liquid_state
        component_count = max(len(liquid_condition) - 2, 1)
        if len(set(deposited_phases)) >= component_count:
            queue = []
            continue
        if outcome.deposited_increment <= 1e-12 and liquid_unchanged:
            queue = []
            continue
        if not recheck_conditions:
            queue = []
            continue
        facts = _screen_candidate_driving_forces(
            liquid_phase_name,
            database,
            liquid_condition,
            state,
            unit=unit,
            warnings=res.warnings,
        )
        candidates = _get_same_temperature_nuclei_candidates(
            liquid_phase_name,
            database,
            liquid_condition,
            critical_undercooling,
            unit=unit,
            state=state,
            facts=facts,
            warnings=res.warnings,
        )
        if bare_newcomer_commotion:
            for phase in state.eligible_universe:
                fact = facts.get(phase)
                if (
                    phase not in state.memory
                    and fact is not None
                    and fact.value is not None
                    and fact.value <= _GROWTH_MANIFOLD_TOLERANCE
                ):
                    candidates.setdefault(phase, float(liquid_condition["T"]))
        queue = _constant_temperature_candidates(
            state,
            previous_facts,
            facts,
            candidates,
            newcomer,
            max(len(liquid_condition) - 2, 1),
            len(set(deposited_phases)),
        )[:1]
        previous_facts = facts
    component_count = max(len(liquid_condition) - 2, 1)
    invariant_phases = list(dict.fromkeys(deposited_phases))
    if len(invariant_phases) >= component_count:
        try:
            _equilib_single(
                database,
                liquid_condition,
                unit=unit,
                phases=[liquid_phase_name, *invariant_phases],
            )
            res.append_output()
            terminal_point = res.equilib_result.points[-1]
            terminal_liquid = terminal_point.phase_map.get(liquid_phase_name)
            terminal = (
                terminal_liquid is None
                or float(terminal_liquid.amount_n) <= 1e-12
            )
            if terminal:
                terminal_point = _balance_deposition_transaction(
                    res,
                    terminal_point,
                    liquid_phase_name,
                    invariant_phases,
                    liquid_condition,
                    None,
                )
                terminal_solids = []
                for phase in invariant_phases:
                    terminal_phase = terminal_point.phase_map.get(phase)
                    if terminal_phase is None:
                        continue
                    terminal_solids.append(terminal_phase)
                    mass = float(terminal_phase.amount_w)
                    deposited_masses[phase] = deposited_masses.get(phase, 0.0) + mass
                    element_masses = deposited_element_masses.setdefault(phase, {})
                    for element, fraction in terminal_phase.elements.w_i.items():
                        element_masses[element] = (
                            element_masses.get(element, 0.0)
                            + mass * float(fraction)
                        )
                _record_deposition_transaction(
                    state,
                    "invariant_close",
                    liquid_condition,
                    None,
                    terminal_solids,
                )
            if terminal and terminal_liquid is not None and start_index > 0:
                reference_liquid = res.equilib_result.points[
                    start_index - 1
                ].phase_map.get(liquid_phase_name)
                if reference_liquid is not None:
                    phase_map = dict(terminal_point.phase_map)
                    phase_map[liquid_phase_name] = terminal_liquid.model_copy(
                        update={
                            "elements_x": dict(reference_liquid.elements_x),
                            "elements_w": dict(reference_liquid.elements_w),
                        }
                    )
                    terminal_point = terminal_point.model_copy(
                        update={"phase_map": phase_map}
                    )
                    equilib_points = res.equilib_result.points
                    equilib_points[-1] = terminal_point
                    _set_result_points(res.equilib_result, equilib_points)
        except Exception:
            terminal = False
        if terminal:
            _collapse_constant_temperature_rows(
                res,
                start_index,
                invariant_phases,
                liquid_phase_name,
                deposited_masses,
                deposited_element_masses,
            )
            return res, True
    _collapse_constant_temperature_rows(
        res,
        start_index,
        deposited_phases,
        liquid_phase_name,
        deposited_masses,
        deposited_element_masses,
    )
    condition.clear()
    condition.update(liquid_condition)
    return res, False


def _direct_endmember_constitution(
    endmembers: dict[str, float],
    elements: dict[str, float],
) -> dict[str, float] | None:
    """Map elemental fractions onto simple one-element endmembers."""
    mapped: dict[str, float] = {}
    matched_elements: set[str] = set()
    for endmember in endmembers:
        element = endmember.split(":", 1)[0]
        if element not in elements or element in matched_elements:
            return None
        matched_elements.add(element)
        mapped[endmember] = elements[element]
    if matched_elements != set(elements):
        return None
    return mapped


def _assemble_cumulative_phase_constitutions(
    result: ScheilResult,
    liquid_phase_name: str,
) -> None:
    """Align solid constitution columns with cumulative Scheil phase amounts."""
    cumulative_masses: dict[str, float] = {}
    weighted_sums: dict[str, dict[str, dict[str, float]]] = {}
    last_committed: dict[str, dict[str, dict[str, float]]] = {}
    last_committed_phase = {}
    previous_amounts: dict[str, float] = {}
    constitution_overrides = []

    for scheil_point, equilib_point in zip(
        result.points,
        result.equilib_result.points,
        strict=True,
    ):
        phase_map = dict(equilib_point.phase_map)
        row_overrides = {}
        terminal_row = scheil_point.cumulative_phases_w.get(
            liquid_phase_name, 0.0
        ) <= 1e-12
        for phase, cumulative_mass in scheil_point.cumulative_phases_w.items():
            if phase == liquid_phase_name or cumulative_mass <= 0.0:
                continue
            deposited_mass = max(
                float(cumulative_mass) - previous_amounts.get(phase, 0.0),
                0.0,
            )
            phase_result = phase_map.get(phase)
            if deposited_mass > 1e-12 and phase_result is not None:
                constitution = _phase_constitution(phase_result)
                last_committed[phase] = constitution
                last_committed_phase[phase] = phase_result
                phase_sums = weighted_sums.setdefault(
                    phase,
                    {
                        attribute: {}
                        for attribute in _PHASE_CONSTITUTION_ATTRIBUTES
                    },
                )
                for attribute, fractions in constitution.items():
                    attribute_sums = phase_sums[attribute]
                    for component, fraction in fractions.items():
                        attribute_sums[component] = (
                            attribute_sums.get(component, 0.0)
                            + deposited_mass * float(fraction)
                        )
                cumulative_masses[phase] = (
                    cumulative_masses.get(phase, 0.0) + deposited_mass
                )

            if phase_result is None and phase in last_committed_phase:
                phase_result = last_committed_phase[phase].model_copy(
                    update={"amount_n": 0.0, "amount_w": 0.0, "stability": 0.0}
                )
            if phase_result is None or phase not in last_committed:
                continue
            cumulative_constitution = _mass_weighted_phase_constitution(
                weighted_sums[phase],
                cumulative_masses[phase],
            )
            assembled_constitution = _phase_constitution(phase_result)
            assembled_constitution["elements_x"] = cumulative_constitution[
                "elements_x"
            ]
            assembled_constitution["elements_w"] = cumulative_constitution[
                "elements_w"
            ]
            if deposited_mass <= 1e-12:
                assembled_constitution["endmembers_x"] = last_committed[phase][
                    "endmembers_x"
                ]
                assembled_constitution["endmembers_w"] = last_committed[phase][
                    "endmembers_w"
                ]
            if terminal_row:
                terminal_endmembers_x = _direct_endmember_constitution(
                    cumulative_constitution["endmembers_x"],
                    cumulative_constitution["elements_x"],
                )
                terminal_endmembers_w = _direct_endmember_constitution(
                    cumulative_constitution["endmembers_w"],
                    cumulative_constitution["elements_w"],
                )
                if terminal_endmembers_x is not None:
                    assembled_constitution["endmembers_x"] = terminal_endmembers_x
                if terminal_endmembers_w is not None:
                    assembled_constitution["endmembers_w"] = terminal_endmembers_w
            row_overrides[phase] = _phase_with_constitution(
                phase_result,
                assembled_constitution,
            )

        previous_amounts = dict(scheil_point.cumulative_phases_w)
        constitution_overrides.append(row_overrides)

    if result.points and result.points[-1].fl_w <= 1e-12:
        initial_liquid = result.equilib_result.points[0].phase_map.get(
            liquid_phase_name
        )
        terminal_overrides = constitution_overrides[-1]
        terminal_amounts = {
            phase: float(amount)
            for phase, amount in result.points[-1].cumulative_phases_w.items()
            if phase != liquid_phase_name and amount > 0.0
        }
        if initial_liquid is not None and terminal_amounts:
            mass_units = {
                "grams", "kilograms", "pounds", "g", "kg", "lbs",
                "mass fraction", "weight fraction", "wt%", "wt.%",
            }
            if result.context.unit[2] in mass_units:
                target = {
                    element: float(result.context.input_condition[element])
                    for element in initial_liquid.elements.w_i
                }
            else:
                target = {
                    element: float(initial_liquid.amount_w) * float(fraction)
                    for element, fraction in initial_liquid.elements.w_i.items()
                }
            closure_phase = max(terminal_amounts, key=terminal_amounts.get)
            accounted = {element: 0.0 for element in target}
            prior_amounts: dict[str, float] = {}
            for scheil_point, equilib_point in zip(
                result.points[:-1],
                result.equilib_result.points[:-1],
                strict=True,
            ):
                for phase, cumulative in scheil_point.cumulative_phases_w.items():
                    if phase == liquid_phase_name:
                        continue
                    increment = max(
                        float(cumulative) - prior_amounts.get(phase, 0.0),
                        0.0,
                    )
                    phase_result = equilib_point.phase_map.get(phase)
                    if increment > 0.0 and phase_result is not None:
                        for element in target:
                            accounted[element] += increment * float(
                                phase_result.elements.w_i.get(element, 0.0)
                            )
                prior_amounts = dict(scheil_point.cumulative_phases_w)
            terminal_increments = {
                phase: max(amount - prior_amounts.get(phase, 0.0), 0.0)
                for phase, amount in terminal_amounts.items()
            }
            for phase, increment in terminal_increments.items():
                if phase == closure_phase or increment <= 0.0:
                    continue
                phase_result = result.equilib_result.points[-1].phase_map.get(phase)
                if phase_result is None:
                    continue
                for element in target:
                    accounted[element] += increment * float(
                        phase_result.elements.w_i.get(element, 0.0)
                    )
            closure_result = result.equilib_result.points[-1].phase_map.get(
                closure_phase
            )
            if closure_result is not None:
                closure_mass = terminal_increments[closure_phase]
                constitution = _phase_constitution(closure_result)
                constitution["elements_w"] = {
                    element: (target[element] - accounted[element]) / closure_mass
                    for element in target
                }
                corrected_phase = _phase_with_constitution(
                    closure_result,
                    constitution,
                )
                terminal_point = result.equilib_result.points[-1]
                phase_map = dict(terminal_point.phase_map)
                phase_map[closure_phase] = corrected_phase
                terminal_point = terminal_point.model_copy(update={"phase_map": phase_map})
                points = result.equilib_result.points
                points[-1] = terminal_point
                _set_result_points(result.equilib_result, points)
                override_accounted = {element: 0.0 for element in target}
                for phase, amount in terminal_amounts.items():
                    phase_result = terminal_overrides.get(
                        phase,
                        result.equilib_result.points[-1].phase_map.get(phase),
                    )
                    if phase_result is None:
                        continue
                    for element in target:
                        override_accounted[element] += amount * float(
                            phase_result.elements.w_i.get(element, 0.0)
                        )
                override_result = terminal_overrides.get(
                    closure_phase,
                    corrected_phase,
                )
                override_constitution = _phase_constitution(override_result)
                override_constitution["elements_w"] = {
                    element: float(override_result.elements.w_i.get(element, 0.0))
                    + (target[element] - override_accounted[element])
                    / terminal_amounts[closure_phase]
                    for element in target
                }
                direct_endmembers_w = _direct_endmember_constitution(
                    override_constitution["endmembers_w"],
                    override_constitution["elements_w"],
                )
                if direct_endmembers_w is not None:
                    override_constitution["endmembers_w"] = direct_endmembers_w
                terminal_overrides[closure_phase] = _phase_with_constitution(
                    override_result,
                    override_constitution,
                )
    result.equilib_result._phase_constitution_overrides = constitution_overrides


def nucleoscheil_cooling(
    liquid_phase_name: str,
    database: Dict,
    condition: Dict[str, Union[float, str]],
    critical_undercooling: Optional[Dict[str, float]] = None,
    delta_T: float = 5.0,
    unit: Optional[List[str]] = None,
    phases: Optional[List[str]] = None,
    progress_callback=None,
    show_progress: bool = True,
) -> ScheilResult:
    """Run nucleation-aware Scheil cooling for one condition."""
    if unit is None:
        unit = ["K", "atm", "moles"]
    if delta_T <= 0:
        raise EquilibError("delta_T must be positive for Nucleoscheil cooling.")
    if critical_undercooling is None:
        critical_undercooling = _default_critical_undercooling(
            liquid_phase_name,
            database,
            condition,
            unit,
            phases,
        )
    normalized_condition = expand_condition_species(database, condition, unit[2])
    critical_undercooling = _critical_undercooling_to_kelvin(
        critical_undercooling,
        unit,
    )
    hidden_helpers = {
        str(name).strip() for name in database.get("cOrderDisorderHelperPhaseNames", [])
    }
    critical_undercooling = {
        str(phase).strip(): value
        for phase, value in critical_undercooling.items()
        if str(phase).strip()
        and str(phase).strip() != liquid_phase_name
        and str(phase).strip() not in hidden_helpers
        and not _is_reference_ser_phase(str(phase).strip())
    }

    # Get phases
    solid_phases = list(critical_undercooling.keys())
    if not solid_phases:
        raise EquilibError(
            "Nucleoscheil cooling requires at least one solid phase in "
            "critical_undercooling or phases."
        )
    result_context = capture_result_context(
        database,
        normalized_condition,
        unit,
        [liquid_phase_name] + solid_phases,
    )
    NPT_local = normalized_condition.copy()
    elements = list(normalized_condition.keys())[2:]
    T_current = NPT_local["T"]
    unit_local = unit.copy()

    # if it is mass based, convert to moles
    is_mass_based = unit[2] in [
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
    ]
    if is_mass_based:
        unit_local[2] = "moles"
        for el in elements:
            NPT_local[el] = NPT_local[el] / cPeriodicTable[el.strip()][1]

    # Step 1: initialize Nucleoscheil and get the liquidus.
    liquidus_warning = None
    try:
        T_init = find_liquidus_transition(
            database,
            NPT_local,
            liquid_phase_name,
            T_current,
            max(T_current * 0.1, _transition_search_floor(unit_local)),
            unit=unit_local,
            phases=[liquid_phase_name] + solid_phases,
        )
    except TransitionError as exc:
        T_init = T_current
        liquidus_warning = liquidus_search_failure_warning(exc)

    T_current = T_init + 0.1
    NPT_local["T"] = T_current
    _equilib_single(database, NPT_local, unit=unit_local, phases=[liquid_phase_name])
    res = ScheilResult(
        liquid_phase_name=liquid_phase_name,
        input_unit=unit,
        context=result_context,
    )
    if liquidus_warning is not None:
        _emit_cooling_warning(
            res,
            liquidus_warning,
            progress_callback=progress_callback,
        )
    # now all units are converted to K
    unit_local[0] = "K"
    T_current = res.T
    NPT_local["T"] = T_current

    state = NucleationState(tuple(solid_phases))
    truncated_due_to_failed_continuation = False
    iteration_count = int(max(T_current - 295.0, 0.0) / delta_T) + len(solid_phases)
    if progress_callback is not None:
        progress_callback(0, iteration_count, "Nucleoscheil Cooling")
    for step_index in tqdm(
        range(iteration_count),
        desc="Nucleoscheil Cooling",
        colour="#8060ff",
        ascii="░▒▓",
        bar_format="{l_bar}{bar:50}{r_bar}",
        disable=not show_progress or progress_callback is not None,
    ):
        if progress_callback is not None:
            progress_callback(step_index + 1, iteration_count, "Nucleoscheil Cooling")
        target_temperature = T_current - delta_T
        if target_temperature < 295.0:
            break
        target_liquid = NPT_local.copy()
        target_liquid["T"] = target_temperature
        current_facts = _screen_candidate_driving_forces(
            liquid_phase_name,
            database,
            NPT_local,
            state,
            unit=unit_local,
            warnings=res.warnings,
        )
        facts = _screen_candidate_driving_forces(
            liquid_phase_name,
            database,
            target_liquid,
            state,
            unit=unit_local,
            warnings=res.warnings,
        )
        grow = _derive_growth_set(state, current_facts)
        events = _criterion_satisfied_candidates(
            liquid_phase_name,
            database,
            target_liquid,
            critical_undercooling,
            state,
            facts,
            unit_local,
            res.warnings,
            excluded_phases=grow,
        )
        reachable_events = [event for event in events if event[0] <= T_current + 1e-8]
        if reachable_events:
            event_temperature, newcomer = reachable_events[0]
            first_nucleation = newcomer not in state.memory
            event_liquid = NPT_local.copy()
            event_liquid["T"] = event_temperature
            start_index = len(res.points)
            try:
                anchor_result = equilib_single(
                    database,
                    event_liquid,
                    unit=unit_local,
                    phases=[liquid_phase_name],
                    include_heat_capacity=False,
                )
                res, terminal = _advance_constant_temperature(
                    liquid_phase_name,
                    res,
                    [newcomer],
                    database,
                    event_liquid,
                    critical_undercooling,
                    unit=unit_local,
                    state=state,
                )
            except Exception as exc:
                res.warnings.append(
                    "NucleoScheil newcomer event failed at "
                    f"T={event_temperature:g}: {type(exc).__name__}: {exc}"
                )
                truncated_due_to_failed_continuation = True
                break
            if first_nucleation and state.event_ledger[-1].status == "accepted":
                _insert_nucleation_anchor(
                    res,
                    start_index,
                    newcomer,
                    anchor_result.point,
                )
            NPT_local = event_liquid
            T_current = event_temperature
            if state.event_ledger[-1].status != "accepted":
                truncated_due_to_failed_continuation = True
                break
            if terminal:
                break
            continue

        T_current = target_temperature
        NPT_local["T"] = T_current
        if not grow:
            if not state.memory:
                try:
                    _equilib_single(
                        database,
                        target_liquid,
                        unit=unit_local,
                        phases=[liquid_phase_name],
                    )
                    res.append_output()
                except Exception as exc:
                    res.warnings.append(
                        "NucleoScheil liquid-head solve failed at "
                        f"T={target_temperature:g}: {type(exc).__name__}: {exc}"
                    )
                    truncated_due_to_failed_continuation = True
                    break
            continue

        phase_selection = [
            liquid_phase_name,
            *(phase for phase in solid_phases if phase in grow),
        ]
        try:
            _equilib_single(
                database,
                target_liquid,
                unit=unit_local,
                phases=phase_selection,
            )
            res.append_output()
        except Exception as exc:
            res.warnings.append(
                "NucleoScheil continuing-growth solve failed at "
                f"T={target_temperature:g}: {type(exc).__name__}: {exc}"
            )
            truncated_due_to_failed_continuation = True
            break

        point = res.equilib_result.points[-1]
        stable_names = set(point.stable_phases.names)
        for phase in grow:
            if phase in stable_names:
                increment = max(float(point.phase(phase).amount_n), 0.0)
                state.deposited_phases[phase] = (
                    state.deposited_phases.get(phase, 0.0) + increment
                )
        if liquid_phase_name not in stable_names:
            break
        NPT_local = _update_liquid_composition_from_equilib(
            liquid_phase_name,
            point,
            target_liquid,
        )
        point = _balance_deposition_transaction(
            res,
            point,
            liquid_phase_name,
            list(grow),
            target_liquid,
            NPT_local,
        )
        _record_deposition_transaction(
            state,
            "growth",
            target_liquid,
            NPT_local,
            [point.phase(phase) for phase in grow if phase in stable_names],
        )
        state.liquid_update_count += 1
        post_growth_facts = _screen_candidate_driving_forces(
            liquid_phase_name,
            database,
            NPT_local,
            state,
            unit=unit_local,
            warnings=res.warnings,
        )
        same_temperature = _get_same_temperature_nuclei_candidates(
            liquid_phase_name,
            database,
            NPT_local,
            critical_undercooling,
            unit=unit_local,
            state=state,
            facts=post_growth_facts,
            warnings=res.warnings,
        )
        bare_newcomers = [
            phase
            for phase in solid_phases
            if phase not in state.memory
            and post_growth_facts.get(phase) is not None
            and post_growth_facts[phase].value is not None
            and post_growth_facts[phase].value <= _GROWTH_MANIFOLD_TOLERANCE
        ]
        bare_newcomers.sort(
            key=lambda phase: (
                post_growth_facts[phase].value,
                solid_phases.index(phase),
            )
        )
        # v0.3.2's multinewcomer commotion is a simultaneous crossing: two
        # previously unseen phases become pairwise stable after one resident's
        # lever step.  Treat the group at the same temperature; a lone phase
        # still pays its ordinary kinetic barrier.
        unseen = [phase for phase in solid_phases if phase not in state.memory]
        minimum_barrier = min(
            (critical_undercooling[phase] for phase in unseen),
            default=np.inf,
        )
        minimum_barrier_crossed = any(
            critical_undercooling[phase] <= minimum_barrier + 1e-12
            for phase in bare_newcomers
        )
        component_count = max(len(NPT_local) - 2, 1)
        allow_bare_newcomers = minimum_barrier_crossed and (
            len(state.memory) + 1 < component_count
        )
        if allow_bare_newcomers:
            same_temperature = {
                phase: float(NPT_local["T"]) for phase in bare_newcomers
            }
        if same_temperature:
            newcomers = list(same_temperature)
            newcomer = newcomers[0]
            first_nucleation = newcomer not in state.memory
            start_index = len(res.points)
            ledger_start = len(state.event_ledger)
            try:
                anchor_result = equilib_single(
                    database,
                    NPT_local,
                    unit=unit_local,
                    phases=[liquid_phase_name],
                    include_heat_capacity=False,
                )
                res, terminal = _advance_constant_temperature(
                    liquid_phase_name,
                    res,
                    newcomers,
                    database,
                    NPT_local,
                    critical_undercooling,
                    unit=unit_local,
                    state=state,
                    allow_bare_newcomers=allow_bare_newcomers,
                )
            except Exception as exc:
                res.warnings.append(
                    "NucleoScheil same-temperature event failed at "
                    f"T={T_current:g}: {type(exc).__name__}: {exc}"
                )
                truncated_due_to_failed_continuation = True
                break
            if allow_bare_newcomers:
                event_slice = state.event_ledger[ledger_start:]
                initial = next(
                    (event for event in event_slice if event.newcomer == newcomer),
                    None,
                )
                if initial is not None and event_slice[-1] is not initial:
                    state.event_ledger[ledger_start:] = [
                        *(event for event in event_slice if event is not initial),
                        initial,
                    ]
            if first_nucleation and state.event_ledger[-1].status == "accepted":
                _insert_nucleation_anchor(
                    res,
                    start_index,
                    newcomer,
                    anchor_result.point,
                )
            if state.event_ledger[-1].status != "accepted":
                truncated_due_to_failed_continuation = True
                break
            if terminal:
                break
    if truncated_due_to_failed_continuation:
        _warn_if_truncated_with_liquid(
            res,
            liquid_phase_name=liquid_phase_name,
            cooling_name="Nucleoscheil cooling",
            progress_callback=progress_callback,
        )
    _collapse_same_temperature_output(res, liquid_phase_name)
    _assemble_cumulative_phase_constitutions(res, liquid_phase_name)
    res.update_scheil_constituents()
    # Reorder
    header = ["task_id"]
    elements = list(normalized_condition.keys())[2:]
    header.append("T_Delta [K]")
    header.extend([f"{x} [sp-mol]" for x in elements])
    header.extend([f"{x} [g]" for x in elements])

    scheil_constituent_names = [
        key.split(" [sp-mol]")[0]
        for key in res.scheil_constituents.keys()
        if key not in header and "[sp-mol]" in key
    ]
    scheil_constituent_names = [key for key in scheil_constituent_names if key.strip()]

    primaries = [
        constituent
        for constituent in scheil_constituent_names
        if "+" not in constituent and constituent != ""
    ]
    eutectics = [
        constituent for constituent in scheil_constituent_names if "+" in constituent
    ]

    header.extend([f"{x} [sp-mol]" for x in primaries])
    header.extend([f"{x} [sp-mol]" for x in eutectics])
    header.extend([f"{x} [g]" for x in primaries])
    header.extend([f"{x} [g]" for x in eutectics])
    res.scheil_constituents = {
        key: res.scheil_constituents[key]
        for key in header
        if key in res.scheil_constituents
    }
    # Private diagnostics are intentionally excluded from result serialization.
    res.nucleation_state = state
    return res


def _nucleoscheil_batch_input(
    liquid_phase_name: str,
    database: Dict,
    condition: Dict[str, Union[float, str]],
    critical_undercooling: Dict[str, Union[float, str]],
    delta_T: float = 5.0,
    unit: Optional[List[str]] = None,
    n_per_batch: int = 1,
) -> List:
    """
    Split conditions into smaller batch dictionaries.

    The batches cover the same full input range.
    """
    if unit is None:
        unit = ["K", "atm", "moles"]

    res = []
    expanded_condition = expand_condition_species(database, condition, unit[2])
    header, NPT = _dict2np(expanded_condition)
    header_uc = list(critical_undercooling)
    undercooling_columns = []
    for phase_name in header_uc:
        value = critical_undercooling[phase_name]
        try:
            len(value)  # type: ignore[arg-type]
        except TypeError:
            value = [value]
        column = np.asarray(value, dtype=float)
        if column.ndim == 0:
            column = column.reshape(1)
        undercooling_columns.append(column)
    uc_target_length = max(len(column) for column in undercooling_columns)
    if any(len(column) not in {1, uc_target_length} for column in undercooling_columns):
        raise ValueError(
            "critical_undercooling columns must have equal length or scalar values."
        )
    UC = np.column_stack(
        [
            np.full(uc_target_length, column[0]) if len(column) == 1 else column
            for column in undercooling_columns
        ]
    )
    n, _m = NPT.shape
    n_uc, _m_uc = UC.shape
    if n_uc not in {1, n}:
        raise ValueError(
            "critical_undercooling columns must have one row or match the "
            "condition row count."
        )
    fullrange = [0, n]

    # Nucleoscheil workers run one solidification path per task; keep n_per_batch
    # for API compatibility but pass scalar row conditions to the solver.
    for i in range(fullrange[0], fullrange[1]):
        subcondition = dict({})
        subundercooling = dict({})
        for j, head in enumerate(header):
            subcondition[head] = NPT[i, j].item()
        uc_row = 0 if n_uc == 1 else i
        for j, head in enumerate(header_uc):
            subundercooling[head] = UC[uc_row, j].item()
        res.append(
            [liquid_phase_name, database, subcondition, subundercooling, delta_T, unit]
        )
    return res


def _nucleoscheil_batch(
    liquid_phase_name: str,
    database: Dict,
    condition: Dict[str, Union[float, str]],
    critical_undercooling: Dict[str, Union[float, str]],
    delta_T: float = 5.0,
    unit: Optional[List[str]] = None,
) -> Dict:
    if unit is None:
        unit = ["K", "atm", "moles"]

    try:
        res = nucleoscheil_cooling(
            liquid_phase_name,
            database,
            condition,
            critical_undercooling,
            delta_T,
            unit,
            show_progress=False,
        )
        return res.to_dict()
    except Exception as exc:
        return {"Error": f"Nucleoscheil cooling failed: {type(exc).__name__}: {exc}"}


def _nucleoscheil_constituent_batch(
    liquid_phase_name: str,
    database: Dict,
    condition: Dict[str, Union[float, str]],
    critical_undercooling: Dict[str, Union[float, str]],
    delta_T: float = 5.0,
    unit: Optional[List[str]] = None,
) -> Dict:
    if unit is None:
        unit = ["K", "atm", "moles"]

    try:
        res = nucleoscheil_cooling(
            liquid_phase_name,
            database,
            condition,
            critical_undercooling,
            delta_T,
            unit,
            show_progress=False,
        )
        return res.scheil_constituents.copy()
    except Exception as exc:
        return {
            "Error": (
                "Nucleoscheil constituent calculation failed: "
                f"{type(exc).__name__}: {exc}"
            )
        }


def _nucleoscheil_singlenode(arg, n_cpu, progress_callback=None) -> List[Dict]:
    return starmap_joblib(
        _nucleoscheil_batch,
        arg,
        n_cpu,
        progress=True,
        progress_callback=progress_callback,
        colour="#8060ff",
        ascii="░▒▓",
        bar_format="{l_bar}{bar:50}{r_bar}",
        desc="Nucleoscheil Batch",
    )


def _nucleoscheil_constituent_singlenode(
    arg,
    n_cpu,
    progress_callback=None,
) -> List[Dict]:
    return starmap_joblib(
        _nucleoscheil_constituent_batch,
        arg,
        n_cpu,
        progress=True,
        progress_callback=progress_callback,
        colour="#8060ff",
        ascii="░▒▓",
        bar_format="{l_bar}{bar:50}{r_bar}",
        desc="Nucleoscheil Constituent Batch",
    )


def nucleoscheil_batch(
    liquid_phase_name: str,
    database: Dict,
    condition: Dict[str, Union[float, str]],
    critical_undercooling: Optional[Dict[str, Union[float, str]]] = None,
    delta_T: float = 5.0,
    unit: Optional[List[str]] = None,
    n_cpu: int | None = None,
    n_per_batch: int = 1,
    phases: Optional[List[str]] = None,
    progress_callback=None,
) -> Dict:
    """Run nucleation-aware Scheil cooling for multiple conditions."""
    if unit is None:
        unit = ["K", "atm", "moles"]
    if n_cpu is None:
        n_cpu = os.cpu_count() or 1
    if critical_undercooling is None:
        critical_undercooling = _default_critical_undercooling(
            liquid_phase_name,
            database,
            condition,
            unit,
            phases,
        )
    # Get a batch input
    arg = _nucleoscheil_batch_input(
        liquid_phase_name,
        database,
        condition,
        critical_undercooling,
        delta_T,
        unit,
        n_per_batch,
    )
    res_mpi = _nucleoscheil_singlenode(arg, n_cpu, progress_callback)
    res = combine_dicts(res_mpi)

    phase_names = [liquid_phase_name] + list(critical_undercooling.keys())
    cAmount_n = [f"{phase}_amount_n [sp-mol]" for phase in phase_names]
    cAmount_w = [f"{phase}_amount_w [g]" for phase in phase_names]
    cEndmember_x = [f"{phase}_endmembers_x_" for phase in phase_names]
    cEndmember_w = [f"{phase}_endmembers_w_" for phase in phase_names]
    cElement_x = [f"{phase}_elements_x_" for phase in phase_names]
    cElement_w = [f"{phase}_elements_w_" for phase in phase_names]
    cEndmemberProperties = [
        f"{phase}{marker}"
        for phase in phase_names
        for marker in (
            "_partial_gibbs_",
            "_standard_gibbs_energy_",
            "_activity_",
            "_partial_enthalpy_",
            "_partial_entropy_",
            "_partial_heat_capacity_",
        )
    ]

    header = ["task_id"]
    header.append("T [K]")
    header.append("P [atm]")
    expanded_condition = expand_condition_species(database, condition, unit[2])
    elements = list(expanded_condition.keys())[2:]

    header.extend([f"{x} [sp-mol]" for x in elements])
    header.extend([f"{x} [g]" for x in elements])
    header.extend(
        [
            "label",
            "fl",
            "fs",
            "fl_w",
            "fs_w",
            "G [J]",
            "H [J]",
            "S [J/K]",
            "Cp [J/K]",
        ]
    )
    header.extend(cAmount_n)
    header.extend(cAmount_w)

    compositions = []
    compositions.extend(cEndmember_x)
    compositions.extend(cEndmember_w)
    compositions.extend(cElement_x)
    compositions.extend(cElement_w)
    compositions.extend(cEndmemberProperties)

    # Re-arrange the dictionary according to the header
    reordered_res = {key: res[key] for key in header if key in res}

    # Add the remaining items from the original dictionary
    for composition in compositions:
        for key, value in res.items():
            if key.startswith(composition):
                reordered_res[key] = value

    if "Error" in res:
        reordered_res["Error"] = res["Error"]

    return reordered_res


def nucleoscheil_constituent_batch(
    liquid_phase_name: str,
    database: Dict,
    condition: Dict[str, Union[float, str]],
    critical_undercooling: Optional[Dict[str, Union[float, str]]] = None,
    delta_T: float = 5.0,
    unit: Optional[List[str]] = None,
    n_cpu: int | None = None,
    n_per_batch: int = 1,
    phases: Optional[List[str]] = None,
    progress_callback=None,
) -> Dict:
    """Run nucleation-aware Scheil constituent calculations in batch."""
    if unit is None:
        unit = ["K", "atm", "moles"]
    if n_cpu is None:
        n_cpu = os.cpu_count() or 1
    if critical_undercooling is None:
        critical_undercooling = _default_critical_undercooling(
            liquid_phase_name,
            database,
            condition,
            unit,
            phases,
        )
    # Get a batch input
    arg = _nucleoscheil_batch_input(
        liquid_phase_name,
        database,
        condition,
        critical_undercooling,
        delta_T,
        unit,
        n_per_batch,
    )
    res_mpi = _nucleoscheil_constituent_singlenode(arg, n_cpu, progress_callback)
    res = combine_dicts(res_mpi)

    header = ["task_id"]
    expanded_condition = expand_condition_species(database, condition, unit[2])
    elements = list(expanded_condition.keys())[2:]
    header.append("T_Delta [K]")
    header.extend([f"{x} [sp-mol]" for x in elements])
    header.extend([f"{x} [g]" for x in elements])

    # Re-arrange the dictionary according to the header

    scheil_constituent_names = [
        key.split(" [sp-mol]")[0]
        for key in res.keys()
        if key not in header and "[sp-mol]" in key
    ]
    scheil_constituent_names = [key for key in scheil_constituent_names if key.strip()]

    primaries = [
        constituent
        for constituent in scheil_constituent_names
        if "+" not in constituent and constituent != ""
    ]
    eutectics = [
        constituent for constituent in scheil_constituent_names if "+" in constituent
    ]

    header.extend([f"{x} [sp-mol]" for x in primaries])
    header.extend([f"{x} [sp-mol]" for x in eutectics])
    header.extend([f"{x} [g]" for x in primaries])
    header.extend([f"{x} [g]" for x in eutectics])
    reordered_res = {key: res[key] for key in header if key in res}
    if "Error" in res:
        reordered_res["Error"] = res["Error"]

    # # Add the remaining items from the original dictionary
    # for key, value in res.items():
    #     if key not in reordered_res:
    #         reordered_res[key] = value

    return reordered_res
