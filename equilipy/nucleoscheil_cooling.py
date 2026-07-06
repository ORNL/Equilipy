"""Nucleoscheil solidification routines."""

from __future__ import annotations

import copy
import os
from typing import Dict, List, Optional, Union

import numpy as np
from tqdm import tqdm

from ._parallel import starmap_joblib
from .composition import expand_condition_species
from .equilib_single import _equilib_single, equilib_single
from .exceptions import EquilibError, InputConditionError, TransitionError
from .find_transition import find_first_transition
from .list_phases import list_phases
from .results import capture_result_context
from .results.common import combine_dicts
from .results.scheil import ScheilResult
from .scheil_cooling import _warn_if_truncated_with_liquid
from .utils import _dict2np
from .variables import cPeriodicTable


# ============================================================================
# Nucleoscheil functions
def _default_critical_undercooling(
    liquid_phase_name: str,
    database: Dict,
    condition: Dict[str, Union[float, str]],
    unit: Optional[List[str]] = None,
    phases: Optional[List[str]] = None,
) -> Dict[str, float]:
    """Return zero nucleation undercooling for all eligible solid phases."""
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
    return {
        str(phase): 0.0
        for phase in phase_names
        if (
            str(phase).strip()
            and str(phase).strip() != liquid_phase_name
            and not _is_reference_ser_phase(str(phase).strip())
        )
    }


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
) -> (float, str):
    if unit is None:
        unit = ["K", "atm", "moles"]

    res = []
    solid_phases = list(critical_undercooling.keys())

    for phase in solid_phases:
        # 1 Get liquidus temperature for each solid phase
        try:
            T_liquidus = find_first_transition(
                database,
                condition,
                condition["T"],
                max(float(condition["T"]) * 0.1, _transition_search_floor(unit)),
                unit=unit,
                phases=[liquid_phase_name, phase],
            )
        except TransitionError:
            continue
        # Append a tuple of (float, str)
        res.append((T_liquidus - critical_undercooling[phase], phase))

    if not res:
        raise EquilibError(
            "No nucleating phase could be initialized. Check phase selection "
            "and nucleation undercooling values."
        )

    # Sort tuples in descending order based on the first element
    # (the temperature).
    res.sort(key=lambda x: x[0], reverse=True)

    # The highest temperature nucleation event is now the first item
    return res[0]


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
        res[element] = np.round(float(fraction) * liquid_amount, 10)
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
        if equilib_result.phase(phase).amount_n > 1e-12:
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
            solid_phases: List[str] = []
            for point in points[start_index:end_index]:
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
) -> Dict[str, float]:
    if unit is None:
        unit = ["K", "atm", "moles"]

    res = {}
    NPT_local = condition.copy()
    T_current = condition["T"]
    solid_phases = list(critical_undercooling.keys())
    # 3.2 Identify which phase candidates nucleate next
    for phase in solid_phases:
        NPT_local["T"] = T_current + critical_undercooling[phase]
        try:
            res_temp = equilib_single(
                database,
                NPT_local,
                unit=unit,
                phases=[liquid_phase_name, phase],
                include_heat_capacity=False,
            )
        except (EquilibError, InputConditionError, TransitionError, ValueError):
            continue
        if phase in res_temp.stable_phases.names:
            res[phase] = res_temp.phase(phase).amount_n
    return res


def _get_same_temperature_nuclei_candidates(
    liquid_phase_name: str,
    database: Dict,
    condition: Dict[str, Union[float, str]],
    critical_undercooling: Dict[str, float],
    unit: Optional[List[str]] = None,
    excluded_phases: Optional[List[str]] = None,
) -> Dict[str, float]:
    if unit is None:
        unit = ["K", "atm", "moles"]
    excluded = set(excluded_phases or [])

    res = {}
    for phase in critical_undercooling:
        if phase in excluded:
            continue
        try:
            res_temp = equilib_single(
                database,
                condition,
                unit=unit,
                phases=[liquid_phase_name, phase],
                include_heat_capacity=False,
            )
        except (EquilibError, InputConditionError, TransitionError, ValueError):
            continue
        if phase not in res_temp.stable_phases.names:
            continue
        amount = res_temp.phase(phase).amount_n
        if amount > 1e-12:
            res[phase] = amount
    return res


def _check_eutectic(
    liquid_phase_name: str,
    latest_scheil_result: ScheilResult,
    nuclei_current: List[str],
    database: Dict,
    condition: Dict[str, Union[float, str]],
    critical_undercooling: Dict[str, float],
    unit: Optional[List[str]] = None,
) -> tuple[ScheilResult, bool]:
    if unit is None:
        unit = ["K", "atm", "moles"]

    # First make a copy of the result
    res = copy.deepcopy(latest_scheil_result)
    lfinal = False

    # Start from the liquid state before the current possible invariant event.
    NPT_start = _update_liquid_composition(liquid_phase_name, res, condition)
    nuclei_final = []
    candidates = list(nuclei_current)
    elements = list(condition.keys())[2:]
    NPT_liq = NPT_start.copy()
    if not candidates:
        return res, lfinal

    # At a fixed temperature, a solid precipitate changes the liquid
    # composition. Recheck the new liquid at the same temperature before
    # falling back to undercooling-controlled nucleation. This catches
    # peritectic/invariant cascades without changing the normal cooling loop.
    for _ in range(max(len(critical_undercooling), len(elements)) + 1):
        phase_selection = [liquid_phase_name] + candidates
        eqres = equilib_single(
            database,
            NPT_liq,
            unit=unit,
            phases=phase_selection,
            include_heat_capacity=False,
        )
        stable_solids = _stable_solid_candidates(
            eqres,
            candidates,
            liquid_phase_name,
        )
        if not stable_solids:
            return res, lfinal

        _append_unique_phases(nuclei_final, stable_solids)
        liquid_phase = eqres.phase(liquid_phase_name)
        if float(liquid_phase.amount_n) < 1e-5:
            lfinal = True
            break
        if len(nuclei_final) >= len(elements):
            lfinal = True
            break

        NPT_liq = _update_liquid_composition_from_equilib(
            liquid_phase_name,
            eqres,
            NPT_liq,
        )
        # The exact-temperature cascade is needed for multicomponent invariant
        # paths. Binary paths keep the original undercooling-controlled
        # selection so high-undercooling competing phases cannot bypass their
        # nucleation barrier.
        nuclei_more_candidates = {}
        if len(elements) > 2:
            nuclei_more_candidates = _get_same_temperature_nuclei_candidates(
                liquid_phase_name,
                database,
                NPT_liq,
                critical_undercooling,
                unit=unit,
                excluded_phases=nuclei_final,
            )
        if not nuclei_more_candidates:
            nuclei_more_candidates = _get_nuclei_candiates(
                liquid_phase_name,
                database,
                NPT_liq,
                critical_undercooling,
                unit=unit,
            )
            nuclei_more_candidates = {
                phase: amount
                for phase, amount in nuclei_more_candidates.items()
                if phase not in nuclei_final
            }
        candidates = [
            phase
            for phase in nuclei_more_candidates
            if phase not in nuclei_final
        ]
        if not candidates:
            return res, lfinal

    if lfinal:
        _equilib_single(
            database,
            NPT_start,
            unit=unit,
            phases=[liquid_phase_name] + nuclei_final,
        )
        res.append_output()
    return res, lfinal


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
    try:
        T_init = find_first_transition(
            database,
            NPT_local,
            T_current,
            max(T_current * 0.1, _transition_search_floor(unit_local)),
            unit=unit_local,
            phases=[liquid_phase_name] + solid_phases,
        )
    except TransitionError:
        T_init = T_current

    T_current = T_init + 0.1
    NPT_local["T"] = T_current
    _equilib_single(database, NPT_local, unit=unit_local, phases=[liquid_phase_name])
    res = ScheilResult(
        liquid_phase_name=liquid_phase_name,
        input_unit=unit,
        context=result_context,
    )
    # now all units are converted to K
    unit_local[0] = "K"
    T_current = res.T
    NPT_local["T"] = T_current

    # Step2: Get the first nucleating phase and temperature
    _equilib_single(database, NPT_local, unit=unit_local, phases=[liquid_phase_name])
    res.append_output()

    T_current, nuclei_init = _init_nuclei(
        liquid_phase_name, database, NPT_local, critical_undercooling, unit=unit_local
    )
    NPT_local["T"] = T_current
    res.data[-1].T = float(T_current)

    _equilib_single(
        database,
        NPT_local,
        unit=unit_local,
        phases=[liquid_phase_name, nuclei_init],
    )
    res.append_output()
    nuclei_prev = [nuclei_init]
    phase_selection = [liquid_phase_name] + nuclei_prev
    truncated_due_to_failed_continuation = False

    # Step2 Loop: Proceed to next temperature step
    iteration_count = int(T_init / delta_T)
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
        # Loop Step 1: Update temperature
        T_current = T_current - delta_T
        if T_current < 295.0:
            break
        NPT_local["T"] = T_current

        # Loop Step 2: Get nuclei
        NPT_liq = _update_liquid_composition(liquid_phase_name, res, NPT_local)
        nuclei_candidates = _get_nuclei_candiates(
            liquid_phase_name, database, NPT_liq, critical_undercooling, unit=unit_local
        )
        nuclei_candidates = list(nuclei_candidates.keys())
        nuclei_number = len(nuclei_candidates)
        if (
            nuclei_number == 0
            or len(nuclei_prev) == 0
            or set(nuclei_candidates) == set(nuclei_prev)
        ):
            # condition 1: if nuclei_number is zero or nuclei_candidates are
            # same as previous, continue with previous nuclei.
            if len(nuclei_prev) != 0:
                phase_selection = [liquid_phase_name] + nuclei_prev
            nuclei_current = nuclei_prev.copy()
        else:
            # condition 2: if new nuclei are identified, conduct metastable
            # equilib calculation with new nuclei added.
            phase_selection = [liquid_phase_name] + nuclei_candidates
            nuclei_current = nuclei_candidates.copy()

        # Check if it is final temperature
        try:
            res, lfinal = _check_eutectic(
                liquid_phase_name,
                res,
                nuclei_current,
                database,
                NPT_liq,
                critical_undercooling,
                unit=unit_local,
            )
        except Exception:
            truncated_due_to_failed_continuation = True
            break
        if lfinal and liquid_phase_name not in res.data[-1].label.split("+"):
            break
        else:
            try:
                _equilib_single(
                    database,
                    NPT_liq,
                    unit=unit_local,
                    phases=phase_selection,
                )
            except Exception:
                truncated_due_to_failed_continuation = True
                break
            res.append_output()
            nuclei_prev = [x for x in phase_selection if x != liquid_phase_name]
            NPT_local = NPT_liq
            T_current = NPT_liq["T"]

        # # if not move on to the next
        # _equilib_single(
        #     database, NPT_liq, unit=unit_local, phases=phase_selection
        # )
        # res.append_output()
        # nuclei_prev = [ x for x in phase_selection if x !=liquid_phase_name]

        # # Check if there is any other phases nucleating at this temperature
        # if liquid_phase_name in res.data[-1].label.split('+'):
        #         NPT_liq = _update_liquid_composition(liquid_phase_name, res, NPT_liq)
        #         res, lfinal = _check_eutectic(
        #             liquid_phase_name,
        #             res,
        #             nuclei_current,
        #             database,
        #             NPT_liq,
        #             critical_undercooling,
        #             unit=unit_local,
        #         )
        # else:
        #     lfinal = True

        # if lfinal: break
        # else: NPT_local = NPT_liq
    if truncated_due_to_failed_continuation:
        _warn_if_truncated_with_liquid(
            res,
            liquid_phase_name=liquid_phase_name,
            cooling_name="Nucleoscheil cooling",
            progress_callback=progress_callback,
        )
    _collapse_same_temperature_output(res, liquid_phase_name)
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
            np.full(uc_target_length, column[0])
            if len(column) == 1
            else column
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
