#!/usr/bin/env python3
"""Scheil-Gulliver cooling and batch calculation routines."""

from __future__ import annotations

import os
from typing import Dict, List, Optional, Union

import numpy as np
from tqdm import tqdm

import equilipy.equilifort as fort

from ._parallel import starmap_joblib
from .composition import expand_condition_species
from .equilib_single import _equilib_single
from .exceptions import (
    EquilibError,
    TransitionError,
    liquidus_search_failure_warning,
)
from .find_transition import (
    _PhaseState,
    capture_phase_state,
    find_first_transition,
    find_liquidus_transition,
    find_transitions,
)
from .list_phases import list_phases
from .phase_selection import (
    default_scheil_phase_selection,
    scheil_ordered_phase_exclusion_notice,
)
from .results import capture_result_context
from .results.common import combine_dicts
from .results.equilib import EquilibPoint, get_assemblage_name
from .results.scheil import ScheilResult
from .utils import _dict2np
from .variables import cPeriodicTable

# from mpi4py.futures import MPIPoolExecutor
# from mpi4py import MPI

_TRUNCATED_LIQUID_FRACTION_WARNING_THRESHOLD = 1e-3


def _emit_cooling_warning(
    result: ScheilResult,
    message: str,
    progress_callback=None,
) -> None:
    """Record and surface a solidification warning."""
    warnings = getattr(result, "warnings", None)
    if warnings is None:
        result.warnings = []
        warnings = result.warnings
    warnings.append(message)
    print(message)
    if progress_callback is not None:
        progress_callback(0, 0, message)


def _warn_if_truncated_with_liquid(
    result: ScheilResult,
    *,
    liquid_phase_name: str,
    cooling_name: str,
    progress_callback=None,
) -> None:
    """Warn when a failed continuation leaves substantial liquid unprocessed."""
    points = result.points
    if not points:
        return
    final_point = points[-1]
    final_liquid_fraction = float(final_point.fl)
    if final_liquid_fraction <= _TRUNCATED_LIQUID_FRACTION_WARNING_THRESHOLD:
        return
    if liquid_phase_name not in final_point.label.split("+"):
        return
    _emit_cooling_warning(
        result,
        (
            f"Warning: {cooling_name} terminated at T={final_point.T:.6g} K "
            f"with liquid fraction {final_liquid_fraction:.6g} because "
            "equilibrium below this temperature failed to converge. "
            "Results are truncated."
        ),
        progress_callback=progress_callback,
    )


def _stable_assemblage_ids(assemblage) -> set[int]:
    """Return real stable phase ids from a Fortran assemblage array."""
    return {int(phase_id) for phase_id in np.asarray(assemblage).ravel() if phase_id}


def _current_stable_phase_names() -> set[str]:
    """Return stable phase names from the current Fortran assemblage."""
    assemblage_ids = _stable_assemblage_ids(fort.modulethermo.iassemblage)
    return {
        str(name).strip()
        for name in get_assemblage_name(assemblage_ids)
        if str(name).strip()
    }


def _is_pure_liquid_state(liquid_phase_name: str) -> bool:
    """Return True when the current equilibrium state is only liquid."""
    return _current_stable_phase_names() == {liquid_phase_name}


def _latest_stable_phase_names(result: ScheilResult) -> list[str]:
    """Return the most recent equilibrium stable phase names."""
    try:
        names = result.equilib_result.points[-1].stable_phases.names
    except (AttributeError, IndexError):
        return []
    return [str(name).strip() for name in names if str(name).strip()]


def _latest_phase_state(result: ScheilResult) -> _PhaseState:
    """Return transition phase state from the latest stored equilibrium point."""
    try:
        phases = result.equilib_result.points[-1].stable_phases
    except (AttributeError, IndexError):
        return _PhaseState(set())
    return _PhaseState(
        {int(phase.id) for phase in phases.values() if int(phase.id) != 0}
    )


def _capture_or_latest_phase_state(result: ScheilResult) -> _PhaseState:
    """Return live Fortran phase state, falling back to stored result state."""
    try:
        return capture_phase_state()
    except AttributeError:
        return _latest_phase_state(result)


def _previous_stable_phase_retry(
    result: ScheilResult,
    liquid_phase_name: str,
    phases: Optional[List[str]],
) -> list[str] | None:
    """Return a restricted retry phase set from the latest Scheil point."""
    requested = set(phases) if phases is not None else None
    retry_phases: list[str] = []
    for phase_name in [liquid_phase_name, *_latest_stable_phase_names(result)]:
        if not phase_name:
            continue
        if requested is not None and phase_name not in requested:
            continue
        if phase_name not in retry_phases:
            retry_phases.append(phase_name)

    if len(retry_phases) <= 1:
        return None
    return retry_phases


def _transition_search_floor(unit: List[str]) -> float:
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


def _transition_temperature_pairs(
    temperatures: np.ndarray,
    tolerance: float = 1e-8,
) -> list[tuple[float, float]]:
    """Return validated high/low transition-temperature pairs."""
    values = np.asarray(temperatures, dtype=float).ravel()
    pairs: list[tuple[float, float]] = []
    for index in range(0, len(values) - 1, 2):
        high_temperature = max(float(values[index]), float(values[index + 1]))
        low_temperature = min(float(values[index]), float(values[index + 1]))
        if any(
            abs(high_temperature - high) <= tolerance
            and abs(low_temperature - low) <= tolerance
            for high, low in pairs
        ):
            continue
        pairs.append((high_temperature, low_temperature))
    return pairs


def _reset_equilib_points(result: ScheilResult, points: list[EquilibPoint]) -> None:
    """Replace the equilibrium point history while preserving scalar storage."""
    result.equilib_result.data = points[0] if len(points) == 1 else list(points)


def _rebuild_scheil_points(result: ScheilResult) -> None:
    """Rebuild Scheil cumulative rows from the equilibrium point history."""
    rows = []
    previous_row = None
    for equilib_point in result.equilib_result.points:
        row = result._scheil_point_from_equilib_step(
            equilib_point,
            previous_point=previous_row,
        )
        rows.append(row)
        previous_row = row
    result.data = rows[0] if len(rows) == 1 else rows


def _insert_transition_hot_side_point(
    result: ScheilResult,
    hot_side_point: EquilibPoint,
) -> None:
    """Insert a hot-side transition point before the latest cold-side row."""
    points = result.equilib_result.points
    if not points:
        return
    points.insert(max(len(points) - 1, 0), hot_side_point)
    _reset_equilib_points(result, points)
    _rebuild_scheil_points(result)


def _update_condition_from_latest_liquid(
    condition: Dict[str, Union[float, str]],
    result: ScheilResult,
    liquid_phase_name: str,
) -> bool:
    """Update condition amounts from the latest liquid composition."""
    latest_point = result.equilib_result.points[-1]
    if liquid_phase_name not in latest_point.stable_phases.names:
        return False

    liquid_amount = result.phase_amount_n(liquid_phase_name)[-1]
    liquid_phase_obj = latest_point.phase(liquid_phase_name)
    for component, fraction in liquid_phase_obj.elements.x_i.items():
        condition[component] = float(fraction) * liquid_amount
    return True


def scheil_cooling(
    liquid_phase_name: str,
    database: Dict,
    condition: Dict[str, Union[float, str]],
    delta_T: float = 5.0,
    unit: Optional[List[str]] = None,
    phases: Optional[List[str]] = None,
    include_ordered: bool = False,
    start_from_liquidus=True,
    progress_callback=None,
    show_progress: bool = True,
) -> ScheilResult:
    """
    Perform a Scheil-Gulliver cooling simulation.

    This function simulates the solidification of a liquid phase under the assumptions
    of the Scheil-Gulliver model:
    1. No diffusion in solid phases.
    2. Infinitely fast diffusion in the liquid phase.
    3. Equilibrium is maintained at the solid-liquid interface.

    The simulation proceeds by cooling the system in discrete temperature steps,
    calculating the equilibrium at each step, removing the newly formed solid phases,
    and then continuing the cooling with the remaining liquid.

    Args:
        liquid_phase_name (str): The name of the liquid phase in the database.
        database (Dict): The thermodynamic database dictionary.
        condition (Dict[str, Union[float, str]]): A dictionary defining the initial
            conditions (temperature, pressure, and composition). This dictionary
            is not modified by the function.
        unit (List[str], optional): A list specifying the units for temperature,
            pressure, and amount. Defaults to ['K', 'atm', 'moles'].
        delta_T (float, optional): The temperature step for each cooling iteration.
            Defaults to 5.0.
        IterMax (int, optional): The maximum number of iterations to prevent
            infinite loops. Defaults to 5000.
        phases (Optional[List[str]], optional): A list of phase names to
            consider in the equilibrium calculation. If None, Scheil uses the
            database phase list but excludes ordered SUBOM phases by default.
            Defaults to None.
        include_ordered (bool, optional): Include ordered SUBOM phases in the
            default phase list when ``phases`` is None. Explicit ``phases`` are
            never filtered. Defaults to False.
        show_progress (bool, optional): Whether to show the terminal tqdm progress
            bar when no GUI progress callback is provided. Defaults to True.

    Returns
    -------
        ScheilResult: An object containing the detailed results of the Scheil
            simulation, including temperature, phase amounts, and compositions
            at each step.

    Raises
    ------
        EquilibError: If the initial condition does not contain a stable liquid phase
            or if the temperature drops below 0 K during the simulation.
    """
    if unit is None:
        unit = ["K", "atm", "moles"]
    if delta_T <= 0:
        raise EquilibError("delta_T must be positive for Scheil cooling.")
    normalized_condition = expand_condition_species(database, condition, unit[2])
    elements = list(normalized_condition.keys())[2:]
    default_excluded_ordered_phases: list[str] = []
    if phases is None:
        phases = list(list_phases(database, elements))
        if not include_ordered:
            phases, default_excluded_ordered_phases = default_scheil_phase_selection(
                database,
                phases,
            )
    else:
        phases = list(phases)
    result_context = capture_result_context(
        database,
        normalized_condition,
        unit,
        phases,
    )

    # Create a copy of the initial conditions to avoid modifying the user's dictionary
    current_condition = normalized_condition.copy()

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
        # atomic_masses = np.array([cPeriodicTable[el.strip()][1] for el in elements])
        for el in elements:
            current_condition[el] = (
                current_condition[el] / cPeriodicTable[el.strip()][1]
            )
    unit_local = unit.copy()
    unit_local[2] = "moles"
    # Start from liquidus
    liquidus_warning = None
    if start_from_liquidus:
        if phases is not None and not any(
            phase != liquid_phase_name for phase in phases
        ):
            raise TransitionError(
                "No phase transition found because only the liquid phase was selected."
            )
        Tmax = current_condition["T"]
        Tmin = max(current_condition["T"] * 0.1, _transition_search_floor(unit_local))
        try:
            Ts = find_liquidus_transition(
                database,
                current_condition,
                liquid_phase_name,
                Tmax,
                Tmin,
                unit=unit_local,
                phases=phases,
            )
        except TransitionError as exc:
            liquidus_warning = liquidus_search_failure_warning(exc)
        else:
            current_condition["T"] = Ts + 0.1

    # --- Initial Equilibrium Calculation ---
    # Perform a single equilibrium calculation to establish the starting point
    # Scheil evolves phase amounts and the remaining-liquid composition; it does
    # not consume an equilibrium heat-capacity perturbation.  At a solved phase
    # boundary, the fixed-assemblage T-dT auxiliary state may not exist.  Keep
    # every path solve on the requested temperature and use the fixed-state Cp
    # assembled there.
    _equilib_single(
        database,
        current_condition,
        unit=unit_local,
        phases=phases,
        include_heat_capacity=False,
    )
    state_old = capture_phase_state()
    assemblage_old = set(state_old.ids)

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
    if default_excluded_ordered_phases:
        _emit_cooling_warning(
            res,
            scheil_ordered_phase_exclusion_notice(default_excluded_ordered_phases),
            progress_callback=progress_callback,
        )
    stable_phase_names = list(res.phase_amounts_n.keys())

    # Check if the liquid phase is stable at the start
    if liquid_phase_name not in stable_phase_names:
        raise EquilibError(
            f"The specified liquid phase '{liquid_phase_name}' is not stable "
            "at the initial conditions."
        )
    initial_liquid_amount = res.phase_amount_n(liquid_phase_name)[-1]
    if initial_liquid_amount <= 0:
        raise EquilibError(
            f"The specified liquid phase '{liquid_phase_name}' has zero amount "
            "at the initial conditions."
        )
    solidification_started = False
    truncated_due_to_failed_continuation = False

    IterMax = int(current_condition["T"] // delta_T)
    if progress_callback is not None:
        progress_callback(0, IterMax, "Scheil Cooling")

    # --- Main Cooling Loop ---
    # The for loop iterates up to IterMax times, with a progress bar.
    for step_index in tqdm(
        range(IterMax),
        desc="Scheil Cooling",
        colour="#8060ff",
        ascii="░▒▓",
        bar_format="{l_bar}{bar:50}{r_bar}",
        disable=not show_progress or progress_callback is not None,
    ):
        if progress_callback is not None:
            progress_callback(step_index + 1, IterMax, "Scheil Cooling")
        # Get liquid composition and amount from the previous step.
        liquid_amount = res.phase_amount_n(liquid_phase_name)[-1]
        if step_index == 0:
            current_condition["T"] = res.T  # This converts temperature unit to Kelvin
            unit_local[0] = "K"

            liquid_phase_obj = res.equilib_result.phase(liquid_phase_name)

            if liquid_phase_obj.elements.x_i:
                components_tuple, values_tuple = zip(
                    *liquid_phase_obj.elements.x_i.items(),
                    strict=False,
                )
                liquid_components = list(components_tuple)
                liquid_composition = list(values_tuple)
            else:
                liquid_components = []
                liquid_composition = []

        else:
            liquid_phase_obj = res.equilib_result.points[-1].phase(liquid_phase_name)

            if liquid_phase_obj.elements.x_i:
                components_tuple, values_tuple = zip(
                    *liquid_phase_obj.elements.x_i.items(),
                    strict=False,
                )
                liquid_components = list(components_tuple)
                liquid_composition = list(values_tuple)
            else:
                liquid_components = []
                liquid_composition = []

        # --- Termination Conditions ---
        if liquid_amount / initial_liquid_amount < 1e-5:
            # Solidification is considered complete
            break

        current_temp = float(current_condition["T"])
        if current_temp <= 0 and "K" in unit_local[0]:
            raise EquilibError(
                "Scheil-Gulliver solidification results in melting point "
                "below 0 K. Check phase selections."
            )

        # --- Update Conditions for Next Step ---
        # Decrease temperature
        current_condition["T"] = current_temp - delta_T

        # Update composition based on remaining liquid
        new_ni = np.array(liquid_composition) * liquid_amount
        for component_index, el in enumerate(liquid_components):
            current_condition[el] = new_ni[component_index]

        # --- Perform Equilibrium Calculation for the Current Step ---
        step_equilibrium_failed = False
        used_stable_set_retry = False
        try:
            _equilib_single(
                database,
                current_condition,
                unit=unit_local,
                phases=phases,
                include_heat_capacity=False,
            )
            if solidification_started and _is_pure_liquid_state(liquid_phase_name):
                retry_phases = _previous_stable_phase_retry(
                    res,
                    liquid_phase_name,
                    phases,
                )
                if retry_phases is None:
                    step_equilibrium_failed = True
                else:
                    fort.resetthermo()
                    try:
                        _equilib_single(
                            database,
                            current_condition,
                            unit=unit_local,
                            phases=retry_phases,
                            include_heat_capacity=False,
                        )
                        used_stable_set_retry = True
                    except Exception:
                        step_equilibrium_failed = True
        except Exception:
            retry_phases = _previous_stable_phase_retry(
                res,
                liquid_phase_name,
                phases,
            )
            if retry_phases is None:
                step_equilibrium_failed = True
            else:
                try:
                    _equilib_single(
                        database,
                        current_condition,
                        unit=unit_local,
                        phases=retry_phases,
                        include_heat_capacity=False,
                    )
                    used_stable_set_retry = True
                except Exception:
                    step_equilibrium_failed = True

        # Before updating the result, calculate any transition point.
        if step_equilibrium_failed:
            state_new = None
        elif used_stable_set_retry:
            state_new = state_old
        else:
            state_new = capture_phase_state()
        assemblage_new = set() if state_new is None else set(state_new.ids)
        reached_terminal_liquidus = False
        step_output_appended = False
        if step_equilibrium_failed or assemblage_new != assemblage_old:
            # # Transition occured, Reset variables
            Tmax = float(current_temp)
            Tmin = float(current_condition["T"])
            current_condition["T"] = Tmax
            transition_search_failed = False
            try:
                Ts = find_transitions(
                    database,
                    current_condition,
                    Tmax,
                    Tmin,
                    unit=unit_local,
                    phases=phases,
                    state_at_T_max=state_old,
                    state_at_T_min=state_new,
                )
            except TransitionError:
                transition_search_failed = True
                Ts = np.asarray([])
                current_condition["T"] = Tmin
                try:
                    _equilib_single(
                        database,
                        current_condition,
                        unit=unit_local,
                        phases=phases,
                        include_heat_capacity=False,
                    )
                except Exception:
                    current_condition["T"] = Tmax
                    reached_terminal_liquidus = True
                    truncated_due_to_failed_continuation = True
            if reached_terminal_liquidus:
                break
            if not transition_search_failed:
                if len(Ts) == 0:
                    raise EquilibError(
                        "A phase change was detected during Scheil cooling, but no "
                        "transition temperature could be resolved. Check the selected "
                        "phases and cooling step size."
                    )
                transition_pairs = _transition_temperature_pairs(Ts)
                last_liquid_temperature = float(Tmax)

                for high_transition, low_transition in transition_pairs:
                    # Preserve the original Scheil evolution by solving and
                    # updating the cold side first.  The hot-side bracket point
                    # is captured afterward and inserted before the cold-side
                    # row for reporting.
                    hot_side_condition = current_condition.copy()
                    hot_side_condition["T"] = high_transition

                    current_condition["T"] = low_transition
                    _equilib_single(
                        database,
                        current_condition,
                        unit=unit_local,
                        phases=phases,
                        include_heat_capacity=False,
                    )
                    if solidification_started and _is_pure_liquid_state(
                        liquid_phase_name
                    ):
                        continue
                    res.append_output()
                    step_output_appended = True
                    solidification_started = True
                    if not _update_condition_from_latest_liquid(
                        current_condition,
                        res,
                        liquid_phase_name,
                    ):
                        reached_terminal_liquidus = True
                    last_liquid_temperature = float(low_transition)

                    try:
                        _equilib_single(
                            database,
                            hot_side_condition,
                            unit=unit_local,
                            phases=phases,
                            include_heat_capacity=False,
                        )
                        if not _is_pure_liquid_state(liquid_phase_name):
                            _insert_transition_hot_side_point(
                                res,
                                EquilibPoint.from_fortran(),
                            )
                    except Exception:
                        pass
                    finally:
                        fort.resetthermo()

                    if reached_terminal_liquidus:
                        break

                if reached_terminal_liquidus:
                    break

                current_condition["T"] = float(Tmin)
                terminal_search_needed = False
                try:
                    _equilib_single(
                        database,
                        current_condition,
                        unit=unit_local,
                        phases=phases,
                        include_heat_capacity=False,
                    )
                except Exception:
                    terminal_search_needed = True
                else:
                    terminal_search_needed = (
                        liquid_phase_name not in _current_stable_phase_names()
                    )

                if terminal_search_needed:
                    current_condition["T"] = last_liquid_temperature
                    try:
                        terminal_temperature = find_first_transition(
                            database,
                            current_condition,
                            last_liquid_temperature,
                            Tmin,
                            unit=unit_local,
                            phases=phases,
                        )
                        current_condition["T"] = terminal_temperature
                        _equilib_single(
                            database,
                            current_condition,
                            unit=unit_local,
                            phases=phases,
                            include_heat_capacity=False,
                        )
                        reached_terminal_liquidus = True
                    except TransitionError:
                        current_condition["T"] = last_liquid_temperature
                        reached_terminal_liquidus = True
                        truncated_due_to_failed_continuation = True
        if not used_stable_set_retry:
            state_old = _capture_or_latest_phase_state(res)
        assemblage_old = set(state_old.ids)
        if step_output_appended:
            if liquid_phase_name not in _latest_stable_phase_names(res):
                break
            continue
        if solidification_started and _is_pure_liquid_state(liquid_phase_name):
            continue
        res.append_output()
        if any(
            phase_name != liquid_phase_name
            for phase_name in res.equilib_result.points[-1].stable_phases.names
        ):
            solidification_started = True
        if reached_terminal_liquidus:
            break

        # Check if liquid phase disappeared in the last step
        if liquid_phase_name not in res.equilib_result.points[-1].stable_phases.names:
            break
    else:
        # This block runs only if the for loop completes without a 'break'
        print(
            f"Warning: Reached maximum iterations ({IterMax}) before "
            "solidification completed."
        )

    if truncated_due_to_failed_continuation:
        _warn_if_truncated_with_liquid(
            res,
            liquid_phase_name=liquid_phase_name,
            cooling_name="Scheil cooling",
            progress_callback=progress_callback,
        )

    # Finalize results
    res.update_scheil_constituents()
    # Reorder
    header = ["task_id"]
    elements = list(current_condition.keys())[2:]
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


def _scheil_batch_input(
    liquid_phase_name: str,
    database: Dict,
    condition: Dict[str, Union[float, str]],
    delta_T: float = 5.0,
    unit: Optional[List[str]] = None,
    phases: Optional[List[str]] = None,
    include_ordered: bool = False,
    start_from_liquidus=True,
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
    n, m = NPT.shape
    fullrange = [0, n]

    # Scheil workers run one solidification path per task; keep n_per_batch for
    # API compatibility but pass scalar row conditions to scheil_cooling.
    for i in range(fullrange[0], fullrange[1]):
        subcondition = dict({})
        for j, head in enumerate(header):
            subcondition[head] = NPT[i, j].item()
        res.append(
            [
                liquid_phase_name,
                database,
                subcondition,
                delta_T,
                unit,
                phases,
                include_ordered,
                start_from_liquidus,
            ]
        )
    return res


def _scheil_batch(
    liquid_phase_name: str,
    database: Dict,
    condition: Dict[str, Union[float, str]],
    delta_T: float = 5.0,
    unit: Optional[List[str]] = None,
    phases: Optional[List[str]] = None,
    include_ordered: bool = False,
    start_from_liquidus=True,
) -> Dict:
    if unit is None:
        unit = ["K", "atm", "moles"]

    try:
        res = scheil_cooling(
            liquid_phase_name,
            database,
            condition,
            delta_T=delta_T,
            unit=unit,
            phases=phases,
            include_ordered=include_ordered,
            start_from_liquidus=start_from_liquidus,
            show_progress=False,
        )
        return res.to_dict()
    except Exception as exc:
        return {"Error": f"Scheil cooling failed: {type(exc).__name__}: {exc}"}


def _scheil_constituent_batch(
    liquid_phase_name: str,
    database: Dict,
    condition: Dict[str, Union[float, str]],
    delta_T: float = 5.0,
    unit: Optional[List[str]] = None,
    phases: Optional[List[str]] = None,
    include_ordered: bool = False,
    start_from_liquidus=True,
) -> Dict:
    if unit is None:
        unit = ["K", "atm", "moles"]

    try:
        res = scheil_cooling(
            liquid_phase_name,
            database,
            condition,
            delta_T=delta_T,
            unit=unit,
            phases=phases,
            include_ordered=include_ordered,
            start_from_liquidus=start_from_liquidus,
            show_progress=False,
        )
        return res.scheil_constituents.copy()
    except Exception as exc:
        return {
            "Error": (
                f"Scheil constituent calculation failed: {type(exc).__name__}: {exc}"
            )
        }


def _scheil_singlenode(arg, n_cpu, progress_callback=None) -> List[Dict]:
    return starmap_joblib(
        _scheil_batch,
        arg,
        n_cpu,
        progress=True,
        progress_callback=progress_callback,
        colour="#8060ff",
        ascii="░▒▓",
        bar_format="{l_bar}{bar:50}{r_bar}",
        desc="Scheil Batch",
    )


def _scheil_constituent_singlenode(arg, n_cpu, progress_callback=None) -> List[Dict]:
    return starmap_joblib(
        _scheil_constituent_batch,
        arg,
        n_cpu,
        progress=True,
        progress_callback=progress_callback,
        colour="#8060ff",
        ascii="░▒▓",
        bar_format="{l_bar}{bar:50}{r_bar}",
        desc="Scheil Constituent Batch",
    )


# def _scheil_multinodes(arg,n_cpu)-> List[Dict]:
#     print(f'Process with mpi: {n_cpu} processors')
#     with MPIPoolExecutor(max_workers=n_cpu) as pool:
#         res = pool.starmap(_scheil_batch, arg)
#     return res


def scheil_batch(
    liquid_phase_name: str,
    database: Dict,
    condition: Dict[str, Union[float, str]],
    delta_T: float = 5.0,
    unit: Optional[List[str]] = None,
    phases: Optional[List[str]] = None,
    include_ordered: bool = False,
    start_from_liquidus=True,
    n_cpu: int | None = None,
    n_per_batch: int = 1,
    progress_callback=None,
) -> Dict:
    """Run Scheil-Gulliver cooling for multiple conditions."""
    if unit is None:
        unit = ["K", "atm", "moles"]
    if n_cpu is None:
        n_cpu = os.cpu_count() or 1

    # Get a batch input
    arg = _scheil_batch_input(
        liquid_phase_name,
        database,
        condition,
        delta_T,
        unit,
        phases,
        include_ordered,
        start_from_liquidus,
        n_per_batch,
    )
    res_mpi = _scheil_singlenode(arg, n_cpu, progress_callback)
    res = combine_dicts(res_mpi)

    res_keys = list(res.keys())
    phase_names = [
        x.split("_amount_n [sp-mol]")[0]
        for x in res_keys
        if "_amount_n [sp-mol]" in x
    ]
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
    expanded_condition = expand_condition_species(database, condition, unit[2])
    elements = list(expanded_condition.keys())[2:]
    header.append("T [K]")
    header.append("P [atm]")
    header.extend([f"n_{x} [sp-mol]" for x in elements])
    header.extend([f"w_{x} [g]" for x in elements])
    header.extend(
        [
            "label",
            "fl",
            "fs",
            "fl_w",
            "fs_w",
            "G [J]",
            "H [J]",
            "Q [J]",
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


def scheil_constituent_batch(
    liquid_phase_name: str,
    database: Dict,
    condition: Dict[str, Union[float, str]],
    delta_T: float = 5.0,
    unit: Optional[List[str]] = None,
    phases: Optional[List[str]] = None,
    include_ordered: bool = False,
    start_from_liquidus=True,
    n_cpu: int | None = None,
    n_per_batch: int = 1,
    progress_callback=None,
) -> Dict:
    """Run Scheil constituent calculations for multiple conditions."""
    if unit is None:
        unit = ["K", "atm", "moles"]
    if n_cpu is None:
        n_cpu = os.cpu_count() or 1

    # Get a batch input
    arg = _scheil_batch_input(
        liquid_phase_name,
        database,
        condition,
        delta_T,
        unit,
        phases,
        include_ordered,
        start_from_liquidus,
        n_per_batch,
    )
    res_mpi = _scheil_constituent_singlenode(arg, n_cpu, progress_callback)
    res = combine_dicts(res_mpi)
    header = ["task_id"]
    expanded_condition = expand_condition_species(database, condition, unit[2])
    elements = list(expanded_condition.keys())[2:]
    header.append("T_Delta [K]")
    header.extend([f"{x} [sp-mol]" for x in elements])
    header.extend([f"{x} [g]" for x in elements])

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
