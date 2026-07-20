"""Equilibrium cooling path calculation."""

from __future__ import annotations

from typing import Any

import numpy as np
from tqdm import tqdm

import equilipy.equilifort as fort

from .equilib_single import _equilib_single
from .exceptions import (
    EquilibError,
    TransitionError,
    liquidus_search_failure_warning,
)
from .find_transition import find_liquidus_transition, find_transitions
from .results import capture_result_context
from .results.equilib import EquilibPoint, EquilibResult
from .results.scheil import ScheilPoint, ScheilResult


def _stable_assemblage_ids(assemblage: Any) -> set[int]:
    """Return real stable phase ids from a Fortran assemblage array."""
    return {int(phase_id) for phase_id in np.asarray(assemblage).ravel() if phase_id}


def _transition_search_floor(unit: list[str]) -> float:
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


def _append_current_equilibrium(result: EquilibResult) -> None:
    """Copy current Fortran equilibrium state into a result and clear thermo state."""
    try:
        result.append_output()
    finally:
        fort.resetthermo()


def _run_equilibrium_step(
    database: dict,
    condition: dict,
    unit: list[str],
    phases: list[str] | None,
) -> set[int]:
    """Run one equilibrium step and return its stable assemblage ids."""
    _equilib_single(database, condition, unit=unit, phases=phases)
    return _stable_assemblage_ids(fort.modulethermo.iassemblage)


def _equilib_point_to_cooling_point(
    point: EquilibPoint,
    liquid_phase_name: str,
) -> ScheilPoint:
    """Convert one equilibrium point into the cooling-path result shape."""
    stable_summary = point.stable_phase_summary
    phase_names = stable_summary.get("name", np.array([]))
    amounts_n = stable_summary.get("amount_n", np.array([]))
    amounts_w = stable_summary.get("amount_w", np.array([]))

    phase_amounts_n: dict[str, float] = {}
    phase_amounts_w: dict[str, float] = {}
    liquid_n = 0.0
    liquid_w = 0.0
    total_n = 0.0
    total_w = 0.0
    label_phases: list[str] = []

    for index, raw_name in enumerate(phase_names):
        phase_name = str(raw_name).strip()
        if not phase_name or phase_name.lower() == "nan":
            continue

        amount_n = float(np.round(amounts_n[index], 10))
        amount_w = float(np.round(amounts_w[index], 10))
        phase_amounts_n[phase_name] = amount_n
        phase_amounts_w[phase_name] = amount_w
        label_phases.append(phase_name)
        total_n += amount_n
        total_w += amount_w
        if phase_name == liquid_phase_name:
            liquid_n += amount_n
            liquid_w += amount_w

    if not label_phases:
        raise EquilibError("No phase appears stable during equilibrium cooling.")

    fl = liquid_n / total_n if total_n > 0.0 else 0.0
    fl_w = liquid_w / total_w if total_w > 0.0 else 0.0
    return ScheilPoint(
        T=point.T,
        P=point.P,
        G=point.G,
        H=point.H,
        S=point.S,
        Cp=point.Cp,
        fl=fl,
        fs=1.0 - fl,
        fl_w=fl_w,
        fs_w=1.0 - fl_w,
        label="+".join(sorted(label_phases)),
        cumulative_phases_n=phase_amounts_n,
        cumulative_phases_w=phase_amounts_w,
    )


def _cooling_result_from_equilibrium_path(
    equilib_result: EquilibResult,
    liquid_phase_name: str,
    unit: list[str],
) -> ScheilResult:
    """Build a Scheil-compatible cooling result from an equilibrium path."""
    result = ScheilResult.__new__(ScheilResult)
    result.context = equilib_result.context
    result.equilib_result = equilib_result
    result.liquid_phase_name = liquid_phase_name
    result.input_unit = list(unit)
    result.warnings = list(equilib_result.warnings)
    first_point = equilib_result.points[0]
    result.n_i = first_point.n_i
    result.w_i = first_point.w_i
    result.scheil_constituents = {}
    points = [
        _equilib_point_to_cooling_point(point, liquid_phase_name)
        for point in equilib_result.points
    ]
    result.data = points[0] if len(points) == 1 else points
    return result


def equilib_cooling(
    liquid_phase_name: str,
    database: dict,
    condition: dict[str, float | str],
    delta_T: float = 5.0,
    unit: list[str] | None = None,
    phases: list[str] | None = None,
    start_from_liquidus: bool = True,
    progress_callback=None,
) -> ScheilResult:
    """
    Calculate an equilibrium cooling path.

    The input shape intentionally mirrors :func:`scheil_cooling`: a liquid phase
    name, parsed database, starting condition, temperature step, units, optional
    phase selection, and an option to start just above the liquidus. Unlike
    Scheil-Gulliver cooling, the bulk composition is held fixed at every
    temperature step. The returned object uses the same cooling-path output
    shape as :func:`scheil_cooling`, while ``result.equilib_result`` stores the
    full equilibrium history.
    """
    if unit is None:
        unit = ["K", "atm", "moles"]
    if delta_T <= 0:
        raise EquilibError("delta_T must be positive for equilibrium cooling.")

    unit_local = list(unit)
    current_condition = dict(condition)
    equilib_path = EquilibResult(
        context=capture_result_context(database, condition, unit, phases)
    )

    liquidus_warning = None
    if start_from_liquidus:
        t_max = float(current_condition["T"])
        t_min = max(t_max * 0.1, _transition_search_floor(unit_local))
        try:
            liquidus = find_liquidus_transition(
                database,
                current_condition,
                liquid_phase_name,
                t_max,
                t_min,
                unit=unit_local,
                phases=phases,
            )
        except TransitionError as exc:
            liquidus_warning = liquidus_search_failure_warning(exc)
        else:
            current_condition["T"] = liquidus + 0.1

    assemblage_old = _run_equilibrium_step(
        database,
        current_condition,
        unit_local,
        phases,
    )
    _append_current_equilibrium(equilib_path)
    if liquidus_warning is not None:
        equilib_path.points[0].warnings.append(liquidus_warning)

    current_temperature = float(equilib_path.T)
    current_condition["T"] = current_temperature
    unit_local[0] = "K"

    if liquid_phase_name not in equilib_path.stable_phases.names:
        raise EquilibError(
            f"The specified liquid phase '{liquid_phase_name}' is not stable "
            "at the initial equilibrium cooling condition."
        )

    iteration_count = int(current_temperature // delta_T)
    if progress_callback is not None:
        progress_callback(0, iteration_count, "Equilibrium Cooling")
    for step_index in tqdm(
        range(iteration_count),
        desc="Equilibrium Cooling",
        colour="#8060ff",
        ascii="░▒▓",
        bar_format="{l_bar}{bar:50}{r_bar}",
        disable=progress_callback is not None,
    ):
        if progress_callback is not None:
            progress_callback(step_index + 1, iteration_count, "Equilibrium Cooling")
        next_temperature = current_temperature - delta_T
        if next_temperature <= 0:
            break

        current_condition["T"] = next_temperature
        assemblage_new = _run_equilibrium_step(
            database,
            current_condition,
            unit_local,
            phases,
        )
        fort.resetthermo()

        if assemblage_new != assemblage_old:
            transition_condition = dict(current_condition)
            transition_condition["T"] = current_temperature
            transition_temperatures = find_transitions(
                database,
                transition_condition,
                current_temperature,
                next_temperature,
                unit=unit_local,
                phases=phases,
            )
            if len(transition_temperatures) and abs(
                current_temperature - transition_temperatures[0]
            ) < 0.1:
                transition_temperatures = transition_temperatures[1:]

            for transition_temperature in transition_temperatures:
                if not (
                    next_temperature <= transition_temperature <= current_temperature
                ):
                    continue
                current_condition["T"] = float(transition_temperature)
                _run_equilibrium_step(
                    database,
                    current_condition,
                    unit_local,
                    phases,
                )
                _append_current_equilibrium(equilib_path)

        current_condition["T"] = next_temperature
        assemblage_new = _run_equilibrium_step(
            database,
            current_condition,
            unit_local,
            phases,
        )
        _append_current_equilibrium(equilib_path)

        assemblage_old = set(assemblage_new)
        current_temperature = next_temperature
        if liquid_phase_name not in equilib_path.points[-1].stable_phases.names:
            break

    return _cooling_result_from_equilibrium_path(equilib_path, liquid_phase_name, unit)


__all__ = ["equilib_cooling"]
