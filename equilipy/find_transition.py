"""Phase transition temperature search routines."""

from __future__ import annotations

from typing import List, NamedTuple, Optional, Set

import numpy as np

import equilipy.equilifort as fort

from .equilib_single import _equilib_single
from .exceptions import TransitionError


class _PhaseState(NamedTuple):
    """Stable phase ids for one transition endpoint."""

    ids: Set[int]


def _phase_state_from_assemblage(assemblage: np.ndarray) -> _PhaseState:
    """Build a transition phase state from a Fortran assemblage array."""
    stable_set = {int(phase_id) for phase_id in assemblage if int(phase_id) != 0}
    return _PhaseState(stable_set)


def capture_phase_state() -> _PhaseState:
    """Capture the current stable phase state without rerunning equilibrium."""
    return _phase_state_from_assemblage(fort.modulethermo.iassemblage.copy())


def _reset_thermo_state() -> None:
    """Clear transient Fortran thermo output when the backend provides it."""
    reset = getattr(fort, "resetthermo", None)
    if callable(reset):
        reset()


def _phase_state_is_valid(phase_state: _PhaseState | None) -> bool:
    """Return True when an equilibrium probe produced a usable phase set."""
    return phase_state is not None and bool(phase_state.ids)


def find_transitions(
    database: dict,
    condition: dict,
    T_max: float,
    T_min: float,
    unit: Optional[List[str]] = None,
    phases: Optional[List[str]] = None,
    T_tol: float = 1e-1,
    max_depth: int = 15,
    state_at_T_max: Optional[_PhaseState] = None,
    state_at_T_min: Optional[_PhaseState] = None,
) -> np.ndarray:
    """
    Find phase transition intervals with recursive bisection.

    The search stays inside the requested temperature interval. Empty
    assemblages are treated as failed probes rather than real phase states, and
    simple liquid/solid two-phase searches seed the first liquidus bracket with
    a directed cooling scan before using recursive bisection for lower
    transitions.

    Parameters
    ----------
    database : dict
        The database object.
    condition : dict
        A dictionary with condition conditions (e.g., {'P': 1, 'N': 1}).
        Temperature will be set internally.
    T_max : float
        The upper bound of the temperature search range.
    T_min : float
        The lower bound of the temperature search range.
    phases : Optional[List[str]], optional
        A list of all possible phase names. If None, all phases in the
        database are considered. By default None.
    T_tol : float, optional
        The temperature tolerance for identifying a transition interval.
        When the search interval is smaller than this, the interval boundaries
        are recorded. By default 1E-1.
    max_depth : int, optional
        The maximum recursion depth to prevent infinite loops. By default 15.

    Returns
    -------
    np.ndarray
        A sorted NumPy array containing high/low temperature pairs that bracket
        each transition. For example:
        [T1_high, T1_low, T2_high, T2_low, ...].
    """
    if unit is None:
        unit = ["K", "atm", "moles"]

    phase_state_cache: dict[float, _PhaseState] = {}

    def _known_state_for_temperature(T: float) -> Optional[_PhaseState]:
        """Return a caller-provided endpoint state when T matches it."""
        endpoint_tolerance = max(T_tol * 1e-9, 1e-10)
        if (
            state_at_T_min is not None
            and abs(float(T) - float(T_min)) <= endpoint_tolerance
            and _phase_state_is_valid(state_at_T_min)
        ):
            return state_at_T_min
        if (
            state_at_T_max is not None
            and abs(float(T) - float(T_max)) <= endpoint_tolerance
            and _phase_state_is_valid(state_at_T_max)
        ):
            return state_at_T_max
        return None

    def get_phase_state(T: float) -> _PhaseState:
        """Get a valid stable phase state at a given temperature."""
        temperature = float(T)
        known_state = _known_state_for_temperature(T)
        if known_state is not None:
            return known_state
        cached_state = phase_state_cache.get(temperature)
        if cached_state is not None:
            return cached_state

        ntp_local = condition.copy()
        ntp_local["T"] = temperature
        _reset_thermo_state()
        try:
            _equilib_single(
                database,
                ntp_local,
                unit=unit,
                phases=phases,
                include_heat_capacity=False,
            )
            assemblage = fort.modulethermo.iassemblage.copy()
            phase_state = _phase_state_from_assemblage(assemblage)
            if not phase_state.ids:
                raise TransitionError(
                    "Equilibrium calculation produced no stable phases while "
                    f"checking transition temperature T={temperature:g}."
                )
            phase_state_cache[temperature] = phase_state
            return phase_state
        finally:
            _reset_thermo_state()

    def get_valid_phase_state_near(
        T: float,
        neighbor_T: float,
    ) -> tuple[float, _PhaseState]:
        """Return a valid state at T or slightly inside the search interval."""
        temperature = float(T)
        neighbor_temperature = float(neighbor_T)
        try:
            return temperature, get_phase_state(temperature)
        except Exception as first_error:
            span = abs(neighbor_temperature - temperature)
            if span <= 0:
                raise first_error
            direction = 1.0 if neighbor_temperature > temperature else -1.0
            max_offset = max(span * 0.45, min(span, T_tol))
            offsets = [
                max(span * 1e-8, 1e-8),
                max(T_tol * 0.01, span * 1e-7, 1e-7),
                max(T_tol * 0.05, span * 1e-6, 1e-6),
                max(T_tol * 0.10, span * 1e-5, 1e-5),
                max(T_tol * 0.25, span * 1e-4, 1e-4),
                max(T_tol * 0.50, span * 1e-3, 1e-3),
                min(max(T_tol, span * 0.01), max_offset),
                min(max(T_tol * 2.0, span * 0.05), max_offset),
                min(span * 0.25, max_offset),
                max_offset,
            ]
            tried_temperatures: set[float] = {temperature}
            for offset in offsets:
                if offset <= 0:
                    continue
                trial_temperature = temperature + direction * min(offset, max_offset)
                lower = min(temperature, neighbor_temperature)
                upper = max(temperature, neighbor_temperature)
                if not lower <= trial_temperature <= upper:
                    continue
                if trial_temperature in tried_temperatures:
                    continue
                tried_temperatures.add(trial_temperature)
                try:
                    return trial_temperature, get_phase_state(trial_temperature)
                except Exception:
                    continue
            raise first_error

    def _find_recursive(
        t_low: float,
        t_high: float,
        depth: int,
        transitions_list: list[tuple[float, float]],
    ) -> None:
        """Recursively search for transitions in the interval [t_low, t_high]."""
        if depth > max_depth:
            raise TransitionError(
                "Maximum transition recursion depth reached in interval "
                f"[{t_low:g}, {t_high:g}]."
            )

        try:
            valid_low, state_low = get_valid_phase_state_near(t_low, t_high)
            valid_high, state_high = get_valid_phase_state_near(t_high, t_low)
        except Exception as exc:
            raise TransitionError(
                "Equilibrium calculation failed while searching transition range "
                f"[{t_low:g}, {t_high:g}] with phases={phases}."
            ) from exc

        set_low = state_low.ids
        set_high = state_high.ids
        if set_low == set_high:
            return

        if abs(valid_high - valid_low) < T_tol:
            transitions_list.append((valid_low, valid_high))
            return

        t_mid = (valid_low + valid_high) / 2
        try:
            valid_mid, state_mid = get_valid_phase_state_near(t_mid, valid_high)
        except Exception:
            try:
                valid_mid, state_mid = get_valid_phase_state_near(t_mid, valid_low)
            except Exception:
                transitions_list.append((valid_low, valid_high))
                return

        if (
            abs(valid_mid - valid_low) < max(T_tol * 1e-3, 1e-10)
            or abs(valid_high - valid_mid) < max(T_tol * 1e-3, 1e-10)
        ):
            transitions_list.append((valid_low, valid_high))
            return

        if state_mid.ids != set_low:
            _find_recursive(valid_low, valid_mid, depth + 1, transitions_list)
        if state_mid.ids != set_high:
            _find_recursive(valid_mid, valid_high, depth + 1, transitions_list)

    def _validated_transition_pairs(
        transition_pairs: list[tuple[float, float]],
    ) -> list[tuple[float, float]]:
        """Keep only brackets whose lower and upper phase sets differ."""
        validated_pairs: list[tuple[float, float]] = []
        for t_low, t_high in transition_pairs:
            try:
                valid_low, state_low = get_valid_phase_state_near(t_low, t_high)
                valid_high, state_high = get_valid_phase_state_near(t_high, t_low)
            except Exception as exc:
                raise TransitionError(
                    "Equilibrium calculation failed while validating transition "
                    f"range [{t_low:g}, {t_high:g}] with phases={phases}."
                ) from exc
            if state_low.ids != state_high.ids:
                validated_pairs.append((valid_low, valid_high))
        return _deduplicate_transition_pairs(validated_pairs)

    def _deduplicate_transition_pairs(
        transition_pairs: list[tuple[float, float]],
    ) -> list[tuple[float, float]]:
        """Sort transition brackets high-to-low without losing pair identity."""
        sorted_pairs = sorted(
            transition_pairs,
            key=lambda pair: max(pair),
            reverse=True,
        )
        unique_pairs: list[tuple[float, float]] = []
        for pair in sorted_pairs:
            if any(
                np.isclose(pair[0], existing[0])
                and np.isclose(pair[1], existing[1])
                for existing in unique_pairs
            ):
                continue
            unique_pairs.append(pair)
        return unique_pairs

    def _flatten_transition_pairs(
        transition_pairs: list[tuple[float, float]],
    ) -> np.ndarray:
        """Return [high, low, high, low, ...] for each validated bracket."""
        values: list[float] = []
        for t_low, t_high in transition_pairs:
            values.extend([max(t_low, t_high), min(t_low, t_high)])
        return np.asarray(values)

    def _cold_side_after_transition(
        transition_high: float,
        high_state: _PhaseState,
    ) -> tuple[float, _PhaseState] | None:
        """Return a valid cold-side state immediately below a transition."""
        span = max(0.0, float(transition_high) - float(T_min))
        if span <= 0:
            return None
        offsets = [
            T_tol * 0.5,
            T_tol,
            T_tol * 2.0,
            T_tol * 5.0,
            min(max(T_tol * 10.0, 1.0), span),
            min(max(T_tol * 20.0, 2.0), span),
        ]
        tried_temperatures: set[float] = set()
        for offset in offsets:
            if offset <= 0:
                continue
            trial_temperature = max(float(T_min), float(transition_high) - offset)
            if trial_temperature in tried_temperatures:
                continue
            tried_temperatures.add(trial_temperature)
            try:
                trial_temperature, trial_state = get_valid_phase_state_near(
                    trial_temperature,
                    transition_high,
                )
            except Exception:
                continue
            if trial_state.ids != high_state.ids:
                return trial_temperature, trial_state
        return None

    def _should_seed_first_transition() -> bool:
        """Return True for simple liquidus searches that need directed seeding."""
        if phases is None or len(phases) != 2:
            return False
        return any("LIQUID" in str(phase).upper() for phase in phases)

    try:
        _, state_min = get_valid_phase_state_near(T_min, T_max)
        _, state_max = get_valid_phase_state_near(T_max, T_min)
    except Exception as exc:
        raise TransitionError(
            "Equilibrium calculation failed while checking transition "
            f"boundaries T_min={T_min:g}, T_max={T_max:g} with phases={phases}."
        ) from exc

    if state_min.ids == state_max.ids:
        raise TransitionError(
            "No phase transitions were found in the requested interval for "
            f"T_max={T_max:g}, T_min={T_min:g}, unit={unit}, phases={phases}."
        )

    transitions: list[tuple[float, float]] = []
    first_cold_side = None
    first_transition_high: float | None = None
    if _should_seed_first_transition():
        try:
            first_transition_high = find_first_transition(
                database,
                condition,
                T_max,
                T_min,
                unit=unit,
                phases=phases,
                T_tol=T_tol,
                max_depth=max_depth,
            )
            first_high_state = get_phase_state(first_transition_high)
            first_cold_side = _cold_side_after_transition(
                first_transition_high,
                first_high_state,
            )
        except TransitionError:
            first_cold_side = None
    else:
        first_cold_side = None

    if first_cold_side is not None:
        first_transition_low, _first_low_state = first_cold_side
        transitions.append((first_transition_low, first_transition_high))

    recursive_transitions: list[tuple[float, float]] = []
    _find_recursive(T_min, T_max, 0, recursive_transitions)
    if first_cold_side is None:
        transitions.extend(recursive_transitions)
    else:
        first_transition_low, _first_low_state = first_cold_side
        reentrant_skip_width = max(T_tol * 5.0, 5.0)
        transitions.extend(
            pair
            for pair in recursive_transitions
            if max(pair) < first_transition_low - reentrant_skip_width
        )

    if transitions:
        validated_transitions = _validated_transition_pairs(transitions)
        if validated_transitions:
            return _flatten_transition_pairs(validated_transitions)

    raise TransitionError(
        "No phase transitions were found in the requested interval for "
        f"T_max={T_max:g}, T_min={T_min:g}, unit={unit}, phases={phases}."
    )


def find_first_transition(
    database: dict,
    condition: dict,
    T_max: float,
    T_min: float,
    unit: Optional[List[str]] = None,
    phases: Optional[List[str]] = None,
    T_tol: float = 1e-2,
    max_depth: int = 20,
) -> float:
    """
    Find the first phase transition point below T_max.

    This function is optimized to find only the highest temperature transition.
    A transition is defined as any change in the set of stable phase IDs.

    Args:
        database (dict): The database object.
        condition (dict): Conditions (e.g., {'P': 1, 'N': 1}). 'T' will be set.
        T_max (float): The upper bound of the temperature search.
        T_min (float): The lower bound of the temperature search.
        unit (List[str], optional): Input units. Defaults to ['K', 'atm', 'moles'].
        phases (Optional[List[str]], optional): Phases to consider.
            Defaults to None.
        T_tol (float, optional): Temperature tolerance to stop the search.
            Defaults to 1E-1.
        max_depth (int, optional): Maximum search depth to prevent errors.
            Defaults to 20.

    Returns
    -------
        float: The highest transition temperature found.
    """
    if unit is None:
        unit = ["K", "atm", "moles"]

    phase_state_cache: dict[float, _PhaseState] = {}

    def get_phase_state(T: float) -> _PhaseState:
        """Return a valid phase state for a directed liquidus probe."""
        temperature = float(T)
        cached_state = phase_state_cache.get(temperature)
        if cached_state is not None:
            return cached_state
        ntp_local = condition.copy()
        ntp_local["T"] = temperature
        _reset_thermo_state()
        try:
            _equilib_single(
                database,
                ntp_local,
                unit=unit,
                phases=phases,
                include_heat_capacity=False,
            )
            phase_state = _phase_state_from_assemblage(
                fort.modulethermo.iassemblage.copy()
            )
            if not phase_state.ids:
                raise TransitionError(
                    "Equilibrium calculation produced no stable phases while "
                    f"checking transition temperature T={temperature:g}."
                )
            phase_state_cache[temperature] = phase_state
            return phase_state
        finally:
            _reset_thermo_state()

    try:
        hot_state = get_phase_state(T_max)
    except Exception as exc:
        raise TransitionError(
            "Equilibrium calculation failed while checking transition "
            f"temperature T={T_max:g} with phases={phases}."
        ) from exc

    hot_temperature = float(T_max)
    cold_temperature = float(T_min)
    samples_per_pass = 64

    for _ in range(max_depth):
        if hot_temperature - cold_temperature <= T_tol:
            return hot_temperature

        previous_same_temperature = hot_temperature
        found_bracket = False
        for trial_temperature in np.linspace(
            hot_temperature,
            cold_temperature,
            samples_per_pass + 1,
        )[1:]:
            try:
                trial_state = get_phase_state(float(trial_temperature))
            except Exception:
                continue

            if trial_state.ids == hot_state.ids:
                previous_same_temperature = float(trial_temperature)
                continue

            hot_temperature = previous_same_temperature
            cold_temperature = float(trial_temperature)
            found_bracket = True
            break

        if found_bracket:
            continue

        raise TransitionError(
            f"No phase transition found between T_max={T_max:g} and "
            f"T_min={T_min:g} for unit={unit}, phases={phases}."
        )

    raise TransitionError(
        "Maximum transition recursion depth reached while searching "
        f"T_max={T_max:g}, T_min={T_min:g}, unit={unit}, phases={phases}."
    )


def find_liquidus_transition(
    database: dict,
    condition: dict,
    liquid_phase_name: str,
    T_max: float,
    T_min: float,
    unit: Optional[List[str]] = None,
    phases: Optional[List[str]] = None,
    T_tol: float = 1e-2,
    max_iterations: int = 32,
) -> float:
    """Compatibility wrapper for the old all-phase liquidus bisection path."""
    return find_first_transition(
        database,
        condition,
        T_max,
        T_min,
        unit=unit,
        phases=phases,
        T_tol=T_tol,
        max_depth=max_iterations,
    )
