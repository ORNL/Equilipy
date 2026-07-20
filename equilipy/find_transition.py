"""Phase transition temperature search routines."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import List, NamedTuple, Optional, Sequence, Set

import numpy as np

import equilipy.equilifort as fort
import equilipy.variables as var

from .equilib_single import _equilib_single
from .exceptions import EquilibError, SolverFailureReason, TransitionError
from .results.equilib import get_assemblage_name


class TransitionArray(np.ndarray):
    """Transition brackets with locator diagnostics attached."""

    solver_call_count: int
    solver_calls_per_transition: tuple[int, ...]
    engine: str
    locator_iterations: tuple[int, ...]
    locator_statuses: tuple[str, ...]
    marginal_phases: tuple[tuple[str, ...], ...]

    def __new__(
        cls,
        values,
        *,
        solver_call_count: int = 0,
        solver_calls_per_transition: Sequence[int] = (),
        engine: str = "driving_force",
        locator_iterations: Sequence[int] = (),
        locator_statuses: Sequence[str] = (),
        marginal_phases: Sequence[Sequence[str]] = (),
    ):
        obj = np.asarray(values, dtype=float).view(cls)
        obj.solver_call_count = int(solver_call_count)
        obj.solver_calls_per_transition = tuple(
            int(value) for value in solver_calls_per_transition
        )
        obj.engine = str(engine)
        obj.locator_iterations = tuple(int(value) for value in locator_iterations)
        obj.locator_statuses = tuple(str(value) for value in locator_statuses)
        obj.marginal_phases = tuple(
            tuple(str(phase) for phase in phases) for phases in marginal_phases
        )
        return obj

    def __array_finalize__(self, obj) -> None:
        if obj is None:
            return
        self.solver_call_count = getattr(obj, "solver_call_count", 0)
        self.solver_calls_per_transition = getattr(
            obj,
            "solver_calls_per_transition",
            (),
        )
        self.engine = getattr(obj, "engine", "driving_force")
        self.locator_iterations = getattr(obj, "locator_iterations", ())
        self.locator_statuses = getattr(obj, "locator_statuses", ())
        self.marginal_phases = getattr(obj, "marginal_phases", ())


class _PhaseState(NamedTuple):
    """Stable phase ids for one transition endpoint."""

    ids: Set[int]


@dataclass(frozen=True)
class _TransitionProbe:
    """Fortran-backed transition facts at one temperature."""

    temperature: float
    temperature_K: float
    state: _PhaseState
    phase_names: tuple[str, ...]
    phase_amounts: dict[str, float]
    phase_potentials: dict[str, float]
    phase_potential_derivatives: dict[str, float]
    component_count: int
    solver_calls: int = 1

    @property
    def phase_set(self) -> frozenset[str]:
        return frozenset(self.phase_names)

    def amount(self, phase_name: str) -> float:
        return float(self.phase_amounts.get(phase_name, 0.0))


@dataclass(frozen=True)
class _EventPrediction:
    """One predicted transition event."""

    phase_name: str
    kind: str
    temperature: float
    value: float
    derivative: float
    solver_calls: int = 0


@dataclass(frozen=True)
class _LocatedEvent:
    """One located transition bracket."""

    high: float
    low: float
    phase_name: str
    kind: str
    solver_calls: int
    iterations: int
    status: str


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
    Find phase transition intervals with driving-force prediction.

    The search stays inside the requested temperature interval.  A converged
    Fortran equilibrium supplies stable phases, phase amounts, Leveling-row
    driving forces, and entropy facts.  Candidate crossings are predicted from
    ``dDF/dT`` when the entropy plane is determined, or from one same-side exact
    driving-force re-solve when the active assemblage is underdetermined.  The
    prediction is corrected by bounded Fortran re-solves. Recursive bisection
    and directed scan retry are not production fallbacks: discontinuous
    swap/cascade events raise with a diagnostic trail for later algorithm work.

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
        The maximum correction-iteration count to prevent infinite loops. By
        default 15.

    Returns
    -------
    np.ndarray
        A sorted NumPy array containing high/low temperature pairs that bracket
        each transition. For example:
        [T1_high, T1_low, T2_high, T2_low, ...].
    """
    if unit is None:
        unit = ["K", "atm", "moles"]
    if T_tol <= 0:
        raise TransitionError("T_tol must be positive for transition search.")

    del state_at_T_max, state_at_T_min
    probe_cache: dict[float, _TransitionProbe] = {}
    total_calls = 0

    def probe(T: float) -> _TransitionProbe:
        nonlocal total_calls
        temperature = float(T)
        cached = probe_cache.get(temperature)
        if cached is not None:
            return cached
        result = _probe_transition_state(
            database,
            condition,
            temperature,
            unit=unit,
            phases=phases,
        )
        total_calls += result.solver_calls
        probe_cache[temperature] = result
        return result

    try:
        hot_probe = probe(float(T_max))
        cold_probe = probe(float(T_min))
    except Exception as exc:
        raise TransitionError(
            "Equilibrium calculation failed while checking transition "
            f"boundaries T_min={T_min:g}, T_max={T_max:g} with phases={phases}."
        ) from exc

    if hot_probe.state.ids == cold_probe.state.ids:
        raise TransitionError(
            "No phase transitions were found in the requested interval for "
            f"T_max={T_max:g}, T_min={T_min:g}, unit={unit}, phases={phases}."
        )

    located_events: list[_LocatedEvent] = []
    current_probe = hot_probe
    max_events = max(1, int(max_depth) * max(1, len(phases or current_probe.phase_names)))
    errors: list[str] = []

    for _event_index in range(max_events):
        if current_probe.temperature <= float(T_min) + max(T_tol * 1e-6, 1e-12):
            break
        if current_probe.state.ids == cold_probe.state.ids:
            break
        calls_before_event = total_calls
        try:
            event = _march_to_next_transition(
                current_probe,
                cold_probe,
                probe,
                T_tol=float(T_tol),
                max_iterations=max(3, int(max_depth)),
            )
        except TransitionError as exc:
            errors.append(str(exc))
            break
        event = _LocatedEvent(
            high=event.high,
            low=event.low,
            phase_name=event.phase_name,
            kind=event.kind,
            solver_calls=total_calls - calls_before_event,
            iterations=event.iterations,
            status=event.status,
        )
        located_events.append(event)
        if event.low <= float(T_min):
            break
        current_probe = probe(event.low)

    if not located_events:
        detail = f" Details: {' | '.join(errors)}" if errors else ""
        raise TransitionError(
            "No phase transitions were found by the driving-force predictor in "
            f"the requested interval T_max={T_max:g}, T_min={T_min:g}, "
            f"unit={unit}, phases={phases}.{detail} "
            f"{_discontinuous_event_note()}"
        )

    values: list[float] = []
    calls_per_transition: list[int] = []
    iterations: list[int] = []
    statuses: list[str] = []
    marginal_phases: list[tuple[str, ...]] = []
    for event in _deduplicate_located_events(located_events):
        values.extend([float(event.high), float(event.low)])
        calls_per_transition.append(int(event.solver_calls))
        iterations.append(int(event.iterations))
        statuses.append(str(event.status))
        marginal_phases.append((event.phase_name,))
    return TransitionArray(
        values,
        solver_call_count=total_calls,
        solver_calls_per_transition=calls_per_transition,
        engine="driving_force",
        locator_iterations=iterations,
        locator_statuses=statuses,
        marginal_phases=marginal_phases,
    )


def _directed_phase_change_bracket(
    hot_probe: _TransitionProbe,
    cold_probe: _TransitionProbe,
    probe,
    *,
    samples: int,
) -> tuple[_TransitionProbe, _TransitionProbe] | None:
    """Diagnostic-only first high-to-low phase-set change from a directed scan."""
    if samples <= 0 or hot_probe.temperature <= cold_probe.temperature:
        return None
    previous = hot_probe
    temperatures = np.linspace(
        hot_probe.temperature,
        cold_probe.temperature,
        int(samples) + 1,
    )[1:]
    for temperature in temperatures:
        current = probe(float(temperature))
        if current.state.ids != previous.state.ids:
            return previous, current
        previous = current
    return None


def diagnose_directed_phase_change_bracket(
    database: dict,
    condition: dict,
    T_max: float,
    T_min: float,
    unit: Optional[List[str]] = None,
    phases: Optional[List[str]] = None,
    *,
    samples: int = 64,
) -> dict:
    """Run the retired directed scan as an explicit diagnostic only."""
    if unit is None:
        unit = ["K", "atm", "moles"]
    probe_cache: dict[float, _TransitionProbe] = {}

    def probe(T: float) -> _TransitionProbe:
        temperature = float(T)
        cached = probe_cache.get(temperature)
        if cached is not None:
            return cached
        result = _probe_transition_state(
            database,
            condition,
            temperature,
            unit=unit,
            phases=phases,
        )
        probe_cache[temperature] = result
        return result

    hot_probe = probe(float(T_max))
    cold_probe = probe(float(T_min))
    bracket = _directed_phase_change_bracket(
        hot_probe,
        cold_probe,
        probe,
        samples=int(samples),
    )
    payload = {
        "engine": "directed_scan_diagnostic",
        "solver_call_count": len(probe_cache),
        "hot": _format_probe_state("hot", hot_probe),
        "cold": _format_probe_state("cold", cold_probe),
        "bracket": None,
    }
    if bracket is not None:
        payload["bracket"] = (
            float(bracket[0].temperature),
            float(bracket[1].temperature),
        )
        payload["bracket_trail"] = _format_transition_trail(bracket[0], bracket[1])
    return payload


def _probe_transition_state(
    database: dict,
    condition: dict,
    temperature: float,
    *,
    unit: Sequence[str],
    phases: Sequence[str] | None,
) -> _TransitionProbe:
    """Run one Fortran equilibrium and capture transition-search facts."""
    ntp_local = condition.copy()
    ntp_local["T"] = float(temperature)
    _reset_thermo_state()
    try:
        _equilib_single(
            database,
            ntp_local,
            unit=list(unit),
            phases=list(phases) if phases is not None else None,
            include_heat_capacity=False,
        )
        assemblage = fort.modulethermo.iassemblage.copy()
        phase_state = _phase_state_from_assemblage(assemblage)
        if not phase_state.ids:
            raise TransitionError(
                "Equilibrium calculation produced no stable phases while "
                f"checking transition temperature T={float(temperature):g}."
            )
        raw_phase_names = tuple(str(name).strip() for name in get_assemblage_name(assemblage))
        phase_names = tuple(name for name in raw_phase_names if name)
        phase_amounts: dict[str, float] = {}
        for phase_name, amount, entry in zip(
            raw_phase_names,
            np.asarray(fort.modulethermo.dmolesphase, dtype=float),
            np.asarray(assemblage, dtype=int),
            strict=False,
        ):
            if not phase_name or int(entry) == 0:
                continue
            phase_amounts[phase_name] = phase_amounts.get(phase_name, 0.0) + float(amount)

        phase_potentials, phase_derivatives = _phase_driving_force_facts(
            phase_state,
            assemblage,
        )
        component_count = len(getattr(var, "iElementDBIndex", ()))
        return _TransitionProbe(
            temperature=float(temperature),
            temperature_K=float(fort.modulethermoio.dtemperature),
            state=phase_state,
            phase_names=phase_names,
            phase_amounts=phase_amounts,
            phase_potentials=phase_potentials,
            phase_potential_derivatives=phase_derivatives,
            component_count=int(component_count),
        )
    finally:
        _reset_thermo_state()


def _phase_driving_force_facts(
    phase_state: _PhaseState,
    assemblage: np.ndarray,
) -> tuple[dict[str, float], dict[str, float]]:
    """Return per-phase driving-force values and entropy slopes in SI units."""
    del phase_state
    groups = _phase_row_groups()
    if not groups:
        return {}, {}

    ideal_constant = float(fort.modulethermo.didealconstant)
    temperature_K = float(fort.modulethermoio.dtemperature)
    raw_potential = np.asarray(fort.modulegemsolver.dphasepotential, dtype=float)
    raw_entropy = np.asarray(fort.modulethermo.dpartialentropy, dtype=float)
    raw_composition = np.asarray(fort.modulethermo.dlevelingcompositionspecies, dtype=float)
    energy_factor = ideal_constant * temperature_K

    entropy_plane = _active_entropy_plane(assemblage, groups, ideal_constant)
    potentials: dict[str, float] = {}
    derivatives: dict[str, float] = {}
    for phase_name, rows in groups.items():
        usable_rows = [row for row in rows if 0 <= row < raw_potential.size]
        if not usable_rows:
            continue
        row_values = raw_potential[usable_rows] * energy_factor
        local_index = int(np.argmin(row_values))
        row = int(usable_rows[local_index])
        potentials[phase_name] = float(row_values[local_index])
        if entropy_plane is None or row >= raw_entropy.size or row >= raw_composition.shape[0]:
            continue
        phase_entropy = float(raw_entropy[row] * ideal_constant)
        plane_entropy = float(np.dot(raw_composition[row], entropy_plane))
        derivatives[phase_name] = -(phase_entropy - plane_entropy)
    return potentials, derivatives


def _phase_row_groups() -> dict[str, list[int]]:
    """Map public phase names to selected Fortran species rows."""
    phase_names = [str(name).strip() for name in getattr(var, "cPhaseNameSys", ())]
    if not phase_names:
        return {}
    n_solution = len(getattr(var, "iSys2DBSoln", ()))
    boundaries = np.asarray(fort.modulethermo.nspeciesphase, dtype=int)
    groups: dict[str, list[int]] = {}
    for phase_index, phase_name in enumerate(phase_names[:n_solution]):
        if not phase_name or phase_index + 1 >= boundaries.size:
            continue
        first = int(boundaries[phase_index])
        last = int(boundaries[phase_index + 1])
        groups[phase_name] = list(range(first, last))

    selected_species = [int(index) for index in getattr(var, "iSys2DBSpecies", ())]
    compound_indices = [int(index) for index in getattr(var, "iSys2DBComp", ())]
    for offset, phase_name in enumerate(phase_names[n_solution:]):
        if not phase_name or offset >= len(compound_indices):
            continue
        try:
            row = selected_species.index(compound_indices[offset])
        except ValueError:
            continue
        groups[phase_name] = [int(row)]
    return groups


def _active_entropy_plane(
    assemblage: np.ndarray,
    groups: dict[str, list[int]],
    ideal_constant: float,
) -> np.ndarray | None:
    """Solve the active Gibbs-plane entropy vector from stable phase rows."""
    compositions: list[np.ndarray] = []
    entropies: list[float] = []
    raw_composition = np.asarray(fort.modulethermo.dlevelingcompositionspecies, dtype=float)
    raw_entropy = np.asarray(fort.modulethermo.dpartialentropy, dtype=float)
    mole_fractions = np.asarray(fort.modulethermo.dmolfraction, dtype=float)
    phase_names = [str(name).strip() for name in getattr(var, "cPhaseNameSys", ())]

    for entry in np.asarray(assemblage, dtype=int).ravel():
        if entry == 0:
            continue
        if entry < 0:
            phase_index = int(-entry - 1)
            if not (0 <= phase_index < len(phase_names)):
                continue
            phase_name = phase_names[phase_index]
            rows = groups.get(phase_name, [])
            if not rows:
                continue
            weights = np.asarray([mole_fractions[row] for row in rows], dtype=float)
            total = float(np.sum(weights))
            if total > 0.0:
                normalized_weights = weights / total
                composition = np.dot(normalized_weights, raw_composition[rows])
                entropy = float(
                    np.dot(normalized_weights, raw_entropy[rows]) * ideal_constant
                )
            else:
                composition = np.mean(raw_composition[rows], axis=0)
                entropy = float(np.mean(raw_entropy[rows]) * ideal_constant)
        else:
            phase_name = _compound_name_from_assemblage_entry(int(entry))
            rows = groups.get(phase_name, [])
            if not rows:
                continue
            row = rows[0]
            composition = np.asarray(raw_composition[row], dtype=float)
            entropy = float(raw_entropy[row] * ideal_constant)
        if composition.size == 0 or not np.all(np.isfinite(composition)):
            continue
        compositions.append(composition)
        entropies.append(entropy)

    if not compositions:
        return None
    matrix = np.vstack(compositions)
    rhs = np.asarray(entropies, dtype=float)
    if np.linalg.matrix_rank(matrix) < matrix.shape[1]:
        # A liquid-only multicomponent state does not determine the elemental
        # entropy plane.  A least-squares vector is finite but its projection
        # onto an inactive phase is arbitrary and can reverse dDF/dT.  Let the
        # transition marcher obtain an exact same-assemblage secant instead.
        return None
    try:
        solution, *_ = np.linalg.lstsq(matrix, rhs, rcond=None)
    except np.linalg.LinAlgError:
        return None
    if not np.all(np.isfinite(solution)):
        return None
    return np.asarray(solution, dtype=float)


def _compound_name_from_assemblage_entry(entry: int) -> str:
    """Return the public compound phase name for a positive assemblage entry."""
    try:
        names = get_assemblage_name(np.asarray([entry], dtype=int))
    except Exception:
        return ""
    return str(names[0]).strip() if names else ""


def _format_float(value: float) -> str:
    """Return a compact transition-trail float."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not np.isfinite(number):
        return str(number)
    return f"{number:.12g}"


def _format_phase_map(values: dict[str, float], *, limit: int = 10) -> str:
    """Format phase-keyed facts for failure trails."""
    if not values:
        return "{}"
    parts = [
        f"{name}={_format_float(value)}"
        for name, value in sorted(values.items())[:limit]
    ]
    if len(values) > limit:
        parts.append(f"...+{len(values) - limit}")
    return "{" + ", ".join(parts) + "}"


def _format_probe_state(label: str, probe: _TransitionProbe) -> str:
    """Format one bracket endpoint for a transition failure trail."""
    return (
        f"{label}(T={_format_float(probe.temperature)}, "
        f"phases={list(probe.phase_names)}, "
        f"amounts={_format_phase_map(probe.phase_amounts)}, "
        f"df={_format_phase_map(probe.phase_potentials)}, "
        f"ddf_dT={_format_phase_map(probe.phase_potential_derivatives)})"
    )


def _format_candidate(candidate: _EventPrediction) -> str:
    """Format one predicted event candidate."""
    return (
        f"{candidate.phase_name}:{candidate.kind}"
        f"(T_pred={_format_float(candidate.temperature)}, "
        f"value={_format_float(candidate.value)}, "
        f"derivative={_format_float(candidate.derivative)}, "
        f"calls={candidate.solver_calls})"
    )


def _format_transition_trail(
    hot_probe: _TransitionProbe,
    cold_probe: _TransitionProbe,
    *,
    candidates: Sequence[_EventPrediction] = (),
    history: Sequence[str] = (),
) -> str:
    """Return endpoint, candidate, and correction facts for TransitionError."""
    candidate_text = (
        "[" + "; ".join(_format_candidate(candidate) for candidate in candidates) + "]"
        if candidates
        else "[]"
    )
    history_text = "[" + "; ".join(history) + "]" if history else "[]"
    return (
        f"trail: {_format_probe_state('hot', hot_probe)}; "
        f"{_format_probe_state('cold', cold_probe)}; "
        f"candidates={candidate_text}; history={history_text}"
    )


def _discontinuous_event_note() -> str:
    """Return the named limitation note for unsolved discontinuous events."""
    return (
        "named_issue=discontinuous_swap_or_cascade; "
        "needs swap/cascade transition treatment from the v0.4 taxonomy; "
        "evidence_commit=02db1a37."
    )


def _format_march_history(history: Sequence[str]) -> str:
    """Return a compact Newton-march history string."""
    return "march=[" + "; ".join(history) + "]"


def _march_to_next_transition(
    hot_probe: _TransitionProbe,
    cold_probe: _TransitionProbe,
    probe,
    *,
    T_tol: float,
    max_iterations: int,
) -> _LocatedEvent:
    """March by Newton-predicted event temperatures before local correction."""
    previous: _TransitionProbe | None = None
    current = hot_probe
    history: list[str] = []
    max_marches = max(1, int(max_iterations)) * max(
        1,
        len(current.phase_potentials),
    )
    progress_floor = max(T_tol * 1.0e-9, 1.0e-12)

    for march_index in range(max_marches):
        candidates = _march_predictions(current, cold_probe, previous_probe=previous)
        if not candidates:
            if previous is None:
                span = current.temperature - cold_probe.temperature
                offset = min(max(span * 0.05, T_tol), span * 0.45)
                secant_temperature = current.temperature - offset
                if cold_probe.temperature < secant_temperature < current.temperature:
                    secant_probe = probe(secant_temperature)
                    history.append(
                        f"step={march_index + 1}:current="
                        f"{_format_float(current.temperature)};next="
                        f"{_format_float(secant_temperature)};"
                        "candidate=same_assemblage_secant"
                    )
                    if secant_probe.state.ids != current.state.ids:
                        return _refine_phase_change_bracket(
                            current,
                            secant_probe,
                            probe,
                            T_tol=T_tol,
                            max_iterations=max_iterations,
                        )
                    previous = current
                    current = secant_probe
                    continue
            if current.state.ids != cold_probe.state.ids:
                return _refine_phase_change_bracket(
                    current,
                    cold_probe,
                    probe,
                    T_tol=T_tol,
                    max_iterations=max_iterations,
                )
            raise TransitionError(
                "Newton march could not predict a candidate event. "
                f"{_format_transition_trail(current, cold_probe)} "
                f"{_format_march_history(history)} {_discontinuous_event_note()}"
            )
        candidate = candidates[0]
        next_temperature = float(candidate.temperature)
        history.append(
            f"step={march_index + 1}:current={_format_float(current.temperature)};"
            f"next={_format_float(next_temperature)};"
            f"candidate={_format_candidate(candidate)}"
        )
        if next_temperature >= current.temperature - progress_floor:
            raise TransitionError(
                "Newton march stopped making monotone descent progress. "
                f"{_format_transition_trail(current, cold_probe, candidates=candidates)} "
                f"{_format_march_history(history)} {_discontinuous_event_note()}"
            )
        if next_temperature <= cold_probe.temperature:
            if current.state.ids != cold_probe.state.ids:
                return _refine_phase_change_bracket(
                    current,
                    cold_probe,
                    probe,
                    T_tol=T_tol,
                    max_iterations=max_iterations,
                )
            raise TransitionError(
                "Newton march predicted outside the remaining interval. "
                f"{_format_transition_trail(current, cold_probe, candidates=candidates)} "
                f"{_format_march_history(history)} {_discontinuous_event_note()}"
            )

        try:
            next_probe = probe(next_temperature)
        except EquilibError as exc:
            if exc.reason is not SolverFailureReason.LAGRANGIAN_UNCONVERGED:
                raise
            if current.state.ids != cold_probe.state.ids:
                return _refine_phase_change_bracket(
                    current,
                    cold_probe,
                    probe,
                    T_tol=T_tol,
                    max_iterations=max_iterations,
                )
            retry_temperature = 0.5 * (current.temperature + next_temperature)
            history.append(
                f"step={march_index + 1}:failed="
                f"{_format_float(next_temperature)};retry="
                f"{_format_float(retry_temperature)};"
                "reason=lagrangian_unconverged"
            )
            next_probe = probe(retry_temperature)
        if next_probe.state.ids != current.state.ids:
            return _refine_phase_change_bracket(
                current,
                next_probe,
                probe,
                T_tol=T_tol,
                max_iterations=max_iterations,
            )
        previous = current
        current = next_probe

    raise TransitionError(
        "Newton march exceeded its structural iteration budget without "
        "bracketing a phase change. "
        f"{_format_transition_trail(current, cold_probe)} "
        f"{_format_march_history(history)} {_discontinuous_event_note()}"
    )


def _refine_phase_change_bracket(
    high_probe: _TransitionProbe,
    low_probe: _TransitionProbe,
    probe,
    *,
    T_tol: float,
    max_iterations: int,
) -> _LocatedEvent:
    """Refine an already bracketed phase-set change to the requested width."""
    high = high_probe
    low = low_probe
    calls = 0
    for iteration in range(max(1, int(max_iterations))):
        if high.temperature - low.temperature <= T_tol * (1.0 + 1.0e-12):
            phase_name, kind = _phase_change_label(high, low)
            return _LocatedEvent(
                high=float(high.temperature),
                low=float(low.temperature),
                phase_name=phase_name,
                kind=kind,
                solver_calls=calls,
                iterations=iteration,
                status="converged",
            )
        midpoint = 0.5 * (high.temperature + low.temperature)
        try:
            trial = probe(midpoint)
        except EquilibError as exc:
            if exc.reason is not SolverFailureReason.LAGRANGIAN_UNCONVERGED:
                raise
            subdivision_temperatures = (
                0.5 * (high.temperature + midpoint),
                0.5 * (midpoint + low.temperature),
            )
            successful_probes = []
            for subdivision_temperature in subdivision_temperatures:
                try:
                    successful_probes.append(probe(subdivision_temperature))
                except EquilibError as subdivision_exc:
                    if (
                        subdivision_exc.reason
                        is not SolverFailureReason.LAGRANGIAN_UNCONVERGED
                    ):
                        raise
            if not successful_probes:
                raise TransitionError(
                    "Transition refinement could not step around a "
                    "lagrangian_unconverged boundary. "
                    f"{_format_transition_trail(high, low)}"
                ) from exc
            calls += sum(item.solver_calls for item in successful_probes)
            ordered = [
                high,
                *sorted(
                    successful_probes,
                    key=lambda item: item.temperature,
                    reverse=True,
                ),
                low,
            ]
            for upper, lower in zip(ordered, ordered[1:], strict=False):
                if upper.state.ids != lower.state.ids:
                    high, low = upper, lower
                    break
            continue
        calls += trial.solver_calls
        if trial.state.ids == high.state.ids:
            high = trial
        else:
            low = trial
    if high.temperature - low.temperature <= T_tol * (1.0 + 1.0e-12):
        phase_name, kind = _phase_change_label(high, low)
        return _LocatedEvent(
            high=float(high.temperature),
            low=float(low.temperature),
            phase_name=phase_name,
            kind=kind,
            solver_calls=calls,
            iterations=max(1, int(max_iterations)),
            status="converged",
        )
    raise TransitionError(
        "Bracketed phase-set correction did not reach T_tol. "
        f"{_format_transition_trail(high_probe, low_probe)}"
    )


def _phase_change_label(
    high_probe: _TransitionProbe,
    low_probe: _TransitionProbe,
) -> tuple[str, str]:
    """Return the visible phase and kind for one phase-set bracket."""
    entrants = sorted(low_probe.phase_set - high_probe.phase_set)
    if entrants:
        return entrants[0], "phase_bracket_entry"
    exits = sorted(high_probe.phase_set - low_probe.phase_set)
    if exits:
        return exits[0], "phase_bracket_exit"
    return "phase_set_change", "phase_bracket_change"


def _march_predictions(
    hot_probe: _TransitionProbe,
    cold_probe: _TransitionProbe,
    *,
    previous_probe: _TransitionProbe | None = None,
) -> list[_EventPrediction]:
    """Predict candidate event temperatures from current-state facts only."""
    predictions: list[_EventPrediction] = []
    hot_phases = hot_probe.phase_set
    for phase_name in sorted(set(hot_probe.phase_potentials) - hot_phases):
        prediction = _predict_march_entry(
            hot_probe,
            cold_probe,
            phase_name,
            previous_probe=previous_probe,
        )
        if prediction is not None:
            predictions.append(prediction)
    for phase_name in sorted(hot_phases - cold_probe.phase_set):
        prediction = _predict_exit_event(hot_probe, cold_probe, phase_name)
        if prediction is not None:
            predictions.append(prediction)
    return sorted(predictions, key=lambda item: item.temperature, reverse=True)


def _predict_march_entry(
    hot_probe: _TransitionProbe,
    cold_probe: _TransitionProbe,
    phase_name: str,
    *,
    previous_probe: _TransitionProbe | None = None,
) -> _EventPrediction | None:
    """Predict an inactive phase entry from current driving-force facts."""
    value = float(hot_probe.phase_potentials.get(phase_name, np.nan))
    if previous_probe is not None and previous_probe.state.ids == hot_probe.state.ids:
        previous_value = float(previous_probe.phase_potentials.get(phase_name, np.nan))
        derivative = (value - previous_value) / (
            hot_probe.temperature - previous_probe.temperature
        )
        if _valid_entry_prediction(hot_probe, cold_probe, value, derivative):
            temperature = hot_probe.temperature - value / derivative
            return _EventPrediction(
                phase_name=phase_name,
                kind="march_df_secant_entry",
                temperature=float(temperature),
                value=value,
                derivative=float(derivative),
            )
        return None
    derivative = float(hot_probe.phase_potential_derivatives.get(phase_name, np.nan))
    if not _valid_entry_prediction(hot_probe, cold_probe, value, derivative):
        return None
    temperature = hot_probe.temperature - value / derivative
    return _EventPrediction(
        phase_name=phase_name,
        kind="march_driving_force_entry",
        temperature=float(temperature),
        value=value,
        derivative=derivative,
    )


def _locate_next_transition(
    hot_probe: _TransitionProbe,
    cold_probe: _TransitionProbe,
    probe,
    *,
    T_tol: float,
    max_iterations: int,
) -> _LocatedEvent:
    """Locate the next high-to-low transition without recursive bisection."""
    errors: list[str] = []
    candidates = _event_predictions(hot_probe, cold_probe, probe, T_tol=T_tol)
    if not candidates:
        raise TransitionError(
            "No driving-force or active-amount event candidate could be predicted "
            f"from T={hot_probe.temperature:g} to T={cold_probe.temperature:g}. "
            f"{_format_transition_trail(hot_probe, cold_probe)} "
            f"{_discontinuous_event_note()}"
        )
    for candidate in candidates:
        try:
            return _correct_predicted_event(
                hot_probe,
                cold_probe,
                candidate,
                probe,
                T_tol=T_tol,
                max_iterations=max_iterations,
            )
        except TransitionError as exc:
            errors.append(f"{_format_candidate(candidate)}: {exc}")
    raise TransitionError(
        "Driving-force transition correction failed for all candidates. "
        + " | ".join(errors)
        + " "
        + _format_transition_trail(hot_probe, cold_probe, candidates=candidates)
        + " "
        + _discontinuous_event_note()
    )


def _event_predictions(
    hot_probe: _TransitionProbe,
    cold_probe: _TransitionProbe,
    probe,
    *,
    T_tol: float,
) -> list[_EventPrediction]:
    """Predict candidate entry/exit temperatures in descending order."""
    predictions: list[_EventPrediction] = []
    hot_phases = hot_probe.phase_set
    cold_phases = cold_probe.phase_set
    entry_phase_names = sorted(set(hot_probe.phase_potentials) - hot_phases)
    for phase_name in entry_phase_names:
        prediction = _predict_entry_event(
            hot_probe,
            cold_probe,
            phase_name,
            probe,
            T_tol=T_tol,
        )
        if prediction is None and phase_name in cold_phases:
            amount_temperature, amount_calls = _entry_amount_prediction(
                hot_probe,
                cold_probe,
                phase_name,
                probe,
                T_tol=T_tol,
            )
            if amount_temperature is not None:
                prediction = _EventPrediction(
                    phase_name=phase_name,
                    kind="amount_entry",
                    temperature=float(amount_temperature),
                    value=cold_probe.amount(phase_name),
                    derivative=np.nan,
                    solver_calls=amount_calls,
                )
        if prediction is not None:
            if _predicted_entry_is_active(hot_probe, prediction, probe) or (
                phase_name in cold_phases
            ):
                predictions.append(prediction)
            continue
    for phase_name in sorted(hot_phases - cold_phases):
        prediction = _predict_exit_event(hot_probe, cold_probe, phase_name)
        if prediction is not None:
            predictions.append(prediction)
    return sorted(predictions, key=lambda item: item.temperature, reverse=True)


def _predicted_entry_is_active(
    hot_probe: _TransitionProbe,
    prediction: _EventPrediction,
    probe,
) -> bool:
    """Return True when a predicted entry changes the local Fortran state."""
    if prediction.kind == "amount_entry":
        return True
    try:
        trial = probe(prediction.temperature)
    except Exception:
        return False
    return trial.state.ids != hot_probe.state.ids or prediction.phase_name in trial.phase_set


def _predict_entry_event(
    hot_probe: _TransitionProbe,
    cold_probe: _TransitionProbe,
    phase_name: str,
    probe,
    *,
    T_tol: float,
) -> _EventPrediction | None:
    """Predict an inactive phase entry from driving-force slope facts."""
    value = float(hot_probe.phase_potentials.get(phase_name, np.nan))
    if not np.isfinite(value) or value <= 0.0:
        return None
    calls = 0
    derivative = np.nan
    slope = _same_side_driving_force_slope(
        hot_probe,
        cold_probe,
        phase_name,
        probe,
        T_tol=T_tol,
    )
    if slope is not None:
        derivative, slope_calls = slope
        calls += slope_calls
    else:
        derivative = float(hot_probe.phase_potential_derivatives.get(phase_name, np.nan))
    if not _valid_entry_prediction(hot_probe, cold_probe, value, derivative):
        derivative = float(hot_probe.phase_potential_derivatives.get(phase_name, np.nan))
    if not _valid_entry_prediction(hot_probe, cold_probe, value, derivative):
        return None
    temperature = hot_probe.temperature - value / derivative
    if not (cold_probe.temperature < temperature < hot_probe.temperature):
        return None
    return _EventPrediction(
        phase_name=phase_name,
        kind="driving_force_entry",
        temperature=float(temperature),
        value=value,
        derivative=float(derivative),
        solver_calls=calls,
    )


def _valid_entry_prediction(
    hot_probe: _TransitionProbe,
    cold_probe: _TransitionProbe,
    value: float,
    derivative: float,
) -> bool:
    if not np.isfinite(derivative) or abs(derivative) <= 1e-14:
        return False
    temperature = hot_probe.temperature - value / derivative
    return bool(cold_probe.temperature < temperature < hot_probe.temperature)


def _same_side_driving_force_slope(
    hot_probe: _TransitionProbe,
    cold_probe: _TransitionProbe,
    phase_name: str,
    probe,
    *,
    T_tol: float,
) -> tuple[float, int] | None:
    """Estimate dDF/dT from one same-assemblage exact Fortran re-solve."""
    span = hot_probe.temperature - cold_probe.temperature
    if span <= 0.0:
        return None
    offsets = [
        min(max(span * 0.05, T_tol), span * 0.45),
        min(max(span * 0.02, T_tol * 0.5), span * 0.35),
        min(max(span * 0.01, T_tol * 0.25), span * 0.25),
        min(max(span * 0.005, T_tol * 0.10), span * 0.10),
    ]
    base_value = float(hot_probe.phase_potentials.get(phase_name, np.nan))
    calls = 0
    for offset in offsets:
        if offset <= 0.0:
            continue
        trial_temperature = hot_probe.temperature - offset
        if not (cold_probe.temperature < trial_temperature < hot_probe.temperature):
            continue
        trial_probe = probe(trial_temperature)
        calls += trial_probe.solver_calls
        if trial_probe.state.ids != hot_probe.state.ids:
            continue
        trial_value = float(trial_probe.phase_potentials.get(phase_name, np.nan))
        if not np.isfinite(trial_value) or trial_value == base_value:
            continue
        derivative = (trial_value - base_value) / (
            trial_probe.temperature - hot_probe.temperature
        )
        if np.isfinite(derivative) and abs(derivative) > 1e-14:
            return float(derivative), calls
    return None


def _predict_exit_event(
    hot_probe: _TransitionProbe,
    cold_probe: _TransitionProbe,
    phase_name: str,
) -> _EventPrediction | None:
    """Predict a present phase exit from the Fortran lever amounts."""
    value = hot_probe.amount(phase_name)
    if not np.isfinite(value) or value <= 0.0:
        return None
    cold_value = cold_probe.amount(phase_name)
    derivative = (cold_value - value) / (cold_probe.temperature - hot_probe.temperature)
    if not np.isfinite(derivative) or abs(derivative) <= 1e-14:
        return None
    temperature = hot_probe.temperature - value / derivative
    if not (cold_probe.temperature < temperature < hot_probe.temperature):
        return None
    kind = (
        "rigid_lever_exit"
        if len(hot_probe.phase_names) == hot_probe.component_count
        else "amount_secant_exit"
    )
    return _EventPrediction(
        phase_name=phase_name,
        kind=kind,
        temperature=float(temperature),
        value=float(value),
        derivative=float(derivative),
    )


def _correct_predicted_event(
    hot_probe: _TransitionProbe,
    cold_probe: _TransitionProbe,
    candidate: _EventPrediction,
    probe,
    *,
    T_tol: float,
    max_iterations: int,
) -> _LocatedEvent:
    """Correct a predicted event with bounded Fortran re-solves."""
    high = hot_probe
    low = cold_probe
    predicted = float(candidate.temperature)
    calls = int(candidate.solver_calls)
    last_status = "predicted"
    history: list[str] = []

    for iteration in range(max_iterations):
        bracket = _local_phase_bracket(
            predicted,
            high,
            low,
            probe,
            T_tol=T_tol,
        )
        calls += bracket[2]
        history.append(
            f"iter={iteration + 1}:pred={_format_float(predicted)};"
            f"local=({_format_float(bracket[0])},{_format_float(bracket[1])});"
            f"local_calls={bracket[2]};"
            f"bracket=[{_format_float(low.temperature)},{_format_float(high.temperature)}]"
        )
        if bracket[0] is not None and bracket[1] is not None:
            return _LocatedEvent(
                high=float(bracket[0]),
                low=float(bracket[1]),
                phase_name=candidate.phase_name,
                kind=candidate.kind,
                solver_calls=calls,
                iterations=iteration + 1,
                status="converged",
            )

        wide = _wide_phase_bracket(
            predicted,
            high,
            low,
            probe,
            T_tol=T_tol,
        )
        calls += wide[2]
        history.append(
            f"iter={iteration + 1}:wide=({getattr(wide[0], 'temperature', None)},"
            f"{getattr(wide[1], 'temperature', None)});wide_calls={wide[2]}"
        )
        if wide[0] is not None and wide[1] is not None:
            high, low = wide[0], wide[1]
            endpoint = _hot_endpoint_entry_bracket(
                high,
                low,
                candidate,
                probe,
                T_tol=T_tol,
            )
            calls += endpoint[2]
            history.append(
                f"iter={iteration + 1}:endpoint="
                f"({_format_float(endpoint[0])},{_format_float(endpoint[1])});"
                f"endpoint_calls={endpoint[2]}"
            )
            if endpoint[0] is not None and endpoint[1] is not None:
                return _LocatedEvent(
                    high=float(endpoint[0]),
                    low=float(endpoint[1]),
                    phase_name=candidate.phase_name,
                    kind=candidate.kind,
                    solver_calls=calls,
                    iterations=iteration + 1,
                    status="converged",
                )
            materialized = _pairwise_materialized_entry_bracket(
                high,
                low,
                candidate,
                probe,
                T_tol=T_tol,
                max_iterations=max_iterations,
            )
            calls += materialized[2]
            history.append(
                f"iter={iteration + 1}:materialized="
                f"({_format_float(materialized[0])},"
                f"{_format_float(materialized[1])});"
                f"materialized_calls={materialized[2]}"
            )
            if materialized[0] is not None and materialized[1] is not None:
                return _LocatedEvent(
                    high=float(materialized[0]),
                    low=float(materialized[1]),
                    phase_name=candidate.phase_name,
                    kind="materialized_entry",
                    solver_calls=calls,
                    iterations=iteration + 1,
                    status="converged",
                )
            refreshed = _refresh_prediction(
                high,
                low,
                candidate.phase_name,
                candidate.kind,
                probe,
                T_tol,
            )
            calls += refreshed[1]
            if refreshed[0] is None:
                last_status = "wide_bracket_prediction_failed"
                history.append(f"iter={iteration + 1}:refresh=None")
                break
            predicted = refreshed[0]
            last_status = "wide_bracketed"
            history.append(
                f"iter={iteration + 1}:refresh={_format_float(predicted)};"
                f"refresh_calls={refreshed[1]}"
            )
            continue

        trial_temperature = _bounded_temperature(predicted, high.temperature, low.temperature)
        if trial_temperature in {high.temperature, low.temperature}:
            last_status = "bounded_temperature_on_endpoint"
            history.append(f"iter={iteration + 1}:trial_on_endpoint")
            break
        trial = probe(trial_temperature)
        calls += trial.solver_calls
        if trial.state.ids == high.state.ids:
            high = trial
        else:
            low = trial
        refreshed = _refresh_prediction(high, low, candidate.phase_name, candidate.kind, probe, T_tol)
        calls += refreshed[1]
        if refreshed[0] is None:
            last_status = "prediction_refresh_failed"
            history.append(
                f"iter={iteration + 1}:trial={_format_float(trial_temperature)};"
                f"trial_phases={list(trial.phase_names)};refresh=None"
            )
            break
        predicted = refreshed[0]
        last_status = "iterating"
        history.append(
            f"iter={iteration + 1}:trial={_format_float(trial_temperature)};"
            f"trial_phases={list(trial.phase_names)};"
            f"refresh={_format_float(predicted)};refresh_calls={refreshed[1]}"
        )

    raise TransitionError(
        f"{candidate.kind} for {candidate.phase_name} did not converge after "
        f"{max_iterations} correction iterations; status={last_status}; "
        f"last bracket=[{low.temperature:g}, {high.temperature:g}]. "
        f"{_format_transition_trail(hot_probe, cold_probe, candidates=[candidate], history=history)}"
    )


def _local_phase_bracket(
    predicted: float,
    high_probe: _TransitionProbe,
    low_probe: _TransitionProbe,
    probe,
    *,
    T_tol: float,
) -> tuple[float | None, float | None, int]:
    """Return a phase-separating local bracket around a prediction."""
    calls = 0
    center = _bounded_temperature(predicted, high_probe.temperature, low_probe.temperature)
    widths = [T_tol, T_tol * 0.5, T_tol * 2.0, T_tol * 4.0]
    for width in widths:
        if width <= 0.0:
            continue
        high_temperature = min(high_probe.temperature, center + width * 0.5)
        low_temperature = max(low_probe.temperature, center - width * 0.5)
        if high_temperature <= low_temperature:
            continue
        high = probe(high_temperature)
        low = probe(low_temperature)
        calls += high.solver_calls + low.solver_calls
        if high.state.ids != low.state.ids and high_temperature - low_temperature <= T_tol * 1.000001:
            return float(high_temperature), float(low_temperature), calls
    return None, None, calls


def _wide_phase_bracket(
    predicted: float,
    high_probe: _TransitionProbe,
    low_probe: _TransitionProbe,
    probe,
    *,
    T_tol: float,
) -> tuple[_TransitionProbe | None, _TransitionProbe | None, int]:
    """Return a directed phase-separating bracket around a rough prediction."""
    calls = 0
    center = _bounded_temperature(predicted, high_probe.temperature, low_probe.temperature)
    span = high_probe.temperature - low_probe.temperature
    if span <= 0.0:
        return None, None, calls
    widths: list[float] = []
    width = max(T_tol * 2.0, span * 1e-6)
    while width < span:
        widths.append(width)
        width *= 2.0
    widths.append(span)
    for width in widths:
        high_temperature = min(high_probe.temperature, center + width * 0.5)
        low_temperature = max(low_probe.temperature, center - width * 0.5)
        if high_temperature <= low_temperature:
            continue
        high = probe(high_temperature)
        low = probe(low_temperature)
        calls += high.solver_calls + low.solver_calls
        high_matches = high.state.ids == high_probe.state.ids
        low_matches = low.state.ids == high_probe.state.ids
        if high_matches and not low_matches:
            return high, low, calls
    return None, None, calls


def _hot_endpoint_entry_bracket(
    high_probe: _TransitionProbe,
    low_probe: _TransitionProbe,
    candidate: _EventPrediction,
    probe,
    *,
    T_tol: float,
) -> tuple[float | None, float | None, int]:
    """Return an endpoint-adjacent entry bracket from an existing phase bracket."""
    if "entry" not in candidate.kind:
        return None, None, 0
    if candidate.phase_name in high_probe.phase_set:
        return None, None, 0
    if candidate.phase_name not in low_probe.phase_set:
        return None, None, 0
    span = high_probe.temperature - low_probe.temperature
    if span <= 0.0:
        return None, None, 0
    offset = min(float(T_tol), span * 0.5)
    if offset <= 0.0:
        return None, None, 0
    trial_temperature = high_probe.temperature - offset
    if not (low_probe.temperature < trial_temperature < high_probe.temperature):
        return None, None, 0
    trial = probe(trial_temperature)
    entrant_active = (
        trial.state.ids != high_probe.state.ids
        and candidate.phase_name in trial.phase_set
    )
    if not entrant_active and trial.amount(candidate.phase_name) <= 0.0:
        return None, None, trial.solver_calls
    return float(high_probe.temperature), float(trial.temperature), trial.solver_calls


def _pairwise_materialized_entry_bracket(
    high_probe: _TransitionProbe,
    low_probe: _TransitionProbe,
    candidate: _EventPrediction,
    probe,
    *,
    T_tol: float,
    max_iterations: int,
) -> tuple[float | None, float | None, int]:
    """Refine a two-phase entry bracket when the entrant materializes abruptly."""
    if "entry" not in candidate.kind:
        return None, None, 0
    if len(high_probe.phase_potentials) != 2:
        return None, None, 0
    if candidate.phase_name in high_probe.phase_set:
        return None, None, 0
    if candidate.phase_name not in low_probe.phase_set:
        return None, None, 0
    span = high_probe.temperature - low_probe.temperature
    if span <= 0.0:
        return None, None, 0

    calls = 0
    inactive = high_probe
    active: _TransitionProbe | None = None
    offset = float(T_tol)
    for _attempt in range(max(1, int(max_iterations))):
        if offset >= span:
            break
        trial_temperature = high_probe.temperature - offset
        if not (low_probe.temperature < trial_temperature < high_probe.temperature):
            break
        trial = probe(trial_temperature)
        calls += trial.solver_calls
        if _probe_materializes_entry(trial, high_probe, candidate.phase_name):
            active = trial
            break
        inactive = trial
        offset *= 2.0

    if active is None:
        return None, None, calls

    for _attempt in range(max(1, int(max_iterations))):
        if inactive.temperature - active.temperature <= T_tol * (1.0 + 1.0e-12):
            return float(inactive.temperature), float(active.temperature), calls
        midpoint = 0.5 * (inactive.temperature + active.temperature)
        trial = probe(midpoint)
        calls += trial.solver_calls
        if _probe_materializes_entry(trial, high_probe, candidate.phase_name):
            active = trial
        else:
            inactive = trial
    if inactive.temperature - active.temperature <= T_tol * (1.0 + 1.0e-12):
        return float(inactive.temperature), float(active.temperature), calls
    return None, None, calls


def _probe_materializes_entry(
    trial: _TransitionProbe,
    high_probe: _TransitionProbe,
    phase_name: str,
) -> bool:
    """Return True when one probe has the candidate entrant active."""
    return (
        trial.state.ids != high_probe.state.ids
        and (
            phase_name in trial.phase_set
            or trial.amount(phase_name) > 0.0
        )
    )


def _bounded_temperature(temperature: float, high: float, low: float) -> float:
    """Clamp a trial temperature strictly inside a high-to-low interval."""
    span = float(high) - float(low)
    if span <= 0.0:
        return float(high)
    margin = max(span * 1e-10, 1e-12)
    return float(min(float(high) - margin, max(float(low) + margin, temperature)))


def _refresh_prediction(
    high_probe: _TransitionProbe,
    low_probe: _TransitionProbe,
    phase_name: str,
    kind: str,
    probe,
    T_tol: float,
) -> tuple[float | None, int]:
    """Recompute an event prediction after one correction probe."""
    if "entry" in kind:
        amount_prediction = _entry_amount_prediction(
            high_probe,
            low_probe,
            phase_name,
            probe,
            T_tol=T_tol,
        )
        if amount_prediction[0] is not None:
            return amount_prediction
        prediction = _predict_entry_event(
            high_probe,
            low_probe,
            phase_name,
            probe,
            T_tol=T_tol,
        )
    else:
        prediction = _predict_exit_event(high_probe, low_probe, phase_name)
    if prediction is None:
        return None, 0
    return float(prediction.temperature), int(prediction.solver_calls)


def _entry_amount_prediction(
    high_probe: _TransitionProbe,
    low_probe: _TransitionProbe,
    phase_name: str,
    probe,
    *,
    T_tol: float,
) -> tuple[float | None, int]:
    """Predict an entry boundary from cold-side entrant amounts."""
    if low_probe.amount(phase_name) <= 0.0:
        return None, 0
    span = high_probe.temperature - low_probe.temperature
    if span <= 0.0:
        return None, 0
    calls = 0
    active_probes: list[_TransitionProbe] = []
    fractions = (
        0.02,
        0.05,
        0.10,
        0.20,
        0.35,
        0.50,
        0.70,
        0.85,
        0.92,
        0.96,
        0.98,
        0.99,
    )
    for fraction in fractions:
        trial_temperature = high_probe.temperature - span * fraction
        if not (low_probe.temperature < trial_temperature < high_probe.temperature):
            continue
        trial = probe(trial_temperature)
        calls += trial.solver_calls
        if trial.state.ids == high_probe.state.ids:
            continue
        if trial.amount(phase_name) > 0.0:
            active_probes.append(trial)
            if len(active_probes) >= 2:
                break
    active_probes.append(low_probe)
    unique: list[_TransitionProbe] = []
    seen: set[float] = set()
    for item in sorted(active_probes, key=lambda probe_item: probe_item.temperature, reverse=True):
        if item.temperature in seen:
            continue
        seen.add(item.temperature)
        unique.append(item)
    if len(unique) < 2:
        endpoint_prediction = low_probe.temperature + min(T_tol * 0.5, span * 0.5)
        if high_probe.temperature > endpoint_prediction > low_probe.temperature:
            return float(endpoint_prediction), calls
        return None, calls
    near, far = unique[0], unique[1]
    near_amount = near.amount(phase_name)
    far_amount = far.amount(phase_name)
    derivative = (near_amount - far_amount) / (near.temperature - far.temperature)
    if not np.isfinite(derivative) or abs(derivative) <= 1e-14:
        endpoint_prediction = low_probe.temperature + min(T_tol * 0.5, span * 0.5)
        if high_probe.temperature > endpoint_prediction > low_probe.temperature:
            return float(endpoint_prediction), calls
        return None, calls
    predicted = near.temperature - near_amount / derivative
    if high_probe.temperature > predicted > low_probe.temperature:
        return float(predicted), calls
    endpoint_prediction = low_probe.temperature + min(T_tol * 0.5, span * 0.5)
    if high_probe.temperature > endpoint_prediction > low_probe.temperature:
        return float(endpoint_prediction), calls
    return None, calls


def _deduplicate_located_events(events: list[_LocatedEvent]) -> list[_LocatedEvent]:
    """Sort transition brackets high-to-low without losing event diagnostics."""
    sorted_events = sorted(events, key=lambda event: event.high, reverse=True)
    unique: list[_LocatedEvent] = []
    for event in sorted_events:
        if any(
            np.isclose(event.high, existing.high)
            and np.isclose(event.low, existing.low)
            for existing in unique
        ):
            continue
        unique.append(event)
    return unique


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

    This function returns the hot-side temperature of the highest transition
    bracket produced by :func:`find_transitions`.

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

    probe_cache: dict[float, _TransitionProbe] = {}

    def probe(T: float) -> _TransitionProbe:
        temperature = float(T)
        cached = probe_cache.get(temperature)
        if cached is not None:
            return cached
        result = _probe_transition_state(
            database,
            condition,
            temperature,
            unit=unit,
            phases=phases,
        )
        probe_cache[temperature] = result
        return result

    try:
        hot_probe = probe(float(T_max))
        cold_probe = probe(float(T_min))
    except Exception as exc:
        raise TransitionError(
            "Equilibrium calculation failed while checking transition "
            f"boundaries T_min={T_min:g}, T_max={T_max:g} with phases={phases}."
        ) from exc
    if hot_probe.state.ids == cold_probe.state.ids:
        raise TransitionError(
            f"No phase transition found between T_max={T_max:g} and "
            f"T_min={T_min:g} for unit={unit}, phases={phases}."
        )
    event = _march_to_next_transition(
        hot_probe,
        cold_probe,
        probe,
        T_tol=float(T_tol),
        max_iterations=max(3, int(max_depth)),
    )
    return float(event.high)


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
    """Find the first transition by marching down from a liquid hot state.

    A liquidus search does not need a solved equilibrium at an arbitrary cold
    search bound.  In multicomponent systems that cold assemblage can be much
    harder to minimize than the liquidus neighborhood and used to abort an
    otherwise well-posed search.  ``T_min`` is therefore only a lower bound for
    hot-side driving-force predictions; every physical probe starts at
    ``T_max`` and marches toward the first phase entry.
    """
    if unit is None:
        unit = ["K", "atm", "moles"]
    if T_min >= T_max:
        raise TransitionError(
            f"Liquidus search requires T_min < T_max; got {T_min:g} and {T_max:g}."
        )

    probe_cache: dict[float, _TransitionProbe] = {}

    def probe(T: float) -> _TransitionProbe:
        temperature = float(T)
        cached = probe_cache.get(temperature)
        if cached is not None:
            return cached
        result = _probe_transition_state(
            database,
            condition,
            temperature,
            unit=unit,
            phases=phases,
        )
        probe_cache[temperature] = result
        return result

    try:
        hot_probe = probe(float(T_max))
    except Exception as exc:
        raise TransitionError(
            "Equilibrium calculation failed at the hot liquidus boundary "
            f"T_max={T_max:g} with phases={phases}."
        ) from exc

    # From-liquidus semantics: when the caller's hot bound is not yet a
    # liquid-only state, ascend until one is found and anchor the march
    # there.  The liquidus may legitimately sit above the requested T_max.
    hot_T = float(T_max)
    ascent_step = max(50.0, 0.05 * abs(hot_T))
    for _ in range(12):
        if hot_probe.phase_set == {liquid_phase_name}:
            break
        hot_T += ascent_step
        ascent_step *= 1.6
        try:
            hot_probe = probe(hot_T)
        except Exception as exc:
            raise TransitionError(
                "Liquidus ascent failed while heating toward a liquid-only "
                f"state at T={hot_T:g}; last stable phases="
                f"{sorted(hot_probe.phase_set)}."
            ) from exc
    hot_phases = hot_probe.phase_set
    if hot_phases != {liquid_phase_name}:
        raise TransitionError(
            "No liquid-only state found while heating from "
            f"T_max={T_max:g} up to T={hot_T:g}; stable phases at the top="
            f"{sorted(hot_phases)}."
        )

    # The march only needs a temperature floor.  Copying the hot facts avoids
    # solving the unrelated, potentially ill-conditioned cold assemblage.
    lower_bound = replace(
        hot_probe,
        temperature=float(T_min),
        solver_calls=0,
    )
    try:
        event = _march_to_next_transition(
            hot_probe,
            lower_bound,
            probe,
            T_tol=float(T_tol),
            max_iterations=max(3, int(max_iterations)),
        )
    except Exception as exc:
        if isinstance(exc, TransitionError):
            detail = str(exc)
        else:
            detail = f"{type(exc).__name__}: {exc}"
        raise TransitionError(
            "Hot-side liquidus march failed between "
            f"T_max={T_max:g} and T_min={T_min:g} with phases={phases}. "
            f"Cause: {detail}"
        ) from exc
    return float(event.high)
