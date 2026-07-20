"""Transition-finding algorithm regressions.

Author:
    S.Y. Kwon
Date:
    07/15/2026
Tests:
    1. test_driving_force_locator_matches_v032_bisection_digest:
       The replacement transition locator matches the v0.3.2 bisection digest
       centers on narrow Al-Si windows while using fewer Fortran solves.
    2. test_find_first_transition_uses_predictor_without_scan_fallback:
       Broad Scheil liquidus startup uses the predictor/corrector directly.
    3. test_directed_scan_survives_only_as_explicit_diagnostic:
       The retired directed scan is available only through an explicit
       diagnostic helper, not production transition routing.
    4. test_failure_trail_names_discontinuous_swap_cascade_issue:
       Unsmooth/discontinuous failures carry endpoint, candidate, and v0.4
       taxonomy evidence in their TransitionError.
    5. test_hallstedt_long_range_newton_marching_benchmark:
       Long-range Newton marching closes the Hallstedt Al-2.5wt%Fe falsifier
       where single-shot predict/correct lost to legacy scan.
"""

from __future__ import annotations

import time
from pathlib import Path

import pytest

import equilipy as eq
import equilipy.find_transition as transition_module


REPO_ROOT = Path(__file__).resolve().parents[1]
HALLSTEDT_TDB = (
    REPO_ROOT
    / "_ongoing_development"
    / "archive"
    / "closed_tracks_2026-07"
    / "devop_tdb_parser"
    / "TDB_parsers"
    / "CompatabilityTest"
    / "Compatible"
    / "2023Hallstedt_AlCoCrFeMnNiVC.tdb"
)


def _dev_fixture(path: Path) -> Path:
    """Development baseline fixture; released trees skip tests that need it."""
    if not path.exists():
        pytest.skip(f"dev baseline fixture not present: {path.name}")
    return path



@pytest.mark.parametrize(
    (
        "label",
        "condition",
        "phases",
        "T_max",
        "T_min",
        "unit",
        "T_tol",
        "legacy_center",
        "legacy_call_count",
        "max_call_count",
    ),
    [
        (
            "alsi_binary_liquidus_narrow_window",
            {"T": 950.0, "P": 1.0, "Al": 0.92, "Si": 0.08},
            ["LIQUID", "FCC_A1"],
            890.0,
            870.0,
            ["K", "atm", "mole fraction"],
            0.25,
            880.6162109375,
            96,
            30,
        ),
        (
            "alsi_scheil_spot_first_event_narrow_window",
            {"T": 700.0, "P": 1.0, "Al": 94.0, "Si": 6.0},
            ["LIQUID", "FCC_A1", "DIAMOND_A4"],
            630.0,
            618.0,
            ["C", "atm", "grams"],
            0.25,
            623.15625,
            8,
            30,
        ),
        (
            "alsi_scheil_spot_second_event_narrow_window",
            {"T": 700.0, "P": 1.0, "Al": 94.0, "Si": 6.0},
            ["LIQUID", "FCC_A1", "DIAMOND_A4"],
            582.0,
            572.0,
            ["C", "atm", "grams"],
            0.25,
            577.078125,
            8,
            10,
        ),
    ],
)
def test_driving_force_locator_matches_v032_bisection_digest(
    liu_alfesi_dat_database: dict,
    label: str,
    condition: dict[str, float],
    phases: list[str],
    T_max: float,
    T_min: float,
    unit: list[str],
    T_tol: float,
    legacy_center: float,
    legacy_call_count: int,
    max_call_count: int,
) -> None:
    """DF-predict locator parity against the v0.3.2 bisection digest."""
    transitions = eq.find_transitions(
        liu_alfesi_dat_database,
        condition,
        T_max,
        T_min,
        unit=unit,
        phases=phases,
        T_tol=T_tol,
        max_depth=12,
    )

    assert getattr(transitions, "engine") == "driving_force"
    assert getattr(transitions, "locator_statuses") == ("converged",)
    assert len(transitions) == 2, label
    high, low = (float(transitions[0]), float(transitions[1]))
    assert high > low
    assert high - low <= T_tol * (1.0 + 1.0e-12)
    assert (high + low) / 2.0 == pytest.approx(legacy_center, abs=0.5)
    solver_call_count = getattr(transitions, "solver_call_count")
    if legacy_call_count > max_call_count:
        assert solver_call_count < legacy_call_count
    assert solver_call_count <= max_call_count


def _all_alfesi_phases(database: dict) -> list[str]:
    """Return all real Al-Fe-Si phases for public transition probes."""
    return [
        phase
        for phase in eq.list_phases(database, ["Al", "Fe", "Si"])
        if not phase.endswith("_SER(s)")
    ]


def test_find_first_transition_uses_predictor_without_scan_fallback(
    liu_alfesi_dat_database: dict,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Broad Scheil liquidus startup does not call the retired scan fallback."""
    phases = _all_alfesi_phases(liu_alfesi_dat_database)
    condition = {"T": 1600.0, "P": 1.0, "Al": 91.0, "Fe": 2.5, "Si": 6.5}
    call_count = 0
    original = transition_module._equilib_single

    def counted_equilib_single(*args, **kwargs):
        nonlocal call_count
        call_count += 1
        return original(*args, **kwargs)

    monkeypatch.setattr(transition_module, "_equilib_single", counted_equilib_single)
    monkeypatch.setattr(
        transition_module,
        "_directed_phase_change_bracket",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("directed scan is diagnostic-only")
        ),
    )

    liquidus = eq.find_first_transition(
        liu_alfesi_dat_database,
        condition,
        1600.0,
        160.0,
        unit=["C", "atm", "grams"],
        phases=phases,
        T_tol=1.0e-2,
        max_depth=20,
    )

    assert liquidus == pytest.approx(657.6265, abs=0.5)
    assert call_count <= 300


def test_directed_scan_survives_only_as_explicit_diagnostic(
    liu_alfesi_dat_database: dict,
) -> None:
    """The retired scan is diagnostic-only and reports its own engine."""
    diagnostic = transition_module.diagnose_directed_phase_change_bracket(
        liu_alfesi_dat_database,
        {"T": 1600.0, "P": 1.0, "Al": 91.0, "Fe": 2.5, "Si": 6.5},
        1600.0,
        160.0,
        unit=["C", "atm", "grams"],
        phases=_all_alfesi_phases(liu_alfesi_dat_database),
        samples=64,
    )

    assert diagnostic["engine"] == "directed_scan_diagnostic"
    assert diagnostic["bracket"] is not None


def test_failure_trail_names_discontinuous_swap_cascade_issue(
    liu_alfesi_dat_database: dict,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Unsmooth predictor failure reports the v0.4 swap/cascade limitation."""
    monkeypatch.setattr(transition_module, "_event_predictions", lambda *_args, **_kwargs: [])
    monkeypatch.setattr(transition_module, "_march_predictions", lambda *_args, **_kwargs: [])

    def blocked_refiner(high_probe, low_probe, *_args, **_kwargs):
        raise transition_module.TransitionError(
            transition_module._format_transition_trail(high_probe, low_probe)
        )

    monkeypatch.setattr(
        transition_module,
        "_refine_phase_change_bracket",
        blocked_refiner,
    )

    with pytest.raises(transition_module.TransitionError) as raised:
        eq.find_transitions(
            liu_alfesi_dat_database,
            {"T": 950.0, "P": 1.0, "Al": 0.92, "Si": 0.08},
            890.0,
            870.0,
            unit=["K", "atm", "mole fraction"],
            phases=["LIQUID", "FCC_A1"],
            T_tol=0.25,
            max_depth=12,
        )

    message = str(raised.value)
    assert "named_issue=discontinuous_swap_or_cascade" in message
    assert "evidence_commit=02db1a37" in message
    assert "trail: hot(" in message
    assert "df={" in message


def test_hallstedt_long_range_newton_marching_benchmark(record_property) -> None:
    """Long-range descent falsifier for the retired universal-speed hypothesis.

    Claude's 07/16/2026 benchmark falsified "predict-correct is universally
    faster": near-event DF prediction beat bisection, but far-above-liquidus
    single-shot correction lost to the legacy directed scan. Newton marching
    must find the Hallstedt Al-2.5wt%Fe liquidus and eutectic accurately while
    recording elapsed time and solve-count facts.
    """
    database = eq.read_tdb(_dev_fixture(HALLSTEDT_TDB))
    phases = [
        phase
        for phase in eq.list_phases(database, ["Al", "Fe"])
        if not phase.endswith("_SER(s)")
    ]
    condition = {"T": 2273.15, "P": 1.0, "Al": 97.5, "Fe": 2.5}

    started = time.perf_counter()
    transitions = eq.find_transitions(
        database,
        condition,
        2273.15,
        800.0,
        unit=["K", "atm", "wt%"],
        phases=phases,
        T_tol=0.1,
        max_depth=20,
    )
    elapsed_s = time.perf_counter() - started

    centers = [
        (float(transitions[index]) + float(transitions[index + 1])) / 2.0
        for index in range(0, len(transitions), 2)
    ]
    record_property("elapsed_s", elapsed_s)
    record_property("solver_call_count", getattr(transitions, "solver_call_count"))
    record_property(
        "solver_calls_per_transition",
        getattr(transitions, "solver_calls_per_transition"),
    )

    assert centers == pytest.approx([951.63, 926.44], abs=0.5)
    assert getattr(transitions, "solver_call_count") <= 30
    assert elapsed_s < 4.0
