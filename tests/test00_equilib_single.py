"""Public equilibrium and phase-selection smoke tests.

Author:
    S.Y. Kwon
Date:
    06/30/2026
Tests:
    1. test_equilib_single_al_cu_si: Convergence check for G and public result
       context on a small Al-Cu-Si equilibrium point.
    2. test_equilib_single_prunes_selected_phases_for_zero_components: Selected
       phase pruning when a component amount is zero.
    3. test_equilib_batch_rejects_unknown_phase: Clear public error handling
       for unknown selected phases.
    4. test_property_solution_liquid_matches_single_phase_equilibrium:
       Consistency between ``property_solution`` and one-phase equilibrium.
"""

import numpy as np
import pytest

import equilipy as eq


def test_equilib_single_al_cu_si(alcumgsi_fs73_database: dict) -> None:
    """Solve one Al-Cu-Si equilibrium and check the stable phases."""
    condition = {
        "T": 700.0,
        "P": 1.0,
        "Al": 0.060606061,
        "Cu": 0.42424242,
        "Si": 0.515151515,
    }

    result = eq.equilib_single(
        alcumgsi_fs73_database,
        condition,
        unit=["K", "atm", "moles"],
    )

    assert result.G == pytest.approx(-35245.83605808, abs=1e-2)
    assert result.T == pytest.approx(700.0)
    assert result.context.component_names == ["Al", "Cu", "Si"]
    assert result.stable_phases.names


def test_equilib_single_prunes_selected_phases_for_zero_components(
    alcumgsi_fs73_database: dict,
) -> None:
    """Zero-amount components prune the selected phase list."""
    selected_phases = list(eq.list_phases(alcumgsi_fs73_database, ["Al", "Cu", "Si"]))
    active_phases = set(eq.list_phases(alcumgsi_fs73_database, ["Al", "Cu"]))

    result = eq.equilib_single(
        alcumgsi_fs73_database,
        {"T": 700.0, "P": 1.0, "Al": 0.5, "Cu": 0.5, "Si": 0.0},
        unit=["K", "atm", "moles"],
        phases=selected_phases,
    )

    assert np.isfinite(result.G)
    assert set(result.stable_phases.names).issubset(active_phases)


def test_equilib_batch_rejects_unknown_phase(alcumgsi_fs73_database: dict) -> None:
    """equilib_batch raises for a phase name not in the database."""
    condition = {
        "T": [700.0],
        "P": [1.0],
        "Al": [0.060606061],
        "Cu": [0.42424242],
        "Si": [0.515151515],
    }

    with pytest.raises(eq.InputConditionError, match="NOT_A_PHASE"):
        eq.equilib_batch(
            alcumgsi_fs73_database,
            condition,
            unit=["K", "atm", "moles"],
            phases=["LIQUID", "NOT_A_PHASE"],
            n_cpu=1,
        )


def test_property_solution_liquid_matches_single_phase_equilibrium(
    liu_alfesi_dat_database: dict,
) -> None:
    """property_solution on LIQUID matches the single-phase equilibrium."""
    condition = {"T": 1000.0, "P": 1.0, "Al": 0.7, "Fe": 0.2, "Si": 0.1}

    direct = eq.property_solution(
        liu_alfesi_dat_database,
        "LIQUID",
        condition,
        unit=["K", "atm", "moles"],
    )
    equilibrium = eq.equilib_single(
        liu_alfesi_dat_database,
        condition,
        unit=["K", "atm", "moles"],
        phases=["LIQUID"],
    )
    phase = equilibrium.phase("LIQUID")

    assert direct.G == pytest.approx(phase.G)
    assert direct.H == pytest.approx(phase.H)
    assert direct.S == pytest.approx(phase.S)
    assert direct.Cp == pytest.approx(phase.Cp)
    assert direct.phase("LIQUID").elements_x == pytest.approx(
        {"Al": 0.7, "Fe": 0.2, "Si": 0.1}
    )
