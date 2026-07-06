"""Public solution-model and active-component smoke tests.

Author:
    S.Y. Kwon
Date:
    06/30/2026
Tests:
    1. test_tdb_order_disorder_single_phase_equilibrium: ThermoCalc reference
       check for FCC_4SL and BCC_B2 G/H/S/Cp values.
    2. test_tdb_dispart_pure_ordered_equilibrium_uses_disordered_boundary:
       DIS_PART pure ordered endmember inheritance check.
    3. test_cas_subq_only_single_phase_equilibrium: FactSage reference check
       for three CAS_SUBQ-only G values.
    4. test_cas_all_phase_equilibrium_matches_factsage_row_224: FactSage row
       check for CAS all-phase stable phases and G.
    5. test_saved_databases_expose_rkmp_and_qkto_models: Confirm public
       databases still expose RKMP and QKTO model coverage.
    6. test_expand_condition_species_accepts_formula_and_oxide_components:
       Formula and oxide pseudo-component expansion checks.
    7. test_active_component_driving_force_uses_component_plane: Active
       component driving-force plane check.
    8. test_active_component_rejects_invalid_basis_or_negative_coefficients:
       Active-component projection rejection checks.
"""

import numpy as np
import pytest

import equilipy as eq
from equilipy.active_component import ActiveComponentBasis
from equilipy.composition import (
    expand_condition_species,
    parse_formula_stoichiometry,
)

FACTSAGE_SUBQ_G_ABS_TOL = 2e-1
ORDER_DISORDER_G_ABS_TOL = 2e-1


# ThermoCalc TC-Python 2025b reference:
# _ongoing_development/debug_backend/liu_alfe_tcpython_single_point_refs.json
@pytest.mark.parametrize(
    (
        "phase_name",
        "condition",
        "expected_g",
        "expected_h",
        "expected_s",
        "expected_cp",
    ),
    [
        (
            "FCC_4SL",
            {"T": 900.0, "P": 1.0, "Al": 0.25, "Fe": 0.75},
            -55434.49399284888,
            -229.6510163415005,
            61.33871441834154,
            38.706322,
        ),
        (
            "BCC_B2",
            {"T": 900.0, "P": 1.0, "Al": 0.5, "Fe": 0.5},
            -63052.23133623364,
            -9986.691666617415,
            58.961710744018035,
            34.812987,
        ),
    ],
)
def test_tdb_order_disorder_single_phase_equilibrium(
    liu_alfesi_tdb_database: dict,
    phase_name: str,
    condition: dict[str, float],
    expected_g: float,
    expected_h: float,
    expected_s: float,
    expected_cp: float,
) -> None:
    """TDB order/disorder single-phase equilibrium matches references."""
    result = eq.equilib_single(
        liu_alfesi_tdb_database,
        condition,
        phases=[phase_name],
        unit=["K", "atm", "moles"],
        include_heat_capacity=True,
    )

    assert result.stable_phases.names == [phase_name]
    assert result.G == pytest.approx(expected_g, abs=ORDER_DISORDER_G_ABS_TOL)
    assert result.H == pytest.approx(expected_h, abs=5e-1)
    assert result.S == pytest.approx(expected_s, abs=1e-3)
    assert result.Cp == pytest.approx(expected_cp, abs=5e-1)


def test_tdb_liquid_al_keeps_high_power_gibbs_term(
    liu_alfesi_tdb_database: dict,
) -> None:
    """Check TDB custom Gibbs powers are preserved for liquid Al.

    :author: S.Y. Kwon
    :date: 06/30/2026
    :task:
        Validate the ``GLIQAL`` branch just below 933.47 K, where the
        ``7.9337e-20*T**7`` contribution is about 49 J/mol. This guards
        against dropping small custom coefficients whose high powers are
        thermodynamically significant.
    """
    result = eq.equilib_single(
        liu_alfesi_tdb_database,
        {"T": 933.4269272672689, "P": 1.0, "Al": 1.0},
        phases=["LIQUID"],
        unit=["K", "atm", "moles"],
        include_heat_capacity=True,
    )

    assert result.stable_phases.names == ["LIQUID"]
    assert result.G == pytest.approx(-37844.93762983438, abs=1e-4)
    assert result.H == pytest.approx(28838.637016324494, abs=1e-4)
    assert result.S == pytest.approx(71.4395232215808, abs=1e-6)
    assert result.Cp == pytest.approx(31.747612, abs=1e-4)


def test_tdb_dispart_pure_ordered_equilibrium_uses_disordered_boundary(
    liu_alfesi_tdb_database: dict,
) -> None:
    """Check pure ordered equilibrium uses the DIS_PART boundary.

    :author: S.Y. Kwon
    :date: 07/01/2026
    :task:
        Validate that a pure-Al ``FCC_4SL`` equilibrium with its ``DIS_PART``
        partner present does not fall back to the raw implicit zero ordered
        endmember. The ordered corner should use the matching ``FCC_A1``
        ``Al:VA`` boundary thermodynamics.
    """
    condition = {"T": 773.15, "P": 1.0, "Al": 1.0}

    disordered = eq.property_solution(
        liu_alfesi_tdb_database,
        "FCC_A1",
        condition,
        unit=["K", "atm", "moles"],
        include_heat_capacity=True,
    )
    result = eq.equilib_single(
        liu_alfesi_tdb_database,
        condition,
        phases=["FCC_A1", "FCC_4SL"],
        unit=["K", "atm", "moles"],
        include_heat_capacity=True,
    )

    assert set(result.stable_phases.names).issubset({"FCC_A1", "FCC_4SL"})
    assert result.G == pytest.approx(disordered.G, abs=1e-6)
    assert result.H == pytest.approx(disordered.H, abs=1e-6)
    assert result.S == pytest.approx(disordered.S, abs=1e-8)
    assert result.Cp == pytest.approx(disordered.Cp, abs=1e-4)


# FactSage CAS_SUBQ-only workbook reference rows:
# _ongoing_development/cas_subq_only_comparison.csv
@pytest.mark.parametrize(
    ("condition", "factsage_row", "expected_g"),
    [
        (
            {"T": 1600.0, "P": 1.0, "SiO2": 1.0},
            134,
            -1105341.9,
        ),
        (
            {"T": 1600.0, "P": 1.0, "Al2O3": 0.5, "SiO2": 0.5},
            139,
            -1539360.4,
        ),
        (
            {"T": 1600.0, "P": 1.0, "CaO": 0.2, "Al2O3": 0.3, "SiO2": 0.5},
            158,
            -1336548.4,
        ),
    ],
)
def test_cas_subq_only_single_phase_equilibrium(
    cas_fs83_database: dict,
    condition: dict[str, float],
    factsage_row: int,
    expected_g: float,
) -> None:
    """CAS SUBQ-only single-phase equilibrium solves correctly."""
    result = eq.equilib_single(
        cas_fs83_database,
        condition,
        phases=["Slag-liq"],
        unit=["C", "atm", "moles"],
        include_heat_capacity=False,
    )

    assert result.stable_phases.names == ["Slag-liq"]
    assert result.G == pytest.approx(
        expected_g,
        abs=FACTSAGE_SUBQ_G_ABS_TOL,
    ), f"FactSage row {factsage_row}"


def test_cas_all_phase_equilibrium_matches_factsage_row_224(
    cas_fs83_database: dict,
) -> None:
    """CAS all-phase equilibrium matches FactSage row 224."""
    phases = [
        phase
        for phase in eq.list_phases(cas_fs83_database, ["Ca", "Al", "Si", "O"])
        if not phase.endswith("_SER(s)")
    ]

    result = eq.equilib_single(
        cas_fs83_database,
        {"T": 1400.0, "P": 1.0, "CaO": 0.2, "Al2O3": 0.3, "SiO2": 0.5},
        phases=phases,
        unit=["C", "atm", "moles"],
        include_heat_capacity=False,
    )

    assert result.stable_phases.names == [
        "Slag-liq",
        "Al2O3_corundum(alpha)(s4)",
    ]
    assert result.G == pytest.approx(-1294910.4, abs=FACTSAGE_SUBQ_G_ABS_TOL)


def test_saved_databases_expose_rkmp_and_qkto_models(
    liu_alfesi_dat_database: dict,
    cas_fs83_database: dict,
) -> None:
    """Saved databases expose RKMP and QKTO model names."""
    assert any(
        name.strip() == "LIQUID" and phase_type.strip() == "RKMP"
        for name, phase_type in zip(
            liu_alfesi_dat_database["cSolnPhaseNameCS"],
            liu_alfesi_dat_database["cSolnPhaseTypeCS"],
            strict=False,
        )
    )
    assert any(
        phase_type.strip() == "QKTO"
        for phase_type in cas_fs83_database["cSolnPhaseTypeCS"]
    )


def test_expand_condition_species_accepts_formula_and_oxide_components(
    liu_alfesi_dat_database: dict,
    cas_mqm_fs73_database: dict,
) -> None:
    """expand_condition_species accepts formula and oxide components."""
    alfesi = expand_condition_species(
        liu_alfesi_dat_database,
        {"T": 1600.0, "P": 1.0, "Al13Fe4": 1.0},
    )
    cas = expand_condition_species(
        cas_mqm_fs73_database,
        {"T": 1600.0, "P": 1.0, "CaO": 1.0, "Al2O3": 1.0, "SiO2": 1.0},
    )

    assert alfesi["Al"] == pytest.approx(13.0)
    assert alfesi["Fe"] == pytest.approx(4.0)
    assert list(cas) == ["T", "P", "Ca", "O", "Al", "Si"]
    assert cas["O"] == pytest.approx(6.0)


def test_active_component_driving_force_uses_component_plane() -> None:
    """Active-component driving force uses the component plane."""
    basis = _formula_oxide_basis()
    coefficients = np.asarray([0.30, 0.20, 0.50])
    component_potentials = _oxide_component_potentials()
    candidate_gibbs = (
        float(coefficients @ np.asarray(list(component_potentials.values())))
        - 12_345.0
    )

    driving_force = basis.driving_force(
        candidate_gibbs,
        coefficients,
        component_potentials,
    )

    assert driving_force == pytest.approx(-12_345.0)


def test_active_component_rejects_invalid_basis_or_negative_coefficients() -> None:
    """Invalid basis or negative coefficients are rejected."""
    oxide_basis = _formula_oxide_basis()
    bad_oxide_projection = oxide_basis.project(
        {"Ca": 0.3, "Al": 0.4, "Si": 0.5, "O": 0.0}
    )
    signed_basis = ActiveComponentBasis.from_stoichiometries(
        {"AB": {"A": 1.0, "B": 1.0}, "A": {"A": 1.0}},
        ["A", "B"],
    )
    signed_projection = signed_basis.project(
        {"A": 0.0, "B": 1.0},
        require_nonnegative=True,
    )

    assert not bad_oxide_projection.compatible
    assert bad_oxide_projection.residual_norm > 1e-3
    assert not signed_projection.compatible
    assert signed_projection.nonnegative


def _formula_oxide_basis() -> ActiveComponentBasis:
    return ActiveComponentBasis.from_stoichiometries(
        {
            "CaO": parse_formula_stoichiometry("CaO"),
            "Al2O3": parse_formula_stoichiometry("Al2O3"),
            "SiO2": parse_formula_stoichiometry("SiO2"),
        },
        ["Ca", "Al", "Si", "O"],
    )


def _oxide_component_potentials() -> dict[str, float]:
    return {"CaO": -610_000.0, "Al2O3": -1_580_000.0, "SiO2": -850_000.0}
