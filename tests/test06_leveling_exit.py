"""Compounds-only classical Leveling fast-exit certification regressions.

Author:
    S.Y. Kwon
Date:
    07/16/2026
Tests:
    1. test_solution_interior_falsifies_compounds_only_leveling_plane:
       A synthetic ideal solution lies below a compound plane even though
       both solution endmembers lie above it.  The fast-exit certificate must
       retain the solution witness and the production equilibrium must contain
       the solution phase.
    2. test_compounds_only_leveling_plane_records_standard_certificate:
       A high-energy ideal solution remains above the compound plane and records
       one fresh certified-settled sweep with the standard reason.
    3. test_system_without_solution_phases_keeps_trivial_fast_exit:
       The no-solution branch returns compounds-only without paying a sweep.
    4. test_public_equilibrium_uses_certified_compounds_only_fast_exit:
       The public equilibrium and PostProcess path records one honest certified
       sweep and returns the compound equilibrium with coherent terminal state.
"""

from __future__ import annotations

import hashlib
import tempfile
from pathlib import Path

import equilipy.equilifort as fort
import numpy as np

import equilipy as eq
from equilipy.equilib_single import _preprocess_single
from equilipy.minimize import (
    check_phase_assemblage,
    prepare_minimization,
    run_leveling,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
SYNTHETIC_TDB_TEXT = """$ ========================================================================
$ Synthetic compounds-only Leveling fast-exit falsifier
$
$ Author: S.Y. Kwon
$ Date:   07/16/2026
$
$ The two ideal-solution endmembers are 100 J/mol above the AL_LINE/CU_LINE
$ plane.  At 300 K, ideal mixing lowers the equimolar solution by RT*ln(2),
$ so the solution interior is stable even though classical endmember Leveling
$ sees only positive solution-endmember potentials.
$ ========================================================================

ELEMENT AL FCC_A1 26.9815384 0 0 !
ELEMENT CU FCC_A1 63.5460000 0 0 !

TYPE_DEFINITION % SEQ * !
DEFINE_SYSTEM_DEFAULT ELEMENT 2 !

PHASE IDEAL_SOLN % 1 1 !
CONSTITUENT IDEAL_SOLN :AL,CU: !
PARAMETER G(IDEAL_SOLN,AL;0) 1 100; 6000 N !
PARAMETER G(IDEAL_SOLN,CU;0) 1 100; 6000 N !

PHASE HIGH_SOLN % 1 1 !
CONSTITUENT HIGH_SOLN :AL,CU: !
PARAMETER G(HIGH_SOLN,AL;0) 1 3000; 6000 N !
PARAMETER G(HIGH_SOLN,CU;0) 1 3000; 6000 N !

PHASE AL_LINE % 1 1 !
CONSTITUENT AL_LINE :AL: !
PARAMETER G(AL_LINE,AL;0) 1 0; 6000 N !

PHASE CU_LINE % 1 1 !
CONSTITUENT CU_LINE :CU: !
PARAMETER G(CU_LINE,CU;0) 1 0; 6000 N !
"""
SYNTHETIC_TDB = Path(tempfile.gettempdir()) / "leveling_exit_ideal_solution.tdb"
SYNTHETIC_TDB.write_text(SYNTHETIC_TDB_TEXT)
SYNTHETIC_TDB_SHA256 = (
    "cfe9f9db732c29486911d86d9453736fa888f45ea0347bc6b9fafec140e7b798"
)
UNIT = ["K", "atm", "moles"]


def _fixture_database() -> dict:
    """Return the hash-pinned synthetic compounds/ideal-solution database."""
    actual_hash = hashlib.sha256(SYNTHETIC_TDB.read_bytes()).hexdigest()
    assert actual_hash == SYNTHETIC_TDB_SHA256
    return eq.read_tdb(SYNTHETIC_TDB)


def _prepare_leveling(
    database: dict,
    temperature_k: float,
    *,
    phases: list[str] | None = None,
) -> None:
    """Prepare and run classical Leveling for the equimolar synthetic binary."""
    _preprocess_single(
        database,
        {"T": temperature_k, "P": 1.0, "Al": 0.5, "Cu": 0.5},
        UNIT,
        phases=phases,
    )
    prepare_minimization()
    run_leveling()


def _exercise_provisional_fast_exit() -> tuple[int, float, float]:
    """Exercise the post-Leveling certificate boundary and return its state."""
    gem = fort.modulegemsolver
    gem.dgemfunctionnorm = 0.0
    gem.lcompbdonly = True
    sweep_count_before = int(gem.ngemcertificationsweep)

    fort.certifylevelingcompoundsonly()

    assert int(fort.modulethermoio.infothermo) == 0
    return (
        int(gem.ngemcertificationsweep) - sweep_count_before,
        float(gem.dpeaexitminphasepotential),
        float(gem.dpeaexittolerance),
    )


def test_solution_interior_falsifies_compounds_only_leveling_plane() -> None:
    """A below-plane solution interior must cancel the provisional fast exit."""
    database = _fixture_database()

    _prepare_leveling(database, 300.0)
    try:
        gem = fort.modulegemsolver
        thermo = fort.modulethermo
        n_species = int(thermo.nspecies)

        # Classical Leveling sees positive solution-endmember potentials and
        # selects only the two zero-energy line compounds.
        endmember_potentials = np.asarray(gem.dphasepotential, dtype=float)[:2]
        assert np.all(endmember_potentials > 0.0)
        assert np.all(np.asarray(thermo.iassemblage, dtype=int) > 0)

        sweep_delta, min_phase_potential, tolerance = _exercise_provisional_fast_exit()

        assert sweep_delta == 1
        assert not bool(gem.lcompbdonly)
        assert int(gem.ipeaexitfreshminpointsweep) == 1
        assert min_phase_potential < tolerance

        # The standard solution pseudo-compound rows retain the below-plane
        # witness for the normal PEA path.  Candidate identities are checked
        # numerically rather than through partition-sensitive display labels.
        solution_rows = np.asarray(gem.dphasepotential, dtype=float)[n_species:]
        assert np.min(solution_rows) == min_phase_potential

        check_phase_assemblage()
        assert np.any(np.asarray(thermo.iassemblage, dtype=int) < 0)
    finally:
        fort.resetthermo()

    result = eq.equilib_single(
        database,
        {"T": 300.0, "P": 1.0, "Al": 0.5, "Cu": 0.5},
        unit=UNIT,
    )
    assert set(result.stable_phases.names) == {"IDEAL_SOLN"}
    assert result.G < -1000.0


def test_compounds_only_leveling_plane_records_standard_certificate() -> None:
    """An above-plane solution minimum yields one honest fast-exit certificate."""
    database = _fixture_database()

    phases = ["HIGH_SOLN", "AL_LINE", "CU_LINE"]
    _prepare_leveling(database, 300.0, phases=phases)
    try:
        gem = fort.modulegemsolver
        sweep_delta, min_phase_potential, tolerance = _exercise_provisional_fast_exit()

        assert sweep_delta == 1
        assert bool(gem.lcompbdonly)
        assert int(gem.ipeaexitstatus) == int(gem.pea_exit_status_certified_settled)
        assert int(gem.ipeaexitreason) == int(gem.pea_exit_reason_compounds_only)
        assert int(gem.ipeaexitfreshminpointsweep) == 1
        assert min_phase_potential >= tolerance
        assert float(gem.dgemtimingcertification) >= 0.0
    finally:
        fort.resetthermo()

    result = eq.equilib_single(
        database,
        {"T": 300.0, "P": 1.0, "Al": 0.5, "Cu": 0.5},
        unit=UNIT,
        phases=phases,
    )
    assert set(result.stable_phases.names) == {"AL_LINE", "CU_LINE"}
    assert result.G == 0.0


def test_system_without_solution_phases_keeps_trivial_fast_exit() -> None:
    """A true compound-only system returns without a certificate sweep."""
    database = _fixture_database()
    phases = ["AL_LINE", "CU_LINE"]
    _preprocess_single(
        database,
        {"T": 300.0, "P": 1.0, "Al": 0.5, "Cu": 0.5},
        UNIT,
        phases=phases,
    )

    try:
        prepare_minimization()
        fort.gemsolver()
        gem = fort.modulegemsolver

        assert int(fort.modulethermo.nsolnphasessys) == 0
        assert bool(gem.lcompbdonly)
        assert int(gem.ngemcertificationsweep) == 0
        assert int(gem.ipeaexitstatus) == int(gem.pea_exit_status_unknown)
    finally:
        fort.resetthermo()


def test_public_equilibrium_uses_certified_compounds_only_fast_exit(
    monkeypatch,
) -> None:
    """The public minimize/PostProcess path must expose the certified fast exit."""
    database = _fixture_database()
    phases = ["HIGH_SOLN", "AL_LINE", "CU_LINE"]

    try:
        # equilib_single normally resets Fortran state after capturing its public
        # result.  Defer that reset so this test can assert the same call's exit
        # certificate and terminal fields after PostProcess has completed.
        with monkeypatch.context() as patch:
            patch.setattr(fort, "resetthermo", lambda: None)
            result = eq.equilib_single(
                database,
                {"T": 300.0, "P": 1.0, "Al": 0.5, "Cu": 0.5},
                unit=UNIT,
                phases=phases,
                include_heat_capacity=False,
            )

            gem = fort.modulegemsolver
            assert set(result.stable_phases.names) == {"AL_LINE", "CU_LINE"}
            assert result.G == 0.0
            assert int(gem.ngemcertificationsweep) == 1
            assert int(gem.ipeaexitstatus) == int(
                gem.pea_exit_status_certified_settled
            )
            assert int(gem.ipeaexitreason) == int(gem.pea_exit_reason_compounds_only)
            assert bool(gem.lconverged)
            assert not bool(gem.lphasechange)
            assert int(gem.iphasechangereason) == int(gem.phase_change_reason_none)
    finally:
        fort.resetthermo()
