"""Public Scheil solidification smoke tests.

Author:
    S.Y. Kwon
Date:
    06/30/2026
Tests:
    1. test_alfesi_scheil_al25fe65si_matches_tcpython_markers: TC-Python marker
       check for Al-2.5Fe-6.5Si Scheil labels and liquidus behavior.
    2. test_alsi_scheil_terminal_transition_records_refined_hot_side:
       Terminal Al-Si transition check for refined hot-side bracket recording.
    3. test_alfesi_nucleoscheil_al25fe_does_not_reenter_liquid_only:
       NucleoScheil regression for liquid-only row re-entry after nucleation.
    4. test_nucleoscheil_records_failed_df_screen_warning:
       A failed batched driving-force screen remains an honest UNKNOWN.
    5. test_fs73_no_fcc_scheil_completes_from_found_liquidus:
       Restricted Example06 regression for liquidus startup and completion.
"""

import importlib
import math
from pathlib import Path

import pytest

import equilipy as eq
from equilipy.nucleoscheil_cooling import _init_nuclei

nucleoscheil_module = importlib.import_module("equilipy.nucleoscheil_cooling")
scheil_module = importlib.import_module("equilipy.scheil_cooling")

TCPYTHON_2025B_LIQUIDUS_C = 657.6265
HALLSTEDT_TDB = (
    Path(__file__).resolve().parents[1]
    / "tests"
    / "test_devop"
    / "fixtures"
    / "databases"
    / "2023Hallstedt_AlCoCrFeMnNiVC.tdb"
)


def test_hallstedt_316l_scheil_starts_near_liquidus() -> None:
    """The 316L path does not fall back to its 2000 C input temperature."""
    if not HALLSTEDT_TDB.is_file():
        pytest.skip("Hallstedt public database is unavailable")
    database = eq.read_tdb(HALLSTEDT_TDB)

    result = eq.scheil_cooling(
        "LIQUID",
        database,
        {"T": 2000.0, "P": 1.0, "Fe": 71.97, "Cr": 17.0, "Ni": 11.0, "C": 0.03},
        delta_T=5.0,
        unit=["C", "atm", "mass fraction"],
        start_from_liquidus=True,
        show_progress=False,
    )
    table = result.to_dict()

    assert table["T [C]"][0] == pytest.approx(1470.5, abs=2.0)
    assert table["T [C]"][0] < 1500.0
    assert table["label"][0] == "LIQUID"
    assert not any("LIQUIDUS_SEARCH_FAILED" in warning for warning in result.warnings)


def test_fs73_no_fcc_scheil_completes_from_found_liquidus(
    alcumgsi_fs73_database: dict,
) -> None:
    """Restricted Example06 completes without reverting to its input T."""
    condition = {
        "T": 900.0,
        "P": 1.0,
        "Al": 0.89260,
        "Cu": 0.01745,
        "Mg": 0.00114,
        "Si": 0.0881,
    }
    phases = [
        phase
        for phase in eq.list_phases(alcumgsi_fs73_database, list(condition)[2:])
        if "FCC_A1" not in phase
    ]

    result = eq.scheil_cooling(
        "LIQUID",
        alcumgsi_fs73_database,
        condition,
        delta_T=10.0,
        unit=["C", "atm", "g"],
        phases=phases,
        show_progress=False,
    )
    table = result.to_dict()

    assert table["T [C]"][0] == pytest.approx(1076.4, abs=2.0)
    assert table["T [C]"][0] > condition["T"]
    assert table["label"][0] == "LIQUID"
    assert table["label"][-1] == "DIAMOND_A4+GAMMA+HCP_A3+MG2SI(s)"
    assert table["fl_w"][-1] == pytest.approx(0.0)
    assert table["fs_w"][-1] == pytest.approx(1.0)
    assert result.warnings == []


def test_scheil_failed_liquidus_probe_records_typed_warning(
    liu_alfesi_dat_database: dict,
    monkeypatch,
) -> None:
    """An unsolvable liquidus probe is visible in the result and GUI report."""
    phases = [
        phase
        for phase in eq.list_phases(liu_alfesi_dat_database, ["Al", "Si"])
        if not phase.endswith("_SER(s)")
    ]

    def fail_liquidus(*args, **kwargs):
        raise eq.TransitionError("synthetic unsolvable liquidus probe")

    monkeypatch.setattr(scheil_module, "find_liquidus_transition", fail_liquidus)
    result = eq.scheil_cooling(
        "LIQUID",
        liu_alfesi_dat_database,
        {"T": 660.0, "P": 1.0, "Al": 94.0, "Si": 6.0},
        delta_T=100.0,
        unit=["C", "atm", "grams"],
        phases=phases,
        start_from_liquidus=True,
        show_progress=False,
    )

    assert result.to_dict()["T [C]"][0] == pytest.approx(660.0)
    assert any(
        warning.startswith("Warning [LIQUIDUS_SEARCH_FAILED]")
        and "synthetic unsolvable liquidus probe" in warning
        for warning in result.warnings
    )


def test_alfesi_scheil_al25fe65si_matches_tcpython_markers(
    liu_alfesi_dat_database: dict,
) -> None:
    """The Al-Fe-Si Scheil path matches TC-Python reference markers."""
    phases = [
        phase
        for phase in eq.list_phases(liu_alfesi_dat_database, ["Al", "Fe", "Si"])
        if not phase.endswith("_SER(s)")
    ]
    condition = {"T": 1600.0, "P": 1.0, "Al": 91.0, "Fe": 2.5, "Si": 6.5}

    result = eq.scheil_cooling(
        "LIQUID",
        liu_alfesi_dat_database,
        condition,
        delta_T=5.0,
        unit=["C", "atm", "grams"],
        phases=phases,
        start_from_liquidus=True,
        show_progress=False,
    )
    table = result.to_dict()
    labels = table["label"]

    assert table["T [C]"][0] == pytest.approx(
        TCPYTHON_2025B_LIQUIDUS_C,
        abs=20.0,
    )
    assert labels[0] == "LIQUID"
    assert labels[1].startswith("ALFESI_ALPHA+")
    assert labels[-1] == "ALFESI_BETA+DIAMOND_A4+FCC_A1"
    assert table["fl_w"][-1] == pytest.approx(0.0)
    assert table["fs_w"][-1] == pytest.approx(1.0)
    assert all(math.isfinite(value) for value in table["G [J]"])
    assert all(math.isfinite(value) for value in table["Cp [J/K]"])

    solid_start = next(index for index, label in enumerate(labels) if label != "LIQUID")
    assert all(label != "LIQUID" for label in labels[solid_start:])


def test_alsi_scheil_terminal_transition_records_refined_hot_side(
    liu_alfesi_dat_database: dict,
) -> None:
    """Check terminal eutectic transition rows keep find_transitions tolerance.

    :author: S.Y. Kwon
    :date: 06/30/2026
    :task:
        Reproduce the Al-6Si Scheil path where the terminal transition from
        ``FCC_A1+LIQUID`` to ``DIAMOND_A4+FCC_A1`` was bracketed internally but
        only the cold-side row was recorded. The row immediately before the
        terminal solid assemblage should be the refined hot-side bracket, not
        the previous coarse cooling step.
    """
    phases = [
        phase
        for phase in eq.list_phases(liu_alfesi_dat_database, ["Al", "Si"])
        if not phase.endswith("_SER(s)")
    ]

    result = eq.scheil_cooling(
        "LIQUID",
        liu_alfesi_dat_database,
        {"T": 1600.0, "P": 1.0, "Al": 94.0, "Si": 6.0},
        delta_T=5.0,
        unit=["C", "atm", "grams"],
        phases=phases,
        start_from_liquidus=True,
        show_progress=False,
    )
    table = result.to_dict()
    labels = table["label"]

    terminal_index = labels.index("DIAMOND_A4+FCC_A1")
    assert labels[terminal_index - 1] == "FCC_A1+LIQUID"
    assert (table["T [C]"][terminal_index - 1] - table["T [C]"][terminal_index]) < 0.1
    assert table["fl_w"][terminal_index] == pytest.approx(0.0)
    assert table["fs_w"][terminal_index] == pytest.approx(1.0)


def test_alfesi_nucleoscheil_al25fe_does_not_reenter_liquid_only(
    liu_alfesi_tdb_database: dict,
) -> None:
    """Check NucleoScheil keeps accepted nuclei through cooling.

    :author: S.Y. Kwon
    :date: 06/30/2026
    :task:
        Reproduce the Al-2.5Fe NucleoScheil condition where ``AL13FE4``
        nucleates with 13.5 K undercooling but the report temporarily returned
        pure ``LIQUID`` rows even though accumulated ``AL13FE4`` remained
        nonzero. Once solidification starts, later reported rows must not
        re-enter a liquid-only label.
    """
    result = eq.nucleoscheil_cooling(
        "LIQUID",
        liu_alfesi_tdb_database,
        {"T": 2000.0, "P": 1.0, "Al": 97.5, "Fe": 2.5},
        {"AL13FE4": 13.5, "FCC_A1": 0.5},
        delta_T=0.1,
        unit=["C", "atm", "wt%"],
        phases=["AL13FE4", "FCC_A1", "LIQUID"],
        show_progress=False,
    )
    table = result.to_dict()
    labels = table["label"]

    assert "AL13FE4+LIQUID" in labels
    solid_start = next(index for index, label in enumerate(labels) if label != "LIQUID")
    liquid_reentry = [
        (index, table["T [C]"][index])
        for index, label in enumerate(labels[solid_start:], start=solid_start)
        if label == "LIQUID"
    ]
    assert liquid_reentry == []


def test_nucleoscheil_records_failed_df_screen_warning(monkeypatch) -> None:
    """A typed failed DF screen cannot be interpreted as nucleation."""
    calls = []

    def fake_probe(liquid, database, condition, candidates, unit):
        calls.append(((liquid, *candidates), condition["T"]))
        raise eq.EquilibError(
            "forced batched-screen failure",
            reason=eq.SolverFailureReason.DF_SWEEP_UNKNOWN,
        )

    monkeypatch.setattr(
        nucleoscheil_module,
        "_batched_liquid_driving_force_values",
        fake_probe,
    )
    warnings = []
    with pytest.raises(eq.EquilibError, match="No phase"):
        _init_nuclei(
            "LIQUID",
            {},
            {"T": 2000.0, "P": 1.0, "Fe": 1.0},
            {"V3C2": 0.0, "FCC_A1": 0.0},
            warnings=warnings,
        )

    assert calls == [(("LIQUID", "V3C2", "FCC_A1"), 2000.0)]
    assert len(warnings) == 1
    assert "DF screen UNKNOWN" in warnings[0]
    assert "reason=df_sweep_unknown" in warnings[0]
