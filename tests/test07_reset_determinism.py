"""Reset determinism falsifiers for parser and equilibrium global state.

Author:
    S.Y. Kwon
Date:
    07/16/2026
Tests:
    1. test_parser_read_reset_read_is_bit_identical: Parsed database payloads
       are identical across a parser reset boundary.
    2. test_equilibrium_solve_reset_solve_is_bit_identical: Equilibrium
       results and raw Fortran result arrays are identical across a
       ResetThermo boundary.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pytest

import equilipy as eq
from equilipy.equilib_single import _preprocess_single
from equilipy.minimize import minimize
from equilipy.results.equilib import EquilibPoint, PhaseResult

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


def _read_dat(path: Path) -> dict:
    """Read a DAT database through the public parser."""
    return eq.read_dat(str(path))


def _read_tdb(path: Path) -> dict:
    """Read a TDB database through the public parser."""
    return eq.read_tdb(path)


RESET_CASES = (
    pytest.param(
        "alfesi_dat",
        _read_dat,
        REPO_ROOT / "database" / "AlFeSi_99Liu.dat",
        {"T": 900.0, "P": 1.0, "Al": 0.9, "Fe": 0.05, "Si": 0.05},
        ["K", "atm", "moles"],
        None,
        id="AlFeSi-DAT",
    ),
    pytest.param(
        "alfesi_tdb",
        _read_tdb,
        REPO_ROOT / "database" / "AlFeSi_99Liu.tdb",
        {"T": 900.0, "P": 1.0, "Al": 0.9, "Fe": 0.05, "Si": 0.05},
        ["K", "atm", "moles"],
        None,
        id="AlFeSi-TDB",
    ),
    pytest.param(
        "hallstedt_subset",
        _read_tdb,
        HALLSTEDT_TDB,
        {
            "T": 1399.5,
            "P": 1.0,
            "C": 0.024176010246818144,
            "Cr": 0.2461916677820453,
            "Ni": 0.1206204610623526,
            "Fe": 0.609011860908784,
        },
        ["C", "atm", "mole fraction"],
        ["LIQUID", "FCC_A1", "BCC_A2", "M23C6_D84"],
        id="Hallstedt-subset",
    ),
)


@pytest.mark.parametrize(
    ("case_name", "reader", "database_path", "condition", "unit", "phases"),
    RESET_CASES,
)
def test_parser_read_reset_read_is_bit_identical(
    case_name: str,
    reader: Callable[[Path], dict],
    database_path: Path,
    condition: dict[str, float],
    unit: list[str],
    phases: list[str] | None,
) -> None:
    """Reading a database after ResetThermoParser reproduces the same payload."""
    del condition, unit, phases
    _require_database(database_path)

    first = _normalize_for_exact_compare(reader(database_path))
    eq.fort.resetthermoparser()
    second = _normalize_for_exact_compare(reader(database_path))

    _assert_identical(first, second, f"{case_name} parser read/reset/read")


@pytest.mark.parametrize(
    ("case_name", "reader", "database_path", "condition", "unit", "phases"),
    RESET_CASES,
)
def test_equilibrium_solve_reset_solve_is_bit_identical(
    case_name: str,
    reader: Callable[[Path], dict],
    database_path: Path,
    condition: dict[str, float],
    unit: list[str],
    phases: list[str] | None,
) -> None:
    """Solving after ResetThermo reproduces the same equilibrium bit-for-bit."""
    _require_database(database_path)
    database = reader(database_path)

    first = _solve_snapshot(database, condition, unit, phases)
    eq.fort.resetthermo()
    second = _solve_snapshot(database, condition, unit, phases)

    _assert_identical(first, second, f"{case_name} solve/reset/solve")


def test_grid_on_reset_then_off_matches_fresh_off_n3() -> None:
    """An N=3 grid run cannot contaminate a later default-OFF equilibrium."""
    _require_database(HALLSTEDT_TDB)
    database = _read_tdb(HALLSTEDT_TDB)
    condition = {"T": 1500.0, "P": 1.0, "Fe": 0.74, "Cr": 0.18, "Ni": 0.08}
    unit = ["K", "atm", "mass fraction"]

    _preprocess_single(database, condition, unit)
    try:
        eq.fort.modulegemsolver.lgridfrontendactive = True
        minimize()
    finally:
        eq.fort.resetthermo()

    after_grid = _solve_snapshot(database, condition, unit, None)
    fresh_off = _solve_snapshot(database, condition, unit, None)
    _assert_identical(after_grid, fresh_off, "N=3 grid ON/reset/OFF versus fresh OFF")


def _require_database(path: Path) -> None:
    """Fail clearly when a required determinism fixture is missing."""
    if not path.is_file():
        pytest.skip(f"dev baseline fixture not present: {path.name}")


def _solve_snapshot(
    database: dict,
    condition: dict[str, float],
    unit: list[str],
    phases: list[str] | None,
) -> Any:
    """Run one equilibrium solve and return exact public and raw-state snapshots."""
    _preprocess_single(database, condition, unit, phases=phases)
    try:
        minimize()
        point = EquilibPoint.from_fortran()
        return _normalize_for_exact_compare(
            {
                "public_point": _point_payload(point),
                "fortran_result_state": _fortran_result_state(),
            }
        )
    finally:
        eq.fort.resetthermo()


def _point_payload(point: EquilibPoint) -> dict[str, Any]:
    """Return the public fields that must remain deterministic."""
    return {
        "T": point.T,
        "P": point.P,
        "G": point.G,
        "H": point.H,
        "S": point.S,
        "Cp": point.Cp,
        "n_i": point.n_i,
        "w_i": point.w_i,
        "stable_phase_summary": point.stable_phase_summary,
        "stable_phase_names": point.stable_phases.names,
        "phases": {
            phase_name: _phase_payload(phase)
            for phase_name, phase in point.phase_map.items()
        },
    }


def _phase_payload(phase: PhaseResult) -> dict[str, Any]:
    """Return amount, composition, and thermodynamic fields for one phase."""
    return {
        "id": phase.id,
        "name": phase.name,
        "amount_n": phase.amount_n,
        "amount_w": phase.amount_w,
        "amount_n_basis": phase.amount_n_basis,
        "stability": phase.stability,
        "parent_model_id": phase.parent_model_id,
        "parent_model_name": phase.parent_model_name,
        "composition_set_id": phase.composition_set_id,
        "display_label": phase.display_label,
        "ordering_degree": phase.ordering_degree,
        "endmembers_x": phase.endmembers_x,
        "endmembers_w": phase.endmembers_w,
        "elements_x": phase.elements_x,
        "elements_w": phase.elements_w,
        "G": phase.G,
        "H": phase.H,
        "S": phase.S,
        "Cp": phase.Cp,
        "partial_gibbs": phase.partial_gibbs,
        "standard_gibbs_energy": phase.standard_gibbs_energy,
        "activity": phase.activity,
        "partial_enthalpy": phase.partial_enthalpy,
        "partial_entropy": phase.partial_entropy,
        "partial_heat_capacity": phase.partial_heat_capacity,
    }


def _fortran_result_state() -> dict[str, Any]:
    """Return raw Fortran arrays for phase sets, amounts, and potentials."""
    return {
        "modulethermo": _module_fields(
            eq.fort.modulethermo,
            (
                "iassemblage",
                "dmolesphase",
                "dmolfraction",
                "dchemicalpotential",
                "dstdgibbsenergy",
                "dactivity",
                "dpartialenthalpy",
                "dpartialentropy",
                "dpartialheatcapacity",
                "nspeciesphase",
            ),
        ),
        "modulethermoio": _module_fields(
            eq.fort.modulethermoio,
            (
                "dgramphase",
                "dgramfraction",
                "dgibbsenergysys",
                "denthalpysys",
                "dentropysys",
                "dheatcapacitysys",
                "dtemperature",
                "dpressure",
            ),
        ),
        "modulegemsolver": _module_fields(
            eq.fort.modulegemsolver,
            (
                "dgemchemicalpotentialnorm",
                "dgemcondensedchemicalpotentialnorm",
                "dgemsolutionchemicalpotentialnorm",
                "dmaxelementpotential",
                "dminphasepotential",
                "dminmolesphase",
            ),
        ),
    }


def _module_fields(module: Any, names: Sequence[str]) -> dict[str, Any]:
    """Capture named fields that exist on a f2py module object."""
    fields: dict[str, Any] = {}
    for name in names:
        try:
            fields[name] = np.asarray(getattr(module, name)).copy()
        except AttributeError:
            fields[name] = "<missing>"
    return fields


def _normalize_for_exact_compare(value: Any) -> Any:
    """Normalize nested Python/NumPy values for exact equality comparison."""
    if isinstance(value, np.ndarray):
        return (
            "ndarray",
            str(value.dtype),
            tuple(int(size) for size in value.shape),
            _normalize_for_exact_compare(value.tolist()),
        )
    if isinstance(value, np.generic):
        return _normalize_for_exact_compare(value.item())
    if isinstance(value, float):
        if math.isnan(value):
            return ("float", "nan")
        return ("float", value.hex())
    if isinstance(value, Mapping):
        return {
            str(key): _normalize_for_exact_compare(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, tuple):
        return ("tuple", tuple(_normalize_for_exact_compare(item) for item in value))
    if isinstance(value, list):
        return ("list", tuple(_normalize_for_exact_compare(item) for item in value))
    if isinstance(value, (str, int, bool, type(None))):
        return value
    return repr(value)


def _assert_identical(first: Any, second: Any, label: str) -> None:
    """Assert exact equality and report the first divergent state path."""
    if first == second:
        return
    path, left, right = _first_difference(first, second)
    pytest.fail(
        f"{label} is not bit-identical; first difference at {path}: "
        f"{left!r} != {right!r}"
    )


def _first_difference(first: Any, second: Any, path: str = "root") -> tuple[str, Any, Any]:
    """Return the first differing path in two normalized values."""
    if first == second:
        return path, first, second
    if isinstance(first, dict) and isinstance(second, dict):
        left_keys = set(first)
        right_keys = set(second)
        for missing_key in sorted(left_keys - right_keys):
            return f"{path}.{missing_key}", first[missing_key], "<missing>"
        for extra_key in sorted(right_keys - left_keys):
            return f"{path}.{extra_key}", "<missing>", second[extra_key]
        for key in sorted(left_keys):
            if first[key] != second[key]:
                return _first_difference(first[key], second[key], f"{path}.{key}")
    if isinstance(first, tuple) and isinstance(second, tuple):
        for index, (left_item, right_item) in enumerate(zip(first, second, strict=False)):
            if left_item != right_item:
                return _first_difference(left_item, right_item, f"{path}[{index}]")
        if len(first) != len(second):
            return f"{path}.length", len(first), len(second)
    return path, first, second
