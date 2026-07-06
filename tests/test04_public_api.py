"""Public helper and result-table API smoke tests.

Author:
    S.Y. Kwon
Date:
    06/30/2026
Tests:
    1. test_gibbs_cp_custom_ln_t_roundtrips_as_custom_power_99: G/H/S/Cp helper
       roundtrip for the custom ln(T) power term.
    2. test_hscp2g_keeps_tiny_high_power_custom_cp_coefficients: G/H/S/Cp
       helper keeps tiny coefficients with thermodynamically significant powers.
    3. test_result_table_preserves_column_order_and_broadcasts_scalars:
       ResultTable column ordering and scalar broadcasting.
    4. test_result_table_converts_to_polars: ResultTable Polars conversion.
    5. test_debug_mode_toggle_updates_fortran_flag: public debug switch updates
       the Fortran minimizer flag.
"""

import numpy as np
import pytest

import equilipy as eq
from equilipy.results import ResultTable


def test_gibbs_cp_custom_ln_t_roundtrips_as_custom_power_99() -> None:
    """Custom ln(T) Cp terms roundtrip as custom power 99."""
    gibbs_energies = np.zeros((1, 16))
    gibbs_energies[0, 0] = 298.15
    gibbs_energies[0, 1] = 6000.0
    gibbs_energies[0, 8] = 5.0
    gibbs_energies[0, 9] = 99.0

    h298, s298, cp = eq.G2HSCp(gibbs_energies=gibbs_energies)
    reconstructed = eq.HSCp2G(h298, s298, cp)

    assert cp[0, 6] == pytest.approx(5.0)
    assert cp[0, 7] == pytest.approx(99.0)
    assert reconstructed[0, 8] == pytest.approx(5.0)
    assert reconstructed[0, 9] == pytest.approx(99.0)


def test_hscp2g_keeps_tiny_high_power_custom_cp_coefficients() -> None:
    """Check tiny custom Cp coefficients are not dropped by raw magnitude."""
    coefficient = 7.9337e-20
    power = 6.0
    t_ref = 298.15
    t_break = 933.4269272672689
    cp = np.array(
        [
            [
                t_ref,
                t_break,
                0.0,
                0.0,
                0.0,
                0.0,
                coefficient,
                power,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ],
            [
                t_break,
                1200.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ],
        ],
    )

    reconstructed = eq.HSCp2G(0.0, 0.0, cp)
    expected_h_offset = coefficient / (power + 1.0) * (
        t_break ** (power + 1.0) - t_ref ** (power + 1.0)
    )
    expected_s_offset = coefficient / power * (t_break**power - t_ref**power)

    assert reconstructed[0, 8] == pytest.approx(
        coefficient / (power + 1.0) - coefficient / power,
    )
    assert reconstructed[0, 9] == pytest.approx(power + 1.0)
    assert reconstructed[1, 2] == pytest.approx(expected_h_offset)
    assert reconstructed[1, 3] == pytest.approx(-expected_s_offset)


def test_result_table_preserves_column_order_and_broadcasts_scalars() -> None:
    """ResultTable preserves column order and broadcasts scalars."""
    table = ResultTable.from_dict(
        {
            "T [K]": [1000.0, 1100.0],
            "P [atm]": 1.0,
            "Label": np.array(["A", "B"]),
        }
    )

    assert table.available_columns() == ["T [K]", "P [atm]", "Label"]
    assert table.to_rows() == [
        {"T [K]": 1000.0, "P [atm]": 1.0, "Label": "A"},
        {"T [K]": 1100.0, "P [atm]": 1.0, "Label": "B"},
    ]


def test_result_table_converts_to_polars() -> None:
    """ResultTable converts to a polars DataFrame."""
    table = ResultTable.from_dict({"T [K]": [1000.0], "G [J]": [-1.0]})

    frame = table.to_polars()

    assert frame.columns == ["T [K]", "G [J]"]
    assert frame.shape == (1, 2)


def test_debug_mode_toggle_updates_fortran_flag() -> None:
    """The debug-mode toggle updates the Fortran flag."""
    try:
        eq.set_debug_mode(True)
        assert eq.debug_mode_enabled()
        assert bool(eq._f2py_extension.modulegemsolver.ldebugmode)
    finally:
        eq.set_debug_mode(False)

    assert not eq.debug_mode_enabled()
    assert not bool(eq._f2py_extension.modulegemsolver.ldebugmode)
