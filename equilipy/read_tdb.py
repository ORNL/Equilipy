"""Read Thermo-Calc TDB database files."""

from __future__ import annotations

from collections.abc import Collection
from pathlib import Path

from .database_ir import read_tdb as _read_database_ir_tdb
from .database_ir.runtime import to_chemsage_compounds


def read_tdb(
    file_name: str | Path,
    *,
    strict: bool = True,
    auto_correct: bool = False,
    editable: bool = False,
    as_ir: bool = False,
    remove_redundant_phases: bool = False,
    include_stoichiometric_compounds: bool = True,
    include_endmember_compounds: bool = False,
    endmember_compound_phases: Collection[str] | str | None = None,
):
    """Read a TDB file.

    By default this returns the same ChemSage-style calculation dictionary shape
    as :func:`read_dat`, populated from executable ``PHASE`` records for
    stoichiometric compounds, Bragg-Williams/RKMP solutions, and CEF/SUBL
    solutions. Pass ``editable=True`` to get the editable ``DatabaseIR`` object
    used by the GUI and TDB writer. ``as_ir=True`` is kept as a compatibility
    alias for the same behavior. Known disordered helper aliases such as
    ``A1_FCC`` and ``A2_BCC`` are preserved for calculation by default because
    ordered phases may need their DIS_PART reference parameters even when a
    canonical physical phase is also present. Pass ``remove_redundant_phases=True``
    only when intentionally inspecting the older canonicalized import behavior.
    """
    if not isinstance(file_name, (str, Path)):
        raise AssertionError("Error: File name must be a string or Path")
    database = _read_database_ir_tdb(
        file_name,
        strict=strict,
        auto_correct=auto_correct,
        remove_redundant_phases=remove_redundant_phases,
    )
    if editable or as_ir:
        return database
    return to_chemsage_compounds(
        database,
        include_stoichiometric_compounds=include_stoichiometric_compounds,
        include_endmember_compounds=include_endmember_compounds,
        endmember_compound_phases=endmember_compound_phases,
    )
