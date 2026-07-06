"""Split a TDB database onto an element subset."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

from .database_ir import DatabaseIR
from .database_ir import read_tdb as _read_tdb_ir
from .database_ir import write_tdb as _write_tdb
from .database_ir.ops import split_database


def split_tdb(
    database: DatabaseIR | str | Path,
    elements: Sequence[str],
    out: str | Path | None = None,
) -> DatabaseIR:
    """Split a TDB database onto an element subset.

    Input
    -----
    database : Path to a ``.tdb`` file, or an editable ``DatabaseIR`` from
               ``eq.read_tdb(path, editable=True)``.
    elements : Element symbols to keep, e.g. ``["Al", "Fe"]``.  ``VA`` and
               the electron pseudo-element are kept automatically.
    out      : Optional path; when given, the split database is also written
               there as a TDB file.

    Output
    ------
    The split ``DatabaseIR``.  For equilibrium calculations, read the
    written file back with ``eq.read_tdb(out)``.

    Example
    -------
    >>> import equilipy as eq
    >>> eq.split_tdb("2023Hallstedt_AlCoCrFeMnNiVC.tdb", ["Al", "Fe"],
    ...              out="AlFe_Hallstedt.tdb")
    """
    if isinstance(database, (str, Path)):
        database = _read_tdb_ir(database)
    split_ir = split_database(database, elements)
    if out is not None:
        _write_tdb(split_ir, out)
    return split_ir
