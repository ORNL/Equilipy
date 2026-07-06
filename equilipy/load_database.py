"""Load parsed database dictionaries into module variables."""

from __future__ import annotations

import numpy as np

from . import variables as var
from .exceptions import DatabaseLoadError

_REQUIRED_DATABASE_KEYS = (
    "nElementsCS",
    "nSpeciesCS",
    "nSolnPhasesSysCS",
    "nCountSublatticeCS",
    "nMaxSpeciesPhaseCS",
    "nParamCS",
    "nMagParamCS",
    "nMaxSublatticeCS",
    "nSolnPhasesSysMax",
    "nGibbsCoeff",
    "nMaxGibbsEqs",
    "nParamMax",
    "iPhaseCS",
    "iParticlesPerMoleCS",
    "iPhaseSublatticeCS",
    "iDisorderedPhaseCS",
    "nGibbsEqSpecies",
    "nSublatticePhaseCS",
    "nParamPhaseCS",
    "nSpeciesPhaseCS",
    "nMagParamPhaseCS",
    "dAtomicMass",
    "cRegularParamCS",
    "cElementNameCS",
    "cSolnPhaseTypeCS",
    "cSolnPhaseNameCS",
    "cSpeciesNameCS",
    "nPairsSROCS",
    "nConstituentSublatticeCS",
    "nSublatticeElementsCS",
    "dGibbsMagneticCS",
    "dStoichSublatticeCS",
    "dZetaSpeciesCS",
    "dGibbsCoeffSpeciesTemp",
    "dStoichSpeciesCS",
    "iMagneticParamCS",
    "dMagneticParamCS",
    "iRegularParamCS",
    "dRegularParamCS",
    "cPairNameCS",
    "iConstituentSublatticeCS",
    "iPairIDCS",
    "iChemicalGroupCS",
    "dSublatticeChargeCS",
    "dStoichPairsCS",
    "dConstituentCoefficientsCS",
    "dCoordinationNumberCS",
    "cConstituentNameSUBCS",
    "cPhaseNames",
    "cEndmemberNameCS",
    "nPureSpeciesCS",
    "nSolnPhaseCS",
    "indx",
    "INFO",
    "DataBase",
    "iCounterGibbsEqn",
    "cPeriodicTable",
)

_FORTRAN_DATABASE_KEY = None


def load_database(db):
    """Load a parsed database dictionary into shared module variables."""
    if not isinstance(db, dict):
        raise DatabaseLoadError(
            f"Expected parsed database dictionary, got {type(db).__name__}."
        )

    missing = [key for key in _REQUIRED_DATABASE_KEYS if key not in db]
    if missing:
        missing_text = ", ".join(missing)
        raise DatabaseLoadError(f"Database is missing required key(s): {missing_text}")

    for key in _REQUIRED_DATABASE_KEYS:
        setattr(var, key, db[key])
    var.cOrderDisorderHelperPhaseNames = db.get(
        "cOrderDisorderHelperPhaseNames",
        [],
    )

    return None


def mark_fortran_database_stale() -> None:
    """Force the next Fortran database load to recopy all parser arrays."""
    global _FORTRAN_DATABASE_KEY
    _FORTRAN_DATABASE_KEY = None


def load_fortran_database(db, *, force: bool = False) -> None:
    """Load a parsed database into Python globals and, if needed, Fortran.

    Passing the full database dictionary into the Fortran module is expensive for
    large TDB-derived databases.  Cache the source object that was copied so
    repeated calculations can reset only the small per-calculation selection
    arrays instead of recopying thousands of species and parameters.
    """
    global _FORTRAN_DATABASE_KEY

    load_database(db)
    key = _database_cache_key(db)
    if force or _FORTRAN_DATABASE_KEY != key:
        import equilipy.equilifort as fort

        from .utils import _pyvar2fvar

        fort.resetthermo()
        _pyvar2fvar(var)
        _FORTRAN_DATABASE_KEY = key

    _reset_fortran_calculation_scope()


def _database_cache_key(db):
    """Return a cheap identity key for the parsed database payload."""
    return (
        id(db),
        id(db.get("iPhaseCS")),
        id(db.get("dStoichSpeciesCS")),
        id(db.get("dGibbsCoeffSpeciesTemp")),
        id(db.get("iRegularParamCS")),
        id(db.get("dRegularParamCS")),
        int(db.get("nSpeciesCS", -1)),
        int(db.get("nParamCS", -1)),
        int(db.get("nMagParamCS", -1)),
        int(db.get("nSolnPhasesSysCS", -1)),
    )


def _reset_fortran_calculation_scope() -> None:
    """Reset cheap Fortran arrays mutated by phase selection/checksystem."""
    import equilipy.equilifort as fort

    fort.moduleparsecs.info = 0
    fort.modulethermoio.infothermo = 0
    fort.moduleparsecs.iphasecs = var.iPhaseCS
    fort.modulethermo.isolnps = np.arange(
        1,
        int(var.nSolnPhasesSysCS) + 1,
        dtype=int,
    )
    fort.moduleparsecs.iparampasscs = np.zeros(max(int(var.nParamCS), 1), dtype=int)
    fort.moduleparsecs.imagparampasscs = np.zeros(
        max(int(var.nMagParamCS), 1),
        dtype=int,
    )
