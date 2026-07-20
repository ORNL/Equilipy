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

_OD_TOPOLOGY_HELPER_STANDALONE = 1
_OD_TOPOLOGY_DIRECT_TARGET = 2
_OD_TOPOLOGY_HELPER_ONLY = 3
_OD_TOPOLOGY_AMBIGUOUS = 4
_OD_TOPOLOGY_LEGACY_NO_METADATA = 5


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
    (
        var.iOrderDisorderTopologyCS,
        var.iOrderDisorderStandalonePhaseCS,
    ) = _order_disorder_structure_metadata(db)

    return None


def _order_disorder_structure_metadata(
    database: dict,
) -> tuple[np.ndarray, np.ndarray]:
    """Classify DIS_PART identity from parser metadata and phase topology.

    DAT databases do not expose the TDB helper-role metadata and therefore
    receive the explicit legacy-fallback class.  For TDB databases, the
    DIS_PART target role and exact declared sublattice signatures distinguish
    helper+standalone, direct-target, helper-only, and ambiguous graphs.
    """
    phase_count = int(database.get("nSolnPhasesSysCS", 0))
    if (
        "iOrderDisorderTopologyCS" in database
        and "iOrderDisorderStandalonePhaseCS" in database
    ):
        topology = np.asarray(database["iOrderDisorderTopologyCS"], dtype=int)
        standalone = np.asarray(
            database["iOrderDisorderStandalonePhaseCS"],
            dtype=int,
        )
        if topology.shape != (phase_count,) or standalone.shape != (phase_count,):
            raise DatabaseLoadError(
                "Typed order/disorder identity metadata must have one entry "
                "per solution phase."
            )
        return topology.copy(), standalone.copy()

    topology = np.full(
        phase_count,
        _OD_TOPOLOGY_LEGACY_NO_METADATA,
        dtype=int,
    )
    standalone = np.zeros(phase_count, dtype=int)
    if "cOrderDisorderHelperPhaseNames" not in database:
        return topology, standalone

    topology.fill(0)
    names = [
        str(name).strip().upper()
        for name in database.get("cSolnPhaseNameCS", [])
    ]
    models = [
        str(model).strip().upper()
        for model in database.get("cSolnPhaseTypeCS", [])
    ]
    disordered = np.asarray(
        database.get("iDisorderedPhaseCS", np.zeros(phase_count)),
        dtype=int,
    )
    helper_names = {
        str(name).strip().upper()
        for name in database.get("cOrderDisorderHelperPhaseNames", [])
    }

    for ordered_index in range(min(phase_count, len(disordered))):
        target_id = int(disordered[ordered_index])
        if target_id <= 0:
            continue
        target_index = target_id - 1
        if target_index < 0 or target_index >= phase_count:
            topology[ordered_index] = _OD_TOPOLOGY_AMBIGUOUS
            continue
        if target_index >= len(names):
            topology[ordered_index] = _OD_TOPOLOGY_AMBIGUOUS
            continue

        if names[target_index] not in helper_names:
            topology[ordered_index] = _OD_TOPOLOGY_DIRECT_TARGET
            standalone[ordered_index] = target_id
            continue

        target_signature = _solution_topology_signature(database, target_index)
        matches = []
        for candidate_index in range(phase_count):
            if candidate_index in {ordered_index, target_index}:
                continue
            if candidate_index >= len(models) or target_index >= len(models):
                continue
            if models[candidate_index] != models[target_index]:
                continue
            if _solution_topology_signature(database, candidate_index) == target_signature:
                matches.append(candidate_index + 1)

        if len(matches) == 1:
            topology[ordered_index] = _OD_TOPOLOGY_HELPER_STANDALONE
            standalone[ordered_index] = matches[0]
        elif not matches:
            topology[ordered_index] = _OD_TOPOLOGY_HELPER_ONLY
        else:
            topology[ordered_index] = _OD_TOPOLOGY_AMBIGUOUS

    return topology, standalone


def _solution_topology_signature(database: dict, phase_index: int) -> tuple:
    """Return one phase's exact declared sublattice signature."""
    n_sublattice = int(np.asarray(database["nSublatticePhaseCS"])[phase_index])
    counts = np.asarray(database["nConstituentSublatticeCS"], dtype=int)
    ratios = np.asarray(database["dStoichSublatticeCS"], dtype=float)
    constituent_names = np.asarray(database["cConstituentNameSUBCS"])
    sublattices = []
    for sublattice in range(n_sublattice):
        count = int(counts[phase_index, sublattice])
        names = tuple(
            _fixed_width_name(
                constituent_names[phase_index, sublattice, constituent]
            )
            for constituent in range(count)
        )
        sublattices.append((float(ratios[phase_index, sublattice]), names))
    return tuple(sublattices)


def _fixed_width_name(value) -> str:
    """Return an uppercase name from one fixed-width parser character row."""
    array = np.asarray(value)
    if array.dtype.kind in {"S", "a"}:
        text = array.tobytes().decode(errors="ignore")
    else:
        text = "".join(str(item) for item in array.ravel())
    return text.strip().upper()


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
