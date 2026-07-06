"""Public Equilipy API."""

from __future__ import annotations

import equilipy.equilifort as fort
import equilipy.variables as variables

from .equilib_batch import equilib_batch
from .equilib_cooling import equilib_cooling
from .equilib_single import equilib_single
from .exceptions import (
    DatabaseLoadError,
    DatabaseParsingError,
    EquilibError,
    InputConditionError,
    PostProcessError,
    TransitionError,
)
from .find_transition import find_first_transition, find_transitions
from .input_condition import input_condition
from .list_phases import list_phases
from .load_database import load_database
from .minimize import (
    check_phase_assemblage,
    debug_mode_enabled,
    minimize,
    prepare_minimization,
    run_lagrangian_gem,
    run_leveling,
    set_debug_mode,
)
from .nucleoscheil_cooling import (
    nucleoscheil_batch,
    nucleoscheil_constituent_batch,
    nucleoscheil_cooling,
)
from .phase_selection import phase_selection
from .property_compound import property_compound
from .property_solution import property_solution
from .read_dat import read_dat
from .read_tdb import read_tdb
from .result_io import load_result, save_result
from .results import (
    CompositionFractions,
    PhaseCollection,
    ResultColumn,
    ResultContext,
    ResultTable,
    SlnPropertyPhase,
    SlnPropertyPoint,
    SlnPropertyResult,
    SpeciesProperty,
    capture_result_context,
)
from .results.common import combine_dicts
from .results.equilib import (
    EquilibPoint,
    EquilibResult,
    PhaseResult,
    create_phase_from_sys,
    endmembers2elements,
)
from .results.scheil import (
    ScheilConstituentSegment,
    ScheilPoint,
    ScheilResult,
)
from .scheil_cooling import scheil_batch, scheil_constituent_batch, scheil_cooling
from .simplex import simplex_count, simplex_grid, simplex_grid_shift
from .split_tdb import split_tdb
from .system_check import system_check
from .utils import G2HSCp, G2HSCp_single, HSCp2G, NeumanKoppHSCp
from .variables import (
    cPeriodicTable,
    nGibbsCoeff,
    nMaxGibbsEqs,
    nMaxSublatticeCS,
    nParamMax,
    nSolnPhasesSysMax,
    to_dict,
)

_f2py_extension = fort
_global_state_module = variables
_internal_database_loader = load_database

__all__ = [
    "DatabaseParsingError",
    "DatabaseLoadError",
    "EquilibError",
    "EquilibPoint",
    "EquilibResult",
    "CompositionFractions",
    "G2HSCp",
    "G2HSCp_single",
    "HSCp2G",
    "InputConditionError",
    "NeumanKoppHSCp",
    "PostProcessError",
    "PhaseCollection",
    "PhaseResult",
    "TransitionError",
    "ResultColumn",
    "ResultContext",
    "ResultTable",
    "ScheilConstituentSegment",
    "ScheilPoint",
    "ScheilResult",
    "SlnPropertyPhase",
    "SlnPropertyPoint",
    "SlnPropertyResult",
    "SpeciesProperty",
    "cPeriodicTable",
    "capture_result_context",
    "check_phase_assemblage",
    "combine_dicts",
    "create_phase_from_sys",
    "debug_mode_enabled",
    "endmembers2elements",
    "equilib_batch",
    "equilib_cooling",
    "equilib_single",
    "find_first_transition",
    "find_transitions",
    "input_condition",
    "list_phases",
    "load_result",
    "minimize",
    "nGibbsCoeff",
    "nMaxGibbsEqs",
    "nMaxSublatticeCS",
    "nParamMax",
    "nSolnPhasesSysMax",
    "nucleoscheil_batch",
    "nucleoscheil_constituent_batch",
    "nucleoscheil_cooling",
    "phase_selection",
    "prepare_minimization",
    "property_compound",
    "property_solution",
    "read_dat",
    "read_tdb",
    "run_lagrangian_gem",
    "run_leveling",
    "scheil_batch",
    "scheil_constituent_batch",
    "scheil_cooling",
    "set_debug_mode",
    "save_result",
    "simplex_count",
    "simplex_grid",
    "simplex_grid_shift",
    "split_tdb",
    "system_check",
    "to_dict",
]
