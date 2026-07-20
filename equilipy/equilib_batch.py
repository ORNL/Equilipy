"""Batch equilibrium calculation helpers."""

from __future__ import annotations

import os

import numpy as np

import equilipy.equilifort as fort
import equilipy.variables as var

from ._parallel import starmap_joblib
from .composition import expand_condition_species
from .equilib_single import _active_phase_selection, _validated_phase_selection
from .exceptions import EquilibError, InputConditionError
from .input_condition import input_condition
from .list_phases import list_phases
from .load_database import load_fortran_database
from .minimize import minimize
from .phase_selection import phase_selection
from .results import capture_result_context
from .results.equilib import EquilibPoint, EquilibResult
from .utils import _dict2np

# from mpi4py.futures import MPIPoolExecutor
# from mpi4py import MPI


def _equilib_batch(
    database: dict,
    condition: dict,
    unit: list = None,
    phases: list = None,
):
    if unit is None:
        unit = ["K", "atm", "moles"]

    # Get info from condition dictionary.
    condition = expand_condition_species(database, condition, unit[2])
    NPTheader, NPTvals = _dict2np(condition)
    L, _ = NPTvals.shape
    res = EquilibResult(
        context=capture_result_context(
            database,
            condition,
            unit,
            phases,
        )
    )
    points = []
    phase_cache: dict[tuple[str, ...], list[str]] = {}

    def phases_for(elements) -> list[str]:
        key = tuple(str(element) for element in elements)
        if key not in phase_cache:
            phase_cache[key] = list(list_phases(database, list(key)))
        return phase_cache[key]

    for i in range(L):
        load_fortran_database(database)
        row_condition = NPTvals[i, :]
        comp = np.array(row_condition[2:])
        phase_elements = list(NPTheader[2:])
        phase_selector = _validated_phase_selection
        if np.any(np.isclose(comp, 0.0)):
            elements = np.array(NPTheader[2:])
            active_mask = ~np.isclose(comp, 0.0)
            phase_elements = list(elements[active_mask])
            if not phase_elements:
                raise InputConditionError(
                    "At least one component amount must be non-zero."
                )
            row_condition = list(comp[active_mask])
            row_condition = np.array(list(NPTvals[i, :2]) + row_condition)
            phase_selector = _active_phase_selection

        all_phases = phases_for(phase_elements)
        if phases is not None:
            phase_selection(phase_selector(phases, all_phases))

        var.dConditionSys = row_condition
        input_condition(unit, row_condition)

        try:
            minimize()
            points.append(EquilibPoint.from_fortran())

        except EquilibError as error:
            warning_result = EquilibResult(context=res.context)
            if warning_result.append_postprocess_warning(error):
                points.append(warning_result.point)
            else:
                points.append(EquilibPoint.for_error(error))

        fort.resetthermo()

    if points:
        res.data = points[0] if len(points) == 1 else points

    return res


def _equilib_batch_input(database, condition, unit, phases, n_per_batch):
    """
    Split conditions into smaller batch dictionaries.

    The batches cover the same full input range.
    """
    res = []
    header, NPT = _dict2np(condition)
    n, m = NPT.shape
    fullrange = [0, n]

    for i in range(fullrange[0], fullrange[1], n_per_batch):
        subcondition = dict({})
        for j, head in enumerate(header):
            subcondition[head] = NPT[i : min(i + n_per_batch, fullrange[1]), j]
        res.append([database, subcondition, unit, phases])
    return res


def _equilib_singlenode(arg, n_cpu, progress_callback=None, progress=True):
    return starmap_joblib(
        _equilib_batch,
        arg,
        n_cpu,
        progress=progress,
        progress_callback=progress_callback,
        colour="#8060ff",
        ascii="░▒▓",
        bar_format="{l_bar}{bar:50}{r_bar}",
        desc="Equilibrium Batch",
    )


def equilib_batch(
    database: dict,
    condition: dict,
    unit: list = None,
    phases: list = None,
    n_cpu: int | None = None,
    n_per_batch: int = 1,
    progress_callback=None,
    progress: bool = True,
):
    """
    Calculate phase equlibria for multiple NPT conditions.

    Input
    -----
    database : A string of 'Directory/databasename.dat'

    unit     : A list of units for temperature, pressure, mass e.g.['K','atm','moles'].
               Temperature units, 'K'/'C'/'F'/'R'.
               Pressure units, 'atm'/'psi'/'bar'/'Pa'/'kPa'.
               Mass units, 'mass fraction'/'kilograms'/'grams'/'pounds',
               'mole fraction'/'atom fraction'/'atoms'/'moles'.

    condition: Dictionary of T, P, element1, element2 ...

    progress : Set False to suppress the tqdm progress bar, e.g. for
               non-root MPI ranks or non-interactive logs.

    Output
    ------
    EquilibResult: Equilipy result object that includes system and phase
        properties.

    System properties
    EquilibResult.n_i        : Dictionary {element: amount on mole basis, ...}
    EquilibResult.w_i        : Dictionary {element: amount on mass basis, ...}
    EquilibResult.T          : Temperature
    EquilibResult.P          : Pressure
    EquilibResult.G          : System Gibbs energy (total)
    EquilibResult.H          : System enthalpy (total)
    EquilibResult.S          : System entropy (total)
    EquilibResult.Cp         : System heat capacity (total)


    Phase properties
    EquilibResult.stable_phases : PhaseCollection of stable phases for a
                                  one-point result.
    EquilibResult.phases        : PhaseCollection of all phases for a
                                  one-point result.
    EquilibResult.points[i]     : EquilibPoint access for batch results.

    """
    if unit is None:
        unit = ["K", "atm", "moles"]
    if n_cpu is None:
        n_cpu = os.cpu_count() or 1
    expanded_condition = expand_condition_species(database, condition, unit[2])

    # Get a batch input
    arg = _equilib_batch_input(database, expanded_condition, unit, phases, n_per_batch)

    res_mpi = list(_equilib_singlenode(arg, n_cpu, progress_callback, progress))

    # Concatenate all worker points and materialize batch views once.
    points = []
    for result in res_mpi:
        points.extend(result.points)
    res = EquilibResult()
    if points:
        res.data = points[0] if len(points) == 1 else points
    res.context = capture_result_context(
        database,
        expanded_condition,
        unit,
        phases,
    )

    return res
