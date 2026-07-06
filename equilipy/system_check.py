"""System composition validation helpers."""

from __future__ import annotations

import sys

import numpy as np

import equilipy.equilifort as fort
import equilipy.variables as var


def system_check(system_element_names):
    """Check validity of component names in an input table.

    Description
    ===========
    This checks the validity of component names given in a table. Errors occur
    when 1. a component has an invalid element name, and 2. a component is not
    part of the loaded database.


    Revisions
    =========

     Date            Programmer      Description of change
     ----            ----------      ---------------------
     11/19/2020      S.Y. KWON       Original code


    Variables
    =========

    Input
    NPT : pandas dataframe of T, P, element1, element2 ...

    Output
    res : an integer array of atomic number
    """
    # Info of name, mass, number for all elements
    # Consider using python class for storing the information

    Elements = var.cPeriodicTable
    n = len(system_element_names)
    var.cComponentNameSys = system_element_names
    var.iElementSys = np.ones(n)
    var.iElementSysIndex = np.ones(n, dtype=int)
    var.iElementDBIndex = np.ones(n, dtype=int)
    Sym = [x.strip() for x in system_element_names]

    iElementPass = np.ones(var.nElementsCS, dtype=int)
    for i in range(n):
        try:
            var.iElementSys[i] = int(Elements[Sym[i]][0])
        except KeyError:
            print("Error: Invalid element name")
        if "{:<3}".format(system_element_names[i]) not in var.cElementNameCS:
            print("Error:Some input elements are not in the database")
            sys.exit()
    j = 0

    for i in range(var.nElementsCS):
        if var.cElementNameCS[i].strip() in system_element_names:
            k = system_element_names.index(var.cElementNameCS[i].strip())
            var.iElementSysIndex[k] = j
            var.iElementDBIndex[k] = i
            iElementPass[i] = 0
            j += 1
    for i in getattr(var, "iPseudoComponentDependentElementDBIndex", []):
        if 0 <= int(i) < len(iElementPass):
            iElementPass[int(i)] = 0
    var.iElementSys = var.iElementSys[var.iElementSysIndex]

    lSpeciesPass = _get_endmembers(iElementPass)
    fort.moduleparsecs.lendmembers2species = True
    lCompPass = np.logical_and(var.iPhaseCS == 0, lSpeciesPass)

    # Get Species index and name for the system
    IndexSpecies = np.arange(len(lSpeciesPass))
    var.iSys2DBSpecies = IndexSpecies[lSpeciesPass]
    var.cSpeciesNameSys = np.array(var.cSpeciesNameCS)[var.iSys2DBSpecies]

    # Get Solution index and name for the system
    IndexSoln = var.iPhaseCS[lSpeciesPass]
    var.iSys2DBSoln = np.unique(IndexSoln[IndexSoln > 0])
    var.iSys2DBComp = IndexSpecies[lCompPass]

    var.cPhaseNameSys = np.array(var.cPhaseNames)[var.iSys2DBSoln - 1]
    var.cPhaseNameSys = np.append(
        var.cPhaseNameSys, np.array(var.cSpeciesNameCS)[lCompPass]
    )

    counts = {}
    phase_names = []
    for phase_name in var.cPhaseNameSys:
        name = str(phase_name).strip()
        counts[name] = counts.get(name, 0) + 1
        if counts[name] > 1:
            phase_names.append(f"{name}#{counts[name]}")
        else:
            phase_names.append(name)
    var.cPhaseNameSys = phase_names

    # Additonal variables for postprocess
    var.iSys2DBSolnDefault = var.iSys2DBSoln
    var.cPhaseNameSysDefault = var.cPhaseNameSys
    var.iSys2DBCompDefault = var.iSys2DBComp
    var.iSys2DBSpeciesDefault = var.iSys2DBSpecies

    return None


def _get_endmembers(element_pass_flags):
    # This function identifies which species to be considered for a given system

    # Get elements not considered in the system
    iElementNotPass = element_pass_flags.nonzero()[0]

    # First guess the species that need to be considered
    res = np.sum(var.dStoichSpeciesCS[:, iElementNotPass], axis=1) < 1e-15

    # Refine the species based on solution phases
    # Solution phase should have more than one endmember, otherwise remove the
    # endmembers.

    # Get solution index
    IndexSoln = var.iPhaseCS[res]
    n = len(res)
    k = 0

    for i in range(n):
        # Among the species that need to be passed in database
        if res[i]:
            if IndexSoln[k] > 0:
                # If more than two endmembers are involved add the phase
                res[i] = True
            elif IndexSoln[k] < 0:
                res[i] = False
            k += 1

    return res
