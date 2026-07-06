"""Phase listing helpers for parsed databases."""

from __future__ import annotations

import equilipy.variables as var

from .load_database import load_database
from .system_check import system_check


def list_phases(database, elements):
    """List phases available for the selected system elements.

    Description
    ===========
    This function lists all phases


    Revisions
    =========

     Date            Programmer      Description of change
     ----            ----------      ---------------------
     02/05/2024      S.Y. KWON       Original code


    Variables
    =========

    Input
    database : Global variable object that includes parsed database
    elements : A list of system elements e.g. ['Al','Fe','Si']

    Output
    list of phases
    """
    load_database(database)
    system_check(elements)

    # Store the dictionary for PhaseNames, types(soln or compd), and id
    PhaseTypeSys = ["soln"] * len(var.iSys2DBSoln) + ["compd"] * len(var.iSys2DBComp)
    iSys2DB = list(var.iSys2DBSoln) + list(var.iSys2DBComp)
    PhaseTypeIDSys = [[PhaseTypeSys[i], iSys2DB[i]] for i in range(len(iSys2DB))]
    var.PhaseNameSys = dict(zip(var.cPhaseNameSys, PhaseTypeIDSys, strict=False))
    hidden_helpers = {
        str(name).strip()
        for name in database.get("cOrderDisorderHelperPhaseNames", [])
    }
    return [phase for phase in var.cPhaseNameSys if phase not in hidden_helpers]
