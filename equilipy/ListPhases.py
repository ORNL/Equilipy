
from .utils import _pyvar2fvar
from .SystemCheck import system_check
from .ReadDict import read_dict
import equilipy.variables as var

def list_phases(database,elements):
    '''----------------------------------------------------------------------------------------------------------------
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
    ----------------------------------------------------------------------------------------------------------------'''
    read_dict(database)
    _pyvar2fvar(var)
    system_check(elements)
    
    # Store the dictionary for PhaseNames, types(soln or compd), and id
    PhaseTypeSys=['soln']*len(var.iSys2DBSoln)+['compd']*len(var.iSys2DBComp)
    iSys2DB= list(var.iSys2DBSoln)+list(var.iSys2DBComp)
    PhaseTypeIDSys=[[PhaseTypeSys[i],iSys2DB[i]] for i in range(len(iSys2DB))]
    var.PhaseNameSys=dict(zip(var.cPhaseNameSys,PhaseTypeIDSys))
    return var.cPhaseNameSys