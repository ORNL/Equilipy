
import numpy as np
import equilipy.equilifort as fort
from .PhaseSelection import phase_selection
from .InternalFunctions import _dict2np
from .InputCondition import input_condition
from .Minimize import minimize
from .PostProcess import Result
from .ListPhases import list_phases
import equilipy.variables as var
from .ReadDict import read_dict
from .Errors import *

def _preprocess_single(Database,Condition,Unit:list=['K','atm','moles'],ListOfPhases=None):
    read_dict(Database)
    # Get info from input condition
    NTPheader,NTPvals=_dict2np(Condition)
    
    NTPvals = np.squeeze(NTPvals)
    
    try: 
        x=var.PhaseNameSys
        # If this is successfull, systems are defined previously
        if ListOfPhases!=None: phase_selection(ListOfPhases)
        else: KeyError('Error: add a list of phases')
    except AttributeError:
        list_phases(var,NTPheader[2:])
    
    condition=NTPvals
    comp=np.array(NTPvals[2:])
    
    
    if 0 in comp:
        elements=np.array(NTPheader[2:])
        elements=list(elements[comp!=0]) 
        condition=list(comp[comp!=0])
        condition=np.array(list(NTPvals[:2])+condition)
        list_phases(var,elements)
        if ListOfPhases!=None: 
                phase_selection(ListOfPhases)
    else:
        list_phases(var,NTPheader[2:])
        if ListOfPhases!=None: phase_selection(ListOfPhases)
    
    var.dConditionSys=condition
    input_condition(Unit,condition)        
    
    return None

def _equilib_single(Database,Condition,Unit:list=['K','atm','moles'],ListOfPhases=None):
    _preprocess_single(Database,Condition,Unit,ListOfPhases=ListOfPhases)
    
    try: minimize()  
    except EquilibError: pass
        
    
    return None

def equilib_single(Database,Condition,Unit:list=['K','atm','moles'],ListOfPhases=None):
    '''----------------------------------------------------------------------------------------------------------------
    Description
    ===========
    This function conducts equilibrium calculations for one condition.


    Revisions
    =========

     Date            Programmer      Description of change
     ----            ----------      ---------------------
     12/29/2023      S.Y. KWON       Original code


    Variables
    =========

    Input
    datafile : A string of 'Directory/databasename.dat'
    units    : A list of units for temperature, pressure, mass e.g.['K','atm','moles'].
               Temperature units, 'K'/'C'/'F'/'R'.
               Pressure units, 'atm'/'psi'/'bar'/'Pa'/'kPa'.
               Mass units, 'mass fraction'/'kilograms'/'grams'/'pounds'/'mole fraction'/'atom fraction'/'atoms'/'moles'.
    Components: list of strings containing element name[element1, element2 ...]
    NTP       : [T, P, element1 amount, element1 amount ...]

    Output
    Results dataclass
    ----------------------------------------------------------------------------------------------------------------'''
    res = Result()
    _preprocess_single(Database,Condition,Unit,ListOfPhases)
    try: 
        minimize()  
        res.append_output()
    except EquilibError:
        res.append_error()
    
    fort.resetthermoall() 
    
    return res


