
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

def _equilib_single(database,units,Condition,ListOfPhases=None):
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
    Condition: Dictionary containing NTP ensemble {'T':[val], 'P': [val], 'Al': [val], ...}
    
    Output
    Results dataclass
    ----------------------------------------------------------------------------------------------------------------'''
    
    read_dict(database)
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
    else:
        list_phases(var,NTPheader[2:])
    
    var.dConditionSys=condition
    input_condition(units,condition)
    minimize()  
    
    return None


def equilib_single(database,units,Condition,ListOfPhases=None):
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
    
    read_dict(database)
    # Get info from input condition
    NTPheader,NTPvals=_dict2np(Condition)
    
    NTPvals = np.squeeze(NTPvals)
    
    res = Result()
    
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
    else:
        list_phases(var,NTPheader[2:])
    
    var.dConditionSys=condition
    input_condition(units,condition)
    minimize()
    res.append_output()
    fort.resetthermo()        
    
    return res


