
import numpy as np
from tqdm import tqdm
import equilipy.equilifort as fort
from .InternalFunctions import _dict2np
from .InputCondition import input_condition
from .Minimize import minimize
from .PostProcess import Result
from .PhaseSelection import phase_selection
from .ListPhases import list_phases
from .Errors import *
from .ReadDict import read_dict
import equilipy.variables as var


def equilib_batch(Database:str,Unit:list,NTP:dict,ListOfPhases:list=None):
    '''
    Calculate phase equlibria for multiple NTP conditions.
    ----------------------------------------------------------------------------------------------------------------

    Input
    -----
    Database : A string of 'Directory/databasename.dat'
    
    Unit     : A list of units for temperature, pressure, mass e.g.['K','atm','moles'].
               Temperature units, 'K'/'C'/'F'/'R'.
               Pressure units, 'atm'/'psi'/'bar'/'Pa'/'kPa'.
               Mass units, 'mass fraction'/'kilograms'/'grams'/'pounds'/'mole fraction'/'atom fraction'/'atoms'/'moles'.
               
    NTP      : Dictionary of T, P, element1, element2 ...

    Output
    ------
    Result: Equilipy result object that includes system properties and phase properties.
    
    System properties
    Result.N          : Dictionary {element: Amount,...}
    Result.T          : Temperature
    Result.P          : Pressure
    Result.G          : System Gibbs energy (total)
    Result.H          : System enthalpy (total)
    Result.S          : System entropy (total)
    Result.Cp         : System heat capacity (total)
    

    Phase properties
    Result.StablePhases: Dictionary {'ID':[phase_id1,..],'Name':[phase_name1,...],'Amount':[phase_amount1,...]'}
    Result.Phases      : Dictionary {phase_name: phase_object,...}
    
    '''
    # Read Database saved in a dictionary
    read_dict(Database)
    
    # Get info from NTP dictionary
    NTPheader,NTPvals=_dict2np(NTP)
    L,_=NTPvals.shape
    res = Result()
    
    
    try: 
        # If this is successfull, systems are defined previously
        if ListOfPhases!=None: phase_selection(ListOfPhases)
        else: KeyError('Error: add a list of phases')
    except AttributeError:
        list_phases(Database,NTPheader[2:])
        
    for i in tqdm(range(L),colour='#8060ff',ascii="░▒▓",bar_format="{l_bar}{bar:30}{r_bar}"):
        
        # Revise input condition when some of the elements are zero
        condition=NTPvals[i,:]
        comp=np.array(NTPvals[i,2:])
        
        if 0 in comp:
            elements=np.array(NTPheader[2:])
            elements=list(elements[comp!=0]) 
            condition=list(comp[comp!=0])
            condition=np.array(list(NTPvals[i,:2])+condition)
            
            list_phases(Database,elements)
            if ListOfPhases!=None: 
                phase_selection(ListOfPhases)
            
        else:
            list_phases(Database,NTPheader[2:])
            if ListOfPhases!=None: phase_selection(ListOfPhases)
        
        var.dConditionSys=condition
        input_condition(Unit,condition)
        try: 
            minimize()
            res.append_output()
            
        except EquilibError:
            res.append_error()
        
        fort.resetthermo() 
        fort.resetthermoparser()
        
    
    return res


