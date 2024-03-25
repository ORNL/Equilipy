
import numpy as np, os
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
from multiprocessing import Pool
from numba import njit


def _equilib_batch(Database:str,Unit:list,NTP:dict,ListOfPhases:list=None):
    
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
        
    for i in range(L):
        
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


def _batch_input(database,units,conditions,nPerBatch):
    """
    Break fullrange up into smaller sets of ranges that cover all
    the same numbers.
    """
    res = []
    header,NTP=_dict2np(conditions)
    n,m=NTP.shape
    fullrange=[0,n]

    for i in range(fullrange[0], fullrange[1], nPerBatch):
        subcondition=dict({})
        for j, head in enumerate(header):
            subcondition[head] = NTP[i:min(i+nPerBatch, fullrange[1]),j]
        res.append( [database,units,subcondition] ) 
    return res

def _equilib_mpi(arg,nCPU):
    pool = Pool(processes=nCPU)
    
    res = pool.starmap(_equilib_batch, tqdm(arg,colour='#8060ff',ascii="░▒▓",bar_format="{l_bar}{bar:50}{r_bar}"))
    
    pool.close()
    pool.join()
    return res

def equilib_batch(Database:str,Unit:list,NTP:dict,ListOfPhases:list=None,nCPU:int=os.cpu_count(),nPerBatch:int=1):
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
    # Get a batch input
    arg =_batch_input(Database,Unit,NTP,nPerBatch)
    res_mpi=list(_equilib_mpi(arg,nCPU))
    res=res_mpi[0]
    for i,r in enumerate(res_mpi):
        if i==0: continue
        res.append(r)
    return res