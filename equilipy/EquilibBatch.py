
import numpy as np, os
from tqdm import tqdm
import equilipy.equilifort as fort
from .utils import _dict2np
from .InputCondition import input_condition
from .Minimize import minimize
from .PostProcess import Result, check_output_units
from .PhaseSelection import phase_selection
from .ListPhases import list_phases
from .Errors import *
from .ReadDict import read_dict
import equilipy.variables as var
# from mpi4py.futures import MPIPoolExecutor
# from mpi4py import MPI
from multiprocessing import Pool
from numba import njit


def _equilib_batch(Database:dict,NTP:dict,UnitIn:list=['K','atm','moles'],UnitOut:list=None,ListOfPhases:list=None):
    
    if UnitOut==None:
        UnitOut = UnitIn.copy() 
    
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
        input_condition(UnitIn,condition)
        
        try: 
            minimize()
            check_output_units(UnitOut)
            res.append_output()
            
        except EquilibError:
            res.append_error()
        fort.resetthermo() 
        fort.resetthermoparser()
        
    
    return res


def _batch_input(database,conditions,UnitIn,UnitOut,ListOfPhases,nPerBatch):
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
        res.append( [database,subcondition,UnitIn,UnitOut,ListOfPhases] ) 
    return res


def _equilib_singlenode(arg,nCPU):
    with Pool(processes=nCPU) as pool:
        res = pool.starmap(_equilib_batch, tqdm(arg,colour='#8060ff',ascii="░▒▓",bar_format="{l_bar}{bar:50}{r_bar}"))
    return res

# def _equilib_multinodes(arg,nCPU):
#     print(f'Process with mpi: {nCPU} processors')
#     with MPIPoolExecutor(max_workers=nCPU) as pool:
#         res = pool.starmap(_equilib_batch, arg)
#     return res

def equilib_batch(Database:dict,NTP:dict,UnitIn:list=['K','atm','moles'],UnitOut:list=None,ListOfPhases:list=None,nCPU:int=os.cpu_count(),nPerBatch:int=1):
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
    if UnitOut==None:
        UnitOut = UnitIn.copy() 
    # Get a batch input
    arg =_batch_input(Database,NTP,UnitIn,UnitOut,ListOfPhases,nPerBatch)
    
    # # Check if it is on a single node with multiple processors
    # if(MPI.COMM_WORLD.Get_size()==1):
    #     res_mpi=list(_equilib_singlenode(arg,nCPU))
    # else:
    #     res_mpi=list(_equilib_multinodes(arg,MPI.COMM_WORLD.Get_size()))
        
    res_mpi=list(_equilib_singlenode(arg,nCPU))
    
    #Concatenate all results
    res=res_mpi[0]
    for i,r in enumerate(res_mpi):
        if i==0: continue
        res.append(r)
        
    # Convert the datatype into numpy
    res.T = np.array(res.T)
    res.P = np.array(res.P)
    res.G = np.array(res.G)
    res.H = np.array(res.H)
    res.S = np.array(res.S)
    res.Cp = np.array(res.Cp)
    
    
    return res