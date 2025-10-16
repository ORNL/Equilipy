
import numpy as np
from tqdm import tqdm
import equilipy.equilifort as fort
from .utils import _dict2np
from .InputCondition import input_condition
from .Minimize import *
from .PhaseSelection import phase_selection
from .ListPhases import list_phases
from .Errors import *
from .ReadDict import read_dict
import equilipy.variables as var


def single_phase_property(Database:str,Unit:list,NTP:dict,PhaseName:str):
    '''
    Calculate thermodynamic properties for multiple NTP conditions.
    ----------------------------------------------------------------------------------------------------------------

    Input
    -----
    Database : A string of 'Directory/databasename.dat'
    
    Unit     : A list of units for temperature, pressure, mass e.g.['K','atm','moles'].
               Temperature units, 'K'/'C'/'F'/'R'.
               Pressure units, 'atm'/'psi'/'bar'/'Pa'/'kPa'.
               Mass units, 'mass fraction'/'kilograms'/'grams'/'pounds'/'mole fraction'/'atom fraction'/'atoms'/'moles'.
               
    NTP      : Dictionary of T, P, element1, element2 ...
    
    PhaseName: A string of phase name

    Output
    ------
    Result: Equilipy object for phase properties.    

    Phase properties
    Result.T      : list of temperature conditions
    Result.P      : list of pressure conditions
    Result.N      : list of composition conditions
    Result.G      : list of Gibbs energies of the phase
    Result.H      : list of enthalpies of the phase
    Result.S      : list of entropies of the phase
    Result.Cp     : list of heat capacities of the phase
    Result.gi     : dictionary of partial molar Gibbs energy of all endmembers 
    Result.hi     : dictionary of partial molar enthalpy of all endmembers
    Result.si     : dictionary of partial molar entropy of all endmembers
    Result.ai     : dictionary of partial molar activities of all endmembers
    Result.Hess   : Embeded dictionary of Hessian (dlnr_i/dlnx_j)
    '''
    read_dict(Database)
    
    # Get info from NTP dictionary
    NTPheader,NTPvals=_dict2np(NTP)
    L,_=NTPvals.shape
    
    # Check if the phase name is in the database
    AllPhases=list_phases(var,NTPheader[2:])
    if PhaseName not in AllPhases: raise EquilibError(f'Database does not include {PhaseName}')
    else: phase_selection([PhaseName])
    
   
    for i in tqdm(range(L),colour='#8060ff',ascii="░▒▓",bar_format="{l_bar}{bar:50}{r_bar}"):
        
        # check if mass blance is correct
        condition=NTPvals[i,:]
        comp=np.array(NTPvals[i,2:])
        lPass=_check_massbalance(Database,comp,PhaseName)
        if lPass:
            var.dConditionSys=condition
            input_condition(Unit,condition)
            minimize_single_phase()
    return None

def _check_massbalance(Database,Bulk,PhaseName):
    res= True
    var = Database
    p_type,p_id = var.PhaseNameSys[PhaseName]
    if p_type=='compd':
        compd=var.dStoichSpeciesCS[p_id,var.iElementDBIndex]
        for i in range(1,len(Bulk)):
            if compd[i]<1E-10: continue
            else: 
                scaler= abs(Bulk[i]-compd[i])/compd[i]
                if compd[i-1]<1E-10: scaler =0
                else:
                    scaler = abs(scaler - abs(Bulk[i-1]-compd[i-1])/compd[i-1])
            if scaler>1E-10: 
                res = False
                raise EquilibError(f'{PhaseName} cannot represent the input composition')
    else:
        res = False
        raise NameError(f'Error: PhaseSelection cannot identify type of {PhaseName}')
    return res