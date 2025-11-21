#!/usr/bin/env python3
import equilipy.equilifort as fort
from typing import Dict, List, Optional, Union
import numpy as np, os
from tqdm import tqdm
from .PostProcess import ResultScheil, combine_dicts
from .EquilibSingle import _equilib_single
from .Errors import EquilibError
from .FindTransition import find_transitions, find_first_transition
from .utils import _dict2np
from .variables import cPeriodicTable
# from mpi4py.futures import MPIPoolExecutor
# from mpi4py import MPI
from multiprocessing import Pool, set_start_method

def scheil_cooling(
    LiquidPhaseName: str,
    Database: Dict,
    Condition: Dict[str, Union[float, str]],
    dT: float = 5.0,
    Unit: List[str] = ['K', 'atm', 'moles'],
    ListOfPhases: Optional[List[str]] = None,
    StartFromLiquidus = True
) -> ResultScheil:
    """
    Performs a Scheil-Gulliver cooling simulation.

    This function simulates the solidification of a liquid phase under the assumptions
    of the Scheil-Gulliver model:
    1. No diffusion in solid phases.
    2. Infinitely fast diffusion in the liquid phase.
    3. Equilibrium is maintained at the solid-liquid interface.

    The simulation proceeds by cooling the system in discrete temperature steps (dT),
    calculating the equilibrium at each step, removing the newly formed solid phases,
    and then continuing the cooling with the remaining liquid.

    Args:
        LiquidPhaseName (str): The name of the liquid phase in the database.
        Database (Dict): The thermodynamic database dictionary.
        Condition (Dict[str, Union[float, str]]): A dictionary defining the initial
            conditions (temperature, pressure, and composition). This dictionary
            is not modified by the function.
        Unit (List[str], optional): A list specifying the units for temperature,
            pressure, and amount. Defaults to ['K', 'atm', 'moles'].
        dT (float, optional): The temperature step for each cooling iteration.
            Defaults to 5.0.
        IterMax (int, optional): The maximum number of iterations to prevent
            infinite loops. Defaults to 5000.
        ListOfPhases (Optional[List[str]], optional): A list of phase names to
            consider in the equilibrium calculation. If None, all phases in the
            database are considered. Defaults to None.

    Returns:
        ResultScheil: An object containing the detailed results of the Scheil
            simulation, including temperature, phase amounts, and compositions
            at each step.

    Raises:
        EquilibError: If the initial condition does not contain a stable liquid phase
            or if the temperature drops below 0 K during the simulation.
    """
    # Create a copy of the initial conditions to avoid modifying the user's dictionary
    current_condition = Condition.copy()
    Elements = list(Condition.keys())[2:]


    is_mass_based = Unit[2] in ['grams','kilograms','pounds','g','kg','lbs','mass fraction','weight fraction','wt%','wt.%']
    if is_mass_based:
        # atomic_masses = np.array([cPeriodicTable[el.strip()][1] for el in Elements])
        for el in Elements:
            current_condition[el] = current_condition[el]/cPeriodicTable[el.strip()][1]
    Unit_local = Unit.copy()
    Unit_local[2] = 'moles'
    # Start from liquidus
    if StartFromLiquidus:
        Tmax = current_condition['T']
        Tmin = current_condition['T']*0.1
        Ts=find_first_transition(Database, current_condition, Tmax,Tmin,Unit=Unit_local)
        current_condition['T'] = Ts + 0.1

    # --- Initial Equilibrium Calculation ---
    # Perform a single equilibrium calculation to establish the starting point
    _equilib_single(Database, current_condition, Unit=Unit_local, ListOfPhases=ListOfPhases)
    assemblage_old = fort.modulethermo.iassemblage.copy()
    
    res = ResultScheil()
    StablePhaseNames = [x.split('_mole')[0] for x in list(res.ScheilPhases_mole.keys())]

    # Check if the liquid phase is stable at the start
    if LiquidPhaseName not in StablePhaseNames:
        raise EquilibError(f"The specified liquid phase '{LiquidPhaseName}' is not stable at the initial conditions.")
    
    IterMax= int(current_condition['T']//dT)
    
    # --- Main Cooling Loop ---
    # The for loop iterates up to IterMax times, with a progress bar.
    for i in tqdm(range(IterMax), desc="Scheil Cooling", colour='#8060ff', ascii="░▒▓", bar_format="{l_bar}{bar:50}{r_bar}"):
        # Get the composition and amount of the liquid from the previous step based on the OutputUnit
        liquid_amount = res.ScheilPhases_mole[f'{LiquidPhaseName}_mole'][-1]
        if i==0:
            current_condition['T'] = res.T # This converts temperature unit to Kelvin
            Unit_local[0] = 'K'
            
            # --- UPDATED BLOCK ---
            # res.EquilibResult.Phases is a Dict[str, Phase] at i==0
            liquid_phase_obj = res.EquilibResult.Phases[LiquidPhaseName]

            # Use zip(*.items()) to unzip the dict into keys and values
            # This guarantees the order of components and values will match.
            if liquid_phase_obj.Elements_Xmole: # Handle non-empty dict
                components_tuple, values_tuple = zip(*liquid_phase_obj.Elements_Xmole.items())
                liquid_components = list(components_tuple)
                liquid_composition = list(values_tuple) # This is now the list of values
            else: # Handle empty dict
                liquid_components = []
                liquid_composition = []
            
        else:
            # --- UPDATED BLOCK ---
            # res.EquilibResult.Phases is a List[Dict[str, Phase]] for i > 0
            # Get the phases dictionary from the *last* iteration
            last_iteration_phases = res.EquilibResult.Phases[-1]
            liquid_phase_obj = last_iteration_phases[LiquidPhaseName]

            # Use the same zip(*.items()) method
            if liquid_phase_obj.Elements_Xmole: # Handle non-empty dict
                components_tuple, values_tuple = zip(*liquid_phase_obj.Elements_Xmole.items())
                liquid_components = list(components_tuple)
                liquid_composition = list(values_tuple) # This is the list of values
            else: # Handle empty dict
                liquid_components = []
                liquid_composition = []

        # --- Termination Conditions ---
        if liquid_amount < 1e-5:
            # Solidification is considered complete
            break
        
        current_temp = float(current_condition['T'])
        if current_temp <= 0 and 'K' in Unit_local[0]:
            raise EquilibError('Scheil-Gulliver solidification results in melting point below 0 K. Check phase selections.')

        # --- Update Conditions for Next Step ---
        # Decrease temperature
        current_condition['T'] = current_temp - dT
        
        # Update composition based on remaining liquid
        new_ni = np.array(liquid_composition) * liquid_amount
        for i, el in enumerate(liquid_components):
            current_condition[el] = new_ni[i]

        # --- Perform Equilibrium Calculation for the Current Step ---
        _equilib_single(Database, current_condition, Unit=Unit_local, ListOfPhases=ListOfPhases)
        
        # --- Before updating the result, if there is any transition occurs, calculate the transition point
        assemblage_new = fort.modulethermo.iassemblage.copy()
        if set(assemblage_new)!=set(assemblage_old):
            # # Transition occured, Reset variables
            Tmax = float(current_temp)
            Tmin = float(current_condition['T'])
            current_condition['T'] = Tmax
            Ts=find_transitions(Database, current_condition, Tmax,Tmin,Unit=Unit_local)
            if abs(Tmax-Ts[0])<0.1: Ts = Ts[1:]

            for j, T_trans in enumerate(Ts):
                #Update transition point
                current_condition['T'] = T_trans
                _equilib_single(Database, current_condition, Unit=Unit_local,ListOfPhases=ListOfPhases)
                if j<len(Ts)-1:
                    res.append_output()
                    # Get the composition and amount of the liquid from the previous step
                    liquid_amount = res.ScheilPhases_mole[f'{LiquidPhaseName}_mole'][-1]
                    last_iteration_phases = res.EquilibResult.Phases[-1]
                    liquid_phase_obj = last_iteration_phases[LiquidPhaseName]

                    # Use the same zip(*.items()) method
                    if liquid_phase_obj.Elements_Xmole: # Handle non-empty dict
                        components_tuple, values_tuple = zip(*liquid_phase_obj.Elements_Xmole.items())
                        liquid_components = list(components_tuple)
                        liquid_composition = list(values_tuple) # This is the list of values
                    else: # Handle empty dict
                        liquid_components = []
                        liquid_composition = []
                    
                    # Update composition based on remaining liquid
                    new_ni = np.array(liquid_composition) * liquid_amount
                    for k, el in enumerate(liquid_components):
                        current_condition[el] = new_ni[k]
            current_condition['T'] =float(Tmin)
        assemblage_new = fort.modulethermo.iassemblage.copy()
        assemblage_old=assemblage_new.copy()
        res.append_output()
        
        # Check if liquid phase disappeared in the last step
        if LiquidPhaseName not in res.EquilibResult.StablePhases[-1]['Name']:
            break
    else:
        # This block runs only if the for loop completes without a 'break'
        print(f'Warning: Reached maximum iterations ({IterMax}) before solidification completed.')

    # Finalize results
    res.update_scheilconstituents()
    # Reorder
    header = ['task_id']
    elements = list(Condition.keys())[2:]
    header.append('T_Delta [K]')
    header.extend([f'{x} [sp-mol]' for x in elements])
    header.extend([f'{x} [g]' for x in elements])

    ScheilConstituents = [key.split(' [sp-mol]')[0] for key in res.ScheilConstituents.keys() if key not in header and '[sp-mol]' in key]
    ScheilConstituents = [key for key in ScheilConstituents if key.strip()]

    primaries = [constituent for constituent in ScheilConstituents if '+' not in constituent and constituent!='']
    eutectics = [constituent for constituent in ScheilConstituents if '+' in constituent ]

    header.extend([f'{x} [sp-mol]' for x in primaries])
    header.extend([f'{x} [sp-mol]' for x in eutectics])
    header.extend([f'{x} [g]' for x in primaries])
    header.extend([f'{x} [g]' for x in eutectics])
    res.ScheilConstituents = {key: res.ScheilConstituents[key] for key in header if key in res.ScheilConstituents}
    return res


def _scheil_batch_input(
    LiquidPhaseName: str,
    Database: Dict,
    Condition: Dict[str, Union[float, str]],
    dT: float = 5.0,
    Unit: List[str] = ['K', 'atm', 'moles'],
    ListOfPhases: Optional[List[str]] = None,
    StartFromLiquidus = True,
    nPerBatch: int=1)-> List:
    """
    Break fullrange up into smaller sets of ranges that cover all
    the same numbers.
    """
    res = []
    header,NTP=_dict2np(Condition)
    n,m=NTP.shape
    fullrange=[0,n]

    for i in range(fullrange[0], fullrange[1], nPerBatch):
        subcondition=dict({})
        for j, head in enumerate(header):
            subcondition[head] = NTP[i:min(i+nPerBatch, fullrange[1]),j]
        res.append( [LiquidPhaseName,Database,subcondition,dT,Unit,ListOfPhases,StartFromLiquidus] ) 
    return res

def _scheil_batch(
    LiquidPhaseName: str,
    Database: Dict,
    Condition: Dict[str, Union[float, str]],
    dT: float = 5.0,
    Unit: List[str] = ['K', 'atm', 'moles'],
    ListOfPhases: Optional[List[str]] = None,
    StartFromLiquidus = True
) -> Dict:
    try: 
        res = scheil_cooling(
            LiquidPhaseName,
            Database,
            Condition,
            dT,
            Unit,
            ListOfPhases,
            StartFromLiquidus
        )
        return res.to_dict()
    except EquilibError:
        return {'Error': 'Equilibrium calculation failed during Scheil cooling.'}
    
def _scheil_constituent_batch(
    LiquidPhaseName: str,
    Database: Dict,
    Condition: Dict[str, Union[float, str]],
    dT: float = 5.0,
    Unit: List[str] = ['K', 'atm', 'moles'],
    ListOfPhases: Optional[List[str]] = None,
    StartFromLiquidus = True
) -> Dict:
    try: 
        res = scheil_cooling(
            LiquidPhaseName,
            Database,
            Condition,
            dT,
            Unit,
            ListOfPhases,
            StartFromLiquidus
        )
        return res.ScheilConstituents.copy()
    except EquilibError:
        return {'Error': 'Equilibrium calculation failed during Scheil cooling.'}
    
def _scheil_singlenode(arg,nCPU)-> List[Dict]:
    set_start_method('fork')
    with Pool(processes=nCPU) as pool:
        res = pool.starmap(_scheil_batch, arg)
    return res

def _scheil_constituent_singlenode(arg,nCPU)-> List[Dict]:
    set_start_method('fork')
    with Pool(processes=nCPU) as pool:
        res = pool.starmap(_scheil_constituent_batch, arg)
    return res
# def _scheil_multinodes(arg,nCPU)-> List[Dict]:
#     print(f'Process with mpi: {nCPU} processors')
#     with MPIPoolExecutor(max_workers=nCPU) as pool:
#         res = pool.starmap(_scheil_batch, arg)
#     return res

def scheil_batch(
    LiquidPhaseName: str,
    Database: Dict,
    Condition: Dict[str, Union[float, str]],
    dT: float = 5.0,
    Unit: List[str] = ['K', 'atm', 'moles'],
    ListOfPhases: Optional[List[str]] = None,
    StartFromLiquidus = True,
    nCPU: int=os.cpu_count(),
    nPerBatch: int=1) -> Dict:

    # Get a batch input
    arg =_scheil_batch_input(LiquidPhaseName,Database,Condition,dT,Unit,ListOfPhases,StartFromLiquidus,nPerBatch)
    res_mpi=_scheil_singlenode(arg,nCPU)
    res = combine_dicts(res_mpi)

    res_keys = list(res.keys())
    phases = [x.split('_Amount_mole [sp-mol]')[0] for x in res_keys if '_Amount_mole [sp-mol]' in x]
    cAmount_n = [f'{phase}_Amount_mole [sp-mol]' for phase in phases]
    cAmount_w = [f'{phase}_Amount_mass [g]' for phase in phases]
    cEndmember_x = [f'{phase}_Endmembers_Xmole_' for phase in phases]
    cEndmember_w = [f'{phase}_Endmembers_Xmass_' for phase in phases]
    cElement_x = [f'{phase}_Elements_Xmole_' for phase in phases]
    cElement_w = [f'{phase}_Elements_Xmass_' for phase in phases]

    header = ['task_id']
    elements = list(Condition.keys())[2:]
    header.append('T [K]')
    header.append('P [atm]')
    header.extend([f'{x} [sp-mol]' for x in elements])
    header.extend([f'{x} [g]' for x in elements])
    header.extend(['Label','fl [sp-mol/mol]','fs [sp-mol/mol]','fl [g/g]','fs [g/g]','G [J]','H [J]','S [J/K]','Cp [J/K]'])
    header.extend(cAmount_n)
    header.extend(cAmount_w)

    compositions= []
    compositions.extend(cEndmember_x)
    compositions.extend(cEndmember_w)
    compositions.extend(cElement_x)
    compositions.extend(cElement_w)
    # Re-arrange the dictionary according to the header
    reordered_res = {key: res[key] for key in header if key in res}
    
    # Add the remaining items from the original dictionary
    for composition in compositions:
        for key, value in res.items():
            if key.startswith(composition):
                reordered_res[key] = value

    return reordered_res

def scheil_constituent_batch(
    LiquidPhaseName: str,
    Database: Dict,
    Condition: Dict[str, Union[float, str]],
    dT: float = 5.0,
    Unit: List[str] = ['K', 'atm', 'moles'],
    ListOfPhases: Optional[List[str]] = None,
    StartFromLiquidus = True,
    nCPU: int=os.cpu_count(),
    nPerBatch: int=1) -> Dict:

    # Get a batch input
    arg =_scheil_batch_input(LiquidPhaseName,Database,Condition,dT,Unit,ListOfPhases,StartFromLiquidus,nPerBatch)
    res_mpi=_scheil_constituent_singlenode(arg,nCPU)
    res = combine_dicts(res_mpi)
    header = ['task_id']
    elements = list(Condition.keys())[2:]
    header.append('T_Delta [K]')
    header.extend([f'{x} [sp-mol]' for x in elements])
    header.extend([f'{x} [g]' for x in elements])

    ScheilConstituents = [key.split(' [sp-mol]')[0] for key in res.keys() if key not in header and '[sp-mol]' in key]
    ScheilConstituents = [key for key in ScheilConstituents if key.strip()]

    primaries = [constituent for constituent in ScheilConstituents if '+' not in constituent and constituent!='']
    eutectics = [constituent for constituent in ScheilConstituents if '+' in constituent ]

    header.extend([f'{x} [sp-mol]' for x in primaries])
    header.extend([f'{x} [sp-mol]' for x in eutectics])
    header.extend([f'{x} [g]' for x in primaries])
    header.extend([f'{x} [g]' for x in eutectics])
    reordered_res = {key: res[key] for key in header if key in res}
    
    
    # # Add the remaining items from the original dictionary
    # for key, value in res.items():
    #     if key not in reordered_res:
    #         reordered_res[key] = value

    return reordered_res