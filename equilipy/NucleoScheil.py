
from typing import Dict, List, Union, Optional
from .PostProcess import ResultScheil,endmembers2elements, combine_dicts
from .EquilibSingle import _equilib_single,equilib_single
from .FindTransition import find_first_transition, find_transitions
from .utils import _dict2np
from .Errors import EquilibError
from tqdm import tqdm
from multiprocessing import Pool, set_start_method
import copy, os
import numpy as np
from .variables import cPeriodicTable


#============================================================================
# NucleoScheil functions
def _init_nuclei(LiquidPhaseName: str,
    Database: Dict,
    Condition: Dict[str, Union[float, str]],
    CriticalUndercooling: Dict[str, float],
    Unit: List[str] = ['K', 'atm', 'moles']
    ) -> (float, str):

    res = []
    solid_phases = list(CriticalUndercooling.keys())

    for phase in solid_phases:
        # 1 Get liquidus temperature for each solid phase
        T_liquidus = find_first_transition(Database, Condition, Condition['T'], 10, Unit=Unit, ListOfPhases=[LiquidPhaseName, phase])
        # Append a tuple of (float, str)
        res.append((T_liquidus - CriticalUndercooling[phase], phase))

    # Sort the list of tuples in descending order based on the first element (the temperature)
    res.sort(key=lambda x: x[0], reverse=True)

    # The highest temperature nucleation event is now the first item
    return res[0]


def _update_liquid_composition(LiquidPhaseName: str,
    LatestScheilResult: ResultScheil,
    Condition: Dict[str, Union[float, str]]) -> Dict[str, float]:

    res = Condition.copy()
    LiqProp=LatestScheilResult.EquilibResult.Phases[-1][LiquidPhaseName] # Input
    LiqComposition = LiqProp.Elements_Xmole
    # LiqAmount = LatestScheilResult.EquilibResult.Phases[-1][LiquidPhaseName].Amount_mole
    LiqAmount=LatestScheilResult.ScheilPhases_mole[f'{LiquidPhaseName}_mole'][-1]
    if LiqAmount<=0: LiqAmount = 0.0
    for j,element in enumerate(LiqComposition.keys()):
        res[element]= np.round(float(LiqComposition[element]) *LiqAmount,10)
    return res

def _get_nuclei_candiates(LiquidPhaseName: str,
    Database: Dict,
    Condition: Dict[str, Union[float, str]],
    CriticalUndercooling: Dict[str, float],
    Unit: List[str] = ['K', 'atm', 'moles'])-> Dict[str, float]:
    res = {}
    NTP_local = Condition.copy()
    T_current = Condition['T']
    solid_phases= list(CriticalUndercooling.keys())
    # 3.2 Identify which phase candidates nucleate next
    for phase in solid_phases:
        NTP_local['T'] =  T_current + CriticalUndercooling[phase]
        res_temp=equilib_single(Database, NTP_local, Unit=Unit, ListOfPhases=[LiquidPhaseName,phase],CalcHeatCapacity=False)
        if phase in res_temp.StablePhases['Name']: res[phase] = res_temp.Phases[phase].Amount_mole
    return res


def _check_eutectic(
    LiquidPhaseName: str,
    LatestScheilResult: ResultScheil,
    nuclei_current: List[str],
    Database: Dict,
    Condition: Dict[str, Union[float, str]],
    CriticalUndercooling: Dict[str, float],
    Unit: List[str] = ['K', 'atm', 'moles'])-> (ResultScheil,bool):

    # First make a copy of the result
    res = copy.deepcopy(LatestScheilResult)
    lfinal = False


    # Get liquid composition and nuclei candidate
    NTP_local = _update_liquid_composition(LiquidPhaseName, res, Condition)
    nuclei_more_candidates=_get_nuclei_candiates(LiquidPhaseName,Database,NTP_local,CriticalUndercooling,Unit=Unit)
    candidates = list(nuclei_more_candidates.keys())
    phase_selection = [LiquidPhaseName]+candidates

    nuclei_final = nuclei_current.copy()
    nuclei_prev = nuclei_current.copy()
    Unit_local = Unit.copy()
    Elements = list(Condition.keys())[2:]
    NTP_liq = NTP_local.copy()
    if len(candidates)==0: 
        return res, lfinal
    else:
        #Check if it is the final point
        for j in range(len(Elements)):
            eqres=equilib_single(Database, NTP_liq, Unit=Unit, ListOfPhases=phase_selection,CalcHeatCapacity=False)
            LiqFraction = eqres.Phases[LiquidPhaseName].Elements_Xmole
            LiqAmount =float( eqres.Phases[LiquidPhaseName].Amount_mole)
            if LiqAmount<1E-5:
                nuclei_final.extend([x for x in candidates if x not in nuclei_final])
                lfinal = True
                break
            
            for e in Elements:
                NTP_liq[e] = float(LiqFraction[e]*LiqAmount)

            nuclei_more_candidates=_get_nuclei_candiates(LiquidPhaseName,Database,NTP_liq,CriticalUndercooling,Unit=Unit_local)
            candidates = list(nuclei_more_candidates.keys())
            
            if set(candidates)==set(nuclei_prev):
                # Condition 1: If the nuclei_candidates are same as previous, continue with previous nuclei
                phase_selection = [LiquidPhaseName]+nuclei_prev
            else:
                # Condition 2: If new nuclei are identified, conduct metastable equilib calculation with new nuclei added
                phase_selection = [LiquidPhaseName]+candidates

            nuclei_final.extend([x for x in candidates if x not in nuclei_final])
            nuclei_prev = [x for x in phase_selection if x !=LiquidPhaseName]
            if len(nuclei_final)>= len(Elements):
                #Check if final nuclei is correct
                eqres=equilib_single(Database, NTP_liq, Unit=Unit, ListOfPhases=nuclei_final,CalcHeatCapacity=False)
                if not np.isnan(eqres.G): 
                    lfinal = True
                    break
        if lfinal:
            _equilib_single(Database, NTP_local, Unit=Unit,ListOfPhases=[LiquidPhaseName]+nuclei_final)
            res.append_output()
            return res, lfinal
        else:
            return res, lfinal
        


def nucleoscheil_cooling(
    LiquidPhaseName: str,
    Database: Dict,
    Condition: Dict[str, Union[float, str]],
    CriticalUndercooling: Dict[str, float],
    dT: float = 5.0,
    Unit: List[str] = ['K', 'atm', 'moles']) -> ResultScheil:

    
    # Get phases
    solid_phases= list(CriticalUndercooling.keys())
    NTP_local = Condition.copy()
    Elements = list(Condition.keys())[2:]
    T_current = NTP_local['T']
    Unit_local = Unit.copy()

    # if it is mass based, convert to moles
    is_mass_based = Unit[2] in ['grams','kilograms','pounds','g','kg','lbs','mass fraction','weight fraction','wt%','wt.%']
    if is_mass_based:
        Unit_local[2] = 'moles'
        for el in Elements:
            NTP_local[el] = NTP_local[el]/cPeriodicTable[el.strip()][1]

    # Step1: Initialize NucleoScheil and get the liquidus
    T_init=find_first_transition(Database, NTP_local, T_current,T_current*0.1,Unit=Unit_local,ListOfPhases=[LiquidPhaseName]+solid_phases)

    T_current = T_init +0.1
    NTP_local['T'] = T_current
    _equilib_single(Database, NTP_local, Unit=Unit_local, ListOfPhases=[LiquidPhaseName])
    res = ResultScheil()
    # now all units are converted to K
    Unit_local[0] = 'K'
    T_current = res.T
    NTP_local['T'] = T_current


    # Step2: Get the first nucleating phase and temperature
    _equilib_single(Database, NTP_local, Unit=Unit_local, ListOfPhases=[LiquidPhaseName])
    res.append_output()

    T_current, nuclei_init = _init_nuclei(LiquidPhaseName,Database, NTP_local, CriticalUndercooling, Unit=Unit_local)
    NTP_local['T'] = T_current
    res.data[-1].T = float(T_current)

    _equilib_single(Database, NTP_local, Unit=Unit_local, ListOfPhases=[LiquidPhaseName,nuclei_init])
    res.append_output()
    nuclei_prev = [nuclei_init]
    phase_selection = [LiquidPhaseName]+nuclei_prev

    # Step2 Loop: Proceed to next temperature step
    for i in tqdm(range(int(T_init/dT)), desc="NucleoScheil Cooling", colour='#8060ff', ascii="░▒▓", bar_format="{l_bar}{bar:50}{r_bar}"):

        # Loop Step 1: Update temperature
        T_current = T_current - dT
        NTP_local['T'] = T_current

        # Loop Step 2: Get nuclei
        NTP_liq = _update_liquid_composition(LiquidPhaseName, res, NTP_local)
        nuclei_candidates=_get_nuclei_candiates(LiquidPhaseName,Database,NTP_liq,CriticalUndercooling,Unit=Unit_local)
        nuclei_candidates=list(nuclei_candidates.keys())
        nuclei_number=len(nuclei_candidates)
        if nuclei_number==0 or len(nuclei_prev)==0 or set(nuclei_candidates)==set(nuclei_prev):
            # Condition 1: If the nuclei_number is zero or if the nuclei_candidates are same as previous, continue with previous nuclei
            if len(nuclei_prev)!=0: phase_selection = [LiquidPhaseName]+nuclei_prev
            nuclei_current = nuclei_prev.copy()
        else:
            # Condition 2: If new nuclei are identified, conduct metastable equilib calculation with new nuclei added
            phase_selection = [LiquidPhaseName]+nuclei_candidates
            nuclei_current = nuclei_candidates.copy()

        # Check if it is final temperature
        res, lfinal= _check_eutectic(LiquidPhaseName,res,nuclei_current,Database,NTP_liq,CriticalUndercooling,Unit=Unit_local)
        if lfinal and LiquidPhaseName not in res.data[-1].label.split('+'):
            break
        else:
            _equilib_single(Database, NTP_liq, Unit=Unit_local, ListOfPhases=phase_selection)
            res.append_output()
            nuclei_prev = [ x for x in phase_selection if x !=LiquidPhaseName]
            NTP_local = NTP_liq
            T_current = NTP_liq['T']

        # # if not move on to the next
        # _equilib_single(Database, NTP_liq, Unit=Unit_local, ListOfPhases=phase_selection)
        # res.append_output()
        # nuclei_prev = [ x for x in phase_selection if x !=LiquidPhaseName]

        # # Check if there is any other phases nucleating at this temperature
        # if LiquidPhaseName in res.data[-1].label.split('+'):
        #         NTP_liq = _update_liquid_composition(LiquidPhaseName, res, NTP_liq)
        #         res, lfinal= _check_eutectic(LiquidPhaseName,res,nuclei_current,Database,NTP_liq,CriticalUndercooling,Unit=Unit_local)
        # else:
        #     lfinal = True

        # if lfinal: break
        # else: NTP_local = NTP_liq
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


def _nucleoscheil_batch_input(
    LiquidPhaseName: str,
    Database: Dict,
    Condition: Dict[str, Union[float, str]],
    CriticalUndercooling: Dict[str, Union[float, str]],
    dT: float = 5.0,
    Unit: List[str] = ['K', 'atm', 'moles'],
    nPerBatch: int=1
    )-> List:
    """
    Break fullrange up into smaller sets of ranges that cover all
    the same numbers.
    """
    res = []
    header,NTP=_dict2np(Condition)
    header_uc, UC = _dict2np(CriticalUndercooling)
    n,m=NTP.shape
    fullrange=[0,n]

    for i in range(fullrange[0], fullrange[1], nPerBatch):
        subcondition=dict({})
        subundercooling = dict({})
        for j, head in enumerate(header):
            subcondition[head] = NTP[i:min(i+nPerBatch, fullrange[1]),j]
        for j, head in enumerate(header_uc):
            subundercooling[head] = UC[i:min(i+nPerBatch, fullrange[1]),j]
        res.append( [LiquidPhaseName,Database,subcondition,subundercooling,dT,Unit] ) 
    return res

def _nucleoscheil_batch(
    LiquidPhaseName: str,
    Database: Dict,
    Condition: Dict[str, Union[float, str]],
    CriticalUndercooling: Dict[str, Union[float, str]],
    dT: float = 5.0,
    Unit: List[str] = ['K', 'atm', 'moles']) -> Dict:
    try: 
        res = nucleoscheil_cooling(LiquidPhaseName,Database,Condition,CriticalUndercooling,dT,Unit)
        return res.to_dict()
    except EquilibError:
        return {'Error': 'Equilibrium calculation failed during NucleoScheil cooling.'}
    
def _nucleoscheil_constituent_batch(
    LiquidPhaseName: str,
    Database: Dict,
    Condition: Dict[str, Union[float, str]],
    CriticalUndercooling: Dict[str, Union[float, str]],
    dT: float = 5.0,
    Unit: List[str] = ['K', 'atm', 'moles']) -> Dict:
    try: 
        res = nucleoscheil_cooling(LiquidPhaseName,Database,Condition,CriticalUndercooling,dT,Unit)
        return res.ScheilConstituents.copy()
    except EquilibError:
        return {'Error': 'Equilibrium calculation failed during NucleoScheil cooling.'}
    
def _nucleoscheil_singlenode(arg,nCPU)-> List[Dict]:
    set_start_method('fork')
    with Pool(processes=nCPU) as pool:
        res = pool.starmap(_nucleoscheil_batch, arg)
    return res

def _nucleoscheil_constituent_singlenode(arg,nCPU)-> List[Dict]:
    set_start_method('fork')
    with Pool(processes=nCPU) as pool:
        res = pool.starmap(_nucleoscheil_constituent_batch, arg)
    return res

def nucleoscheil_batch(
    LiquidPhaseName: str,
    Database: Dict,
    Condition: Dict[str, Union[float, str]],
    CriticalUndercooling: Dict[str, Union[float, str]],
    dT: float = 5.0,
    Unit: List[str] = ['K', 'atm', 'moles'],
    nCPU: int=os.cpu_count(),
    nPerBatch: int=1) -> Dict:

    # Get a batch input
    arg =_nucleoscheil_batch_input(LiquidPhaseName,Database,Condition,CriticalUndercooling,dT,Unit,nPerBatch)
    res_mpi=_nucleoscheil_singlenode(arg,nCPU)
    res = combine_dicts(res_mpi)

    phases = [LiquidPhaseName] + list(CriticalUndercooling.keys())
    cAmount_n = [f'{phase}_Amount_mole [sp-mol]' for phase in phases]
    cAmount_w = [f'{phase}_Amount_mass [g]' for phase in phases]
    cEndmember_x = [f'{phase}_Endmembers_Xmole_' for phase in phases]
    cEndmember_w = [f'{phase}_Endmembers_Xmass_' for phase in phases]
    cElement_x = [f'{phase}_Elements_Xmole_' for phase in phases]
    cElement_w = [f'{phase}_Elements_Xmass_' for phase in phases]

    header = ['task_id']
    header.append('T [K]')
    header.append('P [atm]')
    elements = list(Condition.keys())[2:]
    
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

def nucleoscheil_constituent_batch(
    LiquidPhaseName: str,
    Database: Dict,
    Condition: Dict[str, Union[float, str]],
    CriticalUndercooling: Dict[str, Union[float, str]],
    dT: float = 5.0,
    Unit: List[str] = ['K', 'atm', 'moles'],
    nCPU: int=os.cpu_count(),
    nPerBatch: int=1) -> Dict:

    # Get a batch input
    arg =_nucleoscheil_batch_input(LiquidPhaseName,Database,Condition,CriticalUndercooling,dT,Unit,nPerBatch)
    res_mpi=_nucleoscheil_constituent_singlenode(arg,nCPU)
    res = combine_dicts(res_mpi)
    
    header = ['task_id']
    elements = list(Condition.keys())[2:]
    header.append('T_Delta [K]')
    header.extend([f'{x} [sp-mol]' for x in elements])
    header.extend([f'{x} [g]' for x in elements])
    
    # Re-arrange the dictionary according to the header
    
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