import numpy as np
import equilipy.equilifort as fort
import equilipy.variables as var
from .Errors import *
from pydantic import BaseModel, Field, ConfigDict
from typing import List, Dict, Any, Optional, Union, TypedDict
from collections import defaultdict


def _count_unstable_compounds(i):    
    n=0
    iphase=fort.modulethermo.iphase
    iLastSoln=fort.modulethermo.nspeciesphase[len(var.iSys2DBSoln)]
    for k in range(iLastSoln,i):
        if iphase[k]<0:
            n = n+1
    return n

def _get_assemblage_name(AssemblageIDs):
    '''
    An internal function that spits out phase name based on Assemblage IDs

    Variables
    =========

    Input
    AssemblageIDs: list of phase assemblage id defined after selecting phases based on input elements

    Output
    AssemblageNames
    '''
    AssemblageNames=list([])
    for i in AssemblageIDs:
        if i ==0:
            # Empty phase: place holder
            AssemblageNames.append("{:<1}".format(''))
        elif i<0:
            # Solution Phase
            # Note that the solution phase id changes with phase selection
            AssemblageNames.append(var.cPhaseNameSys[int(-(i+1))].strip())
        else:
            #Compound Phase
            
            # Note that the compund id doesn't change with phase selection
            # Correction is made by counting the number of ignored compound phases
            nSoln=len(var.iSys2DBSoln)
            iLastSoln=fort.modulethermo.nspeciesphase[nSoln]
            id=-iLastSoln+i+nSoln-1
            n = _count_unstable_compounds(i)
            id = int(id-n)
            AssemblageNames.append(var.cPhaseNameSys[id].strip())
    return AssemblageNames

def combine_dicts(
    dict_list: List[Dict[str, float]]
) -> Dict[str, np.ndarray]:
    """
    Combines a list of "row" dictionaries into a single "columnar"
    dictionary of numpy arrays.
    
    This is the high-performance, dependency-free alternative to
    `pd.DataFrame(list_of_dicts).fillna(0.0)`.

    Args:
        dict_list: A list of dictionaries, where each dict is a "row".
                   e.g., [{'A': 1, 'B': 2}, {'B': 3, 'C': 4}]

    Returns:
        A single dictionary where keys are all unique keys found
        and values are numpy arrays, padded with 0.0.
        e.g., {'A': np.array([1., 0.]),
               'B': np.array([2., 3.]),
               'C': np.array([0., 4.])}
    """
    if not dict_list:
        return {}
        
    all_keys = set().union(*dict_list)
    col_dict = defaultdict(list)

    # Iterate through each row dictionary one at a time
    for i, row_dict in enumerate(dict_list):
        # Step 1: Determine the size 'n' for the current row
        n = 1
        for key, value in row_dict.items():
            if key.startswith('T'):
                try:
                    if hasattr(value, '__iter__') and not isinstance(value, str):
                        n = len(value)
                except TypeError:
                    pass # n remains 1
                break 
        
        # Step 2: For the current row, iterate through ALL possible keys
        for key in all_keys:
            if key in row_dict:
                # Key exists in the current row
                value = row_dict[key]
                # Broadcast scalar values to the size of the row
                if not (hasattr(value, '__iter__') and not isinstance(value, str)):
                    value = [value] * n
                col_dict[key].extend(np.array(value))
            else:
                # Key is missing in the current row, so add padding of size 'n'
                col_dict[key].extend(np.zeros(n))

        # Handle 'task_id' for the current row
        col_dict['task_id'].extend(np.ones(n, dtype=int) * (i + 1))
    
    final_dict ={}
    # Convert the flat lists to numpy arrays
    for key, values in col_dict.items():
        final_dict[key] = np.array(values)
    # final_dict = {key: np.array(values) for key, values in col_dict.items()}
    
    return final_dict

def endmembers2elements(
    PhaseName: str,
    Endmembers: Dict[str, float],
    Elements: list,
    UnitOut: List[str] = ['K', 'atm', 'moles'],
    IsSolution: bool=True ) -> Dict[str, float]:
   
    #1. If the elements in endmembers and input condition are same, return
    EndmembersName=list(Endmembers.keys())

    if set(EndmembersName)==set(Elements): return Endmembers

    endmemeber_amounts = np.array([Endmembers.get(name) for name in EndmembersName])

    is_mass_based = UnitOut[2] in ['grams','kilograms','pounds','g','kg','lbs','mass fraction','weight fraction','wt%','wt.%']
    if is_mass_based:
        atomic_masses = np.array([var.cPeriodicTable[el.strip()][1] for el in Elements])
    
    
    #2. If not, calulate element fraction from database

    iElementDBIndex = var.iElementDBIndex.copy()
    dStoichSpeciesPhase = var.dStoichSpeciesCS[:,iElementDBIndex]

    if IsSolution:
        # 2.1 Get the phase index
        DB_soln_names=var.cSolnPhaseNameCS.copy()
        DB_soln_names = [x.strip() for x in DB_soln_names]

        counts = defaultdict(int)
        all_soln_names = []
        for name in DB_soln_names:
            counts[name] += 1
            count = counts[name]
            if count > 1:
                all_soln_names.append(f'{name}#{count}')
            else:
                all_soln_names.append(name)

        # iSys=var.cSolnPhaseNameCS.index("{:<25}".format(PhaseName))
        
        
        iSys= all_soln_names.index(PhaseName)
        nSpeciesIndex=var.nSpeciesPhaseCS.copy()
        nSpeciesIndex = np.append([0], nSpeciesIndex)
        iFirstSys = nSpeciesIndex[iSys]
        iLastSys = nSpeciesIndex[iSys+1]
        Temp=[str(x).strip() for x in np.array(var.cEndmemberNameCS)[iFirstSys:iLastSys]]

        # iSpeciesDBIndex = np.where(np.isin(var.iSys2DBSpeciesDefault,np.arange(iFirstSys,iLastSys)))[0]
        iSpeciesDBIndex = np.array([iFirstSys+Temp.index(x) for x in EndmembersName])
        dStoichSpeciesPhase = dStoichSpeciesPhase[iSpeciesDBIndex,:]
        ElementValues = np.matmul(endmemeber_amounts,dStoichSpeciesPhase)
        
        
    else:
        stripped_db_names = [name.strip() for name in var.cEndmemberNameCS]
        iSpeciesDBIndex = stripped_db_names.index(EndmembersName[0])
        dStoichSpeciesPhase = dStoichSpeciesPhase[iSpeciesDBIndex,:]
        ElementValues = endmemeber_amounts*dStoichSpeciesPhase
     
    
    # 2.3 Refine it for the system
    
    ElementValues = ElementValues / np.sum(ElementValues)
    
    if is_mass_based:
        ElementValues_mass = ElementValues*atomic_masses
        ElementValues = ElementValues_mass / np.sum(ElementValues_mass)
    
    output = dict(zip(Elements, ElementValues))
    return output
    


    
class Phase(BaseModel):
    """
    Pydantic model for a single phase (solution or compound).
    Holds the thermochemical state of the phase.
    
    UPDATED: Endmembers and xi are combined into a single xi dictionary.
    """
    ID: int
    Name: str
    Amount_mole: float = Field(default=0.0, ge=0.0)
    Amount_mass: float = Field(default=0.0, ge=0.0)
    Stability: float = Field(default=0.0, ge=0.0, le=1.0)
    Endmembers_Xmole: Dict[str, float] = Field(default_factory=dict) 
    Endmembers_Xmass: Dict[str, float] = Field(default_factory=dict) 
    Elements_Xmole: Dict[str, float] = Field(default_factory=dict) 
    Elements_Xmass: Dict[str, float] = Field(default_factory=dict)

    model_config = ConfigDict(arbitrary_types_allowed=True)

# --------------------------------------------------------------------------
# 2. Factory Function to Create a Phase from Fortran
# --------------------------------------------------------------------------

def create_phase_from_sys(iSys: int) -> Phase:
    """
    Factory function to create a Phase object by querying
    the equilipy Fortran and variable modules.
    
    UPDATED: Creates a dictionary for xi.
    """
    Name = str(var.cPhaseNameSys[iSys])
    Elements = np.array(var.cElementNameCS)[var.iElementDBIndex]
    Elements = [x.strip() for x in Elements]
    
    if iSys < len(var.iSys2DBSoln):
        # Solution Phase
        ID = -(iSys + 1)
        if ID in list(fort.modulethermo.iassemblage):
            Amount_mole = float(fort.modulethermo.dmolesphase.copy()[list(fort.modulethermo.iassemblage).index(ID)])
            Amount_mass = float(fort.modulethermoio.dgramphase.copy()[list(fort.modulethermo.iassemblage).index(ID)])
            Stability = 1.0
        else:
            Amount_mole = 0.0
            Amount_mass = 0.0
            Stability = 0.0
        
        iFirstSys = fort.modulethermo.nspeciesphase[iSys]
        iLastSys = fort.modulethermo.nspeciesphase[iSys + 1]

        idx_species = np.array(var.iSys2DBSpecies[iFirstSys:iLastSys])
        
        # --- MODIFIED ---
        # Create the xi dict by zipping names and values
        endmember_names = [str(x).strip() for x in list(np.array(var.cEndmemberNameCS)[idx_species])]
        xi_values = list(fort.modulethermo.dmolfraction[iFirstSys:iLastSys])
        wi_values = list(fort.modulethermoio.dgramfraction[iFirstSys:iLastSys])
        Endmembers_Xmole = dict(zip(endmember_names, xi_values))
        Endmembers_Xmass = dict(zip(endmember_names, wi_values))
        Elements_Xmole =  endmembers2elements(Name, Endmembers_Xmole,Elements)
        Elements_Xmass =  endmembers2elements(Name, Endmembers_Xmole,Elements,UnitOut = ['K', 'atm', 'g'])

    else:
        # Compound phase
        ID = fort.modulethermo.nspeciesphase[len(var.iSys2DBSoln)] + iSys - len(var.iSys2DBSoln) + 1
        if ID in list(fort.modulethermo.iassemblage):
            Amount_mole = float(fort.modulethermo.dmolesphase[list(fort.modulethermo.iassemblage).index(ID)])
            Amount_mass = float(fort.modulethermoio.dgramphase[list(fort.modulethermo.iassemblage).index(ID)])
            Stability = 1.0
        else:
            Amount_mole = 0.0
            Amount_mass = 0.0
            Stability = 0.0
        
        # --- MODIFIED ---
        # Compound phase has only itself as the "endmember" with fraction 1.0
        Endmembers_Xmole = {Name.strip(): 1.0} 
        Endmembers_Xmass = {Name.strip(): 1.0} 
        Elements_Xmole =  endmembers2elements(Name, Endmembers_Xmole,Elements,IsSolution=False)
        Elements_Xmass =  endmembers2elements(Name, Endmembers_Xmole,Elements,UnitOut = ['K', 'atm', 'g'],IsSolution=False)
        
        
    return Phase(
        ID=ID,
        Name=Name.strip(),
        Amount_mole=Amount_mole,
        Amount_mass=Amount_mass,
        Stability=Stability,
        Endmembers_Xmole=Endmembers_Xmole,
        Endmembers_Xmass=Endmembers_Xmass,
        Elements_Xmole = Elements_Xmole,
        Elements_Xmass = Elements_Xmass
    )

# --------------------------------------------------------------------------
# 3. Pydantic Model for a Single "Row" of Results
# --------------------------------------------------------------------------
class StablePhasesDict(TypedDict):
    ID: np.ndarray
    Name: np.ndarray
    Amount_mole: np.ndarray
    Amount_mass: np.ndarray

# Type alias for the dictionary of Phases.
# Keys are dynamic phase names (e.g. 'LIQUID', 'FCC_A1').
PhaseDict = Dict[str, Phase]

class SinglePointResult(BaseModel):
    """
    Pydantic model for a single equilibrium calculation point.
    This represents one "row" of data.
    
    UPDATED: StablePhases dictionary now holds numpy arrays.
    """
    T: float
    P: float
    G: Optional[float] = None
    H: Optional[float] = None
    S: Optional[float] = None
    Cp: Optional[float] = None
    N_x: Dict[str, float] = Field(default_factory=dict) # Element composition
    N_w: Dict[str, float] = Field(default_factory=dict) # Element composition
    
    # UPDATED: The values in this dict will be numpy arrays
    StablePhases: StablePhasesDict = Field(default_factory=TypedDict)
    
    # Holds data for ALL phases (stable and unstable)
    Phases: Dict[str, Phase] = Field(default_factory=dict) 

    model_config = ConfigDict(arbitrary_types_allowed=True)
        
    @classmethod
    def from_fortran(cls) -> "SinglePointResult":
        """
        Factory classmethod to create a new instance by reading
        from the equilipy fortran/variable modules.
        """
        try:
            G = float(fort.modulethermoio.dgibbsenergysys)
            H = float(fort.modulethermoio.denthalpysys)
            S = float(fort.modulethermoio.dentropysys)
            Cp = float(fort.modulethermoio.dheatcapacitysys)
            
            # --- MODIFIED ---
            # Create the StablePhases dictionary with numpy arrays
            stable_phases_dict = StablePhasesDict( {
                'ID': np.array(list(fort.modulethermo.iassemblage)),
                'Name': np.array(list(_get_assemblage_name(fort.modulethermo.iassemblage))),
                'Amount_mole': np.array(list(fort.modulethermo.dmolesphase)),
                'Amount_mass': np.array(list(fort.modulethermoio.dgramphase))
            })
            

            # Assign phase properties
            var.iPhaseSys = [] 
            all_phases = {}
            for i in range(len(var.iSys2DBSoln) + len(var.iSys2DBComp)):
                new_phase = create_phase_from_sys(i) # This now uses np.array
                all_phases[new_phase.Name] = new_phase
                var.iPhaseSys.append(new_phase.ID)
                
            return cls(
                T=float(fort.modulethermoio.dtemperature),
                P=float(fort.modulethermoio.dpressure),
                N_x=dict(zip([f'{x} [sp-mol]' for x in var.cComponentNameSys], fort.modulethermo.dmoleselement.copy()[var.iElementSysIndex])),
                N_w=dict(zip([f'{x} [g]' for x in var.cComponentNameSys], fort.modulethermoio.dgramelement.copy()[var.iElementSysIndex])),
                G=G, H=H, S=S, Cp=Cp,
                StablePhases=stable_phases_dict,
                Phases=all_phases
            )
        except Exception as e:
            print(f"Error reading from Fortran, creating error row: {e}")
            return cls.for_error()

    @classmethod
    def for_error(cls) -> "SinglePointResult":
        """
        Factory classmethod to create an instance representing a
        failed calculation, populated with NaNs.
        """
        # --- MODIFIED ---
        # Create the error StablePhases dictionary with numpy arrays
        stable_phases_dict = {
            'ID': np.array([np.nan]),
            'Name': np.array(['nan']),
            'Amount_mole': np.array([np.nan]),
            'Amount_mass':  np.array([np.nan])
        }

        return cls(
            T=float(fort.modulethermoio.dtemperature),
            P=float(fort.modulethermoio.dpressure),
            N_x=dict(zip([f'{x} [sp-mol]' for x in var.cComponentNameSys], fort.modulethermo.dmoleselement.copy()[var.iElementSysIndex])),
            N_w=dict(zip([f'{x} [g]' for x in var.cComponentNameSys], fort.modulethermoio.dgramelement.copy()[var.iElementSysIndex])),
            G=np.nan, H=np.nan, S=np.nan, Cp=np.nan,
            StablePhases=stable_phases_dict,
            Phases={}
        )

# --------------------------------------------------------------------------
# 4. The Revised "Result" Class (with reusable to_dict helpers)
# --------------------------------------------------------------------------

class Result:
    """
    A class to hold calculation results.
    
    UPDATED: to_dict now uses reusable helper methods for flattening
    dictionary attributes like 'xi', 'ai', 'mui', etc.
    """
    def __init__(self):
        self.data: Union[None, SinglePointResult, List[SinglePointResult]] = None

    # ... (append, append_output, append_error methods are unchanged) ...
    def append(self, other: "Result"):
        if other.data is None: return
        if isinstance(other.data, SinglePointResult):
            to_append = [other.data]
        elif isinstance(other.data, list):
            to_append = other.data 
        else: return
        if self.data is None:
            self.data = other.data
        elif isinstance(self.data, SinglePointResult):
            self.data = [self.data] + to_append
        elif isinstance(self.data, list):
            self.data.extend(to_append)

    def append_output(self):
        new_point = SinglePointResult.from_fortran()
        if self.data is None:
            self.data = new_point
        elif isinstance(self.data, SinglePointResult):
            self.data = [self.data, new_point]
        elif isinstance(self.data, list):
            self.data.append(new_point)
        
    def append_error(self):
        error_point = SinglePointResult.for_error()
        if self.data is None:
            self.data = error_point
        elif isinstance(self.data, SinglePointResult):
            self.data = [self.data, error_point]
        elif isinstance(self.data, list):
            self.data.append(error_point)
        
    # --- Property Accessors (Unchanged) ---
    @property
    def N_x(self) -> Union[Dict[str, float], List[Dict[str, float]]]:
        """
        Returns:
        - Dict[str, float] if single point.
        - List[Dict[str, float]] if multiple points.
        """
        if isinstance(self.data, SinglePointResult):
            return self.data.N_x
        if isinstance(self.data, list):
            return [iter.N_x for iter in self.data]
        return {}
    @property
    def N_w(self) -> Union[Dict[str, float], List[Dict[str, float]]]:
        """
        Returns:
        - Dict[str, float] if single point.
        - List[Dict[str, float]] if multiple points.
        """
        if isinstance(self.data, SinglePointResult):
            return self.data.N_w
        if isinstance(self.data, list):
            return [iter.N_w for iter in self.data]
        return {}
    @property
    def T(self) -> Union[float, List[float]]:
        if isinstance(self.data, SinglePointResult): return self.data.T
        if isinstance(self.data, list): return [iter.T for iter in self.data]
        return None
        
    @property
    def P(self) -> Union[float, List[float]]:
        if isinstance(self.data, SinglePointResult): return self.data.P
        if isinstance(self.data, list): return [iter.P for iter in self.data]
        return None

    # (G, H, S, Cp, StablePhases, Phases properties are also unchanged)
    # ... (omitted for brevity) ...

    @property
    def G(self) -> Union[Optional[float], List[Optional[float]]]:
        if isinstance(self.data, SinglePointResult): return self.data.G
        if isinstance(self.data, list): return [iter.G for iter in self.data]
        return None

    @property
    def H(self) -> Union[Optional[float], List[Optional[float]]]:
        if isinstance(self.data, SinglePointResult): return self.data.H
        if isinstance(self.data, list): return [iter.H for iter in self.data]
        return None

    @property
    def S(self) -> Union[Optional[float], List[Optional[float]]]:
        if isinstance(self.data, SinglePointResult): return self.data.S
        if isinstance(self.data, list): return [iter.S for iter in self.data]
        return None

    @property
    def Cp(self) -> Union[Optional[float], List[Optional[float]]]:
        if isinstance(self.data, SinglePointResult): return self.data.Cp
        if isinstance(self.data, list): return [iter.Cp for iter in self.data]
        return None

    @property
    def StablePhases(self) -> Union[StablePhasesDict, List[StablePhasesDict]]:
        if isinstance(self.data, SinglePointResult):
            return self.data.StablePhases
        if isinstance(self.data, list):
            return [iter.StablePhases for iter in self.data]
        return {} 

    @property
    def Phases(self) -> Union[Dict[str, Phase], List[Dict[str, Phase]]]:
        if isinstance(self.data, SinglePointResult):
            return self.data.Phases
        if isinstance(self.data, list):
            return [iter.Phases for iter in self.data]
        return {}
        
    # --- NEW: Reusable Helper Methods for to_dict ---
    
    def _discover_phase_sub_keys(self, iterations: List[SinglePointResult], attribute_name: str) -> Dict[str, set]:
        """
        Discovers all unique sub-keys for a given phase attribute dict.
        e.g., for 'Endmembers_Xmole', returns {'LIQUID': {'AL', 'FE'}, 'FCC_A1': {'AL', 'VA'}}
        """
        all_sub_keys = defaultdict(set)
        for iter_data in iterations:
            for phase_name, phase_obj in iter_data.Phases.items():
                attribute_dict = getattr(phase_obj, attribute_name, {})
                if attribute_dict:
                    all_sub_keys[phase_name].update(attribute_dict.keys())
        return all_sub_keys

    def _populate_flattened_phase_attribute(self, 
                                           dict_of_lists: defaultdict,
                                           phase_name: str,
                                           phase_obj: Optional[Phase],
                                           sub_keys: set,
                                           attribute_name: str,
                                           column_suffix: str):
        """
        Populates the dict_of_lists for one phase's flattened attribute.
        e.g., for phase 'LIQUID', attribute 'xi', suffix '_xi'.
        """
        for component in sub_keys: # e.g., component = 'AL'
            # e.g., 'LIQUID_xi_AL'
            output_key = f'{phase_name}{column_suffix}_{component}'
            
            if phase_obj:
                attribute_dict = getattr(phase_obj, attribute_name, {})
                value = attribute_dict.get(component, 0.0)
            else:
                value = 0.0
            
            dict_of_lists[output_key].append(value)

    # --- to_dict Method (Refactored) ---
    
    def to_dict(self) -> Dict[str, Any]:
        """
        Exports all results to a flattened dictionary, using reusable
        helpers to flatten phase attributes.
        """
        if self.data is None: return {}
            
        if isinstance(self.data, SinglePointResult):
            iterations = [self.data]
            is_single = True
        else:
            iterations = self.data
            is_single = False
            if not iterations: return {}

        # for iter in iterations:
        #     if iter.G is None or np.isnan(iter.G):
        #         raise EquilibError(f"Calculation failed {iter.G}")

        dict_of_lists = defaultdict(list)
        
        # --- STEP 1: DISCOVERY ---
        all_comp_keys = set()
        all_phase_keys = set()
        
        # Call the reusable discovery helper for 'xi'
        Endmembers_Xmole = self._discover_phase_sub_keys(iterations, "Endmembers_Xmole")
        Endmembers_Xmass = self._discover_phase_sub_keys(iterations, "Endmembers_Xmass")
        Elements_Xmole = self._discover_phase_sub_keys(iterations, "Elements_Xmole")
        Elements_Xmass = self._discover_phase_sub_keys(iterations, "Elements_Xmass")
        
        # --- (Future) Add discovery for other attributes here ---
        # all_ai_components = self._discover_phase_sub_keys(iterations, "ai")
        # all_mui_components = self._discover_phase_sub_keys(iterations, "mui")

        for iter_data in iterations:
            all_comp_keys.update([x.split(' ')[0] for x in list(iter_data.N_x.keys())])
            all_phase_keys.update(iter_data.Phases.keys())
        
        # --- STEP 2: POPULATION ---
        for iter_data in iterations:
            # (T, P, G, H, S, Cp, StablePhases*, N keys are all the same)
            dict_of_lists['T [K]'].append(iter_data.T)
            dict_of_lists['P [atm]'].append(iter_data.P)
            for key in all_comp_keys:
                dict_of_lists[f'{key} [sp-mol]'].append(iter_data.N_x.get(f'{key} [sp-mol]', 0.0))
            for key in all_comp_keys:
                dict_of_lists[f'{key} [g]'].append(iter_data.N_w.get(f'{key} [g]', 0.0))
            dict_of_lists['G [J]'].append(iter_data.G)
            dict_of_lists['H [J]'].append(iter_data.H)
            dict_of_lists['S [J/K]'].append(iter_data.S)
            dict_of_lists['Cp [J/K]'].append(iter_data.Cp)
            
            dict_of_lists['StablePhaseNames'].append(str(iter_data.StablePhases.get('Name', np.array([]))))
            dict_of_lists['StablePhaseIDs'].append(str(iter_data.StablePhases.get('ID', np.array([]))))
            dict_of_lists['StablePhaseAmountMole [sp-mol]'].append(str(iter_data.StablePhases.get('Amount_mole', np.array([]))))
            dict_of_lists['StablePhaseAmountMass [g]'].append(str(iter_data.StablePhases.get('Amount_mass', np.array([]))))

            
                
            for key in all_phase_keys: # e.g., key = 'LIQUID'
                phase = iter_data.Phases.get(key)
                
                # Handle Amount and Stability
                if phase:
                    dict_of_lists[f'{key}_Amount_mole [sp-mol]'].append(phase.Amount_mole)
                    dict_of_lists[f'{key}_Amount_mass [g]'].append(phase.Amount_mass)
                    dict_of_lists[f'{key}_Stability'].append(phase.Stability)
                else:
                    dict_of_lists[f'{key}_Amount_mole [sp-mol]'].append(0.0)
                    dict_of_lists[f'{key}_Amount_mass [g]'].append(0.0)
                    dict_of_lists[f'{key}_Stability'].append(0.0)

                # --- Call the reusable population helper for 'xi' ---
                self._populate_flattened_phase_attribute(dict_of_lists,
                                                         key, phase,
                                                         Endmembers_Xmole[key],
                                                         "Endmembers_Xmole", "_Endmembers_Xmole")
                self._populate_flattened_phase_attribute(dict_of_lists,
                                                         key, phase,
                                                         Endmembers_Xmass[key],
                                                         "Endmembers_Xmass", "_Endmembers_Xmass")
                self._populate_flattened_phase_attribute(dict_of_lists,
                                                         key, phase,
                                                         Elements_Xmole[key],
                                                         "Elements_Xmole", "_Elements_Xmole")
                self._populate_flattened_phase_attribute(dict_of_lists,
                                                         key, phase,
                                                         Elements_Xmass[key],
                                                         "Elements_Xmass", "_Elements_Xmass")
                
                # --- (Future) Add calls for other attributes here ---
                # self._populate_flattened_phase_attribute(dict_of_lists,
                #                                          key, phase,
                #                                          all_ai_components[key],
                #                                          "ai", "_activity")
                
                # self._populate_flattened_phase_attribute(dict_of_lists,
                #                                          key, phase,
                #                                          all_mui_components[key],
                #                                          "mui", "_mui")

        # --- (End of main loops) ---
        
        final_dict = dict(dict_of_lists)

        # Un-list the values if it was a single point calculation
        if is_single:
            for key, value_list in final_dict.items():
                final_dict[key] = value_list[0] if value_list else None
                    
        return final_dict

class ScheilDataRow(BaseModel):
    """
    Pydantic model for a single step in a Scheil simulation.
    This represents one "row" of Scheil-specific data.
    """
    T: float
    P: float
    G: Optional[float]
    H: Optional[float]
    S: Optional[float]
    Cp: Optional[float]
    fl_mole: float = Field(ge=0.0)
    fs_mole: float = Field(ge=0.0)
    fl_mass: float = Field(ge=0.0)
    fs_mass: float = Field(ge=0.0)
    label: str
    
    # This stores the CUMULATIVE phase amounts *at this step*
    cumulative_phases_mole: Dict[str, float] = Field(default_factory=dict)
    cumulative_phases_mass: Dict[str, float] = Field(default_factory=dict)

    model_config = ConfigDict(arbitrary_types_allowed=True)

class ResultScheil:
    """
    Refactored Scheil-Gulliver cooling result object.
    
    This class now mimics the dual-state (single vs. list) behavior
    of the main 'Result' class for full API consistency.
    """
    
    def __init__(self):
        # --- Core Attributes ---
        # self.EquilibResult holds the detailed SinglePointResult for each step
        self.EquilibResult = Result()
        
        # This holds the Scheil-specific data (T, fl, fs, cumulative phases)
        self.data: Union[None, ScheilDataRow, List[ScheilDataRow]] = None
        
        # --- System Properties (set once) ---
        self.N_x: Dict[str, float] = {}
        self.N_w: Dict[str, float] = {}
        
        # Stores final micro-constituent amounts
        self.ScheilConstituents: Dict[str, float] = {}

        # --- Initial Step (called from scheil_cooling) ---
        try:
            # 1. Capture the initial equilibrium state from Fortran
            self.EquilibResult.append_output()
            
            # 2. Get the data for this first step
            # At this point, EquilibResult properties return single values
            self.N_x = self.EquilibResult.N_x
            self.N_w = self.EquilibResult.N_w
            
            T = self.EquilibResult.T
            P = self.EquilibResult.P
            G = self.EquilibResult.G
            H = self.EquilibResult.H
            S = self.EquilibResult.S
            Cp = self.EquilibResult.Cp

            # 3. Process initial Scheil phases
            # EquilibResult.StablePhases is a single Dict[str, np.ndarray]
            CurrentPhases_Names = self.EquilibResult.StablePhases.get('Name', np.array([]))
            CurrentPhases_Amts_mole = self.EquilibResult.StablePhases.get('Amount_mole', np.array([]))
            CurrentPhases_Amts_mass = self.EquilibResult.StablePhases.get('Amount_mass', np.array([]))
            
            current_fl_mole = 0.0
            current_fs_mole = 0.0
            current_fl_mass = 0.0
            current_fs_mass = 0.0
            phase_labels = []
            cumulative_phases_mole = {}
            cumulative_phases_mass = {}
            
            for i, name in enumerate(CurrentPhases_Names):
                if name == ' ' or name == 'nan': 
                    continue
                
                amount_mole = np.round(CurrentPhases_Amts_mole[i],10)
                amount_mass = np.round(CurrentPhases_Amts_mass[i],10)
                cumulative_phases_mole[f'{name}_mole'] = amount_mole # Initial amount
                cumulative_phases_mass[f'{name}_mass'] = amount_mass # Initial amount
                phase_labels.append(name)
                
                if 'LIQ' in name.upper():
                    current_fl_mole += amount_mole
                    current_fl_mass += amount_mass
                else:
                    current_fs_mole += amount_mole
                    current_fs_mass += amount_mass
            
            fl_mole = current_fl_mole/(current_fl_mole+current_fs_mole) if (current_fl_mole+current_fs_mole)>0 else 0.0
            fs_mole = 1.0 - fl_mole

            fl_mass = current_fl_mass/(current_fl_mass+current_fs_mass) if (current_fl_mass+current_fs_mass)>0 else 0.0
            fs_mass = 1.0 - fl_mass
            
            # 4. Update PhaseLabel
            if not phase_labels:
                raise EquilibError('No phase appears stable during Scheil simulation')
            label = '+'.join(sorted(phase_labels))
            
            # 5. Create and store the first ScheilDataRow
            self.data = ScheilDataRow(
                T=T, P=P, G=G, H=H, S=S, Cp=Cp,
                fl_mole=fl_mole,
                fs_mole=fs_mole,
                fl_mass=fl_mass,
                fs_mass=fs_mass,
                label=label,
                cumulative_phases_mole = cumulative_phases_mole,
                cumulative_phases_mass = cumulative_phases_mass
            )
            
        except Exception as e:
            print(f"Error during ResultScheil initialization: {e}")
            raise
        finally:
            fort.resetthermo()
        
    def append_output(self):
        """
        Appends the next Scheil step's results from the current Fortran state.
        Handles the transition from a single data point to a list.
        """
        try:
            # 1. Capture the new equilibrium state
            self.EquilibResult.append_output()
            
            # 2. Get data for this new step
            # Now, EquilibResult properties return LISTS
            last_equilib_step = self.EquilibResult.data[-1] # Get the new SinglePointResult
            
            T = last_equilib_step.T
            P = last_equilib_step.P
            G = last_equilib_step.G
            H = last_equilib_step.H
            S = last_equilib_step.S
            Cp = last_equilib_step.Cp

            # 3. Get the *previous* Scheil step's cumulative phases
            if isinstance(self.data, ScheilDataRow):
                last_scheil_row = self.data
            elif isinstance(self.data, list):
                last_scheil_row = self.data[-1]
            else:
                raise EquilibError("Cannot append_output, ResultScheil is not initialized.")
            
            last_cumulative_phases_mole = last_scheil_row.cumulative_phases_mole
            last_cumulative_phases_mass = last_scheil_row.cumulative_phases_mass

            # 4. Get current stable phases from the *last* equilib step
            CurrentPhases_Names = last_equilib_step.StablePhases.get('Name', np.array([]))
            CurrentPhases_Amts_mole = last_equilib_step.StablePhases.get('Amount_mole', np.array([]))
            CurrentPhases_Amts_mass = last_equilib_step.StablePhases.get('Amount_mass', np.array([]))

            current_fl_mole = 0.0
            current_fs_mole = 0.0
            current_fl_mass = 0.0
            current_fs_mass = 0.0

            phase_labels = []
            new_cumulative_phases_mole = {}
            new_cumulative_phases_mass = {}

            # 5. Update ScheilPhases (cumulative)
            for i, name in enumerate(CurrentPhases_Names):
                if name == ' ' or name == 'nan': 
                    continue
                
                amoun_mole = float(CurrentPhases_Amts_mole[i])
                amoun_mass = float(CurrentPhases_Amts_mass[i])
                phase_labels.append(name)
                
                if 'LIQ' in name.upper():
                    new_cumulative_phases_mole[f'{name}_mole'] = amoun_mole
                    new_cumulative_phases_mass[f'{name}_mass'] = amoun_mass
                    current_fl_mole += amoun_mole
                    current_fl_mass += amoun_mass
                else:
                    # Solid phase: new_cumulative = last_cumulative + new_solid
                    
                    last_cumulative_amount_mole = last_cumulative_phases_mole.get(f'{name}_mole', 0.0)
                    last_cumulative_amount_mass = last_cumulative_phases_mass.get(f'{name}_mass', 0.0)
                    current_fs_mole += amoun_mole +last_cumulative_amount_mole
                    current_fs_mass += amoun_mass +last_cumulative_amount_mass
                    new_cumulative_phases_mole[f'{name}_mole'] = last_cumulative_amount_mole + amoun_mole
                    new_cumulative_phases_mass[f'{name}_mass'] = last_cumulative_amount_mass + amoun_mass
            #Add last cumulative amount for phases that disappeared
            for name_mole, last_amount in last_cumulative_phases_mole.items():
                name = name_mole.split('_mole')[0]

                if name not in CurrentPhases_Names:
                    last_cumulative_amount_mole = last_cumulative_phases_mole.get(f'{name}_mole', 0.0)
                    last_cumulative_amount_mass = last_cumulative_phases_mass.get(f'{name}_mass', 0.0)
                    current_fs_mole += last_cumulative_amount_mole
                    current_fs_mass += last_cumulative_amount_mass
                    
                    if 'LIQ' in name.upper():
                        new_cumulative_phases_mole[f'{name}_mole'] = 0.0
                        new_cumulative_phases_mass[f'{name}_mass'] = 0.0
                    else:
                        new_cumulative_phases_mole[f'{name}_mole'] = last_cumulative_amount_mole
                        new_cumulative_phases_mass[f'{name}_mass'] = last_cumulative_amount_mass


            fl_mole = current_fl_mole/(current_fl_mole+current_fs_mole) if (current_fl_mole+current_fs_mole)>0 else 0.0
            fs_mole = 1.0 - fl_mole

            fl_mass = current_fl_mass/(current_fl_mass+current_fs_mass) if (current_fl_mass+current_fs_mass)>0 else 0.0
            fs_mass = 1.0 - fl_mass

            
            # 6. Update PhaseLabel
            if not phase_labels:
                raise EquilibError('No phase appears stable during Scheil simulation')
            label = '+'.join(sorted(phase_labels))
            
            # 7. Add padding for phases that disappeared
            # (Carry over last cumulative amount)
            for name, last_amount in last_cumulative_phases_mole.items():
                name = name_mole.split('_mole')[0]
                if name_mole not in new_cumulative_phases_mole:
                    new_cumulative_phases_mole[f'{name}_mole'] = last_amount
                    new_cumulative_phases_mass[f'{name}_mass'] = last_cumulative_phases_mass[f'{name}_mass']
                    if 'LIQ' in name.upper():
                        new_cumulative_phases_mole[f'{name}_mole'] = 0.0
                        new_cumulative_phases_mass[f'{name}_mass'] = 0.0
                        

            # 8. Create the new ScheilDataRow
            new_row = ScheilDataRow(
                T=T, P=P, G=G, H=H, S=S, Cp=Cp,
                fl_mole=fl_mole,
                fs_mole=fs_mole,
                fl_mass=fl_mass,
                fs_mass=fs_mass,
                label=label,
                cumulative_phases_mole=new_cumulative_phases_mole,
                cumulative_phases_mass=new_cumulative_phases_mass
            )

            # 9. Handle the dual-state transition
            if isinstance(self.data, ScheilDataRow):
                self.data = [self.data, new_row] # Convert to list
            elif isinstance(self.data, list):
                self.data.append(new_row)
                    
        except Exception as e:
            print(f"Error during ResultScheil.append_output: {e}")
            raise
        finally:
            fort.resetthermo()

    # --- Property Accessors (Dual-State) ---

    @property
    def T(self) -> Union[float, List[float]]:
        if isinstance(self.data, ScheilDataRow): return self.data.T
        if isinstance(self.data, list): return [iter.T for iter in self.data]
        return None
        
    @property
    def P(self) -> Union[float, List[float]]:
        if isinstance(self.data, ScheilDataRow): return self.data.P
        if isinstance(self.data, list): return [iter.P for iter in self.data]
        return None

    @property
    def G(self) -> Union[Optional[float], List[Optional[float]]]:
        if isinstance(self.data, ScheilDataRow): return self.data.G
        if isinstance(self.data, list): return [iter.G for iter in self.data]
        return None

    @property
    def H(self) -> Union[Optional[float], List[Optional[float]]]:
        if isinstance(self.data, ScheilDataRow): return self.data.H
        if isinstance(self.data, list): return [iter.H for iter in self.data]
        return None

    @property
    def S(self) -> Union[Optional[float], List[Optional[float]]]:
        if isinstance(self.data, ScheilDataRow): return self.data.S
        if isinstance(self.data, list): return [iter.S for iter in self.data]
        return None
    @property
    def Cp(self) -> Union[Optional[float], List[Optional[float]]]:
        if isinstance(self.data, ScheilDataRow): return self.data.Cp
        if isinstance(self.data, list): return [iter.Cp for iter in self.data]
        return None
    @property
    def fl_mole(self) -> Union[float, List[float]]:
        if isinstance(self.data, ScheilDataRow): return self.data.fl_mole
        if isinstance(self.data, list): return [iter.fl_mole for iter in self.data]
        return None

    @property
    def fs_mole(self) -> Union[float, List[float]]:
        if isinstance(self.data, ScheilDataRow): return self.data.fs_mole
        if isinstance(self.data, list): return [iter.fs_mole for iter in self.data]
        return None
    
    @property
    def fl_mass(self) -> Union[float, List[float]]:
        if isinstance(self.data, ScheilDataRow): return self.data.fl_mass
        if isinstance(self.data, list): return [iter.fl_mass for iter in self.data]
        return None

    @property
    def fs_mass(self) -> Union[float, List[float]]:
        if isinstance(self.data, ScheilDataRow): return self.data.fs_mass
        if isinstance(self.data, list): return [iter.fs_mass for iter in self.data]

    @property
    def PhaseLabel(self) -> Union[str, List[str]]:
        if isinstance(self.data, ScheilDataRow): return self.data.label
        if isinstance(self.data, list): return [iter.label for iter in self.data]
        return []

    @property
    def ScheilPhases_mole(self) -> Dict[str, List[float]]:
        """
        Dynamically builds the dict-of-lists for cumulative phase amounts
        to maintain compatibility with the 'scheil_cooling' function.
        """
        if self.data is None: return {}

        # Normalize data to a list for iteration
        if isinstance(self.data, ScheilDataRow):
            iterations = [self.data]
        else: # is a list
            iterations = self.data
            
        all_phase_keys = set()
        for row in iterations:
            all_phase_keys.update(row.cumulative_phases_mole.keys())
        
        dict_of_lists = defaultdict(list)
        for row in iterations:
            for key in all_phase_keys:
                amount = row.cumulative_phases_mole.get(key, 0.0)
                dict_of_lists[key].append(amount)
        return dict(dict_of_lists)
    
    @property
    def ScheilPhases_mass(self) -> Dict[str, List[float]]:
        """
        Dynamically builds the dict-of-lists for cumulative phase amounts
        to maintain compatibility with the 'scheil_cooling' function.
        """
        if self.data is None: return {}

        # Normalize data to a list for iteration
        if isinstance(self.data, ScheilDataRow):
            iterations = [self.data]
        else: # is a list
            iterations = self.data
            
        all_phase_keys = set()
        for row in iterations:
            all_phase_keys.update(row.cumulative_phases_mass.keys())
        
        dict_of_lists = defaultdict(list)
        for row in iterations:
            for key in all_phase_keys:
                amount = row.cumulative_phases_mass.get(key, 0.0)
                dict_of_lists[key].append(amount)
        return dict(dict_of_lists)

    # --- Methods ---

    def update_scheilconstituents(self):
        """
        Post-processing calculation for micro-constituents.
        Refactored to work with the internal 'data' attribute.
        """
        k = 0
        self.ScheilConstituents = {} # Reset
        
        # Normalize data to lists for safe iteration
        if self.data is None: return None
        if isinstance(self.data, ScheilDataRow):
            iterations = [self.data]
        else:
            iterations = self.data

        if len(iterations) < 2:
            return None # Not enough data
        
        T_list = [iter.T for iter in iterations]
        label_list = [iter.label for iter in iterations]
        fs_mole_list = [iter.fs_mole for iter in iterations]
        fs_mass_list = [iter.fs_mass for iter in iterations]

        #Update component, solidification range
        self.ScheilConstituents['T_Delta [K]'] = T_list[0] - T_list[-1]
        for component in var.cComponentNameSys:
            self.ScheilConstituents[f'{component} [sp-mol]'] = self.N_x.get(f'{component} [sp-mol]', 0.0)
        for component in var.cComponentNameSys:
            self.ScheilConstituents[f'{component} [g]'] = self.N_w.get(f'{component} [g]', 0.0)
            
        cphase = ''
        
        for i, T in enumerate(T_list):
            label = label_list[i]
            fs_mole = fs_mole_list[i]
            fs_mass = fs_mass_list[i]
            
            constituents = label.split('+')
            constituents = [x for x in constituents if 'LIQ' not in x.upper() and x != '']
            
            if i == 0:
                cphase = '+'.join(constituents)
                k = 0
                continue
            
            if i == len(T_list) - 1:
                if k == 0: 
                    if cphase:
                        self.ScheilConstituents[f'{cphase} [sp-mol]'] = fs_mole_list[k]
                        self.ScheilConstituents[f'{cphase} [g]']   = fs_mass_list[k]
                    return None
                    
                fphase = '+'.join(constituents)
                prev_fs_mole_k = fs_mole_list[k]
                prev_fs_mass_k = fs_mass_list[k]
                
                if cphase != fphase:
                    self.ScheilConstituents[f'{fphase} [sp-mol]'] = self.ScheilConstituents.get(f'{fphase} [sp-mol]', 0.0) + (fs_mole - prev_fs_mole_k)
                    self.ScheilConstituents[f'{fphase} [g]']   = self.ScheilConstituents.get(f'{fphase} [g]', 0.0)   + (fs_mass - prev_fs_mass_k)
                else:
                    self.ScheilConstituents[f'{cphase} [sp-mol]'] = self.ScheilConstituents.get(f'{cphase} [sp-mol]', 0.0) + (fs_mole - prev_fs_mole_k)
                    self.ScheilConstituents[f'{cphase} [g]']   = self.ScheilConstituents.get(f'{cphase} [g]', 0.0)   + (fs_mass - prev_fs_mass_k)
                return None
            
            next_label = label_list[i+1]
            
            if label != next_label and k == 0:
                cphase = '+'.join(constituents)
                # print('Start',i,label,next_label,k,prev_fs,fs)
                self.ScheilConstituents[f'{cphase} [sp-mol]'] = (fs_mole)
                self.ScheilConstituents[f'{cphase} [g]']   = (fs_mass)
                k = i
            elif label != next_label and k > 0:
                prev_fs_mole = fs_mole_list[k]
                prev_fs_mass = fs_mass_list[k]
                # print('Middle',i,label,next_label,k,prev_fs,fs)
                cphase = '+'.join(constituents)
                self.ScheilConstituents[f'{cphase} [sp-mol]'] = self.ScheilConstituents.get(f'{cphase} [sp-mol]', 0.0) + (fs_mole-prev_fs_mole)
                self.ScheilConstituents[f'{cphase} [g]'] = self.ScheilConstituents.get(f'{cphase} [g]', 0.0) + (fs_mass-prev_fs_mass)
                k = i
            
        return None
    
    def to_dict(self) -> Dict[str, Any]:
        """
        Exports all results to a flattened dictionary, handling
        both single-point and multi-point (list) states.
        """
        if self.data is None: return {}

        # 1. Get the detailed, flattened dict from the EquilibResult
        df = dict()
        df_eq = self.EquilibResult.to_dict()
        

        
        # 2. Normalize Scheil data to lists for processing
        if isinstance(self.data, ScheilDataRow):
            iterations = [self.data]
            is_single = True
        else: # is a list
            iterations = self.data
            is_single = False

        # for iter in iterations:
        #     if iter.G is None or np.isnan(iter.G):
        #         raise EquilibError(f"Calculation failed {fort.modulegemsolver.dgemfunctionnorm}")


        # 3. Add/Overwrite with Scheil-specific lists
        df['T [K]'] = [iter.T for iter in iterations]
        df['P [atm]'] = [iter.P for iter in iterations]
        for component in var.cComponentNameSys:
            df[f'{component} [sp-mol]'] = [self.N_x.get(f'{component} [sp-mol]', 0.0) for _ in iterations]
        for component in var.cComponentNameSys:
            df[f'{component} [g]'] = [self.N_w.get(f'{component} [g]', 0.0) for _ in iterations]
        df['Label'] = [iter.label for iter in iterations]
        df['fl [sp-mol/sp-mol]'] = [iter.fl_mole for iter in iterations]
        df['fs [sp-mol/sp-mol]'] = [iter.fs_mole for iter in iterations]
        df['fl [g/g]'] = [iter.fl_mass for iter in iterations]
        df['fs [g/g]'] = [iter.fs_mass for iter in iterations]
        df['G [J]'] = [iter.G for iter in iterations]
        df['H [J]'] = [iter.H for iter in iterations]
        df['S [J/K]'] = [iter.S for iter in iterations]
        df['Cp [J/K]'] = [iter.Cp for iter in iterations]
        
        
        # 4. Add the cumulative ScheilPhases (overwriting amounts)
        cumulative_phases_dict_mole = self.ScheilPhases_mole
        for phase_name_mole, amounts_list in cumulative_phases_dict_mole.items():
            phase_name = phase_name_mole.split('_mole')[0]
            df[f'{phase_name}_Amount_mole [sp-mol]'] = amounts_list
            for key, value in df_eq.items():
                if key.startswith(f'{phase_name}_Endmembers_Xmole_'):  # Check if the key matches the current phase
                    df[key] = value
        
        cumulative_phases_dict_mass = self.ScheilPhases_mass
        for phase_name_mass, amounts_list in cumulative_phases_dict_mass.items():
            phase_name = phase_name_mass.split('_mass')[0]
            df[f'{phase_name}_Amount_mass [g]'] = amounts_list
            for key, value in df_eq.items():
                if key.startswith(f'{phase_name}_Endmembers_Xmass_'):  # Check if the key matches the current phase
                    df[key] = value
            
        # # 5. Un-list all values if it was a single point
        # if is_single:
        #     for key, value_list in df.items():
        #         if isinstance(value_list, list) and len(value_list) == 1:
        #             df[key] = value_list[0]
        #         elif not isinstance(value_list, (dict, str)):
        #             # This handles keys from EquilibResult.to_dict
        #             # that were *already* single values
        #             pass 
        
        return df