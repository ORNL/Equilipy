#!/usr/bin/env python3
import equilipy.equilifort as fort
from typing import Dict, List, Optional, Union
import numpy as np
from tqdm import tqdm
from .PostProcess import ResultScheil
from .EquilibSingle import _equilib_single
from .Errors import EquilibError
from .FindTransition import find_transitions

def scheil_cooling(
    LiquidPhaseName: str,
    Database: Dict,
    Condition: Dict[str, Union[float, str]],
    dT: float = 5.0,
    UnitIn: List[str] = ['K', 'atm', 'moles'],
    UnitOut: List[str] = None,
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
    if UnitOut==None:
        UnitOut = UnitIn.copy()

    # Start from liquidus
    if StartFromLiquidus:
        Tmax = current_condition['T']
        Tmin = current_condition['T']*0.1
        Ts=find_transitions(Database, current_condition, Tmax,Tmin,UnitIn=UnitIn)
        current_condition['T'] = Ts[0] + 0.1

    # --- Initial Equilibrium Calculation ---
    # Perform a single equilibrium calculation to establish the starting point
    _equilib_single(Database, current_condition, UnitIn=UnitIn,UnitOut=UnitOut, ListOfPhases=ListOfPhases)
    assemblage_old = fort.modulethermo.iassemblage.copy()
    res = ResultScheil()
    StablePhaseNames = list(res.ScheilPhases.keys())

    # Check if the liquid phase is stable at the start
    if LiquidPhaseName not in StablePhaseNames:
        raise EquilibError(f"The specified liquid phase '{LiquidPhaseName}' is not stable at the initial conditions.")
    
    IterMax= int(Condition['T']//dT)
    
    # --- Main Cooling Loop ---
    # The for loop iterates up to IterMax times, with a progress bar.
    for i in tqdm(range(IterMax), desc="Scheil Cooling", colour='#8060ff', ascii="░▒▓", bar_format="{l_bar}{bar:50}{r_bar}"):
        # Get the composition and amount of the liquid from the previous step based on the OutputUnit
        liquid_amount = res.ScheilPhases[LiquidPhaseName][-1]
        if i==0:
            current_condition['T'] = res.T #temperature baed on output unit
            liquid_composition = res.EquilibResult.Phases[LiquidPhaseName].xi
            liquid_components = res.EquilibResult.Phases[LiquidPhaseName].Endmembers
        else:
            liquid_composition = res.EquilibResult.Phases[LiquidPhaseName].xi[-1]
            liquid_components = res.EquilibResult.Phases[LiquidPhaseName].Endmembers[-1]

        # --- Termination Conditions ---
        if liquid_amount < 1e-5:
            # Solidification is considered complete
            break
        
        current_temp = float(current_condition['T'])
        if current_temp <= 0 and 'K' in UnitIn[0]:
            raise EquilibError('Scheil-Gulliver solidification results in melting point below 0 K. Check phase selections.')

        # --- Update Conditions for Next Step ---
        # Decrease temperature
        current_condition['T'] = current_temp - dT
        
        # Update composition based on remaining liquid
        new_ni = np.array(liquid_composition) * liquid_amount
        for i, el in enumerate(liquid_components):
            current_condition[el] = new_ni[i]

        # --- Perform Equilibrium Calculation for the Current Step ---
        _equilib_single(Database, current_condition, UnitIn=UnitOut,UnitOut=UnitOut, ListOfPhases=ListOfPhases)
        
        # --- Before updating the result, if there is any transition occurs, calculate the transition point
        assemblage_new = fort.modulethermo.iassemblage.copy()
        if set(assemblage_new)!=set(assemblage_old):
            # # Transition occured, Reset variables
            Tmax = float(current_temp)
            Tmin = float(current_condition['T'])
            current_condition['T'] = Tmax
            Ts=find_transitions(Database, current_condition, Tmax,Tmin,UnitIn=UnitOut)
            if abs(Tmax-Ts[0])<0.1: Ts = Ts[1:]

            for j, T_trans in enumerate(Ts):
                #Update transition point
                current_condition['T'] = T_trans
                _equilib_single(Database, current_condition, UnitIn=UnitOut,UnitOut=UnitOut, ListOfPhases=ListOfPhases)
                if j<len(Ts)-1:
                    res.append_output()
                    # Get the composition and amount of the liquid from the previous step
                    liquid_amount = res.ScheilPhases[LiquidPhaseName][-1]
                    liquid_composition = res.EquilibResult.Phases[LiquidPhaseName].xi[-1]
                    liquid_components = res.EquilibResult.Phases[LiquidPhaseName].Endmembers[-1]
                    
                    # Update composition based on remaining liquid
                    new_ni = np.array(liquid_composition) * liquid_amount
                    for k, el in enumerate(liquid_components):
                        current_condition[el] = new_ni[k]
            current_condition['T'] =float(Tmin)
        assemblage_new = fort.modulethermo.iassemblage.copy()
        assemblage_old=assemblage_new.copy()
        res.append_output()
        
        # Check if liquid phase disappeared in the last step
        if LiquidPhaseName not in res.EquilibResult.StablePhases['Name'][-1]:
            break
    else:
        # This block runs only if the for loop completes without a 'break'
        print(f'Warning: Reached maximum iterations ({IterMax}) before solidification completed.')

    # Finalize results
    res.update_scheilconstituents()
    return res
