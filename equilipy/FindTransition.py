import numpy as np
from .EquilibSingle import _equilib_single
import equilipy.equilifort as fort
from typing import Dict, List, Tuple, Optional, Set

def find_transitions(
    Database: dict,
    NTP: dict,
    Tmax: float,
    Tmin: float,
    UnitIn: List[str] = ['K', 'atm', 'moles'],
    ListOfPhases: Optional[List[str]] = None,
    T_tol: float = 1E-1,
    max_depth: int = 15
) -> np.ndarray:
    """
    Finds all phase transition temperature intervals within a given temperature range using a recursive bisection method.

    This function includes two robustness features:
    1. A pre-check to ensure the initial temperature range brackets a phase change. If not, it slightly expands the range.
    2. A retry mechanism that lowers the minimum temperature if no transitions are found in the initial search.

    Parameters
    ----------
    Database : dict
        The database object.
    NTP : dict
        A dictionary with NTP conditions (e.g., {'P': 1, 'N': 1}). Temperature will be set internally.
    Tmax : float
        The upper bound of the temperature search range.
    Tmin : float
        The lower bound of the temperature search range.
    ListOfPhases : Optional[List[str]], optional
        A list of all possible phase names. If None, all phases in the database are considered. By default None.
    T_tol : float, optional
        The temperature tolerance for identifying a transition interval. When the search interval is smaller
        than this, the interval boundaries are recorded. By default 1E-1.
    max_depth : int, optional
        The maximum recursion depth to prevent infinite loops. By default 15.

    Returns
    -------
    np.ndarray
        A sorted NumPy array containing pairs of temperatures that bracket each transition.
        For example: [T1_low, T1_high, T2_low, T2_high, ...]. Returns an empty array if no transitions are found after all retries.
    """
    
    def get_stable_phases_set(T: float) -> Set[int]:
        """Helper function to get the set of stable phase IDs at a given temperature."""
        ntp_local = NTP.copy()
        ntp_local['T'] = T
        _equilib_single(Database, ntp_local, UnitIn=UnitIn,ListOfPhases=ListOfPhases,CalcHeatCapacity=False)
        assemblage=fort.modulethermo.iassemblage.copy()
        return set(assemblage)

    def _find_recursive(t_low: float, t_high: float, depth: int, transitions_list: list):
        """Recursively search for transitions in the interval [t_low, t_high]."""
        if depth > max_depth:
            print(f"Warning: Max recursion depth reached in interval [{t_low}, {t_high}].")
            return
            
        try:
            set_low = get_stable_phases_set(t_low)
            set_high = get_stable_phases_set(t_high)
        except Exception as e:
            print(f"Warning: Equilibrium calculation failed in range [{t_low}, {t_high}]. Error: {e}")
            return

        if set_low == set_high:
            return

        if (t_high - t_low) < T_tol:
            transitions_list.extend([t_low, t_high])
            return

        t_mid = (t_low + t_high) / 2
        _find_recursive(t_low, t_mid, depth + 1, transitions_list)
        _find_recursive(t_mid, t_high, depth + 1, transitions_list)

    # --- Main logic with retry mechanisms ---
    # Use local copies of Tmin/Tmax for the search
    current_Tmin = Tmin
    current_Tmax = Tmax
    
    for retry_attempt in range(5): # Allow up to 5 major retries
        
        # --- 1. Pre-check to ensure boundaries have different phase sets ---
        search_Tmin, search_Tmax = current_Tmin, current_Tmax
        for adjustment_attempt in range(10):
            set_min = get_stable_phases_set(search_Tmin)
            set_max = get_stable_phases_set(search_Tmax)
            
            if set_min != set_max:
                break  # Found a valid range with different phase sets
            
            # If sets are the same, slightly expand the range and re-check
            print(f"Info: Boundary phase sets are identical. Expanding range: Tmin={search_Tmin-0.1:.2f}, Tmax={search_Tmax+0.1:.2f}")
            search_Tmin -= 0.1
            search_Tmax += 0.1
        else: # This 'else' belongs to the 'for' loop
            print("Warning: Could not find different phase sets at boundaries after 10 adjustments. Proceeding with original range.")
            search_Tmin, search_Tmax = current_Tmin, current_Tmax

        # --- 2. Perform the recursive search ---
        transitions = []
        _find_recursive(search_Tmin, search_Tmax, 0, transitions)
        
        # --- 3. Check results and decide whether to retry ---
        if transitions:
            # Success: found transitions, return the result
            return np.unique(np.asarray(transitions))[::-1]
        
        # Failure: no transitions found, prepare for next retry
        print(f"Warning: No transitions found in range [{current_Tmin:.2f}, {current_Tmax:.2f}]. Retrying with lower Tmin.")
        current_Tmin -= 100.0 # Lower Tmin for the next attempt
        
    # If all retries fail, return an empty array
    print("Error: No transitions found after all retry attempts.")
    return np.array([])
