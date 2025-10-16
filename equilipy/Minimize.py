import numpy as np
import equilipy.equilifort as fort
import equilipy.variables as var
from .Errors import *

def minimize():
    info=fort.modulethermoio.infothermo
    
    # Check the input variables:
    try:
        fort.checkthermoinput()
    except Exception as e: 
        raise EquilibError(f'Equilifort checkthermoinput failed: infothermo={info}, Error: {e}')
    
    #Initialize equilifort:
    try:
        fort.initthermo()
    except Exception as e:  
        raise EquilibError(f'Equilifort initthermo failed: infothermo={info}, Error: {e}')
    
    #Check if the input condition is consistent with database:
    try:
        fort.checksystem()
    except Exception as e:  
        raise EquilibError(f'Equilifort checksystem failed: infothermo={info}, Error: {e}')
    
    
    #Compute thermodynamic data:
    try:
        fort.compthermodata()
    except Exception as e:  
        raise EquilibError(f'Equilifort compthermodata failed: infothermo={info}, Error: {e}')
    
    # Check the thermodynamic data:
    try:
        fort.checkthermodata()
    except Exception as e:  
        raise EquilibError(f'Equilifort checkthermodata failed: infothermo={info}, Error: {e}')
    
    #Estimate the equilibrium phase assemblage and other important properties
    #using the Leveling algorithm:
    try:
        fort.gemsolver()
    except Exception as e:  
        raise EquilibError(f'Equilifort gemsolvernew failed: infothermo={info}, Error: {e}')
    
    # Perform post-processing calculations of results:
    try:
        fort.postprocess()
    except Exception as e:  
        raise EquilibError(f'Equilifort postprocess failed: infothermo={info}, Error: {e}')
    
    if fort.modulegemsolver.dgemfunctionnorm>1E-3:
        raise EquilibError(f'Equilibrium calculation failed: dGEMFunctionNorm={fort.modulegemsolver.dgemfunctionnorm}')

    return None


def minimize_single_phase():
    
    # info=fort.modulethermoio.infothermo
    info = fort.modulethermoio.dtemperaturediff
    # Check the input variables:
    try:
        fort.checkthermoinput()
    except Exception as e: 
        raise EquilibError(f'Equilifort checkthermoinput failed: infothermo={info}, Error: {e}')
    
    #Initialize equilifort:
    try:
        fort.initthermo()
    except Exception as e:  
        raise EquilibError(f'Equilifort initthermo failed: infothermo={info}, Error: {e}')
    
    #Check if the input condition is consistent with database:
    try:
        fort.checksystem()
    except Exception as e:  
        raise EquilibError(f'Equilifort checksystem failed: infothermo={info}, Error: {e}')
        
    #Compute thermodynamic data:
    try:
        fort.compsinglephaseprop()
    except Exception as e:  
        raise EquilibError(f'Equilifort compthermodata failed: infothermo={info}, Error: {e}')
    
    if fort.modulegemsolver.dgemfunctionnorm>1E-3:
        raise EquilibError(f'Equilibrium calculation failed: dGEMFunctionNorm={fort.modulegemsolver.dgemfunctionnorm}')
    return None