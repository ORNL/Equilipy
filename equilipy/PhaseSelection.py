import numpy as np
import equilipy.equilifort as fort
import equilipy.variables as var
from .Errors import InputConditionError
def phase_selection(ListOfPhases):
    '''----------------------------------------------------------------------------------------------------------------
    Description
    ===========
    This function revise variables to pass to equilifort module. Not that this function is active when system elements
    are defined.

    Revisions
    =========

     Date            Programmer      Description of change
     ----            ----------      ---------------------
     02/06/2024      S.Y. KWON       Original code


    Variables
    =========

    Input
    DBobject : Global variable object that includes parsed database
    ListOfPhases : A list of phases selected by the user

    Output
    Revised DBobject
    ----------------------------------------------------------------------------------------------------------------'''
    iPhasePS= np.ones(len(var.iPhaseCS),dtype=int)*(-1)
    iSolnPS = np.ones(var.nSolnPhasesSysCS,dtype=int)*(-1)
    
    # Check if ListOfPhases are list
    if not isinstance(ListOfPhases, list): raise InputConditionError('Make sure if the data type of ListOfPhases is list.')
    
    # Revise var.iSys2DBSoln, var.iSys2DBComp, iSys2DBSpecies, var.PhaseNameSys, var.cPhaseNameSys
    iSys2DBSoln = list([])
    iSys2DBComp = list([])
    iSys2DBSpecies = list([])
    cPhaseNameSys = list([])
    PhaseTypeIDSys = list([])
    NewPhaseList = list([])
    
	#1. Add each phase information into PS variables
    for phase in ListOfPhases:
        # Check the type of phase name
        if not isinstance(phase, str): raise InputConditionError('Make sure if all elements in ListOfPhases are string.')
        
        # For immiscible phases, the second and third phase must involve the first phase:
        if '#' in phase:
            fphase= phase.split('#')[0]
            if fphase not in ListOfPhases:
                raise InputConditionError(f'Immiscibile phase {phase} was selected without involving {fphase}.')
                
        #Get the type of the phase and corresponding id in database
        try: 
            p_type,p_id = var.PhaseNameSys[phase]
            NewPhaseList.append(phase)
        except KeyError:
            continue
        PhaseTypeIDSys.append([p_type,p_id])
        cPhaseNameSys.append(phase)
        
        if p_type=='soln':
            if p_id ==1: iFirst =0
            else:iFirst=var.nSpeciesPhaseCS[int(p_id-2)]
            iLast=var.nSpeciesPhaseCS[int(p_id-1)]
            iPhasePS[iFirst:iLast]=p_id
            iSolnPS[int(p_id-1)]=p_id
            iSys2DBSoln.append(p_id)
            for i in range(iFirst,iLast):
                if sum(var.dStoichSpeciesCS[i,var.iElementDBIndex])>0:
                    iSys2DBSpecies.append(i)
        elif p_type=='compd':
            iPhasePS[p_id]=0
            iSys2DBComp.append(p_id)
            iSys2DBSpecies.append(p_id)
        else:
            raise NameError(f'Error: PhaseSelection cannot identify type of {phase}')
    # print(iSys2DBSpecies)
	#Allocate python variables:
    var.iSys2DBSoln 	= iSys2DBSoln
    var.iSys2DBComp 	= iSys2DBComp
    var.iSys2DBSpecies  = iSys2DBSpecies
    
    var.cPhaseNameSys   = cPhaseNameSys
    var.PhaseNameSys=dict(zip(cPhaseNameSys,PhaseTypeIDSys))
    var.iPhaseCS = iPhasePS
	
    #Allocate to fortran variables
    fort.moduleparsecs.iphasecs= iPhasePS
    fort.modulethermo.isolnps= iSolnPS
    return NewPhaseList

