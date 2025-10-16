import numpy as np, sys
from numba import njit
import equilipy.variables as var
import equilipy.equilifort as fort


def system_check(cElementSys):

    '''----------------------------------------------------------------------------------------------------------------
    Description
    ===========
    This checks the validity of the component names given in a table. Errors occur when 1. A component has an invalid
    element name, and 2. A component is not part of the loaded database.


    Revisions
    =========

     Date            Programmer      Description of change
     ----            ----------      ---------------------
     11/19/2020      S.Y. KWON       Original code


    Variables
    =========

    Input
    NTP : pandas dataframe of T, P, element1, element2 ...

    Output
    res : an integer array of atomic number
    ----------------------------------------------------------------------------------------------------------------'''
    #Info of name, mass, number for all elements
    #Consider using python class for storing the information
    
    Elements=var.cPeriodicTable
    n=len(cElementSys)
    var.cComponentNameSys=cElementSys
    var.iElementSys=np.ones(n)
    var.iElementSysIndex=np.ones(n,dtype=int)
    var.iElementDBIndex=np.ones(n,dtype=int)
    Sym=[x.strip() for x in cElementSys]
    
    iElementPass= np.ones(var.nElementsCS,dtype=int)
    for i in range(n):
        try: var.iElementSys[i]=int(Elements[Sym[i]][0])
        except: print('Error: Invalid element name')
        if "{:<3}".format(cElementSys[i]) not in var.cElementNameCS:
            print('Error:Some input elements are not in the database')
            sys.exit()
    j=0
    
    for i in range(var.nElementsCS):
        if var.cElementNameCS[i].strip() in cElementSys:
            k=cElementSys.index(var.cElementNameCS[i].strip())
            var.iElementSysIndex[k]=j
            var.iElementDBIndex[k]=i
            iElementPass[i]=0
            j+=1
    var.iElementSys=var.iElementSys[var.iElementSysIndex]
    
    # #Identify species that can pass for the system    
    # lSpeciesPass = _get_species_pass(iElementPass)
    # if sum(lSpeciesPass) == 0: 
    #     lSpeciesPass = _get_endmembers(iElementPass)
    #     fort.moduleparsecs.lendmembers2species= True
    # else:
    #     fort.moduleparsecs.lendmembers2species= False

    lSpeciesPass = _get_endmembers(iElementPass)
    fort.moduleparsecs.lendmembers2species= True
    lCompPass    = np.logical_and(var.iPhaseCS==0,lSpeciesPass)
    
    # Get Species index and name for the system
    IndexSpecies=np.arange(len(lSpeciesPass))
    var.iSys2DBSpecies = IndexSpecies[lSpeciesPass]  
    var.cSpeciesNameSys= np.array(var.cSpeciesNameCS)[var.iSys2DBSpecies]
    
    # Get Solution index and name for the system 
    IndexSoln=var.iPhaseCS[lSpeciesPass]
    var.iSys2DBSoln=np.unique(IndexSoln[IndexSoln>0])
    var.iSys2DBComp = IndexSpecies[lCompPass]
    
    var.cPhaseNameSys= np.array(var.cPhaseNames)[var.iSys2DBSoln-1]
    var.cPhaseNameSys=np.append(var.cPhaseNameSys,np.array(var.cSpeciesNameCS)[lCompPass])
    
    # Add # for immiscible phases
    if len(dict.fromkeys(var.cPhaseNameSys))!=len(var.cPhaseNameSys):
        item_count = {}
        result_list = []

        for item in var.cPhaseNameSys:
            if item in item_count:
                item_count[item] += 1
                result_list.append(f'{item}#{item_count[item]}')
            else:
                item_count[item] = 1
                result_list.append(item)
        
        var.cPhaseNameSys=result_list
    var.cPhaseNameSys=[x.strip() for x in var.cPhaseNameSys]
    
    # Additonal variables for postprocess    
    var.iSys2DBSolnDefault =var.iSys2DBSoln
    var.cPhaseNameSysDefault=var.cPhaseNameSys
    var.iSys2DBCompDefault = var.iSys2DBComp
    var.iSys2DBSpeciesDefault= var.iSys2DBSpecies
 
    return None

@njit
def _get_species_pass(iElementPass):
    #This function identifies which species to be considered for a given system
    
    #Get elements not considered in the system
    iElementNotPass= iElementPass.nonzero()[0]
    
    #First guess the species that need to be considered
    res=np.sum(var.dStoichSpeciesCS[:,iElementNotPass],axis=1)<1E-15
    
    # Refine the species based on solution phases
    # Note that solution phase should have more than one endmember, otherwise remove the endmembers
    
    #Get solution index
    IndexSoln=var.iPhaseCS[res]
    n=len(res)
    m = len(IndexSoln)
    k=0
    
    for i in range(n):
        # Among the species that need to be passed in database
        if res[i]:
            if (k==0):
                # Condition for the first endmembers
                
                if (IndexSoln[k]>0) and (IndexSoln[k]==IndexSoln[k+1]):
                    #If more than two endmembers are involved add the phase
                    res[i] = True
                else:
                    res[i] = False
                k+=1
            elif k==m-1:
                # Condition for the last endmembers
                if ((IndexSoln[k]>0) and (IndexSoln[k]==IndexSoln[k-1])):
                    # If more than two endmembers are involved add the phase
                    res[i] = True
                else:
                    res[i] = False
                k+=1
            elif k<m-1:
                # Condition for endmembers that are inbetween
                if (IndexSoln[k]!=IndexSoln[k-1]) and (IndexSoln[k]!=IndexSoln[k+1]) and (IndexSoln[k]>0):
                    res[i]=False
                elif IndexSoln[k]==-1:
                    res[i]=False
                k+=1

    return res

@njit
def _get_endmembers(iElementPass):
    #This function identifies which species to be considered for a given system
    
    #Get elements not considered in the system
    iElementNotPass= iElementPass.nonzero()[0]
    
    #First guess the species that need to be considered
    res=np.sum(var.dStoichSpeciesCS[:,iElementNotPass],axis=1)<1E-15
    
    # Refine the species based on solution phases
    # Note that solution phase should have more than one endmember, otherwise remove the endmembers
    
    #Get solution index
    IndexSoln=var.iPhaseCS[res]
    n=len(res)
    m = len(IndexSoln)
    k=0
    
    for i in range(n):
        # Among the species that need to be passed in database
        if res[i]:
            if (IndexSoln[k]>0):
                #If more than two endmembers are involved add the phase
                res[i] = True
            elif (IndexSoln[k]<0):
                res[i] = False
            k+=1

    return res