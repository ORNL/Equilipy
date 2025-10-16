import numpy as np
import equilipy.equilifort as fort

def _dict2np(dictionary):
    header = list(dictionary.keys())
    vals = list([])
    for i, h in enumerate(header):
        val = dictionary[h]
        try: len(val)
        except TypeError: val=[val]
            
        if i ==0:
            vals.append(val)
            l=len(val)
        elif len(val)==l:
            vals.append(val)
        else:
            print('Error: different length of colums ')
            raise ValueError()
    res = np.squeeze(np.asarray(vals).T)
    
    try: L,n=res.shape
    except ValueError:
        # When only one condition is given
        res = res.reshape((1,len(res)))
	
    return header, res


def _pyvar2fvar(var):       
 
	#Pass it to fortran	
	fort.moduleparsecs.info              = 0
	fort.modulethermoio.infothermo       = 0
	fort.moduleparsecs.nelementscs       = var.nElementsCS
	fort.moduleparsecs.nspeciescs        = var.nSpeciesCS
	fort.moduleparsecs.nsolnphasessyscs  = var.nSolnPhasesSysCS
	fort.moduleparsecs.ncountsublatticecs= var.nCountSublatticeCS
	fort.moduleparsecs.nmaxspeciesphasecs= var.nMaxSpeciesPhaseCS
	# fort.moduleparsecs.nparamcs          = var.nParamCS		#Could not find the use of this in fortran scripts
	# fort.moduleparsecs.nmagparamcs       = var.nMagParamCS	#Could not find the use of this in fortran scripts
	
	# Dimension(:)
	fort.moduleparsecs.iphasecs		      = var.iPhaseCS  # This might cause an issue
	fort.moduleparsecs.iparticlespermolecs = var.iParticlesPerMoleCS #
	fort.moduleparsecs.ngibbseqspecies	  = var.nGibbsEqSpecies
	fort.moduleparsecs.nsublatticephasecs  = var.nSublatticePhaseCS #
	fort.moduleparsecs.iphasesublatticecs  = var.iPhaseSublatticeCS #
	fort.moduleparsecs.datomicmasscs		  = var.dAtomicMass
	
	if var.nParamCS==0:
		#In Fortran parsing, this is an empty array
		fort.moduleparsecs.iparampasscs=np.zeros(1,dtype=int)
	else:
		fort.moduleparsecs.iparampasscs=np.zeros(var.nParamCS,dtype=int)
	if var.nMagParamCS==0:
		#In Fortran parsing, this is an empty array
		fort.moduleparsecs.imagparampasscs=np.zeros(1,dtype=int)
	else:
		fort.moduleparsecs.imagparampasscs=np.zeros(var.nMagParamCS,dtype=int)
	
	fort.moduleparsecs.cregularparamcs=np.asarray(var.cRegularParamCS,np.dtype('c'))
	fort.str1d('cElementNameCS',np.asarray(var.cElementNameCS,np.dtype('c')))
	fort.str1d('cSolnPhaseTypeCS',np.asarray(var.cSolnPhaseTypeCS,np.dtype('c')))
	fort.str1d('cSolnPhaseNameCS',np.asarray(var.cSolnPhaseNameCS,np.dtype('c')))
	fort.str1d('cSpeciesNameCS',np.asarray(var.cSpeciesNameCS,np.dtype('c')))
	fort.str2d('cPairNameCS',np.asarray(var.cPairNameCS,np.dtype('c')))
	fort.str3d('cConstituentNameSUBCS',np.asarray(var.cConstituentNameSUBCS,np.dtype('c')))
	

	#Appending 0th index to 1st index
	fort.moduleparsecs.nparamphasecs	= np.append([0],var.nParamPhaseCS) #CheckSystemExcess,Compthermodata
	fort.moduleparsecs.nspeciesphasecs = np.append([0],var.nSpeciesPhaseCS) #Checksystem,CheckSystemExcess,Compthermodata
	fort.moduleparsecs.nmagparamphasecs= np.append([0],var.nMagParamPhaseCS) #CheckSystemExcess,Compthermodata

	# Dimension(:,:)
	fort.moduleparsecs.npairssrocs				= var.nPairsSROCS
	fort.moduleparsecs.nconstituentsublatticecs  = var.nConstituentSublatticeCS
	fort.moduleparsecs.nsublatticeelementscs	    = var.nSublatticeElementsCS
	fort.moduleparsecs.dgibbsmagneticcs		    = var.dGibbsMagneticCS
	fort.moduleparsecs.dstoichsublatticecs		= var.dStoichSublatticeCS
	fort.moduleparsecs.dzetaspeciescs			= var.dZetaSpeciesCS
	fort.moduleparsecs.dgibbscoeffspeciestemp	= var.dGibbsCoeffSpeciesTemp
	fort.moduleparsecs.dstoichspeciescs	   	    = var.dStoichSpeciesCS
	fort.moduleparsecs.imagneticparamcs		    = var.iMagneticParamCS #Modify
	fort.moduleparsecs.dmagneticparamcs			= var.dMagneticParamCS #Modify
	fort.moduleparsecs.iregularparamcs			= var.iRegularParamCS #Modify
	fort.moduleparsecs.dregularparamcs			= var.dRegularParamCS #Modify


	 # Dimension(:,:,:)
	fort.moduleparsecs.iconstituentsublatticecs	 = var.iConstituentSublatticeCS
	fort.moduleparsecs.ipairidcs					 = var.iPairIDCS
	fort.moduleparsecs.ichemicalgroupcs			 = var.iChemicalGroupCS
	fort.moduleparsecs.dsublatticechargecs	     = var.dSublatticeChargeCS
	fort.moduleparsecs.dstoichpairscs			 = var.dStoichPairsCS
	fort.moduleparsecs.dconstituentcoefficientscs = var.dConstituentCoefficientsCS
	fort.moduleparsecs.dcoordinationnumbercs		 = var.dCoordinationNumberCS

	return None