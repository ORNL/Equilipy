import numpy as np
import equilipy.variables as var

def ParseCSHeader():
	# Part 2: nElement,nSln+1(nGas),nGas species, nSln species,...,nPureCompunds
	var.nElementsCS=int(var.DataBase.pop(0))
	var.nSolnPhasesSysCS=int(var.DataBase.pop(0))
	var.nSpeciesPhaseCS=np.asarray(var.DataBase[:var.nSolnPhasesSysCS],dtype=int)

	del var.DataBase[:var.nSolnPhasesSysCS]

	if var.nSolnPhasesSysCS!=1 and var.nSpeciesPhaseCS[0]==0:
		#No gas species
		var.nSpeciesPhaseCS=np.delete(var.nSpeciesPhaseCS,0)
		var.nSolnPhasesSysCS = var.nSolnPhasesSysCS - 1
	var.nSolnPhaseCS=var.nSpeciesPhaseCS.copy()
	var.nPureSpeciesCS=int(var.DataBase.pop(0))

	if var.nSolnPhasesSysCS==1 and var.nSpeciesPhaseCS[0]==0: var.nSolnPhasesSysCS=0
	
	# Check the number of solution phases are under max number
	assert var.nSolnPhasesSysCS<=var.nSolnPhasesSysMax, \
		'Error: Number of solution phases exceeds the maximum number: '+str(var.nSolnPhasesSysMax)

		# Store the maximum number of species per solution phase in the system:
	var.nMaxSpeciesPhaseCS = max(var.nSpeciesPhaseCS)

	#Compute the number of species per solution phase
	#Why are we doing this? This method should be corrected to simplify code
	for i in range(1,var.nSolnPhasesSysCS):
		var.nSpeciesPhaseCS[i]=var.nSpeciesPhaseCS[i]+var.nSpeciesPhaseCS[i-1]
	#Compute the total number of species in the system
	var.nSpeciesCS = var.nPureSpeciesCS + var.nSpeciesPhaseCS[var.nSolnPhasesSysCS-1]

	#Memory allocation
	j = max(1,var.nSolnPhasesSysCS)
	var.nParamPhaseCS               = np.zeros(j,dtype=int)
	var.nMagParamPhaseCS            = np.zeros(j,dtype=int)
	var.iPhaseSublatticeCS          = np.zeros(j,dtype=int)
	var.nPairsSROCS                 = np.zeros((j,2),dtype=int)
	var.iPairIDCS                   = np.zeros((var.nSolnPhasesSysCS,var.nMaxSpeciesPhaseCS,4),dtype=int)
	var.nGibbsEqSpecies             = np.zeros(var.nSpeciesCS,dtype=int)
	var.iPhaseCS                    = np.zeros(var.nSpeciesCS,dtype=int)
	var.iParticlesPerMoleCS         = np.ones(var.nSpeciesCS, dtype=int)
	var.iMagneticParamCS            = np.zeros((1000,var.nParamMax*2+3),dtype=int)
	var.iRegularParamCS             = np.zeros((1000,var.nParamMax*2+3),dtype=int)
	var.nSublatticePhaseCS          = np.zeros(var.nSolnPhasesSysCS,dtype=int)
	var.nConstituentSublatticeCS    = np.zeros((var.nSolnPhasesSysCS,var.nMaxSublatticeCS),dtype=int)
	var.iConstituentSublatticeCS    = np.zeros((var.nSolnPhasesSysCS,var.nMaxSublatticeCS,var.nMaxSpeciesPhaseCS),dtype=int)
	var.nSublatticeElementsCS       = np.zeros((var.nSolnPhasesSysCS,var.nMaxSublatticeCS),dtype=int)
	var.iChemicalGroupCS            = np.zeros((var.nSolnPhasesSysCS,var.nMaxSublatticeCS,var.nMaxSpeciesPhaseCS),dtype=int)

	#Memory allocation
	var.dAtomicMass                 = np.zeros(var.nElementsCS,dtype=float)
	var.dCoordinationNumberCS       = np.zeros((var.nSolnPhasesSysCS,var.nMaxSpeciesPhaseCS,4),dtype=float)
	var.dStoichSpeciesCS            = np.zeros((var.nSpeciesCS,var.nElementsCS),dtype=float)
	var.dGibbsCoeffSpeciesTemp      = np.zeros((var.nGibbsCoeff,var.nSpeciesCS*var.nMaxGibbsEqs),dtype=float)
	var.dRegularParamCS             = np.zeros((1000,6),dtype=float)
	var.dGibbsMagneticCS            = np.zeros((var.nSpeciesCS,4),dtype=float)
	var.dMagneticParamCS            = np.zeros((1000,2),dtype=float)
	var.dStoichSublatticeCS         = np.zeros((var.nSolnPhasesSysCS,var.nMaxSublatticeCS),dtype=float)
	var.dZetaSpeciesCS              = np.zeros((var.nSolnPhasesSysCS,var.nMaxSpeciesPhaseCS),dtype=float)
	var.dConstituentCoefficientsCS  = np.zeros((var.nSolnPhasesSysCS,var.nMaxSpeciesPhaseCS,5),dtype=float)
	var.dSublatticeChargeCS         = np.zeros((var.nSolnPhasesSysCS,var.nMaxSublatticeCS,var.nMaxSpeciesPhaseCS),dtype=float)
	var.dStoichPairsCS              = np.zeros((var.nSolnPhasesSysCS,var.nMaxSpeciesPhaseCS,var.nElementsCS),dtype=float)

	var.cSolnPhaseNameCS      =["{:<25}".format('')]*j
	var.cSolnPhaseTypeCS      =["{:<8}".format('')]*j
	var.cSpeciesNameCS        =["{:<30}".format('')]*var.nSpeciesCS
	var.cEndmemberNameCS        =["{:<30}".format('')]*var.nSpeciesCS
	var.cRegularParamCS       =["{:<1}".format('')]*1000
	var.cPairNameCS           =np.empty((var.nSolnPhasesSysCS,var.nMaxSpeciesPhaseCS,30),dtype='c')
	var.cConstituentNameSUBCS =np.empty((var.nSolnPhasesSysCS,var.nMaxSublatticeCS,var.nMaxSpeciesPhaseCS,8),dtype='c')

	# Part 3: list of system components
	var.cElementNameCS=['{:<3}'.format(x) for x in var.DataBase[:var.nElementsCS]]
	del var.DataBase[:var.nElementsCS]
	for i in range(len(var.cElementNameCS)):
		if var.cElementNameCS[i][0]=='e': var.cElementNameCS[i]='e- '

	# Part 4: list of atomic mass
	var.dAtomicMass=np.asarray(var.DataBase[:var.nElementsCS],dtype=float)
	del var.DataBase[:var.nElementsCS]

	# Part 5: Definition of the temperature dependence terms
	iGequation=np.asarray(var.DataBase[:7],dtype=int)
	assert sum(iGequation-[6,1,2,3,4,5,6])==0, 'Themperature dependent G is not implemented'
	del var.DataBase[:14]
	return var.DataBase
