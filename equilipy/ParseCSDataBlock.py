import numpy as np
import equilipy.variables as var
from .ParseCSDataBlockGibbs import *
from .ParseCSDataBlockQKTO import *
from .ParseCSDataBlockRKMP import *
from .ParseCSDataBlockSUBL import *
from .ParseCSDataBlockSUBG import *
from .utils import *


def ParseCSDataBlock():
	
	# Intialize variables:
	var.nCountSublatticeCS  = int(0)
	var.iCounterGibbsEqn    = int(0)
	var.nParamCS            = int(0)
	var.nMagParamCS         = int(0)
	var.dTempVec            = np.zeros(15,dtype=float)
	var.dGibbsMagneticCS    = np.zeros((var.nSpeciesCS,4),dtype=float)

	#Loop through all solution phases
	for i in range(1,var.nSolnPhasesSysCS+1):
		#Read phase name
		var.cSolnPhaseNameCS[i-1]="{:<25}".format(var.DataBase.pop(0))

		#Read model name
		var.cSolnPhaseTypeCS[i-1]="{:<8}".format(var.DataBase.pop(0))

		#Check if the model is supported
		if var.cSolnPhaseTypeCS[i-1] in var.cSolnPhaseTypeSupport:
			var.INFO = 0

		else:
			var.INFO = 17

		#Dummy parameter
		if i<=1: var.indx = int(0)+1
		else: var.indx=var.nSpeciesPhaseCS[i-2]+1
		
  
		#Count sublattice phases
		if var.cSolnPhaseTypeCS[i-1] in var.cSolnPhaseTypeSupport[5:]:
			var.nCountSublatticeCS = var.nCountSublatticeCS + 1
			var.iPhaseSublatticeCS[i-1] = var.nCountSublatticeCS

		# Read magnetic terms if necessary
		if var.cSolnPhaseTypeCS[i-1] in var.cSolnPhaseTypeSupport[3:6]:
			#Magnetic ordering is considered
			var.dTempVec[:2]=np.asarray(var.DataBase[:2],dtype=float)
			del var.DataBase[:2]
			for j in range(var.indx,var.nSpeciesPhaseCS[i-1]+1):
				var.dGibbsMagneticCS[j-1,2:4] = var.dTempVec[:2]
		elif var.cSolnPhaseTypeCS[i-1]==var.cSolnPhaseTypeSupport[7]:
			#Quadruplet model
			#Read zeta (FNN/SNN ratio)
			# In SUBG zeta is the same for every FNN pair
			# This method for keeping track of zeta is not very safe,
			# relies on having a predictable FNN pair order in the dat file
			dummy=float(var.DataBase.pop(0))
			var.dZetaSpeciesCS[var.nCountSublatticeCS-1,:var.nMaxSpeciesPhaseCS]=dummy*np.ones(var.nMaxSpeciesPhaseCS)
			# Read in two integers representing the number of species and the number of pairs:
			var.nPairsSROCS[var.nCountSublatticeCS-1,:2]=np.asarray(var.DataBase[:2],dtype=int)
			del var.DataBase[:2]
		elif var.cSolnPhaseTypeCS[i-1]==var.cSolnPhaseTypeSupport[8]:
			#MQM1 model
			#The SUBQ phase data files seems to not have the magnetic term so skipping this part.
			# Read in two integers representing the number of species and the number of pairs:
			var.nPairsSROCS[var.nCountSublatticeCS-1,:2]=np.asarray(var.DataBase[:2],dtype=int)
			del var.DataBase[:2]

		# Loop through species in solution phase:
		for j in range(var.indx,var.nSpeciesPhaseCS[i-1]+1):
			# SUBG and SUBQ phases contain a certain number of species, which are necessarily less
			# than the number of pair fractions. The # of species indicated in the
			# header file actually represents the number of pairs. Therefore, there are
			# fewer species listed than what has been allocated.
			if var.cSolnPhaseTypeCS[i-1] in var.cSolnPhaseTypeSupport[7:]:
				if j>=var.indx + var.nPairsSROCS[var.nCountSublatticeCS-1,0]:
					#Assign ID for all pairs
					var.iPhaseCS[j-1: var.nSpeciesPhaseCS[i-1]+1] = i
					break
			# Store the magnetic ordering terms for each solution:
			k = var.indx
			
			
			var.dGibbsMagneticCS[j-1,2:4] = var.dGibbsMagneticCS[k-1,2:4]

			# Store the phase index corresponding to the current species:
			var.iPhaseCS[j-1] = i

			if var.cSolnPhaseTypeCS[i-1] in var.cSolnPhaseTypeSupport[7:]:
				#The following subroutine parses the Gibbs energy equations (entries 3-5):
				ParseCSDataBlockGibbs(i,j)
				
				# Get pair stoichiometry in terms of constituents
				var.dConstituentCoefficientsCS[var.nCountSublatticeCS-1,j - var.indx,:5]=np.asarray(var.DataBase[:5],dtype=float)
				del var.DataBase[:5]
				if var.cSolnPhaseTypeCS[i-1] == var.cSolnPhaseTypeSupport[8]:
					var.dZetaSpeciesCS[var.nCountSublatticeCS-1,j - var.indx] = float(var.DataBase.pop(0))
			else:
				#The following subroutine parses the Gibbs energy equations (entries 3-5):				
				ParseCSDataBlockGibbs(i,j)
			if 'QKTO' in var.cSolnPhaseTypeCS[i-1]:
				var.dTempVec=np.asarray(var.DataBase[:2],dtype=float)
				del var.DataBase[:2]
				if 'QKTOM' in var.cSolnPhaseTypeCS[i-1]:
					if sum(var.dTempVec-[1,2])>1E-10:
						var.INFO = 1600 + i
						print('Error: QKTOM')


		# Check the type of solution phase to interpret mixing parameters:
		
		# Ideal mixing (IDMX)
		if var.cSolnPhaseTypeCS[i-1]==var.cSolnPhaseTypeSupport[0]:
			pass
		# Quasichemical Kohler-Toop model
		elif 'QKTO' in var.cSolnPhaseTypeCS[i-1]:
			ParseCSDataBlockQKTO(i)
		# Redlich-Kister-Muggiano-Polynomial model
		elif 'RKMP' in var.cSolnPhaseTypeCS[i-1]:
			ParseCSDataBlockRKMP(i)
		# Compound Energy Formalism (sublattice) model
		elif 'SUBL' in var.cSolnPhaseTypeCS[i-1]:
			ParseCSDataBlockSUBL(i)
		# Quadruplet quasichemical model:
		elif 'SUBG' in var.cSolnPhaseTypeCS[i-1]:
			ParseCSDataBlockSUBG(i)
		# Modified quasichemical model:
		elif 'SUBQ' in var.cSolnPhaseTypeCS[i-1]:
			ParseCSDataBlockSUBG(i)
		else:
			# The solution phase type is not supported. Report an error.
			var.INFO = 17
			return
		# Record the index of the mixing parameter for this phase:
		var.nParamPhaseCS[i-1] = var.nParamCS
		var.nMagParamPhaseCS[i-1] = var.nMagParamCS
		
	# END of Loop_SolnPhases
	
	# Begin parsing pure condensed phases:
	for j in range(var.nSpeciesPhaseCS[var.nSolnPhasesSysCS-1]+1,var.nSpeciesCS+1):
		# The phase index of a pure separate phase is set to zero:
		var.iPhaseCS[j-1] = 0
		ParseCSDataBlockGibbs(var.iPhaseCS[j-1],j)

	return None
		#Nest functions
