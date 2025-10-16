import numpy as np, regex as re
import equilipy.variables as var
from .utils import *


def ParseCSDataBlockGibbs(i,j):
    
	if var.FactSage8Plus:
		ChangeBlockGibbs()
	
	#Consider the case when the species name includes some space
	cdum=''
	while True:
		temp=0
		try:
			int(var.DataBase[temp])
			break
		except ValueError:
			try:
				#I should find out the meaning of these parameters
				float(var.DataBase[temp])
				var.DataBase.pop(0)
			except ValueError:
				cdum+=var.DataBase.pop(0)

	var.cSpeciesNameCS[j-1]="{:<30}".format(str(cdum))
	var.cEndmemberNameCS[j-1]="{:<30}".format(str(cdum))

	#Check to see if there is more than one particle per constituent formula mass.
	#Record the number of particles per constituent formula mass. (stoichio of constituents)
	if '):' in var.cSpeciesNameCS[j-1]:
		cDummy=var.cSpeciesNameCS[j-1]
		var.iParticlesPerMoleCS[j-1]=re.search(r'\d+', cDummy[cDummy.find(':'):]).group()


	# Check for a dummy species (a description of a "dummy species" is given in the header).
	if i==0 and '#' in var.cSpeciesNameCS[j-1]:
		#The phase index of a dummy species is set to -1:
		var.iPhaseCS[j-1] = -1
		i = -1

	# Entry 4: Read thermodynamic data for constituent species:
	iGibbsEqType=int(var.DataBase.pop(0))
	var.nGibbsEqSpecies[j-1]=int(var.DataBase.pop(0))
	var.dTempVec=np.asarray(var.DataBase[:var.nElementsCS],dtype=float)
	del var.DataBase[:var.nElementsCS]

	if iGibbsEqType not in [1,2,3,4,5,6,13,16]:
		# The type of Gibbs energy equation is not supported.
		var.INFO = 1400 + i
		return print('The type of Gibbs energy equation is not supported.')

	# Note: ChemSage data-files store the stoichiometry coefficients as real variables, which may
	# be less than one when there are more than 1 particles/mole.  I am going to convert the stoichiometry
	# to an integer to minimize numerical error while still recording the number of particles per mole.
	dTemp = float(var.iParticlesPerMoleCS[j-1])
	var.dStoichSpeciesCS[j-1,:var.nElementsCS] = var.dTempVec * dTemp

	#Loop_GibbsEquations
	for m in range(var.nGibbsEqSpecies[j-1]):
		# Entry 5: Read upper temperature limit and the coefficients for the Gibbs energy equation for this species:
		if var.iCounterGibbsEqn > len(var.dGibbsCoeffSpeciesTemp[0,:]):
			var.INFO = 5
			return print('Error: Entry 5')
		# Read in the standard Gibbs energy terms:

		var.dGibbsCoeffSpeciesTemp[:7,var.iCounterGibbsEqn]=np.asarray(var.DataBase[:7],dtype=float)
		del var.DataBase[:7]

		l=0


		if iGibbsEqType in [4,16,5,6]:
			# Note that Gibbs energy equations designated type 4 and 16 have an additional
			# line, but type 1 does not.
			# Read in the number of additional terms to the standard Gibbs energy equation:
			l=int(var.DataBase.pop(0))
		if l==0:
			pass
		elif l==1:
			var.dGibbsCoeffSpeciesTemp[7:9, var.iCounterGibbsEqn]=np.asarray(var.DataBase[:2],dtype=float)
			del var.DataBase[:2]
		elif l==2:
			var.dGibbsCoeffSpeciesTemp[7:11, var.iCounterGibbsEqn]=np.asarray(var.DataBase[:4],dtype=float)
			del var.DataBase[:4]
		elif l==3:
			var.dGibbsCoeffSpeciesTemp[7:13, var.iCounterGibbsEqn]=np.asarray(var.DataBase[:6],dtype=float)
			del var.DataBase[:6]
		else:
			var.INFO = 1500 + i
			return print('Error: Assigning Gibbs energy equations')
			#SY: revised
		var.iCounterGibbsEqn = var.iCounterGibbsEqn + 1
	#END Loop_GibbsEquations
 
	# Check if the equation has magnetic contributions:
	if iGibbsEqType in [13,16]:
		# This species has magnetic contributions to the Gibbs energy equation:

		# Check if this is a pure condensed phase or solution phase:
		if i==0:
			# Read all four coefficients for a pure condensed phase:
			
			var.dGibbsMagneticCS[j-1,:4]=np.asarray(var.DataBase[:4],dtype=float)
			del var.DataBase[:4]
		else:
		# Read the first two coefficients for a solution phase:
			var.dGibbsMagneticCS[j-1,:2]=np.asarray(var.DataBase[:2],dtype=float)
			del var.DataBase[:2]
	if iGibbsEqType in [2,5]:
			#Molare volume calculation is not implemented yet
			var.dMolarVolume=float(var.DataBase.pop(0))
	if iGibbsEqType in [3,6]:
			#Molare volume calculation is not implemented yet
			var.dMolarVolume=np.asarray(var.DataBase[:11],dtype=float)
			del var.DataBase[:11]
	#Assign an arbitrarily large positive value for the Gibbs energy equations for dummy species:
	if var.iPhaseCS[j-1] == -1:
		var.dGibbsCoeffSpeciesTemp[1,var.iCounterGibbsEqn-1] = 1E6

	return None



def ChangeBlockGibbs():
	
	res=[]
	# add the name of endmember
	EndmemberName=var.DataBase.pop(0)
	res.append(EndmemberName)
	phasename=[x for x in var.cSolnPhaseNameCS if x !="{:<25}".format(' ')]
	
	# Delete dummy integers
	if float(var.DataBase[0])==0.0 and float(var.DataBase[1])==0.0:
		del var.DataBase[:2]
	
	# Revise iGibbsEqType
	res.append(str(int(var.DataBase.pop(0))-6))
	
	# Read the number of Cp intervals for the endmember
	nGibbsInterval=int(var.DataBase.pop(0))
	res.append(str(nGibbsInterval))
	
	res=res+var.DataBase[:var.nElementsCS]
	del var.DataBase[:var.nElementsCS]
	
	
	# Obtain the number of functions to combine
	nFuncions=int(var.DataBase.pop(0))
	
	
	# Not sure what this parameter is. It is always zero
	dummy = int(var.DataBase.pop(0))
	if dummy !=0: print('non-zero dummy')

	#Obtain function names and their coefficients
	if nFuncions==0:
		#Stoichiometric compounds
		H298 = float(var.DataBase.pop(0))
		S298 = float(var.DataBase.pop(0))
		Cp = np.zeros((nGibbsInterval,14),dtype=float)
		for i in range(nGibbsInterval):
			# allocate T min
			if i==0: 
				Cp[i,0]=298.15
			else:
				Cp[i,0]=Cp[i-1,1]
			# allocate Tmax
			Cp[i,1]=float(var.DataBase.pop(0))
			# In ChemSage, Cp = a + b*T + c*T^2 +d*T(-2)
			Cp[i,2:4]=var.DataBase[:2]
			del var.DataBase[:2]
			Cp[i,5]=float(var.DataBase.pop(0))
			Cp[i,4]=float(var.DataBase.pop(0))
			
			
			try: 
				nCustomCp=int(var.DataBase.pop(0))
			except ValueError:
				nCustomCp= int(var.DataBase.pop(0))
				
			
			for j in range(nCustomCp):
				# coefficient
				Cp[i,6+2*j]=float(var.DataBase.pop(0))
				# power
				Cp[i,7+2*j]=float(var.DataBase.pop(0))
		if Cp[-1,1]<6000:Cp[-1,1]=6000
		G=HSCp2G(H298,S298,Cp)
		var.CompoundGibbs[f'{phasename[-1].strip()}_{EndmemberName}']=G
	elif nFuncions==1:
		del var.DataBase[:3]
		Multiplier=var.DataBase.pop(0)
		FunctionName=var.DataBase.pop(0)
		# Convert H S Cp into Gibbs
		H298 = var.FunctionsHSCp[FunctionName]['H298']
		S298 = var.FunctionsHSCp[FunctionName]['S298']
		Cp	 = var.FunctionsHSCp[FunctionName]['Cp']
		G=HSCp2G(H298,S298,Cp)
		H, S,cp_t = G2HSCp(G)
		var.CompoundGibbs[f'{phasename[-1].strip()}_{EndmemberName}']=G
		
	else:
		# Calculate Neuman-Kopp of functions given in H, S, Cp
		FunctionNames=[]
		Multipliers=[]
		for i in range(nFuncions):
			# Ignore dummy parameters
			del var.DataBase[:3]

			Multipliers.append(float(var.DataBase.pop(0)))
			FunctionNames.append(str(var.DataBase.pop(0)))
		# print(FunctionNames)
		# print(Multipliers)
		NKData=NeumanKoppHSCp(var.FunctionsHSCp,FunctionNames,Multipliers)
		G=HSCp2G(NKData['H298'],NKData['S298'],NKData['Cp'])
		var.CompoundGibbs[f'{phasename[-1].strip()}_{EndmemberName}']=G
	
	nCpIntervals,_=G.shape
	
	

	if nGibbsInterval!=nCpIntervals: 
		raise ValueError
	G=list(G[:,1:].flatten())
	for i in range(nCpIntervals):
		Ginterval= G[:15]
		del G[:15]
		# Save default parameters
		res=res + Ginterval[:7]
		del Ginterval[:7]

		# count the number of non-zero coefficients
		CustomPara=np.array(Ginterval)
		nCustomPara=np.count_nonzero(CustomPara[[0,2,4,6]])
		
		if nCustomPara==0:
			res = res + [1, 0.0, 0.0]
		else:
			for j in range(4):
				CustomCoeff=Ginterval[0]
				if CustomCoeff==0.0: 
					continue
				
				CustomPower=Ginterval[1]
				if j==0: 
					res.append(nCustomPara)
				res.append(CustomCoeff)
				res.append(CustomPower)
				del Ginterval[:2]
	
	var.DataBase[:0]=res
	return None