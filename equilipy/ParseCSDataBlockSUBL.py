import numpy as np
import equilipy.variables as var

def ParseCSDataBlockSUBL(i):
	n = 0

	#Read in the number of sublattices per phase:
	var.nSublatticePhaseCS[var.nCountSublatticeCS-1]=int(var.DataBase.pop(0))

	# Report an error if the number of sublattices is out of range:
	if var.nSublatticePhaseCS[var.nCountSublatticeCS-1]<=0 or var.nSublatticePhaseCS[var.nCountSublatticeCS-1]>var.nMaxSublatticeCS:
		var.INFO=33
		return print('Error: The number of sublattices is out of range')

	#Read in the stoichiometry coefficients for each sublattice:
	var.dStoichSublatticeCS[var.nCountSublatticeCS-1,:var.nSublatticePhaseCS[var.nCountSublatticeCS-1]]=\
		np.asarray(var.DataBase[:var.nSublatticePhaseCS[var.nCountSublatticeCS-1]],dtype=float)
	del var.DataBase[:var.nSublatticePhaseCS[var.nCountSublatticeCS-1]]

	#Read in the number of constituents for each sublattice:
	var.nConstituentSublatticeCS[var.nCountSublatticeCS-1,:var.nSublatticePhaseCS[var.nCountSublatticeCS-1]]=\
		np.asarray(var.DataBase[:var.nSublatticePhaseCS[var.nCountSublatticeCS-1]],dtype=int)
	del var.DataBase[:var.nSublatticePhaseCS[var.nCountSublatticeCS-1]]

	#Read in the name of each constituent for each sublattice:
	#LOOP_SUBL_CONST_NAME
	for s in range(var.nSublatticePhaseCS[var.nCountSublatticeCS-1]):
		#Number of constituents on last line (per sublattice):
		for m in range(var.nConstituentSublatticeCS[var.nCountSublatticeCS-1,s]):
			var.cConstituentNameSUBCS[var.nCountSublatticeCS-1, s, m,:]=np.asarray("{:<8}".format(str(var.DataBase.pop(0))),np.dtype('c'))


	# Read in the constituent indices for each component on each sublattice
	#LOOP_SUBL_CONST_ID
	for s in range(var.nSublatticePhaseCS[var.nCountSublatticeCS-1]):
		n=var.nSolnPhaseCS[i-1]
		for m in range(n):
			var.iConstituentSublatticeCS[var.nCountSublatticeCS-1, s, m]=int(var.DataBase.pop(0))
	# SUBLM phases include magnetic mixing terms right before non-ideal mixing terms.
	# The end of the list of mixing terms is indicated by a "0".
	if 'SUBLM' in var.cSolnPhaseTypeCS[i-1]:
		#LOOP_MagneticMixingSUBL
		while True:
			# Read in number of constituents involved in parameter:
			var.iMagneticParamCS[var.nMagParamCS,0]=int(var.DataBase.pop(0))
			if var.iMagneticParamCS[var.nMagParamCS,0]==0:
				# The end of the section of mixing terms is labelled "0".
				break
			else:
				# Read in the list of constituent indices involved in this parameter:
				var.iMagneticParamCS[var.nMagParamCS,1:var.iMagneticParamCS[var.nMagParamCS,0]+2]=\
					np.asarray(var.DataBase[:var.iMagneticParamCS[var.nMagParamCS,0]+1],dtype=int)
				del var.DataBase[:var.iMagneticParamCS[var.nMagParamCS,0]+1]

				j=var.iMagneticParamCS[var.nMagParamCS,0]
				# Read in the mixing parameter:
				var.dMagneticParamCS[var.nMagParamCS,:2]=np.asarray(var.DataBase[:2],dtype=float)
				del var.DataBase[:2]

				# Store the number of equations to temporary memory
				n = var.iMagneticParamCS[var.nMagParamCS, j+1]

				# Correct the indexing scheme (order of mixing) for the first parameter:
				var.iMagneticParamCS[var.nMagParamCS, j+1] = 0

				# Loop through number of mixing terms per component array:
				for k in range(2,n+1):
					# Update counter of the number of parameters:
					var.nMagParamCS += 1

					# Record mixing ID's:
					var.iMagneticParamCS[var.nMagParamCS,:var.nParamMax*2+1] = var.iMagneticParamCS[var.nMagParamCS-1,:var.nParamMax*2+1]

					# Correct the indexing scheme (order of mixing parameter):
					var.iMagneticParamCS[var.nMagParamCS, j+1] = k - 1

					# Read mixing terms:
					var.dMagneticParamCS[var.nMagParamCS,:2] = np.asarray(var.DataBase[:2],dtype=float)
					del var.DataBase[:2]
				# Update counter of the number of parameters:
				var.nMagParamCS += 1

	#Loop through excess mixing parameters:
	#LOOP_ExcessMixingSUBL
	while True:
		# Read in number of constituents involved in parameter:
		var.iRegularParamCS[var.nParamCS,0]=int(var.DataBase.pop(0))
		
		if var.iRegularParamCS[var.nParamCS,0]==0:
			# The end of the section of mixing terms is labelled "0".
			break
		else:
			# Read in the list of constituent indices involved in this parameter:
			j=var.iRegularParamCS[var.nParamCS,0]

			# Read in the list of constituent indices involved in this parameter:
			var.iRegularParamCS[var.nParamCS,1:j+2]=np.asarray(var.DataBase[:j+1],dtype=int)
			del var.DataBase[:j+1]
			
			# Read in the mixing parameter:
			var.dRegularParamCS[var.nParamCS,:6]=np.asarray(var.DataBase[:6],dtype=float)
			del var.DataBase[:6]

			# Store the number of equations to temporary memory
			n = var.iRegularParamCS[var.nParamCS, j+1]

			# Correct the indexing scheme (order of mixing) for the first parameter:
			var.iRegularParamCS[var.nParamCS, j+1] = 0
			# Loop through number of mixing terms per component array:
			#LOOP_ParamArray:
			for k in range(2,n+1):
				# Update counter of the number of parameters:
				var.nParamCS += 1

				# Record mixing ID's:
				var.iRegularParamCS[var.nParamCS,:var.nParamMax*2+1] = var.iRegularParamCS[var.nParamCS-1,:var.nParamMax*2+1]

				# Correct the indexing scheme (order of mixing parameter):
				var.iRegularParamCS[var.nParamCS, j+1] = k - 1

				# Read mixing terms:
				# var.dRegularParamCS[var.nMagParamCS,:6] = np.asarray(var.DataBase[:6],dtype=float)
				var.dRegularParamCS[var.nParamCS,:6] = np.asarray(var.DataBase[:6],dtype=float)
				del var.DataBase[:6]
		# Update counter of the number of parameters:
		var.nParamCS += 1
	return None
