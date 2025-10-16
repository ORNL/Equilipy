import numpy as np
import equilipy.variables as var

def ParseCSDataBlockQKTO(i):
	# RKMPM phases include magnetic mixing terms right before non-ideal mixing terms.
	# The end of the list of mixing terms is indicated by a "0".
	# Loop through excess parameters:
	if 'QKTOM' in var.cSolnPhaseTypeCS[i-1]:
		#Parsing magnetic terms
		while True:
			# Read in number of constituents involved in parameter:
			var.iMagneticParamCS[var.nMagParamCS,0]=int(var.DataBase.pop(0))

			if var.iMagneticParamCS[var.nMagParamCS,0] == 0:
				# The end of the parameter listing is marked by "0":
				break
			elif var.iMagneticParamCS[var.nMagParamCS,0] == 2:
				#Check if the parameter is binary
				var.iMagneticParamCS[var.nMagParamCS,1:5]=np.asarray(var.DataBase[:4],dtype=int)
				del var.DataBase[:4]
				var.dMagneticParamCS[var.nParamCS,:6]=np.asarray(var.DataBase[:6],dtype=float)
				del var.DataBase[:6]
			elif var.iMagneticParamCS[var.nMagParamCS,0] == 3:
				#Ternary mixing
				var.iMagneticParamCS[var.nMagParamCS,1:7]=np.asarray(var.DataBase[:6],dtype=int)
				del var.DataBase[:4]
				var.dMagneticParamCS[var.nParamCS,:6]=np.asarray(var.DataBase[:6],dtype=float)
				del var.DataBase[:6]
			else:
				# Not recognized
				var.INFO = 43
				return print('Error: the magnetic model is not recognized')
			var.nMagParamCS+=1

	while True:
		var.iRegularParamCS[var.nParamCS,0]=int(var.DataBase.pop(0))
		if var.iRegularParamCS[var.nParamCS,0]==0:
			break
		elif var.iRegularParamCS[var.nParamCS,0]==2:
			#Binary mixing terms
			var.iRegularParamCS[var.nParamCS,1:5]=np.asarray(var.DataBase[:4],dtype=int)
			del var.DataBase[:4]
			var.dRegularParamCS[var.nParamCS,:6]=np.asarray(var.DataBase[:6],dtype=float)
			del var.DataBase[:6]
		elif var.iRegularParamCS[var.nParamCS,0]==3:
			#Ternary mixing term
			var.iRegularParamCS[var.nParamCS,1:7]=np.asarray(var.DataBase[:6],dtype=int)
			del var.DataBase[:6]
			var.dRegularParamCS[var.nParamCS,:6]=np.asarray(var.DataBase[:6],dtype=float)
			del var.DataBase[:6]
		else:
			#Not recognized
			var.INFO = 1600 + i
			return print('Error: The Gibbs energy model is not recognized')
		var.nParamCS+=1
	return None
