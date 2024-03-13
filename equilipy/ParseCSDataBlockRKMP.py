import numpy as np
import equilipy.variables as var

def ParseCSDataBlockRKMP(i):
	# RKMPM phases include magnetic mixing terms right before non-ideal mixing terms.
	# The end of the list of mixing terms is indicated by a "0".
	
	if 'RKMPM' in var.cSolnPhaseTypeCS[i-1]:
		#Parsing magnetic terms
		while True:
			# Read in number of constituents involved in parameter:
			var.iMagneticParamCS[var.nMagParamCS,0]=int(var.DataBase.pop(0))

			if var.iMagneticParamCS[var.nMagParamCS,0] == 0:
				# The end of the parameter listing is marked by "0":
				break
			elif var.iMagneticParamCS[var.nMagParamCS,0] == 2:
				#Check if the parameter is binary
				var.iMagneticParamCS[var.nMagParamCS,1:4]=np.asarray(var.DataBase[:3],dtype=int)
				del var.DataBase[:3]
				#Store the number of equations to temporary memory
				k = var.iMagneticParamCS[var.nMagParamCS,3]
				# Loop through the number of terms per constituent matrix
				for j in range(1,k+1):
					# Record mixing ID's:
					var.iMagneticParamCS[var.nMagParamCS,:3] = var.iMagneticParamCS[var.nMagParamCS-j+1,:3]

					# Set the exponent:
					var.iMagneticParamCS[var.nMagParamCS,3]   = j - 1
					var.dMagneticParamCS[var.nMagParamCS,:2] = np.asarray(var.DataBase[:2],dtype=float)
					del var.DataBase[:2]

					# Count the number of parameters:
					var.nMagParamCS+=1
			elif var.iMagneticParamCS[var.nMagParamCS,0] == 3:
				#Ternary mixing
				var.iMagneticParamCS[var.nMagParamCS,1:5]=np.asarray(var.DataBase[:4],dtype=int)
				del var.DataBase[:4]

				#Store the number of equations to temporary memory
				k = var.iMagneticParamCS[var.nMagParamCS,4]

				# Loop through the number of terms per constituent matrix
				for j in range(1,k+1):
					# Record mixing ID's:
					var.iMagneticParamCS[var.nMagParamCS,:5] = var.iMagneticParamCS[var.nMagParamCS-j+1,:5]

					# Set the exponent:
					var.iMagneticParamCS[var.nMagParamCS,4]   = var.iMagneticParamCS[var.nMagParamCS,j]
					var.dMagneticParamCS[var.nMagParamCS,:2] = np.asarray(var.DataBase[:2],dtype=float)
					del var.DataBase[:2]

					# Count the number of parameters:
					var.nMagParamCS+=1
			elif var.iMagneticParamCS[var.nMagParamCS,0] == 4:
				#Quaternary mixing
				var.iMagneticParamCS[var.nMagParamCS,1:6]=np.asarray(var.DataBase[:5],dtype=int)
				del var.DataBase[:5]

				#Store the number of equations to temporary memory
				k = var.iMagneticParamCS[var.nMagParamCS,5]

				# Loop through the number of terms per constituent matrix
				for j in range(1,k+1):
					# Record mixing ID's:
					var.iMagneticParamCS[var.nMagParamCS,:6] = var.iMagneticParamCS[var.nMagParamCS-j+1,:6]

					# Set the exponent:
					var.iMagneticParamCS[var.nMagParamCS,5]   = var.iMagneticParamCS[var.nMagParamCS,j]
					var.dMagneticParamCS[var.nMagParamCS,:2] = np.asarray(var.DataBase[:2],dtype=float)
					del var.DataBase[:2]

					# Count the number of parameters:
					var.nMagParamCS+=1
			else:
				# Not recognized
				var.INFO = 43
				return print('Error: The magnetic model is not recognized')

	while True:
		#Loop through excess mixing parameters:

		# Read in number of constituents involved in parameter
		var.iRegularParamCS[var.nParamCS,0]= int(var.DataBase.pop(0))
		
		if var.iRegularParamCS[var.nParamCS,0]==0:
			break
		elif var.iRegularParamCS[var.nParamCS,0]==2:
			#Binary mixing
			var.iRegularParamCS[var.nParamCS,1:4]=np.asarray(var.DataBase[:3],dtype=int)
			del var.DataBase[:3]
			#Store the number of equations to temporary memory
			k = var.iRegularParamCS[var.nParamCS,3]
			#Loopt thorugh the number of terms per constituent matrix
			for j in range(1,k+1):
				#Record mixing ID's
				var.iRegularParamCS[var.nParamCS,:4]=var.iRegularParamCS[var.nParamCS-j+1,:4]
				#Set the exponent for RKMP
				var.iRegularParamCS[var.nParamCS,3]=j - 1
				var.dRegularParamCS[var.nParamCS,:6]=np.asarray(var.DataBase[:6],dtype=float)

				
				del var.DataBase[:6]
				#Count the number of parameters:
				var.nParamCS+=1


		elif var.iRegularParamCS[var.nParamCS,0]==3:
			#Ternary mixing
			var.iRegularParamCS[var.nParamCS,1:5]=np.asarray(var.DataBase[:4],dtype=int)
			del var.DataBase[:4]

			#Store the number of equations to temporary memory
			k = var.iRegularParamCS[var.nParamCS,4]

			#Loopt thorugh the number of terms per constituent matrix
			for j in range(1,k+1):
				#Record mixing ID's
				var.iRegularParamCS[var.nParamCS,:5]=var.iRegularParamCS[var.nParamCS-j+1,:5]

				#Set the exponent for RKMP
				var.iRegularParamCS[var.nParamCS,4]=var.iRegularParamCS[var.nParamCS,j]

				var.dRegularParamCS[var.nParamCS,:6]=np.asarray(var.DataBase[:6],dtype=float)
				del var.DataBase[:6]

				#Count the number of parameters:
				var.nParamCS+=1

		elif var.iRegularParamCS[var.nParamCS,0]==4:
			#Quaternary mixing
			var.iRegularParamCS[var.nParamCS,1:6]=np.asarray(var.DataBase[:5],dtype=int)
			del var.DataBase[:5]

			#Store the number of equations to temporary memory
			k = var.iRegularParamCS[var.nParamCS,5]

			#Loopt thorugh the number of terms per constituent matrix
			for j in range(1,k+1):
				#Record mixing ID's
				var.iRegularParamCS[var.nParamCS,:6]=var.iRegularParamCS[var.nParamCS-j+1,:6]

				#Set the exponent for RKMP
				var.iRegularParamCS[var.nParamCS,5]=var.iRegularParamCS[var.nParamCS,j]

				var.dRegularParamCS[var.nParamCS,:6]=np.asarray(var.DataBase[:6],dtype=float)
				del var.DataBase[:6]

				#Count the number of parameters:
				var.nParamCS+=1

		else:
			#Not recognized
			var.INFO = 1600 + i
			return print('Error: the Gibbs energy model is not recognized')
		
	return None
