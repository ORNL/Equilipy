import numpy as np
import equilipy.variables as var

def ParseHSCpFunctions():
	funcs = dict()
	nFunctions=int(var.DataBase.pop(0))
	del var.DataBase[:3]
	for k in range(nFunctions):
		# Process function one by one (Read_func_HSCp)
		FuncsName=str(var.DataBase.pop(0))
		nCpIntervals=int(var.DataBase.pop(0))
		Cp = np.zeros((nCpIntervals,14),dtype=float)
		funcs[FuncsName] = {}
		funcs[FuncsName]['H298']=float(var.DataBase.pop(0))
		funcs[FuncsName]['S298']=float(var.DataBase.pop(0))
		for i in range(nCpIntervals):
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
			
			nCustomCp= int(var.DataBase.pop(0))
			
			for j in range(nCustomCp):
				# coefficient
				Cp[i,6+2*j]=float(var.DataBase.pop(0))
				# power
				Cp[i,7+2*j]=float(var.DataBase.pop(0))
		if Cp[-1,1]<6001:Cp[-1,1]=6001
		funcs[FuncsName]['Cp'] = Cp	
	var.FunctionsHSCp= funcs
	var.CompoundGibbs={
		'Descriptor':['Tmin','Tmax','T**0', 'T**1', 'T*ln(T)', 'T**2' ,'T**3', 'T**(-1)', 'T**(p1)','p1', 'T**(p2)','p2', 'T**(p3)', 'p3','T**(p4)', 'p4']
	}
	return None