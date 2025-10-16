import numpy as np
from typing import List, Tuple, Union, Dict, Any
import numpy as np
import equilipy.variables as var
import equilipy.equilifort as fort

def G2HSCp_single(GibbsEnegy: List[float]) -> Tuple[float, float, np.ndarray]:
	'''
	Convert Gibbs energy [J] of one temperature interval to H298 [J], S298 [J/(mol K)], Cp(T) [J/(mol K)]

	Input
	-----
	GibbsEnegy : A list of 14 floats [a, b, c, d, e, f, c1, p1, c2, p2, c3, p3, c4, p4]
				a-f are coefficients of predefined function, ci and pi are coefficient and power of custom functions
				G = a +b*T +c*T*s.ln(T) +d*T**2 +e*T**3 +f*T**(-1) +c1*T**(p1) +c2*T**(p2) +c3*T**(p3) +c4*T**(p4) 

	Output
	------
	H298 : Enthalpy of formation at 298.15 K
	S298 : Vibration entropy at 298.15 K
	Cp(T): A list of 12 floats [a, b, c, d, c1, p1, c2, p2, c3, p3, c4, p4]
		  Cp(T) = a +b*T +c*T**(-2) +d*T**(2) +c1*T**(p1) +c2*T**(p2) +c3*T**(p3) +c4*T**(p4) 
	
	'''
	
	# assert isinstance(Gibbs,list), "Error: input must be a list with length 14"
	assert len(GibbsEnegy)==14, "Error: input must have length of 14"
 
	a, b, c, d, e, f, c1, p1, c2, p2, c3, p3, c4, p4 = GibbsEnegy
	custom= [[c1,p1],[c2,p2],[c3,p3],[c4,p4]]
	
	# Allocate Cp
	Cp = np.zeros(12,dtype=float)
	Cp[0]= -c       # a term
	Cp[1]= -2*d     # +b*T term
	Cp[2]= -2*f     #+c*T**(-2) term
	Cp[3]= -6*e     #+d*T**(2) term
	
	# Allocate H298 and S298   
	T=298.15
	S298= - (b+c) - 2*d*T + f/T**2 -3*e*T**2 - c*np.log(T)  
	H298= + a - c *T - d*T**2 -2*e*T**3+ 2*f/T
	
	# Manage customized functions
	for i, param in enumerate(custom):
		coeff, power= param
		if abs(power)<1E-5: continue
		
		if abs(power+1)<1E-5: 
			print('G2HSCp warning: customized coefficient and power are not considered as it appears as -1')
			continue
		
		S298= S298- coeff*power*T**(power-1)
		H298= H298 +coeff*(1-power)* T**(power)
		
		if abs(power-1)<1E-5: #power-1=0
			Cp[0]= Cp[0]+coeff*power*(1-power)
		elif abs(power-2)<1E-5:#power-1=1
			Cp[1]= Cp[1]+coeff*power*(1-power)
		elif abs(power+1)<1E-5: #power-1=-2
			Cp[2]= Cp[2]+coeff*power*(1-power)
		elif abs(power-3)<1E-5: #power-1=2
			Cp[3]= Cp[3]+coeff*power*(1-power)
		else:
			j = int(2*i+4)
			Cp[j]= +coeff*power*(1-power)
			Cp[j+1]= power-1
	return H298,S298,Cp

def HSCp2G(H298: float, S298: float, Cp: Union[List[List[float]], np.ndarray]) -> np.ndarray:
	'''
	Convert H298 [J], S298 [J/(mol K)], Cp(T) [J/(mol K)] to Gibbs energies [J]
	Cp can have multiple temperature intervals
	
	Input
	-----
	H298 : Enthalpy of formation at 298.15 K
	S298 : Vibration entropy at 298.15 K
	Cp(T): A list of 14 floats [[Tmin, Tmax, a, b, c, d, c1, p1, c2, p2, c3, p3, c4, p4],
								[Tmin, Tmax, a, b, c, d, c1, p1, c2, p2, c3, p3, c4, p4],
								...]
		  Cp(T) = a +b*T +c*T**(-2) +d*T**(2) +c1*T**(p1) +c2*T**(p2) +c3*T**(p3) +c4*T**(p4) 

	Output
	------
	GibbsEnegies : A list of 16 floats [[Tmin, Tmax, a, b, c, d, e, f, c1, p1, c2, p2, c3, p3, c4, p4],
										[Tmin, Tmax, a, b, c, d, e, f, c1, p1, c2, p2, c3, p3, c4, p4],
										...]
				a-f are coefficients of predefined function, ci and pi are coefficient and power of custom functions
				G = a +b*T +c*T*s.ln(T) +d*T**2 +e*T**3 +f*T**(-1) +c1*T**(p1) +c2*T**(p2) +c3*T**(p3) +c4*T**(p4) 
	
	'''
	Cp=np.asarray(Cp)
	n,l = Cp.shape # number of temperature intervals
	
	assert l==14, "Error: input must have length of 14"
	
	GibbsEnergies= np.zeros((n,16))
	

	for i, cp in enumerate(Cp):
		GibbsEnergies[i,0] = cp[0] #minimum temperature
		GibbsEnergies[i,1] = cp[1] #maximumm temperature
		
		# Allocate temperary variables
		a, b, c, d, c1, p1, c2, p2, c3, p3, c4, p4 = cp[2:]
		custom= [[c1,p1],[c2,p2],[c3,p3],[c4,p4]]
				
		
		# IntCp_H298= (a*T + 0.5*b*T**2 +(d/3)*T**3 - c/T)
		# IntCp_S298= (a*np.log(T) +b*T + 0.5*d*T**2- 0.5*c/(T**2))
		IntCp_H= lambda T: _IntCp_H(cp[2:],T)
		IntCp_S= lambda T: _IntCp_S(cp[2:],T)
		
		for k, param in enumerate(custom):
			coeff, power= param
			
			if abs(power)<1E-50: continue
		
			if abs(power+1)<1E-25:
				print('HSCp2G warning: customized coefficient and power are not considered as it appears as -1')
				continue
			# IntCp_H298 = IntCp_H298+(coeff/(power + 1))*T**(power + 1)
			# IntCp_S298 = IntCp_S298 + coeff/(power)*T**(power)
			
			CoeffCustom= (coeff/(power + 1)-coeff/(power))
			if  abs(power-1)<1E-25: #power+1=2
				GibbsEnergies[i,5]= CoeffCustom
			elif abs(power-2)<1E-25: #power+1=3
				GibbsEnergies[i,6]= CoeffCustom
			elif abs(power+2)<1E-25: #power+1=-1
				GibbsEnergies[i,7]= CoeffCustom
			else:
				j = int(2*k+8)
				GibbsEnergies[i,j]= CoeffCustom
				GibbsEnergies[i,j+1]= power+1
		
  
		# allocate Gibbs
		if i!=0:
			Tmin=GibbsEnergies[i,0]
			H298= H298 + IntCp_Hprev-IntCp_H(Tmin)+IntCp_H(298.15)
			S298= S298 + IntCp_Sprev-IntCp_S(Tmin)+IntCp_S(298.15)

		Tmax=GibbsEnergies[i,1]
		
		GibbsEnergies[i,2]= GibbsEnergies[i,2]+H298-IntCp_H(298.15)   # a term
		GibbsEnergies[i,3]= GibbsEnergies[i,3]-S298+a+IntCp_S(298.15) # b*T term
		GibbsEnergies[i,4]= GibbsEnergies[i,4]-a                 #TlnT, 
		GibbsEnergies[i,5]= GibbsEnergies[i,5]-0.5*b             #T**2, 
		GibbsEnergies[i,6]= GibbsEnergies[i,6]+ (1/3-0.5)*d      #T**3, 
		GibbsEnergies[i,7]= GibbsEnergies[i,7]-0.5*c             #T**(-1)
		
		IntCp_Hprev=IntCp_H(Tmax)-IntCp_H(298.15)
		IntCp_Sprev=IntCp_S(Tmax)-IntCp_S(298.15)
	
	# Ensure that the final temperature is 6000 K
	if GibbsEnergies[-1,1] <=6000: GibbsEnergies[-1,1]=6000
		
	return GibbsEnergies

def G2HSCp(GibbsEnegies: Union[List[List[float]], np.ndarray]) -> Tuple[float, float, np.ndarray]:
	'''
	Convert Gibbs energy coefficients to H298, S298, and Cp(T) over temperature intervals.

	This function takes Gibbs energy coefficients for one or more temperature intervals
	and calculates the standard enthalpy of formation (H298), standard entropy (S298),
	and the coefficients for the heat capacity (Cp) polynomial for each interval.

	The Gibbs energy is expressed as:
	G(T) = a +b*T +c*T*s.ln(T) +d*T**2 +e*T**3 +f*T**(-1) +c1*T**(p1) +c2*T**(p2) +c3*T**(p3) +c4*T**(p4) 


	The resulting heat capacity is expressed as:
	Cp(T) = a +b*T +c*T**(-2) +d*T**(2) +c1*T**(p1) +c2*T**(p2) +c3*T**(p3) +c4*T**(p4) 

	Parameters
	----------
	GibbsEnegies : Union[List[List[float]], np.ndarray]
		A list or array of Gibbs energy coefficients. Each row should contain
		[Tmin, Tmax, a, b, c, d, e, f, c1, p1, c2, p2, c3, p3, c4, p4].

	Returns
	-------
	H298 : float
		Enthalpy of formation at 298.15 K [J/mol].
	S298 : float
		Standard entropy at 298.15 K [J/(mol K)].
	Cp : np.ndarray
		An array representing the heat capacity over various temperature intervals.
		Each row contains [Tmin, Tmax, a', b', c', d', c1', p1', ...].
	'''
	
	# assert isinstance(Gibbs,list), "Error: input must be a list with length 14"
	GibbsEnegies=np.asarray(GibbsEnegies)
	n,l=GibbsEnegies.shape
	
	assert l==16, "Error: input must have length of 16"
	
	#Sort Gibbs energy by Tmin
	Tmins=GibbsEnegies[:,0]
	idx=np.argsort(Tmins)
	GibbsEnegies=GibbsEnegies[idx,:]
	
	Cp= np.zeros((n,14))
	Cp[:,0]=GibbsEnegies[:,0]
	Cp[:,1]=GibbsEnegies[:,1]
	for i in range(n):
		H,S,cp=G2HSCp_single(GibbsEnegies[i,2:])
		if i ==0:
			H298=H
			S298=S
		Cp[i,2:]=cp
		
	return H298, S298, Cp

def NeumanKoppHSCp(dataset: Dict[str, Dict[str, Union[float, np.ndarray]]], names: List[str], multipliers: Union[List[float], np.ndarray]) -> Dict[str, Union[float, np.ndarray]]:
	'''
	Calculates the thermodynamic properties of a mixture using the Neumann-Kopp rule.

	This function computes the enthalpy of formation (H298), standard entropy (S298),
	and heat capacity (Cp) for a mixture of substances. It combines the properties
	of individual components, weighted by their respective multipliers.

	Input
	-----
	dataset : dict
		A dictionary containing the thermodynamic data for various functions/species.
		Each entry should have 'H298', 'S298', and 'Cp' data.
	names : list of str
		A list of function/species names to be included in the calculation. These names
		must be keys in the `dataset` dictionary.
	multipliers : list or array of float
		The corresponding multipliers (e.g., mole fractions) for each species in `names`.

	Output
	------
	res : dict
		A dictionary containing the calculated thermodynamic properties for the mixture:
		- 'H298': Enthalpy of formation at 298.15 K [J/mol].
		- 'S298': Standard entropy at 298.15 K [J/(mol K)].
		- 'Cp': A numpy array representing the heat capacity over various temperature
		  intervals. Each row corresponds to an interval and contains coefficients
		  for the Cp polynomial.
	'''
	
	# Data must be dictionary
	headers=list(dataset.keys())
	multipliers=np.asarray(multipliers)
	nmultipliers=len(multipliers)
	
	#Memmory allocation
	T_intervals=[]
	H298s=[]
	S298s=[]
	pi_index=[5,7,9,11]
	ci_index=[4,6,8,10]
	#Ensure that the given names are within the dataset
	for name in names:
		assert name in headers, "Error: given function name is not part of the dataset"
		T_intervals=T_intervals+list(dataset[name]['Cp'][:,:2].flatten())
		H298s.append(dataset[name]['H298'])
		S298s.append(dataset[name]['S298'])
	T_intervals=np.asarray(T_intervals)
	T_intervals=np.unique(T_intervals.round(decimals=3))
	dT_intervals=T_intervals[1:]-T_intervals[:-1]
	nintervals = len(dT_intervals)
	H298_new= np.dot(H298s,multipliers)
	S298_new= np.dot(S298s,multipliers)
	
	
	for i in range(nintervals):
		Cp_new=np.zeros(14)
		Cp_temp=np.zeros((nmultipliers,12))
		custom_power=[]
		Tmax= T_intervals[i+1]
		Cp_new[0]=T_intervals[i]
		Cp_new[1]=Tmax
		
		for j, name in enumerate(names):
			Cps=dataset[name]['Cp']
			
			n_interval_temp=len(Cps)
			for k in range(len(Cps)):
				if Tmax>Cps[k,0] and Tmax<=Cps[k,1]:
					Cp_temp[j,:]=_sort_CustomCp(Cps[k,2:])
					custom_power.append(Cp_temp[j,pi_index])
					break
		Cp_new[2:6]=np.matmul(multipliers,Cp_temp[:,:4])
		custom_power = np.unique(custom_power)
		custom_power = custom_power[custom_power != 0]
		if len(custom_power)!=0:
			powers= Cp_temp[:,pi_index]
			for k,p in enumerate(custom_power):
				imultiplier,ipower=np.asarray(np.where(powers==p))
				lmultiplier=len(imultiplier)
				imultiplier=np.squeeze(np.unique(imultiplier))
				ipower=np.squeeze(np.unique(ipower))
				
				if lmultiplier==1:
					
					Cp_new[int(7+2*k)]=p
					Cp_new[int(6+2*k)]=multipliers[imultiplier]*Cp_temp[imultiplier,ci_index[ipower]]
				else:
					Cp_new[int(7+2*k)]=p
					Cp_new[int(6+2*k)]=np.matmul(multipliers[imultiplier],Cp_temp[imultiplier,ci_index[ipower]])
		if i==0:
			Cps_new=np.zeros((1,14))
			Cps_new[0,:]=Cp_new
		else:
			Cps_new=np.vstack((Cps_new,Cp_new))
				
			
	res = {
		'H298':H298_new,
		'S298':S298_new,
		'Cp': Cps_new
	}
	return res


# Internal functions >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
def _IntCp_H(cp: np.ndarray, T: float) -> float:
	'''
	Calculates the integral of heat capacity (Cp) to determine the change in enthalpy (H).

	This function analytically integrates the heat capacity polynomial for a given temperature T,
	which is used to compute the enthalpy change.

	Parameters
	----------
	cp : np.ndarray
		An array of coefficients for the heat capacity polynomial.
	T : float
		The temperature (in Kelvin) at which to calculate the enthalpy.

	Returns
	-------
	float
		The calculated enthalpy change at the given temperature.
	'''
	a, b, c, d, c1, p1, c2, p2, c3, p3, c4, p4 = cp.copy()
	custom= [[c1,p1],[c2,p2],[c3,p3],[c4,p4]]
	res=(a*T + 0.5*b*T**2 +(d/3)*T**3 - c/T)

	for k, param in enumerate(custom):
			coeff, power= param
			if abs(coeff)<1E-50: continue
			
			if abs(power)<1E-25: continue
		
			if abs(power+1)<1E-25:
				print('_IntCp_H warning: customized coefficient and power are not considered as it appears as -1')
				continue
			res = res+(coeff/(power + 1))*np.power(T,power + 1)
	return res

def _IntCp_S(cp: np.ndarray, T: float) -> float:
	'''
	Calculates the integral of heat capacity over temperature (Cp/T) to determine the change in entropy (S).

	This function analytically integrates the heat capacity polynomial divided by temperature
	for a given temperature T, which is used to compute the entropy change.

	Parameters
	----------
	cp : np.ndarray
		An array of coefficients for the heat capacity polynomial.
	T : float
		The temperature (in Kelvin) at which to calculate the entropy.

	Returns
	-------
	float
		The calculated entropy change at the given temperature.
	'''
	a, b, c, d, c1, p1, c2, p2, c3, p3, c4, p4 = cp.copy()
	custom= [[c1,p1],[c2,p2],[c3,p3],[c4,p4]]
	res=(a*np.log(T) +b*T + 0.5*d*T**2- 0.5*c/(T**2))
	for k, param in enumerate(custom):
			coeff, power= param
			if abs(coeff)<1E-50: continue
			
			if abs(power)<1E-25: continue
		
			if abs(power+1)<1E-25:
				print('_IntCp_S warning: customized coefficient and power are not considered as it appears as -1')
				continue
			res = res+coeff/(power)*T**(power)
	return res

def _sort_CustomCp(cp: np.ndarray) -> np.ndarray:
	'''
	Sorts the custom heat capacity (Cp) parameters by their power values.

	This function rearranges the coefficients and powers of the custom Cp terms
	in ascending order of their powers. This standardization is useful for
	consistent processing and combining of thermodynamic data.

	Parameters
	----------
	cp : np.ndarray
		An array of heat capacity coefficients, including custom terms.

	Returns
	-------
	np.ndarray
		The cp array with custom terms sorted by power.
	'''
	a, b, c, d, c1, p1, c2, p2, c3, p3, c4, p4 = cp
	custom= np.asarray([[c1,p1],[c2,p2],[c3,p3],[c4,p4]])
	
	idx=np.argsort(custom[:,1])
	res= cp
	res[[4,6,8,9]]=custom[idx,0]
	res[[5,7,9,11]]=custom[idx,1]
	
	return res

def _dict2np(dictionary: Dict[str, Union[Any, List[Any]]]) -> Tuple[List[str], np.ndarray]:
    '''
    Converts a dictionary to a NumPy array, ensuring consistent column lengths.

    This function takes a dictionary where keys are headers and values are data columns.
    It converts the dictionary into a NumPy array, handling cases where values are
    single items by wrapping them in a list. It raises a ValueError if the columns
    have different lengths.

    Parameters
    ----------
    dictionary : Dict[str, Union[Any, List[Any]]]
        The input dictionary to convert.

    Returns
    -------
    header : List[str]
        A list of the dictionary keys, representing the headers of the columns.
    res : np.ndarray
        A NumPy array containing the dictionary values, with shape (n_rows, n_columns).
    
    Raises
    ------
    ValueError
        If the columns in the dictionary have different lengths.
    '''
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


def _pyvar2fvar(var: object) -> None:
    '''
    Passes Python variables to Fortran modules.

    This function takes a Python object containing thermodynamic data and assigns its
    attributes to the corresponding variables in the Fortran modules. This is a
    critical step for interfacing with the Fortran backend.

    Parameters
    ----------
    var : object
        A Python object (typically a module or class instance) containing the variables
        to be passed to Fortran.
    '''
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
# Internal functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<