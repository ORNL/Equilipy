from dataclasses import dataclass, field, replace
import numpy as np
import equilipy.equilifort as fort
import equilipy.variables as var
from equilipy.InternalFunctions import _dict2np
from .Errors import *

def _count_unstable_compounds(i):    
    n=0
    iphase=fort.modulethermo.iphase
    iLastSoln=fort.modulethermo.nspeciesphase[len(var.iSys2DBSoln)]
    for k in range(iLastSoln,i):
        if iphase[k]<0:
            n = n+1
    return n

def _get_assemblage_name(AssemblageIDs):
    '''
    An internal function that spits out phase name based on Assemblage IDs

    Variables
    =========

    Input
    AssemblageIDs: list of phase assemblage id defined after selecting phases based on input elements

    Output
    AssemblageNames
    '''
    AssemblageNames=list([])
    for i in AssemblageIDs:
        if i ==0:
            # Empty phase: place holder
            AssemblageNames.append("{:<1}".format(''))
        elif i<0:
            # Solution Phase
            # Note that the solution phase id changes with phase selection
            AssemblageNames.append(var.cPhaseNameSys[int(-(i+1))].strip())
        else:
            #Compound Phase
            
            # Note that the compund id doesn't change with phase selection
            # Correction is made by counting the number of ignored compound phases
            nSoln=len(var.iSys2DBSoln)
            iLastSoln=fort.modulethermo.nspeciesphase[nSoln]
            id=-iLastSoln+i+nSoln-1
            n = _count_unstable_compounds(i)
            id = int(id-n)
            AssemblageNames.append(var.cPhaseNameSys[id].strip())
    return AssemblageNames

def check_output_units(UnitOut:list=['K','atm','moles']):
    
    # Convert Temperature Unit
    if UnitOut[0] == 'C':
        var.dTOut=fort.modulethermoio.dtemperature - 273.15
    elif UnitOut[0] == 'F':
        var.dTOut= fort.modulethermoio.dtemperature*9/5 -459.67
    elif UnitOut[0] == 'R':
        var.dTOut=fort.modulethermoio.dtemperature*9/5
    else:
        var.dTOut=fort.modulethermoio.dtemperature

    # Convert Pressure Unit
    if UnitOut[1] == 'psi':
        var.dPOut=fort.modulethermoio.dpressure/0.068045957064
    elif UnitOut[1] == 'bar':
        var.dPOut= fort.modulethermoio.dpressure/0.98692316931
    elif UnitOut[1] == 'Pa':
        var.dPOut=fort.modulethermoio.dpressure*1E5/0.98692316931
    elif UnitOut[1] == 'kPa':
        var.dPOut=fort.modulethermoio.dpressure*1E2/0.98692316931
    else: #atm
        var.dPOut=fort.modulethermoio.dpressure
       
    if UnitOut[2] in ['grams','kilograms','pounds','g','kg','lbs',
                    'mass fraction','weight fraction','wt%','wt.%']:
        var.dSpeciesFractionOut = fort.modulethermoio.dgramfraction.copy()
        var.dPhaseAmountOut = fort.modulethermoio.dgramphase.copy()
        
    else:
        var.dSpeciesFractionOut = fort.modulethermo.dmolfraction.copy()
        var.dPhaseAmountOut = fort.modulethermo.dmolesphase.copy()
    
    if UnitOut[2] in ['kilograms','kg']:
        var.dPhaseAmountOut = 1E-3*var.dPhaseAmountOut
    elif UnitOut[2] in ['pounds','lbs']:
        var.dPhaseAmountOut = var.dPhaseAmountOut/453.59237
    return

@dataclass
class Phase:    
    '''
    Phase object takes 'iSys' as a variable and assign partinent thermochemical properties. 

    Variables
    =========

    Input
    -----
    iSys : Element in the phase index size of (len(iSys2DBSoln)+len(iSys2DBComp)) that is defined after selecting phases
    based on input list of elements, so the species/solution numbers differ from what is defined in the database.
    
    The related variables are listed below
    -cPhaseNameSys: solution and compound phase names defined by iSys index (Soln index:-int(iSys+1),Compd index: iSys)
    -iSys2DBSoln: List that relates Soln index defined in the system to database index
    -iSys2DBSpecies:List that relates Species index defined in the system to database index
    -iSys2DBComp: List that relates Compd index defined in the system to database index

    Output
    ------
    Phase properties: ID, Name, Amount, Stability (1:stable, 0: not stable)
    '''
    def __init__(self,iSys):
        self.ID: int =field(default_factory=int)
        self.Name: str =field(default_factory=str)
        self.Amount: float = field(default_factory=float)
        self.Stability: float = field(default_factory=float)
        self.Endmembers: list[str] = field(default_factory=list)
        self.xi: list[float] = field(default_factory=list)
        
        
        # Assign phase name
        self.Name = str(var.cPhaseNameSys[iSys])
        
        if iSys<len(var.iSys2DBSoln):
            # Solution Phase:
            self.ID =-(iSys+1)
            if self.ID in list(fort.modulethermo.iassemblage):
                # This phase is stable. Assign the relevant amount.
                self.Amount = float(var.dPhaseAmountOut[list(fort.modulethermo.iassemblage).index(self.ID)])
                self.Stability = 1.0
            else:
                # This phase is not stable. Assign the amount to zero.
                self.Amount = 0.0
                self.Stability = 0.0
            
            # Obtain endmembers index
            iFirstSys = fort.modulethermo.nspeciesphase[iSys]
            iLastSys= fort.modulethermo.nspeciesphase[iSys+1]
            
            idx_species=np.array(var.iSys2DBSpecies[iFirstSys:iLastSys])
            
            # Assign endmemebr properties
            self.Endmembers=list([str(x).strip() for x in list(np.array(var.cSpeciesNameCS)[idx_species])])
            self.xi=list(var.dSpeciesFractionOut[iFirstSys:iLastSys])
            #     ai: list[float] = field(default_factory=list)
            #     gi: list[float] = field(default_factory=list)
            #     hi: list[float] = field(default_factory=list)
            #     si: list[float] = field(default_factory=list)
        else:
            # Compound phase (iSys>len(var.iSys2DBSoln))
            self.ID = fort.modulethermo.nspeciesphase[len(var.iSys2DBSoln)] +iSys-len(var.iSys2DBSoln)+1
            if self.ID in list(fort.modulethermo.iassemblage):
                # This phase is stable. Assign the relevant amount.
                self.Amount = float(var.dPhaseAmountOut[list(fort.modulethermo.iassemblage).index(self.ID)])
                self.Stability = 1.0
            else:
                # This phase is not stable. Assign the amount to zero.
                self.Amount = 0.0
                self.Stability = 0.0
            
            # Compound phase do not have endmembers. Assign the species name instead
            self.Endmembers=list([var.cPhaseNameSys[iSys]])
            self.xi=list([1])
        return

@dataclass
class EmptyPhase:    
    def __init__(self):
        self.ID: int =field(default_factory=int)
        self.Name: str =field(default_factory=str)
        self.Amount: float = field(default_factory=float)
        self.Stability: float = field(default_factory=float)
        self.Endmembers: list[str] = field(default_factory=list)
        self.xi: list[float] = field(default_factory=list)
        return

@dataclass
class Result:
    '''
    Result object that includes system and phase properties

    Output
    ======

    System properties
    -----------------
    N          : Dictionary {element: Amount,...}
    T          : Temperature
    P          : Pressure
    G          : System Gibbs energy (total)
    H          : System enthalpy (total)
    S          : System entropy (total)
    Cp         : System heat capacity (total)
    

    Phase properties
    ----------------
    StablePhases: Dictionary {'ID':[phase_id1,..],'Name':[phase_name1,...],'Amount':[phase_amount1,...]'}
    Phases      : Dictionary {phase_name: phase_object,...}
    '''
   # System properties
    N: dict = field(default_factory=dict)
    T: float = None
    P: float = None
    G: float = field(default_factory=float)
    H: float = field(default_factory=float)
    S: float = field(default_factory=float)
    Cp: float = field(default_factory=float)
    
    # Phase properties
    StablePhases:dict = field(default_factory=dict)
    Phases :dict = field(default_factory=dict)
    
    def append(self,other:object):
        '''
        A method to append calculation results
        ''' 
        # When there is more than two conditions, change results type to list for the first calcuation
        if type(other.T)==float: lSinglePoint=True
        else: lSinglePoint=False
        
        if type(self.T)==float:
            for el in list(self.N.keys()):
                self.N[el]=list([self.N[el]])
            self.T = list([self.T])
            self.P = list([self.P])
            self.G = list([self.G])
            self.H = list([self.H])
            self.S = list([self.S])
            self.Cp = list([self.Cp])
            
            self.StablePhases['ID'] = list([self.StablePhases['ID']])
            self.StablePhases['Name'] = list([self.StablePhases['Name']])
            self.StablePhases['Amount'] = list([self.StablePhases['Amount']])
            
            for name in list(self.Phases.keys()):
                # self.Phases[name].Name=name
                self.Phases[name].Amount=list([self.Phases[name].Amount])
                self.Phases[name].Stability=list([self.Phases[name].Stability])
                self.Phases[name].xi=list([self.Phases[name].xi])
                self.Phases[name].Endmembers=list([self.Phases[name].Endmembers])
        
        # Append system and phase properties
        ElementsBefore= list(self.N.keys())
        ln_before= int(len(self.N[ElementsBefore[0]]))
        for i,el in enumerate(list(other.N.keys())):
            if el not in ElementsBefore: 
                #Add padding
                self.N[el] = list([float(0.0)]*ln_before)
            if lSinglePoint: self.N[el].append(other.N[el])
            else: self.N[el]=self.N[el]+other.N[el]
            ln_after = int(len(self.N[el]))
        
        
        if lSinglePoint:
            self.T.append(other.T)
            self.P.append(other.P)
            self.G.append(other.G)
            self.H.append(other.H)
            self.S.append(other.S)
            self.Cp.append(other.Cp)
            self.StablePhases['ID'].append(other.StablePhases['ID'])
            self.StablePhases['Name'].append(other.StablePhases['Name'])
            self.StablePhases['Amount'].append(other.StablePhases['Amount'])
        else:
            self.T = self.T + other.T
            self.P = self.P + other.P
            self.G = self.G + other.G
            self.H = self.H + other.H
            self.S = self.S + other.S
            self.Cp = self.Cp + other.Cp
            self.StablePhases['ID']=self.StablePhases['ID'] + other.StablePhases['ID']
            self.StablePhases['Name']=self.StablePhases['Name'] + other.StablePhases['Name']
            self.StablePhases['Amount']=self.StablePhases['Amount'] + other.StablePhases['Amount']
        
        # Check through elements and see if we need to add paddings
        for el in list(self.N.keys()):
            ln_differ= ln_after-len(self.N[el])
            if ln_differ>0: self.N[el]=self.N[el]+list([float(0.0)]*ln_differ)
        
        PhasesBefore=list(self.Phases.keys())
        
        lp_after= int(len(self.T))
        if lSinglePoint: 
            lp_before = lp_after-1
        else:
            lp_before = lp_after-len(other.T)
        
        for Pname in list(other.Phases.keys()):
            new=other.Phases[Pname]
            
            # In case it has not appeared previously add padding
            if Pname in PhasesBefore:
                if lSinglePoint:
                    self.Phases[Pname].Amount.append(new.Amount)
                    self.Phases[Pname].Stability.append(new.Stability)
                    self.Phases[Pname].xi.append(new.xi)
                    self.Phases[Pname].Endmembers.append(new.Endmembers)
                else:
                    self.Phases[Pname].Amount=self.Phases[Pname].Amount + new.Amount
                    self.Phases[Pname].Stability=self.Phases[Pname].Stability + new.Stability
                    self.Phases[Pname].xi=self.Phases[Pname].xi + new.xi
                    self.Phases[Pname].Endmembers=self.Phases[Pname].Endmembers + new.Endmembers
            else:
                # Assign Nan for all previous calculations
                self.Phases[Pname] = EmptyPhase()
                self.Phases[Pname].ID=new.ID #Initialize
                self.Phases[Pname].Name = Pname
                self.Phases[Pname].Amount=list([float(0.0)]*lp_before)
                self.Phases[Pname].Stability=list([float(0.0)]*lp_before)
                self.Phases[Pname].xi=list([[int(0) for i in range(len(new.Endmembers))]])*lp_before
                self.Phases[Pname].Endmembers=list([[str('nan') for i in range(len(new.Endmembers))]])*lp_before
                
                # Assign values for current condition
                if lSinglePoint:
                    self.Phases[Pname].Amount.append(new.Amount)
                    self.Phases[Pname].Stability.append(new.Stability)
                    self.Phases[Pname].xi.append(new.xi)
                    self.Phases[Pname].Endmembers.append(new.Endmembers)
                else:
                    self.Phases[Pname].Amount=self.Phases[Pname].Amount + new.Amount
                    self.Phases[Pname].Stability=self.Phases[Pname].Stability + new.Stability
                    self.Phases[Pname].xi=self.Phases[Pname].xi + new.xi
                    self.Phases[Pname].Endmembers=self.Phases[Pname].Endmembers + new.Endmembers

        # Check through all phases and add paddings
        for name in list(self.Phases.keys()):
            lp_differ= lp_after-len(self.Phases[name].Amount)
            
            if lp_differ>0:
                nEndmembers=len(self.Phases[name].Endmembers[-1])
                self.Phases[name].Amount=self.Phases[name].Amount+list([float(0.0)]*lp_differ)
                self.Phases[name].Stability=self.Phases[name].Stability+list([float(0.0)]*lp_differ)
                self.Phases[name].xi=self.Phases[name].xi+([[int(0) for i in range(nEndmembers)]])*lp_differ
                self.Phases[name].Endmembers=self.Phases[name].Endmembers+list([[str('nan') for i in range(nEndmembers)]])*lp_differ

        return None
        
    
    def append_output(self):
        '''
        A method to append calculation results
        '''
        # 1 Check if it is first time assigning variables
        if (self.T==None) and (self.P == None):            
            # Assign system properties
            self.N = dict(zip(var.cComponentNameSys,var.dConditionSys[2:]))
            self.T = float(var.dTOut)
            self.P = float(var.dPOut)
            self.G = float(fort.modulethermoio.dgibbsenergysys)
            # # To be modified
            self.H = float(fort.modulethermoio.denthalpysys)
            self.S = float(fort.modulethermoio.dentropysys)
            self.Cp = float(fort.modulethermoio.dheatcapacitysys)
            
            # Assign stable phase properties
            self.StablePhases['ID'] = list(fort.modulethermo.iassemblage)
            self.StablePhases['Name'] = list(_get_assemblage_name(fort.modulethermo.iassemblage))
            self.StablePhases['Amount'] = list(var.dPhaseAmountOut)
            
            # Assign phase properties
            var.iPhaseSys = list([])
            
            for i in range(len(var.iSys2DBSoln)+len(var.iSys2DBComp)):
                if i<len(var.iSys2DBSoln):
                    id=int(-(i+1))
                else:
                    id=fort.modulethermo.nspeciesphase[len(var.iSys2DBSoln)]\
                    +i-len(var.iSys2DBSoln)+1
                new=Phase(i)
                name=new.Name.strip()
                self.Phases[name]=new
                var.iPhaseSys.append(id)
                
        else:
            # When there is more than two conditions, change results type to list for the first calcuation
            if (type(self.T)==float) and (type(self.P)==float):
                for el in list(self.N.keys()):
                    self.N[el]=list([self.N[el]])
                self.T = list([self.T])
                self.P = list([self.P])
                self.G = list([self.G])
                self.H = list([self.H])
                self.S = list([self.S])
                self.Cp = list([self.Cp])
                
                self.StablePhases['ID'] = list([self.StablePhases['ID']])
                self.StablePhases['Name'] = list([self.StablePhases['Name']])
                self.StablePhases['Amount'] = list([self.StablePhases['Amount']])
                
                for name in list(self.Phases.keys()):
                    # self.Phases[name].Name=name
                    self.Phases[name].Amount=list([self.Phases[name].Amount])
                    self.Phases[name].Stability=list([self.Phases[name].Stability])
                    self.Phases[name].xi=list([self.Phases[name].xi])
                    self.Phases[name].Endmembers=list([self.Phases[name].Endmembers])
            
            # Append system and phase properties
            ElementsBefore= list(self.N.keys())
            ln_before= int(len(self.N[ElementsBefore[0]]))
            for i,el in enumerate(var.cComponentNameSys):
                if el not in ElementsBefore: 
                    #Add padding
                    self.N[el] = list([float(0.0)]*ln_before)
                self.N[el].append(var.dConditionSys[2+i])
                ln_after = int(len(self.N[el]))
            self.T.append(float(var.dTOut))
            self.P.append(float(var.dPOut))
            self.G.append(float(fort.modulethermoio.dgibbsenergysys))
            self.H.append(float(fort.modulethermoio.denthalpysys))
            self.S.append(float(fort.modulethermoio.dentropysys))
            self.Cp.append(float(fort.modulethermoio.dheatcapacitysys))
            self.StablePhases['ID'].append(list(fort.modulethermo.iassemblage))
            self.StablePhases['Name'].append(list(_get_assemblage_name(fort.modulethermo.iassemblage)))
            self.StablePhases['Amount'].append(list(var.dPhaseAmountOut))
            
            # Check through elements and see if we need to add paddings
            for el in list(self.N.keys()):
                ln_differ= ln_after-len(self.N[el])
                if ln_differ>0: self.N[el]=self.N[el]+list([float(0.0)]*ln_differ)
            
            PhasesBefore=list(self.Phases.keys())
            # print('PhasesBefore',PhasesBefore)
            # if len(PhasesBefore)==0: lp_before = 0
            # else: lp_before= int(len(self.Phases[PhasesBefore[0]].Amount))
            lp_before = int(len(self.T))-1
            for i in range(len(var.iSys2DBSoln)+len(var.iSys2DBComp)):
                new=Phase(i)
                name=new.Name.strip()
                # In case it has not appeared previously add padding
                if name in PhasesBefore:
                    self.Phases[name].Amount.append(new.Amount)
                    self.Phases[name].Stability.append(new.Stability)
                    self.Phases[name].xi.append(new.xi)
                    self.Phases[name].Endmembers.append(new.Endmembers)
                    
                else:
                    # Assign Nan for all previous calculations
                    self.Phases[name]=EmptyPhase()
                    self.Phases[name].ID=new.ID
                    self.Phases[name].Name = name
                    self.Phases[name].Amount=list([float(0.0)]*lp_before)
                    self.Phases[name].Stability=list([float(0.0)]*lp_before)
                    self.Phases[name].xi=list([[int(0) for i in range(len(new.Endmembers))]])*lp_before
                    self.Phases[name].Endmembers=list([[str('nan') for i in range(len(new.Endmembers))]])*lp_before
                    
                    # Assign values for current condition
                    self.Phases[name].Amount.append(new.Amount)
                    self.Phases[name].Stability.append(new.Stability)
                    self.Phases[name].xi.append(new.xi)
                    self.Phases[name].Endmembers.append(new.Endmembers)
                lp_after= int(len(self.Phases[name].Amount))
                
            # Check through all phases and add paddings
            for name in list(self.Phases.keys()):
                lp_differ= lp_after-len(self.Phases[name].Amount)
                
                if lp_differ>0:
                    nEndmembers=len(self.Phases[name].Endmembers[-1])
                    self.Phases[name].Amount=self.Phases[name].Amount+list([float(0.0)]*lp_differ)
                    self.Phases[name].Stability=self.Phases[name].Stability+list([float(0.0)]*lp_differ)
                    self.Phases[name].xi=self.Phases[name].xi+([[int(0) for i in range(nEndmembers)]])*lp_differ
                    self.Phases[name].Endmembers=self.Phases[name].Endmembers+list([[str('nan') for i in range(nEndmembers)]])*lp_differ
        return None
    
    def append_error(self):
        '''
        A method to append calculation results
        '''
        # 1 Check if it is first time assigning variables
        if (self.T==None) and (self.P == None):            
            # Assign system properties
            self.N = dict(zip(var.cComponentNameSys,var.dConditionSys[2:]))
            self.T = float(var.dTOut)
            self.P = float(var.dPOut)
            self.G = np.nan
            self.H = np.nan
            self.S = np.nan
            self.Cp = np.nan
            
            # Assign stable phase properties
            self.StablePhases['ID'] = list([np.nan])
            self.StablePhases['Name'] = list(['nan'])
            self.StablePhases['Amount'] = list([np.nan])
            
            # Assign phase properties
            var.iPhaseSys = list([np.nan])
                
        else:
            # When there is more than two conditions, change results type to list for the first calcuation
            if (type(self.T)==float) and (type(self.P)==float):
                for el in list(self.N.keys()):
                    self.N[el]=list([self.N[el]])
                self.T = list([self.T])
                self.P = list([self.P])
                self.G = list([np.nan])
                self.H = list([np.nan])
                self.S = list([np.nan])
                self.Cp = list([np.nan])
                
                self.StablePhases['ID'] = list([np.nan])
                self.StablePhases['Name'] = list(['nan'])
                self.StablePhases['Amount'] = list([np.nan])
                
            
            # Append system and phase properties
            ElementsBefore= list(self.N.keys())
            ln_before= int(len(self.N[ElementsBefore[0]]))
            for i,el in enumerate(var.cComponentNameSys):
                if el not in ElementsBefore: 
                    #Add padding
                    self.N[el] = list([float(0.0)]*ln_before)
                self.N[el].append(var.dConditionSys[2+i])
                ln_after = int(len(self.N[el]))
            self.T.append(float(var.dTOut))
            self.P.append(float(var.dPOut))
            self.G.append(np.nan)
            self.H.append(np.nan)
            self.S.append(np.nan)
            self.Cp.append(np.nan)
            
            self.StablePhases['ID'].append(np.nan)
            self.StablePhases['Name'].append('nan')
            self.StablePhases['Amount'].append(np.nan)
            
            # Check through elements and see if we need to add paddings
            for el in list(self.N.keys()):
                ln_differ= ln_after-len(self.N[el])
                if ln_differ>0: self.N[el]=self.N[el]+list([float(0.0)]*ln_differ)
            
            
            # Check through all phases and add paddings
            for name in list(self.Phases.keys()):
                lp_differ= int(len(self.T))-len(self.Phases[name].Amount)
                
                if lp_differ>0:
                    nEndmembers=len(self.Phases[name].Endmembers[-1])
                    self.Phases[name].Amount=self.Phases[name].Amount+list([float(0.0)]*lp_differ)
                    self.Phases[name].Stability=self.Phases[name].Stability+list([float(0.0)]*lp_differ)
                    self.Phases[name].xi=self.Phases[name].xi+([[int(0) for i in range(nEndmembers)]])*lp_differ
                    self.Phases[name].Endmembers=self.Phases[name].Endmembers+list([[str('nan') for i in range(nEndmembers)]])*lp_differ
        return None
    
    def to_dict(self):
        df = dict({})
        df['T']=self.T
        df['P']=self.P
        df.update(self.N)
        df['G J']=self.G
        df['H J']=self.H
        df['S J/K']=self.S
        df['Cp J/K']=self.Cp
        df['StablePhaseNames']=[str(x) for x in self.StablePhases['Name']]
        df['StablePhaseIDs']=[str(x) for x in self.StablePhases['ID']]
        df['StablePhaseAmount']=[str(x) for x in self.StablePhases['Amount']]
        for name in list(self.Phases.keys()):
            df[f'{name}_Amount']=list(self.Phases[name].Amount)
            df[f'{name}_Stability']=list(self.Phases[name].Stability)
            df[f'{name}_Endmemebrs']=list([f'{x}' for x in self.Phases[name].Endmembers])
            df[f'{name}_xi']=list([f'{x}' for x in self.Phases[name].xi])
        return df
    
@dataclass
class ResultScheil():
    '''
    Scheil-Gulliver cooling result object that includes system and phase properties

    Output
    ======

    System properties
    -----------------
    Components : System Elements
    N          : Dictionary {element: Amount,...}
    T          : Temperature
    P          : Pressure
    G          : System Gibbs energy (total)
    H          : System enthalpy (total)
    S          : System entropy (total)
    Cp         : System heat capacity (total)
    
    
    Scheil system properties (function of T)
    ----------------------------
    T, fl, fs, Precipitating Phases, Cumulative phases, Phases      : Dictionary {phase_name: phase_object,...}
    
    Scheil constituent properties (Final cooling)
    ----------------------------
    ScheilConstituents: dict({phase_config:amount})
    Final amount of each phases: dict({phase_name: amount})
    
    '''
    # # Condition variables
    # T: float = None
    # PrimaryPhases: list[str] = field(default_factory=list)
    # Components: list[str] = field(default_factory=list)
    # PhaseName: str = field(default_factory=str)
    # PrecipitatePhases: list[str] = field(default_factory=list)
    
     # System properties
    N: dict = field(default_factory=dict)
    T: list[float] = field(default_factory=float)
    P: list[float] = field(default_factory=float)
    G: list[float] = field(default_factory=float)
    fl: list[float] = field(default_factory=float)
    fs: list[float] = field(default_factory=float)
    H: list[float] = field(default_factory=list)
    S: list[float] = field(default_factory=list)
    # Cp: list[float] = field(default_factory=list)
    
    # Phase properties
    PhaseLabel= list()
    ScheilPhases= dict({})
    ScheilConstituents= dict({})
    EquilibResult=Result()
    
    
    def __init__(self):
        # Store system properties
        self.EquilibResult.append_output()
        self.N = self.EquilibResult.N
        self.T = self.EquilibResult.T
        self.P = self.EquilibResult.P
        self.G = self.EquilibResult.G
        self.H = self.EquilibResult.H
        self.S = self.EquilibResult.S
        
        # Store predicted stable phases
        CurrentPhases= self.EquilibResult.StablePhases['Name']
        for i, name in enumerate(CurrentPhases):
            if name == ' ': 
                continue
            if 'LIQ' in name.upper():
                self.ScheilPhases[name]= list([float(self.EquilibResult.StablePhases['Amount'][i])])
                self.fl=list([float(self.EquilibResult.StablePhases['Amount'][i])])
                self.fs=list([1-float(self.EquilibResult.StablePhases['Amount'][i])])
        
        # Update PhaseLabel (precipitating phases)
        CurrentPhases= [x for x in CurrentPhases if x != ' ']
        if len(CurrentPhases)==0: raise EquilibError('No phase appears stable during Scheil simulation')
        self.PhaseLabel.append('+'.join(sorted(CurrentPhases)))
        fort.resetthermo() 
        return None
    
    def append_output(self):
        self.EquilibResult.append_output()
        self.N = self.EquilibResult.N
        self.T = self.EquilibResult.T
        self.P = self.EquilibResult.P
        self.G = self.EquilibResult.G
        self.H = self.EquilibResult.H
        self.S = self.EquilibResult.S
        
        nrows = len(self.EquilibResult.T)
        
        # Update ScheilPhases based on stable phases
        CurrentPhases= self.EquilibResult.StablePhases['Name'][-1]
        
        PreviousPhases= list(self.ScheilPhases.keys())
        for i, name in enumerate(CurrentPhases):
            if name == ' ': 
                continue
            elif name in PreviousPhases:
                if 'LIQ' in name.upper():
                    self.ScheilPhases[name].append(float(self.EquilibResult.StablePhases['Amount'][-1][i]))
                    self.fl.append(float(self.EquilibResult.StablePhases['Amount'][-1][i]))
                    self.fs.append(1-float(self.EquilibResult.StablePhases['Amount'][-1][i]))
                else:
                    vals= self.ScheilPhases[name]
                    self.ScheilPhases[name].append(vals[-1]+float(self.EquilibResult.StablePhases['Amount'][-1][i]))
            else:
                self.ScheilPhases[name]= list([float(0.0)]*int(nrows-1))
                self.ScheilPhases[name].append(float(self.EquilibResult.StablePhases['Amount'][-1][i]))
        
        # Update PhaseLabel (precipitating phases)
        CurrentPhases= [x for x in CurrentPhases if x != ' ']
        if len(CurrentPhases)==0: raise EquilibError('No phase appears stable during Scheil simulation')
        self.PhaseLabel.append('+'.join(sorted(CurrentPhases)))
        
        # Add padding if a phase appeared stable previously but not at current condition     
        for i, name in enumerate(PreviousPhases):
            if name not in CurrentPhases:
                if 'LIQ' in name.upper():
                    # End of Scheil
                    self.ScheilPhases[name].append(0.0)
                    self.fl.append(float(self.EquilibResult.StablePhases['Amount'][-1][i]))
                    self.fs.append(1-float(self.EquilibResult.StablePhases['Amount'][-1][i]))
                else:
                    vals= self.ScheilPhases[name]
                    self.ScheilPhases[name].append(vals[-1])
        fort.resetthermo()
    
    def update_scheilconstituents(self):
        k=0
        for i,T in enumerate(self.T):
            constituents=self.PhaseLabel[i].split('+')
            constituents=[x for x in constituents if 'LIQ' not in x.upper() and x != '']
            
            if self.PhaseLabel[i]!=self.PhaseLabel[i-1] and k==0:
                # 1. first time non liquid phase starts appearing
                cphase='+'.join(constituents)
                k=i
            elif self.PhaseLabel[i]!=self.PhaseLabel[i-1] and k>0:
                # 2. Store the data whenever different scheilconstituents appears
                scheilconstituents=list(self.ScheilConstituents.keys())
                if cphase in scheilconstituents:
                    self.ScheilConstituents[cphase]=self.ScheilConstituents[cphase]+self.fl[k]-self.fl[i]
                else:
                    self.ScheilConstituents[cphase]=self.fl[k]-self.fl[i]
                k=i
                if len(constituents)>0:
                    cphase='+'.join(constituents)
                else:
                    # Only liquid phase appears, which is weird. Store into previous assemblage
                    constituents=self.PhaseLabel[i-1].split('+')
                    constituents=[x for x in constituents if 'LIQ' not in x.upper() and x != '']
                    cphase='+'.join(constituents)
            if i==len(self.T)-1:
                fphase='+'.join(constituents)
                scheilconstituents=list(self.ScheilConstituents.keys())
                if cphase != fphase:
                    # Since the scheilconstituents are different it must have went through 2.
                    # Only save the last part.
                    if fphase !='':
                        if fphase in scheilconstituents:
                            self.ScheilConstituents[fphase]=self.ScheilConstituents[fphase]+self.fl[k]
                        else:
                            self.ScheilConstituents[fphase]=self.fl[k]
                    
                else:
                
                    if cphase in scheilconstituents:
                        self.ScheilConstituents[cphase]=self.ScheilConstituents[cphase]+self.fl[k]
                    else:
                        self.ScheilConstituents[cphase]=self.fl[k]
        return None
    
    def to_dict(self):
        df = dict({})
        df['T']=self.T
        df['P']=self.P
        df.update(self.N)
        df['G J']=self.G
        df['H J']=self.H
        df['S J']=self.S
        df['fl mol']=self.fl
        df['fs mol']=self.fs
        df['Label']=self.PhaseLabel
        df.update(self.ScheilPhases)
        return df