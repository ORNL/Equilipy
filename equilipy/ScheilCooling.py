#!/usr/bin/env python3
import numpy as np, sys
from tqdm import tqdm
import equilipy.equilifort as fort
from .PostProcess import *
from .EquilibSingle import _equilib_single
from .Errors import *
def scheil_cooling(LiquidPhaseName:str,Database,Condition,Unit:list=['K','atm','moles'],dT=20.0,IterMax=5000,ListOfPhases=None,Progress=None,Output=None):
    # Calculate equilib
    _equilib_single(Database,Condition,Unit,ListOfPhases=ListOfPhases)
    
    
    if Progress is None and Output is None:
        # Initialize the progress bar on the first call
        IterStart=int(np.min([int(1000/dT),IterMax]))
        Progress = tqdm(total=IterStart,colour='#8060ff',ascii="░▒▓",bar_format="{l_bar}{bar:50}{r_bar}")
        
        # Set recursion limit
        sys.setrecursionlimit(IterStart*2)
        
        # Initialize Scheil cooling result object
        res = ResultScheil()
        StablePhaseNames = list(res.ScheilPhases.keys())
        
        # Check if the target phase (liquid) is in the stable phase
        if LiquidPhaseName in StablePhaseNames:
            #Get liquid composition
            ni=np.array(res.EquilibResult.Phases[LiquidPhaseName].xi)*res.ScheilPhases[LiquidPhaseName]
            components=res.EquilibResult.Phases[LiquidPhaseName].Endmembers
            
            Condition['T']=np.squeeze(Condition['T'])-dT
            for i, el in enumerate(components):
                Condition[el] = ni[i]
                Progress.update(1)
            scheil_cooling(LiquidPhaseName,Database,Condition,Unit=Unit,ListOfPhases=ListOfPhases,Output=res,dT=dT,Progress=Progress)
    else:
        res=Output
        res.append_output()
        # Update current results
        CurrentPhases=res.EquilibResult.StablePhases['Name'][-1]
        Liq_diff=res.EquilibResult.Phases[LiquidPhaseName].Amount[-2]-res.EquilibResult.Phases[LiquidPhaseName].Amount[-1]
        
        if LiquidPhaseName not in CurrentPhases or res.EquilibResult.Phases[LiquidPhaseName].Amount[-1]<1E-5:
            # When liquid is no longer stable or the amount is close to zero,
            # Base 1 of recursive: Break recursive when liquid phase is not stable
            Progress.total = Progress.n
            Progress.close()
            res.update_scheilconstituents()
            return res
        elif res.EquilibResult.Phases[LiquidPhaseName].Amount[-1]<1E-2 and (Liq_diff<1E-10 or Condition['T']<298.15):
            #Base 2 of recursive: Liquid amount sufficiently small and it's no longer reducing
            Progress.total = Progress.n
            Progress.close()
            res.update_scheilconstituents()
            return res
        elif Condition['T']<=0 and 'K' in Unit[0]:
            #Base 2 of recursive: Liquid amount sufficiently small and it's no longer reducing
            
            Progress.total = Progress.n
            Progress.close()
            raise EquilibError('Scheil-Gulliver solidification results in melting point below 0K. Check phase selections.')
        else:
            #Get liquid composition
            ni=np.array(res.EquilibResult.Phases[LiquidPhaseName].xi[-1])*res.ScheilPhases[LiquidPhaseName][-1]
            components=res.EquilibResult.Phases[LiquidPhaseName].Endmembers[-1]
            Condition['T']=np.squeeze(Condition['T'])-dT
            
            for i, el in enumerate(components):
                Condition[el] = ni[i]
            Progress.update(1)
            scheil_cooling(LiquidPhaseName,Database,Condition,Unit=Unit,ListOfPhases=ListOfPhases,Output=res,dT=dT,Progress=Progress)
        
        if Progress.n>Progress.total:
            print('More iterations...')
            Progress.total = IterMax
            sys.setrecursionlimit(IterMax)
        elif  Progress.n==IterMax:
            print(f'Takes more than {IterMax} iterations, aborting the calculation...')
            Progress.close()
    return res