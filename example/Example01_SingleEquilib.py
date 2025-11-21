#!/usr/bin/env python3
import equilipy as eq, os

if __name__ == "__main__":
    # Step 1: Parse database
    fpath=os.path.dirname(os.path.abspath(__file__))
    path ='/'.join(fpath.split('/')[:-1])
    datafile=f'{path}/database/AlCuMgSi_ORNL_FS73'
    DB=eq.read_dat(datafile+'.dat',FactSage8Plus=False)
    
    # Step 2: Parse input data
    NTP = dict({
        'T':500,
        'P': 1,
        'Al':0.96,
        'Mg':0.03,
        'Si':0.01})

    # Step 3: Calculate equilibrium
    res=eq.equilib_single(DB,NTP,Unit=['K','atm','g'])

    # Step 4: Post process
    
    # 4.1: print all stable phases
    print(res.StablePhases['Name'])
    print(res.StablePhases['Amount_mass'])
    print(res.StablePhases['Amount_mole'])

    # 4.2: print all phases
    PhasesAll=list(res.Phases.keys())

    for i,ph in enumerate(PhasesAll):
        print('--------------------------------------------------------------')
        print(f'Amount of {PhasesAll[i]}:',res.Phases[PhasesAll[i]].Amount_mole)
        print(f'Endmembers of {PhasesAll[i]}:',res.Phases[PhasesAll[i]].Elements_Xmole)
        print(f'Composition of {PhasesAll[i]}:',res.Phases[PhasesAll[i]].Endmembers_Xmole)
        print(' ')
