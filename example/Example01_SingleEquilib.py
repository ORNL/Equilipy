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
        'T':700,
        'P': 1,
        'Al':0.060606061,
        'Cu':0.42424242,
        'Si':0.515151515})

    # Step 3: Calculate equilibrium
    res=eq.equilib_single(DB,NTP)

    # Step 4: Post process
    
    # 4.1: print all stable phases
    print(res.StablePhases['Name'])
    print(res.StablePhases['Amount'])

    # 4.2: print all phases
    PhasesAll=list(res.Phases.keys())

    for i,ph in enumerate(PhasesAll):
        print('--------------------------------------------------------------')
        print(f'Amount of {PhasesAll[i]}:',res.Phases[PhasesAll[i]].Amount)
        print(f'Endmembers of {PhasesAll[i]}:',res.Phases[PhasesAll[i]].Endmembers)
        print(f'Composition of {PhasesAll[i]}:',res.Phases[PhasesAll[i]].xi)
        print(' ')
