#!/usr/bin/env python3
import equilipy as eq

if __name__ == "__main__":
    # Step 1: Parse database
    datafile=f'./Database/AlCuMgSi_ORNL'
    DB=eq.read_dat(datafile+'.dat')

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
