#!/usr/bin/env python3
import equilipy as eq

if __name__ == "__main__":
    #Parse database
    datafile=f'./Database/AlCuMgSi_ORNL'
    DB=eq.read_dat(datafile+'.dat')

    # Parse input data
    units=['K','atm','moles']
    NTP = dict({
        'T':700,
        'P': 1,
        'Al':0.060606061,
        'Cu':0.42424242,
        'Si':0.515151515})


    #Get all phases for the input system
    PhasesAll=eq.list_phases(DB,list(NTP.keys())[2:])


    # List of phases to be considered
    phases = PhasesAll[:7]
    print(f'Selected phases: {phases}')
    res=eq.equilib_single(DB,units,NTP, ListOfPhases=phases)


    #print all stable phases
    print(res.StablePhases['Name'])
    print(res.StablePhases['Amount'])

    #print all phases
    PhasesAll=list(res.Phases.keys())

    for i,ph in enumerate(PhasesAll):
        print('--------------------------------------------------------------')
        print(f'Amount of {PhasesAll[i]}:',res.Phases[PhasesAll[i]].Amount)
        print(f'Endmembers of {PhasesAll[i]}:',res.Phases[PhasesAll[i]].Endmembers)
        print(f'Composition of {PhasesAll[i]}:',res.Phases[PhasesAll[i]].xi)
        print(' ')