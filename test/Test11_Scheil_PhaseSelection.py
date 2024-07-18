#!/usr/bin/env python3
import numpy as np, matplotlib.pyplot as plt, polars as pl
import equilipy as eq

if __name__ == "__main__":
    #Parse database
    datafile= './Database/AlCuMgSi_Scheil_SK'
    DB=eq.read_dat(datafile+'.dat')

    #Set input data
    system=['Al','Cu','Mg','Si']
    NTP = dict({
        'T':1000,
        'P': 1,
        'Al':0.75,
        'Cu': 0.05,
        'Mg':0.1,
        'Si': 0.1})
    #Phase selection
    PhasesAll=eq.list_phases(DB,list(NTP.keys())[2:])
    # phases = ['LIQUID','FCC_A1','MG2SI(s)']
    phases=PhasesAll[:2]+PhasesAll[9:]

    TargetPhase='LIQUID'

    # Calculate Scheil cooling
    res=eq.scheil_cooling(TargetPhase,DB,NTP,dT=10,ListOfPhases=phases)
    
    # # Plot Phase amount as function of temperature
    Scheilresult={
        'T': np.array(res.T)
    }
    
    phases=list(res.ScheilPhases.keys())
    for phase in phases:
        Scheilresult[phase]=res.ScheilPhases[phase]
        

    df=pl.DataFrame(Scheilresult)
    print(df)
