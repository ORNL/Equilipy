#!/usr/bin/env python3
import numpy as np, matplotlib.pyplot as plt, polars as pl
import equilipy as eq
if __name__ == "__main__":
    #Parse database
    datafile= './Database/AlCuMgSi_ORNL'
    DB=eq.read_dat(datafile+'.dat')

    #Set input data
    units=['K','atm','moles']

    NTP = dict({
        'T':1000,
        'P': 1,
        'Al':0.8,
        'Cu': 0.05,
        'Mg':0.1,
        'Si': 0.05})
    #Phase selection
    PhasesAll=eq.list_phases(DB,list(NTP.keys())[2:])
    phases = ['LIQUID','FCC_A1','MG2SI(s)']
    # phases=PhasesAll[:2]+PhasesAll[10:]
    # Calculate Scheil cooling
    res=eq.scheil_cooling('LIQUID',DB,units,NTP,dT=10,ListOfPhases=phases)
    res.update_microconstituents()

    # Plot Phase amount as function of temperature
    T= np.array(res.T)
    phases=list(res.ScheilPhases.keys())
    fig, ax = plt.subplots(figsize=(5,4))

    for phase in phases:
        plt.plot(T,res.ScheilPhases[phase],'-',linewidth=3,label=phase)
        
    ax.legend(fontsize=14)
    ax.set_xlabel('Temperature, K', fontsize=16)

    plt.tight_layout()
    plt.show()
