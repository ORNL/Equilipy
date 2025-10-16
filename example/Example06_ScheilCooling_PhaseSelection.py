#!/usr/bin/env python3
import numpy as np, matplotlib.pyplot as plt, polars as pl
import equilipy as eq
import os
if __name__ == "__main__":
    # Step 1: Parse database
    fpath=os.path.dirname(os.path.abspath(__file__))
    path ='/'.join(fpath.split('/')[:-1])
    datafile=f'{path}/database/AlCuMgSi_ORNL_FS73'
    DB=eq.read_dat(datafile+'.dat',FactSage8Plus=False)

    #Set input data

    NTP = dict({
        'T':900,
        'P': 1,
        'Al':0.89260,
        'Cu': 0.01745,
        'Mg':0.00114,
        'Si': 0.0881})
    #Phase selection
    PhasesAll=eq.list_phases(DB,list(NTP.keys())[2:])
    phases = [x for x in PhasesAll if 'FCC_A1' not in x]
    # phases=PhasesAll
    
    # Calculate Scheil cooling
    res=eq.scheil_cooling('LIQUID',DB,NTP,dT=5,ListOfPhases=phases,UnitIn=['C','atm','g'],UnitOut=['K','atm','mol'])
    print('Scheil Constituent information, mol. fr.:', res.ScheilConstituents)
    
    # Save data
    df=pl.DataFrame(res.to_dict())  
    df.write_csv(f'{fpath}/Result_Ex06_ACMS.csv')
    
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
