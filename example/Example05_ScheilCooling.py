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

    # Step 2: Set input data
    system=['Al','Cu','Mg','Si']
    NTP = dict({
        'T':900,
        'P': 1,
        'Al':0.89260,
        'Cu': 0.01745,
        'Mg':0.00114,
        'Si': 0.0881})

    # Step 3: Calculate Scheil cooling based on LIQUID as target phase
    TargetPhase='LIQUID'
    res=eq.scheil_cooling(TargetPhase,DB,NTP,dT=1,Unit=['C','atm','g'])
    
    # Step 4: Post processing
    print('Scheil Constituent information, mol. fr.:', res.ScheilConstituents)
    df=pl.DataFrame(res.to_dict())  
    df.write_csv(f'{fpath}/Result_Ex05_ACMS.csv')

    # Plot Phase amount as function of temperature
    T= np.array(res.T)
    phases=list(res.ScheilPhases_mass.keys())
    fig, ax = plt.subplots(figsize=(5,4))

    for phase in phases:
        plt.plot(T,res.ScheilPhases_mass[phase],'-',linewidth=3,label=phase)
        
    ax.legend(fontsize=14)
    ax.set_xlabel('Temperature, K', fontsize=16)
    ax.set_ylabel('Phase amount, mol', fontsize=16)
    

    plt.tight_layout()
    plt.show()