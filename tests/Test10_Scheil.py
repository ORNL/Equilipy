#!/usr/bin/env python3
import numpy as np, matplotlib.pyplot as plt, polars as pl, os
import equilipy as eq

if __name__ == "__main__":
    #Parse database
    fpath=os.path.dirname(os.path.abspath(__file__))
    path ='/'.join(fpath.split('/')[:-1])
    datafile=f'{path}/database/AlCuMgSi_ORNL_FS83'
    DB=eq.read_dat(datafile+'.dat',FactSage8Plus=True)

    #Set input data
    system=['Al','Cu','Mg','Si']
    NTP = dict({
        'T':1000-273.15,
        'P': 1,
        'Al':0.8,
        'Cu': 0.05,
        'Mg':0.1,
        'Si': 0.05})

    TargetPhase='LIQUID'

    # Calculate Scheil cooling
    res=eq.scheil_cooling(TargetPhase,DB,NTP,UnitIn=['C', 'atm', 'g'],dT=5)
        

    df=pl.DataFrame(res.to_dict())
    sc = pl.DataFrame(res.ScheilConstituents)
    print(df)
    print(sc)