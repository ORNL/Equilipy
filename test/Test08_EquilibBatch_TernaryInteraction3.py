#!/usr/bin/env python3
import polars as pl, time
import numpy as np
from datetime import timedelta
import equilipy as eq

if __name__ == "__main__":
    #Parse database
    datafile= './Database/AlCuMgSi_SK'
    DB=eq.read_dat(datafile+'.dat')


    #Input data
    grid=eq.simplex_grid(3,10)
    l,_=grid.shape
    NTP={
        'T': 700*np.ones(l),
        'P': np.ones(l),
        'Al':grid[:,0],
        'Mg':grid[:,1],
        'Si':grid[:,2]
    }
    starttime=time.time()

    res=eq.equilib_batch(DB,NTP)

    duration= time.time()-starttime
    dftime=pl.DataFrame({'Time, s':duration})

    print('Total processing time:',timedelta(seconds=duration))

    df=pl.DataFrame(res.to_dict())
    print(df.select(list(NTP.keys())+['G J', 'H J', 'S J/K', 'Cp J/K', 'StablePhaseNames']))
