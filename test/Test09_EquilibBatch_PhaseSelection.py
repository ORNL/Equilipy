#!/usr/bin/env python3
import polars as pl, time
import numpy as np
from datetime import timedelta
import equilipy as eq

if __name__ == "__main__":
    system = 'AlCuMgSi'


    #Parse database
    datafile= './Database/AlCuMgSi_SK'
    DB=eq.read_dat(datafile+'.dat')


    #Input data
    grid=eq.simplex_grid(4,5)
    l,_=grid.shape
    NTP={
        'T': 700*np.ones(l),
        'P': np.ones(l),
        'Al': grid[:,0],
        'Cu': grid[:,1],
        'Mg': grid[:,2],
        'Si': grid[:,3],
    }

    starttime=time.time()
    PhasesAll=eq.list_phases(DB,list(NTP.keys())[2:])

    #Print All phases
    print(PhasesAll)

    #Custom selected phases
    phases = ['LIQUID', 'FCC_A1', 'FCC_A1#2', 'HCP_A3', 'HCP_A3#2',\
            'BCC_A2', 'BCC_A2#2', 'BCT_A5', 'BCT_A5#2', 'DIAMOND_A4',\
            'ETA_CU19SI6','GAMMA_CU56SI11']

    res=eq.equilib_batch(DB,NTP,ListOfPhases=phases)

    duration= time.time()-starttime
    dftime=pl.DataFrame({'Time, s':duration})

    print('Total processing time:',timedelta(seconds=duration))

    df=pl.DataFrame(res.to_dict())
    print(df.select(list(NTP.keys())+['G J', 'H J', 'S J/K', 'Cp J/K', 'StablePhaseNames']))