#!/usr/bin/env python3
import polars as pl, time
import numpy as np
from datetime import timedelta
import equilipy as eq


if __name__ == "__main__":
    datafile=f'./Database/AlCuMgSi_SK'
    NTP ={
        'T':700*np.ones(10),
        'P': 1*np.ones(10),
        'Al':1-np.linspace(0.1,0.9,10),
        'Cu':np.linspace(0.1,0.9,10)
    }
    starttime=time.time()
    #Parse database
    DB=eq.read_dat(datafile+'.dat')
    PhasesAll=eq.list_phases(DB,list(NTP.keys())[2:])
    phases = [phase for phase in PhasesAll if '#' not in phase]
    res=eq.equilib_batch(DB,NTP,ListOfPhases=PhasesAll,nCPU=10)

    duration= time.time()-starttime
    dftime=pl.DataFrame({'Time, s':duration})
    print('Total processing time:',timedelta(seconds=duration))

    df=pl.DataFrame(res.to_dict())
    print(df.select(list(NTP.keys())+['G J', 'H J', 'S J/K', 'Cp J/K', 'StablePhaseNames']))