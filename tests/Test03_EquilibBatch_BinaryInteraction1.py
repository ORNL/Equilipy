#!/usr/bin/env python3
import polars as pl, time, os
import numpy as np
from datetime import timedelta
import equilipy as eq


if __name__ == "__main__":
#Parse database
    fpath = os.path.dirname(os.path.abspath(__file__))
    path = os.path.dirname(fpath)
    datafile = os.path.join(path, 'database', 'AlCuMgSi_ORNL_FS83')
    DB=eq.read_dat(datafile+'.dat',FactSage8Plus=True)
    
    NTP ={
        'T':700*np.ones(10),
        'P': 1*np.ones(10),
        'Al':1-np.linspace(0.1,0.9,10),
        'Cu':np.linspace(0.1,0.9,10)
    }
    starttime=time.time()
    #Parse database
    DB=eq.read_dat(datafile+'.dat',FactSage8Plus=True)
    PhasesAll=eq.list_phases(DB,list(NTP.keys())[2:])
    phases = [phase for phase in PhasesAll if '#' not in phase]
    res=eq.equilib_batch(DB,NTP,ListOfPhases=PhasesAll,nCPU=10)

    duration= time.time()-starttime
    dftime=pl.DataFrame({'Time, s':duration})
    print('Total processing time:',timedelta(seconds=duration))

    df=pl.DataFrame(res.to_dict())
    df.write_csv(os.path.join(fpath,'Test03.csv'))
    print(df)