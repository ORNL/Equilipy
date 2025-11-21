#!/usr/bin/env python3
import polars as pl, time, os
import numpy as np
from datetime import timedelta
import equilipy as eq
pl.Config.set_tbl_cols(20)
pl.Config.set_tbl_rows(100)
if __name__ == "__main__":


    #Parse database
    fpath=os.path.dirname(os.path.abspath(__file__))
    path ='/'.join(fpath.split('/')[:-1])
    datafile=f'{path}/database/AlCuMgSi_ORNL_FS83'
    DB=eq.read_dat(datafile+'.dat',FactSage8Plus=True)


    #Input data
    grid=eq.simplex_grid(4,5)
    # grid=grid[30:40,:]
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
    print(df.select(['G [J]', 'H [J]', 'S [J/K]', 'Cp [J/K]', 'StablePhaseNames']))