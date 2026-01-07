#!/usr/bin/env python3
import polars as pl, time, os
import numpy as np
from datetime import timedelta
import equilipy as eq
pl.Config.set_tbl_cols(50)
pl.Config.set_tbl_rows(100)

if __name__ == "__main__":

    #Parse database
    fpath = os.path.dirname(os.path.abspath(__file__))
    path = os.path.dirname(fpath)
    datafile = os.path.join(path, 'database', 'AlCuMgSi_ORNL_FS83')
    DB=eq.read_dat(datafile+'.dat')


    #Input data
    grid=eq.simplex_grid(3,10)
    # grid=grid[40:50]
    l,_=grid.shape
    NTP={
        'T': 700*np.ones(l),
        'P': np.ones(l),
        'Al':grid[:,0],
        'Cu':grid[:,1],
        'Mg':grid[:,2]
    }
    starttime=time.time()

    res=eq.equilib_batch(DB,NTP,nCPU=1)

    duration= time.time()-starttime
    dftime=pl.DataFrame({'Time, s':duration})

    print('Total processing time:',timedelta(seconds=duration))

    df=pl.DataFrame(res.to_dict())
    print(df)
