#!/usr/bin/env python3
import polars as pl, time
from datetime import timedelta
import equilipy as eq
import os

if __name__ == "__main__":
    system = 'AlCuMgSi'
    # system = 'AlCuMg'
    # system = 'AlCuSi'
    # system = 'AlMgSi'
    # system = 'CuMgSi'
    
    # Step 1: Parse database
    fpath=os.path.dirname(os.path.abspath(__file__))
    path ='/'.join(fpath.split('/')[:-1])
    datafile=f'{path}/database/AlCuMgSi_ORNL_FS73'
    DB=eq.read_dat(datafile+'.dat',FactSage8Plus=False)


    # Step 2: Read input data using polars
    df_name= 'Input_ACMS.xlsx'
    NTP=pl.read_excel(f'{fpath}/{df_name}',sheet_name=system).to_dict()

    # Step 3: Calculate batch equilibrium
    starttime=time.time()
    res=eq.equilib_batch(DB,NTP,UnitIn=['K','atm','moles'],UnitOut=['C','atm','g'])
    duration= time.time()-starttime
    dftime=pl.DataFrame({'Time, s':duration})

    # Step 4: Post processing
    print('Total processing time:',timedelta(seconds=duration))
    df=pl.DataFrame(res.to_dict())  
    df.write_csv(f'{fpath}/Result_Ex03_{system}.csv')
