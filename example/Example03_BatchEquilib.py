#!/usr/bin/env python3
import polars as pl, time
from datetime import timedelta
import equilipy as eq

if __name__ == "__main__":
    system = 'AlCuMgSi'
    # system = 'AlCuMg'
    # system = 'AlCuSi'
    # system = 'AlMgSi'
    # system = 'CuMgSi'


    # Step 1: Parse database
    datafile= './Database/AlCuMgSi_ORNL'
    DB=eq.read_dat(datafile+'.dat')


    # Step 2: Read input data using polars
    df_name= 'Input_ACMS.xlsx'
    NTP=pl.read_excel(f'{df_name}',sheet_name=system).to_dict()

    # Step 3: Calculate batch equilibrium
    starttime=time.time()
    res=eq.equilib_batch(DB,NTP)
    duration= time.time()-starttime
    dftime=pl.DataFrame({'Time, s':duration})

    # Step 4: Post processing
    print('Total processing time:',timedelta(seconds=duration))
    df=pl.DataFrame(res.to_dict())  
    df.write_csv(f'Result_Ex03_{system}.csv')
