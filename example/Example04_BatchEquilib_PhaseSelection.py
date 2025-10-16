#!/usr/bin/env python3
import polars as pl, time
from datetime import timedelta
import equilipy as eq
import os

if __name__ == "__main__":
    system = 'AlCuMgSi'
    
    # Step 1: Parse database
    fpath=os.path.dirname(os.path.abspath(__file__))
    path ='/'.join(fpath.split('/')[:-1])
    datafile=f'{path}/database/AlCuMgSi_ORNL_FS73'
    DB=eq.read_dat(datafile+'.dat',FactSage8Plus=False)

    # Step 2: Input data
    df_name= 'Input_ACMS.xlsx'
    NTP=pl.read_excel(f'{fpath}/{df_name}',sheet_name=system).to_dict()

    # Phase selection
    PhasesAll=eq.list_phases(DB,list(NTP.keys())[2:])
    
    # Custom selected phases
    phases = ['LIQUID', 'FCC_A1', 'FCC_A1#2', 'HCP_A3', 'HCP_A3#2',\
            'BCC_A2', 'BCC_A2#2', 'BCT_A5', 'BCT_A5#2', 'DIAMOND_A4',\
            'ETA_CU19SI6','GAMMA_CU56SI11']
    print(f'Following phases are selected: {phases} from {PhasesAll}')
    
    # Step 3: Calculate equilibrium
    starttime=time.time()
    res=eq.equilib_batch(DB,NTP,ListOfPhases=phases)
    duration= time.time()-starttime
    dftime=pl.DataFrame({'Time, s':duration})
    print('Total processing time:',timedelta(seconds=duration))

    # Step 4: Post processing
    df=pl.DataFrame(res.to_dict()) 
    df.write_csv(f'{fpath}/Result_Ex04_{system}.csv')
