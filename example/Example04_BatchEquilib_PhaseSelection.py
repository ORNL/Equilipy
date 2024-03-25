#!/usr/bin/env python3
import polars as pl, time
from datetime import timedelta
import equilipy as eq

if __name__ == "__main__":
    system = 'AlCuMgSi'

    #Parse database
    datafile= './Database/AlCuMgSi_ORNL'
    DB=eq.read_dat(datafile+'.dat')


    #Input data
    df_name= 'Input_ACMS.xlsx'
    units=['K','atm','moles']
    NTP=pl.read_excel(f'{df_name}',sheet_name=system).to_dict()

    starttime=time.time()
    PhasesAll=eq.list_phases(DB,list(NTP.keys())[2:])



    #Custom selected phases
    phases = ['LIQUID', 'FCC_A1', 'FCC_A1#2', 'HCP_A3', 'HCP_A3#2',\
            'BCC_A2', 'BCC_A2#2', 'BCT_A5', 'BCT_A5#2', 'DIAMOND_A4',\
            'ETA_CU19SI6','GAMMA_CU56SI11']

    #Print All phases
    print(f'Following phases are selected: {phases} from {PhasesAll}')
    res=eq.equilib_batch(DB,units,NTP,ListOfPhases=phases)

    duration= time.time()-starttime
    dftime=pl.DataFrame({'Time, s':duration})

    print('Total processing time:',timedelta(seconds=duration))

    df=pl.DataFrame(res.to_dict())
    # dftime.write_csv(f'Result/T9_{system}_t.csv')    
    df.write_csv(f'Result/Ex04_{system}.csv')
