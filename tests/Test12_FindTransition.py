import numpy as np
import equilipy as eq
import polars as pl
if __name__ == "__main__":
    # NucleoScheil input
    datafile = '../database/AlCuMgSi_ORNL_FS83'
    NTP = {
        'T': 2000,  # Initial temperature, will be updated by the function
        'P': 1,
        'Al': 0.955,
        'Cu': 0.005,
        'Mg':0.03,
        'Si':0.01
    }
    Elements= list(NTP.keys())[2:]

    # Parse database
    DB = eq.read_dat(datafile + '.dat')

    # Get phases
    # all_phases = eq.list_phases(DB, list(NTP.keys())[2:])
    # liquid_phase = ['LIQUID']
    # solid_phases= list(T_undercool.keys())
    # phases = liquid_phase + solid_phases
    # phases = ['LIQUID','FCC_A1','AL13FE4']


    # Step1: Initialize NucleoScheil and get the liquidus
    Ts=eq.find_transitions(DB, NTP, 1000,500,Unit=['C','atm','g'])

    NTP1 =NTP.copy()
    n = len(Ts)
    ones= np.ones(n)
    NTP1['T']=Ts
    NTP1['P']=ones
    NTP1['Al']=ones*(NTP['Al'])
    NTP1['Cu']=ones*(NTP['Cu'])
    NTP1['Mg']=ones*(NTP['Mg'])
    NTP1['Si']=ones*(NTP['Si'])
    df=pl.DataFrame(eq.equilib_batch(DB,NTP1).to_dict())
    print(df)