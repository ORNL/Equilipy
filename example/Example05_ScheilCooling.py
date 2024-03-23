#!/usr/bin/env python3
import numpy as np, matplotlib.pyplot as plt, polars as pl
import equilipy as eq

#Parse database
datafile= './Database/AlCuMgSi_ORNL'
DB=eq.read_dat(datafile+'.dat')

#Set input data
units=['K','atm','moles']
system=['Al','Cu','Mg','Si']
# NTP = [1000,1,0.75,0.05,0.1,0.1]
NTP = dict({
    'T':1000,
    'P': 1,
    'Al':0.75,
    'Cu': 0.05,
    'Mg':0.1,
    'Si': 0.1})

TargetPhase='LIQUID'

# Calculate Scheil cooling
res=eq.scheil_cooling(TargetPhase,DB,units,NTP,dT=10)

# Plot Phase amount as function of temperature
T= np.array(res.T)
phases=list(res.ScheilPhases.keys())
fig, ax = plt.subplots(figsize=(5,4))

for phase in phases:
    plt.plot(T,res.ScheilPhases[phase],'-',linewidth=3,label=phase)
    
ax.legend(fontsize=14)
ax.set_xlabel('Temperature, K', fontsize=16)
ax.set_ylabel('Phase amount, mol', fontsize=16)

plt.tight_layout()
plt.show()