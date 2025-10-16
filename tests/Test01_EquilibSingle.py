#!/usr/bin/env python3
import equilipy as eq
import os

#Parse database
fpath=os.path.dirname(os.path.abspath(__file__))
path ='/'.join(fpath.split('/')[:-1])
datafile=f'{path}/database/AlCuMgSi_ORNL_FS73'
DB=eq.read_dat(datafile+'.dat',FactSage8Plus=False)

# Parse input data
# Parse input data
NTP = dict({
    'T':700,
    'P': 1,
    'Al':0.060606061,
    'Cu':0.42424242,
    'Si':0.515151515})
elements=list(NTP.keys())[2:]
# print(eq.list_phases(DB,elements))
res=eq.equilib_single(DB,NTP,UnitIn=['K', 'atm', 'mol'])


#Check all stable phases
print(res.StablePhases['Name'])
print(res.StablePhases['Amount'])
    
#Check properties all phases
PhasesAll=list(res.Phases.keys())
for i,ph in enumerate(PhasesAll):
    print('--------------------------------------------------------------')
    print(f'Amount of {PhasesAll[i]}:',res.Phases[PhasesAll[i]].Amount)
    print(f'Endmembers of {PhasesAll[i]}:',res.Phases[PhasesAll[i]].Endmembers)
    print(f'Composition of {PhasesAll[i]}:',res.Phases[PhasesAll[i]].xi)
    print(' ')

#Check the results conmpared to FactSage
dG=abs(res.G+28103.14808)
if dG>1E-2: 
    raise ValueError
else:
    print(f'Gibbs energy differences: {dG}')
