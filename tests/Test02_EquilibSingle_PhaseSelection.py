#!/usr/bin/env python3
import equilipy as eq
import os

#Parse database
#Parse database
fpath=os.path.dirname(os.path.abspath(__file__))
path ='/'.join(fpath.split('/')[:-1])
datafile=f'{path}/database/AlCuMgSi_ORNL_FS73'
DB=eq.read_dat(datafile+'.dat',FactSage8Plus=False)

# Parse input data
NTP = dict({
    'T':700,
    'P': 1,
    'Al':0.060606061,
    'Cu':0.42424242,
    'Si':0.515151515})


#Get all phases for the input system
PhasesAll=eq.list_phases(DB,list(NTP.keys())[2:])


# List of phases to be considered
phases = ['LIQUID', 'FCC_A1', 'HCP_A3',  'BCC_A2',  'BCT_A5']
res=eq.equilib_single(DB,NTP, ListOfPhases=phases)


#print all stable phases
print(res.StablePhases['Name'])
print(res.StablePhases['Amount_mole'])
print(res.StablePhases['Amount_mass'])

#print all phases
PhasesAll=list(res.Phases.keys())

for i,ph in enumerate(PhasesAll):
    print('--------------------------------------------------------------')
    print(f'Amount of {PhasesAll[i]}:',res.Phases[PhasesAll[i]].Amount_mole)
    print(f'Amount of {PhasesAll[i]}:',res.Phases[PhasesAll[i]].Amount_mass)
    print(f'Endmembers of {PhasesAll[i]}:',res.Phases[PhasesAll[i]].Endmembers_Xmole)
    print(f'Composition of {PhasesAll[i]}:',res.Phases[PhasesAll[i]].Endmembers_Xmass)
    print(f'Elements of {PhasesAll[i]}:',res.Phases[PhasesAll[i]].Elements_Xmole)
    print(f'Composition of {PhasesAll[i]}:',res.Phases[PhasesAll[i]].Elements_Xmass)
    print(' ')
    
#Check the results conmpared to FactSage
dG=abs(res.G+17557.30849)
print(dG)
if dG>1E-2: 
    raise ValueError
else:
    print(f'Gibbs energy differences: {dG}')