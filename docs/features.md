---
layout: default
title: Features and Examples
nav_enabled: true
nav_order: 4
---

# Features and examples
The following features are currently available.
- Single condition equilibrium calculations
- Batch equilibrium calculations
- Scheil-Gulliver solidification
- Phase selection
Examples of using each features are demonstraed here. The corresponding python script files are accessible in [Examples][examples].

## Simple equilibrium (Single condition)
Equilibrium calculation based on the CALPHAD approach requires a thermochemical database and NTP ensemble (N: elemental composition, T: temperature, P: pressure). In general, calculating equilibrium by the CALPHAD approach is a four-step process. An example of calculating a single NTP condition is given in [Example01][example01].

### Step 1: Parse thermochemical database
First, a thermochemical database must be provided. Currently, `equilipy` only supports ChemSage data format `.dat`, available from FactSage 7.3.

To read a `.dat` database in a python script:

```
import equilipy as eq

# Step 1: Parse database
datafile=f'./Database/AlCuMgSi_ORNL'
DB=eq.read_dat(datafile+'.dat')
```

### Step 2: Define NPT ensemble
Once providing a database, users should define a NTP condition for input variable. The input condition in `equilipy` is defined as a python dictionary:

```
NTP = dict({
    'T':700,
    'P': 1,
    'Al':0.060606061,
    'Cu':0.42424242,
    'Si':0.515151515})
```

### Step 3: Calculate equilibrium
For a single NPT condition, `equilib_sinlge()` function is used to calculate phase equilibria. Note that both a database `DB` and NTP condition `NTP` must be given as arguments:

```
res=eq.equilib_single(DB,NTP)
```

### Step 4: Post process
The calculated results are than stored in a result object `res`. Users can access to the information through class methods:

```
# 4.1: To print all stable phases
print(res.StablePhases['Name'])
print(res.StablePhases['Amount'])

# 4.2: To print all relevant phases
PhasesAll=list(res.Phases.keys())

for i,ph in enumerate(PhasesAll):
    print('--------------------------------------------------------------')
    print(f'Amount of {PhasesAll[i]}:',res.Phases[PhasesAll[i]].Amount)
    print(f'Endmembers of {PhasesAll[i]}:',res.Phases[PhasesAll[i]].Endmembers)
    print(f'Composition of {PhasesAll[i]}:',res.Phases[PhasesAll[i]].xi)
    print(' ')
```

## Phase selection
`equilipy` provides meta-stable phase equilibrium calculation through pre-selecting phases. The list of all available phase names can be obtained by `list_phases()`. Users should construct a custom list of phase names and parse it to `ListOfPhases` argument in `equilib_single` at Step 3.  An example of such calculation is given in [Example02][example02].


```
import equilipy as eq

# Step 1: Parse database
datafile=f'./Database/AlCuMgSi_ORNL'
DB=eq.read_dat(datafile+'.dat')

# Step 2: Parse input data
NTP = dict({
    'T':700,
    'P': 1,
    'Al':0.060606061,
    'Cu':0.42424242,
    'Si':0.515151515})

# Step 2: Select phases
# 2.1: Get all available phases of the given system
PhasesAll=eq.list_phases(DB,list(NTP.keys())[2:])
# 2.2: Select phases
phases = PhasesAll[:7]
print(f'Selected phases: {phases}')


# Step 3: Calculate equilibrium
res=eq.equilib_single(DB,NTP, ListOfPhases=phases)


# Step 4: Post process
# 4.1: print all stable phases
print(res.StablePhases['Name'])
print(res.StablePhases['Amount'])

# 4.2 print all phases
PhasesAll=list(res.Phases.keys())

for i,ph in enumerate(PhasesAll):
    print('--------------------------------------------------------------')
    print(f'Amount of {PhasesAll[i]}:',res.Phases[PhasesAll[i]].Amount)
    print(f'Endmembers of {PhasesAll[i]}:',res.Phases[PhasesAll[i]].Endmembers)
    print(f'Composition of {PhasesAll[i]}:',res.Phases[PhasesAll[i]].xi)
    print(' ')
```

{: .note }
Phase selection can also be used to calculation batch equilibrium and Scheil-Gulliver solidification.

## Batch equilibrium
Calculating multiple NTP conditions are also available using a batch process. By default, `equilipy` uses all available processors in the computing node via `multiprocessing`. 

{: .warning }
`multiprocessing` calls python script multiple times. Users should ensure using `if __name__ == "__main__":` in their main script.

An example of calculating batch equilibrium is given in [Example03][example03]. 

{: .note }
[Polars][polars] is used for reading large input data from an `Excel` file, which requires `fastexcel` as the optional dependancy. Install `fastexcel` via 
```pip install fastexcel```


```
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
```

## Scheil-Gulliver solidification
`equilipy` offers phase stability calculations during Scheil-Gulliver solidification. The calculation procedure is similar to that of phase equilibrium calculation except for Step 3. The example is given in [Example05][example05]

```
import numpy as np, matplotlib.pyplot as plt, polars as pl
import equilipy as eq

if __name__ == "__main__":
    # Step 1: Parse database
    datafile= './Database/AlCuMgSi_ORNL'
    DB=eq.read_dat(datafile+'.dat')

    # Step 2: Set input data
    system=['Al','Cu','Mg','Si']
    NTP = dict({
        'T':1000,
        'P': 1,
        'Al':0.75,
        'Cu': 0.05,
        'Mg':0.1,
        'Si': 0.1})

    # Step 3: Calculate Scheil cooling based on LIQUID as target phase
    TargetPhase='LIQUID'
    res=eq.scheil_cooling(TargetPhase,DB,NTP,dT=10)
    
    # Step 4: Post processing
    print('Scheil Constituent information, mol. fr.:', res.ScheilConstituents)
    df=pl.DataFrame(res.to_dict())  
    df.write_csv(f'Result_Ex05_ACMS.csv')

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
```



[examples]: https://github.com/ORNL/Equilipy/blob/main/example
[example01]: https://github.com/ORNL/Equilipy/blob/main/example/Example01_SingleEquilib.py
[example02]: https://github.com/ORNL/Equilipy/blob/main/example/Example02_SingleEquilib_PhaseSelection.py
[example03]: https://github.com/ORNL/Equilipy/blob/main/example/Example03_BatchEquilib.py
[example05]: https://github.com/ORNL/Equilipy/blob/main/example/Example05_ScheilCooling.py
[polars]: https://docs.pola.rs/