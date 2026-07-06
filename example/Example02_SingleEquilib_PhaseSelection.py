#!/usr/bin/env python3
import os
import equilipy as eq


# Step 1: Parse database
fpath = os.path.dirname(os.path.abspath(__file__))
path = os.path.dirname(fpath)
datafile = os.path.join(path, "database", "AlCuMgSi_ORNL_FS73")
DB = eq.read_dat(datafile + ".dat", factsage_8_plus=False)

# Step 2: Parse input data
NPT = dict(
    {"T": 700, "P": 1, "Al": 0.060606061, "Cu": 0.42424242, "Si": 0.515151515}
)

# Step 2: Select phases
# 2.1: Get all available phases of the given system
phases_all = eq.list_phases(DB, list(NPT.keys())[2:])
# 2.2: Select phases
phases = phases_all[:7]
print(f"Selected phases: {phases}")

# Step 3: Calculate equilibrium
res = eq.equilib_single(DB, NPT, phases=phases)

# Step 4: Post process
# 4.1: print all stable phases
print(res.stable_phases.names)
print([phase.amount_w for phase in res.stable_phases.values()])
print([phase.amount_n for phase in res.stable_phases.values()])

# 4.2: print all phases
phases_all = res.phases.names

for phase_name in phases_all    :
    phase = res.phase(phase_name)
    print("--------------------------------------------------------------")
    print(f"Amount of {phase_name}:", phase.amount_n)
    print(f"Elements of {phase_name}:", phase.elements.x_i)
    print(f"Endmembers of {phase_name}:", phase.endmembers.x_i)
    print(" ")
