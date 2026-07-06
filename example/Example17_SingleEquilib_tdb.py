#!/usr/bin/env python3
import os
import sys

fpath = os.path.dirname(os.path.abspath(__file__))
path = os.path.dirname(fpath)
sys.path.insert(0, path)

import equilipy as eq

# Step 1: Parse database
datafile = os.path.join(path, "database", "AlFeSi_99Liu.tdb")
DB = eq.read_tdb(datafile, strict=True, auto_correct=False)

# Step 2: Parse input data
NPT = dict({"T": 600, "P": 1, "Al": 0.93, "Fe": 0.01, "Si": 0.06})

# Step 3: Calculate equilibrium
res = eq.equilib_single(DB, NPT, unit=["C", "atm", "g"])


# Step 4: Post process

# 4.1: print all stable phases
print(res.stable_phases.names)
print(res.G, res.H,res.S, res.G-res.H+res.T*res.S)
print([phase.amount_w for phase in res.stable_phases.values()])
# print([phase.amount_n for phase in res.stable_phases.values()])

# # 4.2: print all phases
# PhasesAll = res.phases.names

# for phase_name in PhasesAll:
#     phase = res.phase(phase_name)
#     print("--------------------------------------------------------------")
#     print(f"Amount of {phase_name}:", phase.amount_n)
#     print(f"Elements of {phase_name}:", phase.elements.x_i)
#     print(f"Endmembers of {phase_name}:", phase.endmembers.x_i)
#     print(" ")
