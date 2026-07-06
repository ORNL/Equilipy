#!/usr/bin/env python3
import os
import time
from datetime import timedelta

import polars as pl

import equilipy as eq

system = "AlCuMgSi"

# Step 1: Parse database
fpath = os.path.dirname(os.path.abspath(__file__))
path = os.path.dirname(fpath)
datafile = os.path.join(path, "database", "AlCuMgSi_ORNL_FS73")
DB = eq.read_dat(datafile + ".dat", factsage_8_plus=False)

# Step 2: Input data
df_name = "Input_ACMS.xlsx"
NPT = pl.read_excel(os.path.join(fpath, df_name), sheet_name=system).to_dict()

# Phase selection
PhasesAll = eq.list_phases(DB, list(NPT.keys())[2:])

# Custom selected phases
phases = [
    "LIQUID",
    "FCC_A1",
    "FCC_A1#2",
    "HCP_A3",
    "HCP_A3#2",
    "BCC_A2",
    "BCC_A2#2",
    "BCT_A5",
    "BCT_A5#2",
    "DIAMOND_A4",
    "ETA_CU19SI6",
    "GAMMA_CU56SI11",
]
print(f"Following phases are selected: {phases} from {PhasesAll}")

# Step 3: Calculate equilibrium
starttime = time.time()
res = eq.equilib_batch(DB, NPT, phases=phases)
duration = time.time() - starttime
dftime = pl.DataFrame({"Time, s": duration})
print("Total processing time:", timedelta(seconds=duration))

# Step 4: Post processing
df = pl.DataFrame(res.to_dict())
print(df)
# df.write_csv(f"{fpath}/Result_Ex04_{system}.csv")
