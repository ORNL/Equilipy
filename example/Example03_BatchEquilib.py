#!/usr/bin/env python3
import os
import time
from datetime import timedelta
import polars as pl
import equilipy as eq



system = 'AlCuMg'
# system = 'AlCuSi'
# system = 'AlMgSi'
# system = 'CuMgSi'
system = "AlCuMgSi"

# Step 1: Parse database
fpath = os.path.dirname(os.path.abspath(__file__))
path = os.path.dirname(fpath)
# datafile = os.path.join(path, "database", "AlCuMgSi_SK.tdb")
# DB = eq.read_tdb(datafile, strict=False, auto_correct=False)

datafile = os.path.join(path, "database", "AlCuMgSi_ORNL_FS83.dat")
DB = eq.read_dat(datafile)


# Step 2: Read input data using polars
df_name = "Input_ACMS.xlsx"
NPT = pl.read_excel(os.path.join(fpath, df_name), sheet_name=system).to_dict()
# Step 3: Calculate batch equilibrium
starttime = time.time()
res = eq.equilib_batch(DB, NPT, unit=["K", "atm", "moles"])
duration = time.time() - starttime
dftime = pl.DataFrame({"Time, s": duration})

# Step 4: Post processing
print("Total processing time:", timedelta(seconds=duration))
df = pl.DataFrame(res.to_dict())
print(df)
df.write_csv(f"{fpath}/Result_Ex03_{system}.csv")
