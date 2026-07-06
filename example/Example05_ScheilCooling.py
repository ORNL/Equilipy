#!/usr/bin/env python3
import os
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import equilipy as eq


# Step 1: Parse database
fpath = os.path.dirname(os.path.abspath(__file__))
path = os.path.dirname(fpath)
datafile = os.path.join(path, "database", "AlCuMgSi_ORNL_FS73")
DB = eq.read_dat(datafile + ".dat", factsage_8_plus=False)

# Step 2: Set input data
system = ["Al", "Cu", "Mg", "Si"]
NPT = dict(
    {"T": 900, "P": 1, "Al": 0.89260, "Cu": 0.01745, "Mg": 0.00114, "Si": 0.0881}
)

# Step 3: Calculate Scheil cooling based on LIQUID as target phase
TargetPhase = "LIQUID"
res = eq.scheil_cooling(TargetPhase, DB, NPT, delta_T=1, unit=["C", "atm", "g"])

# Step 4: Post processing
print("Scheil Constituent information, mol. fr.:", res.scheil_constituents)
result_dict = res.to_dict()
df = pl.DataFrame(result_dict)
preview_columns = [
    column
    for column in df.columns
    if column.startswith(("T [", "fs_w", "label", "LIQUID_endmembers_w_"))
]
print(df.select(preview_columns))
df.write_csv(f"{fpath}/Result_Ex05_ACMS.csv")

# Plot Phase amount as function of temperature
T = np.array(res.T)
phases = list(res.phase_amounts_w.keys())
fig, ax = plt.subplots(figsize=(5, 4))

for phase in phases:
    plt.plot(T, res.phase_amount_w(phase), "-", linewidth=3, label=phase)

ax.legend(fontsize=14)
ax.set_xlabel("Temperature, K", fontsize=16)
ax.set_ylabel("Phase amount, g", fontsize=16)

plt.tight_layout()
plt.show()
