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

# Set input data

NPT = dict(
    {"T": 900, "P": 1, "Al": 0.89260, "Cu": 0.01745, "Mg": 0.00114, "Si": 0.0881}
)
# Phase selection
PhasesAll = eq.list_phases(DB, list(NPT.keys())[2:])
phases = [x for x in PhasesAll if "FCC_A1" not in x]
# phases=PhasesAll

# Calculate Scheil cooling
res = eq.scheil_cooling(
    "LIQUID",
    DB,
    NPT,
    delta_T=5,
    phases=phases,
    unit=["C", "atm", "g"],
)
print("Scheil Constituent information, mol. fr.:", res.scheil_constituents)

# Save data
df = pl.DataFrame(res.to_dict())
df.write_csv(f"{fpath}/Result_Ex06_ACMS.csv")

# Plot Phase amount as function of temperature
T = np.array(res.T)
phases = list(res.phase_amounts_n.keys())
fig, ax = plt.subplots(figsize=(5, 4))

for phase in phases:
    plt.plot(T, res.phase_amount_n(phase), "-", linewidth=3, label=phase)

ax.legend(fontsize=14)
ax.set_xlabel("Temperature, K", fontsize=16)

plt.tight_layout()
plt.show()
