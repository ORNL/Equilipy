#!/usr/bin/env python3
import os

import matplotlib.pyplot as plt
import numpy as np

import equilipy as eq

# Step 1: Parse database
fpath = os.path.dirname(os.path.abspath(__file__))
path = os.path.dirname(fpath)
datafile = os.path.join(path, "database", "AlFeSi_99Liu_v84")
DB = eq.read_dat(datafile + ".dat")


# Step 2: Input data
x = np.linspace(0.0001, 1 - 1e-4, 101)
NPT = {
    "T": 900 * np.ones_like(x),
    "P": np.ones_like(x),
    "Al": 1 - x,
    "Fe": x,
}

# Step 3: Parse available phases and select LIQUID
phases_all = eq.list_phases(DB, list(NPT.keys())[2:])
phases = ["LIQUID"]
missing_phases = [phase for phase in phases if phase not in phases_all]
if missing_phases:
    raise ValueError(f"Selected phase(s) not available: {missing_phases}")

print(f"Parsed {len(phases_all)} phase(s).")
print(f"Selected phases: {phases}")

# Step 4: Calculate equilibrium
res = eq.equilib_batch(DB, NPT, phases=phases, n_cpu=1)
print(f"Result phases: {res.phases.names}")

# Step 5: Plot system thermodynamic properties
fig, axes = plt.subplots(3, 1, figsize=(7, 8), sharex=True)
system_properties = [
    ("G", res.G, "G [J]"),
    ("H", res.H, "H [J]"),
    ("S", res.S, "S [J/K]"),
]
for ax, (title, values, ylabel) in zip(axes, system_properties, strict=False):
    ax.plot(x, values, "-o", ms=3)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
axes[-1].set_xlabel("x(Fe)")
fig.suptitle("Al-Fe LIQUID System Properties at 900 K")
fig.tight_layout()
ghs_path = os.path.join(fpath, "Result_Ex16_GHS.png")
fig.savefig(ghs_path, dpi=200)
print(f"Saved {ghs_path}")

# Step 6: Plot phase partial thermodynamic properties
phase_name = "LIQUID"
species_names = list(res.phases[phase_name].partial_enthalpy)
partial_properties = [
    ("Partial Gibbs energy", res.partial_gibbs, "g [J]"),
    ("Partial enthalpy", res.partial_enthalpy, "h [J]"),
    ("Partial entropy", res.partial_entropy, "s [J/K]"),
    ("Partial heat capacity", res.partial_heat_capacity, "cp [J/K]"),
]

fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
for ax, (title, property_view, ylabel) in zip(
    axes.ravel(),
    partial_properties,
    strict=False,
):
    for species_name in species_names:
        ax.plot(
            x,
            property_view[f"{species_name}@{phase_name}"],
            "-o",
            ms=3,
            label=species_name,
        )
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
    ax.legend()

for ax in axes[-1, :]:
    ax.set_xlabel("x(Fe)")

fig.suptitle("Al-Fe LIQUID Partial Properties at 900 K")
fig.tight_layout()
partial_path = os.path.join(fpath, "Result_Ex16_LiquidPartialProperties.png")
fig.savefig(partial_path, dpi=200)
print(f"Saved {partial_path}")

# Step 7: Plot activities, computed from partial Gibbs energy minus
# reference Gibbs energy.
fig, ax = plt.subplots(figsize=(7, 5))
for species_name in species_names:
    key = f"{species_name}@{phase_name}"
    ax.plot(x, res.activity[key], "-o", ms=3, label=species_name)
ax.set_xlabel("x(Fe)")
ax.set_ylabel("activity [-]")
ax.set_title("Al-Fe LIQUID Activity at 900 K")
ax.grid(True, alpha=0.3)
ax.legend()
fig.tight_layout()
activity_path = os.path.join(fpath, "Result_Ex16_LiquidActivity.png")
fig.savefig(activity_path, dpi=200)
print(f"Saved {activity_path}")

# Step 8: Check that the total properties agree with the partial properties.
property_checks = [
    ("G", res.G, res.partial_gibbs),
    ("H", res.H, res.partial_enthalpy),
    ("S", res.S, res.partial_entropy),
    ("Cp", res.Cp, res.partial_heat_capacity),
]
for name, total_values, partial_view in property_checks:
    weighted_sum = (
        (1 - x) * np.asarray(partial_view[f"Al@{phase_name}"])
        + x * np.asarray(partial_view[f"Fe@{phase_name}"])
    )
    max_error = np.nanmax(np.abs(np.asarray(total_values) - weighted_sum))
    print(f"{name} partial-property balance max abs error: {max_error:.6g}")

for species_name in species_names:
    key = f"{species_name}@{phase_name}"
    print(
        f"{key} activity range:",
        f"{np.nanmin(res.activity[key]):.6g}",
        "to",
        f"{np.nanmax(res.activity[key]):.6g}",
    )

plt.show()
