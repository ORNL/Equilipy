# Python scripting

Everything the GUI does — and more — is available from Python. Every
script follows the same workflow:

```text
Load database  ->  NPT condition  ->  Calculation (Equilibrium / Solidification / ...)  ->  Result classes
```

```python
import equilipy as eq

# 1. Load database
DB = eq.read_dat("database/AlCuMgSi_ORNL_FS83.dat")

# 2. NPT condition: temperature, pressure, and element amounts
NPT = {"T": 900, "P": 1, "Al": 0.75, "Cu": 0.05, "Mg": 0.10, "Si": 0.10}

# 3. Calculation: equilibrium here; solidification and batch variants
#    take the same database and condition
res = eq.equilib_single(DB, NPT)

# 4. Result classes: phases, properties, tables
print(res.stable_phases.names)        # ['LIQUID', 'HCP_A3']
print(res.G, res.H, res.S)            # system properties, J and J/K
print(res.phase("LIQUID").amount_n)   # phase amount, mole basis
```

Each step has its own page:

- **[Database](database)** — read `.dat`/`.tdb` files, split element
  subsystems, write TDB.
- **[Functions](functions)** — equilibrium, batch, phase selection, and
  Scheil-Gulliver solidification.
- **[Result classes](results)** — access calculated properties, tables,
  and save/reload results.
- **[HPC](hpc)** — scale batch calculations across cluster nodes with MPI.

```{toctree}
:maxdepth: 2

database
functions
results
result_api_reference
hpc
```
