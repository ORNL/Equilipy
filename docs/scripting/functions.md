# Functions

All calculations take a database object and an NPT condition: a plain
dictionary with temperature, pressure, and element amounts. Units default to
`["K", "atm", "moles"]` and can be changed with the `unit` argument.

## Single equilibrium

```python
import equilipy as eq

DB = eq.read_dat("AlCuMgSi_ORNL_FS83.dat")
NPT = {"T": 700, "P": 1, "Al": 0.06, "Cu": 0.42, "Si": 0.52}

res = eq.equilib_single(DB, NPT)
print(res.stable_phases.names)
```

Heat capacity `Cp` is computed by default with a perturbed solve; pass
`include_heat_capacity=False` to skip it and save time.

## Phase selection (metastable equilibria)

Restrict the calculation to a phase subset:

```python
all_phases = eq.list_phases(DB, ["Al", "Cu", "Si"])
res = eq.equilib_single(DB, NPT, phases=all_phases[:7])
```

`phases` works the same way in `equilib_batch` and the Scheil functions.

## Batch equilibrium

Each key of the condition dictionary holds one column of values; each row is
one condition. Batch runs execute in parallel with `n_cpu` worker processes.

```python
import polars as pl

NPT = pl.read_excel("Input_ACMS.xlsx", sheet_name="AlCuMgSi").to_dict()  # fastexcel engine, installed with equilipy
res = eq.equilib_batch(DB, NPT, unit=["K", "atm", "moles"])

pl.DataFrame(res.to_dict()).write_csv("results.csv")
```

Useful options: `n_cpu` (worker count; defaults to all cores),
`progress=False` (silence the progress bar in scripts and logs).

For multi-node clusters, see [HPC](hpc).

## Scheil-Gulliver solidification

```python
NPT = {"T": 1000, "P": 1, "Al": 0.75, "Cu": 0.05, "Mg": 0.1, "Si": 0.1}

res = eq.scheil_cooling("LIQUID", DB, NPT, delta_T=10)

print(res.scheil_constituents)               # constituents at the end of solidification
temperature = res.T
amounts = res.phase_amounts_n                # phase amounts along the path
```

When `phases` is omitted, `scheil_cooling()` excludes ordered SUBOM phases
from the default phase list for performance.  The result warning list and
console notice name any excluded phases.  To include them, pass
`include_ordered=True` or pass an explicit `phases=[...]` list.  NucleoScheil
is unaffected because its active phases are controlled by the selected
nucleation/undercooling phases.

Batch variants (`scheil_batch`, `scheil_constituent_batch`) run many
compositions in parallel. `nucleoscheil_cooling` and `nucleoscheil_batch`
add nucleation-undercooling control per phase.

## Saving and loading results

Save any result to a JSON-based `.eqres` bundle and reload it in a fresh
session — no recalculation, no pickle:

```python
eq.save_result(res, "run1.eqres")
loaded = eq.load_result("run1.eqres")
loaded.to_dict()                              # same post-processing as before
```
