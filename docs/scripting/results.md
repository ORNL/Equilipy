# Result classes

Calculations return typed result objects:

- `EquilibResult` — from `equilib_single()` and `equilib_batch()`
- `ScheilResult` — from `scheil_cooling()` and related functions

A single-condition result holds one point; batch and path results hold many.
`result.points` is always a list; `result.point` is a shortcut valid only
for one-point results.

## System properties

```python
res.T, res.P                  # conditions
res.G, res.H, res.S, res.Cp   # J and J/K
res.n_i["Al"]                 # component amount, mole basis
res.w_i["Al"]                 # component amount, mass basis
```

## Phases

```python
res.stable_phases.names            # stable phase names
phase = res.phase("LIQUID")
phase.amount_n                     # phase amount, mole basis
phase.amount_n_basis               # basis for amount_n, e.g. formula_moles
phase.elements.x_i["Al"]           # element mole fraction in the phase
phase.endmembers.x_i               # endmember fractions
```

For stoichiometric condensed phases in pseudo-component systems, `amount_n`
is reported as formula moles.  This matches FactSage-style phase amount
semantics for oxide components such as `CaO`, `Al2O3`, and `SiO2`.
Solution phases usually keep their model-native amount basis, but
rank-reduced pseudo-component systems may report `pseudo_formula_moles` when
the solution endmembers project cleanly into the user component basis.  The
exact basis is exposed through `phase.amount_n_basis`.

Per-species activities use `"species@PHASE"` keys:

```python
res.activity["Al@LIQUID"]
```

## Scheil results

```python
res.scheil_constituents            # constituents at end of solidification
res.fl, res.fs                     # liquid/solid fractions along the path
res.phase_amounts_n                # {phase: amounts along the path}
res.phase_amount_n("FCC_A1")       # one phase's amounts
```

## Tables and export

```python
res.to_dict()                      # flattened table as a dictionary
pl.DataFrame(res.to_dict())        # -> polars dataframe
res.to_table()                     # ResultTable: column selection, CSV export
res.available_columns()            # list of exportable columns
```

Naming conventions used throughout: `n_*`/`amount_n` are mole basis,
`amount_n_basis` states the phase amount basis, `w_*`/`amount_w` mass basis;
`x_i` mole fractions, `w_i` mass fractions; `fl`/`fs` liquid/solid fractions.

For the complete class-level API, see the
[Result API reference](result_api_reference).
