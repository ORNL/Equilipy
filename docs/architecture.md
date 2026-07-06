# Program architecture

This page documents the current Equilipy workflow as implemented in the Python
and Fortran code. It is a living current-state document for Phase 1 debugging.

The purpose is to make the data path explicit:

```text
user input
-> Python database/condition preparation
-> f2py transfer into Fortran modules
-> Fortran thermodynamic/minimization workflow
-> Fortran post-processing state
-> Python result objects
```

This page is intentionally descriptive. It records how the code works now, even
where the current workflow has cleanup or state-isolation risks.

## Main public entry points

| Public API | Primary file | Main result |
|---|---|---|
| `read_dat` | `equilipy/read_dat.py` | parsed database dictionary |
| `list_phases` | `equilipy/list_phases.py` | list of phase names for selected elements |
| `equilib_single` | `equilipy/equilib_single.py` | `EquilibResult` |
| `equilib_batch` | `equilipy/equilib_batch.py` | `EquilibResult` |
| `scheil_cooling` | `equilipy/scheil_cooling.py` | `ScheilResult` |
| `scheil_batch` | `equilipy/scheil_cooling.py` | flattened dictionary |
| `scheil_constituent_batch` | `equilipy/scheil_cooling.py` | flattened constituent dictionary |
| `equilib_cooling` | `equilipy/equilib_cooling.py` | Scheil-compatible cooling path with embedded equilibrium history |
| `nucleoscheil_cooling` | `equilipy/nucleoscheil_cooling.py` | Scheil-like result |
| `find_transitions` | `equilipy/find_transition.py` | transition temperature array |
| `find_first_transition` | `equilipy/find_transition.py` | highest transition temperature |
| `property_solution` | `equilipy/property_solution.py` | phase-local `EquilibResult` |
| `property_compound` | `equilipy/property_compound.py` | one-species `EquilibResult` |

## State ownership overview

Equilipy currently uses three major state layers.

| Layer | Owner | Examples | Main purpose |
|---|---|---|---|
| Parsed Python database state | `equilipy.variables` | `nElementsCS`, `cElementNameCS`, `dStoichSpeciesCS`, `cPhaseNames` | Stores ChemSage parser output and selected-system maps |
| f2py/Fortran input and solver state | `equilipy.equilifort` modules | `moduleparsecs`, `modulethermoio`, `modulethermo`, `modulegemsolver`, `modulesubmin` | Holds database arrays, input condition, minimizer state, and output state |
| Python result state | `equilipy.results` and `equilipy.post_process` | `EquilibPoint`, `EquilibResult`, `PhaseResult`, `ScheilResult` | Copies selected Fortran state into user-facing Python objects |

The important practical point is that the Fortran modules are global mutable
state inside each Python process. A process can run many calculations, so reset
coverage matters.

## Calculation workflow map

This is the working Phase 1 flowchart for the equilibrium calculation path. Each
box is a documentation target: we will identify the Python globals, Fortran
globals, allocation points, mutation points, and reset/deallocation points for
that step.

```text
1. Parse database
   read and parse ChemSage/FactSage database
        |
        v
2. Receive user input
   user provides NPT condition as a Python dictionary
        |
        v
3. Check selected system
   confirm input elements are a subset of the database system
   define phases related to the selected element system
        |
        v
4. Calculate Gibbs energies and run leveling
   compute compound/endmember Gibbs energies
   estimate an initial assemblage with leveling
        |
        v
5. Calculate solution-phase driving forces
   use elemental potentials from leveling
   run subminimization for each candidate solution phase
        |
        v
6. Iterate assemblage search
   repeat leveling and driving-force calculations
   update the candidate phase assemblage in the GEM solver loop
        |
        v
7. Run Lagrangian multiphase minimization
   solve phase amounts, species fractions, and element potentials
        |
        v
8. Check convergence
   if not converged, change the phase assemblage and continue iterating
        |
        v
9. Post-process converged result
   calculate output phase/species amounts, thermodynamic properties,
   and Python result objects
```

The same workflow as a tracking table:

| Step | Current code path | Global-state mapping target |
|---|---|---|
| 1. Parse database | `read_dat()`, parser helpers | Python parser globals in `equilipy.variables`; later copied to `ModuleParseCS`. See [](developer/database_parser_contract.md). |
| 2. Receive input | public API condition dictionary, `_dict2np()` | `var.dConditionSys`; `ModuleThermoIO` input fields |
| 3. Check selected system | `list_phases()`, `_pyvar2fvar()`, `system_check()`, optional `phase_selection()` | selected-system maps in `variables.py`; active parser arrays in `ModuleParseCS`; selected phase arrays in `ModuleThermo` |
| 4. Gibbs energies and leveling | `compthermodata()`, `gemsolver()`, `InitGEMSolver`, `LevelingSolver` | standard Gibbs arrays, element potentials, initial assemblage, leveling history |
| 5. Solution driving forces | subminimization routines called from the solver | `ModuleSubMin`; solution Gibbs/driving-force arrays; candidate phase state |
| 6. GEM solver iteration | `GEMSolver` / GEM solver loop | `ModuleGEMSolver`; assemblage updates; convergence flags; iteration counters |
| 7. Lagrangian minimization | `MultiPhaseMinimizer` and related Lagrangian routines | phase amounts, species mole fractions, update vectors, element potentials |
| 8. Convergence check | solver convergence logic and Python `dgemfunctionnorm` check | `dGEMFunctionNorm`, convergence flags, phase-change history |
| 9. Post-process | `postprocess()`, `EquilibResult.append_output()`, `EquilibPoint.from_fortran()` | output arrays in `ModuleThermoIO`; phase/species arrays in `ModuleThermo`; Python result objects |

## Database parsing workflow

`read_dat(file_name, factsage_8_plus=True)` parses a ChemSage `.dat` file and
returns a plain Python dictionary.

```text
read_dat
-> set supported solution phase types in variables.py
-> read file text into var.DataBase
-> ParseHSCpFunctions, when FactSage 8+ mode is enabled
-> ParseCSHeader
-> ParseCSDataBlock
-> assemble var.cPhaseNames
-> return var.to_dict()
```

Key files:

- `equilipy/read_dat.py`
- `equilipy/variables.py`
- `equilipy/parse_HSCp_functions.py`
- `equilipy/parse_chemsage_header.py`
- `equilipy/parse_chemsage_data_block.py`

Important state:

- Parser output is first written into `equilipy.variables`.
- `variables.to_dict()` exports a subset of that module state as the database
  dictionary.
- The returned dictionary is the object users pass into calculation functions.

Current-state note:

- `read_dat` does not directly populate Fortran solver modules for a
  calculation. Calculation entrypoints call `load_fortran_database()`, which
  copies the parsed database to Fortran when the cached payload is stale.

## Loading a parsed database for calculation

`load_database(database)` copies a parsed database dictionary back into
`equilipy.variables`.

```text
load_database(database)
-> assign database dictionary fields into equilipy.variables
-> return None
```

`load_fortran_database(database)` owns the Python-to-Fortran database transfer:

```text
load_fortran_database(database)
-> load_database(database)
-> _pyvar2fvar(var) if the cached Fortran payload is stale
-> reset cheap per-calculation Fortran masks
-> return None
```

Then `list_phases(database, elements)` prepares the selected system in Python
and updates the Fortran selected-system flag:

```text
list_phases
-> load_database(database)
-> system_check(elements)
-> build var.PhaseNameSys mapping
-> return var.cPhaseNameSys
```

`_pyvar2fvar(var)` is the Python-to-Fortran database transfer helper in this
runtime path. It lives in `equilipy/utils.py` and is called by
`load_fortran_database()` when the cached Fortran database payload is stale.

The helper resets parser/calculation status:

- `fort.moduleparsecs.info = 0`
- `fort.modulethermoio.infothermo = 0`

Then it writes parser arrays from `equilipy.variables` into
`fort.moduleparsecs` and uses the f2py string helpers `fort.str1d`,
`fort.str2d`, and `fort.str3d` for Fortran character arrays.

`system_check(elements)` builds the current selected-system maps:

- `var.cComponentNameSys`
- `var.iElementSys`
- `var.iElementSysIndex`
- `var.iElementDBIndex`
- `var.iSys2DBSpecies`
- `var.iSys2DBSoln`
- `var.iSys2DBComp`
- `var.cSpeciesNameSys`
- `var.cPhaseNameSys`
- default copies used by post-processing

It also sets:

- `fort.moduleparsecs.lendmembers2species = True`

Current-state note:

- `list_phases()` takes a `database` argument, but several call sites pass
  `var` after `load_database()` has already populated `equilipy.variables`.
  This works because `load_database()` catches incompatible inputs and because
  `var` already holds the needed state, but it is a stateful pattern worth
  cleaning up or documenting more tightly.

## Input condition transfer

Input conditions enter as dictionaries such as:

```python
condition = {
    "T": 700,
    "P": 1,
    "Al": 0.75,
    "Cu": 0.05,
    "Mg": 0.10,
    "Si": 0.10,
}
```

`_dict2np()` converts the dictionary to:

- a header list
- a 2D NumPy array of rows

`input_condition(unit, composition, include_heat_capacity=True)` writes the
active row into `fort.modulethermoio`:

| Fortran field | Source |
|---|---|
| `dTemperatureDiff` | `0.01` if heat capacity is enabled, otherwise `0.00` |
| `cInputUnitTemperature` | `unit[0]` |
| `cInputUnitPressure` | `unit[1]` |
| `cInputUnitMass` | `unit[2]` |
| `dTemperature` | first condition value |
| `dPressure` | second condition value |
| `dElementMass` | element amounts placed by atomic number |
| `iPrintResultsMode` | fixed to `2` |

Before writing `dElementMass`, the composition vector is reordered with
`var.iElementSysIndex` so that Python input order matches the database/system
order expected by Fortran. The reordered values are then placed into
`dElementMass` by atomic number using `var.iElementSys`.

## Optional phase selection

`phase_selection(phases)` narrows the active phase set after the selected system
has been prepared.

It updates Python state:

- `var.iSys2DBSoln`
- `var.iSys2DBComp`
- `var.iSys2DBSpecies`
- `var.cPhaseNameSys`
- `var.PhaseNameSys`
- `var.iPhaseCS`

It also updates Fortran state:

- `fort.moduleparsecs.iphasecs`
- `fort.modulethermo.isolnps`

Current-state note:

- Phase selection depends on `var.PhaseNameSys`, so it must run after
  `list_phases()`.
- Immiscible phases are represented by the exact selected-system phase names,
  including `#N` suffixes such as `FCC_A1#2`.

## Single equilibrium workflow

The public `equilib_single()` path is:

```text
equilib_single
-> EquilibResult()
-> _preprocess_single
   -> load_database
   -> _dict2np
   -> list_phases
      -> load_database
      -> _pyvar2fvar
      -> system_check
   -> optional phase_selection
   -> remove zero-composition elements, if needed
   -> repeat list_phases/phase_selection for the reduced element set
   -> var.dConditionSys = composition_condition
   -> input_condition
-> minimize
-> EquilibResult.append_output or EquilibResult.append_error
-> fort.resetthermoall()
-> return EquilibResult
```

`_preprocess_single()` removes zero-composition elements from the selected
system before minimization. If any input element amount is zero, the selected
elements and composition are reduced, and `list_phases()` is run again for the
reduced system.

The internal `_equilib_single()` uses the same preprocessing and minimization
path, but it does not create an `EquilibResult` and does not call `fort.resetthermoall()`.
It is used by transition search and Scheil calculations when the caller needs to
inspect Fortran state directly after minimization.

Current-state reset behavior:

- `equilib_single()` calls `fort.resetthermoall()` after copying output into the
  Python result.
- `_equilib_single()` does not reset; the caller owns the Fortran state after the
  call.

## Fortran minimization workflow

`minimize()` is the Python wrapper around the Fortran calculation sequence:

```text
minimize
-> fort.checkthermoinput()
-> fort.initthermo()
-> fort.checksystem()
-> fort.compthermodata()
-> fort.checkthermodata()
-> fort.gemsolver()
-> fort.postprocess()
-> check fort.modulegemsolver.dgemfunctionnorm
```

Key Fortran files:

| Stage | Python call | Fortran file |
|---|---|---|
| input validation | `checkthermoinput` | `fsrc/3_SystemDef/CheckThermoInput.f90` |
| allocation/initialization | `initthermo` | `fsrc/3_SystemDef/InitThermo.f90` |
| system consistency | `checksystem` | `fsrc/3_SystemDef/CheckSystem.f90` |
| thermodynamic data | `compthermodata` | `fsrc/3_SystemDef/CompThermoData.f90` |
| data validation | `checkthermodata` | `fsrc/3_SystemDef/CheckThermoData.f90` |
| minimization | `gemsolver` | `fsrc/7_MultiPhaseMinimizer/GEMSolver.f90` |
| output preparation | `postprocess` | `fsrc/3_SystemDef/PostProcess.f90` |

The current `GEMSolver` path calls:

```text
GEMSolver
-> InitGEMSolver
-> LevelingSolver
-> MultiPhaseMinimizer
```

Related Fortran solver modules:

- `ModuleParseCS`: database/parser arrays transferred from Python.
- `ModuleThermoIO`: user input and output arrays visible through f2py.
- `ModuleThermo`: internal thermodynamic and phase-assemblage state.
- `ModuleGEMSolver`: global solver iteration, convergence, and active-set state.
- `ModuleSubMin`: solution-phase subminimization state.

Current-state convergence and error behavior:

- Python raises `EquilibError` if a wrapped Fortran call raises an exception.
- `fort.modulethermoio.infothermo` is captured for error-message context, but
  it is not currently checked as a failure condition.
- Python also raises `EquilibError` after `postprocess()` if
  `fort.modulegemsolver.dgemfunctionnorm > 1e-3`.
- The public `equilib_single()` catches `EquilibError` and appends an error row
  to the result.
- The public `_equilib_single()` catches `EquilibError` and suppresses it, so
  callers that inspect Fortran state need to be careful.

## EquilibResult construction

`EquilibResult.append_output()` copies the current Fortran state into
`EquilibPoint.from_fortran()`.

`EquilibPoint.from_fortran()` reads:

- system properties from `fort.modulethermoio`
- phase assemblage from `fort.modulethermo`
- phase names and database maps from `equilipy.variables`

Main copied fields:

| EquilibResult field | Source |
|---|---|
| `T`, `P` | `fort.modulethermoio` |
| `G`, `H`, `S`, `Cp` | `fort.modulethermoio` |
| `n_i` | `fort.modulethermo.dmoleselement` plus `var.cComponentNameSys` |
| `w_i` | `fort.modulethermoio.dgramelement` plus `var.cComponentNameSys` |
| `stable_phase_summary` | `fort.modulethermo.iassemblage`, `fort.modulethermo.dmolesphase`, `fort.modulethermoio.dgramphase` |
| `phase_map` | `create_phase_from_sys()` for each selected phase |

`create_phase_from_sys()` constructs `PhaseResult` objects for both solution and
compound phases. It reads species mole/mass fractions from Fortran and converts
endmember fractions back to element fractions using database stoichiometry in
`equilipy.variables`.

Current-state note:

- EquilibResult construction must happen before reset because it reads global Fortran
  output arrays.
- Some post-processing helpers still depend on `equilipy.variables` global maps,
  especially selected-system phase/species indices.

## Batch equilibrium workflow

`equilib_batch()` splits input rows into sub-batches and runs them with joblib.

```text
equilib_batch
-> _equilib_batch_input
   -> split condition dictionary into chunks of n_per_batch rows
-> _equilib_singlenode
   -> starmap_joblib(_equilib_batch, args, n_cpu)
-> combine returned EquilibResult objects
-> return EquilibResult
```

Each worker task runs `_equilib_batch()`:

```text
_equilib_batch
-> load_database
-> _dict2np
-> for each row:
   -> reduce zero-composition elements, if needed
   -> list_phases
   -> optional phase_selection
   -> var.dConditionSys = row_condition
   -> input_condition
   -> minimize
   -> EquilibResult.append_output or EquilibResult.append_error
   -> fort.resetthermo()
   -> fort.resetthermoparser()
-> return EquilibResult
```

Parallel execution details:

- `equilipy/_parallel.py` uses `joblib.Parallel`.
- The backend is `loky`, so normal parallel execution uses separate worker
  processes.
- `max_nbytes=None` disables joblib memmapping for large NumPy objects.
- If `n_cpu == 1`, the helper runs in the current process without joblib.

Current-state reset behavior:

- Each row inside `_equilib_batch()` calls `fort.resetthermo()` and
  `fort.resetthermoparser()`.
- Public `equilib_batch()` does not call `fort.resetthermoall()` after combining
  worker results.

Known Phase 1 concern:

- `tests/equilib_test04.py` has shown intermittent first-worker mismatch. The
  batch path is therefore a priority target for reset/load isolation review.

## Transition search workflow

`find_transitions()` and `find_first_transition()` repeatedly call internal
`_equilib_single()` at different temperatures, then read:

```python
fort.modulethermo.iassemblage.copy()
fort.modulethermo.dmolesphase.copy()
```

to detect changes in the stable phase set and, when available, the phase
amounts associated with each phase id.

Current-state notes:

- Transition search intentionally uses `_equilib_single()` instead of public
  `equilib_single()` so Fortran state remains available for inspection.
- The transition helper resets transient Fortran thermodynamic state around
  each temperature probe through `fort.resetthermo()` when that routine is
  available.
- If the hot endpoint is a single liquid phase, `find_transitions()` first tries
  a liquidus-specific driving-force root.  It keeps the liquid metastable,
  evaluates non-liquid phase driving forces against the liquid elemental
  potentials, and then validates the returned temperature with ordinary stable
  phase-set checks.
- For ordinary transition brackets, if the endpoint phase sets differ by
  exactly one phase, the search treats the event as a zero-phase-fraction
  appearance or disappearance.  It first solves the metastable endpoint
  assemblage, evaluates the missing phase's driving force, and uses a
  safeguarded secant/false-position estimate for the transition temperature.
  The estimated point is then validated with ordinary stable phase-set checks;
  if that validation bracket is still too wide, only the local bracket is
  refined by phase membership.
- Coupled events where one phase appears while another disappears still use the
  recursive stable-set split.

## Scheil workflow

`scheil_cooling()` simulates solidification by repeatedly calculating
equilibrium, reading the liquid composition from the previous result, lowering
temperature, and updating the remaining liquid composition.

High-level flow:

```text
scheil_cooling
-> copy user condition
-> convert mass input to mole-like basis, if needed
-> optionally find liquidus using find_liquidus_transition
-> _equilib_single for initial state
-> ScheilResult() reads initial Fortran output
-> for each cooling step:
   -> read liquid amount/composition from previous ScheilResult state
   -> lower temperature by delta_T
   -> update current_condition element amounts
   -> _equilib_single
   -> compare fort.modulethermo.iassemblage with previous assemblage
   -> if assemblage changed, refine transition temperatures with find_transitions
   -> append output to ScheilResult
   -> stop when liquid disappears or liquid amount is small
-> update_scheil_constituents
-> return ScheilResult
```

Batch Scheil functions:

- `scheil_batch()` runs `scheil_cooling()` for each row and combines
  `ScheilResult.to_dict()` outputs.
- `scheil_constituent_batch()` runs `scheil_cooling()` for each row and combines
  `ScheilResult.scheil_constituents`.
- Both use joblib through `_parallel.starmap_joblib()`.

Current-state notes:

- Scheil uses `_equilib_single()` internally, so Fortran state is intentionally
  left available after each equilibrium step.
- `scheil_cooling()` does not make a direct reset call. Instead,
  `ScheilResult.__init__()` and `ScheilResult.append_output()` call
  `fort.resetthermo()` in `finally` blocks after copying the current Fortran
  state.
- Transition refinement inside Scheil also calls `_equilib_single()` repeatedly.

## Solution-phase property workflow

`property_solution()` prepares one selected solution phase and evaluates
phase-local thermodynamic properties without running the global phase
assemblage minimizer.

```text
property_solution
-> expand_condition_species
-> _dict2np
-> for each row:
   -> _preprocess_single(..., phases=[phase_name])
   -> checkthermoinput
   -> initthermo
   -> checksystem
   -> compthermodata
   -> checkthermodata
   -> initgemsolver
   -> infer explicit endmember fractions, or read endmember_fractions={...}
   -> compexcessgibbsenergy(selected_phase)
   -> build EquilibResult / PhaseResult
   -> resetthermo
```

Current-state notes:

- The result uses the same `EquilibResult` / `PhaseResult` fields as equilibrium
  calculations, including phase `G/H/S/Cp`, endmember mole fractions,
  chemical potentials, activities, and partial `H/S/Cp`.
- If the elemental composition uniquely determines endmember fractions, the
  function evaluates the phase directly.
- For phases with more endmembers than independent elements, callers must pass
  `endmember_fractions={...}` until a true phase-local constrained
  subminimizer is added.

## Compound property workflow

`property_compound()` evaluates a stoichiometric compound species or a solution
endmember as one mole of that database species.  It builds the element
composition directly from the species stoichiometry, so `n_i` reports raw
element mole amounts rather than normalized fractions.

```text
property_compound
-> resolve compound/endmember name to one database species
-> build T/P plus one-mole species stoichiometry
-> _preprocess_single(..., phases=None)
-> checkthermoinput
-> initthermo
-> checksystem
-> compthermodata
-> checkthermodata
-> build EquilibResult / PhaseResult
-> resetthermo
```

## Reset routines and current usage

Relevant Fortran reset routines:

| Routine | Fortran file | Current purpose |
|---|---|---|
| `resetthermo` | `fsrc/3_SystemDef/ResetThermo.f90` | deallocates `ModuleThermo`, `ModuleThermoIO`, `ModuleGEMSolver`, and `ModuleSubMin` arrays |
| `resetthermoparser` | `fsrc/3_SystemDef/ResetThermoParser.f90` | deallocates `ModuleParseCS` arrays |
| `resetthermoall` | `fsrc/3_SystemDef/ResetThermoAll.f90` | calls `ResetThermo` and `ResetThermoParser` |

Observed Python usage:

| Caller | Reset behavior |
|---|---|
| `equilib_single()` | calls `fort.resetthermoall()` after result copy |
| `_equilib_single()` | no reset |
| `_equilib_batch()` | calls `fort.resetthermo()` and `fort.resetthermoparser()` after each row |
| `equilib_batch()` | no final reset after result combination |
| `find_transitions()` | resets transient thermo state around each temperature probe |
| `scheil_cooling()` | no direct reset; `ScheilResult` methods call `fort.resetthermo()` after copying state |
| `property_solution()` | calls `fort.resetthermo()` after copying each point |

This table is one of the first Phase 1 targets to verify and refine.

## Known current-state risk points

These are not conclusions yet; they are places the Phase 1 investigation should
inspect first.

| Risk area | Why it matters |
|---|---|
| Internal `_equilib_single()` leaves Fortran state live | Needed by transition/Scheil, but can leak state if callers do not reset deliberately |
| Batch worker row reset uses `resetthermo` + `resetthermoparser` | Need to confirm this is equivalent to the required clean state for every row/task |
| Public and internal workflows use different reset policies | Makes bugs hard to reproduce across single, batch, transition, and Scheil paths |
| Python `variables.py` is reused across calculations in each process | Selected-system maps can remain valid-looking even when stale |
| `EquilibResult` construction depends on Python globals plus Fortran globals | EquilibResult copy must happen before reset, and stale maps can affect phase names/compositions |
| `property_solution()` only handles explicit or caller-supplied endmember fractions today | A true phase-local constrained subminimizer is still needed for underdetermined phases |
| Transition and Scheil loops repeatedly call `_equilib_single()` | These loops intentionally inspect Fortran state, so they need explicit state ownership rules |

## Next documentation tasks

This first architecture page should be expanded into three focused pages during
Phase 1:

- `docs/developer/minimization_workflow.md`
- `docs/developer/fortran_global_state.md`
- `docs/developer/database_parser_contract.md`
- `docs/developer/python_state_and_results.md`

Those pages should replace broad descriptions with variable-level tables:

- variable owner
- allocation point
- mutation point
- reset/deallocation point
- Python writer/reader
- Fortran writer/reader
- result-copy behavior
- state-leak risk
