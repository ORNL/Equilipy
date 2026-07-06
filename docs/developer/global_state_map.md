# Global State Map

This page maps the current Equilipy calculation state across Python, f2py, and
Fortran. It is written for debugging `equilib_batch()` flakiness and process
reuse. It describes the current implementation, not a proposed reset strategy.

## State Layers

Equilipy currently has three mutable state layers:

```text
Python parsed database and selected-system state
-> f2py-visible Fortran module state
-> Python result objects copied from Python/Fortran state
```

| Layer | Owner | Lifetime | Main contents | Reset / rebuild owner |
|---|---|---|---|---|
| Parsed database state | `equilipy.variables` | Python process-global | Arrays parsed from `.dat`, including database species, phases, stoichiometry, Gibbs parameters, and model metadata. | Reassigned by `read_dat()` / `load_database()`; no full Python-state reset helper exists. |
| Selected-system state | `equilipy.variables` | Python process-global | Active elements, system species, system phases, phase maps, species maps, and selected phase maps. | Rebuilt by `list_phases()` and optionally narrowed by `phase_selection()`. |
| Input-condition state | `equilipy.variables` + `fort.modulethermoio` | Python and Fortran process-global | Current condition vector, units, temperature, pressure, element masses, heat-capacity flag. | Python caller writes `var.dConditionSys`; `input_condition()` writes Fortran input fields. |
| Fortran parser state | `fort.moduleparsecs` | Process-global inside each Python process | Database arrays copied from `equilipy.variables` through f2py. | `_pyvar2fvar()` writes; `fort.resetthermoparser()` deallocates parser arrays. |
| Fortran thermodynamic state | `fort.modulethermo` | Process-global inside each Python process | Active thermodynamic system, species/phase arrays, assemblage, mole amounts, chemical potentials, Gibbs arrays. | Allocated/reallocated by `checksystem`, `initthermo`, `compthermodata`, GEM/PEA routines; mostly cleared by `fort.resetthermo()`. |
| Fortran solver state | `fort.modulegemsolver`, `fort.modulesubmin` | Process-global inside each Python process | Iteration counters, convergence norms, line-search state, subminimization work arrays. | Initializers deallocate/reallocate many arrays; `fort.resetthermo()` clears major arrays but not every submin local-work array. |
| Result state | `equilipy.results` | User-owned Python object | `EquilibResult`, `EquilibPoint`, `PhaseResult`, `ScheilResult`, and captured result context. | No reset. Result objects must copy what they need before Fortran state changes. |

## Public Calculation Workflows

### `read_dat()`

| Step | State read | State written | Fortran state |
|---|---|---|---|
| Parse database file | `.dat` file | `equilipy.variables` parser/database arrays | None in the normal Python parser path |
| Return dictionary | `equilipy.variables` | returned database dictionary | None |

Risk: the returned dictionary is the safe portable object. `equilipy.variables`
remains mutable process-global state after parsing.

### `load_database(database)`

| Step | State read | State written | Fortran state |
|---|---|---|---|
| Copy dictionary fields | input database dictionary | `equilipy.variables` attributes | None directly |

Risk: this rebuilds Python parser/database state only. It does not reset or
clear Fortran module arrays.

### `list_phases(database, elements)`

Current flow:

```text
load_database(database)
utils._pyvar2fvar(var)
system_check(elements)
build var.PhaseNameSys
return var.cPhaseNameSys
```

| State | Effect |
|---|---|
| Python parser/database state | Reassigned through `load_database()` |
| `fort.moduleparsecs` | Written by `_pyvar2fvar()` |
| Python selected-system state | Written by `system_check()` and the phase-map construction |
| Fortran selected-system state | `system_check()` writes selected-system arrays used by later Fortran routines |

Risk: this is the main bridge from database/elements to selected system. If a
previous worker leaves incompatible Fortran state, `_pyvar2fvar()` and
`system_check()` must overwrite enough state before `minimize()`.

Phase 6 note: `list_phases()` builds Python-side system maps, not the final
active Fortran `ModuleThermo` system. The final compact active system is rebuilt
inside `CheckSystem()` and `CompThermoData()` during `minimize()`.

### `phase_selection(phases)`

Current flow:

```text
validate against var.PhaseNameSys
build narrowed iSys2DBSoln / iSys2DBComp / iSys2DBSpecies
rewrite var.cPhaseNameSys and var.PhaseNameSys
write fort.moduleparsecs.iphasecs
write fort.modulethermo.isolnps
```

Risk: this function assumes `list_phases()` has already built
`var.PhaseNameSys` for the active element system. `equilib_single()` validates
and reorders requested phases through `_validated_phase_selection()`;
`_equilib_batch()` currently calls `phase_selection()` directly.

Phase 6 note: `phase_selection()` also mutates `var.iPhaseCS` from the raw
database phase-ownership vector into a selected-phase mask. The next
`load_database()` restores raw parser state from the database dictionary.

### `input_condition(unit, composition)`

| Python state read | Fortran state written |
|---|---|
| `var.iElementSys`, `var.iElementSysIndex`, caller-provided composition | `modulethermoio.dtemperaturediff`, unit strings, `dtemperature`, `dpressure`, `delementmass`, `iprintresultsmode` |

Risk: `input_condition()` resets `delementmass` to zeros before writing current
composition. It depends on current selected-system element ordering.

### `minimize()`

Current Fortran sequence:

```text
checkthermoinput()
initthermo()
checksystem()
compthermodata()
checkthermodata()
gemsolver()
postprocess()
```

Python checks `fort.modulegemsolver.dgemfunctionnorm` after `postprocess()` and
raises `EquilibError` if it is too large.
Python also checks `fort.modulethermoio.infothermo` after each Fortran stage and
raises `EquilibError` immediately if a stage sets a nonzero status without
throwing through f2py.

Risk: `minimize()` does not own reset. It assumes the active database, selected
system, phase selection, and input condition are already coherent.

Phase 6 note: `CheckSystem()` and `CompThermoData()` are the key Fortran
system-setting routines. A future TDB/IR backend can initially target the
`ModuleParseCS` fields consumed by those routines; bypassing them means
replacing the current system-setting contract.

## Equilibrium Entry Points

### `equilib_single()`

Current flow:

```text
expand_condition_species()
capture_result_context()
_preprocess_single()
    expand_condition_species()
    load_database()
    _dict2np()
    list_phases()
    _validated_phase_selection() if phases are provided
    phase_selection() if phases are provided
    zero-composition branch may call list_phases() again
    input_condition()
minimize()
EquilibResult.append_output() or append_error()
fort.resetthermoall()
```

Important state ownership:

- setup is owned by `_preprocess_single()`;
- result capture happens before `resetthermoall()`;
- after return, Fortran parser, thermodynamic, and reinit arrays should be
  deallocated, but Python `equilipy.variables` remains populated.

### `_equilib_single()`

Current flow:

```text
_preprocess_single()
minimize()
return None
```

Risk: no result object and no reset. The caller owns result capture and cleanup.

### `equilib_batch()`

Current flow:

```text
expand_condition_species()
_equilib_batch_input(..., n_per_batch)
_equilib_singlenode()
    starmap_joblib(_equilib_batch, args, n_cpu)
append worker results into one EquilibResult
recapture result context
```

`n_cpu=1` uses the same `_equilib_batch()` worker function in the main process.
`n_cpu>1` uses joblib/loky worker processes.

### `_equilib_batch()`

Current flow:

```text
load_database()
expand_condition_species()
_dict2np()
for each row:
    if a component is zero:
        reduce active elements
        list_phases(database, reduced_elements)
        phase_selection(phases) if phases are provided
    else:
        list_phases(database, all_elements)
        phase_selection(phases) if phases are provided
    var.dConditionSys = row_condition
    input_condition()
    minimize()
    append_output() or append_error()
    fort.resetthermo()
    fort.resetthermoparser()
return worker EquilibResult
```

Important current differences from `equilib_single()`:

- direct `phase_selection(phases)` instead of `_validated_phase_selection()`;
- no `resetthermoall()` after the worker result is complete;
- row cleanup is `resetthermo()` + `resetthermoparser()`;
- Python `equilipy.variables` remains populated after each row and worker task.

These differences are investigation targets. They are not proof that any one
reset is required.

## Cooling and Transition Workflows

| Workflow | Equilibrium dependency | State risk |
|---|---|---|
| `find_transitions()` | Repeated equilibrium calls over a temperature interval | Failed trial calculations must not leave stale successful state visible to later trials. |
| `scheil_cooling()` | Repeated equilibrium-like steps while carrying liquid state | Intentionally reuses selected database/system state; broad parser reset inside the path may be wrong. |
| `equilib_cooling()` | Repeated equilibrium calculations over a cooling grid | Should eventually follow the same isolation rules discovered for batch equilibrium. |
| `nucleoscheil_cooling()` | Scheil/equilibrium path plus nucleation decisions | Adds undercooling state and transition calls; inherits equilibrium state risks. |

## Reset Coverage

| Reset call | Source | Clears | Does not clear | Current Python use |
|---|---|---|---|---|
| `fort.resetthermo()` | `ResetThermo.f90` | Most `ModuleThermo`, major `ModuleSubMin`, and `ModuleThermoIO` gram output arrays; `ModuleThermoIO.INFOThermo`; `ModuleGEMSolver` scalar counters, norms, indexes, and logical flags. | `ModuleParseCS`; Python `equilipy.variables`; some locally managed submin arrays; `ModuleThermoIO.dPair/cPair`. | `_equilib_batch()` after each row; Scheil result paths. |
| `fort.resetthermoparser()` | `ResetThermoParser.f90` | Main `ModuleParseCS` parser/database arrays; parser-side `INFO`; `ModuleThermoIO.INFOThermo`. | Thermodynamic solver arrays; Python `equilipy.variables`. | `_equilib_batch()` after each row. |
| `fort.resetthermoall()` | `ResetThermoAll.f90` | Calls `ResetThermo` and `ResetThermoParser`. | Python `equilipy.variables`; any Fortran array not covered by those routines. | `equilib_single()` after result capture. |

## High-Risk State for Batch Flakiness

| State | Why it matters | Current cleanup path | Investigation question |
|---|---|---|---|
| `var.PhaseNameSys`, `var.cPhaseNameSys` | `phase_selection()` validates against current Python phase map. | Rebuilt by `list_phases()` and narrowed by `phase_selection()`. | Does batch ever validate against a stale or differently ordered map? |
| `var.iSys2DBSoln`, `var.iSys2DBComp`, `var.iSys2DBSpecies` | Result capture and phase setup rely on these maps. | Rebuilt by `system_check()` / `phase_selection()`. | Are maps identical for duplicate rows and across workers? |
| `fort.moduleparsecs.iphasecs` | Fortran parser phase pass mask. | Written by `_pyvar2fvar()` and `phase_selection()`; reset by parser reset. | Does each row write the intended mask before `minimize()`? |
| `fort.modulethermo.isolnps` | Fortran solution phase selection. | Written by `phase_selection()`; deallocated by `resetthermo()`. | Does stale `isolnps` survive when no phases are provided? |
| `fort.modulethermo.iassemblage` | Stable assemblage read by post-processing and result capture. | Allocated by GEM init paths; deallocated by `resetthermo()`. | Is bad output already present immediately after `minimize()`? |
| `fort.modulegemsolver.dgemfunctionnorm` | Success/failure criterion read by Python. | Reset to zero by `resetthermo()` and overwritten during solver paths. | Is it nonzero only because the current solve produced it? |
| `ModuleSubMin` local arrays | Some are cleaned by `SubMinInit` but not clearly by `resetthermo()`. | Mostly initializer-owned; partial global reset coverage. | Can a worker reuse stale submin vectors if a path exits before `SubMinInit` cleanup? |
| `ModuleThermoIO.dPair/cPair` | SUBG pair output arrays are not cleared by `resetthermo()`. | Deallocated before next `CalculateCompositionSUBG` allocation. | Can stale pair output affect result capture or later calculations? |

## Phase 3 Reproduction Requirement

The backend comparison should record:

- execution mode;
- worker PID;
- task index;
- number of rows in the task;
- input case name;
- `G`, `H`, `S`, `Cp`;
- stable phase names per point;
- whether the task returned an error point;
- whether duplicate rows disagree.

The comparison must include joblib/loky and multiprocessing-style process
execution. If joblib differs but multiprocessing does not, the next step is to
identify which process-global variable survives in a reused joblib worker.
