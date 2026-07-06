# Single Equilibrium Workflow

This page maps one `equilib_single()` call from Python input to Fortran
minimization and back to the Python result object. It is meant as a compact
workflow guide; detailed variable ownership is documented in
[Global State Map](global_state_map) and [Fortran Global State](fortran_global_state).

## High-Level Path

```{mermaid}
flowchart TD
    A["equilib_single(database, condition, unit, phases)"] --> B["expand_condition_species()"]
    B --> C["EquilibResult(context=capture_result_context(...))"]
    C --> D["_preprocess_single()"]
    D --> E["minimize()"]
    E --> F{"minimize succeeded?"}
    F -->|Yes| G["res.append_output()"]
    F -->|No: EquilibError| H["res.append_error()"]
    G --> I["fort.resetthermoall()"]
    H --> I
    I --> J["return EquilibResult"]
```

`equilib_single()` always returns an `EquilibResult`. A failed internal
calculation appends an error point instead of returning stale successful solver
state.

## Python Preprocessing

```{mermaid}
flowchart TD
    subgraph PREPROCESS ["_preprocess_single()"]
        P0["expand_condition_species(database, condition, unit[2])"] --> P1["load_database(database)"]
        P1 --> P2["_dict2np(condition)"]
        P2 --> P3["list_phases(database, active_elements)"]
        P3 --> P4{"phases provided?"}
        P4 -->|Yes| P5["_validated_phase_selection(phases, all_phases)"]
        P5 --> P6["phase_selection(validated_phases)"]
        P4 -->|No| P7["use all available phases"]
        P6 --> P8{"zero-composition elements present?"}
        P7 --> P8
        P8 -->|Yes| P9["drop zero elements and rerun list_phases()"]
        P9 --> P10["revalidate requested phases"]
        P8 -->|No| P11["keep original active elements"]
        P10 --> P12["var.dConditionSys = composition_condition"]
        P11 --> P12
        P12 --> P13["input_condition(unit, composition_condition, include_heat_capacity)"]
    end
```

| Step | Python function | Main state touched |
|---|---|---|
| Formula expansion | `expand_condition_species()` | Rewrites compound-formula inputs such as `Al13Fe4` into elemental amounts. |
| Database load | `load_database()` | Copies parsed database arrays into `equilipy.variables` and Fortran parser modules. |
| Phase discovery | `list_phases()` | Builds the system phase list for the active elements. |
| Phase filtering | `_validated_phase_selection()` and `phase_selection()` | Validates user-requested phases against available phases, then writes selected phase maps. |
| Condition transfer | `input_condition()` | Writes temperature, pressure, composition, units, and heat-capacity mode into Fortran input modules. |

## Fortran Minimization Sequence

`minimize()` is a Python wrapper around the Fortran minimization stages. Each
stage is followed by an `INFOThermo` check so failures become Python
`EquilibError` exceptions.

```{mermaid}
flowchart TD
    M0["minimize()"] --> M1["fort.checkthermoinput()"]
    M1 --> M2["fort.initthermo()"]
    M2 --> M3["fort.checksystem()"]
    M3 --> M4["fort.compthermodata()"]
    M4 --> M5["fort.checkthermodata()"]
    M5 --> M6["fort.gemsolver()"]
    M6 --> M7["fort.postprocess()"]
    M7 --> M8{"dGEMFunctionNorm <= 1e-3?"}
    M8 -->|Yes| M9["return to Python result capture"]
    M8 -->|No| M10["raise EquilibError"]
```

| Stage | Purpose | Important output state |
|---|---|---|
| `checkthermoinput()` | Validate and convert units. | Temperature, pressure, and composition in solver units; `INFOThermo`. |
| `initthermo()` | Initialize constants, counters, and tolerances. | Solver constants and tolerance arrays. |
| `checksystem()` | Match user elements/species/phases against the database. | System element/species/phase maps and normalized element amounts. |
| `compthermodata()` | Evaluate species and phase thermodynamic data at `T`, `P`. | Chemical potentials, Gibbs functions, stoichiometry, excess parameters. |
| `checkthermodata()` | Verify thermodynamic data and mass-balance feasibility. | `INFOThermo` success/error state. |
| `gemsolver()` | Solve the Gibbs energy minimization problem. | Stable assemblage, phase amounts, species fractions, element potentials. |
| `postprocess()` | Convert converged solver state into output arrays. | System `G`, `H`, `S`, `Cp`, phase masses, species amounts, mass fractions. |

## Result Capture and Reset

```{mermaid}
sequenceDiagram
    participant User
    participant Py as Python
    participant Var as Python globals
    participant Fort as Fortran modules
    participant Res as EquilibResult

    User->>Py: equilib_single(database, condition)
    Py->>Py: expand_condition_species()
    Py->>Res: create result context
    Py->>Var: load_database()
    Py->>Fort: input_condition()
    Py->>Fort: minimize()
    Fort-->>Py: solver state in modules
    Py->>Res: append_output() or append_error()
    Py->>Fort: resetthermoall()
    Py-->>User: EquilibResult
```

The order matters:

- Result context is captured before preprocessing mutates global state for the
  next calculation.
- Solver output is copied into `EquilibResult` before `resetthermoall()`.
- Reset happens after both successful and failed minimization attempts.

## Error Ownership

| Error source | How it should surface |
|---|---|
| Invalid input type or unknown phase | `InputConditionError` before minimization. |
| Fortran stage sets `INFOThermo` | `EquilibError` from `minimize()`, then `res.append_error()`. |
| Poor GEM convergence norm | `EquilibError` with `dGEMFunctionNorm`, then `res.append_error()`. |
| Unexpected Python/native exception | Should be fixed at the source; do not silently reuse previous solver state. |

## Developer Checklist

When changing this path, verify:

- Formula species expansion still happens before phase discovery and result
  context capture.
- Requested phases are validated after zero-composition elements are removed.
- `input_condition()` receives the same composition vector that was used to
  build the active phase system.
- `append_output()` reads Fortran state before any reset.
- Failed calculations produce error points instead of stale success points.
- `fort.resetthermoall()` still runs before returning to the caller.
