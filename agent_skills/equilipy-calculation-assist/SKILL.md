---
name: equilipy-calculation-assist
description: Guide AI assistants that operate inside or alongside the Equilipy GUI to set up, populate, and run calculation workflows through the whitelisted GUI action bridge. Use when a user asks to create calculation sessions, load .dat or .tdb databases, add equilibrium or solidification modules, set temperature/pressure/composition conditions, configure batch grids, choose phases, set units, run calculations, or explain what GUI action should happen next.
---

# Equilipy Calculation Assist

## Core Rule

Treat Equilipy as the executor. The assistant may plan, explain, and request
whitelisted GUI actions, but it must not claim that a GUI action has already
happened. Equilipy executes requested actions and reports the result.

Use the action bridge when the user asks to create, set up, populate, or run a
calculation. For conversation-only requests, explain using the current GUI
context and do not emit an action block.

## Workflow

1. Read the GUI context first: workspace, active calculation directory, project
   database files, loaded sessions, selected module, units, and diagnostics.
2. Identify the calculation target:
   - equilibrium module for equilibrium at fixed T/P/composition
   - solidification module for Scheil, nucleation Scheil, or equilibrium cooling
   - batch mode when the user gives ranges, grids, CSV-like columns, or multiple
     points
3. Request only enough actions to satisfy the user's intent. Do not invent a
   database file name if the context lists multiple candidates and the user did
   not choose one.
4. Load a database before setting batch grids or selecting phases.
5. Set units before composition when the user specifies wt%, mole fraction, or
   another amount basis.
6. Select phases only when the user asks for explicit phase selection or when a
   solidification workflow requires a known liquid phase.
7. Run `calculation.calculate` only when the user explicitly asks to run now.

## Action Output

For text-only providers, put requested actions at the end of the reply inside
the exact action tags:

```text
<EQUILIPY_ACTIONS>
{"actions":[{"name":"calculation.create_session"}]}
</EQUILIPY_ACTIONS>
```

For providers with structured output, return an object with `reply` and
`actions` fields. Use the schema in
`references/gui-action-contract.md` when the integration supports JSON schema.

The visible reply should be short and honest:

```text
I can set that up. Equilipy will create the session, load the database, and add
an equilibrium module.
```

## Action Order

Use this default order for new calculations:

1. `calculation.create_session`
2. `calculation.load_database`
3. `calculation.add_module`
4. `calculation.set_units`
5. `calculation.set_type`
6. `calculation.set_condition` or `calculation.set_batch_condition` /
   `calculation.set_batch_grid`
7. `calculation.select_phases`
8. `calculation.set_nucleation_undercooling`
9. `calculation.set_result_columns`
10. `calculation.calculate`

Skip steps that are already satisfied by the current GUI context or not relevant
to the user's request.

## Common Recipes

### Single equilibrium setup

Use when the user asks for an equilibrium calculation at one condition:

```json
{
  "actions": [
    {"name": "calculation.create_session"},
    {"name": "calculation.load_database", "args": {"file": "AlFe.dat"}},
    {"name": "calculation.add_module", "args": {"kind": "equilibrium"}},
    {
      "name": "calculation.set_condition",
      "args": {
        "T": 900,
        "P": 1,
        "composition": {"Al": 0.8, "Fe": 0.2}
      }
    }
  ]
}
```

### Solidification setup

Use `kind: solidification`, set the liquid phase when known, and set
`from_liquidus` only when requested or implied by a liquidus-start workflow:

```json
{
  "actions": [
    {"name": "calculation.create_session"},
    {"name": "calculation.load_database", "args": {"file": "AlCuMgSi.tdb"}},
    {"name": "calculation.add_module", "args": {"kind": "solidification"}},
    {
      "name": "calculation.set_type",
      "args": {"mode": "single", "model": "scheil"}
    },
    {
      "name": "calculation.set_condition",
      "args": {
        "T": 700,
        "P": 1,
        "delta_t": 0.25,
        "liquid_phase": "LIQUID",
        "from_liquidus": true,
        "composition": {"Al": 90, "Cu": 4, "Mg": 1, "Si": 5}
      }
    }
  ]
}
```

### Batch grid

Use `calculation.set_batch_grid` when the user describes composition ranges.
The balance species receives the residual to `total`.

```json
{
  "actions": [
    {"name": "calculation.set_units", "args": {"amount": "wt%"}},
    {
      "name": "calculation.set_batch_grid",
      "args": {
        "T": 700,
        "P": 1,
        "balance": "Al",
        "axes": {
          "Cu": [0, 10, 1],
          "Mg": [0, 3, 0.5],
          "Si": [0, 5, 0.5]
        },
        "total": 100
      }
    }
  ]
}
```

## Guardrails

- Use only action names listed in `references/gui-action-contract.md`.
- Do not ask the model to manipulate Qt widgets directly.
- Do not call file-system tools to mutate `.eq`, `.dat`, `.tdb`, or session
  files when the GUI action bridge can do the task.
- If the database folder is missing or no database file exists, request
  `calculation.load_database`; Equilipy will prompt or create the folder.
- If the user asks for database editing, TDB parsing, Gibbs implementation, or
  TDB merging, do not use this calculation skill beyond loading a database for a
  calculation. Those tasks need separate skills.
- If the request is ambiguous about database file, calculation type, units, or
  whether to run immediately, ask a short clarification instead of guessing.

## References

- Read `references/gui-action-contract.md` when implementing a provider adapter,
  validating action JSON, or copy-pasting the action schema into another AI
  assistant.
