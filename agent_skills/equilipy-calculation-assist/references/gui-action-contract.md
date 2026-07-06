# GUI Action Contract

This reference is the copy-paste contract for AI assistants that need to request
Equilipy GUI calculation actions. Equilipy remains the executor; the assistant
only emits requested actions.

## Text Action Block

Use this format when the provider only supports text output:

```text
<EQUILIPY_ACTIONS>
{"actions":[{"name":"calculation.create_session"}]}
</EQUILIPY_ACTIONS>
```

The block must contain valid JSON. Comments, trailing commas, and prose inside
the JSON are not allowed.

## Structured Output Shape

Use this object shape when the provider supports JSON or schema output:

```json
{
  "reply": "I can set that up. Equilipy will execute the requested GUI actions.",
  "actions": [
    {
      "name": "calculation.create_session",
      "args": {}
    }
  ]
}
```

## JSON Schema

```json
{
  "type": "object",
  "properties": {
    "reply": {
      "type": "string",
      "description": "Brief user-visible response. Do not claim actions ran."
    },
    "actions": {
      "type": "array",
      "description": "Whitelisted Equilipy GUI actions to request.",
      "items": {
        "type": "object",
        "properties": {
          "name": {
            "type": "string",
            "enum": [
              "calculation.create_session",
              "calculation.load_database",
              "calculation.add_module",
              "calculation.set_condition",
              "calculation.set_units",
              "calculation.set_type",
              "calculation.set_batch_condition",
              "calculation.set_batch_grid",
              "solidification.set_batch_grid",
              "calculation.set_nucleation_undercooling",
              "calculation.set_result_columns",
              "calculation.select_phases",
              "calculation.calculate",
              "calculation.calculate_all"
            ]
          },
          "args": {
            "type": "object",
            "additionalProperties": true
          }
        },
        "required": ["name"],
        "additionalProperties": false
      }
    }
  },
  "required": ["reply", "actions"],
  "additionalProperties": false
}
```

## Action Catalog

### `calculation.create_session`

Create a calculation session.

```json
{"name": "calculation.create_session", "args": {"name": "AlFe study"}}
```

The `name` argument is optional.

### `calculation.load_database`

Load a `.dat` or `.tdb` database from the active calculation directory's
`database` or `databases` folder.

```json
{"name": "calculation.load_database", "args": {"file": "AlFe.dat"}}
```

If `file` is omitted and exactly one database file exists, Equilipy loads it. If
no active directory is set, Equilipy prompts the user.

### `calculation.add_module`

Add a calculation module.

```json
{"name": "calculation.add_module", "args": {"kind": "equilibrium"}}
```

Allowed `kind` values:

- `equilibrium`
- `solidification`

### `calculation.set_units`

Set module or session units.

```json
{
  "name": "calculation.set_units",
  "args": {
    "temperature": "K",
    "pressure": "atm",
    "amount": "wt%"
  }
}
```

Allowed temperature units: `K`, `C`, `F`, `R`.
Allowed pressure units: `atm`, `psi`, `bar`, `Pa`, `kPa`.
Use the GUI amount-unit label for amount units, such as `wt%`.
Set `"scope": "session"` only when the user asks to change session defaults.

### `calculation.set_type`

Set single or batch mode, and solidification model when relevant.

```json
{
  "name": "calculation.set_type",
  "args": {
    "mode": "batch",
    "model": "nucleoscheil",
    "cpu_count": 2
  }
}
```

Allowed `mode` values: `single`, `batch`.
Allowed solidification `model` values: `scheil`, `nucleoscheil`,
`equilibrium_cooling`.

### `calculation.set_condition`

Set a single-condition module input.

```json
{
  "name": "calculation.set_condition",
  "args": {
    "T": 900,
    "P": 1,
    "composition": {
      "Al": 0.8,
      "Fe": 0.2
    }
  }
}
```

For solidification modules, these optional arguments are also supported:

```json
{
  "delta_t": 0.25,
  "liquid_phase": "LIQUID",
  "from_liquidus": true
}
```

Composition may be an object or a list of entries with `species` and `amount`.

### `calculation.set_batch_condition`

Set explicit batch columns.

```json
{
  "name": "calculation.set_batch_condition",
  "args": {
    "condition": {
      "T": [700, 800],
      "P": [1, 1],
      "Al": [90, 80],
      "Fe": [10, 20]
    },
    "label": "assistant batch"
  }
}
```

All condition columns must have compatible lengths.

### `calculation.set_batch_grid`

Generate a composition grid. The balance species gets `total - sum(other
species)` for each grid point.

```json
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
    "fixed": {
      "Fe": 0.2
    },
    "total": 100
  }
}
```

Axis arrays are `[start, stop, step]`.

### `solidification.set_batch_grid`

Alias for `calculation.set_batch_grid`. Prefer `calculation.set_batch_grid`
unless the provider wants to emphasize a solidification workflow.

### `calculation.select_phases`

Select phases by name.

```json
{
  "name": "calculation.select_phases",
  "args": {
    "phases": ["LIQUID", "FCC_A1"]
  }
}
```

Only select phases that are present in the active database or clearly requested
by the user.

### `calculation.set_nucleation_undercooling`

Set critical undercooling values for nucleation Scheil workflows.

```json
{
  "name": "calculation.set_nucleation_undercooling",
  "args": {
    "undercooling": {
      "FCC_A1": 5,
      "AL13FE4": 12.5
    }
  }
}
```

Use only for solidification modules.

### `calculation.set_result_columns`

Set result table columns.

```json
{
  "name": "calculation.set_result_columns",
  "args": {
    "columns": ["T", "fs", "stable_phases"]
  }
}
```

### `calculation.calculate`

Run the selected/current module.

```json
{"name": "calculation.calculate"}
```

Only emit this when the user asks to run or calculate now.

### `calculation.calculate_all`

Run all calculation modules.

```json
{"name": "calculation.calculate_all"}
```

Only emit this when the user asks to run all modules.

## Refusal and Clarification Cases

Ask a short clarification instead of emitting actions when:

- no database file is specified and multiple database files are available
- the user gives composition without units and the expected basis is unclear
- the user asks to merge, edit, or export TDB databases
- the user asks to run calculations but required database or composition details
  are missing
- the requested action is not in the action catalog
