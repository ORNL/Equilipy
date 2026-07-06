# Database

The Database workspace inspects and edits thermochemical databases before
they are used in calculations.

## Open and inspect

**Open** loads a `.tdb`, `.eqdb`, or `.json` file. The tree shows elements,
functions, compounds, and solution phases; functions display their
piecewise Gibbs-energy ranges together with converted `G`, `Cp`, `H`, and
`S` forms.

## Edit and validate

- Edit function names and piecewise Gibbs-energy ranges in place.
- **Validate** checks the database and reports diagnostics, e.g. adjacent
  Gibbs-energy ranges that do not match at a transition temperature, or
  function names too long for portable TDB export.

Unsupported or non-portable content is reported as a warning instead of
being silently treated as calculation-ready.

## Export

**Export** writes the database as Equilipy `.eqdb`, JSON, or TDB
(normal or FactSage-style).

Before the file dialog opens, a periodic-table window lists every element
in the database (all selected by default). Deselecting elements exports an
**exact subsystem split**: only the phases, parameters, and functions
representable in the chosen subsystem are written — for example, an Al-Fe
database extracted from a 10-element assessment. The same operation is
available in scripts as [`eq.split_tdb`](../scripting/database).

Order/disorder databases are exported with the Thermo-Calc-compatible
`DIS_PART` helper convention, so the written TDB loads in Thermo-Calc,
Pandat, and FactSage-family tools as well as Equilipy.
