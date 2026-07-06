# Database

## Reading databases

```python
import equilipy as eq

DB = eq.read_dat("AlCuMgSi_ORNL_FS83.dat")   # ChemSage/FactSage .dat
DB = eq.read_tdb("AlFeSi_99Liu.tdb")         # Thermo-Calc style .tdb
```

Both return the calculation database object accepted by every calculation
function. `read_dat` auto-detects legacy (7.3) and modern (8.3+) ChemSage
dialects, including validated MQM/SUBG/SUBQ solution-model canaries.
`read_tdb` supports CEF solution models, magnetic contributions, and
order/disorder (`DIS_PART`) descriptions.

To list the phases available for a set of elements:

```python
phases = eq.list_phases(DB, ["Al", "Cu", "Si"])
```

## Splitting a database

`split_tdb` extracts an element subsystem from a multicomponent TDB — for
example the Al-Fe binary from the Al-Fe-Si assessment. The split is
thermodynamically exact: every parameter representable in the subsystem is
kept, everything else is removed.

```python
eq.split_tdb("AlFeSi_99Liu.tdb", ["Al", "Fe"],
             out="AlFe_99Liu.tdb")

DB = eq.read_tdb("AlFe_99Liu.tdb")           # ready for calculations
```

`VA` and the electron pseudo-element are kept automatically. Order/disorder
phase pairs are kept or dropped together.

## Editing and writing TDB

For database editing, parse to the editable representation and write back:

```python
from equilipy.database_ir import write_tdb

ir = eq.read_tdb("AlFeSi_99Liu.tdb", editable=True)
# ... inspect or modify ir.elements / ir.functions / ir.phases / ir.parameters
write_tdb(ir, "AlFeSi_99Liu_out.tdb")        # or export_style="factsage"
```

The writer normalizes order/disorder databases to the Thermo-Calc-compatible
`DIS_PART` helper convention so exported files load in Thermo-Calc, Pandat,
and FactSage-family tools.
