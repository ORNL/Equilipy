# TDB Format Developer Notes

This page documents the subset of the Thermo-Calc Database (TDB) syntax that is
important for Equilipy development.  It is not intended to replace commercial
software manuals.  Its purpose is to make Equilipy's reader, writer, validator,
GUI database editor, and runtime converter behave consistently.

The central rule is:

> A TDB file is a sequence of command records ending in `!`.  Each record
> defines either system metadata, symbols/functions, phase models, constituent
> arrays, or thermodynamic parameters.

Equilipy parses those records into `DatabaseIR`, validates/normalizes that IR,
and can export a stricter TDB form for round-trip use.

## Record Structure

TDB records are free-form text records terminated by `!`:

```tdb
FUNCTION GHSERAL 298.15 -7976.15+137.093038*T-24.3671976*T*LN(T);
    700 Y -11276.23976923516+223.048446*T-38.5844296*T*LN(T);
    933.47 Y -11278.37838094356+188.684153*T-31.748192*T*LN(T);
    2900 N !
```

Physical lines are only formatting.  The command above is one record because it
has one terminating `!`.

Comments start with `$` and are not executable TDB commands:

```tdb
$ This is a comment.
$PARAMETER G(FCC_A1,AL:VA;0) 298.15 GHSERAL#; 6000 N !
```

Equilipy must not parse commented-out commands as real phases, functions, or
parameters.

## Main Command Families

### `ELEMENT`

Declares a database element, its standard element reference (SER) phase, molar
mass, H298, and S298.

```tdb
ELEMENT Al FCC_A1 26.982 4577.3 28.322 !
ELEMENT Va VACUUM 0 0 0 !
ELEMENT /- ELECTRON_GAS 0 0 0 !
```

Role in Equilipy:

- creates element/species identity
- defines available components for phase filtering and calculations
- supplies atomic masses for unit conversion

Writer policy:

- element symbols are exported in conventional form, for example `Al`, `Fe`,
  `Si`, and `Va`
- `Va` should remain `Va` in the element table
- the default-system command may still use the traditional uppercase `VA`

### `SPECIES`

Declares a species that is not simply a single element.  Many alloy TDBs omit
explicit `SPECIES` records and rely only on `ELEMENT`.

```tdb
SPECIES AL+3 AL1/+3 !
```

Role in Equilipy:

- used for charged species or non-element species
- relevant for ionic or complex solution models

Current limitation:

- most current Equilipy TDB work focuses on metallic CEF/SUBL-style phases

### `FUNCTION`

Defines a named expression over one or more temperature intervals.

```tdb
FUNCTION GHSERFE 298.15 1225.7+124.134*T-23.5143*T*LN(T);
    1811 Y -25383.5748952178+299.31255*T-46*T*LN(T);
    6000 N !
```

Fields:

- function name: `GHSERFE`
- first lower temperature: `298.15`
- Gibbs expression: `1225.7+...`
- upper temperature: `1811`
- continuation flag: `Y` means another range follows; `N` means final range

Function references use `#`:

```tdb
PARAMETER G(BCC_A2,FE:VA;0) 298.15 GHSERFE#; 6000 N !
```

Role in Equilipy:

- stored as `FunctionDefinition`
- expanded/evaluated during runtime conversion
- reused for compact writer output

Important writer rule:

- numeric `D` exponents are normalized to `E`, but only in numeric literals
- function names such as `LDF0MNNI` must not be rewritten to `LEF0MNNI`

### System Defaults: `TYPE_DEF`, `DEFINE_SYSTEM_DEFAULT`, `DEFAULT_COMMAND`

These commands configure phase behavior and default species.

Common required defaults:

```tdb
TYPE_DEF % SEQ * !
DEFINE_SYSTEM_DEFAULT SPECIE 2 !
DEFAULT_COMMAND DEF_SYS_ELEMENT VA /-!
```

`TYPE_DEF` can also assign magnetic models or order-disorder relationships.

Magnetic example:

```tdb
TYPE_DEF A GES AMEND_PHASE_DESCRIPTION @ MAGNETIC -1 0.4 !
TYPE_DEF B GES AMEND_PHASE_DESCRIPTION @ MAGNETIC -3 0.28 !
```

The `@` placeholder allows one magnetic type definition to be reused by phases
with the same magnetic model parameters.

Order-disorder example:

```tdb
TYPE_DEF C GES AMEND_PHASE_DESCRIPTION B2_BCC DIS_PART BCC_A2 !
TYPE_DEF D GES AMEND_PHASE_DESCRIPTION FCC_4SL DIS_PART FCC_A1 !
```

Role in Equilipy:

- magnetic type definitions are used to attach Curie/Neel and magnetic moment
  behavior
- `DIS_PART` records link ordered phases to their disordered part
- runtime import canonicalizes known duplicate disordered aliases for
  calculation
- database export supports two order/disorder alias modes: `preserve` and
  `canonical`

Alias policy:

```text
A2_BCC == BCC_A2
A1_FCC == FCC_A1
```

Equilipy treats these as the same physical disordered phases.  The calculation
runtime folds them to the canonical phase names by default.

Writer modes:

- `preserve` keeps concrete helper aliases used by `DIS_PART`.  If the source
  only has canonical `BCC_A2` or `FCC_A1`, this mode synthesizes `A2_BCC` or
  `A1_FCC` by copying the canonical phase and its parameters, then points
  `DIS_PART` at the helper alias.
- `canonical` collapses helper aliases and points `DIS_PART` at `BCC_A2` or
  `FCC_A1`.

### `PHASE`

Declares a phase, type codes, sublattice count, and site ratios.

```tdb
PHASE BCC_A2 %B 2 1 3 !
PHASE B2_BCC %BC 3 0.5 0.5 3 !
PHASE FCC_4SL:F %AD 5 0.25 0.25 0.25 0.25 1 !
```

Fields:

- phase name: `BCC_A2`
- optional phase option suffix: `:F`
- type codes: `%B`
- number of sublattices: `2`
- site ratios: `1 3`

Role in Equilipy:

- creates a `Phase`
- determines CEF/SUBL runtime shape
- determines endmember enumeration

### `CONSTITUENT` / `CONSTITUENTS` / `CONST`

Declares allowed constituents for each sublattice.

```tdb
CONSTITUENT BCC_A2 :Al,Fe,Ni:Va: !
CONSTITUENT B2_BCC :Al,Fe,Ni:Al,Fe,Ni:Va: !
CONSTITUENT FCC_4SL :Al,Fe,Ni:Al,Fe,Ni:Al,Fe,Ni:Al,Fe,Ni:Va: !
```

Role in Equilipy:

- defines allowed species per sublattice
- determines possible endmember configurations
- determines which interactions are valid

### `PARAMETER` / `PARAM` / `PAR`

Defines thermodynamic parameters.  Equilipy accepts common abbreviations:

```tdb
PARAMETER G(BCC_A2,FE:VA;0) 298.15 GHSERFE#; 6000 N !
PARAM L(BCC_A2,AL,FE:VA;0) 298.15 -122360+31.6*T; 6000 N !
PAR TC(BCC_A2,FE:VA;0) 298.15 1043; 6000 N !
```

Parameter family roles:

| Family | Role |
| --- | --- |
| `G` | Endmember Gibbs energy or ordered-phase Gibbs/order parameter |
| `L` | Excess interaction parameter |
| `TC` | Magnetic Curie/Neel temperature parameter |
| `BMAGN` | Magnetic moment parameter |
| `BMAG`, `BM` | Accepted aliases; writer exports as `BMAGN` |

The target section is phase-specific:

```tdb
G(B2_BCC,Al:Fe:Va;0)
```

means sublattice 1 has `Al`, sublattice 2 has `Fe`, sublattice 3 has `Va`.

The final `;0` is the parameter order.  For Redlich-Kister interactions, `;1`,
`;2`, etc. indicate polynomial order.

### Ternary CEF Interaction Orders

For CEF/SUBL phases, a ternary interaction on one mixed sublattice has two
different TDB meanings depending on whether the target is supplied as a lone
record or as a complete order set.

A lone ternary `;0` record:

```tdb
PAR L(TEST,Al,Co,Cr:Va;0) 298.15 3000; 6000 N !
```

is evaluated as:

```text
y_Al * y_Co * y_Cr * L0
```

It must remain one runtime interaction row.  It is **not** equivalent to a
synthetic `;0`, `;1`, `;2` set.

A complete ternary set:

```tdb
PAR L(TEST,Al,Co,Cr:Va;0) 298.15 L0; 6000 N !
PAR L(TEST,Al,Co,Cr:Va;1) 298.15 L1; 6000 N !
PAR L(TEST,Al,Co,Cr:Va;2) 298.15 L2; 6000 N !
```

is evaluated in TC-style composition-dependent form.  For a strictly ternary
mixed sublattice this is:

```text
y_Al * y_Co * y_Cr * (y_Al*L0 + y_Co*L1 + y_Cr*L2)
```

In a higher-component sublattice, the complete ternary set uses Muggianu
redistribution over the inactive constituents:

```text
y_Al * y_Co * y_Cr *
(
  (y_Al + (1 - y_Al - y_Co - y_Cr)/3) * L0
  + (y_Co + (1 - y_Al - y_Co - y_Cr)/3) * L1
  + (y_Cr + (1 - y_Al - y_Co - y_Cr)/3) * L2
)
```

Equilipy therefore carries an internal marker for complete ternary CEF sets.
The marker is a runtime encoding detail only; the TDB writer should preserve
the original `;0`/`;1`/`;2` records rather than writing an artificial syntax.

## Endmembers

An endmember is one concrete choice of one constituent per sublattice.

Example:

```tdb
PHASE B2_BCC %BC 3 0.5 0.5 3 !
CONSTITUENT B2_BCC :Al,Fe:Al,Fe:Va: !
```

Concrete endmembers include:

```text
Al:Al:Va
Al:Fe:Va
Fe:Al:Va
Fe:Fe:Va
```

Each concrete endmember normally needs a `G(...;0)` record, directly or through
a valid wildcard shorthand.

Equilipy validation behavior:

- missing endmembers are diagnosed
- all-vacancy endmembers, such as `VA:VA:VA`, are generated as zero-energy
  compatibility records with a warning
- order-disorder pairs may need generated zero compatibility records for
  ordered endmembers that commercial software silently tolerates

## Wildcards and Comma Lists

Wildcards and comma lists are not equivalent.

### Concrete Endmember

```tdb
PAR G(FCC_4SL,Al:Ni:Ni:Ni:Va;0) 298.15 -1000; 6000 N !
```

This covers only:

```text
Al:Ni:Ni:Ni:Va
```

### Wildcard Endmember Shorthand

```tdb
PAR G(FCC_4SL,*:Ni:Ni:Ni:Va;0) 298.15 -1000; 6000 N !
```

This can cover concrete endmembers where the first sublattice is any allowed
species:

```text
Al:Ni:Ni:Ni:Va
Fe:Ni:Ni:Ni:Va
Co:Ni:Ni:Ni:Va
...
```

Equilipy treats single-site wildcard `G` records as valid endmember shorthand.
For export, wildcard records may be expanded or symmetry-reduced depending on
writer options and phase symmetry.

### Comma-List `G(...)` Records

```tdb
PAR G(FCC_4SL,Mn,Ni:Mn,Ni:*:*:Va;0) 298.15 SFMNNI#; 6000 N !
```

This is not one concrete endmember.  `Mn,Ni` means a mixed constituent set on
that sublattice.  In ordered CEF phases, these records behave like
interaction/order terms.

Equilipy rule:

- `G(A:B:C;0)` covers a concrete endmember
- `G(*:B:C;0)` can cover matching concrete endmembers
- `G(A,B:B:C;0)` does not cover concrete endmembers; it is treated as an
  interaction/order parameter

This distinction matters.  Treating comma-list `G(...)` as concrete endmember
coverage can suppress necessary generated endmembers and shift phase
energetics.

## `:F` Phase Option and Symmetric Sublattices

The `:F` option marks equivalent sublattices in a phase.

Example:

```tdb
PHASE FCC_4SL:F %AD 5 0.25 0.25 0.25 0.25 1 !
CONSTITUENT FCC_4SL :Al,Fe,Ni:Al,Fe,Ni:Al,Fe,Ni:Al,Fe,Ni:Va: !
```

The first four sublattices are equivalent FCC sublattices.

If a TDB has:

```tdb
PAR L(FCC_4SL,Al,Ni:*:*:*:Va;0) 298.15 10000; 6000 N !
```

the runtime must apply that interaction to equivalent positions:

```text
Al,Ni:*:*:*:Va
*:Al,Ni:*:*:Va
*:*:Al,Ni:*:Va
*:*:*:Al,Ni:Va
```

Important: this is interaction expansion, not endmember population.

The original parameter value is preserved for each equivalent interaction.
Equilipy must not divide it by the number of equivalent positions and must not
site-ratio-scale it to `10000*0.25`.  Commercial behavior is that the expanded
equivalent terms keep the original value.

## Order-Disorder Phases

Order-disorder pairs are declared through `DIS_PART`.  Some databases use the
canonical disordered phase name:

```tdb
TYPE_DEF C GES AMEND_PHASE_DESCRIPTION B2_BCC DIS_PART BCC_A2 !
TYPE_DEF D GES AMEND_PHASE_DESCRIPTION FCC_4SL DIS_PART FCC_A1 !
```

Others use helper aliases:

```tdb
TYPE_DEF C GES AMEND_PHASE_DESCRIPTION B2_BCC DIS_PART A2_BCC !
TYPE_DEF D GES AMEND_PHASE_DESCRIPTION FCC_4SL DIS_PART A1_FCC !
```

Meaning:

- `B2_BCC` is ordered
- `BCC_A2` is the associated disordered part
- `FCC_4SL` is ordered
- `FCC_A1` is the associated disordered part

Equilipy runtime canonical alias policy:

```text
A2_BCC -> BCC_A2
A1_FCC -> FCC_A1
```

Runtime import folds helper aliases to canonical names for calculation, so
`A2_BCC` and `BCC_A2` are not separate physical phases in Equilipy.  Writer
export is controlled by `order_disorder_alias_mode`.

In `preserve` mode, if the source TDB declared and referenced `A2_BCC` or
`A1_FCC`, the writer keeps those records and keeps `DIS_PART` pointing at them:

```tdb
TYPE_DEF C GES AMEND_PHASE_DESCRIPTION B2_BCC DIS_PART A2_BCC !
TYPE_DEF D GES AMEND_PHASE_DESCRIPTION FCC_4SL DIS_PART A1_FCC !
PHASE A2_BCC %B 2 1 3 !
PHASE A1_FCC %A 2 1 1 !
```

If the source only has canonical `BCC_A2` or `FCC_A1`, `preserve` mode
synthesizes the missing helper phase and copies the canonical parameters:

```tdb
TYPE_DEF C GES AMEND_PHASE_DESCRIPTION B2_BCC DIS_PART A2_BCC !
PHASE BCC_A2 %B 2 1 3 !
PHASE A2_BCC %B 2 1 3 !
PAR G(BCC_A2,Al:Va;0) 298.15 GHSERAL#; 6000 N !
PAR G(A2_BCC,Al:Va;0) 298.15 GHSERAL#; 6000 N !
```

In `canonical` mode, the same order-disorder pair is emitted without the helper
alias:

```tdb
TYPE_DEF C GES AMEND_PHASE_DESCRIPTION B2_BCC DIS_PART BCC_A2 !
PHASE BCC_A2 %B 2 1 3 !
PAR G(BCC_A2,Al:Va;0) 298.15 GHSERAL#; 6000 N !
```

### Partitioned CEF Equation Used by `DIS_PART`

For an ordered phase `O` with disordered partner `D`, the ordered sublattice
site fractions are collapsed onto the disordered sublattices before the
disordered reference is evaluated.  For each disordered sublattice group `d`:

```text
ybar[d,k] = sum(s -> d, a[s] * y[s,k]) / sum(s -> d, a[s])
```

where `a[s]` is the ordered sublattice site ratio, `y[s,k]` is the ordered
site fraction, and `s -> d` means ordered sublattice `s` maps to disordered
sublattice `d`.

The partitioned ordered Gibbs energy is:

```text
G_part^O(y) = G^D(ybar) + G_corr^O(y) - G_corr^O(y = ybar)
```

Equivalently, when the raw ordered CEF expression is used as the correction
function:

```text
G_part^O(y) = G^O_raw(y) + G^D(ybar) - G^O_raw(y = ybar)
```

This is the Hillert CEF order-disorder construction.  The value that must
vanish in the random/disordered state is the ordering contribution:

```text
G_corr^O(y) - G_corr^O(y = ybar)
```

not the disordered-minus-random bracket by itself.  At the random state,
`G_part^O(y = ybar) = G^D(ybar)`.

For phases with vacancy sublattices, the same collapse is applied per mapped
disordered sublattice group.  For example, in a B2/A2 pair:

```tdb
PHASE BCC_A2 %B 2 1 3 !
PHASE BCC_B2 %BC 3 0.5 0.5 3 !
```

the first two ordered B2 sublattices collapse into the substitutional A2
sublattice, while the pure vacancy sublattice maps to the vacancy sublattice.

Magnetic order-disorder needs extra care.  Hillert-Jarl magnetic Gibbs energy
is evaluated from model variables such as the Curie/Neel temperature (`TC`) and
magnetic moment (`BMAGN`) and is nonlinear in those variables.  Equilipy
therefore partitions the magnetic model variables:

```text
TC_part   = TC^O(y)   + TC^D(ybar)   - TC^O(y = ybar)
BMAG_part = BMAG^O(y) + BMAG^D(ybar) - BMAG^O(y = ybar)
```

and then evaluates the Hillert-Jarl magnetic Gibbs contribution from those
partitioned variables.  Partitioning already-evaluated magnetic Gibbs values is
not generally equivalent because the magnetic Gibbs function is nonlinear in
`TC` and `BMAGN`.  The corresponding partial-molar `G`, `H`, `S`, and `Cp`
terms are analytical, not numerical finite differences.

## Writer Normalization Examples

### Repeated Range Lower Bounds

Bad old development export:

```tdb
PAR G(ALCRFE_E,Al:Cr:Fe;0) 298.15 expr1; 700 Y 700 expr2;
    933.47 Y 933.47 expr3; 6000 N !
```

Correct export:

```tdb
PAR G(ALCRFE_E,Al:Cr:Fe;0) 298.15 expr1; 700 Y expr2;
    933.47 Y expr3; 6000 N !
```

Only the first range carries an explicit lower bound.  Later lower bounds are
implied by the previous upper bound.

Reader repair policy:

- if a segment body starts with the exact previous lower bound as a standalone
  token, strip it
- do not strip legitimate coefficients such as `700*T`

### Legacy `,,` Lower Bound and Reference Tags

Legacy input:

```tdb
PAR G(FCC_A1,Al:Va),, +GHSERAL; 2900 N 91Din !
```

Normalized internal/export form:

```tdb
PAR G(FCC_A1,Al:Va;0) 298.15 GHSERAL#; 2900 N !
```

Rules:

- leading `,,` means default lower temperature `298.15`
- missing upper bound forms such as `;,, N` become `; 6000 N`
- trailing literature tags such as `91Din` are removed during strict export

### Runtime Alias Folding vs Editable Export

Input may contain:

```tdb
PHASE A1_FCC %A 2 1 1 !
PHASE FCC_A1 %A 2 1 1 !
TYPE_DEF Y GES AMEND_PHASE_DESCRIPTION FCC_4SL DIS_PART A1_FCC !
```

Runtime import for calculation folds the helper alias to the canonical
disordered phase:

```tdb
PHASE FCC_A1 %A 2 1 1 !
TYPE_DEF D GES AMEND_PHASE_DESCRIPTION FCC_4SL DIS_PART FCC_A1 !
```

`preserve` export keeps the source helper alias when the source provided it:

```tdb
PHASE A1_FCC %A 2 1 1 !
PHASE FCC_A1 %A 2 1 1 !
TYPE_DEF D GES AMEND_PHASE_DESCRIPTION FCC_4SL DIS_PART A1_FCC !
```

Similarly at runtime:

```text
A2_BCC -> BCC_A2
```

## Developer Checklist

When changing the TDB parser/writer/runtime path, check the following:

1. Comments: commented-out records must not enter IR.
2. Ranges: multi-range expressions must not duplicate lower bounds after the
   first range.
3. Numeric exponents: convert `1.D-4` to `1.0E-4` without corrupting function
   names.
4. Endmembers: concrete `G` and wildcard `G` can cover endmembers; comma-list
   `G` cannot.
5. All-vacancy: generate zero compatibility records with a warning.
6. Aliases: fold `A2_BCC` to `BCC_A2` and `A1_FCC` to `FCC_A1` for runtime
   calculation.  In writer `preserve` mode, keep or synthesize helper aliases;
   in `canonical` mode, remove helper aliases.
7. `:F`: expand equivalent interactions without coefficient scaling.
8. Ternary CEF: preserve lone `;0` as unweighted; mark complete `;0`/`;1`/`;2`
   sets for composition-dependent runtime evaluation.
9. Magnetic aliases: accept `BM`, `BMAG`, and `BMAGN`; write `BMAGN`.
10. External checks: compare selected strict exports with TC-Python when possible.
11. Runtime checks: compare not only stable phases but also G, H, S, Cp and
    important endmember site fractions for order-disorder phases.
