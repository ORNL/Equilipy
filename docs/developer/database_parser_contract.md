# Database parser contract

This page identifies the Python parser variables and Fortran parser variables
involved in Step 1, "Parse database", of the Equilipy workflow.

The immediate purpose is to define what a future parser, for example a `.tdb`
parser, must produce so the existing Fortran backend can run without knowing
which database format was used.

```{note}
This contract is based on `utils._pyvar2fvar()`, `variables.py`,
`parse_chemsage_header.py`, and `ModuleParseCS.f90`. It describes the current
runtime interface, not an ideal future API.
```

## Dimension symbols

| Symbol | Meaning | Current maximum or source |
|---|---|---|
| `E` | number of database elements, `nElementsCS` | practical max `119` if electron `e-` plus atomic numbers `1..118` are supported |
| `S` | total species count, `nSpeciesCS` | dynamic, limited by memory |
| `P` | solution-phase count, `nSolnPhasesSysCS` | `<= nSolnPhasesSysMax = 500` |
| `SPmax` | maximum species count in any solution phase, `nMaxSpeciesPhaseCS` | dynamic, `<= S` |
| `C` | number of sublattice/SRO phases, `nCountSublatticeCS` | `<= P` |
| `SLmax` | maximum sublattices per phase, `nMaxSublatticeCS` | `5` |
| `Gcoeff` | Gibbs coefficient row count, `nGibbsCoeff` | `13` |
| `Geqmax` | maximum Gibbs equations per species, `nMaxGibbsEqs` | `6` |
| `Nparam` | excess parameter count, `nParamCS` | current ChemSage parser preallocates `1000` rows |
| `Nmag` | magnetic parameter count, `nMagParamCS` | current ChemSage parser preallocates `1000` rows |
| `NparamMax` | max components in a mixing parameter, `nParamMax` | `4` |
| `ParamWidth` | integer parameter row width | `nParamMax * 2 + 3 = 11` |

Current supported solution phase type strings are:

```text
IDMX, QKTO, SUBL, RKMP, RKMPM, SUBLM, SUBG, SUBQ, SUBI, SUBM
```

## Required parser conventions

The parser output is stored in `equilipy.variables` and then copied into
`fort.moduleparsecs` by `utils._pyvar2fvar()`.

Important conventions:

- Python arrays are zero-based, but several stored IDs are one-based because
  the Fortran code uses one-based phase/species/constituent indexing.
- `nSpeciesPhaseCS`, `nParamPhaseCS`, and `nMagParamPhaseCS` are cumulative
  counts per solution phase. `_pyvar2fvar()` prepends a zero before sending
  them to Fortran.
- Pure condensed compounds have `iPhaseCS = 0`.
- Dummy species have `iPhaseCS = -1`.
- Solution species have `iPhaseCS = 1..P`.
- Unused padded entries should be zero or blank. Downstream code relies on the
  count arrays to decide which entries are active.
- Character arrays must be padded or convertible to the fixed lengths used by
  the Fortran string helpers.
- When `P = 0` (pure-compound-only database), the current parser sets
  `j = max(1, P)` and uses `j` for the first dimension of `nParamPhaseCS`,
  `nMagParamPhaseCS`, `iPhaseSublatticeCS`, `nPairsSROCS`,
  `cSolnPhaseNameCS`, and `cSolnPhaseTypeCS`. Other multi-dimensional arrays
  use `P` directly, so they become zero-length when `P = 0`. A new parser
  must preserve this split.
- In the same `P = 0` path, `nSpeciesPhaseCS` remains a length-1 cumulative
  sentinel with value `0`, because current parsing code indexes
  `nSpeciesPhaseCS[-1]` for pure species.

## Scalar and generated fields

| Python variable | Fortran target | Dimension | Maximum size | Required contents |
|---|---|---|---|---|
| `nElementsCS` | `moduleparsecs.nelementscs` | scalar integer | practical max `119` | Number of element entries in the database system. Includes electron as `e-` if present. |
| `nSpeciesCS` | `moduleparsecs.nspeciescs` | scalar integer | dynamic | Total species count: all solution-phase species plus pure condensed species. |
| `nSolnPhasesSysCS` | `moduleparsecs.nsolnphasessyscs` | scalar integer | `500` | Number of solution phases in the database after removing an empty gas phase, if applicable. |
| `nCountSublatticeCS` | `moduleparsecs.ncountsublatticecs` | scalar integer | `P` | Count of phases that use sublattice/SRO-style data structures. |
| `nMaxSpeciesPhaseCS` | `moduleparsecs.nmaxspeciesphasecs` | scalar integer | `S` | Maximum species or pair-record count among solution phases. Used to size many padded phase arrays. |
| `nParamCS` | sizes `moduleparsecs.iparampasscs` | scalar integer | current parser capacity `1000` | Total number of excess mixing-parameter rows that are active. `_pyvar2fvar()` does not copy it to `moduleparsecs.nparamcs`; Fortran can derive the total from the final `nParamPhaseCS` boundary. |
| `nMagParamCS` | sizes `moduleparsecs.imagparampasscs` | scalar integer | current parser capacity `1000` | Total number of magnetic mixing-parameter rows that are active. `_pyvar2fvar()` does not copy it to `moduleparsecs.nmagparamcs`; Fortran can derive the total from the final `nMagParamPhaseCS` boundary. |
| generated zeros | `moduleparsecs.iparampasscs` | `(max(Nparam, 1),)` | `Nparam` active rows | Fortran pass flags for excess parameters. `_pyvar2fvar()` initializes this to zeros; `CheckSystemExcess.f90` marks active parameters. |
| generated zeros | `moduleparsecs.imagparampasscs` | `(max(Nmag, 1),)` | `Nmag` active rows | Fortran pass flags for magnetic parameters. `_pyvar2fvar()` initializes this to zeros; `CheckSystemExcess.f90` marks active parameters. |

Parser constants that should match `ModuleParseCS.f90`:

| Python variable | Fortran constant | Value | Meaning |
|---|---|---|---|
| `nMaxSublatticeCS` | `nMaxSublatticeCS` | `5` | Maximum number of sublattices in a phase. |
| `nSolnPhasesSysMax` | `nSolnPhasesSysMax` | `500` | Maximum solution-phase count accepted by the parser. |
| `nGibbsCoeff` | `nGibbsCoeff` | `13` | Rows in the Gibbs coefficient matrix. |
| `nMaxGibbsEqs` | `nMaxGibbsEqs` | `6` | Maximum Gibbs temperature intervals per species. |
| `nParamMax` | `nParamMax` | `4` | Maximum order/components for mixing parameters. |

## One-dimensional arrays

| Python variable | Fortran target | Python dimension | Fortran dimension after transfer | Maximum size | Required contents |
|---|---|---|---|---|---|
| `iPhaseCS` | `moduleparsecs.iphasecs` | `(S,)` | `(S,)` | `S` | Phase membership for every species. `1..P` means solution phase, `0` means pure condensed compound, `-1` means dummy species. |
| `iParticlesPerMoleCS` | `moduleparsecs.iparticlespermolecs` | `(S,)` | `(S,)` | `S` | Particles per constituent formula mass. Defaults to `1`; used to scale species stoichiometry. |
| `nGibbsEqSpecies` | `moduleparsecs.ngibbseqspecies` | `(S,)` | `(S,)` | `S`; each value `<= Geqmax` | Number of Gibbs-energy temperature intervals for each species. |
| `nSublatticePhaseCS` | `moduleparsecs.nsublatticephasecs` | `(P,)` | `(P,)` | `P`; each value `<= SLmax` | Number of sublattices for each sublattice phase. Active entries are indexed through `iPhaseSublatticeCS`. |
| `iPhaseSublatticeCS` | `moduleparsecs.iphasesublatticecs` | `(max(1, P),)` | `(max(1, P),)` | `max(1, P)` | Map from solution phase to sublattice-phase counter. `0` for non-sublattice phases; length-1 zero sentinel when `P = 0`. |
| `dAtomicMass` | `moduleparsecs.datomicmasscs` | `(E,)` | `(E,)` | practical max `119` | Atomic mass for each database element in the same order as `cElementNameCS`. |
| `nParamPhaseCS` | `moduleparsecs.nparamphasecs` | `(max(1, P),)` cumulative | `(max(1, P) + 1,)` with prepended `0` | `max(1, P) + 1` | Cumulative excess-parameter count after each solution phase. Fortran uses differences between adjacent entries to find phase-local parameters. Length-1 zero sentinel when `P = 0`. |
| `nSpeciesPhaseCS` | `moduleparsecs.nspeciesphasecs` | normally `(P,)` cumulative; `(1,)` zero sentinel when `P = 0` | `(len(nSpeciesPhaseCS) + 1,)` with prepended `0` | `max(1, P) + 1` after transfer | Cumulative solution-species count after each solution phase. Pure condensed species start after the last entry. |
| `nMagParamPhaseCS` | `moduleparsecs.nmagparamphasecs` | `(max(1, P),)` cumulative | `(max(1, P) + 1,)` with prepended `0` | `max(1, P) + 1` | Cumulative magnetic-parameter count after each solution phase. Length-1 zero sentinel when `P = 0`. |
| `cRegularParamCS` | `moduleparsecs.cregularparamcs` | `(Nparam_capacity,)` char length `1` | `(Nparam_capacity,)` char length `1` | current parser capacity `1000` | Model-specific one-character flag for regular/SRO parameters. Blank for models that do not need a flag. |
| `cElementNameCS` | `str1d("cElementNameCS", ...)` | `(E,)` char length `3` | `(E,)` char length `3` | practical max `119` | Element symbols padded to 3 characters. Electron is normalized to `e- `. |
| `cSolnPhaseTypeCS` | `str1d("cSolnPhaseTypeCS", ...)` | `(max(1, P),)` char length `8` | `(max(1, P),)` char length `8` | `max(1, P)` with `P <= 500` | Solution model type string such as `IDMX`, `RKMP`, `SUBL`, `SUBG`, padded to 8 characters. Blank sentinel when `P = 0`. |
| `cSolnPhaseNameCS` | `str1d("cSolnPhaseNameCS", ...)` | `(max(1, P),)` char length `25` | `(max(1, P),)` char length `25` | `max(1, P)` with `P <= 500` | Solution phase names padded to 25 characters. Blank sentinel when `P = 0`. |
| `cSpeciesNameCS` | `str1d("cSpeciesNameCS", ...)` | `(S,)` char length currently padded to `30` | `ModuleParseCS` declares char length `25` | `S` | Species/endmember/compound names. The current Python parser pads to 30, but `ModuleParseCS` declares 25; verify before formalizing a new parser. |

## Two-dimensional arrays

| Python variable | Fortran target | Python dimension | Maximum size | Required contents |
|---|---|---|---|---|
| `nPairsSROCS` | `moduleparsecs.npairssrocs` | `(max(P, 1), 2)` | `(max(1, P), 2)` | SRO/SUBG/SUBQ pair-count metadata. Active rows correspond to sublattice/SRO phase counter. Column meanings are model-specific in the current parser; column 2 is used as a pair-record loop bound in several paths. |
| `nConstituentSublatticeCS` | `moduleparsecs.nconstituentsublatticecs` | `(P, SLmax)` | `(P, 5)` | Number of constituents on each sublattice for each sublattice phase. Active rows are sublattice/SRO phase rows. |
| `nSublatticeElementsCS` | `moduleparsecs.nsublatticeelementscs` | `(P, SLmax)` | `(P, 5)` | Number of element/constituent entries represented on each sublattice for SRO/SUBG/SUBQ handling. |
| `dGibbsMagneticCS` | `moduleparsecs.dgibbsmagneticcs` | `(S, 4)` | `(S, 4)` | Magnetic Gibbs-energy coefficients for species and solution phases. Nonmagnetic entries are zero. |
| `dStoichSublatticeCS` | `moduleparsecs.dstoichsublatticecs` | `(P, SLmax)` | `(P, 5)` | Stoichiometric coefficient of each sublattice in a sublattice phase. |
| `dZetaSpeciesCS` | `moduleparsecs.dzetaspeciescs` | `(P, SPmax)` | `(P, SPmax)` | Zeta/FNN-SNN ratio or related SRO species coefficient for SUBG/SUBQ-style models. |
| `dGibbsCoeffSpeciesTemp` | `moduleparsecs.dgibbscoeffspeciestemp` | `(Gcoeff, S * Geqmax)` | `(13, S * 6)` | Flattened Gibbs-energy coefficient table. Columns are **contiguously packed**: all intervals for species 0, then all intervals for species 1, etc., with no gaps between species that use fewer than `Geqmax` intervals. The array is allocated `(13, S * 6)` to guarantee space, but only `sum(nGibbsEqSpecies)` columns are active. Fortran locates each species' columns via cumulative sums of `nGibbsEqSpecies`. Row 0 (Python) / row 1 (Fortran) stores the upper temperature limit; rows 1–6 store the standard Gibbs coefficients; rows 7–12 store extra custom-power term pairs. |
| `dStoichSpeciesCS` | `moduleparsecs.dstoichspeciescs` | `(S, E)` | `(S, E)` | Species stoichiometry matrix in database element order. Current parser scales by `iParticlesPerMoleCS`. |
| `iMagneticParamCS` | `moduleparsecs.imagneticparamcs` | `(Nmag_capacity, ParamWidth)` | current parser `(1000, 11)` | Integer metadata for magnetic mixing parameters: phase-local constituent IDs, parameter order, and model-specific indexing. Active rows are `0:nMagParamCS-1`. |
| `dMagneticParamCS` | `moduleparsecs.dmagneticparamcs` | `(Nmag_capacity, 2)` | current parser `(1000, 2)` | Numeric coefficients associated with `iMagneticParamCS`. Active rows are `0:nMagParamCS-1`. |
| `iRegularParamCS` | `moduleparsecs.iregularparamcs` | `(Nparam_capacity, ParamWidth)` | current parser `(1000, 11)` | Integer metadata for excess Gibbs mixing parameters. First entry usually stores component/parameter type count; following entries store constituent IDs, exponents, order, or SRO-specific IDs. |
| `dRegularParamCS` | `moduleparsecs.dregularparamcs` | `(Nparam_capacity, 6)` | current parser `(1000, 6)` | Numeric coefficients for each excess Gibbs mixing parameter row. Active rows are `0:nParamCS-1`. |
| `cPairNameCS` | `str2d("cPairNameCS", ...)` | `(P, SPmax, 30)` byte chars | `(P, SPmax)` char length `30` | Human-readable SRO pair names for SUBG/SUBQ-style phases. Unused entries blank. |

## Three-dimensional arrays

| Python variable | Fortran target | Python dimension | Maximum size | Required contents |
|---|---|---|---|---|
| `iConstituentSublatticeCS` | `moduleparsecs.iconstituentsublatticecs` | `(P, SLmax, SPmax)` | `(P, 5, SPmax)` | Constituent/species index map for each sublattice. Values are one-based IDs used by Fortran model routines. |
| `iPairIDCS` | `moduleparsecs.ipairidcs` | `(P, SPmax, 4)` | `(P, SPmax, 4)` | Four integer IDs defining each SRO pair record. Active pair rows are bounded by `nPairsSROCS`. |
| `iChemicalGroupCS` | `moduleparsecs.ichemicalgroupcs` | `(P, SLmax, SPmax)` | `(P, 5, SPmax)` | Chemical group identifiers for constituents in sublattice/SRO models. |
| `dSublatticeChargeCS` | `moduleparsecs.dsublatticechargecs` | `(P, SLmax, SPmax)` | `(P, 5, SPmax)` | Charge assigned to each constituent on each sublattice. Zero when not used. |
| `dStoichPairsCS` | `moduleparsecs.dstoichpairscs` | `(P, SPmax, E)` | `(P, SPmax, E)` | Stoichiometry of each SRO pair record in database element order. |
| `dConstituentCoefficientsCS` | `moduleparsecs.dconstituentcoefficientscs` | `(P, SPmax, 5)` | `(P, SPmax, 5)` | Pair/constituent coefficients used by SUBG/SUBQ-style models. |
| `dCoordinationNumberCS` | `moduleparsecs.dcoordinationnumbercs` | `(P, SPmax, 4)` | `(P, SPmax, 4)` | Coordination numbers associated with SRO pair records. |
| `cConstituentNameSUBCS` | `str3d("cConstituentNameSUBCS", ...)` | `(P, SLmax, SPmax, 8)` byte chars | `(P, SLmax, SPmax)` char length `8` | Constituent names for each sublattice entry. Unused entries blank. |

## Python-only database fields

The parser dictionary returned by `variables.to_dict()` contains additional
fields that are not copied by `_pyvar2fvar()` but are still used by Python
preprocessing and post-processing.

| Python variable | Dimension | Required contents |
|---|---|---|
| `cPhaseNames` | list length `P + nPureSpeciesCS` | Public phase-name list used by phase selection and user-facing phase lookup. |
| `nPureSpeciesCS` | scalar | Number of pure condensed species. Equals `S - nSpeciesPhaseCS[-1]`. |
| `nSolnPhaseCS` | `(P,)` | Per-phase species counts (non-cumulative), saved before `nSpeciesPhaseCS` is converted to cumulative. Used by Python post-processing to reconstruct individual phase sizes. |
| `cEndmemberNameCS` | `(S,)` char length `30` | Endmember names padded to 30 characters. Used by Python result post-processing. |
| `INFO` | scalar | Parser status code. |
| `DataBase` | list | Remaining parser token list after parsing. Usually empty or diagnostic only after a successful parse. |
| `iCounterGibbsEqn` | scalar | Number of Gibbs equation columns populated in `dGibbsCoeffSpeciesTemp`. Equals `sum(nGibbsEqSpecies)` after a successful parse. |
| `cPeriodicTable` | dict | Element symbol to atomic-number/mass map used by Python preprocessing. |
| `indx` | scalar | Parser-internal 1-based species index tracking current position during parsing. Exported for `load_database()` roundtrip compatibility. A new parser should include it, even though `_pyvar2fvar()` does not use it directly. |

## Minimum checklist for a new parser

Before connecting a new parser to `load_database()` and `_pyvar2fvar()`, verify:

1. Count fields are internally consistent:
   `S = nPureSpeciesCS + nSpeciesPhaseCS[-1]`.
2. Cumulative arrays are truly cumulative:
   `nSpeciesPhaseCS`, `nParamPhaseCS`, and `nMagParamPhaseCS`.
3. Gibbs columns are flattened in the exact order expected by
   `nGibbsEqSpecies`.
4. Species, phase, constituent, and pair IDs use the current one-based Fortran
   conventions where required.
5. Character fields are padded to the expected lengths.
6. Unused padded rows are zero or blank.
7. The parser dictionary includes both f2py-transfer fields and Python-only
   post-processing fields.
8. SUBL/SUBG/SUBQ/SRO arrays are populated only when the corresponding model
   needs them, but their shapes are still valid.
9. A pure-compound-only database is tested separately because several current
   arrays use `max(P, 1)` while others may be zero-length.
