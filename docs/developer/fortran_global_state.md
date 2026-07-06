# Fortran global-state inventory

This page records the current Fortran module state exposed through
`equilipy.equilifort`. It is a developer contract for memory ownership,
state isolation, and solver/result invariants.

```{note}
This inventory is based on static source inspection. A variable listed as
"weak/no observed use" is not automatically dead code. f2py access, generated
bindings, conditional Fortran paths, and old parser paths can hide real use.
Treat those entries as review targets, not deletion candidates.
```

## Scope

The calculation path uses Fortran modules as process-global mutable state:

```text
Python input
-> f2py writes into ModuleParseCS and ModuleThermoIO
-> Fortran allocates ModuleThermo, ModuleGEMSolver, and ModuleSubMin state
-> Fortran postprocess writes ModuleThermoIO output arrays
-> Python copies selected values into result objects
-> reset routines deallocate global arrays
```

The main modules are:

| Module | File | Role | Main reset owner |
|---|---|---|---|
| `ModuleParseCS` | `fsrc/1_Modules/ModuleParseCS.f90` | Parsed ChemSage database arrays | `ResetThermoParser` |
| `ModuleThermoIO` | `fsrc/1_Modules/ModuleThermoIO.f90` | User input, error status, public output arrays | `ResetThermo` for output arrays |
| `ModuleThermo` | `fsrc/1_Modules/ModuleThermo.f90` | Internal thermodynamic system and phase state | `ResetThermo` |
| `ModuleGEMSolver` | `fsrc/1_Modules/ModuleGEMSolver.f90` | Global equilibrium minimizer state | `ResetThermo`, plus local initializer cleanup |
| `ModuleSubMin` | `fsrc/1_Modules/ModuleSubMin.f90` | Solution-phase subminimization state | Partial `ResetThermo`, plus `SubMinInit` cleanup |

## Database-To-GEM State Bridge

The current database/system-setting path is DAT/ChemSage shaped. Python parses a
database into `equilipy.variables`, then `_pyvar2fvar()` copies that parser state
into `ModuleParseCS`. The Fortran calculation does not consume a high-level
database object; it consumes the global parser arrays plus the global input
condition.

```text
Python database dict
-> equilipy.variables parser arrays
-> _pyvar2fvar()
-> ModuleParseCS
-> CheckThermoInput()
-> InitThermo()
-> CheckSystem()
-> CheckSystemExcess()
-> CompThermoData()
-> CheckThermoData()
-> GEMSolver()
-> PostProcess()
```

Important ownership points:

- `load_database()` rewrites Python parser globals only. It does not clear or
  rewrite Fortran arrays until `_pyvar2fvar()` is called.
- `system_check()` is a Python system-filtering helper. It builds
  `iElementSys`, `iElementSysIndex`, `iElementDBIndex`, `iSys2DBSoln`,
  `iSys2DBComp`, `iSys2DBSpecies`, `cPhaseNameSys`, and `PhaseNameSys`.
- `phase_selection()` narrows those Python maps and writes the Fortran masks
  `ModuleParseCS.iPhaseCS` and `ModuleThermo.iSolnPS`. It also mutates
  `var.iPhaseCS` into a selected-phase mask, so reload the database before
  treating `var.iPhaseCS` as the raw parser phase-ownership vector.
- `CheckSystem()` converts input amounts to moles, normalizes element amounts,
  filters species/phases for the active system, allocates active
  `ModuleThermo` arrays, and writes `nSpeciesPhase`, `iSpeciesPass`,
  `cElementName`, `dMolesElement`, `dAtomicMass`, `cSolnPhaseName`, and
  `cSolnPhaseType`.
- `CheckSystemExcess()` is not only validation. It sets active regular and
  magnetic parameter pass masks and fills active sublattice/MQM phase metadata.
- `CompThermoData()` compacts database species into active system species,
  evaluates standard properties, maps excess and magnetic parameters, and fills
  `dStoichSpecies`, `dAtomFractionSpecies`, `dStdGibbsEnergy`,
  `dStdEnthalpy`, `dStdEntropy`, `dStdHeatCapacity`, `dChemicalPotential`,
  `dPartialEnthalpy`, `dPartialEntropy`, and `dPartialHeatCapacity`.
- `PostProcess()` computes system properties, mass-basis output arrays, and
  species activities. Python result capture must copy both Fortran arrays and
  Python maps before either layer is reset or reused.

For Phase 6 planning, the working source map lives in
`_ongoing_development/source_var.md`. If that map is promoted into public
developer documentation, keep this section as the short orientation and avoid
duplicating every variable table here.

## Correctness-sensitive subroutine expectations

These routines are central to the PEA/GEM setup, minimization, postprocess, and
reset lifecycle. The descriptions below are the expected current behavior, not
change history.

| Routine | State read | State written | Expected invariant after return |
|---|---|---|---|
| `InitGEMSolverNew.f90` | current GEM assemblage, species/phase dimensions, old chemical potentials, standard Gibbs state | `nSpeciesLevel`, expanded `dAtomFractionSpecies`, `dStoichSpeciesLevel`, `dChemicalPotential`, `dPhasePotential`, `iPhaseLevel`, PEA/GEM history arrays | The expanded species-level arrays have deterministic defaults before any solution minimum point is written: atom/stoichiometry arrays are zeroed, chemical/phase potentials are set to `5D9`, and `iPhaseLevel` is zero. Unfilled pseudo-solution slots must therefore be energetically unselectable and must not carry stale memory. |
| `CompInitMinSolnPoint.f90` | `nElements`, `nSolnPhasesSys`, `nSpeciesPhase`, solution stoichiometry, `dStdGibbsEnergy`, `dElementPotential`, miscibility flags | `dMolFraction`, `dAtomFractionSpecies`, `dStoichSpeciesLevel`, `dChemicalPotential`, `iPhaseLevel`, `dGEMFunctionNorm` | Each solution phase receives a bounded set of candidate minimum points. The number of endpoint trials is controlled by local `nElementOrConstituent = MIN(nElements, nConstituents)`, so temporary arrays and loops never index more endmember extrema than exist and never create more independent minimum-point candidates than the element/constituent dimension permits. Immiscible-copy propagation must prove `k <= nSolnPhasesSys` before reading `lMiscibility(k)`, `cSolnPhaseName(k)`, or `nSpeciesPhase(k)`, because Fortran `.AND.` evaluation is not guaranteed to short-circuit. |
| `CompDrivingForceAll.f90` | `nElements`, `nSolnPhasesSys`, `nSpeciesPhase`, `dStdGibbsEnergy`, `dElementPotential`, miscibility flags, current solution mole fractions | `dPhasePotential`, `dDrivingForceSoln`, `dChemicalPotential`, `dMolFraction` | Stoichiometric and solution-phase driving forces are refreshed for the active system. Immiscible-copy propagation must prove `k <= nSolnPhasesSys` before any `k`-indexed array access, and only same-name miscibility copies may receive alternate local minima. |
| `GetNewAssemblage.f90` | `dPhasePotential`, current `iAssemblage`, current GEM stoichiometry/potential/mole-fraction arrays, `dMolesElement` | candidate `iAssemblage`, `dMolesPhase`, `dChemicalPotentialGEM`, `dStoichSpeciesGEM`, `dAtomFractionSpeciesGEM`, `dMolFractionGEM`, `iPhaseGEM` | When no strictly nonnegative candidate is found, the fallback replacement index is initialized from `MAXLOC(dMinValuedMolesPhase)` before tolerance checks. The fallback must never use a stale loop index to replace a phase. |
| `PostProcessPEA.f90` | PEA/GEM assemblage history, `iPhaseGEM`, `dMolesPhaseLast`, `dMolFractionGEM`, `lMiscibility` | final `iAssemblage`, `dMolesPhase`, `dMolesSpecies`, `dMolFraction`, `nConPhases`, `nSolnPhases`, `lSolnPhases`, `dGEMFunctionNormLast` | Final assemblage slots remain within `nElements`. PEA minimum-point pseudo species with `iAssemblageLast(i) > nSpecies` are mapped back to their owning solution phase when `iPhaseGEM(i)` is unavailable, instead of being treated as pure compounds. |
| `GEMNewton.f90` | active assemblage, species stoichiometry, `dMolesSpecies`, `dChemicalPotential`, `iParticlesPerMole` | Newton matrix/vector work arrays, `dUpdateVar`, local `INFOLocal` | Solution-phase element-balance Hessian and RHS terms are scaled per species by `iParticlesPerMole(m:n)`. This is required when constituents in the same phase do not share one particles-per-mole value. |
| `ResetThermo.f90` | allocated state from `ModuleThermo`, `ModuleThermoIO`, `ModuleGEMSolver`, `ModuleSubMin` | deallocates active arrays; resets `INFOThermo`, `dElementMass`, GEM iteration counters, swap/index scalars, driving-force scalars, convergence flags, and postprocess flags | After `resetthermo()`, transient thermo/GEM scalar status must not leak into the next calculation in the same process. Parser arrays are not reset here; use `resetthermoall()` when that layer must also be cleared. |
| `ResetThermoParser.f90` | allocated parser/database arrays from `ModuleParseCS` | deallocates parser arrays; resets `ModuleParseCS.INFO` and `ModuleThermoIO.INFOThermo` | Parser reset clears parser-side failure status as well as arrays. Thermodynamic, GEM, and reinit state are outside this routine's ownership. |

## Python access map

Python currently touches a small subset of the Fortran globals directly.

| Python file | Fortran module fields |
|---|---|
| `load_database.py` | calls `utils._pyvar2fvar()` through `load_fortran_database()` when the cached Fortran payload is stale |
| `list_phases.py` | loads Python database globals and calls `system_check()` for the selected elements |
| `utils.py` | writes most `ModuleParseCS` database arrays through `_pyvar2fvar()` and resets `moduleparsecs.info` / `modulethermoio.infothermo` |
| `system_check.py` | writes `moduleparsecs.lendmembers2species` |
| `input_condition.py` | writes `modulethermoio` temperature, pressure, units, element masses, output mode |
| `phase_selection.py` | writes `moduleparsecs.iphasecs` and `modulethermo.isolnps` |
| `minimize.py` | calls the Fortran workflow and reads `modulethermoio.infothermo` and `modulegemsolver.dgemfunctionnorm` |
| `results/equilib.py` and `results/scheil.py` | read `ModuleThermo` phase/species arrays and `ModuleThermoIO` system/output arrays; `post_process.py` is now a compatibility facade |
| `find_transition.py` | reads `modulethermo.iassemblage` |
| `scheil_cooling.py` | reads `modulethermo.iassemblage` |

Current reset calls from Python:

| Python path | Reset call |
|---|---|
| `equilib_single()` | `fort.resetthermoall()` |
| `_equilib_batch()` row loop | `fort.resetthermo()` then `fort.resetthermoparser()` |
| `ScheilResult` initialization and append paths | `fort.resetthermo()` |

`resetthermoall()` calls `ResetThermo` and `ResetThermoParser`.

## ModuleParseCS

`ModuleParseCS` stores parsed ChemSage database state. In the current Python
workflow, most of this state is first parsed into `equilipy.variables`, then
copied into `fort.moduleparsecs` by `utils._pyvar2fvar()`.

### Variables

| Group | Variables |
|---|---|
| Scalar counts/status | `nElementsCS`, `nSpeciesCS`, `nSolnPhasesSysCS`, `INFO`, `iMiscSUBI`, `nParamCS`, `nCountSublatticeCS`, `nMaxSpeciesPhaseCS`, `nMagParamCS` |
| Integer parameters | `nSolnPhasesSysMax`, `nMaxSublatticeCS`, `nSolnTypeSupport`, `nGibbsCoeff`, `nMaxGibbsEqs`, `nParamMax` |
| Integer arrays | `nSpeciesPhaseCS`, `nGibbsEqSpecies`, `iPhaseCS`, `iParticlesPerMoleCS`, `nParamPhaseCS`, `iParamPassCS`, `nSublatticePhaseCS`, `iPhaseSublatticeCS`, `iMagParamPassCS`, `nMagParamPhaseCS`, `iSUBIMixTypeCS`, `nInterpolationOverrideCS`, `iRegularParamCS`, `nConstituentSublatticeCS`, `nPairsSROCS`, `iMagneticParamCS`, `iSUBIParamDataCS`, `nSublatticeElementsCS`, `iInterpolationOverrideCS`, `iConstituentSublatticeCS`, `iPairIDCS`, `iChemicalGroupCS` |
| Real arrays | `dAtomicMassCS`, `dGibbsCoeffSpeciesTemp`, `dRegularParamCS`, `dGibbsMagneticCS`, `dMagneticParamCS`, `dStoichSublatticeCS`, `dStoichSpeciesCS`, `dZetaSpeciesCS`, `dStoichConstituentCS`, `dQKTOParamsCS`, `dSublatticeChargeCS`, `dStoichPairsCS`, `dConstituentCoefficientsCS`, `dCoordinationNumberCS` |
| Logical | `lEndmembers2Species` |
| Character arrays | `cElementNameCS`, `cSolnPhaseTypeCS`, `cSolnPhaseNameCS`, `cSpeciesNameCS`, `cPairNameCS`, `cConstituentNameSUBCS`, `cRegularParamCS` |
| Character parameters | `cSolnPhaseTypeSupport` |

### Observed use

Active bridge variables written by Python:

- Counts/status: `INFO`, `nElementsCS`, `nSpeciesCS`, `nSolnPhasesSysCS`,
  `nCountSublatticeCS`, `nMaxSpeciesPhaseCS`
- Phase/species arrays: `iPhaseCS`, `iParticlesPerMoleCS`,
  `nGibbsEqSpecies`, `nSpeciesPhaseCS`, `nParamPhaseCS`,
  `nMagParamPhaseCS`
- Sublattice arrays: `nSublatticePhaseCS`, `iPhaseSublatticeCS`,
  `nPairsSROCS`, `nConstituentSublatticeCS`, `nSublatticeElementsCS`,
  `iConstituentSublatticeCS`, `iPairIDCS`, `iChemicalGroupCS`
- Thermodynamic arrays: `dAtomicMassCS`, `dGibbsMagneticCS`,
  `dStoichSublatticeCS`, `dZetaSpeciesCS`, `dGibbsCoeffSpeciesTemp`,
  `dStoichSpeciesCS`, `iMagneticParamCS`, `dMagneticParamCS`,
  `iRegularParamCS`, `dRegularParamCS`, `dSublatticeChargeCS`,
  `dStoichPairsCS`, `dConstituentCoefficientsCS`,
  `dCoordinationNumberCS`
- Character arrays: `cRegularParamCS`, `cElementNameCS`,
  `cSolnPhaseTypeCS`, `cSolnPhaseNameCS`, `cSpeciesNameCS`,
  `cPairNameCS`, `cConstituentNameSUBCS`

Weak/no observed use in the current static pass:

- `iMiscSUBI`
- `iSUBIMixTypeCS`
- `nInterpolationOverrideCS`
- `iSUBIParamDataCS`
- `iInterpolationOverrideCS`
- `dQKTOParamsCS`
- `dStoichConstituentCS` from Fortran source inspection, although Python
  parser state still contains related data

### Allocation and deallocation

In the current Python-driven workflow, most `ModuleParseCS` allocatable arrays
are allocated by f2py assignment from NumPy arrays in `utils._pyvar2fvar()`. The
original Fortran parser path also has parser-owned allocation logic.

`ResetThermoParser` deallocates the main parser groups:

- Core integer/species group:
  `nSpeciesPhaseCS`, `nGibbsEqSpecies`, `iPhaseCS`,
  `iParticlesPerMoleCS`, `nParamPhaseCS`, `iParamPassCS`,
  `dStoichSpeciesCS`, `iRegularParamCS`, `cRegularParamCS`,
  `nMagParamPhaseCS`, `iMagneticParamCS`, `dMagneticParamCS`,
  `iMagParamPassCS`
- Core real group:
  `dAtomicMassCS`, `dGibbsCoeffSpeciesTemp`, `dRegularParamCS`,
  `dGibbsMagneticCS`
- Character group:
  `cElementNameCS`, `cSolnPhaseTypeCS`, `cSolnPhaseNameCS`,
  `cSpeciesNameCS`
- Sublattice group:
  `cConstituentNameSUBCS`, `iPhaseSublatticeCS`,
  `nConstituentSublatticeCS`, `iConstituentSublatticeCS`,
  `dStoichSublatticeCS`, `nSublatticePhaseCS`,
  `nSublatticeElementsCS`, `dZetaSpeciesCS`,
  `dConstituentCoefficientsCS`, `dSublatticeChargeCS`,
  `dStoichPairsCS`, `iChemicalGroupCS`, `cPairNameCS`
- SRO pair group:
  `nPairsSROCS`, `iPairIDCS`, `dCoordinationNumberCS`

Review note:

- `ResetThermoParser` has a second `dGibbsMagneticCS` deallocation guard after
  the core real group. This is usually harmless because the second guard should
  be false after the first deallocation.
- `ResetThermoParser` also resets `ModuleThermoIO.INFOThermo` so parser reset
  does not leave stale Fortran error status visible to the next parser/system
  setup call.
- Some declared parser arrays have no observed allocation/deallocation in the
  current path. They should be reviewed before removing because they may belong
  to older or not-yet-exercised model types.

## ModuleThermoIO

`ModuleThermoIO` is the user input and public Fortran output boundary.

### Variables

| Group | Variables |
|---|---|
| Input scalars | `iPrintResultsMode`, `dTemperature`, `dPressure`, `dTemperatureDiff` |
| Input arrays | `dElementMass` |
| Input strings | `cInputUnitTemperature`, `cInputUnitPressure`, `cInputUnitMass` |
| Output status/scalars | `INFOThermo`, `dGibbsEnergySys`, `dEnthalpySys`, `dEntropySys`, `dHeatCapacitySys` |
| Postprocess diagnostics | `dPostProcessGEMNormAtT`, `dPostProcessGEMNormPerturbed`, `iPostProcessIterAtT`, `iPostProcessIterPerturbed`, `lPostProcessAssemblageChanged` |
| Output real arrays | `dPair`, `dGramFraction`, `dGramSpecies`, `dGramPhase`, `dGramElement` |
| Output character arrays | `cPair` |

### Observed use

Active Python input fields:

- `dTemperatureDiff`
- `cInputUnitTemperature`
- `cInputUnitPressure`
- `cInputUnitMass`
- `dTemperature`
- `dPressure`
- `dElementMass`
- `iPrintResultsMode`

Active Python output fields:

- `INFOThermo`
- `dGibbsEnergySys`
- `dEnthalpySys`
- `dEntropySys`
- `dHeatCapacitySys`
- `dGramPhase`
- `dGramFraction`
- `dGramElement`
- `dTemperature`
- `dPressure`

Removed in the 06/23/2026 cleanup after a static active-use audit:

- Dormant reinit flags: `lReinitAvailable`, `lReinitLoaded`,
  `lReinitRequested`
- Unused legacy input/output fields: `iCounter`, `lPreset`,
  `cThermoFileName`, compound input buffers, legacy phase/species output
  buffers, and `lSpeciesStable`

### Allocation and deallocation

- `dGramFraction`, `dGramSpecies`, `dGramPhase`, and `dGramElement` are
  allocated in `PostProcess.f90` after deallocating previous contents.
- Those gram arrays are also deallocated by `ResetThermo`.
- `ModuleThermoIO.dPair` and `ModuleThermoIO.cPair` are allocated in
  `CalculateCompositionSUBG.f90` after deallocating previous contents.
- Nonallocatable input fields are overwritten by each Python call to
  `input_condition()`. `ResetThermo` also sets `dElementMass = 0.0d0`.

Review note:

- `ScheilResult` resets only `ResetThermo`, not `ResetThermoParser`. That is
  intentional for iterative reuse of the selected database/parser state, but it
  means parser state remains live through a Scheil sequence.
- `ModuleThermoIO.dPair` and `ModuleThermoIO.cPair` are not currently
  deallocated by `ResetThermo`. They are cleaned up when
  `CalculateCompositionSUBG.f90` is reached again, but global reset coverage
  should include them.
- These output pair arrays are distinct from `ModuleThermo.cPairName`, which is
  an internal pair-name array and is deallocated by `ResetThermo`.

## ModuleThermo

`ModuleThermo` is the largest internal state module. It stores selected-system
dimensions, species/phase maps, standard Gibbs-property arrays, phase
assemblage state, Lagrangian state, and solution-model arrays.

### Variables

| Group | Variables |
|---|---|
| Scalar counts | `nElements`, `nSpecies`, `nParam`, `nMaxParam`, `nDummySpecies`, `nElemOrComp`, `nMagParam`, `nSpeciesLevel`, `nConPhases`, `nSolnPhases`, `nSolnPhasesSys`, `nChargedConstraints`, `nElementsSys`, `nMaxSublatticeSys`, `nMaxConstituentSys`, `nCountSublattice` |
| Parameters | `iTolNum`, `nElementsPT`, `nMaxCompounds` |
| Integer arrays | `iPhase`, `nSpeciesPhase`, `iParticlesPerMole`, `iPhaseLevel`, `iCandidate`, `iSolnPS`, `iAssemblage`, `nParamPhase`, `iElementSystem`, `iSpeciesPass`, `nMagParamPhase`, `nSublatticePhase`, `iPhaseSublattice`, `iPhaseElectronID`, `iPhaseGEM`, `iPhaseGEMinit`, `iShuffled` |
| Integer matrix arrays | `iRegularParam`, `iterHistoryLevel`, `nConstituentSublattice`, `nPairsSRO`, `iMagneticParam`, `iSUBLParamData`, `nSublatticeElements` |
| Integer 3D arrays | `iConstituentPass`, `iConstituentSublattice`, `iPairID`, `iChemicalGroup` |
| Real scalars | `dIdealConstant`, `dNormalizeSum`, `dNormalizeInput`, `dMassScale`, `dMinPhasePotential`, `dToleranceLevel`, `dMaxElementPotential`, `dPEATol`, `dMinMolesPhase` |
| Real static array | `dTolerance` |
| Real arrays | `dGibbsSolnPhase`, `dMolesSpecies`, `dMagGibbsEnergy`, `dMagEnthalpy`, `dMagEntropy`, `dMagHeatCapacity`, `dStdGibbsEnergy`, `dStdEnthalpy`, `dStdEntropy`, `dStdHeatCapacity`, `dChemicalPotential`, `dExcessGibbsParam`, `dSpeciesTotalAtoms`, `dPhasePotential`, `dExcessHParam`, `dExcessSParam`, `dExcessCpParam`, `dPartialEnthalpy`, `dPartialEntropy`, `dPartialHeatCapacity`, `dPartialEnthalpyXS`, `dPartialEntropyXS`, `dPartialHeatCapacityXS`, `dElementPotential`, `dMolesPhase`, `dMolesElement`, `dMolFraction`, `dAtomicMass`, `dChemicalPotentialOld`, `dChemicalPotentialGEM`, `dMolFractionOld`, `dChemicalPotentialGEMinit` |
| Real matrix arrays | `dAtomFractionSpecies`, `dStoichSublattice`, `dStoichSpecies`, `dCoeffGibbsMagnetic`, `dZetaSpecies`, `dMagneticParam`, `dAtomFractionSpeciesOld`, `dStoichSpeciesLevel`, `dStoichSpeciesGEM`, `dAtomFractionSpeciesGEM`, `dMolesPhaseHistory`, `dEffStoichSolnPhase`, `dMolFractionGEM`, `dStoichSpeciesGEMinit`, `dAtomFractionSpeciesGEMinit`, `dMolFractionGEMinit` |
| Real 3D arrays | `dSiteFraction`, `dCoordinationNumber`, `dSublatticeCharge`, `dStoichPairs`, `dConstituentCoefficients` |
| Character arrays | `cElementName`, `cSpeciesName`, `cSolnPhaseType`, `cSolnPhaseName`, `cConstituentNameSUB`, `cRegularParam`, `cPairName` |
| Logical scalars | `lpseudo`, `lSkipLagrange` |

### Observed use

Important Python-facing result fields:

- `iPhase`
- `nSpeciesPhase`
- `iAssemblage`
- `dMolesPhase`
- `dMolFraction`
- `dMolesElement`

Important Fortran workflow fields:

- System dimensions and maps are prepared in `CheckSystem.f90`.
- Gibbs-property arrays are populated in `CompThermoData.f90`.
- Phase assemblage and Lagrangian state are initialized in
  `InitGEMSolver.f90` or `InitGEMSolverNew.f90`.
- `iShuffled` is allocated in `ShuffleAssemblage.f90`.
- Post-processing reads phase/species/mole arrays and writes `ModuleThermoIO`
  output arrays.

Weak/no observed use in the current static pass:

- `nElementsSys`
- `nMaxCompounds`
- `dMassScale`
- `dMinMolesPhase`
- `lpseudo`

These should be validated with targeted searches and runtime instrumentation
before removal.

### Allocation and deallocation

Major allocation owners:

| Routine | Main allocations |
|---|---|
| `CheckSystem.f90` | selected-system arrays, species/phase arrays, standard property arrays, sublattice arrays |
| `CheckSystemExcess.f90` | excess Gibbs/H/S/Cp parameter arrays and magnetic parameter arrays |
| `CompThermoData.f90` | Gibbs-property values and temporary/reallocated stoichiometry arrays |
| `InitGEMSolver.f90` | `dMolesPhase`, `dMolesSpecies`, `dElementPotential`, `iAssemblage`, GEM/history/update arrays |
| `InitGEMSolverNew.f90` | newer PEA/GEM initialization arrays; deterministic defaults for expanded pseudo-solution species-level slots |
| `CompInitMinSolnPoint.f90` | minimum-point candidate work arrays for solution-phase leveling; local `nElementOrConstituent` bounds candidate count and temporary dimensions |
| `GetNewAssemblage.f90` | trial assemblage replacement work arrays; fallback replacement index must be initialized before use |
| `ShuffleAssemblage.f90` | `iShuffled` |
| `PostProcess.f90` and `PostProcessPEA.f90` | post-processing working arrays and result-related state |

`ResetThermo` deallocates most active `ModuleThermo` arrays, including:

- species/phase maps and names
- standard property arrays
- chemical potentials
- moles and mole fractions
- sublattice, constituent, pair, and magnetic arrays
- GEM/PEA arrays
- postprocess gram arrays from `ModuleThermoIO`

`PostProcessPEA` must translate both ordinary stable species and PEA
minimum-point pseudo species into the final phase assemblage. For solution
minimum points, a positive assemblage value above `nSpecies` identifies a
minimum point, not a stoichiometric compound; it should be mapped back to the
owning solution phase through `iPhaseGEM` or by subtracting `nSpecies`.

Review note:

- Many initializer routines deallocate before allocating. This protects normal
  repeated calls if the initializer is reached.
- Failure paths remain important. If an exception exits Python after Fortran has
  allocated arrays but before the expected reset call, process-global state can
  survive into the next calculation in the same worker process.

## ModuleGEMSolver

`ModuleGEMSolver` stores global equilibrium-minimizer control state.

### Variables

| Group | Variables |
|---|---|
| Iteration counters | `iterLast`, `iterStep`, `iterRevert`, `iterGlobal`, `iterPEA`, `iterLG`, `iterLastCon`, `iterLastSoln`, `iterSwap`, `iterLastMiscGapCheck` |
| Last/swap indexes | `iConPhaseLast`, `iSolnPhaseLast`, `iSolnSwap`, `iPureConSwap` |
| Driving-force indexes | `iMinDrivingForceStoich`, `iMinDrivingForceSoln`, `iSpeciesRemove` |
| Parameter | `iterGlobalMax` |
| Integer arrays | `iAssemblageGEM`, `iAssemblageGEMinit`, `iterHistory` |
| Real scalars | `dGEMFunctionNorm`, `dGEMFunctionNormLast`, `dMaxSpeciesChange`, `dMinGibbs`, `dMinDrivingForceStoich`, `dMinDrivingForceSoln`, `dSpeciesRemove`, `dPlateau`, `xT`, `dGParam`, `dHParam`, `dSParam`, `dCpParam`, `dMaxPotentialTol` |
| Real arrays | `dSumMolFractionSoln`, `dMolesPhaseLast`, `dUpdateVar`, `dDrivingForceSoln`, `dPartialExcessGibbs`, `dPartialExcessGibbsLast`, `dUpdateVarLast`, `dMolesSpeciesLast`, `dElementPotentialLast`, `dPartialGParam`, `dPartialHParam`, `dPartialSParam`, `dPartialCpParam` |
| Logical scalars | `lDebugMode`, `lRevertSystem`, `lConverged`, `lSubConverged`, `lNegativeMolesPhase`, `lGibbsMinCheck`, `lPhaseChange`, `lCompbdOnly`, `lPostProcess` |
| Logical arrays | `lSolnPhases`, `lMiscibility`, `lPhaseChangeHistory` |

### Observed use

Active fields include:

- `iterGlobal`, `iterPEA`, and convergence-related flags in the solver path
- `dGEMFunctionNorm`, read by Python after minimization
- update/history arrays in the Lagrangian workflow
- excess partial arrays during Gibbs energy and chemical-potential calculations
- `iAssemblageGEM` and `dMolesPhaseLast` in leveling/postprocess paths

`ResetThermo` resets the scalar solver-status fields that can otherwise survive
worker-process reuse:

- iteration counters: `iterLast`, `iterStep`, `iterRevert`, `iterGlobal`,
  `iterPEA`, `iterLG`, `iterLastCon`, `iterLastSoln`, `iterSwap`,
  `iterLastMiscGapCheck`
- swap/index state: `iConPhaseLast`, `iSolnPhaseLast`, `iSolnSwap`,
  `iPureConSwap`, `iMinDrivingForceStoich`, `iMinDrivingForceSoln`,
  `iSpeciesRemove`
- convergence/driving-force scalars: `dGEMFunctionNorm`,
  `dGEMFunctionNormLast`, `dMaxSpeciesChange`, `dMinGibbs`,
  `dMinDrivingForceStoich`, `dMinDrivingForceSoln`, `dSpeciesRemove`,
  `dPlateau`, `dMaxPotentialTol`
- thermodynamic parameter scratch scalars: `xT`, `dGParam`, `dHParam`,
  `dSParam`, `dCpParam`
- logical flags: `lDebugMode`, `lRevertSystem`, `lConverged`,
  `lSubConverged`, `lNegativeMolesPhase`, `lGibbsMinCheck`, `lPhaseChange`,
  `lCompbdOnly`, `lPostProcess`

Weak/no observed use in the current static pass:

- `iterStep`
- `iterRevert`
- `iterLG`
- `iterSwap`
- `iSolnSwap`
- `iPureConSwap`
- `iMinDrivingForceStoich`
- `iMinDrivingForceSoln`
- `iSpeciesRemove`
- `dMaxSpeciesChange`
- `dMinGibbs`
- `dMinDrivingForceStoich`
- `dMinDrivingForceSoln`
- `dSpeciesRemove`
- `dPlateau`
- `lSubConverged`
- `lNegativeMolesPhase`
- `lGibbsMinCheck`

### Allocation and deallocation

Major allocation owners:

- `InitGEMSolver.f90`
- `InitGEMLineSearch.f90`
- `LevelingSolver.f90`
- `PostProcessPEA.f90`
- Gibbs-model routines that compute partial excess parameters

`ResetThermo` deallocates the main allocatable arrays:

- `iterHistory`
- `dSumMolFractionSoln`
- `dDrivingForceSoln`
- `dUpdateVar`
- `dPartialExcessGibbs`
- `dPartialExcessGibbsLast`
- `lSolnPhases`
- `lMiscibility`
- `dMolesPhaseLast`
- `iAssemblageGEM`
- `lPhaseChangeHistory`
- `iAssemblageGEMinit`

Review note:

- `dUpdateVarLast`, `dMolesSpeciesLast`, `dElementPotentialLast`,
  `dPartialGParam`, `dPartialHParam`, `dPartialSParam`, and `dPartialCpParam`
  should be checked against all allocation paths. Some are deallocated by local
  initializers rather than clearly by `ResetThermo`.
- Solver scalar status is now reset by `ResetThermo`; debugging should treat a
  nonzero `dGEMFunctionNorm` after reset as an error.

## ModuleSubMin

`ModuleSubMin` stores the solution-phase subminimization state. This is
especially relevant for Phase 1 because the current plan includes debugging
zeroing endmembers, solution Gibbs energies, chemical potentials, and
Lagrangian max-iteration cases.

### Variables

| Group | Variables |
|---|---|
| Integer scalars | `nVar`, `iFirstSUB`, `iLastSUB`, `iSolnPhaseIndexOther`, `iterSub`, `iterSubLg`, `iterSubAdam`, `iAdamNeg` |
| Integer arrays | `iHessian` |
| Real scalars | `dDrivingForce`, `dDrivingForceLast`, `dSubMinFunctionNorm`, `dConverge`, `dSumPairs`, `dSubminGibbsEst`, `dMaxPotentialVector` |
| Real parameters | `dSubMinTolerance`, `dMinMoleFraction`, `dTolEuclideanNorm`, `dTolDrivingForceChange`, `alpha`, `beta1`, `beta2`, `epsilon` |
| Real arrays | `dChemicalPotentialStar`, `dRHS`, `dRHSLast`, `dPiDot`, `dPiLast`, `dPi`, `dLambda`, `dLambdaLast`, `dvAdam`, `dsAdam`, `dvAdamLast`, `dsAdamLast`, `dxDot`, `dPotentialVector`, `dHessian` |
| Logical scalars | `lSubMinConverged`, `lNegativeFraction` |

### Observed use

Active fields include:

- `nVar`, `iFirstSUB`, `iLastSUB`
- `dDrivingForce`, `dDrivingForceLast`, `dSubMinFunctionNorm`
- `dChemicalPotentialStar`, `dRHS`, `dRHSLast`
- `dPi`, `dPiLast`, `dPiDot`
- `dLambda`, `dLambdaLast`
- Adam/line-search vectors
- `dPotentialVector`
- `dHessian`, `iHessian`
- `lSubMinConverged`, `lNegativeFraction`

Weak/no observed use in the current static pass:

- `dConverge`
- `dSumPairs`

### Allocation and deallocation

`SubMinInit.f90` deallocates and reallocates most subminimization arrays:

- `dChemicalPotentialStar`
- `dRHS`
- `dRHSLast`
- `dHessian`
- `iHessian`
- `dPi`
- `dPiLast`
- `dPiDot`
- `dLambda`
- `dLambdaLast`
- `dvAdam`
- `dsAdam`
- `dvAdamLast`
- `dsAdamLast`
- `dxDot`
- `dPotentialVector`

`Subminimization.f90` deallocates only:

- `dChemicalPotentialStar`
- `dRHS`
- `dHessian`
- `iHessian`

`ResetThermo` deallocates:

- `iHessian`
- `dChemicalPotentialStar`
- `dRHS`
- `dRHSLast`
- `dHessian`
- `dPi`
- `dPiLast`
- `dPiDot`
- `dLambda`
- `dLambdaLast`
- `dvAdam`
- `dsAdam`
- `dvAdamLast`
- `dsAdamLast`
- `dxDot`
- `dPotentialVector`

## Reset coverage summary

| Area | Current coverage | Concern |
|---|---|---|
| Parser/database arrays | `ResetThermoParser` covers the active parser groups | Some declared legacy/model arrays have no observed current allocation or reset |
| Main thermodynamic arrays | `ResetThermo` covers most active `ModuleThermo` arrays | Failure paths before reset can keep process-global state alive |
| Output gram arrays | Allocated in `PostProcess`, deallocated in `PostProcess` and `ResetThermo` | Looks reasonably covered |
| Output pair arrays | `ModuleThermoIO.dPair`/`cPair` are deallocated before the next `CalculateCompositionSUBG` allocation | Not globally reset today; distinct from `ModuleThermo.cPairName`, which is reset |
| GEM solver arrays and scalar status | Mostly covered by `ResetThermo` and initializer cleanup; core scalar counters, norms, indexes, and logical flags are reset by `ResetThermo` | Some line-search/partial arrays need exact allocation-path review |
| Subminimization arrays | Partly covered by `ResetThermo`; mostly covered by `SubMinInit` on next entry | Several allocated arrays are not globally reset |

## Phase 1 follow-up checklist

1. Add a small Fortran or Python-accessible diagnostic routine that reports
   `allocated()` status for each major allocatable group after a calculation and
   after each reset routine.
2. Expand `ResetThermo` or add dedicated reset helpers so every
   `ModuleSubMin` allocatable array and the `ModuleThermoIO` pair arrays have
   global reset coverage.
3. Review `ModuleGEMSolver` arrays that are locally deallocated but not clearly
   handled by `ResetThermo`.
4. Remove or mark dormant parser/model variables only after targeted tests for
   SUBL, SUBG, SUBQ, SUBI, magnetic, and SRO cases.
5. Instrument the flaky parallel batch path to log worker process id, reset
   calls, selected system dimensions, allocated-state diagnostics, and final
   `dGEMFunctionNorm`.
6. Decide where `_pyvar2fvar()` should live and keep only one implementation to
   avoid duplicate f2py transfer logic.

## Phase 2 Batch-Flakiness Allocation Map

This section narrows the static source inventory to state that can plausibly
explain `equilib_batch()` differences between serial execution, joblib/loky
workers, and multiprocessing-style workers.

### Source Files Reviewed

| Area | Files / routines |
|---|---|
| Module declarations | `ModuleParseCS.f90`, `ModuleThermoIO.f90`, `ModuleThermo.f90`, `ModuleGEMSolver.f90`, `ModuleSubMin.f90` |
| Reset routines | `ResetThermo.f90`, `ResetThermoParser.f90`, `ResetThermoAll.f90` |
| System setup | `CheckSystem.f90`, `CheckSystemExcess.f90`, `InitThermo.f90` |
| PEA / GEM setup | `CompThermoData.f90`, `InitGEMSolver.f90`, `InitGEMSolverNew.f90`, `LevelingSolver.f90`, `PostProcessPEA.f90`, `PostProcess.f90` |
| Solver / submin | `GEMSolver.f90`, `MultiPhaseMinimizer.f90`, `GEMNewton.f90`, `GEMLineSearch.f90`, `SubMinInit.f90`, `Subminimization.f90` |
| Gibbs model outputs | `CalculateCompositionSUBG.f90`, `CompExcessGibbsEnergySUBG.f90`, other `4_GibbsModels/*.f90` routines |

### Current Reset Coverage by Owner

| Owner | Allocated / written in | Reset / deallocated in | Batch-flakiness risk |
|---|---|---|---|
| `ModuleParseCS` parser arrays | f2py assignment in `utils._pyvar2fvar()`; original Fortran parser paths | `ResetThermoParser.f90`; overwritten by later f2py assignment | High if a reused worker starts with incompatible parser dimensions and `_pyvar2fvar()` does not overwrite every active field. |
| `ModuleThermoIO` input scalars and `dElementMass` | `input_condition.py` | `input_condition.py` overwrites active fields; `ResetThermo.f90` zeros `dElementMass` | Medium. Input appears actively overwritten every row, but stale scalar fields not written by Python could persist. |
| `ModuleThermoIO` gram output arrays | `PostProcess.f90` | `PostProcess.f90` before allocation; `ResetThermo.f90` | Medium. These are result-capture inputs, so stale values matter if `postprocess()` fails or is skipped. |
| `ModuleThermoIO.dPair/cPair` | `CalculateCompositionSUBG.f90` | Deallocated before next `CalculateCompositionSUBG` allocation | Medium. Not currently covered by `ResetThermo`; likely output-only, but stale SUBG pair output should be checked. |
| `ModuleThermo` selected-system arrays | `CheckSystem.f90`, `CheckSystemExcess.f90`, `CompThermoData.f90` | `ResetThermo.f90`; many setup routines deallocate before reallocating | High. Phase/species dimensions and maps define the solver system and result capture. |
| `ModuleThermo.iAssemblage`, `dMolesPhase`, `dMolesSpecies`, `dMolFraction` | `InitGEMSolver.f90`, `InitGEMSolverNew.f90`, GEM/PEA routines | `ResetThermo.f90`; some initializers deallocate before reallocating | High. These are core solver and result-capture arrays. |
| `ModuleGEMSolver` history/update arrays | `InitGEMSolver.f90`, `InitGEMLineSearch.f90`, `LevelingSolver.f90`, `PostProcessPEA.f90` | Mostly `ResetThermo.f90` and local initializers; scalar counters and flags are reset by `ResetThermo.f90` | High. Worker reuse can expose stale convergence/history array state if a path exits early. |
| `ModuleGEMSolver.dGEMFunctionNorm` | `CompFunctionNorm.f90`, `CompInitMinSolnPoint.f90`, line-search/GEM paths | Reset to zero by `ResetThermo.f90`; expected to be overwritten by solver | High. Python uses it as a failure criterion after minimization. Instrument before/after solve. |
| `ModuleSubMin` work arrays | `SubMinInit.f90` | Mostly `SubMinInit.f90` on next entry; partial `ResetThermo.f90`; partial `Subminimization.f90` cleanup | High. Several arrays allocated in `SubMinInit` are not globally reset. This is a strong candidate for allocated-state diagnostics. |

### High-Risk Variables to Instrument First

| Variable | Module | Why first |
|---|---|---|
| `iPhaseCS` | `ModuleParseCS` | Phase pass mask transferred from Python and narrowed by phase selection. |
| `iSolnPS` | `ModuleThermo` | Solution phase-selection mask written by `phase_selection()`. |
| `iPhase`, `nSpeciesPhase` | `ModuleThermo` | Active phase/species structure used by all solver paths. |
| `iAssemblage` | `ModuleThermo` | Stable phase IDs read by result capture and transition logic. |
| `dMolesPhase`, `dMolesSpecies`, `dMolFraction` | `ModuleThermo` | Core phase/species amounts read by post-processing. |
| `dChemicalPotential`, `dElementPotential` | `ModuleThermo` | Solver state that can affect convergence and output. |
| `dGEMFunctionNorm` | `ModuleGEMSolver` | Python success/failure gate. |
| `iterGlobal`, `iterPEA`, `lConverged`, `lPhaseChange` | `ModuleGEMSolver` | Confirm solver counters/flags are initialized for each row. |
| `dRHS`, `dHessian`, `dRHSLast`, `dPi`, `dLambda`, `dPotentialVector` | `ModuleSubMin` | Some submin arrays have initializer-owned cleanup but incomplete global reset coverage. |
| `dGramPhase`, `dGramFraction`, `dGramElement` | `ModuleThermoIO` | Result-capture outputs; should be fresh after `postprocess()`. |
| `dPair`, `cPair` | `ModuleThermoIO` | SUBG output arrays not covered by `ResetThermo`. |

### Static Findings Relevant to Joblib vs Multiprocessing

- Fortran module variables are `SAVE`d process globals. A joblib/loky worker
  process can keep them between tasks if the worker is reused.
- Multiprocessing can appear safer if the old execution pattern created or
  destroyed processes more often. That would mask, not fix, missing lifecycle
  ownership.
- `ResetThermo` is broad but not complete for every allocatable work array. Some
  arrays rely on the next initializer to deallocate before allocation.
- `ResetThermoParser` clears parser arrays but does not clear thermodynamic or
  solver arrays.
- `ResetThermoAll` is a composition of `ResetThermo` and
  `ResetThermoParser`. It still does not reset Python `equilipy.variables`.
- Therefore, if joblib differs from multiprocessing, the first question is
  which process-global Python/Fortran variable survives in the reused worker and
  is not deterministically overwritten before `minimize()`.

### Phase 2 Exit Notes

The source map does not justify adding broad per-row `resetthermoall()`.
It does justify Phase 3 reproduction and Phase 4 instrumentation around:

- phase-selection maps before `minimize()`;
- parser and selected-system dimensions after `list_phases()`;
- `ModuleSubMin` allocated-state coverage;
- `dGEMFunctionNorm` and `iAssemblage` immediately before and after `minimize()`.
