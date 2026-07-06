!-------------------------------------------------------------------------------------------------------------
!
!> \file    ModuleThermoIO.f90
!> \brief   Shared input, output, and status state for Equilipy/Thermochimica calls.
!
!
! Revisions:
! ==========
!
!   Date            Programmer          Description of change
!   ----            ----------          ---------------------
!   10/17/2011      M.H.A. Piro         Original code
!   04/23/2012      M.H.A. Piro         Added dOxygen documentation
!   06/23/2026      S.Y. Kwon           Removed unused legacy Thermochimica I/O buffers
!   06/23/2026      S.Y. Kwon           Organized declarations by rank and role
!   06/26/2026      S.Y. Kwon           Added dependent element flags for pseudo-component bases
!
!
! Purpose:
! ========
!
!> \details This module stores the process-global input and output state used
!! by the Fortran runtime and exposed through the f2py module.  Python writes
!! the input condition state before minimization.  Fortran updates the status,
!! thermodynamic properties, postprocess diagnostics, and selected result
!! arrays after minimization.
!
!
! Required input variables:
! =========================
!
! dTemperature              Absolute temperature.  CheckThermoInput normalizes this to Kelvin.
! dPressure                 Hydrostatic pressure.  CheckThermoInput normalizes this to atm.
! dTemperatureDiff          Temperature perturbation used by PostProcess heat-capacity calculation.
! dElementMass              Input amount for each periodic-table element index.
! iDependentElementInput    Marks elements admitted for species screening but excluded from active mass balance.
! cInputUnitTemperature     Temperature unit: K, C, F, or R.
! cInputUnitPressure        Pressure unit: atm, psi, bar, Pa, or kPa.
! cInputUnitMass            Composition unit, normalized to moles by CheckSystem.
! iPrintResultsMode         PrintResults mode.  Zero suppresses Fortran console output.
!
!
! Output/updated variables:
! =========================
!
! INFOThermo                    Zero for success, otherwise a ThermoDEBUG error code.
! dGibbsEnergySys               Integral system Gibbs energy after PostProcess.
! dEnthalpySys                  Integral system enthalpy after PostProcess.
! dEntropySys                   Integral system entropy after PostProcess.
! dHeatCapacitySys              Integral system heat capacity after PostProcess.
! dPostProcessGEMNormAtT        GEM norm at the requested temperature.
! dPostProcessGEMNormPerturbed  GEM norm at the T-dT heat-capacity perturbation.
! iPostProcessIterAtT           Lagrangian iteration count at the requested temperature.
! iPostProcessIterPerturbed     Lagrangian iteration count at the T-dT perturbation.
! lPostProcessAssemblageChanged True when the T-dT postprocess solve changes assemblage.
! dGramFraction                 Species mass fractions within each solution phase.
! dGramSpecies                  Species masses after PostProcess.
! dGramPhase                    Stable phase masses after PostProcess.
! dGramElement                  Element masses after PostProcess.
! cPair                         SUBG/SUBQ pair names from CalculateCompositionSUBG.
! dPair                         SUBG/SUBQ pair fractions from CalculateCompositionSUBG.
!
!
! Primary callers/users:
! ======================
!
! input_condition.py        Writes temperature, pressure, units, and input element amounts.
! CheckThermoInput          Validates and normalizes temperature/pressure units.
! CheckSystem               Validates and normalizes element amounts.
! PostProcess               Writes integral properties, diagnostics, and gram outputs.
! CalculateCompositionSUBG  Writes pair output arrays for SUBG/SUBQ phases.
! PrintResults              Uses status and print mode for optional console output.
! ThermoDEBUG               Interprets INFOThermo error codes.
! Python result capture     Reads system properties and gram outputs through f2py.
!
!
! Numerical assumptions:
! ======================
!
! - dElementMass is indexed by periodic-table element number, with index zero
!   available for compatibility with legacy Thermochimica conventions.
! - iDependentElementInput uses the same periodic-table index as dElementMass.
! - CheckThermoInput and CheckSystem convert user units to the internal units
!   before thermodynamic data and minimization routines consume the values.
! - PostProcess rescales output properties by dNormalizeInput before Python
!   result capture reads this module.
!
!
!-------------------------------------------------------------------------------------------------------------
module ModuleThermoIO

    implicit none

    SAVE

    ! Integer scalars:
    integer                          :: iPrintResultsMode
    integer                          :: INFOThermo
    integer                          :: iPostProcessIterAtT, iPostProcessIterPerturbed

    ! Integer 1D arrays:
    integer, dimension(0:118)         :: iDependentElementInput

    ! Real scalars:
    real(8)                          :: dTemperature, dPressure, dTemperatureDiff
    real(8)                          :: dGibbsEnergySys, dEnthalpySys, dEntropySys, dHeatCapacitySys
    real(8)                          :: dPostProcessGEMNormAtT, dPostProcessGEMNormPerturbed

    ! Logical scalars:
    logical                          :: lPostProcessAssemblageChanged = .FALSE.

    ! Character scalars:
    character(15)                    :: cInputUnitTemperature, cInputUnitPressure, cInputUnitMass


    ! Real 1D arrays:
    real(8),       dimension(0:118)          :: dElementMass
    real(8), dimension(:), allocatable       :: dPair
    real(8), dimension(:), allocatable       :: dGramFraction, dGramSpecies, dGramPhase, dGramElement

    ! Character 1D arrays:
    character(30), dimension(:), allocatable :: cPair

end module ModuleThermoIO
