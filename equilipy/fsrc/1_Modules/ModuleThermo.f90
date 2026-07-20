
module ModuleThermo
!
!-------------------------------------------------------------------------------------------------------------
!
!> \file    ModuleThermo.f90
!> \brief   Shared thermodynamic system, phase, species, and property state.
!> \author  M.H.A. Piro
!> \sa      ModuleGEMSolver.f90
!> \sa      ModuleParseCS.f90
!> \sa      InitThermo.f90
!> \sa      CheckSystem.f90
!> \sa      CompThermoData.f90
!> \sa      PostProcess.f90
!
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   04/25/2012      M.H.A. Piro         Original shared thermodynamic state module
    !   07/20/2026      S.Y. Kwon           Separated thermodynamic derivative storage and structural order/disorder phase identities.
    !
    !
! Purpose:
! ========
!
!> \details ModuleThermo stores the parsed thermodynamic system definition,
!! current thermodynamic state, phase/species amounts, and property work arrays.
!! Runtime minimizer controls, active-set repair state, and PEA candidate pools
!! are owned by ModuleGEMSolver.
!
!
! Required input variables:
! =========================
!
! nElements, nSpecies              Active element and species counts supplied by InitThermo/CheckSystem.
! nSpeciesPhase, iPhase            Phase/species topology used by Leveling, PEA, and Lagrangian GEM.
! iParticlesPerMole                Formula-particle scaling for species Gibbs energy and chemical potentials.
! dStoichSpecies                   Elemental stoichiometry for each species.
! dLevelingSpeciesFormulaAtoms       Formula-atom scaling used when Leveling compares ordered/disordered endmembers.
! dStdGibbsEnergy                  Species standard Gibbs energies after database and model evaluation.
! dExcessGibbsParam                Excess Gibbs parameters for solution models.
! nSublatticePhase                 Sublattice counts for CEF/order-disorder phases.
! iDisorderedPhase                 Thermodynamic order/disorder DIS_PART helper mapping.
! iODCompanionPhase                Solver-identity map to the parser-declared DIS_PART helper.
! iODStandalonePhase               Independent standalone disordered phase, when separately declared.
! iODTopologyClass                 Typed structural DIS_PART graph classification by ordered phase.
! iODProjectionTopology            Supported ordering-projector topology by ordered phase.
!
!
! Output/updated variables:
! =========================
!
! dElementPotential                Element potentials from Leveling/Lagrangian solves.
! dChemicalPotential               Species residual potentials used by minimizer convergence checks.
! dLevelingChemicalPotential       Leveling chemical potentials corrected for order/disorder and formula scaling.
! dMolesElement                    Bulk element amounts after input normalization.
! dMolesPhase, dMolesSpecies       Current phase and species amounts.
! dMolFraction                     Current phase-local species mole fractions.
! dAtomFractionSpecies             Species atom fractions used by Leveling and mass-balance equations.
! dPartialExcessGibbs              Partial excess Gibbs energy of solution species.
! dPartialEnthalpy/Entropy/HeatCapacity
!                                  Partial thermodynamic properties for post-processing.
! dSiteFraction                    Sublattice site fractions for CEF/order-disorder phases.
!
!
! Primary callers/users:
! ======================
!
! InitThermo                       Initializes counts, constants, and tolerances.
! CheckSystem                      Converts Python/database inputs into normalized Fortran state.
! CompThermoData                   Evaluates Gibbs, enthalpy, entropy, heat capacity, and magnetic terms.
! LevelingSolver                   Uses species Gibbs energies and compositions for initial active-set search.
! Level2Lagrange                   Converts Leveling/PEA rows into phase/species amounts.
! RunLagrangianGEM                 Updates thermodynamic state for the active assemblage.
! PostProcess                      Converts normalized internal state back to user-facing properties.
! ResetThermo                      Deallocates and resets thermodynamic/module state between calculations.
!
!
! Numerical assumptions:
! ======================
!
! - dMolesElement is normalized for numerical conditioning; PostProcess rescales
!   output properties using dNormalizeInput.
! - Vacancy is represented in species/endmember topology, but conserved-element
!   mass balance uses the active thermodynamic elements supplied to CheckSystem.
! - Ordered/disordered species Gibbs energies are stored per formula as required
!   by the phase model; Leveling-specific comparison scaling is tracked through
!   dLevelingSpeciesFormulaAtoms and dLevelingChemicalPotential.
!
!-------------------------------------------------------------------------------------------------------------
    implicit none
!
    SAVE
!
    integer,       parameter                         :: iTolNum = 15, nElementsPT = 118

    integer::   nElements, nSpecies, nParam, nMaxParam, nDummySpecies, nElemOrComp, nMagParam
    integer::   nSpeciesLevel
    integer::   nConPhases, nSolnPhases, nSolnPhasesSys, nChargedConstraints, nElementsSys !Added
    integer::   nMaxSublatticeSys, nMaxConstituentSys, nCountSublattice

    integer,       dimension(:),   allocatable::  iPhase, nSpeciesPhase, iParticlesPerMole, iPhaseLevel, iCandidate,iSolnPS
    integer,       dimension(:),   allocatable::  iAssemblage, nParamPhase, iElementSystem, iSpeciesPass, nMagParamPhase
    integer,       dimension(:),   allocatable::  nSublatticePhase, iPhaseSublattice, iDisorderedPhase
    integer,       dimension(:),   allocatable::  iODCompanionPhase, iODStandalonePhase
    integer,       dimension(:),   allocatable::  iODTopologyClass, iODProjectionTopology
    integer,       dimension(:),   allocatable::  iPhaseElectronID
    integer,       dimension(:),   allocatable::  iShuffled
    ! integer,       dimension(:),   allocatable::  iSub2SysSoln, iSub2SysComp, iSub2SysSpecies! the whole line is added
    integer,       dimension(:,:), allocatable::  iRegularParam, nConstituentSublattice, nPairsSRO, iMagneticParam
    integer,       dimension(:,:), allocatable::  iSUBLParamData
    integer,       dimension(:,:), allocatable::  nSublatticeElements
    integer,       dimension(:,:,:),allocatable:: iConstituentPass, iConstituentSublattice, iPairID
    integer,       dimension(:,:,:),allocatable:: iChemicalGroup
!
    real(8),       parameter                  ::  dDefaultPostLevelingZeroEndmemberFloor = 1D-30
    real(8)                                    ::  dIdealConstant, dNormalizeSum, dNormalizeInput
    real(8)                                    ::  dPostLevelingZeroEndmemberFloor = dDefaultPostLevelingZeroEndmemberFloor
    real(8)                                    ::  xT, dGParam, dHParam, dSParam, dCpParam
    real(8),       dimension(iTolNum)          ::  dTolerance
    real(8),       dimension(:),   allocatable ::  dGibbsSolnPhase, dMolesSpecies, dMagGibbsEnergy, dMagEnthalpy
    real(8),       dimension(:),   allocatable ::  dMagEntropy,dMagHeatCapacity
    real(8),       dimension(:),   allocatable ::  dStdGibbsEnergy, dStdEnthalpy, dStdEntropy, dStdHeatCapacity
    real(8),       dimension(:),   allocatable ::  dChemicalPotential, dActivity, dExcessGibbsParam
    real(8),       dimension(:),   allocatable ::  dSpeciesTotalAtoms, dLevelingSpeciesTotalAtoms
    real(8),       dimension(:),   allocatable ::  dLevelingSpeciesFormulaAtoms
    real(8),       dimension(:),   allocatable ::  dExcessHParam, dExcessSParam, dExcessCpParam
    real(8),       dimension(:),   allocatable ::  dPartialEnthalpy, dPartialEntropy, dPartialHeatCapacity
    real(8),       dimension(:),   allocatable ::  dPartialExcessGibbs
    real(8),       dimension(:),   allocatable ::  dPartialEnthalpyXS, dPartialEntropyXS, dPartialHeatCapacityXS
    real(8),       dimension(:),   allocatable ::  dElementPotential, dMolesPhase, dMolesElement, dMolFraction, dAtomicMass
    real(8),       dimension(:),   allocatable ::  dChemicalPotentialOld, dMolFractionOld
    real(8),       dimension(:),   allocatable ::  dLevelingChemicalPotential, dLevelingChemicalPotentialOld
    real(8),       dimension(:),   allocatable ::  dPartialGParam, dPartialHParam
    real(8),       dimension(:),   allocatable ::  dPartialSParam, dPartialCpParam

    real(8),       dimension(:,:), allocatable ::  dAtomFractionSpecies, dLevelingCompositionSpecies
    real(8),       dimension(:,:), allocatable ::  dStoichSublattice, dStoichSpecies
    real(8),       dimension(:,:), allocatable ::  dCoeffGibbsMagnetic, dZetaSpecies, dMagneticParam, dAtomFractionSpeciesOld
    real(8),       dimension(:,:), allocatable ::  dLevelingCompositionSpeciesOld
    real(8),       dimension(:,:), allocatable ::  dStoichSpeciesLevel

    real(8),      dimension(:,:,:),allocatable ::  dSiteFraction, dCoordinationNumber, dSublatticeCharge, dStoichPairs
    real(8),      dimension(:,:,:),allocatable ::  dConstituentCoefficients

    character(12), dimension(:),   allocatable  ::  cElementName
    character(30), dimension(:),   allocatable  ::  cSpeciesName
    character(8),  dimension(:),   allocatable   ::  cSolnPhaseType
    character(25), dimension(:),   allocatable  ::  cSolnPhaseName
    character(8),  dimension(:,:,:),allocatable   :: cConstituentNameSUB
    character,     dimension(:),   allocatable    :: cRegularParam
    character(30),  dimension(:,:), allocatable :: cPairName

    logical                                                       :: lSkipLagrange
    logical                                                       :: lOrderDisorderEvaluation = .FALSE.

end module ModuleThermo
