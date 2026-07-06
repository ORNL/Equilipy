
!-------------------------------------------------------------------------------------------------------------
!
!> \file    ModuleParseCS.f90
!> \brief   Parser-side ChemSage/TDB database state before transfer to runtime thermodynamic modules.
!
!
! Revisions:
! ==========
!
!   Date            Programmer          Description of change
!   ----            ----------          ---------------------
!   04/24/2012      M.H.A Piro          Original parser-state module
!   03/25/2019      P. Bajpai           Added SUBQ to supported models
!   06/23/2026      S.Y. Kwon           Added SUBOM to supported models
!   06/23/2026      S.Y. Kwon           Documented parser-to-runtime state contract
!   06/23/2026      S.Y. Kwon           Organized declarations by rank and role
!
!
! Purpose:
! ========
!
!> \details This module stores thermodynamic database data in the ChemSage
!! parser namespace before it is copied into the runtime modules used by
!! equilibrium, property, and minimizer routines.  Python parser/runtime code
!! fills these `*CS` arrays for ChemSage `.dat` and TDB-derived databases.
!! Fortran validation and setup routines then consume the parser state to build
!! `ModuleThermo` and related runtime arrays.
!
!
! Required input variables:
! =========================
!
! INFO                       Parser/setup status propagated to validation routines.
! nElementsCS                Number of elements in the parsed database system.
! nSpeciesCS                 Number of parsed species, including solution endmembers and pure phases.
! nSolnPhasesSysCS           Number of parsed solution phases.
! nSpeciesPhaseCS            Cumulative species/endmember offsets for each solution phase.
! iPhaseCS                   Phase ownership for each parsed species.
! iParticlesPerMoleCS        Particles per mole of formula for each species.
! dAtomicMassCS              Atomic mass by selected-system element.
! dStoichSpeciesCS           Species stoichiometry matrix.
! cElementNameCS             Element names.
! cSpeciesNameCS             Species/endmember names.
! cSolnPhaseNameCS           Solution phase names.
! cSolnPhaseTypeCS           Solution model type for each parsed solution phase.
!
!
! Output/updated variables:
! =========================
!
! nParamPhaseCS              Number of excess Gibbs parameters per solution phase.
! nParamCS                   Total excess Gibbs parameter count derived from nParamPhaseCS.
! iRegularParamCS            Excess Gibbs parameter constituent indexes and exponents.
! dRegularParamCS            Excess Gibbs parameter values.
! nMagParamPhaseCS           Number of magnetic parameters per solution phase.
! nMagParamCS                Total magnetic parameter count derived from nMagParamPhaseCS.
! iMagneticParamCS           Magnetic parameter metadata.
! dMagneticParamCS           Magnetic parameter values.
! nSublatticePhaseCS         Number of sublattices per CEF/SUBG/SUBQ/SUBOM phase.
! nConstituentSublatticeCS   Number of constituents on each sublattice.
! dStoichSublatticeCS        Sublattice stoichiometric coefficients.
! iConstituentSublatticeCS   Constituent indexes for each phase/sublattice.
! iDisorderedPhaseCS         Order/disorder phase-pair mapping for ordered CEF phases.
! dSublatticeChargeCS        Sublattice charges for ionic phases.
! nPairsSROCS, iPairIDCS     SUBG/SUBQ pair-count and pair-index data.
! dCoordinationNumberCS      SUBG/SUBQ coordination-number data.
!
!
! Primary callers/users:
! ======================
!
! Python parser/runtime loaders  Populate this module from `.dat` or TDB-derived database state.
! CheckThermoData                Validates parsed model types, parameters, and species metadata.
! CompThermoData                 Transfers parser state into runtime thermodynamic arrays.
! ResetThermoParser              Deallocates parser-side allocatable arrays.
!
!
! Numerical assumptions:
! ======================
!
! - `*CS` arrays are parser-side storage.  Runtime routines should use
!   `ModuleThermo` after setup rather than reading parser state directly.
! - `nParamCS` and `nMagParamCS` are derived from phase-level parameter counts
!   during Fortran setup, so Python does not need to map them directly.
! - `nSolnPhasesSysMax`, `nMaxSublatticeCS`, `nGibbsCoeff`, `nMaxGibbsEqs`,
!   and `nParamMax` define parser allocation/validation limits.
! - `SUBOM` identifies ordered/disordered CEF phases and must remain in
!   `cSolnPhaseTypeSupport` for order/disorder TDB databases.
!
!-------------------------------------------------------------------------------------------------------------


module ModuleParseCS

    implicit none

    SAVE

    ! Integer parameters:
    integer,        parameter                   :: nSolnPhasesSysMax = 500, nMaxSublatticeCS = 5
    integer,        parameter                   :: nSolnTypeSupport = 11
    integer,        parameter                   :: nGibbsCoeff = 13, nMaxGibbsEqs = 6, nParamMax = 4

    ! Character parameters:
    character(8),   dimension(nSolnTypeSupport), parameter :: cSolnPhaseTypeSupport = &
                                                    ['IDMX    ','QKTO    ','SUBL    ','RKMP    ','RKMPM   ','SUBLM   ','SUBG    ', &
                                                    'SUBQ    ','SUBI    ','SUBM    ','SUBOM   ']


    ! Integer scalars:
    integer                                     :: nElementsCS, nSpeciesCS, nSolnPhasesSysCS, INFO, iMiscSUBI
    integer                                     :: nParamCS, nCountSublatticeCS, nMaxSpeciesPhaseCS, nMagParamCS

    ! Logical scalars:
    logical                                     :: lEndmembers2Species


    ! Integer 1D arrays:
    integer,        dimension(:),   allocatable :: nSpeciesPhaseCS, nGibbsEqSpecies
    integer,        dimension(:),   allocatable :: iPhaseCS, iParticlesPerMoleCS
    integer,        dimension(:),   allocatable :: nParamPhaseCS, iParamPassCS
    integer,        dimension(:),   allocatable :: nSublatticePhaseCS, iPhaseSublatticeCS
    integer,        dimension(:),   allocatable :: iDisorderedPhaseCS
    integer,        dimension(:),   allocatable :: iMagParamPassCS, nMagParamPhaseCS, iSUBIMixTypeCS
    integer,        dimension(:),   allocatable :: nInterpolationOverrideCS

    ! Real 1D arrays:
    real(8),        dimension(:),   allocatable :: dAtomicMassCS

    ! Character 1D arrays:
    character(3),   dimension(:),   allocatable :: cElementNameCS
    character(8),   dimension(:),   allocatable :: cSolnPhaseTypeCS
    character(25),  dimension(:),   allocatable :: cSolnPhaseNameCS
    character(25),  dimension(:),   allocatable :: cSpeciesNameCS
    character,      dimension(:),   allocatable :: cRegularParamCS


    ! Integer 2D arrays:
    integer,        dimension(:,:), allocatable :: iRegularParamCS, nConstituentSublatticeCS
    integer,        dimension(:,:), allocatable :: nPairsSROCS, iMagneticParamCS
    integer,        dimension(:,:), allocatable :: iSUBIParamDataCS, nSublatticeElementsCS

    ! Real 2D arrays:
    real(8),        dimension(:,:), allocatable :: dGibbsCoeffSpeciesTemp, dRegularParamCS
    real(8),        dimension(:,:), allocatable :: dGibbsMagneticCS, dMagneticParamCS
    real(8),        dimension(:,:), allocatable :: dStoichSublatticeCS, dStoichSpeciesCS
    real(8),        dimension(:,:), allocatable :: dZetaSpeciesCS, dStoichConstituentCS
    real(8),        dimension(:,:), allocatable :: dQKTOParamsCS

    ! Character 2D arrays:
    character(30),  dimension(:,:), allocatable :: cPairNameCS


    ! Integer 3D arrays:
    integer,        dimension(:,:,:), allocatable :: iInterpolationOverrideCS
    integer,        dimension(:,:,:), allocatable :: iConstituentSublatticeCS, iPairIDCS, iChemicalGroupCS

    ! Real 3D arrays:
    real(8),        dimension(:,:,:), allocatable :: dSublatticeChargeCS, dStoichPairsCS
    real(8),        dimension(:,:,:), allocatable :: dConstituentCoefficientsCS
    real(8),        dimension(:,:,:), allocatable :: dCoordinationNumberCS

    ! Character 3D arrays:
    character(8),   dimension(:,:,:), allocatable :: cConstituentNameSUBCS

end module ModuleParseCS
