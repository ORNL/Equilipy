


subroutine ResetThermoParser

    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ResetThermoParser.f90
    !> \brief   Deallocate allocatable variables used by the ModuleParseCS.f90 module.
    !> \author  M.H.A. Piro
    !> \date    Apr. 26, 2012
    !> \sa      ModuleParseCS.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   02/17/2012      M.H.A. Piro         Original code
    !   07/20/2026      S.Y. Kwon           Reset parser state including structural order/disorder partition metadata.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to deallocate all allocatable arrays from the
    !! ModuleParseCS.f90 module, which are allocated in the ParseCS*.f90 subroutines and then used in several
    !! subroutines by Thermochimica.  A value of INFOThermo = 18 is returned if an error has occured during
    !! deallocation.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! INFO                  An error is returned if deallocation is unsuccessful.
    ! INFOThermo            An integer scalar identifying whether the program exits successfully or if
    !                       it encounters an error.  A description for each error is given in ThermoDebug.f90.
    !
    !-------------------------------------------------------------------------------------------------------------
    USE ModuleParseCS
    USE ModuleThermoIO, ONLY: INFOThermo

    implicit none

    integer::   i


    ! Initialize variables:
    i    = 0
    INFO = 0
    INFOThermo = 0
    nElementsCS = 0
    nSpeciesCS = 0
    nSolnPhasesSysCS = 0
    iMiscSUBI = 0
    nParamCS = 0
    nCountSublatticeCS = 0
    nMaxSpeciesPhaseCS = 0
    nMagParamCS = 0
    lEndmembers2Species = .FALSE.

    ! Deallocate each parser array under its own guard; Python-side database
    ! transfer paths can allocate selected arrays without allocating the
    ! historical group sentinel first.
    INFO = 0
    if (allocated(nSpeciesPhaseCS)) deallocate(nSpeciesPhaseCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(nGibbsEqSpecies)) deallocate(nGibbsEqSpecies, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(iPhaseCS)) deallocate(iPhaseCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(iParticlesPerMoleCS)) deallocate(iParticlesPerMoleCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(nParamPhaseCS)) deallocate(nParamPhaseCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(iParamPassCS)) deallocate(iParamPassCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(nSublatticePhaseCS)) deallocate(nSublatticePhaseCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(iPhaseSublatticeCS)) deallocate(iPhaseSublatticeCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(iDisorderedPhaseCS)) deallocate(iDisorderedPhaseCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(iOrderDisorderTopologyCS)) deallocate(iOrderDisorderTopologyCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(iOrderDisorderStandalonePhaseCS)) then
        deallocate(iOrderDisorderStandalonePhaseCS, STAT = INFO)
    end if
    i = i + INFO
    INFO = 0
    if (allocated(iMagParamPassCS)) deallocate(iMagParamPassCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(nMagParamPhaseCS)) deallocate(nMagParamPhaseCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(iSUBIMixTypeCS)) deallocate(iSUBIMixTypeCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(nInterpolationOverrideCS)) deallocate(nInterpolationOverrideCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(dAtomicMassCS)) deallocate(dAtomicMassCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(cElementNameCS)) deallocate(cElementNameCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(cSolnPhaseTypeCS)) deallocate(cSolnPhaseTypeCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(cSolnPhaseNameCS)) deallocate(cSolnPhaseNameCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(cSpeciesNameCS)) deallocate(cSpeciesNameCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(cRegularParamCS)) deallocate(cRegularParamCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(iRegularParamCS)) deallocate(iRegularParamCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(nConstituentSublatticeCS)) deallocate(nConstituentSublatticeCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(nPairsSROCS)) deallocate(nPairsSROCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(iMagneticParamCS)) deallocate(iMagneticParamCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(iSUBIParamDataCS)) deallocate(iSUBIParamDataCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(nSublatticeElementsCS)) deallocate(nSublatticeElementsCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(dGibbsCoeffSpeciesTemp)) deallocate(dGibbsCoeffSpeciesTemp, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(dRegularParamCS)) deallocate(dRegularParamCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(dGibbsMagneticCS)) deallocate(dGibbsMagneticCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(dMagneticParamCS)) deallocate(dMagneticParamCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(dStoichSublatticeCS)) deallocate(dStoichSublatticeCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(dStoichSpeciesCS)) deallocate(dStoichSpeciesCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(dZetaSpeciesCS)) deallocate(dZetaSpeciesCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(dStoichConstituentCS)) deallocate(dStoichConstituentCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(dQKTOParamsCS)) deallocate(dQKTOParamsCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(cPairNameCS)) deallocate(cPairNameCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(iInterpolationOverrideCS)) deallocate(iInterpolationOverrideCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(iConstituentSublatticeCS)) deallocate(iConstituentSublatticeCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(iPairIDCS)) deallocate(iPairIDCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(iChemicalGroupCS)) deallocate(iChemicalGroupCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(dSublatticeChargeCS)) deallocate(dSublatticeChargeCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(dStoichPairsCS)) deallocate(dStoichPairsCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(dConstituentCoefficientsCS)) deallocate(dConstituentCoefficientsCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(dCoordinationNumberCS)) deallocate(dCoordinationNumberCS, STAT = INFO)
    i = i + INFO
    INFO = 0
    if (allocated(cConstituentNameSUBCS)) deallocate(cConstituentNameSUBCS, STAT = INFO)
    i = i + INFO

    ! Return an INFOThermo if deallocation of any of the allocatable variables failed:
    if (i > 0) then
        INFOThermo = 18
    end if

    return

end subroutine ResetThermoParser
