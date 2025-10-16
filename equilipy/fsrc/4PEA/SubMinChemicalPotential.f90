!
subroutine SubMinChemicalPotential(iSolnPhaseIndex)
    !---------------------------------------------------------------------------
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this subroutine is to compute the chemical potentials of
    ! only the constituents of this solution phase.
    !
    ! Pertinent variables:
    ! ====================
    !
    ! dChemicalPotential    A double real vector representing the chemical
    !                        potential of every species in the system.
    ! dPartialExcessGibbs   A double real vector representing the partial
    !                        molar excess Gibbs energy of every species in
    !                        the system.
    !
    !---------------------------------------------------------------------------
    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleSubMin
!
    implicit none
!
    integer::  iSolnPhaseIndex
!
    ! Initialize variables:
    dChemicalPotential(iFirstSUB:iLastSUB)  = 0D0
    dPartialEnthalpy(iFirstSUB:iLastSUB)  = 0D0
    dPartialEntropy(iFirstSUB:iLastSUB)  = 0D0
    dPartialHeatCapacity(iFirstSUB:iLastSUB)  = 0D0
!
!
    ! Compute excess terms based on solution phase type:
    call CompExcessGibbsEnergy(iSolnPhaseIndex)
!
end subroutine SubMinChemicalPotential
!
!
