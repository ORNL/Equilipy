!
subroutine SubMinDrivingForce
    !---------------------------------------------------------------------------
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this subroutine is to compute the driving force of a
    ! particular solution phase within the Subminimization routine.
    !
    ! Pertinent variables:
    ! ====================
    !
    ! dDrivingForce     A double real scalar representing the driving force of
    !                    a solution phase.  The driving force represents the
    !                    difference between the molar Gibbs energy of that
    !                    particular phase and the corresponding value calculated
    !                    from the chemical potentials of the system components.
    !                    When negative, this term indicates that a phase should
    !                    be added to the system.  When positive, the phase
    !                    should not be added to hte system.
    !
    !---------------------------------------------------------------------------
    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleSubMin
!
    implicit none
!
    integer::  i, j
    real(8):: dTemp
!
!
    ! Initialize variables:
    dDrivingForce = 0D0
    dTemp         = 0D0
!
    ! Loop through all species in this phase:
    do j = 1, nVar
!
        ! Absolute species index:
        i = iFirstSUB + j - 1
!
        ! Update the driving force:
        dDrivingForce = dDrivingForce + dMolFraction(i) * (dChemicalPotential(i) - dChemicalPotentialStar(j))
!
        dTemp = dTemp + sum(dMolfraction(i)*dStoichSpecies(i,:))/DFLOAT(iParticlesPerMole(i))
    end do
!
    ! Normalize the driving force as to one mole of atom
    dDrivingForce=dDrivingForce/dTemp
    
!
end subroutine SubMinDrivingForce
!
!
