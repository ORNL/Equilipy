!
subroutine SubMinFunctionNorm(iSolnPhaseIndex)
    !---------------------------------------------------------------------------
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this subroutine is to compute the functional norm to
    ! test for convergence.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! INFO          An integer scalar used by LAPACK to identify a successful
    !                exit (INFO = 0) or an error (INFO /= 0).
    !
    !---------------------------------------------------------------------------
    USE ModuleThermo!, ONLY: dMolFraction, dStoichSpecies, iPhaseElectronID
    USE ModuleSubMin
!
    implicit none
!
    integer :: i, j
    integer :: iSolnPhaseIndex
    real(8)  :: dTemp, dDrivingForceTemp
!
!
    ! Initialize variables:
    dSubMinFunctionNorm = -1D0
    dTemp               = 0D0
    dDrivingForceTemp   = 0D0
!
    ! Compute residual of mole fractions:
    do i = iFirstSUB, iLastSUB
        dSubMinFunctionNorm = dSubMinFunctionNorm + dMolFraction(i)
    end do
!
    dSubMinFunctionNorm = DABS(dSubMinFunctionNorm)
    ! Compute residual of charge neutrality constraints:
    if (iPhaseElectronID(iSolnPhaseIndex) /= 0) then
!
        j = iPhaseElectronID(iSolnPhaseIndex)
!
        do i = iFirstSUB, iLastSUB
            dSubMinFunctionNorm = dSubMinFunctionNorm + dMolFraction(i) * dStoichSpecies(i,j)
        end do
    end if
!
    dSubMinFunctionNorm = DABS(dSubMinFunctionNorm)
    !Compute residual chemical potential (Driving force)

    do j = 1, nVar
!
        ! Absolute species index:
        i = iFirstSUB + j - 1
!
        ! Update the driving force:
        dDrivingForceTemp = dDrivingForceTemp + dMolFraction(i) * (dChemicalPotential(i) - dChemicalPotentialStar(j))
!
        dTemp = dTemp + sum(dMolfraction(i)*dStoichSpecies(i,:))/DFLOAT(iParticlesPerMole(i))
    end do
    
    ! Compute absolute quantity:
    ! Normalize the driving force as to one mole of atom
    dSubMinFunctionNorm = dSubMinFunctionNorm +DABS(dDrivingForceTemp/dTemp-dDrivingForce)
!
end subroutine SubMinFunctionNorm
!
!
