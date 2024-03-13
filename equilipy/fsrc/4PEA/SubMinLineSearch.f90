!
subroutine SubMinLineSearch(iSolnPhaseIndex)
    !---------------------------------------------------------------------------
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this subroutine is to perform a line search on the
    ! direction vector computed by the SubMinNewton subroutine.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! dStepLength   A double real scalar representing the steplength that is
    !                applied to updating the mole fractions.
    ! dTemp         A double real scalar representing a temporary variable.
    ! dMaxChange    A double real scalar representing the maximum change in the
    !                mole fraction of a solution phase constituent.
    !
    !---------------------------------------------------------------------------
    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleSubMin
!
    implicit none
!
    integer :: i, j, k, iSolnPhaseIndex
    real(8) :: dStepLength, dTemp, dMaxChange
!
!
    ! Initialize variables:
    dStepLength = 1D0
    dMaxChange  = 0.2D0
    ! print*,'dRHS',dRHS
    ! Initialize steplength (determine a steplength that prevents negative mole fractions
    ! and constrains maximum change):
    do j = 1, nVar
!
        ! Absolute species index:
        i = iFirstSUB + j - 1
!
        ! Check if the mole fraction of this constituent is driven to be negative:
        if (dMolFraction(i) + dRHS(j) <= 0D0) then
!
            ! Determine a step length that reduces the mole fraction by a factor of 100:
            dTemp = -0.99D0 * dMolFraction(i) / dRHS(j)
!
            ! Update the step length:
            dStepLength = DMIN1(dStepLength, dTemp)
!
        end if
!
        ! Compute a step length that constrains the maximum change by dMaxChange:
        dTemp=dStepLength
        if(dRHS(j) /= 0D0)then
           dTemp       = DABS(dMaxChange / dRHS(j))
        end if
        dStepLength = DMIN1(dStepLength, dTemp)
!
    end do
!
    ! Update the Species amount:
    dTemp = 0D0
!
    do j = 1, nVar
!
        ! Absolute species index:
        i = iFirstSUB + j - 1
!
        ! Apply step length:
        dMolFraction(i) = dMolFraction(i) + dStepLength * dRHS(j)
!
        ! Store maximum change to the mole fraction:
        dTemp = DMAX1(dTemp, DABS(dRHS(j)))
!
    end do
!
    ! Iterate to satisfy Wolfe conditions:
    LOOP_WOLFE: do k = 1, 5
!
        ! Exit if the minimum mole fraction of this phase is below a specified tolerance:
        if (MINVAL(dMolFraction(iFirstSUB:iLastSUB)) < dMinMoleFraction) exit LOOP_WOLFE
!
        ! Compute the chemical potentials of solution phase constituents:
        call SubMinChemicalPotential(iSolnPhaseIndex)
!
        ! Compute the driving force of this solution phase:
        call SubMinDrivingForce
!
        ! Check for a converging/diverging solution:
        if (dDrivingForce < dDrivingForceLast) then
!
            ! The driving force is less than the last driving force:
            exit LOOP_WOLFE
!
        else
            ! Divergence has been detected.  Dampen the sub-system:
            !SYMF
            dStepLength = dStepLength * 0.5D0
            dTemp       = 0D0
!
            ! Adjust the mole fractions of species:
            do j = 1, nVar
!
                ! Absolute species index:
                i = iFirstSUB + j - 1
!
                ! Apply step length:
                dMolFraction(i) = dMolFraction(i) - dStepLength * dRHS(j)
!
                ! Constrain the minimum mole fraction to an arbitrarily small value:
                dMolFraction(i) = DMAX1(dMolFraction(i), 1D-100)
!
                dTemp = DMAX1(dTemp, dStepLength * dRHS(j))
!
            end do
!
        end if
!
    end do LOOP_WOLFE
!
    ! Check convergence (maximum change to any mole fraction):
    if (dTemp <= dSubMinTolerance) lSubMinConverged = .TRUE.
!
    ! Store the driving force from the last iteration:
    dDrivingForceLast = dDrivingForce
!
end subroutine SubMinLineSearch
!
!
