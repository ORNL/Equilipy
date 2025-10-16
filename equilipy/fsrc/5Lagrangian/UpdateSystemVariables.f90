subroutine UpdateSystemVariables(dStepLength)
    !---------------------------------------------------------------------------
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this subroutine is to update the system variables using
    ! a specified steplength.
    !
    !---------------------------------------------------------------------------
    USE ModuleThermo
    USE ModuleGEMSolver
!
    integer                       :: i, j, k
    real(8)                       :: dTemp, dStepLength
!
!
    ! Loop through all solution phases expected to be stable:
    do j = 1, nSolnPhases
        k     = -iAssemblage(nElements - j + 1)
        dTemp = 0D0
        do i = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
            dMolesSpecies(i) = dStepLength * dMolesSpecies(i) + (1D0 - dStepLength) * dMolesSpeciesLast(i)
            dMolesSpecies(i) = DMAX1(dMolesSpecies(i), dTolerance(8))
            dTemp            = dTemp + dMolesSpecies(i)
        end do
        ! Compute the total number of moles of each solution phase:
        dMolesPhase(nElements-j+1) = dTemp
    end do
!
    ! Dampen the number of moles of pure condensed phases:
    do i = 1, nConPhases
        k = nSpecies - nConPhases + i
        dMolesPhase(i) = dStepLength * dMolesPhase(i) + (1D0 - dStepLength) * dMolesPhaseLast(i)
        dMolesSpecies(k) = dMolesPhase(i)
        dMolFraction(k)  = 1D0
    end do
!
    ! Dampen the element potentials:
    LOOP_Gamma: do i = 1, nElements
        if (dElementPotential(i) == 0D0) then
            dElementPotential(i) = dElementPotentialLast(i)
            cycle LOOP_Gamma
        end if
        dElementPotential(i) = dStepLength * dElementPotential(i) + (1D0 - dStepLength) * dElementPotentialLast(i)
    end do LOOP_Gamma
!
end subroutine UpdateSystemVariables
!
!
