!> \brief Try a mass-conserving binary SUBG/SUBQ pair-exchange correction.
!!
!! \details When a single binary MQM phase is active, the three pair variables
!! have one internal exchange degree of freedom after mass balance is fixed.
!! This routine searches that coordinate for the pair-exchange stationarity
!! root and accepts it only if the standard GEM residual decreases.
subroutine TryBinarySUBGQPairExchange(lAccepted)
!
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TryBinarySUBGQPairExchange.f90
    !> \brief   Rescue no-descent binary SUBG/SUBQ Lagrangian line searches.
    !> \author  S.Y. Kwon
    !> \date    Jun. 29, 2026
    !> \sa      GEMLineSearch.f90
    !> \sa      CompFunctionNorm.f90
    !> \sa      CompChemicalPotential.f90
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/20/2026      S.Y. Kwon           Added a binary pair-exchange correction for SUBG and SUBQ solution phases.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details For a binary SUBG/SUBQ phase with three pair/quadruplet species,
    !! the mass-balance equations leave one pair-exchange coordinate.  Newton
    !! directions can be too large near trace pair fractions and the ordinary
    !! damped line search may miss the bracketed root.  This routine moves only
    !! along the null vector of the effective pair stoichiometry, solves
    !! ``sum_i v_i mu_i = 0`` by bisection, recomputes the two-element tangent
    !! plane from two independent pair equations, and keeps the trial only if
    !! ``CompFunctionNorm`` improves.
    !
    !
    ! Required input variables:
    ! =========================
    !
    ! iAssemblage              Active phase assemblage.
    ! dMolesSpecies            Current species amounts.
    ! dMolesPhase              Current active phase amounts.
    ! dElementPotential        Current elemental potentials.
    !
    !
    ! Output/updated variables:
    ! =========================
    !
    !> \param[out] lAccepted   True when the pair-exchange trial reduces the GEM norm.
    !
    ! dMolesSpecies            Updated pair species amounts when accepted.
    ! dMolFraction             Updated pair fractions when accepted.
    ! dMolesPhase              Updated active solution phase amount when accepted.
    ! dElementPotential        Updated two-element tangent plane when accepted.
    ! dGEMFunctionNorm         Recomputed residual norm after accepted or restored state.
    !
    !
    ! Called subroutines/functions:
    ! =============================
    !
    ! CompChemicalPotential    Recomputes SUBG/SUBQ pair chemical potentials.
    ! CompFunctionNorm         Evaluates the standard Lagrangian GEM residual.
    !
    !
    ! Primary callers:
    ! ================
    !
    ! GEMLineSearch            Calls this after a normal line search finds no descent.
    !
    !
    ! Numerical assumptions:
    ! ======================
    !
    ! - Applies only to one active solution phase, no condensed phases, two active
    !   elements, and exactly three SUBG/SUBQ pair species.
    ! - The pair-exchange null vector is formed from effective stoichiometry
    !   ``dStoichSpecies / iParticlesPerMole``, matching the GEM mass and plane
    !   basis.
    ! - No state is kept unless the ordinary ``CompFunctionNorm`` decreases.
    !
    !-------------------------------------------------------------------------------------------------------------

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    logical, intent(out) :: lAccepted

    integer :: i, j, iSolnPhase, iFirst, iLast, iPhaseSlot
    integer :: iIter, iPairA, iPairB, iBestA, iBestB
    real(8) :: dNormSave, dMoleFloor, dAlphaLow, dAlphaHigh, dAlphaMid
    real(8) :: dFlow, dFhigh, dFmid, dSpan, dTotalMoles, dDet, dBestDet
    real(8) :: dMuA, dMuB
    real(8), dimension(3) :: dMolesLocal, dNullVector, dTrialMoles
    real(8), dimension(3) :: dPairPotential
    real(8), dimension(3,2) :: dEffectiveStoich
    real(8), dimension(:), allocatable :: dMolesSpeciesSave, dMolFractionSave
    real(8), dimension(:), allocatable :: dMolesPhaseSave, dElementPotentialSave
    logical :: lCompEverything, lInfo

    lAccepted = .FALSE.
    lCompEverything = .FALSE.

    if (nElements /= 2) return
    if (nSolnPhases /= 1) return
    if (nConPhases /= 0) return

    iPhaseSlot = nElements
    iSolnPhase = -iAssemblage(iPhaseSlot)
    if ((iSolnPhase <= 0).OR.(iSolnPhase > nSolnPhasesSys)) return
    if ((TRIM(cSolnPhaseType(iSolnPhase)) /= 'SUBG').AND.&
        (TRIM(cSolnPhaseType(iSolnPhase)) /= 'SUBQ')) return

    iFirst = nSpeciesPhase(iSolnPhase-1) + 1
    iLast = nSpeciesPhase(iSolnPhase)
    if ((iLast - iFirst + 1) /= 3) return

    allocate(dMolesSpeciesSave(nSpecies), dMolFractionSave(nSpecies), &
        dMolesPhaseSave(nElements), dElementPotentialSave(nElements))
    dMolesSpeciesSave = dMolesSpecies
    dMolFractionSave = dMolFraction
    dMolesPhaseSave = dMolesPhase
    dElementPotentialSave = dElementPotential
    dNormSave = dGEMFunctionNorm

    do i = 1, 3
        dMolesLocal(i) = dMolesSpecies(iFirst + i - 1)
        do j = 1, 2
            dEffectiveStoich(i,j) = dStoichSpecies(iFirst + i - 1,j) / &
                DBLE(iParticlesPerMole(iFirst + i - 1))
        end do
    end do

    dNullVector(1) = dEffectiveStoich(2,1) * dEffectiveStoich(3,2) - &
        dEffectiveStoich(3,1) * dEffectiveStoich(2,2)
    dNullVector(2) = dEffectiveStoich(3,1) * dEffectiveStoich(1,2) - &
        dEffectiveStoich(1,1) * dEffectiveStoich(3,2)
    dNullVector(3) = dEffectiveStoich(1,1) * dEffectiveStoich(2,2) - &
        dEffectiveStoich(2,1) * dEffectiveStoich(1,2)

    if (MAXVAL(DABS(dNullVector)) <= dTolerance(8)) goto 900
    dNullVector = dNullVector / MAXVAL(DABS(dNullVector))

    dMoleFloor = DMAX1(dTolerance(8), 1D-300)
    dAlphaLow = -HUGE(1D0)
    dAlphaHigh = HUGE(1D0)
    do i = 1, 3
        if (dNullVector(i) > 0D0) then
            dAlphaLow = DMAX1(dAlphaLow, (dMoleFloor - dMolesLocal(i)) / dNullVector(i))
        else if (dNullVector(i) < 0D0) then
            dAlphaHigh = DMIN1(dAlphaHigh, (dMoleFloor - dMolesLocal(i)) / dNullVector(i))
        end if
    end do

    if ((dAlphaLow <= -HUGE(1D0)/2D0).OR.(dAlphaHigh >= HUGE(1D0)/2D0)) goto 900
    if (dAlphaHigh <= dAlphaLow) goto 900
    dSpan = dAlphaHigh - dAlphaLow
    dAlphaLow = dAlphaLow + 1D-12 * dSpan
    dAlphaHigh = dAlphaHigh - 1D-12 * dSpan

    call EvaluateExchangeResidual(dAlphaLow, dFlow, lInfo)
    if (.NOT.lInfo) goto 900
    call EvaluateExchangeResidual(dAlphaHigh, dFhigh, lInfo)
    if (.NOT.lInfo) goto 900
    if (dFlow == 0D0) then
        dAlphaMid = dAlphaLow
    else if (dFhigh == 0D0) then
        dAlphaMid = dAlphaHigh
    else
        if (dFlow * dFhigh > 0D0) goto 900
        dAlphaMid = 0.5D0 * (dAlphaLow + dAlphaHigh)
        do iIter = 1, 80
            call EvaluateExchangeResidual(dAlphaMid, dFmid, lInfo)
            if (.NOT.lInfo) goto 900
            if (DABS(dFmid) < 1D-12) exit
            if (dFlow * dFmid <= 0D0) then
                dAlphaHigh = dAlphaMid
                dFhigh = dFmid
            else
                dAlphaLow = dAlphaMid
                dFlow = dFmid
            end if
            dAlphaMid = 0.5D0 * (dAlphaLow + dAlphaHigh)
            if (DABS(dAlphaHigh - dAlphaLow) <= 1D-14 * DMAX1(1D0, DABS(dAlphaMid))) exit
        end do
    end if

    dTrialMoles = dMolesLocal + dAlphaMid * dNullVector
    if (MINVAL(dTrialMoles) <= dMoleFloor) goto 900
    dTotalMoles = SUM(dTrialMoles)
    if (dTotalMoles <= dMoleFloor) goto 900

    dMolesSpecies(iFirst:iLast) = dTrialMoles
    dMolesPhase(iPhaseSlot) = dTotalMoles
    dMolFraction(iFirst:iLast) = dTrialMoles / dTotalMoles
    call CompChemicalPotential(lCompEverything)
    dPairPotential = dChemicalPotential(iFirst:iLast)

    dBestDet = 0D0
    iBestA = 0
    iBestB = 0
    do iPairA = 1, 2
        do iPairB = iPairA + 1, 3
            dDet = dEffectiveStoich(iPairA,1) * dEffectiveStoich(iPairB,2) - &
                dEffectiveStoich(iPairA,2) * dEffectiveStoich(iPairB,1)
            if (DABS(dDet) > DABS(dBestDet)) then
                dBestDet = dDet
                iBestA = iPairA
                iBestB = iPairB
            end if
        end do
    end do
    if (DABS(dBestDet) <= dTolerance(8)) goto 900

    dMuA = dPairPotential(iBestA)
    dMuB = dPairPotential(iBestB)
    dElementPotential(1) = (dMuA * dEffectiveStoich(iBestB,2) - &
        dMuB * dEffectiveStoich(iBestA,2)) / dBestDet
    dElementPotential(2) = (dEffectiveStoich(iBestA,1) * dMuB - &
        dEffectiveStoich(iBestB,1) * dMuA) / dBestDet

    call CompFunctionNorm
    if (dGEMFunctionNorm < dNormSave) then
        lAccepted = .TRUE.
        goto 910
    end if

900 continue
    dMolesSpecies = dMolesSpeciesSave
    dMolFraction = dMolFractionSave
    dMolesPhase = dMolesPhaseSave
    dElementPotential = dElementPotentialSave
    call CompChemicalPotential(lCompEverything)
    call CompFunctionNorm

910 continue
    if (allocated(dMolesSpeciesSave)) deallocate(dMolesSpeciesSave)
    if (allocated(dMolFractionSave)) deallocate(dMolFractionSave)
    if (allocated(dMolesPhaseSave)) deallocate(dMolesPhaseSave)
    if (allocated(dElementPotentialSave)) deallocate(dElementPotentialSave)

    return

contains

    subroutine EvaluateExchangeResidual(dAlpha, dResidual, lValid)

        implicit none

        real(8), intent(in) :: dAlpha
        real(8), intent(out) :: dResidual
        logical, intent(out) :: lValid

        dTrialMoles = dMolesLocal + dAlpha * dNullVector
        lValid = .FALSE.
        dResidual = 0D0
        if (MINVAL(dTrialMoles) <= dMoleFloor) return
        dTotalMoles = SUM(dTrialMoles)
        if (dTotalMoles <= dMoleFloor) return

        dMolesSpecies(iFirst:iLast) = dTrialMoles
        dMolesPhase(iPhaseSlot) = dTotalMoles
        dMolFraction(iFirst:iLast) = dTrialMoles / dTotalMoles
        call CompChemicalPotential(lCompEverything)
        dResidual = SUM(dNullVector * dChemicalPotential(iFirst:iLast))
        lValid = .TRUE.

        return

    end subroutine EvaluateExchangeResidual

end subroutine TryBinarySUBGQPairExchange
!
!
