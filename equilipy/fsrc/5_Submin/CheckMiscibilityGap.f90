!> \brief Check whether a miscibility-gap solution phase should be added.
!!
!! \details Starts Subminimization from the lowest endmember driving-force
!! extremums, ranked against the current elemental-potential plane. This keeps
!! Equilipy's targeted miscibility-gap check and does not sweep every endmember
!! the way Thermochimica-style exhaustive checks can.
!!
!-------------------------------------------------------------------------------------------------------------
!
!> \file    CheckMiscibilityGap.f90
!> \brief   Check whether a miscibility-gap solution phase should be added.
!> \author  M.H.A. Piro
!> \date    Aug. 30, 2012
!> \sa      Subminimization.f90
!> \sa      CheckSolnPhaseAdd.f90
!> \sa      CheckConvergence.f90
!
! Revisions:
! ==========
!
!   Date            Programmer          Description of change
!   ----            ----------          ---------------------
!   08/30/2012      M.H.A. Piro         Original code
!   06/25/2026      S.Y. Kwon           Ranked miscibility starts with corrected Leveling potentials.
!   06/27/2026      S.Y. Kwon           Documented the targeted ranked-start rule and bounded starts by
!                                       the number of available constituents.
!
!
! Purpose:
! ========
!
!> \details The purpose of this subroutine is to check whether a non-ideal
!! solution phase with known miscibility-gap copies should be added. The
!! routine computes each endmember's driving force relative to the current
!! elemental-potential plane, sorts those values, and only launches
!! Subminimization from the best ranked extremums. The number of starts is
!! limited by both nElements and the number of phase-local constituents.
!
!
! Required input variables:
! =========================
!
!> \param[in]   iSolnPhaseIndex  Absolute solution phase index being checked.
!
! dElementPotential             Current elemental-potential plane.
! dLevelingChemicalPotential    Leveling Gibbs value for each phase-local endmember.
! dLevelingCompositionSpecies   Leveling composition vector for each phase-local endmember.
!
!
! Output/updated variables:
! =========================
!
!> \param[out]  lAddPhase        True when the subminimized phase has negative driving force.
!
! dMolFraction                  Updated with the ranked-start composition before each Subminimization call.
! iterLastMiscGapCheck          Records the global iteration of the miscibility-gap check.
!
!
! Called subroutines/functions:
! =============================
!
! Qsort                         Ranks endmember driving-force starts from low to high.
! Subminimization               Minimizes the phase against the current elemental potentials.
! CompStoichSolnPhase           Updates effective phase stoichiometry when the phase is accepted.
!
!
! Primary callers:
! ================
!
! CheckConvergence              Checks whether a known miscibility-gap copy should enter the active set.
!
!
! Numerical assumptions:
! ======================
!
! - Starts use a dominant endmember fraction with all other constituents set to dMinMoleFraction.
! - Only the lowest MIN(nElements,nConstituents) starts are tested; the routine is not an exhaustive
!   every-endmember miscibility-gap search.
!
!-------------------------------------------------------------------------------------------------------------



subroutine CheckMiscibilityGap(iSolnPhaseIndex,lAddPhase)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer::  i, iFirst, iLast, nConstituents, iSolnPhaseIndex, nStartCandidates
    real(8)::  dMinMoleFraction, dMaxMoleFraction
    real(8), dimension(:), allocatable:: dEndmemberPotential
    integer, dimension(:), allocatable:: iEndmemberPotential
    logical::  lAddPhase

    iterLastMiscGapCheck = iterGlobal
    iFirst           = nSpeciesPhase(iSolnPhaseIndex-1) + 1
    iLast            = nSpeciesPhase(iSolnPhaseIndex)
    nConstituents    = iLast - iFirst + 1
    nStartCandidates = MIN(nElements, nConstituents)
    dMinMoleFraction = 1D-5
    dMaxMoleFraction = 1D0 - dMinMoleFraction * DFLOAT(nConstituents-1)
    dMaxMoleFraction = DMAX1(dMaxMoleFraction, 0.9D0)
    lAddPhase        = .FALSE.

    if(allocated(dEndmemberPotential)) deallocate(dEndmemberPotential)
    if(allocated(iEndmemberPotential)) deallocate(iEndmemberPotential)
    allocate(dEndmemberPotential(nConstituents),iEndmemberPotential(nConstituents))
    dEndmemberPotential = dLevelingChemicalPotential(iFirst:iLast) - &
        MATMUL(dLevelingCompositionSpecies(iFirst:iLast,:), dElementPotential)

    call Qsort(dEndmemberPotential,iEndmemberPotential,nConstituents)

    LOOP_Constituents: do i = 1, nStartCandidates

        dMolFraction(iFirst:iLast) = dMinMoleFraction
        dMolFraction(iFirst+iEndmemberPotential(i)-1) = dMaxMoleFraction

        call Subminimization(iSolnPhaseIndex, lAddPhase)

        if (lAddPhase) exit LOOP_Constituents

    end do LOOP_Constituents

    if(lAddPhase) call CompStoichSolnPhase(iSolnPhaseIndex)

    return

end subroutine CheckMiscibilityGap
