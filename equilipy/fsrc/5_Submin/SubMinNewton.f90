!-------------------------------------------------------------------------------------------------------------
!
!> \file    SubMinNewton.f90
!> \brief   Build and solve the Subminimization Newton system.
!> \author  M.H.A. Piro
!> \date    Aug. 21, 2012
!> \sa      Subminimization.f90
!> \sa      SubMinLineSearch.f90
!> \sa      SubMinTraceSpeciesControl.f90
!> \sa      CompHessianSUBG.f90
!> \sa      CompHessianRKMP.f90
!> \sa      CompHessianQKTO.f90
!
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   08/21/2012      M.H.A. Piro         Original classic Subminimization Newton routine
    !   07/20/2026      S.Y. Kwon           Added guarded analytical Newton curvature while retaining the stable ideal-log path for SUBG and SUBQ phases.
    !
    !
! Purpose:
! ========
!
!> \details This routine computes the Newton direction for the non-CEF
!! Subminimization path.  By default it builds the classic ideal-mixing
!! endmember KKT system.  If the solution model provides an analytical Hessian
!! on its native variables, this routine uses that Hessian to build the same
!! Newton update before handing the direction to SubMinLineSearch.  RKMP/QKTO
!! use ordinary one-sublattice species fractions.  The redundant Hessian is projected to the
!! sum(x)=1 simplex.  If the analytical Hessian is
!! indefinite and returns an ascent direction for the local stationarity
!! residual, the routine rebuilds the original positive diagonal Newton
!! approximation.  Ionic phases add one charge-neutrality row in the classic
!! branch.
!
!
! Required input variables:
! =========================
!
!> \param[in] iSolnPhaseIndex  Absolute solution phase index being subminimized.
!
! dMolFraction             Current active solution-variable fractions.
! dChemicalPotential       Current solution-variable chemical potentials.
! dChemicalPotentialStar   Element-potential projection for each solution variable.
! lSubMinTraceInactive     Phase-local trace-active mask.
!
!
! Output/updated variables:
! =========================
!
! dHessian                     Symmetric Newton/KKT matrix sent to DSYSV.
! dRHS                         Newton update direction on successful return.
! dPotentialVector             Updated active chemical-potential residual.
! iSubMinNewtonDSYSVInfo       DSYSV status code from the latest solve.
! dSubMinNewtonSymmetryResidual
!                              Largest pre-solve matrix symmetry mismatch.
!
!
! Called subroutines/functions:
! =============================
!
! CompHessianSUBG              Builds the analytical SUBG/SUBQ chemical-potential Hessian.
! CompHessianRKMP              Builds the analytical RKMP chemical-potential Hessian.
! CompHessianQKTO              Builds the analytical QKTO chemical-potential Hessian.
! UpdateSubMinPotentialVector  Updates the residual used to form the Newton RHS.
! DSYSV                        Solves the symmetric indefinite KKT system.
!
!
! Primary callers:
! ================
!
! Subminimization              Calls this in the non-CEF endmember/solution-variable path.
!
!
! Numerical assumptions:
! ======================
!
! - This routine is for the non-CEF path.  CEF phases that converge in
!   SubMinSiteFractionCEF do not call this routine.
! - SUBG/SUBQ phases use formal pair/quadruplet fractions and intentionally
!   remain on the classic ideal-log Newton path.  RKMP/QKTO phases use ordinary
!   one-sublattice species fractions.  Their analytical Hessians are redundant
!   on the full variable basis, so they are projected by eliminating the largest
!   current fraction as the reference variable.
! - If phase-local trace handling has temporarily removed a RKMP/QKTO
!   endmember, the routine uses the classic trace-aware branch for that reduced
!   active set.
! - Analytical Hessian directions must be descent directions for the projected
!   stationarity residual.  Otherwise the solve is globalized by rebuilding the
!   classic diagonal KKT matrix rather than sending an ascent direction to the
!   line search.
! - The matrix should be symmetric before DSYSV.  A symmetry guard records the
!   largest mismatch and turns the solve into a controlled error if the mismatch
!   exceeds roundoff-scaled tolerance.
! - Trace-inactive endmembers have zero update and do not participate in the
!   normalization or charge rows.
!
!-------------------------------------------------------------------------------------------------------------



subroutine SubMinNewton(iSolnPhaseIndex)
    USE ModuleThermoIO, ONLY: INFOThermo
    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleSubMin
!
    implicit none
!
    integer :: i, j, k, iSolnPhaseIndex, INFO, nEqn, m, LWORK
    integer :: iRef, nIndependent, nVariableOut, iInfoHessian
    integer, dimension(nVar) :: iIndependent
    real(8) :: dMatrixScale, dSymmetryTolerance, dSymmetryDiff
    real(8) :: dDirectionalDerivative
    real(8) :: dProjectedHessian
    real(8), dimension(1) :: dWorkQuery
    real(8), dimension(nVar) :: dSolutionChemicalPotential, dDirection
    real(8), dimension(nVar,nVar) :: dSolutionHessian
    real(8), dimension(:), allocatable :: WORK
    character(len=8) :: cPhaseTypeLocal
    logical :: lUseSolutionHessian
!
!
    ! Initialize variables:

    nEqn         = nVar + 1
    INFO         = 0
    LWORK        = 0
    nIndependent = 0
    nVariableOut = 0
    iInfoHessian = 0
    iRef         = 0
    iHessian     = 0
    dHessian     = 0D0
    dRHS         = -1D0
    m            = 1
    dSolutionChemicalPotential = 0D0
    dSolutionHessian = 0D0
    dDirection   = 0D0
    dSubminGibbsEst =0D0
    dDirectionalDerivative = 0D0
    iSubMinNewtonDSYSVInfo = 0
    dSubMinNewtonSymmetryResidual = 0D0
    lSubMinNewtonAnalyticalHessianAccepted = .FALSE.
    cPhaseTypeLocal = TRIM(cSolnPhaseType(iSolnPhaseIndex))
    lUseSolutionHessian = ((TRIM(cPhaseTypeLocal) == 'RKMP').OR.&
        (TRIM(cPhaseTypeLocal) == 'QKTO')).AND.&
        (iPhaseElectronID(iSolnPhaseIndex) == 0).AND.(nVar > 1)

    if (lUseSolutionHessian.AND.allocated(lSubMinTraceInactive)) then
        if (ANY(lSubMinTraceInactive(1:nVar))) lUseSolutionHessian = .FALSE.
    end if
    
!
    if (lUseSolutionHessian) then
        lSubMinNewtonAnalyticalHessianAttempted = .TRUE.
        select case (TRIM(cPhaseTypeLocal))
        case ('RKMP')
            call CompHessianRKMP(iSolnPhaseIndex, nVar, dSolutionChemicalPotential, &
                dSolutionHessian, nVariableOut, iInfoHessian)
        case ('QKTO')
            call CompHessianQKTO(iSolnPhaseIndex, nVar, dSolutionChemicalPotential, &
                dSolutionHessian, nVariableOut, iInfoHessian)
        end select

        if ((iInfoHessian /= 0).OR.(nVariableOut /= nVar)) then
            INFOThermo = 28
            iSubMinNewtonDSYSVInfo = -1000 - iInfoHessian
            dRHS = 0D0
            return
        end if

        call UpdateSubMinPotentialVector

        dRHS = 0D0
        iRef = MAXLOC(dMolFraction(iFirstSUB:iLastSUB), DIM=1)
        do i = 1, nVar
            if (i == iRef) cycle
            nIndependent = nIndependent + 1
            iIndependent(nIndependent) = i
        end do
        nEqn = nIndependent

        do j = 1, nIndependent
            i = iIndependent(j)
            dRHS(j) = dPotentialVector(i) - dPotentialVector(iRef)
            do k = 1, nIndependent
                dProjectedHessian = dSolutionHessian(i,iIndependent(k)) - dSolutionHessian(i,iRef) - &
                    dSolutionHessian(iRef,iIndependent(k)) + dSolutionHessian(iRef,iRef)
                dHessian(j,k) = dProjectedHessian
            end do
        end do
    else
        call UpdateSubMinPotentialVector
!
        ! Construct diagonal and part of the arrow head:
        ! SY: Note that dRHS is not xi but (xi-yi)
        do j = 1, nVar
            i                 = nSpeciesPhase(iSolnPhaseIndex-1) + j
            if (allocated(lSubMinTraceInactive)) then
                if (lSubMinTraceInactive(j)) then
                    dHessian(j,j) = 1D0
                    dRHS(j) = 0D0
                    cycle
                end if
            end if
            dHessian(j,j)     = 1D0 / dMolFraction(i)
            dHessian(nEqn,j)  = -1D0
            dHessian(j,nEqn)  = -1D0
            dRHS(j)           = dDrivingForce - (dChemicalPotential(i) + 1 - dChemicalPotentialStar(j))
            dRHS(nEqn)        = dRHS(nEqn) + dMolFraction(i)
        end do
!
        ! Apply an additional row/column if the phase is ionic:
        if (iPhaseElectronID(iSolnPhaseIndex) /= 0) then
            nEqn = nEqn + 1
            k    = iPhaseElectronID(iSolnPhaseIndex)
            m    = 2
!
            dRHS(nEqn) = 0D0
            do j = 1, nVar
                i                 = nSpeciesPhase(iSolnPhaseIndex-1) + j
                if (allocated(lSubMinTraceInactive)) then
                    if (lSubMinTraceInactive(j)) cycle
                end if
                dHessian(nEqn,j)  = -dStoichSpecies(i,k)
                dHessian(j,nEqn)  = dHessian(nEqn,j)
                dRHS(nEqn)        = dRHS(nEqn) + dStoichSpecies(i,k) * dMolFraction(i)
            end do
        end if
    end if
!
    ! Call the symmetric indefinite linear equation solver:
    dMatrixScale = 1D0
    dSubMinNewtonSymmetryResidual = 0D0
    do j = 1, nEqn
        do i = 1, nEqn
            dMatrixScale = MAX(dMatrixScale, ABS(dHessian(i,j)))
            if (i > j) then
                dSymmetryDiff = ABS(dHessian(i,j) - dHessian(j,i))
                if (dSymmetryDiff > dSubMinNewtonSymmetryResidual) then
                    dSubMinNewtonSymmetryResidual = dSymmetryDiff
                end if
            end if
        end do
    end do
    dSymmetryTolerance = 100D0 * 2.220446049250313D-16 * dMatrixScale

    if (dSubMinNewtonSymmetryResidual > dSymmetryTolerance) then
        INFO = -999
    else
        LWORK = -1
        dWorkQuery = 0D0
        call DSYSV('U', nEqn, 1, dHessian, SIZE(dHessian,1), iHessian, &
            dRHS, nEqn, dWorkQuery, LWORK, INFO)

        if (INFO == 0) then
            LWORK = MAX(1,INT(dWorkQuery(1)))
            allocate(WORK(LWORK))
            call DSYSV('U', nEqn, 1, dHessian, SIZE(dHessian,1), iHessian, &
                dRHS, nEqn, WORK, LWORK, INFO)
            if (allocated(WORK)) deallocate(WORK)
        end if
    end if

    if (lUseSolutionHessian.AND.(INFO == 0)) then
        dDirectionalDerivative = 0D0
        do j = 1, nIndependent
            i = iIndependent(j)
            dDirectionalDerivative = dDirectionalDerivative - &
                (dPotentialVector(i) - dPotentialVector(iRef)) * dRHS(j)
        end do

        if (dDirectionalDerivative >= 0D0) then
            lUseSolutionHessian = .FALSE.
            INFO = 0
            LWORK = 0
            nEqn = nVar + 1
            iHessian = 0
            dHessian = 0D0
            dRHS = -1D0

            do j = 1, nVar
                i = nSpeciesPhase(iSolnPhaseIndex-1) + j
                if (allocated(lSubMinTraceInactive)) then
                    if (lSubMinTraceInactive(j)) then
                        dHessian(j,j) = 1D0
                        dRHS(j) = 0D0
                        cycle
                    end if
                end if
                dHessian(j,j) = 1D0 / dMolFraction(i)
                dHessian(nEqn,j) = -1D0
                dHessian(j,nEqn) = -1D0
                dRHS(j) = dDrivingForce - (dChemicalPotential(i) + 1 - dChemicalPotentialStar(j))
                dRHS(nEqn) = dRHS(nEqn) + dMolFraction(i)
            end do

            dMatrixScale = 1D0
            dSubMinNewtonSymmetryResidual = 0D0
            do j = 1, nEqn
                do i = 1, nEqn
                    dMatrixScale = MAX(dMatrixScale, ABS(dHessian(i,j)))
                    if (i > j) then
                        dSymmetryDiff = ABS(dHessian(i,j) - dHessian(j,i))
                        if (dSymmetryDiff > dSubMinNewtonSymmetryResidual) then
                            dSubMinNewtonSymmetryResidual = dSymmetryDiff
                        end if
                    end if
                end do
            end do
            dSymmetryTolerance = 100D0 * 2.220446049250313D-16 * dMatrixScale

            if (dSubMinNewtonSymmetryResidual > dSymmetryTolerance) then
                INFO = -999
            else
                LWORK = -1
                dWorkQuery = 0D0
                call DSYSV('U', nEqn, 1, dHessian, SIZE(dHessian,1), iHessian, &
                    dRHS, nEqn, dWorkQuery, LWORK, INFO)

                if (INFO == 0) then
                    LWORK = MAX(1,INT(dWorkQuery(1)))
                    if (allocated(WORK)) deallocate(WORK)
                    allocate(WORK(LWORK))
                    call DSYSV('U', nEqn, 1, dHessian, SIZE(dHessian,1), iHessian, &
                        dRHS, nEqn, WORK, LWORK, INFO)
                    if (allocated(WORK)) deallocate(WORK)
                end if
            end if
        end if
    end if
    iSubMinNewtonDSYSVInfo = INFO
!
!
    if (INFO /= 0) then
        ! Return an error and reset dRHS:
        INFOThermo = 28
        dRHS       = 0D0
    else if (lUseSolutionHessian) then
        lSubMinNewtonAnalyticalHessianAccepted = .TRUE.
        dDirection = 0D0
        do j = 1, nIndependent
            dDirection(iIndependent(j)) = dRHS(j)
        end do
        dDirection(iRef) = -SUM(dDirection(1:nVar))
        dRHS = 0D0
        dRHS(1:nVar) = dDirection(1:nVar)
    end if
    if (allocated(WORK)) deallocate(WORK)
!
end subroutine SubMinNewton
!
!
