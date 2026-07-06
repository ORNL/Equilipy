!-------------------------------------------------------------------------------------------------------------
!
!> \file    GEMNewtonCEF.f90
!> \brief   Compute a mixed phase-amount/site-fraction Newton direction for active CEF phases.
!> \author  S.Y. Kwon
!> \date    Jun. 25, 2026
!> \sa      RunLagrangianGEM.f90
!> \sa      GEMLineSearchCEF.f90
!> \sa      CompGradientSUBL.f90
!> \sa      CompHessianSUBL.f90
!
! Revisions:
! ==========
!
!   Date            Programmer          Description of change
!   ----            ----------          ---------------------
!   06/25/2026      S.Y. Kwon           Added mixed-coordinate CEF Lagrangian Newton direction.
!   06/25/2026      S.Y. Kwon           Computed CEF composition and mass Jacobian terms from numeric
!                                        species stoichiometry.
!   06/28/2026      S.Y. Kwon           Generalized the guarded site-fraction path to mixed active sets
!                                        with CEF site variables and supported non-CEF fraction variables.
!   06/28/2026      S.Y. Kwon           Corrected non-CEF phase-row residuals to use molar solution Gibbs
!                                        energy on the same basis as the phase-composition plane.
!   07/01/2026      S.Y. Kwon           Documented the CEF handoff-to-site-fraction conversion contract.
!   07/01/2026      S.Y. Kwon           Excluded lower-bound CEF trace site variables from the Newton map
!                                        unless their exchange residual indicates growth.
!   07/01/2026      S.Y. Kwon           Added an optional residual Levenberg direction for CEF KKT
!                                        stagnation after the analytical site Hessian gives no descent.
!   07/02/2026      S.Y. Kwon           Added passive KKT pivot and direction-norm diagnostics for
!                                        main and residual-LM CEF solves.
!   07/02/2026      S.Y. Kwon           Reported 2x2 pivot scale from eigenvalues so near-null
!                                        Bunch-Kaufman blocks are not hidden by large off-diagonal entries.
!   07/02/2026      S.Y. Kwon           Counted SUBG/SUBQ ride-along activation in mixed CEF active sets.
!   07/03/2026      S.Y. Kwon           Stored phase stationarity residuals for passive lower-bound
!                                        complementarity diagnostics.
!   07/03/2026      S.Y. Kwon           Added a bound-active phase-amount KKT retry for
!                                        complementarity-valid raw negative CEF phase directions.
!   07/04/2026      S.Y. Kwon           Added C3-c1 primal-block Levenberg retries when KKT inertia is wrong.
!   07/05/2026      S.Y. Kwon           Seeded lower-bound CEF site constituents with growth residuals
!                                        so trace components can leave the numerical floor.
!
!
! Purpose:
! ========
!
!> \details This file contains the guarded mixed-coordinate Lagrangian path.
!! If a fixed active set contains CEF family phases (`SUBL`, `SUBLM`, or
!! `SUBOM`) together with supported non-CEF solution models, the Newton
!! variables are phase amounts, CEF site fractions for CEF phases, independent
!! phase-local fractions for non-CEF phases, and elemental potentials.
!! Leveling, PEA, and Level2Lagrange may hand off CEF constitutions as product
!! endmember fractions, but this routine first maps that representation to
!! sublattice site fractions.  The CEF site-gradient/site-Hessian terms are
!! then evaluated from `CompGradientSUBL` and `CompHessianSUBL`; supported
!! non-CEF phases use their analytical phase-fraction Hessians.
!
!
! Required input variables:
! =========================
!
! iAssemblage            Current fixed active phase assemblage.
! dMolesPhase            Current phase amounts.
! dMolFraction           Current solution endmember fractions used only to rebuild CEF sites
!                        at the start of the mixed-coordinate solve.
! dElementPotential      Current elemental potentials.
!
!
! Output/updated variables:
! =========================
!
! dGEMCEFPhaseDirection     Additive phase-amount Newton direction.
! dGEMCEFSiteDirection      Additive independent site-fraction Newton direction.
! dGEMCEFElementDirection   Additive elemental-potential Newton direction.
! iGEMCEF*                  Maps between mixed-coordinate variables and active phases/sites.
!
!
! Called subroutines/functions:
! =============================
!
! CompGradientSUBL          Computes CEF scalar Gibbs energy and site-gradient.
! CompHessianSUBL           Computes CEF site Hessian.
! DSYSV                     Solves the symmetric indefinite KKT system.
!
!
! Primary callers:
! ================
!
! RunLagrangianGEM          Dispatches to this routine for CEF fixed assemblages.
! GEMLineSearchCEF          Applies the direction computed here.
!
!
! Numerical assumptions:
! ======================
!
! - CEF phase compositions are linear in site fractions, so the composition
!   Hessian term in the Lagrangian site Hessian is zero.
! - Per-sublattice normalization is eliminated by subtracting one reference
!   constituent from each sublattice.
! - CEF phases use site-fraction variables; supported non-CEF phases use
!   independent phase-local fraction variables with one reference constituent.
! - Product endmember fractions for CEF phases are rebuilt only so the existing
!   Gibbs model routines can evaluate the current site-fraction state.  They do
!   not define the Newton residual or direction.
! - Trace site constituents satisfy complementarity.  If a constituent is in
!   the trace region and its residual points toward a smaller site fraction, it
!   is excluded from the KKT variable map; if the residual points toward growth,
!   it remains active.
! - Unsupported non-CEF phases keep the caller on the legacy GEM path.
!
!-------------------------------------------------------------------------------------------------------------



!> \brief Decide whether to use the mixed solution-coordinate Lagrangian path.
!!
!! \details The historical name is CEF-oriented, but the implementation now
!! supports SUBL/SUBLM/SUBOM site fractions and ordinary RKMP/QKTO solution
!! mole fractions.  Standalone SUBG/SUBQ active sets use the production
!! mixed-coordinate route after the reviewed C1 validation.  IDMX phases may
!! join that switched SUBG/SUBQ route as ordinary ideal species fractions.

subroutine UseCEFLagrangian(lUseCEF)
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    GEMNewtonCEF.f90
    !> \brief   Select the mixed solution-coordinate Lagrangian Newton path.
    !> \author  S.Y. Kwon
    !> \date    Jul. 2, 2026
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/02/2026      S.Y. Kwon           Activated the mixed-coordinate path for RKMP/QKTO-only active
    !                                        sets so ordinary solution fractions are solved explicitly.
    !   07/04/2026      S.Y. Kwon           Added default-off standalone SUBG/SUBQ routing switch for
    !                                        C1-b0 plumbing.
    !   07/04/2026      S.Y. Kwon           Counted SUBG/SUBQ ride-along whenever any supported
    !                                        mixed-coordinate route is active.
    !   07/04/2026      S.Y. Kwon           Allowed IDMX gas phases to join the switched SUBG/SUBQ
    !                                        mixed-coordinate route.
    !   07/04/2026      S.Y. Kwon           Promoted standalone SUBG/SUBQ mixed-coordinate routing to
    !                                        the production default.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details Decide whether the fixed-active-set Lagrangian solve should use
    !! the mixed coordinate Newton system in GEMNewtonCEF.  Despite the legacy
    !! CEF-oriented routine name, this path supports CEF site fractions and
    !! ordinary RKMP/QKTO solution fractions.  SUBG/SUBQ activate this path
    !! when the standalone switch is enabled, which is now the production
    !! default.  They may also ride along when another supported model already
    !! activated the mixed-coordinate path.  IDMX is switch-gated and may only
    !! join this path when a SUBG/SUBQ phase is present, so ordinary gas-only
    !! problems do not silently change route.
    !
    !
    ! Required input variables:
    ! =========================
    !
    ! nSolnPhases                 Number of active solution phases.
    ! iAssemblage                 Active phase assemblage.
    ! cSolnPhaseType              Solution model type for each active solution.
    ! iPhaseSublattice            Sublattice metadata used to identify CEF phases.
    !
    !
    ! Output/updated variables:
    ! =========================
    !
    !> \param[out] lUseCEF        True when GEMNewtonCEF/GEMLineSearchCEF should
    !!                             solve the current active set.
    !
    !
    ! Called subroutines/functions:
    ! =============================
    !
    ! None.
    !
    !
    ! Primary callers:
    ! ================
    !
    ! RunLagrangianGEM            Chooses GEMNewtonCEF/GEMLineSearchCEF when this returns true.
    !
    !
    ! Numerical assumptions:
    ! ======================
    !
    ! - RKMP/QKTO use ordinary one-sublattice mole fractions as Newton
    !   variables.
    ! - SUBG/SUBQ pair/quadruplet variables use the mixed-coordinate path when
    !   they ride along with another supported active model or the standalone
    !   switch is enabled.
    ! - IDMX phases use ordinary ideal species fractions.  They are admitted
    !   only under the same SUBG/SUBQ switch and only when SUBG/SUBQ is active
    !   in the same fixed assemblage.
    !
    !-------------------------------------------------------------------------------------------------------------

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    logical, intent(out) :: lUseCEF

    integer :: l, k, nSUBGQRideAlong, nIDMXRideAlong
    character(len=8) :: cTypeLocal
    logical :: lAnyCEF, lAnySupportedNonCEF, lStandaloneSUBGQActive

    lUseCEF = .FALSE.
    lAnyCEF = .FALSE.
    lAnySupportedNonCEF = .FALSE.
    lStandaloneSUBGQActive = .FALSE.
    nSUBGQRideAlong = 0
    nIDMXRideAlong = 0
    if (.NOT.lGEMCEFSiteLagrangianEnabled) return
    if (nSolnPhases <= 0) return
    if (.NOT.allocated(iPhaseSublattice)) return

    do l = 1, nSolnPhases
        k = -iAssemblage(nElements - l + 1)
        if (k <= 0) return
        if (k > SIZE(iPhaseSublattice)) return

        cTypeLocal = TRIM(cSolnPhaseType(k))
        if (((TRIM(cTypeLocal) == 'SUBL').OR.&
             (TRIM(cTypeLocal) == 'SUBLM').OR.&
             (TRIM(cTypeLocal) == 'SUBOM')).AND.(iPhaseSublattice(k) > 0)) then
            lAnyCEF = .TRUE.
        else if ((TRIM(cTypeLocal) == 'RKMP').OR.&
                 (TRIM(cTypeLocal) == 'QKTO')) then
            lAnySupportedNonCEF = .TRUE.
        else if ((TRIM(cTypeLocal) == 'SUBG').OR.&
                 (TRIM(cTypeLocal) == 'SUBQ')) then
            nSUBGQRideAlong = nSUBGQRideAlong + 1
            cycle
        else if (TRIM(cTypeLocal) == 'IDMX') then
            if (.NOT.lSUBQStandaloneEnabled) return
            nIDMXRideAlong = nIDMXRideAlong + 1
            cycle
        else
            return
        end if
    end do

    if ((nIDMXRideAlong > 0).AND.(.NOT.(lSUBQStandaloneEnabled.AND.(nSUBGQRideAlong > 0)))) return
    lStandaloneSUBGQActive = lSUBQStandaloneEnabled.AND.(nSUBGQRideAlong > 0)
    lUseCEF = lAnyCEF .OR. lAnySupportedNonCEF .OR. lStandaloneSUBGQActive
    if (lUseCEF.AND.(nSUBGQRideAlong > 0)) iGEMSUBGQRideAlongUsed = nSUBGQRideAlong

    return

end subroutine UseCEFLagrangian



subroutine GEMNewtonCEF(INFOLocal)

    USE ModuleThermo
    USE ModuleThermoIO, ONLY: INFOThermo
    USE ModuleGEMSolver
    USE ModuleSubMin, ONLY: dSubMinHenrianTraceThreshold

    implicit none

    integer, intent(out) :: INFOLocal

    integer :: i, j, k, l, p, q, s, c, r, e
    integer :: iPhaseVar, iSiteVar, iPhaseID, iSlot, iSpecies
    integer :: iFirstLocal, iLastLocal, iSpeciesRef, nVariableOut
    integer :: nVarLocal, nKKT, nSiteCapacity, nSpeciesCapacity
    integer :: nSiteOut, nSpeciesOut, iInfo, iPhaseRow, iSiteRow, iLambdaRow
    integer :: INFO_SY, INFO_BASE, nPrimalLocal, nConstraintLocal, iRegSteps
    integer, allocatable :: IPIV(:)
    integer, allocatable :: IPIVBase(:)
    real(8) :: dPhaseAmount, dScalarGibbs
    real(8) :: dScalarEnthalpy, dScalarEntropy, dScalarHeatCapacity
    real(8) :: dPhaseResidual, dSiteResidual, dMatrixScale, dSymDiff
    real(8), allocatable :: A(:,:), B(:), ABase(:,:), BBase(:)
    real(8), allocatable :: AFactorBase(:,:), BSolutionBase(:)
    real(8), allocatable :: dSiteGradient(:), dSiteGradientH(:)
    real(8), allocatable :: dSiteGradientS(:), dSiteGradientCp(:)
    real(8), allocatable :: dSiteHessian(:,:), dAmountHessian(:,:), dJacobian(:,:)
    real(8), allocatable :: dPhaseComp(:), dSiteCompDeriv(:)
    real(8), allocatable :: dSiteResidualByVar(:)
    logical :: lCEFPhaseLocal, lInitialInertiaWrong, lRegularizationAccepted
    logical :: lInertiaRegularized

    INFOLocal = 0
    INFO_SY = 0
    nSiteCapacity = nMaxSublatticeSys * nMaxConstituentSys
    nSpeciesCapacity = MAX(1,nSpecies)

    call BuildCEFLagrangianMaps

    if ((nGEMCEFPhaseVariables <= 0).OR.(nGEMCEFSiteVariables < 0)) then
        INFOLocal = 1
        return
    end if

    nVarLocal = nGEMCEFPhaseVariables + nGEMCEFSiteVariables + nElements
    nKKT = nVarLocal

    allocate(A(nKKT,nKKT), B(nKKT), IPIV(nKKT))
    allocate(ABase(nKKT,nKKT), BBase(nKKT), AFactorBase(nKKT,nKKT), &
        BSolutionBase(nKKT), IPIVBase(nKKT))
    allocate(dSiteGradient(nSiteCapacity), dSiteGradientH(nSiteCapacity), &
        dSiteGradientS(nSiteCapacity), dSiteGradientCp(nSiteCapacity))
    allocate(dSiteHessian(nSiteCapacity,nSiteCapacity), &
        dAmountHessian(nSpeciesCapacity,nSpeciesCapacity), &
        dJacobian(nSpeciesCapacity,nSpeciesCapacity))
    allocate(dPhaseComp(nElements), dSiteCompDeriv(nElements))
    allocate(dSiteResidualByVar(MAX(1,nGEMCEFSiteVariables)))

    A = 0D0
    B = 0D0
    IPIV = 0
    ABase = 0D0
    BBase = 0D0
    AFactorBase = 0D0
    BSolutionBase = 0D0
    IPIVBase = 0
    dSiteResidualByVar = 0D0
    dGEMCEFPhaseDirection = 0D0
    dGEMCEFSiteDirection = 0D0
    dGEMCEFElementDirection = 0D0
    dUpdateVar = 0D0
    iGEMNewtonSolver = 0
    iGEMNewtonDSYSVInfo = 0
    iGEMNewtonKKTSize = nKKT
    iGEMNewtonPivot1x1Count = 0
    iGEMNewtonPivot2x2Count = 0
    iGEMNewtonPivotPositiveCount = 0
    iGEMNewtonPivotNegativeCount = 0
    iGEMNewtonPivotZeroCount = 0
    dGEMNewtonSymmetryResidual = 0D0
    dGEMNewtonMinPivotScale = 0D0
    dGEMNewtonMaxPivotScale = 0D0
    dGEMNewtonDirectionNorm = 0D0
    lInertiaRegularized = .FALSE.

    if (lGEMCEFResidualLMDirection) then
        call BuildCEFResidualLMDirection(INFOLocal)
        call CleanupGEMNewtonCEF
        return
    end if

    ! Phase-amount stationarity rows and lambda coupling.
    do p = 1, nGEMCEFPhaseVariables
        iPhaseRow = p
        iSlot = iGEMCEFPhaseSlot(p)
        dPhaseComp = 0D0

        if (iGEMCEFPhaseSoln(p) > 0) then
            k = iGEMCEFPhaseSoln(p)
            if (IsCEFPhase(k)) then
                iPhaseID = iPhaseSublattice(k)
                call SetCEFMolFractionFromSite(k, iPhaseID, dGEMCEFPhaseSiteLast(p,:,:))
                call ComputeCEFPhaseComposition(k, iPhaseID, dGEMCEFPhaseSiteLast(p,:,:), dPhaseComp)
                call CompGradientSUBL(k, nSiteCapacity, dSiteGradient, dScalarGibbs, &
                    dSiteGradientH, dSiteGradientS, dSiteGradientCp, &
                    dScalarEnthalpy, dScalarEntropy, dScalarHeatCapacity, nSiteOut, iInfo)
                if (iInfo /= 0) then
                    INFOLocal = iInfo
                    exit
                end if
                dPhaseResidual = dScalarGibbs - SUM(dElementPotential(1:nElements) * dPhaseComp(1:nElements))
            else
                call ComputeNonCEFPhaseComposition(k, dPhaseComp)
                if (dMolesPhase(iSlot) <= dTolerance(8)) then
                    INFOLocal = 1
                    exit
                end if
                dPhaseResidual = dGibbsSolnPhase(k) / dMolesPhase(iSlot) - &
                    SUM(dElementPotential(1:nElements) * dPhaseComp(1:nElements))
            end if
        else
            iSpecies = iGEMCEFPhaseSpecies(p)
            dPhaseComp = dStoichSpecies(iSpecies,1:nElements)
            dPhaseResidual = dStdGibbsEnergy(iSpecies) - &
                SUM(dElementPotential(1:nElements) * dPhaseComp(1:nElements))
        end if

        if (allocated(dGEMCEFPhaseResidual)) dGEMCEFPhaseResidual(p) = dPhaseResidual
        B(iPhaseRow) = -dPhaseResidual
        do e = 1, nElements
            iLambdaRow = nGEMCEFPhaseVariables + nGEMCEFSiteVariables + e
            A(iPhaseRow,iLambdaRow) = -dPhaseComp(e)
            A(iLambdaRow,iPhaseRow) = -dPhaseComp(e)
        end do
    end do

    if (INFOLocal /= 0) then
        call CleanupGEMNewtonCEF
        return
    end if

    ! Site-fraction stationarity rows and CEF Hessian blocks.
    do p = 1, nGEMCEFPhaseVariables
        k = iGEMCEFPhaseSoln(p)
        if (k <= 0) cycle
        iSlot = iGEMCEFPhaseSlot(p)
        dPhaseAmount = dMolesPhase(iSlot)
        lCEFPhaseLocal = IsCEFPhase(k)

        if (lCEFPhaseLocal) then
            iPhaseID = iPhaseSublattice(k)
            call SetCEFMolFractionFromSite(k, iPhaseID, dGEMCEFPhaseSiteLast(p,:,:))
            call CompGradientSUBL(k, nSiteCapacity, dSiteGradient, dScalarGibbs, &
                dSiteGradientH, dSiteGradientS, dSiteGradientCp, &
                dScalarEnthalpy, dScalarEntropy, dScalarHeatCapacity, nSiteOut, iInfo)
            if (iInfo /= 0) then
                INFOLocal = iInfo
                exit
            end if
            call CompHessianSUBL(k, nSiteCapacity, nSpeciesCapacity, dSiteHessian, &
                dAmountHessian, dJacobian, nSiteOut, nSpeciesOut, iInfo)
            if (iInfo /= 0) then
                INFOLocal = iInfo
                exit
            end if
        else
            call CompNonCEFHessian(k, nSpeciesCapacity, dAmountHessian, nSpeciesOut, iInfo)
            if (iInfo /= 0) then
                INFOLocal = iInfo
                exit
            end if
        end if

        do iSiteVar = 1, nGEMCEFSiteVariables
            if (iGEMCEFVarPhaseVar(iSiteVar) /= p) cycle
            iSiteRow = nGEMCEFPhaseVariables + iSiteVar
            c = iGEMCEFVarCon(iSiteVar)
            r = iGEMCEFVarRef(iSiteVar)

            if (lCEFPhaseLocal) then
                iPhaseID = iGEMCEFVarPhaseID(iSiteVar)
                s = iGEMCEFVarSub(iSiteVar)
                call ComputeCEFCompositionDerivative(k, iPhaseID, s, c, r, &
                    dGEMCEFPhaseSiteLast(p,:,:), dSiteCompDeriv)
                dSiteResidual = dSiteGradient(LocalCEFIndex(iPhaseID,s,c)) - &
                    dSiteGradient(LocalCEFIndex(iPhaseID,s,r)) - &
                    SUM(dElementPotential(1:nElements) * dSiteCompDeriv(1:nElements))
            else
                call ComputeNonCEFCompositionDerivative(k, c, r, dSiteCompDeriv)
                iFirstLocal = nSpeciesPhase(k-1) + 1
                iSpecies = iFirstLocal + c - 1
                iSpeciesRef = iFirstLocal + r - 1
                dSiteResidual = (dChemicalPotential(iSpecies) - &
                    SUM(dElementPotential(1:nElements) * &
                    dStoichSpecies(iSpecies,1:nElements) / DBLE(iParticlesPerMole(iSpecies)))) - &
                    (dChemicalPotential(iSpeciesRef) - &
                    SUM(dElementPotential(1:nElements) * &
                    dStoichSpecies(iSpeciesRef,1:nElements) / DBLE(iParticlesPerMole(iSpeciesRef))))
            end if
            dSiteResidualByVar(iSiteVar) = dSiteResidual

            B(iSiteRow) = -dPhaseAmount * dSiteResidual
            do e = 1, nElements
                iLambdaRow = nGEMCEFPhaseVariables + nGEMCEFSiteVariables + e
                A(iSiteRow,iLambdaRow) = -dPhaseAmount * dSiteCompDeriv(e)
                A(iLambdaRow,iSiteRow) = -dPhaseAmount * dSiteCompDeriv(e)
            end do

            do j = 1, nGEMCEFSiteVariables
                if (iGEMCEFVarPhaseVar(j) /= p) cycle
                q = nGEMCEFPhaseVariables + j
                if (lCEFPhaseLocal) then
                    A(iSiteRow,q) = dPhaseAmount * ProjectCEFHessian(dSiteHessian, iPhaseID, &
                        s, c, r, iGEMCEFVarSub(j), iGEMCEFVarCon(j), iGEMCEFVarRef(j))
                else
                    A(iSiteRow,q) = dPhaseAmount * ProjectNonCEFHessian(dAmountHessian, &
                        c, r, iGEMCEFVarCon(j), iGEMCEFVarRef(j))
                end if
                A(q,iSiteRow) = A(iSiteRow,q)
            end do
        end do
    end do

    if (INFOLocal /= 0) then
        call CleanupGEMNewtonCEF
        return
    end if

    ! Mass-balance residual rows.
    do e = 1, nElements
        iLambdaRow = nGEMCEFPhaseVariables + nGEMCEFSiteVariables + e
        B(iLambdaRow) = -dMolesElement(e)
    end do

    do p = 1, nGEMCEFPhaseVariables
        iSlot = iGEMCEFPhaseSlot(p)
        if (iGEMCEFPhaseSoln(p) > 0) then
            k = iGEMCEFPhaseSoln(p)
            if (IsCEFPhase(k)) then
                iPhaseID = iPhaseSublattice(k)
                call ComputeCEFPhaseComposition(k, iPhaseID, dGEMCEFPhaseSiteLast(p,:,:), dPhaseComp)
            else
                call ComputeNonCEFPhaseComposition(k, dPhaseComp)
            end if
        else
            iSpecies = iGEMCEFPhaseSpecies(p)
            dPhaseComp = dStoichSpecies(iSpecies,1:nElements)
        end if

        do e = 1, nElements
            iLambdaRow = nGEMCEFPhaseVariables + nGEMCEFSiteVariables + e
            B(iLambdaRow) = B(iLambdaRow) + dMolesPhase(iSlot) * dPhaseComp(e)
        end do
    end do

    ! Site contribution to the mass-balance Jacobian.
    do iSiteVar = 1, nGEMCEFSiteVariables
        p = iGEMCEFVarPhaseVar(iSiteVar)
        iSlot = iGEMCEFPhaseSlot(p)
        iPhaseID = iGEMCEFVarPhaseID(iSiteVar)
        c = iGEMCEFVarCon(iSiteVar)
        r = iGEMCEFVarRef(iSiteVar)
        if (iPhaseID > 0) then
            s = iGEMCEFVarSub(iSiteVar)
            call ComputeCEFCompositionDerivative(iGEMCEFVarSolnPhase(iSiteVar), iPhaseID, s, c, r, &
                dGEMCEFPhaseSiteLast(p,:,:), dSiteCompDeriv)
        else
            call ComputeNonCEFCompositionDerivative(iGEMCEFVarSolnPhase(iSiteVar), c, r, dSiteCompDeriv)
        end if
        do e = 1, nElements
            iLambdaRow = nGEMCEFPhaseVariables + nGEMCEFSiteVariables + e
            iSiteRow = nGEMCEFPhaseVariables + iSiteVar
            A(iLambdaRow,iSiteRow) = -dMolesPhase(iSlot) * dSiteCompDeriv(e)
            A(iSiteRow,iLambdaRow) = A(iLambdaRow,iSiteRow)
        end do
    end do

    if (lGEMCEFBndPhaseActive) then
        call ApplyCEFBoundPhaseConstraint(A, B, nKKT, INFOLocal)
        if (INFOLocal /= 0) then
            call CleanupGEMNewtonCEF
            return
        end if
    end if

    do i = 1, nKKT
        do j = i + 1, nKKT
            dSymDiff = DABS(A(i,j) - A(j,i))
            if (dSymDiff > dGEMNewtonSymmetryResidual) dGEMNewtonSymmetryResidual = dSymDiff
        end do
    end do
    dMatrixScale = MAX(1D0,MAXVAL(DABS(A)))
    if (dGEMNewtonSymmetryResidual > 100D0 * 2.220446049250313D-16 * dMatrixScale) then
        INFOLocal = -999
        iGEMNewtonDSYSVInfo = INFOLocal
        call CleanupGEMNewtonCEF
        return
    end if

    nPrimalLocal = nGEMCEFPhaseVariables + nGEMCEFSiteVariables
    nConstraintLocal = nElements
    ABase = A
    BBase = B

    call SolveCEFKKTSystem(A, B, IPIV, nKKT, INFO_SY)
    INFO_BASE = INFO_SY
    AFactorBase = A
    BSolutionBase = B
    IPIVBase = IPIV

    INFOLocal = INFO_SY
    iGEMNewtonDSYSVInfo = INFO_SY
    call StoreCEFPivotDiagnostics(nKKT, A, IPIV)
    lInitialInertiaWrong = CEFKKTInertiaWrong(INFO_SY, nPrimalLocal, nConstraintLocal)
    if (lGEMCEFInertiaRegularizationActive.AND.lInitialInertiaWrong) then
        call TryCEFPrimalRegularization(ABase, BBase, A, B, IPIV, nKKT, &
            nPrimalLocal, nConstraintLocal, INFO_SY, lRegularizationAccepted, iRegSteps)
        if (lRegularizationAccepted) then
            INFOLocal = INFO_SY
            iGEMNewtonDSYSVInfo = INFO_SY
            lInertiaRegularized = .TRUE.
        else
            A = AFactorBase
            B = BSolutionBase
            IPIV = IPIVBase
            INFO_SY = INFO_BASE
            INFOLocal = INFO_BASE
            iGEMNewtonDSYSVInfo = INFO_BASE
            call StoreCEFPivotDiagnostics(nKKT, A, IPIV)
        end if
    end if

    if (INFO_SY == 0) then
        do i = 1, nKKT
            if (B(i) /= B(i)) INFOLocal = 1
        end do
    end if

    if (INFOLocal == 0) then
        dGEMNewtonDirectionNorm = DSQRT(SUM(B(1:nKKT)**2))
        dGEMCEFPhaseDirection(1:nGEMCEFPhaseVariables) = B(1:nGEMCEFPhaseVariables)
        if (nGEMCEFSiteVariables > 0) then
            dGEMCEFSiteDirection(1:nGEMCEFSiteVariables) = &
                B(nGEMCEFPhaseVariables+1:nGEMCEFPhaseVariables+nGEMCEFSiteVariables)
        end if
        dGEMCEFElementDirection(1:nElements) = &
            B(nGEMCEFPhaseVariables+nGEMCEFSiteVariables+1:nKKT)
        lGEMCEFSiteDirectionActive = .TRUE.
        lRevertSystem = .FALSE.
        if (lInertiaRegularized) then
            iGEMNewtonSolver = 7
        else if (lGEMCEFBndPhaseActive) then
            iGEMNewtonSolver = 6
        else
            iGEMNewtonSolver = 4
        end if
    else
        lGEMCEFSiteDirectionActive = .FALSE.
        lRevertSystem = .TRUE.
        dGEMCEFPhaseDirection = 0D0
        dGEMCEFSiteDirection = 0D0
        dGEMCEFElementDirection = 0D0
    end if

    call CleanupGEMNewtonCEF

    return

contains

    subroutine SolveCEFKKTSystem(AMatrix, BRhs, IPIVLocal, nLocal, iInfoOut)

        implicit none

        integer, intent(in) :: nLocal
        real(8), intent(inout) :: AMatrix(nLocal,nLocal), BRhs(nLocal)
        integer, intent(inout) :: IPIVLocal(nLocal)
        integer, intent(out) :: iInfoOut

        integer :: LWORKLocal
        real(8) :: dWorkQueryLocal(1)
        real(8), allocatable :: WORKLocal(:)

        iInfoOut = 0
        LWORKLocal = -1
        dWorkQueryLocal = 0D0
        call DSYSV('U', nLocal, 1, AMatrix, nLocal, IPIVLocal, BRhs, nLocal, &
            dWorkQueryLocal, LWORKLocal, iInfoOut)
        if (iInfoOut == 0) then
            LWORKLocal = MAX(1,INT(dWorkQueryLocal(1)))
            allocate(WORKLocal(LWORKLocal))
            call DSYSV('U', nLocal, 1, AMatrix, nLocal, IPIVLocal, BRhs, nLocal, &
                WORKLocal, LWORKLocal, iInfoOut)
            if (allocated(WORKLocal)) deallocate(WORKLocal)
        end if

        return

    end subroutine SolveCEFKKTSystem

    logical function CEFKKTInertiaWrong(iInfoIn, nPrimalIn, nConstraintIn)

        implicit none

        integer, intent(in) :: iInfoIn, nPrimalIn, nConstraintIn

        CEFKKTInertiaWrong = (iInfoIn /= 0).OR.&
            (iGEMNewtonPivotZeroCount /= 0).OR.&
            (iGEMNewtonPivotPositiveCount /= nPrimalIn).OR.&
            (iGEMNewtonPivotNegativeCount /= nConstraintIn)

        return

    end function CEFKKTInertiaWrong

    subroutine TryCEFPrimalRegularization(ABaseIn, BBaseIn, AMatrix, BRhs, IPIVLocal, &
        nLocal, nPrimalIn, nConstraintIn, iInfoOut, lAcceptedOut, iStepsOut)

        implicit none

        integer, intent(in) :: nLocal, nPrimalIn, nConstraintIn
        real(8), intent(in) :: ABaseIn(nLocal,nLocal), BBaseIn(nLocal)
        real(8), intent(inout) :: AMatrix(nLocal,nLocal), BRhs(nLocal)
        integer, intent(inout) :: IPIVLocal(nLocal)
        integer, intent(out) :: iInfoOut, iStepsOut
        logical, intent(out) :: lAcceptedOut

        integer, parameter :: iMaxRegularizationStep = 10
        integer :: iStepLocal
        real(8) :: dMuLocal

        lAcceptedOut = .FALSE.
        iInfoOut = 0
        iStepsOut = 0
        dMuLocal = InitialCEFRegularizationMu(ABaseIn, nLocal, nPrimalIn)
        iGEMInertiaRegularizationAttemptedUsed = iGEMInertiaRegularizationAttemptedUsed + 1
        iGEMInertiaRegularizationAttemptedTotal = iGEMInertiaRegularizationAttemptedTotal + 1

        do iStepLocal = 1, iMaxRegularizationStep
            AMatrix = ABaseIn
            BRhs = BBaseIn
            IPIVLocal = 0
            call ApplyCEFPrimalShift(AMatrix, nLocal, nPrimalIn, dMuLocal)
            call SolveCEFKKTSystem(AMatrix, BRhs, IPIVLocal, nLocal, iInfoOut)
            call StoreCEFPivotDiagnostics(nLocal, AMatrix, IPIVLocal)
            iStepsOut = iStepLocal

            if (.NOT.CEFKKTInertiaWrong(iInfoOut, nPrimalIn, nConstraintIn)) then
                lAcceptedOut = .TRUE.
                iGEMInertiaRegularizationAcceptedUsed = iGEMInertiaRegularizationAcceptedUsed + 1
                iGEMInertiaRegularizationAcceptedTotal = iGEMInertiaRegularizationAcceptedTotal + 1
                exit
            end if

            dMuLocal = 10D0 * dMuLocal
        end do

        iGEMInertiaRegularizationStepLast = iStepsOut
        iGEMInertiaRegularizationStepTotal = iGEMInertiaRegularizationStepTotal + iStepsOut
        if (.NOT.lAcceptedOut) then
            iGEMInertiaRegularizationFailedUsed = iGEMInertiaRegularizationFailedUsed + 1
            iGEMInertiaRegularizationFailedTotal = iGEMInertiaRegularizationFailedTotal + 1
        end if

        return

    end subroutine TryCEFPrimalRegularization

    real(8) function InitialCEFRegularizationMu(AMatrix, nLocal, nPrimalIn)

        implicit none

        integer, intent(in) :: nLocal, nPrimalIn
        real(8), intent(in) :: AMatrix(nLocal,nLocal)

        real(8) :: dHScaleLocal, dPivotScaleLocal

        dHScaleLocal = 1D0
        if (nPrimalIn > 0) dHScaleLocal = DMAX1(1D0, &
            MAXVAL(DABS(AMatrix(1:nPrimalIn,1:nPrimalIn))))
        dPivotScaleLocal = dGEMNewtonMinPivotScale
        if ((dPivotScaleLocal /= dPivotScaleLocal).OR.(dPivotScaleLocal < 0D0)) then
            dPivotScaleLocal = 0D0
        end if
        InitialCEFRegularizationMu = DMAX1(1D-14, &
            DMAX1(1D-12 * dHScaleLocal, 10D0 * dPivotScaleLocal))

        return

    end function InitialCEFRegularizationMu

    subroutine ApplyCEFPrimalShift(AMatrix, nLocal, nPrimalIn, dMuIn)

        implicit none

        integer, intent(in) :: nLocal, nPrimalIn
        real(8), intent(inout) :: AMatrix(nLocal,nLocal)
        real(8), intent(in) :: dMuIn

        integer :: iLocal, iBoundVar

        iBoundVar = BoundPhaseVariableIndex()
        do iLocal = 1, nPrimalIn
            if (iLocal == iBoundVar) cycle
            AMatrix(iLocal,iLocal) = AMatrix(iLocal,iLocal) + dMuIn
        end do

        return

    end subroutine ApplyCEFPrimalShift

    integer function BoundPhaseVariableIndex()

        implicit none

        integer :: iRowLocal

        BoundPhaseVariableIndex = 0
        if (.NOT.lGEMCEFBndPhaseActive) return
        do iRowLocal = 1, nGEMCEFPhaseVariables
            if (iGEMCEFPhaseSlot(iRowLocal) == iGEMCEFBndPhaseSlot) then
                BoundPhaseVariableIndex = iRowLocal
                exit
            end if
        end do

        return

    end function BoundPhaseVariableIndex

    subroutine ApplyCEFBoundPhaseConstraint(AMatrix, BRhs, nLocal, iInfoOut)

        implicit none

        integer, intent(in) :: nLocal
        real(8), intent(inout) :: AMatrix(nLocal,nLocal), BRhs(nLocal)
        integer, intent(out) :: iInfoOut

        integer :: iBoundVar, iRowLocal
        real(8) :: dFixedStep

        iInfoOut = 0
        iBoundVar = BoundPhaseVariableIndex()

        if (iBoundVar <= 0) then
            iInfoOut = 1
            return
        end if

        dFixedStep = dGEMCEFBndPhaseStep
        do iRowLocal = 1, nLocal
            if (iRowLocal == iBoundVar) cycle
            BRhs(iRowLocal) = BRhs(iRowLocal) - AMatrix(iRowLocal,iBoundVar) * dFixedStep
        end do

        AMatrix(:,iBoundVar) = 0D0
        AMatrix(iBoundVar,:) = 0D0
        AMatrix(iBoundVar,iBoundVar) = 1D0
        BRhs(iBoundVar) = dFixedStep

        return

    end subroutine ApplyCEFBoundPhaseConstraint

    subroutine StoreCEFPivotDiagnostics(nLocal, AFactor, IPIVLocal)

        implicit none

        integer, intent(in) :: nLocal
        integer, intent(in) :: IPIVLocal(nLocal)
        real(8), intent(in) :: AFactor(nLocal,nLocal)

        integer :: iLocal
        real(8) :: dTolLocal, dPivotLocal
        real(8) :: dBlockA, dBlockB, dBlockC, dTraceLocal, dDiscLocal, dEig1, dEig2

        if (nLocal <= 0) return

        dTolLocal = 1D-12 * DMAX1(1D0,MAXVAL(DABS(AFactor)))
        iGEMNewtonPivot1x1Count = 0
        iGEMNewtonPivot2x2Count = 0
        iGEMNewtonPivotPositiveCount = 0
        iGEMNewtonPivotNegativeCount = 0
        iGEMNewtonPivotZeroCount = 0
        dGEMNewtonMinPivotScale = HUGE(1D0)
        dGEMNewtonMaxPivotScale = 0D0

        iLocal = 1
        do while (iLocal <= nLocal)
            if ((IPIVLocal(iLocal) < 0).AND.(iLocal < nLocal)) then
                iGEMNewtonPivot2x2Count = iGEMNewtonPivot2x2Count + 1
                dBlockA = AFactor(iLocal,iLocal)
                dBlockB = AFactor(iLocal,iLocal+1)
                dBlockC = AFactor(iLocal+1,iLocal+1)
                dTraceLocal = 0.5D0 * (dBlockA + dBlockC)
                dDiscLocal = DSQRT(MAX(0D0,0.25D0 * (dBlockA - dBlockC)**2 + dBlockB**2))
                dEig1 = dTraceLocal + dDiscLocal
                dEig2 = dTraceLocal - dDiscLocal
                dPivotLocal = DMIN1(DABS(dEig1),DABS(dEig2))
                call CountCEFPivotSign(dEig1, dTolLocal)
                call CountCEFPivotSign(dEig2, dTolLocal)
                iLocal = iLocal + 2
            else
                iGEMNewtonPivot1x1Count = iGEMNewtonPivot1x1Count + 1
                dPivotLocal = DABS(AFactor(iLocal,iLocal))
                call CountCEFPivotSign(AFactor(iLocal,iLocal), dTolLocal)
                iLocal = iLocal + 1
            end if

            dGEMNewtonMinPivotScale = DMIN1(dGEMNewtonMinPivotScale, dPivotLocal)
            dGEMNewtonMaxPivotScale = DMAX1(dGEMNewtonMaxPivotScale, dPivotLocal)
        end do

        if (dGEMNewtonMinPivotScale == HUGE(1D0)) dGEMNewtonMinPivotScale = 0D0

        return

    end subroutine StoreCEFPivotDiagnostics

    subroutine CountCEFPivotSign(dValue, dTolLocal)

        implicit none

        real(8), intent(in) :: dValue, dTolLocal

        if (dValue > dTolLocal) then
            iGEMNewtonPivotPositiveCount = iGEMNewtonPivotPositiveCount + 1
        else if (dValue < -dTolLocal) then
            iGEMNewtonPivotNegativeCount = iGEMNewtonPivotNegativeCount + 1
        else
            iGEMNewtonPivotZeroCount = iGEMNewtonPivotZeroCount + 1
        end if

        return

    end subroutine CountCEFPivotSign

    subroutine BuildCEFResidualLMDirection(iInfoOut)

        implicit none

        integer, intent(out) :: iInfoOut

        integer :: iVar, jVar, iRowLocal, INFO_LM, LWORK_LM
        integer, allocatable :: IPIV_LM(:)
        real(8) :: dStepLocal, dDirectionScale, dDirectionCap, dDamping
        real(8) :: dWorkQueryLM(1)
        real(8), allocatable :: dResidual0(:), dResidualPlus(:), dResidualMinus(:)
        real(8), allocatable :: dJacobianLM(:,:), dNormalLM(:,:), dRhsLM(:), WORK_LM(:)
        real(8), allocatable :: dMolesPhaseBase(:), dMolesSpeciesBase(:)
        real(8), allocatable :: dMolFractionBase(:), dElementPotentialBase(:)
        real(8), allocatable :: dSiteBase(:,:,:)
        logical :: lPlusFeasible, lMinusFeasible

        iInfoOut = 0
        INFO_LM = 0
        dDamping = 1D-3

        allocate(dResidual0(nKKT), dResidualPlus(nKKT), dResidualMinus(nKKT))
        allocate(dJacobianLM(nKKT,nKKT), dNormalLM(nKKT,nKKT), dRhsLM(nKKT))
        allocate(IPIV_LM(nKKT))
        allocate(dMolesPhaseBase(nElements), dMolesSpeciesBase(nSpecies), &
            dMolFractionBase(nSpecies), dElementPotentialBase(nElements))
        allocate(dSiteBase(MAX(1,nGEMCEFPhaseVariables),nMaxSublatticeSys,nMaxConstituentSys))

        dMolesPhaseBase = dMolesPhase
        dMolesSpeciesBase = dMolesSpecies
        dMolFractionBase = dMolFraction
        dElementPotentialBase = dElementPotential
        dSiteBase = dGEMCEFPhaseSiteLast

        call BuildCEFResidualVector(dResidual0, iInfoOut)
        if (iInfoOut /= 0) goto 900

        dJacobianLM = 0D0
        do iVar = 1, nKKT
            dStepLocal = CEFResidualLMStep(iVar)

            call RestoreCEFResidualLMState(dMolesPhaseBase, dMolesSpeciesBase, &
                dMolFractionBase, dElementPotentialBase, dSiteBase)
            call PerturbCEFResidualLMVariable(iVar, dStepLocal, lPlusFeasible)
            if (lPlusFeasible) then
                call BuildCEFResidualVector(dResidualPlus, iInfoOut)
                if (iInfoOut /= 0) lPlusFeasible = .FALSE.
            end if

            call RestoreCEFResidualLMState(dMolesPhaseBase, dMolesSpeciesBase, &
                dMolFractionBase, dElementPotentialBase, dSiteBase)
            call PerturbCEFResidualLMVariable(iVar, -dStepLocal, lMinusFeasible)
            if (lMinusFeasible) then
                call BuildCEFResidualVector(dResidualMinus, iInfoOut)
                if (iInfoOut /= 0) lMinusFeasible = .FALSE.
            end if

            if (lPlusFeasible.AND.lMinusFeasible) then
                dJacobianLM(:,iVar) = (dResidualPlus - dResidualMinus) / (2D0 * dStepLocal)
            else if (lPlusFeasible) then
                dJacobianLM(:,iVar) = (dResidualPlus - dResidual0) / dStepLocal
            else if (lMinusFeasible) then
                dJacobianLM(:,iVar) = (dResidual0 - dResidualMinus) / dStepLocal
            else
                dJacobianLM(:,iVar) = 0D0
            end if
        end do

        call RestoreCEFResidualLMState(dMolesPhaseBase, dMolesSpeciesBase, &
            dMolFractionBase, dElementPotentialBase, dSiteBase)

        dNormalLM = MATMUL(TRANSPOSE(dJacobianLM), dJacobianLM)
        dRhsLM = -MATMUL(TRANSPOSE(dJacobianLM), dResidual0)
        do iVar = 1, nKKT
            dNormalLM(iVar,iVar) = dNormalLM(iVar,iVar) + dDamping
        end do

        LWORK_LM = -1
        dWorkQueryLM = 0D0
        call DSYSV('U', nKKT, 1, dNormalLM, nKKT, IPIV_LM, dRhsLM, nKKT, &
            dWorkQueryLM, LWORK_LM, INFO_LM)
        if (INFO_LM == 0) then
            LWORK_LM = MAX(1,INT(dWorkQueryLM(1)))
            allocate(WORK_LM(LWORK_LM))
            call DSYSV('U', nKKT, 1, dNormalLM, nKKT, IPIV_LM, dRhsLM, nKKT, &
                WORK_LM, LWORK_LM, INFO_LM)
            if (allocated(WORK_LM)) deallocate(WORK_LM)
        end if
        call StoreCEFPivotDiagnostics(nKKT, dNormalLM, IPIV_LM)

        if (INFO_LM /= 0) then
            iInfoOut = INFO_LM
            goto 900
        end if

        do iVar = 1, nKKT
            if (dRhsLM(iVar) /= dRhsLM(iVar)) then
                iInfoOut = 1
                goto 900
            end if
        end do

        dDirectionScale = 1D0
        do iVar = 1, nKKT
            dDirectionCap = CEFResidualLMDirectionCap(iVar)
            if (dDirectionCap > 0D0) then
                dDirectionScale = DMAX1(dDirectionScale, DABS(dRhsLM(iVar)) / dDirectionCap)
            end if
        end do
        dRhsLM = dRhsLM / dDirectionScale
        dGEMNewtonDirectionNorm = DSQRT(SUM(dRhsLM(1:nKKT)**2))

        dGEMCEFPhaseDirection(1:nGEMCEFPhaseVariables) = dRhsLM(1:nGEMCEFPhaseVariables)
        if (nGEMCEFSiteVariables > 0) then
            dGEMCEFSiteDirection(1:nGEMCEFSiteVariables) = &
                dRhsLM(nGEMCEFPhaseVariables+1:nGEMCEFPhaseVariables+nGEMCEFSiteVariables)
        end if
        do iRowLocal = 1, nElements
            dGEMCEFElementDirection(iRowLocal) = &
                dRhsLM(nGEMCEFPhaseVariables+nGEMCEFSiteVariables+iRowLocal)
        end do

        lGEMCEFSiteDirectionActive = .TRUE.
        lRevertSystem = .FALSE.
        iGEMNewtonSolver = 5
        iGEMNewtonDSYSVInfo = INFO_LM

    900 continue
        call RestoreCEFResidualLMState(dMolesPhaseBase, dMolesSpeciesBase, &
            dMolFractionBase, dElementPotentialBase, dSiteBase)

        if (iInfoOut /= 0) then
            lGEMCEFSiteDirectionActive = .FALSE.
            lRevertSystem = .TRUE.
            dGEMCEFPhaseDirection = 0D0
            dGEMCEFSiteDirection = 0D0
            dGEMCEFElementDirection = 0D0
            iGEMNewtonDSYSVInfo = iInfoOut
        end if

        if (allocated(dResidual0)) deallocate(dResidual0)
        if (allocated(dResidualPlus)) deallocate(dResidualPlus)
        if (allocated(dResidualMinus)) deallocate(dResidualMinus)
        if (allocated(dJacobianLM)) deallocate(dJacobianLM)
        if (allocated(dNormalLM)) deallocate(dNormalLM)
        if (allocated(dRhsLM)) deallocate(dRhsLM)
        if (allocated(IPIV_LM)) deallocate(IPIV_LM)
        if (allocated(WORK_LM)) deallocate(WORK_LM)
        if (allocated(dMolesPhaseBase)) deallocate(dMolesPhaseBase)
        if (allocated(dMolesSpeciesBase)) deallocate(dMolesSpeciesBase)
        if (allocated(dMolFractionBase)) deallocate(dMolFractionBase)
        if (allocated(dElementPotentialBase)) deallocate(dElementPotentialBase)
        if (allocated(dSiteBase)) deallocate(dSiteBase)

        return

    end subroutine BuildCEFResidualLMDirection

    subroutine BuildCEFResidualVector(dResidualOut, iInfoOut)

        implicit none

        real(8), intent(out) :: dResidualOut(nKKT)
        integer, intent(out) :: iInfoOut

        integer :: iRowLocal, iPhaseVar, iSiteVar, iSlotLocal, iSolnPhaseLocal
        integer :: iSpeciesLocal, iPhaseIDLocal, iSubLocal, iConLocal, iRefLocal
        integer :: iElemLocal, iInfoLocal, nSiteOutLocal
        real(8) :: dScalarGibbsLocal, dScalarHLocal, dScalarSLocal, dScalarCpLocal
        real(8), allocatable :: dGradientLocal(:), dGradientHLocal(:)
        real(8), allocatable :: dGradientSLocal(:), dGradientCpLocal(:)
        real(8), allocatable :: dPhaseCompositionLocal(:), dSiteDerivativeLocal(:)
        logical :: lCompEverythingLocal

        iInfoOut = 0
        dResidualOut = 0D0
        lCompEverythingLocal = .FALSE.
        call CompChemicalPotential(lCompEverythingLocal)

        allocate(dGradientLocal(nSiteCapacity), dGradientHLocal(nSiteCapacity), &
            dGradientSLocal(nSiteCapacity), dGradientCpLocal(nSiteCapacity))
        allocate(dPhaseCompositionLocal(nElements), dSiteDerivativeLocal(nElements))

        iRowLocal = 0

        do iPhaseVar = 1, nGEMCEFPhaseVariables
            iRowLocal = iRowLocal + 1
            iSlotLocal = iGEMCEFPhaseSlot(iPhaseVar)
            dPhaseCompositionLocal = 0D0

            if (iGEMCEFPhaseSoln(iPhaseVar) > 0) then
                iSolnPhaseLocal = iGEMCEFPhaseSoln(iPhaseVar)
                if (IsCEFPhase(iSolnPhaseLocal)) then
                    iPhaseIDLocal = iPhaseSublattice(iSolnPhaseLocal)
                    call SetCEFMolFractionFromSite(iSolnPhaseLocal, iPhaseIDLocal, &
                        dGEMCEFPhaseSiteLast(iPhaseVar,:,:))
                    call ComputeCEFPhaseComposition(iSolnPhaseLocal, iPhaseIDLocal, &
                        dGEMCEFPhaseSiteLast(iPhaseVar,:,:), dPhaseCompositionLocal)
                    call CompGradientSUBL(iSolnPhaseLocal, nSiteCapacity, dGradientLocal, &
                        dScalarGibbsLocal, dGradientHLocal, dGradientSLocal, dGradientCpLocal, &
                        dScalarHLocal, dScalarSLocal, dScalarCpLocal, nSiteOutLocal, iInfoLocal)
                    if (iInfoLocal /= 0) then
                        iInfoOut = iInfoLocal
                        exit
                    end if
                    dResidualOut(iRowLocal) = dScalarGibbsLocal - &
                        SUM(dElementPotential(1:nElements) * dPhaseCompositionLocal(1:nElements))
                else
                    call ComputeNonCEFPhaseComposition(iSolnPhaseLocal, dPhaseCompositionLocal)
                    if (dMolesPhase(iSlotLocal) <= dTolerance(8)) then
                        iInfoOut = 1
                        exit
                    end if
                    dResidualOut(iRowLocal) = dGibbsSolnPhase(iSolnPhaseLocal) / &
                        dMolesPhase(iSlotLocal) - &
                        SUM(dElementPotential(1:nElements) * dPhaseCompositionLocal(1:nElements))
                end if
            else
                iSpeciesLocal = iGEMCEFPhaseSpecies(iPhaseVar)
                dPhaseCompositionLocal = dStoichSpecies(iSpeciesLocal,1:nElements)
                dResidualOut(iRowLocal) = dStdGibbsEnergy(iSpeciesLocal) - &
                    SUM(dElementPotential(1:nElements) * dPhaseCompositionLocal(1:nElements))
            end if
        end do

        if (iInfoOut /= 0) goto 900

        do iSiteVar = 1, nGEMCEFSiteVariables
            iRowLocal = nGEMCEFPhaseVariables + iSiteVar
            iSolnPhaseLocal = iGEMCEFVarSolnPhase(iSiteVar)
            iConLocal = iGEMCEFVarCon(iSiteVar)
            iRefLocal = iGEMCEFVarRef(iSiteVar)
            if (iGEMCEFVarPhaseID(iSiteVar) > 0) then
                iPhaseIDLocal = iGEMCEFVarPhaseID(iSiteVar)
                iSubLocal = iGEMCEFVarSub(iSiteVar)
                call SetCEFMolFractionFromSite(iSolnPhaseLocal, iPhaseIDLocal, &
                    dGEMCEFPhaseSiteLast(iGEMCEFVarPhaseVar(iSiteVar),:,:))
                call CompGradientSUBL(iSolnPhaseLocal, nSiteCapacity, dGradientLocal, &
                    dScalarGibbsLocal, dGradientHLocal, dGradientSLocal, dGradientCpLocal, &
                    dScalarHLocal, dScalarSLocal, dScalarCpLocal, nSiteOutLocal, iInfoLocal)
                if (iInfoLocal /= 0) then
                    iInfoOut = iInfoLocal
                    exit
                end if
                call ComputeCEFCompositionDerivative(iSolnPhaseLocal, iPhaseIDLocal, iSubLocal, &
                    iConLocal, iRefLocal, &
                    dGEMCEFPhaseSiteLast(iGEMCEFVarPhaseVar(iSiteVar),:,:), dSiteDerivativeLocal)
                dResidualOut(iRowLocal) = dGradientLocal(LocalCEFIndex(iPhaseIDLocal,iSubLocal,iConLocal)) - &
                    dGradientLocal(LocalCEFIndex(iPhaseIDLocal,iSubLocal,iRefLocal)) - &
                    SUM(dElementPotential(1:nElements) * dSiteDerivativeLocal(1:nElements))
            else
                call ComputeNonCEFCompositionDerivative(iSolnPhaseLocal, iConLocal, iRefLocal, &
                    dSiteDerivativeLocal)
                iFirstLocal = nSpeciesPhase(iSolnPhaseLocal-1) + 1
                iSpeciesLocal = iFirstLocal + iConLocal - 1
                iSpeciesRef = iFirstLocal + iRefLocal - 1
                dResidualOut(iRowLocal) = (dChemicalPotential(iSpeciesLocal) - &
                    SUM(dElementPotential(1:nElements) * &
                    dStoichSpecies(iSpeciesLocal,1:nElements) / DBLE(iParticlesPerMole(iSpeciesLocal)))) - &
                    (dChemicalPotential(iSpeciesRef) - &
                    SUM(dElementPotential(1:nElements) * &
                    dStoichSpecies(iSpeciesRef,1:nElements) / DBLE(iParticlesPerMole(iSpeciesRef))))
            end if
        end do

        if (iInfoOut /= 0) goto 900

        do iElemLocal = 1, nElements
            iRowLocal = nGEMCEFPhaseVariables + nGEMCEFSiteVariables + iElemLocal
            dResidualOut(iRowLocal) = -dMolesElement(iElemLocal)
        end do

        do iPhaseVar = 1, nGEMCEFPhaseVariables
            iSlotLocal = iGEMCEFPhaseSlot(iPhaseVar)
            dPhaseCompositionLocal = 0D0
            if (iGEMCEFPhaseSoln(iPhaseVar) > 0) then
                iSolnPhaseLocal = iGEMCEFPhaseSoln(iPhaseVar)
                if (IsCEFPhase(iSolnPhaseLocal)) then
                    iPhaseIDLocal = iPhaseSublattice(iSolnPhaseLocal)
                    call ComputeCEFPhaseComposition(iSolnPhaseLocal, iPhaseIDLocal, &
                        dGEMCEFPhaseSiteLast(iPhaseVar,:,:), dPhaseCompositionLocal)
                else
                    call ComputeNonCEFPhaseComposition(iSolnPhaseLocal, dPhaseCompositionLocal)
                end if
            else
                iSpeciesLocal = iGEMCEFPhaseSpecies(iPhaseVar)
                dPhaseCompositionLocal = dStoichSpecies(iSpeciesLocal,1:nElements)
            end if
            do iElemLocal = 1, nElements
                iRowLocal = nGEMCEFPhaseVariables + nGEMCEFSiteVariables + iElemLocal
                dResidualOut(iRowLocal) = dResidualOut(iRowLocal) + &
                    dMolesPhase(iSlotLocal) * dPhaseCompositionLocal(iElemLocal)
            end do
        end do

    900 continue
        if (allocated(dGradientLocal)) deallocate(dGradientLocal)
        if (allocated(dGradientHLocal)) deallocate(dGradientHLocal)
        if (allocated(dGradientSLocal)) deallocate(dGradientSLocal)
        if (allocated(dGradientCpLocal)) deallocate(dGradientCpLocal)
        if (allocated(dPhaseCompositionLocal)) deallocate(dPhaseCompositionLocal)
        if (allocated(dSiteDerivativeLocal)) deallocate(dSiteDerivativeLocal)

        return

    end subroutine BuildCEFResidualVector

    subroutine RestoreCEFResidualLMState(dMolesPhaseBase, dMolesSpeciesBase, &
        dMolFractionBase, dElementPotentialBase, dSiteBase)

        implicit none

        real(8), intent(in) :: dMolesPhaseBase(nElements)
        real(8), intent(in) :: dMolesSpeciesBase(nSpecies)
        real(8), intent(in) :: dMolFractionBase(nSpecies)
        real(8), intent(in) :: dElementPotentialBase(nElements)
        real(8), intent(in) :: dSiteBase(:,:,:)

        logical :: lCompEverythingLocal

        lCompEverythingLocal = .FALSE.
        dMolesPhase = dMolesPhaseBase
        dMolesSpecies = dMolesSpeciesBase
        dMolFraction = dMolFractionBase
        dElementPotential = dElementPotentialBase
        dGEMCEFPhaseSiteLast = dSiteBase
        call CompChemicalPotential(lCompEverythingLocal)

        return

    end subroutine RestoreCEFResidualLMState

    subroutine PerturbCEFResidualLMVariable(iVariableIn, dStepIn, lFeasibleOut)

        implicit none

        integer, intent(in) :: iVariableIn
        real(8), intent(in) :: dStepIn
        logical, intent(out) :: lFeasibleOut

        integer :: iSiteVariableLocal, iElementLocal, iPhaseVariableLocal
        integer :: iSlotLocal, iSolnPhaseLocal, iPhaseIDLocal, iSubLocal, iConLocal, iRefLocal
        integer :: iFirstPerturb, iLastPerturb, iSpeciesLocal, iSpeciesRefLocal
        real(8) :: dTrialValue, dSiteFloorLocal
        logical :: lCompEverythingLocal

        lFeasibleOut = .TRUE.
        lCompEverythingLocal = .FALSE.
        dSiteFloorLocal = DMAX1(dTraceSpeciesRemoveFraction, 1D-30)

        if (iVariableIn <= nGEMCEFPhaseVariables) then
            iPhaseVariableLocal = iVariableIn
            iSlotLocal = iGEMCEFPhaseSlot(iPhaseVariableLocal)
            dTrialValue = dMolesPhase(iSlotLocal) + dStepIn
            if (dTrialValue <= dTolerance(8)) then
                lFeasibleOut = .FALSE.
                return
            end if
            dMolesPhase(iSlotLocal) = dTrialValue
            call RebuildCEFResidualLMPhaseMoles(iPhaseVariableLocal)
        else if (iVariableIn <= nGEMCEFPhaseVariables + nGEMCEFSiteVariables) then
            iSiteVariableLocal = iVariableIn - nGEMCEFPhaseVariables
            iSolnPhaseLocal = iGEMCEFVarSolnPhase(iSiteVariableLocal)
            iPhaseIDLocal = iGEMCEFVarPhaseID(iSiteVariableLocal)
            iConLocal = iGEMCEFVarCon(iSiteVariableLocal)
            iRefLocal = iGEMCEFVarRef(iSiteVariableLocal)
            if (iPhaseIDLocal > 0) then
                iSubLocal = iGEMCEFVarSub(iSiteVariableLocal)
                iPhaseVariableLocal = iGEMCEFVarPhaseVar(iSiteVariableLocal)
                dGEMCEFPhaseSiteLast(iPhaseVariableLocal,iSubLocal,iConLocal) = &
                    dGEMCEFPhaseSiteLast(iPhaseVariableLocal,iSubLocal,iConLocal) + dStepIn
                dGEMCEFPhaseSiteLast(iPhaseVariableLocal,iSubLocal,iRefLocal) = &
                    dGEMCEFPhaseSiteLast(iPhaseVariableLocal,iSubLocal,iRefLocal) - dStepIn
                if ((dGEMCEFPhaseSiteLast(iPhaseVariableLocal,iSubLocal,iConLocal) <= dSiteFloorLocal).OR.&
                    (dGEMCEFPhaseSiteLast(iPhaseVariableLocal,iSubLocal,iRefLocal) <= dSiteFloorLocal)) then
                    lFeasibleOut = .FALSE.
                    return
                end if
                call RebuildCEFResidualLMPhaseMoles(iPhaseVariableLocal)
            else
                iFirstPerturb = nSpeciesPhase(iSolnPhaseLocal-1) + 1
                iLastPerturb = nSpeciesPhase(iSolnPhaseLocal)
                iSpeciesLocal = iFirstPerturb + iConLocal - 1
                iSpeciesRefLocal = iFirstPerturb + iRefLocal - 1
                dMolFraction(iSpeciesLocal) = dMolFraction(iSpeciesLocal) + dStepIn
                dMolFraction(iSpeciesRefLocal) = dMolFraction(iSpeciesRefLocal) - dStepIn
                if ((dMolFraction(iSpeciesLocal) <= dSiteFloorLocal).OR.&
                    (dMolFraction(iSpeciesRefLocal) <= dSiteFloorLocal)) then
                    lFeasibleOut = .FALSE.
                    return
                end if
                iSlotLocal = iGEMCEFPhaseSlot(iGEMCEFVarPhaseVar(iSiteVariableLocal))
                dMolesSpecies(iFirstPerturb:iLastPerturb) = &
                    dMolesPhase(iSlotLocal) * dMolFraction(iFirstPerturb:iLastPerturb)
            end if
        else
            iElementLocal = iVariableIn - nGEMCEFPhaseVariables - nGEMCEFSiteVariables
            dElementPotential(iElementLocal) = dElementPotential(iElementLocal) + dStepIn
        end if

        call CompChemicalPotential(lCompEverythingLocal)

        return

    end subroutine PerturbCEFResidualLMVariable

    subroutine RebuildCEFResidualLMPhaseMoles(iPhaseVariableIn)

        implicit none

        integer, intent(in) :: iPhaseVariableIn

        integer :: iSolnPhaseLocal, iPhaseIDLocal, iSlotLocal
        integer :: iFirstRebuild, iLastRebuild

        iSolnPhaseLocal = iGEMCEFPhaseSoln(iPhaseVariableIn)
        if (iSolnPhaseLocal <= 0) then
            dMolesSpecies(iGEMCEFPhaseSpecies(iPhaseVariableIn)) = &
                dMolesPhase(iGEMCEFPhaseSlot(iPhaseVariableIn))
            return
        end if

        iSlotLocal = iGEMCEFPhaseSlot(iPhaseVariableIn)
        if (IsCEFPhase(iSolnPhaseLocal)) then
            iPhaseIDLocal = iPhaseSublattice(iSolnPhaseLocal)
            call SetCEFMolesFromSiteLocal(iSolnPhaseLocal, iPhaseIDLocal, &
                dGEMCEFPhaseSiteLast(iPhaseVariableIn,:,:), dMolesPhase(iSlotLocal))
        else
            iFirstRebuild = nSpeciesPhase(iSolnPhaseLocal-1) + 1
            iLastRebuild = nSpeciesPhase(iSolnPhaseLocal)
            dMolesSpecies(iFirstRebuild:iLastRebuild) = &
                dMolesPhase(iSlotLocal) * dMolFraction(iFirstRebuild:iLastRebuild)
        end if

        return

    end subroutine RebuildCEFResidualLMPhaseMoles

    subroutine SetCEFMolesFromSiteLocal(iSolnPhaseIndexIn, iPhaseIDIn, dSiteIn, dPhaseAmountIn)

        implicit none

        integer, intent(in) :: iSolnPhaseIndexIn, iPhaseIDIn
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)
        real(8), intent(in) :: dPhaseAmountIn

        integer :: iFirstLocalSet, iLastLocalSet, iLocalSpeciesSet, iLocalSet
        integer :: iSubSet, iConSet
        real(8) :: dProductSet, dProductSumSet

        iFirstLocalSet = nSpeciesPhase(iSolnPhaseIndexIn-1) + 1
        iLastLocalSet = nSpeciesPhase(iSolnPhaseIndexIn)
        dProductSumSet = 0D0
        do iLocalSpeciesSet = iFirstLocalSet, iLastLocalSet
            iLocalSet = iLocalSpeciesSet - iFirstLocalSet + 1
            dProductSet = 1D0
            do iSubSet = 1, nSublatticePhase(iPhaseIDIn)
                iConSet = iConstituentSublattice(iPhaseIDIn,iSubSet,iLocalSet)
                if (iConSet > 0) dProductSet = dProductSet * DMAX1(dSiteIn(iSubSet,iConSet), 1D-75)
            end do
            dMolFraction(iLocalSpeciesSet) = dProductSet
            dProductSumSet = dProductSumSet + dProductSet
        end do

        if (dProductSumSet > 0D0) then
            do iLocalSpeciesSet = iFirstLocalSet, iLastLocalSet
                dMolFraction(iLocalSpeciesSet) = dMolFraction(iLocalSpeciesSet) / dProductSumSet
                dMolesSpecies(iLocalSpeciesSet) = dPhaseAmountIn * dMolFraction(iLocalSpeciesSet)
            end do
        end if

        return

    end subroutine SetCEFMolesFromSiteLocal

    real(8) function CEFResidualLMStep(iVariableIn)

        implicit none

        integer, intent(in) :: iVariableIn

        integer :: iSlotLocal, iElementLocal

        if (iVariableIn <= nGEMCEFPhaseVariables) then
            iSlotLocal = iGEMCEFPhaseSlot(iVariableIn)
            CEFResidualLMStep = 1D-6 * DMAX1(1D0, DABS(dMolesPhase(iSlotLocal)))
        else if (iVariableIn <= nGEMCEFPhaseVariables + nGEMCEFSiteVariables) then
            CEFResidualLMStep = 1D-6
        else
            iElementLocal = iVariableIn - nGEMCEFPhaseVariables - nGEMCEFSiteVariables
            CEFResidualLMStep = 1D-6 * DMAX1(1D0, DABS(dElementPotential(iElementLocal)))
        end if

        return

    end function CEFResidualLMStep

    real(8) function CEFResidualLMDirectionCap(iVariableIn)

        implicit none

        integer, intent(in) :: iVariableIn

        integer :: iSlotLocal, iElementLocal

        if (iVariableIn <= nGEMCEFPhaseVariables) then
            iSlotLocal = iGEMCEFPhaseSlot(iVariableIn)
            CEFResidualLMDirectionCap = 2D-1 * DMAX1(1D0, DABS(dMolesPhase(iSlotLocal)))
        else if (iVariableIn <= nGEMCEFPhaseVariables + nGEMCEFSiteVariables) then
            CEFResidualLMDirectionCap = 1D-2
        else
            iElementLocal = iVariableIn - nGEMCEFPhaseVariables - nGEMCEFSiteVariables
            CEFResidualLMDirectionCap = 5D-1 * DMAX1(1D0, DABS(dElementPotential(iElementLocal)))
        end if

        return

    end function CEFResidualLMDirectionCap

    subroutine BuildCEFLagrangianMaps

        implicit none

        integer :: nPhaseMax, nSiteMax, iCursor, iRef, cSeed
        integer :: iFirst, iLast, iLocal, nLocal, nMapSiteCapacity, nMapSiteOut, iMapInfo
        real(8) :: dRefValue, dSiteTraceThreshold, dMapResidual, dTraceSeed
        real(8), parameter :: dPassiveTraceSiteSeed = 1D-8
        real(8) :: dMapScalarGibbs, dMapScalarH, dMapScalarS, dMapScalarCp
        real(8), allocatable :: dMapSiteGradient(:), dMapSiteGradientH(:)
        real(8), allocatable :: dMapSiteGradientS(:), dMapSiteGradientCp(:)
        real(8), allocatable :: dMapSiteCompDeriv(:)
        logical :: lIncludeSiteVariable, lActiveSetContainsOrderDisorderMapped
        logical :: lTraceElementEdgeSystem

        nPhaseMax = nConPhases + nSolnPhases
        nSiteMax = MAX(1,nSpecies + nSolnPhases * nMaxSublatticeSys * nMaxConstituentSys)
        nMapSiteCapacity = nMaxSublatticeSys * nMaxConstituentSys
        dSiteTraceThreshold = DMAX1(DMAX1(dTraceSpeciesRemoveFraction, 1D-30), dSubMinHenrianTraceThreshold)
        lActiveSetContainsOrderDisorderMapped = ActiveSetContainsOrderDisorderMapped()
        lTraceElementEdgeSystem = TraceElementEdgeSystem()

        if (allocated(iGEMCEFPhaseSlot)) deallocate(iGEMCEFPhaseSlot)
        if (allocated(iGEMCEFPhaseSoln)) deallocate(iGEMCEFPhaseSoln)
        if (allocated(iGEMCEFPhaseSpecies)) deallocate(iGEMCEFPhaseSpecies)
        if (allocated(iGEMCEFVarPhaseVar)) deallocate(iGEMCEFVarPhaseVar)
        if (allocated(iGEMCEFVarSolnPhase)) deallocate(iGEMCEFVarSolnPhase)
        if (allocated(iGEMCEFVarPhaseID)) deallocate(iGEMCEFVarPhaseID)
        if (allocated(iGEMCEFVarSub)) deallocate(iGEMCEFVarSub)
        if (allocated(iGEMCEFVarCon)) deallocate(iGEMCEFVarCon)
        if (allocated(iGEMCEFVarRef)) deallocate(iGEMCEFVarRef)
        if (allocated(dGEMCEFPhaseDirection)) deallocate(dGEMCEFPhaseDirection)
        if (allocated(dGEMCEFSiteDirection)) deallocate(dGEMCEFSiteDirection)
        if (allocated(dGEMCEFElementDirection)) deallocate(dGEMCEFElementDirection)
        if (allocated(dGEMCEFPhaseResidual)) deallocate(dGEMCEFPhaseResidual)
        if (allocated(dGEMCEFSiteLast)) deallocate(dGEMCEFSiteLast)
        if (allocated(dGEMCEFPhaseSiteLast)) deallocate(dGEMCEFPhaseSiteLast)

        allocate(iGEMCEFPhaseSlot(nPhaseMax), iGEMCEFPhaseSoln(nPhaseMax), &
            iGEMCEFPhaseSpecies(nPhaseMax))
        allocate(iGEMCEFVarPhaseVar(nSiteMax), iGEMCEFVarSolnPhase(nSiteMax), &
            iGEMCEFVarPhaseID(nSiteMax), iGEMCEFVarSub(nSiteMax), &
            iGEMCEFVarCon(nSiteMax), iGEMCEFVarRef(nSiteMax))
        allocate(dGEMCEFPhaseDirection(nPhaseMax), dGEMCEFSiteDirection(nSiteMax), &
            dGEMCEFElementDirection(nElements))
        allocate(dGEMCEFPhaseResidual(nPhaseMax))
        allocate(dGEMCEFSiteLast(MAX(1,nCountSublattice),nMaxSublatticeSys,nMaxConstituentSys))
        allocate(dGEMCEFPhaseSiteLast(nPhaseMax,nMaxSublatticeSys,nMaxConstituentSys))
        allocate(dMapSiteGradient(nMapSiteCapacity), dMapSiteGradientH(nMapSiteCapacity), &
            dMapSiteGradientS(nMapSiteCapacity), dMapSiteGradientCp(nMapSiteCapacity), &
            dMapSiteCompDeriv(nElements))

        iGEMCEFPhaseSlot = 0
        iGEMCEFPhaseSoln = 0
        iGEMCEFPhaseSpecies = 0
        iGEMCEFVarPhaseVar = 0
        iGEMCEFVarSolnPhase = 0
        iGEMCEFVarPhaseID = 0
        iGEMCEFVarSub = 0
        iGEMCEFVarCon = 0
        iGEMCEFVarRef = 0
        dGEMCEFPhaseDirection = 0D0
        dGEMCEFSiteDirection = 0D0
        dGEMCEFElementDirection = 0D0
        dGEMCEFPhaseResidual = 0D0
        dGEMCEFSiteLast = 0D0
        dGEMCEFPhaseSiteLast = 0D0
        nGEMCEFPhaseVariables = 0
        nGEMCEFSiteVariables = 0
        lGEMCEFSiteDirectionActive = .FALSE.

        do i = 1, nConPhases
            nGEMCEFPhaseVariables = nGEMCEFPhaseVariables + 1
            iGEMCEFPhaseSlot(nGEMCEFPhaseVariables) = i
            iGEMCEFPhaseSpecies(nGEMCEFPhaseVariables) = iAssemblage(i)
        end do

        do l = 1, nSolnPhases
            iSlot = nElements - l + 1
            k = -iAssemblage(iSlot)
            if (allocated(iActiveSlotThermoPhase)) then
                if ((iSlot > 0).AND.(iSlot <= SIZE(iActiveSlotThermoPhase))) then
                    if (iActiveSlotThermoPhase(iSlot) > 0) k = iActiveSlotThermoPhase(iSlot)
                end if
            end if
            nGEMCEFPhaseVariables = nGEMCEFPhaseVariables + 1
            iGEMCEFPhaseSlot(nGEMCEFPhaseVariables) = iSlot
            iGEMCEFPhaseSoln(nGEMCEFPhaseVariables) = k

            if (IsCEFPhase(k)) then
                iPhaseID = iPhaseSublattice(k)
                if (allocated(dActiveSlotSiteFraction)) then
                    if ((iSlot > 0).AND.(iSlot <= SIZE(dActiveSlotSiteFraction,1))) then
                        if (SUM(dActiveSlotSiteFraction(iSlot,:,:)) > 0D0) then
                            dGEMCEFPhaseSiteLast(nGEMCEFPhaseVariables,:,:) = dActiveSlotSiteFraction(iSlot,:,:)
                        else
                            call BuildCEFSiteFromMolFraction(k, iPhaseID, &
                                dGEMCEFPhaseSiteLast(nGEMCEFPhaseVariables,:,:))
                        end if
                    else
                        call BuildCEFSiteFromMolFraction(k, iPhaseID, &
                            dGEMCEFPhaseSiteLast(nGEMCEFPhaseVariables,:,:))
                    end if
                else
                    call BuildCEFSiteFromMolFraction(k, iPhaseID, &
                        dGEMCEFPhaseSiteLast(nGEMCEFPhaseVariables,:,:))
                end if
                call NormalizeCEFSite(iPhaseID, dGEMCEFPhaseSiteLast(nGEMCEFPhaseVariables,:,:))
                dGEMCEFSiteLast(iPhaseID,:,:) = dGEMCEFPhaseSiteLast(nGEMCEFPhaseVariables,:,:)
                call SetCEFMolFractionFromSite(k, iPhaseID, dGEMCEFPhaseSiteLast(nGEMCEFPhaseVariables,:,:))
                call CompGradientSUBL(k, nMapSiteCapacity, dMapSiteGradient, dMapScalarGibbs, &
                    dMapSiteGradientH, dMapSiteGradientS, dMapSiteGradientCp, &
                    dMapScalarH, dMapScalarS, dMapScalarCp, nMapSiteOut, iMapInfo)

                do s = 1, nSublatticePhase(iPhaseID)
                    iRef = 1
                    dRefValue = dGEMCEFPhaseSiteLast(nGEMCEFPhaseVariables,s,1)
                    do c = 2, nConstituentSublattice(iPhaseID,s)
                        if (dGEMCEFPhaseSiteLast(nGEMCEFPhaseVariables,s,c) > dRefValue) then
                            iRef = c
                            dRefValue = dGEMCEFPhaseSiteLast(nGEMCEFPhaseVariables,s,c)
                        end if
                    end do

                    do c = 1, nConstituentSublattice(iPhaseID,s)
                        if (c == iRef) cycle
                        lIncludeSiteVariable = .TRUE.
                        if ((iMapInfo == 0).AND.&
                            (dGEMCEFPhaseSiteLast(nGEMCEFPhaseVariables,s,c) <= dSiteTraceThreshold)) then
                            call ComputeCEFCompositionDerivative(k, iPhaseID, s, c, iRef, &
                                dGEMCEFPhaseSiteLast(nGEMCEFPhaseVariables,:,:), dMapSiteCompDeriv)
                            dMapResidual = dMapSiteGradient(LocalCEFIndex(iPhaseID,s,c)) - &
                                dMapSiteGradient(LocalCEFIndex(iPhaseID,s,iRef)) - &
                                SUM(dElementPotential(1:nElements) * dMapSiteCompDeriv(1:nElements))
                            if (lTraceElementEdgeSystem.AND.&
                                (.NOT.lActiveSetContainsOrderDisorderMapped).AND.&
                                (dMapResidual < -dTraceSpeciesResidualTolerance)) then
                                dTraceSeed = DMIN1(3D-4, &
                                    0.5D0 * dGEMCEFPhaseSiteLast(nGEMCEFPhaseVariables,s,iRef))
                                if (dTraceSeed > dGEMCEFPhaseSiteLast(nGEMCEFPhaseVariables,s,c)) then
                                    dGEMCEFPhaseSiteLast(nGEMCEFPhaseVariables,s,iRef) = &
                                        dGEMCEFPhaseSiteLast(nGEMCEFPhaseVariables,s,iRef) - &
                                        (dTraceSeed - dGEMCEFPhaseSiteLast(nGEMCEFPhaseVariables,s,c))
                                    dGEMCEFPhaseSiteLast(nGEMCEFPhaseVariables,s,c) = dTraceSeed
                                    do cSeed = 1, nConstituentSublattice(iPhaseID,s)
                                        if ((cSeed == iRef).OR.(cSeed == c)) cycle
                                        if (dGEMCEFPhaseSiteLast(nGEMCEFPhaseVariables,s,cSeed) < &
                                            dPassiveTraceSiteSeed) then
                                            dGEMCEFPhaseSiteLast(nGEMCEFPhaseVariables,s,iRef) = &
                                                dGEMCEFPhaseSiteLast(nGEMCEFPhaseVariables,s,iRef) - &
                                                (dPassiveTraceSiteSeed - &
                                                dGEMCEFPhaseSiteLast(nGEMCEFPhaseVariables,s,cSeed))
                                            dGEMCEFPhaseSiteLast(nGEMCEFPhaseVariables,s,cSeed) = &
                                                dPassiveTraceSiteSeed
                                        end if
                                    end do
                                    call NormalizeCEFSite(iPhaseID, &
                                        dGEMCEFPhaseSiteLast(nGEMCEFPhaseVariables,:,:))
                                    dGEMCEFSiteLast(iPhaseID,:,:) = &
                                        dGEMCEFPhaseSiteLast(nGEMCEFPhaseVariables,:,:)
                                    call SetCEFMolFractionFromSite(k, iPhaseID, &
                                        dGEMCEFPhaseSiteLast(nGEMCEFPhaseVariables,:,:))
                                    call CompGradientSUBL(k, nMapSiteCapacity, dMapSiteGradient, &
                                        dMapScalarGibbs, dMapSiteGradientH, dMapSiteGradientS, &
                                        dMapSiteGradientCp, dMapScalarH, dMapScalarS, dMapScalarCp, &
                                        nMapSiteOut, iMapInfo)
                                end if
                            end if
                            lIncludeSiteVariable = dMapResidual < -dTraceSpeciesResidualTolerance
                        end if
                        if (.NOT.lIncludeSiteVariable) cycle
                        nGEMCEFSiteVariables = nGEMCEFSiteVariables + 1
                        iGEMCEFVarPhaseVar(nGEMCEFSiteVariables) = nGEMCEFPhaseVariables
                        iGEMCEFVarSolnPhase(nGEMCEFSiteVariables) = k
                        iGEMCEFVarPhaseID(nGEMCEFSiteVariables) = iPhaseID
                        iGEMCEFVarSub(nGEMCEFSiteVariables) = s
                        iGEMCEFVarCon(nGEMCEFSiteVariables) = c
                        iGEMCEFVarRef(nGEMCEFSiteVariables) = iRef
                    end do
                end do
            else
                iFirst = nSpeciesPhase(k-1) + 1
                iLast = nSpeciesPhase(k)
                nLocal = iLast - iFirst + 1
                if (nLocal <= 1) cycle
                iRef = 1
                dRefValue = dMolFraction(iFirst)
                do iLocal = 2, nLocal
                    if (dMolFraction(iFirst+iLocal-1) > dRefValue) then
                        iRef = iLocal
                        dRefValue = dMolFraction(iFirst+iLocal-1)
                    end if
                end do

                do iLocal = 1, nLocal
                    if (iLocal == iRef) cycle
                    nGEMCEFSiteVariables = nGEMCEFSiteVariables + 1
                    iGEMCEFVarPhaseVar(nGEMCEFSiteVariables) = nGEMCEFPhaseVariables
                    iGEMCEFVarSolnPhase(nGEMCEFSiteVariables) = k
                    iGEMCEFVarPhaseID(nGEMCEFSiteVariables) = 0
                    iGEMCEFVarSub(nGEMCEFSiteVariables) = 0
                    iGEMCEFVarCon(nGEMCEFSiteVariables) = iLocal
                    iGEMCEFVarRef(nGEMCEFSiteVariables) = iRef
                end do
            end if
        end do

        if (allocated(dMapSiteGradient)) deallocate(dMapSiteGradient)
        if (allocated(dMapSiteGradientH)) deallocate(dMapSiteGradientH)
        if (allocated(dMapSiteGradientS)) deallocate(dMapSiteGradientS)
        if (allocated(dMapSiteGradientCp)) deallocate(dMapSiteGradientCp)
        if (allocated(dMapSiteCompDeriv)) deallocate(dMapSiteCompDeriv)

        return

    end subroutine BuildCEFLagrangianMaps

    logical function PhaseIsOrderDisorderMapped(iPhaseIn)

        implicit none

        integer, intent(in) :: iPhaseIn
        integer :: iOrderedPhase

        PhaseIsOrderDisorderMapped = .FALSE.
        if (.NOT.allocated(iDisorderedPhase)) return
        if (iPhaseIn <= 0) return

        if (iPhaseIn <= SIZE(iDisorderedPhase)) then
            if (iDisorderedPhase(iPhaseIn) > 0) then
                PhaseIsOrderDisorderMapped = .TRUE.
                return
            end if
        end if

        do iOrderedPhase = 1, MIN(nSolnPhasesSys, SIZE(iDisorderedPhase))
            if (iDisorderedPhase(iOrderedPhase) == iPhaseIn) then
                PhaseIsOrderDisorderMapped = .TRUE.
                return
            end if
        end do

        return

    end function PhaseIsOrderDisorderMapped

    logical function ActiveSetContainsOrderDisorderMapped()

        implicit none

        integer :: iSlotLocal, iPhaseLocal

        ActiveSetContainsOrderDisorderMapped = .FALSE.
        do iSlotLocal = 1, nElements
            if (iAssemblage(iSlotLocal) >= 0) cycle
            iPhaseLocal = -iAssemblage(iSlotLocal)
            if (allocated(iActiveSlotThermoPhase)) then
                if ((iSlotLocal > 0).AND.(iSlotLocal <= SIZE(iActiveSlotThermoPhase))) then
                    if (iActiveSlotThermoPhase(iSlotLocal) > 0) iPhaseLocal = &
                        iActiveSlotThermoPhase(iSlotLocal)
                end if
            end if
            if (PhaseIsOrderDisorderMapped(iPhaseLocal)) then
                ActiveSetContainsOrderDisorderMapped = .TRUE.
                return
            end if
        end do

        return

    end function ActiveSetContainsOrderDisorderMapped

    logical function TraceElementEdgeSystem()

        implicit none

        integer :: iElementLocal
        real(8) :: dTotalElementMoles, dMinPositiveElementFraction

        TraceElementEdgeSystem = .FALSE.
        dTotalElementMoles = SUM(DMAX1(dMolesElement(1:nElements), 0D0))
        if (dTotalElementMoles <= 0D0) return

        dMinPositiveElementFraction = HUGE(1D0)
        do iElementLocal = 1, nElements
            if (dMolesElement(iElementLocal) <= 0D0) cycle
            dMinPositiveElementFraction = DMIN1( &
                dMinPositiveElementFraction, &
                dMolesElement(iElementLocal) / dTotalElementMoles)
        end do

        if (dMinPositiveElementFraction < 1D-4) TraceElementEdgeSystem = .TRUE.

        return

    end function TraceElementEdgeSystem

    subroutine BuildCEFSiteFromMolFraction(iSolnPhaseIndexIn, iPhaseIDIn, dSiteOut)

        implicit none

        integer, intent(in) :: iSolnPhaseIndexIn, iPhaseIDIn
        real(8), intent(out) :: dSiteOut(nMaxSublatticeSys,nMaxConstituentSys)

        integer :: iFirstLocal, iLastLocal, iLocalSpecies, iLocal, iSub, iCon

        ! Collapse the Leveling/PEA product-fraction handoff into sublattice
        ! marginals.  The following Newton system is built in these site
        ! variables, not in the original product fractions.
        dSiteOut = 0D0
        iFirstLocal = nSpeciesPhase(iSolnPhaseIndexIn-1) + 1
        iLastLocal = nSpeciesPhase(iSolnPhaseIndexIn)
        do iLocalSpecies = iFirstLocal, iLastLocal
            iLocal = iLocalSpecies - iFirstLocal + 1
            do iSub = 1, nSublatticePhase(iPhaseIDIn)
                iCon = iConstituentSublattice(iPhaseIDIn,iSub,iLocal)
                if (iCon > 0) dSiteOut(iSub,iCon) = dSiteOut(iSub,iCon) + &
                    DMAX1(dMolFraction(iLocalSpecies), 0D0)
            end do
        end do

        return

    end subroutine BuildCEFSiteFromMolFraction

    subroutine NormalizeCEFSite(iPhaseIDIn, dSiteInOut)

        implicit none

        integer, intent(in) :: iPhaseIDIn
        real(8), intent(inout) :: dSiteInOut(nMaxSublatticeSys,nMaxConstituentSys)

        integer :: iSub, iCon
        real(8) :: dSumLocal, dFloor

        dFloor = DMAX1(dTraceSpeciesRemoveFraction, 1D-30)
        do iSub = 1, nSublatticePhase(iPhaseIDIn)
            dSumLocal = 0D0
            do iCon = 1, nConstituentSublattice(iPhaseIDIn,iSub)
                dSiteInOut(iSub,iCon) = DMAX1(dSiteInOut(iSub,iCon), dFloor)
                dSumLocal = dSumLocal + dSiteInOut(iSub,iCon)
            end do
            if (dSumLocal > 0D0) then
                do iCon = 1, nConstituentSublattice(iPhaseIDIn,iSub)
                    dSiteInOut(iSub,iCon) = dSiteInOut(iSub,iCon) / dSumLocal
                end do
            end if
        end do

        return

    end subroutine NormalizeCEFSite

    subroutine SetCEFMolFractionFromStoredSite(iSolnPhaseIndexIn, iPhaseIDIn)

        implicit none

        integer, intent(in) :: iSolnPhaseIndexIn, iPhaseIDIn

        call SetCEFMolFractionFromSite(iSolnPhaseIndexIn, iPhaseIDIn, dGEMCEFSiteLast(iPhaseIDIn,:,:))

        return

    end subroutine SetCEFMolFractionFromStoredSite

    subroutine SetCEFMolFractionFromSite(iSolnPhaseIndexIn, iPhaseIDIn, dSiteIn)

        implicit none

        integer, intent(in) :: iSolnPhaseIndexIn, iPhaseIDIn
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)

        integer :: iFirstLocal, iLastLocal, iLocalSpecies, iLocal, iSub, iCon
        real(8) :: dProduct, dSumLocal

        ! Rebuild product fractions from accepted site fractions only as the
        ! input representation expected by CompGradientSUBL/CompHessianSUBL
        ! and the shared thermodynamic arrays.
        iFirstLocal = nSpeciesPhase(iSolnPhaseIndexIn-1) + 1
        iLastLocal = nSpeciesPhase(iSolnPhaseIndexIn)
        dSumLocal = 0D0
        do iLocalSpecies = iFirstLocal, iLastLocal
            iLocal = iLocalSpecies - iFirstLocal + 1
            dProduct = 1D0
            do iSub = 1, nSublatticePhase(iPhaseIDIn)
                iCon = iConstituentSublattice(iPhaseIDIn,iSub,iLocal)
                if (iCon > 0) dProduct = dProduct * DMAX1(dSiteIn(iSub,iCon), 1D-75)
            end do
            dMolFraction(iLocalSpecies) = dProduct
            dSumLocal = dSumLocal + dProduct
        end do
        if (dSumLocal > 0D0) then
            do iLocalSpecies = iFirstLocal, iLastLocal
                dMolFraction(iLocalSpecies) = dMolFraction(iLocalSpecies) / dSumLocal
            end do
        end if

        return

    end subroutine SetCEFMolFractionFromSite

    subroutine ComputeCEFPhaseComposition(iSolnPhaseIndexIn, iPhaseIDIn, dSiteIn, dCompositionOut)

        implicit none

        integer, intent(in) :: iSolnPhaseIndexIn, iPhaseIDIn
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)
        real(8), intent(out) :: dCompositionOut(nElements)

        integer :: iFirstLocal, iLastLocal, iLocalSpecies, iLocal, iSub, iCon
        real(8) :: dProduct, dProductSum
        real(8), allocatable :: dProductSpecies(:)

        dCompositionOut = 0D0
        iFirstLocal = nSpeciesPhase(iSolnPhaseIndexIn-1) + 1
        iLastLocal = nSpeciesPhase(iSolnPhaseIndexIn)
        allocate(dProductSpecies(iLastLocal-iFirstLocal+1))
        dProductSpecies = 0D0
        dProductSum = 0D0

        do iLocalSpecies = iFirstLocal, iLastLocal
            iLocal = iLocalSpecies - iFirstLocal + 1
            dProduct = 1D0
            do iSub = 1, nSublatticePhase(iPhaseIDIn)
                iCon = iConstituentSublattice(iPhaseIDIn,iSub,iLocal)
                if (iCon > 0) dProduct = dProduct * DMAX1(dSiteIn(iSub,iCon), 1D-75)
            end do
            dProductSpecies(iLocal) = dProduct
            dProductSum = dProductSum + dProduct
        end do

        if (dProductSum > 0D0) then
            do iLocalSpecies = iFirstLocal, iLastLocal
                iLocal = iLocalSpecies - iFirstLocal + 1
                dCompositionOut = dCompositionOut + &
                    dProductSpecies(iLocal) * dStoichSpecies(iLocalSpecies,1:nElements) / dProductSum
            end do
        end if

        if (allocated(dProductSpecies)) deallocate(dProductSpecies)

        return

    end subroutine ComputeCEFPhaseComposition

    subroutine ComputeCEFCompositionDerivative(iSolnPhaseIndexIn, iPhaseIDIn, iSubIn, iConIn, &
        iRefIn, dSiteIn, dDerivativeOut)

        implicit none

        integer, intent(in) :: iSolnPhaseIndexIn, iPhaseIDIn, iSubIn, iConIn, iRefIn
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)
        real(8), intent(out) :: dDerivativeOut(nElements)

        integer :: iFirstLocal, iLastLocal, iLocalSpecies, iLocal, iSub, iCon
        real(8) :: dProduct, dProductSum, dProductDerivative, dProductDerivativeSum
        real(8) :: dYCon, dYRef, dDerivativeFraction
        real(8), allocatable :: dProductSpecies(:), dProductDerivativeSpecies(:)

        dDerivativeOut = 0D0

        iFirstLocal = nSpeciesPhase(iSolnPhaseIndexIn-1) + 1
        iLastLocal = nSpeciesPhase(iSolnPhaseIndexIn)
        allocate(dProductSpecies(iLastLocal-iFirstLocal+1), &
            dProductDerivativeSpecies(iLastLocal-iFirstLocal+1))
        dProductSpecies = 0D0
        dProductDerivativeSpecies = 0D0
        dProductSum = 0D0
        dProductDerivativeSum = 0D0
        dYCon = DMAX1(dSiteIn(iSubIn,iConIn), 1D-75)
        dYRef = DMAX1(dSiteIn(iSubIn,iRefIn), 1D-75)

        do iLocalSpecies = iFirstLocal, iLastLocal
            iLocal = iLocalSpecies - iFirstLocal + 1
            dProduct = 1D0
            do iSub = 1, nSublatticePhase(iPhaseIDIn)
                iCon = iConstituentSublattice(iPhaseIDIn,iSub,iLocal)
                if (iCon > 0) dProduct = dProduct * DMAX1(dSiteIn(iSub,iCon), 1D-75)
            end do

            dProductDerivative = 0D0
            iCon = iConstituentSublattice(iPhaseIDIn,iSubIn,iLocal)
            if (iCon == iConIn) dProductDerivative = dProductDerivative + dProduct / dYCon
            if (iCon == iRefIn) dProductDerivative = dProductDerivative - dProduct / dYRef

            dProductSpecies(iLocal) = dProduct
            dProductDerivativeSpecies(iLocal) = dProductDerivative
            dProductSum = dProductSum + dProduct
            dProductDerivativeSum = dProductDerivativeSum + dProductDerivative
        end do

        if (dProductSum > 0D0) then
            do iLocalSpecies = iFirstLocal, iLastLocal
                iLocal = iLocalSpecies - iFirstLocal + 1
                dDerivativeFraction = (dProductDerivativeSpecies(iLocal) * dProductSum - &
                    dProductSpecies(iLocal) * dProductDerivativeSum) / dProductSum**2
                dDerivativeOut = dDerivativeOut + &
                    dDerivativeFraction * dStoichSpecies(iLocalSpecies,1:nElements)
            end do
        end if

        if (allocated(dProductSpecies)) deallocate(dProductSpecies)
        if (allocated(dProductDerivativeSpecies)) deallocate(dProductDerivativeSpecies)

        return

    end subroutine ComputeCEFCompositionDerivative

    integer function LocalCEFIndex(iPhaseIDIn, iSubIn, iConIn)

        implicit none

        integer, intent(in) :: iPhaseIDIn, iSubIn, iConIn
        integer :: iSubLocal

        LocalCEFIndex = iConIn
        do iSubLocal = 1, iSubIn - 1
            LocalCEFIndex = LocalCEFIndex + nConstituentSublattice(iPhaseIDIn,iSubLocal)
        end do

        return

    end function LocalCEFIndex

    real(8) function ProjectCEFHessian(dHessianIn, iPhaseIDIn, iSubA, iConA, iRefA, &
        iSubB, iConB, iRefB)

        implicit none

        real(8), intent(in) :: dHessianIn(nSiteCapacity,nSiteCapacity)
        integer, intent(in) :: iPhaseIDIn, iSubA, iConA, iRefA, iSubB, iConB, iRefB
        integer :: iA, iRA, iB, iRB

        iA = LocalCEFIndex(iPhaseIDIn, iSubA, iConA)
        iRA = LocalCEFIndex(iPhaseIDIn, iSubA, iRefA)
        iB = LocalCEFIndex(iPhaseIDIn, iSubB, iConB)
        iRB = LocalCEFIndex(iPhaseIDIn, iSubB, iRefB)

        ProjectCEFHessian = dHessianIn(iA,iB) - dHessianIn(iA,iRB) - &
            dHessianIn(iRA,iB) + dHessianIn(iRA,iRB)

        return

    end function ProjectCEFHessian

    logical function IsCEFPhase(iSolnPhaseIndexIn)

        implicit none

        integer, intent(in) :: iSolnPhaseIndexIn
        character(len=8) :: cTypeLocal

        IsCEFPhase = .FALSE.
        if ((iSolnPhaseIndexIn <= 0).OR.(iSolnPhaseIndexIn > nSolnPhasesSys)) return
        if (.NOT.allocated(iPhaseSublattice)) return
        if (iSolnPhaseIndexIn > SIZE(iPhaseSublattice)) return
        if (iPhaseSublattice(iSolnPhaseIndexIn) <= 0) return

        cTypeLocal = TRIM(cSolnPhaseType(iSolnPhaseIndexIn))
        IsCEFPhase = (TRIM(cTypeLocal) == 'SUBL').OR.&
            (TRIM(cTypeLocal) == 'SUBLM').OR.&
            (TRIM(cTypeLocal) == 'SUBOM')

        return

    end function IsCEFPhase

    subroutine ComputeNonCEFPhaseComposition(iSolnPhaseIndexIn, dCompositionOut)

        implicit none

        integer, intent(in) :: iSolnPhaseIndexIn
        real(8), intent(out) :: dCompositionOut(nElements)

        integer :: iFirstLocal, iLastLocal, iSpeciesLocal

        dCompositionOut = 0D0
        iFirstLocal = nSpeciesPhase(iSolnPhaseIndexIn-1) + 1
        iLastLocal = nSpeciesPhase(iSolnPhaseIndexIn)

        do iSpeciesLocal = iFirstLocal, iLastLocal
            dCompositionOut = dCompositionOut + dMolFraction(iSpeciesLocal) * &
                dStoichSpecies(iSpeciesLocal,1:nElements) / DBLE(iParticlesPerMole(iSpeciesLocal))
        end do

        return

    end subroutine ComputeNonCEFPhaseComposition

    subroutine ComputeNonCEFCompositionDerivative(iSolnPhaseIndexIn, iConIn, iRefIn, dDerivativeOut)

        implicit none

        integer, intent(in) :: iSolnPhaseIndexIn, iConIn, iRefIn
        real(8), intent(out) :: dDerivativeOut(nElements)

        integer :: iFirstLocal, iSpeciesLocal, iSpeciesRefLocal

        iFirstLocal = nSpeciesPhase(iSolnPhaseIndexIn-1) + 1
        iSpeciesLocal = iFirstLocal + iConIn - 1
        iSpeciesRefLocal = iFirstLocal + iRefIn - 1

        dDerivativeOut = dStoichSpecies(iSpeciesLocal,1:nElements) / &
            DBLE(iParticlesPerMole(iSpeciesLocal)) - &
            dStoichSpecies(iSpeciesRefLocal,1:nElements) / DBLE(iParticlesPerMole(iSpeciesRefLocal))

        return

    end subroutine ComputeNonCEFCompositionDerivative

    subroutine CompNonCEFHessian(iSolnPhaseIndexIn, nSpeciesDimIn, dHessianOut, nSpeciesOutLocal, iInfoOut)

        implicit none

        integer, intent(in) :: iSolnPhaseIndexIn, nSpeciesDimIn
        real(8), intent(out) :: dHessianOut(nSpeciesDimIn,nSpeciesDimIn)
        integer, intent(out) :: nSpeciesOutLocal, iInfoOut

        real(8), allocatable :: dChemicalPotentialScratch(:)
        integer :: iFirstLocal, iLastLocal, iLocal
        character(len=8) :: cTypeLocal

        dHessianOut = 0D0
        nSpeciesOutLocal = 0
        iInfoOut = 0
        allocate(dChemicalPotentialScratch(nSpeciesDimIn))
        dChemicalPotentialScratch = 0D0

        cTypeLocal = TRIM(cSolnPhaseType(iSolnPhaseIndexIn))
        select case (TRIM(cTypeLocal))
        case ('SUBG', 'SUBQ')
            call CompHessianSUBG(iSolnPhaseIndexIn, nSpeciesDimIn, dChemicalPotentialScratch, &
                dHessianOut, nSpeciesOutLocal, iInfoOut)
        case ('RKMP')
            call CompHessianRKMP(iSolnPhaseIndexIn, nSpeciesDimIn, dChemicalPotentialScratch, &
                dHessianOut, nSpeciesOutLocal, iInfoOut)
        case ('QKTO')
            call CompHessianQKTO(iSolnPhaseIndexIn, nSpeciesDimIn, dChemicalPotentialScratch, &
                dHessianOut, nSpeciesOutLocal, iInfoOut)
        case ('IDMX')
            iFirstLocal = nSpeciesPhase(iSolnPhaseIndexIn-1) + 1
            iLastLocal = nSpeciesPhase(iSolnPhaseIndexIn)
            nSpeciesOutLocal = iLastLocal - iFirstLocal + 1
            if ((nSpeciesOutLocal <= 0).OR.(nSpeciesOutLocal > nSpeciesDimIn)) then
                iInfoOut = 1
            else
                do iLocal = 1, nSpeciesOutLocal
                    dHessianOut(iLocal,iLocal) = 1D0 / &
                        DMAX1(dMolFraction(iFirstLocal+iLocal-1), 1D-75)
                end do
            end if
        case default
            iInfoOut = 1
        end select

        if (allocated(dChemicalPotentialScratch)) deallocate(dChemicalPotentialScratch)

        return

    end subroutine CompNonCEFHessian

    real(8) function ProjectNonCEFHessian(dHessianIn, iConA, iRefA, iConB, iRefB)

        implicit none

        real(8), intent(in) :: dHessianIn(nSpeciesCapacity,nSpeciesCapacity)
        integer, intent(in) :: iConA, iRefA, iConB, iRefB

        ProjectNonCEFHessian = dHessianIn(iConA,iConB) - dHessianIn(iConA,iRefB) - &
            dHessianIn(iRefA,iConB) + dHessianIn(iRefA,iRefB)

        return

    end function ProjectNonCEFHessian

    subroutine CleanupGEMNewtonCEF

        implicit none

        integer :: iStat

        iStat = 0
        if (allocated(A)) deallocate(A)
        if (allocated(B)) deallocate(B)
        if (allocated(ABase)) deallocate(ABase)
        if (allocated(BBase)) deallocate(BBase)
        if (allocated(AFactorBase)) deallocate(AFactorBase)
        if (allocated(BSolutionBase)) deallocate(BSolutionBase)
        if (allocated(IPIV)) deallocate(IPIV)
        if (allocated(IPIVBase)) deallocate(IPIVBase)
        if (allocated(dSiteGradient)) deallocate(dSiteGradient)
        if (allocated(dSiteGradientH)) deallocate(dSiteGradientH)
        if (allocated(dSiteGradientS)) deallocate(dSiteGradientS)
        if (allocated(dSiteGradientCp)) deallocate(dSiteGradientCp)
        if (allocated(dSiteHessian)) deallocate(dSiteHessian)
        if (allocated(dAmountHessian)) deallocate(dAmountHessian)
        if (allocated(dJacobian)) deallocate(dJacobian)
        if (allocated(dPhaseComp)) deallocate(dPhaseComp)
        if (allocated(dSiteCompDeriv)) deallocate(dSiteCompDeriv)
        if (allocated(dSiteResidualByVar)) deallocate(dSiteResidualByVar)
        if (iStat /= 0) INFOThermo = 24

        return

    end subroutine CleanupGEMNewtonCEF

end subroutine GEMNewtonCEF
