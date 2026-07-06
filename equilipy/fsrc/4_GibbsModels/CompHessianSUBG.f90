!> \brief Compute analytical SUBG/SUBQ chemical-potential Jacobians.
!!
!! \details This diagnostic/minimizer routine evaluates the SUBG/SUBQ
!! model-species chemical potentials and their analytical derivatives with
!! respect to formal pair/quadruplet mole fractions.  The returned Hessian is
!! the redundant composition-space Jacobian \f$\partial\mu_i/\partial x_j\f$;
!! callers should project it before solving constrained Newton systems.
!!
!! \param[in]  iSolnIndex                Absolute solution phase index.
!! \param[in]  nPairDim                  Caller-provided pair/quad capacity.
!! \param[out] dPairChemicalPotential    Chemical potentials in internal G/(RT) units.
!! \param[out] dPairHessian              Composition Jacobian in internal G/(RT) units.
!! \param[out] nPairOut                  Number of active pair/quad variables written.
!! \param[out] iInfo                     Zero on success; nonzero for unsupported/undersized calls.



subroutine CompHessianSUBG(iSolnIndex, nPairDim, dPairChemicalPotential, dPairHessian, nPairOut, iInfo)
!
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CompHessianSUBG.f90
    !> \brief   Compute analytical SUBG/SUBQ chemical-potential Jacobians.
    !> \author  S.Y. Kwon
    !> \date    Jun. 27, 2026
    !> \sa      CompExcessGibbsEnergySUBG.f90
    !> \sa      solution_models.md
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   06/27/2026      S.Y. Kwon           Added analytical SUBG/SUBQ base and one-anion Q-term Jacobian.
    !   06/27/2026      S.Y. Kwon           Removed in-place quotient aliasing in analytical derivative setup.
    !   06/27/2026      S.Y. Kwon           Added one-anion G-type SUBQ excess Jacobian support.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details This routine evaluates \f$\mu_q\f$ and
    !! \f$J_{qr}=\partial \mu_q/\partial x_r\f$ for formal SUBG/SUBQ
    !! pair/quadruplet fractions.  It follows the same variable convention as
    !! CompExcessGibbsEnergySUBG: dMolFraction is the formal pair/quadruplet
    !! fraction, dChemicalPotential is the derivative with respect to that
    !! variable, and dStoichSpecies stores the formal composition used by the
    !! elemental-potential plane.
    !
    !
    ! Required input variables:
    ! =========================
    !
    !> \param[in] iSolnIndex  Absolute solution phase index.
    !> \param[in] nPairDim    First/second dimension of output arrays.
    !
    ! cSolnPhaseType             Solution model type. SUBG and SUBQ are accepted.
    ! dMolFraction               Current formal pair/quadruplet fractions.
    ! dStdGibbsEnergy            Dimensionless standard pair/quadruplet Gibbs energies.
    ! iPairID                    Pair/quadruplet topology.
    ! dCoordinationNumber        Coordination-number data for formal-complex stoichiometry.
    ! dZetaSpecies               SRO pair normalization.
    ! dExcessGibbsParam          Dimensionless SUBG/SUBQ excess parameters.
    !
    !
    ! Output/updated variables:
    ! =========================
    !
    !> \param[out] dPairChemicalPotential  Pair/quad chemical potentials in G/(RT) units.
    !> \param[out] dPairHessian            Redundant composition Jacobian dmu_i/dx_j in G/(RT) units.
    !> \param[out] nPairOut                Number of active pair/quad variables written.
    !> \param[out] iInfo                   Zero on success. 1=unsupported phase, 2=bad topology,
    !!                                      3=caller arrays too small, 4=unsupported excess topology.
    !
    !
    ! Called subroutines/functions:
    ! =============================
    !
    ! BuildMQMCompositionDerivatives  Builds normalized SUBG/SUBQ auxiliary variables and derivatives.
    ! AccumulateBaseJacobian          Adds standard and configurational chemical potentials/Jacobian.
    ! AccumulateQExcessJacobian       Adds supported G/Q-type excess chemical-potential derivatives.
    ! CompExcessGibbsEnergy           Recomputes runtime chemical potentials for the output vector.
    !
    !
    ! Primary callers:
    ! ================
    !
    ! Python diagnostics/devop tests      Validate analytical curvature against finite differences.
    ! Future SUBG/SUBQ Subminimization   Candidate source for phase-local Newton directions.
    !
    !
    ! Numerical assumptions:
    ! ======================
    !
    ! - Returned values use internal dimensionless Gibbs-energy convention, G/(RT).
    ! - dPairHessian is the redundant composition-space Jacobian.  It is singular
    !   under the sum(x)=1 constraint and should be projected before inversion.
    ! - The base SUBG/SUBQ terms are fully analytical.  The first excess support
    !   covers G/Q-type one-anion terms used by current MQM fixtures.  Other
    !   excess topologies return iInfo=4 until their chain-rule derivatives are
    !   implemented explicitly.
    !
    !-------------------------------------------------------------------------------------------------------------
    USE ModuleThermo
    USE ModuleThermoIO

    implicit none

    integer, intent(in)  :: iSolnIndex, nPairDim
    real(8), intent(out) :: dPairChemicalPotential(nPairDim)
    real(8), intent(out) :: dPairHessian(nPairDim,nPairDim)
    integer, intent(out) :: nPairOut, iInfo

    integer :: iSPI, iFirst, iLast
    integer :: nSub1, nSub2, nPhaseElements
    real(8), allocatable :: dNi(:), dYi(:), dFi(:)
    real(8), allocatable :: dNiGrad(:,:), dYiGrad(:,:), dFiGrad(:,:)
    real(8), allocatable :: dNij(:,:), dNsij(:,:), dXij(:,:), dXsij(:,:)
    real(8), allocatable :: dNijGrad(:,:,:), dNsijGrad(:,:,:)
    real(8), allocatable :: dXijGrad(:,:,:), dXsijGrad(:,:,:)
    real(8), allocatable :: dBaseMu(:), dBaseJac(:,:)
    real(8), allocatable :: dExcessMu(:), dExcessJac(:,:)

    dPairChemicalPotential = 0D0
    dPairHessian = 0D0
    nPairOut = 0
    iInfo = 0

    if ((iSolnIndex <= 0).OR.(iSolnIndex > nSolnPhasesSys)) then
        iInfo = 1
        return
    end if
    if ((TRIM(cSolnPhaseType(iSolnIndex)) /= 'SUBG').AND.&
        (TRIM(cSolnPhaseType(iSolnIndex)) /= 'SUBQ')) then
        iInfo = 1
        return
    end if

    iSPI = iPhaseSublattice(iSolnIndex)
    if (iSPI <= 0) then
        iInfo = 2
        return
    end if

    iFirst = nSpeciesPhase(iSolnIndex-1) + 1
    iLast = nSpeciesPhase(iSolnIndex)
    nPairOut = iLast - iFirst + 1
    if ((nPairOut <= 0).OR.(nPairOut /= nPairsSRO(iSPI,2))) then
        iInfo = 2
        return
    end if
    if (nPairOut > nPairDim) then
        iInfo = 3
        return
    end if

    nSub1 = nSublatticeElements(iSPI,1)
    nSub2 = nSublatticeElements(iSPI,2)
    nPhaseElements = nSub1 + nSub2

    allocate(dNi(nPhaseElements), dYi(nPhaseElements), dFi(nPhaseElements))
    allocate(dNiGrad(nPhaseElements,nPairOut), dYiGrad(nPhaseElements,nPairOut))
    allocate(dFiGrad(nPhaseElements,nPairOut))
    allocate(dNij(nSub1,nSub2), dNsij(nSub1,nSub2), dXij(nSub1,nSub2), dXsij(nSub1,nSub2))
    allocate(dNijGrad(nSub1,nSub2,nPairOut), dNsijGrad(nSub1,nSub2,nPairOut))
    allocate(dXijGrad(nSub1,nSub2,nPairOut), dXsijGrad(nSub1,nSub2,nPairOut))
    allocate(dBaseMu(nPairOut), dBaseJac(nPairOut,nPairOut))
    allocate(dExcessMu(nPairOut), dExcessJac(nPairOut,nPairOut))

    call BuildMQMCompositionDerivatives(iSPI, iFirst, nPairOut, nSub1, nSub2, nPhaseElements, &
        dNi, dNiGrad, dYi, dYiGrad, dNij, dNijGrad, dNsij, dNsijGrad, &
        dXij, dXijGrad, dXsij, dXsijGrad, dFi, dFiGrad, iInfo)
    if (iInfo /= 0) goto 900

    call AccumulateBaseJacobian(iSolnIndex, iSPI, iFirst, nPairOut, nSub1, nSub2, nPhaseElements, &
        dNi, dNiGrad, dYi, dYiGrad, dXij, dXijGrad, dXsij, dXsijGrad, &
        dFi, dFiGrad, dBaseMu, dBaseJac)

    dExcessMu = 0D0
    dExcessJac = 0D0
    call AccumulateQExcessJacobian(iSolnIndex, iSPI, iFirst, nPairOut, nSub1, nSub2, &
        dXij, dXijGrad, dYi, dYiGrad, dExcessMu, dExcessJac, iInfo)
    if (iInfo /= 0) goto 900

    dPairHessian(1:nPairOut,1:nPairOut) = dBaseJac + dExcessJac

    call CompExcessGibbsEnergy(iSolnIndex)
    dPairChemicalPotential(1:nPairOut) = dChemicalPotential(iFirst:iLast)

900 continue
    if (allocated(dNi)) deallocate(dNi)
    if (allocated(dYi)) deallocate(dYi)
    if (allocated(dFi)) deallocate(dFi)
    if (allocated(dNiGrad)) deallocate(dNiGrad)
    if (allocated(dYiGrad)) deallocate(dYiGrad)
    if (allocated(dFiGrad)) deallocate(dFiGrad)
    if (allocated(dNij)) deallocate(dNij)
    if (allocated(dNsij)) deallocate(dNsij)
    if (allocated(dXij)) deallocate(dXij)
    if (allocated(dXsij)) deallocate(dXsij)
    if (allocated(dNijGrad)) deallocate(dNijGrad)
    if (allocated(dNsijGrad)) deallocate(dNsijGrad)
    if (allocated(dXijGrad)) deallocate(dXijGrad)
    if (allocated(dXsijGrad)) deallocate(dXsijGrad)
    if (allocated(dBaseMu)) deallocate(dBaseMu)
    if (allocated(dBaseJac)) deallocate(dBaseJac)
    if (allocated(dExcessMu)) deallocate(dExcessMu)
    if (allocated(dExcessJac)) deallocate(dExcessJac)

    return

contains

    subroutine BuildMQMCompositionDerivatives(iSPIIn, iFirstIn, nPairLocal, nSub1In, nSub2In, nPhaseElementsIn, &
        dNiOut, dNiGradOut, dYiOut, dYiGradOut, dNijOut, dNijGradOut, dNsijOut, dNsijGradOut, &
        dXijOut, dXijGradOut, dXsijOut, dXsijGradOut, dFiOut, dFiGradOut, iInfoOut)

        implicit none

        integer, intent(in) :: iSPIIn, iFirstIn, nPairLocal, nSub1In, nSub2In, nPhaseElementsIn
        real(8), intent(out) :: dNiOut(nPhaseElementsIn), dNiGradOut(nPhaseElementsIn,nPairLocal)
        real(8), intent(out) :: dYiOut(nPhaseElementsIn), dYiGradOut(nPhaseElementsIn,nPairLocal)
        real(8), intent(out) :: dNijOut(nSub1In,nSub2In), dNijGradOut(nSub1In,nSub2In,nPairLocal)
        real(8), intent(out) :: dNsijOut(nSub1In,nSub2In), dNsijGradOut(nSub1In,nSub2In,nPairLocal)
        real(8), intent(out) :: dXijOut(nSub1In,nSub2In), dXijGradOut(nSub1In,nSub2In,nPairLocal)
        real(8), intent(out) :: dXsijOut(nSub1In,nSub2In), dXsijGradOut(nSub1In,nSub2In,nPairLocal)
        real(8), intent(out) :: dFiOut(nPhaseElementsIn), dFiGradOut(nPhaseElementsIn,nPairLocal)
        integer, intent(out) :: iInfoOut

        integer :: i, j, k, l, m, q, nA, nX
        real(8) :: dZa, dZb, dZx, dZy
        real(8) :: dSumCation, dSumAnion, dSumNij, dSumNsij, dTempValue
        real(8) :: dSumCationGrad(nPairLocal), dSumAnionGrad(nPairLocal)
        real(8) :: dSumNijGrad(nPairLocal), dSumNsijGrad(nPairLocal), dTempGrad(nPairLocal)

        iInfoOut = 0
        dNiOut = 0D0
        dYiOut = 0D0
        dFiOut = 0D0
        dNijOut = 0D0
        dNsijOut = 0D0
        dXijOut = 0D0
        dXsijOut = 0D0
        dNiGradOut = 0D0
        dYiGradOut = 0D0
        dFiGradOut = 0D0
        dNijGradOut = 0D0
        dNsijGradOut = 0D0
        dXijGradOut = 0D0
        dXsijGradOut = 0D0

        dSumCation = 0D0
        dSumAnion = 0D0
        dSumCationGrad = 0D0
        dSumAnionGrad = 0D0

        do i = 1, nSub1In
            do k = 1, nPairLocal
                l = iFirstIn + k - 1
                dZa = dCoordinationNumber(iSPIIn,k,1)
                dZb = dCoordinationNumber(iSPIIn,k,2)
                if (dZa <= 0D0 .OR. dZb <= 0D0) then
                    iInfoOut = 2
                    return
                end if
                if (i == iPairID(iSPIIn,k,1)) then
                    dNiOut(i) = dNiOut(i) + dMolFraction(l) / dZa
                    dYiOut(i) = dYiOut(i) + dMolFraction(l) / 2D0
                    dNiGradOut(i,k) = dNiGradOut(i,k) + 1D0 / dZa
                    dYiGradOut(i,k) = dYiGradOut(i,k) + 0.5D0
                end if
                if (i == iPairID(iSPIIn,k,2)) then
                    dNiOut(i) = dNiOut(i) + dMolFraction(l) / dZb
                    dYiOut(i) = dYiOut(i) + dMolFraction(l) / 2D0
                    dNiGradOut(i,k) = dNiGradOut(i,k) + 1D0 / dZb
                    dYiGradOut(i,k) = dYiGradOut(i,k) + 0.5D0
                end if
            end do
            dSumCation = dSumCation + dNiOut(i)
            dSumCationGrad = dSumCationGrad + dNiGradOut(i,:)
        end do

        do i = 1, nSub2In
            j = i + nSub1In
            do k = 1, nPairLocal
                l = iFirstIn + k - 1
                dZx = dCoordinationNumber(iSPIIn,k,3)
                dZy = dCoordinationNumber(iSPIIn,k,4)
                if (dZx <= 0D0 .OR. dZy <= 0D0) then
                    iInfoOut = 2
                    return
                end if
                if (j == iPairID(iSPIIn,k,3)) then
                    dNiOut(j) = dNiOut(j) + dMolFraction(l) / dZx
                    dYiOut(j) = dYiOut(j) + dMolFraction(l) / 2D0
                    dNiGradOut(j,k) = dNiGradOut(j,k) + 1D0 / dZx
                    dYiGradOut(j,k) = dYiGradOut(j,k) + 0.5D0
                end if
                if (j == iPairID(iSPIIn,k,4)) then
                    dNiOut(j) = dNiOut(j) + dMolFraction(l) / dZy
                    dYiOut(j) = dYiOut(j) + dMolFraction(l) / 2D0
                    dNiGradOut(j,k) = dNiGradOut(j,k) + 1D0 / dZy
                    dYiGradOut(j,k) = dYiGradOut(j,k) + 0.5D0
                end if
            end do
            dSumAnion = dSumAnion + dNiOut(j)
            dSumAnionGrad = dSumAnionGrad + dNiGradOut(j,:)
        end do

        if (dSumCation <= 0D0 .OR. dSumAnion <= 0D0) then
            iInfoOut = 2
            return
        end if

        do i = 1, nSub1In
            call QuotientGradient(nPairLocal, dNiOut(i), dNiGradOut(i,:), dSumCation, dSumCationGrad, &
                dTempValue, dTempGrad)
            dNiOut(i) = dTempValue
            dNiGradOut(i,:) = dTempGrad
        end do
        do i = 1, nSub2In
            j = i + nSub1In
            call QuotientGradient(nPairLocal, dNiOut(j), dNiGradOut(j,:), dSumAnion, dSumAnionGrad, &
                dTempValue, dTempGrad)
            dNiOut(j) = dTempValue
            dNiGradOut(j,:) = dTempGrad
        end do

        dSumNij = 0D0
        dSumNsij = 0D0
        dSumNijGrad = 0D0
        dSumNsijGrad = 0D0
        do i = 1, nSub1In
            do j = 1, nSub2In
                m = iConstituentSublattice(iSPIIn,1,i) + &
                    ((iConstituentSublattice(iSPIIn,2,j) - 1) * nSub1In)
                if (dZetaSpecies(iSPIIn,m) <= 0D0) then
                    iInfoOut = 2
                    return
                end if
                do k = 1, nPairLocal
                    l = iFirstIn + k - 1
                    nA = 0
                    if (i == iPairID(iSPIIn,k,1)) nA = nA + 1
                    if (i == iPairID(iSPIIn,k,2)) nA = nA + 1
                    nX = 0
                    if ((j + nSub1In) == iPairID(iSPIIn,k,3)) nX = nX + 1
                    if ((j + nSub1In) == iPairID(iSPIIn,k,4)) nX = nX + 1
                    dNijOut(i,j) = dNijOut(i,j) + dMolFraction(l) * DBLE(nA * nX)
                    dNsijOut(i,j) = dNsijOut(i,j) + dMolFraction(l) * DBLE(nA * nX) / dZetaSpecies(iSPIIn,m)
                    dNijGradOut(i,j,k) = dNijGradOut(i,j,k) + DBLE(nA * nX)
                    dNsijGradOut(i,j,k) = dNsijGradOut(i,j,k) + DBLE(nA * nX) / dZetaSpecies(iSPIIn,m)
                end do
                dSumNij = dSumNij + dNijOut(i,j)
                dSumNsij = dSumNsij + dNsijOut(i,j)
                dSumNijGrad = dSumNijGrad + dNijGradOut(i,j,:)
                dSumNsijGrad = dSumNsijGrad + dNsijGradOut(i,j,:)
            end do
        end do

        if (dSumNij <= 0D0 .OR. dSumNsij <= 0D0) then
            iInfoOut = 2
            return
        end if

        do i = 1, nSub1In
            do j = 1, nSub2In
                call QuotientGradient(nPairLocal, dNijOut(i,j), dNijGradOut(i,j,:), dSumNij, dSumNijGrad, &
                    dXijOut(i,j), dXijGradOut(i,j,:))
                call QuotientGradient(nPairLocal, dNsijOut(i,j), dNsijGradOut(i,j,:), dSumNsij, dSumNsijGrad, &
                    dXsijOut(i,j), dXsijGradOut(i,j,:))
            end do
        end do

        do i = 1, nSub1In
            do j = 1, nSub2In
                dFiOut(i) = dFiOut(i) + dXsijOut(i,j)
                dFiGradOut(i,:) = dFiGradOut(i,:) + dXsijGradOut(i,j,:)
                k = j + nSub1In
                dFiOut(k) = dFiOut(k) + dXsijOut(i,j)
                dFiGradOut(k,:) = dFiGradOut(k,:) + dXsijGradOut(i,j,:)
            end do
        end do

        return

    end subroutine BuildMQMCompositionDerivatives

    subroutine AccumulateBaseJacobian(iSolnIndexIn, iSPIIn, iFirstIn, nPairLocal, nSub1In, nSub2In, &
        nPhaseElementsIn, dXiIn, dXiGradIn, dYiIn, dYiGradIn, dXijIn, dXijGradIn, dXsijIn, dXsijGradIn, &
        dFiIn, dFiGradIn, dMuOut, dJacOut)

        implicit none

        integer, intent(in) :: iSolnIndexIn, iSPIIn, iFirstIn, nPairLocal, nSub1In, nSub2In, nPhaseElementsIn
        real(8), intent(in) :: dXiIn(nPhaseElementsIn), dXiGradIn(nPhaseElementsIn,nPairLocal)
        real(8), intent(in) :: dYiIn(nPhaseElementsIn), dYiGradIn(nPhaseElementsIn,nPairLocal)
        real(8), intent(in) :: dXijIn(nSub1In,nSub2In), dXijGradIn(nSub1In,nSub2In,nPairLocal)
        real(8), intent(in) :: dXsijIn(nSub1In,nSub2In), dXsijGradIn(nSub1In,nSub2In,nPairLocal)
        real(8), intent(in) :: dFiIn(nPhaseElementsIn), dFiGradIn(nPhaseElementsIn,nPairLocal)
        real(8), intent(out) :: dMuOut(nPairLocal), dJacOut(nPairLocal,nPairLocal)

        integer :: i, j, k, l, m, ii, jj, kk, ll, ka, la, nA, nX, iWeight
        real(8) :: dZa, dZb, dZx, dZy, dPowXij, dPowYi, dSum

        dMuOut = 0D0
        dJacOut = 0D0

        do k = 1, nPairLocal
            l = iFirstIn + k - 1
            dMuOut(k) = dStdGibbsEnergy(l)
            dZa = dCoordinationNumber(iSPIIn,k,1)
            dZb = dCoordinationNumber(iSPIIn,k,2)
            dZx = dCoordinationNumber(iSPIIn,k,3)
            dZy = dCoordinationNumber(iSPIIn,k,4)

            if (iPairID(iSPIIn,k,1) > 0) then
                i = iPairID(iSPIIn,k,1)
                call AddLogTerm(nPairLocal, 1D0 / dZa, dXiIn(i), dXiGradIn(i,:), dMuOut(k), dJacOut(k,:))
            end if
            if (iPairID(iSPIIn,k,2) > 0) then
                i = iPairID(iSPIIn,k,2)
                call AddLogTerm(nPairLocal, 1D0 / dZb, dXiIn(i), dXiGradIn(i,:), dMuOut(k), dJacOut(k,:))
            end if
            if (iPairID(iSPIIn,k,3) > 0) then
                i = iPairID(iSPIIn,k,3)
                call AddLogTerm(nPairLocal, 1D0 / dZx, dXiIn(i), dXiGradIn(i,:), dMuOut(k), dJacOut(k,:))
            end if
            if (iPairID(iSPIIn,k,4) > 0) then
                i = iPairID(iSPIIn,k,4)
                call AddLogTerm(nPairLocal, 1D0 / dZy, dXiIn(i), dXiGradIn(i,:), dMuOut(k), dJacOut(k,:))
            end if

            do i = 1, nSub1In
                do j = 1, nSub2In
                    m = iConstituentSublattice(iSPIIn,1,i) + &
                        ((iConstituentSublattice(iSPIIn,2,j) - 1) * nSub1In)
                    nA = 0
                    if (i == iPairID(iSPIIn,k,1)) nA = nA + 1
                    if (i == iPairID(iSPIIn,k,2)) nA = nA + 1
                    nX = 0
                    if ((j + nSub1In) == iPairID(iSPIIn,k,3)) nX = nX + 1
                    if ((j + nSub1In) == iPairID(iSPIIn,k,4)) nX = nX + 1
                    if (nA * nX > 0) then
                        call AddLogTerm(nPairLocal, DBLE(nA * nX) / dZetaSpecies(iSPIIn,m), &
                            dXsijIn(i,j), dXsijGradIn(i,j,:), dMuOut(k), dJacOut(k,:))
                        call AddLogTerm(nPairLocal, -DBLE(nA * nX) / dZetaSpecies(iSPIIn,m), &
                            dFiIn(i), dFiGradIn(i,:), dMuOut(k), dJacOut(k,:))
                        call AddLogTerm(nPairLocal, -DBLE(nA * nX) / dZetaSpecies(iSPIIn,m), &
                            dFiIn(j+nSub1In), dFiGradIn(j+nSub1In,:), dMuOut(k), dJacOut(k,:))
                    end if
                end do
            end do

            ii = iPairID(iSPIIn,k,1)
            jj = iPairID(iSPIIn,k,2)
            kk = iPairID(iSPIIn,k,3)
            ll = iPairID(iSPIIn,k,4)
            ka = kk - nSub1In
            la = ll - nSub1In
            iWeight = 1
            if (ii /= jj) iWeight = iWeight * 2
            if (kk /= ll) iWeight = iWeight * 2

            if (TRIM(cSolnPhaseType(iSolnIndexIn)) == 'SUBG') then
                dPowXij = 1D0
                dPowYi = 1D0
            else
                dPowXij = 0.75D0
                dPowYi = 0.5D0
            end if

            dSum = dMolFraction(l)
            if (dSum <= 0D0) dSum = 1D-75
            dMuOut(k) = dMuOut(k) + DLOG(dSum) - DLOG(DBLE(iWeight))
            dJacOut(k,k) = dJacOut(k,k) + 1D0 / dSum
            call AddLogTerm(nPairLocal, -dPowXij, dXijIn(ii,ka), dXijGradIn(ii,ka,:), dMuOut(k), dJacOut(k,:))
            call AddLogTerm(nPairLocal, -dPowXij, dXijIn(ii,la), dXijGradIn(ii,la,:), dMuOut(k), dJacOut(k,:))
            call AddLogTerm(nPairLocal, -dPowXij, dXijIn(jj,ka), dXijGradIn(jj,ka,:), dMuOut(k), dJacOut(k,:))
            call AddLogTerm(nPairLocal, -dPowXij, dXijIn(jj,la), dXijGradIn(jj,la,:), dMuOut(k), dJacOut(k,:))
            call AddLogTerm(nPairLocal, dPowYi, dYiIn(ii), dYiGradIn(ii,:), dMuOut(k), dJacOut(k,:))
            call AddLogTerm(nPairLocal, dPowYi, dYiIn(jj), dYiGradIn(jj,:), dMuOut(k), dJacOut(k,:))
            call AddLogTerm(nPairLocal, dPowYi, dYiIn(kk), dYiGradIn(kk,:), dMuOut(k), dJacOut(k,:))
            call AddLogTerm(nPairLocal, dPowYi, dYiIn(ll), dYiGradIn(ll,:), dMuOut(k), dJacOut(k,:))
        end do

        return

    end subroutine AccumulateBaseJacobian

    subroutine AccumulateQExcessJacobian(iSolnIndexIn, iSPIIn, iFirstIn, nPairLocal, nSub1In, nSub2In, &
        dXijIn, dXijGradIn, dYiIn, dYiGradIn, dMuOut, dJacOut, iInfoOut)

        implicit none

        integer, intent(in) :: iSolnIndexIn, iSPIIn, iFirstIn, nPairLocal, nSub1In, nSub2In
        real(8), intent(in) :: dXijIn(nSub1In,nSub2In), dXijGradIn(nSub1In,nSub2In,nPairLocal)
        real(8), intent(in) :: dYiIn(nSub1In+nSub2In), dYiGradIn(nSub1In+nSub2In,nPairLocal)
        real(8), intent(inout) :: dMuOut(nPairLocal), dJacOut(nPairLocal,nPairLocal)
        integer, intent(out) :: iInfoOut

        integer :: iParam, a, b, x, y, p, q, r, d, w, xx, yy, iBlock
        integer :: iBlockLocal, i, j, k, l, iQuad2
        logical :: lAsymmetric1(nSub1In), lAsymmetric2(nSub1In)
        real(8) :: dXi1, dXi2, dXiDen, dXi1Grad(nPairLocal), dXi2Grad(nPairLocal), dXiDenGrad(nPairLocal)
        real(8) :: dChi1, dChi2, dChiDen, dChi1Grad(nPairLocal), dChi2Grad(nPairLocal), dChiDenGrad(nPairLocal)
        real(8) :: dYdk, dYik, dYjk, dYdkGrad(nPairLocal), dYikGrad(nPairLocal), dYjkGrad(nPairLocal)
        real(8) :: dFactor, dFactorGrad(nPairLocal), dGex, dGexGrad(nPairLocal)
        real(8) :: dDgexBase, dDgexBaseGrad(nPairLocal), dDgex, dDgexGrad(nPairLocal)
        real(8) :: dTemp, dTempGrad(nPairLocal)
        real(8) :: dTFDG, dTFDGGrad(nPairLocal), dParameter
        real(8) :: dQuotientTemp, dQuotientGradTemp(nPairLocal)

        iInfoOut = 0

        do iParam = nParamPhase(iSolnIndexIn-1) + 1, nParamPhase(iSolnIndexIn)
            if (dExcessGibbsParam(iParam) == 0D0) cycle
            if ((TRIM(cRegularParam(iParam)) /= 'Q').AND.(TRIM(cRegularParam(iParam)) /= 'G')) then
                iInfoOut = 4
                return
            end if

            a = iRegularParam(iParam,2)
            b = iRegularParam(iParam,3)
            xx = iRegularParam(iParam,4)
            yy = iRegularParam(iParam,5)
            x = xx - nSub1In
            y = yy - nSub1In
            p = iRegularParam(iParam,6)
            q = iRegularParam(iParam,7)
            r = iRegularParam(iParam,8)
            d = iRegularParam(iParam,10)
            w = iRegularParam(iParam,11)
            if ((w /= 0).OR.(x /= y).OR.(a == b).OR.(x < 1).OR.(x > nSub2In)) then
                iInfoOut = 4
                return
            end if

            iBlock = (x - 1) * (nSub1In * (nSub1In + 1) / 2)
            if (a < b) then
                iBlock = iBlock + nSub1In + a + ((b - 2) * (b - 1) / 2)
            else
                iInfoOut = 4
                return
            end if
            iBlockLocal = iBlock
            if ((iBlockLocal < 1).OR.(iBlockLocal > nPairLocal)) then
                iInfoOut = 2
                return
            end if

            lAsymmetric1 = .FALSE.
            lAsymmetric2 = .FALSE.
            lAsymmetric1(a) = .TRUE.
            lAsymmetric2(b) = .TRUE.
            if (iChemicalGroup(iSPIIn,1,a) /= iChemicalGroup(iSPIIn,1,b)) then
                do i = 1, nSub1In
                    if (iChemicalGroup(iSPIIn,1,i) == iChemicalGroup(iSPIIn,1,a)) then
                        lAsymmetric1(i) = .TRUE.
                    else if (iChemicalGroup(iSPIIn,1,i) == iChemicalGroup(iSPIIn,1,b)) then
                        lAsymmetric2(i) = .TRUE.
                    end if
                end do
            end if

            call BuildAsymmetricFractions(iSPIIn, iFirstIn, nPairLocal, nSub1In, x, lAsymmetric1, lAsymmetric2, &
                dXi1, dXi1Grad, dXi2, dXi2Grad, dXiDen, dXiDenGrad, &
                dChi1, dChi1Grad, dChi2, dChi2Grad, dChiDen, dChiDenGrad)
            call QuotientGradient(nPairLocal, dChi1, dChi1Grad, dChiDen, dChiDenGrad, &
                dQuotientTemp, dQuotientGradTemp)
            dChi1 = dQuotientTemp
            dChi1Grad = dQuotientGradTemp
            call QuotientGradient(nPairLocal, dChi2, dChi2Grad, dChiDen, dChiDenGrad, &
                dQuotientTemp, dQuotientGradTemp)
            dChi2 = dQuotientTemp
            dChi2Grad = dQuotientGradTemp

            call BuildTernaryFactor(iSPIIn, iFirstIn, nPairLocal, nSub1In, x, a, b, d, r, &
                lAsymmetric1, lAsymmetric2, dXi1, dXi1Grad, dXi2, dXi2Grad, &
                dFactor, dFactorGrad, dYdk, dYdkGrad, dYik, dYikGrad, dYjk, dYjkGrad)

            dParameter = dExcessGibbsParam(iParam)
            if (TRIM(cRegularParam(iParam)) == 'Q') then
                call PowerProductQ(nPairLocal, dParameter, p, q, dXi1, dXi1Grad, dXi2, dXi2Grad, &
                    dXiDen, dXiDenGrad, dFactor, dFactorGrad, dGex, dGexGrad)
                call DivideByLinear(nPairLocal, -dGex * DBLE(p + q), &
                    -DBLE(p + q) * dGexGrad, dXiDen, dXiDenGrad, dDgexBase, dDgexBaseGrad)
            else
                call PowerProductG(nPairLocal, dParameter, p, q, dChi1, dChi1Grad, dChi2, dChi2Grad, &
                    dFactor, dFactorGrad, dGex, dGexGrad)
                call DivideByLinear(nPairLocal, -dGex * DBLE(p + q), &
                    -DBLE(p + q) * dGexGrad, dChiDen, dChiDenGrad, dDgexBase, dDgexBaseGrad)
            end if

            dMuOut(iBlockLocal) = dMuOut(iBlockLocal) + 0.5D0 * dGex
            dJacOut(iBlockLocal,:) = dJacOut(iBlockLocal,:) + 0.5D0 * dGexGrad

            do iQuad2 = 1, nPairLocal
                i = iPairID(iSPIIn,iQuad2,1)
                j = iPairID(iSPIIn,iQuad2,2)
                k = iPairID(iSPIIn,iQuad2,3) - nSub1In
                l = iPairID(iSPIIn,iQuad2,4) - nSub1In

                call BuildTernaryFactorDerivative(iSPIIn, nPairLocal, nSub1In, x, a, b, d, r, &
                    i, j, k, l, lAsymmetric1, lAsymmetric2, dXi1, dXi1Grad, dXi2, dXi2Grad, &
                    dYdk, dYdkGrad, dYik, dYikGrad, dYjk, dYjkGrad, dTFDG, dTFDGGrad)

                if (TRIM(cRegularParam(iParam)) == 'Q') then
                    call AccumulateQPairDerivative(nPairLocal, nSub1In, i, j, k, l, x, p, q, &
                        lAsymmetric1, lAsymmetric2, dXi1, dXi1Grad, dXi2, dXi2Grad, &
                        dGex, dGexGrad, dDgexBase, dDgexBaseGrad, dDgex, dDgexGrad)
                else
                    call AccumulateGPairDerivative(nPairLocal, i, j, k, l, x, p, q, &
                        TRIM(cSolnPhaseType(iSolnIndexIn)) == 'SUBQ', lAsymmetric1, lAsymmetric2, &
                        dChi1, dChi1Grad, dChi2, dChi2Grad, dChiDen, dChiDenGrad, &
                        dGex, dGexGrad, dDgexBase, dDgexBaseGrad, dDgex, dDgexGrad)
                end if

                call MultiplyGradient(nPairLocal, dGex, dGexGrad, dTFDG, dTFDGGrad, dTemp, dTempGrad)
                dDgex = dDgex + dTemp
                dDgexGrad = dDgexGrad + dTempGrad

                dMuOut(iQuad2) = dMuOut(iQuad2) + 0.5D0 * dMolFraction(iFirstIn + iBlockLocal - 1) * dDgex
                dJacOut(iQuad2,:) = dJacOut(iQuad2,:) + &
                    0.5D0 * dMolFraction(iFirstIn + iBlockLocal - 1) * dDgexGrad
                dJacOut(iQuad2,iBlockLocal) = dJacOut(iQuad2,iBlockLocal) + 0.5D0 * dDgex
            end do
        end do

        return

    end subroutine AccumulateQExcessJacobian

    subroutine AddLogTerm(nLocal, dCoeff, dValue, dGradient, dMuInOut, dRowInOut)

        implicit none

        integer, intent(in) :: nLocal
        real(8), intent(in) :: dCoeff, dValue, dGradient(nLocal)
        real(8), intent(inout) :: dMuInOut, dRowInOut(nLocal)
        real(8) :: dSafe

        dSafe = DMAX1(dValue, 1D-75)
        dMuInOut = dMuInOut + dCoeff * DLOG(dSafe)
        if (dValue > 1D-75) dRowInOut = dRowInOut + dCoeff * dGradient / dSafe

        return

    end subroutine AddLogTerm

    subroutine QuotientGradient(nLocal, dA, dAGrad, dB, dBGrad, dOut, dOutGrad)

        implicit none

        integer, intent(in) :: nLocal
        real(8), intent(in) :: dA, dAGrad(nLocal), dB, dBGrad(nLocal)
        real(8), intent(out) :: dOut, dOutGrad(nLocal)

        dOut = dA / dB
        dOutGrad = (dAGrad * dB - dA * dBGrad) / (dB * dB)

        return

    end subroutine QuotientGradient

    subroutine MultiplyGradient(nLocal, dA, dAGrad, dB, dBGrad, dOut, dOutGrad)

        implicit none

        integer, intent(in) :: nLocal
        real(8), intent(in) :: dA, dAGrad(nLocal), dB, dBGrad(nLocal)
        real(8), intent(out) :: dOut, dOutGrad(nLocal)

        dOut = dA * dB
        dOutGrad = dAGrad * dB + dA * dBGrad

        return

    end subroutine MultiplyGradient

    subroutine DivideByLinear(nLocal, dA, dAGrad, dB, dBGrad, dOut, dOutGrad)

        implicit none

        integer, intent(in) :: nLocal
        real(8), intent(in) :: dA, dAGrad(nLocal), dB, dBGrad(nLocal)
        real(8), intent(out) :: dOut, dOutGrad(nLocal)

        dOut = dA / dB
        dOutGrad = (dAGrad * dB - dA * dBGrad) / (dB * dB)

        return

    end subroutine DivideByLinear

    subroutine PowGradient(nLocal, dA, dAGrad, iPower, dOut, dOutGrad)

        implicit none

        integer, intent(in) :: nLocal, iPower
        real(8), intent(in) :: dA, dAGrad(nLocal)
        real(8), intent(out) :: dOut, dOutGrad(nLocal)

        if (iPower == 0) then
            dOut = 1D0
            dOutGrad = 0D0
        else
            dOut = dA**iPower
            dOutGrad = DBLE(iPower) * dA**(iPower - 1) * dAGrad
        end if

        return

    end subroutine PowGradient

    subroutine PowerProductQ(nLocal, dParameter, iP, iQ, dXi1, dXi1Grad, dXi2, dXi2Grad, &
        dXiDen, dXiDenGrad, dFactor, dFactorGrad, dOut, dOutGrad)

        implicit none

        integer, intent(in) :: nLocal, iP, iQ
        real(8), intent(in) :: dParameter, dXi1, dXi1Grad(nLocal), dXi2, dXi2Grad(nLocal)
        real(8), intent(in) :: dXiDen, dXiDenGrad(nLocal), dFactor, dFactorGrad(nLocal)
        real(8), intent(out) :: dOut, dOutGrad(nLocal)
        real(8) :: dA, dB, dC, dAB, dABC
        real(8) :: dAGrad(nLocal), dBGrad(nLocal), dCGrad(nLocal), dABGrad(nLocal), dABCGrad(nLocal)

        call PowGradient(nLocal, dXi1, dXi1Grad, iP, dA, dAGrad)
        call PowGradient(nLocal, dXi2, dXi2Grad, iQ, dB, dBGrad)
        call PowGradient(nLocal, dXiDen, dXiDenGrad, iP + iQ, dC, dCGrad)
        call MultiplyGradient(nLocal, dA, dAGrad, dB, dBGrad, dAB, dABGrad)
        call DivideByLinear(nLocal, dAB, dABGrad, dC, dCGrad, dABC, dABCGrad)
        call MultiplyGradient(nLocal, dABC, dABCGrad, dFactor, dFactorGrad, dOut, dOutGrad)
        dOut = dParameter * dOut
        dOutGrad = dParameter * dOutGrad

        return

    end subroutine PowerProductQ

    subroutine PowerProductG(nLocal, dParameter, iP, iQ, dChi1, dChi1Grad, dChi2, dChi2Grad, &
        dFactor, dFactorGrad, dOut, dOutGrad)

        implicit none

        integer, intent(in) :: nLocal, iP, iQ
        real(8), intent(in) :: dParameter, dChi1, dChi1Grad(nLocal), dChi2, dChi2Grad(nLocal)
        real(8), intent(in) :: dFactor, dFactorGrad(nLocal)
        real(8), intent(out) :: dOut, dOutGrad(nLocal)
        real(8) :: dA, dB, dAB
        real(8) :: dAGrad(nLocal), dBGrad(nLocal), dABGrad(nLocal)

        call PowGradient(nLocal, dChi1, dChi1Grad, iP, dA, dAGrad)
        call PowGradient(nLocal, dChi2, dChi2Grad, iQ, dB, dBGrad)
        call MultiplyGradient(nLocal, dA, dAGrad, dB, dBGrad, dAB, dABGrad)
        call MultiplyGradient(nLocal, dAB, dABGrad, dFactor, dFactorGrad, dOut, dOutGrad)
        dOut = dParameter * dOut
        dOutGrad = dParameter * dOutGrad

        return

    end subroutine PowerProductG

    subroutine AccumulateQPairDerivative(nLocal, nSub1In, iQuadA, iQuadB, iQuadX, iQuadY, iAnion, iP, iQ, &
        lAsym1, lAsym2, dXi1, dXi1Grad, dXi2, dXi2Grad, dGex, dGexGrad, &
        dDgexBase, dDgexBaseGrad, dDgexOut, dDgexGradOut)

        implicit none

        integer, intent(in) :: nLocal, nSub1In, iQuadA, iQuadB, iQuadX, iQuadY, iAnion, iP, iQ
        logical, intent(in) :: lAsym1(nSub1In), lAsym2(nSub1In)
        real(8), intent(in) :: dXi1, dXi1Grad(nLocal), dXi2, dXi2Grad(nLocal)
        real(8), intent(in) :: dGex, dGexGrad(nLocal), dDgexBase, dDgexBaseGrad(nLocal)
        real(8), intent(out) :: dDgexOut, dDgexGradOut(nLocal)

        integer :: ii, nMatch1, nMatch2
        real(8) :: dTerm, dTermGrad(nLocal)

        dDgexOut = 0D0
        dDgexGradOut = 0D0
        nMatch1 = 0
        nMatch2 = 0

        do ii = 1, nSub1In
            if (lAsym1(ii)) then
                if (ii == iQuadA .AND. iAnion == iQuadX) nMatch1 = nMatch1 + 1
                if (ii == iQuadA .AND. iAnion == iQuadY) nMatch1 = nMatch1 + 1
                if (ii == iQuadB .AND. iAnion == iQuadX) nMatch1 = nMatch1 + 1
                if (ii == iQuadB .AND. iAnion == iQuadY) nMatch1 = nMatch1 + 1
            end if
            if (lAsym2(ii)) then
                if (ii == iQuadA .AND. iAnion == iQuadX) nMatch2 = nMatch2 + 1
                if (ii == iQuadA .AND. iAnion == iQuadY) nMatch2 = nMatch2 + 1
                if (ii == iQuadB .AND. iAnion == iQuadX) nMatch2 = nMatch2 + 1
                if (ii == iQuadB .AND. iAnion == iQuadY) nMatch2 = nMatch2 + 1
            end if
        end do

        if (nMatch1 > 0) then
            call DivideByLinear(nLocal, dGex * DBLE(iP), dGexGrad * DBLE(iP), &
                4D0 * dXi1, 4D0 * dXi1Grad, dTerm, dTermGrad)
            dDgexOut = dDgexOut + DBLE(nMatch1) * (dDgexBase / 4D0 + dTerm)
            dDgexGradOut = dDgexGradOut + DBLE(nMatch1) * (dDgexBaseGrad / 4D0 + dTermGrad)
        end if
        if (nMatch2 > 0) then
            call DivideByLinear(nLocal, dGex * DBLE(iQ), dGexGrad * DBLE(iQ), &
                4D0 * dXi2, 4D0 * dXi2Grad, dTerm, dTermGrad)
            dDgexOut = dDgexOut + DBLE(nMatch2) * (dDgexBase / 4D0 + dTerm)
            dDgexGradOut = dDgexGradOut + DBLE(nMatch2) * (dDgexBaseGrad / 4D0 + dTermGrad)
        end if

        return

    end subroutine AccumulateQPairDerivative

    subroutine AccumulateGPairDerivative(nLocal, iQuadA, iQuadB, iQuadX, iQuadY, iAnion, iP, iQ, &
        lUseSubQHalf, lAsym1, lAsym2, dChi1, dChi1Grad, dChi2, dChi2Grad, dChiDen, dChiDenGrad, &
        dGex, dGexGrad, dDgexBase, dDgexBaseGrad, dDgexOut, dDgexGradOut)

        implicit none

        integer, intent(in) :: nLocal, iQuadA, iQuadB, iQuadX, iQuadY, iAnion, iP, iQ
        logical, intent(in) :: lUseSubQHalf
        logical, intent(in) :: lAsym1(:), lAsym2(:)
        real(8), intent(in) :: dChi1, dChi1Grad(nLocal), dChi2, dChi2Grad(nLocal)
        real(8), intent(in) :: dChiDen, dChiDenGrad(nLocal), dGex, dGexGrad(nLocal)
        real(8), intent(in) :: dDgexBase, dDgexBaseGrad(nLocal)
        real(8), intent(out) :: dDgexOut, dDgexGradOut(nLocal)

        real(8) :: dChiFactor, dDen, dDenGrad(nLocal), dTerm, dTermGrad(nLocal)

        dDgexOut = 0D0
        dDgexGradOut = 0D0

        dChiFactor = 0D0
        if ((iAnion == iQuadX).AND.(iAnion == iQuadY)) then
            dChiFactor = 1D0
        else if (lUseSubQHalf.AND.((iAnion == iQuadX).OR.(iAnion == iQuadY))) then
            dChiFactor = 0.5D0
        end if
        if (dChiFactor == 0D0) return

        if (lAsym1(iQuadA).AND.lAsym1(iQuadB)) then
            call MultiplyGradient(nLocal, dChi1, dChi1Grad, dChiDen, dChiDenGrad, dDen, dDenGrad)
            call DivideByLinear(nLocal, dChiFactor * dGex * DBLE(iP), &
                dChiFactor * dGexGrad * DBLE(iP), dDen, dDenGrad, dTerm, dTermGrad)
            dDgexOut = dDgexOut + dTerm
            dDgexGradOut = dDgexGradOut + dTermGrad
        end if

        if (lAsym2(iQuadA).AND.lAsym2(iQuadB)) then
            call MultiplyGradient(nLocal, dChi2, dChi2Grad, dChiDen, dChiDenGrad, dDen, dDenGrad)
            call DivideByLinear(nLocal, dChiFactor * dGex * DBLE(iQ), &
                dChiFactor * dGexGrad * DBLE(iQ), dDen, dDenGrad, dTerm, dTermGrad)
            dDgexOut = dDgexOut + dTerm
            dDgexGradOut = dDgexGradOut + dTermGrad
        end if

        if ((lAsym1(iQuadA).OR.lAsym2(iQuadA)).AND.(lAsym1(iQuadB).OR.lAsym2(iQuadB))) then
            dDgexOut = dDgexOut + dChiFactor * dDgexBase
            dDgexGradOut = dDgexGradOut + dChiFactor * dDgexBaseGrad
        end if

        return

    end subroutine AccumulateGPairDerivative

    subroutine BuildAsymmetricFractions(iSPIIn, iFirstIn, nPairLocal, nSub1In, iAnion, lAsym1, lAsym2, &
        dXi1Out, dXi1GradOut, dXi2Out, dXi2GradOut, dXiDenOut, dXiDenGradOut, &
        dChi1Out, dChi1GradOut, dChi2Out, dChi2GradOut, dChiDenOut, dChiDenGradOut)

        implicit none

        integer, intent(in) :: iSPIIn, iFirstIn, nPairLocal, nSub1In, iAnion
        logical, intent(in) :: lAsym1(nSub1In), lAsym2(nSub1In)
        real(8), intent(out) :: dXi1Out, dXi1GradOut(nPairLocal), dXi2Out, dXi2GradOut(nPairLocal)
        real(8), intent(out) :: dXiDenOut, dXiDenGradOut(nPairLocal)
        real(8), intent(out) :: dChi1Out, dChi1GradOut(nPairLocal), dChi2Out, dChi2GradOut(nPairLocal)
        real(8), intent(out) :: dChiDenOut, dChiDenGradOut(nPairLocal)

        integer :: i, j, k, l, iQuad, kk, ll
        real(8) :: dWeight

        dXi1Out = 0D0
        dXi2Out = 0D0
        dXiDenOut = 0D0
        dChi1Out = 0D0
        dChi2Out = 0D0
        dChiDenOut = 0D0
        dXi1GradOut = 0D0
        dXi2GradOut = 0D0
        dXiDenGradOut = 0D0
        dChi1GradOut = 0D0
        dChi2GradOut = 0D0
        dChiDenGradOut = 0D0

        do iQuad = 1, nPairLocal
            i = iPairID(iSPIIn,iQuad,1)
            j = iPairID(iSPIIn,iQuad,2)
            kk = iPairID(iSPIIn,iQuad,3) - nSub1In
            ll = iPairID(iSPIIn,iQuad,4) - nSub1In
            l = iFirstIn + iQuad - 1
            if (lAsym1(i).AND.lAsym1(j)) then
                if ((iAnion == kk).AND.(iAnion == ll)) then
                    dChi1Out = dChi1Out + dMolFraction(l)
                    dChi1GradOut(iQuad) = dChi1GradOut(iQuad) + 1D0
                else if ((iAnion == kk).OR.(iAnion == ll)) then
                    dChi1Out = dChi1Out + 0.5D0 * dMolFraction(l)
                    dChi1GradOut(iQuad) = dChi1GradOut(iQuad) + 0.5D0
                end if
            end if
            if (lAsym2(i).AND.lAsym2(j)) then
                if ((iAnion == kk).AND.(iAnion == ll)) then
                    dChi2Out = dChi2Out + dMolFraction(l)
                    dChi2GradOut(iQuad) = dChi2GradOut(iQuad) + 1D0
                else if ((iAnion == kk).OR.(iAnion == ll)) then
                    dChi2Out = dChi2Out + 0.5D0 * dMolFraction(l)
                    dChi2GradOut(iQuad) = dChi2GradOut(iQuad) + 0.5D0
                end if
            end if
            if ((lAsym1(i).OR.lAsym2(i)).AND.(lAsym1(j).OR.lAsym2(j))) then
                if ((iAnion == kk).AND.(iAnion == ll)) then
                    dChiDenOut = dChiDenOut + dMolFraction(l)
                    dChiDenGradOut(iQuad) = dChiDenGradOut(iQuad) + 1D0
                else if ((iAnion == kk).OR.(iAnion == ll)) then
                    dChiDenOut = dChiDenOut + 0.5D0 * dMolFraction(l)
                    dChiDenGradOut(iQuad) = dChiDenGradOut(iQuad) + 0.5D0
                end if
            end if

            do k = 1, nSub1In
                dWeight = 0D0
                if (k == i .AND. iAnion == kk) dWeight = dWeight + 0.25D0
                if (k == i .AND. iAnion == ll) dWeight = dWeight + 0.25D0
                if (k == j .AND. iAnion == kk) dWeight = dWeight + 0.25D0
                if (k == j .AND. iAnion == ll) dWeight = dWeight + 0.25D0
                if (lAsym1(k)) then
                    dXi1Out = dXi1Out + dWeight * dMolFraction(l)
                    dXi1GradOut(iQuad) = dXi1GradOut(iQuad) + dWeight
                end if
                if (lAsym2(k)) then
                    dXi2Out = dXi2Out + dWeight * dMolFraction(l)
                    dXi2GradOut(iQuad) = dXi2GradOut(iQuad) + dWeight
                end if
            end do
        end do
        dXiDenOut = dXi1Out + dXi2Out
        dXiDenGradOut = dXi1GradOut + dXi2GradOut

        return

    end subroutine BuildAsymmetricFractions

    subroutine BuildTernaryFactor(iSPIIn, iFirstIn, nPairLocal, nSub1In, iAnion, iA, iB, iD, iR, &
        lAsym1, lAsym2, dXi1, dXi1Grad, dXi2, dXi2Grad, dFactorOut, dFactorGradOut, &
        dYdkOut, dYdkGradOut, dYikOut, dYikGradOut, dYjkOut, dYjkGradOut)

        implicit none

        integer, intent(in) :: iSPIIn, iFirstIn, nPairLocal, nSub1In, iAnion, iA, iB, iD, iR
        logical, intent(in) :: lAsym1(nSub1In), lAsym2(nSub1In)
        real(8), intent(in) :: dXi1, dXi1Grad(nPairLocal), dXi2, dXi2Grad(nPairLocal)
        real(8), intent(out) :: dFactorOut, dFactorGradOut(nPairLocal)
        real(8), intent(out) :: dYdkOut, dYdkGradOut(nPairLocal), dYikOut, dYikGradOut(nPairLocal)
        real(8), intent(out) :: dYjkOut, dYjkGradOut(nPairLocal)
        real(8) :: dRatio, dRatioGrad(nPairLocal), dOneMinus, dOneMinusGrad(nPairLocal)
        real(8) :: dPower, dPowerGrad(nPairLocal)

        call BuildYFraction(iSPIIn, iFirstIn, nPairLocal, nSub1In, iD, iAnion, dYdkOut, dYdkGradOut)
        call BuildYFraction(iSPIIn, iFirstIn, nPairLocal, nSub1In, iA, iAnion, dYikOut, dYikGradOut)
        call BuildYFraction(iSPIIn, iFirstIn, nPairLocal, nSub1In, iB, iAnion, dYjkOut, dYjkGradOut)

        dFactorOut = 1D0
        dFactorGradOut = 0D0
        if (iD <= 0) return

        if (lAsym2(iD)) then
            call DivideByLinear(nPairLocal, dYdkOut, dYdkGradOut, dXi2, dXi2Grad, dRatio, dRatioGrad)
            call DivideByLinear(nPairLocal, dYjkOut, dYjkGradOut, dXi2, dXi2Grad, dOneMinus, dOneMinusGrad)
            dOneMinus = 1D0 - dOneMinus
            dOneMinusGrad = -dOneMinusGrad
            call PowGradient(nPairLocal, dOneMinus, dOneMinusGrad, iR - 1, dPower, dPowerGrad)
            call MultiplyGradient(nPairLocal, dRatio, dRatioGrad, dPower, dPowerGrad, dFactorOut, dFactorGradOut)
        else if (lAsym1(iD)) then
            call DivideByLinear(nPairLocal, dYdkOut, dYdkGradOut, dXi1, dXi1Grad, dRatio, dRatioGrad)
            call DivideByLinear(nPairLocal, dYikOut, dYikGradOut, dXi1, dXi1Grad, dOneMinus, dOneMinusGrad)
            dOneMinus = 1D0 - dOneMinus
            dOneMinusGrad = -dOneMinusGrad
            call PowGradient(nPairLocal, dOneMinus, dOneMinusGrad, iR - 1, dPower, dPowerGrad)
            call MultiplyGradient(nPairLocal, dRatio, dRatioGrad, dPower, dPowerGrad, dFactorOut, dFactorGradOut)
        else
            dOneMinus = 1D0 - dXi1 - dXi2
            dOneMinusGrad = -dXi1Grad - dXi2Grad
            call PowGradient(nPairLocal, dOneMinus, dOneMinusGrad, iR - 1, dPower, dPowerGrad)
            call MultiplyGradient(nPairLocal, dYdkOut, dYdkGradOut, dPower, dPowerGrad, dFactorOut, dFactorGradOut)
        end if

        return

    end subroutine BuildTernaryFactor

    subroutine BuildYFraction(iSPIIn, iFirstIn, nPairLocal, nSub1In, iCation, iAnion, dValueOut, dGradOut)

        implicit none

        integer, intent(in) :: iSPIIn, iFirstIn, nPairLocal, nSub1In, iCation, iAnion
        real(8), intent(out) :: dValueOut, dGradOut(nPairLocal)
        integer :: k, l, i, j, x, y
        real(8) :: dWeight

        dValueOut = 0D0
        dGradOut = 0D0
        if (iCation <= 0) return
        do k = 1, nPairLocal
            i = iPairID(iSPIIn,k,1)
            j = iPairID(iSPIIn,k,2)
            x = iPairID(iSPIIn,k,3) - nSub1In
            y = iPairID(iSPIIn,k,4) - nSub1In
            l = iFirstIn + k - 1
            dWeight = 0D0
            if (iCation == i .AND. iAnion == x) dWeight = dWeight + 0.25D0
            if (iCation == i .AND. iAnion == y) dWeight = dWeight + 0.25D0
            if (iCation == j .AND. iAnion == x) dWeight = dWeight + 0.25D0
            if (iCation == j .AND. iAnion == y) dWeight = dWeight + 0.25D0
            dValueOut = dValueOut + dWeight * dMolFraction(l)
            dGradOut(k) = dGradOut(k) + dWeight
        end do

        return

    end subroutine BuildYFraction

    subroutine BuildTernaryFactorDerivative(iSPIIn, nPairLocal, nSub1In, iAnion, iA, iB, iD, iR, &
        iQuadA, iQuadB, iQuadX, iQuadY, lAsym1, lAsym2, dXi1, dXi1Grad, dXi2, dXi2Grad, &
        dYdk, dYdkGrad, dYik, dYikGrad, dYjk, dYjkGrad, dValueOut, dGradOut)

        implicit none

        integer, intent(in) :: iSPIIn, nPairLocal, nSub1In, iAnion, iA, iB, iD, iR
        integer, intent(in) :: iQuadA, iQuadB, iQuadX, iQuadY
        logical, intent(in) :: lAsym1(nSub1In), lAsym2(nSub1In)
        real(8), intent(in) :: dXi1, dXi1Grad(nPairLocal), dXi2, dXi2Grad(nPairLocal)
        real(8), intent(in) :: dYdk, dYdkGrad(nPairLocal), dYik, dYikGrad(nPairLocal), dYjk, dYjkGrad(nPairLocal)
        real(8), intent(out) :: dValueOut, dGradOut(nPairLocal)
        integer :: e, nA, nX
        real(8) :: dConst1, dConst2, dConst3, dSum1, dSum2
        real(8) :: dTerm, dTermGrad(nPairLocal), dRatio, dRatioGrad(nPairLocal)
        real(8) :: dDen, dDenGrad(nPairLocal), dNum, dNumGrad(nPairLocal)

        dValueOut = 0D0
        dGradOut = 0D0
        if (iD <= 0) return

        nX = 0
        if (iQuadX == iAnion) nX = nX + 1
        if (iQuadY == iAnion) nX = nX + 1

        if (lAsym2(iD)) then
            nA = 0
            if (iQuadA == iD) nA = nA + 1
            if (iQuadB == iD) nA = nA + 1
            dConst1 = DBLE(nA * nX) / 4D0
            call DivideByLinear(nPairLocal, dConst1, 0D0 * dXi2Grad, dYdk, dYdkGrad, dTerm, dTermGrad)
            dValueOut = dValueOut + dTerm
            dGradOut = dGradOut + dTermGrad

            dSum2 = 0D0
            do e = 1, nSub1In
                nA = 0
                if (iQuadA == e) nA = nA + 1
                if (iQuadB == e) nA = nA + 1
                if (lAsym2(e)) dSum2 = dSum2 + DBLE(nA * nX) / (4D0 * dXi2)
            end do
            dValueOut = dValueOut - dSum2
            dGradOut = dGradOut + dSum2 * dXi2Grad / dXi2

            nA = 0
            if (iQuadA == iB) nA = nA + 1
            if (iQuadB == iB) nA = nA + 1
            dConst2 = DBLE(nA * nX) / 4D0
            dNum = dConst2 - dYjk * dSum2
            dNumGrad = -dYjkGrad * dSum2 + dYjk * dSum2 * dXi2Grad / dXi2
            dDen = dXi2 * (1D0 - dYjk / dXi2)
            dDenGrad = dXi2Grad * (1D0 - dYjk / dXi2) - dXi2 * &
                ((dYjkGrad * dXi2 - dYjk * dXi2Grad) / (dXi2 * dXi2))
            call DivideByLinear(nPairLocal, dNum, dNumGrad, dDen, dDenGrad, dRatio, dRatioGrad)
            dValueOut = dValueOut - DBLE(iR - 1) * dRatio
            dGradOut = dGradOut - DBLE(iR - 1) * dRatioGrad
        else if (lAsym1(iD)) then
            nA = 0
            if (iQuadA == iD) nA = nA + 1
            if (iQuadB == iD) nA = nA + 1
            dConst1 = DBLE(nA * nX) / 4D0
            call DivideByLinear(nPairLocal, dConst1, 0D0 * dXi1Grad, dYdk, dYdkGrad, dTerm, dTermGrad)
            dValueOut = dValueOut + dTerm
            dGradOut = dGradOut + dTermGrad

            dSum1 = 0D0
            do e = 1, nSub1In
                nA = 0
                if (iQuadA == e) nA = nA + 1
                if (iQuadB == e) nA = nA + 1
                if (lAsym1(e)) dSum1 = dSum1 + DBLE(nA * nX) / (4D0 * dXi1)
            end do
            dValueOut = dValueOut - dSum1
            dGradOut = dGradOut + dSum1 * dXi1Grad / dXi1

            nA = 0
            if (iQuadA == iA) nA = nA + 1
            if (iQuadB == iA) nA = nA + 1
            dConst2 = DBLE(nA * nX) / 4D0
            dNum = dConst2 - dYik * dSum1
            dNumGrad = -dYikGrad * dSum1 + dYik * dSum1 * dXi1Grad / dXi1
            dDen = dXi1 * (1D0 - dYik / dXi1)
            dDenGrad = dXi1Grad * (1D0 - dYik / dXi1) - dXi1 * &
                ((dYikGrad * dXi1 - dYik * dXi1Grad) / (dXi1 * dXi1))
            call DivideByLinear(nPairLocal, dNum, dNumGrad, dDen, dDenGrad, dRatio, dRatioGrad)
            dValueOut = dValueOut - DBLE(iR - 1) * dRatio
            dGradOut = dGradOut - DBLE(iR - 1) * dRatioGrad
        else
            dValueOut = -DBLE(iR)
            nA = 0
            if (iQuadA == iD) nA = nA + 1
            if (iQuadB == iD) nA = nA + 1
            dConst1 = DBLE(nA * nX) / 4D0
            call DivideByLinear(nPairLocal, dConst1, 0D0 * dXi1Grad, dYdk, dYdkGrad, dTerm, dTermGrad)
            dValueOut = dValueOut + dTerm
            dGradOut = dGradOut + dTermGrad

            dSum1 = 0D0
            dSum2 = 0D0
            do e = 1, nSub1In
                nA = 0
                if (iQuadA == e) nA = nA + 1
                if (iQuadB == e) nA = nA + 1
                dConst3 = DBLE(nA * nX) / 4D0
                if (lAsym1(e)) dSum1 = dSum1 + dConst3
                if (lAsym2(e)) dSum2 = dSum2 + dConst3
            end do
            dNum = 1D0 - dSum1 - dSum2
            dNumGrad = 0D0
            dDen = 1D0 - dXi1 - dXi2
            dDenGrad = -dXi1Grad - dXi2Grad
            call DivideByLinear(nPairLocal, dNum, dNumGrad, dDen, dDenGrad, dRatio, dRatioGrad)
            dValueOut = dValueOut + DBLE(iR - 1) * dRatio
            dGradOut = dGradOut + DBLE(iR - 1) * dRatioGrad
        end if

        return

    end subroutine BuildTernaryFactorDerivative

end subroutine CompHessianSUBG
