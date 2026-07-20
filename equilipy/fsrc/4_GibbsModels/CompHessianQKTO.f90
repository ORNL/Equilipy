!> \brief Compute analytical QKTO chemical-potential Jacobians.
!!
!! \details This diagnostic/minimizer routine evaluates one-sublattice QKTO
!! model-species chemical potentials and their analytical derivatives with
!! respect to phase-local mole fractions.  The returned Hessian is the
!! redundant composition-space Jacobian \f$\partial\mu_i/\partial x_j\f$;
!! callers should project it before solving constrained Newton systems.
!!
!! \param[in]  iSolnIndex                 Absolute solution phase index.
!! \param[in]  nSpeciesDim                Caller-provided species capacity.
!! \param[out] dSpeciesChemicalPotential  Chemical potentials in internal G/(RT) units.
!! \param[out] dSpeciesHessian            Composition Jacobian in internal G/(RT) units.
!! \param[out] nSpeciesOut                Number of active phase species written.
!! \param[out] iInfo                      Zero on success; nonzero for unsupported/undersized calls.



subroutine CompHessianQKTO(iSolnIndex, nSpeciesDim, dSpeciesChemicalPotential, dSpeciesHessian, &
    nSpeciesOut, iInfo)
!
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CompHessianQKTO.f90
    !> \brief   Compute analytical QKTO chemical-potential Jacobians.
    !> \author  S.Y. Kwon
    !> \date    Jun. 28, 2026
    !> \sa      CompExcessGibbsEnergyQKTO.f90
    !> \sa      PolyRegularQKTO.f90
    !> \sa      KohlerInterpolate.f90
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/20/2026      S.Y. Kwon           Added the analytical ideal and excess chemical-potential Jacobian for QKTO phases.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details This routine evaluates \f$\mu_i\f$ and
    !! \f$J_{ij}=\partial \mu_i/\partial x_j\f$ for a nonmagnetic QKTO
    !! one-sublattice phase.  Each Kohler-Toop excess term is treated as the
    !! scalar
    !! \f$g=L x_T^{k-\sum a_m}\prod_m x_m^{a_m}\f$,
    !! then scalar curvature is projected to partial chemical-potential
    !! curvature:
    !! \f$J_{ij}=K_{ij}-\sum_k x_k K_{kj}+\delta_{ij}/x_i\f$.
    !
    !
    ! Required input variables:
    ! =========================
    !
    !> \param[in] iSolnIndex   Absolute solution phase index.
    !> \param[in] nSpeciesDim  First/second dimension of output arrays.
    !
    ! cSolnPhaseType             Solution model type.  QKTO is accepted.
    ! dMolFraction               Current phase-local species mole fractions.
    ! dExcessGibbsParam          Dimensionless QKTO excess parameters.
    ! iRegularParam              QKTO parameter topology and exponents.
    !
    !
    ! Output/updated variables:
    ! =========================
    !
    !> \param[out] dSpeciesChemicalPotential  Species chemical potentials in G/(RT) units.
    !> \param[out] dSpeciesHessian            Redundant composition Jacobian dmu_i/dx_j in G/(RT) units.
    !> \param[out] nSpeciesOut                Number of active phase species written.
    !> \param[out] iInfo                      Zero on success. 1=unsupported phase, 2=bad topology,
    !!                                         3=caller arrays too small, 4=unsupported magnetic variant.
    !
    !
    ! Called subroutines/functions:
    ! =============================
    !
    ! AccumulateQKTOParameterScalar  Adds scalar excess value, gradient, and Hessian for one QKTO parameter.
    ! MultiplyLinearPower            Differentiates a linear factor raised to an integer power.
    ! ProjectScalarHessian           Converts scalar excess curvature to partial chemical-potential curvature.
    ! CompExcessGibbsEnergy          Recomputes runtime chemical potentials for the output vector.
    !
    !
    ! Primary callers:
    ! ================
    !
    ! Python diagnostics/devop tests  Validate analytical curvature against finite differences.
    ! SubMinNewton                    Candidate source for phase-local Newton directions.
    !
    !
    ! Numerical assumptions:
    ! ======================
    !
    ! - Returned values use the internal dimensionless Gibbs-energy convention, G/(RT).
    ! - dSpeciesHessian is redundant under the sum(x)=1 constraint and should be
    !   projected before inversion.
    ! - This routine supports the nonmagnetic QKTO excess curvature.  QKTOM returns
    !   iInfo=4 because magnetic curvature is not included here.
    ! - Validation points should stay inside the simplex.  The ideal term uses the
    !   same trace floor as the runtime Gibbs evaluator for numerical protection.
    !
    !-------------------------------------------------------------------------------------------------------------
    USE ModuleThermo

    implicit none

    integer, intent(in)  :: iSolnIndex, nSpeciesDim
    real(8), intent(out) :: dSpeciesChemicalPotential(nSpeciesDim)
    real(8), intent(out) :: dSpeciesHessian(nSpeciesDim,nSpeciesDim)
    integer, intent(out) :: nSpeciesOut, iInfo

    integer :: iFirst, iLast, iParam
    real(8), allocatable :: dScalarGradient(:), dScalarHessian(:,:)

    dSpeciesChemicalPotential = 0D0
    dSpeciesHessian = 0D0
    nSpeciesOut = 0
    iInfo = 0

    if ((iSolnIndex <= 0).OR.(iSolnIndex > nSolnPhasesSys)) then
        iInfo = 1
        return
    end if
    if (TRIM(cSolnPhaseType(iSolnIndex)) == 'QKTOM') then
        iInfo = 4
        return
    end if
    if (TRIM(cSolnPhaseType(iSolnIndex)) /= 'QKTO') then
        iInfo = 1
        return
    end if

    iFirst = nSpeciesPhase(iSolnIndex-1) + 1
    iLast = nSpeciesPhase(iSolnIndex)
    nSpeciesOut = iLast - iFirst + 1
    if (nSpeciesOut <= 0) then
        iInfo = 2
        return
    end if
    if (nSpeciesOut > nSpeciesDim) then
        iInfo = 3
        return
    end if

    allocate(dScalarGradient(nSpeciesOut), dScalarHessian(nSpeciesOut,nSpeciesOut))
    dScalarGradient = 0D0
    dScalarHessian = 0D0

    do iParam = nParamPhase(iSolnIndex-1) + 1, nParamPhase(iSolnIndex)
        call AccumulateQKTOParameterScalar(iFirst, nSpeciesOut, iParam, dScalarGradient, &
            dScalarHessian, iInfo)
        if (iInfo /= 0) goto 900
    end do

    call ProjectScalarHessian(iFirst, nSpeciesOut, dScalarHessian, dSpeciesHessian)

    call CompExcessGibbsEnergy(iSolnIndex)
    dSpeciesChemicalPotential(1:nSpeciesOut) = dChemicalPotential(iFirst:iLast)

900 continue
    if (allocated(dScalarGradient)) deallocate(dScalarGradient)
    if (allocated(dScalarHessian)) deallocate(dScalarHessian)

    return

contains

    subroutine AccumulateQKTOParameterScalar(iFirstIn, nLocal, iParamIn, dGradientOut, &
        dHessianOut, iInfoOut)

        implicit none

        integer, intent(in) :: iFirstIn, nLocal, iParamIn
        real(8), intent(inout) :: dGradientOut(nLocal), dHessianOut(nLocal,nLocal)
        integer, intent(out) :: iInfoOut

        integer :: i, iMixType, iSpecies, iExponent, iExponentSum, iKohlerExponent
        real(8) :: dTermValue, dMoleSum
        real(8) :: dTermGradient(nLocal), dTermHessian(nLocal,nLocal), dFactorGradient(nLocal)

        iInfoOut = 0
        dTermValue = dExcessGibbsParam(iParamIn)
        dTermGradient = 0D0
        dTermHessian = 0D0

        iMixType = iRegularParam(iParamIn,1)
        if ((iMixType < 2).OR.(iMixType > 4)) then
            iInfoOut = 2
            return
        end if

        iExponentSum = 0
        dMoleSum = 0D0
        dFactorGradient = 0D0
        do i = 1, iMixType
            iSpecies = iRegularParam(iParamIn,i+1)
            if ((iSpecies < 1).OR.(iSpecies > nLocal)) then
                iInfoOut = 2
                return
            end if
            iExponent = iRegularParam(iParamIn,iMixType+1+i)
            iExponentSum = iExponentSum + iExponent
            dMoleSum = dMoleSum + dMolFraction(iFirstIn+iSpecies-1)
            dFactorGradient(iSpecies) = dFactorGradient(iSpecies) + 1D0
        end do

        iKohlerExponent = iMixType - iExponentSum
        call MultiplyLinearPower(nLocal, dMoleSum, dFactorGradient, iKohlerExponent, &
            dTermValue, dTermGradient, dTermHessian)

        do i = 1, iMixType
            iSpecies = iRegularParam(iParamIn,i+1)
            iExponent = iRegularParam(iParamIn,iMixType+1+i)
            dFactorGradient = 0D0
            dFactorGradient(iSpecies) = 1D0
            call MultiplyLinearPower(nLocal, dMolFraction(iFirstIn+iSpecies-1), dFactorGradient, iExponent, &
                dTermValue, dTermGradient, dTermHessian)
        end do

        dGradientOut = dGradientOut + dTermGradient
        dHessianOut = dHessianOut + dTermHessian

        return

    end subroutine AccumulateQKTOParameterScalar

    subroutine MultiplyLinearPower(nLocal, dFactor, dFactorGradient, iExponent, dValue, &
        dGradient, dHessian)

        implicit none

        integer, intent(in) :: nLocal, iExponent
        real(8), intent(in) :: dFactor, dFactorGradient(nLocal)
        real(8), intent(inout) :: dValue, dGradient(nLocal), dHessian(nLocal,nLocal)

        integer :: i, j
        real(8) :: dOldValue, dPowerValue
        real(8) :: dOldGradient(nLocal), dPowerGradient(nLocal)
        real(8) :: dOldHessian(nLocal,nLocal), dPowerHessian(nLocal,nLocal)

        if (iExponent == 0) return

        dOldValue = dValue
        dOldGradient = dGradient
        dOldHessian = dHessian

        dPowerValue = dFactor ** iExponent
        dPowerGradient = DFLOAT(iExponent) * (dFactor ** (iExponent - 1)) * dFactorGradient
        dPowerHessian = 0D0
        if (iExponent /= 1) then
            do j = 1, nLocal
                do i = 1, nLocal
                    dPowerHessian(i,j) = DFLOAT(iExponent) * DFLOAT(iExponent - 1) * &
                        (dFactor ** (iExponent - 2)) * dFactorGradient(i) * dFactorGradient(j)
                end do
            end do
        end if

        dValue = dOldValue * dPowerValue
        dGradient = dOldGradient * dPowerValue + dOldValue * dPowerGradient
        do j = 1, nLocal
            do i = 1, nLocal
                dHessian(i,j) = dOldHessian(i,j) * dPowerValue + &
                    dOldGradient(i) * dPowerGradient(j) + dPowerGradient(i) * dOldGradient(j) + &
                    dOldValue * dPowerHessian(i,j)
            end do
        end do

        return

    end subroutine MultiplyLinearPower

    subroutine ProjectScalarHessian(iFirstIn, nLocal, dScalarHessianIn, dSpeciesHessianOut)

        implicit none

        integer, intent(in) :: iFirstIn, nLocal
        real(8), intent(in) :: dScalarHessianIn(nLocal,nLocal)
        real(8), intent(out) :: dSpeciesHessianOut(nSpeciesDim,nSpeciesDim)

        integer :: i, j, k
        real(8) :: dWeightedColumn

        dSpeciesHessianOut = 0D0
        do j = 1, nLocal
            dWeightedColumn = 0D0
            do k = 1, nLocal
                dWeightedColumn = dWeightedColumn + dMolFraction(iFirstIn+k-1) * dScalarHessianIn(k,j)
            end do
            do i = 1, nLocal
                dSpeciesHessianOut(i,j) = dScalarHessianIn(i,j) - dWeightedColumn
                if (i == j) dSpeciesHessianOut(i,j) = dSpeciesHessianOut(i,j) + &
                    1D0 / DMAX1(dMolFraction(iFirstIn+i-1), 1D-75)
            end do
        end do

        return

    end subroutine ProjectScalarHessian

end subroutine CompHessianQKTO
