!
!
!
subroutine GEMNewton(INFOLocal)
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    GEMNewton.f90
    !> \brief   Compute the direction vector for the GEMSolver using Newton's method.
    !> \author  M.H.A. Piro
    !> \date    Apr. 25, 2012
    !> \sa      GEMSolver.f90
    !> \sa      GEMLineSearch.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   04/25/2012      M.H.A. Piro         Original code (new GEM solver)
    !   05/25/2012      M.H.A. Piro         Check for a NAN immediately after call to DGESV.
    !   01/31/2013      M.H.A. Piro         Check if a charged phase is contained in the database, but is
    !                                       not represented by the current phase assemblage.
    !   03/04/2013      M.H.A. Piro         Fix bug in correction process when dealing with ionic phases
    !                                       the loop should count back from the number of constraints,
    !                                       not the number of charged phases).
    !   07/20/2026      S.Y. Kwon           Added analytical solution-phase Newton directions with projected constitution and descent safeguards.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \brief The purpose of this subroutine is to compute the direction vector for the Gibbs energy
    !! minimization (GEM) solver using Newton's method.  The Hessian matrix and its corresponding constraint
    !! vector are first constructed and then the direction vector representing the system parameters is solved
    !! with the DSYSV symmetric indefinite driver routine from LAPACK.  The updated element potentials, adjustments to the number of
    !! moles of solution phases and the number of moles of pure condensed phases are applied in the
    !! GEMLineSearch.f90 subroutine.
    !!
    !! Thermochimica is capable of handling ionic phases, which have an additional charge neutrality
    !! constraint imposed for each ionic phase.  Thus, an electron is added as a system component for every
    !! charged phase in the system.  It may be possible that an ionic phase is not predicted to be stable at
    !! a particular iteration and, thus, there aren't any stable species in the system representing that electron.
    !! To prevent a numerical singularity in the Hessian matrix, a check is performed after the Hessian matrix
    !! has been constructed ensuring that the Hessian does not contain a zero row.  In the event that the
    !! Hessian matrix contains all zeroes in the jth row (and necessarily, the jth column), a unit value is
    !! assigned to A(j,j).  Since the total balance of an electron is necessarily zero (i.e., ensuring charge
    !! neutrality) and there aren't any species for this solution phase, the corresponding value on the b vector
    !! will also be zero.  This procedure effectively ignores the jth row while preventing a numerical
    !! singularity.
    !
    !
    ! References:
    ! ===========
    !
    !> \details For further information regarding this methodology, refer to the following material:
    !! <ul>
    !! <li>  W.B. White, S.M. Johnson, G.B. Dantzig, "Chemical Equilibrium in Complex Mixtures," Journal of
    !!        Chemical Physics, V. 28, N. 5, 1958.
    !!
    !! <li>  G. Eriksson, "Thermodynamic Studies of High Temperature Equilibria," Acta Chemica Scandinavica,
    !!        25, 1971.
    !!
    !! <li>  G. Eriksson, E. Rosen, "General Equations for the Calculation of Equilibria in Multiphase Systems,"
    !!        Chemica Scripta, 4, 1973.
    !! </ul>
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[out]  INFOLocal        An integer scalar used by LAPACK indicating a successful exit or an error.
    !
    ! nVarLocal                      An integer scalar representing the total number of unknowns/linear equations.
    ! nElements                 An integer scalar representing the total number of elements in the system.
    ! nSpeciesPhase             An integer vector representing the number of species in a particular solution
    !                            phase (accumulative indexing)
    ! dStoichSpecies            A double real matrix representing stoichiometry coefficients.
    ! dMolesSpecies             A double real vector representing the number of moles of each species.
    ! dMolesPhase               A double real vector representing the number of moles of each phase.
    ! dMolesElement             A double real vector representing the number of moles of each element.
    ! JacobianLong              A double real matrix representing part of the Jacobian matrix that involves the
    !                            stoichiometry coefficients of solution species.
    ! JacobianShort             A double real vector that incorporates the JacobianLong matrix along with the
    !                            updated number of moles of each solution species.
    ! A                         Hessian matrix
    ! B                         Constraint vector (before call to LAPACK); unknown vector (after call to LAPACK)
    ! dEffStoichSolnPhase       A double real matrix representing the effective stoichiometry of a solution phase.
    ! dUpdateVar                A double real vector represending the updated system variables.
    !
    !-------------------------------------------------------------------------------------------------------------
    USE ModuleThermo
    USE ModuleThermoIO, ONLY: INFOThermo
    USE ModuleGEMSolver
!
    implicit none
!
    integer                              :: i, j, k, l, m,n, INFOLocal, nVarLocal
    integer                              :: INFO_SY, LWORK
    integer, dimension(:),   allocatable :: IPIV
    real(8)                              :: dTemp, dMatrixScale, dSymmetryTolerance, dSymmetryDiff
    real(8), dimension(1)                :: dWorkQuery
    real(8), dimension(:),   allocatable :: B, WORK
    real(8), dimension(:,:), allocatable :: A
!
!
    ! Determine the number of unknowns/linear equations:
    nVarLocal = nElements + nConPhases + nSolnPhases
!
    ! Allocate memory:
    allocate(A(nVarLocal, nVarLocal))
    allocate(B(nVarLocal))
    allocate(IPIV(nVarLocal))
!
    ! Initialize variables:
    IPIV                = 0
    INFOLocal           = 0
    INFO_SY             = 0
    A                   = 0D0
    B                   = 0D0
    dUpdateVar          = 0D0
    dEffStoichSolnPhase = 0D0
    iGEMNewtonSolver    = 0
    iGEMNewtonDSYSVInfo = 0
    iGEMNewtonKKTSize   = nVarLocal
    iGEMNewtonPivot1x1Count = 0
    iGEMNewtonPivot2x2Count = 0
    iGEMNewtonPivotPositiveCount = 0
    iGEMNewtonPivotNegativeCount = 0
    iGEMNewtonPivotZeroCount = 0
    iGEMAnalyticalHessianFallbackCount = 0
    dGEMNewtonSymmetryResidual = 0D0
    dGEMNewtonMinPivotScale = 0D0
    dGEMNewtonMaxPivotScale = 0D0
    dGEMNewtonDirectionNorm = 0D0
    if (allocated(dGEMAnalyticalSpeciesDirection)) dGEMAnalyticalSpeciesDirection = 0D0
    if (allocated(lGEMAnalyticalSpeciesDirection)) lGEMAnalyticalSpeciesDirection = .FALSE.
!
    ! Construct the Hessian matrix (elements):
    do j = 1, nElements
        do i = j, nElements
            do k = 1, nSolnPhases
                ! ! Absolute solution phase index:
                ! m = -iAssemblage(nElements - k + 1)
                ! ! Loop through species in phase:
                ! do l = nSpeciesPhase(m-1) + 1, nSpeciesPhase(m)
                !     dTemp  = dStoichSpecies(l,i) * dStoichSpecies(l,j) * dMolesSpecies(l)
                !     A(i,j) = A(i,j) + dTemp / (DFLOAT(iParticlesPerMole(l))**2)
                ! end do

                ! Absolute solution phase index:
                l = -iAssemblage(nElements - k + 1)
                m = nSpeciesPhase(l-1) + 1
                n = nSpeciesPhase(l)
                dTemp = sum(dStoichSpecies(m:n,i) * dStoichSpecies(m:n,j) * dMolesSpecies(m:n) / &
                    (DBLE(iParticlesPerMole(m:n))**2))
                A(i,j) = A(i,j)+ dTemp
            end do
            ! Apply symmetry:
            A(j,i) = A(i,j)
        end do
    end do
!
    ! Compute the constraint vector (elements):
    do j = 1, nElements
        B(j) = dMolesElement(j)
        ! do l = 1, nSolnPhases
        !     k = -iAssemblage(nElements - l + 1)
        !     do i = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
        !         dTemp = dStoichSpecies(i,j) * dMolesSpecies(i) * (dChemicalPotential(i) - 1D0)
        !         B(j)  = B(j) + dTemp / DFLOAT(iParticlesPerMole(i))
        !     end do
        ! end do

        do k = 1, nSolnPhases
            l = -iAssemblage(nElements - k + 1)
            m = nSpeciesPhase(l-1) + 1
            n = nSpeciesPhase(l)
            dTemp = sum(dStoichSpecies(m:n,j) * dMolesSpecies(m:n) * &
                (dChemicalPotential(m:n) - 1D0) / DBLE(iParticlesPerMole(m:n)))
            B(j)  = B(j) + dTemp
        end do
    end do
!
    ! Construct the Hessian matrix and constraint vector (contribution from solution phases):
    do j = nElements + 1, nElements + nSolnPhases
        l = 2 * nElements - j + 1       ! Relative solution phase index (in iAssemblage vector).
        k = -iAssemblage(l)             ! Absolute solution phase index.
!
        ! Compute the stoichiometry of this phase:
        call CompStoichSolnPhase(k)

        do i = 1,nElements
            A(i,j) = dEffStoichSolnPhase(k,i) * dMolesPhase(l)
            A(j,i) = A(i,j)
        end do
        B(j) = dGibbsSolnPhase(k)
    end do
    
!
    ! Construct the Hessian matrix and constraint vector (contribution from pure condensed phases):
    do j = nElements + nSolnPhases + 1, nVarLocal
        k = j - nElements - nSolnPhases
        do i = 1, nElements
            A(i,j) = dStoichSpecies(iAssemblage(k),i)
            A(j,i) = A(i,j)
        end do
        B(j) = dStdGibbsEnergy(iAssemblage(k))
    end do
!
    ! Check if the Hessian is properly structured if the system contains any charged phases:
    if (nCountSublattice > 0) then
        ! Loop through elements
        LOOP_SUB: do j = nElements, nElements - nChargedConstraints + 1, -1
            dTemp = 0D0
            ! Loop through coefficients along column:
            do i = 1, nElements
                dTemp = dTemp + DABS(A(i,j))
                if (dTemp > 0D0) cycle LOOP_SUB
            end do
            ! The phase corresponding to this electron is not stable.
            A(j,j) = 1D0
        end do LOOP_SUB
    end if
!
    ! Call the linear equation solver:
    if ((nConPhases > 1) .OR. (nSolnPhases > 0)) then
        dMatrixScale = 1D0
        dGEMNewtonSymmetryResidual = 0D0
        do j = 1, nVarLocal
            do i = 1, nVarLocal
                dMatrixScale = MAX(dMatrixScale, ABS(A(i,j)))
                if (i > j) then
                    dSymmetryDiff = ABS(A(i,j) - A(j,i))
                    if (dSymmetryDiff > dGEMNewtonSymmetryResidual) then
                        dGEMNewtonSymmetryResidual = dSymmetryDiff
                    end if
                end if
            end do
        end do
        dSymmetryTolerance = 100D0 * 2.220446049250313D-16 * dMatrixScale

        if (dGEMNewtonSymmetryResidual > dSymmetryTolerance) then
            INFO_SY = -999
        else
            ! The Newton matrix is symmetric indefinite. Use the symmetric
            ! Bunch-Kaufman solver directly; a singular structured solve is a
            ! Newton failure, not a reason to fall back to general LU.
            LWORK = -1
            dWorkQuery = 0D0
            call DSYSV('U', nVarLocal, 1, A, nVarLocal, IPIV, B, nVarLocal, &
                dWorkQuery, LWORK, INFO_SY)

            if (INFO_SY == 0) then
                LWORK = MAX(1,INT(dWorkQuery(1)))
                allocate(WORK(LWORK))
                call DSYSV('U', nVarLocal, 1, A, nVarLocal, IPIV, B, nVarLocal, &
                    WORK, LWORK, INFO_SY)
                if (allocated(WORK)) deallocate(WORK)
            end if
        end if

        iGEMNewtonDSYSVInfo = INFO_SY
        INFOLocal = INFO_SY
        call StoreNewtonPivotDiagnostics(nVarLocal, A, IPIV)
        if (INFO_SY == 0) iGEMNewtonSolver = 1
    else
        do i = 1, nElements
            B(i) = dElementPotential(i)
        end do
        B(nElements + 1) = dMolesPhase(1)
        iGEMNewtonSolver = 3
        dGEMNewtonDirectionNorm = DSQRT(SUM(B(1:nVarLocal)**2))
    end if

    ! ! Check solver status:
    ! do i = 1, nElements
    !     j = iAssemblage(i)
    !     if (j==-3) then
    !         m = nSpeciesPhase(-j-1) + 1
    !         n = nSpeciesPhase(-j)
    !         print*, dABS(dChemicalPotential(m:n)&
    !          - MATMUL(dStoichSpecies(m:n,:),B(:nElements))/iParticlesPerMole(m:n))
    !         print*, dStoichSpecies(m,:)
    !         print*, dStoichSpecies(m+1,:)
    !         print*, dStoichSpecies(m+2,:)
    !         print*, dStoichSpecies(m+3,:)
    !     end if
    ! end do
!
!
    ! Check for a NAN:
    LOOP_CheckNan: do i = 1, nVarLocal
        if (B(i) /= B(i)) then
            INFOLocal = 1
            exit LOOP_CheckNan
        end if
    end do LOOP_CheckNan
!
    ! Store the updated variables if LAPACK is successful:
    if (INFOLocal == 0) then
        dGEMNewtonDirectionNorm = DSQRT(SUM(B(1:nVarLocal)**2))
        do j = 1, nVarLocal
            dUpdateVar(j) = B(j)
        end do
        call BuildAnalyticalSpeciesDirections
!
        ! Reset:
        lRevertSystem = .FALSE.
    else
        ! The system failed.  Revert to a previous assemblage.
        lRevertSystem = .TRUE.
        dUpdateVar    = 0D0
    end if
!
    ! Deallocate memory of local variables:
    i = 0
    if (allocated(WORK)) deallocate(WORK)
    deallocate(A, B, IPIV, STAT = i)
    if (i /= 0) INFOThermo = 24
!
    return
!
contains

    subroutine StoreNewtonPivotDiagnostics(nLocal, AFactor, IPIVLocal)

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
                call CountNewtonPivotSign(dEig1, dTolLocal)
                call CountNewtonPivotSign(dEig2, dTolLocal)
                iLocal = iLocal + 2
            else
                iGEMNewtonPivot1x1Count = iGEMNewtonPivot1x1Count + 1
                dPivotLocal = DABS(AFactor(iLocal,iLocal))
                call CountNewtonPivotSign(AFactor(iLocal,iLocal), dTolLocal)
                iLocal = iLocal + 1
            end if

            dGEMNewtonMinPivotScale = DMIN1(dGEMNewtonMinPivotScale, dPivotLocal)
            dGEMNewtonMaxPivotScale = DMAX1(dGEMNewtonMaxPivotScale, dPivotLocal)
        end do

        if (dGEMNewtonMinPivotScale == HUGE(1D0)) dGEMNewtonMinPivotScale = 0D0

        return

    end subroutine StoreNewtonPivotDiagnostics

    subroutine CountNewtonPivotSign(dValue, dTolLocal)

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

    end subroutine CountNewtonPivotSign

    subroutine BuildAnalyticalSpeciesDirections

        implicit none

        integer :: iSolnSlot, iSolnPhase, iFirstLocal, iLastLocal, nLocal
        integer :: nVariableOut, iInfoHessian, nSiteCapacityLocal, nSiteOutScratch
        real(8), allocatable :: dSolutionChemicalPotential(:), dSolutionHessian(:,:)
        real(8), allocatable :: dSolutionHessianScratch(:,:), dSolutionJacobianScratch(:,:)
        real(8), allocatable :: dSpeciesDirectionLocal(:)
        character(len=8) :: cPhaseTypeLocal
        logical :: lUseSolutionHessian

        if (.NOT.allocated(dGEMAnalyticalSpeciesDirection)) return
        if (.NOT.allocated(lGEMAnalyticalSpeciesDirection)) return

        do iSolnSlot = 1, nSolnPhases
            iSolnPhase = -iAssemblage(nElements - iSolnSlot + 1)
            if ((iSolnPhase <= 0).OR.(iSolnPhase > nSolnPhasesSys)) cycle

            cPhaseTypeLocal = TRIM(cSolnPhaseType(iSolnPhase))
            lUseSolutionHessian = ((TRIM(cPhaseTypeLocal) == 'SUBG').OR.&
                (TRIM(cPhaseTypeLocal) == 'SUBQ').OR.&
                (TRIM(cPhaseTypeLocal) == 'RKMP').OR.&
                (TRIM(cPhaseTypeLocal) == 'QKTO').OR.&
                (TRIM(cPhaseTypeLocal) == 'SUBL').OR.&
                (TRIM(cPhaseTypeLocal) == 'SUBLM').OR.&
                (TRIM(cPhaseTypeLocal) == 'SUBOM')).AND.&
                (iPhaseElectronID(iSolnPhase) == 0)
            if (.NOT.lUseSolutionHessian) cycle

            iFirstLocal = nSpeciesPhase(iSolnPhase-1) + 1
            iLastLocal = nSpeciesPhase(iSolnPhase)
            nLocal = iLastLocal - iFirstLocal + 1
            if (nLocal <= 1) cycle

            if (allocated(lTraceSpeciesInactive)) then
                if (ANY(lTraceSpeciesInactive(iFirstLocal:iLastLocal))) cycle
            end if

            allocate(dSolutionChemicalPotential(nLocal), dSolutionHessian(nLocal,nLocal))
            allocate(dSpeciesDirectionLocal(nLocal))
            dSolutionChemicalPotential = 0D0
            dSolutionHessian = 0D0
            dSpeciesDirectionLocal = 0D0
            nVariableOut = 0
            iInfoHessian = 0
            nSiteOutScratch = 0

            select case (TRIM(cPhaseTypeLocal))
            case ('SUBG', 'SUBQ')
                call CompHessianSUBG(iSolnPhase, nLocal, dSolutionChemicalPotential, &
                    dSolutionHessian, nVariableOut, iInfoHessian)
            case ('RKMP')
                call CompHessianRKMP(iSolnPhase, nLocal, dSolutionChemicalPotential, &
                    dSolutionHessian, nVariableOut, iInfoHessian)
            case ('QKTO')
                call CompHessianQKTO(iSolnPhase, nLocal, dSolutionChemicalPotential, &
                    dSolutionHessian, nVariableOut, iInfoHessian)
            case ('SUBL', 'SUBLM', 'SUBOM')
                nSiteCapacityLocal = nMaxSublatticeSys * nMaxConstituentSys
                allocate(dSolutionHessianScratch(nSiteCapacityLocal,nSiteCapacityLocal), &
                    dSolutionJacobianScratch(nLocal,nLocal))
                call CompHessianSUBL(iSolnPhase, nSiteCapacityLocal, &
                    nLocal, dSolutionHessianScratch, dSolutionHessian, dSolutionJacobianScratch, &
                    nSiteOutScratch, nVariableOut, iInfoHessian)
            end select

            if ((iInfoHessian == 0).AND.(nVariableOut == nLocal)) then
                call SolveAnalyticalSpeciesDirection(iSolnSlot, iSolnPhase, nLocal, &
                    dSolutionHessian, dSpeciesDirectionLocal, iInfoHessian)
                if (iInfoHessian == 0) then
                    dGEMAnalyticalSpeciesDirection(iFirstLocal:iLastLocal) = dSpeciesDirectionLocal(1:nLocal)
                    lGEMAnalyticalSpeciesDirection(iSolnPhase) = .TRUE.
                end if
            end if

            if (allocated(dSolutionChemicalPotential)) deallocate(dSolutionChemicalPotential)
            if (allocated(dSolutionHessian)) deallocate(dSolutionHessian)
            if (allocated(dSolutionHessianScratch)) deallocate(dSolutionHessianScratch)
            if (allocated(dSolutionJacobianScratch)) deallocate(dSolutionJacobianScratch)
            if (allocated(dSpeciesDirectionLocal)) deallocate(dSpeciesDirectionLocal)
        end do

        call ProjectAnalyticalDirectionsToMassBalance

        return

    end subroutine BuildAnalyticalSpeciesDirections


    subroutine ProjectAnalyticalDirectionsToMassBalance

        implicit none

        integer :: i, j, l, k, m, n, iCol, INFO_LU
        integer, allocatable :: IPIVLocal(:)
        real(8) :: dTrialPhaseMoles, dTrialSpeciesMoles, dCompositionTerm
        real(8), allocatable :: ALocal(:,:), BLocal(:), dTrialFraction(:)
        logical :: lProjectDirections

        if (nConPhases + nSolnPhases /= nElements) return
        if (.NOT.allocated(lGEMAnalyticalSpeciesDirection)) return
        if (.NOT.allocated(dGEMAnalyticalSpeciesDirection)) return

        lProjectDirections = .TRUE.
        do l = 1, nSolnPhases
            k = -iAssemblage(nElements - l + 1)
            if ((k <= 0).OR.(k > SIZE(lGEMAnalyticalSpeciesDirection))) then
                lProjectDirections = .FALSE.
                exit
            end if
            if (.NOT.lGEMAnalyticalSpeciesDirection(k)) then
                lProjectDirections = .FALSE.
                exit
            end if
        end do
        if (.NOT.lProjectDirections) return

        allocate(ALocal(nElements,nElements), BLocal(nElements), IPIVLocal(nElements))
        allocate(dTrialFraction(nSpecies))
        ALocal = 0D0
        BLocal = dMolesElement(1:nElements)
        IPIVLocal = 0
        dTrialFraction = 0D0

        iCol = 0
        do i = 1, nConPhases
            iCol = iCol + 1
            k = iAssemblage(i)
            do j = 1, nElements
                ALocal(j,iCol) = dStoichSpecies(k,j)
            end do
        end do

        do l = 1, nSolnPhases
            iCol = iCol + 1
            k = -iAssemblage(nElements - l + 1)
            m = nSpeciesPhase(k-1) + 1
            n = nSpeciesPhase(k)
            dTrialPhaseMoles = 0D0
            do i = m, n
                dTrialSpeciesMoles = dMolesSpecies(i) + dGEMAnalyticalSpeciesDirection(i)
                if (dTrialSpeciesMoles <= dTolerance(8)) then
                    lProjectDirections = .FALSE.
                    exit
                end if
                dTrialPhaseMoles = dTrialPhaseMoles + dTrialSpeciesMoles
                dTrialFraction(i) = dTrialSpeciesMoles
            end do
            if (.NOT.lProjectDirections) exit
            if (dTrialPhaseMoles <= dTolerance(8)) then
                lProjectDirections = .FALSE.
                exit
            end if

            do i = m, n
                dTrialFraction(i) = dTrialFraction(i) / dTrialPhaseMoles
            end do

            do j = 1, nElements
                dCompositionTerm = 0D0
                do i = m, n
                    dCompositionTerm = dCompositionTerm + dTrialFraction(i) * &
                        dStoichSpecies(i,j) / DBLE(iParticlesPerMole(i))
                end do
                ALocal(j,iCol) = dCompositionTerm
            end do
        end do

        if (lProjectDirections) then
            call DGESV(nElements, 1, ALocal, nElements, IPIVLocal, BLocal, nElements, INFO_LU)
            if (INFO_LU /= 0) lProjectDirections = .FALSE.
        end if

        if (lProjectDirections) then
            iCol = 0
            do i = 1, nConPhases
                iCol = iCol + 1
                if (BLocal(iCol) <= dTolerance(8)) then
                    lProjectDirections = .FALSE.
                    exit
                end if
                dUpdateVar(nElements + nSolnPhases + i) = BLocal(iCol)
            end do
        end if

        if (lProjectDirections) then
            do l = 1, nSolnPhases
                iCol = nConPhases + l
                if (BLocal(iCol) <= dTolerance(8)) then
                    lProjectDirections = .FALSE.
                    exit
                end if
                k = -iAssemblage(nElements - l + 1)
                m = nSpeciesPhase(k-1) + 1
                n = nSpeciesPhase(k)
                do i = m, n
                    dGEMAnalyticalSpeciesDirection(i) = BLocal(iCol) * dTrialFraction(i) - dMolesSpecies(i)
                end do
            end do
        end if

        if (allocated(ALocal)) deallocate(ALocal)
        if (allocated(BLocal)) deallocate(BLocal)
        if (allocated(IPIVLocal)) deallocate(IPIVLocal)
        if (allocated(dTrialFraction)) deallocate(dTrialFraction)

        return

    end subroutine ProjectAnalyticalDirectionsToMassBalance


    subroutine SolveAnalyticalSpeciesDirection(iSolnSlotIn, iSolnPhaseIn, nLocalIn, dSolutionHessianIn, &
        dSpeciesDirectionOut, iInfoOut)

        implicit none

        integer, intent(in) :: iSolnSlotIn, iSolnPhaseIn, nLocalIn
        real(8), intent(in) :: dSolutionHessianIn(nLocalIn,nLocalIn)
        real(8), intent(out) :: dSpeciesDirectionOut(nLocalIn)
        integer, intent(out) :: iInfoOut

        integer :: iLocal, jLocal, iRow, jRow, iRef, nIndependent
        integer :: iSpecies, iRefSpecies, iElement, INFO_SY, LWORK, iSolveAttempt
        integer, allocatable :: iIndependent(:), IPIVLocal(:)
        real(8) :: dPlaneLocal, dPlaneRef, dPhaseDirection, dPhaseAmount
        real(8) :: dDeltaNIdeal, dProjectedHessian, dIdealReferenceCurvature
        real(8) :: dWorkQueryLocal(1)
        real(8), allocatable :: dProjectedMatrix(:,:), dProjectedMatrixModel(:,:)
        real(8), allocatable :: dProjectedMatrixIdealLog(:,:)
        real(8), allocatable :: dProjectedRHS(:), dProjectedRHSBase(:), dDeltaFraction(:)
        real(8), allocatable :: WORKLocal(:)
        logical :: lUseIdealLogCurvature

        iInfoOut = 0
        dSpeciesDirectionOut = 0D0
        if (nLocalIn <= 1) then
            iInfoOut = 1
            return
        end if

        iRef = MAXLOC(dMolFraction(nSpeciesPhase(iSolnPhaseIn-1)+1:nSpeciesPhase(iSolnPhaseIn)), DIM=1)
        nIndependent = nLocalIn - 1
        allocate(iIndependent(nIndependent), dProjectedMatrix(nIndependent,nIndependent), &
            dProjectedMatrixModel(nIndependent,nIndependent), &
            dProjectedMatrixIdealLog(nIndependent,nIndependent), &
            dProjectedRHS(nIndependent), dProjectedRHSBase(nIndependent), &
            dDeltaFraction(nLocalIn), IPIVLocal(nIndependent))

        iRow = 0
        do iLocal = 1, nLocalIn
            if (iLocal == iRef) cycle
            iRow = iRow + 1
            iIndependent(iRow) = iLocal
        end do

        dProjectedMatrix = 0D0
        dProjectedMatrixModel = 0D0
        dProjectedMatrixIdealLog = 0D0
        dProjectedRHS = 0D0
        dProjectedRHSBase = 0D0
        iRefSpecies = nSpeciesPhase(iSolnPhaseIn-1) + iRef
        dPlaneRef = 0D0
        do iElement = 1, nElements
            dPlaneRef = dPlaneRef + dUpdateVar(iElement) * dStoichSpecies(iRefSpecies,iElement) / &
                DBLE(iParticlesPerMole(iRefSpecies))
        end do

        do iRow = 1, nIndependent
            iLocal = iIndependent(iRow)
            iSpecies = nSpeciesPhase(iSolnPhaseIn-1) + iLocal
            dPlaneLocal = 0D0
            do iElement = 1, nElements
                dPlaneLocal = dPlaneLocal + dUpdateVar(iElement) * dStoichSpecies(iSpecies,iElement) / &
                    DBLE(iParticlesPerMole(iSpecies))
            end do
            dProjectedRHSBase(iRow) = (dPlaneLocal - dPlaneRef) - &
                (dChemicalPotential(iSpecies) - dChemicalPotential(iRefSpecies))

            do jRow = 1, nIndependent
                jLocal = iIndependent(jRow)
                dProjectedHessian = dSolutionHessianIn(iLocal,jLocal) - &
                    dSolutionHessianIn(iLocal,iRef) - &
                    dSolutionHessianIn(iRef,jLocal) + dSolutionHessianIn(iRef,iRef)
                dProjectedMatrixModel(iRow,jRow) = dProjectedHessian
            end do
        end do

        dIdealReferenceCurvature = 1D0 / DMAX1(dMolFraction(iRefSpecies), dTolerance(8))
        do iRow = 1, nIndependent
            iLocal = iIndependent(iRow)
            iSpecies = nSpeciesPhase(iSolnPhaseIn-1) + iLocal
            do jRow = 1, nIndependent
                dProjectedMatrixIdealLog(iRow,jRow) = dIdealReferenceCurvature
                if (iRow == jRow) then
                    dProjectedMatrixIdealLog(iRow,jRow) = dProjectedMatrixIdealLog(iRow,jRow) + &
                        1D0 / DMAX1(dMolFraction(iSpecies), dTolerance(8))
                end if
            end do
        end do

        dProjectedMatrix = dProjectedMatrixModel
        lUseIdealLogCurvature = .FALSE.
        ! Use the analytical solution-model Hessian as the first active-set
        ! direction.  The ideal-log matrix is the Gunnar-style fallback for a
        ! failed projected solve.

        INFO_SY = 0
        LOOP_HessianSolve: do iSolveAttempt = 1, 2
            dProjectedRHS = dProjectedRHSBase
            LWORK = -1
            dWorkQueryLocal = 0D0
            call DSYSV('U', nIndependent, 1, dProjectedMatrix, nIndependent, IPIVLocal, &
                dProjectedRHS, nIndependent, dWorkQueryLocal, LWORK, INFO_SY)

            if (INFO_SY == 0) then
                LWORK = MAX(1,INT(dWorkQueryLocal(1)))
                allocate(WORKLocal(LWORK))
                call DSYSV('U', nIndependent, 1, dProjectedMatrix, nIndependent, IPIVLocal, &
                    dProjectedRHS, nIndependent, WORKLocal, LWORK, INFO_SY)
                if (allocated(WORKLocal)) deallocate(WORKLocal)
            end if

            if (INFO_SY == 0) exit LOOP_HessianSolve
            if (lUseIdealLogCurvature) exit LOOP_HessianSolve

            lUseIdealLogCurvature = .TRUE.
            dProjectedMatrix = dProjectedMatrixIdealLog
        end do LOOP_HessianSolve

        if (INFO_SY /= 0) then
            iInfoOut = INFO_SY
            goto 900
        end if

        if (lUseIdealLogCurvature) iGEMAnalyticalHessianFallbackCount = iGEMAnalyticalHessianFallbackCount + 1

        dDeltaFraction = 0D0
        do iRow = 1, nIndependent
            dDeltaFraction(iIndependent(iRow)) = dProjectedRHS(iRow)
        end do
        dDeltaFraction(iRef) = -SUM(dDeltaFraction(1:nLocalIn))

        dPhaseDirection = 0D0
        do iLocal = 1, nLocalIn
            iSpecies = nSpeciesPhase(iSolnPhaseIn-1) + iLocal
            dPlaneLocal = 0D0
            do iElement = 1, nElements
                dPlaneLocal = dPlaneLocal + dUpdateVar(iElement) * dStoichSpecies(iSpecies,iElement) / &
                    DBLE(iParticlesPerMole(iSpecies))
            end do
            dDeltaNIdeal = dMolesSpecies(iSpecies) * &
                (dUpdateVar(nElements + iSolnSlotIn) + dPlaneLocal - dChemicalPotential(iSpecies))
            dPhaseDirection = dPhaseDirection + dDeltaNIdeal
        end do

        dPhaseAmount = SUM(dMolesSpecies(nSpeciesPhase(iSolnPhaseIn-1)+1:nSpeciesPhase(iSolnPhaseIn)))
        do iLocal = 1, nLocalIn
            iSpecies = nSpeciesPhase(iSolnPhaseIn-1) + iLocal
            dSpeciesDirectionOut(iLocal) = dMolFraction(iSpecies) * dPhaseDirection + &
                dPhaseAmount * dDeltaFraction(iLocal)
            if (dSpeciesDirectionOut(iLocal) /= dSpeciesDirectionOut(iLocal)) then
                iInfoOut = 2
                exit
            end if
        end do

    900 continue
        if (allocated(iIndependent)) deallocate(iIndependent)
        if (allocated(dProjectedMatrix)) deallocate(dProjectedMatrix)
        if (allocated(dProjectedMatrixModel)) deallocate(dProjectedMatrixModel)
        if (allocated(dProjectedMatrixIdealLog)) deallocate(dProjectedMatrixIdealLog)
        if (allocated(dProjectedRHS)) deallocate(dProjectedRHS)
        if (allocated(dProjectedRHSBase)) deallocate(dProjectedRHSBase)
        if (allocated(dDeltaFraction)) deallocate(dDeltaFraction)
        if (allocated(IPIVLocal)) deallocate(IPIVLocal)
        if (allocated(WORKLocal)) deallocate(WORKLocal)

        return

    end subroutine SolveAnalyticalSpeciesDirection

end subroutine GEMNewton
!
!
