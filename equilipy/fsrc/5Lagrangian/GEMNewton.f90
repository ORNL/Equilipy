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
    !                                        not represented by the current phase assemblage.
    !   03/04/2013      M.H.A. Piro         Fix bug in correction process when dealing with ionic phases
    !                                        the loop should count back from the number of constraints,
    !                                        not the number of charged phases).
    !
    !
    ! Purpose:
    ! ========
    !
    !> \brief The purpose of this subroutine is to compute the direction vector for the Gibbs energy
    !! minimization (GEM) solver using Newton's method.  The Hessian matrix and its corresponding constraint
    !! vector are first constructed and then the direction vector representing the system parameters is solved
    !! with the DGESV driver routine from LAPACK.  The updated element potentials, adjustments to the number of
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
    integer, dimension(:),   allocatable :: IPIV
    real(8)                              :: dTemp
    real(8), dimension(:),   allocatable :: B
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
    A                   = 0D0
    B                   = 0D0
    dUpdateVar          = 0D0
    dEffStoichSolnPhase = 0D0
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
                dTemp = sum(dStoichSpecies(m:n,i) * dStoichSpecies(m:n,j) * dMolesSpecies(m:n))
                A(i,j) = A(i,j)+ dTemp/ (DFLOAT(iParticlesPerMole(l))**2)
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
            dTemp = sum(dStoichSpecies(m:n,j) * dMolesSpecies(m:n) * (dChemicalPotential(m:n) - 1D0))
            dTemp = dTemp / DFLOAT(iParticlesPerMole(i))
            B(j)  = B(j) + dTemp
        end do
    end do
!
    ! Construct the Hessian matrix and constraint vector (contribution from solution phases):
    do j = nElements + 1, nElements + nSolnPhases
        l = 2 * nElements - j + 1       ! Relative solution phase index (in iAssemblage vector).
        k = -iAssemblage(l)             ! Absolute solution phase index.
        m = nSpeciesPhase(k-1) + 1
        n = nSpeciesPhase(k)
!
        ! Compute the stoichiometry of this phase:
        ! call CompStoichSolnPhase(k)

        do i = 1,nElements
            ! A(i,j) = dEffStoichSolnPhase(k,i) * dMolesPhase(l)
            A(i,j) = sum(dStoichSpecies(m:n,i) * dMolesSpecies(m:n))
            A(j,i) = A(i,j)
        end do
        ! B(j) = sum(dChemicalPotential(m:n) * dMolFraction(m:n))* dMolesPhase(l)
        B(j) = sum(dChemicalPotential(m:n) * dMolesSpecies(m:n))
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
        call dgesv( nVarLocal, 1, A, nVarLocal, IPIV, B, nVarLocal, INFOLocal )
    else
        do i = 1, nElements
            B(i) = dElementPotential(i)
        end do
        B(nElements + 1) = dMolesPhase(1)
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
        do j = 1, nVarLocal
            dUpdateVar(j) = B(j)
        end do
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
    deallocate(A, B, IPIV, STAT = i)
    if (i /= 0) INFOThermo = 24
!
    return
!
end subroutine GEMNewton
!
!
