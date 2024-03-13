!
subroutine GetNewAssemblage(iter)
!
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    GetNewAssemblage.f90
    !> \brief   Determine the next phase assemblage to be considered in Leveling.
    !> \author  M.H.A. Piro
    !> \date    Apr. 24, 2012
    !> \sa      ShuffleAssemblage.f90
    !> \sa      LevelingSolver.f90
    !> \sa      PostLevelingSolver.f90
    !
    !
    ! References:
    ! ===========
    !
    ! For further information regarding this method, refer to the following material:
    !
    !        M.H.A. Piro, "Computation of Thermodynamic Equilibria Pertinent to Nuclear Materials
    !        in Multi-Physics Codes," PhD Dissertation, Royal Military College of Canada, 2011.
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   03/31/2011      M.H.A. Piro         Original code
    !   07/31/2011      M.H.A. Piro         Clean up code: remove unnecessary variables, update variable names
    !   10/21/2011      M.H.A. Piro         Clean up code: modules, simplify iteration history check.
    !   01/21/2013      M.H.A. Piro         Improved the iteration history check and created the
    !                                        CheckLevelingIterHistory subroutine.  Also, a previous check for
    !                                        the mass balance constraints computed the sum of coefficients of
    !                                        A along the 2nd dimension.  This was changed from an integer
    !                                        vector to a double vector and the absolute quantity of each
    !                                        coefficient is taken to appropriately handle anions.
    !   02/14/2013      M.H.A. Piro         The ShuffleDummySpecies subroutine was created to shuffle dummy
    !                                        species to the front of the iAssemblage vector.
    !   03/24/2021      S.Y. Kwon           dLeveling is removed and dPhasePotential is introduced.
    !    
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to provide a new estimated phase assemblage to be tested in
    !! the Leveling subroutine.  The phase with the most negative relative Gibbs energy will be introduced
    !! into the previous estimated phase assemblage.  In order to avoid violating the Gibbs Phase Rule,
    !! this new phase will replace a phase from the previous assemblage.  The conditions for considering
    !! a new phase assemblage are:
    !!
    !! <ol>
    !! <li> The number of moles of each phase must be non-negative and real, </li>
    !! <li> The phase assemblage has not been previously tested (this is only
    !!      perfomed after x iterations), </li>
    !! <li> The numerical adjustments to the Gibbs Plane are real.  This last check does
    !!      not add any additional expense because this would have to be computed anyways. </li>
    !! </ol>
    !!
    !! There is an important, yet subtle, verification performed in the third (3) condition.  The Phase Rule
    !! dictates that the maximum number of phases that can coexist in an isobaric-isothermal closed system
    !! cannot exceed the number of system components (elements).  An additional condition to the Phase Rule,
    !! which is normally implied in thermodynamics texts but not explicity stated, is that only one pure
    !! separate phase of one component can exist at equilibrium.  For example, U(BCC) and U(FCC) cannot coexist.
    !! If two pure separate phases of the same component (same X) are included in the phase assemblage, then the
    !! Gibbs Plane is not uniquely defined (A is not a unique matrix) resulting in non-real adjustments to the
    !! Gibbs Plane.  By ensuring that the adjustments applied to the Gibbs Plane are real not only avoids obvious
    !! numerical problems, but also guarentees that the Phase Rule is explicitly satisfied.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! INFO                      An integer scalar used internally by LAPACK that returns 0 for a successful exit
    !                            and a non-zero value when an error has occurred.
    ! iNewPhase                 Integer index representing the new phase that is to be introduced to the system.
    ! iterHistoryLevel          Integer matrix storing the history of the phase assemblages that were tested.
    ! iAssemblage               Integer vector containing the indices of phases in the assemblage
    ! iSpeciesAtoms             Stoichiometry coefficient of a particular species (e.g., UO2: 1 U, 2 O).
    ! nElements                 An integer scalar representing the number of elements in the system.
    ! dMolesPhase               Number of moles of a phase.  These are directly mapped to phases in iAssemblage
    ! dMolesElement             Number of moles of a particular element in the system.
    ! dChemicalPotential        The Relative Gibbs Energy is defined as the difference between the chemical
    !                            potential and standard Gibbs energy of the pure species.
    ! dLevel                    The adjustment applied to the element potentials.
    ! dAtomFractionSpecies      Atomic fraction of a particular element (e.g., UO2: 0.333 U, 0.667 O).
    !
    !-------------------------------------------------------------------------------------------------------------
    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver
!
    implicit none
!
    integer                                :: i, j, k, m, n, INFO, iNewPhase, iter, iPhaseTypeOut, ipriori
    integer,dimension(nElements)           :: iTempAssemblage, IPIV,idxpriori
    real(8),dimension(nElements)           :: dTempChemicalPotentialGEM, dMinValuedMolesPhase, dpriori
    real(8),dimension(nElements,nElements) :: A, dTempStoichSpeciesGEM, dTempAtomFractionSpeciesGEM
    real(8),dimension(nElements,nSpecies)  :: dTempMolFractionGEM
    logical                                :: lPhasePass
!
!
    ! Initialize variables:
    IPIV   = 0
    A      = 0D0
    dpriori= 0

    iTempAssemblage             = 0
    dTempChemicalPotentialGEM   = 0D0
    dTempStoichSpeciesGEM       = 0D0
    dTempAtomFractionSpeciesGEM = 0D0
    dTempMolFractionGEM         = 0D0
    dMinValuedMolesPhase        = -1D3
    
    ! Determine index of the species with the most negative relative Gibbs energy:
    iNewPhase = MAXVAL(MINLOC(dPhasePotential))
    
    ! Shuffle the phase assemblage to make the best candidate phase to be tested first:
    call ShuffleAssemblage(iNewPhase,iPhaseTypeOut)
!
    ! Shuffle variables according to the result of Shuffle Assemblage
    dTempChemicalPotentialGEM   = dChemicalPotentialGEM
    dTempStoichSpeciesGEM       = dStoichSpeciesGEM
    dTempAtomFractionSpeciesGEM = dAtomFractionSpeciesGEM
    dTempMolFractionGEM         = dMolFractionGEM
   
    do i = 1, nElements
        j = iShuffled(i)
        dChemicalPotentialGEM(i)     = dTempChemicalPotentialGEM(j)
        dStoichSpeciesGEM(i,:)       = dTempStoichSpeciesGEM(j,:)
        dAtomFractionSpeciesGEM(i,:) = dTempAtomFractionSpeciesGEM(j,:)
        dMolFractionGEM(i,:)         = dTempMolFractionGEM(j,:)
    end do

    ! ! Store the indices of the estimated phase assemblage at each iteration
    iterHistoryLevel(:,iter) = iAssemblage
    dMolesPhaseHistory(:,iter)= dMolesPhase
!
    ! Store the sorted version of previous assemblage variables
    iTempAssemblage = iAssemblage
    dTempChemicalPotentialGEM   = dChemicalPotentialGEM
    dTempStoichSpeciesGEM       = dStoichSpeciesGEM
    dTempAtomFractionSpeciesGEM = dAtomFractionSpeciesGEM
    dTempMolFractionGEM         = dMolFractionGEM

!
    ! Loop through all "phases" in the current phase assemblage to determine which one should be
    ! substituted for the new "phase":
    LOOP_NewAssemblage: do k = 1, nElements

        !1. Adjust variables to the new assemblage
        iAssemblage(k)               = iNewPhase
        dChemicalPotentialGEM(k)     = dChemicalPotential(iNewPhase)
        dStoichSpeciesGEM(k,:)       = dStoichSpeciesLevel(iNewPhase,:)
        dAtomFractionSpeciesGEM(k,:) = dAtomFractionSpecies(iNewPhase,:)

    

        j = iPhaseLevel(iNewPhase)
        if (j>0) then
            m = nSpeciesPhase(j-1) + 1      ! First constituent in phase.
            n = nSpeciesPhase(j)            ! Last  constituent in phase.
            dMolFractionGEM(k,m:n)      = dMolFraction(m:n)
        end if

        ! 2. Calculate phase amount
        do i = 1, nElements
            do j = 1, nElements
                A(j,i) = dStoichSpeciesGEM(i,j)
            end do
            dMolesPhase(i) = dMolesElement(i)
        end do

        ! ! Check if this phase (vertex) yield the bulk composition to be in 
        ! call InSimplex(dStoichSpeciesGEM, nElements, dMolesElement, lVertexPass)
        ! if (.NOT.lVertexPass) cycle

        ! 2.1 Call the linear equation solver to compute molar quantities of the phase assemblage:
        INFO = 0
        IPIV = 0
        call DGESV( nElements, 1, A, nElements, IPIV, dMolesPhase, nElements, INFO )

        if(INFO==0) dMinValuedMolesPhase(k)=MINVAL(dMolesPhase)

        if ((INFO==0).AND.(dMinValuedMolesPhase(k) > 0D0)) exit LOOP_NewAssemblage

        if (k<nElements) then 
            ! If negative phase amount appears, revert the system and move onto the next
            ! Revert to previous assemblage
            iAssemblage             = iTempAssemblage
            dStoichSpeciesGEM       = dTempStoichSpeciesGEM
            dChemicalPotentialGEM   = dTempChemicalPotentialGEM
            dAtomFractionSpeciesGEM = dTempAtomFractionSpeciesGEM
            dMolFractionGEM         = dTempMolFractionGEM
        else
            ! if k=nElements and we didn't find the species to remove without resulting positive amount
            ! Revert to previous assemblage
            iAssemblage             = iTempAssemblage
            dStoichSpeciesGEM       = dTempStoichSpeciesGEM
            dChemicalPotentialGEM   = dTempChemicalPotentialGEM
            dAtomFractionSpeciesGEM = dTempAtomFractionSpeciesGEM
            dMolFractionGEM         = dTempMolFractionGEM

            !Check if the most possible assemblage has tolerable negative amount
            if((MINVAL(DABS(dMinValuedMolesPhase))<1D-3)) then
                !Check if there is zero in dMinValuedMolesPhase
                i = MAXLOC(dMinValuedMolesPhase,DIM=1)
                if(.NOT.(lPhasePass).AND.(count(DABS(dMinValuedMolesPhase)<1D-3)>1)) then
                    LOOP_i:do m =1,nElements
                        if (ABS(dMinValuedMolesPhase(m))<1D-3) then
                            i = m
                            exit LOOP_i
                        end if
                    end do LOOP_i
                end if
            endif

            iAssemblage(i)               = iNewPhase
            dChemicalPotentialGEM(i)     = dChemicalPotential(iNewPhase)
            dStoichSpeciesGEM(i,:)       = dStoichSpeciesLevel(iNewPhase,:)
            dAtomFractionSpeciesGEM(i,:) = dAtomFractionSpecies(iNewPhase,:)

            j = iPhaseLevel(iNewPhase)
            if (j>0) then
                m = nSpeciesPhase(j-1) + 1      ! First constituent in phase.
                n = nSpeciesPhase(j)            ! Last  constituent in phase.
                dMolFractionGEM(i,m:n)      = dMolFraction(m:n)
            end if

            ! 2. Calculate phase amount
            do i = 1, nElements
                do j = 1, nElements
                    A(j,i) = dStoichSpeciesGEM(i,j)
                end do
                dMolesPhase(i) = dMolesElement(i)
            end do

            ! 2.1 Call the linear equation solver to compute molar quantities of the phase assemblage:
            INFO = 0
            IPIV = 0
            call DGESV( nElements, 1, A, nElements, IPIV, dMolesPhase, nElements, INFO )
            if (INFO /= 0) INFOThermo = 10
        end if
                
    end do LOOP_NewAssemblage

    ! Re-Calculate elemental potentiall based on new assemblage
    do j = 1,nElements
        do i = 1,nElements
            A(i,j) = dAtomFractionSpeciesGEM(i,j)
        end do
        dElementPotential(j) = dChemicalPotentialGEM(j)
    end do

    ! Reinitialize variables:
    INFO = 0
    IPIV = 0
!
    ! Call linear equation solver to solve the adjustments applied to the Gibbs Plane:
    call DGESV( nElements, 1, A, nElements, IPIV, dElementPotential, nElements, INFO ) !SY: dLeveling ->dElementalPotential
!
    if (INFO /= 0) INFOThermo = 10

    ! ! Calculate phase potential of each species
    ! dPhasePotential = dChemicalPotential - MATMUL(dAtomFractionSpecies,dElementPotential)
    ! dMinPhasePotential = MINVAL(dPhasePotential)

    ! Check if this phase assemblage has been previously cosidered
    call CheckLevelingIterHistory(iter,lPhasePass)

    if(.NOT.(lPhasePass)) then
        ! Relax the criteria for phase potential
        dToleranceLevel = -1D-5
        ! print*,'Repeated'
    end if
    return
     
    
!
!
!
!
end subroutine GetNewAssemblage
!
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
!
!
!
!
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
!
!
!
!
    !---------------------------------------------------------------------------
    !                       END - GetNewAssemblage.f90
    !---------------------------------------------------------------------------
