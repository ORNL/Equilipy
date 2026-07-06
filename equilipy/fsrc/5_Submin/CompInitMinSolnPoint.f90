subroutine CompInitMinSolnPoint
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CompMinSolnPoint.f90
    !> \brief   Calculate the minimum point of a solution phase for leveling purpose
    !> \author  S.Y. Kwon
    !> \date    Nov. 01, 2021
    !
    !
    ! Revisions:
    ! ==========
    !
!   Date            Programmer          Description of change
!   ----            ----------          ---------------------
!
!   11/01/2021      S.Y. Kwon           Original code
!   06/25/2026      S.Y. Kwon           Used corrected Leveling potentials for endmember starts and
!                                       recomputed phase potentials after restoring stored local minima
!   06/25/2026      S.Y. Kwon           Limited SUBOM initial PEA starts to the nElements lowest endmember
!                                       driving forces after subtracting the elemental-potential plane
!   06/26/2026      S.Y. Kwon           Rejected unknown submin starts from valid initial PEA pseudo-compound rows.
!   06/27/2026      S.Y. Kwon           Accepted general negative driving-force witnesses as valid
!                                       initial PEA pseudo-compound rows.
!   06/28/2026      S.Y. Kwon           Included all Leveling-degenerate endpoint starts within 1D-12.
!   06/28/2026      S.Y. Kwon           Preferred converged starts over early-exit witnesses when ranking
!                                       initial PEA pseudo-compound candidates.
!   07/01/2026      S.Y. Kwon           Skipped mapped ordered phases without active ordering degrees of
!                                       freedom during initial PEA solution-candidate generation.
!   07/01/2026      S.Y. Kwon           Ranked converged and negative-witness SUBOM starts together by
!                                       driving force so ordered below-plane witnesses survive PEA.
!   07/03/2026      S.Y. Kwon           Registered switch-gated second SUBOM composition-set rows in the
!                                       PEA candidate pool after initial solution minimization.
!
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to calculate the global minimum point of a solution phase
    !  for the leveling purpose. All possible local minimums are calculated at once.
    !
    !-------------------------------------------------------------------------------------------------------------
    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver
!
    implicit none
!
    integer :: i, j, k, l, n, m, nConstituents, nElementOrConstituent
    integer :: nStartCount, nEndpointStartCount
    real(8), parameter :: dEndpointTieTolerance = 1D-12
    real(8) :: dMinMoleFraction, dMaxMoleFraction, dNormComponent
    real(8) :: dEndpointCutoffPotential
    logical :: lAddPhase, lSelectedCandidateValid, lSelectedCandidateDistinct
    logical :: lHasConvergedStart, lHasNegativeWitnessStart
    logical :: lUseSortedEndmemberStartsOnly
    logical :: OrderDisorderPhaseIsEligible
    integer, dimension(:), allocatable:: iIndex, iEndmemberPotential, iCandidateStatusTemp
    real(8), dimension(:), allocatable :: dDrivingForceTemp, dEndmemberPotential
    real(8), dimension(:,:), allocatable :: dAtomFractionTemp,dStoichSpeciesTemp, dMolFractionTemp
!
!
    l = MAX(1,nSolnPhasesSys)
!
    if (allocated(dMolFractionOld)) deallocate(dMolFractionOld)
    if (allocated(dEffStoichSolnPhase)) deallocate(dEffStoichSolnPhase)
    if (allocated(dSumMolFractionSoln)) deallocate(dSumMolFractionSoln)
    if (allocated(dDrivingForceSoln)) deallocate(dDrivingForceSoln)
    if (allocated(iSubMinCandidateStatusSoln)) deallocate(iSubMinCandidateStatusSoln)
    if (allocated(dMolesSpecies)) deallocate(dMolesSpecies)
    if (allocated(dPartialExcessGibbs)) then
        deallocate(dPartialExcessGibbs,dPartialEnthalpyXS,dPartialEntropyXS,dPartialHeatCapacityXS)
    end if

    allocate(&
        dMolesSpecies(nSpecies),&
        dMolFractionOld(nSpecies),&
        dEffStoichSolnPhase(l,nElements),&
        dSumMolFractionSoln(l),&
        dDrivingForceSoln(l),&
        iSubMinCandidateStatusSoln(l),&
        dPartialExcessGibbs(nSpecies),&
        dPartialEnthalpyXS(nSpecies),&
        dPartialEntropyXS(nSpecies),&
        dPartialHeatCapacityXS(nSpecies))
    dPartialExcessGibbs=0D0
    dPartialEnthalpyXS=0D0
    dPartialEntropyXS=0D0
    dPartialHeatCapacityXS=0D0
    dDrivingForceSoln = 9D5
    iSubMinCandidateStatusSoln = SUBMIN_CANDIDATE_UNKNOWN
    dMinMoleFraction = 1D-5
    
!
!
!   
  
    dMolFractionOld=dMolFraction
    dChemicalPotentialOld = dChemicalPotential(:nSpecies)
    LOOP_Soln: do i = 1, nSolnPhasesSys
        !initialize variables
        m                = nSpeciesPhase(i-1) + 1      ! First constituent in phase.
        n                = nSpeciesPhase(i)            ! Last  constituent in phase.
        nConstituents    = n - m + 1
        nElementOrConstituent = MIN(nElements,nConstituents)
        lUseSortedEndmemberStartsOnly = .FALSE.
        if (.NOT.OrderDisorderPhaseIsEligible(i)) then
            iPhaseLevel(nSpecies+i) = i
            dChemicalPotential(nSpecies+i) = 5D9
            dPhasePotential(nSpecies+i) = 5D9
            call SetLevelingSolutionCandidateRow(nSpecies+i, i, dMolFraction(m:n), .FALSE.)
            cycle LOOP_Soln
        end if
        if (TRIM(cSolnPhaseType(i)) == 'SUBOM') then
            nElementOrConstituent = MIN(nConstituents, nElements)
            lUseSortedEndmemberStartsOnly = .TRUE.
        end if
        dMaxMoleFraction = 1D0 - dMinMoleFraction * DFLOAT(nConstituents-1)
        dMaxMoleFraction = DMAX1(dMaxMoleFraction, 0.99D0)
        lAddPhase        = .FALSE.       
        
!
        if (.NOT.lMiscibility(i)) then
            ! Normal soln or first soln phase when immiscibility is considered
            ! Global minimum point of the solution phase is computed
            if(allocated(iIndex)) deallocate(iIndex)
            if(allocated(dDrivingForceTemp)) deallocate(dDrivingForceTemp)
            if(allocated(dAtomFractionTemp)) deallocate(dAtomFractionTemp)
            if(allocated(dStoichSpeciesTemp)) deallocate(dStoichSpeciesTemp)
            if(allocated(dMolFractionTemp)) deallocate(dMolFractionTemp)
            if(allocated(iCandidateStatusTemp)) deallocate(iCandidateStatusTemp)

            if(allocated(dEndmemberPotential)) deallocate(dEndmemberPotential)
            if(allocated(iEndmemberPotential)) deallocate(iEndmemberPotential)
            allocate(dEndmemberPotential(nConstituents),iEndmemberPotential(nConstituents))

            dEndmemberPotential = dLevelingChemicalPotential(m:n) - &
                MATMUL(dLevelingCompositionSpecies(m:n,:), dElementPotential)
            call Qsort(dEndmemberPotential,iEndmemberPotential,nConstituents)

            nEndpointStartCount = nElementOrConstituent
            dEndpointCutoffPotential = dEndmemberPotential(nElementOrConstituent)
            do while (nEndpointStartCount < nConstituents)
                if (DABS(dEndmemberPotential(nEndpointStartCount+1) - &
                    dEndpointCutoffPotential) > dEndpointTieTolerance) exit
                nEndpointStartCount = nEndpointStartCount + 1
            end do

            nStartCount = nEndpointStartCount
            if (.NOT.lUseSortedEndmemberStartsOnly) nStartCount = nEndpointStartCount + 1
            allocate(&
            iIndex(nStartCount),&
            iCandidateStatusTemp(nStartCount),&
            dDrivingForceTemp(nStartCount),&
            dAtomFractionTemp(nStartCount,nElements),&
            dStoichSpeciesTemp(nStartCount,nElements),&
            dMolFractionTemp(nStartCount, nConstituents))
!
            iIndex             = 1
            dDrivingForceTemp  = 9E5
            iCandidateStatusTemp = SUBMIN_CANDIDATE_UNKNOWN
            dAtomFractionTemp  = 0D0
            dStoichSpeciesTemp = 0D0
            dMolFractionTemp   = 0D0
!
            ! Perform subminimization multiple times by initializing from all extremums of the domain space:
            LOOP_Constituents1: do j = 1, nStartCount
            
                if (lUseSortedEndmemberStartsOnly) then
                    dMolFraction(m:n) = dMinMoleFraction
                    dMolFraction(m+iEndmemberPotential(j)-1) = dMaxMoleFraction
                else if(j>1) then
                    dMolFraction(m:n)   = dMinMoleFraction
                    dMolFraction(m+iEndmemberPotential(j-1)-1) = dMaxMoleFraction
                end if
                ! Perform subminimization:
                call Subminimization(i, lAddPhase)
!
                ! Store info of all local minima
                call CompStoichSolnPhase(i)
                dDrivingForceTemp(j)   = dDrivingForceSoln(i)
                iCandidateStatusTemp(j) = iSubMinCandidateStatusSoln(i)
                dMolFractionTemp(j,:)  = dMolFraction(m:n)
                dAtomFractionTemp(j,:) = dEffStoichSolnPhase(i,:)/sum(dEffStoichSolnPhase(i,:))
                dStoichSpeciesTemp(j,:)= dEffStoichSolnPhase(i,:)
            end do LOOP_Constituents1

            ! Multiple starts should escape bad endpoints.  A negative witness is a
            ! thermodynamic phase-entry proof, so rank it with converged starts by
            ! driving force instead of letting a converged random state hide it.
            lHasConvergedStart = ANY(iCandidateStatusTemp == SUBMIN_CANDIDATE_CONVERGED)
            lHasNegativeWitnessStart = ANY(iCandidateStatusTemp == SUBMIN_CANDIDATE_NEGATIVE_WITNESS)
            if (lHasConvergedStart.OR.lHasNegativeWitnessStart) then
                do j = 1, nStartCount
                    if ((iCandidateStatusTemp(j) /= SUBMIN_CANDIDATE_CONVERGED).AND.&
                        (iCandidateStatusTemp(j) /= SUBMIN_CANDIDATE_NEGATIVE_WITNESS)) then
                        dDrivingForceTemp(j) = 9D5
                    end if
                end do
            end if
!
            ! Sort the results according to the driving force (ascending order)
            ! note that output dDrivingForceTemp is the sorted version
            call Qsort(dDrivingForceTemp,iIndex,nStartCount)
            ! print*, cSolnPhaseName(i),iIndex, dMolFractionTemp(iIndex(1),:)
!
            ! Update variables for the current phase 
!
            dAtomFractionSpecies(nSpecies+i,:) = dAtomFractionTemp(iIndex(1),:)
            dStoichSpeciesLevel(nSpecies+i,:)  = dStoichSpeciesTemp(iIndex(1),:)
            dMolFraction(m:n) = dMolFractionTemp(iIndex(1),:)
            dGibbsSolnPhase(i) = 0D0
            call CompExcessGibbsEnergy(i)
            dChemicalPotential(nSpecies+i)     = dDrivingForceTemp(1) + &
            dot_product(dAtomFractionSpecies(nSpecies+i,:),dElementPotential(:))
            dDrivingForceSoln(i)               = dDrivingForceTemp(1)
            iSubMinCandidateStatusSoln(i)      = iCandidateStatusTemp(iIndex(1))
            iPhaseLevel(nSpecies+i)            = i
            lSelectedCandidateValid = (iSubMinCandidateStatusSoln(i) == SUBMIN_CANDIDATE_CONVERGED).OR.&
                (iSubMinCandidateStatusSoln(i) == SUBMIN_CANDIDATE_NEGATIVE_WITNESS)
            call SetLevelingSolutionCandidateRow(nSpecies+i, i, dMolFraction(m:n), lSelectedCandidateValid)
! !
            ! Update variables for immiscible phases
            LOOP_immiscibleSoln: do j = 2, nStartCount
!
                k = i + j - 1
                if (k > nSolnPhasesSys) exit LOOP_immiscibleSoln
                if (.NOT.lMiscibility(k)) cycle LOOP_immiscibleSoln
                if (cSolnPhaseName(i) /= cSolnPhaseName(k)) cycle LOOP_immiscibleSoln

                m  = nSpeciesPhase(k-1) + 1      ! First constituent in phase.
                n  = nSpeciesPhase(k)            ! Last  constituent in phase.
!
                dAtomFractionSpecies(nSpecies+k,:) = dAtomFractionTemp(iIndex(j),:)
                dStoichSpeciesLevel(nSpecies+k,:)  = dStoichSpeciesTemp(iIndex(j),:)
                dMolFraction(m:n) = dMolFractionTemp(iIndex(j),:)
                dGibbsSolnPhase(k) = 0D0
                call CompExcessGibbsEnergy(k)
                iPhaseLevel(nSpecies+k)            = k
                iSubMinCandidateStatusSoln(k)      = iCandidateStatusTemp(iIndex(j))
                lSelectedCandidateDistinct          = .TRUE.
                ! Store the local minima if it is not equal to global minimum
                if (sum(abs(dAtomFractionTemp(iIndex(1),:)-dAtomFractionTemp(iIndex(j),:)))/DFLOAT(nElements)>1D-4) then
                    dChemicalPotential(nSpecies+k) = dDrivingForceTemp(j) + &
                    dot_product(dAtomFractionSpecies(nSpecies+k,:),dElementPotential(:))
                else
                    dChemicalPotential(nSpecies+k) = 5D9
                    lSelectedCandidateDistinct = .FALSE.
                end if
                lSelectedCandidateValid = (iSubMinCandidateStatusSoln(k) == SUBMIN_CANDIDATE_CONVERGED).OR.&
                    (iSubMinCandidateStatusSoln(k) == SUBMIN_CANDIDATE_NEGATIVE_WITNESS)
                lSelectedCandidateValid = lSelectedCandidateValid.AND.lSelectedCandidateDistinct
                call SetLevelingSolutionCandidateRow(nSpecies+k, k, dMolFraction(m:n), lSelectedCandidateValid)
            end do LOOP_immiscibleSoln
        end if
!
!
!
    end do LOOP_Soln

    call RegisterSUBOMTwoSetCandidateRows
!
        ! Calculate functional norm: Mass balance 
    dGEMFunctionNormLast = dGEMFunctionNorm
    dGEMFunctionNorm=0D0
    do j = 1, nElements
        dNormComponent = dMolesElement(j)
        do i = 1,nElements
            dNormComponent = dNormComponent - dMolesPhase(i) * dStoichSpeciesGEM(i,j)
        end do
        dGEMFunctionNorm = dGEMFunctionNorm + (dNormComponent)**(2)
    end do

    ! Calculate chemical potential balance
    do i = 1, nElements
        k     = iAssemblage(i)
        dNormComponent = 0D0
        do j = 1, nElements
            dNormComponent = dNormComponent + dElementPotential(j) * dAtomFractionSpeciesGEM(i,j)
        end do
        dNormComponent            = dNormComponent - dChemicalPotentialGEM(i)
        dGEMFunctionNorm = dGEMFunctionNorm + (dNormComponent)**(2)
    end do
    dGEMFunctionNorm = dGEMFunctionNorm**(0.5)
    
!
    !After Subminimization, the chemical potential of solution components changes.
    !These need to be reverted back to pure
    
    dChemicalPotential(:nSpecies)     = (dStdGibbsEnergy) &
    *DFLOAT(iParticlesPerMole)/ dSpeciesTotalAtoms
!
    ! deallocate(dEffStoichSolnPhase,dSumMolFractionSoln,&
    ! dDrivingForceSoln,dPartialExcessGibbs)
!
end subroutine CompInitMinSolnPoint
!
!
!
