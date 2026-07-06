!> Calculate updated minimum solution points during PEA and keep active slots synchronized.
subroutine CompMinSolnPoint
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CompMinSolnPoint.f90
    !> \brief   Calculate the minimum point of a solution phase for GEM purpose
    !> \author  S.Y. Kwon
    !> \date    May. 04, 2022
    !
    !
    ! Revisions:
    ! ==========
    !
!   Date            Programmer          Description of change
!   ----            ----------          ---------------------
!
!   05/04/2022      S.Y. Kwon           Original code
!   06/25/2026      S.Y. Kwon           Added corrected symmetry-broken SUBOM starts during PEA solution
!                                       minimization so ordered phases do not remain trapped at random states
!   06/25/2026      S.Y. Kwon           Limited SUBOM PEA starts to the nElements lowest endmember driving
!                                       forces after subtracting the elemental-potential plane
!   06/25/2026      S.Y. Kwon           Refresh active PEA solution-slot mole fractions after subminimization
!                                       updates a selected solution candidate row
!   06/26/2026      S.Y. Kwon           Rejected unknown submin starts from refreshed PEA pseudo-compound rows.
!   06/27/2026      S.Y. Kwon           Accepted general negative driving-force witnesses as valid
!                                       refreshed PEA pseudo-compound rows.
!   06/28/2026      S.Y. Kwon           Preferred converged starts over early-exit witnesses when ranking
!                                       refreshed PEA pseudo-compound candidates.
!   06/28/2026      S.Y. Kwon           Stopped overwriting active PEA slot fractions without refreshing
!                                       the matching pseudo-compound stoichiometry.
!   07/01/2026      S.Y. Kwon           Skipped mapped ordered phases without active ordering degrees of
!                                       freedom during refreshed PEA solution-candidate generation.
!   07/01/2026      S.Y. Kwon           Ranked converged and negative-witness SUBOM starts together by
!                                       driving force so ordered below-plane witnesses survive PEA.
!   07/01/2026      S.Y. Kwon           Included all refreshed SUBOM endpoint starts tied within the
!                                       Leveling-potential cutoff.
!   07/03/2026      S.Y. Kwon           Registered switch-gated second SUBOM composition-set rows in the
!                                       PEA candidate pool after refreshed solution minimization.
!
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to calculate the minimum point of a solution during GEM.
    !  dMoleFraction is carried over from the previous iterationg, assuming that it is close to
    !  the minimum point in the next iteration.
    !
    !-------------------------------------------------------------------------------------------------------------
    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver
!
    implicit none
!
    integer :: i,j,k,l,m,n,nConstituents,nStartConstituents
    real(8):: dNormComponent, dStartMinMoleFraction, dStartMaxMoleFraction
    real(8), parameter :: dEndpointTieTolerance = 1D-12
    real(8) :: dEndpointCutoffPotential
    logical :: lAddPhase, lSelectedCandidateValid, lSelectedCandidateDistinct
    logical :: lHasConvergedStart, lHasNegativeWitnessStart
    logical :: OrderDisorderPhaseIsEligible
    integer, dimension(:), allocatable :: iIndex, iEndmemberPotential, iCandidateStatusTemp
    real(8), dimension(:), allocatable :: dDrivingForceTemp, dEndmemberPotential
    real(8), dimension(:,:), allocatable :: dAtomFractionTemp, dStoichSpeciesTemp, dMolFractionTemp
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
    ! dChemicalPotential(nSpecies+1:)=5D9
    dPartialExcessGibbs=0D0
    dPartialEnthalpyXS=0D0
    dPartialEntropyXS=0D0
    dPartialHeatCapacityXS=0D0
    dDrivingForceSoln = 9D5
    iSubMinCandidateStatusSoln = SUBMIN_CANDIDATE_UNKNOWN
    dStartMinMoleFraction = 1D-5
!
!
    dMolFractionOld=dMolFraction
    dChemicalPotentialOld = dChemicalPotential(:nSpecies)
    k = 1
    LOOP_Soln: do i = 1, nSolnPhasesSys
        !initialize variables
        m                = nSpeciesPhase(i-1) + 1      ! First constituent in phase.
        n                = nSpeciesPhase(i)            ! Last  constituent in phase.
        nConstituents    = n - m + 1
        lAddPhase        = .False.
        if (.NOT.OrderDisorderPhaseIsEligible(i)) then
            iPhaseLevel(nSpecies+i) = i
            dChemicalPotential(nSpecies+i) = 5D9
            dPhasePotential(nSpecies+i) = 5D9
            call SetLevelingSolutionCandidateRow(nSpecies+i, i, dMolFraction(m:n), .FALSE.)
            cycle LOOP_Soln
        end if

        if (TRIM(cSolnPhaseType(i)) == 'SUBOM') then
            nStartConstituents = MIN(nConstituents, nElements)

            if (allocated(iIndex)) deallocate(iIndex)
            if (allocated(iEndmemberPotential)) deallocate(iEndmemberPotential)
            if (allocated(dDrivingForceTemp)) deallocate(dDrivingForceTemp)
            if (allocated(dEndmemberPotential)) deallocate(dEndmemberPotential)
            if (allocated(dAtomFractionTemp)) deallocate(dAtomFractionTemp)
            if (allocated(dStoichSpeciesTemp)) deallocate(dStoichSpeciesTemp)
            if (allocated(dMolFractionTemp)) deallocate(dMolFractionTemp)
            if (allocated(iCandidateStatusTemp)) deallocate(iCandidateStatusTemp)

            allocate(iEndmemberPotential(nConstituents))
            allocate(dEndmemberPotential(nConstituents))

            dEndmemberPotential = dLevelingChemicalPotential(m:n) - &
                MATMUL(dLevelingCompositionSpecies(m:n,:), dElementPotential)
            call Qsort(dEndmemberPotential, iEndmemberPotential, nConstituents)

            dEndpointCutoffPotential = dEndmemberPotential(nStartConstituents)
            do while (nStartConstituents < nConstituents)
                if (DABS(dEndmemberPotential(nStartConstituents+1) - &
                    dEndpointCutoffPotential) > dEndpointTieTolerance) exit
                nStartConstituents = nStartConstituents + 1
            end do

            allocate(iIndex(nStartConstituents), iCandidateStatusTemp(nStartConstituents))
            allocate(dDrivingForceTemp(nStartConstituents))
            allocate(dAtomFractionTemp(nStartConstituents,nElements))
            allocate(dStoichSpeciesTemp(nStartConstituents,nElements))
            allocate(dMolFractionTemp(nStartConstituents,nConstituents))

            iIndex = 1
            iCandidateStatusTemp = SUBMIN_CANDIDATE_UNKNOWN
            dDrivingForceTemp = 9D5
            dAtomFractionTemp = 0D0
            dStoichSpeciesTemp = 0D0
            dMolFractionTemp = 0D0

            do j = 1, nStartConstituents
                dMolFraction(m:n) = dStartMinMoleFraction
                dStartMaxMoleFraction = 1D0 - dStartMinMoleFraction * DFLOAT(nConstituents-1)
                dMolFraction(m+iEndmemberPotential(j)-1) = &
                    DMAX1(dStartMaxMoleFraction, 0.99D0)

                call Subminimization(i, lAddPhase)
                call CompStoichSolnPhase(i)

                dDrivingForceTemp(j) = dDrivingForceSoln(i)
                iCandidateStatusTemp(j) = iSubMinCandidateStatusSoln(i)
                dMolFractionTemp(j,:) = dMolFraction(m:n)
                dAtomFractionTemp(j,:) = dEffStoichSolnPhase(i,:) / SUM(dEffStoichSolnPhase(i,:))
                dStoichSpeciesTemp(j,:) = dEffStoichSolnPhase(i,:)
            end do

            ! Multiple starts should escape bad endpoints.  A negative witness is a
            ! thermodynamic phase-entry proof, so rank it with converged starts by
            ! driving force instead of letting a converged random state hide it.
            lHasConvergedStart = ANY(iCandidateStatusTemp == SUBMIN_CANDIDATE_CONVERGED)
            lHasNegativeWitnessStart = ANY(iCandidateStatusTemp == SUBMIN_CANDIDATE_NEGATIVE_WITNESS)
            if (lHasConvergedStart.OR.lHasNegativeWitnessStart) then
                do j = 1, nStartConstituents
                    if ((iCandidateStatusTemp(j) /= SUBMIN_CANDIDATE_CONVERGED).AND.&
                        (iCandidateStatusTemp(j) /= SUBMIN_CANDIDATE_NEGATIVE_WITNESS)) then
                        dDrivingForceTemp(j) = 9D5
                    end if
                end do
            end if

            call Qsort(dDrivingForceTemp, iIndex, nStartConstituents)

            dMolFraction(m:n) = dMolFractionTemp(iIndex(1),:)
            dGibbsSolnPhase(i) = 0D0
            call CompExcessGibbsEnergy(i)
            call CompStoichSolnPhase(i)
            dDrivingForceSoln(i) = dDrivingForceTemp(1)
            iSubMinCandidateStatusSoln(i) = iCandidateStatusTemp(iIndex(1))
        else
            call Subminimization(i, lAddPhase)
        end if

        call CompStoichSolnPhase(i)
        ! print*,'Composition',i,m,n,dMolFraction(m:n),dDrivingForceSoln(i)
        dAtomFractionSpecies(nSpecies+i,:) = dEffStoichSolnPhase(i,:)/sum(dEffStoichSolnPhase(i,:))
        dStoichSpeciesLevel(nSpecies+i,:)  = dEffStoichSolnPhase(i,:)
        iPhaseLevel(nSpecies+i)            = i
        lSelectedCandidateDistinct          = .TRUE.
        if (lMiscibility(i)) then
            if (sum(abs(dAtomFractionSpecies(nSpecies+i,:)-dAtomFractionSpecies(nSpecies+k,:)))/DFLOAT(nElements)>1D-4) then
                dChemicalPotential(nSpecies+i)     = dDrivingForceSoln(i) + &
                dot_product(dAtomFractionSpecies(nSpecies+i,:),dElementPotential(:))
            else
                dChemicalPotential(nSpecies+i) = 5D9
                lSelectedCandidateDistinct = .FALSE.
            end if
        else
            k=i
            dChemicalPotential(nSpecies+i)     = dDrivingForceSoln(i) + &
            dot_product(dAtomFractionSpecies(nSpecies+i,:),dElementPotential(:))  
        end if
        lSelectedCandidateValid = (iSubMinCandidateStatusSoln(i) == SUBMIN_CANDIDATE_CONVERGED).OR.&
            (iSubMinCandidateStatusSoln(i) == SUBMIN_CANDIDATE_NEGATIVE_WITNESS)
        lSelectedCandidateValid = lSelectedCandidateValid.AND.lSelectedCandidateDistinct
        call SetLevelingSolutionCandidateRow(nSpecies+i, i, dMolFraction(m:n), lSelectedCandidateValid)
        
        
    end do LOOP_Soln

    call RegisterSUBOMTwoSetCandidateRows

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
!
end subroutine CompMinSolnPoint
!
!
