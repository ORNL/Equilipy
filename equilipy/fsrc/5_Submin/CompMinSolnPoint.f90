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
    !   07/20/2026      S.Y. Kwon           Refreshed typed solution candidates from model-appropriate starts and used one best grid point per phase in grid-mode certification sweeps.
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
    USE GridDiscovery, ONLY: GridTangentRowIndex
!
    implicit none
!
    integer :: i,j,k,l,m,n,iLoop,iLevelRow,iBaseLevelRow,nConstituents,nElementOrConstituent
    integer :: iPoint, iBestPoint, iGridRow, iGridOffset
    integer :: nStartCount, nEndpointStartCount
    real(8):: dNormComponent, dStartMinMoleFraction, dStartMaxMoleFraction
    real(8):: dSubMinElapsed
    real(8), parameter :: dEndpointTieTolerance = 1D-12
    real(8) :: dEndpointCutoffPotential, dGridPotential, dBestGridPotential
    logical :: lAddPhase, lSelectedCandidateValid, lSelectedCandidateDistinct
    logical :: lHasConvergedStart, lHasNegativeWitnessStart
    logical :: lUseSortedEndmemberStartsOnly
    logical :: OrderDisorderPhaseIsEligible
    integer, dimension(:), allocatable :: iIndex, iEndmemberPotential, iCandidateStatusTemp
    real(8), dimension(:), allocatable :: dDrivingForceTemp, dEndmemberPotential, dSubMinTimeTemp
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
    LOOP_Soln: do iLoop = 1, nSolnPhasesSys
        i = iLoop
        if (lPEADFSweepReversePhaseOrder) i = nSolnPhasesSys + 1 - iLoop
        !initialize variables
        m                = nSpeciesPhase(i-1) + 1      ! First constituent in phase.
        n                = nSpeciesPhase(i)            ! Last  constituent in phase.
        iLevelRow        = GridTangentRowIndex(i)
        nConstituents    = n - m + 1
        lUseSortedEndmemberStartsOnly = .FALSE.
        lAddPhase        = .False.
        if (.NOT.OrderDisorderPhaseIsEligible(i)) then
            iPhaseLevel(iLevelRow) = i
            dChemicalPotential(iLevelRow) = 5D9
            dPhasePotential(iLevelRow) = 5D9
            call SetLevelingSolutionCandidateRow(iLevelRow, i, dMolFraction(m:n), .FALSE.)
            cycle LOOP_Soln
        end if

        if (lGridFrontEndActive.AND.allocated(iGridPointPhase).AND.&
            allocated(iGridPointLevelRow).AND.allocated(iGridPointFractionOffset).AND.&
            allocated(dGridPointFraction)) then
            iBestPoint = 0
            dBestGridPotential = HUGE(1D0)
            do iPoint = 1, nGridPoint
                if (iGridPointPhase(iPoint) /= i) cycle
                iGridRow = iGridPointLevelRow(iPoint)
                if ((iGridRow < 1).OR.(iGridRow > nSpeciesLevel)) cycle
                dGridPotential = dLevelingChemicalPotential(iGridRow) - &
                    dot_product(dLevelingCompositionSpecies(iGridRow,:),dElementPotential)
                if (dGridPotential < dBestGridPotential) then
                    dBestGridPotential = dGridPotential
                    iBestPoint = iPoint
                end if
            end do
            if (iBestPoint > 0) then
                iGridOffset = iGridPointFractionOffset(iBestPoint)
                dMolFraction(m:n) = dGridPointFraction(iGridOffset:iGridOffset+nConstituents-1)
            end if
        else if (TRIM(cSolnPhaseType(i)) == 'SUBOM') then
            lUseSortedEndmemberStartsOnly = .TRUE.
        end if

        if (allocated(iIndex)) deallocate(iIndex)
        if (allocated(iEndmemberPotential)) deallocate(iEndmemberPotential)
        if (allocated(dDrivingForceTemp)) deallocate(dDrivingForceTemp)
        if (allocated(dEndmemberPotential)) deallocate(dEndmemberPotential)
        if (allocated(dSubMinTimeTemp)) deallocate(dSubMinTimeTemp)
        if (allocated(dAtomFractionTemp)) deallocate(dAtomFractionTemp)
        if (allocated(dStoichSpeciesTemp)) deallocate(dStoichSpeciesTemp)
        if (allocated(dMolFractionTemp)) deallocate(dMolFractionTemp)
        if (allocated(iCandidateStatusTemp)) deallocate(iCandidateStatusTemp)

        if (lUseSortedEndmemberStartsOnly) then
            allocate(iEndmemberPotential(nConstituents))
            allocate(dEndmemberPotential(nConstituents))

            dEndmemberPotential = dLevelingChemicalPotential(m:n) - &
                MATMUL(dLevelingCompositionSpecies(m:n,:), dElementPotential)
            call Qsort(dEndmemberPotential, iEndmemberPotential, nConstituents)

            nElementOrConstituent = MIN(nElements,nConstituents)
            nEndpointStartCount = nElementOrConstituent
            dEndpointCutoffPotential = dEndmemberPotential(nElementOrConstituent)
            do while (nEndpointStartCount < nConstituents)
                if (DABS(dEndmemberPotential(nEndpointStartCount+1) - &
                    dEndpointCutoffPotential) > dEndpointTieTolerance) exit
                nEndpointStartCount = nEndpointStartCount + 1
            end do
            nStartCount = nEndpointStartCount
        else
            nStartCount = 1
        end if

        allocate(iIndex(nStartCount), iCandidateStatusTemp(nStartCount))
        allocate(dDrivingForceTemp(nStartCount), dSubMinTimeTemp(nStartCount))
        allocate(dAtomFractionTemp(nStartCount,nElements))
        allocate(dStoichSpeciesTemp(nStartCount,nElements))
        allocate(dMolFractionTemp(nStartCount,nConstituents))

        iIndex = 1
        iCandidateStatusTemp = SUBMIN_CANDIDATE_UNKNOWN
        dDrivingForceTemp = 9D5
        dSubMinTimeTemp = 0D0
        dAtomFractionTemp = 0D0
        dStoichSpeciesTemp = 0D0
        dMolFractionTemp = 0D0
        if (lUseSortedEndmemberStartsOnly) then
            dStartMaxMoleFraction = 1D0 - dStartMinMoleFraction * DFLOAT(nConstituents-1)
            dStartMaxMoleFraction = DMAX1(dStartMaxMoleFraction, 0.99D0)
        end if

        do j = 1, nStartCount
            if (lUseSortedEndmemberStartsOnly) then
                dMolFraction(m:n) = dStartMinMoleFraction
                dMolFraction(m+iEndmemberPotential(j)-1) = dStartMaxMoleFraction
            else
                if (SUM(dMolFraction(m:n)) > 0D0) then
                    dMolFraction(m:n) = dMolFraction(m:n) / SUM(dMolFraction(m:n))
                else
                    dMolFraction(m:n) = 1D0 / DFLOAT(nConstituents)
                end if
            end if

            call TimedCandidateSubminimization(i, lAddPhase, dSubMinElapsed)
            call CompStoichSolnPhase(i)

            dDrivingForceTemp(j) = dDrivingForceSoln(i)
            dSubMinTimeTemp(j) = dSubMinElapsed
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
            do j = 1, nStartCount
                if ((iCandidateStatusTemp(j) /= SUBMIN_CANDIDATE_CONVERGED).AND.&
                    (iCandidateStatusTemp(j) /= SUBMIN_CANDIDATE_NEGATIVE_WITNESS)) then
                    dDrivingForceTemp(j) = 9D5
                end if
            end do
        end if

        call Qsort(dDrivingForceTemp, iIndex, nStartCount)
        do j = 1, nStartCount
            call RecordCandidateCertificate(i, j, iLevelRow, &
                iCandidateStatusTemp(iIndex(j)), &
                dDrivingForceTemp(j), dSubMinTimeTemp(iIndex(j)), m, n, j, &
                dMolFractionTemp(iIndex(j),:))
        end do

        dMolFraction(m:n) = dMolFractionTemp(iIndex(1),:)
        dGibbsSolnPhase(i) = 0D0
        call CompExcessGibbsEnergy(i)
        call CompStoichSolnPhase(i)
        dDrivingForceSoln(i) = dDrivingForceTemp(1)
        iSubMinCandidateStatusSoln(i) = iCandidateStatusTemp(iIndex(1))

        call CompStoichSolnPhase(i)
        dAtomFractionSpecies(iLevelRow,:) = &
            dEffStoichSolnPhase(i,:)/sum(dEffStoichSolnPhase(i,:))
        dStoichSpeciesLevel(iLevelRow,:) = dEffStoichSolnPhase(i,:)
        iPhaseLevel(iLevelRow) = i
        lSelectedCandidateDistinct          = .TRUE.
        if (lMiscibility(i)) then
            iBaseLevelRow = GridTangentRowIndex(k)
            if (sum(abs(dAtomFractionSpecies(iLevelRow,:)-&
                dAtomFractionSpecies(iBaseLevelRow,:)))/DFLOAT(nElements)>1D-4) then
                dChemicalPotential(iLevelRow) = dDrivingForceSoln(i) + &
                    dot_product(dAtomFractionSpecies(iLevelRow,:),dElementPotential(:))
            else
                dChemicalPotential(iLevelRow) = 5D9
                lSelectedCandidateDistinct = .FALSE.
            end if
        else
            k=i
            dChemicalPotential(iLevelRow) = dDrivingForceSoln(i) + &
                dot_product(dAtomFractionSpecies(iLevelRow,:),dElementPotential(:))
        end if
        lSelectedCandidateValid = (iSubMinCandidateStatusSoln(i) == SUBMIN_CANDIDATE_CONVERGED).OR.&
            (iSubMinCandidateStatusSoln(i) == SUBMIN_CANDIDATE_NEGATIVE_WITNESS)
        lSelectedCandidateValid = lSelectedCandidateValid.AND.lSelectedCandidateDistinct
        call SetLevelingSolutionCandidateRow(iLevelRow, i, dMolFraction(m:n), lSelectedCandidateValid)
        
        
    end do LOOP_Soln

    call ReconcileOrderDisorderCandidateRows
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
contains

    subroutine RecordCandidateCertificate(iSolnPhaseIn, iCopyIn, iLevelRowIn, &
        iSubMinStatusIn, dDrivingForceIn, dTimingIn, iFirstSpeciesIn, iLastSpeciesIn, iRankIn, &
        dFractionIn)

        implicit none

        integer, intent(in) :: iSolnPhaseIn, iCopyIn, iLevelRowIn
        integer, intent(in) :: iSubMinStatusIn, iFirstSpeciesIn, iLastSpeciesIn, iRankIn
        real(8), intent(in) :: dDrivingForceIn, dTimingIn
        real(8), dimension(:), intent(in) :: dFractionIn

        integer :: iCopy, iRecord, nCopyCapacity, nPlane

        nGEMCertEmissionCount = nGEMCertEmissionCount + 1
        call StagePEADFCandidate(iSolnPhaseIn,iLevelRowIn,CertificateBasis(iSolnPhaseIn),&
            CertificateStatus(iSubMinStatusIn,dDrivingForceIn),CertificateProof(iSubMinStatusIn),&
            iRankIn,dDrivingForceIn,dTimingIn,iFirstSpeciesIn,iLastSpeciesIn,dFractionIn)
        if (.NOT.allocated(iGEMCertPhase)) return
        if (nGEMCertCapacity <= 0) return

        nCopyCapacity = MAX(1,nElements)
        if ((iCopyIn < 1).OR.(iCopyIn > nCopyCapacity)) then
            nGEMCertDropped = nGEMCertDropped + 1
            return
        end if
        iCopy = iCopyIn
        iRecord = (iSolnPhaseIn - 1) * nCopyCapacity + iCopy
        if ((iRecord < 1).OR.(iRecord > nGEMCertCapacity)) then
            nGEMCertDropped = nGEMCertDropped + 1
            return
        end if

        if (iGEMCertPhase(iRecord) == 0) then
            nGEMCertCount = nGEMCertCount + 1
        end if

        iGEMCertPhase(iRecord) = iSolnPhaseIn
        iGEMCertCopy(iRecord) = iCopy
        iGEMCertLevelRow(iRecord) = iLevelRowIn
        iGEMCertBasis(iRecord) = CertificateBasis(iSolnPhaseIn)
        iGEMCertNorm(iRecord) = GEM_CERT_NORMALIZATION_PER_MOLE_ATOMS
        iGEMCertStatus(iRecord) = CertificateStatus(iSubMinStatusIn, dDrivingForceIn)
        iGEMCertProof(iRecord) = CertificateProof(iSubMinStatusIn)
        iGEMCertSubMinStatus(iRecord) = iSubMinStatusIn
        iGEMCertRank(iRecord) = iRankIn
        iGEMCertFirstSpecies(iRecord) = iFirstSpeciesIn
        iGEMCertLastSpecies(iRecord) = iLastSpeciesIn
        dGEMCertDrivingForce(iRecord) = dDrivingForceIn
        ! Timing records only the Subminimization call; emission overhead is not included.
        dGEMCertTiming(iRecord) = dTimingIn
        if (allocated(dGEMCertPlane)) then
            dGEMCertPlane(iRecord,:) = 0D0
            nPlane = MIN(nElements,SIZE(dGEMCertPlane,2))
            if (nPlane > 0) dGEMCertPlane(iRecord,:nPlane) = dElementPotential(:nPlane)
        end if

        return

    end subroutine RecordCandidateCertificate

    integer function CertificateBasis(iSolnPhaseIn)

        implicit none

        integer, intent(in) :: iSolnPhaseIn

        select case (TRIM(cSolnPhaseType(iSolnPhaseIn)))
        case ('SUBL','SUBLM','SUBOM')
            CertificateBasis = GEM_CERT_BASIS_SITE_FRACTION
        case ('SUBG')
            CertificateBasis = GEM_CERT_BASIS_PAIR_FRACTION
        case ('SUBQ')
            CertificateBasis = GEM_CERT_BASIS_QUADRUPLET_FRACTION
        case default
            CertificateBasis = GEM_CERT_BASIS_ENDMEMBER_FRACTION
        end select

        return

    end function CertificateBasis

    integer function CertificateStatus(iSubMinStatusIn, dDrivingForceIn)

        implicit none

        integer, intent(in) :: iSubMinStatusIn
        real(8), intent(in) :: dDrivingForceIn

        select case (iSubMinStatusIn)
        case (SUBMIN_CANDIDATE_CONVERGED, SUBMIN_CANDIDATE_NEGATIVE_WITNESS)
            if (dDrivingForceIn < 0D0) then
                CertificateStatus = GEM_CERT_STATUS_FAVORABLE
            else
                CertificateStatus = GEM_CERT_STATUS_FRESH_NEGATIVE
            end if
        case (SUBMIN_CANDIDATE_MAX_ITER, SUBMIN_CANDIDATE_REJECTED)
            CertificateStatus = GEM_CERT_STATUS_FAILED
        case (SUBMIN_CANDIDATE_DUPLICATE)
            CertificateStatus = GEM_CERT_STATUS_DUPLICATE
        case default
            CertificateStatus = GEM_CERT_STATUS_UNKNOWN
        end select

        return

    end function CertificateStatus

    integer function CertificateProof(iSubMinStatusIn)

        implicit none

        integer, intent(in) :: iSubMinStatusIn

        select case (iSubMinStatusIn)
        case (SUBMIN_CANDIDATE_CONVERGED)
            CertificateProof = GEM_CERT_PROOF_SUBMIN_CONVERGED
        case (SUBMIN_CANDIDATE_NEGATIVE_WITNESS)
            CertificateProof = GEM_CERT_PROOF_NEGATIVE_WITNESS
        case (SUBMIN_CANDIDATE_MAX_ITER)
            CertificateProof = GEM_CERT_PROOF_MAX_ITER
        case (SUBMIN_CANDIDATE_REJECTED)
            CertificateProof = GEM_CERT_PROOF_REJECTED
        case (SUBMIN_CANDIDATE_DUPLICATE)
            CertificateProof = GEM_CERT_PROOF_DUPLICATE
        case default
            CertificateProof = GEM_CERT_PROOF_UNKNOWN
        end select

        return

    end function CertificateProof

    subroutine TimedCandidateSubminimization(iSolnPhaseIn, lAddPhaseOut, dElapsedOut)

        implicit none

        integer, intent(in) :: iSolnPhaseIn
        logical, intent(out) :: lAddPhaseOut
        real(8), intent(out) :: dElapsedOut

        real(8) :: dStartTime, dStopTime

        call cpu_time(dStartTime)
        call Subminimization(iSolnPhaseIn, lAddPhaseOut)
        call cpu_time(dStopTime)
        dElapsedOut = DMAX1(0D0,dStopTime - dStartTime)

        return

    end subroutine TimedCandidateSubminimization

end subroutine CompMinSolnPoint
!
!
