!
subroutine SubMinInit(iSolnPhaseIndex,iterSubMax)
!
    !---------------------------------------------------------------------------
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this subroutine is to initialize the subminimiation routine.
    ! Variables are initialized, allocatable arrays are allocated and the
    ! chemical potential terms of the solution phase constituents defined by
    ! the element potentials are computed (i.e., dChemicalPotentialStar).
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! iFirstSUB            Absolute index of first species in phase.
    ! iLastSUB             Absolute index of last species in phase.
    ! nVar              Total number of species in phase.
    ! iterSubMax        An integer scalar representing the maximum number
    !                    of iterations allowed in subminimization.
    ! lSubMinConverged  A logical scalar indicating convergence (true).
    ! dDrivingForceLast A double real scalar representing the driving force
    !                    from the last subminimization iteration.
    !
    !---------------------------------------------------------------------------
!
    USE ModuleThermo
    USE ModuleSubMin
    USE ModuleGEMSolver!, ONLY: lMiscibility, dDrivingForceSoln
!
    implicit none
!
    integer::  i, j, k, iSolnPhaseIndex, iterSubMax
!
!
    ! Initialize variables:
    dDrivingForceLast = 10D0
    ! dMaxPotentialVectorLast = 1D3
    lSubMinConverged  = .FALSE.
!
    dDrivingForceSoln(iSolnPhaseIndex) = 0D0
!
!
    ! Check if allocatable arrays are already allocated.  Deallocate if necessary:
    if (allocated(dChemicalPotentialStar)) deallocate(dChemicalPotentialStar)
    if (allocated(dRHS))                   deallocate(dRHS)
    if (allocated(dRHSLast))               deallocate(dRHSLast)
    if (allocated(dHessian))               deallocate(dHessian)
    if (allocated(iHessian))               deallocate(iHessian)
    if (allocated(dPi)) deallocate(dPi)
    if (allocated(dPiLast)) deallocate(dPiLast)
    if (allocated(dPiDot)) deallocate(dPiDot)
    if (allocated(dLambda)) deallocate(dLambda)
    if (allocated(dLambdaLast)) deallocate(dLambdaLast)
    if (allocated(dvAdam)) deallocate(dvAdam)
    if (allocated(dsAdam)) deallocate(dsAdam)
    if (allocated(dvAdamLast)) deallocate(dvAdamLast)
    if (allocated(dsAdamLast)) deallocate(dsAdamLast)
    if (allocated(dxDot)) deallocate(dxDot)
    if (allocated(dPotentialVector)) deallocate(dPotentialVector)

    ! Determine prefactor for allocating arrays (depends on whether the phase is ionic):
    if (iPhaseElectronID(iSolnPhaseIndex) == 0) then
        i = 1
    else
        i = 2
    end if
    
    !We will revise this part with nVarSub to remove zeroing species
    ! All variable related to Submin will be revised including dChemicalPotentialStar and variables for Adam
    ! and dRHS,iHessian and dHessian

    ! Allocate allocatable arrays:
    allocate(dChemicalPotentialStar(nVar),dPi(nVar+i),dPiDot(nVar+i),&
    dLambda(nVar+i),dLambdaLast(nVar+i),dvAdam(nVar), dsAdam(nVar), dvAdamLast(nVar),&
    dsAdamLast(nVar),dxDot(nVar),dPotentialVector(nVar))
    allocate(dRHS(nVar+i), iHessian(nVar+i), dHessian(nVar+i,nVar+i))
    allocate(dRHSLast(nVar+i))
!
    ! Initialize variables:
    dRHS                   = 0D0
    dRHSLast               = 0D0
    dChemicalPotentialStar = 0D0
    dPi         = 0D0
    dPiDot      = 0D0
    dLambda     = 0D0
    dLambdaLast = 0D0
    dvAdam      = 0D0
    dsAdam      = 0D0
    dvAdamLast  = 0D0
    dsAdamLast  = 0D0
    dxDot       = 0D0
    dSubminGibbsEst = 0D0
    dPotentialVector= 0D0
    iAdamNeg = 0
    lNegativeFraction = .FALSE.
!
    ! Set a default value for iterSubMax if it is not specified:
    iterSub = 1
    if(dMaxElementPotential>1D-1) then
        dMaxPotentialTol = 1D-1
    elseif(dMaxElementPotential>1D-2) then 
        dMaxPotentialTol = 1D-3
    else
        dMaxPotentialTol = 1D-4
    end if
    iterSubMax = 1000
    ! dMaxPotentialTol = 1D-3
    ! Loop through all constituents in this solution phase:
    do k = 1, nVar
!
        ! Absolute species index:
        i = iFirstSUB + k -1 
!
        ! Compute the chemical potentials of all constituents defined by the element potentials:
        dChemicalPotentialStar(k) = 0D0
!
        do j = 1, nElements
            dChemicalPotentialStar(k) = dChemicalPotentialStar(k) + dElementPotential(j) &
                * dStoichSpecies(i,j)
        end do
        dChemicalPotentialStar(k) = dChemicalPotentialStar(k) / DFLOAT(iParticlesPerMole(i))
!
        

        ! Initialize the mole fractions:
        dMolFraction(i) = DMIN1(dMolFraction(i),1D0)
        dMolFraction(i) = DMAX1(dMolFraction(i), 1D-50)
!
!
    end do

    

    dMolFraction(iFirstSUB:iLastSUB)=dMolFraction(iFirstSUB:iLastSUB)/sum(dMolFraction(iFirstSUB:iLastSUB))
    
    ! Compute the chemical potentials of solution phase constituents:
    call SubMinChemicalPotential(iSolnPhaseIndex)

!
    ! Compute the driving force of this solution phase:
    call SubMinDrivingForce
!
    
    do j = 1, nVar
        i                 = iFirstSUB + j - 1 
        dSubminGibbsEst = dSubminGibbsEst+(dChemicalPotential(i)-dChemicalPotentialStar(j))
    end do
    dSubminGibbsEst = dSubminGibbsEst/nVar

    do j = 1, nVar
        i                 = iFirstSUB + j - 1
        dPotentialVector(j) = dSubminGibbsEst-(dChemicalPotential(i)-dChemicalPotentialStar(j))
    end do
    dMaxPotentialVector = MAXVAL(DABS(dPotentialVector))

    ! If this phase contains a miscibility gap, determine the absolute index of the corresponding
    ! solution phase with the miscibility gap:
    if (lMiscibility(iSolnPhaseIndex)) then
!
        ! Check the name of the solution phases:
        if (cSolnPhaseName(iSolnPhaseIndex) == cSolnPhaseName(iSolnPhaseIndex-1)) then
            ! The last solution phase in the indexing scheme corresponds to the miscibility gap:
            iSolnPhaseIndexOther = iSolnPhaseIndex - 1
        else
            ! Loop through all solution phases in the data-file system:
            LOOP_SolnSys: do i = 1, nSolnPhasesSys
!
                ! Cycle if this is same phase:
                if (i == iSolnPhaseIndex) cycle LOOP_SolnSys
!
                ! Check if these phases are the same:
                if (cSolnPhaseName(iSolnPhaseIndex) == cSolnPhaseName(i)) then
                    iSolnPhaseIndexOther = i
                    exit LOOP_SolnSys
                end if
            end do LOOP_SolnSys
        end if
    end if
!
end subroutine SubMinInit
!
!
