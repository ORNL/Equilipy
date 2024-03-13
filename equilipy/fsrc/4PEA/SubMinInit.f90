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
    USE ModuleGEMSolver, ONLY: lMiscibility, dDrivingForceSoln
!
    implicit none
!
    integer::  i, j, k, iSolnPhaseIndex, iterSubMax
!
!
    ! Initialize variables:
    dDrivingForceLast = 10D0
    lSubMinConverged  = .FALSE.
!
    dDrivingForceSoln(iSolnPhaseIndex) = 0D0
!
!
    ! Check if allocatable arrays are already allocated.  Deallocate if necessary:
    if (allocated(dChemicalPotentialStar)) deallocate(dChemicalPotentialStar)
    if (allocated(dRHS))                   deallocate(dRHS)
    if (allocated(dHessian))               deallocate(dHessian)
    if (allocated(iHessian))               deallocate(iHessian)
!
    ! Allocate allocatable arrays:
    allocate(dChemicalPotentialStar(nVar))
!
    ! Determine prefactor for allocating arrays (depends on whether the phase is ionic):
    if (iPhaseElectronID(iSolnPhaseIndex) == 0) then
        i = 1
    else
        i = 2
    end if
!
    ! Allocate allocatable arrays:
    allocate(dRHS(nVar+i), iHessian(nVar+i), dHessian(nVar+i,nVar+i))
!
    ! Initialize variables:
    dRHS                   = 0D0
    dChemicalPotentialStar = 0D0
!
    ! Set a default value for iterSubMax if it is not specified:
    iterSubMax = 100
!
    ! Loop through all constituents in this solution phase:
    do k = 1, nVar
!
        ! Absolute species index:
        i = nSpeciesPhase(iSolnPhaseIndex-1) + k
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
        ! ! Mole fractions are often given before proceeding Subminimization.
        ! ! Define the mole fractions in case they are not given or given as zeros.
        ! if((MINVAL(DABS(dElementPotential))<1E-15).AND.&
        ! (.NOT.lMiscibility(iSolnPhaseIndex)).AND.&
        ! (MAXVAL(dMolfraction(iFirstSUB:iLastSUB))<1E-15)) then
        !     !Subminimization during Leveling solver
        !     dMolFraction(i) = DEXP(dChemicalPotentialStar(k) - dStdGibbsEnergy(i))
        !     dMolFraction(i) = DMIN1(dMolFraction(i),1D0)
        !     dMolFraction(i) = dMAX1(dMolFraction(i),1E-15)
        ! else
        !     ! Take the MolFraction as it is given
        !     dMolFraction(i) = DMIN1(dMolFraction(i),1D0)
        !     dMolFraction(i) = dMAX1(dMolFraction(i),1E-15)
        ! end if
        ! Initialize the mole fractions:
        dMolFraction(i) = DMAX1(dMolFraction(i), 1D-15)
!
!
    end do
    ! dMolFraction(iFirstSUB:iLastSUB)=dMolFraction(iFirstSUB:iLastSUB)/sum(dMolFraction(iFirstSUB:iLastSUB))
!
    ! Compute the chemical potentials of solution phase constituents:
    call SubMinChemicalPotential(iSolnPhaseIndex)
!
    ! Compute the driving force of this solution phase:
    call SubMinDrivingForce
!
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
