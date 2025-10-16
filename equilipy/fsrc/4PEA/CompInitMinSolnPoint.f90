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
    integer :: i, j, k, l, n, m, nConstituents
    real(8) :: dMinMoleFraction, dMaxMoleFraction, dNormComponent
    logical :: lAddPhase
    integer, dimension(:), allocatable:: iIndex, iEndmemberPotential
    real(8), dimension(:), allocatable :: dDrivingForceTemp, dEndmemberPotential
    real(8), dimension(:,:), allocatable :: dAtomFractionTemp,dStoichSpeciesTemp, dMolFractionTemp
!
!
    l = MAX(1,nSolnPhasesSys)
!
    if (allocated(dMolFractionOld))        deallocate(dMolFractionOld)
    if (allocated(dEffStoichSolnPhase)) deallocate(dEffStoichSolnPhase)
    if (allocated(dSumMolFractionSoln)) deallocate(dSumMolFractionSoln)
    if (allocated(dDrivingForceSoln))   deallocate(dDrivingForceSoln)
    if (allocated(dMolesSpecies))       deallocate(dMolesSpecies)
    if (allocated(dPartialExcessGibbs)) then
        deallocate(dPartialExcessGibbs,dPartialEnthalpyXS,dPartialEntropyXS,dPartialHeatCapacityXS)
    end if
    allocate(&
        dMolesSpecies(nSpecies),&
        dMolFractionOld(nSpecies),&
        dEffStoichSolnPhase(l,nElements),&
        dSumMolFractionSoln(l),&
        dDrivingForceSoln(l),&
        dPartialExcessGibbs(nSpecies),&
        dPartialEnthalpyXS(nSpecies),&
        dPartialEntropyXS(nSpecies),&
        dPartialHeatCapacityXS(nSpecies))
    dPartialExcessGibbs=0D0
    dPartialEnthalpyXS=0D0
    dPartialEntropyXS=0D0
    dPartialHeatCapacityXS=0D0
    dMinMoleFraction = 1D-5
    
!
!
!   
  
    dMolFractionOld=dMolFraction
    LOOP_Soln: do i = 1, nSolnPhasesSys
        !initialize variables
        m                = nSpeciesPhase(i-1) + 1      ! First constituent in phase.
        n                = nSpeciesPhase(i)            ! Last  constituent in phase.
        nConstituents    = n - m + 1
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

            if(allocated(dEndmemberPotential)) deallocate(dEndmemberPotential)
            if(allocated(iEndmemberPotential)) deallocate(iEndmemberPotential)
            allocate(dEndmemberPotential(nConstituents),iEndmemberPotential(nConstituents))
            allocate(&
            iIndex(nElements+1),&
            dDrivingForceTemp(nElements+1),&
            dAtomFractionTemp(nElements+1,nElements),&
            dStoichSpeciesTemp(nElements+1,nElements),&
            dMolFractionTemp(nElements+1, nConstituents))


            dEndmemberPotential= dStdGibbsEnergy(m:n)- MATMUL(dAtomFractionSpecies(m:n,:),dElementPotential)
            call Qsort(dEndmemberPotential,iEndmemberPotential,nConstituents)
!
            iIndex             = 1
            dDrivingForceTemp  = 9E5
            dAtomFractionTemp  = 0D0
            dStoichSpeciesTemp = 0D0
            dMolFractionTemp   = 0D0
!
            ! Perform subminimization multiple times by initializing from all extremums of the domain space:
            LOOP_Constituents1: do j = 1, nElements+1
            
                if(j>1) then
                    dMolFraction(m:n)   = dMinMoleFraction
                    dMolFraction(m+iEndmemberPotential(j-1)-1) = dMaxMoleFraction
                end if
                ! Perform subminimization:
                call Subminimization(i, lAddPhase)
!
                ! Store info of all local minima
                call CompStoichSolnPhase(i)
                dDrivingForceTemp(j)   = dDrivingForceSoln(i)
                dMolFractionTemp(j,:)  = dMolFraction(m:n)
                dAtomFractionTemp(j,:) = dEffStoichSolnPhase(i,:)/sum(dEffStoichSolnPhase(i,:))
                dStoichSpeciesTemp(j,:)= dEffStoichSolnPhase(i,:)
            end do LOOP_Constituents1
!
            ! Sort the results according to the driving force (ascending order)
            ! note that output dDrivingForceTemp is the sorted version
            call Qsort(dDrivingForceTemp,iIndex,nElements+1)
            ! print*, cSolnPhaseName(i),iIndex, dMolFractionTemp(iIndex(1),:)
!
            ! Update variables for the current phase 
!
            dAtomFractionSpecies(nSpecies+i,:) = dAtomFractionTemp(iIndex(1),:)
            dStoichSpeciesLevel(nSpecies+i,:)  = dStoichSpeciesTemp(iIndex(1),:)
            dMolFraction(m:n) = dMolFractionTemp(iIndex(1),:)
            dChemicalPotential(nSpecies+i)     = dDrivingForceTemp(1) + &
            dot_product(dAtomFractionSpecies(nSpecies+i,:),dElementPotential(:))
            iPhaseLevel(nSpecies+i)            = i
! !
            ! Update variables for immiscible phases
            LOOP_immiscibleSoln: do j = 2, nConstituents
!
                k = i + j - 1
                m  = nSpeciesPhase(k-1) + 1      ! First constituent in phase.
                n  = nSpeciesPhase(k)            ! Last  constituent in phase.
!
                !The following solution phases must have True value for lMiscibility
                if (lMiscibility(k).AND.(k<=nSolnPhasesSys).AND.(cSolnPhaseName(i)==cSolnPhaseName(k))) then
!
                    ! Store the local minima if it is not equal to global minimum
!
                    dAtomFractionSpecies(nSpecies+k,:) = dAtomFractionTemp(iIndex(j),:)
                    dStoichSpeciesLevel(nSpecies+k,:)  = dStoichSpeciesTemp(iIndex(j),:)
                    dMolFraction(m:n) = dMolFractionTemp(iIndex(j),:)
                    dChemicalPotential(nSpecies+k) = dDrivingForceTemp(j) + &
                    dot_product(dAtomFractionSpecies(nSpecies+k,:),dElementPotential(:))
                    iPhaseLevel(nSpecies+k)            = k
                end if
!
            end do LOOP_immiscibleSoln
        end if
!
!
!
    end do LOOP_Soln
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
