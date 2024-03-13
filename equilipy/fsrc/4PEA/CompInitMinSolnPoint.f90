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
    real(8) :: dMinMoleFraction, dMaxMoleFraction
    logical :: lAddPhase
    integer, dimension(:), allocatable:: iIndex
    real(8), dimension(:), allocatable :: dDrivingForceTemp
    real(8), dimension(:,:), allocatable :: dAtomFractionTemp,dStoichSpeciesTemp
!
!
    l = MAX(1,nSolnPhasesSys)
!
    if (allocated(dMolFraction))        deallocate(dMolFraction)
    if (allocated(dEffStoichSolnPhase)) deallocate(dEffStoichSolnPhase)
    if (allocated(dSumMolFractionSoln)) deallocate(dSumMolFractionSoln)
    if (allocated(dDrivingForceSoln))   deallocate(dDrivingForceSoln)
    if (allocated(dMolesSpecies))       deallocate(dMolesSpecies)
    if (allocated(dPartialExcessGibbs)) then
        deallocate(dPartialExcessGibbs,dPartialEnthalpyXS,dPartialEntropyXS,dPartialHeatCapacityXS)
    end if
    allocate(&
        dMolesSpecies(nSpecies),&
        dMolFraction(nSpecies),&
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
!
!
!
    LOOP_Soln: do i = 1, nSolnPhasesSys
        !initialize variables
        m                = nSpeciesPhase(i-1) + 1      ! First constituent in phase.
        n                = nSpeciesPhase(i)            ! Last  constituent in phase.
        nConstituents    = n - m + 2
        dMinMoleFraction = 1D-15
        dMaxMoleFraction = 1D0 - dMinMoleFraction
!
        lAddPhase        = .False.        
!
        if (.NOT.lMiscibility(i)) then
            ! Normal soln or first soln phase when immiscibility is considered
            ! Global minimum point of the solution phase is computed
            allocate(&
                iIndex(nConstituents),&
                dDrivingForceTemp(nConstituents),&
                dAtomFractionTemp(nConstituents,nElements),&
                dStoichSpeciesTemp(nConstituents,nElements)&
            )         
!
            iIndex             = 1
            dDrivingForceTemp  = 9E5
            dAtomFractionTemp  = 0D0
            dStoichSpeciesTemp = 0D0
!
            ! Perform subminimization multiple times by initializing from all extremums of the domain space:
            LOOP_Constituents1: do j = 1, nConstituents
!
                ! Initialize the mole fractions:
                if (j==1)then
                    dMolFraction(m:n)   = 0D0
                else
                    dMolFraction(m:n)   = dMinMoleFraction
                    dMolFraction(m+j-1) = dMaxMoleFraction
                end if
!
!
                ! Perform subminimization:
                call Subminimization(i, lAddPhase)
!
                ! Store info of all local minima
                call CompStoichSolnPhase(i)
                dDrivingForceTemp(j)   = dDrivingForceSoln(i)
                dAtomFractionTemp(j,:) = dEffStoichSolnPhase(i,:)/sum(dEffStoichSolnPhase(i,:))
                dStoichSpeciesTemp(j,:)= dEffStoichSolnPhase(i,:)
            end do LOOP_Constituents1
!
            ! Sort the results according to the driving force (ascending order)
            ! note that output dDrivingForceTemp is the sorted version
            call Qsort(dDrivingForceTemp,iIndex,nConstituents)
!
            ! Update variables for the current phase 
!
            dAtomFractionSpecies(nSpecies+i,:) = dAtomFractionTemp(iIndex(1),:)
            dStoichSpeciesLevel(nSpecies+i,:)  = dStoichSpeciesTemp(iIndex(1),:)
            ! SY: Had an issue with Compound Energy Formalism, so just set to zero for now
            dChemicalPotential(nSpecies+i)     = 0
            ! dChemicalPotential(nSpecies+i)     = dDrivingForceTemp(1) + &
            ! dot_product(dAtomFractionSpecies(nSpecies+i,:),dElementPotential(:))
            iPhaseLevel(nSpecies+i)            = i
!
            ! Update variables for immiscible phases
            LOOP_immiscibleSoln: do j = 2, nConstituents
!
                k = i + j - 1
!
                !The following solution phases must have True value for lMiscibility
                if (lMiscibility(k).AND.(k<=nSolnPhasesSys)) then
!
                    ! Store the local minima if it is not equal to global minimum
!
                    dAtomFractionSpecies(nSpecies+k,:) = dAtomFractionTemp(iIndex(j),:)
                    dStoichSpeciesLevel(nSpecies+k,:)  = dStoichSpeciesTemp(iIndex(j),:)
                    ! SY: Had an issue with Compound Energy Formalism, so just set to zero for now
                    dChemicalPotential(nSpecies+i)     = 0
                    ! dChemicalPotential(nSpecies+k) = dDrivingForceTemp(j) + &
                    ! dot_product(dAtomFractionSpecies(nSpecies+k,:),dElementPotential(:))
                    iPhaseLevel(nSpecies+k)            = k
                else
                    exit LOOP_immiscibleSoln
                end if
!
            end do LOOP_immiscibleSoln
!
            deallocate(iIndex,dDrivingForceTemp,dAtomFractionTemp,dStoichSpeciesTemp)
        else
            !Second phase when immiscibility is considered
            !All regarding values were set in the previous condition.
            cycle LOOP_Soln
!
        end if
!
!
!
    end do LOOP_Soln
!
    !After Subminimization, the chemical potential of solution components changes.
    !These need to be reverted back to pure
    dChemicalPotential(:nSpecies)     = dChemicalPotentialOld
!
    ! deallocate(dEffStoichSolnPhase,dSumMolFractionSoln,&
    ! dDrivingForceSoln,dPartialExcessGibbs)
!
end subroutine CompInitMinSolnPoint
!
!
!
