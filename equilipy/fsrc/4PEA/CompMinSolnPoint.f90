!
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
    integer :: i,j,k,l,m,n,nConstituents, iFloor
    real(8):: log2iter,dNormComponent
    logical :: lAddPhase
!
!
    l = MAX(1,nSolnPhasesSys)
!
    if (allocated(dMolFractionOld)) deallocate(dMolFractionOld)
    if (allocated(dEffStoichSolnPhase)) deallocate(dEffStoichSolnPhase)
    if (allocated(dSumMolFractionSoln)) deallocate(dSumMolFractionSoln)
    if (allocated(dDrivingForceSoln)) deallocate(dDrivingForceSoln)
    if (allocated(dMolesSpecies)) deallocate(dMolesSpecies)
    if (allocated(dPartialExcessGibbs)) then
        deallocate(dPartialExcessGibbs,dPartialEnthalpyXS,dPartialEntropyXS,dPartialHeatCapacityXS)
    end if
    
    allocate(&
        dMolesSpecies(nSpecies),&
        dEffStoichSolnPhase(l,nElements),&
        dSumMolFractionSoln(l),&
        dDrivingForceSoln(l),&
        dMolFractionOld(nSpecies),&
        dPartialExcessGibbs(nSpecies),&
        dPartialEnthalpyXS(nSpecies),&
        dPartialEntropyXS(nSpecies),&
        dPartialHeatCapacityXS(nSpecies))
    ! dChemicalPotential(nSpecies+1:)=5D9
    dPartialExcessGibbs=0D0
    dPartialEnthalpyXS=0D0
    dPartialEntropyXS=0D0
    dPartialHeatCapacityXS=0D0
!
!
    dMolFractionOld=dMolFraction
    dChemicalPotentialOld = dChemicalPotential(:nSpecies)
    LOOP_Soln: do i = 1, nSolnPhasesSys
        !initialize variables
        m                = nSpeciesPhase(i-1) + 1      ! First constituent in phase.
        n                = nSpeciesPhase(i)            ! Last  constituent in phase.
        nConstituents    = n - m + 1
        lAddPhase        = .False.        
!
        if (.NOT.lMiscibility(i)) then
            ! Normal soln or first soln phase when immiscibility is considered
            ! Initial mole fraction is given from previous iteration
!
!
            ! Perform subminimization:
            call Subminimization(i, lAddPhase)
            ! print*,'Composition',i,m,n,dMolFraction(m:n),dDrivingForceSoln(i)
!
            call CompStoichSolnPhase(i)
!
!
            dAtomFractionSpecies(nSpecies+i,:) = dEffStoichSolnPhase(i,:)/sum(dEffStoichSolnPhase(i,:))
            dStoichSpeciesLevel(nSpecies+i,:)  = dEffStoichSolnPhase(i,:)
            iPhaseLevel(nSpecies+i)            = i
            dChemicalPotential(nSpecies+i)     = dDrivingForceSoln(i) + &
            dot_product(dAtomFractionSpecies(nSpecies+i,:),dElementPotential(:))

            do j = 1,nElements
                if(iAssemblage(j)-nSpecies+1==i) then
                    dMolFractionGEM(j,m:n) = dMolFraction(m:n)
                end if
            end do
!
        else
            !Second phase when immiscibility is considered
            ! imod = INT(10*log(float(iterGlobal))/log(2.0))
            log2iter=log(float(iterGlobal))/log(2.0)
            iFloor = FLOOR(log2iter)
!
            ! if((iterGlobal<5).OR.(imod==5).OR.(dMaxElementPotential<5D-4)) then
            if((iterGlobal<16).OR.((log2iter-iFloor<1D-20).AND.(iFloor<4))) then
                ! print*, 'iterGlobal',iterGlobal
                call CheckMiscibilityGap(i,lAddPhase)
                dAtomFractionSpecies(nSpecies+i,:) = dEffStoichSolnPhase(i,:)/sum(dEffStoichSolnPhase(i,:))
                dStoichSpeciesLevel(nSpecies+i,:)  = dEffStoichSolnPhase(i,:)
                iPhaseLevel(nSpecies+i)            = i
                
                if(lAddPhase) then
                    dChemicalPotential(nSpecies+i)     = dDrivingForceSoln(i) + &
                    dot_product(dAtomFractionSpecies(nSpecies+i,:),dElementPotential(:))
                    
                else
                    dChemicalPotential(nSpecies+i)     = 9D5
                end if
            else
                ! Perform subminimization:
                call Subminimization(i, lAddPhase)
                ! print*,'Composition',i,m,n,dMolFraction(m:n),dDrivingForceSoln(i)
                
    !
                call CompStoichSolnPhase(i)
    !
    !
                dAtomFractionSpecies(nSpecies+i,:) = dEffStoichSolnPhase(i,:)/sum(dEffStoichSolnPhase(i,:))
                dStoichSpeciesLevel(nSpecies+i,:)  = dEffStoichSolnPhase(i,:)
                iPhaseLevel(nSpecies+i)            = i
                dChemicalPotential(nSpecies+i)     = dDrivingForceSoln(i) + &
                dot_product(dAtomFractionSpecies(nSpecies+i,:),dElementPotential(:))

                do j = 1,nElements
                    if(iAssemblage(j)-nSpecies+1==i) then
                        dMolFractionGEM(j,m:n) = dMolFraction(m:n)
                    end if
                end do
            end if
!
            cycle LOOP_Soln
!
        end if
    end do LOOP_Soln

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
