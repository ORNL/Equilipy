!
!
subroutine GEMSolverNew
!
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    GEMSolver.f90
    !> \brief   Gibbs Energy Minimization solver.
    !> \author  S.Y. Kwon
    !> \date    May 2, 2022
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   05/02/2022      S.Y. Kwon            Original code
    !   10/21/2023      S.Y. Kwon            Cleaning the code
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to compute the quantities of species and phases at thermodynamic
    !! equilibrium using the Gibbs Energy Minimization (GEM) method.  This subroutine uses the theory developed by
    !! Capitani 1983
    !! Issue: Immiscibility distinguising problem, Subminimization issue
    !
    ! Pertinent variables:
    ! ====================
    !
    ! nConPhases            The number of pure condensed phases in the assemblage
    ! nSolnPhases           The number of solution phases in the assemblage
    ! nSolnPhasesSys        The number of solution phases in the system
    ! iAssemblage           Integer vector containing the indices of phases in the assemblage
    !                        (1:nConphases represent pure condensed phases and (nElements-nSolnPhases:nSolnPhases)
    !                        represent solution phases.
    ! INFOThermo            An integer scalar identifying whether the program exits successfully or if
    !                        it encounters an error.
    !-------------------------------------------------------------------------------------------------------------
    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver
!
    implicit none
!
    integer                             :: i,j,k,l,m,n,o,p,q,iterLevel, nTempSolnPhases
    integer,dimension(nElements)        :: iAssemblageLast
    real(8)                             :: dGibbsGlobalLast, dGibbsGlobal, dTemp
    real(8),dimension(32)               :: dGprogress
    logical                             :: lStableAssemblage

    ! Initialize variables
    iterPEA         = 0
    iAssemblageLast = 0
    dGibbsGlobal    = 0
    dGibbsGlobalLast= 0
    dMaxElementPotential = 10
    dGprogress      = 1D0
    k               = 1
    
    
!
!
    ! Initializatrion: Allocate variables and calculate minimum points of each solution
    call InitGEMSolverNew
!
!
!
    LOOP_GEM: do iterPEA = 1,100
        ! Initialize variables related to history tracking
        iterHistoryLevel = 0
        dMolesPhaseHistory = 0d0
        dElementPotentialLast=dElementPotential
        dGibbsGlobalLast = dGibbsGlobal
        dToleranceLevel = -1D-8
        
!
        ! Part 1: Leveling
        LOOP_Leveling: do iterLevel = 1, 100
            
            ! Calculate phase potential of each species
            dPhasePotential = dChemicalPotential - MATMUL(dAtomFractionSpecies,dElementPotential)
            dMinPhasePotential = MINVAL(dPhasePotential)
            ! Check global minimum: if all elements in phase potential are positive, the system is in global minimum
            if ((dMinPhasePotential> dToleranceLevel)) exit LOOP_Leveling
!
            ! Determine the next phase assemblage to be tested:
            call GetNewAssemblage(iterLevel)
!
            ! Exit if an error has been encountered:
            if (INFOThermo /= 0) then

                call CompMinSolnPoint
                dPhasePotential = dChemicalPotential - MATMUL(dAtomFractionSpecies,dElementPotential)
                dMinPhasePotential = MINVAL(dPhasePotential)
                exit LOOP_GEM
            end if
!
        end do LOOP_Leveling
       
        
!       
        ! Part2: Calculate the minimum point of each solution phase based on ElementPotential
        call CompMinSolnPoint
!
        
        ! Part3: Check convergence of Capitani algorithm
        dGibbsGlobal = dot_product(dMolesElement,dElementPotential)
        dPhasePotential = dChemicalPotential - MATMUL(dAtomFractionSpecies,dElementPotential)
        dMinPhasePotential = MINVAL(dPhasePotential)
        dMaxElementPotential = MAXVAL(DABS(dElementPotential-dElementPotentialLast))

        ! When stable phases are only composed of stoichiomatric compounds
        CheckCompdOnly: if ((dMinPhasePotential>=-1D-7).AND.(dGEMFunctionNorm<1D-7)) then 
            do i = 1, nElements
                if(iAssemblage(i)>nSpecies) exit CheckCompdOnly
            end do 
            lCompbdOnly = .True.
            exit LOOP_GEM
        end if CheckCompdOnly

        ! When a single phase appears and the amount converges to unity, exit the loop
        if((MINVAL(iAssemblage)==MAXVAL(iAssemblage)).AND.(MAXVAL(dMolesPhase)>=1000)) exit LOOP_GEM

        ! Global minimum
        if ((dMinPhasePotential >=-dPEATol).AND.&
        (dMaxElementPotential<dPEATol).AND.(dGEMFunctionNorm<dPEATol)) exit LOOP_GEM
!
    end do LOOP_GEM
    ! print*, 'GEMSolverNew:',dGEMFunctionNorm, iterPEA, iterGlobal, dMinPhasePotential,dMaxElementPotential,INFOThermo
    ! print*, 'iAssembalge', iAssemblage
    ! print*, 'dMolesPhase', dMolesPhase

    ! do i =1,nElements
    !     l=iPhaseGEM(i)
    !     m = nSpeciesPhase(l-1) + 1      ! First constituent in phase.
    !     n = nSpeciesPhase(l)            ! Last  constituent in phase.
    !     dMolesSpecies(m:n)= dMolFractionGEM(i,m:n)* dMolesPhase(i)
    !     print*, 'dMolFractionGEM(i,m:n)', l, dMolFractionGEM(i,m:n)
    !     ! print*, 'dMolesSpecies(m:n)',l,dMolesSpecies(m:n)
    !     ! print*, 'dStoichSpecies(m:n,1)',dStoichSpecies(m:n,1)
    !     ! print*, 'dStoichSpecies(m:n,2)',dStoichSpecies(m:n,2)
    !     ! print*, 'dStoichSpecies(m:n,3)',dStoichSpecies(m:n,3)
    !     ! print*, 'dStoichSpecies(m:n,4)',dStoichSpecies(m:n,4)
    ! end do


    
    call PostProcessPEA
    ! print*, '============================================='
    ! print*, 'GEMSolverNew After:',dGEMFunctionNorm
    ! print*, 'iAssembalge', iAssemblage
    ! print*, 'dMolesPhase', dMolesPhase
    ! do i = 1, nSolnPhases
    !     k = -iAssemblage(nElements - i + 1)       ! Absolute solution phase index
    !     m = nSpeciesPhase(k-1) + 1      ! First constituent in phase.
    !     n = nSpeciesPhase(k)            ! Last  constituent in phase.
    !     ! dTemp = dTemp - sum(dStoichSpecies(m:n,j) * dMolesSpecies(m:n))
    !     print*,k,cSolnPhaseName(k)
    !     print*, 'dMolFraction(m:n)', dMolFraction(m:n)
    !     ! print*, 'dStoichSpecies(m:n,1)',dStoichSpecies(m:n,1)
    !     ! print*, 'dStoichSpecies(m:n,2)',dStoichSpecies(m:n,2)
    !     ! print*, 'dStoichSpecies(m:n,3)',dStoichSpecies(m:n,3)
    !     ! print*, 'dStoichSpecies(m:n,4)',dStoichSpecies(m:n,4)
    ! !     print*,k,cSolnPhaseName(k),dot_product(dMolFraction(m:n),dChemicalPotential(m:n))/&
    ! !     sum(matmul(dMolFraction(m:n)/DFLOAT(iParticlesPerMole(m:n)), dStoichSpecies(m:n,:))), dMolFraction(m:n)
    ! end do
    

    ! if (dGEMFunctionNorm>1D0) print*, iAssemblage
    
!
!
    return
!
end subroutine GEMSolverNew
!
!
