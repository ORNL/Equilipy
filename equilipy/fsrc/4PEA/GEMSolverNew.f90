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
    real(8)                             :: dGibbsGlobalLast, dGibbsGlobal
    real(8),dimension(32)               :: dGprogress
    logical                             :: lStableAssemblage

    ! Initialize variables
    iterGlobal      = 0
    iAssemblageLast = 0
    iterUBC         = 0
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
    LOOP_GEM: do iterGlobal = 1,1025
        ! Initialize variables related to history tracking
        iterHistoryLevel = 0
        dMolesPhaseHistory = 0d0
        dElementPotentialLast=dElementPotential
        dGibbsGlobalLast = dGibbsGlobal
        dToleranceLevel = -1D-8
        
!
        ! Part 1: Leveling
        LOOP_Leveling: do iterLevel = 1, 100
            
!
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
        dMaxElementPotential = MAXVAL(DABS(dElementPotential-dElementPotentialLast))
        dGprogress(k) = dGibbsGlobal-dGibbsGlobalLast
        
        if(MOD(k,32)==0) k=0
        k = k+1
        

        if((sum(ABS(dGprogress))<1D-15).AND.(iterGlobal>256)) then
            !Relax the tolerance when minimizer doesn't progress well
            dToleranceLevel = -FLOAT(iterGlobal)*1D-3
        end if
!       
        ! When a single phase appears and the amount converges to unity, exit the loop
        if((MINVAL(iAssemblage)==MAXVAL(iAssemblage)).AND.(MAXVAL(dMolesPhase)>=1000)) exit LOOP_GEM

        ! Global minimum
        if ((MINVAL(dPhasePotential) >dToleranceLevel).AND.&
        (dMaxElementPotential<1D-4)) exit LOOP_GEM
!
    end do LOOP_GEM
!
    ! POST PROCESSING
    call PostProcessPEA
    
    call CompFunctionNorm

    ! Define Candidate phases to be added 
    do i =1,nSpecies
        if(dPhasePotential(i)<500) iCandidate(i)=1
    end do
!
    
!
!
    return
!
end subroutine GEMSolverNew
!
!
