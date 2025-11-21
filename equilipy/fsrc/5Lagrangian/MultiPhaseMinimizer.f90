


subroutine MultiPhaseMinimizer

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none
    integer::   INFO, k, j, m, n
    
    iterGlobal = 0

    LOOP_GEMSolver: do iterGlobal = 1, iterGlobalMax
    ! LOOP_GEMSolver: do iterGlobal = 1, 50
        ! 1. Correct phase assemblage according to UBC algorithm if necessary
        if(lPhaseChange) then
            call GEMSolverNew
            ! print*,'dGEMFunctionNorm',dGEMFunctionNorm
            ! if ((dGEMFunctionNorm<1D-7.and.MINVAL(dPhasePotential) >=-1D-7).or.&
            if ((nSolnPhases==0 .and.dGEMFunctionNorm<1D-7 ).or.lCompbdOnly) return
        end if
        ! print*, iterGlobal, dMolesPhase
        
        ! 2. Gradien decent using Lagrangian multiplier method
        ! 2.1 Construct the Hessian matrix and compute the direction vector:
        call GEMNewton(INFO)

        ! 2.2 Perform a line search using the direction vector:
        call GEMLineSearch

        ! Check convergence: 
        call CheckConvergence
        ! if(lPhaseChange) print*, iterGlobal, iAssemblage

        if (lConverged) exit LOOP_GEMSolver 
             

    end do LOOP_GEMSolver
    ! print*,'iterGlobal',iterGlobal
   
    
    
    ! if (iterGlobal>=1000) then
    !     print*, 'Maxing out'
    ! end if
    ! if(dGEMFunctionNorm>1D-7) print*, iterGlobal, dTemperature, cElementName,dMolesElement
    ! if(dGEMFunctionNorm>1D-7) print*, iterGlobal, dGEMFunctionNorm
    ! print*, iterGlobal,iterPEA, dGEMFunctionNorm, dMolesPhase
    ! print*, iAssemblage
    ! if(iterGlobal>500) print*, iterGlobal, dGEMFunctionNorm
    DEALLOCATE(dElementPotentialLast,dMolesSpeciesLast,dMolesPhaseLast)
    
    

    return

end subroutine MultiPhaseMinimizer
