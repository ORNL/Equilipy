
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    GEMSolver.f90
    !> \brief   Gibbs Energy Minimization solver.
    !> \author  M.H.A. Piro
    !> \date    Apr. 25, 2012
    !> \sa      Thermochimica.f90
    !> \sa      InitGEMSolver.f90
    !> \sa      GEMNewton.f90
    !> \sa      GEMLineSearch.f90
    !> \sa      CheckPhaseAssemblage.f90
    !> \sa      CheckConvergence.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   04/25/2012      M.H.A. Piro         Original code
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to compute the quantities of species and phases at thermodynamic
    !! equilibrium using the Gibbs Energy Minimization (GEM) method.  This subroutine uses values of
    !! dMolesPhase, dChemicalPotential and iAssemblage from the Leveling and PostLeveling subroutines as initial
    !! estimates for computation.
    !!
    !! The main subroutines used by this solver are summarized below:
    !! <table border="1" width="800">
    !! <tr>
    !!    <td> <b> File name </td> <td> Description </b> </td>
    !! </tr>
    !! <tr>
    !!    <td> InitGEMSolver.f90 </td>
    !!    <td> Initialize the GEMSolver by establishing the initial phase assemblage and composition.  </td>
    !! </tr>
    !! <tr>
    !!    <td> CheckSysOnlyPureConPhases.f90 </td>
    !!    <td> Check the system if there are only pure condensed phases. The system may already be converged.</td>
    !! </tr>
    !! <tr>
    !!    <td> GEMNewton.f90 </td>
    !!    <td> Compute the direction vector using Newton's method.  </td>
    !! </tr>
    !! <tr>
    !!    <td> GEMLineSearch.f90 </td>
    !!    <td> Perform a line search along the direction vector.  </td>
    !! </tr>
    !! <tr>
    !!    <td> CheckPhaseAssemblage.f90 </td>
    !!    <td> Check if the phase assemblage needs to be adjusted.  </td>
    !! </tr>
    !! <tr>
    !!    <td> CheckConvergence.f90 </td>
    !!    <td> Check if the system has converged.  </td>
    !! </tr>
    !! </table>
    !
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
    ! INFO                  An integer scalar identifying an error from LAPACK.  This is used by the GEMNewton
    !                        subroutine to indicate whether there is a singularity in the Hessian matrix.
    ! lConverged            A logical variable indicating whether the code has convered (.TRUE.) or not (.FALSE.).
    ! lRevertSystem         A logical scalar indicating whether the system should be reverted to a previously
    !                        successful phase assemblage.
    ! dTolerance            A double real vector representing numerical tolerances (defined in InitThermo.f90).
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine GEMSolver

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    lPhaseChange = .True.
    lCompbdOnly = .False.
    
    call InitGEMSolver

    ! Leveling solely based on stoichiometric compounds and solution endmembers
    call LevelingSolver
    
    if(nSolnPhasesSys==0.or.nElements==1.or. (MINVAL(dPhasePotential)>-1D-13.and.dGEMFunctionNorm<1D-10)) then
        ! If there is no solution phases in the system, return leveling resul
        lCompbdOnly = .True.
        return   
    end if

    ! If there exist solution phases in the system, make sure global minimum
    call MultiPhaseMinimizer

    return

end subroutine GEMSolver
