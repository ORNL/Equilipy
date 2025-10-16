
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CheckConvergence.f90
    !> \brief   Check convergence in the non-linear solver.
    !> \author  M.H.A. Piro
    !> \date    Apr. 26, 2012
    !> \sa      GEMSolver.f90
    !> \todo    Consider removing check for residuals of chemical potential terms...this may be redundant.
    !
    !
    ! References:
    ! ===========
    !
    !        M.H.A. Piro, "Computation of Thermodynamic Equilibria Pertinent to Nuclear Materials
    !        in Multi-Physics Codes," PhD Dissertation, Royal Military College of Canada, 2011.
    !
    !        M.H.A. Piro, T.M. Besmann, S. Simunovic, B.J. Lewis and W.T. Thompson, "Numerical
    !        Verification of Equilibrium Thermodynamic Computations in Nuclear Fuel Performance
    !        Codes," Journal of Nuclear Materials, 414 (2011) 399-407.
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   03/31/2011      M.H.A. Piro         Original code
    !   07/31/2011      M.H.A. Piro         Clean up code: remove unnecessary variables, update variable names
    !   09/12/2011      M.H.A. Piro         Added feature: check if dSumMolFractionSoln or dRelError is a NAN
    !   10/25/2011      M.H.A. Piro         Clean up code: modules, simplify programming
    !   10/27/2011      M.H.A. Piro         Modified check for Gibbs' Criteria: exclude dummy species
    !   04/05/2012      M.H.A. Piro         Clean up code: simplify relative error calculation of mass balance
    !                                        equations by using the dFunction vector.
    !   04/26/2012      M.H.A. Piro         Implementing Gibbs energy Minimization algorithm.
    !   08/22/2012      M.H.A. Piro         Add check for miscibility gap.
    !   01/31/2013      M.H.A. Piro         Apply appropriate check for residual/relative errors of the mass
    !                                        balance equations when the system component is an electron (i.e.,
    !                                        dMolesElement = 0D0).
    !   02/04/2013      M.H.A. Piro         Add a check to ensure that the Phase Rule is satisfied (this was
    !                                        necessarily satisfied in the previous version because charged
    !                                        phases were not considered.
    !   05/12/2014      M.H.A. Piro         Change the 9th check for a global minimum so that all metastable
    !                                        solution phases are checked for a global minimum, not just phases
    !                                        with a known miscibility gap.
    !   08/21/2015      M.H.A. Piro         Change tolerance for residuals of chemical potentials of individual
    !                                        species to only apply when x > 1.66D-24.  Practically, who cares
    !                                        if the residual of the chemical potential is large for a species
    !                                        whose mole fraction is 1D-50?!?  1.66D-24 is the fraction of a
    !                                        single atom in one mole.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to check if convergence has been achieved.  The conditions for
    !! thermodynamic equilibrium are (refer to above references for more details):
    !! <ol>
    !! <li> None of the phases in the estimated phase assemblage are "dummy species", which were introduced
    !!      by the ChemSage data-file, </li>
    !! <li> The number of moles of all species and phases are positive and non-zero. </li>
    !! <li> Gibbs' Phase Rule has been satisfied. </li>
    !! <li> The standard Gibbs energy of a pure condensed phase is not below the Gibbs Plane.
    !!      If so, this phase should be added to the assemblage. </li>
    !! <li> The residuals of the chemical potential terms are below a specified tolerance. </li>
    !! <li> The sum of site fractions on each sublattice must equal unity within tolerance. </li>
    !! <li> The driving force for every solution phase must be non-negative. </li>
    !! <li> The relative errors of the mass balance equations are within tolerance. </li>
    !! <li> The integral Gibbs energy of the system must be at a true global minimum and not
    !!       a local minima. </li>
    !! </ol>
    !!
    !! Each criterion listed above is sequentially tested and the system is considered converged
    !! when all are satisfied.  Control is returned to the PGESolver subroutine when any of the
    !! criterions has not been satisfied.  Note that the order of testing is done in a fashion
    !! that progressively increases computational expense.  For example, testing the mass balance
    !! constraints is the most computationally expensive task, which is why it is performed last.
    !!
    !! Note that the Gibbs Phase Rule is necessarily satisfied because iAssemblage is dimensioned by
    !! nElements.  Therefore, it is impossible for nPhases > nElements and the Gibbs Phase Rule is
    !! implicitly satisfied at all times.
    !!
    !
    ! Pertinent variables:
    ! ====================
    !
    ! lConverged                Logical variable: true if convedRelGibbsEnergyd, false if not converged.
    ! dStoichSpecies            Number of atoms of a particular element per formula mass of a particular
    !                           species (double matrix).
    ! iAssemblage               Integer vector storing the indices of phases contributing to the equilibrium
    !                           phase assemblage.
    ! iElement                  Integer vector storing the atomic numbers of elements in the system.
    ! dRelError                 Relative error of mass balance equations.
    ! dChemicalPotential        A double real vector representing the difference between the chemical potential
    !                           (defined by the element potentials) and the standard molar Gibbs energy of the
    !                           pure species divided by the total number of atoms per formula mass.
    ! dMolesPhase               The number of moles of a phase.
    ! dMolesElement             Total number of moles of an element.
    ! dEffStoichSolnPhase       The "effective" stoichiometry of a solution phase.
    ! dSumMolFractionSoln       Sum of mole fractions in a solution phase.
    ! dTolerance                Acceptable numerical tolerance (defined in Init subroutine).
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine CheckConvergence

    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver

    implicit none

    integer :: iMinMolesPhase, i, j, k, m, n
    logical :: lPhasePass


    ! Initialize variables
    dMaxElementPotential = MAXVAL(dabs(dElementPotential - dElementPotentialLast))
    dUpdateVarLast = dUpdateVar
    lPhaseChangeHistory(iterGlobal) = lPhaseChange

    
  
    ! lDebugMode = .True.



    ! Reset lSolnPhases to ensure it is correct
    lSolnPhases     = .FALSE.
    do i = 1, nSolnPhases
        if (iAssemblage(nElements+1-i) < 0) lSolnPhases(-iAssemblage(nElements+1-i)) = .TRUE.
    end do


    ! TEST #1: Check convergence
    ! -------------------------------------------------------------------------
    if (MINVAL(dPhasePotential) >=-1D-10) then
        if (dGEMFunctionNorm<1D-7) then
            ! Global minimum
            lPhaseChange    = .FALSE.
            lConverged = .TRUE.
            return
        else if((.not.lPhaseChange.and.&
            dMaxElementPotential<1D-6) .and.&
            (dGEMFunctionNorm<1D-6)) then
            !All positive and sufficiently progressed
            lPhaseChange    = .FALSE.
            lConverged = .TRUE.
            return
        end if
    end if

    ! Initialize variables:
    lConverged      = .FALSE.
    lPhaseChange    = .FALSE.

    ! TEST #2: Check if any of the phases in the assemblage are "dummy" phases:
    ! -------------------------------------------------------------------------
    if (lDebugMode) print *, "Test 2: Check dummy phases"
    LOOP_TEST1: do i = 1, nElements
        if (iAssemblage(i) > 0) then
            if (iPhase(iAssemblage(i)) < 0) then
                ! If there exist dummy phase, go through phase change
                iterHistory(1:nElements,iterGlobal) = iAssemblage(1:nElements)
                lConverged      = .FALSE.
                lPhaseChange    = .True.
                return
            end if
        end if
    end do LOOP_TEST1

    ! TEST #3: Check to make sure that the number of moles of all phases are non-negative.
    ! ------------------------------------------------------------------------------------
    dMolesPhase(nConPhases + 1 : nElements - nSolnPhases) = 0D0
    if (lDebugMode) print *, "Test 3: Check phase with negative mole ", MINVAL(dMolesPhase)
    if (MINVAL(dMolesPhase) < 0D0) then
        iterHistory(1:nElements,iterGlobal) = iAssemblage(1:nElements)
        lConverged      = .FALSE.
        lPhaseChange = .True.
        return
    end if

    ! TEST #4: Check if phase assemblaged needs to be changed
    ! -----------------------------------------------------------------------------------
    if (lDebugMode) print *, "Test 4: Negative driving force", dMaxElementPotential,dGEMFunctionNorm
    if(MINVAL(dPhasePotential)<-1D-10) then
        lPhaseChange    = .True.
        lConverged      = .FALSE.
        return
    end if


    ! TEST #5: Check if system progresses slowly (Need to come up with better implementation)
    ! Lagrangian some cases do not progress very well even though there is no need for phase change.
    ! One identified case is when one of the endmembers become small ~1D-7.
    ! There are some other cases as well, but could not identified the pattern yet.
    ! Implementing Entropy based minimization might help resolving issue for the first case
    ! -----------------------------------------------------------------------------------    
    if(dMaxElementPotential<1D-6.and.dGEMFunctionNorm-dGEMFunctionNormLast>=-1D-8) then
        if (lPostProcess.and. iterGlobal<3) then
            !if system diverges when calculating ProstProces, check phase assemblage again
            lPhaseChange = .True.
        else if(iterGlobal>=100.and.ALL(.not.lPhaseChangeHistory(iterGlobal-50:iterGlobal))) then
            lPhaseChange = .True.
        else if(iterGlobal>300.and.ALL(.not.lPhaseChangeHistory(iterGlobal-100:iterGlobal))) then
            lPhaseChange = .True.
        else if(iterGlobal>500.and.ALL(.not.lPhaseChangeHistory(iterGlobal-200:iterGlobal))) then
            lPhaseChange = .True.
        end if
    end if

    return

end subroutine CheckConvergence
