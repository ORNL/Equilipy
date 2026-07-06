
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
    !   06/24/2026      S.Y. Kwon           Evaluate phase-potential active-set triggers only after the fixed
    !                                        Lagrangian residual satisfies the convergence tolerance.
    !   06/24/2026      S.Y. Kwon           Ignored internal disordered helpers of active ordered phases in
    !                                        phase-potential active-set checks.
    !   06/24/2026      S.Y. Kwon           Tested a tighter convergence residual for active SUBOM phases.
!   06/25/2026      S.Y. Kwon           Removed the SUBOM-only pseudo-endmember tolerance; CEF order
!                                        coordinates are handled by site-fraction minimization.
!   06/25/2026      S.Y. Kwon           Tightened fixed-assemblage postprocess convergence for derivative
!                                        properties.
!   06/25/2026      S.Y. Kwon           Refreshed phase potentials from the current elemental-potential plane
!                                        before testing active-set phase-potential triggers.
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

    integer :: iMinMolesPhase, i, j, k, m, n, l, nTraceResidualCheck
    integer :: iPotential, iPotentialPhase, iActiveSoln
    real(8) :: dTraceResidual, dMinPotentialForPhaseCheck, dConvergenceNormTolerance
    logical :: lPhasePass, lSkipPhasePotential


    ! Initialize variables
    dMaxElementPotential = MAXVAL(dabs(dElementPotential - dElementPotentialLast))
    dUpdateVarLast = dUpdateVar
    lPhaseChangeHistory(iterGlobal) = .FALSE.
    iPhaseChangeReasonHistory(iterGlobal) = PHASE_CHANGE_REASON_NONE

    
  
    ! lDebugMode = .True.

    ! Initialize variables:
    lConverged      = .FALSE.
    lPhaseChange    = .FALSE.
    iPhaseChangeReason = PHASE_CHANGE_REASON_NONE
    dConvergenceNormTolerance = 1D-5
    if (lPostProcess) dConvergenceNormTolerance = 1D-10

    ! Reset lSolnPhases to ensure it is correct
    lSolnPhases     = .FALSE.
    do i = 1, nSolnPhases
        if (iAssemblage(nElements+1-i) < 0) then
            j = -iAssemblage(nElements+1-i)
            lSolnPhases(j) = .TRUE.
        end if
    end do


    ! TEST #1: Check to make sure that the number of moles of all phases are non-negative.
    ! ------------------------------------------------------------------------------------
    dMolesPhase(nConPhases + 1 : nElements - nSolnPhases) = 0D0
    if (lDebugMode) print *, "Test 1: Check phase with negative mole ", MINVAL(dMolesPhase)
    if (MINVAL(dMolesPhase) < 0D0) then
        iterHistory(1:nElements,iterGlobal) = iAssemblage(1:nElements)
        lConverged      = .FALSE.
        lPhaseChange = .True.
        iPhaseChangeReason = PHASE_CHANGE_REASON_NEGATIVE_PHASE_AMOUNT
        lPhaseChangeHistory(iterGlobal) = .TRUE.
        iPhaseChangeReasonHistory(iterGlobal) = iPhaseChangeReason
        return
    end if

    ! TEST #2: GEM stagnation is handled in RunLagrangianGEM by monitoring
    ! relative progress in dGEMFunctionNorm.  CheckConvergence only handles
    ! thermodynamic convergence and explicit phase-assemblage triggers below.
    ! ------------------------------------------------------------------------------------

    ! TEST #3: Check Lagrangian Convergence
    ! -------------------------------------------------------------------------
    if(dMaxElementPotential>1D-6) return

    ! TEST #5: Check if any of the phases in the assemblage are "dummy" phases:
    ! -------------------------------------------------------------------------
    if (lDebugMode) print *, "Test 5: Check dummy phases"
    LOOP_TEST1: do i = 1, nElements
        if (iAssemblage(i) > 0) then
            if (iPhase(iAssemblage(i)) < 0) then
                ! If there exist dummy phase, go through phase change
                iterHistory(1:nElements,iterGlobal) = iAssemblage(1:nElements)
                lConverged      = .FALSE.
                lPhaseChange    = .True.
                iPhaseChangeReason = PHASE_CHANGE_REASON_DUMMY_PHASE
                lPhaseChangeHistory(iterGlobal) = .TRUE.
                iPhaseChangeReasonHistory(iterGlobal) = iPhaseChangeReason
                return
            end if
        end if
    end do LOOP_TEST1

    ! TEST #6: Check Phase Potentials: Global Minimum
    ! -------------------------------------------------------------------------
    if (dGEMFunctionNorm<dConvergenceNormTolerance) then
        if (lTraceSpeciesControlEnabled .AND. (.NOT.lPostProcess)) then
            nTraceResidualCheck = 0
            do l = 1, nSolnPhases
                k = -iAssemblage(nElements - l + 1)
                m = nSpeciesPhase(k-1) + 1
                n = nSpeciesPhase(k)
                do i = m, n
                    if (allocated(lTraceSpeciesInactive)) then
                        if (lTraceSpeciesInactive(i)) cycle
                    end if
                    if (.NOT.allocated(lTraceSpeciesReinjected)) cycle
                    if (.NOT.lTraceSpeciesReinjected(i)) cycle
                    if ((dMolFraction(i) > dTraceSpeciesRemoveFraction).AND.&
                        (dMolFraction(i) < 1D-5)) then
                        nTraceResidualCheck = nTraceResidualCheck + 1
                    end if
                end do
            end do

            if (nTraceResidualCheck == 1) then
                do l = 1, nSolnPhases
                    k = -iAssemblage(nElements - l + 1)
                    m = nSpeciesPhase(k-1) + 1
                    n = nSpeciesPhase(k)
                    do i = m, n
                        if (allocated(lTraceSpeciesInactive)) then
                            if (lTraceSpeciesInactive(i)) cycle
                        end if
                        if (.NOT.allocated(lTraceSpeciesReinjected)) cycle
                        if (.NOT.lTraceSpeciesReinjected(i)) cycle
                        if ((dMolFraction(i) > dTraceSpeciesRemoveFraction).AND.&
                            (dMolFraction(i) < 1D-5)) then
                            dTraceResidual = 0D0
                            do j = 1, nElements
                                dTraceResidual = dTraceResidual + dElementPotential(j) * dStoichSpecies(i,j)
                            end do
                            dTraceResidual = dTraceResidual / DFLOAT(iParticlesPerMole(i))
                            dTraceResidual = DABS(dChemicalPotential(i) - dTraceResidual)
                            if (dTraceResidual > dTraceSpeciesResidualTolerance) return
                        end if
                    end do
                end do
            end if
        end if
        ! The fixed active set has converged.  Only now interpret the stored
        ! phase potentials as an active-set trigger.  While the Lagrangian
        ! residual is still large, a negative phase potential can be stale
        ! Leveling/PEA information and should not mask a line-search failure.
        dPhasePotential = dLevelingChemicalPotential - &
            MATMUL(dLevelingCompositionSpecies,dElementPotential)
        dMinPotentialForPhaseCheck = HUGE(1D0)
        do iPotential = 1, nSpecies
            iPotentialPhase = iPhase(iPotential)
            lSkipPhasePotential = .FALSE.

            if ((iPotentialPhase > 0).AND.allocated(iDisorderedPhase)) then
                do l = 1, nSolnPhases
                    iActiveSoln = -iAssemblage(nElements - l + 1)
                    if ((iActiveSoln > 0).AND.(iActiveSoln <= SIZE(iDisorderedPhase))) then
                        if (iDisorderedPhase(iActiveSoln) == iPotentialPhase) then
                            lSkipPhasePotential = .TRUE.
                            exit
                        end if
                    end if
                end do
            end if

            if (lSkipPhasePotential) cycle
            dMinPotentialForPhaseCheck = DMIN1(dMinPotentialForPhaseCheck, dPhasePotential(iPotential))
        end do

        if ((.not. lPostProcess) .AND. (dMinPotentialForPhaseCheck <-1D-10)) then
            lPhaseChange    = .TRUE.
            lConverged      = .FALSE.
            iPhaseChangeReason = PHASE_CHANGE_REASON_PHASE_POTENTIAL
            lPhaseChangeHistory(iterGlobal) = .TRUE.
            iPhaseChangeReasonHistory(iterGlobal) = iPhaseChangeReason
            return
        end if

        ! Global minimum
        lPhaseChange    = .FALSE.
        lConverged = .TRUE.
        return
    else
        ! print*,'Lagrangian issue:','dGEMFunctionNorm',dGEMFunctionNorm
    end if

    return

end subroutine CheckConvergence
