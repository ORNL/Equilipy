!
!
subroutine LevelingSolver
!
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    LevelingSolver.f90
    !> \brief   A linear solver that estimates thermodynamic equilibrium.
    !> \author  M.H.A. Piro
    !> \date    Apr. 24, 2012
    !> \sa      Thermochimica.f90
    !> \sa      GetFirstAssemblage.f90
    !> \sa      GetNewAssemblage.f90
    !> \sa      PostLevelingSolver.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !
    !   03/31/2011      M.H.A. Piro         Original code
    !   07/29/2011      M.H.A. Piro         Cleaned up code
    !   07/31/2011      M.H.A. Piro         Clean up code: remove unnecessary variables, update variable names
    !   10/20/2011      M.H.A. Piro         Clean up code: Update modules, simplifying programming.
    !   01/03/2012      M.H.A. Piro         Call the GetNewAssemblage subroutine at the end of Leveling if
    !                                       the GetFirstAssemblage calculates the correct phase assemblage.
    !   02/17/2012      M.H.A. Piro         Check if allocatable arrays have been allocated and if so, have they
    !                                       changed in dimension.
    !   03/24/2021      S.Y. Kwon           Introducing dPhasePotential for checking global minimum.
    !   12/07/2021      S.Y. Kwon           Considering minimum point of each solution in Leveling
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to provide initial estimates of the equilibrium phase
    !! assemblage for the GEMSolver using the "Leveling" technique of Eriksson and Thompson.  The fundamental
    !! principle of the leveling technique is to temporarily assume that all species may be treated as pure
    !! stoichiometric phases.  In more mathematical terms, the mixing terms in the chemical potential
    !! function of a solution species becomes zero and the Gibbs energy of the system becomes a linear function.
    !! The Leveling subroutine can therefore be interpreted as a linear optimizer that provides initial
    !! estimates for the GEMSolver subroutine.
    !!
    !! There are several advantages in employing the Leveling technique:
    !! <ol>
    !! <li> Close estimates are provided for dElementPotential and dChemicalPotential, </li>
    !! <li> An estimate is provided for the phase assemblage, which is often fairly close to the equilibrium
    !!      assemblage (i.e., which phases are expected to be most stable), </li>
    !! <li> The estimated phase assemblage is determined with relatively few iterations, and </li>
    !! <li> Estimates are provided for the mole fractions of all constituents in all solution phases. </li>
    !! </ol>
    !!
    !
    ! References:
    ! ===========
    !
    !> \details For further information regarding the "Leveling" method, refer to the following literature:
    !!
    !! <ul>
    !!
    !! <li>     G. Eriksson and W.T. Thompson, "A Procedure to Estimate Equilibrium
    !!          Concentrations in Multicomponent Systems and Related Applications,"
    !!          CALPHAD, V. 13, N. 4, pp. 389-400 (1989).
    !! </li>
    !! <li>
    !!          M.H.A. Piro, "Computation of Thermodynamic Equilibria Pertinent to Nuclear Materials
    !!          in Multi-Physics Codes," PhD Dissertation, Royal Military College of Canada, 2011.
    !! </li>
    !! <li>
    !!          M.H.A. Piro, M.J. Welland, B.J. Lewis, W.T. Thompson and D.R. Olander, "Development
    !!          of a Self-Standing Numerical Tool to Compute Chemical Equilibria in Nuclear Materials,"
    !!          Top Fuel Conference, Paris, France (2009).
    !! </li>
    !! <li>
    !!          M.H.A. Piro and S. Simunovic, "Performance Enhancing Algorithms for Computing
    !!          Thermodynamic Equilibria," CALPHAD, 39 (2012) 104-110.
    !! </li>
    !! </ul>
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! INFOThermo            An integer scalar indicating a successful exit (0) or an error that is
    !                        encountered.
    ! nElements             The number of elements in the system.
    ! iAssemblage           Integer vector containing the indices of species (temporarily treated as pure
    !                        stoichiometric phases) estimated to be part of the equilibrium phase
    !                        assemblage.
    ! iSpeciesTotalAtoms    An integer vector representing the total number of atoms per formula mass of
    !                        a particular species.
    ! dMolesPhase           Number of moles of a phase at equilibrium
    ! dMolesElement         Total number of moles of an element (normalized for numerical purposes).
    ! dLevel                The adjustment applied to the element potentials as a result of leveling.
    ! dChemicalPotential    A double real vector representing the chemical potential of each species and
    !                        pure condensed phase.  This variable is numerical adjusted and represented
    !                        with respect to the element potentials. (per mole of atoms)
    ! dElementPotential     The chemical potentals of the pure elements.
    ! dTolerance            A double real vector representing numerical tolerances (defined in
    !                        InitThermo.f90).
    !
    !-------------------------------------------------------------------------------------------------------------
    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver
!
    implicit none
!
    integer :: iter, i, j, k
    integer,dimension(nElements)::          iAssemblageLast
    real(8),dimension(-1:nSolnPhasesSys)::  dTempVec
    real(8),dimension(nSpecies):: dTempVec1
!
!
    call InitGEMSolverNew
!
    ! Consider the minimum point of each solution phase (Eric Method)
    ! call CompInitMinSolnPoint
!
!
!
!
!
    ! START LEVELING:
    LOOP_Leveling: do iter = 1, 500
!
        !SY: Check global minimum: if all elements in phase potential are positive, the system is in global minimum
        dPhasePotential = dChemicalPotential - MATMUL(dAtomFractionSpecies,dElementPotential)
!
        ! print*,'dPhasePotential', MINVAL(dPhasePotential)
        if (MINVAL(dPhasePotential) > -1D-10) exit LOOP_Leveling
!
        ! Determine the next phase assemblage to be tested:
        call GetNewAssemblage(iter)
!
!
        ! Exit if an error has been encountered:
        if (INFOThermo /= 0) exit LOOP_Leveling
!
    end do LOOP_Leveling
!
!
    ! If the GetFirstAssemblage subroutine determined the correct phase assemblage, then the
    ! GetNewAssemblage subroutine was not called and the number of moles of each phase was not computed and
    ! the number of moles of each constituent was not computed.
    if ((iter == 1).AND.(INFOThermo == 0)) then
!
        ! Since the phase assemblage is comprised of only pure species, the number of moles of each phase
        ! can be easily computed (the stoichiometry matrix is a diagonal matrix) by the following:
        do i = 1, nElements
            dMolesPhase(i) = dMolesElement(i) / sum(dStoichSpeciesLevel(iAssemblage(i),:))
        end do
!
    end if
!
!
    ! If minimum point of a solution appears to be stable, distribute the dMolesPhase into soln components
    if (allocated(dMolesPhaseLast)) deallocate(dMolesPhaseLast)
    allocate(dMolesPhaseLast(nElements))
    iAssemblageLast         = iAssemblage
    dMolesPhaseLast         = dMolesPhase
    dTempVec                = 0D0
    dTempVec1               = 0D0
    nConPhases              = 0
    nSolnPhases             = 0
    iAssemblage             = 0
    dMolesPhase             = 0
!
!
!
    ! Calculate the total number of moles for each phase:
    ! All the species amount were added when more than 
    ! two species were predicted for a solution phase
    do i = 1, nElements
        j  = iPhaseLevel(iAssemblageLast(i))
        dTempVec(j) = dTempVec(j) + dMolesPhaseLast(i)
    end do
!
    ! Count the number of pure condensed phases and solution phases are assumed to be part of the phase
    ! assemblage and establish iAssemblage based on the results of Leveling:
!
    LOOP_AddPhase: do i = 1, nElements
!
        if ((iPhaseLevel(iAssemblageLast(i)) <= 0).AND.(nConPhases + nSolnPhases < nElements)) then
            nConPhases              = nConPhases + 1
            iAssemblage(nConPhases) = iAssemblageLast(i)
            dMolesPhase(nConPhases) = dMolesPhaseLast(i)
        elseif ((iPhaseLevel(iAssemblageLast(i)) > 0).AND.(nConPhases + nSolnPhases < nElements)) then
!
            do j = 1,nSolnPhases
                k = nElements - j + 1
                ! Ensure that this solution phase is not already stored:
                if (iAssemblage(k) == -iPhaseLevel(iAssemblageLast(i))) cycle LOOP_AddPhase
            end do
!
            nSolnPhases = nSolnPhases + 1
!
            j = nElements - nSolnPhases + 1
            k = iPhaseLevel(iAssemblageLast(i))
!
            iAssemblage(j) = -k
            dMolesPhase(j) = DMAX1(dTempVec(k),dTolerance(9))
            lSolnPhases(k) = .TRUE.
        end if
    end do LOOP_AddPhase
!
!
    !Revert back modified variables
    dTempVec1 = dPhasePotential(:nSpecies)
!
!
    deallocate(dAtomFractionSpecies,dChemicalPotential,dPhasePotential)
    allocate(dAtomFractionSpecies(nSpecies,nElements),&
    dChemicalPotential(nSpecies),dPhasePotential(nSpecies))
!
    dAtomFractionSpecies = dAtomFractionSpeciesOld
    dChemicalPotential   = dChemicalPotentialOld
    dPhasePotential      = dTempVec1
!
    return
!
end subroutine LevelingSolver
!
!
