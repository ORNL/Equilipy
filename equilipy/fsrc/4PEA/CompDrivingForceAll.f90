!
!
subroutine CompDrivingForceAll
!
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CompDrivingForceAll.f90
    !> \brief   Compute the driving force of all phases.
    !> \author  S.Y. Kwon
    !> \date    Dec. 20, 2021
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   12/20/2021      S.Y. Kwon         Original code
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to compute the driving force of all pure condensed phases in
    !! the database.  The driving force is defined as the difference between the standard molar Gibbs energy of
    !! a pure condensed phase and the corresponding value computed from the element potentials.  This value is
    !! used to determine whether a pure condensed phase should be added to the system.  For a more thorough
    !! explanation of the chemical significance of the driving force, refer to the following literature:
    !!
    !! H.L. Lukas, S.G. Fries, B. Sundman, Computational Thermodynamics - The Calphad Method, Cambridge
    !! University Press, New York, 2007.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[out]  iMaxDrivingForce    An integer scalar representing the index of the pure condensed phase
    !!                                   with the maximum driving force.  A value of zero is returned by default.
    !> \param[out]  dMaxDrivingForce    A double real scalar representing the maximum driving force of all pure
    !!                                   condensed phases.  A value of zero is returned by default.
    !
    !-------------------------------------------------------------------------------------------------------------
    USE ModuleThermo
    USE ModuleGEMSolver
!
    implicit none
!
    integer :: i, j, idummy(1), m, n
    real(8) :: dTemp
    logical :: ldummy
!
    dPhasePotential = 0d0
!
    m = nSpeciesPhase(nSolnPhasesSys) + 1
    n = nSpecies - nDummySpecies
!
    ! Calculate driving force of all stoichiometric phases:
    do i = m, n
!
        ! Compute the chemical potential of this phase as defined by the element potentials:
        dTemp = 0D0
        do j = 1, nElements
            dTemp = dTemp + dElementPotential(j) * dStoichSpecies(i,j)
        end do
!
!
        ! Compute the driving force:
        dTemp = dStdGibbsEnergy(i) - dTemp
!
        ! Normalize per gram-atom:
        dTemp = dTemp / dSpeciesTotalAtoms(i)
!
        ! !SY: Update phase potential of condensted phase
        dPhasePotential(i) = dTemp
!
!
    end do
!
    !Store the most stable stoichiometric phase to be added
    idummy                 = MINLOC(dPhasePotential)
    iMinDrivingForceStoich = idummy(1)
    dMinDrivingForceStoich = dPhasePotential(iMinDrivingForceStoich)
!
    ! Calculate driving force of all solution phases:
    LOOP_SolnPhaseSys: do i = 1, nSolnPhasesSys
!
        ! Skip this phase if it is already predicted to be stable:
        if (lSolnPhases(i)) cycle LOOP_SolnPhaseSys
!
!
        ! Compute the mole fractions of all constituents in this solution phase:
        call CompMolFraction(i)
!
        ! Skip this phase if it is not the first "phase" in a phase with a miscibility gap:
        if (lMiscibility(i)) then
            call CheckMiscibilityGap(i,ldummy)
            ! print*,'miscibility',i,ldummy
            cycle LOOP_SolnPhaseSys
        end if
!
    end do LOOP_SolnPhaseSys
!
    !Store the most stable solution phase to be added
    idummy                 = MINLOC(dPhasePotential(:m))
    iMinDrivingForceSoln   = idummy(1)
    dMinDrivingForceSoln   = dPhasePotential(iMinDrivingForceSoln)
!
    return
!
end subroutine CompDrivingForceAll
!
!
!
