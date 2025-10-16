!
subroutine InitGEMLineSearch(dStepLength)
!
    !---------------------------------------------------------------------------
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this subroutine is to initialize the line search algorithm.
    ! Specifically, the initial step length needs to be determined before the
    ! line search loop starts.
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer      Description of change
    !   ----            ----------      ---------------------
    !
    !   05/08/2012      M.H.A. Piro     Original Code
    !   07/04/2012      M.H.A. Piro     If the number of moles of any soln
    !                                    phase is to be significantly reduced
    !                                    and the system is not stagnant, then
    !                                    then further dampen the system.
    !   09/29/2012      M.H.A. Piro     Revert the system if the phase assemblage
    !                                    has not changed in 50 iterations and
    !                                    the maximum change to the system
    !                                    variables is extremely large.
    !   12/22/2021      S.Y Kwon        Revised the code to be compatible with
    !                                    a new minimization method i.e., separating
    !                                    gradient descent and phase assemblage check.
    !                                    There are stil many redundant codes carried
    !                                    over from previous codes. Require cleaning.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! dStepLength       A double real scalar representing the step length.
    ! dMaxIncrease      A double real scalar representing the maximum
    !                    increase in the number of moles of a solution phase.
    ! dMaxDecrease      A double real scalar representing the maximum
    !                    decrease in the number of moles of a solution phase.
    ! lCompEverything   A logical scalar indicating whether everything should
    !                    be computed in a particular subroutine (true) or
    !                    not (false).
    !
    !---------------------------------------------------------------------------
    USE ModuleThermo
    USE ModuleGEMSolver
!
    implicit none
!
    integer                       :: i, j, k, l, nMisciblePhases
    real(8)                       :: dStepLength, dTemp, dMaxIncrease, dMaxDecrease, dMaxChange, dMaxGamma
    logical                       :: lCompEverything
!x
!
    ! Initialize variables:
    if (allocated(dMolesSpeciesLast)) deallocate(dMolesSpeciesLast)
    if (allocated(dMolesPhaseLast)) deallocate(dMolesPhaseLast)
    if (allocated(dElementPotentialLast)) deallocate(dElementPotentialLast)
    allocate(dMolesSpeciesLast(nSpecies),dMolesPhaseLast(nElements),dElementPotentialLast(nElements))
    dMolesSpeciesLast =dMolesSpecies
    dMolesPhaseLast   =dMolesPhase
    dElementPotentialLast =dElementPotential
    lCompEverything  = .FALSE.
    dStepLength = 1D0
    nMisciblePhases  = 0
    dMaxGamma = 0D0
!
    ! Count the number of stable miscible phases:
    do j = 1, nSolnPhases
        k = -iAssemblage(nElements - j + 1)
        if (lMiscibility(k)) nMisciblePhases = nMisciblePhases + 1
    end do
!
    ! Update the number of moles of pure condensed phases:
    do i = 1, nConPhases
        j = nElements + nSolnPhases + i
        dMolesPhase(i) = dUpdateVar(j)
    end do
    
    ! Update the number of moles of solution phases:
    do l = 1, nSolnPhases
        k = -iAssemblage(nElements - l + 1)     ! Absolute solution phase index.
        
        do i = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
            dTemp = 0D0
            do j = 1, nElements
                dTemp = dTemp + dUpdateVar(j) * dStoichSpecies(i,j)
            end do
!
            dTemp = dTemp / DFLOAT(iParticlesPerMole(i))
!
            ! NOTE: The variable dMolFraction is used temporarily to represent
            ! the fractional update (i.e., x_i/y_i)to dMolesSpecies, but it does not replace
            ! dMolesSpecies in the event that further dampening is required. 
            dMolFraction(i) = 1D0 + dUpdateVar(nElements + l) &
            +(dTemp - dChemicalPotential(i))

            !Start of the lambda correciton----------------------------------------------------
            !where dTemp stands for lambda at each iteration
            if (dMolFraction(i) /= 1D0) dTemp = 1D0 / (1D0 - dMolFraction(i))
!
            if ((dMolFraction(i) < 0D0).AND.(dTemp < dStepLength)) then
                ! if (dMolFraction(i) < 0D0) print*, 'Negative species', dMolFraction(i), iterGlobal, dGEMFunctionNorm, dMolesPhase

                if (dMolesSpecies(i) < 1D-50) then
                !     ! This species is very small and the system is trying to make it negative.
                !     ! Reduce the mass of this species by an arbitrary positive value less than
                !     ! unity:
                !     dMolFraction(i) = 1D-3
                ! else
                    ! If the species is being driven negative particularly aggressively,
                    ! then reduce the step size by a larger margin to avoid float troubles.
                    dStepLength = dTemp * 0.99D0
                end if
            end if

            
!             ! End of the lambda correciton----------------------------------------------------
            dMolesSpecies(i) = dMolesSpecies(i) * dMolFraction(i)
!
        end do
    end do
!
	! Update the element potentials:
    do i = 1, nElements
        dElementPotential(i) = dUpdateVar(i)
    end do

    ! Update the system variables:
    call UpdateSystemVariables(dStepLength)

    ! Compute the chemical potentials of solution species:
    call CompChemicalPotential(lCompEverything)
!
    ! Compute the functional norm:
    call CompFunctionNorm
!
end subroutine InitGEMLineSearch
!
!
