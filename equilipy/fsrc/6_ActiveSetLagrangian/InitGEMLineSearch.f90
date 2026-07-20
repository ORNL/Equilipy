!> \brief Initialize Lagrangian GEM line search from the Newton update.
!!
!! \details Applies the raw Newton proposal to phase amounts, species amounts,
!! and elemental potentials, then computes the first lambda-damped line-search
!! state.  Raw full-step phase moles are stored before damping so active-set
!! diagnostics can see when a phase wants to leave the assemblage.
!!
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    InitGEMLineSearch.f90
    !> \brief   Initialize the line search for Lagrangian GEM.
    !> \author  M.H.A. Piro
    !> \date    May 8, 2012
    !> \sa      GEMLineSearch.f90
    !> \sa      GEMNewton.f90
    !> \sa      UpdateSystemVariables.f90
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   05/08/2012      M.H.A. Piro         Original code
    !   07/04/2012      M.H.A. Piro         Dampened large solution phase decreases
    !   09/29/2012      M.H.A. Piro         Added stale-assemblage revert logic
    !   12/22/2021      S.Y. Kwon           Split gradient descent from phase assemblage checks
    !   07/20/2026      S.Y. Kwon           Initialized line-search species directions from analytical solution curvature when available.
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to initialize the Lagrangian
    !! GEM line search.  It saves the previous active-set state, applies the
    !! full Newton proposal, computes the lambda damping needed to maintain
    !! positive species amounts, and evaluates the first damped trial state.
    !
    ! Required input variables:
    ! =========================
    !
    ! dUpdateVar          Newton unknown vector from GEMNewton.
    ! dMolesSpecies       Current active-set species moles before applying the update.
    ! dMolesPhase         Current active-set phase amounts before applying the update.
    ! dElementPotential   Current elemental potentials before applying the update.
    !
    !
    ! Output/updated variables:
    ! =========================
    !
    !> \param[out]  dStepLength   Initial lambda-damped line-search step.
    !
    ! dMolesSpeciesLast                 Previous species moles.
    ! dMolesPhaseLast                   Previous phase amounts.
    ! dElementPotentialLast             Previous elemental potentials.
    ! dGEMLSRawPhaseMoles       Raw full-step phase amounts before lambda damping.
    ! iGEMLineSearchNegativeFactorCount Count of raw negative species update factors.
    ! iGEMLineSearchNegativePhaseCount  Count of raw negative full-step phase amounts.
    !
    !
    ! Called subroutines/functions:
    ! =============================
    !
    ! UpdateSystemVariables   Applies the damped line-search step.
    ! CompChemicalPotential   Recomputes species potentials for the damped state.
    ! CompFunctionNorm        Recomputes the combined and split GEM residual norms.
    !
    !
    ! Primary callers:
    ! ================
    !
    ! GEMLineSearch   Starts each Lagrangian line search from the raw Newton proposal.
    !
    !
    ! Numerical assumptions:
    ! ======================
    !
    ! - Raw full-step phase amounts are diagnostics only.  Positivity is still
    !   enforced by the existing lambda damping and mole floor.
    ! - A negative raw phase amount is evidence for active-set review, not a
    !   direct phase removal inside this routine.
    !
    !-------------------------------------------------------------------------------------------------------------



subroutine InitGEMLineSearch(dStepLength)
    USE ModuleThermo
    USE ModuleGEMSolver
!
    implicit none
!
    integer                       :: i, j, k, l, nMisciblePhases
    real(8)                       :: dStepLength, dTemp, dMaxIncrease, dMaxDecrease, dMaxChange, dMaxGamma
    logical                       :: lCompEverything
    logical                       :: lUseAnalyticalSpeciesDirection
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
    iGEMLineSearchNegativePhaseCount = 0
    dGEMLineSearchMinRawPhaseMoles = 0D0
    dGEMLineSearchMinFinalPhaseMoles = 0D0
    if (allocated(dGEMLSRawPhaseMoles)) then
        dGEMLSRawPhaseMoles = dMolesPhaseLast
    end if
    if (allocated(dGEMLSFinalPhaseMoles)) then
        dGEMLSFinalPhaseMoles = dMolesPhaseLast
    end if
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
        if (allocated(dGEMLSRawPhaseMoles)) then
            dGEMLSRawPhaseMoles(i) = dMolesPhase(i)
        end if
    end do
    
    ! Update the number of moles of solution phases:
    do l = 1, nSolnPhases
        k = -iAssemblage(nElements - l + 1)     ! Absolute solution phase index.
        lUseAnalyticalSpeciesDirection = .FALSE.
        if (allocated(lGEMAnalyticalSpeciesDirection)) then
            if ((k > 0).AND.(k <= SIZE(lGEMAnalyticalSpeciesDirection))) then
                lUseAnalyticalSpeciesDirection = lGEMAnalyticalSpeciesDirection(k)
            end if
        end if
        
        do i = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
            if (lUseAnalyticalSpeciesDirection) then
                dMolFraction(i) = 1D0 + dGEMAnalyticalSpeciesDirection(i) / DMAX1(dMolesSpecies(i), dTolerance(8))
            else
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
            end if

            if (dMolFraction(i) < 0D0) then
                iGEMLineSearchNegativeFactorCount = iGEMLineSearchNegativeFactorCount + 1
            end if

            if (dMolesSpecies(i) * dMolFraction(i) < dTolerance(8)) then
                iGEMLineSearchFloorCount = iGEMLineSearchFloorCount + 1
            end if

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

        if (allocated(dGEMLSRawPhaseMoles)) then
            dGEMLSRawPhaseMoles(nElements - l + 1) = &
                SUM(dMolesSpecies(nSpeciesPhase(k-1) + 1:nSpeciesPhase(k)))
        end if
    end do

    if (allocated(dGEMLSRawPhaseMoles)) then
        dGEMLineSearchMinRawPhaseMoles = MINVAL(dGEMLSRawPhaseMoles)
        iGEMLineSearchNegativePhaseCount = COUNT(dGEMLSRawPhaseMoles < 0D0)
    end if
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
