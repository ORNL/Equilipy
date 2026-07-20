!
!
subroutine GEMLineSearch
!
    !---------------------------------------------------------------------------
    !
    !> \file    GEMLineSearch.f90
    !> \brief   Perform a line search for the GEMSolver.f90
    !> \author  M.H.A. Piro
    !> \date    Apr. 25, 2012
    !> \sa      GEMSolver.f90
    !> \sa      GEMNewton.f90
    !> \sa      CompChemicalPotential.f90
    !> \sa      CompFunctionNorm.f90
    !> \todo    Figure out a long term solution for the condition to allow the solution
    !!           to only perform 1 iteration if the functional norm is less than 1E-6.
    !! \todo    Figure out a more permanent solution to how dMaxGamma shoudl be handled.
    !!           I think that I need to get away from that and focus on ensuring that the
    !!           Wolfe conditions have been satisfied.
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   04/25/2012      M.H.A. Piro         Original code
    !   05/02/2012      M.H.A. Piro         Improved calculation of initial step
    !                                       length (constrain changes to the Element
    !                                       potentials).
    !   05/08/2012      M.H.A. Piro         Improved calculation of initial step length
    !                                       (constrain the maximum change to the total
    !                                       number of moles of a solutoin phase).
    !   04/11/2013      M.H.A. Piro         When computing a steplength that constrains
    !                                       the maximum change to the element potential
    !                                       to less than or equal to unity, exclude elements
    !                                       with zero moles.  In other words, exclude
    !                                       electrons corresponding to ionic phases
    !                                       that are not currently stable.  This is also
    !                                       done when updating the element potentials.
    !   05/08/2013      M.H.A. Piro         Changed tolerance for calculating the steplength
    !                                       when the number of moles of a species tends to zero
    !                                       to 1D-50 from 1D-100.
    !   06/08/2013      M.H.A. Piro         In determining an initial steplength, do not constrain
    !                                       the maximum decrease in the number of moles of a
    !                                       solution phase by a certain increment (e.g., 50%) if
    !                                       the number of moles of that phase is below a certain
    !                                       value (e.g., 10**(-9)).  The motivation for doing this
    !                                       is that a solution phase may be driving out of the system
    !                                       but it may be inhibited if this condition is not made.
    !   06/08/2013      M.H.A. Piro         Exit the Wolfe loop if the functional norm is below a
    !                                       certain tolerance (e.g., 10**(-6)).
    !   04/02/2014      M.H.A. Piro         I changed one of the conditions to satisfy the line search.
    !                                       Specifically, if the relative change of the functional
    !                                       norm is 0.95 < F_norm < 1.0 to 0.97 < F_norm < 1.0.
    !                                       Consider the following scenario: the functional norm
    !                                       is 1 at global iteration 5 and currently we are at iter 6.
    !                                       The f-norm of the first line search loop is 6.1 and then
    !                                       its 5.86 (dTemp = 0.961).  Before, the system would then
    !                                       exit, but the f-norm on the global scale is not converging.
    !   08/20/2015      M.H.A. Piro         Previously, the maximum change in the element potentials was
    !                                       constrained to 1.  Now, it is a linear function that varies
    !                                       with respect to the functional norm.  The motivation for this
    !                                       is that an initially poor guess may be far from equilibrium,
    !                                       and laxing this constraint accelerates convergence.  Another
    !                                       change is the way that dStepLength is initialized when
    !                                       dMolesSpecies tends to zero.  Previously, this would compute
    !                                       the minimum dStepLength that corresponds to reducing dMolesSpecies
    !                                       by 100.  Now, this is NOT applied when dMolesSpecies is incredibly
    !                                       small (e.g., less than numerical tolerance).
    !   12/18/2018      M.H.A. Piro         I effectively removed the constraint applied to the maximum
    !                                       change to the element potentials. I think that this is inefficient
    !                                       because it's probably best to leave it to ensuring that the
    !                                       Wolfe conditions have been satisfied. Furthermore, it seems to
    !                                       be an issue for SUBG phases, which are effectively on a knife edge.
    !   07/20/2026      S.Y. Kwon           Added model-aware Gibbs and merit line searches with diagnostics for active solution phases.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to perform a line search using
    !! the direction vector computed by the Newton/Broyden solver.  The system
    !! is updated using an appropriate step length that satisfies the Wolfe
    !! conditions.  Specifically, values of dChemicalPotential and dMolesPhase
    !! are updated.  It is possible for the system of equations to be ill-
    !! behaved and yield inappropriate results.  An initial step-length is
    !! computed by normalizing the largest change of the system variables by a
    !! pre-defined value.  The maximum change to the element potentials is 1 and
    !! the maximum change to the number of moles of a solution phase is twice of
    !! the previous value.  For more information, refer to Chapter 6 of the
    !! above reference.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! nElements             The number of elements in the system.
    ! nConPhases            The number of pure condensed phases in the system.
    ! nSolnPhases           The number of solution phases in the system.
    ! dUpdateVar            A double real vector that contains updates to the
    !                        system variables (element potentials, moles of
    !                        solution phases and moles of pure condensed phases).
    ! dStepLength           Step length applied to the direction vector
    !                        (dUpdateVar)
    ! dMolesPhase           A double real vector representing the number of
    !                        moles of phases predicted to be stable.
    ! dLevel                The adjustment applied to the chemical potentials of
    !                        the elements
    ! iPhaseDampen          Integer vector that counts the number of times that
    !                        the number of moles
    !                        of a solution phase had to be dampened.
    ! dChemicalPotential    A double real vector representing the chemical
    !                        potential of each species and pure condensed phase.
    ! dGEMFunctionNorm      A double real scalar representing the norm of the
    !                        functional vector in the PGESolver.
    ! dGEMFunctionNormLast  A double real scalar representing the norm of the
    !                        functional vector from the previous iteration.
    !
    !---------------------------------------------------------------------------
!
    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver
!
    implicit none
!
    integer                       :: iterWolfe, iWolfeCount, iSolnPhase
    real(8)                       :: dStepLength
    real(8)                       :: dBestFunctionNorm, dBestStepLength
    real(8)                       :: dLineSearchNormTolerance
    real(8)                       :: dTrialGibbs, dTrialMerit, dBestGibbs, dBestMerit
    real(8), dimension(:), allocatable :: dBestMolesSpecies, dBestMolesPhase, dBestElementPotential
    logical                       :: lCompEverything, lPairExchangeAccepted, lSlopeRecorded
!
!
    ! Initialize variable:
    dGEMFunctionNormLast    = dGEMFunctionNorm
    lCompEverything         = .FALSE.
    iWolfeCount             = 0
    iGEMLineSearchIterationCount = 0
    iGEMLineSearchNegativeFactorCount = 0
    iGEMLineSearchFloorCount = 0
    iGEMLineSearchNoDescent = 0
    iGEMLineSearchNoDescentClass = 0
    dGEMLineSearchInitialNorm = dGEMFunctionNorm
    dGEMLineSearchBestNorm = dGEMFunctionNorm
    dGEMLineSearchFinalNorm = dGEMFunctionNorm
    dGEMLineSearchInitialGibbs = CurrentActiveGibbs()
    dGEMLineSearchInitialMerit = ActiveMassMerit(dGEMLineSearchInitialGibbs)
    dGEMLineSearchBestGibbs = dGEMLineSearchInitialGibbs
    dGEMLineSearchFinalGibbs = dGEMLineSearchInitialGibbs
    dGEMLineSearchBestMerit = dGEMLineSearchInitialMerit
    dGEMLineSearchFinalMerit = dGEMLineSearchInitialMerit
    dGEMNewtonDirNormSlope = 0D0
    dGEMNewtonDirGibbsSlope = 0D0
    dGEMNewtonDirMeritSlope = 0D0
    dGEMLineSearchInitialStep = 0D0
    dGEMLineSearchBestStep = 0D0
    dGEMLineSearchFinalStep = 0D0
    dLineSearchNormTolerance = 1D-5
    if (lPostProcess) dLineSearchNormTolerance = 1D-10
    do iSolnPhase = 1, nSolnPhases
        if ((-iAssemblage(nElements - iSolnPhase + 1) > 0).AND.&
            (-iAssemblage(nElements - iSolnPhase + 1) <= SIZE(cSolnPhaseType))) then
            if (cSolnPhaseType(-iAssemblage(nElements - iSolnPhase + 1)) == 'SUBOM') then
                dLineSearchNormTolerance = MIN(dLineSearchNormTolerance, 1D-9)
            end if
        end if
    end do
!

    ! Initialize the line search method:
    call InitGEMLineSearch(dStepLength)
    dGEMLineSearchInitialStep = dStepLength
    dGEMLineSearchFinalStep = dStepLength
    allocate(dBestMolesSpecies(nSpecies), dBestMolesPhase(nElements), dBestElementPotential(nElements))
    dBestFunctionNorm   = dGEMFunctionNormLast
    dBestStepLength     = 0D0
    dBestGibbs          = dGEMLineSearchInitialGibbs
    dBestMerit          = dGEMLineSearchInitialMerit
    lSlopeRecorded      = .FALSE.
    dBestMolesSpecies   = dMolesSpeciesLast
    dBestMolesPhase     = dMolesPhaseLast
    dBestElementPotential = dElementPotentialLast
    dTrialGibbs = CurrentActiveGibbs()
    dTrialMerit = ActiveMassMerit(dTrialGibbs)
    if (dStepLength > 0D0) then
        call RecordDirectionalSlope(dStepLength, dTrialGibbs, dTrialMerit)
        lSlopeRecorded = .TRUE.
    end if
    if (dTrialGibbs < dBestGibbs) then
        dBestGibbs = dTrialGibbs
        dBestMerit = dTrialMerit
    end if
    if (dGEMFunctionNorm < dBestFunctionNorm) then
        dBestFunctionNorm     = dGEMFunctionNorm
        dBestStepLength       = dStepLength
        dBestMolesSpecies     = dMolesSpecies
        dBestMolesPhase       = dMolesPhase
        dBestElementPotential = dElementPotential
    end if
    
    
!     ! Commence line search:
    LOOP_WOLFE: do iterWolfe = 1, 8
        iWolfeCount = iterWolfe
        iGEMLineSearchIterationCount = iWolfeCount
!
        ! If the functional norm is already small, call it a day:
        if (dGEMFunctionNorm < dLineSearchNormTolerance) exit LOOP_WOLFE
!
        ! Check if the system has sufficiently progressed, otherwise dampen:
        if (dGEMFunctionNorm < 0.999D0*dGEMFunctionNormLast) then
!           
            ! The system has sufficiently progressed, exit:
            exit LOOP_WOLFE
            
        else
            ! Dampen the system variables:
            dStepLength = dStepLength*0.5 ! 0.8 makes dStepLength vary from 0.01 to 1
            dGEMLineSearchFinalStep = dStepLength
            call UpdateSystemVariables(dStepLength)
!
            ! Compute the chemical potentials of solution species:
            call CompChemicalPotential(lCompEverything)
            ! Compute the functional norm:
            call CompFunctionNorm
            dTrialGibbs = CurrentActiveGibbs()
            dTrialMerit = ActiveMassMerit(dTrialGibbs)
            if ((.NOT.lSlopeRecorded).AND.(dStepLength > 0D0)) then
                call RecordDirectionalSlope(dStepLength, dTrialGibbs, dTrialMerit)
                lSlopeRecorded = .TRUE.
            end if
            if (dTrialGibbs < dBestGibbs) then
                dBestGibbs = dTrialGibbs
                dBestMerit = dTrialMerit
            end if
            if (dGEMFunctionNorm < dBestFunctionNorm) then
                dBestFunctionNorm     = dGEMFunctionNorm
                dBestStepLength       = dStepLength
                dBestMolesSpecies     = dMolesSpecies
                dBestMolesPhase       = dMolesPhase
                dBestElementPotential = dElementPotential
            end if

            if (dStepLength < 1D-12) exit LOOP_WOLFE
            
!
            ! Reiterate:
            cycle LOOP_WOLFE
!
        end if
! !
        end do LOOP_WOLFE

        if (dBestFunctionNorm < dGEMFunctionNormLast) then
            dMolesSpecies    = dBestMolesSpecies
            dMolesPhase      = dBestMolesPhase
            dElementPotential = dBestElementPotential
            dGEMLineSearchFinalStep = dBestStepLength
        else
            dMolesSpecies    = dMolesSpeciesLast
            dMolesPhase      = dMolesPhaseLast
            dElementPotential = dElementPotentialLast
            dGEMLineSearchFinalStep = 0D0
            iGEMLineSearchNoDescent = 1
        end if
        dGEMLineSearchBestNorm = dBestFunctionNorm
        dGEMLineSearchBestStep = dBestStepLength
        dGEMLineSearchBestGibbs = dBestGibbs
        dGEMLineSearchBestMerit = dBestMerit
        if (dBestFunctionNorm >= dGEMLineSearchInitialNorm) iGEMLineSearchNoDescent = 1
        
    call CompChemicalPotential(lCompEverything)
    call CompFunctionNorm

    if (iGEMLineSearchNoDescent > 0) then
        call TryBinarySUBGQPairExchange(lPairExchangeAccepted)
        if (lPairExchangeAccepted) then
            iGEMLineSearchNoDescent = 0
            iGEMLineSearchNoDescentClass = 0
            dGEMLineSearchBestNorm = dGEMFunctionNorm
            dGEMLineSearchBestStep = 0D0
            dGEMLineSearchFinalStep = 0D0
            dGEMLineSearchBestGibbs = CurrentActiveGibbs()
            dGEMLineSearchBestMerit = ActiveMassMerit(dGEMLineSearchBestGibbs)
        end if
    end if
    dGEMLineSearchFinalNorm = dGEMFunctionNorm
    dGEMLineSearchFinalGibbs = CurrentActiveGibbs()
    dGEMLineSearchFinalMerit = ActiveMassMerit(dGEMLineSearchFinalGibbs)
    call ClassifyNoDescent
    ! print*,'iterWolfe',iterWolfe
!
    if (allocated(dBestMolesSpecies)) deallocate(dBestMolesSpecies)
    if (allocated(dBestMolesPhase)) deallocate(dBestMolesPhase)
    if (allocated(dBestElementPotential)) deallocate(dBestElementPotential)
    

    return
!
contains

    real(8) function CurrentActiveGibbs()

        implicit none

        integer :: iSlot, iSolnPhase, iSpecies

        CurrentActiveGibbs = 0D0
        do iSlot = 1, nConPhases
            iSpecies = iAssemblage(iSlot)
            if (iSpecies > 0) then
                CurrentActiveGibbs = CurrentActiveGibbs + dMolesPhase(iSlot) * dStdGibbsEnergy(iSpecies)
            end if
        end do

        do iSlot = 1, nSolnPhases
            iSolnPhase = -iAssemblage(nElements - iSlot + 1)
            if ((iSolnPhase > 0).AND.(iSolnPhase <= nSolnPhasesSys)) then
                CurrentActiveGibbs = CurrentActiveGibbs + dGibbsSolnPhase(iSolnPhase)
            end if
        end do

        return

    end function CurrentActiveGibbs

    real(8) function ActiveMassMerit(dGibbsIn)

        implicit none

        real(8), intent(in) :: dGibbsIn

        ActiveMassMerit = dGibbsIn + 0.5D0 * dGEMMassBalanceNorm * dGEMMassBalanceNorm

        return

    end function ActiveMassMerit

    subroutine RecordDirectionalSlope(dStepIn, dTrialGibbsIn, dTrialMeritIn)

        implicit none

        real(8), intent(in) :: dStepIn, dTrialGibbsIn, dTrialMeritIn

        dGEMNewtonDirNormSlope = (dGEMFunctionNorm - dGEMLineSearchInitialNorm) / dStepIn
        dGEMNewtonDirGibbsSlope = (dTrialGibbsIn - dGEMLineSearchInitialGibbs) / dStepIn
        dGEMNewtonDirMeritSlope = (dTrialMeritIn - dGEMLineSearchInitialMerit) / dStepIn

        return

    end subroutine RecordDirectionalSlope

    subroutine ClassifyNoDescent

        implicit none

        real(8) :: dGibbsTol

        if (iGEMLineSearchNoDescent <= 0) return

        dGibbsTol = 1D-12 * DMAX1(1D0,DABS(dGEMLineSearchInitialGibbs))
        if (dGEMLineSearchBestGibbs < dGEMLineSearchInitialGibbs - dGibbsTol) then
            iGEMLineSearchNoDescentClass = 1
        else
            iGEMLineSearchNoDescentClass = 2
        end if

        return

    end subroutine ClassifyNoDescent

end subroutine GEMLineSearch

!
!
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
!
!
!
!
!
!
!
    !---------------------------------------------------------------------------
    !                            END - GEMLineSearch.f90
    !---------------------------------------------------------------------------
