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
    !   Date            Programmer      Description of change
    !   ----            ----------      ---------------------
    !   04/25/2012      M.H.A. Piro     Original code
    !   05/02/2012      M.H.A. Piro     Improved calculation of initial step
    !                                    length (constrain changes to the Element
    !                                    potentials).
    !   05/08/2012      M.H.A. Piro     Improved calculation of initial step length
    !                                    (constrain the maximum change to the total
    !                                    number of moles of a solutoin phase).
    !   04/11/2013      M.H.A. Piro     When computing a steplength that constrains
    !                                    the maximum change to the element potential
    !                                    to less than or equal to unity, exclude elements
    !                                    with zero moles.  In other words, exclude
    !                                    electrons corresponding to ionic phases
    !                                    that are not currently stable.  This is also
    !                                    done when updating the element potentials.
    !   05/08/2013      M.H.A. Piro     Changed tolerance for calculating the steplength
    !                                    when the number of moles of a species tends to zero
    !                                    to 1D-50 from 1D-100.
    !   06/08/2013      M.H.A. Piro     In determining an initial steplength, do not constrain
    !                                    the maximum decrease in the number of moles of a
    !                                    solution phase by a certain increment (e.g., 50%) if
    !                                    the number of moles of that phase is below a certain
    !                                    value (e.g., 10**(-9)).  The motivation for doing this
    !                                    is that a solution phase may be driving out of the system
    !                                    but it may be inhibited if this condition is not made.
    !   06/08/2013      M.H.A. Piro     Exit the Wolfe loop if the functional norm is below a
    !                                    certain tolerance (e.g., 10**(-6)).
    !   04/02/2014      M.H.A. Piro     I changed one of the conditions to satisfy the line search.
    !                                    Specifically, if the relative change of the functional
    !                                    norm is 0.95 < F_norm < 1.0 to 0.97 < F_norm < 1.0.
    !                                    Consider the following scenario: the functional norm
    !                                    is 1 at global iteration 5 and currently we are at iter 6.
    !                                    The f-norm of the first line search loop is 6.1 and then
    !                                    its 5.86 (dTemp = 0.961).  Before, the system would then
    !                                    exit, but the f-norm on the global scale is not converging.
    !   08/20/2015      M.H.A. Piro     Previously, the maximum change in the element potentials was
    !                                    constrained to 1.  Now, it is a linear function that varies
    !                                    with respect to the functional norm.  The motivation for this
    !                                    is that an initially poor guess may be far from equilibrium,
    !                                    and laxing this constraint accelerates convergence.  Another
    !                                    change is the way that dStepLength is initialized when
    !                                    dMolesSpecies tends to zero.  Previously, this would compute
    !                                    the minimum dStepLength that corresponds to reducing dMolesSpecies
    !                                    by 100.  Now, this is NOT applied when dMolesSpecies is incredibly
    !                                    small (e.g., less than numerical tolerance).
    !   12/18/2018      M.H.A. Piro     I effectively removed the constraint applied to the maximum
    !                                    change to the element potentials. I think that this is inefficient
    !                                    because it's probably best to leave it to ensuring that the
    !                                    Wolfe conditions have been satisfied. Furthermore, it seems to
    !                                    be an issue for SUBG phases, which are effectively on a knife edge.
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
    integer                       :: iterWolfe
    real(8)                       :: dStepLength, dTemp, dWolfeFunctionNormLast
    logical                       :: lCompEverything
!
!
    ! Initialize variable:
    dWolfeFunctionNormLast  = 1D-12
    dGEMFunctionNormLast    = dGEMFunctionNorm
    lCompEverything         = .FALSE.
!

    ! Initialize the line search method:
    call InitGEMLineSearch(dStepLength)
    
    
!     ! Commence line search:
    LOOP_WOLFE: do iterWolfe = 1, 20
!
        ! Compute the fractional change in the functional norm:
        dTemp = MIN1(dGEMFunctionNorm,dGEMFunctionNormLast) / dWolfeFunctionNormLast
!
        ! If the functional norm is already small, call it a day:
        if (dGEMFunctionNorm < 1D-7) exit LOOP_WOLFE
!
        ! Check if the system is diverging or has sufficiently progressed, otherwise dampen:
        if (dGEMFunctionNorm < 0.999D0*dGEMFunctionNormLast) then
!           
            ! The system has sufficiently progressed, exit:
            exit LOOP_WOLFE
            
        elseif ((iterWolfe > 1).AND.(dTemp >= 0.9D0).AND.(dTemp <= 1D0)) then
            ! The system has sufficiently progressed, exit:
            exit LOOP_WOLFE

        else
            ! Dampen the system variables:
            ! If the functional norm has only increased by a nominal amount (i.e., 1%), then exit:
            dTemp = MIN1(dGEMFunctionNorm,dGEMFunctionNormLast) / dWolfeFunctionNormLast
            if ((dTemp > 1D0).AND.(dTemp < 1.01D0)) exit LOOP_WOLFE
!
            dStepLength = dStepLength*0.8 ! 0.8 makes dStepLength vary from 0.01 to 1
            call UpdateSystemVariables(dStepLength)
!
            ! Compute the chemical potentials of solution species:
            call CompChemicalPotential(lCompEverything)
!
            dWolfeFunctionNormLast = dGEMFunctionNorm
!
            ! Compute the functional norm:
            call CompFunctionNorm

            if(dGEMFunctionNorm >= dWolfeFunctionNormLast .and.(iterWolfe > 2)) then
                ! dStepLength = 2

                ! call UpdateSystemVariables(dStepLength)

                exit LOOP_WOLFE
                
            end if
            
!
            ! Reiterate:
            cycle LOOP_WOLFE
!
        end if
! !
    end do LOOP_WOLFE
    
    call CompChemicalPotential(lCompEverything)
    call CompFunctionNorm
    ! print*,'iterWolfe',iterWolfe
!

    

    return
!
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
