!
!
subroutine Subminimization(iSolnPhaseIndex,lPhasePass)
!
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    Subminimization.f90
    !> \brief   Determine whether a particular non-ideal solution phase should be added to the system
    !!           by performing a subminimization routine.
    !> \author  M.H.A. Piro
    !> \date    Aug. 21, 2012
    !> \sa      CheckSolnPhaseAdd.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   08/21/2012      M.H.A. Piro         Original code
    !   08/30/2012      M.H.A. Piro         Enforce the Wolfe conditions in the line search algorithm.  Also,
    !                                        a check is performed to ensure that the solution phase is not the
    !                                        same as the other "solution phase."
    !   08/31/2012      M.H.A. Piro         Enforce a minimum mole fraction to prevent the system from
    !                                        continuing to push the mole fraction of a particular constituent
    !                                        to zero.  The tolerance is defined as the quotient of machine
    !                                        precision and the relative error tolerance of the mass balance
    !                                        (i.e., 10**(-15) / 10**(-5)).
    !   09/03/2012      M.H.A. Piro         Revert the last change...do not test for a minimum mole fraction.
    !   09/19/2012      M.H.A. Piro         Replaced the DGESV solver with the ArrowSolver.  This solver
    !                                        exploits the fact that the Hessian is a symmetric arrow matrix.
    !   09/24/2012      M.H.A. Piro         Refine convergence criteria: if the driving force doesn't change
    !                                        by a predefined amount (e.g., 1%).
    !   10/22/2012      M.H.A. Piro         Set a minimum constraint for dMolFraction of 1D-10 when initializing
    !                                        Subminimization.  This prevents extremely small values from being
    !                                        considered.
    !   02/28/2013      M.H.A. Piro         The routine was upgraded to allow for ionic phases (i.e., charge
    !                                        neutrality constraints must be considered).
    !   03/05/2013      M.H.A. Piro         Fix bug in SubminNewton in constructing the Hessian matrix: the
    !                                        stoichiometry coefficients of ionic species should be multiplied
    !                                        by (-1).
    !   03/12/2013      M.H.A. Piro         Compute the functional norm (currently computed by the residual
    !                                        of the mass fractions and charge neutrality constraints) and
    !                                        only proceed if this is below tolerance.
    !   03/14/2013      M.H.A. Piro         Constrain the minimum mole fraction to an arbitrarily small value
    !                                        (e.g., 1D-100) to avoid an floating point underflow error.
    !   06/04/2013      M.H.A. Piro         Lax tolerance for dSubMinFunctionNorm.
    !   11/26/2021      S.Y. Kwon           Lagrangian was based on Joule per mole of species. This method has
    !                                       an issue of calculating correct driving force and the corresponding
    !                                       composition of the minimum point. I modified Lagrangian so that it 
    !                                       is based on Joule per mole of atoms.
    !                                       
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to determine whether a particular non-ideal solution phase
    !! should be added to the system by performing a subminimization.  The criteria for adding any type of phase
    !! (regardless of whether it is a pure stoichiometric phase, ideal or non-ideal solution phase) is based on
    !! whether the driving force is positive (it should not be added) or negative (it should be added to the
    !! system).  The driving force is defined as the difference between the molar Gibbs energy of a particular
    !! phase and the corresponding value defined by the element potentials.  A graphical interpretation of the
    !! driving force is the difference between the molar Gibbs energy of a particular phase (represented by a
    !! line for a stoichiometric phase or a point on a curve for a solution phase) and the corresponding value
    !! on the Gibbs hyperplane.
    !!
    !! For further information regarding the term "driving force", refer to:
    !!
    !!      H.L. Lukas, S.G. Fries and B. Sundman, "Computational Thermodynamics: The Calphad Method,"
    !!      Cambridge University Press, New York (2007).
    !!
    !! The premise of the subminimization method is to determine a unique combination of mole fractions for a
    !! particular non-ideal solution phase that minimizes the driving force of that phase.  As an additional
    !! constraint, the sum of the mole fractions of this phase must equal unity.  In the subminimization approach,
    !! a new Lagrangian function is defined as:
    !!
    !! \f$ L_{\lambda} = \sum_{i=1}^{N_{\lambda}} x_{i(\lambda)}^n \left( \mu_{i(\lambda)}^n -
    !! \sum_{j=1}^E a_{i,j} \Gamma_j^m\right) - \pi_{\lambda}^n \left( \sum_{i=1}^{N_{\lambda}}
    !! x_{i(\lambda)}^n - 1 \right)  \f$
    !!
    !! which solves for \f$ N_{\lambda} + 1 \f$ variables corresponding to \f$ x_{i(\lambda)}^n \f$ and
    !! \f$ \pi_{\lambda}^n  \f$.  For more information regarding this approach, refer to the following literature:
    !!
    !!      C.E. Harvie, J.P. Greenberg and J.H. Weare, "A Chemical Equilibrium Algorithm for Highly Non-Ideal
    !!      Multiphase Systems: Free Energy Minimization," Geochimica et Cosmochimica Acta, V. 51 (1987)
    !!      1045-1057.
    !!
    !
    !  In a previous version of Thermochimica, a solution phase was added to the system based on a different
    !  approach.  The previous method would compute the mole fractions of solution phase constituents as a
    !  function of the element potentials and that phase would be added to the system if the sum of mole
    !  fractions of all constituents in that phase was greater than unity.   The biggest issue with doing
    !  this is that the mole fraction of a constituent in a non-ideal solution phase IS NOT a function of
    !  the element potentials due to the inclusion of the partial molar excess Gibbs energy of mixing term.
    !  The previous approach relied heavily on a convoluted technique of dampening excess terms due to the highly
    !  non-linear nature of the equations.  Essentially, the mole fraction is written in that aproach as an
    !  expontential function of mole fractions.
    !
    !  The Subminimization approach achieves the same net effect, but due to the simplicity of the approach, the
    !  numerics are much more robust and efficient.  One of the main issues with the previous approach is that
    !  the mole fractions may not be unique to the current estimated element potentials.  In the subminimization
    !  approach, the mole fractions are unique to the current estimated element potentials.  Furthermore, the
    !  computational expense associated with each global iteration is far less.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in]   iSolnPhaseIndex  An integer scalar representing the absolute index of the solution phase that
    !!                                is being considered.
    !> \param[out]  lPhasePass       A logical scalar indicating whether the phase should be added (i.e., TRUE)
    !                                 to the system or not (i.e., FALSE).
    !
    ! lDuplicate                    A logical scalar indicating whether the mole fractions for this "phase" are
    !                                a duplicate of the mole fractions of the corresponding solution phase.
    !
    !
    !-------------------------------------------------------------------------------------------------------------
!
    USE ModuleThermo
    USE ModuleSubMin
    uSE ModuleGEMSolver!, ONLY: lMiscibility, dDrivingForceSoln
    USE ModuleThermoIO
!
    implicit none
!
    integer :: iSolnPhaseIndex, iterSubMax, i,j, k, l
    logical :: lPhasePass, lDuplicate
!
    ! Initialize local variables:
    lPhasePass = .FALSE.
    lDuplicate = .FALSE.
!
!
    ! Initializa Global variables
    iFirstSUB            = nSpeciesPhase(iSolnPhaseIndex-1) + 1
    iLastSUB             = nSpeciesPhase(iSolnPhaseIndex)
    nVar                 = iLastSUB - iFirstSUB + 1
    iterSubLg =0
    iterSubAdam =0
    
!
    

    ! Initialize the subminimization subroutine:
    call SubMinInit(iSolnPhaseIndex,iterSubMax)

    ! Subminimization iteration loop:
    LOOP_IterSub: do iterSub = 1, iterSubMax
    ! LOOP_IterSub: do iterSub = 1, 200
        ! Hyrbid algorithm
        ! if((nVar>1)&
        ! .and.(.NOT.lNegativeFraction)&
        ! .and.(dMaxElementPotential>1D-2)) then
        !     ! Update next iteration
        !     call SubMinAdam(iSolnPhaseIndex)
        !     iterSubAdam = iterSubAdam +1
        ! else
        !     ! Compute the direction vector:
        !     call SubMinNewton(iSolnPhaseIndex)

        !     ! Compute an appropriate step length:
        !     call SubMinLineSearch(iSolnPhaseIndex)
        !     iterSubLg = iterSubLg+1
        ! end if

        ! ! Adam only
        ! ! Update next iteration
        ! call SubMinAdam(iSolnPhaseIndex)

        ! Lagrangian only
        ! Compute the direction vector:
        call SubMinNewton(iSolnPhaseIndex)

        ! Compute an appropriate step length:
        call SubMinLineSearch(iSolnPhaseIndex)

        ! Compute the functional norm:
        call SubMinFunctionNorm(iSolnPhaseIndex)

        ! In case when handling zeroing species fails, exit the loop:
        if (MINVAL(dMolFraction(iFirstSUB:iLastSUB)) <= 0D0) exit LOOP_IterSub
        

        ! Only check for convergence if the functional norm is below tolerance:
        if ((dSubMinFunctionNorm < dTolerance(1)) &
        .and.(dMaxPotentialVector<dMaxPotentialTol)&
        ) lSubMinConverged = .TRUE.

        ! Check if the solution phases represening the miscibility gap duplicate one another:
        if (lMiscibility(iSolnPhaseIndex)) call SubMinCheckDuplicate(lDuplicate)
!
        ! Exit if the subminimization has converged:
        if (lDuplicate) exit LOOP_IterSub
!
        ! Exit if the subminimization has converged:
        if ((lSubMinConverged).OR.(INFOThermo /= 0)) exit LOOP_IterSub

        

    end do LOOP_IterSub
    ! if(iSolnPhaseIndex==11 .or. iSolnPhaseIndex==12) then
    !     print*, 'Subminimization did not converge ',iterSub,lSubMinConverged,iSolnPhaseIndex, cSolnPhaseName(iSolnPhaseIndex), dSubMinFunctionNorm, dMaxPotentialVector, dMolFraction(iFirstSUB:iLastSUB)
    ! end if

    
    ! If the composition of phases representing a miscibility gap duplicate one another (i.e., they have virtually
    ! the same composition), then set the driving force to zero to prevent this phase from being added to the system.
    if (lDuplicate) dDrivingForce = 9D5
!
    ! If the driving force is less than a specified tolerance and the system has converged,
    ! add this solution phase to the system:
    if ((dDrivingForce < dTolerance(4)).AND.(lSubMinConverged)) lPhasePass = .TRUE.
!
    ! If the mole fraction and charge neutrality constraints were not satisfied, then set the driving force to zero:

    if (dSubMinFunctionNorm > 10D0*dTolerance(1)) dDrivingForce = 9D5
!
    ! Update the solution driving force:
    dDrivingForceSoln(iSolnPhaseIndex)  = dDrivingForce
    dPhasePotential(iFirstSUB:iLastSUB) = dDrivingForce
!
    ! Deallocate allocatable arrays:
    deallocate(dChemicalPotentialStar, dRHS)
!
    deallocate(dHessian, iHessian)
!
    return
!
end subroutine Subminimization
!
    !---------------------------------------------------------------------------
    !                       END - Subminimization.f90
    !---------------------------------------------------------------------------
