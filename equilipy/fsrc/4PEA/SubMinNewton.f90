!
subroutine SubMinNewton(iSolnPhaseIndex)
    !---------------------------------------------------------------------------
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this subroutine is to compute the direction vector using
    ! Lagrangian multiplier method. Here the Gibbs energy is based on one mole
    ! atomsl
    !
    ! Pertinent variables:
    ! ====================
    !
    ! INFO          An integer scalar used by LAPACK to identify a successful
    !                exit (INFO = 0) or an error (INFO /= 0).
    ! nEqn          The number of equations in the Hessian matrix.
    ! dRHS          A double real vector representing the right hand side (RHS)
    !                of the Hessian matrix and returns the direction vector.
    !
    !---------------------------------------------------------------------------
    USE ModuleThermoIO, ONLY: INFOThermo
    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleSubMin
!
    implicit none
!
    integer :: i, j, k, iSolnPhaseIndex, INFO, nEqn, m
!
!
    ! Initialize variables:
    nEqn     = nVar + 1
    INFO     = 0
    iHessian = 0
    dHessian = 0D0
    dRHS     = -1D0
    m        = 1
!
!
    ! Construct diagonal and part of the arrow head:
    ! SY: Note that dRHS is not xi but (xi-yi)
    do j = 1, nVar
        i                 = nSpeciesPhase(iSolnPhaseIndex-1) + j
        dHessian(j,j)     = 1D0 / dMolFraction(i)
        dHessian(nEqn,j)  = -1D0
        dHessian(j,nEqn)  = -1D0
        dRHS(j)           = dDrivingForce - (dChemicalPotential(i) + 1 - dChemicalPotentialStar(j))
        dRHS(nEqn)        = dRHS(nEqn) + dMolFraction(i)
    end do
!
    ! Apply an additional row/column if the phase is ionic:
    if (iPhaseElectronID(iSolnPhaseIndex) /= 0) then
        nEqn = nEqn + 1
        k    = iPhaseElectronID(iSolnPhaseIndex)
        m    = 2
!
        dRHS(nEqn) = 0D0
        do j = 1, nVar
            i                 = nSpeciesPhase(iSolnPhaseIndex-1) + j
            dHessian(nEqn,j)  = -dStoichSpecies(i,k)
            dHessian(j,nEqn)  = dHessian(nEqn,j)
            dRHS(nEqn)        = dRHS(nEqn) + dStoichSpecies(i,k) * dMolFraction(i)
        end do
    end if
!
    ! Call the linear equation solver:
    call DGESV( nEqn, 1, dHessian, nEqn, iHessian, dRHS, nEqn, INFO )
!
!
    if (INFO /= 0) then
        ! Return an error and reset dRHS:
        INFOThermo = 28
        dRHS       = 0D0
    end if
!
end subroutine SubMinNewton
!
!
