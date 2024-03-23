!
!
!
subroutine PostProcess
!
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    PostProcess.f90
    !> \brief   Perform post-processing of results.
    !> \author  M.H.A. Piro
    !> \date    January 14, 2013
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   01/14/2013      M.H.A. Piro         Original code.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to perform post-procesing of results.
    !
    !
    ! Pertinent variables:
    ! ====================
!
    !
    !-------------------------------------------------------------------------------------------------------------
    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver
!
    implicit none
    integer:: i
    real(8):: dSpeciesTotalMole, dTemp
!
    !1. Recalculate Elemental potential if nElement/=nElementSys
    !2. Calculate dMolesPhase, dMoleFraction, dPhaseActivity, dBondFraction, dTotalBonds
    !3. Calcuate Thermochemical variables for the system: H, G, V, S, Cp, T, P
!
    ! Initiate variables:
    dGibbsEnergySys  = 0d0
    dEnthalpySys     = 0d0
    dEntropySys      = 0d0
!
!
    ! Multiply the number of moles of all phases by the normalizing constant:
    dNormalizeInput = 1D0 / dNormalizeInput
    dMolesPhase     = dMolesPhase   * dNormalizeInput
    dMolesElement   = dMolesElement * dNormalizeInput
    dMolesSpecies   = dMolesSpecies * dNormalizeInput
    dElementMass    = dElementMass  * dNormalizeInput
!
    ! Compute the integral Gibbs energy of the system:
    ! dGibbsEnergySys = dot_product(dElementPotential, dMolesElement) * dIdealConstant * dTemperature
    do i =1,nElements
        dGibbsEnergySys = dGibbsEnergySys +dElementPotential(i)*dMolesElement(i)
    end do
!
    do i =1,nSpecies
        dTemp= dMolesSpecies(i)
        if (dTemp>0D0) then
            dEntropySys = dEntropySys +dPartialEntropy(i)*dTemp
            dEnthalpySys = dEnthalpySys +dPartialEnthalpy(i)*dTemp
        end if
    end do
    dGibbsEnergySys = dGibbsEnergySys* dIdealConstant * dTemperature
    dEnthalpySys = dEnthalpySys *dIdealConstant * dTemperature
    dEntropySys = dEntropySys *dIdealConstant 

    return
!
end subroutine PostProcess
!
!
