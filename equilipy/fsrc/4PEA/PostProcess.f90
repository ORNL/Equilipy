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
    real(8):: dSpeciesTotalMole
!
    !1. Recalculate Elemental potential if nElement/=nElementSys
    !2. Calculate dMolesPhase, dMoleFraction, dPhaseActivity, dBondFraction, dTotalBonds
    !3. Calcuate Thermochemical variables for the system: H, G, V, S, Cp, T, P
!
    ! Initiate variables:
    dGibbsEnergySys  = 0d0
    dEnthalpySys     = 0d0
    dEntropySys      = 0d0
    dHeatCapacitySys = 0d0
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
    dGibbsEnergySys = dot_product(dElementPotential, dMolesElement) * dIdealConstant * dTemperature
    ! dEnthalpySys =  dot_product(dPartialEnthalpy, dMolesSpecies) * dIdealConstant * dTemperature
    ! dEntropySys =  dot_product(dPartialEntropy, dMolesSpecies)* dIdealConstant
    ! dHeatCapacitySys =  dot_product(dPartialHeatCapacity, dMolesSpecies) * dIdealConstant
!
    do i =1,nSpecies
        if (dMolesSpecies(i)>0D0) then
            ! dSpeciesTotalMole = sum(dStoichSpecies(i,:))
            ! dGibbsEnergySys= dGibbsEnergySys +dChemicalPotential(i)*dMolesSpecies(i)
            dEnthalpySys = dEnthalpySys +dPartialEnthalpy(i)*dMolesSpecies(i)
            dEntropySys = dEntropySys +dPartialEntropy(i)*dMolesSpecies(i)
            dHeatCapacitySys = dHeatCapacitySys +dPartialHeatCapacity(i)*dMolesSpecies(i)
        end if
    end do
    ! dGibbsEnergySys = dGibbsEnergySys* dIdealConstant * dTemperature
    dEnthalpySys = dEnthalpySys *dIdealConstant * dTemperature
    dEntropySys = dEntropySys *dIdealConstant 
    dHeatCapacitySys = dHeatCapacitySys *dIdealConstant

    return
!
end subroutine PostProcess
!
!
