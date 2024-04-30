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
    integer:: i, j, k, iFirstPP, iLastPP
    real(8):: dSpeciesTotalMole, dTemp, dEnthalpySysTemp
!
    !1. Recalculate Elemental potential if nElement/=nElementSys
    !2. Calculate dMolesPhase, dMoleFraction, dPhaseActivity, dBondFraction, dTotalBonds
    !3. Calcuate Thermochemical variables for the system: H, G, V, S, Cp, T, P
!
    ! Initiate variables:
    dEnthalpySysTemp = 0d0
    

    ! Store enthalpy value for Heat capacity calculation
    do i =1,nElements
        j=iAssemblage(i)
        if(j>0) then
            ! Stoichimetric compound phase
            dTemp= dMolesPhase(i)
            dMolesSpecies(j) = dMolesPhase(i)
            if (dTemp>0D0) then
                dEnthalpySysTemp    = dEnthalpySysTemp +dPartialEnthalpy(j)*dTemp
            end if
        else if(j<0) then
            ! Solution phase
            dTemp= dMolesPhase(i)
            iFirstPP = nSpeciesPhase(-j-1) + 1
            iLastPP = nSpeciesPhase(-j)
            dMolesSpecies(iFirstPP:iLastPP)= dTemp*dMolFraction(iFirstPP:iLastPP)
            do k =iFirstPP, iLastPP
                if(dTemp>0D0.and.dMolesSpecies(k)>0D0) then
                    dEnthalpySysTemp    = dEnthalpySysTemp +dPartialEnthalpy(k)*dMolesSpecies(k)
                end if
            end do  
        end if
    end do
    dEnthalpySysTemp = dEnthalpySysTemp*dIdealConstant*dTemperature/dNormalizeInput
    
    ! Recalculate phase equilibria by correcting temperature
    dTemperature = dTemperature + dTemperatureDiff
    call CompThermoData
    call CheckThermoData
    call GEMSolverNew

    dGibbsEnergySys  = 0d0
    dEnthalpySys     = 0d0
    dEntropySys      = 0d0
    dMolesSpecies    = 0D0
    dHeatCapacitySys = 0d0

    ! Multiply the number of moles of all phases by the normalizing constant:
    dNormalizeInput = 1D0 / dNormalizeInput
    dMolesPhase     = dMolesPhase   * dNormalizeInput
    dMolesElement   = dMolesElement * dNormalizeInput
    dMolesSpecies   = dMolesSpecies * dNormalizeInput
    dElementMass    = dElementMass  * dNormalizeInput
!
    ! Compute the integral Gibbs energy of the system:
    do i =1,nElements
        j=iAssemblage(i)
        if(j>0) then
            ! Stoichimetric compound phase
            dTemp= dMolesPhase(i)
            dMolesSpecies(j) = dMolesPhase(i)
            if (dTemp>0D0) then
                dEntropySys     = dEntropySys +dPartialEntropy(j)*dTemp
                dEnthalpySys    = dEnthalpySys +dPartialEnthalpy(j)*dTemp
                dGibbsEnergySys = dGibbsEnergySys +dChemicalPotential(j)*dTemp
            end if
            
        else if(j<0) then
            ! Solution phase
            dTemp= dMolesPhase(i)
            iFirstPP = nSpeciesPhase(-j-1) + 1
            iLastPP = nSpeciesPhase(-j)
            dMolesSpecies(iFirstPP:iLastPP)= dTemp*dMolFraction(iFirstPP:iLastPP)
            
            do k =iFirstPP, iLastPP
                if(dTemp>0D0.and.dMolesSpecies(k)>0D0) then
                    dEntropySys     = dEntropySys +dPartialEntropy(k)*dMolesSpecies(k)
                    dEnthalpySys    = dEnthalpySys +dPartialEnthalpy(k)*dMolesSpecies(k)
                    dGibbsEnergySys = dGibbsEnergySys +dChemicalPotential(k)*dMolesSpecies(k)
                end if
            end do  
        end if
    end do
    dGibbsEnergySys = dGibbsEnergySys* dIdealConstant * dTemperature
    dEnthalpySys = dEnthalpySys *dIdealConstant * dTemperature
    dEntropySys = dEntropySys *dIdealConstant 
    dHeatCapacitySys = (dEnthalpySys - dEnthalpySysTemp)/dTemperatureDiff
    
    return
!
end subroutine PostProcess
!
!
