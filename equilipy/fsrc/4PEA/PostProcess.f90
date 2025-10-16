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
    !   12/14/2024      S.Y. Kwon           Revised code.
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
    real(8):: dTemp, dEnthalpySysTemp, dEntropySysTemp, dGibbsEnergySysTemp, dPostNorm,dHeatCapacitySysTemp
    logical:: lCompEveryPhases

!
    !1. Recalculate Elemental potential if nElement/=nElementSys
    !2. Calculate dMolesPhase, dMoleFraction, dPhaseActivity, dBondFraction, dTotalBonds
    !3. Calcuate Thermochemical variables for the system: H, G, V, S, Cp, T, P
!
    ! Initiate variables:
    dGibbsEnergySysTemp = 0d0
    dEnthalpySysTemp    = 0d0
    dEntropySysTemp     = 0d0
    dHeatCapacitySysTemp= 0d0
    
    lCompEveryPhases = .True.

    ! Store enthalpy value for Heat capacity calculation
    do i =1,nElements
        j=iAssemblage(i)
        if(j>0) then
            ! Stoichimetric compound phase
            dTemp= dMolesPhase(i)
            dMolesSpecies(j) = dMolesPhase(i)
            if (dTemp>0D0) then
                dGibbsEnergySysTemp = dGibbsEnergySysTemp +dChemicalPotential(j)*dTemp
                dEnthalpySysTemp    = dEnthalpySysTemp +dPartialEnthalpy(j)*dTemp
                dEntropySysTemp     = dEntropySysTemp +dPartialEntropy(j)*dTemp
                dHeatCapacitySysTemp= dHeatCapacitySysTemp + dPartialHeatCapacity(j)*dTemp
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
                    dEntropySysTemp     = dEntropySysTemp +dPartialEntropy(k)*dMolesSpecies(k)
                    dGibbsEnergySysTemp = dGibbsEnergySysTemp +dChemicalPotential(k)*dMolesSpecies(k)
                    dHeatCapacitySysTemp= dHeatCapacitySysTemp + dPartialHeatCapacity(k)*dMolesSpecies(k)
                end if
            end do
        end if
    end do
    dNormalizeInput     = 1D0 / dNormalizeInput
    dGibbsEnergySysTemp = dGibbsEnergySysTemp* dIdealConstant * dTemperature*dNormalizeInput
    dEnthalpySysTemp    = dEnthalpySysTemp*dIdealConstant*dTemperature* dNormalizeInput
    dEntropySysTemp     =  dEntropySysTemp*dIdealConstant *dNormalizeInput
    dHeatCapacitySysTemp= dHeatCapacitySysTemp*dIdealConstant* dNormalizeInput


    if (dTemperatureDiff<1D-5) then
        
        dMolesPhase     = dMolesPhase   * dNormalizeInput
        dMolesElement   = dMolesElement * dNormalizeInput
        dElementMass    = dElementMass  * dNormalizeInput
        dGibbsEnergySys = dGibbsEnergySysTemp
        dEnthalpySys    = dEnthalpySysTemp
        dEntropySys     =  dEntropySysTemp
        dHeatCapacitySys= dHeatCapacitySysTemp
        
        
    else
        ! Revert the dNormalize Input
        dNormalizeInput = 1D0 / dNormalizeInput

        ! Recalculate phase equilibria by correcting temperature
        dTemperature = dTemperature + dTemperatureDiff
        call CompThermoData
        
        ! print*, 'PostProcess'
        
        ! print*, 'PostProcess'
        !1. Re-Calculate elemental potential
        !2. Re-Calculate driving force
        if (lCompbdOnly) then
            !Replace this part with Leveling
            call GEMSolver
        else
            
            ! Replace this part with lagrangian solver
            lPhaseChange = .False.
            lPostProcess = .True.
            call MultiPhaseMinimizer
            
        end if
        

        dGibbsEnergySys  = 0d0
        dEnthalpySys     = 0d0
        dEntropySys      = 0d0
        dMolesSpecies    = 0D0
        dHeatCapacitySys = 0d0

        ! Multiply the number of moles of all phases by the normalizing constant:
        dNormalizeInput = 1D0 / dNormalizeInput
        dMolesPhase     = dMolesPhase   * dNormalizeInput
        dElementMass    = dElementMass  * dNormalizeInput
        dMolesElement   = dMolesElement * dNormalizeInput
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
        dHeatCapacitySys = ((dEnthalpySys - dEnthalpySysTemp)+dTemperature*(dEntropySys - dEntropySysTemp))/(2*dTemperatureDiff)

        ! !Numerical convergence  check
        ! dPostNorm           = 0d0
        ! ! 1. G-H+TS=0
        ! dPostNorm = ABS(dGibbsEnergySys-dEnthalpySys+dTemperature*dEntropySys)
        
        ! ! 2. dG/dT+S=0
        ! dPostNorm = dPostNorm + ABS(dEntropySys+(dGibbsEnergySys-dGibbsEnergySysTemp)/dTemperatureDiff)

        ! ! 3. d(G/T)/d(T)+H/T^2=0
        ! dTemp = dTemperature-dTemperatureDiff
        ! dPostNorm = dPostNorm+ABS((dGibbsEnergySys/dTemperature-dGibbsEnergySysTemp/(dTemp))/(dTemperature-dTemp)-dEnthalpySys/(dTemperature**2))
        ! if (dPostNorm>1D-4) print*, 'Not converging?', dPostNorm
    end if

    ! Update mass fraction
    if(allocated(dGramFraction)) deallocate(dGramFraction)
    if(allocated(dGramSpecies)) deallocate(dGramSpecies)
    if(allocated(dGramPhase)) deallocate(dGramPhase)
    if(allocated(dGramElement)) deallocate(dGramElement)
    allocate(dGramElement(nElements),dGramPhase(nElements),dGramFraction(nSpecies),dGramSpecies(nSpecies))
    dGramFraction = 1d0
    dGramSpecies = 0d0
    dGramPhase = 0d0
    dGramElement = 0d0

    dGramElement = dMolesElement * dAtomicMass
    dGramSpecies = dMolesSpecies*matmul(dStoichSpecies,dAtomicMass)

    do i =1,nElements
        j=iAssemblage(i)
        if(j>0) then
            ! Stoichimetric compound phase
            dGramPhase(i)=dGramSpecies(j)
        else if(j<0) then
            ! Solution phase
            iFirstPP = nSpeciesPhase(-j-1) + 1
            iLastPP = nSpeciesPhase(-j)
            dGramPhase(i)=sum(dGramSpecies(iFirstPP:iLastPP))
        end if
    end do
    
    do i = 1,nSolnPhasesSys
        iFirstPP = nSpeciesPhase(i-1) + 1
        iLastPP = nSpeciesPhase(i)
        dGramFraction(iFirstPP:iLastPP) = dMolFraction(iFirstPP:iLastPP)*matmul(dStoichSpecies(iFirstPP:iLastPP,:),dAtomicMass)
        dGramFraction(iFirstPP:iLastPP) = dGramFraction(iFirstPP:iLastPP)/sum(dGramFraction(iFirstPP:iLastPP))
    end do




   

    return
!
end subroutine PostProcess
!
!
