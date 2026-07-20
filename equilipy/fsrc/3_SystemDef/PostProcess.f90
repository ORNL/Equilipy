!> \brief Perform equilibrium postprocessing and heat-capacity endpoint solves.
!!
!! \details Converts the converged internal minimizer state into reported
!! system, phase, species, and mass properties.  When heat capacity is
!! requested, the fixed phase assemblage is first solved at T-dT and then
!! solved again at the requested T so the reported state remains the final
!! requested-temperature Lagrangian state.
!
!-------------------------------------------------------------------------------------------------------------
!
!> \file    PostProcess.f90
!> \brief   Perform post-processing of equilibrium results.
!> \author  M.H.A. Piro
!> \date    January 14, 2013
!> \sa      RunLagrangianGEM.f90
!
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   01/14/2013      M.H.A. Piro         Original code.
    !   12/14/2024      S.Y. Kwon           Revised code.
    !   07/20/2026      S.Y. Kwon           Accumulated same-parent composition sets from slot-local properties and preserved terminal solver failures.
    !
    !
! Purpose:
! ========
!
!> \details The purpose of this subroutine is to perform post-processing of
!! the selected equilibrium assemblage.  Heat capacity is evaluated from a
!! fixed-assemblage finite-temperature perturbation, but the last Lagrangian
!! solve and the reported thermochemical state must remain at the requested
!! temperature.
!
!
! Required input variables:
! =========================
!
! iAssemblage              Active phase assemblage from MultiPhaseMinimizer.
! dTemperature             Requested temperature on entry.
! dTemperatureDiff         Temperature perturbation used for heat capacity.
! dNormalizeInput          Internal input normalization factor.
! dMolesPhase/dMolFraction Active phase amounts and solution fractions.
!
!
! Output/updated variables:
! =========================
!
! dGibbsEnergySys          Reported system Gibbs energy at requested T.
! dEnthalpySys             Reported system enthalpy at requested T.
! dEntropySys              Reported system entropy at requested T.
! dHeatCapacitySys         Reported fixed-assemblage heat capacity.
! dPostProcessGEMNormAtT   Fixed-assemblage GEM norm at requested T.
! dPostProcessGEMNormPerturbed
!                           Fixed-assemblage GEM norm at T-dT.
! dMolesPhase/dMolesSpecies Reported requested-T phase and species amounts.
!
!
! Called subroutines/functions:
! =============================
!
! CompThermoData            Recomputes standard thermodynamic data after a temperature change.
! CompChemicalPotential     Updates species chemical potentials and partial thermodynamic properties.
! CompFunctionNorm          Updates GEM residuals before fixed-assemblage postprocess solves.
! RunLagrangianGEM          Refines the fixed active assemblage without PEA while lPostProcess is true.
!
!
! Primary callers:
! ================
!
! minimize                 Calls this after MultiPhaseMinimizer/GEMSolver.
!
!
! Numerical assumptions:
! ======================
!
! - Heat-capacity perturbations keep the active phase assemblage fixed.
! - The perturbed T-dT solve is performed first; the requested-T solve is
!   performed last so reported arrays are not restored from a stale snapshot.
! - Solution species chemical potentials and partial properties are already
!   on the same formula-stoichiometry basis used in the mass-balance equations;
!   PostProcess must not divide their phase contribution by formula atom count.
! - Phase and species amounts are normalized to the original user input only
!   after the final requested-T state has been accumulated.
! - Multiple composition sets of one ordered parent share the legacy phase
!   arrays, so their integral and mass properties use slot-local snapshots.
!
!-------------------------------------------------------------------------------------------------------------



subroutine PostProcess

    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver
!
    implicit none
    integer:: i, j, k, iFirstPP, iLastPP
    integer:: iPostProcessFailureStatus, iPostProcessFailureReason
    integer, dimension(:), allocatable :: iAssemblageBeforePostProcess
    real(8):: dTemp, dOutputScale, dTemperatureAtT, dTemperatureStepPP
    real(8):: dPostProcessNormTolerance
    real(8):: dEnthalpySysTemp, dEntropySysTemp, dGibbsEnergySysTemp
    real(8):: dEnthalpySysAtT, dEntropySysAtT, dGibbsEnergySysAtT, dHeatCapacitySysAtT
    logical:: lCompEveryPhases, lUseSlotProperties, lPostProcessFailure

!
    !1. Recalculate Elemental potential if nElement/=nElementSys
    !2. Calculate dMolesPhase, dMoleFraction, dPhaseActivity, dBondFraction, dTotalBonds
    !3. Calcuate Thermochemical variables for the system: H, G, V, S, Cp, T, P
!
    ! Initiate variables:
    dGibbsEnergySysAtT  = 0d0
    dEnthalpySysAtT     = 0d0
    dEntropySysAtT      = 0d0
    dHeatCapacitySysAtT = 0d0
    dGibbsEnergySysTemp = 0d0
    dEnthalpySysTemp    = 0d0
    dEntropySysTemp     = 0d0
    dOutputScale        = 1D0 / dNormalizeInput
    dTemperatureAtT     = dTemperature
    dTemperatureStepPP  = dTemperatureDiff
    dPostProcessNormTolerance = 1D-10
    dPostProcessGEMNormAtT = dGEMFunctionNorm
    iPostProcessIterAtT    = iterGlobal
    dPostProcessGEMNormPerturbed = dGEMFunctionNorm
    iPostProcessIterPerturbed    = iterGlobal
    lCompEveryPhases = .False.
    ! Postprocess may refine a fixed assemblage for property evaluation, but it
    ! cannot turn an unsettled minimizer state into a successful equilibrium.
    ! Preserve the incoming raw failure before either endpoint solve resets the
    ! Lagrangian status fields.
    lPostProcessFailure = (.NOT.lConverged).OR.&
        (iGEMExitStatus /= GEM_EXIT_STATUS_OK)
    if (lPostProcessFailure) then
        iPostProcessFailureStatus = iGEMExitStatus
        if (iPostProcessFailureStatus == GEM_EXIT_STATUS_OK) &
            iPostProcessFailureStatus = GEM_EXIT_STATUS_LAGRANGIAN_UNCONVERGED
        iPostProcessFailureReason = iPhaseChangeReason
    else
        iPostProcessFailureStatus = GEM_EXIT_STATUS_OK
        iPostProcessFailureReason = PHASE_CHANGE_REASON_NONE
    end if

    lPostProcessAssemblageChanged = .FALSE.

    if (dTemperatureStepPP >= 1D-5) then
        allocate(iAssemblageBeforePostProcess(nElements))
        iAssemblageBeforePostProcess = iAssemblage

        ! Solve the fixed assemblage at T-dT first.  This gives the perturbed
        ! endpoint for heat capacity without leaving the reported state there.
        dTemperature = dTemperatureAtT - dTemperatureStepPP
        call CompThermoData
        call CompChemicalPotential(lCompEveryPhases)
        call FitCompoundOnlyPostProcessPlane
        call CompFunctionNorm
        lPhaseChange = .FALSE.
        if (dGEMFunctionNorm > dPostProcessNormTolerance) then
            lPostProcess = .TRUE.
            iGEMLagrangianCallSite = 3
            call RunLagrangianGEM
            if ((.NOT.lConverged).OR.&
                (iGEMExitStatus /= GEM_EXIT_STATUS_OK)) then
                lPostProcessFailure = .TRUE.
                iPostProcessFailureStatus = iGEMExitStatus
                iPostProcessFailureReason = iPhaseChangeReason
            end if
            iGEMLagrangianCallSite = 1
            call CompChemicalPotential(lCompEveryPhases)
            call CompFunctionNorm
            lPostProcess = .FALSE.
        else
            iterGlobal = 0
        end if
        call CompChemicalPotential(lCompEveryPhases)
        dPostProcessGEMNormPerturbed = dGEMFunctionNorm
        iPostProcessIterPerturbed    = iterGlobal

        lPostProcessAssemblageChanged = ANY(iAssemblage /= iAssemblageBeforePostProcess)

        dGibbsEnergySysTemp  = 0d0
        dEnthalpySysTemp     = 0d0
        dEntropySysTemp      = 0d0
        dMolesSpecies        = 0D0

        ! Compute the perturbed T-dT integral properties:
        do i =1,nElements
            j=iAssemblage(i)
            if(j>0) then
                ! Stoichimetric compound phase
                dTemp= dMolesPhase(i)
                dMolesSpecies(j) = dMolesPhase(i)
                if (dTemp>0D0) then
                    dEntropySysTemp     = dEntropySysTemp +dPartialEntropy(j)*dTemp
                    dEnthalpySysTemp    = dEnthalpySysTemp +dPartialEnthalpy(j)*dTemp
                    dGibbsEnergySysTemp = dGibbsEnergySysTemp +dChemicalPotential(j)*dTemp
                end if

            else if(j<0) then
                ! Solution phase
                dTemp= dMolesPhase(i)
                iFirstPP = nSpeciesPhase(-j-1) + 1
                iLastPP = nSpeciesPhase(-j)
                if (.NOT.PartitionPhaseActive(-j)) then
                    dMolesSpecies(iFirstPP:iLastPP)= dTemp*dMolFraction(iFirstPP:iLastPP)
                    do k =iFirstPP, iLastPP
                        if(dTemp>0D0.and.dMolesSpecies(k)>0D0) then
                            dEntropySysTemp     = dEntropySysTemp + &
                                dPartialEntropy(k)*dMolesSpecies(k)
                            dEnthalpySysTemp    = dEnthalpySysTemp + &
                                dPartialEnthalpy(k)*dMolesSpecies(k)
                            dGibbsEnergySysTemp = dGibbsEnergySysTemp + &
                                dChemicalPotential(k)*dMolesSpecies(k)
                        end if
                    end do
                else
                    lUseSlotProperties = UseActiveSlotProperties(i, -j)
                    if (lUseSlotProperties) then
                        dMolesSpecies(iFirstPP:iLastPP) = dMolesSpecies(iFirstPP:iLastPP) + &
                            dTemp*dActiveSlotMolFraction(i,iFirstPP:iLastPP)
                    else
                        dMolesSpecies(iFirstPP:iLastPP) = dTemp*dMolFraction(iFirstPP:iLastPP)
                    end if

                    do k =iFirstPP, iLastPP
                        if(dTemp>0D0.and.dMolesSpecies(k)>0D0) then
                            if (lUseSlotProperties) then
                                dEntropySysTemp = dEntropySysTemp + dActiveSlotPartialS(i,k) * &
                                    dTemp*dActiveSlotMolFraction(i,k)
                                dEnthalpySysTemp = dEnthalpySysTemp + dActiveSlotPartialH(i,k) * &
                                    dTemp*dActiveSlotMolFraction(i,k)
                                dGibbsEnergySysTemp = dGibbsEnergySysTemp + dActiveSlotChemPot(i,k) * &
                                    dTemp*dActiveSlotMolFraction(i,k)
                            else
                                dEntropySysTemp = dEntropySysTemp + dPartialEntropy(k)*dTemp*dMolFraction(k)
                                dEnthalpySysTemp = dEnthalpySysTemp + dPartialEnthalpy(k)*dTemp*dMolFraction(k)
                                dGibbsEnergySysTemp = dGibbsEnergySysTemp + &
                                    dChemicalPotential(k)*dTemp*dMolFraction(k)
                            end if
                        end if
                    end do
                end if
            end if
        end do
        dGibbsEnergySysTemp = dGibbsEnergySysTemp* dIdealConstant * dTemperature*dOutputScale
        dEnthalpySysTemp    = dEnthalpySysTemp *dIdealConstant * dTemperature*dOutputScale
        dEntropySysTemp     = dEntropySysTemp *dIdealConstant*dOutputScale

        ! Solve the requested T endpoint last so reported thermodynamic
        ! properties and arrays are computed from the final requested-T state.
        dTemperature = dTemperatureAtT
        call CompThermoData
        call CompChemicalPotential(lCompEveryPhases)
        call FitCompoundOnlyPostProcessPlane
        call CompFunctionNorm
        lPhaseChange = .FALSE.
        if (dGEMFunctionNorm > dPostProcessNormTolerance) then
            lPostProcess = .TRUE.
            iGEMLagrangianCallSite = 4
            call RunLagrangianGEM
            if (((.NOT.lConverged).OR.&
                (iGEMExitStatus /= GEM_EXIT_STATUS_OK)).AND.&
                (.NOT.lPostProcessFailure)) then
                lPostProcessFailure = .TRUE.
                iPostProcessFailureStatus = iGEMExitStatus
                iPostProcessFailureReason = iPhaseChangeReason
            end if
            iGEMLagrangianCallSite = 1
            call CompChemicalPotential(lCompEveryPhases)
            call CompFunctionNorm
            lPostProcess = .FALSE.
        else
            iterGlobal = 0
        end if
        call CompChemicalPotential(lCompEveryPhases)
        dPostProcessGEMNormAtT = dGEMFunctionNorm
        iPostProcessIterAtT    = iterGlobal
        lPostProcessAssemblageChanged = lPostProcessAssemblageChanged .OR. &
            ANY(iAssemblage /= iAssemblageBeforePostProcess)

        deallocate(iAssemblageBeforePostProcess)
    else
        call CompChemicalPotential(lCompEveryPhases)
        call FitCompoundOnlyPostProcessPlane
        call CompFunctionNorm
        dPostProcessGEMNormPerturbed  = dPostProcessGEMNormAtT
        iPostProcessIterPerturbed     = iPostProcessIterAtT
    end if

    ! Compute requested-T integral properties from the final requested-T state.
    dGibbsEnergySysAtT  = 0d0
    dEnthalpySysAtT     = 0d0
    dEntropySysAtT      = 0d0
    dHeatCapacitySysAtT = 0d0
    dMolesSpecies       = 0D0

    do i =1,nElements
        j=iAssemblage(i)
        if(j>0) then
            ! Stoichimetric compound phase
            dTemp= dMolesPhase(i)
            dMolesSpecies(j) = dMolesPhase(i)
            if (dTemp>0D0) then
                dGibbsEnergySysAtT  = dGibbsEnergySysAtT +dChemicalPotential(j)*dTemp
                dEnthalpySysAtT     = dEnthalpySysAtT +dPartialEnthalpy(j)*dTemp
                dEntropySysAtT      = dEntropySysAtT +dPartialEntropy(j)*dTemp
                dHeatCapacitySysAtT = dHeatCapacitySysAtT + dPartialHeatCapacity(j)*dTemp
            end if
        else if(j<0) then
            ! Solution phase
            dTemp= dMolesPhase(i)
            iFirstPP = nSpeciesPhase(-j-1) + 1
            iLastPP = nSpeciesPhase(-j)
            if (.NOT.PartitionPhaseActive(-j)) then
                dMolesSpecies(iFirstPP:iLastPP)= dTemp*dMolFraction(iFirstPP:iLastPP)
                do k =iFirstPP, iLastPP
                    if(dTemp>0D0.and.dMolesSpecies(k)>0D0) then
                        dEnthalpySysAtT     = dEnthalpySysAtT + &
                            dPartialEnthalpy(k)*dMolesSpecies(k)
                        dEntropySysAtT      = dEntropySysAtT + &
                            dPartialEntropy(k)*dMolesSpecies(k)
                        dGibbsEnergySysAtT  = dGibbsEnergySysAtT + &
                            dChemicalPotential(k)*dMolesSpecies(k)
                        dHeatCapacitySysAtT = dHeatCapacitySysAtT + &
                            dPartialHeatCapacity(k)*dMolesSpecies(k)
                    end if
                end do
            else
                lUseSlotProperties = UseActiveSlotProperties(i, -j)
                if (lUseSlotProperties) then
                    dMolesSpecies(iFirstPP:iLastPP) = dMolesSpecies(iFirstPP:iLastPP) + &
                        dTemp*dActiveSlotMolFraction(i,iFirstPP:iLastPP)
                else
                    dMolesSpecies(iFirstPP:iLastPP) = dTemp*dMolFraction(iFirstPP:iLastPP)
                end if
                do k =iFirstPP, iLastPP
                    if(dTemp>0D0.and.dMolesSpecies(k)>0D0) then
                        if (lUseSlotProperties) then
                            dEnthalpySysAtT = dEnthalpySysAtT + dActiveSlotPartialH(i,k) * &
                                dTemp*dActiveSlotMolFraction(i,k)
                            dEntropySysAtT = dEntropySysAtT + dActiveSlotPartialS(i,k) * &
                                dTemp*dActiveSlotMolFraction(i,k)
                            dGibbsEnergySysAtT = dGibbsEnergySysAtT + dActiveSlotChemPot(i,k) * &
                                dTemp*dActiveSlotMolFraction(i,k)
                            dHeatCapacitySysAtT = dHeatCapacitySysAtT + dActiveSlotPartialCp(i,k) * &
                                dTemp*dActiveSlotMolFraction(i,k)
                        else
                            dEnthalpySysAtT = dEnthalpySysAtT + dPartialEnthalpy(k)*dTemp*dMolFraction(k)
                            dEntropySysAtT = dEntropySysAtT + dPartialEntropy(k)*dTemp*dMolFraction(k)
                            dGibbsEnergySysAtT = dGibbsEnergySysAtT + &
                                dChemicalPotential(k)*dTemp*dMolFraction(k)
                            dHeatCapacitySysAtT = dHeatCapacitySysAtT + &
                                dPartialHeatCapacity(k)*dTemp*dMolFraction(k)
                        end if
                    end if
                end do
            end if
        end if
    end do
    dGibbsEnergySysAtT  = dGibbsEnergySysAtT* dIdealConstant * dTemperature*dOutputScale
    dEnthalpySysAtT     = dEnthalpySysAtT*dIdealConstant*dTemperature* dOutputScale
    dEntropySysAtT      = dEntropySysAtT*dIdealConstant*dOutputScale
    dHeatCapacitySysAtT = dHeatCapacitySysAtT*dIdealConstant*dOutputScale

    if (dTemperatureStepPP >= 1D-5) then
        dHeatCapacitySys = ((dEnthalpySysAtT - dEnthalpySysTemp) + &
            dTemperatureAtT*(dEntropySysAtT - dEntropySysTemp))/ &
            (2D0*dTemperatureStepPP)
    else
        dHeatCapacitySys= dHeatCapacitySysAtT
    end if

    dNormalizeInput = dOutputScale
    dMolesPhase     = dMolesPhase   * dNormalizeInput
    dMolesElement   = dMolesElement * dNormalizeInput
    dElementMass    = dElementMass  * dNormalizeInput
    dMolesSpecies   = dMolesSpecies * dNormalizeInput
    dGibbsEnergySys = dGibbsEnergySysAtT
    dEnthalpySys    = dEnthalpySysAtT
    dEntropySys     = dEntropySysAtT

    ! Compute species activities from the dimensionless chemical potentials.
    ! Python result capture copies this array directly.
    if (.not. allocated(dActivity)) allocate(dActivity(nSpecies))
    dActivity = 0D0
    do i = 1,nSpecies
        dTemp = dChemicalPotential(i) - dStdGibbsEnergy(i)
        dTemp = DMAX1(DMIN1(dTemp, 709D0), -745D0)
        dActivity(i) = DEXP(dTemp)
    end do

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
            if (.NOT.PartitionPhaseActive(-j)) then
                dGramPhase(i)=sum(dGramSpecies(iFirstPP:iLastPP))
            else if (UseActiveSlotProperties(i, -j)) then
                dGramPhase(i) = dMolesPhase(i) * SUM(&
                    dActiveSlotMolFraction(i,iFirstPP:iLastPP) * &
                    MATMUL(dStoichSpecies(iFirstPP:iLastPP,:),dAtomicMass))
            else
                dGramPhase(i)=sum(dGramSpecies(iFirstPP:iLastPP))
            end if
        end if
    end do
    
    do i = 1,nSolnPhasesSys
        iFirstPP = nSpeciesPhase(i-1) + 1
        iLastPP = nSpeciesPhase(i)
        dGramFraction(iFirstPP:iLastPP) = dMolFraction(iFirstPP:iLastPP)*matmul(dStoichSpecies(iFirstPP:iLastPP,:),dAtomicMass)
        dGramFraction(iFirstPP:iLastPP) = dGramFraction(iFirstPP:iLastPP)/sum(dGramFraction(iFirstPP:iLastPP))
    end do

    ! The requested-T solve must not erase a failure at T-dT.  Postprocess
    ! properties may have been populated for diagnostics, but the raw terminal
    ! status remains failed and Python exposes the stage and reason.
    if (lPostProcessFailure) then
        lConverged = .FALSE.
        iGEMExitStatus = iPostProcessFailureStatus
        if (iGEMExitStatus == GEM_EXIT_STATUS_OK) &
            iGEMExitStatus = GEM_EXIT_STATUS_LAGRANGIAN_UNCONVERGED
        iPhaseChangeReason = iPostProcessFailureReason
    end if




   

    return
!
contains

    logical function UseActiveSlotProperties(iSlotIn, iPhaseIn)

        implicit none

        integer, intent(in) :: iSlotIn, iPhaseIn
        integer             :: iOtherSlot

        UseActiveSlotProperties = .FALSE.
        if (.NOT.PartitionPhaseActive(iPhaseIn)) return
        if (.NOT.allocated(lActiveSlotPropValid)) return
        if (.NOT.allocated(iActiveSlotThermoPhase)) return
        if (.NOT.allocated(dActiveSlotMolFraction)) return
        if (.NOT.allocated(dActiveSlotChemPot)) return
        if (.NOT.allocated(dActiveSlotPartialH)) return
        if (.NOT.allocated(dActiveSlotPartialS)) return
        if (.NOT.allocated(dActiveSlotPartialCp)) return
        if ((iSlotIn < 1).OR.(iSlotIn > SIZE(lActiveSlotPropValid))) return
        if (iActiveSlotThermoPhase(iSlotIn) /= iPhaseIn) return
        if (.NOT.lActiveSlotPropValid(iSlotIn)) return
!
!       Global phase arrays remain the authoritative legacy path for one
!       composition set.  Slot-local properties are required only when two
!       active slots share one thermodynamic parent and would otherwise
!       overwrite each other's constitution and partial properties.
        do iOtherSlot = 1, MIN(nElements,SIZE(iActiveSlotThermoPhase))
            if (iOtherSlot == iSlotIn) cycle
            if (iAssemblage(iOtherSlot) /= -iPhaseIn) cycle
            if (iActiveSlotThermoPhase(iOtherSlot) /= iPhaseIn) cycle
            UseActiveSlotProperties = .TRUE.
            return
        end do

        return
    end function UseActiveSlotProperties


    logical function PartitionPhaseActive(iPhaseIn)

        implicit none

        integer, intent(in) :: iPhaseIn

        PartitionPhaseActive = .FALSE.
        if (.NOT.lODPartitionUnifiedActive) return
        if (.NOT.allocated(iODTopologyClass)) return
        if ((iPhaseIn < 1).OR.(iPhaseIn > SIZE(iODTopologyClass))) return
        PartitionPhaseActive = &
            (iODTopologyClass(iPhaseIn) >= OD_TOPOLOGY_HELPER_STANDALONE).AND.&
            (iODTopologyClass(iPhaseIn) <= OD_TOPOLOGY_HELPER_ONLY)

        return
    end function PartitionPhaseActive

    subroutine FitCompoundOnlyPostProcessPlane

        implicit none

        integer :: iLocal, jLocal, iSpeciesLocal
        integer :: M, N, NRHS, LDA, LDB, LWORK, INFO
        real(8) :: dWorkQuery(1), dResidualMax, dResidualLocal
        real(8), allocatable :: ALocal(:,:), BLocal(:), WORK(:)

        if (nSolnPhases /= 0) return
        if (nConPhases <= 0) return
        if (nElements <= 0) return

        M = nConPhases
        N = nElements
        NRHS = 1
        LDA = MAX(1,M)
        LDB = MAX(1,M,N)
        LWORK = -1

        allocate(ALocal(LDA,N), BLocal(LDB))
        ALocal = 0D0
        BLocal = 0D0

        do iLocal = 1, nConPhases
            iSpeciesLocal = iAssemblage(iLocal)
            if ((iSpeciesLocal <= 0).OR.(iSpeciesLocal > nSpecies)) then
                deallocate(ALocal, BLocal)
                return
            end if
            ALocal(iLocal,1:nElements) = dStoichSpecies(iSpeciesLocal,1:nElements)
            BLocal(iLocal) = dChemicalPotential(iSpeciesLocal)
        end do

        call DGELS('N', M, N, NRHS, ALocal, LDA, BLocal, LDB, &
            dWorkQuery, LWORK, INFO)
        if (INFO == 0) then
            LWORK = MAX(1,INT(dWorkQuery(1)))
            allocate(WORK(LWORK))
            call DGELS('N', M, N, NRHS, ALocal, LDA, BLocal, LDB, &
                WORK, LWORK, INFO)
            if (allocated(WORK)) deallocate(WORK)
        end if

        if (INFO == 0) then
            dResidualMax = 0D0
            do iLocal = 1, nConPhases
                iSpeciesLocal = iAssemblage(iLocal)
                dResidualLocal = -dChemicalPotential(iSpeciesLocal)
                do jLocal = 1, nElements
                    dResidualLocal = dResidualLocal + &
                        BLocal(jLocal) * dStoichSpecies(iSpeciesLocal,jLocal)
                end do
                dResidualMax = MAX(dResidualMax, DABS(dResidualLocal))
            end do

            if (dResidualMax <= 1D-8) then
                dElementPotential(1:nElements) = BLocal(1:nElements)
            end if
        end if

        if (allocated(ALocal)) deallocate(ALocal)
        if (allocated(BLocal)) deallocate(BLocal)

        return

    end subroutine FitCompoundOnlyPostProcessPlane

end subroutine PostProcess
!
!
