
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CompFunctionNorm.f90
    !> \brief   Compute the functional norm for the line search.
    !> \author  M.H.A. Piro
    !> \date    Apr. 25, 2012
    !> \sa      GEMLineSearch.f90
    !> \sa      CompStoichSolnPhase.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   04/25/2012      M.H.A. Piro         Original code
    !   07/20/2026      S.Y. Kwon           Extended the GEM residual norm to CEF site exchange, phase complementarity, and trace-bound conditions.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to compute the functional norm for the line search algorithm to
    !! determine whether the system is converging sufficinetly or diverging. The functional vector is not directly
    !! computed because it is not needed. Thus, the functional norm is computed directly.  This term incorporates
    !! the residuals of the mass balance equations, the average residual between the chemical potential of each
    !! species and the corresponding value computed from the element potentials.  When the active fixed
    !! assemblage is solved by the CEF site-fraction KKT path, solution stationarity for CEF phases is measured
    !! by CompSublatticeExchangeNorm instead of by product-endmember chemical-potential residuals.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
! dGEMFunctionNorm      A double real scalar representing the functional norm of the functional  vector.  The
!                        functional vector represents the relative errors of the mass balance equations,
!                        the departure of chemical potentials of stable phases from the corresponding values
!                        computed from the element potentials.
! dGEMMassBalanceNorm   Square-rooted mass-balance part of dGEMFunctionNorm.
! dGEMChemicalPotentialNorm
!                       Square-rooted chemical-potential part of dGEMFunctionNorm, including solution,
!                       pure condensed, and sublattice exchange residuals.
    !
    !-------------------------------------------------------------------------------------------------------------


subroutine CompFunctionNorm

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer :: i, j, k, l, iSlot, iThermoPhase, iPhaseID
    real(8) :: dNormComponent, dMassNormSquared, dSolutionChemicalNormSquared
    real(8) :: dCondensedChemicalNormSquared
    real(8) :: dExchangeNormLocal
    real(8) :: dResidual, dResidualWeight
    real(8), allocatable :: dPhaseComp(:)
    logical :: lAddPhase, lUseCEFExchangeResidual, lTraceInactive

    ! Initialize variables:
    dGEMFunctionNorm    = 0D0
    dGEMMassBalanceNorm = 0D0
    dGEMChemicalPotentialNorm = 0D0
    dGEMSolutionChemicalPotentialNorm = 0D0
    dGEMCondensedChemicalPotentialNorm = 0D0
    dEffStoichSolnPhase = 0D0
    dMassNormSquared = 0D0
    dSolutionChemicalNormSquared = 0D0
    dCondensedChemicalNormSquared = 0D0
    lAddPhase = .False.
    allocate(dPhaseComp(nElements))

    ! Compute the residuals of the mass balance equations for each element:
    do j = 1, nElements
        dNormComponent = dMolesElement(j)
        do i = 1, nSolnPhases
            iSlot = nElements - i + 1
            k = -iAssemblage(iSlot)       ! Absolute displayed solution phase index
            iThermoPhase = ActiveSlotThermoPhaseForNorm(iSlot, k)
            if (lGEMCEFSiteLagrangianActive.AND.&
                IsNormCEFPhase(iThermoPhase)) then
                iPhaseID = iPhaseSublattice(iThermoPhase)
                call ComputeActiveSlotCEFPhaseComposition(iThermoPhase, iPhaseID, iSlot, dPhaseComp)
                dNormComponent = dNormComponent - dPhaseComp(j) * dMolesPhase(iSlot)
            else
                ! Compute the stoichiometry of this solution phase:
                call CompStoichSolnPhase(iThermoPhase)
                dNormComponent = dNormComponent - dEffStoichSolnPhase(iThermoPhase,j) * dMolesPhase(iSlot)
            end if
        end do

        do i = 1,nConPhases
            dNormComponent = dNormComponent - dMolesPhase(i) * dStoichSpecies(iAssemblage(i),j)
        end do
        dMassNormSquared = dMassNormSquared + (dNormComponent)**(2)
    end do
    ! if (dGEMFunctionNorm**(0.5)>1D-5) print*, 'mass balance', dGEMFunctionNorm**(0.5)
    ! Compute the residuals of the Gibbs energy difference between each solution phase and the element potentials:
    do l = 1, nSolnPhases
        iSlot = nElements - l + 1
        k = -iAssemblage(iSlot)       ! Absolute displayed solution phase index
        iThermoPhase = ActiveSlotThermoPhaseForNorm(iSlot, k)
        lUseCEFExchangeResidual = .FALSE.
        if (allocated(iPhaseSublattice)) then
            if ((iThermoPhase > 0).AND.(iThermoPhase <= SIZE(iPhaseSublattice))) then
                if (lGEMCEFSiteLagrangianEnabled .AND. lGEMCEFSiteLagrangianActive) then
                    lUseCEFExchangeResidual = IsNormCEFPhase(iThermoPhase)
                end if
            end if
        end if
        ! CEF phases solved by GEMNewtonCEF already contribute their phase and
        ! site-exchange residuals through CompSublatticeExchangeNorm.  Do not
        ! also test their product endmember chemical potentials here.
        if (lUseCEFExchangeResidual) cycle

        ! Postprocess is a fixed-assemblage property evaluation.  A solution
        ! phase with zero phase amount is not an active stationarity constraint
        ! for that fixed endpoint; retaining its old composition only pollutes
        ! the diagnostic norm without changing G/H/S/Cp or mass balance.
        if (lPostProcess.AND.(dMolesPhase(iSlot) <= DMAX1(dTolerance(8), 1D-300))) cycle

        do i = nSpeciesPhase(iThermoPhase-1) + 1, nSpeciesPhase(iThermoPhase)
            lTraceInactive = .FALSE.
            if (lTraceSpeciesControlEnabled .AND. allocated(lTraceSpeciesInactive)) then
                lTraceInactive = lTraceSpeciesInactive(i)
            end if
            dNormComponent = 0D0
            do j = 1, nElements
                dNormComponent = dNormComponent + dElementPotential(j) * dStoichSpecies(i,j)
            end do
            
            ! Normalize the residual term by the number of particles per formula mass:
            dNormComponent = dNormComponent / DFLOAT(iParticlesPerMole(i))
            dResidual = dChemicalPotential(i) - dNormComponent
            if (lTraceInactive.AND.(dResidual >= -dTraceSpeciesResidualTolerance)) cycle

            ! Active solution species satisfy complementarity at the lower
            ! bound.  If a near-zero species has a negative grand-potential
            ! residual, it wants to grow and must not be hidden by x_i = 0.
            dResidualWeight = dMolFraction(i)
            if ((dMolFraction(i) <= dTraceSpeciesRemoveFraction).AND.&
                (dResidual < -dTraceSpeciesResidualTolerance)) dResidualWeight = 1D0
            if (lTraceInactive.AND.(dResidual < -dTraceSpeciesResidualTolerance)) dResidualWeight = 1D0

            ! Compute the residual term weighted by the active complementarity measure:
            dNormComponent = DABS(dResidual) * dResidualWeight
            dSolutionChemicalNormSquared = dSolutionChemicalNormSquared + (dNormComponent)**(2)
        end do
        ! if (k==1) then
        !     print*, 'dGEMFunctionNorm', dGEMFunctionNorm**(0.5) ,k
        !     print*, 'dMolFraction(i)',dMolFraction(nSpeciesPhase(k-1) + 1: nSpeciesPhase(k))
        !     print*, 'dChemicalPotential(i)',dChemicalPotential(nSpeciesPhase(k-1) + 1: nSpeciesPhase(k))
        !     print*, 'Element', MATMUL(dStoichSpecies(nSpeciesPhase(k-1) + 1: nSpeciesPhase(k),:),dElementPotential(:))
        !     print*, 'Delta', dChemicalPotential(nSpeciesPhase(k-1) + 1: nSpeciesPhase(k))-MATMUL(dStoichSpecies(nSpeciesPhase(k-1) + 1: nSpeciesPhase(k),:),dElementPotential(:))
        ! end if
    end do
    ! if (dGEMFunctionNorm**(0.5)>1D-5)  print*, 'Chem balance, Sln', dGEMFunctionNorm**(0.5) 

    call CompSublatticeExchangeNorm(dExchangeNormLocal)
    dSublatticeExchangeNorm = dExchangeNormLocal


    ! Compute the residuals of the chemical potentials of pure condensed phases and the element potentials:
    do i = 1, nConPhases
        k     = iAssemblage(i)
        dNormComponent = 0D0
        do j = 1, nElements
            dNormComponent = dNormComponent + dElementPotential(j) * dStoichSpecies(k,j)
        end do
        dNormComponent            = dNormComponent - dChemicalPotential(k)
        dCondensedChemicalNormSquared = dCondensedChemicalNormSquared + (dNormComponent)**(2)
    end do
    ! if (dGEMFunctionNorm**(0.5)>1D-5)  print*, 'Chem balance, cmpd', dGEMFunctionNorm**(0.5)
    ! print*, 'Chem balance, cmpd', dGEMFunctionNorm**(0.5)
    
    ! Finally, store residual families and the combined functional norm:
    dGEMMassBalanceNorm = dMassNormSquared**(0.5)
    dGEMSolutionChemicalPotentialNorm = dSolutionChemicalNormSquared**(0.5)
    dGEMCondensedChemicalPotentialNorm = dCondensedChemicalNormSquared**(0.5)
    dGEMChemicalPotentialNorm = (dSolutionChemicalNormSquared + &
        dExchangeNormLocal**2 + dCondensedChemicalNormSquared)**(0.5)
    dGEMFunctionNorm = (dMassNormSquared + dGEMChemicalPotentialNorm**2)**(0.5)

    if (lDebugMode) print *, 'dGEMFunctionNorm = ', dGEMFunctionNorm

    if (allocated(dPhaseComp)) deallocate(dPhaseComp)

    return

contains

    integer function ActiveSlotThermoPhaseForNorm(iSlotIn, iDisplayPhaseIn)

        implicit none

        integer, intent(in) :: iSlotIn, iDisplayPhaseIn

        ActiveSlotThermoPhaseForNorm = iDisplayPhaseIn
        if (.NOT.allocated(iActiveSlotThermoPhase)) return
        if ((iSlotIn <= 0).OR.(iSlotIn > SIZE(iActiveSlotThermoPhase))) return
        if (iActiveSlotThermoPhase(iSlotIn) > 0) ActiveSlotThermoPhaseForNorm = iActiveSlotThermoPhase(iSlotIn)

        return

    end function ActiveSlotThermoPhaseForNorm

    logical function IsNormCEFPhase(iSolnPhaseIndexIn)

        implicit none

        integer, intent(in) :: iSolnPhaseIndexIn

        IsNormCEFPhase = .FALSE.
        if (iSolnPhaseIndexIn <= 0) return
        if (.NOT.allocated(iPhaseSublattice)) return
        if (iSolnPhaseIndexIn > SIZE(iPhaseSublattice)) return
        if (iPhaseSublattice(iSolnPhaseIndexIn) <= 0) return

        IsNormCEFPhase = (TRIM(cSolnPhaseType(iSolnPhaseIndexIn)) == 'SUBL').OR.&
            (TRIM(cSolnPhaseType(iSolnPhaseIndexIn)) == 'SUBLM').OR.&
            (TRIM(cSolnPhaseType(iSolnPhaseIndexIn)) == 'SUBOM')

        return

    end function IsNormCEFPhase

    subroutine ComputeActiveSlotCEFPhaseComposition(iSolnPhaseIndexIn, iPhaseIDIn, iSlotIn, dCompositionOut)

        implicit none

        integer, intent(in) :: iSolnPhaseIndexIn, iPhaseIDIn, iSlotIn
        real(8), intent(out) :: dCompositionOut(nElements)

        real(8), allocatable :: dSiteLocal(:,:)

        allocate(dSiteLocal(nMaxSublatticeSys,nMaxConstituentSys))
        if (allocated(dActiveSlotSiteFraction)) then
            if ((iSlotIn > 0).AND.(iSlotIn <= SIZE(dActiveSlotSiteFraction,1)).AND.&
                (SUM(dActiveSlotSiteFraction(iSlotIn,:,:)) > 0D0)) then
                dSiteLocal = dActiveSlotSiteFraction(iSlotIn,:,:)
            else
                call BuildCEFNormSiteFromMolFraction(iSolnPhaseIndexIn, iPhaseIDIn, dSiteLocal)
            end if
        else
            call BuildCEFNormSiteFromMolFraction(iSolnPhaseIndexIn, iPhaseIDIn, dSiteLocal)
        end if
        call ComputeCEFNormPhaseComposition(iSolnPhaseIndexIn, iPhaseIDIn, dSiteLocal, dCompositionOut)
        if (allocated(dSiteLocal)) deallocate(dSiteLocal)

        return

    end subroutine ComputeActiveSlotCEFPhaseComposition

    subroutine BuildCEFNormSiteFromMolFraction(iSolnPhaseIndexIn, iPhaseIDIn, dSiteOut)

        implicit none

        integer, intent(in) :: iSolnPhaseIndexIn, iPhaseIDIn
        real(8), intent(out) :: dSiteOut(nMaxSublatticeSys,nMaxConstituentSys)

        integer :: iFirstLocal, iLastLocal, iLocalSpecies, iLocal, iSub, iCon
        real(8) :: dSumLocal

        dSiteOut = 0D0
        iFirstLocal = nSpeciesPhase(iSolnPhaseIndexIn-1) + 1
        iLastLocal = nSpeciesPhase(iSolnPhaseIndexIn)
        do iLocalSpecies = iFirstLocal, iLastLocal
            iLocal = iLocalSpecies - iFirstLocal + 1
            do iSub = 1, nSublatticePhase(iPhaseIDIn)
                iCon = iConstituentSublattice(iPhaseIDIn,iSub,iLocal)
                if (iCon > 0) dSiteOut(iSub,iCon) = dSiteOut(iSub,iCon) + &
                    DMAX1(dMolFraction(iLocalSpecies), 0D0)
            end do
        end do

        do iSub = 1, nSublatticePhase(iPhaseIDIn)
            dSumLocal = SUM(dSiteOut(iSub,1:nConstituentSublattice(iPhaseIDIn,iSub)))
            if (dSumLocal > 0D0) then
                dSiteOut(iSub,1:nConstituentSublattice(iPhaseIDIn,iSub)) = &
                    dSiteOut(iSub,1:nConstituentSublattice(iPhaseIDIn,iSub)) / dSumLocal
            end if
        end do

        return

    end subroutine BuildCEFNormSiteFromMolFraction

    subroutine ComputeCEFNormPhaseComposition(iSolnPhaseIndexIn, iPhaseIDIn, dSiteIn, dCompositionOut)

        implicit none

        integer, intent(in) :: iSolnPhaseIndexIn, iPhaseIDIn
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)
        real(8), intent(out) :: dCompositionOut(nElements)

        integer :: iFirstLocal, iLastLocal, iLocalSpecies, iLocal, iSub, iCon
        real(8) :: dProduct, dProductSum
        real(8), allocatable :: dProductSpecies(:)

        dCompositionOut = 0D0
        iFirstLocal = nSpeciesPhase(iSolnPhaseIndexIn-1) + 1
        iLastLocal = nSpeciesPhase(iSolnPhaseIndexIn)
        allocate(dProductSpecies(iLastLocal-iFirstLocal+1))
        dProductSpecies = 0D0
        dProductSum = 0D0

        do iLocalSpecies = iFirstLocal, iLastLocal
            iLocal = iLocalSpecies - iFirstLocal + 1
            dProduct = 1D0
            do iSub = 1, nSublatticePhase(iPhaseIDIn)
                iCon = iConstituentSublattice(iPhaseIDIn,iSub,iLocal)
                if (iCon > 0) dProduct = dProduct * DMAX1(dSiteIn(iSub,iCon), 1D-75)
            end do
            dProductSpecies(iLocal) = dProduct
            dProductSum = dProductSum + dProduct
        end do

        if (dProductSum > 0D0) then
            do iLocalSpecies = iFirstLocal, iLastLocal
                iLocal = iLocalSpecies - iFirstLocal + 1
                dCompositionOut = dCompositionOut + &
                    dProductSpecies(iLocal) * dStoichSpecies(iLocalSpecies,1:nElements) / dProductSum
            end do
        end if

        if (allocated(dProductSpecies)) deallocate(dProductSpecies)

        return

    end subroutine ComputeCEFNormPhaseComposition

end subroutine CompFunctionNorm
