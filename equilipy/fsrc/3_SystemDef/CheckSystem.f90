!
!
subroutine CheckSystem
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CheckSystem.f90
    !> \brief   Check for consistency between the system and the data-file.
    !> \author  M.H.A. Piro
    !> \date    Apr. 24, 2012
    !> \sa      Thermochimica.f90
    !> \sa      GetElementName.f90
    !> \todo    Do a better job allocating variables iPairID and dCoordinationNumber.
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   10/17/2011      M.H.A. Piro         Original code
    !   05/01/2012      M.H.A. Piro         Corrected storage of mixing component indices when solution species
    !                                       are no longer considered.  This affected the variable iSpeciesPass.
    !                                       This also affected the way that pure condensed phases are represented
    !                                       by iSpeciesPass.
    !   02/17/2012      M.H.A. Piro         Check if allocatable arrays have been allocated and if so, have they
    !                                       changed in dimension.
    !   05/24/2012      M.H.A. Piro         Correct the minimum number of moles of an element allowed =
    !                                        (normalizer) * (machine precision) / (mass balance tolerance)
    !   08/21/2012      M.H.A. Piro         Redefine the tolerance for the minimum number of moles of a solution
    !                                        phase that is introduced to the system by the minimum number of
    !                                        moles of an element in the system.
    !   01/14/2013      M.H.A. Piro         Move code relavent to excess terms into a separate subroutine.
    !   02/05/2013      M.H.A. Piro         The check for the minimum number of system components now includes
    !                                        constraints imposed by charge neutrality.
    !   05/14/2013      M.H.A. Piro         Fixed bug in allocating arrays specific to ionic phases.  This was
    !                                        done when the size of cSolnPhaseName changes. This should be independent.
    !   11/03/2021      S.Y. Kwon           Removed Compound variables
    !   06/26/2026      S.Y. Kwon           Added dependent-element screening for pseudo-component bases.
    !   06/28/2026      S.Y. Kwon           Moved pseudo-component dependent-element detection from input labels
    !                                       to rank analysis of the screened species stoichiometry.
    !   06/28/2026      S.Y. Kwon           Honored phase selection when screening pure condensed species.
!   06/28/2026      S.Y. Kwon           Identified pseudo-component dependent elements from the coupled species
!                                       submatrix and excluded species outside that active-bundle basis.
!   07/02/2026      S.Y. Kwon           Added a separate order/disorder active-companion map for minimizer
!                                       identity diagnostics.
!   07/03/2026      S.Y. Kwon           Used the shared order/disorder helper-alias classifier when resolving
!                                       active companion aliases.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to ensure that the selection of system components in the
    !! parsed ChemSage data-file and the data provided to Thermochimica are consistent.  System components are
    !! always taken to be chemical elements in Thermochimica, or "elements" for sake of brevity.  An element will
    !! only be considered if thermodynamic data is provided by the data-file and if the mass of that particular
    !! element is provided.  If the mass of a particular element is provided but there isn't any data for it
    !! (via the data-file), that element will not be considered.  Similarly, if thermodynamic data (via the
    !! data-file) is provided for a particular element, but the mass of that element is not available, that
    !! element will not be considered.
    !!
    !! This subroutine will select all species and phases that are relevant to this system from the data
    !! parsed from the ChemSage data-file.
    !!
    !! The variable cInputThermo(3) can accept the following values:
    !!
    ! cThermoInputUnits(3)      Description
    ! --------------------      -----------
    !> \details
    !!
    !! <table border="1" width="800">
    !! <tr>
    !!    <td> <b> Units </td> <td> Description </b> </td>
    !! </tr>
    !! <tr>
    !!    <td> "mass fraction" </td>
    !!    <td> Mass fraction in dimensionless units (e.g., gram/gram, kilogram/kilogram, pound/pound, wt%).  </td>
    !! </tr>
    !! <tr>
    !!    <td> "mole fraction" </td>
    !!    <td> Mole fraction in dimensionless units (e.g., mole/mole, mol%).  </td>
    !! </tr>
    !! <tr>
    !!    <td> "atom fraction" </td>
    !!    <td> Atom fraction in dimensionless units (e.g., atom/atom, at%).  </td>
    !! </tr>
    !! <tr>
    !!    <td> "kilograms" </td>
    !!    <td> All quantities are in kilograms.  </td>
    !! </tr>
    !! <tr>
    !!    <td> "grams" </td>
    !!    <td> All quantities are in grams.  </td>
    !! </tr>
    !! <tr>
    !!    <td> "pounds" or "lbs" </td>
    !!    <td> All quantities are in pounds.  </td>
    !! </tr>
    !! <tr>
    !!    <td> "moles" </td>
    !!    <td> All quantities are in moles.  </td>
    !! </tr>
    !! <tr>
    !!    <td> "gram-atoms" </td>
    !!    <td> All quantities are in gram-atoms (same as moles for the pure elements);  </td>
    !! </tr>
    !! <tr>
    !!    <td> "atoms"   </td>
    !!    <td>  All quantities are in atoms.  </td>
    !! </tr>
    !! </table>
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! dNormalizeInput                   A double scalar that normalizes the mass units.  This is performed to
    !                                    minimize numerical error when evaluating the Jacobian/Broyden matrix in
    !                                    the GEMSolver.
    !
    ! dTemperature                      Temperature [K]
    ! dPressure                         Absolute hydrostatic pressure [atm]
    ! dElementMass                      Total mass of each element, where the coefficient corresponds to the
    !                                    atomic number (e.g., dMolesElement(92) refers to uranium).
    ! dElementMoleFractionMin           A minimum allowable value for the number of moles of an element.
    ! INFOThermo                        A scalar integer that indicates a successful exit or identifies an error.
    ! cInputThermo                      A character vector containing the units for temperature, pressure and
    !                                    elemental quantity.
    !
    !-------------------------------------------------------------------------------------------------------------
!
    USE ModuleParseCS
    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver, ONLY: lSolnPhases, lMiscibility
!
    implicit none
!
    integer                                 :: i, j, k, m, n, o, nMaxSpeciesPhase, nCountSublatticeTemp, iCon1, iCon2, iCon3, iCon4
    integer                                 :: mm, nn, nConstituentPass, c, s, nSpeciesCurrentPhase
    integer,parameter                       :: iDependentElement = -2
    integer,dimension(0:nSolnPhasesSysCS+1) :: iTempVec
    integer,dimension(nSolnPhasesSysCS)     :: iSolnPhaseCS2Sys
    real(8)                                 :: dSum, dElementMoleFractionMin
    character(3),dimension(0:nElementsPT)   :: cElementNamePT
    character(12)                           :: cDummy
    logical                                 :: lPos, lNeg, lConsiderEndmember, lActiveStoich, lDependentStoich
    logical                                 :: lHasDependentElement
!
!
    ! Check to see if the allocatable arrays have already been allocated
    if (allocated(iElementSystem)) then
        ! Do nothing.
    else
        ! Allocate memory:
        allocate(iElementSystem(1:nElementsCS))
        allocate(iSpeciesPass(nSpeciesCS))
    end if
!
    ! Check if there are any solution phases with a sublattice:
    if (nCountSublatticeCS > 0) then
        ! Allocate array to check if a constituent passes:
        j = SIZE(iConstituentSublatticeCS, DIM=3)
        if (allocated(iConstituentPass)) deallocate(iConstituentPass)
        allocate(iConstituentPass(nCountSublatticeCS,nMaxSublatticeCS,j))
        iConstituentPass = 0
    end if
!
    ! Initialize variables:
    iElementSystem      = 0
    iTempVec            = 0
    iSolnPhaseCS2Sys    = 0
    iSpeciesPass        = 0
    nCountSublattice    = 0
    nMaxSublatticeSys   = 0
    nMaxConstituentSys  = 0
    nMaxSpeciesPhase    = 0
    nCountSublatticeTemp   = 0
    nChargedConstraints = 0
    n                   = 0
    dSum                = 0D0
    dElementMoleFractionMin = dTolerance(6)
!
    ! Get the name of all the elements on the periodic table:
    call GetElementName(cElementNamePT)
!
    ! Perform a mass conversion:
    LOOP_Small: do j = 1, nElementsCS
        ! Map the character string of this element to its atomic number:
        LOOP_Big: do i = 0, nElementsPT
            if (cElementNameCS(j) == cElementNamePT(i)) then
                if (iDependentElementInput(i) /= 0) then
                    iElementSystem(j) = iDependentElement
                    dElementMass(i)   = 0D0
                    cycle LOOP_Small
                end if
                ! Convert the mass of each element to moles:
                select case (cInputUnitMass)
                case ('grams','g','mass fraction','weight fraction','wt%','wt.%')
                    ! Convert mass unit to moles:
                    dElementMass(i) = dElementMass(i) / dAtomicMassCS(j)
                case ('kilograms','kg')
                    ! Convert mass unit to moles:
                    dElementMass(i) = 1000*dElementMass(i) / dAtomicMassCS(j)
                case ('pounds','lbs')
                    dElementMass(i) = 453.59237*dElementMass(i) / dAtomicMassCS(j)
                case ('moles','mol','atoms','gram-atoms','gram-moles',&
                    'mole fraction','atom fraction','mol%','at%')
                    ! Do nothing
                case default
                    ! The character string representing input units is not recognized.
                    INFOThermo = 4
                    return
                end select
                iElementSystem(j) = i
                dSum = dSum + dElementMass(i)
                cycle LOOP_Small
            elseif (cElementNameCS(j) == 'e-') then
                ! Electron
                iElementSystem(j) = i
                dElementMass(i)   = 0D0
                cycle LOOP_Small
            end if
        end do LOOP_Big
        ! The chemical element stored by cElementNameCS(j) does not correspond to an
        ! element in cElementNamePT.  Report an error and exit:
        INFOThermo = 31
        return
    end do LOOP_Small
!
    cInputUnitMass = 'moles'
!
!
    nElemOrComp = nElementsCS
    ! call CheckSysComponents
    ! call CheckElements
    ! !SY: Compound variables that considers user defined component (e.g.CaO) as system elements, are removed
    ! !This part will be solve automatically by Lagrangian. Only calculate element for now
    ! if (nElemOrComp /= nElementsCS) then
    !     ! ! call CheckCompounds
    !     ! nElemOrComp = nCompounds
    !     dSum = 0
    !     do i = 1, nElemOrComp
    !         dSum = dSum + dElementMass(iElementSystem(i))
    !     end do
    ! end if
!
    ! ! Make sure that the sum of element masses is not zero:
    ! if (dSum == 0D0) then
    !     INFOThermo = 5
    !     return
    ! end if
!
    ! Make dSum multiplicative and normalize:
    dNormalizeInput = dNormalizeInput / dSum
!
    k = 0
    nElements = 0
    ! Normalize dElementMass to mole fraction and establish the elements of the system:
    do j = 1, nElemOrComp
        if (iElementSystem(j) == iDependentElement) cycle
        dElementMass(iElementSystem(j)) = dElementMass(iElementSystem(j)) * dNormalizeInput
        if (dElementMass(iElementSystem(j)) < dElementMoleFractionMin) then
            ! Element j should not be considered.
            dElementMass(iElementSystem(j)) = 0
            iElementSystem(j) = 0
            if (cElementNameCS(j) == 'e-') then
                nElements         = nElements + 1
                iElementSystem(j) = -1
                k = k + 1
            end if
        else
            nElements = nElements + 1
        end if
    end do

    call DetectRankDependentElements(iDependentElement)

    nElements = 0
    k = 0
    lHasDependentElement = .FALSE.
    do j = 1, nElemOrComp
        if (iElementSystem(j) == iDependentElement) cycle
        if ((iElementSystem(j) > 0).OR.(iElementSystem(j) == -1)) then
            nElements = nElements + 1
            if (iElementSystem(j) == -1) k = k + 1
        end if
    end do
    do j = 1, nElemOrComp
        if (iElementSystem(j) == iDependentElement) lHasDependentElement = .TRUE.
    end do



    LOOP_checkElements: do i = 1, nElementsPT
        if (dElementMass(i) > 0) then
            do j = 1, nElementsCS
                if (iElementSystem(j) == i) cycle LOOP_checkElements
            end do
            print *, "WARNING: Element ", cElementNamePT(i), " not in database and therefore omitted from calculation"
        end if
    end do LOOP_checkElements

    ! The system requires a minimum of two elements in order to be considered:
    ! SY: when nElements==1, no result is given. Add a function that print and terminate the calculations in CheckThermoData
    if (nElements < 1 + k) then
        INFOThermo = 5
        return
    end if
    
    ! Check to see if the system has to be re-adjusted:
    ! The system is smaller than what is in the data-file.  Re-compute the system parameters.
    ! Loop through solution phases:
    LOOP_SolnPhases: do o = 1, nSolnPhasesSysCS
        if(allocated(iSolnPS)) then
            if (iSolnPS(o)>0) then
                i=o
            else
                if ((cSolnPhaseTypeCS(o) == 'SUBL').OR.(cSolnPhaseTypeCS(o) == 'SUBLM').OR.&
                    (cSolnPhaseTypeCS(o) == 'SUBOM').OR.&
                    (cSolnPhaseTypeCS(o) == 'SUBG').OR.(cSolnPhaseTypeCS(o) == 'SUBQ')) then
                    nCountSublatticeTemp = nCountSublatticeTemp + 1
                end if
                cycle
            end if
        else
            i=o
        end if
        ! Loop through species in solution phases:
        nSpeciesCurrentPhase = 0
        LOOP_SpeciesInSolnPhase: do j = nSpeciesPhaseCS(i) + 1, nSpeciesPhaseCS(i+1)
            lActiveStoich = .FALSE.
            lDependentStoich = .FALSE.
            do k = 1, nElemOrComp
                if (((iElementSystem(k) > 0).OR.(iElementSystem(k) == -1)).AND.&
                    (DABS(dStoichSpeciesCS(j,k)) > 0D0)) lActiveStoich = .TRUE.
                if ((iElementSystem(k) == iDependentElement).AND.&
                    (DABS(dStoichSpeciesCS(j,k)) > 0D0)) lDependentStoich = .TRUE.
                if ((dStoichSpeciesCS(j,k) > 0).AND.(iElementSystem(k) == 0)) then
                    ! This species should not be considered
                    cycle LOOP_SpeciesInSolnPhase
                end if
            end do
            if ((.NOT.lActiveStoich).AND.lDependentStoich) cycle LOOP_SpeciesInSolnPhase
            if (lHasDependentElement.AND.lActiveStoich.AND.(.NOT.lDependentStoich)) cycle LOOP_SpeciesInSolnPhase
            nSpecies = nSpecies + 1
            nSpeciesCurrentPhase = nSpeciesCurrentPhase + 1
            iSpeciesPass(j) = nSpeciesCurrentPhase
            if (cSolnPhaseTypeCS(i) == 'SUBG' .OR. cSolnPhaseTypeCS(i) == 'SUBQ') then
                k = iPhaseSublatticeCS(i)
                iCon1 = iPairIDCS(k,j-nSpeciesPhaseCS(i),1)
                iCon2 = iPairIDCS(k,j-nSpeciesPhaseCS(i),2)
                iCon3 = iPairIDCS(k,j-nSpeciesPhaseCS(i),3) - nSublatticeElementsCS(k,1)
                iCon4 = iPairIDCS(k,j-nSpeciesPhaseCS(i),4) - nSublatticeElementsCS(k,1)
                if (iCon1 > 0) iConstituentPass(k,1,iCon1) = 1
                if (iCon2 > 0) iConstituentPass(k,1,iCon2) = 1
                if (iCon3 > 0) iConstituentPass(k,2,iCon3) = 1
                if (iCon4 > 0) iConstituentPass(k,2,iCon4) = 1
            end if
        end do LOOP_SpeciesInSolnPhase

        ! For electrons, there have to be species with positive and negative stoichiometries still in the system.
        ! Otherwise, remove the species that use that electron.
        do k = 1, nElementsCS
            if (cElementNameCS(k) == 'e-') then
                lPos = .FALSE.
                lNeg = .FALSE.
                do j = nSpeciesPhaseCS(i) + 1, nSpeciesPhaseCS(i+1)
                    if ((iSpeciesPass(j) > 0) .AND. dStoichSpeciesCS(j,k) > 0D0) lPos = .TRUE.
                    if ((iSpeciesPass(j) > 0) .AND. dStoichSpeciesCS(j,k) < 0D0) lNeg = .TRUE.
                end do
                if (lPos .NEQV. lNeg) then
                    do j = nSpeciesPhaseCS(i) + 1, nSpeciesPhaseCS(i+1)
                        if ((iSpeciesPass(j) > 0) .AND. DABS(dStoichSpeciesCS(j,k)) > 0D0) then
                            iSpeciesPass(j) = 0
                            nSpecies = nSpecies - 1
                            nSpeciesCurrentPhase  = nSpeciesCurrentPhase-1
                        end if
                    end do
                end if
            end if
        end do

        ! Store temporary counter for the number of charged phases from the CS data-file:
        if ((cSolnPhaseTypeCS(i) == 'SUBL').OR.(cSolnPhaseTypeCS(i) == 'SUBLM').OR. &
            (cSolnPhaseTypeCS(i) == 'SUBOM')) then
            nCountSublatticeTemp = nCountSublatticeTemp + 1
            if (nSpeciesCurrentPhase > 0) then
                k = iPhaseSublatticeCS(i)
                ! Loop through species in phase to determine which constituents are stable:
                do j = nSpeciesPhaseCS(i) + 1, nSpeciesPhaseCS(i+1)
                    if (iSpeciesPass(j) > 0) then
                        ! This species has passed.
                        nn = j - nSpeciesPhaseCS(i)
                        ! Loop through sublattices per phase:
                        do s = 1, nSublatticePhaseCS(k)
                            mm = iConstituentSublatticeCS(k, s, nn)
                            iConstituentPass(k, s, mm) = 1
                        end do
                    end if
                end do
                ! If a sublattice has no constituents the whole phase should get removed
                LOOP_CHECK_SUBLATTICES: do s = 1, nSublatticePhaseCS(k)
                    nConstituentPass = 0
                    do c = 1, nConstituentSublatticeCS(k,s)
                        if (iConstituentPass(k,s,c) /= 0) then
                            nConstituentPass = nConstituentPass + 1
                        end if
                    end do
                    if (nConstituentPass < 1) then
                        iConstituentPass(k,:,:) = 0
                        nSpecies = nSpecies - nSpeciesCurrentPhase
                        iSpeciesPass(nSpeciesPhaseCS(i) + 1:nSpeciesPhaseCS(i+1)) = 0
                        exit LOOP_CHECK_SUBLATTICES
                    end if
                end do LOOP_CHECK_SUBLATTICES
            end if
                
            
        else if((cSolnPhaseTypeCS(i) == 'SUBG').OR.(cSolnPhaseTypeCS(i) == 'SUBQ')) then
            nCountSublatticeTemp = nCountSublatticeTemp + 1

        end if
!
        ! Count the number of solution phases in the system:
        iTempVec(nSolnPhasesSys+1) = nSpecies
        if (lEndmembers2Species) then
            lConsiderEndmember=iTempVec(nSolnPhasesSys+1) - iTempVec(nSolnPhasesSys) >= 1
        else
            lConsiderEndmember=iTempVec(nSolnPhasesSys+1) - iTempVec(nSolnPhasesSys) > 1
        end if
        
        if (lConsiderEndmember) then
            nSolnPhasesSys = nSolnPhasesSys + 1
            iSolnPhaseCS2Sys(i) = nSolnPhasesSys
            nMaxSpeciesPhase = MAX(nMaxSpeciesPhase, iTempVec(nSolnPhasesSys) - iTempVec(nSolnPhasesSys-1))

            ! Check if this is a charged phase:
            if ((cSolnPhaseTypeCS(i) == 'SUBL').OR.(cSolnPhaseTypeCS(i) == 'SUBLM').OR. &
                    (cSolnPhaseTypeCS(i) == 'SUBOM').OR. &
                    (cSolnPhaseTypeCS(i) == 'SUBG').OR.(cSolnPhaseTypeCS(i) == 'SUBQ')) then
                ! Count the number of charged phases:
                nCountSublattice = nCountSublattice + 1
                ! Determine the maximum number of sublattice of any stable phase:
                nMaxSublatticeSys  = MAX(nMaxSublatticeSys,nSublatticePhaseCS(nCountSublatticeTemp))
                m = MAXVAL(nConstituentSublatticeCS(nCountSublatticeTemp,1:nMaxSublatticeCS))
                nMaxConstituentSys = MAX(nMaxConstituentSys,m)
            end if
            
        elseif (iTempVec(nSolnPhasesSys+1) - iTempVec(nSolnPhasesSys) == 1) then
            ! There is only one species in this solution phase.  This solution phase should not be considered.
            iSpeciesPass(nSpeciesPhaseCS(i) + 1 : nSpeciesPhaseCS(i+1))          = 0
            nSpecies                 = nSpecies - 1
            iTempVec(nSolnPhasesSys) = nSpecies
        else
            ! Do nothing.  The number of species in this solution phase is zero and this phase will not be
            ! considered.
        end if
    end do LOOP_SolnPhases
!
    ! Loop through pure condensed phases:
    LOOP_PureConPhases: do j = nSpeciesPhaseCS(nSolnPhasesSysCS+1) + 1, nSpeciesCS
        if (iPhaseCS(j) /= 0) cycle LOOP_PureConPhases
        if (SUM(dStoichSpeciesCS(j,1:nElemOrComp)) == 0) cycle LOOP_PureConPhases
        lActiveStoich = .FALSE.
        lDependentStoich = .FALSE.
        do k = 1, nElemOrComp
            if (((iElementSystem(k) > 0).OR.(iElementSystem(k) == -1)).AND.&
                (DABS(dStoichSpeciesCS(j,k)) > 0D0)) lActiveStoich = .TRUE.
            if ((iElementSystem(k) == iDependentElement).AND.&
                (DABS(dStoichSpeciesCS(j,k)) > 0D0)) lDependentStoich = .TRUE.
            if ((dStoichSpeciesCS(j,k) > 0).AND.(iElementSystem(k) == 0)) then
                ! This species should not be considered
                cycle LOOP_PureConPhases
            end if
        end do
        if ((.NOT.lActiveStoich).AND.lDependentStoich) cycle LOOP_PureConPhases
        if (lHasDependentElement.AND.lActiveStoich.AND.(.NOT.lDependentStoich)) cycle LOOP_PureConPhases
        nSpecies        = nSpecies + 1
        iSpeciesPass(j) = 1
    end do LOOP_PureConPhases
        
    ! Re-establish the character vector representing the element names:
    j = 0
    do i = 1, nElemOrComp
        if ((iElementSystem(i) > 0).OR.(iElementSystem(i) == -1)) then
            cDummy           = cElementNameCS(i)
            if (cDummy == 'e-') nChargedConstraints = nChargedConstraints + 1
        end if
    end do
!
    ! Add dummy species representing electrons to the number of species in the system:
    nSpecies = nSpecies + nChargedConstraints
!
    ! Check if these variables have already been allocated:
    if (allocated(dChemicalPotential)) then
        ! Check to see if the number of species has changed:
!
        i = SIZE(dChemicalPotential)
        if (i /= nSpecies) then
            ! The number of species has changed.
            deallocate(dChemicalPotential,dLevelingChemicalPotential,dActivity,iPhase,&
                dSpeciesTotalAtoms,dLevelingSpeciesTotalAtoms,dLevelingSpeciesFormulaAtoms,cSpeciesName,&
                iParticlesPerMole,dCoeffGibbsMagnetic,dMagGibbsEnergy,dMagEnthalpy,&
                dMagEntropy,dMagHeatCapacity,&
                dPartialEnthalpy,dPartialEntropy,dPartialHeatCapacity,&
                dStdGibbsEnergy,dStdEnthalpy,dStdEntropy,dStdHeatCapacity, STAT = n)
            if (n /= 0) then
                INFOThermo = 19
                return
            end if
            ! Allocate memory for variables:
            allocate(dChemicalPotential(nSpecies),dLevelingChemicalPotential(nSpecies),&
                dActivity(nSpecies),iPhase(nSpecies),dSpeciesTotalAtoms(nSpecies),&
                dLevelingSpeciesTotalAtoms(nSpecies),dLevelingSpeciesFormulaAtoms(nSpecies))
            allocate(cSpeciesName(nSpecies),dStdGibbsEnergy(nSpecies),dStdEnthalpy(nSpecies))
            allocate(dStdEntropy(nSpecies),dStdHeatCapacity(nSpecies))
            allocate(iParticlesPerMole(nSpecies), dCoeffGibbsMagnetic(nSpecies,4), dMagGibbsEnergy(nSpecies))
            allocate(dMagEnthalpy(nSpecies),dMagEntropy(nSpecies),dMagHeatCapacity(nSpecies))
            allocate(dPartialEnthalpy(nSpecies),dPartialEntropy(nSpecies),dPartialHeatCapacity(nSpecies))
        end if
!
        ! Check to see if the number of elements has changed:
        j = SIZE(cElementName)
        if (j /= nElements) then
            ! The number of elements has changed.
            deallocate(cElementName,dMolesElement,dAtomicMass, STAT = n)
            if (n /= 0) then
                INFOThermo = 19
                return
            end if
            ! Allocate memory:
            allocate(cElementName(nElements),dMolesElement(nElements),dAtomicMass(nElements))
        end if
!
        ! Check to see if either the number of species or the number of elements has changed:
        if ((i /= nSpecies).OR.(j /= nElements)) then
            deallocate(dAtomFractionSpecies,dLevelingCompositionSpecies,dStoichSpecies, STAT = n)
            if (n /= 0) then
                INFOThermo = 19
                return
            end if
            ! Allocate memory:
            allocate(dAtomFractionSpecies(nSpecies,nElements),&
                dLevelingCompositionSpecies(nSpecies,nElements),dStoichSpecies(nSpecies,nElements))
        end if
!
        ! Check to see if the number of solution phases in the system has changed:
        k = SIZE(cSolnPhaseName)
        if (k /= nSolnPhasesSys) then
            ! The number of solution phases in the system has changed.
            deallocate(nSpeciesPhase,nParamPhase,cSolnPhaseType,cSolnPhaseName,lSolnPhases,dGibbsSolnPhase, &
                iDisorderedPhase, iODCompanionPhase, &
                lMiscibility,nMagParamPhase, STAT = n)
            if (n /= 0) then
                INFOThermo = 19
                return
            end if
            ! Allocate memory:
            allocate(nSpeciesPhase(0:nSolnPhasesSys),nParamPhase(0:nSolnPhasesSys),nMagParamPhase(0:nSolnPhasesSys))
            allocate(cSolnPhaseType(nSolnPhasesSys),cSolnPhaseName(nSolnPhasesSys))
            allocate(lSolnPhases(nSolnPhasesSys),dGibbsSolnPhase(nSolnPhasesSys),lMiscibility(nSolnPhasesSys))
            allocate(iDisorderedPhase(nSolnPhasesSys), iODCompanionPhase(nSolnPhasesSys))
!
        end if
!
        ! Only allocate if there are sublattice phases:
        if (nCountSublattice > 0) then
!
            deallocate(iPhaseSublattice,nSublatticePhase,nConstituentSublattice,dStoichSublattice, &
                    dSiteFraction,cConstituentNameSUB,iConstituentSublattice,nSublatticeElements, &
                     nPairsSRO,iPairID,dCoordinationNumber, dZetaSpecies, &
                     dSublatticeCharge,iChemicalGroup,dStoichPairs,cPairName,dConstituentCoefficients, STAT = n)
!
            allocate(iPhaseSublattice(nSolnPhasesSys),nSublatticePhase(nCountSublattice))
            allocate(nConstituentSublattice(nCountSublattice,nMaxSublatticeSys))
            allocate(dStoichSublattice(nCountSublattice,nMaxSublatticeSys))
            allocate(dSiteFraction(nCountSublattice,nMaxSublatticeSys,nMaxConstituentSys))
            allocate(cConstituentNameSUB(nCountSublattice,nMaxSublatticeSys,MAX(nMaxConstituentSys,nMaxSpeciesPhase)))
            allocate(iConstituentSublattice(nCountSublattice,nMaxSublatticeSys,nMaxSpeciesPhase))
            allocate(nSublatticeElements(nCountSublattice,nMaxSublatticeSys))
            j = MAXVAL(nSublatticeElementsCS)
            allocate(nPairsSRO(nCountSublattice,2))
            allocate(iPairID(nCountSublattice,nMaxSpeciesPhase,4))
            allocate(dCoordinationNumber(nCountSublattice,nMaxSpeciesPhase,4))
            allocate(dZetaSpecies(nCountSublattice,nMaxSpeciesPhase))
            allocate(dConstituentCoefficients(nCountSublattice,nMaxSpeciesPhase,5))
            allocate(dSublatticeCharge(nCountSublattice,nMaxSublatticeSys,j))
            allocate(iChemicalGroup(nCountSublattice,nMaxSublatticeSys,j))
            allocate(dStoichPairs(nCountSublattice,MAXVAL(nPairsSROCS(:,1)),nElements))
            allocate(cPairName(nCountSublattice,MAXVAL(nPairsSROCS(:,1))))
        end if
!
    else
!
        ! Allocate memory for variables:
        allocate(dChemicalPotential(nSpecies),dLevelingChemicalPotential(nSpecies),&
            dActivity(nSpecies),iPhase(nSpecies),dSpeciesTotalAtoms(nSpecies),&
            dLevelingSpeciesTotalAtoms(nSpecies),dLevelingSpeciesFormulaAtoms(nSpecies))
        allocate(dPartialEnthalpy(nSpecies),dPartialEntropy(nSpecies),dPartialHeatCapacity(nSpecies))
        allocate(cSpeciesName(nSpecies),dStdGibbsEnergy(nSpecies),dStdEnthalpy(nSpecies))
        allocate(dStdEntropy(nSpecies),dStdHeatCapacity(nSpecies))
        allocate(iParticlesPerMole(nSpecies),dCoeffGibbsMagnetic(nSpecies,4),dMagGibbsEnergy(nSpecies))
        allocate(dMagEnthalpy(nSpecies),dMagEntropy(nSpecies),dMagHeatCapacity(nSpecies))
        allocate(cElementName(nElements),dMolesElement(nElements),dAtomicMass(nElements))
        allocate(dAtomFractionSpecies(nSpecies,nElements),&
            dLevelingCompositionSpecies(nSpecies,nElements),dStoichSpecies(nSpecies,nElements))
        allocate(nSpeciesPhase(0:nSolnPhasesSys),nParamPhase(0:nSolnPhasesSys),nMagParamPhase(0:nSolnPhasesSys))
        allocate(cSolnPhaseType(nSolnPhasesSys),cSolnPhaseName(nSolnPhasesSys))
        allocate(lSolnPhases(nSolnPhasesSys),dGibbsSolnPhase(nSolnPhasesSys),lMiscibility(nSolnPhasesSys))
        allocate(iDisorderedPhase(nSolnPhasesSys), iODCompanionPhase(nSolnPhasesSys))
!
        ! Only allocate if there are charged phases:
        if (nCountSublattice > 0) then
            allocate(iPhaseSublattice(nSolnPhasesSys),nSublatticePhase(nCountSublattice))
            allocate(nConstituentSublattice(nCountSublattice,nMaxSublatticeSys))
            allocate(dStoichSublattice(nCountSublattice,nMaxSublatticeSys))
            allocate(dSiteFraction(nCountSublattice,nMaxSublatticeSys,nMaxConstituentSys))
            allocate(cConstituentNameSUB(nCountSublattice,nMaxSublatticeSys,MAX(nMaxConstituentSys,nMaxSpeciesPhase)))
            allocate(iConstituentSublattice(nCountSublattice,nMaxSublatticeSys,nMaxSpeciesPhase))
            allocate(nSublatticeElements(nCountSublattice,nMaxSublatticeSys))
            j = MAXVAL(nSublatticeElementsCS)
            allocate(nPairsSRO(nCountSublattice,2))
            allocate(iPairID(nCountSublattice,nMaxSpeciesPhase,4))
            allocate(dCoordinationNumber(nCountSublattice,nMaxSpeciesPhase,4))
            allocate(dZetaSpecies(nCountSublattice,nMaxSpeciesPhase))
            allocate(dConstituentCoefficients(nCountSublattice,nMaxSpeciesPhase,5))
            allocate(dSublatticeCharge(nCountSublattice,nMaxSublatticeSys,j))
            allocate(iChemicalGroup(nCountSublattice,nMaxSublatticeSys,j))
            allocate(dStoichPairs(nCountSublattice,MAXVAL(nPairsSROCS(:,1)),nElements))
            allocate(cPairName(nCountSublattice,MAXVAL(nPairsSROCS(:,1))))
        end if
!
    end if
!
    ! Initialize variables:
    iPhase               = 0
    dSpeciesTotalAtoms   = 0
    dLevelingSpeciesTotalAtoms = 0D0
    dLevelingSpeciesFormulaAtoms = 0D0
    iParticlesPerMole    = 0
    dStoichSpecies       = 0
    nSpeciesPhase        = 0
    nParamPhase          = 0
    nMagParamPhase       = 0
    dChemicalPotential   = 0D0
    dLevelingChemicalPotential = 0D0
    dActivity            = 0D0
    dStdGibbsEnergy      = 0D0
    dStdEnthalpy         = 0D0
    dStdEntropy          = 0D0
    dStdHeatCapacity     = 0D0
    dMolesElement        = 0D0
    dAtomFractionSpecies = 0D0
    dLevelingCompositionSpecies = 0D0
    dGibbsSolnPhase      = 0D0
    iDisorderedPhase     = 0
    iODCompanionPhase    = 0
    dCoeffGibbsMagnetic  = 0D0
    dMagGibbsEnergy      = 0D0
    dMagEnthalpy         = 0D0
    dMagEntropy          = 0D0
    dMagHeatCapacity     = 0D0
    dPartialEnthalpy     = 0D0
    dPartialEntropy      = 0D0
    dPartialHeatCapacity = 0D0
    lSolnPhases          = .FALSE.
    lMiscibility         = .FALSE.
    nElementsSys         = nElements
!
    ! Initialize arrays (if necessary) for sublattice phases:
    if (nCountSublattice > 0) then
        dSiteFraction          = 0D0
        dStoichSublattice      = 0D0
        iConstituentSublattice = 0
        nSublatticePhase       = 0
        nConstituentSublattice = 0
        nSublatticeElements  = 0
        nPairsSRO            = 0
        iPairID              = 0
        dCoordinationNumber  = 0D0
        dZetaSpecies         = 0D0
        dConstituentCoefficients = 0D0
        dSublatticeCharge    = 0D0
        iChemicalGroup       = 0
        dStoichPairs         = 0D0
    end if
!
    ! Re-establish the character vector representing the element names:
    ! SY: Compound variables are removed
    j = 0
    do i = 1, nElemOrComp
        if ((iElementSystem(i) > 0).OR.(iElementSystem(i) == -1)) then
            j = j + 1
            cElementName(j)  = cElementNameCS(i)
            cDummy           = cElementName(j)
            dAtomicMass(j)   = dAtomicMassCS(i)
            if (iElementSystem(i) > 0) dMolesElement(j) = dElementMass(iElementSystem(i))
        end if
    end do
!
    ! Redefine the tolerance for the minimum number of moles of a solution phase that is introduced to the system:
    dTolerance(9) = DMIN1(1000D0 * MINVAL(dMolesElement, MASK = dMolesElement > 0D0),dTolerance(9))
!
    ! Re-establish the nSpeciesPhase vector:
    nSpeciesPhase(0:nSolnPhasesSys) = iTempVec(0:nSolnPhasesSys)
    do i = 1, nSolnPhasesSysCS
        if ((iSolnPhaseCS2Sys(i) > 0).AND.(iDisorderedPhaseCS(i) > 0)) then
            if (iSolnPhaseCS2Sys(iDisorderedPhaseCS(i)) > 0) then
                iDisorderedPhase(iSolnPhaseCS2Sys(i)) = iSolnPhaseCS2Sys(iDisorderedPhaseCS(i))
            end if
        end if
    end do
    iODCompanionPhase = iDisorderedPhase
    call ResolveOrderDisorderCompanionAliases
!
    ! Check the excess terms:
    call CheckSystemExcess
!
    return
!
contains

    subroutine ResolveOrderDisorderCompanionAliases

        implicit none

        integer :: iOrderedPhase, iHelperPhase, iAliasPhase

        if (.NOT.allocated(iODCompanionPhase)) return
        if (.NOT.allocated(iDisorderedPhase)) return

        do iOrderedPhase = 1, nSolnPhasesSys
            if (iOrderedPhase > SIZE(iDisorderedPhase)) exit
            iHelperPhase = iDisorderedPhase(iOrderedPhase)
            if (iHelperPhase <= 0) cycle
            if (iHelperPhase > nSolnPhasesSys) cycle

            iAliasPhase = FindOrderDisorderAliasPhase(iHelperPhase, iOrderedPhase)
            if (iAliasPhase > 0) iODCompanionPhase(iOrderedPhase) = iAliasPhase
        end do

        return

    end subroutine ResolveOrderDisorderCompanionAliases


    integer function FindOrderDisorderAliasPhase(iHelperPhase, iOrderedPhase)

        implicit none

        integer, intent(in) :: iHelperPhase, iOrderedPhase
        integer             :: iCandidatePhase

        FindOrderDisorderAliasPhase = 0

        do iCandidatePhase = 1, nSolnPhasesSys
            if (iCandidatePhase == iOrderedPhase) cycle
            if (iCandidatePhase == iHelperPhase) cycle
            if (OrderDisorderHelperAliasMatch(iHelperPhase, iCandidatePhase)) then
                FindOrderDisorderAliasPhase = iCandidatePhase
                return
            end if
        end do

        return

    end function FindOrderDisorderAliasPhase


    logical function OrderDisorderHelperAliasMatch(iHelperPhase, iCandidatePhase)

        implicit none

        integer, intent(in) :: iHelperPhase, iCandidatePhase
        integer             :: iHelperCSPhase, iCandidateCSPhase
        integer             :: iHelperKind, iCandidateKind
        integer             :: OrderDisorderHelperAliasKind
        character(25)       :: cHelperName, cCandidateName

        OrderDisorderHelperAliasMatch = .FALSE.
        if ((iHelperPhase <= 0).OR.(iHelperPhase > nSolnPhasesSys)) return
        if ((iCandidatePhase <= 0).OR.(iCandidatePhase > nSolnPhasesSys)) return

        iHelperCSPhase = SystemPhaseToCSPhase(iHelperPhase)
        iCandidateCSPhase = SystemPhaseToCSPhase(iCandidatePhase)
        if ((iHelperCSPhase <= 0).OR.(iCandidateCSPhase <= 0)) return

        cHelperName = UpperPhaseName(cSolnPhaseNameCS(iHelperCSPhase))
        cCandidateName = UpperPhaseName(cSolnPhaseNameCS(iCandidateCSPhase))

        if (TRIM(cHelperName) == TRIM(cCandidateName)) then
            OrderDisorderHelperAliasMatch = .TRUE.
            return
        end if

        iHelperKind = OrderDisorderHelperAliasKind(iHelperPhase)
        iCandidateKind = OrderDisorderHelperAliasKind(iCandidatePhase)

        if ((iHelperKind == 1).AND.(TRIM(cCandidateName) == 'BCC_A2')) then
            OrderDisorderHelperAliasMatch = .TRUE.
        else if ((iCandidateKind == 1).AND.(TRIM(cHelperName) == 'BCC_A2')) then
            OrderDisorderHelperAliasMatch = .TRUE.
        else if ((iHelperKind == 2).AND.(TRIM(cCandidateName) == 'FCC_A1')) then
            OrderDisorderHelperAliasMatch = .TRUE.
        else if ((iCandidateKind == 2).AND.(TRIM(cHelperName) == 'FCC_A1')) then
            OrderDisorderHelperAliasMatch = .TRUE.
        end if

        return

    end function OrderDisorderHelperAliasMatch


    integer function SystemPhaseToCSPhase(iSystemPhase)

        implicit none

        integer, intent(in) :: iSystemPhase
        integer             :: iCSPhase

        SystemPhaseToCSPhase = 0
        if (iSystemPhase <= 0) return

        do iCSPhase = 1, nSolnPhasesSysCS
            if (iSolnPhaseCS2Sys(iCSPhase) == iSystemPhase) then
                SystemPhaseToCSPhase = iCSPhase
                return
            end if
        end do

        return

    end function SystemPhaseToCSPhase


    character(25) function UpperPhaseName(cName)

        implicit none

        character(*), intent(in) :: cName
        integer                  :: iChar, iCode, nChar

        UpperPhaseName = ' '
        nChar = MIN(LEN_TRIM(cName), LEN(UpperPhaseName))
        do iChar = 1, nChar
            iCode = IACHAR(cName(iChar:iChar))
            if ((iCode >= IACHAR('a')).AND.(iCode <= IACHAR('z'))) then
                UpperPhaseName(iChar:iChar) = ACHAR(iCode - 32)
            else
                UpperPhaseName(iChar:iChar) = cName(iChar:iChar)
            end if
        end do

        return

    end function UpperPhaseName

    subroutine DetectRankDependentElements(iDependentElementLocal)

        implicit none

        integer, intent(in) :: iDependentElementLocal
        integer :: iAtomicNumber, iCandidate, nActive, nRank, nRows
        integer, dimension(:), allocatable :: iActive
        real(8), dimension(:,:), allocatable :: dRankMatrix

        do
            call BuildRankMatrix(iActive, nActive, dRankMatrix, nRows, iDependentElementLocal)
            if ((nActive <= 1).OR.(nRows <= 0)) exit

            nRank = MatrixRank(dRankMatrix, nRows, nActive)
            if (nRank >= nActive) then
                iCandidate = FindCoupledDependentElement(iActive, nActive, iDependentElementLocal)
            else
                iCandidate = FindDependentElement(iActive, nActive, dRankMatrix, nRows, nRank, &
                    iDependentElementLocal)
            end if
            if (iCandidate <= 0) exit

            iAtomicNumber = iElementSystem(iCandidate)
            iElementSystem(iCandidate) = iDependentElementLocal
            if (iAtomicNumber > 0) dElementMass(iAtomicNumber) = 0D0
        end do

        if (allocated(iActive)) deallocate(iActive)
        if (allocated(dRankMatrix)) deallocate(dRankMatrix)

        return

    end subroutine DetectRankDependentElements


    subroutine BuildRankMatrix(iActive, nActive, dRankMatrix, nRows, iDependentElementLocal)

        implicit none

        integer, intent(in) :: iDependentElementLocal
        integer, intent(out) :: nActive, nRows
        integer, dimension(:), allocatable, intent(out) :: iActive
        real(8), dimension(:,:), allocatable, intent(out) :: dRankMatrix
        integer :: iCol, iRow, iSpecies

        call CollectActiveElements(iActive, nActive, iDependentElementLocal)
        call CountRankRows(nRows, iDependentElementLocal)

        allocate(dRankMatrix(MAX(1,nRows), MAX(1,nActive)))
        dRankMatrix = 0D0
        if ((nActive <= 0).OR.(nRows <= 0)) return

        iRow = 0
        do iSpecies = 1, nSpeciesCS
            if (.NOT.SpeciesVisibleForRank(iSpecies)) cycle
            if (.NOT.SpeciesPassesRankScreen(iSpecies, iDependentElementLocal)) cycle
            iRow = iRow + 1
            do iCol = 1, nActive
                dRankMatrix(iRow,iCol) = dStoichSpeciesCS(iSpecies,iActive(iCol))
            end do
        end do

        return

    end subroutine BuildRankMatrix


    subroutine CollectActiveElements(iActive, nActive, iDependentElementLocal)

        implicit none

        integer, intent(in) :: iDependentElementLocal
        integer, intent(out) :: nActive
        integer, dimension(:), allocatable, intent(out) :: iActive
        integer :: iElement

        nActive = 0
        do iElement = 1, nElemOrComp
            if (iElementSystem(iElement) == iDependentElementLocal) cycle
            if ((iElementSystem(iElement) > 0).OR.(iElementSystem(iElement) == -1)) then
                nActive = nActive + 1
            end if
        end do

        allocate(iActive(MAX(1,nActive)))
        iActive = 0
        nActive = 0
        do iElement = 1, nElemOrComp
            if (iElementSystem(iElement) == iDependentElementLocal) cycle
            if ((iElementSystem(iElement) > 0).OR.(iElementSystem(iElement) == -1)) then
                nActive = nActive + 1
                iActive(nActive) = iElement
            end if
        end do

        return

    end subroutine CollectActiveElements


    subroutine CountRankRows(nRows, iDependentElementLocal)

        implicit none

        integer, intent(in) :: iDependentElementLocal
        integer, intent(out) :: nRows
        integer :: iSpecies

        nRows = 0
        do iSpecies = 1, nSpeciesCS
            if (.NOT.SpeciesVisibleForRank(iSpecies)) cycle
            if (.NOT.SpeciesPassesRankScreen(iSpecies, iDependentElementLocal)) cycle
            nRows = nRows + 1
        end do

        return

    end subroutine CountRankRows


    logical function SpeciesVisibleForRank(iSpecies)

        implicit none

        integer, intent(in) :: iSpecies
        integer :: iPhaseLocal

        SpeciesVisibleForRank = .FALSE.
        if (SUM(DABS(dStoichSpeciesCS(iSpecies,1:nElemOrComp))) <= 0D0) return

        if (iSpecies <= nSpeciesPhaseCS(nSolnPhasesSysCS+1)) then
            iPhaseLocal = iPhaseCS(iSpecies)
            if (iPhaseLocal <= 0) return
            if (allocated(iSolnPS)) then
                if (iPhaseLocal > SIZE(iSolnPS)) return
                if (iSolnPS(iPhaseLocal) <= 0) return
            end if
        else
            if (iPhaseCS(iSpecies) /= 0) return
        end if

        SpeciesVisibleForRank = .TRUE.

        return

    end function SpeciesVisibleForRank


    logical function SpeciesPassesRankScreen(iSpecies, iDependentElementLocal)

        implicit none

        integer, intent(in) :: iSpecies, iDependentElementLocal
        integer :: iElement
        logical :: lHasActiveStoich

        lHasActiveStoich = .FALSE.
        SpeciesPassesRankScreen = .FALSE.
        do iElement = 1, nElemOrComp
            if (((iElementSystem(iElement) > 0).OR.(iElementSystem(iElement) == -1)).AND.&
                (DABS(dStoichSpeciesCS(iSpecies,iElement)) > 0D0)) lHasActiveStoich = .TRUE.
            if ((dStoichSpeciesCS(iSpecies,iElement) > 0D0).AND.(iElementSystem(iElement) == 0)) return
        end do

        SpeciesPassesRankScreen = lHasActiveStoich

        return

    end function SpeciesPassesRankScreen


    integer function FindDependentElement(iActive, nActive, dRankMatrix, nRows, nRank, &
        iDependentElementLocal)

        implicit none

        integer, intent(in) :: nActive, nRows, nRank, iDependentElementLocal
        integer, dimension(nActive), intent(in) :: iActive
        real(8), dimension(nRows,nActive), intent(in) :: dRankMatrix
        integer :: iCandidate, iIndex, iPriority, iPriorityBest

        FindDependentElement = 0
        iPriorityBest = 1000000

        do iIndex = 1, nActive
            iCandidate = iActive(iIndex)
            if (.NOT.CanRemoveDependentColumn(iCandidate, iActive, nActive, dRankMatrix, nRows, &
                nRank, iDependentElementLocal)) cycle

            iPriority = DependentElementPriority(cElementNameCS(iCandidate))
            if (iPriority < iPriorityBest) then
                iPriorityBest = iPriority
                FindDependentElement = iCandidate
            end if
        end do

        return

    end function FindDependentElement


    integer function FindCoupledDependentElement(iActive, nActive, iDependentElementLocal)

        implicit none

        integer, intent(in) :: nActive, iDependentElementLocal
        integer, dimension(nActive), intent(in) :: iActive
        integer :: iCandidate, iIndex, iPriority, iPriorityBest

        FindCoupledDependentElement = 0
        iPriorityBest = 1000000

        do iIndex = 1, nActive
            iCandidate = iActive(iIndex)
            if (.NOT.CanUseCoupledDependentColumn(iCandidate, iActive, nActive, &
                iDependentElementLocal)) cycle

            iPriority = DependentElementPriority(cElementNameCS(iCandidate))
            if (iPriority < iPriorityBest) then
                iPriorityBest = iPriority
                FindCoupledDependentElement = iCandidate
            end if
        end do

        return

    end function FindCoupledDependentElement


    logical function CanUseCoupledDependentColumn(iCandidate, iActive, nActive, &
        iDependentElementLocal)

        implicit none

        integer, intent(in) :: iCandidate, nActive, iDependentElementLocal
        integer, dimension(nActive), intent(in) :: iActive
        integer :: iCol, iNewCol, iSpecies, nCoupledRank, nReducedRank, nRowsCoupled
        real(8), dimension(:,:), allocatable :: dCoupledMatrix, dReducedMatrix

        CanUseCoupledDependentColumn = .FALSE.
        if (iElementSystem(iCandidate) <= 0) return
        if (HasPureActiveSpecies(iCandidate, iActive, nActive, iDependentElementLocal)) return
        if (nActive <= 1) return

        nRowsCoupled = 0
        do iSpecies = 1, nSpeciesCS
            if (.NOT.SpeciesVisibleForRank(iSpecies)) cycle
            if (.NOT.SpeciesPassesRankScreen(iSpecies, iDependentElementLocal)) cycle
            if (DABS(dStoichSpeciesCS(iSpecies,iCandidate)) <= 0D0) cycle
            nRowsCoupled = nRowsCoupled + 1
        end do
        if (nRowsCoupled <= 0) return

        allocate(dCoupledMatrix(nRowsCoupled,nActive))
        allocate(dReducedMatrix(nRowsCoupled,nActive-1))
        dCoupledMatrix = 0D0
        dReducedMatrix = 0D0

        nRowsCoupled = 0
        do iSpecies = 1, nSpeciesCS
            if (.NOT.SpeciesVisibleForRank(iSpecies)) cycle
            if (.NOT.SpeciesPassesRankScreen(iSpecies, iDependentElementLocal)) cycle
            if (DABS(dStoichSpeciesCS(iSpecies,iCandidate)) <= 0D0) cycle

            nRowsCoupled = nRowsCoupled + 1
            do iCol = 1, nActive
                dCoupledMatrix(nRowsCoupled,iCol) = dStoichSpeciesCS(iSpecies,iActive(iCol))
            end do
        end do

        iNewCol = 0
        do iCol = 1, nActive
            if (iActive(iCol) == iCandidate) cycle
            iNewCol = iNewCol + 1
            dReducedMatrix(1:nRowsCoupled,iNewCol) = dCoupledMatrix(1:nRowsCoupled,iCol)
        end do

        nCoupledRank = MatrixRank(dCoupledMatrix, nRowsCoupled, nActive)
        nReducedRank = MatrixRank(dReducedMatrix, nRowsCoupled, nActive-1)
        CanUseCoupledDependentColumn = (nCoupledRank < nActive).AND.(nReducedRank == nCoupledRank)

        if (allocated(dCoupledMatrix)) deallocate(dCoupledMatrix)
        if (allocated(dReducedMatrix)) deallocate(dReducedMatrix)

        return

    end function CanUseCoupledDependentColumn


    logical function CanRemoveDependentColumn(iCandidate, iActive, nActive, dRankMatrix, nRows, &
        nRank, iDependentElementLocal)

        implicit none

        integer, intent(in) :: iCandidate, nActive, nRows, nRank, iDependentElementLocal
        integer, dimension(nActive), intent(in) :: iActive
        real(8), dimension(nRows,nActive), intent(in) :: dRankMatrix
        integer :: iCol, iNewCol, nReducedRank
        real(8), dimension(:,:), allocatable :: dReducedMatrix

        CanRemoveDependentColumn = .FALSE.
        if (iElementSystem(iCandidate) <= 0) return
        if (HasPureActiveSpecies(iCandidate, iActive, nActive, iDependentElementLocal)) return
        if (nActive <= 1) return

        allocate(dReducedMatrix(MAX(1,nRows), MAX(1,nActive-1)))
        dReducedMatrix = 0D0
        iNewCol = 0
        do iCol = 1, nActive
            if (iActive(iCol) == iCandidate) cycle
            iNewCol = iNewCol + 1
            dReducedMatrix(1:nRows,iNewCol) = dRankMatrix(1:nRows,iCol)
        end do

        nReducedRank = MatrixRank(dReducedMatrix, nRows, nActive-1)
        CanRemoveDependentColumn = (nReducedRank == nRank)

        if (allocated(dReducedMatrix)) deallocate(dReducedMatrix)

        return

    end function CanRemoveDependentColumn


    logical function HasPureActiveSpecies(iCandidate, iActive, nActive, iDependentElementLocal)

        implicit none

        integer, intent(in) :: iCandidate, nActive, iDependentElementLocal
        integer, dimension(nActive), intent(in) :: iActive
        integer :: iElement, iSpecies
        logical :: lHasOtherActive

        HasPureActiveSpecies = .FALSE.
        do iSpecies = 1, nSpeciesCS
            if (.NOT.SpeciesVisibleForRank(iSpecies)) cycle
            if (.NOT.SpeciesPassesRankScreen(iSpecies, iDependentElementLocal)) cycle
            if (DABS(dStoichSpeciesCS(iSpecies,iCandidate)) <= 0D0) cycle

            lHasOtherActive = .FALSE.
            do iElement = 1, nActive
                if (iActive(iElement) == iCandidate) cycle
                if (DABS(dStoichSpeciesCS(iSpecies,iActive(iElement))) > 0D0) then
                    lHasOtherActive = .TRUE.
                    exit
                end if
            end do
            if (.NOT.lHasOtherActive) then
                HasPureActiveSpecies = .TRUE.
                return
            end if
        end do

        return

    end function HasPureActiveSpecies


    integer function DependentElementPriority(cName)

        implicit none

        character(*), intent(in) :: cName

        select case (ADJUSTL(cName))
        case ('O')
            DependentElementPriority = 1
        case ('S')
            DependentElementPriority = 2
        case ('Se')
            DependentElementPriority = 3
        case ('Te')
            DependentElementPriority = 4
        case ('N')
            DependentElementPriority = 5
        case ('P')
            DependentElementPriority = 6
        case ('C')
            DependentElementPriority = 7
        case ('F')
            DependentElementPriority = 8
        case ('Cl')
            DependentElementPriority = 9
        case ('Br')
            DependentElementPriority = 10
        case ('I')
            DependentElementPriority = 11
        case ('H')
            DependentElementPriority = 12
        case default
            DependentElementPriority = 1000
        end select

        return

    end function DependentElementPriority


    integer function MatrixRank(dInput, nRows, nCols)

        implicit none

        integer, intent(in) :: nRows, nCols
        real(8), dimension(nRows,nCols), intent(in) :: dInput
        integer :: i, iPivotRow, iRank, j
        real(8) :: dFactor, dPivot, dScale, dTemp, dTol
        real(8), dimension(:,:), allocatable :: dWork

        MatrixRank = 0
        if ((nRows <= 0).OR.(nCols <= 0)) return

        allocate(dWork(nRows,nCols))
        dWork = dInput

        dScale = MAXVAL(DABS(dWork))
        dTol = 1D-7 * DMAX1(1D0, dScale)
        iRank = 0

        do j = 1, nCols
            iPivotRow = 0
            dPivot = dTol
            do i = iRank + 1, nRows
                if (DABS(dWork(i,j)) > dPivot) then
                    dPivot = DABS(dWork(i,j))
                    iPivotRow = i
                end if
            end do

            if (iPivotRow <= 0) cycle

            iRank = iRank + 1
            if (iPivotRow /= iRank) then
                do i = j, nCols
                    dTemp = dWork(iRank,i)
                    dWork(iRank,i) = dWork(iPivotRow,i)
                    dWork(iPivotRow,i) = dTemp
                end do
            end if

            dPivot = dWork(iRank,j)
            if (DABS(dPivot) <= dTol) cycle

            do i = j, nCols
                dWork(iRank,i) = dWork(iRank,i) / dPivot
            end do

            do i = iRank + 1, nRows
                dFactor = dWork(i,j)
                if (DABS(dFactor) <= dTol) cycle
                dWork(i,j:nCols) = dWork(i,j:nCols) - dFactor * dWork(iRank,j:nCols)
            end do
        end do

        MatrixRank = iRank

        if (allocated(dWork)) deallocate(dWork)

        return

    end function MatrixRank

end subroutine CheckSystem
!
!
