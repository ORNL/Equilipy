!
!
subroutine CompDrivingForceAll
!
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CompDrivingForceAll.f90
    !> \brief   Compute the driving force of all phases.
    !> \author  S.Y. Kwon
    !> \date    Dec. 20, 2021
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   12/20/2021      S.Y. Kwon         Original code
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to compute the driving force of all pure condensed phases in
    !! the database.  The driving force is defined as the difference between the standard molar Gibbs energy of
    !! a pure condensed phase and the corresponding value computed from the element potentials.  This value is
    !! used to determine whether a pure condensed phase should be added to the system.  For a more thorough
    !! explanation of the chemical significance of the driving force, refer to the following literature:
    !!
    !! H.L. Lukas, S.G. Fries, B. Sundman, Computational Thermodynamics - The Calphad Method, Cambridge
    !! University Press, New York, 2007.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[out]  iMaxDrivingForce    An integer scalar representing the index of the pure condensed phase
    !!                                   with the maximum driving force.  A value of zero is returned by default.
    !> \param[out]  dMaxDrivingForce    A double real scalar representing the maximum driving force of all pure
    !!                                   condensed phases.  A value of zero is returned by default.
    !
    !-------------------------------------------------------------------------------------------------------------
    USE ModuleThermo
    USE ModuleGEMSolver
!
    implicit none
!
    integer :: i, j,k,  m, n, nConstituents, nElementOrConstituent
    real(8) :: dTemp, dMinMoleFraction,dMaxMoleFraction
    integer, dimension(:), allocatable:: iIndex, iEndmemberPotential
    real(8), dimension(:), allocatable :: dDrivingForceTemp, dEndmemberPotential
    real(8), dimension(:,:), allocatable :: dChemicalPotentialTemp, dMolFractionTemp, dPhasePotentialTemp
    logical :: ldummy, lAddPhase

!
    dPhasePotential = 0d0
    
    

!
    m = nSpeciesPhase(nSolnPhasesSys) + 1
    n = nSpecies - nDummySpecies

!
    ! Calculate driving force of all stoichiometric phases:
    do i = m, n
!
        ! Compute the chemical potential of this phase as defined by the element potentials:
        dTemp = 0D0
        do j = 1, nElements
            dTemp = dTemp + dElementPotential(j) * dStoichSpecies(i,j)
        end do
!
!
        ! Compute the driving force:
        dTemp = dStdGibbsEnergy(i) - dTemp
!
        ! Normalize per gram-atom:
        dTemp = dTemp / dSpeciesTotalAtoms(i)
!
        ! !SY: Update phase potential of condensted phase
        dPhasePotential(i) = dTemp
!
!
    end do
!
    !Store the most stable stoichiometric phase to be added
    dMinMoleFraction = 1D-5
!
    ! Calculate driving force of all solution phases:
    LOOP_SolnPhaseSys: do i = 1, nSolnPhasesSys
        !initialize variables
        m                = nSpeciesPhase(i-1) + 1      ! First constituent in phase.
        n                = nSpeciesPhase(i)            ! Last  constituent in phase.
        nConstituents    = n - m + 1
        dMaxMoleFraction = 1D0 - dMinMoleFraction * DFLOAT(nConstituents-1)
        dMaxMoleFraction = DMAX1(dMaxMoleFraction, 0.9D0)
        nElementOrConstituent = MIN(nElements,nConstituents)
        ! if (nElementOrConstituent<=1) cycle LOOP_SolnPhaseSys
        lAddPhase        = .FALSE.     
        if (.NOT.lMiscibility(i)) then
            if(allocated(iIndex)) deallocate(iIndex)
            if(allocated(dDrivingForceTemp)) deallocate(dDrivingForceTemp)
            if(allocated(dChemicalPotentialTemp)) deallocate(dChemicalPotentialTemp)
            if(allocated(dPhasePotentialTemp)) deallocate(dPhasePotentialTemp)
            if(allocated(dMolFractionTemp)) deallocate(dMolFractionTemp)
            if(allocated(dEndmemberPotential)) deallocate(dEndmemberPotential)
            if(allocated(iEndmemberPotential)) deallocate(iEndmemberPotential)
            
            allocate(dEndmemberPotential(nConstituents),iEndmemberPotential(nConstituents))
            allocate(&
            iIndex(nElementOrConstituent+1),&
            dDrivingForceTemp(nElementOrConstituent+1),&
            dChemicalPotentialTemp(nElementOrConstituent+1,nConstituents),&
            dMolFractionTemp(nElementOrConstituent+1, nConstituents),&
            dPhasePotentialTemp(nElementOrConstituent+1, nConstituents))
            dEndmemberPotential= dStdGibbsEnergy(m:n)- MATMUL(dStoichSpecies(m:n,:),dElementPotential)/iParticlesPerMole(m:n)
            call Qsort(dEndmemberPotential,iEndmemberPotential,nConstituents)

            iIndex             = 1
            dDrivingForceTemp  = 9E5
            dMolFractionTemp   = 0D0
            do j = 1, nElementOrConstituent+1
            
                if(j==1) then
                    dMolFraction(m:n)   = dMolFraction(m:n)/sum(dMolFraction(m:n))
                else
                    dMolFraction(m:n)   = dMinMoleFraction
                    dMolFraction(m+iEndmemberPotential(j-1)-1) = dMaxMoleFraction
                end if
                ! Perform subminimization:
                call Subminimization(i, lAddPhase)
    !
                ! Store info of all local minima
                call CompStoichSolnPhase(i)
                
                dDrivingForceTemp(j)   = dDrivingForceSoln(i)
                dMolFractionTemp(j,:)  = dMolFraction(m:n)
                dChemicalPotentialTemp(j,:) = dChemicalPotential(m:n)
                dPhasePotentialTemp(j,:)    = dPhasePotential(m:n)
            end do

            !Sort the results according to the driving force (ascending order)
            ! note that output dDrivingForceTemp is the sorted version
            call Qsort(dDrivingForceTemp,iIndex,nElementOrConstituent+1)

            dMolFraction(m:n)       = dMolFractionTemp(iIndex(1),:)
            dDrivingForceSoln(i)    = dDrivingForceTemp(1)
            dChemicalPotential(m:n) = dChemicalPotentialTemp(iIndex(1),:)
            dPhasePotential(m:n)    = dPhasePotentialTemp(iIndex(1),:)

            LOOP_immiscibleSoln: do j = 2, nElementOrConstituent
                k = i + j - 1
                m  = nSpeciesPhase(k-1) + 1      ! First constituent in phase.
                n  = nSpeciesPhase(k)            ! Last  constituent in phase.
!
                !The following solution phases must have True value for lMiscibility
                if (lMiscibility(k).AND.(k<=nSolnPhasesSys).AND.(cSolnPhaseName(i)==cSolnPhaseName(k))) then
!
                    ! Store the local minima if it is not equal to global minimum
                    dMolFraction(m:n)       = dMolFractionTemp(iIndex(j),:)
                    dDrivingForceSoln(k)    = dDrivingForceTemp(j)
                    dChemicalPotential(m:n) = dChemicalPotentialTemp(iIndex(j),:)
                    dPhasePotential(m:n)    = dPhasePotentialTemp(iIndex(j),:)
                end if
!
            end do LOOP_immiscibleSoln

        end if
!
    end do LOOP_SolnPhaseSys
    
!
    return
!
end subroutine CompDrivingForceAll
!
!
!
