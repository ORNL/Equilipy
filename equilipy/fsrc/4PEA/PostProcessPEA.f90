subroutine PostProcessPEA
    !
        !-------------------------------------------------------------------------------------------------------------
        !
        !> \file    GEMSolver.f90
        !> \brief   Gibbs Energy Minimization solver.
        !> \author  S.Y. Kwon
        !> \date    May 2, 2022
        !
        !
        ! Revisions:
        ! ==========
        !
        !   Date            Programmer          Description of change
        !   ----            ----------          ---------------------
        !   11/02/2023      S.Y. Kwon            Original code
        !
        ! Purpose:
        ! ========
        !
        !> \details TBA

        ! Pertinent variables:
        ! ====================
        !
        ! nConPhases            The number of pure condensed phases in the assemblage
        ! nSolnPhases           The number of solution phases in the assemblage
        ! nSolnPhasesSys        The number of solution phases in the system
        !-------------------------------------------------------------------------------------------------------------
        USE ModuleThermo
        USE ModuleThermoIO
        USE ModuleGEMSolver
    !
        implicit none
    !
        integer                             :: i,j,k,l,m,n,o,p,q,r, nTempSolnPhases,iMiscIndex
        integer,dimension(nElements)        :: iAssemblageLast, iAssemblageNew
        real(8),dimension(nSpecies)         :: dTempPhasePotential,dMolesPerSpecies,dTempMu
        real(8)                             :: dTemp
        logical                             :: lStableAssemblage, lCompEveryPhases, lDuplicate

        ! 1. Remove small fraction of negative phase amount: set to zero
        if (MINVAL(dMolesPhase)<0D0) then
            dMolesPhase(MINLOC(dMolesPhase,dim=1))=0
        end if

        ! 2. Remove the assemblage that has zero amount
        do i=1,nElements
            if (dMolesPhase(i)<1D-50) then
                k = iAssemblage(i) 
                m=nSpeciesPhase(k-1) + 1
                n= nSpeciesPhase(k) 
                iAssemblage(i)=0
                nSolnPhases=nSolnPhases-1
            end if
        end do

        ! Memory allocation for the post processing
        if (allocated(dMolesPhaseLast)) deallocate(dMolesPhaseLast)
        allocate(dMolesPhaseLast(nElements))
        iAssemblageLast         = iAssemblage
        dMolesPhaseLast         = dMolesPhase
        dTempPhasePotential     = 0D0
        dTempMu                 = 0D0
        nConPhases              = 0
        nSolnPhases             = 0
        iAssemblage             = 0
        dMolesPhase             = 0D0
        dMolesSpecies           = 0D0

        ! 3. Revise phase assemblages, phase amount, number of compounds, and solution phases to Lagrangain vairables
        ! Count the number of pure condensed phases and solution phases are assumed to be part of the phase
        ! assemblage and establish iAssemblage based on the results of Leveling:
        ! print*, 'iAssemblage Last:',iAssemblageLast
        
        LOOP_AddPhase: do i = 1, nElements
            if (iAssemblageLast(i)==0) then
                l = -1
            else
                l = iPhaseLevel(iAssemblageLast(i))
            end if
            

            if ((l == 0).AND.(nConPhases + nSolnPhases <= nElements)) then
                ! Stoichiometric phase
                nConPhases                        = nConPhases + 1
                iAssemblage(nConPhases)           = iAssemblageLast(i)
                dMolesPhase(nConPhases)           = dMolesPhaseLast(i)
                dMolesSpecies(iAssemblageLast(i)) = dMolesPhaseLast(i)

            elseif ((l > 0).AND.(nConPhases + nSolnPhases <= nElements)) then
                ! Solution phases
                m = nSpeciesPhase(l-1) + 1      ! First constituent in phase.
                n = nSpeciesPhase(l)            ! Last  constituent in phase.     
                

                ! Loop through solution phases and check the sln phase is already stored
                do j = 1,nSolnPhases
                    k = nElements - j + 1

                    ! Just add the phase if the current sln phase is micisble phase
                    if(.not.lMiscibility(l)) then
                        ! When the solution phase is already stored:
                        if (iAssemblage(k) == -l) then
                            if (iAssemblageLast(i)>nSpecies) then
                                ! Minimum points
                                dMolesPhase(k) = dMolesPhase(k)+dMolesPhaseLast(i)
                                dMolesSpecies(m:n)=dMolesSpecies(m:n)+dMolFractionGEM(i,m:n)*dMolesPhaseLast(i)
                                
                            else
                                ! Pure species in solution phase
                                dMolesPhase(k)                    = dMolesPhase(k)+dMolesPhaseLast(i)
                                dMolesSpecies(iAssemblageLast(i)) = dMolesSpecies(iAssemblageLast(i))+dMolesPhaseLast(i)
                            end if
                            ! print*,'Stored before and Miscible'
                            ! print*,'processing: ',iAssemblageLast(i),cSolnPhaseName(l),iAssemblage(i)
                            cycle LOOP_AddPhase
                        end if
                    ! If the current phase is immisicble phase, check if corresponding miscible phase is stored before
                    else
                        ! Get Soln index of miscible phase
                        loop_misc1:do p = l, 1, -1
                            if(.Not.lMiscibility(p)) then
                                iMiscIndex =p
                                exit loop_misc1
                            end if
                        end do loop_misc1

                        ! Check if the phase is stored through all corresponding miscible immisible phases
                        do r = iMiscIndex,l

                            !Check if immiscible phase has the same composition as previous phase
                            lDuplicate = .False.
                            ! Solution phases
                            p = nSpeciesPhase(r-1) + 1      ! First constituent in phase.
                            q = nSpeciesPhase(r)            ! Last  constituent in phase.     
                            dTemp = SUM(DABS(dMolFractionGEM(i,m:n)-dMolFraction(p:q)))/float(q-p)
                            if (dTemp < 1D-2) lDuplicate = .TRUE.

                            ! When the relevant phase is already stored:
                            if ((iAssemblage(k) == -r).AND.lDuplicate) then
                                if (iAssemblageLast(i)>nSpecies) then
                                    ! Minimum points
                                    dMolesPhase(k) = dMolesPhase(k)+dMolesPhaseLast(i)
                                    dMolesSpecies(p:q)=dMolesSpecies(p:q)+dMolFractionGEM(i,m:n)*dMolesPhaseLast(i)
                                    dMolesSpecies(m:n) = 0D0
                                else
                                    ! Pure species in solution phase
                                    dMolesPhase(k)                    = dMolesPhase(k)+dMolesPhaseLast(i)
                                    dMolesSpecies(p+iAssemblageLast(i)-m) = dMolesSpecies(iAssemblageLast(i))+dMolesPhaseLast(i)
                                    dMolesSpecies(iAssemblageLast(i)) = 0D0
                                end if
                                ! print*,'Stored before, immiscible, duplicated'
                                ! print*,'processing: ',iAssemblageLast(i),cSolnPhaseName(l),iAssemblage(i)

                                cycle LOOP_AddPhase
                            ! When the immiscible phase actually has different compositions, treat the phase as if it has not been stored.
                            else if ((iAssemblage(k) == -r).AND.(.not.lDuplicate)) then
                                
                                nSolnPhases = nSolnPhases + 1
                                o = nElements - nSolnPhases + 1
                                iAssemblage(o) = -l

                                if (iAssemblageLast(i)>nSpecies) then
                                    ! Minimum points
                                    dMolesPhase(o) = dMolesPhaseLast(i)
                                    dMolesSpecies(m:n)=dMolFractionGEM(i,m:n)*dMolesPhaseLast(i)
                
                                else
                                    ! Pure species in solution phase
                                    dMolesPhase(o)                    = dMolesPhaseLast(i)
                                    dMolesSpecies(iAssemblageLast(i)) = dMolesPhaseLast(i)
                                end if
                
                                lSolnPhases(o) = .TRUE.
                                ! print*,'Stored before, immiscible, not duplicated'
                                ! print*,'processing: ',iAssemblageLast(i),cSolnPhaseName(l),iAssemblage(o),i,o
                                cycle LOOP_AddPhase
                            end if
                        end do

                    end if
                end do

                ! When the solution phase is not stored before.
                ! Just add the phase if it is not a immiscible phase
                nSolnPhases = nSolnPhases + 1
                j = nElements - nSolnPhases + 1
                

                if(.not.lMiscibility(l)) then
                    ! print*,'Not stored before and Miscible'
                    iAssemblage(j) = -l
                    if (iAssemblageLast(i)>nSpecies) then
                        ! Minimum points
                        dMolesPhase(j) = dMolesPhaseLast(i)
                        dMolesSpecies(m:n)=dMolFractionGEM(i,m:n)*dMolesPhaseLast(i)

                    else
                        ! Pure species in solution phase
                        dMolesPhase(j)                    = dMolesPhaseLast(i)
                        dMolesSpecies(iAssemblageLast(i)) = dMolesPhaseLast(i)
                    end if

                    lSolnPhases(l) = .TRUE.
                ! If the current phase is immisicble phase, check if corresponding miscible phase is stored before 
                else
                    ! Get Soln index of miscible phase
                    loop_misc2:do p = l, 1, -1
                        if(.Not.lMiscibility(p)) then
                            iMiscIndex =p
                            exit loop_misc2
                        end if
                    end do loop_misc2
                    p = nSpeciesPhase(iMiscIndex-1) + 1      ! First constituent in phase.
                    q = nSpeciesPhase(iMiscIndex)            ! Last  constituent in phase.     
                    ! print*,'Not stored before and immisicible: Changing from ',l, 'to',iMiscIndex
                    ! Any of relavent phase has not been stored. Store the miscible version of phase
                    !Check if immiscible phase has the same composition as previous phase
                    iAssemblage(j) = -iMiscIndex
                    if (iAssemblageLast(i)>nSpecies) then
                        ! Minimum points
                        dMolesPhase(j) = dMolesPhaseLast(i)
                        dMolesSpecies(p:q)=dMolFractionGEM(i,m:n)*dMolesPhaseLast(i)
                        dMolesSpecies(m:n) = 0D0
                    else
                        ! Pure species in solution phase
                        dMolesPhase(j)                    = dMolesPhaseLast(i)
                        dMolesSpecies(p+iAssemblageLast(i)-m) = dMolesPhaseLast(i)
                        dMolesSpecies(iAssemblageLast(i)) = 0D0
                    end if
                    lSolnPhases(iMiscIndex) = .TRUE.
                    lSolnPhases(iMiscIndex) = .FALSE.
                end if
            end if
            ! print*,'processing: ',iAssemblageLast(i),cSolnPhaseName(l),iAssemblage(j),i,j
        end do LOOP_AddPhase
        ! print*,'iAssemblage new',iAssemblage

        ! 4. Re-calculate mole fraction based on dMolesSpecies, which makes it more accurate
        do i = 1, nSolnPhases
            j = -iAssemblage(nElements - i + 1)       ! Absolute solution phase index
            m = nSpeciesPhase(j-1) + 1
            n = nSpeciesPhase(j) 
            dTemp = sum(dMolesSpecies(m:n))
            if (dTemp>0D0) dMolFraction(m:n)=dMolesSpecies(m:n)/dTemp
        end do
       
    
        !Revert back modified variables
        dTempPhasePotential = dPhasePotential(:nSpecies)
        dTempMu             = dChemicalPotential(:nSpecies)
    !
    !
        deallocate(dAtomFractionSpecies,dPhasePotential,dChemicalPotential)
        allocate(dAtomFractionSpecies(nSpecies,nElements),dPhasePotential(nSpecies),dChemicalPotential(nSpecies))
    !
        dAtomFractionSpecies = dAtomFractionSpeciesOld
        dPhasePotential      = dTempPhasePotential    
        dChemicalPotential   = dTempMu
        iterUBC = iterGlobal

        ! Re-calculate chemical properties
        ! lCompEveryPhases = .True.
        ! call CompDrivingForceAll
        ! print*, 'element potential',dElementPotential
        ! print*,'dMolesElement',dMolesElement
        
        
        ! print*, sum(dMolesSpecies)
        ! Convert the unit from 'per mole of atom' to 'per mole of species'
        ! dMolesPerSpecies =sum(dStoichSpecies,dim=2)
        ! print*,'dChemicalPotential',dChemicalPotential
        ! dChemicalPotential=dChemicalPotential*dMolesPerSpecies
        dChemicalPotential = MATMUL(dStoichSpecies,dElementPotential)/FLOAT(iParticlesPerMole)
        dPartialEnthalpy=dPartialEnthalpy*dSpeciesTotalAtoms
        dPartialEntropy=dPartialEntropy*dSpeciesTotalAtoms
        dPartialHeatCapacity=dPartialHeatCapacity*dSpeciesTotalAtoms

        deallocate(dAtomFractionSpeciesGEM,dChemicalPotentialGEM,dStoichSpeciesGEM,&
        dMolesPhaseHistory,dMolFractionGEM,dMolesPhaseLast,dMolFractionOld,iPhaseLevel,dStoichSpeciesLevel)


    end subroutine PostProcessPEA