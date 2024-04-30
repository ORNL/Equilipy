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
        integer                             :: i,j,k,l,m,n,o,p,q,r,s, nTempSolnPhases,iMiscIndex,iImmiscIndex
        integer,dimension(nElements)        :: iAssemblageLast, iAssemblageNew
        real(8),dimension(nSpecies)         :: dTempPhasePotential,dMolesPerSpecies
        real(8)                             :: dTemp
        logical                             :: lStableAssemblage, lCompEveryPhases, lDuplicate

        lCompEveryPhases =.FALSE.
        lSolnPhases = .FALSE.

        ! 1. Remove small fraction of negative phase amount: set to zero
        if (MINVAL(dMolesPhase)<0D0) then
            dMolesPhase(MINLOC(dMolesPhase,dim=1))=0
        end if

        ! 2. Remove the assemblage that has zero amount
        do i=1,nElements
            if (dMolesPhase(i)<1D-50) then
                iAssemblage(i)=0
            end if
        end do

        ! Memory allocation for the post processing
        if (allocated(dMolesPhaseLast)) deallocate(dMolesPhaseLast)
        allocate(dMolesPhaseLast(nElements))
        iAssemblageLast         = iAssemblage
        dMolesPhaseLast         = dMolesPhase
        dTempPhasePotential     = 0D0
        nConPhases              = 0
        nSolnPhases             = 0
        iAssemblage             = 0
        dMolesPhase             = 0D0
        dMolesSpecies           = 0D0
        lSkipLagrange           = .False.

        ! 3. Revise phase assemblages, phase amount, number of compounds, and solution phases to Lagrangain vairables
        ! Count the number of pure condensed phases and solution phases are assumed to be part of the phase
        ! assemblage and establish iAssemblage based on the results of Leveling:
        
        LOOP_AddPhase: do i = 1, nElements
            if (iAssemblageLast(i)==0) then
                l = -1
            else
                l = iPhaseGEM(i)
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
                
                ! Get Soln index of miscible phase
                loop_misc1:do p = l, 1, -1
                    if(.Not.lMiscibility(p)) then
                        iMiscIndex =p
                        exit loop_misc1
                    end if
                end do loop_misc1

                ! Loop through solution phases and check the sln phase is already stored
                do j = 1,nSolnPhases
                    k = nElements - j + 1
                    r= -iAssemblage(k)

                    !Check next solution phase when the solution phase is not stored
                    if (cSolnPhaseName(r) /= cSolnPhaseName(l)) cycle
                    
                    !Check if the immiscibility is considered in the system
                    ! When the solution phase is already stored, check duplicate:
                    lDuplicate = .TRUE.
                        
                    p = nSpeciesPhase(r-1) + 1      ! First constituent in phase.
                    q = nSpeciesPhase(r)            ! Last  constituent in phase.     

                    dTemp = MAXVAL(DABS(MATMUL(dMolFractionGEM(i,m:n),dStoichSpecies(m:n,:))-&
                    MATMUL(dMolFraction(p:q),dStoichSpecies(p:q,:))))
                    ! print*, 'Processing ', l, cSolnPhaseName(l), MATMUL(dMolFractionGEM(i,m:n),dStoichSpecies(m:n,:))
                    ! print*, 'Comparing with', r, cSolnPhaseName(r), MATMUL(dMolFraction(p:q),dStoichSpecies(p:q,:))

                    Immiscibility: do o = iMiscIndex, nSolnPhasesSys
                        if(lMiscibility(o).and.dTemp > 5D-2) then
                            lDuplicate = .FALSE.
                            exit Immiscibility
                        end if
                    end do Immiscibility

                    if (lDuplicate) then
                        ! Phase are duplicated. Just add them to previous one
                        if (iAssemblageLast(i)>nSpecies) then
                            ! Minimum points
                            dMolesPhase(k) = dMolesPhase(k)+dMolesPhaseLast(i)
                            dMolesSpecies(p:q)=dMolesSpecies(p:q)+dMolFractionGEM(i,m:n)*dMolesPhaseLast(i)
                        else
                            ! Pure species in solution phase
                            dMolesPhase(k)                    = dMolesPhase(k)+dMolesPhaseLast(i)
                            dMolesSpecies(p+iAssemblageLast(i)-m) = dMolesSpecies(iAssemblageLast(i))+dMolesPhaseLast(i)
                        end if

                        cycle LOOP_AddPhase
                    else
                        ! when the phase is not duplicated through out all stored sln phases, create immiscible phase
                        ! When the immiscible phase actually has different compositions, treat the phase as if it has not been stored.
                        
                        iImmiscIndex = iMiscIndex+1

                        do s= 1, nSolnPhases
                            o = nElements - nSolnPhases + 1
                            if(iAssemblage(o)==-iImmiscIndex) iImmiscIndex = iImmiscIndex+1
                        end do
                        nSolnPhases = nSolnPhases + 1
                        o = nElements - nSolnPhases + 1

                        iAssemblage(o) = -(iImmiscIndex)
                        p = nSpeciesPhase(iImmiscIndex-1) + 1      ! First constituent in phase.
                        q = nSpeciesPhase(iImmiscIndex)            ! Last  constituent in phase.     
                        
                        if (iAssemblageLast(i)>nSpecies) then
                            ! Minimum points
                            dMolesPhase(o) = dMolesPhaseLast(i)
                            dMolesSpecies(p:q)=dMolFractionGEM(i,m:n)*dMolesPhaseLast(i)
        
                        else
                            ! Pure species in solution phase
                            dMolesPhase(o)                    = dMolesPhaseLast(i)
                            dMolesSpecies(p) = dMolesPhaseLast(i)
                        end if
        
                        lSolnPhases(o) = .TRUE.
                        cycle LOOP_AddPhase
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
                        dMolFraction(m:n)=dMolFractionGEM(i,m:n)

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
                    
                    ! Any of relavent phase has not been stored. Store the miscible version of phase
                    !Check if immiscible phase has the same composition as previous phase
                    iAssemblage(j) = -iMiscIndex
                    if (iAssemblageLast(i)>nSpecies) then
                        ! Minimum points
                        dMolesPhase(j) = dMolesPhaseLast(i)
                        dMolesSpecies(p:q)=dMolFractionGEM(i,m:n)*dMolesPhaseLast(i)
                        dMolFraction(p:q)=dMolFractionGEM(i,m:n)
                    else
                        ! Pure species in solution phase
                        dMolesPhase(j)                    = dMolesPhaseLast(i)
                        dMolesSpecies(p+iAssemblageLast(i)-m) = dMolesPhaseLast(i)
                        dMolesSpecies(iAssemblageLast(i)) = 0D0
                    end if
                    lSolnPhases(iMiscIndex) = .TRUE.
                end if
            end if
        end do LOOP_AddPhase

        ! 4. Re-calculate mole fraction based on dMolesSpecies, which makes it more accurate
        do i = 1, nSolnPhases
            j = -iAssemblage(nElements - i + 1)       ! Absolute solution phase index
            m = nSpeciesPhase(j-1) + 1
            n = nSpeciesPhase(j) 
            dTemp = sum(dMolesSpecies(m:n))
            if (dTemp>0D0) dMolFraction(m:n)=dMolesSpecies(m:n)/dTemp
            if(MINVAL(dMolesSpecies(m:n))<1D-100) then
                ! print*, 'Skip LG:', cSolnPhaseName(j)
                lSkipLagrange = .True.
            end if
        end do
        
    
        !Revert back modified variables
        dTempPhasePotential = dPhasePotential(:nSpecies)
        deallocate(dAtomFractionSpecies,dPhasePotential,dChemicalPotential)
        allocate(dAtomFractionSpecies(nSpecies,nElements),dPhasePotential(nSpecies),dChemicalPotential(nSpecies))

        dAtomFractionSpecies = dAtomFractionSpeciesOld
        dPhasePotential      = dTempPhasePotential  
        iterUBC = iterGlobal

        !Calculate chemical potentials based on mole fraction
        call CompChemicalPotential(lCompEveryPhases)
        !Double check duplicate and combine them again

        call CompFunctionNorm
        dGEMFunctionNormLast=dGEMFunctionNorm

        ! dChemicalPotential = MATMUL(dStoichSpecies,dElementPotential)/FLOAT(iParticlesPerMole)
        ! dPartialEnthalpy=dPartialEnthalpy*dSpeciesTotalAtoms
        ! dPartialEntropy=dPartialEntropy*dSpeciesTotalAtoms
        ! dPartialHeatCapacity=dPartialHeatCapacity*dSpeciesTotalAtoms

        deallocate(dAtomFractionSpeciesGEM,dChemicalPotentialGEM,dStoichSpeciesGEM,iPhaseGEM,&
        dMolesPhaseHistory,dMolFractionGEM,dMolesPhaseLast,dMolFractionOld,iPhaseLevel,dStoichSpeciesLevel)
        ! if(dGEMFunctionNormLast<1D-5) lSkipLagrange = .True.

    end subroutine PostProcessPEA