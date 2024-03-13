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
        integer                             :: i,j,k,l,m,n,o,p,q,r, nTempSolnPhases
        integer,dimension(nElements)        :: iAssemblageLast, iAssemblageNew
        real(8),dimension(nSpecies)         :: dTempPhasePotential,dMolesPerSpecies
        real(8)                             :: dTemp
        logical                             :: lStableAssemblage, lCompEveryPhases

        ! 1. Remove small fraction of negative phase amount: set to zero
        if (MINVAL(dMolesPhase)<0D0) then
            dMolesPhase(MINLOC(dMolesPhase,dim=1))=0
        end if

        ! 2. Remove the assemblage that has zero amount
        do i=1,nElements
            if (dMolesPhase(i)<1D-20) then
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
        nConPhases              = 0
        nSolnPhases             = 0
        iAssemblage             = 0
        dMolesPhase             = 0D0
        dMolesSpecies           = 0D0

        ! 3. Revise phase assemblages, phase amount, number of compounds, and solution phases to Lagrangain vairables
        ! Count the number of pure condensed phases and solution phases are assumed to be part of the phase
        ! assemblage and establish iAssemblage based on the results of Leveling:
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
                

                do j = 1,nSolnPhases
                    k = nElements - j + 1
                    
                    ! When the solution phase is already stored:
                    if (iAssemblage(k) == -l) then
                        ! o = nElements - nSolnPhases + 1
                        if (iAssemblageLast(i)>nSpecies) then
                            
                            ! Minimum points
                            dMolesPhase(k) = dMolesPhase(k)+dMolesPhaseLast(i)
                            dMolesSpecies(m:n)=dMolesSpecies(m:n)+dMolFractionGEM(i,m:n)*dMolesPhaseLast(i)
                            
                        else
                            ! Pure species in solution phase
                            dMolesPhase(k)                    = dMolesPhase(k)+dMolesPhaseLast(i)
                            dMolesSpecies(iAssemblageLast(i)) = dMolesSpecies(iAssemblageLast(i))+dMolesPhaseLast(i)
                        end if
                        cycle LOOP_AddPhase
                    end if
                end do
                ! When the solution phase is not stored before.

                nSolnPhases = nSolnPhases + 1

                j = nElements - nSolnPhases + 1
                iAssemblage(j) = -l

                if (iAssemblageLast(i)>nSpecies) then
                    ! Minimum points
                    dMolesPhase(j) = dMolesPhaseLast(i)
                    dMolesSpecies(m:n)=dMolesSpecies(m:n)+dMolFractionGEM(i,m:n)*dMolesPhaseLast(i)

                else
                    ! Pure species in solution phase

                    dMolesPhase(j)                    = dMolesPhaseLast(i)
                    dMolesSpecies(iAssemblageLast(i)) = dMolesPhaseLast(i)
                end if

                lSolnPhases(l) = .TRUE.
            end if
        end do LOOP_AddPhase

        ! 4. Re-calculate mole fraction based on dMolesSpecies, which makes it more accurate
        do i = 1, nSolnPhases
            j = -iAssemblage(nElements - i + 1)       ! Absolute solution phase index
            m = nSpeciesPhase(j-1) + 1
            n = nSpeciesPhase(j) 
            dTemp = sum(dMolesSpecies(m:n))
            if (dTemp>0D0) dMolFraction(m:n)=dMolesSpecies(m:n)/dTemp
        end do

        ! 5. Check if an immiscible phase are predicted without miscible phase (Need to revise)
        ! Check also if both phases actually have different compositions
        ! Need to debug this
        nTempSolnPhases= nSolnPhases
        iAssemblageLast = iAssemblage
        iAssemblageNew = iAssemblage

        i = 1
        do r = 1, nTempSolnPhases
            j = -iAssemblage(nElements - i + 1)       ! Absolute solution phase index
    
            if(lMiscibility(j)) then
                ! Check if miscible phase is also predicted to be stable
                lStableAssemblage=.FALSE.
                LOOP_stable:do o = 1,nTempSolnPhases
                    if ((j-1==-iAssemblage(nElements - o + 1)).AND.(cSolnPhaseName(j-1)==cSolnPhaseName(j))) then
                        !Note that the order of miscible phase is directly followed by immiscible phase
                        p=nElements - o + 1
                        q=j-1
                        lStableAssemblage=.True.
                        exit LOOP_stable
                    else if((j-2==-iAssemblage(nElements - o + 1)).AND.(cSolnPhaseName(j-2)==cSolnPhaseName(j))) then
                        p=nElements - o + 1
                        q=j-2
                        lStableAssemblage=.True.
                        exit LOOP_stable
                    else
                        p=nElements - o + 1
                        q=j-1
                        lStableAssemblage=.FALSE.
                    end if
                end do  LOOP_stable
    
                !Index of immiscible phase
                m = nSpeciesPhase(j-1) + 1
                n = nSpeciesPhase(j) 
    
                !Index of stable phase
                k = nSpeciesPhase(q-1) + 1
                l = nSpeciesPhase(q) 
    
                dTemp = MAXVAL(DABS(dMolFraction(m:n)-dMolFraction(k:l)))
    
                ! print*,'lStableAssemblage',lStableAssemblage, j,i,q
    
                if(.NOT.lStableAssemblage) then
                    iAssemblageNEW(nElements - i + 1) = -q
                    dMolFraction(k:l)=dMolFraction(m:n)
                    dMolesSpecies(k:l)=dMolesSpecies(m:n)
                else if(lStableAssemblage.AND. (dTemp < 1D-3)) then
                    
                    dMolesPhase(p)=dMolesPhase(p)+dMolesPhase(nElements - i + 1)
                    dMolesSpecies(k:l)=dMolesSpecies(m:n)+dMolesSpecies(k:l)
                    dMolFraction(k:l)=dMolesSpecies(k:l)/sum(dMolesSpecies(k:l))
    
                    dMolesSpecies(m:n)=0    
                    do p = i, nTempSolnPhases
                        iAssemblageNEW(nElements-p+1) =iAssemblage(nElements-p)
                        dMolesPhase(nElements-p+1) =dMolesPhase(nElements-p)
                    end do
                    iAssemblageNEW(nElements-nTempSolnPhases+1)=0
                    dMolesPhase(nElements-nTempSolnPhases+1)=0
                    iAssemblage = iAssemblageNEW
                    nSolnPhases = nSolnPhases-1
                    i = i-1
                end if

            end if
            i = i+1
        end do

       
    
        !Revert back modified variables
        dTempPhasePotential = dPhasePotential(:nSpecies)
    !
    !
        deallocate(dAtomFractionSpecies,dPhasePotential)
        allocate(dAtomFractionSpecies(nSpecies,nElements),dPhasePotential(nSpecies))
    !
        dAtomFractionSpecies = dAtomFractionSpeciesOld
        dPhasePotential      = dTempPhasePotential    
        iterUBC = iterGlobal

        ! ! Re-calculate chemical properties
        ! lCompEveryPhases = .FALSE.
        ! call CompChemicalPotential(lCompEveryPhases)
        
        

        ! Convert the unit from 'per mole of atom' to 'per mole of species'
        dMolesPerSpecies =sum(dStoichSpecies,dim=2)
        ! dChemicalPotential=dChemicalPotential*dMolesPerSpecies
        dChemicalPotential = MATMUL(dStoichSpecies,dElementPotential)/FLOAT(iParticlesPerMole)
        dPartialEnthalpy=dPartialEnthalpy*dMolesPerSpecies
        dPartialEntropy=dPartialEntropy*dMolesPerSpecies
        dPartialHeatCapacity=dPartialHeatCapacity*dMolesPerSpecies


    end subroutine PostProcessPEA