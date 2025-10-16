!
subroutine CheckLevelingIterHistory(iter,lPhasePass)
    !---------------------------------------------------------------------------
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this subroutine is to check the iteration history for a
    ! leveling iteration.  The iteration history is scanned backwards and then
    ! every phase in the assemblage in question is compared to the full
    ! assemblage at a particular iteration.  Thus, the order of phases in the
    ! integer array does not affect the calculation.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in]   iter        The iteration index in Leveling.
    !> \param[out]  lPhasePass  A logical scalar indicating whether the phase
    !                            assemblaged has passed (TRUE) or not (FALSE).
    !
    !---------------------------------------------------------------------------
    USE ModuleThermo
!
    implicit none
!
    integer:: i, j,k, iter
    logical:: lPhasePass, lSameAssemblage
!
!
    ! Initialize variables:
    lPhasePass = .TRUE.

    ! Loop through iteration history:
    LOOP_Iter: do i = iter, 1, -1 
        lSameAssemblage = .TRUE.       

        LOOP_NEW:do j = 1, nElements
            LOOP_OLD: do k = 1, nElements
                if (iAssemblage(j)==iterHistoryLevel(k,i)) cycle LOOP_NEW
            end do LOOP_OLD
            lSameAssemblage = .FALSE.
        end do LOOP_NEW

        if(lSameAssemblage) then
            lPhasePass = .FALSE.
            exit LOOP_Iter
        end if
        
    end do LOOP_Iter
!

    return
!
end subroutine CheckLevelingIterHistory
!
!
!
