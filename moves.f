        MODULE moves
        USE basin_hopping       ,ONLY: bhop_move
        USE cut_splice          ,ONLY: cut_splice_move
        USE monte               ,ONLY: pre_cut_splice_period
     &    , cut_splice_freq
!DEBUG BEGINS==============================================
!       USE potential           ,ONLY: potential_type
!DEBUG ENDS==============================================
        CONTAINS

        SUBROUTINE generate_config(mc_step)
        INTEGER,INTENT(IN)      :: mc_step

        IF (mc_step > pre_cut_splice_period .AND.
     &    mc_step == cut_splice_freq) THEN

          CALL cut_splice_move

        ELSE

!DEBUG BEGINS==============================================
!       IF(potential_type==2) PRINT *, 'before bhop_move '
!DEBUG ENDS==============================================
          CALL bhop_move
!DEBUG BEGINS==============================================
!       IF(potential_type==2) PRINT *, 'after bhop_move '
!DEBUG ENDS==============================================

        END IF

        END SUBROUTINE generate_config

        END MODULE moves
