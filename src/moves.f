        MODULE moves
        USE basin_hopping       ,ONLY: bhop_move
        USE cut_splice          ,ONLY: cut_splice_move
        USE monte               ,ONLY: is_cut_splice
     &    , pre_cut_splice_period
     &    , cut_splice_freq
!DEBUG BEGINS==============================================
!       USE potential           ,ONLY: potential_type
!DEBUG ENDS==============================================
        PUBLIC          :: generate_config

        CONTAINS

        SUBROUTINE generate_config(mc_step)
        INTEGER,INTENT(IN)      :: mc_step
        INTEGER                 :: remainder

        remainder = MOD(mc_step, cut_splice_freq)

        IF ( is_cut_splice .AND.
     &    mc_step > pre_cut_splice_period .AND.
     &    remainder == 0) THEN

          CALL cut_splice_move
!DEBUG BEGINS==============================================
!         PRINT *, 'after_cut_splice_move'
!DEBUG ENDS==============================================

        ELSE

!DEBUG BEGINS==============================================
!       IF(potential_type==2) PRINT *, 'before bhop_move '
!DEBUG ENDS==============================================
          CALL bhop_move
!DEBUG BEGINS==============================================
!       IF(potential_type==2) PRINT *, 'after bhop_move '
!DEBUG ENDS==============================================

!DEBUG BEGINS==============================================
!         PRINT *, 'after bhop move'
!DEBUG ENDS==============================================

        END IF

        END SUBROUTINE generate_config

        END MODULE moves
