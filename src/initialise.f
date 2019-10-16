        MODULE initialise

        LOGICAL             :: is_resume_calc
        INTEGER             :: total_mc_step
        LOGICAL             :: is_multi_stage_potent
        INTEGER             :: mc_step_1, mc_step_2
        INTEGER             :: mc_save_step 
        INTEGER             :: resume_mc_step
        CHARACTER(LEN=21) :: 
     &    resume_mc_step_dat='01_resume_MC_step.dat'
        CHARACTER(LEN=19) :: 
     &    all_mc_steps='02_all_MC_steps.dat'



        PRIVATE :: printout_session_in


        PUBLIC  :: check_mc_step
        PUBLIC  :: read_session
        PUBLIC  :: read_params
        PUBLIC  :: printout_session_feed
        PUBLIC  :: periodic_save_mc_step

        CONTAINS

        SUBROUTINE read_resume_mc_step
        IMPLICIT NONE
        INTEGER                 :: f_mc
        OPEN(NEWUNIT=f_mc,FILE=resume_mc_step_dat,STATUS='old')
        READ(f_mc,*) resume_mc_step
        CLOSE(f_mc)
!DEBUG STARTS=================================================
!         PRINT *, 'resume mc step is', resume_mc_step
!DEBUG ENDS===================================================
        END SUBROUTINE read_resume_mc_step

        SUBROUTINE periodic_save_mc_step(iteration)
        IMPLICIT NONE
        INTEGER,INTENT(IN)      :: iteration
        INTEGER                 :: f_mc
        INTEGER                 :: remainder

        remainder = MOD(iteration,mc_save_step)
        IF (remainder == 0) THEN

        OPEN(NEWUNIT=f_mc,FILE=resume_mc_step_dat,STATUS='replace')
        WRITE(f_mc,*) iteration 
        CLOSE(f_mc)

        OPEN(NEWUNIT=f_mc,FILE=all_mc_steps,ACCESS='append')
        WRITE(f_mc,*) iteration 
        CLOSE(f_mc)

        END IF

        END SUBROUTINE periodic_save_mc_step

        SUBROUTINE check_mc_step(iteration)
! print MC steps periodically for testing and debugging
        IMPLICIT NONE
        INTEGER,INTENT(IN)      :: iteration
        INTEGER                 :: remainder
        INTEGER                 :: f_in
        REAL                    :: rec_time
        
        CALL CPU_TIME(rec_time)
        remainder = MOD(iteration,mc_save_step)

        IF (remainder == 0) THEN
        OPEN(NEWUNIT=f_in,FILE='mc_steps_taken.dat',ACCESS='append')
        WRITE(f_in,*) 'MC step:', iteration,'time:', rec_time
        CLOSE(f_in)
        END IF

        END SUBROUTINE check_mc_step

        SUBROUTINE read_session
! read session simulation parameters
        USE directory           ,ONLY: session_dir
        USE coord_grad_ene,     ONLY: atoms,material,in_file
        USE potential,          ONLY: potential_type
     &    ,potent_1, potent_2
        IMPLICIT NONE
        CHARACTER(14)   :: session_file='session_in.dat'
        INTEGER         :: f_session
        CHARACTER       :: dummy
        INTEGER         :: iter

        OPEN(NEWUNIT=f_session
     &    ,FILE='./'//session_dir//'/'//session_file,STATUS='old')
        READ(f_session,*) dummy, atoms       
        READ(f_session,*) dummy, material
        READ(f_session,*) dummy, potential_type
        READ(f_session,*) ! skip row
        READ(f_session,*) dummy, is_resume_calc
        READ(f_session,*) dummy, total_mc_step
        READ(f_session,*) dummy, mc_save_step
        READ(f_session,*) dummy, in_file
        READ(f_session,*) ! skip row
        READ(f_session,*) dummy, is_multi_stage_potent
        IF (is_multi_stage_potent) THEN
          READ(f_session,*) dummy, potent_1
          READ(f_session,*) dummy, mc_step_1
          READ(f_session,*) dummy, potent_2
          READ(f_session,*) dummy, mc_step_2
        ELSE
          DO iter=1,4
            READ(f_session,*)
          END DO
        END IF
        CLOSE(f_session)

        END SUBROUTINE read_session

        SUBROUTINE read_params
! read parameters for each simulation components
        USE gupta               ,ONLY: read_gupta_param
        USE random_coord        ,ONLY: read_random_param
        USE monte               ,ONLY: read_mc_param
        USE basin_hopping       ,ONLY: read_bh_param
        IMPLICIT NONE

        CALL read_gupta_param
        CALL read_random_param
        CALL read_mc_param
        CALL read_bh_param

        END SUBROUTINE read_params
        
        SUBROUTINE printout_session_feed
! information are output in a file to allow a review of input /
! parameters to perform the same calculation again when necessary
        USE random_coord        ,ONLY: printout_random_param
        USE monte               ,ONLY: printout_mc_param
        USE basin_hopping       ,ONLY: printout_bh_param
        IMPLICIT NONE
        CALL printout_session_in
        CALL printout_random_param
        CALL printout_mc_param
        CALL printout_bh_param

        END SUBROUTINE printout_session_feed

        SUBROUTINE printout_session_in
! printout session file inputs
        USE coord_grad_ene      ,ONLY: atoms,material,in_file
        USE potential           ,ONLY: potential_type
     &    ,potent_1, potent_2
        USE directory           ,ONLY: saved_session

        IMPLICIT NONE
        CHARACTER(20)   :: session_file='saved_session_in.dat'
        INTEGER         :: f_session

        OPEN(NEWUNIT=f_session
     &    , FILE='./'//saved_session//'/'//session_file
     &    , STATUS='replace')
        WRITE(f_session,'(X, A, T33, I8)') 'number_of_atoms', atoms 
        WRITE(f_session,'(X, A, T33, 20A)') 
     &    'material_or_element', material
        WRITE(f_session,'(X, A, T33, I8)') 
     &    'potential_type',potential_type
        WRITE(f_session,'(60("*"))')
        WRITE(f_session,'(X, A, T33, L)') 'is_resume_calc'
     &  , is_resume_calc
        WRITE(f_session,'(X, A, T33, I8)')
     &    'total_MC_steps', total_mc_step
        WRITE(f_session,'(X, A, T33, I8)') 'mc_save_step', mc_save_step
        WRITE(f_session,'(X, A, T33, 20A)') 'read_in_debug', in_file
        WRITE(f_session,'(60("*"))')
        WRITE(f_session,'(X, A, T33, L)') 'is_multi_stage_poten'
     &  , is_multi_stage_potent
        IF (is_multi_stage_potent) THEN
          WRITE(f_session,'(X, A, T33, I8)') 'potent_1', potent_1
          WRITE(f_session,'(X, A, T33, I8)') 'MC_step_1', mc_step_1
          WRITE(f_session,'(X, A, T33, I8)') 'potent_2', potent_2
          WRITE(f_session,'(X, A, T33, I8)') 'MC_step_2', mc_step_2
        END IF
        CLOSE(f_session)

        END SUBROUTINE printout_session_in


        END MODULE initialise
