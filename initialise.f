        MODULE initialise

        INTEGER             :: total_mc_step
        INTEGER,PARAMETER   :: mc_save_step = 10
        
        CONTAINS

        SUBROUTINE check_mc_step(iteration)
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
        USE coord_grad_ene,     ONLY: atoms,material,in_file
        USE potential,          ONLY: potential_type
        IMPLICIT NONE
        CHARACTER(14)   :: session_file='session_in.dat'
        INTEGER         :: f_session
        CHARACTER       :: dummy

        OPEN(NEWUNIT=f_session,FILE=session_file,STATUS='old')
        READ(f_session,*) dummy, atoms       
        READ(f_session,*) dummy, material
        READ(f_session,*) dummy, potential_type
        READ(f_session,*) dummy, total_mc_step
        READ(f_session,*) dummy, in_file
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
        
        END MODULE initialise
