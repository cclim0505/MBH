        MODULE initialise

        CHARACTER(2)               :: material
        CHARACTER(6)               :: potential_type
        INTEGER                    :: total_mc_step
        
        CONTAINS

        SUBROUTINE read_session
! read session simulation parameters
        USE coord_grad_ene,     ONLY: atoms
        IMPLICIT NONE
        CHARACTER(14)   :: session_file='session_in.dat'
        INTEGER         :: f_session
        CHARACTER       :: dummy

        OPEN(NEWUNIT=f_session,FILE=session_file,STATUS='old')
        READ(f_session,*) dummy, atoms       
        READ(f_session,*) dummy, material
        READ(f_session,*) dummy, potential_type
        READ(f_session,*) dummy, total_mc_step
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
