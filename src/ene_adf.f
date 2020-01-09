        MODULE adf
        USE constants           ,ONLY:DBL
        USE coord_grad_ene      ,ONLY:coord, material, energy

        PUBLIC  :: test_adf
        PUBLIC  :: test_adf_2


        PUBLIC  :: run_adf_calc

        PUBLIC  :: gen_adf_xyz
        PUBLIC  :: read_adf_xyz
        PUBLIC  :: read_adf_energy



        CONTAINS

        SUBROUTINE test_adf
        IMPLICIT NONE

        CALL read_adf_xyz(coord)        ! read given coord
        CALL gen_adf_xyz(coord)         ! gen given coord
        CALL run_adf_calc               ! run ADF calculation

        CALL gen_adf_xyz(coord)         ! gen final_coord

        CALL read_adf_energy(energy)
        PRINT *, 'energy is', energy

        END SUBROUTINE test_adf


        SUBROUTINE test_adf_2
        IMPLICIT NONE

        END SUBROUTINE test_adf_2

        SUBROUTINE run_adf_calc
        IMPLICIT NONE
        CHARACTER(LEN=6)       :: in_file='adf.in'
        CHARACTER(LEN=7)       :: out_file='adf.out'
        CHARACTER(LEN=22)      :: run_command

        run_command = 'adf < '//in_file//' > '//out_file

        CALL EXECUTE_COMMAND_LINE(run_command)
! adf < adf_in > adf_out

        END SUBROUTINE run_adf_calc



        SUBROUTINE gen_adf_xyz(x_coord)
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN)    :: x_coord
        CHARACTER(LEN=13)               :: filename='adf_input.xyz'
        INTEGER :: natoms
        INTEGER :: iter
        INTEGER :: f_out

        natoms = SIZE(x_coord,2)

        OPEN(NEWUNIT=f_out,FILE=TRIM(filename),STATUS='replace')

        DO iter=1,natoms
          WRITE(f_out,*) material,x_coord(1,iter),x_coord(2,iter)
     &      ,x_coord(3,iter)
        END DO

        CLOSE(f_out)
        END SUBROUTINE gen_adf_xyz

        SUBROUTINE read_adf_xyz(x_coord)
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(INOUT)    :: x_coord
        CHARACTER(LEN=14)               :: filename='adf_output.xyz'
        CHARACTER(LEN=1)        :: dummy
        INTEGER :: natoms
        INTEGER :: iter
        INTEGER :: f_in

        natoms = SIZE(x_coord,2)

        OPEN(NEWUNIT=f_in,FILE=TRIM(filename),STATUS='old'
     &       ,POSITION='append')

! Go to the last entry containing the final xyz.
        DO iter=1,natoms
          BACKSPACE(f_in)
        END DO

        DO iter=1,natoms
          READ(f_in,*) dummy,x_coord(1,iter),x_coord(2,iter)
     &      ,x_coord(3,iter)
        END DO

        CLOSE(f_in)
        END SUBROUTINE read_adf_xyz

        SUBROUTINE read_adf_energy(ene)
        IMPLICIT NONE
        REAL(KIND=DBL),INTENT(OUT)      :: ene
        CHARACTER(LEN=7)  :: ori_file='adf.out'
        CHARACTER(LEN=11) :: file_in='adf_tac.out'
        CHARACTER(LEN=25) :: tac_command
        INTEGER           :: f_num
        INTEGER           :: ierr

        CHARACTER(LEN=8)                       :: dum1, dum2, dum3, dum4
        REAL(KIND=DBL)                         :: dum5
        CHARACTER(LEN=5),PARAMETER             :: match1='Total'
        CHARACTER(LEN=7),PARAMETER             :: match2='energy'

! Print output files in reverse, then look for total energy keyword.
        tac_command = 'tac '//ori_file//' > '//file_in
        CALL EXECUTE_COMMAND_LINE(tac_command)

        OPEN(NEWUNIT=f_num, FILE=file_in, STATUS='old')

        DO 
          ierr = 0
          READ(f_num,*,IOSTAT=ierr) dum1, dum2, dum3, dum4, dum5
          IF(ierr /= 0) CYCLE
          IF(dum3 == match1 .AND. dum4 == match2) THEN
            ene = dum5
            EXIT
          END IF
        END DO


        CLOSE(f_num)


        END SUBROUTINE read_adf_energy

        END MODULE adf
