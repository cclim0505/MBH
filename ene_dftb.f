        MODULE dftb
        USE constants           ,ONLY:DBL
        USE coord_grad_ene      ,ONLY:coord,gradient,atoms

        REAL(KIND=DBL),PARAMETER        :: au_r0 = 2.884
        REAL(KIND=DBL),PARAMETER        :: au_e0 = 27.21138624598853
        REAL(KIND=DBL),PARAMETER        :: au_f0 = 51.42208619083232
        CHARACTER(LEN=18)          :: run_dftb = 'dftb+19 > dftb.log'

        CONTAINS

        SUBROUTINE indv_dftb_energy(natoms,ene_array)
        IMPLICIT NONE
        INTEGER,INTENT(IN)                            :: natoms
        REAL(KIND=DBL),DIMENSION(natoms),INTENT(OUT)  :: ene_array

        CALL read_indv_dftb_energy(natoms,ene_array)

        END SUBROUTINE indv_dftb_energy

        SUBROUTINE read_indv_dftb_energy(natoms,ene_array)
        IMPLICIT NONE
        INTEGER,INTENT(IN)                     :: natoms
        REAL(KIND=DBL),DIMENSION(:),INTENT(OUT) :: ene_array

        INTEGER                                :: iter
        INTEGER                                :: f_num
        CHARACTER(LEN=12)                      :: file_in='detailed.out'
        CHARACTER(LEN=8)                       :: dum1, dum2, dum3, dum4
        CHARACTER(LEN=1)                       :: dummy
        CHARACTER(LEN=4),PARAMETER             :: match1='Atom'
        CHARACTER(LEN=8),PARAMETER             :: match2='resolved'
        CHARACTER(LEN=5),PARAMETER             :: match3='total'
        CHARACTER(LEN=8),PARAMETER             :: match4='energies'
        INTEGER                                :: ierr

        OPEN(NEWUNIT=f_num,FILE=file_in,STATUS='old')

        DO
          ierr = 0
          READ(f_num,*,IOSTAT=ierr) dum1, dum2, dum3, dum4
          IF(ierr /= 0) CYCLE
          IF(dum1 == match1 .AND. dum2 == match2
     &      .AND. dum3 == match3 .AND. dum4 == match4) EXIT
        END DO

        DO iter=1,natoms
          READ(f_num,*) dummy, ene_array(iter)
        END DO

        CLOSE(f_num)

!DEBUG BEGINS==============================================
!       DO iter=1,natoms
!         PRINT *, ene_array(iter)
!       END DO
!DEBUG ENDS==============================================

        END SUBROUTINE read_indv_dftb_energy

        SUBROUTINE dftb_both_ene_grad(x_coord,natoms,ene,grad)
        IMPLICIT NONE
        REAl(KIND=DBL),DIMENSION(:,:),INTENT(IN)  :: x_coord
        INTEGER,INTENT(IN)                        :: natoms
        REAL(KIND=DBL),INTENT(OUT)                :: ene
        REAl(KIND=DBL),DIMENSION(:,:),INTENT(OUT) :: grad

        CALL coord_2_gen(x_coord,natoms)
        CALL EXECUTE_COMMAND_LINE(run_dftb)
        CALL read_dftb_energy(ene)
        CALL read_dftb_gradient(natoms,grad)

        END SUBROUTINE dftb_both_ene_grad

        SUBROUTINE dftb_energy(x_coord,natoms,ene)
        IMPLICIT NONE
        REAl(KIND=DBL),DIMENSION(:,:),INTENT(IN) :: x_coord
        INTEGER,INTENT(IN)                       :: natoms
        REAL(KIND=DBL),INTENT(OUT)               :: ene

        CALL coord_2_gen(x_coord,natoms)
        CALL EXECUTE_COMMAND_LINE(run_dftb)
        CALL read_dftb_energy(ene)

        END SUBROUTINE dftb_energy

        SUBROUTINE dftb_gradient(x_coord,natoms,grad)
        IMPLICIT NONE
        REAl(KIND=DBL),DIMENSION(:,:),INTENT(IN)  :: x_coord
        INTEGER,INTENT(IN)                        :: natoms
        REAl(KIND=DBL),DIMENSION(:,:),INTENT(OUT) :: grad

        CALL coord_2_gen(x_coord,natoms)
        CALL EXECUTE_COMMAND_LINE(run_dftb)
        CALL read_dftb_gradient(natoms,grad)

        END SUBROUTINE dftb_gradient

        SUBROUTINE coord_2_gen(x_coord,natoms)
        USE coord_grad_ene          ,ONLY: material
        IMPLICIT NONE
        REAl(KIND=DBL),DIMENSION(:,:),INTENT(IN)  :: x_coord
        INTEGER,INTENT(IN)   :: natoms
        CHARACTER(LEN=14)    :: gen_file = 'dftb.gen'
        INTEGER              :: iter
        INTEGER              :: f_gen
        CHARACTER(LEN=2)     :: geometry = 'C'

        OPEN(NEWUNIT=f_gen, FILE=gen_file, STATUS='replace')

        WRITE (f_gen,*) natoms,"   " , geometry
        WRITE (f_gen,*) material

        DO iter=1,natoms
           WRITE (f_gen,*) iter, 1, x_coord(1,iter), x_coord(2,iter),
     &       x_coord(3,iter)
        END DO
        CLOSE(f_gen)

        END SUBROUTINE coord_2_gen

        SUBROUTINE read_dftb_energy(ene)
        IMPLICIT NONE
        REAL(KIND=DBL),INTENT(OUT)      :: ene
        CHARACTER(LEN=12) :: file_in='detailed.out'
        INTEGER           :: f_num
        INTEGER           :: ierr

        CHARACTER(LEN=8)                       :: dum1, dum2
        REAL(KIND=DBL)                         :: dum3
        CHARACTER(LEN=5),PARAMETER             :: match1='Total'
        CHARACTER(LEN=7),PARAMETER             :: match2='energy:'

        OPEN(NEWUNIT=f_num, FILE=file_in, STATUS='old')

        DO 
          ierr = 0
          READ(f_num,*,IOSTAT=ierr) dum1, dum2, dum3 
          IF(ierr /= 0) CYCLE
          IF(dum1 == match1 .AND. dum2 == match2) THEN
            ene = dum3
            EXIT
          END IF
        END DO


        CLOSE(f_num)

!DEBUG BEGINS==============================================
!       PRINT *, 'dftb energy in hartree', ene
!DEBUG ENDS==============================================

        ene = ene * au_e0

        END SUBROUTINE read_dftb_energy

!       SUBROUTINE read_dftb_energy(ene)
!       IMPLICIT NONE
!       REAL(KIND=DBL),INTENT(OUT)      :: ene
!       CHARACTER(LEN=11) :: res_file='results.tag'
!       INTEGER           :: f_res
!       INTEGER           :: iter
!       INTEGER           :: nline=1

!       OPEN(NEWUNIT=f_res, FILE=res_file, STATUS='old')
!       DO iter=1,nline
!         READ(f_res,*)         !skip
!       END DO
!       READ(f_res,*) ene
!       CLOSE(f_res)

!DEBUG BEGINS==============================================
!       PRINT *, 'dftb energy in hartree', ene
!DEBUG ENDS==============================================

!       ene = ene * au_e0

!       END SUBROUTINE read_dftb_energy

        SUBROUTINE read_dftb_gradient(natoms,grad)
        IMPLICIT NONE
        INTEGER,INTENT(IN)              :: natoms
        REAL(KIND=DBL),DIMENSION(3,natoms),INTENT(OUT)      :: grad

        CHARACTER(LEN=12) :: file_in='detailed.out'
        INTEGER           :: f_num
        INTEGER           :: ierr
        INTEGER           :: iter
        INTEGER           :: dummy

        CHARACTER(LEN=8)                       :: dum1, dum2
        CHARACTER(LEN=5),PARAMETER             :: match1='Total'
        CHARACTER(LEN=6),PARAMETER             :: match2='Forces'

        OPEN(NEWUNIT=f_num, FILE=file_in, STATUS='old')

        DO 
          ierr = 0
          READ(f_num,*,IOSTAT=ierr) dum1, dum2
          IF(ierr /= 0) CYCLE
          IF(dum1 == match1 .AND. dum2 == match2) EXIT
        END DO

        DO iter=1,natoms
          READ(f_num,*) dummy, grad(1,iter), grad(2,iter), grad(3,iter)
        END DO

        CLOSE(f_num)

!DEBUG BEGINS==============================================
!       PRINT *, 'read in force is'
!       DO iter=1,natoms
!         PRINT *, grad(1,iter), grad(2,iter), grad(3,iter)
!       END DO
!DEBUG ENDS==============================================

        grad = -grad * au_f0    ! add negative sign to change from force to gradient

        END SUBROUTINE read_dftb_gradient

!       SUBROUTINE read_dftb_gradient(natoms,grad)
!       IMPLICIT NONE
!       INTEGER,INTENT(IN)              :: natoms
!       REAL(KIND=DBL),DIMENSION(3,natoms),INTENT(OUT)      :: grad
!       CHARACTER(LEN=11) :: res_file='results.tag'
!       INTEGER           :: f_res
!       INTEGER           :: iter
!       INTEGER           :: nline

!       nline = 9               ! defaul guess

!==============================================
! Determining the number of lines to skip in restuls.tag
!       nline = 6 + (natoms/3) 

!       IF(MOD(natoms,3) == 0) THEN
!         nline = nline
!       ELSE
!         nline = nline + 1
!       END IF
!==============================================

!DEBUG BEGINS==============================================
!       PRINT *, 'nline is', nline
!DEBUG ENDS==============================================

!       OPEN(NEWUNIT=f_res, FILE=res_file, STATUS='old')
!       DO iter=1,nline
!         READ(f_res,*)         !skip
!       END DO

!       DO iter=1,natoms
!         READ(f_res,*) grad(1,iter), grad(2,iter), grad(3,iter)
!       END DO

!DEBUG BEGINS==============================================
!       PRINT *, 'read in force is'
!       DO iter=1,natoms
!         PRINT *, grad(1,iter), grad(2,iter), grad(3,iter)
!       END DO
!DEBUG ENDS==============================================

!       CLOSE(f_res)
!       grad = -grad * au_f0    ! add negative sign to change from force to gradient

!       END SUBROUTINE read_dftb_gradient

        END MODULE dftb
