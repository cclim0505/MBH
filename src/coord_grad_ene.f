        MODULE coord_grad_ene
        USE constants   ,ONLY:DBL
        INTEGER                                    :: atoms
        CHARACTER(2)                               :: material

        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE  :: coord
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE  :: old_coord
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE  :: optim_coord
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE  :: lowest_coord
        INTEGER                                    :: lowest_step

        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE  :: gradient

        REAL(KIND=DBL)                             :: energy
        REAL(KIND=DBL)                             :: old_energy
        REAL(KIND=DBL)                             :: lowest_energy

        CHARACTER(LEN=23) :: 
     &    resume_old_coord='01_resume_old_coord.xyz'
        CHARACTER(LEN=26) :: 
     &    resume_lowest_coord='01_resume_lowest_coord.xyz'
        CHARACTER(LEN=24) :: 
     &    resume_lowest_ene='01_resume_lowest_ene.dat'
        CHARACTER(LEN=24) :: local_coords='03_local_coords.xyz'
        CHARACTER(LEN=25) :: local_energies='03_local_energies.dat'

!DEBUG===============================================================
        CHARACTER(LEN=18)        :: in_file 
!DEBUG===============================================================

        PUBLIC :: read_atoms           ! for testing only
        PUBLIC :: read_coord           ! for testing only

        PRIVATE :: calc_centroid
        PRIVATE :: read_coord_ene

        PUBLIC :: allocate_coord_gradient
        PUBLIC :: printout_single_coord
        PUBLIC :: read_single_coord
        PUBLIC :: print_coord
        PUBLIC :: print_lowest_ene
        PUBLIC :: print_ene_dat
        PUBLIC :: printout_xyz
        PUBLIC :: print_update_lowest_coord
        PUBLIC :: print_lowest_coord
        PUBLIC :: print_local_coord
        PUBLIC :: set_coord_to_origin

        PUBLIC :: resume_print_lowest_coord
        PUBLIC :: resume_print_old_coord
        PUBLIC :: read_resume_coord_ene

        CONTAINS


        SUBROUTINE allocate_coord_gradient
! allocate memory for coordinates and gradient
        IMPLICIT NONE

        IF(ALLOCATED(coord)) DEALLOCATE(coord)
        IF(ALLOCATED(old_coord)) DEALLOCATE(old_coord)
        IF(ALLOCATED(optim_coord)) DEALLOCATE(optim_coord)
        IF(ALLOCATED(lowest_coord)) DEALLOCATE(lowest_coord)
        ALLOCATE(coord(3,atoms), old_coord(3,atoms)
     &    , optim_coord(3,atoms), lowest_coord(3,atoms))

        IF(ALLOCATED(gradient)) DEALLOCATE(gradient)
        ALLOCATE(gradient(3,atoms))

        END SUBROUTINE allocate_coord_gradient

        SUBROUTINE read_coord_ene(file_name,x_coord,energy)
! read input coordinates
        IMPLICIT NONE
        CHARACTER(LEN=*),INTENT(IN)                        :: file_name
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(INOUT)        :: x_coord
        REAL(KIND=DBL),INTENT(INOUT)                       :: energy
        INTEGER                 :: iter
        INTEGER                 :: f_coord
        CHARACTER(len=1)        :: dummy

        OPEN(NEWUNIT=f_coord,FILE=TRIM(file_name), STATUS='old')
        READ(f_coord,*)
        READ(f_coord,*) energy
        DO iter=1,atoms
         READ(f_coord,*) dummy
     &     ,x_coord(1,iter),x_coord(2,iter),x_coord(3,iter)
        END DO
        CLOSE(f_coord)

        END SUBROUTINE read_coord_ene

        SUBROUTINE read_resume_coord_ene
! read in old coordinates and energy to resume calculation
        IMPLICIT NONE

        CALL read_coord_ene(resume_old_coord,old_coord,old_energy)
        CALL read_coord_ene(resume_lowest_coord
     &    ,lowest_coord,lowest_energy)
        coord = old_coord

        END SUBROUTINE read_resume_coord_ene

        SUBROUTINE printout_single_coord(f_num,x_coord)
! print out coordinates to a file. Takes in file number for the output
! file
        IMPLICIT NONE
        INTEGER,INTENT(IN)                           :: f_num
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN)     :: x_coord
        INTEGER                                      :: iter

        WRITE(f_num,*) SIZE(x_coord,2)
        WRITE(f_num,*)
        DO iter=1,SIZE(x_coord,2)
          WRITE(f_num,*) material ,x_coord(1,iter),x_coord(2,iter)
     &      ,x_coord(3,iter)
        END DO

        END SUBROUTINE printout_single_coord

        SUBROUTINE read_single_coord(f_num,x_coord)
! reads xyz coordinate file, takes in file number as argument
        IMPLICIT NONE
        INTEGER,INTENT(IN)                           :: f_num
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(INOUT)  :: x_coord
        INTEGER                                      :: iter
        CHARACTER(len=1)                             :: dummy

        READ(f_num,*)
        READ(f_num,*)
        DO iter=1,SIZE(x_coord,2)
          READ(f_num,*) dummy,x_coord(1,iter),x_coord(2,iter)
     &      ,x_coord(3,iter)
        END DO

        END SUBROUTINE read_single_coord

        SUBROUTINE print_coord
! print out coordinates for checking
        IMPLICIT NONE
        INTEGER :: iter
        PRINT *, '%%%%%%%%%%%%%%%%%%'
        PRINT *, 'printing coordinates'
        DO iter=1,atoms
          PRINT *,  coord(1,iter),coord(2,iter),coord(3,iter)
        END DO
        PRINT *, 'end of coordinates'
        PRINT *, '%%%%%%%%%%%%%%%%%%'
        END SUBROUTINE print_coord

        SUBROUTINE print_lowest_ene
! print lowest energy of structure to a file
        IMPLICIT NONE 
        CHARACTER(LEN=20)      :: filename='00_lowest_energy.dat'
        INTEGER :: f_out

        OPEN(NEWUNIT=f_out,FILE=filename,STATUS='replace')
        WRITE(f_out,*) lowest_energy
        CLOSE(f_out)

        END SUBROUTINE print_lowest_ene

        SUBROUTINE print_ene_dat(filename)
! print old energies when updated to a file
        IMPLICIT NONE 
        CHARACTER(LEN=*),INTENT(IN)       :: filename
        INTEGER :: counter=1
        INTEGER :: f_out

        OPEN(NEWUNIT=f_out,FILE=TRIM(filename),ACCESS='append')
        WRITE(f_out,*) old_energy
        CLOSE(f_out)

        counter = counter + 1

        END SUBROUTINE print_ene_dat

        SUBROUTINE printout_xyz(filename,x_coord)
! printout coordinates
        IMPLICIT NONE
        CHARACTER(LEN=*),INTENT(IN)                 :: filename
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN)    :: x_coord
        INTEGER :: natoms
        INTEGER :: iter
        INTEGER :: f_out

        natoms = SIZE(x_coord,2)

        OPEN(NEWUNIT=f_out,FILE=TRIM(filename),ACCESS='append')

        WRITE(f_out,*) natoms
        WRITE(f_out,*) 
        DO iter=1,natoms
          WRITE(f_out,*) material,x_coord(1,iter),x_coord(2,iter)
     $      ,x_coord(3,iter)
        END DO

        CLOSE(f_out)

        END SUBROUTINE printout_xyz

        SUBROUTINE print_update_lowest_coord(step)
! update the lowest energy coordinate from time to time
        IMPLICIT NONE
        INTEGER,INTENT(IN)      :: step
        CHARACTER(LEN=24) :: filename='02_all_lowest_coords.xyz'
        INTEGER           :: iter
        INTEGER           :: f_out

        OPEN(NEWUNIT=f_out,FILE=TRIM(filename),ACCESS='append')

        WRITE(f_out,*) atoms
        WRITE(f_out,*) lowest_energy, "step: ", step
        DO iter=1,atoms
          WRITE(f_out,*) material,lowest_coord(1,iter),
     &      lowest_coord(2,iter) ,lowest_coord(3,iter)
        END DO

        CLOSE(f_out)

        END SUBROUTINE print_update_lowest_coord

        SUBROUTINE print_lowest_coord
! print out the lowest energy coordinate
        IMPLICIT NONE
        CHARACTER(LEN=19) :: filename='00_lowest_coord.xyz'
        INTEGER           :: iter
        INTEGER           :: f_out

        OPEN(NEWUNIT=f_out,FILE=TRIM(filename),STATUS='replace')

        WRITE(f_out,*) atoms
        WRITE(f_out,*) lowest_energy, "step: ", lowest_step
        DO iter=1,atoms
          WRITE(f_out,*) material,lowest_coord(1,iter),
     &      lowest_coord(2,iter) ,lowest_coord(3,iter)
        END DO

        CLOSE(f_out)

        END SUBROUTINE print_lowest_coord

        SUBROUTINE resume_print_lowest_coord(step)
! print out the lowest energy coordinate
        IMPLICIT NONE
        INTEGER,INTENT(IN)      :: step
        CHARACTER(LEN=26) :: filename
        INTEGER           :: iter
        INTEGER           :: f_out

        filename = resume_lowest_coord
        OPEN(NEWUNIT=f_out,FILE=TRIM(filename)
     &    ,STATUS='replace')

        WRITE(f_out,*) atoms
        WRITE(f_out,*) lowest_energy, "step: ", step
        DO iter=1,atoms
          WRITE(f_out,*) material,lowest_coord(1,iter),
     &      lowest_coord(2,iter) ,lowest_coord(3,iter)
        END DO

        CLOSE(f_out)

        filename = resume_lowest_ene
        OPEN(NEWUNIT=f_out,FILE=TRIM(filename)
     &    ,STATUS='replace')
        WRITE(f_out,*) lowest_energy
        CLOSE(f_out)

        END SUBROUTINE resume_print_lowest_coord

        SUBROUTINE resume_print_old_coord
! print out the lowest energy coordinate
        IMPLICIT NONE
        CHARACTER(LEN=26) :: filename
        INTEGER           :: iter
        INTEGER           :: f_out

        filename = resume_old_coord
        OPEN(NEWUNIT=f_out,FILE=TRIM(filename)
     &    ,STATUS='replace')

        WRITE(f_out,*) atoms
        WRITE(f_out,*) old_energy
        DO iter=1,atoms
          WRITE(f_out,*) material,old_coord(1,iter),
     &      old_coord(2,iter) ,old_coord(3,iter)
        END DO

        CLOSE(f_out)

        END SUBROUTINE resume_print_old_coord

        SUBROUTINE print_local_coord
! print local minima coordinates after each MC step
        IMPLICIT NONE
        CHARACTER(LEN=24) :: filename
        INTEGER           :: iter
        INTEGER           :: f_out

        filename=local_coords
        OPEN(NEWUNIT=f_out,FILE=TRIM(filename),ACCESS='append')

        WRITE(f_out,*) atoms
        WRITE(f_out,*) energy
        DO iter=1,atoms
          WRITE(f_out,*) material,coord(1,iter),
     &      coord(2,iter) ,coord(3,iter)
        END DO

        CLOSE(f_out)

        filename=local_energies
        OPEN(NEWUNIT=f_out,FILE=TRIM(filename),ACCESS='append')
        WRITE(f_out,*) energy
        CLOSE(f_out)

        END SUBROUTINE print_local_coord

        SUBROUTINE calc_centroid(x_coord,centroid)
! calculate centroid for set_coord_to_origin
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN)    :: x_coord
        REAL(KIND=DBL),DIMENSION(3),INTENT(OUT)     :: centroid
        INTEGER                         :: iter
        INTEGER                         :: natoms

        natoms = SIZE(x_coord,2)

        centroid (:) = 0.0
        DO iter=1,natoms
          centroid(1) = centroid(1) + x_coord(1,iter)
          centroid(2) = centroid(2) + x_coord(2,iter)
          centroid(3) = centroid(3) + x_coord(3,iter)
        END DO

        centroid(:) = centroid(:) / REAL(atoms)

!DEBUG==============================================
!       PRINT *, 'centroid'
!       PRINT *, centroid(1), centroid(2), centroid(3)
!DEBUG==============================================

        END SUBROUTINE calc_centroid

        SUBROUTINE set_coord_to_origin(x_coord)
! set coordinates' centroid at origin
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(INOUT)    :: x_coord
        INTEGER                         :: iter
        INTEGER                         :: natoms
        REAL(KIND=DBL),DIMENSION(3)     :: centroid

        natoms = SIZE(x_coord,2)

        CALL calc_centroid(x_coord,centroid)

        DO iter=1,natoms
          x_coord(:,iter) = x_coord (:,iter) - centroid(:)
        END DO

        END SUBROUTINE set_coord_to_origin

        SUBROUTINE read_atoms
! determine the number of atoms
        IMPLICIT NONE
        INTEGER         :: f_atoms
        OPEN(NEWUNIT=f_atoms,FILE=TRIM(in_file), STATUS='old')
        READ(f_atoms,*) atoms
        CLOSE(f_atoms)
        END SUBROUTINE read_atoms

        SUBROUTINE read_coord
! read input coordinates
        IMPLICIT NONE
        INTEGER                 :: iter
        INTEGER                 :: f_coord
        CHARACTER(len=1)        :: dummy

        OPEN(NEWUNIT=f_coord,FILE=TRIM(in_file), STATUS='old')
        READ(f_coord,*)
        READ(f_coord,*)
        DO iter=1,atoms
         READ(f_coord,*) dummy,coord(1,iter),coord(2,iter),coord(3,iter)
        END DO
        CLOSE(f_coord)

        END SUBROUTINE read_coord


        END MODULE coord_grad_ene
