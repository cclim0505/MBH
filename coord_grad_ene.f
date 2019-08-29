        MODULE coord_grad_ene
        USE constants   ,ONLY:DBL
        INTEGER                                    :: atoms

        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE  :: coord
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE  :: old_coord
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE  :: optim_coord
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE  :: lowest_coord

        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE  :: gradient

        REAL(KIND=DBL)                             :: energy
        REAL(KIND=DBL)                             :: old_energy
        REAL(KIND=DBL)                             :: lowest_energy

        CHARACTER(LEN=18)        :: in_file = 'test2.xyz'

        CONTAINS

        SUBROUTINE read_atoms
! determine the number of atoms
        IMPLICIT NONE
        INTEGER         :: atoms
        INTEGER         :: f_atoms
        OPEN(NEWUNIT=f_atoms,FILE=TRIM(in_file), STATUS='old')
        READ(f_atoms,*) atoms
        CLOSE(f_atoms)
        END SUBROUTINE read_atoms

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

        SUBROUTINE print_ene_dat(filename)
        IMPLICIT NONE 
        CHARACTER(LEN=*),INTENT(IN)       :: filename
        INTEGER :: counter=1
        INTEGER :: f_out

        OPEN(NEWUNIT=f_out,FILE=TRIM(filename),ACCESS='append')
        WRITE(f_out,*) old_energy
        CLOSE(f_out)

        counter = counter + 1

        END SUBROUTINE print_ene_dat

        SUBROUTINE print_coord_xyz(filename)
        IMPLICIT NONE
        CHARACTER(LEN=*),INTENT(IN)       :: filename
        INTEGER :: iter
        INTEGER :: f_out

        OPEN(NEWUNIT=f_out,FILE=TRIM(filename),ACCESS='append')

        WRITE(f_out,*) atoms
        WRITE(f_out,*) 
        DO iter=1,atoms
          WRITE(f_out,*) 'Au',coord(1,iter),coord(2,iter)
     $      ,coord(3,iter)
        END DO

        CLOSE(f_out)

        END SUBROUTINE print_coord_xyz

        SUBROUTINE calc_centroid(x_coord,centroid)
! calculate centroid for set_coord_to_origin
        IMPLICIT NONE
        INTEGER :: iter
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN)    :: x_coord
        REAL(KIND=DBL),DIMENSION(3),INTENT(OUT)     :: centroid

        centroid (:) = 0.0
        DO iter=1,atoms
          centroid(1) = centroid(1) + x_coord(1,iter)
          centroid(2) = centroid(2) + x_coord(2,iter)
          centroid(3) = centroid(3) + x_coord(3,iter)
        END DO

        centroid(:) = centroid(:) / REAL(atoms)

!DEBUG==============================================
        PRINT *, 'centroid'
        PRINT *, centroid(1), centroid(2), centroid(3)
!DEBUG==============================================

        END SUBROUTINE calc_centroid

        SUBROUTINE set_coord_to_origin
! set coordinates' centroid at origin
        IMPLICIT NONE
        INTEGER                         :: iter
        REAL(KIND=DBL),DIMENSION(3)     :: centroid

        CALL calc_centroid(coord,centroid)

        DO iter=1,atoms
          coord(:,iter) = coord (:,iter) - centroid(:)
        END DO

        END SUBROUTINE set_coord_to_origin


        END MODULE coord_grad_ene
