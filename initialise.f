        MODULE initialise
        USE constants           ,ONLY: DBL
        INTEGER                                    :: atoms
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE  :: coord
        REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE    :: posit
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE  :: gradient
        CHARACTER(LEN=18)        :: in_file = 'test2.xyz'
        CONTAINS

!       SUBROUTINE read_session
!       IMPLICIT NONE
!       CHARACTER(LEN=14)       :: session_in = 'session_in.dat'
!       things to read
!       material type
!       number of atoms
!       END SUBROUTINE read_session

        SUBROUTINE init_gupta
! read gupta potential parameters
        USE gupta
        IMPLICIT NONE
        OPEN (20,FILE='AuAu_parameter1.dat',STATUS='old')
        READ(20,*) a_ij, eta, p_ij, q_ij, r_zero
        CLOSE(20)
        END SUBROUTINE init_gupta

        SUBROUTINE read_atoms
! determine the number of atoms
        IMPLICIT NONE
        OPEN(20,FILE=TRIM(in_file), STATUS='old')
        READ(20,*) atoms
        CLOSE(20)
        END SUBROUTINE read_atoms

        SUBROUTINE read_coord
! read input coordinates
        IMPLICIT NONE
        INTEGER :: iter
        CHARACTER(len=1)        :: dummy
        IF(ALLOCATED(coord)) DEALLOCATE(coord)
        ALLOCATE(coord(3,atoms))
        IF(ALLOCATED(gradient)) DEALLOCATE(gradient)
        ALLOCATE(gradient(3,atoms))
        OPEN(20,FILE=TRIM(in_file), STATUS='old')
        READ(20,*)
        READ(20,*)
        DO iter=1,atoms
          READ(20,*) dummy,coord(1,iter),coord(2,iter),coord(3,iter)
        END DO
        CLOSE(20)
        END SUBROUTINE read_coord

        SUBROUTINE convert_coord_to_x(coord,posit)
! convert coordinates from matrix form to array form
        IMPLICIT NONE
        INTEGER :: iter,jter
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN) :: coord
        REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE,INTENT(OUT) :: posit
        IF (ALLOCATED(posit)) DEALLOCATE(posit)
        ALLOCATE(posit(atoms*3))
        DO iter=1,atoms
          jter = iter*3
          posit(jter-2) = coord(1,iter)
          posit(jter-1) = coord(2,iter)
          posit(jter  ) = coord(3,iter)
        END DO

!       PRINT *,
!       PRINT *, 'POSIT values'
!       DO iter=1,atoms*3
!         PRINT *, posit(iter)
!       END DO
        END SUBROUTINE convert_coord_to_x

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

        SUBROUTINE calc_centroid(coord,centroid)
! calculate centroid for set_coord_to_origin
        IMPLICIT NONE
        INTEGER :: iter
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN)    :: coord
        REAL(KIND=DBL),DIMENSION(3),INTENT(OUT)     :: centroid

        centroid (:) = 0.0
        DO iter=1,atoms
          centroid(1) = centroid(1) + coord(1,iter)
          centroid(2) = centroid(2) + coord(2,iter)
          centroid(3) = centroid(3) + coord(3,iter)
        END DO

        centroid(:) = centroid(:) / REAL(atoms)

        PRINT *, 'centroid'
        PRINT *, centroid(1), centroid(2), centroid(3)

        END SUBROUTINE calc_centroid

        SUBROUTINE set_coord_to_origin
! set coordinates' centroid at origin
        IMPLICIT NONE
        INTEGER :: iter
        REAL(KIND=DBL),DIMENSION(3)     :: centroid
        CALL calc_centroid(coord,centroid)
        DO iter=1,atoms
          coord(:,iter) = coord (:,iter) - centroid(:)
        END DO
        END SUBROUTINE set_coord_to_origin

        END MODULE initialise
