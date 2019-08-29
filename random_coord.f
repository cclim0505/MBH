        MODULE random_coord
        USE constants
        USE coord_grad_ene      ,ONLY:atoms,coord
        REAL(KIND=SGL)  ::      max_radius
        REAL(KIND=SGL)  ::      radius_ratio
        REAL(KIND=SGL)  ::      ref_radius

        CONTAINS

        SUBROUTINE set_random_coord
        IMPLICIT NONE
        CALL set_max_radius
        CALL init_random_coord(atoms,coord)
        END SUBROUTINE set_random_coord


        SUBROUTINE read_random_param
        IMPLICIT NONE
        CHARACTER(LEN=18)    :: random_param_file = 'param_random.dat'
        INTEGER                  :: f_random
        CHARACTER(LEN=2)         :: dummy

        OPEN(NEWUNIT=f_random, FILE=random_param_file, STATUS='old')
        READ(f_random,*) dummy, radius_ratio
        READ(f_random,*) dummy, ref_radius
        CLOSE(f_random)
        
        END SUBROUTINE read_random_param
! add adjust centroid subroutine here
!       SUBROUTINE centroid

        SUBROUTINE set_max_radius
! set maximum radius based on radius parameters
        IMPLICIT NONE
        max_radius = sphere_radius(atoms, ref_radius)
        max_radius = max_radius * radius_ratio
        END SUBROUTINE set_max_radius

        REAL FUNCTION sphere_radius(num,r_zero)
! calculate sphere radius based on formula
        IMPLICIT NONE
        INTEGER,INTENT(IN)        :: num        ! number of atoms
        REAL(KIND=SGL),INTENT(IN) :: r_zero     ! nearest neighbour distance
        REAL(KIND=SGL),PARAMETER  :: one_third = 1.0 / 3.0

        sphere_radius = (3.0 * num) / (4.0 * REAL(PI) * SQRT(2.0))
        sphere_radius = sphere_radius ** one_third
        sphere_radius = r_zero * (1.0 + sphere_radius)

        END FUNCTION sphere_radius

        SUBROUTINE init_random_coord(natoms,x_coord)
! initialise random coordinates
        IMPLICIT NONE
        INTEGER,INTENT(IN)                             :: natoms
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(INOUT)    :: x_coord
        REAL(KIND=SGL),DIMENSION(natoms,3)              :: random_array
        INTEGER         :: iter

        CALL RANDOM_NUMBER(random_array)

        DO iter=1,natoms
          x_coord(1,iter) = max_radius * random_array(iter,1)
          x_coord(2,iter) = PI * random_array(iter,2)
          x_coord(3,iter) = 2.0 * PI * random_array(iter,3)
        END DO

        DO iter=1,natoms
          CALL polar_2_cartesian(x_coord(:,iter))
        END DO

!DEBUG STARTS==============================================
!       PRINT *, 'x_coord after polar_2_cartesian'
!       DO iter=1,natoms
!         PRINT *, x_coord(1,iter), x_coord(2,iter), x_coord(3,iter)
!       END DO
!DEBUG ENDS==============================================

        END SUBROUTINE init_random_coord

        SUBROUTINE polar_2_cartesian(x_coord)
! change polar coordinates to cartesian coordinates
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(3),INTENT(INOUT)  :: x_coord
        REAL(KIND=SGL)    :: r, theta, phi
        REAL(KIND=SGL)    :: x, y, z

        r = REAL(x_coord(1))
        theta = REAL(x_coord(2))
        phi = REAL(x_coord(3))

        x = r*SIN(theta)*COS(phi)
        y = r*SIN(theta)*SIN(phi)
        z = r*COS(theta)
        
        x_coord(1) = x
        x_coord(2) = y
        x_coord(3) = z
        END SUBROUTINE polar_2_cartesian

        END MODULE random_coord
