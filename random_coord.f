        MODULE random_coord
        USE constants
        USE initialise, ONLY:atoms,coord
        REAL(KIND=SGL)  ::      max_radius
        CONTAINS

        SUBROUTINE set_random_coord
        IMPLICIT NONE
        CALL get_radius_param
        CALL init_random_coord(atoms,coord)
        END SUBROUTINE set_random_coord

! add adjust centroid subroutine here
!       SUBROUTINE centroid

        SUBROUTINE get_radius_param
        IMPLICIT NONE
        REAL(KIND=SGL)  ::      radius_ratio
        REAL(KIND=SGL)  ::      ref_radius
        radius_ratio = 0.88
        ref_radius = 4.0782     ! in angstroms for fcc gold
        max_radius = sphere_radius(atoms, ref_radius)
        max_radius = max_radius * radius_ratio
        END SUBROUTINE get_radius_param

        REAL FUNCTION sphere_radius(num,r_zero)
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: num        ! number of atoms
        REAL(KIND=SGL),INTENT(IN) :: r_zero     ! nearest neighbour distance
        REAL(KIND=SGL),PARAMETER  :: one_third = 1.0 / 3.0
        sphere_radius = (3.0 * num) / (4.0 * REAL(PI) * SQRT(2.0))
        sphere_radius = sphere_radius ** one_third
        sphere_radius = r_zero * (1.0 + sphere_radius)
        END FUNCTION sphere_radius

        SUBROUTINE init_random_coord(atoms,coord)
        IMPLICIT NONE
        INTEGER,INTENT(IN)                             :: atoms
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(INOUT)    :: coord
        REAL(KIND=SGL),DIMENSION(atoms,3)              :: random_array
        INTEGER         :: iter
        CALL RANDOM_NUMBER(random_array)
        PRINT *, random_array

        DO iter=1,atoms
          coord(1,iter) = max_radius * random_array(iter,1)
          coord(2,iter) = PI * random_array(iter,2)
          coord(3,iter) = 2.0 * PI * random_array(iter,3)
        END DO

        DO iter=1,atoms
          CALL polar_2_cartesian(coord(:,iter))
        END DO

        PRINT *, 'coord after polar_2_cartesian'
        DO iter=1,atoms
          PRINT *, coord(1,iter), coord(2,iter), coord(3,iter)
        END DO

        END SUBROUTINE init_random_coord

        SUBROUTINE polar_2_cartesian(coord)
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(3),INTENT(INOUT)  :: coord
        REAL(KIND=SGL)    :: r, theta, phi
        REAL(KIND=SGL)    :: x, y, z
        r = REAL(coord(1))
        theta = REAL(coord(2))
        phi = REAL(coord(3))

        x = r*SIN(theta)*COS(phi)
        y = r*SIN(theta)*SIN(phi)
        z = r*COS(theta)
        
        coord(1) = x
        coord(2) = y
        coord(3) = z
        END SUBROUTINE polar_2_cartesian

        END MODULE random_coord
