        MODULE random_coord
        USE constants
        USE coord_grad_ene      ,ONLY:atoms,coord
        REAL(KIND=SGL)  ::      max_radius
        REAL(KIND=SGL)  ::      radius_ratio
        REAL(KIND=SGL)  ::      ref_radius
        REAL(KIND=SGL),PARAMETER  ::      min_atomic_dist = 0.5D0

        CONTAINS

        SUBROUTINE set_random_coord
        IMPLICIT NONE
        CALL set_max_radius
!       CALL init_random_coord(atoms,coord)
        CALL init_random_improved(atoms,coord)
        END SUBROUTINE set_random_coord

        SUBROUTINE get_random3(array)
! get an array of 3 random numbers
        USE mpi_var             ,ONLY: myid
        IMPLICIT NONE
        REAL(KIND=SGL),DIMENSION(3),INTENT(OUT)     :: array
        INTEGER,DIMENSION(1)            :: seed
        REAL                            :: real_seed
        CALL CPU_TIME(real_seed)
        seed = INT(1E8*real_seed)
        CALL RANDOM_SEED(PUT=seed+myid)
        CALL RANDOM_NUMBER(array)
        END SUBROUTINE get_random3


        SUBROUTINE get_min_dist(x_coord,iteration,min_dist)
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN)        :: x_coord
        INTEGER,INTENT(IN)                              :: iteration
        REAL(KIND=DBL),INTENT(OUT)                      :: min_dist
        REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE       :: temp_x
        REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE         :: dist
        INTEGER         :: iter


        IF (ALLOCATED(temp_x)) DEALLOCATE(temp_x)
        ALLOCATE (temp_x(3))
        IF (ALLOCATED(dist)) DEALLOCATE(dist)
        ALLOCATE(dist(iteration-1))

        DO iter=1,iteration-1
          temp_x(:) = x_coord(:,iteration) - x_coord(:,iter)
          dist(iter) = NORM2(temp_x(:))
        END DO

        min_dist = MINVAL(dist)

        END SUBROUTINE get_min_dist

        SUBROUTINE init_random_improved(natoms,x_coord)
! initialise random coordinates, with improved algorith and min_dist
! considered
        IMPLICIT NONE
        INTEGER,INTENT(IN)                             :: natoms
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(INOUT)    :: x_coord
        REAL(KIND=SGL),DIMENSION(3)                    :: random_3
        INTEGER                         :: iter

        REAL(KIND=DBL)          :: min_dist

        CALL get_random3(random_3)

        x_coord(1,1) = max_radius * random_3(1)
        x_coord(2,1) = PI * random_3(2)
        x_coord(3,1) = 2.0 * PI * random_3(3)

        CALL polar_2_cartesian(x_coord(:,1))

        iter = 2
        DO
          IF (iter > natoms) EXIT
          CALL get_random3(random_3)

          x_coord(1,iter) = max_radius * random_3(1)
          x_coord(2,iter) = PI * random_3(2)
          x_coord(3,iter) = 2.0 * PI * random_3(3)
          CALL polar_2_cartesian(x_coord(:,iter))
          CALL get_min_dist(x_coord,iter,min_dist)

          IF (min_dist < min_atomic_dist) CYCLE

          iter = iter + 1
        END DO

        END SUBROUTINE init_random_improved


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
        USE mpi_var             ,ONLY: myid
! initialise random coordinates
        IMPLICIT NONE
        INTEGER,INTENT(IN)                             :: natoms
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(INOUT)    :: x_coord
        REAL(KIND=SGL),DIMENSION(natoms,3)              :: random_array
        INTEGER         :: iter

        INTEGER,DIMENSION(1)            :: seed
        REAL                            :: real_seed

        CALL CPU_TIME(real_seed)
        seed = INT(1E8*real_seed)
        CALL RANDOM_SEED(PUT=seed+myid)
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
