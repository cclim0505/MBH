        MODULE random_coord
        USE constants           ,ONLY:SGL,DBL,PI
        USE coord_grad_ene      ,ONLY:atoms,coord,material

        REAL(KIND=SGL),PARAMETER        :: Au_ref_radius = 4.08
        REAL(KIND=SGL),PARAMETER        :: Ag_ref_radius = 4.09
        REAL(KIND=SGL),PARAMETER        :: Cu_ref_radius = 3.61
        REAL(KIND=SGL),PARAMETER        :: C_ref_radius  = 3.57

        REAL(KIND=SGL)  ::      radius_ratio
        REAL(KIND=SGL)  ::      ref_radius
        LOGICAL         ::      is_fixed_radius
        REAL(KIND=SGL)  ::      fixed_radius
        REAL(KIND=SGL),PARAMETER  ::      min_atomic_dist = 0.5D0
        REAL(KIND=SGL)  ::      max_radius

        PRIVATE          :: get_min_dist
        PRIVATE          :: init_random_improved
        PRIVATE          :: set_max_radius
        PRIVATE          :: init_random_coord
        PRIVATE          :: set_ref_radius

        PUBLIC          :: get_random3
        PUBLIC          :: set_random_coord
        PUBLIC          :: read_random_param
        PUBLIC          :: polar_2_cartesian
        PUBLIC          :: printout_random_param

        CONTAINS


        SUBROUTINE set_ref_radius
        IMPLICIT NONE

        SELECT CASE (material) 
        CASE('Au')
          ref_radius = Au_ref_radius
        CASE('Ag')
          ref_radius = Ag_ref_radius
        CASE('Cu')
          ref_radius = Cu_ref_radius
        CASE('C')
          ref_radius = C_ref_radius
        END SELECT

        END SUBROUTINE set_ref_radius


        SUBROUTINE set_random_coord
! main subroutine to set up random initial coordinates of cluster
        IMPLICIT NONE
        CALL set_max_radius
!DEBUG STARTS==============================================
!       PRINT *, "Max radius is", max_radius
!DEBUG ENDS==============================================
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
! get the minimum value of the calculated distances between newly 
! generated atom and previously generated atoms
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
! initialise random coordinates, with improved algorithm and min_dist
! considered
        IMPLICIT NONE
        INTEGER,INTENT(IN)                             :: natoms
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(INOUT)    :: x_coord
        REAL(KIND=SGL),DIMENSION(3)                    :: random_3
        INTEGER                         :: iter

        REAL(KIND=DBL)          :: min_dist

! Assign first atom
        CALL get_random3(random_3)

        x_coord(1,1) = max_radius * random_3(1)
        x_coord(2,1) = PI * random_3(2)
        x_coord(3,1) = 2.0 * PI * random_3(3)

        CALL polar_2_cartesian(x_coord(:,1))

! Assign the rest
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
! read in parameters to generate random initial coordinates
        USE directory           ,ONLY: session_dir
        IMPLICIT NONE
        CHARACTER(LEN=18)    :: random_param_file = 'param_random.dat'
        INTEGER                  :: f_random
        CHARACTER(LEN=2)         :: dummy

        OPEN(NEWUNIT=f_random
     &    ,FILE='./'//session_dir//'/'//random_param_file, STATUS='old')
        READ(f_random,*) dummy, radius_ratio
        READ(f_random,*) dummy, ref_radius
        READ(f_random,*) dummy, is_fixed_radius
        READ(f_random,*) dummy, fixed_radius
        CLOSE(f_random)
        
        END SUBROUTINE read_random_param

        SUBROUTINE printout_random_param
! printout parameters to generate random initial coordinates
        USE directory           ,ONLY: saved_session
        IMPLICIT NONE
        CHARACTER(LEN=22)    :: 
     &    random_param_file = 'saved_param_random.dat'
        INTEGER                  :: f_random
        CHARACTER(LEN=2)         :: dummy

        OPEN(NEWUNIT=f_random 
     &    , FILE='./'//saved_session//'/'//random_param_file
     &    , STATUS='replace')
        WRITE(f_random,*) 'confining_radius_ratio', radius_ratio
        WRITE(f_random,*) 'ref_radius', ref_radius
        WRITE(f_random,*) 'is_fixed_radius', is_fixed_radius
        WRITE(f_random,*) 'fixed_radius', fixed_radius
        CLOSE(f_random)
        
        END SUBROUTINE printout_random_param

        SUBROUTINE set_max_radius
! set maximum radius based on radius parameters
        IMPLICIT NONE
        !CALL set_ref_radius
        IF (is_fixed_radius) THEN
          max_radius = fixed_radius
        ELSE
          max_radius = sphere_radius(atoms, ref_radius)
          max_radius = max_radius * radius_ratio
        END IF
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

!=====================================================================
! Following functions/subroutines are not used
!=====================================================================

        SUBROUTINE init_random_coord(natoms,x_coord)
! Old random coordinates generator, not used.
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

        END MODULE random_coord
