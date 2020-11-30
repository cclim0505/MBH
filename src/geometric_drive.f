        MODULE geometric_drive
        USE constants           ,ONLY: SGL, DBL, PI
        USE coord_grad_ene      ,ONLY: atoms,coord,set_coord_to_origin
     &    ,in_file, read_atoms, read_coord
        USE random_coord        ,ONLY: get_random3,polar_2_cartesian
        USE basin_hopping       ,ONLY: random_move

        PUBLIC      :: locate_closest2centre
        PUBLIC      :: move_target2peri
        PUBLIC      :: set_ring_geometry

        PUBLIC      :: cage_move
        PUBLIC      :: multi_cage_move
        PUBLIC      :: cage_drive
        PUBLIC      :: multi_cage_drive

        PUBLIC      :: test_ring
        PUBLIC      :: test_cage1
        PUBLIC      :: test_cage2
        PUBLIC      :: test_cage3

        CONTAINS

        SUBROUTINE locate_closest2centre(in_coord 
     &      ,closest_index,largest_radius, mean_radius)
! Locate the index of the atoms that is closest to the origin.
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN)  :: in_coord
        INTEGER,INTENT(OUT)     :: closest_index
        REAL(KIND=SGL),INTENT(OUT)  :: largest_radius, mean_radius
        INTEGER         :: iter
        REAL(KIND=SGL),DIMENSION(:),ALLOCATABLE  :: radii

! Identify target atom
        IF (ALLOCATED(radii)) DEALLOCATE(radii)
        ALLOCATE(radii(atoms))

        DO iter=1,atoms
          radii(iter) = NORM2(in_coord(:,iter))
        END DO
        largest_radius = MAXVAL(radii)
        mean_radius = SUM(radii) / REAL(atoms)
        closest_index = MINLOC(radii,1)

!DEBUG BEGINS==============================================
!        PRINT '("Largest radius:", F14.7)', largest_radius
!        PRINT '("Mean radius:", F14.7)', mean_radius
!        PRINT '("Closest index:", I3)', closest_index
!DEBUG ENDS==============================================

        END SUBROUTINE locate_closest2centre

        SUBROUTINE move_target2peri(target_index,radius,in_coord)
! Move target atom to the periphery of the cluster.
        IMPLICIT NONE
        INTEGER,INTENT(IN)      :: target_index
        REAL(KIND=SGL),INTENT(IN)  :: radius
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(INOUT)  :: in_coord
        REAL(KIND=SGL),DIMENSION(3)       :: random_3
        REAL(KIND=DBL),DIMENSION(3)       :: polar_coord

        CALL get_random3(random_3)
        polar_coord(1) = radius
        polar_coord(2) = PI * random_3(2)
        polar_coord(3) = 2.0 * PI * random_3(3)

!DEBUG BEGINS==============================================
!        PRINT *, 'Random three are:', random_3
!        PRINT *, 'Polar coord are:', polar_coord
!DEBUG ENDS==============================================

        CALL polar_2_cartesian(polar_coord)

!DEBUG BEGINS==============================================
!        PRINT *, 'Polar coords after transformation:', polar_coord
!DEBUG ENDS==============================================

        in_coord(:,target_index) = polar_coord

        END SUBROUTINE move_target2peri

        SUBROUTINE set_ring_geometry(ring_radius, natoms, in_coord)
! Set ring geometry
        IMPLICIT NONE
        REAL(KIND=SGL),INTENT(IN)       :: ring_radius
        INTEGER,INTENT(IN)              :: natoms
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(INOUT)   :: in_coord
        REAL(KIND=DBL)                  :: ring_theta
        INTEGER         :: iter

        DO iter=1,natoms
          ring_theta = (2.0*PI/REAL(natoms)) * REAL(iter - 1)
          in_coord(1,iter) = ring_radius
          in_coord(2,iter) = ring_theta
          in_coord(3,iter) = 0.0D0
          CALL polar_2_cartesian(in_coord(:,iter))
        END DO
        
        END SUBROUTINE set_ring_geometry

        SUBROUTINE cage_drive
        IMPLICIT NONE
        CALL cage_move(coord)
        END SUBROUTINE cage_drive

        SUBROUTINE multi_cage_drive
        IMPLICIT NONE
        CALL multi_cage_move(6,coord)
        END SUBROUTINE multi_cage_drive

        SUBROUTINE multi_cage_move(steps,x_coord)
        IMPLICIT NONE
        INTEGER,INTENT(IN)  ::  steps
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(INOUT)     :: x_coord
        REAL(KIND=SGL)      :: cluster_radius, mean_radius
        INTEGER             :: atom_index
        INTEGER             ::  iter

        DO iter=1,steps
          CALL set_coord_to_origin(x_coord)
          CALL locate_closest2centre(x_coord, atom_index
     &      , cluster_radius, mean_radius)
          CALL move_target2peri(atom_index, cluster_radius, x_coord)

        END DO
        CALL random_move(x_coord)

        END SUBROUTINE multi_cage_move

        SUBROUTINE cage_move(x_coord)
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(INOUT)     :: x_coord
        REAL(KIND=SGL)  :: cluster_radius, mean_radius
        INTEGER         :: atom_index

        CALL set_coord_to_origin(x_coord)
        CALL locate_closest2centre(x_coord, atom_index, cluster_radius
     &    , mean_radius)
        CALL move_target2peri(atom_index, cluster_radius, x_coord)
        CALL random_move(x_coord)


!DEBUG BEGINS==============================================
!        PRINT *, "atom_index:", atom_index
!DEBUG ENDS==============================================

        END SUBROUTINE cage_move

!        SUBROUTINE ring_move(x_coord,ring_radius)
!        IMPLICIT NONE
!        REAL(KIND=DBL),DIMENSION(:,:),INTENT(INOUT)     :: x_coord
!        REAL(KIND=SGL),INTENT(IN)      :: ring_radius
!
!        END SUBROUTINE ring_move

        SUBROUTINE test_ring
! Set coordinate to desired sample.
        IMPLICIT NONE
        REAL(KIND=SGL)      :: radius, increment
        INTEGER             :: iter, jter

        radius = 0.0
        increment = 1.0

        DO iter=1,5
          radius = radius + increment
          CALL set_ring_geometry(radius,atoms,coord)
          PRINT *, 'This is ring with radius:', radius
          PRINT *, ''
          DO jter=1,atoms
            PRINT *, coord(1,jter),coord(2,jter),coord(3,jter)
          END DO
          PRINT *, ''
        END DO


        END SUBROUTINE test_ring

        SUBROUTINE test_cage1
! Set coordinate to desired sample.
! Let's use a linear 8 atom chain cluster as example.
        IMPLICIT NONE
        INTEGER         :: atom_index
        REAL(KIND=SGL)  :: cluster_radius, mean_radius
        INTEGER         :: iter

        coord(1,1) =    -4.55001573000000
        coord(1,2) =    -3.27530796000000
        coord(1,3) =    -1.96398695000000
        coord(1,4) =   -0.685393760000000
        coord(1,5) =    0.617912780000000
        coord(1,6) =     1.89651279000000
        coord(1,7) =     3.20785785000000
        coord(1,8) =     4.48259494000000

        coord(2,1) =      0.304017480000000     
        coord(2,2) =      0.184004180000000     
        coord(2,3) =      6.058044000000000E-002
        coord(2,4) =     -6.006297000000000E-002
        coord(2,5) =     -0.182575340000000     
        coord(2,6) =     -0.303000260000000     
        coord(2,7) =     -0.426463100000000     
        coord(2,8) =     -0.546393510000000     

        coord(3,1) =     -2.358769000000000E-002
        coord(3,2) =     -7.321819000000000E-002
        coord(3,3) =     -0.124204960000000     
        coord(3,4) =     -0.173281820000000     
        coord(3,5) =     -0.223303300000000     
        coord(3,6) =     -0.272300450000000     
        coord(3,7) =     -0.322091520000000     
        coord(3,8) =     -0.371105690000000     

        CALL cage_move(coord)



!        CALL set_coord_to_origin(coord)
!
!        PRINT *, "Set to origin"
!        DO iter=1,atoms
!          PRINT *, coord(1,iter), coord(2,iter), coord(3,iter)
!        END DO
!        PRINT *, " "
!
!        CALL locate_closest2centre(coord, atom_index, cluster_radius,
!     &    mean_radius)
!        CALL move_target2peri(atom_index, cluster_radius, coord)
!
!        DO iter=1,atoms
!          PRINT *, coord(1,iter), coord(2,iter), coord(3,iter)
!        END DO
!
!        CALL random_move(coord)

        DO iter=1,atoms
          PRINT *, coord(1,iter), coord(2,iter), coord(3,iter)
        END DO



        END SUBROUTINE test_cage1

        SUBROUTINE test_cage2
! Set coordinate to desired sample.
        IMPLICIT NONE
        INTEGER         :: iter

        in_file = "in.xyz"
        

        CALL read_atoms
        CALL read_coord
        CALL cage_drive

!        CALL cage_move(coord)

        DO iter=1,atoms
          PRINT *, coord(1,iter), coord(2,iter), coord(3,iter)
        END DO

        END SUBROUTINE test_cage2


        SUBROUTINE test_cage3
        IMPLICIT NONE
        INTEGER     :: iter
        in_file = "in.xyz"
        

        CALL read_atoms
        CALL read_coord
        CALL multi_cage_drive

        DO iter=1,atoms
          PRINT *, coord(1,iter), coord(2,iter), coord(3,iter)
        END DO

        END SUBROUTINE test_cage3



        END MODULE geometric_drive

!   Linear Au8 xyz coordinates
!   -4.55001573000000       0.304017480000000      -2.358769000000000E-002
!   -3.27530796000000       0.184004180000000      -7.321819000000000E-002
!   -1.96398695000000       6.058044000000000E-002 -0.124204960000000     
!  -0.685393760000000      -6.006297000000000E-002 -0.173281820000000     
!   0.617912780000000      -0.182575340000000      -0.223303300000000     
!    1.89651279000000      -0.303000260000000      -0.272300450000000     
!    3.20785785000000      -0.426463100000000      -0.322091520000000     
!    4.48259494000000      -0.546393510000000      -0.371105690000000     
