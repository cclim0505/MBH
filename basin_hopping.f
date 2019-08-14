        MODULE basin_hopping
        USE constants
        USE initialise          ,ONLY: coord,atoms,set_coord_to_origin
        USE random_coord        ,ONLY: max_radius, polar_2_cartesian
        CONTAINS

        SUBROUTINE bhop_move
! basin-hopping step
        IMPLICIT NONE
!       CALL random_move(coord)
        CALL get_new_radius(max_radius)
        END SUBROUTINE bhop_move


        SUBROUTINE get_new_radius(radius)
! get new radius for each structure before angular move
        IMPLICIT NONE
        REAL(KIND=SGL),INTENT(INOUT) :: radius
        REAL(KIND=DBL),DIMENSION(atoms) :: radius_array


        CALL set_coord_to_origin

        CALL calc_all_radius(coord,radius_array)
!       CALL calc_distance(coord)

        radius = REAL(MAXVAL(radius_array))

        END SUBROUTINE get_new_radius

        SUBROUTINE calc_all_radius(coord,radius_array)
! calculate distance of each atom from origin
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN) :: coord
        REAL(KIND=DBL),DIMENSION(:),INTENT(INOUT) :: radius_array
        INTEGER         ::      iter

        DO iter=1,atoms
          radius_array(iter) = NORM2(coord(:,iter))
        END DO

        END SUBROUTINE calc_all_radius

        SUBROUTINE angular_displacement
! angular displacement move for basin-hopping
        IMPLICIT NONE
        INTEGER     :: chosen_index=1
        INTEGER     :: iter
        REAL(KIND=DBL),DIMENSION(3)     :: chosen_coord

        chosen_coord = coord(:,chosen_index)
        CALL displace_angle(chosen_coord)
        coord(:,chosen_index) = chosen_coord

        END SUBROUTINE angular_displacement


        SUBROUTINE displace_angle(coord)
! random angle move for angular displacement
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(3),INTENT(INOUT)          :: coord
        REAL(KIND=SGL)          :: r, theta, phi
        REAL(KIND=SGL),DIMENSION(2)     :: random_array

        r = max_radius
        CALL RANDOM_NUMBER(random_array)
        theta = REAL(PI * random_array(1))
        phi = REAL(2.0 * PI * random_array(2))

        coord(1) = r
        coord(2) = theta
        coord(3) = phi

        CALL polar_2_cartesian(coord)
         
        END SUBROUTINE displace_angle

        SUBROUTINE random_move(coord)
! random move for basin-hopping
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(INOUT)     :: coord
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE       :: temp_coord
        REAL(KIND=SGL),DIMENSION(atoms,3)     :: random_array
        INTEGER         :: iter
        REAL(KIND=SGL)  :: radius=1.0

        IF(ALLOCATED(temp_coord)) DEALLOCATE(temp_coord)
        ALLOCATE(temp_coord(3,atoms))   ! use size function for allocation

        CALL RANDOM_NUMBER(random_array)

        DO iter=1,atoms
          temp_coord(1,iter) = radius * random_array(iter,1)
          temp_coord(2,iter) = PI * random_array(iter,2)
          temp_coord(3,iter) = 2.0 * PI * random_array(iter,3)
        END DO

        DO iter=1,atoms
          CALL polar_2_cartesian(temp_coord(:,iter))
        END DO

        coord = coord + temp_coord

        END SUBROUTINE random_move

        END MODULE basin_hopping
