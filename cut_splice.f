        MODULE cut_splice
        USE constants           ,ONLY:DBL
        USE coord_grad_ene      ,ONLY:set_coord_to_origin
     &    ,printout_xyz
        CONTAINS

        SUBROUTINE rotate_upperlower_cut(upper_cut,lower_cut,phi
     &    ,is_upper)
        USE inertia             ,ONLY:rotate_anticlock
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT) 
     &    :: upper_cut, lower_cut
        REAL(KIND=DBL),INTENT(IN)          :: phi
        LOGICAL,INTENT(IN)                 :: is_upper

        IF (is_upper) THEN
          CALL rotate_anticlock(1,phi,upper_cut, .TRUE.)
        ELSE
          CALL rotate_anticlock(1,phi,lower_cut, .TRUE.)
        END IF

        END SUBROUTINE rotate_upperlower_cut

        SUBROUTINE init_upper_lower_cut(x_coord,cut_point)
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN)        :: x_coord
        INTEGER,INTENT(OUT)     :: cut_point
        INTEGER         :: natoms

        natoms = SIZE(x_coord,2)
        cut_point = natoms / 2

        END SUBROUTINE init_upper_lower_cut

        SUBROUTINE cut_cluster(x_coord,cut_point,upper_cut,lower_cut)
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN)        :: x_coord
        INTEGER,INTENT(IN)         :: cut_point
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT) 
     &    :: upper_cut, lower_cut
        INTEGER         :: natoms
        INTEGER         :: upper_size, lower_size
        INTEGER         :: iter, jter

        natoms = SIZE(x_coord,2)
        upper_size = cut_point
        lower_size = natoms - cut_point


        DO iter=1,upper_size
          upper_cut(:,iter) = x_coord(:,iter)
        END DO

        jter = 1
        DO iter=upper_size + 1, natoms
          lower_cut(:,jter) = x_coord(:,iter)
          jter = jter + 1
        END DO

!DEBUG BEGINS==============================================
        PRINT *, 'upper size and lower size', upper_size, lower_size
        PRINT *, 'upper cut coordinates'
        DO iter=1,upper_size
          PRINT *, iter, upper_cut(1,iter)
     &      , upper_cut(2,iter), upper_cut(3,iter)
        END DO


        PRINT *, 'lower cut coordinates'
        DO iter=1,lower_size
          PRINT *, iter, lower_cut(1,iter)
     &      , lower_cut(2,iter), lower_cut(3,iter)
        END DO

        CALL printout_xyz('xc_coord.xyz',x_coord)
        CALL printout_xyz('uppercut.xyz',upper_cut)
        CALL printout_xyz('lowercut.xyz',lower_cut)
!DEBUG ENDS==============================================
! 
        END SUBROUTINE cut_cluster

        SUBROUTINE calc_distance_from_furthest(furthest,x_coord
     &    ,distances)
        IMPLICIT NONE
        INTEGER,INTENT(IN)                              :: furthest
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(INOUT)        :: x_coord
        REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE,INTENT(OUT) :: distances
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE     :: temp_coord
        REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE         :: temp_dist
        INTEGER,DIMENSION(:),ALLOCATABLE     :: loc_array
        INTEGER,DIMENSION(1)         :: temp
        INTEGER         :: iter,natoms

        natoms = SIZE(x_coord,2)

        IF (ALLOCATED(distances)) DEALLOCATE(distances)
        IF (ALLOCATED(temp_coord)) DEALLOCATE (temp_coord)
        IF (ALLOCATED(temp_dist)) DEALLOCATE (temp_dist)
        IF (ALLOCATED(loc_array)) DEALLOCATE (loc_array)
        ALLOCATE(distances(natoms))
        ALLOCATE(temp_dist(natoms))
        ALLOCATE(loc_array(natoms))
        ALLOCATE(temp_coord(SIZE(x_coord,1),SIZE(x_coord,2)))

        DO iter=1,natoms
          temp_coord(:,iter) = x_coord(:,furthest) - x_coord(:,iter)
          distances(iter) = NORM2(temp_coord(:,iter))
        END DO

        DO iter=1,natoms 
          temp_dist(iter) = MINVAL(distances)
          temp = MINLOC(distances)
          loc_array(iter) = temp(1)
          distances(temp) = 1.0E6
        END DO

!DEBUG BEGINS==============================================
        PRINT *, 'loc_array is'
        DO iter=1,natoms
          PRINT *, loc_array(iter)
        END DO

        PRINT *, 'coordinates before sort'
        DO iter=1,natoms
          PRINT *, iter, x_coord(1,iter)
     &      , x_coord(2,iter), x_coord(3,iter)
        END DO
!DEBUG ENDS==============================================

        DO iter=1,natoms
          temp_coord(:,iter) = x_coord(:,loc_array(iter))
        END DO

        x_coord = temp_coord

!DEBUG BEGINS==============================================
        PRINT *, 'coordinates after sort'
        DO iter=1,natoms
          PRINT *, iter, x_coord(1,iter)
     &      , x_coord(2,iter), x_coord(3,iter)
        END DO

!DEBUG ENDS==============================================

! Sort distances in ascending order

        END SUBROUTINE calc_distance_from_furthest


        SUBROUTINE get_furthest_atom(x_coord,furthest)
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN)        :: x_coord
        INTEGER,INTENT(OUT)                             :: furthest
        REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE         ::  lengths
        REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE         :: temp_len
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE       :: temp_coord

        INTEGER,DIMENSION(:),ALLOCATABLE     :: loc_array
        INTEGER,DIMENSION(1)         :: temp
        INTEGER         :: natoms
        INTEGER         :: iter
        
        natoms = SIZE(x_coord,2)
        IF (ALLOCATED(lengths)) DEALLOCATE (lengths)
        IF (ALLOCATED(temp_len)) DEALLOCATE (temp_len)
        IF (ALLOCATED(loc_array)) DEALLOCATE (loc_array)
        IF (ALLOCATED(temp_coord)) DEALLOCATE (temp_coord)

        ALLOCATE(lengths(natoms),temp_len(natoms),loc_array(natoms))
        ALLOCATE(temp_coord(SIZE(x_coord,1),SIZE(x_coord,2)))

        DO iter=1,natoms
          lengths(iter) = NORM2(x_coord(:,iter))
        END DO

        temp = MAXLOC(lengths) 
        furthest = temp(1)

!DEBUG BEGINS==============================================
        PRINT *, 'furthest is', furthest
!       PRINT *, 'lengths are='
!       DO iter=1,natoms
!         PRINT *, lengths(iter)
!       END DO
!DEBUG ENDS==============================================

! Sort lengths in descending order
!       DO iter=1,natoms 
!         temp_len(iter) = MAXVAL(lengths)
!         temp = MAXLOC(lengths)
!         loc_array(iter) = temp(1)
!         lengths(temp) = 0.0D0
!       END DO

!       lengths = temp_len

!DEBUG BEGINS==============================================
!       PRINT *, 'sorted lengths are='
!       DO iter=1,natoms
!         PRINT *, lengths(iter)
!       END DO

!       PRINT *, 'loc_array is'
!       DO iter=1,natoms
!         PRINT *, loc_array(iter)
!       END DO
!DEBUG ENDS==============================================

        END SUBROUTINE get_furthest_atom


        END MODULE cut_splice
