        MODULE basin_hopping
        USE constants           ,ONLY: SGL,DBL,PI
        USE coord_grad_ene      ,ONLY: atoms,coord,set_coord_to_origin
        USE random_coord        ,ONLY: max_radius, polar_2_cartesian
        USE gupta               

        REAL(KIND=SGL)  :: random_displacement_ratio 
        REAL(KIND=SGL)  :: energy_compare_ratio 
        CONTAINS

        SUBROUTINE read_bh_param
! read basin hopping parameters
        IMPLICIT NONE
        CHARACTER(LEN=18)        :: bh_param_file = 'param_BH.dat'
        INTEGER                  :: f_bh
        CHARACTER(LEN=2)         :: dummy

        OPEN(NEWUNIT=f_bh, FILE=bh_param_file, STATUS='old')
        READ(f_bh,*) dummy, random_displacement_ratio
        READ(f_bh,*) dummy, energy_compare_ratio
        CLOSE(f_bh)

        END SUBROUTINE read_bh_param

        SUBROUTINE bhop_move
! basin-hopping step
        IMPLICIT NONE
        LOGICAL         :: do_angle_disp ! perform angular displacement
        INTEGER         :: atom_index    ! index of atom to be moved

        CALL random_move(coord)

        do_angle_disp = .FALSE.
        CALL sort_energies(do_angle_disp,atom_index)

        IF (do_angle_disp .EQV. .TRUE.) THEN
          CALL get_new_radius(max_radius)
          CALL angular_displacement(atom_index)
        END IF

        END SUBROUTINE bhop_move


        SUBROUTINE get_new_radius(radius)
! get new radius for each structure before angular move
        IMPLICIT NONE
        REAL(KIND=SGL),INTENT(INOUT) :: radius
        REAL(KIND=DBL),DIMENSION(atoms) :: radius_array

        CALL set_coord_to_origin(coord)
        CALL calc_all_radius(coord,radius_array)
!       CALL calc_distance(coord)

        radius = REAL(MAXVAL(radius_array))

        END SUBROUTINE get_new_radius

        SUBROUTINE calc_all_radius(x_coord,radius_array)
! calculate distance of each atom from origin
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN) :: x_coord
        REAL(KIND=DBL),DIMENSION(:),INTENT(INOUT) :: radius_array
        INTEGER         ::      iter

        DO iter=1,atoms
          radius_array(iter) = NORM2(x_coord(:,iter))
        END DO

        END SUBROUTINE calc_all_radius

        SUBROUTINE angular_displacement(chosen_index)
! angular displacement move for basin-hopping
        IMPLICIT NONE
        INTEGER,INTENT(IN)     :: chosen_index  
        REAL(KIND=DBL),DIMENSION(3)     :: chosen_coord

        chosen_coord = coord(:,chosen_index)
        CALL displace_angle(chosen_coord)
        coord(:,chosen_index) = chosen_coord

        END SUBROUTINE angular_displacement


        SUBROUTINE displace_angle(x_coord)
! random angle move for angular displacement
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(3),INTENT(INOUT)          :: x_coord
        REAL(KIND=SGL)          :: r, theta, phi
        REAL(KIND=SGL),DIMENSION(2)     :: random_array

        r = max_radius
        CALL RANDOM_NUMBER(random_array)
        theta = REAL(PI * random_array(1))
        phi = REAL(2.0 * PI * random_array(2))

        x_coord(1) = r
        x_coord(2) = theta
        x_coord(3) = phi

        CALL polar_2_cartesian(x_coord)
         
        END SUBROUTINE displace_angle

        SUBROUTINE random_move(x_coord)
! random move for basin-hopping
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(INOUT)     :: x_coord
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE       :: temp_coord
        REAL(KIND=SGL),DIMENSION(atoms,3)     :: random_array
        INTEGER         :: iter
        REAL(KIND=SGL)  :: radius

        radius = random_displacement_ratio

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

        x_coord = x_coord + temp_coord

        DEALLOCATE(temp_coord)

        END SUBROUTINE random_move


        SUBROUTINE sort_energies(do_angle_disp,highest_index)
! sort energies of each individual atom
        USE potential           ,ONLY: calc_indv_energy 
        IMPLICIT NONE
        LOGICAL,INTENT(OUT)              :: do_angle_disp    ! YES or NO to carry out angular displacement step
        INTEGER,INTENT(OUT)              :: highest_index
        REAL(KIND=DBL),DIMENSION(atoms)  :: energy_array
        REAL(KIND=DBL)                   :: highest_ene,lowest_ene
!       INTEGER                          :: iter              ! DEBUG use only 

        energy_array = 0.0D0
        CALL calc_indv_energy(coord,atoms,energy_array)
!       CALL indv_gupta_energy(coord,atoms,energy_array)

!DEBUG BEGINS==============================================
!       PRINT *, "energy array in sort_energies"
!       DO iter=1,atoms
!         PRINT *, energy_array(iter)
!       END DO
!DEBUG ENDS==============================================

        highest_ene = MAXVAL(energy_array)
        lowest_ene = MINVAL(energy_array)
        highest_index = MAXLOC(energy_array,DIM=1)

!DEBUG BEGINS==============================================
!       PRINT *, "highest_ene =" , highest_ene
!       PRINT *, "lowest_ene =" , lowest_ene
        PRINT *, "highest_index =" , highest_index
!DEBUG ENDS==============================================

        lowest_ene = ABS(lowest_ene) * DBLE(energy_compare_ratio)
        highest_ene = ABS(highest_ene)
        do_angle_disp = .FALSE.
        IF (highest_ene > lowest_ene) do_angle_disp = .TRUE.

!DEBUG BEGINS==============================================
!       PRINT *, "new lowest_ene =" , lowest_ene
!       PRINT *, "do_angle_disp =", do_angle_disp
!DEBUG ENDS==============================================

        END SUBROUTINE sort_energies 

        END MODULE basin_hopping
