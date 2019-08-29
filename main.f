        PROGRAM main
        USE constants           ,ONLY: DBL
        USE initialise          ,ONLY: read_session,read_params
     &    ,total_mc_step
        USE coord_grad_ene      ,ONLY: allocate_coord_gradient 
     &    ,read_atoms,read_coord 
     &    ,coord,optim_coord,old_coord,lowest_coord
     &    ,energy,old_energy,lowest_energy,atoms
     &    , print_coord,print_coord_xyz
        USE gupta               ,ONLY: gupta_energy
        USE random_coord        ,ONLY: set_random_coord
        USE optimization        ,ONLY: local_minim
        USE basin_hopping       ,ONLY: bhop_move
        USE monte               ,ONLY: monte_carlo

        IMPLICIT NONE
        INTEGER                 :: iter
        REAL(KIND=DBL),DIMENSION(3)     :: centre
        REAL(KIND=DBL),DIMENSION(10) :: radius_array


!=========================================
        CALL read_session
        CALL read_params
        CALL allocate_coord_gradient

!DEBUG BEGINS==============================================
        CALL read_coord
        CALL gupta_energy(coord,atoms,energy)
        PRINT *, 'ground_state is ', energy
        CALL print_coord_xyz('3input.xyz')
!DEBUG ENDS==============================================

        CALL set_random_coord


!DEBUG BEGINS==============================================
        CALL print_coord
        CALL print_coord_xyz('1aa.xyz')
!DEBUG ENDS==============================================

        CALL local_minim
        coord = optim_coord
        old_coord = coord
        CALL gupta_energy(coord,atoms,energy)
        old_energy = energy

        lowest_coord = old_coord
        lowest_energy = old_energy

!DEBUG BEGINS==============================================
        CALL print_coord
        CALL print_coord_xyz('1bb.xyz')
!DEBUG ENDS==============================================


        DO iter=1,total_mc_step
          ! IF condition to change MC, and basin hopping parameters  
          CALL bhop_move
          CALL local_minim
          coord = optim_coord
          CALL gupta_energy(coord,atoms,energy)

!DEBUG BEGINS==============================================
        CALL print_coord_xyz('2serial.xyz')
!DEBUG ENDS==============================================

        CALL monte_carlo

        END DO

!DEBUG BEGINS==============================================
        PRINT *, 'lowest energy is ', old_energy
        CALL print_coord_xyz('4MBH_result.xyz')
!DEBUG ENDS==============================================

! Output global coordinates and energy


!=========================================
!TESTING ONLY
!       CALL read_atoms
!       CALL read_coord
!       CALL gupta_energy(coord,atoms,energy)
!       CALL print_coord
!=========================================

!       CALL print_coord
        

!       CALL sort_energies




!       CALL set_random_coord
!       CALL calc_all_radius(coord,radius_array)

!       CALL print_coord

!       CALL set_coord_to_origin
!       CALL calc_all_radius(coord,radius_array)

!       CALL calc_centroid(coord,centre)
!       CALL print_coord
!       CALL gupta_energy(coord,atoms,energy)
!       CALL gupta_gradient(coord,gradient)

!       CALL set_coord_to_origin
!       CALL print_coord
!       CALL bhop_move
!       CALL angular_displacement

!       CALL optim_lbfgs(coord,optim_coord)
!       CALL print_coord

!       CALL open_gupta_parameter

        END PROGRAM main
