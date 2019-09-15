        MODULE simulation
        CONTAINS

        SUBROUTINE set_up_universal
        USE initialise          ,ONLY: read_session,read_params
        USE coord_grad_ene      ,ONLY: allocate_coord_gradient 
        IMPLICIT NONE

        CALL read_session
        CALL read_params
        CALL allocate_coord_gradient

        END SUBROUTINE set_up_universal

        SUBROUTINE simulate_BH
        USE constants           ,ONLY: DBL
        USE initialise          ,ONLY: total_mc_step
        USE coord_grad_ene      ,ONLY: read_atoms,read_coord 
     &    ,coord,optim_coord,old_coord,lowest_coord
     &    ,energy,old_energy,lowest_energy,atoms
     &    , print_coord,print_coord_xyz,print_lowest_coord
     &    , print_lowest_ene 
     &    , print_local_coord
     &    , gradient
        USE random_coord        ,ONLY: set_random_coord
        USE optimization        ,ONLY: local_minim
     &    , optim_ierr
        USE basin_hopping       ,ONLY: bhop_move
        USE monte               ,ONLY: monte_carlo
        USE inertia             ,ONLY: calc_inertia_tensor
     &    ,print_inertia_tensor 
        USE potential           ,ONLY: calc_energy


!DEBUG BEGINS==============================================
!       USE dftb                ,ONLY: coord_2_gen,dftb_energy
!    &    ,dftb_gradient, dftb_both_ene_grad
!DEBUG ENDS==============================================

        IMPLICIT NONE
        INTEGER                 :: iter

!DEBUG BEGINS==============================================
!       CALL read_coord
!       CALL calc_energy(coord,atoms,energy)
!       PRINT *, 'ground_state is ', energy
!       CALL print_coord_xyz('3input.xyz')
!DEBUG ENDS==============================================


        CALL set_random_coord


!DEBUG BEGINS==============================================
!       CALL print_coord
!       CALL print_coord_xyz('1aa.xyz')
!DEBUG ENDS==============================================

        CALL local_minim
        coord = optim_coord
        old_coord = coord
        CALL calc_energy(coord,atoms,energy)
        old_energy = energy

        lowest_coord = old_coord
        lowest_energy = old_energy

!DEBUG BEGINS==============================================
!       CALL print_coord
!       CALL print_coord_xyz('1bb.xyz')
!DEBUG ENDS==============================================


        DO iter=1,total_mc_step
          ! IF condition to change MC, and basin hopping parameters  
          CALL bhop_move
          CALL local_minim
          IF (optim_ierr == 0) THEN
            coord = optim_coord
            CALL calc_energy(coord,atoms,energy)

!DEBUG BEGINS==============================================
!         CALL print_coord_xyz('2serial.xyz')
!DEBUG ENDS==============================================

            CALL print_local_coord
            CALL monte_carlo
          ELSE 
            coord = old_coord 
            energy = old_energy
          END IF

        END DO

        CALL print_lowest_coord         ! output lowest coord and energy
        CALL print_lowest_ene           ! output lowest energy

!DEBUG BEGINS==============================================
        CALL calc_inertia_tensor
        CALL print_inertia_tensor
!DEBUG ENDS==============================================

!DEBUG BEGINS==============================================
!       PRINT *, 'lowest energy is ', lowest_energy
!       CALL print_coord_xyz('4MBH_result.xyz')
!DEBUG ENDS==============================================


!DEBUG BEGINS==============================================
!       CALL coord_2_gen(coord,atoms)
!       CALL dftb_both_ene_grad(coord,atoms,energy,gradient)
!       CALL dftb_energy(coord,atoms,energy) 
!       PRINT *, 'dftb energy is', energy
!       CALL dftb_gradient(coord,atoms,gradient)
!       PRINT *, 'dftb gradient is' 
!       PRINT *, gradient
!DEBUG ENDS==============================================

        END SUBROUTINE simulate_BH

        END MODULE simulation
