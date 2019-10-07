        MODULE simulation
        CONTAINS

        SUBROUTINE set_up_universal
! set up universal initial values across all MPI processes
        USE initialise          ,ONLY: read_session,read_params
        USE coord_grad_ene      ,ONLY: allocate_coord_gradient 
        IMPLICIT NONE

        CALL read_session
        CALL read_params
        CALL allocate_coord_gradient

        END SUBROUTINE set_up_universal

        SUBROUTINE test_align_then_cut
        IMPLICIT NONE
        END SUBROUTINE test_align_then_cut

        SUBROUTINE test_cut_splice
! testing cut and splice routine within process
        USE constants           ,ONLY:DBL,PI
        USE coord_grad_ene      ,ONLY:coord,atoms
     &    ,read_single_coord
     &    ,printout_xyz
        USE cut_splice          ,ONLY:get_furthest_atom
     &    ,calc_distance_from_furthest
     &    ,cut_cluster,init_upper_lower_cut
     &    ,rotate_upperlower_cut
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE :: lengths
        CHARACTER(LEN=22)       :: input_file='eig_222_input.xyz' 
        INTEGER                 :: f_in
        INTEGER                 :: iter, loop_end=15
        INTEGER                 :: furthest
        INTEGER                 :: cut_point
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE
     &    :: upper_cut, lower_cut
        
        IF (ALLOCATED(lengths))  DEALLOCATE(lengths)
        ALLOCATE(lengths(atoms))



        OPEN(NEWUNIT=f_in,FILE=input_file,STATUS='old')
        DO iter=1,loop_end
          CALL read_single_coord(f_in,coord)
          CALL get_furthest_atom(coord,furthest)
          CALL calc_distance_from_furthest(furthest,coord,lengths)
          CALL init_upper_lower_cut(coord,cut_point)
          IF (ALLOCATED(upper_cut)) DEALLOCATE(upper_cut)
          IF (ALLOCATED(lower_cut)) DEALLOCATE(lower_cut)
          ALLOCATE(upper_cut(3,cut_point))
          ALLOCATE(lower_cut(3,atoms-cut_point))
          CALL cut_cluster(coord,cut_point,upper_cut,lower_cut)
          CALL rotate_upperlower_cut(upper_cut,lower_cut,0.5*PI,.TRUE.)

!DEBUG BEGINS==============================================
          CALL printout_xyz('rotxc_coord.xyz',coord)
          CALL printout_xyz('rotuppercut.xyz',upper_cut)
          CALL printout_xyz('rotlowercut.xyz',lower_cut)
!DEBUG ENDS==============================================
        END DO

        END SUBROUTINE test_cut_splice

        SUBROUTINE test_eig_rotate
! testing principal axis alignment by obtaining rotation matrix created
! from eigenfunction calculated
        USE constants           ,ONLY:DBL
        USE coord_grad_ene      ,ONLY:coord,atoms
     &    ,read_single_coord
     &    ,print_coord, printout_single_coord
     &    ,set_coord_to_origin
        USE inertia             ,ONLY: calc_inertia_tensor
     &    ,print_inertia_tensor 
     &    ,inertia_tensor, eig_val, eig_vec
     &    ,calc_tensor_eig, print_eigs
     &    ,diag_tensor, print_matrix
     &    ,eig_rotate, realign_to_zaxes, realign_eig_vec
     &    ,printout_single_eigs
        USE cut_splice          ,ONLY:get_furthest_atom
     &    ,calc_distance_from_furthest
     &    ,cut_cluster,init_upper_lower_cut
     &    ,rotate_upperlower_cut
        IMPLICIT NONE
        CHARACTER(LEN=22)       :: input_file='eig_222_input.xyz' 
        CHARACTER(LEN=22)       :: output_file='rotate_output.xyz' 
        CHARACTER(LEN=22)       :: eigvec_file='eigvec.dat' 
        CHARACTER(LEN=22)       :: orivec_file='orivec.dat' 
        INTEGER                 :: f_in
        INTEGER                 :: f_out
        INTEGER                 :: f_eig
        INTEGER                 :: f_ori
        INTEGER                 :: iter, loop_end=15

        REAL(KIND=DBL),DIMENSION(3,3)   :: temp_matrix
        REAL(KIND=DBL),DIMENSION(3,3)   :: rotation_matrix
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE   :: temp_coord

!===============================================================
! CUT AND SPLICE VARIABLES
        REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE :: lengths
        INTEGER                 :: furthest
        INTEGER                 :: cut_point
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE
     &    :: upper_cut, lower_cut
!===============================================================

        IF (ALLOCATED(lengths))  DEALLOCATE(lengths)
        ALLOCATE(lengths(atoms))


        IF (ALLOCATED(temp_coord)) DEALLOCATE(temp_coord)
        ALLOCATE(temp_coord(SIZE(coord,1),SIZE(coord,2)))

        OPEN(NEWUNIT=f_in,FILE=input_file,STATUS='old')
        OPEN(NEWUNIT=f_out,FILE=output_file,STATUS='replace'
     &    ,ACCESS='append')
        OPEN(NEWUNIT=f_eig,FILE=eigvec_file,STATUS='replace'
     &    ,ACCESS='append')
        OPEN(NEWUNIT=f_ori,FILE=orivec_file,STATUS='replace'
     &    ,ACCESS='append')

        DO iter=1,loop_end
          CALL read_single_coord(f_in,coord)
          CALL set_coord_to_origin(coord)
          CALL calc_inertia_tensor(coord)
          CALL calc_tensor_eig(inertia_tensor,eig_val,eig_vec)
          CALL printout_single_eigs(f_ori)

          CALL print_eigs


!         CALL eig_rotate(coord,eig_vec,temp_coord)
!         CALL realign_to_zaxes
!         CALL realign_eig_vec

!=============================================================================
          rotation_matrix = TRANSPOSE(eig_vec)
          temp_coord = MATMUL(rotation_matrix,coord)
!         eig_vec = MATMUL(rotation_matrix,eig_vec)

          CALL printout_single_coord(f_out,temp_coord)
          coord = temp_coord
!         CALL printout_single_eigs(f_eig)

!         CALL calc_inertia_tensor(temp_coord)
!         CALL calc_tensor_eig(inertia_tensor,eig_val,eig_vec)
!         CALL printout_single_eigs(f_eig)
!=============================================================================

!         CALL diag_tensor(inertia_tensor,eig_vec,temp_matrix)
!         CALL print_matrix('temp matrix',temp_matrix)
!         CALL calc_tensor_eig(temp_matrix,eig_val,eig_vec)
!=============================================================================
!         CUT AND SPLICE TEST HERE
          CALL get_furthest_atom(coord,furthest)
          CALL calc_distance_from_furthest(furthest,coord,lengths)
          CALL init_upper_lower_cut(coord,cut_point)
          IF (ALLOCATED(upper_cut)) DEALLOCATE(upper_cut)
          IF (ALLOCATED(lower_cut)) DEALLOCATE(lower_cut)
          ALLOCATE(upper_cut(3,cut_point))
          ALLOCATE(lower_cut(3,atoms-cut_point))
          CALL cut_cluster(coord,cut_point,upper_cut,lower_cut)


!=============================================================================
          PRINT *, ''
          PRINT *, ''

        END DO

        CLOSE(f_in)


        END SUBROUTINE test_eig_rotate

        SUBROUTINE test_improved_random
! testing improved version for initializing random coordinates
        USE coord_grad_ene        ,ONLY:coord,printout_xyz
        USE random_coord          ,ONLY:set_random_coord
        IMPLICIT NONE

        CALL set_random_coord
        CALL printout_xyz('improved_random.xyz',coord)

        END SUBROUTINE test_improved_random

        SUBROUTINE simulate_BH
! simulate basin-hopping optimization
        USE constants           ,ONLY: DBL
        USE initialise          ,ONLY: total_mc_step,check_mc_step
        USE coord_grad_ene      ,ONLY: read_coord 
     &    ,coord,optim_coord,old_coord,lowest_coord
     &    ,energy,old_energy,lowest_energy,atoms
     &    , print_coord,printout_xyz,print_lowest_coord
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
!       CALL printout_xyz('3input.xyz',coord)
!DEBUG ENDS==============================================


        CALL set_random_coord


!DEBUG BEGINS==============================================
!       CALL print_coord
!       CALL printout_xyz('1aa.xyz',coord)
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
!       CALL printout_xyz('1bb.xyz',coord)
!DEBUG ENDS==============================================


        DO iter=1,total_mc_step
          ! IF condition to change MC, and basin hopping parameters  
          CALL bhop_move
          CALL local_minim
          IF (optim_ierr == 0) THEN
            coord = optim_coord
            CALL calc_energy(coord,atoms,energy)

!DEBUG BEGINS==============================================
!         CALL printout_xyz('2serial.xyz',coord)
!DEBUG ENDS==============================================

            CALL print_local_coord
            CALL monte_carlo
          ELSE 
            coord = old_coord 
            energy = old_energy
          END IF

          CALL check_mc_step(iter)

        END DO

        CALL print_lowest_coord         ! output lowest coord and energy
        CALL print_lowest_ene           ! output lowest energy

!DEBUG BEGINS==============================================
        CALL calc_inertia_tensor(lowest_coord)
        CALL print_inertia_tensor
!DEBUG ENDS==============================================

!DEBUG BEGINS==============================================
!       PRINT *, 'lowest energy is ', lowest_energy
!       CALL printout_xyz('4MBH_result.xyz',coord)
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
