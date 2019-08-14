        PROGRAM main
        USE types
        USE initialise
        USE gupta
        USE param
        USE random_coord
        USE basin_hopping
        USE optimization
        IMPLICIT NONE
        REAL(KIND=DBL)  :: energy 
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE  :: optim_coord
        REAL(KIND=DBL),DIMENSION(3)     :: centre
        REAL(KIND=DBL),DIMENSION(10) :: radius_array
!       TYPE(system)    :: cluster
!       cluster = system('AuAu', 4, 3.0)
!       PRINT *, cluster

        CALL init_gupta
        CALL read_atoms
        CALL read_coord
        IF (ALLOCATED(optim_coord)) DEALLOCATE(optim_coord)
        ALLOCATE (optim_coord(3,atoms))
!       CALL print_coord
        

!       CALL gupta_energy(coord,atoms,energy)


        CALL set_random_coord
        CALL calc_all_radius(coord,radius_array)

        CALL print_coord

        CALL set_coord_to_origin
        CALL calc_all_radius(coord,radius_array)

!       CALL calc_centroid(coord,centre)
!       CALL print_coord
!       CALL gupta_energy(coord,atoms,energy)
!       CALL gupta_gradient(coord,gradient)

!       CALL set_coord_to_origin
!       CALL print_coord

        CALL bhop_move
!       CALL angular_displacement

!       CALL optim_lbfgs(coord,optim_coord)
!       CALL print_coord




!       CALL open_gupta_parameter
!       CALL convert_coord_to_x(coord,posit)
!       CALL agrad(posit,gradient)

        END PROGRAM main
