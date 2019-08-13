        PROGRAM main
        USE types
        USE initialise
        USE gupta
        USE param
        USE random_coord
        USE basin_hopping
        IMPLICIT NONE
        REAL(KIND=DBL)  :: energy 
!       TYPE(system)    :: cluster
!       cluster = system('AuAu', 4, 3.0)
!       PRINT *, cluster

        CALL init_gupta
        CALL read_atoms
        CALL read_coord
!       CALL print_coord

!       CALL gupta_energy(coord,atoms,energy)


        CALL set_random_coord
        CALL print_coord
        CALL gupta_energy(coord,atoms,energy)

        CALL set_coord_to_origin
        CALL print_coord

!       CALL bhop_move
!       CALL angular_displacement


!       CALL gupta_gradient(coord,gradient)

!       CALL open_gupta_parameter
!       CALL convert_coord_to_x(coord,posit)
!       CALL agrad(posit,gradient)

        END PROGRAM main
