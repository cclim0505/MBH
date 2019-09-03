        MODULE potential
        USE coord_grad_ene     ,ONLY: coord, gradient, atoms, energy
        USE gupta              ,ONLY: gupta_energy, gupta_gradient
!       USE dftb              ,ONLY: dftb_energy, dftb_gradient

        INTEGER         :: potential_type

        CONTAINS

        SUBROUTINE calc_energy
        IMPLICIT NONE

        IF (potential_type == 1) THEN
          CALL gupta_energy(coord,atoms,energy)
!       ELSE IF (potential_type == 2) THEN
!         CALL dftb_energy(coord,atoms,energy)
!       ELSE IF (potential_type == 3) THEN
!         CALL density_func_energy(coord,atoms,energy)
        END IF

        END SUBROUTINE calc_energy

        SUBROUTINE calc_grad
        IMPLICIT NONE

        IF (potential_type == 1) THEN
          CALL gupta_gradient(coord,gradient)
!       ELSE IF (potential_type == 2) THEN
!         CALL dftb_gradient(coord,gradient)
!       ELSE IF (potential_type == 3) THEN
!         CALL density_func_gradient(coord,gradient)
        END IF

        END SUBROUTINE calc_grad

        END MODULE potential
