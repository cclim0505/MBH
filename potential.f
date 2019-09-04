        MODULE potential
        USE constants          ,ONLY: DBL
!       USE coord_grad_ene     ,ONLY: gradient, energy
        USE gupta              ,ONLY: gupta_energy, gupta_gradient
!       USE dftb              ,ONLY: dftb_energy, dftb_gradient

        INTEGER         :: potential_type

        CONTAINS

        SUBROUTINE calc_energy(x_coord,natoms,ene_pot)
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN)    :: x_coord
        INTEGER,INTENT(IN)                          :: natoms
        REAL(KIND=DBL),INTENT(OUT)                  :: ene_pot


        IF (potential_type == 1) THEN
          CALL gupta_energy(x_coord,natoms,ene_pot)
!       ELSE IF (potential_type == 2) THEN
!         CALL dftb_energy(x_coord,natoms,ene_pot)
!       ELSE IF (potential_type == 3) THEN
!         CALL density_func_energy(x_coord,natoms,ene_pot)
        END IF

        END SUBROUTINE calc_energy

        SUBROUTINE calc_gradient(x_coord,g_gradient)
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN)     :: x_coord
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(OUT)    :: g_gradient

        IF (potential_type == 1) THEN
          CALL gupta_gradient(x_coord,g_gradient)
!       ELSE IF (potential_type == 2) THEN
!         CALL dftb_gradient(x_coord,g_gradient)
!       ELSE IF (potential_type == 3) THEN
!         CALL density_func_gradient(x_coord,g_gradient)
        END IF

        END SUBROUTINE calc_gradient

        SUBROUTINE calc_indv_energy(x_coord,natoms,ene_pot)
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN)    :: x_coord
        INTEGER,INTENT(IN)                          :: natoms
        REAL(KIND=DBL),INTENT(OUT)                  :: ene_pot


        IF (potential_type == 1) THEN
!         CALL indv_gupta_energy(x_coord,natoms,ene_pot)
!       ELSE IF (potential_type == 2) THEN
!         CALL indv_dftb_energy(x_coord,natoms,ene_pot)
!       ELSE IF (potential_type == 3) THEN
!         CALL indv_density_func_energy(x_coord,natoms,ene_pot)
        END IF

        END SUBROUTINE calc_indv_energy

        END MODULE potential
