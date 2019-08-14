      MODULE param
      IMPLICIT NONE
      REAL:: epsilon_a,zeta_a,p_a,q_a,rzero_a      ! element A-A interaction
      REAL:: epsilon_b,zeta_b,p_b,q_b,rzero_b      ! element B-B interaction
      REAL:: epsilon_ab,zeta_ab,p_ab,q_ab,rzero_ab ! element A-B interaction

      SAVE

      CONTAINS
        SUBROUTINE open_gupta_parameter
! Read parameters for Gupta type many body potential based on type of material.
         OPEN (16,file="AuAu_parameter1.dat")
         READ (16,*) epsilon_a,zeta_a,p_a,q_a,rzero_a
         READ (16,*) epsilon_b,zeta_b,p_b,q_b,rzero_b
         READ (16,*) epsilon_ab,zeta_ab,p_ab,q_ab,rzero_ab
         CLOSE(16)
        END SUBROUTINE open_gupta_parameter
  
      END MODULE param
 
