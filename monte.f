        MODULE monte
        USE constants
        USE coord_grad_ene      ,ONLY: energy,old_energy 
     &    ,lowest_energy
     &    ,coord,old_coord,lowest_coord
     &    ,print_ene_dat,print_lowest_coord
        REAL(KIND=SGL)        :: tempera=0.80000
        REAL(KIND=SGL)        :: acceptance_ratio=0.50000

        CONTAINS

        SUBROUTINE read_mc_param
        IMPLICIT NONE
        CHARACTER(LEN=18)        :: mc_param_file = 'param_MC.dat'
        INTEGER                  :: f_mc
        CHARACTER(LEN=2)         :: dummy

        OPEN(NEWUNIT=f_mc, FILE=mc_param_file, STATUS='old')
        READ(f_mc,*) dummy, tempera
        READ(f_mc,*) dummy, acceptance_ratio
        CLOSE(f_mc)
        
        END SUBROUTINE read_mc_param

        SUBROUTINE monte_carlo
        IMPLICIT NONE
!       REAL(KIND=SGL)  ::      rnum
        LOGICAL         ::      is_accept
        REAL(KIND=DBL)  ::      energy_diff
        REAL(KIND=SGL)  ::      prob_diff

        CALL assign_lowest_energy_coord

        is_accept = .FALSE.
        energy_diff = energy - old_energy

        IF (energy_diff < 0.0 ) THEN
          is_accept = .TRUE.
        ELSE
           CALL calc_prob(energy_diff,prob_diff)
!          CALL RANDOM_NUMBER(rnum)
           IF (prob_diff > acceptance_ratio) THEN
             is_accept = .TRUE.
           END IF
        END IF

        IF (is_accept .EQV. .TRUE.) THEN
          old_energy = energy
          old_coord = coord

!DEBUG BEGINS==============================================
          CALL print_ene_dat('5lowest_ene.dat')
!DEBUG ENDS==============================================

        ELSE
          energy = old_energy
          coord = old_coord
        END IF

        END SUBROUTINE monte_carlo

        SUBROUTINE assign_lowest_energy_coord
        IMPLICIT NONE

        IF (energy < lowest_energy) THEN
          lowest_energy = energy
          lowest_coord = coord
          CALL print_lowest_coord
        END IF

        END SUBROUTINE assign_lowest_energy_coord



        SUBROUTINE calc_prob(ene,prob)
        IMPLICIT NONE
        REAL(KIND=DBL), INTENT(IN)      :: ene
        REAL(KIND=SGL), INTENT(OUT)     :: prob
!       prob = - REAL(ene) / (KB*tempera)
        prob = - REAL(ene) / (tempera)
        prob = EXP(prob)
        END SUBROUTINE calc_prob

        END MODULE monte
