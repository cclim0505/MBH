        MODULE monte
        USE constants           ,ONLY:SGL,DBL
        USE coord_grad_ene      ,ONLY: energy,old_energy 
     &    ,lowest_energy, lowest_step
     &    ,coord,old_coord,lowest_coord
     &    ,print_ene_dat,print_update_lowest_coord
     &    ,resume_print_lowest_coord
     &    ,resume_print_old_coord

        REAL(KIND=SGL)        :: tempera=0.80000
        REAL(KIND=SGL)        :: acceptance_ratio=0.50000
        LOGICAL               :: is_cut_splice
        LOGICAL               :: is_cage_drive
        LOGICAL               :: is_ring_drive
        INTEGER                 :: pre_cut_splice_period
        INTEGER                 :: cut_splice_freq
        INTEGER                 :: cage_drive_freq

        PRIVATE  :: assign_lowest_energy_coord
        PRIVATE  :: calc_prob

        PUBLIC  :: read_mc_param
        PUBLIC  :: printout_mc_param
        PUBLIC  :: monte_carlo

        CONTAINS


        SUBROUTINE read_mc_param
! read parameters for Monte Carlo steps
        USE directory           ,ONLY: session_dir
        IMPLICIT NONE
        CHARACTER(LEN=18)        :: mc_param_file = 'param_MC.dat'
        INTEGER                  :: f_mc
        CHARACTER(LEN=2)         :: dummy

        OPEN(NEWUNIT=f_mc
     &    , FILE='./'//session_dir//'/'//mc_param_file, STATUS='old')
        READ(f_mc,*) dummy, tempera
        READ(f_mc,*) dummy, acceptance_ratio
        READ(f_mc,*) dummy, is_cut_splice
        READ(f_mc,*) dummy, pre_cut_splice_period
        READ(f_mc,*) dummy, cut_splice_freq
        READ(f_mc,*) dummy, is_cage_drive
        READ(f_mc,*) dummy, cage_drive_freq
        READ(f_mc,*) dummy, is_ring_drive
        CLOSE(f_mc)
        
        END SUBROUTINE read_mc_param

        SUBROUTINE printout_mc_param
! printout parameters for Monte Carlo steps
        USE directory           ,ONLY: saved_session
        IMPLICIT NONE
        CHARACTER(LEN=18)        :: mc_param_file = 'saved_param_MC.dat'
        INTEGER                  :: f_mc

        OPEN(NEWUNIT=f_mc
     &    , FILE='./'//saved_session//'/'//mc_param_file
     &    , STATUS='replace')
        WRITE(f_mc,*) 'temperature_ratio', tempera
        WRITE(f_mc,*) 'acceptance_ratio', acceptance_ratio
        WRITE(f_mc,*) 'is_cut_splice', is_cut_splice
        WRITE(f_mc,*) 'pre_cut_splice_period', pre_cut_splice_period
        WRITE(f_mc,*) 'cut_splice_freq', cut_splice_freq
        WRITE(f_mc,*) 'is_cage_drive', is_cage_drive
        WRITE(f_mc,*) 'cage_drive_freq', cage_drive_freq
        WRITE(f_mc,*) 'is_ring_drive', is_ring_drive
        CLOSE(f_mc)
        
        END SUBROUTINE printout_mc_param

        SUBROUTINE monte_carlo(step)
! Monte Carlo simulation to accept or reject new configuration
        IMPLICIT NONE
        INTEGER,INTENT(IN)      :: step
!       REAL(KIND=SGL)  ::      rnum
        LOGICAL         ::      is_accept
        REAL(KIND=DBL)  ::      energy_diff
        REAL(KIND=SGL)  ::      prob_diff

        CALL assign_lowest_energy_coord(step)

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
          CALL resume_print_old_coord

!DEBUG BEGINS==============================================
!         CALL print_ene_dat('03_lowest_energies.dat')
!DEBUG ENDS==============================================

        ELSE
          energy = old_energy
          coord = old_coord
        END IF

        END SUBROUTINE monte_carlo

        SUBROUTINE assign_lowest_energy_coord(step)
! update coordinates with lower energy when condition is met
        IMPLICIT NONE
        INTEGER,INTENT(IN)      :: step

        IF (energy < lowest_energy) THEN
          lowest_energy = energy
          lowest_coord = coord
          lowest_step  = step
          CALL print_update_lowest_coord(step)
          CALL resume_print_lowest_coord(step)
        END IF

        END SUBROUTINE assign_lowest_energy_coord



        SUBROUTINE calc_prob(ene,prob)
! calculate probabilty of exponential function
        IMPLICIT NONE
        REAL(KIND=DBL), INTENT(IN)      :: ene
        REAL(KIND=SGL), INTENT(OUT)     :: prob
!       prob = - REAL(ene) / (KB*tempera)
        prob = - REAL(ene) / (tempera)
        prob = EXP(prob)
        END SUBROUTINE calc_prob

        END MODULE monte
