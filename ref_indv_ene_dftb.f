        PROGRAM main
        IMPLICIT NONE
        INTEGER         :: iter
        INTEGER,PARAMETER         :: natoms=12
        INTEGER         :: f_num
        CHARACTER(LEN=12)      :: file_in='detailed.out'
        CHARACTER(LEN=8)       :: dum1, dum2, dum3, dum4
        CHARACTER(LEN=1)       :: dummy
        CHARACTER(LEN=4),PARAMETER      :: match1='Atom'
        CHARACTER(LEN=8),PARAMETER      :: match2='resolved'
        CHARACTER(LEN=5),PARAMETER      :: match3='total'
        CHARACTER(LEN=8),PARAMETER      :: match4='energies'
        INTEGER         :: ierr

        REAL(KIND= 8),DIMENSION(natoms)  :: ene_array

        OPEN(NEWUNIT=f_num,FILE=file_in,STATUS='old')

        DO
          ierr = 0
          READ(f_num,*,IOSTAT=ierr) dum1, dum2, dum3, dum4

          IF(ierr /= 0) CYCLE
          IF(dum1 == match1 .AND. dum2 == match2
     &      .AND. dum3 == match3 .AND. dum4 == match4) EXIT

        END DO

        DO iter=1,natoms
          READ(f_num,*) dummy, ene_array(iter)
        END DO

        CLOSE(f_num)

        DO iter=1,natoms
          PRINT *, ene_array(iter)
        END DO

        END PROGRAM main
