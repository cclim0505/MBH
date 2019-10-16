        MODULE mpi_var
        INTEGER                 :: myid, numprocs, mpi_ierr
        INTEGER,PARAMETER       :: root = 0

        CONTAINS
        SUBROUTINE ch_work_dir
! change working directory for each worker process
        IMPLICIT NONE
        CHARACTER(20)   :: dir
        CHARACTER(6)    :: prefix="worker"

        WRITE(dir,*) myid

        IF (myid < 10) THEN
          dir = prefix//"00"//ADJUSTL(dir)
        ELSE IF (myid < 100) THEN
          dir = prefix//"0"//ADJUSTL(dir)
        ELSE
          dir = prefix//ADJUSTL(dir)
        END IF

        dir = "./"//dir
        CALL chdir(dir)

        END SUBROUTINE ch_work_dir

        END MODULE mpi_var
