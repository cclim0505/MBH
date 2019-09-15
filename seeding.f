        PROGRAM main
        USE mpi
        IMPLICIT NONE

        INTEGER         :: myid, numprocs, mpi_ierr
        INTEGER,DIMENSION(1)         :: seed
        REAL                         :: real_seed
        REAL,DIMENSION(5)    :: array

        CALL MPI_INIT(mpi_ierr)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, mpi_ierr)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, mpi_ierr)

!       PRINT *, 'I am process', myid
        CALL CPU_TIME(real_seed)
!       PRINT *, 'from id', myid, 'real_seed', real_seed
        seed = INT(1E7*real_seed)
!       PRINT *, 'seed is', seed


!       CALL RANDOM_SEED(GET=seed)
        CALL RANDOM_SEED(PUT=seed+myid)
        CALL RANDOM_NUMBER(array)
        PRINT *, array, 'from process', myid

        CALL MPI_FINALIZE(mpi_ierr)

        END PROGRAM main
