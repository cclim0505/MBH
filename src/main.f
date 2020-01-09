        PROGRAM main
        USE simulation          ,ONLY:simulate_BH
     &    ,simulate_sampling
     &    ,set_up_universal
     &    ,test_eig_rotate
     &    ,test_cut_splice
     &    ,test_improved_random
        USE adf                 ,ONLY:test_adf
!==============================================================
!       MPI OPTIONS
        USE mpi
        USE mpi_var
!==============================================================

        IMPLICIT NONE

!==============================================================
!       MPI OPTIONS
        CALL MPI_INIT(mpi_ierr)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, mpi_ierr)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, mpi_ierr)
!==============================================================

        CALL set_up_universal

!==============================================================
!       MPI OPTIONS
        CALL ch_work_dir
!==============================================================

!       CALL test_eig_rotate

!       CALL test_cut_splice

!       CALL test_improved_random

!       CALL simulate_BH

!       CALL simulate_sampling  ! main subroutine



!==============================================================
!       Testing for ADF
        CALL test_adf
!==============================================================



!==============================================================
!       MPI OPTIONS
        CALL MPI_FINALIZE(mpi_ierr)
!==============================================================

        END PROGRAM main
