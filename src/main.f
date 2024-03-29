        PROGRAM main
        USE banner
        USE simulation          ,ONLY:simulate_sampling
     &    ,set_up_universal
     &    ,test_eig_rotate
     &    ,test_cut_splice
        USE geometric_drive
        USE coord_grad_ene
        USE random_coord        ,ONLY:test_2d_init

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
!       Banner output
        IF (myid == 0) CALL print_version
!==============================================================

        CALL set_up_universal

!==============================================================
!       MPI OPTIONS
        CALL ch_work_dir
!==============================================================

        CALL simulate_sampling     ! actual MBH sampling.

!==============================================================
! Testing section
!        CALL test_eig_rotate

!        CALL test_cut_splice


!        CALL test_cage3

!         CALL test_ring
         
!         CALL test_2d_init

!        CALL test_coord1

!==============================================================
!       MPI OPTIONS
        CALL MPI_FINALIZE(mpi_ierr)
!==============================================================

        END PROGRAM main
