        MODULE array_matrix
        USE constants

        CONTAINS

        SUBROUTINE mat_2_arr(mat,arr)
! transform 2D array into 1D array for LBFGS
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN)                :: mat
        REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE,INTENT(OUT)     :: arr
        INTEGER         :: iter,jter
        INTEGER         :: isize,jsize
        INTEGER         :: counter 

        isize = SIZE(mat,1)
        jsize = SIZE(mat,2)

        IF (ALLOCATED(arr)) DEALLOCATE(arr)
        ALLOCATE(arr(SIZE(mat)))

        counter = 1
        DO iter=1,isize
          DO jter=1,jsize
             arr(counter) = mat(iter,jter)
             counter = counter + 1
          END DO
        END DO

        END SUBROUTINE mat_2_arr

        SUBROUTINE arr_2_mat(arr,mat)
! transform 1D array into 2D array after LBFGS
        USE coord_grad_ene,          ONLY:atoms
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(:),INTENT(IN)                  :: arr
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT)   :: mat
        INTEGER         :: iter,jter
        INTEGER         :: arr_size
        INTEGER         :: struct_count
        INTEGER         :: counter

        arr_size = SIZE(arr)
        struct_count = arr_size / (3*atoms)

        IF (ALLOCATED(mat)) DEALLOCATE(mat)
        ALLOCATE(mat(3,atoms))

        counter = 1
        DO iter=1,3
          DO jter=1,atoms
              mat(iter,jter) = arr(counter)
              counter = counter + 1
          END DO
        END DO

        END SUBROUTINE arr_2_mat


        END MODULE array_matrix
