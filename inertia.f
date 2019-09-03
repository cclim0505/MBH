        MODULE inertia
        USE constants,          ONLY:DBL 
        USE coord_grad_ene,     ONLY:atoms, coord

        REAL(KIND=DBL),DIMENSION(3,3)   :: inertia_tensor

        CONTAINS

        SUBROUTINE calc_inertia_tensor
        IMPLICIT NONE
        REAL(KIND=DBL)          :: diag_xx, diag_yy, diag_zz
        REAL(KIND=DBL)          :: prod_xy, prod_xz, prod_yz

        CALL calc_diag(1,diag_xx)
        CALL calc_diag(2,diag_yy)
        CALL calc_diag(3,diag_zz)

        CALL calc_product(1,2,prod_xy)
        CALL calc_product(1,3,prod_xz)
        CALL calc_product(2,3,prod_yz)
        
        inertia_tensor(1,1) = diag_xx
        inertia_tensor(2,2) = diag_yy
        inertia_tensor(3,3) = diag_zz

        inertia_tensor(1,2) = prod_xy
        inertia_tensor(2,1) = prod_xy

        inertia_tensor(1,3) = prod_xz
        inertia_tensor(3,1) = prod_xz

        inertia_tensor(2,3) = prod_yz
        inertia_tensor(3,2) = prod_yz

        END SUBROUTINE calc_inertia_tensor

        SUBROUTINE print_inertia_tensor
        IMPLICIT NONE
        INTEGER         :: iter
        DO iter=1,3
          PRINT *, inertia_tensor(iter,1), inertia_tensor(iter,2)
     &      ,inertia_tensor(iter,3)
        END DO
        END SUBROUTINE print_inertia_tensor

        SUBROUTINE calc_diag(axis,diag)
        IMPLICIT NONE
        INTEGER,INTENT(IN)              :: axis
        REAL(KIND=DBL),INTENT(OUT)      :: diag
        INTEGER                         :: aindex,bindex
        REAL(KIND=DBL)              :: temp
        REAL(KIND=DBL)              :: temp_sum
        INTEGER         :: iter

        IF (axis == 1) THEN
          aindex = 2
          bindex = 3
        ELSEIF (axis == 2) THEN
          aindex = 1
          bindex = 3
        ELSEIF (axis == 3) THEN
          aindex = 1
          bindex = 2
        END IF

        temp_sum = 0.0D0
        DO iter=1,atoms
          temp = coord(aindex,iter)**2 + coord(bindex,iter)**2
          temp_sum = temp_sum + temp
        END DO

        diag = temp_sum

        END SUBROUTINE calc_diag

        SUBROUTINE calc_product(aindex,bindex,prod)
        IMPLICIT NONE
        INTEGER,INTENT(IN)          :: aindex, bindex
        REAL(KIND=DBL),INTENT(OUT)  :: prod
        REAL(KIND=DBL)              :: temp
        REAL(KIND=DBL)              :: temp_sum
        INTEGER                     :: iter

        temp_sum = 0.0D0
        DO iter=1,atoms
          temp = coord(aindex,iter) * coord(bindex,iter)
          temp_sum = temp_sum + temp
        END DO
        
        prod = -temp_sum

        END SUBROUTINE calc_product

        END MODULE inertia
