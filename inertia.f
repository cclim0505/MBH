        MODULE inertia
        USE constants          ,ONLY:DBL 
        USE coord_grad_ene     ,ONLY:atoms, coord

        REAL(KIND=DBL),DIMENSION(3,3)   :: inertia_tensor
        REAL(KIND=DBL),DIMENSION(3)     :: eig_val
        REAL(KIND=DBL),DIMENSION(3,3)   :: eig_vec

        CONTAINS

        SUBROUTINE rotate_anticlock(rotate_axis,phi,coord_x
     &    ,is_clockwise)
! rotate input coordinates either along x, y or z-axis, and update input coordinates
        IMPLICIT NONE
        INTEGER,INTENT(IN)                           :: rotate_axis ! 1,2,3 are x,y and z respectively
        REAL(KIND=DBL),INTENT(IN)                    :: phi
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(INOUT)  :: coord_x
        LOGICAL,INTENT(IN)                           :: is_clockwise 
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE    :: temp_x
! rotational matrices
        REAL(KIND=DBL),DIMENSION(3,3)                :: rotate_x
        REAL(KIND=DBL),DIMENSION(3,3)                :: rotate_y
        REAL(KIND=DBL),DIMENSION(3,3)                :: rotate_z

        IF (ALLOCATED(temp_x)) DEALLOCATE(temp_x)
        ALLOCATE(temp_x(SIZE(coord_x,1), SIZE(coord_x,2)))
        temp_x = coord_x

!DEBUG BEGINS==============================================
!       PRINT *, temp_x
!       PRINT *, coord_x
!DEBUG ENDS==============================================

        SELECT CASE(rotate_axis) 
        CASE(1)
          rotate_x(1,1) = 1.0D0
          rotate_x(1,2) = 0.0D0
          rotate_x(1,3) = 0.0D0
          rotate_x(2,1) = 0.0D0
          rotate_x(2,2) = COS(phi)
          rotate_x(2,3) = -SIN(phi)
          rotate_x(3,1) = 0.0D0
          rotate_x(3,2) = SIN(phi)
          rotate_x(3,3) = COS(phi)
          IF (is_clockwise) rotate_x = TRANSPOSE(rotate_x)
          temp_x = MATMUL(rotate_x,temp_x)
        CASE(2)
          rotate_y(1,1) = COS(phi)
          rotate_y(1,2) = 0.0D0
          rotate_y(1,3) = SIN(phi)
          rotate_y(2,1) = 0.0D0
          rotate_y(2,2) = 1.0D0
          rotate_y(2,3) = 0.0D0
          rotate_y(3,1) = -SIN(phi)
          rotate_y(3,2) = 0.0D0
          rotate_y(3,3) = COS(phi)
          IF (is_clockwise) rotate_y = TRANSPOSE(rotate_y)
          temp_x = MATMUL(rotate_y,temp_x)
        CASE(3)
          rotate_z(1,1) = COS(phi)
          rotate_z(1,2) = -SIN(phi)
          rotate_z(1,3) = 0.0D0
          rotate_z(2,1) = SIN(phi)
          rotate_z(2,2) = COS(phi)
          rotate_z(2,3) = 0.0D0
          rotate_z(3,1) = 0.0D0
          rotate_z(3,2) = 0.0D0
          rotate_z(3,3) = 1.0D0
          IF (is_clockwise) rotate_z = TRANSPOSE(rotate_z)
          temp_x = MATMUL(rotate_z,temp_x)
        END SELECT
        coord_x = temp_x

!DEBUG BEGINS==============================================
!       PRINT *, '============================================'
!       PRINT *, 'after rotation'
!       PRINT *, temp_x
!       PRINT *, coord_x
!       PRINT *, '============================================'
!DEBUG ENDS==============================================

        END SUBROUTINE rotate_anticlock

        SUBROUTINE get_phi(vec,phi)
        USE constants           ,ONLY:PI
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(3),INTENT(IN)   :: vec
        REAL(KIND=DBL),INTENT(OUT)               :: phi
        REAL(KIND=DBL)          :: vec_x, vec_y
        INTEGER                 :: x_sign, y_sign
        REAL(KIND=DBL)          :: temp

        vec_x = vec(1)
        vec_y = vec(2)
        
        phi = ATAN( vec_y / vec_x )
        phi = ABS(phi)

        temp = SIGN(1.0 , vec(1))
        IF (temp > 0.0) THEN
          x_sign = 1
        ELSE
          x_sign = -1
        END IF


        temp = SIGN(1.0 , vec(2))
        IF (temp > 0.0) THEN
          y_sign = 1
        ELSE
          y_sign = -1
        END IF

        IF (y_sign > 0 .AND. x_sign > 0) THEN           ! Q1
          phi = phi
        ELSE IF (y_sign > 0 .AND. x_sign < 0) THEN      ! Q2
          phi = PI - phi
        ELSE IF (y_sign < 0 .AND. x_sign < 0) THEN      ! Q3
          phi = PI + phi
        ELSE IF (y_sign < 0 .AND. x_sign > 0) THEN      ! Q4
          phi = 2.0D0 * PI - phi
        END IF

        END SUBROUTINE get_phi

        SUBROUTINE get_psi(vec,psi)
        USE constants           ,ONLY:PI
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(3),INTENT(IN)   :: vec
        REAL(KIND=DBL),INTENT(OUT)               :: psi
        REAL(KIND=DBL),DIMENSION(3)              :: vec_x
        INTEGER                                  :: x_sign, y_sign
        REAL(KIND=DBL)                           :: temp

        vec_x =  (/1.0, 0.0 ,0.0/)
        psi = DOT_PRODUCT(vec,vec_x)
        psi = ACOS(psi)

        temp = SIGN(1.0 , vec(1))
        IF (temp > 0.0) THEN
          x_sign = 1
        ELSE
          x_sign = -1
        END IF

        temp = SIGN(1.0 , vec(2))
        IF (temp > 0.0) THEN
          y_sign = 1
        ELSE
          y_sign = -1
        END IF

        IF (y_sign > 0 .AND. x_sign > 0) THEN           ! Q1
          psi = psi
        ELSE IF (y_sign > 0 .AND. x_sign < 0) THEN      ! Q2
          psi = psi
        ELSE IF (y_sign < 0 .AND. x_sign < 0) THEN      ! Q3
          psi = 2.0D0 * PI - psi
        ELSE IF (y_sign < 0 .AND. x_sign > 0) THEN      ! Q4
          psi = 2.0D0 * PI - psi
        END IF

        END SUBROUTINE get_psi

        SUBROUTINE realign_eig_vec
        USE constants           ,ONLY:PI
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(3)    :: vec
        REAL(KIND=DBL),DIMENSION(3)    :: vec_2, vec_3
        REAL(KIND=DBL)          :: theta, phi, psi

        vec = eig_vec(:,1)

        CALL get_phi(vec,phi)
        theta = ACOS(vec(3))

        CALL rotate_anticlock(3,phi,eig_vec,.TRUE.)
        CALL rotate_anticlock(2,theta,eig_vec,.TRUE.)

        vec =  eig_vec(:,2)
        CALL get_psi(vec,psi)
        CALL rotate_anticlock(3,psi,eig_vec,.TRUE.)


        CALL rotate_anticlock(3,phi,coord,.TRUE.)
        CALL rotate_anticlock(2,theta,coord,.TRUE.)
        CALL rotate_anticlock(3,psi,coord,.TRUE.)

        END SUBROUTINE realign_eig_vec

        SUBROUTINE realign_to_zaxes
! find the rotational axis which has the least moment of inertia and
! align it along the z axis
        USE constants           ,ONLY:PI
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(3)    :: vec
        REAL(KIND=DBL)          :: theta, phi

        vec = eig_vec(:,1)

        phi = ATAN( vec(2) / vec(1) )
        theta = ACOS(vec(3))

        phi = (2.0D0 * PI) - phi
        theta = (2.0D0 * PI)- theta

!DEBUG BEGINS==============================================
        PRINT *, 'vec =', vec
        PRINT *, 'phi and theta are =', phi, theta
!DEBUG ENDS==============================================

!DEBUG BEGINS==============================================
!       PRINT *, 'coord before calling rotate counterclock'
!       PRINT *, coord
!DEBUG ENDS==============================================
        CALL rotate_anticlock(3,phi,coord,.FALSE.)
        CALL rotate_anticlock(2,theta,coord,.FALSE.)
!DEBUG BEGINS==============================================
!       PRINT *, 'coord at the end of realing to zaxes'
!       PRINT *, coord
!DEBUG ENDS==============================================

        END SUBROUTINE realign_to_zaxes

        SUBROUTINE eig_rotate(coord_in, rotation, coord_out)
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN)   :: coord_in
        REAL(KIND=DBL),DIMENSION(3,3),INTENT(IN)   :: rotation
        REAL(KIND=DBL),DIMENSION(:,:)
     &    ,ALLOCATABLE,INTENT(OUT)  :: coord_out

        IF (ALLOCATED(coord_out)) DEALLOCATE(coord_out)
        ALLOCATE(coord_out(SIZE(coord_in,1),SIZE(coord_in,2)))

        PRINT *, coord_in
        PRINT *, rotation

        coord_out = MATMUL(rotation, coord_in)


        END SUBROUTINE eig_rotate

        SUBROUTINE print_matrix(text_out,matrix)
        IMPLICIT NONE
        CHARACTER(LEN=*)        :: text_out
        REAL(KIND=DBL),DIMENSION(3,3),INTENT(IN)   :: matrix
        INTEGER         :: iter
        PRINT *, '================================================'
        PRINT *, '*****',text_out,'*****'
        DO iter=1,3
          PRINT *, matrix(iter,1), matrix(iter,2), matrix(iter,3)
        END DO
        PRINT *, '================================================'
        END SUBROUTINE print_matrix

        SUBROUTINE diag_tensor(tensor_in,vec,tensor_out)
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(3,3),INTENT(IN)   :: tensor_in
        REAL(KIND=DBL),DIMENSION(3,3),INTENT(IN)   :: vec
        REAL(KIND=DBL),DIMENSION(3,3),INTENT(OUT)  :: tensor_out

        REAL(KIND=DBL),DIMENSION(3,3)   :: trans_vec    ! vec transposed
        REAL(KIND=DBL),DIMENSION(3,3)   :: temp

        trans_vec =  TRANSPOSE(vec)
        temp = MATMUL(trans_vec, tensor_in)
        temp = MATMUL(temp, vec)

!       temp(1,2) = 0.0D0
!       temp(1,3) = 0.0D0
!       temp(2,1) = 0.0D0
!       temp(2,3) = 0.0D0
!       temp(3,1) = 0.0D0
!       temp(3,2) = 0.0D0

        tensor_out = temp

        END SUBROUTINE diag_tensor

        SUBROUTINE print_eigs
        IMPLICIT NONE
        INTEGER         :: iter
        PRINT *, 'eigenvales are'
        PRINT *, eig_val(1), eig_val(2), eig_val(3)
        PRINT *, 'eigenvectors are'
        DO iter=1,3
          PRINT *, eig_vec(iter,1), eig_vec(iter,2), eig_vec(iter,3)
        END DO
        END SUBROUTINE print_eigs

        SUBROUTINE printout_single_eigs(f_num)
        IMPLICIT NONE
        INTEGER,INTENT(IN)      :: f_num
        INTEGER,SAVE    :: counter = 0
        INTEGER         :: iter

        counter = counter + 1

        WRITE(f_num,*) counter
        WRITE(f_num,*) 
        DO iter=1,3
          WRITE(f_num,*)  eig_vec(iter,1), eig_vec(iter,2)
     &      , eig_vec(iter,3)
        END DO


        END SUBROUTINE printout_single_eigs




        SUBROUTINE calc_tensor_eig(tensor,val,vec)
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(3,3),INTENT(IN)    :: tensor
        REAL(KIND=DBL),DIMENSION(3),INTENT(OUT)     :: val      ! array of eigenvalues
        REAL(KIND=DBL),DIMENSION(3,3),INTENT(OUT)   :: vec      ! matrix of eigenvectors

        INTEGER, PARAMETER             :: lwmax=1000
        INTEGER                        :: lwork
        REAl(KIND=4),DIMENSION(3,3)    :: matrix                ! temporary matrix for ssyev call
        REAl(KIND=4),DIMENSION(3)      :: w
        REAl(KIND=4),DIMENSION(lwmax)  :: work
        INTEGER                        :: ssyev_info

        matrix = tensor

        lwork = -1
        CALL ssyev('V','U', 3, matrix, 
     &        3, w, work, lwork, ssyev_info ) 

        lwork = MIN(lwmax, INT(work(1)))

!       PRINT *, 'lwork = ', lwork

        CALL ssyev('V','U', 3, matrix, 
     &        3, w, work, lwork, ssyev_info ) 

        val = w
        vec = matrix

        END SUBROUTINE calc_tensor_eig

        SUBROUTINE calc_inertia_tensor(x_coord)
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN) :: x_coord
        REAL(KIND=DBL)          :: diag_xx, diag_yy, diag_zz
        REAL(KIND=DBL)          :: prod_xy, prod_xz, prod_yz

        CALL calc_diag(1,x_coord,diag_xx)
        CALL calc_diag(2,x_coord,diag_yy)
        CALL calc_diag(3,x_coord,diag_zz)

        CALL calc_product(1,2,x_coord,prod_xy)
        CALL calc_product(1,3,x_coord,prod_xz)
        CALL calc_product(2,3,x_coord,prod_yz)
        
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
        PRINT *, '=========INERTIA TENSOR========================'
        DO iter=1,3
          PRINT *, inertia_tensor(iter,1), inertia_tensor(iter,2)
     &      ,inertia_tensor(iter,3)
        END DO
        PRINT *, '==============================================='
        END SUBROUTINE print_inertia_tensor

        SUBROUTINE calc_diag(axis,x_coord,diag)
        IMPLICIT NONE
        INTEGER,INTENT(IN)              :: axis
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN) :: x_coord
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
          temp = x_coord(aindex,iter)**2 + x_coord(bindex,iter)**2
          temp_sum = temp_sum + temp
        END DO

        diag = temp_sum

        END SUBROUTINE calc_diag

        SUBROUTINE calc_product(aindex,bindex,x_coord,prod)
        IMPLICIT NONE
        INTEGER,INTENT(IN)          :: aindex, bindex
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN) :: x_coord
        REAL(KIND=DBL),INTENT(OUT)  :: prod
        REAL(KIND=DBL)              :: temp
        REAL(KIND=DBL)              :: temp_sum
        INTEGER                     :: iter

        temp_sum = 0.0D0
        DO iter=1,atoms
          temp = x_coord(aindex,iter) * x_coord(bindex,iter)
          temp_sum = temp_sum + temp
        END DO
        
        prod = -temp_sum

        END SUBROUTINE calc_product

        END MODULE inertia
