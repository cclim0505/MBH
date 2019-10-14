        MODULE gupta
        USE constants           ,ONLY: DBL
        REAL(KIND=DBL)        :: a_ij
        REAL(KIND=DBL)        :: eta
        REAL(KIND=DBL)        :: p_ij
        REAL(KIND=DBL)        :: q_ij
        REAL(KIND=DBL)        :: r_zero
        REAl(KIND=DBl),DIMENSION(:,:),ALLOCATABLE  :: distance !  distances matrix, for gupta band energy calculation

        PRIVATE  :: calc_distance
        PRIVATE  :: gupta_repulsive
        PRIVATE  :: gupta_band
        PRIVATE  :: set_rij_size
        PRIVATE  :: gupta_repulse_gradient
        PRIVATE  :: gupta_band_gradient
        PRIVATE  :: indv_gupta_repulsive
        PRIVATE  :: indv_gupta_band

        PUBLIC  :: read_gupta_param
        PUBLIC  :: gupta_both_ene_grad
        PUBLIC  :: gupta_energy
        PUBLIC  :: gupta_gradient
        PUBLIC  :: indv_gupta_energy


        CONTAINS

        SUBROUTINE read_gupta_param
! read gupta potential parameters
        IMPLICIT NONE
        INTEGER         :: f_gupta

        OPEN (NEWUNIT=f_gupta,FILE='AuAu_parameter1.dat',STATUS='old')
        READ(f_gupta,*) a_ij, eta, p_ij, q_ij, r_zero
        CLOSE(f_gupta)

        END SUBROUTINE read_gupta_param

        SUBROUTINE gupta_both_ene_grad(x_coord,natoms,ene,grad)
! calculate gupta energy and gradient
        IMPLICIT NONE
        REAl(KIND=DBL),DIMENSION(:,:),INTENT(IN)  :: x_coord
        INTEGER,INTENT(IN)                        :: natoms
        REAL(KIND=DBL),INTENT(OUT)                :: ene
        REAl(KIND=DBL),DIMENSION(:,:),INTENT(OUT) :: grad

        CALL gupta_energy(x_coord,natoms,ene)
        CALL gupta_gradient(x_coord,grad)

        END SUBROUTINE gupta_both_ene_grad

        SUBROUTINE calc_distance(x_coord)
! calculate the distances between all atomic pairs
        IMPLICIT NONE
        REAl(KIND=DBL),DIMENSION(:,:),INTENT(IN) :: x_coord
        REAL(KIND=DBL),DIMENSION(3)              :: temp
        INTEGER :: natoms
        INTEGER :: iter,jter 
        natoms = SIZE(x_coord,2)
! distances matrix, for gupta band energy calculation
        IF (ALLOCATED(distance)) DEALLOCATE(distance)
        ALLOCATE(distance(natoms,natoms))
        DO iter=1,natoms
          DO jter=1,natoms
            temp = x_coord(:,iter) -  x_coord(:,jter)
            distance(iter,jter) = NORM2(temp) 
          END DO
        END DO
        END SUBROUTINE calc_distance

        SUBROUTINE gupta_energy(x_coord,natoms,energy)
! compute gupta energy
        IMPLICIT NONE
        REAl(KIND=DBL),DIMENSION(:,:),INTENT(IN) :: x_coord
        INTEGER,INTENT(IN)                       :: natoms
        REAL(KIND=DBL),INTENT(OUT)               :: energy
        REAL(KIND=DBL)                           :: band, repulse
        INTEGER                                  :: iter,jter,counter
        REAl(KIND=DBl),DIMENSION(:),ALLOCATABLE  :: rij
        REAL(KIND=DBL),DIMENSION(3)              :: temp
        INTEGER                                  :: rij_size

!=====================================================================
! distances matrix, for gupta band energy calculation
!       IF (ALLOCATED(distance)) DEALLOCATE(distance)
!       ALLOCATE(distance(natoms,natoms))
!       DO iter=1,natoms
!         DO jter=1,natoms
!           temp = x_coord(:,iter) -  x_coord(:,jter)
!           distance(iter,jter) = NORM2(temp) 
!         END DO
!       END DO
        CALL calc_distance(x_coord)

!=====================================================================
! Single counting for distances(rij), for gupta repulsive energy calculation
!        CALL combination(natoms,2,rij_size)
        CALL set_rij_size(natoms,rij_size)
        IF (ALLOCATED(rij)) DEALLOCATE(rij)
        ALLOCATE(rij(rij_size))
        counter = 1
        DO iter=1 ,natoms-1
          DO jter=iter+1, natoms
            temp = x_coord(:,iter) - x_coord(:,jter) 
            rij(counter) = NORM2(temp)
            counter = counter + 1
          END DO
        END DO

!=====================================================================
! rij checking
!       PRINT *, '' 
!       PRINT *, 'rij' 
!       DO iter=1,rij_size
!         PRINT *, rij(iter)
!       END DO
!       PRINT *, '' 

!=====================================================================
! Calculate new gupta band
        CALL gupta_repulsive(rij,repulse)
        CALL gupta_band(distance,band)
        energy = band + repulse
        energy = energy / REAL(natoms)

!DEBUG==============================================
!       PRINT *, "band energy in eV is", band
!       PRINT *, "repulsive energy in eV is", repulse
!       PRINT *, "energy per atom in eV is", energy
!DEBUG==============================================

        END SUBROUTINE gupta_energy

        SUBROUTINE gupta_repulsive(rij,energy)
! compute gupta repulsive energy
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(:),INTENT(IN)  :: rij
        REAL(KIND=DBL),INTENT(OUT)              :: energy
        REAL(KIND=DBL)                          :: temp, temp_sum
        INTEGER                                 :: iter
        INTEGER                                 :: rij_size
        rij_size = SIZE(rij)
        temp_sum = 0.0D0
        DO iter=1,rij_size
          temp = (rij(iter) / r_zero) - 1.0D0
          temp = temp * (-p_ij)
          temp = a_ij*DEXP(temp)
          temp_sum = temp_sum + temp
        END DO
        energy =  2.0D0 * temp_sum 
        END SUBROUTINE gupta_repulsive

        SUBROUTINE gupta_band(distance,energy)
! compute gupta band energy
        REAL(KIND=DBL),DIMENSION(:,:),INTENT(IN)        :: distance
        REAL(KIND=DBL),INTENT(OUT)                      :: energy
        REAL(KIND=DBL)          :: temp, temp_sum, temp_total
        INTEGER :: iter,jter
        INTEGER :: distance_dim
        distance_dim = SIZE(distance,1)
        temp_total = 0.0D0
        DO iter=1,distance_dim
          temp_sum = 0.0D0
          DO jter=1,distance_dim
             IF(iter==jter) CYCLE 
             temp = (distance(iter,jter) / r_zero) - 1.0D0
             temp = temp * 2.0D0*(-q_ij)
             temp = (eta**2)*DEXP(temp)
             temp_sum = temp_sum + temp
          END DO
          temp_total = temp_total - DSQRT(temp_sum)
        END DO
        energy = temp_total
        END SUBROUTINE gupta_band

        SUBROUTINE set_rij_size(natoms,rij_size)
! determines the array size for array rij
        IMPLICIT NONE
        INTEGER,INTENT(IN)      :: natoms
        INTEGER,INTENT(OUT)     :: rij_size
        INTEGER                 :: temp
        temp = (natoms**2) - natoms
        rij_size = temp / 2
        END SUBROUTINE set_rij_size

        SUBROUTINE gupta_gradient(x_coord,grad)
! calculate gradient of gupta potential
        IMPLICIT NONE
        REAl(KIND=DBL),DIMENSION(:,:),INTENT(IN)  :: x_coord
        REAl(KIND=DBL),DIMENSION(:,:),INTENT(OUT) :: grad
        REAl(KIND=DBL),DIMENSION(:,:),ALLOCATABLE :: repulse
        REAl(KIND=DBL),DIMENSION(:,:),ALLOCATABLE :: band
        INTEGER :: size_1, size_2
!       INTEGER :: iter

        size_1 = SIZE(grad,1)
        size_2 = SIZE(grad,2)

        IF (ALLOCATED(repulse)) DEALLOCATE(repulse)
        IF (ALLOCATED(band)) DEALLOCATE(band)
        ALLOCATE(repulse(size_1,size_2))
        ALLOCATE(band(size_1,size_2))

        CALL gupta_repulse_gradient(x_coord,repulse)
        CALL gupta_band_gradient(x_coord,band)

        grad = repulse + band

!DEBUG==============================================
!       PRINT *, 'repulse and band value'
!       PRINT *, repulse(1,1), band(1,1)

!       PRINT *, ''
!       PRINT *, 'grad values'
!       DO iter=1,size_2
!         PRINT *, grad(1,iter), grad(2,iter), grad(3,iter)
!       END DO
!       PRINT *, ''
!DEBUG==============================================

        DEALLOCATE(repulse, band)

        END SUBROUTINE gupta_gradient


        SUBROUTINE gupta_repulse_gradient(x_coord,grad)
! calculate gradient of gupta repulsive term
        IMPLICIT NONE
        REAl(KIND=DBL),DIMENSION(:,:),INTENT(IN)  :: x_coord
        REAl(KIND=DBL),DIMENSION(:,:),INTENT(OUT) :: grad
        REAL(KIND=DBL)  :: temp_x, temp_y, temp_z
        REAL(KIND=DBL)  :: sum_temp_x, sum_temp_y, sum_temp_z
        REAL(KIND=DBL)  :: expon
        REAL(KIND=DBL)  :: temp
        REAL(KIND=DBL)  :: rij
        INTEGER         :: natoms
        INTEGER         :: iter, jter
        natoms = SIZE(x_coord,2)
        grad = 0.0D0
        DO iter=1,natoms
          sum_temp_x = 0.0D0
          sum_temp_y = 0.0D0
          sum_temp_z = 0.0D0
          DO jter=1,natoms
            IF(iter==jter) CYCLE 
            rij = distance(iter,jter)
            temp = (rij / r_zero) - 1.0D0
            temp = temp * (-p_ij)
            temp = a_ij*DEXP(temp)
            expon = temp * (-p_ij) / r_zero
!=======================
            temp = expon / rij

!DEBUG==============================================
!           PRINT *, 'GRAD_ATT', iter, jter, temp
!DEBUG==============================================

            temp_x = expon * (x_coord(1,iter) - x_coord(1,jter)) / rij
            temp_y = expon * (x_coord(2,iter) - x_coord(2,jter)) / rij
            temp_z = expon * (x_coord(3,iter) - x_coord(3,jter)) / rij

!DEBUG==============================================
!           PRINT *, 'tempxyz',temp_x, temp_y, temp_z
!DEBUG==============================================

            sum_temp_x = sum_temp_x + temp_x
            sum_temp_y = sum_temp_y + temp_y
            sum_temp_z = sum_temp_z + temp_z
          END DO
          grad(1,iter) = 2.0D0* sum_temp_x
          grad(2,iter) = 2.0D0* sum_temp_y
          grad(3,iter) = 2.0D0* sum_temp_z
        END DO

!DEBUG==============================================
!       PRINT *, ''
!       PRINT *, 'repulse gradient'
!       DO iter=1,natoms
!         PRINT *, grad(1,iter), grad(2,iter) ,grad(3,iter)
!       END DO
!       PRINT *, ''
!DEBUG==============================================

        END SUBROUTINE gupta_repulse_gradient

        SUBROUTINE gupta_band_gradient(x_coord,grad)
! calculate gradient of gupta band term
        IMPLICIT NONE
        REAl(KIND=DBL),DIMENSION(:,:),INTENT(IN)  :: x_coord
        REAl(KIND=DBL),DIMENSION(:,:),INTENT(OUT) :: grad
        REAL(KIND=DBL)  :: temp_x, temp_y, temp_z
        REAL(KIND=DBL)  :: sum_temp_x, sum_temp_y, sum_temp_z
        REAL(KIND=DBL)  :: temp
        REAL(KIND=DBL)  :: expon_sum
        REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE  :: band_denom   ! band energy denominator
        REAL(KIND=DBL)  :: rij
        INTEGER         :: natoms
        INTEGER         :: iter, jter
        natoms = SIZE(x_coord,2)
        IF (ALLOCATED(band_denom)) DEALLOCATE(band_denom)
        ALLOCATE(band_denom(natoms))
        grad = 0.0D0
!=======================================================
! Calculate band_denom terms
        DO iter=1,natoms
          expon_sum = 0.0D0
          DO jter=1,natoms
             IF(iter==jter) CYCLE
             rij = distance(iter,jter)
             temp = (rij / r_zero ) - 1.0D0
             temp = temp * 2.0D0 * (-q_ij)
             temp = (eta**2)*DEXP(temp)
             expon_sum = expon_sum + temp
          END DO
          expon_sum = 0.5D0 / DSQRT(expon_sum)
          band_denom(iter) = expon_sum
        END DO
!=======================================================
! Summation for all band terms
        DO iter=1,natoms
          sum_temp_x = 0.0D0
          sum_temp_y = 0.0D0
          sum_temp_z = 0.0D0
          DO jter=1,natoms
            IF(iter==jter) CYCLE
            rij = distance(iter,jter)
            temp = (rij / r_zero ) - 1.0D0
            temp = temp * 2.0D0 * (-q_ij)       ! check here
            temp = (eta**2)*DEXP(temp)
            temp = temp * 2.0D0 * (-q_ij)
            temp = temp / (rij * r_zero)
            temp_x = temp * (x_coord(1,iter) - x_coord(1,jter))
            temp_y = temp * (x_coord(2,iter) - x_coord(2,jter))
            temp_z = temp * (x_coord(3,iter) - x_coord(3,jter))
            sum_temp_x = sum_temp_x + temp_x
            sum_temp_y = sum_temp_y + temp_y
            sum_temp_z = sum_temp_z + temp_z
          END DO
          grad(1,iter) =   sum_temp_x * band_denom(iter) 
          grad(2,iter) =   sum_temp_y * band_denom(iter) 
          grad(3,iter) =   sum_temp_z * band_denom(iter) 
        END DO
!=======================================================
! Summation for remaining band terms
        DO iter=1,natoms
          sum_temp_x = 0.0D0
          sum_temp_y = 0.0D0
          sum_temp_z = 0.0D0
          DO jter=1,natoms
            IF(iter==jter) CYCLE
            rij = distance(iter,jter)
            temp = (rij / r_zero ) - 1.0D0
            temp = temp * 2.0D0 * (-q_ij)       
            temp = (eta**2)*DEXP(temp)
            temp = temp * 2.0D0 * (-q_ij)
            temp = temp / (rij * r_zero)
            temp_x = temp * (x_coord(1,iter) - x_coord(1,jter))
            temp_y = temp * (x_coord(2,iter) - x_coord(2,jter))
            temp_z = temp * (x_coord(3,iter) - x_coord(3,jter))
            sum_temp_x = sum_temp_x + temp_x * band_denom(jter) 
            sum_temp_y = sum_temp_y + temp_y * band_denom(jter) 
            sum_temp_z = sum_temp_z + temp_z * band_denom(jter) 
          END DO
          grad(1,iter) = grad(1,iter) + sum_temp_x
          grad(2,iter) = grad(2,iter) + sum_temp_y
          grad(3,iter) = grad(3,iter) + sum_temp_z
        END DO
!=======================================================
! Add negative sign to band term
        grad = -grad
        END SUBROUTINE gupta_band_gradient

        SUBROUTINE indv_gupta_energy(x_coord,natoms,energy_array)
! calculate gupta energy of individual atoms
        IMPLICIT NONE
        REAl(KIND=DBL),DIMENSION(:,:),INTENT(IN)      :: x_coord
        INTEGER,INTENT(IN)                            :: natoms
        REAL(KIND=DBL),DIMENSION(natoms),INTENT(OUT)  :: energy_array
        REAL(KIND=DBL),DIMENSION(natoms)              :: repulsive_array
        REAL(KIND=DBL),DIMENSION(natoms)              :: band_array

        repulsive_array = 0.0D0
        band_array = 0.0D0

        CALL calc_distance(x_coord)
        CALL indv_gupta_repulsive(natoms,repulsive_array)
        CALL indv_gupta_band(natoms,band_array)
        energy_array = repulsive_array + band_array

        END SUBROUTINE indv_gupta_energy

        SUBROUTINE indv_gupta_repulsive(natoms,energy_array)
! calculate gupta repulsive energy of individual atoms
        IMPLICIT NONE
        INTEGER,INTENT(IN)                      :: natoms
        REAL(KIND=DBL),DIMENSION(:),INTENT(INOUT)  :: energy_array
        REAL(KIND=DBL)                          :: temp, temp_sum
        INTEGER                                 :: iter,jter

        DO iter=1,natoms
          temp_sum = 0.0D0
          DO jter=1,natoms
            IF (iter == jter) CYCLE
            temp = (distance(iter,jter)/ r_zero) - 1.0D0
            temp = temp * (-p_ij)
            temp = a_ij*DEXP(temp)
            temp_sum = temp_sum + temp
          END DO
          energy_array(iter) = temp_sum
        END DO

        END SUBROUTINE indv_gupta_repulsive

        SUBROUTINE indv_gupta_band(natoms,energy_array)
! calculate gupta band energy of individual atoms
        IMPLICIT NONE
        INTEGER,INTENT(IN)                      :: natoms
        REAL(KIND=DBL),DIMENSION(:),INTENT(INOUT)  :: energy_array
        REAL(KIND=DBL)                          :: temp, temp_sum
        INTEGER                                 :: iter,jter

        temp_sum = 0.0D0
        DO iter=1,natoms
          temp_sum = 0.0D0
          DO jter=1,natoms
             IF(iter==jter) CYCLE 
             temp = (distance(iter,jter) / r_zero) - 1.0D0
             temp = temp * 2.0D0*(-q_ij)
             temp = (eta**2)*DEXP(temp)
             temp_sum = temp_sum + temp
          END DO
          energy_array(iter) = -DSQRT(temp_sum)
        END DO

        END SUBROUTINE indv_gupta_band

!       SUBROUTINE old_gupta_band(rij,energy)
!       IMPLICIT NONE
!       REAL(KIND=DBL),DIMENSION(:),INTENT(IN)  :: rij
!       REAL(KIND=DBL),INTENT(OUT)              :: energy
!       REAL(KIND=DBL)                          :: temp, temp_sum
!       INTEGER                                 :: iter
!       INTEGER                                 :: rij_size
!       rij_size = SIZE(rij)
!       temp_sum = 0.0D0
!       DO iter=1,rij_size
!         temp = (rij(iter) / r_zero) - 1.0D0
!         temp = temp * 2.0D0*(-q_ij)
!         temp = (eta**2)*DEXP(temp)
!         temp_sum = temp_sum + temp
!       END DO
!       temp_sum = temp_sum * 2.0D0
!       energy = -DSQRT(temp_sum)
!       END SUBROUTINE old_gupta_band

!       SUBROUTINE combination (n,k,answer)
!       IMPLICIT NONE
!       INTEGER,INTENT(IN)      :: n,k
!       INTEGER,INTENT(OUT)     :: answer
!       INTEGER                 :: temp1, temp2, temp3
!       temp1=factorial(n)
!       temp2=factorial(k)
!       temp3=factorial(n-k)
!       answer = temp1 / (temp2 * temp3)
!       PRINT *, "n and k", n ,k
!       PRINT *, "combination answer" , answer
!       END SUBROUTINE combination 

!       INTEGER FUNCTION factorial(n)
!       IMPLICIT NONE
!       INTEGER,INTENT(IN)      :: n
!       INTEGER                 :: temp1, temp2
!       temp1 = n
!       temp2 = 1
!       DO WHILE( temp1 > 0)
!         temp2 = temp2 * temp1
!         temp1 = temp1 - 1
!       END DO
!       factorial = temp2
!       END FUNCTION factorial

        END MODULE gupta
