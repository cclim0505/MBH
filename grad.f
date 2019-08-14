!
!  Energy and gradient for Many-Body Potential.
!  Developed by senior
      SUBROUTINE agrad(x,grad)
      USE constants,         ONLY:DBL
      USE initialise,             ONLY:atoms
      USE param
      IMPLICIT NONE
      INTEGER        :: J1, J2, J3, J4
      REAL(KIND=DBL) :: ZETA,ESILON,P,Q,RZERO
      REAL(KIND=DBL),DIMENSION(atoms*3) :: X,GRAD
      REAL(KIND=DBL),DIMENSION(atoms)   :: GRAD_ATT
      REAL(KIND=DBL) :: ELJ_ATT,XMUL3
      REAL(KIND=DBL) :: DIST,dist1
      REAL(KIND=DBL),DIMENSION(atoms,atoms):: GRAD_ATTM,GRAD_REPM
      REAL(KIND=DBL) :: DUMMYX, DUMMYY, DUMMYZ, XMUL2, rcut
      INTEGER :: natoms_a
      natoms_a = atoms
                                
        DO J1=1,ATOMS
         J3=3*J1
         DIST=X(J3-2)**2.+X(J3-1)**2.+X(J3)**2.
         GRAD_REPM(J1,J1)=0.0D0
         GRAD_ATTM(J1,J1)=0.0D0
         ELJ_ATT=0.0D0
         DO J2=1,ATOMS
            J4=3*J2
            IF(J1.NE.J2)THEN
              DIST=(X(J3-2)-X(J4-2))**2.+(X(J3-1)-X(J4-1))**2.
     + +(X(J3)-X(J4))**2.
              DIST=DSQRT(DIST)
              IF ((J1<=natoms_a) .and. (J2<=natoms_a) ) THEN
                 P=P_a
                 Q=Q_a
                 ESILON=EPSILON_a
                 ZETA=ZETA_a
                 RZERO=RZERO_a
              ELSE IF ( (J1 > natoms_a ) .and. (J2 > natoms_a) ) THEN
                 P=P_b
                 Q=Q_b
                 ESILON=EPSILON_b
                 ZETA=ZETA_b
                 RZERO=RZERO_b
                ELSE
                  P=P_ab
                  Q=Q_ab
                  eSILON=EPSILON_ab
                  ZETA=ZETA_ab
                  RZERO=RZERO_ab
                ENDIF
              dist1=dist/rzero
              RCUT=3D0*rzero*atoms**(1D0/3D0)
              GRAD_REPM(J1,J2)=  ESILON*(-P/rzero)*DEXP(P*(1.-DIST1))

!ccccccccccccccccccccccc 
              if (dist >rcut ) then 
!               write(*,*)'whoooooooooooooooopppppppppppssssssssss'
             grad_repm(j1,j2)=grad_repm(j1,j2)+esilon*(-p/rzero)*3D0*
     +      (rcut-dist)
!                  wall=wall+1.
              end if
!ccccccccccccccccccccccccccc
              grad_repm(j1,j2)=grad_repm(j1,j2)/DIST
              PRINT *, 'grad_repm of agrad', j1, j2, grad_repm(j1,j2)
!              GRAD_ATTM(J1,J2)=(ZETA**2)*(-2D0*Q/rzero)*DEXP(2.*Q*(1.-DIST1))/DIST
              GRAD_ATTM(J1,J2)=(ZETA**2)*(-2D0*Q/rzero)
              GRAD_ATTM(J1,J2)=GRAD_ATTM(J1,J2)*DEXP(2.*Q*(1.-DIST1))
              GRAD_ATTM(J1,J2)=GRAD_ATTM(J1,J2)/DIST
              ELJ_ATT=ELJ_ATT+GRAD_ATTM(J1,J2)*((DIST)/(-2.*Q/rzero))
            ENDIF
         ENDDO
         GRAD_ATT(J1)=DSQRT(ELJ_ATT)
      ENDDO

      PRINT *, 'GRAD_ATT array of agrad'
      DO j1=1,atoms
        PRINT *, GRAD_ATT(j1)
      END DO

      PRINT *, 'GRAD_ATTM array of agrad'
      DO j1=1,atoms
        DO j2=1,atoms
          PRINT *, j1,j2,GRAD_ATTM(j1,j2)
        END DO
      END DO


      DO J1=1,ATOMS
         J3=3*J1
         DUMMYX=0.0D0
         DUMMYY=0.0D0
         DUMMYZ=0.0D0
         DO J2=1,ATOMS
           J4=3*J2
           XMUL2=GRAD_REPM(J1,J2)-(GRAD_ATTM(J1,J2)/(2.*GRAD_ATT(J1)))
           XMUL3=GRAD_REPM(J2,J1)-(GRAD_ATTM(J2,J1)/(2.*GRAD_ATT(J2)))
           DUMMYX=DUMMYX+(XMUL2+XMUL3)*(X(J3-2)-X(J4-2))
           DUMMYY=DUMMYY+(XMUL2+XMUL3)*(X(J3-1)-X(J4-1))
           DUMMYZ=DUMMYZ+(XMUL2+XMUL3)*(X(J3)  -X(J4))
         ENDDO
         DIST=X(J3-2)**2.+X(J3-1)**2.+X(J3)**2.
         GRAD(J3-2)=DUMMYX
         GRAD(J3-1)=DUMMYY
         GRAD(J3)=DUMMYZ
      ENDDO
!      open(50,file='check_grad1.dat',status='old')
!      write(50,*)grad

      PRINT *, "================="
      PRINT *, "grad values of agrad"
      DO j1=1,atoms*3
         PRINT *, j1, grad(j1)
      END DO
      PRINT *, ""

      END SUBROUTINE agrad
