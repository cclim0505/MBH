        SUBROUTINE DFTB_force(x,force)
        USE mod_constants,    ONLY:DBL
        USE mod_atoms,        ONLY:atoms
        USE mod_material,     ONLY:material

        IMPLICIT NONE
        INTEGER                           :: j1,j2
        INTEGER                           :: nline     ! number of lines to skip
        CHARACTER(LEN=1)                  :: skip      ! dummy 
        REAL(KIND=DBL),DIMENSION(atoms*3) :: x,force
!        REAL(KIND=DBL)                    :: energy
        REAL(KIND=DBL)                    :: au_r0     ! resizing constant
        REAL(KIND=DBL)                    :: au_f0     ! force : a.u. to meV 
        CHARACTER(LEN=5)                  :: geometry  ! DFTB+ geometry input
        CHARACTER(LEN=2)                  :: mat       ! material of Pure cluester(short form)
        CHARACTER(LEN=8)                  :: name_atoms
        CHARACTER(LEN=1)                  :: name_atoms2,name_atoms3


!===================================================
! Initialise values and constants
!===================================================
        mat = material(1:2)
        au_r0 = 2.884
        au_f0 = 51.42208619083232
        geometry = 'C'                     ! c = cluster
        name_atoms2=char(mod (atoms/10,10)+48)
        name_atoms3=char(mod (atoms/1,10)+48)
        name_atoms=name_atoms2//name_atoms3//'.gen'
        
! Determining the number of lines to skip in restuls.tag
        nline = 6 + (atoms/3) 
        IF(MOD(atoms,3) == 0) THEN
          nline = nline
        ELSE
          nline = nline + 1
        END IF
!        WRITE(*,*) mat, name_atoms
!        WRITE(*,*) material
!===================================================
! Convert coordinates into atomic units
!===================================================
!        DO j1=1,atoms*3
           
!        x(j1) = x(j1) / au_r0

!        END DO
!===================================================
! Write coordinates into .gen file
!===================================================
        OPEN(20, file=mat//name_atoms, status='replace')

        WRITE (20,*) atoms,"   " , geometry
        WRITE (20,*) mat

        DO j1=1,atoms
           j2=j1*3
           WRITE (20,*) j1, 1, x(j2-2), x(j2-1), x(j2)
        END DO
        CLOSE(20)

!===================================================
! Run DFTB+
!===================================================
        CALL SYSTEM ('dftb+ > dftb_force.out')

!===================================================
! Read forces from Results.tag
!===================================================
        OPEN(21, file='results.tag', status='old')

        DO j1=1,nline              ! skipping file lines
            READ (21,*) skip
        END DO


        DO j1=1,atoms
           j2=j1*3
            READ (21,*)  force(j2-2), force(j2-1), force(j2)
        END DO

        CLOSE(21)


!Check the forces read from 'results.tag'
!        OPEN(22, file='checkDFTBforce.dat', access='append')
!        DO j1=1,atoms
!           j2=j1*3
!            WRITE (22,'(I3, 3F20.10)')j1, force(j2-2), force(j2-1),
!     1 force(j2)
!            WRITE (22,*)  '    '
!        END DO
!        CLOSE(22)


!===================================================
! Convert forces into eV/Angstrom
!===================================================

        DO j1=1,atoms*3
           force(j1) = force(j1) * au_f0
        END DO


        END SUBROUTINE DFTB_force


             

