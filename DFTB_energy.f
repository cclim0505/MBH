        SUBROUTINE DFTB_energy(x,energy)
        USE mod_constants,    ONLY:DBL
        USE mod_atoms,        ONLY:atoms
        USE mod_material,     ONLY:material
        IMPLICIT NONE
        INTEGER                           :: j1,j2
        INTEGER                           :: nline     ! number of lines to skip
        CHARACTER(LEN=1)                  :: skip      ! dummy 
        REAL(KIND=DBL),DIMENSION(atoms*3) :: x
        REAL(KIND=DBL)                    :: energy
        REAL(KIND=DBL)                    :: au_r0     ! resizing constant
        REAL(KIND=DBL)                    :: au_e0     ! force : a.u. to meV 
        CHARACTER(LEN=5)                  :: geometry  ! DFTB+ geometry input
        CHARACTER(LEN=2)                  :: mat       ! material of Pure cluester(short form)
        CHARACTER(LEN=8)                  :: name_atoms
        CHARACTER(LEN=1)                  :: name_atoms2,name_atoms3

!===================================================
! Initialise values and constants
!===================================================
        mat = material(1:2)
        nline = 1
        au_r0 = 2.884
        au_e0 = 27.2113699175    
        geometry = 'C'                     ! c = cluster
        name_atoms2=char(mod (atoms/10,10)+48)
        name_atoms3=char(mod (atoms/1,10)+48)
        name_atoms=name_atoms2//name_atoms3//'.gen'
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
        CALL SYSTEM ('dftb+ > dftb_energy.out')

!===================================================
! Read forces from Results.tag
!===================================================
        OPEN(21, file='results.tag', status='old')

        DO j1=1,nline              ! skipping file lines
            READ (21,*) skip
        END DO

            READ (21,*) energy 

        CLOSE(21)


!Check the forces read from 'results.tag'
!        OPEN(22, file='checkDFTB_energy.dat', access='append')
!            WRITE (22,'(F20.10)') energy
!        CLOSE(22)


!===================================================
! Convert Energy into meV
!===================================================

        energy = energy * au_e0

        END SUBROUTINE DFTB_energy


             

