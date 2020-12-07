        PROGRAM main
        USE coord_grad_ene          ,ONLY:atoms, material,coord
     &    ,allocate_coord_gradient ,printout_xyz
        USE geometric_drive         ,ONLY:multi_cage_move
        IMPLICIT NONE
        CHARACTER(LEN=30)       :: input_file
        INTEGER                 :: iter
        INTEGER                 :: f_atoms, f_coord
        INTEGER                 :: num_moves
        CHARACTER(len=20)       :: num_moves_text
        CHARACTER(len=1)        :: dummy

        CALL GET_COMMAND_ARGUMENT(1,input_file)
        CALL GET_COMMAND_ARGUMENT(2,num_moves_text)

        READ(num_moves_text,*) num_moves

        OPEN(NEWUNIT=f_atoms,FILE=TRIM(input_file), STATUS='old')
        READ(f_atoms,*) atoms
        READ(f_atoms,*) 
        READ(f_atoms,*) material
        CLOSE(f_atoms)

        PRINT *, "Atoms and material:", atoms, material

        CALL allocate_coord_gradient

        OPEN(NEWUNIT=f_coord,FILE=TRIM(input_file), STATUS='old')
        READ(f_coord,*)
        READ(f_coord,*)
        DO iter=1,atoms
         READ(f_coord,*) dummy,coord(1,iter),coord(2,iter),coord(3,iter)
        END DO
        CLOSE(f_coord)


        CALL multi_cage_move(num_moves,coord)

        CALL printout_xyz('moved_cage.xyz',coord)

        PRINT *, ""
        PRINT *, "Check output structure in 'moved_cage.xyz'."
        PRINT *, ""

        END PROGRAM main
