#MBH
Development project for Modified Basin-Hopping program


## HOW TO USE
This is the source folder for MBH with Gupta and DFTB potential.

## Makefile
Make file to compile when necessary

## Creating or copying into a production folder


Run 20_create_prod_dir.sh to create a production directory/folder.
Run your jobs from there.
The following are copied to the production folder:
executable, run.out
session_param folder
folders that start with "00?"


# Version update and changes

1.0
26-01-2021
Multiple versions previously. This is the version that has the first version changes documentation.
+ Added version number banner subroutine.
+ Multiple changes and updates to the geometric drive (ring & cage) MC moves.
+ Improved random generation of initial configuration.
+ Added features to random generation: adjustable minimumm atomic distance, random reset and its frequency
+ Added 2D initial configuration feature.
+ Added angular-displacement on-off feature for basin hopping.
+ Added script 08_get_all_local_coords.sh
