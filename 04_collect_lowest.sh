#!/bin/bash

natoms=$(awk 'NR==1 {print $2}' ./session_param/session_in.dat)
let natoms=$natoms+2
alldir='worker???'

out_folder='worker_collected'
out_coord='lowest_coord.xyz'
out_ene='lowest_ene.dat'
lowest_coord='00_lowest_coord.xyz'
lowest_ene='00_lowest_energy.dat'

final_folder='worker_final'


# Deletes previous pevious collected folder, updates with new one
rm -r $out_folder
rm -r $final_folder
mkdir $out_folder
mkdir $final_folder


cat $alldir/$lowest_coord > ./$out_folder/$out_coord
cat $alldir/$lowest_ene > ./$out_folder/$out_ene


notdone_ene='01_resume_lowest_ene.dat'
ene_out='01_unfinished_lowest_ene.dat'
notdone_coord='01_resume_lowest_coord.xyz'
coord_out='01_unfinished_last_coord.xyz'

# Copy all lowest coord and energies to target folder.
cat $alldir/$notdone_coord > ./$out_folder/$coord_out
cat $alldir/$notdone_ene > ./$out_folder/$ene_out


# Determine the lowest energy structure with python script.
./04_lowest.py
