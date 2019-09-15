#!/bin/bash

out_folder='worker_collected'
out_coord='lowest_coord.xyz'
out_ene='lowest_ene.dat'
lowest_coord='00_lowest_coord.xyz'
lowest_ene='00_lowest_energy.dat'

rm -r $out_folder
mkdir $out_folder


cat worker???/$lowest_coord > ./$out_folder/$out_coord
cat worker???/$lowest_ene > ./$out_folder/$out_ene
