#!/bin/bash
# get lowest energy structure

dir=$1
lowest='01_resume_lowest_coord.xyz'
target='999_lowest'
final='lowest_coord.xyz'

ene='01_resume_lowest_ene.dat'
fene='lowest_ene.dat'

rm -r $target
mkdir $target
cp ./$dir/$lowest ./$target/$final
cp ./$dir/$ene ./$target/$fene


