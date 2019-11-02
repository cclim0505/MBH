#!/bin/bash

prod='99_prod_MBH'

##gupta='param_Gupta'
##dftb='param_DFTBplus_skf'

mkdir $prod
cp *.sh ./$prod
cp *.out ./$prod
cp README* ./$prod
cp -r 00_worker* ./$prod
cp -r 00_check* ./$prod
cp -r param* ./$prod/
cp -r *session* ./$prod

mv $prod ../
