#!/bin/bash

dir='session_param'

session='session_in.dat'
monte='param_MC.dat'
basin='param_BH.dat'
random='param_random.dat'

echo ">>>>>>>>>>"$session"<<<<<<<<<<"
head -15 $dir/$session
echo

echo ">>>>>>>>>>"$random"<<<<<<<<<<"
head -5 $dir/$random
echo

echo ">>>>>>>>>>"$monte"<<<<<<<<<<"
head -5 $dir/$monte
echo

echo ">>>>>>>>>>"$basin"<<<<<<<<<<"
head -5 $dir/$basin
