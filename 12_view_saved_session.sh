#!/bin/bash

dir='saved_session'

session='saved_session_in.dat'
monte='saved_param_MC.dat'
basin='saved_param_BH.dat'
random='saved_param_random.dat'

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
