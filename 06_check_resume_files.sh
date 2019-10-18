#!/bin/bash

## declare test function
function test_file_exist(){
	if test -s "$1"; then
		echo "OK: $1 "
	else
		echo "NIL: $1				[x]"
		let "counter++"
	fi	
}



## declare test resume function
function test_resume(){
	for iter in $FILE1 $FILE2 $FILE3 $FILE4
	do
		test_file_exist "$iter"
	done
}

### 		    			###
### 		MAIN			###
### 		    			###

counter=0


for num in {000..019}
do
	FILE1=01_resume_MC_step.dat
	FILE2=01_resume_old_coord.xyz
	FILE3=01_resume_lowest_coord.xyz
	FILE4=01_resume_lowest_ene.dat
	dir=worker
	dir=$dir$num
	echo $dir
	FILE1=./$dir/$FILE1
	FILE2=./$dir/$FILE2
	FILE3=./$dir/$FILE3
	FILE4=./$dir/$FILE4
	test_resume		## call function
done

echo
echo "Number of resume file errors: "$counter


### 		    			###
### 		END MAIN		###
### 		    			###








##dir=worker
##worker=000
##dir=$dir$worker

## assign file names to be tested
##FILE1=01_resume_MC_step.dat
##FILE2=01_resume_old_coord.xyz
##FILE3=01_resume_lowest_coord.xyz
##FILE4=01_resume_lowest_energy.dat
##FILE1=./$dir/$FILE1
##FILE2=./$dir/$FILE2
##FILE3=./$dir/$FILE3
##FILE4=./$dir/$FILE4
##test_resume		## call function
