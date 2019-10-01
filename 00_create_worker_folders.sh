#!/bin/bash
# 
#  	Script to create worker folders for MBH optimization
#	in parallel mode
# 
# 

# Input variables

workers=1		# the number of parallel processes intended
tar_folder='00_worker_folder2copy'

#===================================================================

echo 'number of workers is' $workers

#===================================================================


for (( iter=0; iter<$workers; iter++))
do
	if [ $iter -lt 10 ];then
		out_dir='worker00'$iter
	elif [ $iter -lt 100 ];then
		out_dir='worker0'$iter
	else
		out_dir='worker'$iter
	fi
	cp -r $tar_folder $out_dir/
done

