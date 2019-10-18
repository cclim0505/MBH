#!/bin/bash

dir='00_check_MC_plot'
output='MC_overview.dat'
graph='MC_overview.png'
tail worker???/01_resume_MC_step.dat

cat worker???/01_resume_MC_step.dat > $output
gnuplot ./$dir/plot_MC.gp
display $graph
mv $graph ./$dir
mv $output ./$dir
