#Set the output to a png file
set terminal png size 500,500

#Output file nameo
set output 'MC_overview.png'

#Graphic title
set title 'MC Steps taken'

#Plot bar chart
set xrange[-1:]
set yrange[0:600] 
set boxwidth 0.7
set style fill solid
plot 'MC_overview.dat' with boxes
