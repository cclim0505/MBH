#!/usr/bin/python
import numpy as np

f = open('./worker_collected/01_unfinished_lowest_ene.dat','r')
f_array = f.readlines()
f.close()
ene=[]

for row in f_array:
    p = row.split()
    ene.append(float(p[0]))

#print(ene)

print('Lowest energy is: ',np.amin(ene),'eV')
print('from worker no.: ',np.argmin(ene))
