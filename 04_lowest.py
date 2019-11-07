#!/usr/bin/python
"""
Script to find the lowest energy structure

Script finds the lowest energy structure
among all worker folders and outputs
it in a separate file.
"""

import numpy as np

collect_folder = 'worker_collected'
final_folder = 'worker_final'
lowest_file = collect_folder + '/01_unfinished_lowest_ene.dat'

f = open(lowest_file,'r')
f_array = f.readlines()
f.close()
ene=[]

for row in f_array:
    p = row.split()
    ene.append(float(p[0]))

low_ene = np.amin(ene)
low_index = np.argmin(ene)

print('Lowest energy is: ', low_ene, ' eV')
print('from worker no.: ', low_index)

# Save energy, index and coordinates of lowest
# in a file.


if low_index < 10:
    low_folder = "worker00" + str(low_index)
elif low_index < 100:
    low_folder = "worker0" + str(low_index)

#print(low_folder)

lowfile = str(low_folder) + '/01_resume_lowest_coord.xyz'

#print(lowfile)

f = open(lowfile, 'r')
struct = f.readlines()
f.close()


outfile = 'lowest_struct.xyz'
#struct[1] = 'ene (ev) = ' + str(low_ene) \
#	    + '\tfrom ' + low_folder + '\n'

struct[1] = struct[1][:-1]
struct[1] = struct[1] + ' from ' + low_folder + '\n'

f = open(final_folder + '/' + outfile, 'w')
for line in struct:
    f.write(line)
f.close()


print("Final result in " + final_folder)
