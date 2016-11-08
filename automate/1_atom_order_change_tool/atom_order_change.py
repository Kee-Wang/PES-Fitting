#!/usr/bin/env python
import subprocess
import os
import shlex
import numpy as np
import math
def cl(command):
    #ip::string, command line as string input
    #op::string, return value is the output of command line
    #Notice, each time when change directly, cl starts from currect directory.
    #Use three \' if you want to input multiple line
    arg = shlex.split(command)
    p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    (output, err) = p.communicate()
    print(output)
    return output

configuration = list()
#train_x = raw_input('Please input the data file name: ')
train_x = 'testpoint_v2b_co2h2o.dat' #'pts.dat'
f = open(train_x)
i=0


element_1=list()

count_config = 1
for line in f:
    i=i+1
    if i==1:
        natom = int(line)
        height = natom + 2

    lineinconfig = (i-1) % height + 1

    if lineinconfig == 1:
        molecule_coord = list() #Renew molecue_coord each configuration

    if lineinconfig == 2:
        energy = [float(line.split()[0])] #Recond the energy. This can be dipole
        dipole = [float(dip) for dip in line.split()[1:]]

    if lineinconfig == 1:
        if int(line) != natom:
            print('Number of atoms is not consistent for configuration'+str(count_config)+' in line: ' + str(i))


    if count_config==1 and lineinconfig >= 3 and lineinconfig <= natom+2:
        element_1.append(line.split()[0])

    if lineinconfig >= 3 and lineinconfig <= natom+2:
        atom_n = lineinconfig - 2
        element = line.split()[0]
        if element != element_1[atom_n-1]:
            print('Symmetry order not consistent for configuration '+str(count_config)+' in line: ' + str(i))
        coordinate = [float(coor) for coor in line.split()[1:]]#Read cooridnates and turns into float
        atom_coord = [element,coordinate]
        molecule_coord.append(atom_coord)

    if lineinconfig == height: #Reset some of the values
        configuration.append([[natom],[energy, dipole],molecule_coord])

#        if count_config == :
#print('Lines: '+str(i))
#print('Configurations: '+str(count_config))
#        break
        count_config = count_config + 1
f.close()
print('Lines: '+str(i))
print('Configurations: '+str(count_config-1))

i=0
for molecule in configuration: #To print in column
    # print('{:<2d}'.format(molecule[0][0])) #Number of atoms. Align number of atom to the very left
    # print(' '),# To align the colums of energy and coordiante
    # if len(molecule[1][1]) == 0:#Print option for having dipole or not
    #     print('{:14.8f}'.format(molecule[1][0][0]))#Print energy only
    # else:#Print also dipole (if input file has it)
    #     print('{:14.8f}{:14.8f}{:14.8f}{:14.8f}'.format(molecule[1][0][0],molecule[1][1][0],molecule[1][1][1],molecule[1][1][2]))
#    for atom in molecule[2]:#The atom part: coordiante of atoms

#        print('{} {:14.8f}{:14.8f}{:14.8f}'.format(atom[0],atom[1][0],atom[1][1],atom[1][2]))
#    print(molecule[2][2])
    atom1=np.array(molecule[2][2][1])
    atom2=np.array(molecule[2][5][1])
    dis = math.sqrt(np.sum(np.square(atom1-atom2)))
    if dis >=8:
        i=i+1
        print(dis)
print(i)


# f=open('result','w')
# for molecule in configuration: #To print in column
#     f.write('{:<2d}'.format(molecule[0][0])+'\n') #Number of atoms. Align number of atom to the very left
#     f.write('  '),# To align the colums of energy and coordiante
#     if len(molecule[1][1]) == 0:#Print option for having dipole or not
#         f.write('{:14.8f}'.format(molecule[1][0][0])+'\n')#Print energy only
#     else:#Print also dipole (if input file has it)
#         f.write('{:14.8f}{:14.8f}{:14.8f}{:14.8f}'.format(molecule[1][0][0],molecule[1][1][0],molecule[1][1][1],molecule[1][1][2])+'\n')
#     for atom in molecule[2]:#The atom part: coordiante of atoms
#         f.write('{} {:14.8f}{:14.8f}{:14.8f}'.format(atom[0],atom[1][0],atom[1][1],atom[1][2])+'\n')
# f.close()

print('test over')
