#!/usr/bin/env python
import subprocess
import os
import shlex
import numpy as np
import json
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
train_x = 'pts.dat'#'testpoint_v2b_co2h2o.dat' #'pts.dat'
f = open(train_x)
i=0
energy = list()
energy_cout = 0


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
        energy = float(line.split()[0]) #Recond the energy. This can be dipole
        dipole = [float(dip) for dip in line.split()[1:]]

    if lineinconfig >= 3 and lineinconfig <= natom+2:
        atom_n = lineinconfig - 2
        element = line.split()[0]
        coordinate = [float(coor) for coor in line.split()[1:]]#Read cooridnates and turns into float
        atom_coord = [element,coordinate]
        molecule_coord.append(atom_coord)

    if lineinconfig == height: #Reset some of the values
        configuration.append([[[natom],[energy, dipole]],molecule_coord])

        if count_config == 10:
            print('Lines: '+str(i))
            print('Configuration: '+str(count_config))
            break
        count_config = count_config + 1
f.close()

#f=open('result','w')
for a in configuration: #To print in column
    for b in a:
        print(str(b[0][0]))
        print(str(b[1][0]))
        for i in b[1][1]:
            print(i)
a = np.array(configuration)

print (','.join(configuration))
#            f.write(str(c[0])+'\n')
#f.close()

print('test over')
