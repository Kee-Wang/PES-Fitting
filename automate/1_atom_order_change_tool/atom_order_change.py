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

#train_x = raw_input('Please input the data file name: ')
#train_x ='testpoint_v2b_co2h2o.dat' #'pts.dat'
#train_x = 'dep_pts.dat'

class configs:
    '''Read configuraitons and do stuff.

        Input::train_x = file contains configuraitons in format.
                1. Sample:

                    6
                       -0.00476416    -0.90668842     0.00000000     0.00000000
                    O   1.39511118877395       -1.15984255934314       4.127641665696761E-015
                    O   1.39511118877394        1.15984255934314      -1.010148437748930E-014
                    O  -1.37600000000000      -4.445978668838478E-015  2.222989334419239E-015
                    H  -1.99243331768939      -0.734636620951108       6.321962738314389E-015
                    H  -1.99243331768940       0.734636620951091       1.574609909686252E-015
                    C   1.37600000000000       0.000000000000000E+000  2.222989334419239E-015

                2. Mapping to parameters:




                    Breakdown:
                        configs =

                            [
                                [config],    #configs[0], first configuraiton
                                [config],    #configs[1], second configuraiton
                                ...
                                [config]     #configs[-1], last configuraiton
                            ]


                        config = configs[0] =    #Take first configuraiton for example
                            [
                                [molecule_count_total],    #number of atoms in a molecule
                                [[energy], [dipole]]
                                [molecule]
                            ]
                                    #molecule_count_total (int):    configs[0][0] = 6
                                    #energy (float)(Hartree):    configs[0][1][0] = -0.00476416
                                    #dipole (float)(a.u.): configs[0][1][1][0]=-0.90668842


                        molecue = configs[0][2] = config[2] =
                            [
                                [[element],   [x,   y,  z]],   #First atom
                                [[element],   [x,   y,  z]],   #Second atom
                                ...
                                [[element],   [x,   y,  z]]    #Last atom
                            ]
                                    #element  (string) : configs[0][2][0][0]
                                    #x/y/z (float)(Angstrom):   configs[0][2][0][1][0] = '1.39511118877395'

    '''
    def __init__(self,train_x,dip = False,first_n_configs=False):
        '''Read input file into configs with checks

        1. Result:  configs[a][b][c][d][e], type = List
            configs stores the all the input obtained from train_x
                User's guide:
                [a]:  The ([a]+1)_th configuration #Python list count from 0

                    [b]=0: Number of atoms

                    [b]=1:

                        [c]=0: energy
                        [c]=1: dipole

                            [d]=0: dipole_x
                            [d]=1: dipole_y
                            [d]=2: dipole_z

                    [b]=2: molecule

                        [c]: The ([c]+1)_th atom

                            [d]=0: element
                            [d]=1: cartisian coordinate

                                [e]=0: x
                                [e]=1: y
                                [e]=2: z

        2. Check input:
            While reading files, it can check:
            * Number of atoms: existence, type and consistency
            * Energy: existence and type
            * Dipole (if dip = True): existence, type and dimension
            * Element: existence, type and consistency
            * Cartisian coordiante: existence, type and dimension
            Notice: consistency error includes the existence error because not existent is also considered as inconsistent.
        '''
        f = open(train_x)
        configs_count = 0
        line_count = 0
        element_1 = list()
        configs = list()
        dipole = list()
        for line in f:
            line = line.strip() #The delete the newline character
            line_count = line_count + 1
            if line_count == 1: # This is number of atoms.  #Check type
                if line.isdigit():
                    molecule_count_total = int(line)
                    #break
                    config_count_totol = molecule_count_total + 2
                else:
                    print("Type error: molecule number in line "+str(line_count)+" is : "+line)#Check type

            config_count = (line_count - 1) % config_count_totol + 1

            if configs_count == 1 and config_count >= 3 and config_count <= molecule_count_total+2:
                element_1.append(line.split()[0]) #Read elements in first config

            if config_count == 1:
                molecule_coord = list() #Renew molecue_coord each configuration
                configs_count = configs_count + 1
                try:
                    if int(line) != molecule_count_total:
                        print('Consistency error: number of atoms in line: ' + str(line_count)+' is : ' + line) #Check consistency
                except:
                    print("Type error: number of atoms in line "+line_count+" is : "+line) #Check input
                    #break

            if config_count == 2:
                try:
                    energy = [float(line.split()[0])] #Record the energy.
                except:
                    print("Type error: energy in line "+line_count+" is : "+line) #Check energy type
                    #break


                try:
                    if dip == True and len(dipole) !=3:
                        dipole = [float(dip) for dip in line.split()[1:]] #Record dipole.
                        print("Dimensino error: Dipole dimension in line "+str(line_count)+" is : "+str(len(dipole)))
                except:
                    print("Type error: Dipole in line "+str(line_count)+" is: "+line) #Check energy type
                    #break


            if config_count >= 3 and config_count <= molecule_count_total+2:
                molecule_count = config_count - 2
                element = line.split()[0]

                if element != element_1[molecule_count-1] and len(element) != 0: # Element check type and existence
                    print('Consistency error: atom in line: ' + str(line_count)+' is: '+line)

                coordinate = [float(coor) for coor in line.split()[1:]]#Read cooridnates and turns into float
                atom_coord = [element,coordinate]
                molecule_coord.append(atom_coord)

            if config_count == config_count_totol: #Reset some of the values
                configs.append([[molecule_count_total],[energy, dipole],molecule_coord])


                if configs_count == first_n_configs: #Only read first n configs
                    break



        self.configs = configs
        self.dip = dip
        self.configs_count = configs_count



        configs_sorted = self.energy_sort(configs)
        energy_lowest=configs_sorted[0][1][0][0]
        energy_highest=configs_sorted[-1][1][0][0]

        print('''
        Reading summary:

            Inputfile: '''+train_x+'''
            Number of atoms: '''+str(molecule_count_total)+'''
            Number of configurations: '''+str(configs_count)+'''
            Number of lines in file: '''+str(line_count)+'''
            Lowest energy (Hartree): {:14.8f}
            Highest energy (Hartree): {:14.8f}

        Reading finished.
        '''.format(energy_lowest, energy_highest)
            )




        f.close()

    def write_print(self,configs):
        """Print configs into screen.

        """
        for config in configs: #Write config by config
            print('{:<2d}'.format(config[0][0])) #Number of atoms. Align number of atom to the very left
            print(' '),# To align the colums of energy and coordiante
            if self.dip:#Print option for having dipole or not
                print('{:14.8f}{:14.8f}{:14.8f}{:14.8f}'.format(config[1][0][0],config[1][1][0],config[1][1][1],config[1][1][2]))
            else:#Print also dipole (if input file has it)
                print('{:14.8f}'.format(config[1][0][0]))#Print energy only
            for molecule in config[2]:#The atom part: coordiante of atoms
                print('{} {:14.8f}{:14.8f}{:14.8f}'.format(molecule[0],molecule[1][0],molecule[1][1],molecule[1][2]))

    def write_file(self,configs,filename):
        """Write formated configurations into filename.

        """
        f=open(filename,'w')
        for config in configs: #Write config by config
            f.write('{:<2d}'.format(config[0][0])+'\n') #Number of atoms. Align number of atom to the very left
            f.write('  '), #To align the colums of energy and coordiante. Comma to continue on the same line
            if self.dip:#Print energy and dipole
                f.write('{:14.8f}{:14.8f}{:14.8f}{:14.8f}'.format(config[1][0][0],config[1][1][0],config[1][1][1],config[1][1][2])+'\n')
            else:#Print energy only
                f.write('{:14.8f}'.format(config[1][0][0])+'\n')#Print energy only
            for molecule in config[2]:#The atom part: coordiante of atoms
                f.write('{} {:14.8f}{:14.8f}{:14.8f}'.format(molecule[0],molecule[1][0],molecule[1][1],molecule[1][2])+'\n')
        print('Configs are written to the file: '+filename)
        f.close()

    def energy_sort(self,configs,reverse=False): #sort based on energy
        if reverse: #From high to low
            configs.sort(key= lambda item:item[1][0],reverse = True) #In key, it helps iterate though list.
            return configs
        else: #Default, from low to high.
            configs.sort(key= lambda item:item[1][0])
            return configs

    def configs_list(self):
        return self.configs

    def energy_threshold(self,configs,lower=False,upper=False):


        print('Excuting energy_threshold:----------------')
        configs = self.energy_sort(configs)
        configs_new_count = 0
        configs_new = list()
        energy_lowest=configs[0][1][0][0]
        energy_highest=configs[-1][1][0][0]


        if lower is not False:
            energy_lowest = float(lower)
        if upper is not False:
            energy_highest = float(upper)
        if energy_lowest > energy_highest:
            print('Boundary error: lower energy bound is greater then upper energy bound. Returning origianl configurations.')
            self.error_message()
            return configs
        else:
            for config in configs:
                if config[1][0][0] >= energy_lowest and config[1][0][0] <= energy_highest:
                    configs_new.append(config)
                    configs_new_count = configs_new_count + 1
            percentage = round(float(configs_new_count)/float(self.configs_count)*100,2)
            print('There are '+str(configs_new_count)+'/'+str(self.configs_count)+'('+str(percentage)+'%) configs between '+str(energy_lowest)+' and '+ str(energy_highest)+' Hartree')
            print('End of energy_threshold-------------------')
            return configs_new

    def error_message(self):
        """To show error message and waiting for decision

        """
        decision = raw_input('Do you want to continue? (y/n) \n')
        if decision is not 'y':
            print('Exiting program.')
            exit()




#a = configs(train_x)
#list1 = a.configs_list()

#a.write_print(list1)
#a.write_file(list1,'list1test')
#list_sorted = a.energy_sort(list1,reverse=True)

#list_cut = a.energy_threshold(list1,lower=0.0)
#a.write_file(list_cut,'v2bgr0.abs')
        #Calculate the distance betwo atom a and b.
        #         print(molecule[2][2])
        #         atom1=np.array(molecule[2][2][1])
        #         atom2=np.array(molecule[2][5][1])
        #         dis = math.sqrt(np.sum(np.square(atom1-atom2)))
        #     if dis >=8:
        #         i=i+1
        #         print(dis)
        # print(i)






print('test over')
