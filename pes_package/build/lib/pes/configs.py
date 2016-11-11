#!/usr/bin/env python
#import configs

"""TODO:
        1. Add info method
        2. Give hint on type, consistency checking
        3. Add plot distribution method (urgent)
        4. Add random selection method
        5. In __init__, check the same config.
        """

class configs():
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
        self.train_x = train_x
        self.molecule_count_total = molecule_count_total
        self.line_count = line_count

        self.configs_sorted = self.sort(configs)
        self.energy_lowest=self.configs_sorted[0][1][0][0]
        self.energy_highest=self.configs_sorted[-1][1][0][0]

        print('File reading finshed.\n')
        self.info()
        self.help()


        f.close()

    def info(self,configs=False):


        print('============================')
        print('||          INFO          ||')
        print('============================')

        print('\nInputfile:  {}\nNumber of atoms:    {:<2d}\n'.format(self.train_x, self.molecule_count_total))
        self.order()
        print('Number of configurations:    {:<6d}\nNumber of lines in file:    {:<2d}\nLowest energy (Hartree):    {:14.8f}\nHighest energy (Hartree):   {:14.8f}\n'.format(self.configs_count,self.line_count,self.energy_lowest, self.energy_highest))
        print('========END OF INFO=======\n')

    def order(self,configs=False):
        if configs is False:
            configs = self.configs

        print('Ordering: \n')
        molecule_count = 1
        for atom in configs[0][2]:
            print('{} ({})'.format(atom[0],molecule_count))
            molecule_count = molecule_count + 1



    def help(self):
        syntax="""
        =============================
        ||       SYNTAX HELP      ||
        =============================

        Suppose you have:
            a = configs('input.xyz')

        Now you can:

            a.prt(configs=False)
            a.write('output.xyz', configs=False)
            a.sort(configs=False,reverse=False,key ='energy')
            a.list(first_n_configs=False)
            a.threshold_energy(configs,lower=False,upper=False)
            a.distance(config,atom1,atom2)
            a.plot()
            a.order()

        *To repeat info:
            a.info()

        *To repeat syntax help:
            a.help()

        =======END OF HELP SYNTAX====

        """
        #for line in syntax.split('\n'):
        #            print(line.strip())
        print(syntax)

    def prt(self,configs=False):
        """Print configs into screen.

        """
        configs_count=0

        if configs is False:
            print('Using default input configurations')
            configs = self.configs
        try:
            temp = configs[0][0][0]
        except:#If has to configuration, this makes sure it can be safely iterated in the next statemnt.
            configs = [configs]
            configs_count = configs_count + 1
        for config in configs: #Write config by config
            configs_count = configs_count + 1
            print('{:<2d}'.format(config[0][0])) #Number of atoms. Align number of atom to the very left
            print(' '),# To align the colums of energy and coordiante
            if self.dip:#Print option for having dipole or not
                print('{:14.8f}{:14.8f}{:14.8f}{:14.8f}'.format(config[1][0][0],config[1][1][0],config[1][1][1],config[1][1][2]))
            else:#Print also dipole (if input file has it)
                print('{:14.8f}'.format(config[1][0][0]))#Print energy only
            for molecule in config[2]:#The atom part: coordiante of atoms
                print('{} {:14.8f}{:14.8f}{:14.8f}'.format(molecule[0],molecule[1][0],molecule[1][1],molecule[1][2]))
        print('Printed {} configurations.'.format(configs_count))

    def write(self,filename, configs=False): #Default value has to be after non-default value.
        """Write formated configurations into filename.

        """
        f=open(filename,'w')
        configs_count=0

        if configs is False:
            print('Using default input configurations')
            configs = self.configs
        try:
            temp = configs[0][0][0]
        except:#If has to configuration, this makes sure it can be safely iterated in the next statemnt.
            configs = [configs]
            configs_count = configs_count + 1
        for config in configs: #Write config by config
            configs_count = configs_count + 1
            f.write('{:<2d}'.format(config[0][0])+'\n') #Number of atoms. Align number of atom to the very left
            f.write('  '), #To align the colums of energy and coordiante. Comma to continue on the same line
            if self.dip:#Print energy and dipole
                f.write('{:14.8f}{:14.8f}{:14.8f}{:14.8f}'.format(config[1][0][0],config[1][1][0],config[1][1][1],config[1][1][2])+'\n')
            else:#Print energy only
                f.write('{:14.8f}'.format(config[1][0][0])+'\n')#Print energy only
            for molecule in config[2]:#The atom part: coordiante of atoms
                f.write('{} {:14.8f}{:14.8f}{:14.8f}'.format(molecule[0],molecule[1][0],molecule[1][1],molecule[1][2])+'\n')
        print('{} configs are written to the file: {}'.format(configs_count,filename))
        f.close()

    def sort(self,configs=False,reverse=False,key = 'energy'): #sort based on energy
        """Sort given configurations based on the given key and return the list.

        TODO: sort key = distance
        """

        if configs is False:
            print('Using default input configurations')
            configs = self.configs
        if key == 'energy':

            if reverse: #From high to low
                configs.sort(key= lambda item:item[1][0],reverse = True) #In key, it helps iterate though list.
                return configs
            else: #Default, from low to high.
                configs.sort(key= lambda item:item[1][0])
                return configs
        else:
            return

    def sort_distance(self,configs=False,reverse=False):
        return

    def list(self,first_n_configs=False):
        """Return a list of reading file.

        """


        if first_n_configs == False:
            return self.configs
        else:
            return self.configs[0:int(first_n_configs)]

    def threshold_energy(self,configs,lower=False,upper=False):
        """First sort the list then select the configs based on energy threshold

        """

        print('Excuting energy_threshold:----------------')
        configs = self.sort(configs)
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

    def threshold_distance(self,configs,lower=False,upper=False):
        return

    def error_message(self):
        """To show error message and waiting for decision
        Internal
        """
        decision = raw_input('Do you want to continue? (y/n) \n')
        if decision is not 'y':
            print('Exiting program.')
            exit()

    def distance(self,config,atom1,atom2):
        """To calculate distance with given two atoms.
        Internal
        atom = molecue[i] = [[element],[x, y, z]]
        """

        import numpy as np
        import math
        atom1 = int(atom1)-1
        atom2 = int(atom2)-1
        atom1 = np.array(config[2][atom1][1])
        atom2 = np.array(config[2][atom2][1])
        dis = math.sqrt(np.sum(np.square(atom1-atom2)))
        return dis


    def switch(self):

        new atom 1 is the old atom :
        molecule_count = 1
        if molecule_count <= self.molecule_count_total:
            atom1 = raw_input('New atom {:<2d} is the old atom {:<2d}: ',format(molecule_count))
        molecule_count = molecule_count + 1

    def translate(self):

        atom1 = raw_input('Atom 1: ')
        pass


    def plot(self):
        import numpy as np
        import matplotlib.pyplot as plt
        x = np.arange(0,3 * np.pi, 0.1)
        y = np.sin(x)
        plt.plot(x,y)
        plt.show()

train_x = 'testpoint_v2b_co2h2o.dat'

a = configs(train_x, first_n_configs=1000)
b = a.sort(a.list()[0:3])
for config in b:
    print('{:14.8f}'.format(a.distance(config,3,6)))

a.order()
#a.write('test.write')
#a.plot()


# count = 0
# for config in b:
#     dis = a.distance(config[2][2],config[2][5])
#     if dis > 8:
#         count = count + 1
#         print(dis)
#         a.prt(config)
#
# print(count)




#print('Test finished')
