#!/usr/bin/env python
#import configs

"""TODO:
        1. Add info method
        2. Give hint on type, consistency checking
        3. Add plot distribution method (urgent)
        4. Add random selection method
        5. In __init__, check the same config
        6. To speed up, use Numpy whenever possible, for example in sorting.
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
            * Skip blank lines automatically.
            * Number of atoms: existence, type and consistency
            * Energy: existence and type
            * Dipole (if dip = True): existence, type and dimension
            * Element: existence, type and consistency
            * Cartisian coordiante: existence, type and dimension
            Notice: consistency error includes the existence error because not existent is also considered as inconsistent.
        '''
        """Constants"""
        self.hartree_to_cm = 219474.63

        self.logo()
        print('Reading file...\n')
        f = open(train_x)
        configs_count = 0
        line_count = 0
        element_1 = list()
        configs = list()
        dipole = list()
        blank_line_count = 0
        for line in f:
            line = line.strip() #The delete the newline character
            if len(line) == 0:
                blank_line_count += 1
                continue
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
                        print("Dimension error: Dipole dimension in line "+str(line_count)+" is : "+str(len(dipole)))
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


        self.blank_line_count = blank_line_count
        self.configs = configs
        self.dip = dip
        self.configs_count = configs_count
        self.train_x = train_x
        self.molecule_count_total = molecule_count_total
        self.line_count = line_count

        self.configs_sorted = self.sort(configs)
        self.energy_lowest=self.energy_array_sorted[0]
        self.energy_highest=self.energy_array_sorted[-1]
        self.energy_lowest_cm = self.energy_array_sorted_cm[0]
        self.energy_highest_cm = self.energy_array_sorted_cm[-1]
        print('Number of blank lines in file:    {:<2d}'.format(self.blank_line_count))
        print('Configuration check finished!')
        #print('**Status: File reading finshed.\n')
        print('Number of configurations:    {:<6d}'.format(self.configs_count))

        print('Reading finished. You can use self.info() or self.help() to start')
        #self.info()
        #self.help()


        f.close()

    def logo(self):
        print("""
                ######        ######  #######  #####
                #     # #   # #     # #       #     #
                #     #  # #  #     # #       #
                ######    #   ######  #####    #####
                #         #   #       #             #
                #         #   #       #       #     #
                #         #   #       #######  #####
                                            Version 0.0.1

                                    --A Bowman Group Product
                                    """
                )
        """--Developed by Kee"""


    def info(self,configs=False):


        print('============================')
        print('||          INFO          ||')
        print('============================')

        print('Inputfile:  {}'.format(self.train_x))
        print('Number of atoms:    {:<2d}'.format(self.molecule_count_total))
        self.order()
        print('Number of configurations:    {:<6d}'.format(self.configs_count))

        print('\nLowest energy (Hartree):    {:14.8f}'.format(self.energy_lowest))
        print('Lowest energy (wavenumber):    {:14.8f}'.format(self.energy_lowest_cm))
        self.prt(self.configs_sorted[0])
        print('Highest energy (Hartree):   {:14.8f}'.format(self.energy_highest))
        print('Highest energy (wavenumber):   {:14.8f}'.format(self.energy_highest_cm))
        print('Configuration that contains highest energy: ')
        self.prt(self.configs_sorted[-1])
        print('========END OF INFO=======\n')

    def order(self,monomer = False, configs=False):
        """TODO: Add feature to show the group of monomers"""


        configs=self.configs_check(configs)
        print('Atom numbering:')
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
                *Print the configs to the screen with pretty format

            a.write('output.xyz', configs=False)
                *Write the list of configs into output file with pretty format

            a.plot(self,configs = False,binwidth=False)
                *Plot energy distribution with gien binwith (default is 50 cm-1)

            a.sort(configs=False,reverse=False,key ='energy')
                *Sort the configs according to energy and return the list

            a.switch()
                *Switch the order to atoms of all configs as return as list.

            a.list(first_n_configs=False)
                *Return default configurations as a list

            a.threshold_energy(configs,lower=False,upper=False)
                *Select all configurations that satisfy the energy threshold

            a.distance(config,atom_A,atom_B)
                *With given one configuraiton, calculate the distance of two atoms

            a.order()
                *Show the order of atoms

            a.add_expand(self, n_configs=False,configs=False, first_n_configs=False, monomer=False, dis_lower=False, dis_uppder=False, dis_new_lower=False, dis_new_upper=False)
                *Select configs randomly in the given region[dis_lower, dis_upper], then separate them along given two atoms to the given new distance [dis_new_lower, dis_new_upper].


        *To show info:
            a.info()

        *To show syntax help:
            a.help()

        =======END OF HELP SYNTAX====

        """
        #for line in syntax.split('\n'):
        #            print(line.strip())
        print(syntax)

    def prt(self,configs=False):
        """Print configs into screen.

        """

        #print('*Status: printing...')
        configs_count=0
        configs=self.configs_check(configs)

        for config in configs: #Write config by config
            configs_count = configs_count + 1
            print('    {:<2d}'.format(config[0][0])) #Number of atoms. Align number of atom to the very left
            print(' '),# To align the colums of energy and coordiante
            if self.dip:#Print option for having dipole or not
                print('    {:14.8f}{:14.8f}{:14.8f}{:14.8f}'.format(config[1][0][0],config[1][1][0],config[1][1][1],config[1][1][2]))
            else:#Print also dipole (if input file has it)
                print('    {:14.8f}'.format(config[1][0][0]))#Print energy only
            for molecule in config[2]:#The atom part: coordiante of atoms
                print('    {} {:14.8f}{:14.8f}{:14.8f}'.format(molecule[0],molecule[1][0],molecule[1][1],molecule[1][2]))

        print('----Printed {} configuration(s).\n'.format(configs_count))
        #print('*Status: printing finished.\n')

    def write(self,filename, configs=False): #Default value has to be after non-default value.
        """Write formated configurations into filename.

        """

        print('*Status: writting file...\n')
        f=open(filename,'w')
        configs_count=0

        configs=self.configs_check(configs)
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
        print('----{} configs are written to the file: {}\n'.format(configs_count,filename))
        print('*Status: writting finished.\n')
        f.close()

    def sort(self,configs=False,reverse=False,key = 'energy'): #sort based on energy
        """Sort given configurations based on the given key and return the list.

        TODO: sort key = distance
        """
        import numpy as np
        import copy
        print('*Status: sorting...\n')
        configs=self.configs_check(configs)

        if key == 'energy':
            print('\n----Sorting by energy...\n')
            if reverse: #From high to low
                configs.sort(key= lambda item:item[1][0],reverse = True) #In key, it helps iterate though list.
                #print('\nReverse sort finished.\n')
                #return configs
            else: #Default, from low to high.
                configs.sort(key= lambda item:item[1][0])
                energy_array_sorted = list()
                energy_array = list()

                for config in configs:
                    #print(config[1][0][0])
                    energy_array.append(config[1][0][0])
                self.energy_array_sorted = energy_array
                energy_array_sorted_cm = copy.deepcopy(self.energy_array_sorted)
                self.energy_array_sorted_cm = np.array(energy_array_sorted_cm)* self.hartree_to_cm
                #print('\nSort finished.\n')
                #return configs
        else:#Reserved for other key expansion in the future.

            return

        print('\n*Status: sorting finished.\n')
        return configs

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

    def distance(self,config,atom_A,atom_B):
        """To calculate distance with given two atoms.
        Internal
        atom = molecue[i] = [[element],[x, y, z]]
        """

        import numpy as np
        import math
        #print('This is atom_B:',atom_B)
        atom_A = int(atom_A)-1
        atom_B = int(atom_B)-1
        atom_A = np.array(config[2][atom_A][1])
        atom_B = np.array(config[2][atom_B][1])
        dis = math.sqrt(np.sum(np.square(atom_A-atom_B)))
        return dis

    def switch(self,configs=False):
        import copy
        print('Switching order:\nThe origianl numbering of given configs is:')
        self.order(configs)
        molecule_count = 1
        molecule_new_count = 1


        if configs is False:#Default input check module
            configs = self.configs
            print('Switching using orignial configuraitons')

        print('\nWhat new order do you what?\n')

        try:
            temp = configs[0][0][0]
        except:#If has to configuration, this makes sure it can be safely iterated in the next statemnt.
            configs = [configs]

        configs_new = copy.deepcopy(configs)#This is the correct way to create a different list with same value

        while molecule_new_count <= self.molecule_count_total:#Will repeat molecule_count_total times
            configs_new_count = 0

            try:
                molecule_count=int(raw_input('I want new atom ({:d}) to be the old atom number: '.format(molecule_count)))
            except:
                pass
            #molecule_count = self.molecule_count_total - molecule_new_count + 1 #Test arguement(reverse order)
            for config_new in configs_new:
                configs_new[configs_new_count][2][molecule_new_count-1] =configs[configs_new_count][2][molecule_count-1]
                configs_new_count = configs_new_count + 1
            molecule_new_count = molecule_new_count + 1
            #        self.prt(configs_new)
        print('New numbering:')
        self.order(configs_new)
        print('''New configs are returned as list, please use a.write(configs_new,'ouput') to save configs.''')
        return configs_new

    def configs_check(self,configs,silence=True):
        """To check if give the configs arugment.

            Usage: configs=configs_check(configs)
        """

        if configs is False:
            if silence is not True:
                print('Using original list.')
            return self.configs_sorted
        else:
            try:
                check = configs[0][0][0]#Check if the list has only one configuration but one layer smaller
                return configs
            except:#If has to configuration, this makes sure it can be safely iterated in the next statemnt.
                return [configs]

    def translate(self,atom_A=False, atom_B=False, dis_new = False, config = False):
        """Give a dimer and designated atom groups, translate in vector AB from previous distance to new distance

            Input:
                dis_new :: the distance after translation, essiential
                monomer :: contains the info about monomer groups and referenced atoms, omit this if use default.
                config :: take one configuraiton so that it can do the translation. If given multiple configuraitons, return only first one.

            Output:
                configs_new :: none-nested configuraiton, subset of configs.

            Theory:
                For point A and B. Now translate B along AB direction to distance dis_new = |AC| and become C. Parameter function for line AB:
                    xC = xB + m*t
                    yC = yB + n*t
                    zC = zB + p*t
                where:
                    t = (dis_new-dis)/sqrt((m**2+n**2+p**2))
                    (m,n,p) = vector(AB)

                Suppose D is also the same monomer as B, and in ordor to translate D:
                    xD_new = xD + m*t
                    yD_new = yD + n*t
                    zD_new = zD + p*t
                where t and m,n,p take the same value as previous one.


        """
        import numpy as np
        import copy
        #print("-------------Atom_A",atom_A)
        configs = self.configs_check(config)
        #print(type(self.monomers()))
        monomers = self.monomers()
        #print('test')
        #print(monomers)

        monomer_A = monomers[0]
        monomer_B = monomers[1]
        #print ('This is monomer_B:',monomer_B)
        if atom_A is False: #Assign atom a and atom b
            try:
                trail = raw_input()
                self.order()
                atom_A = int(raw_input('Please assign the number of first atom: '))
                atom_B = int(raw_input('Please assign the number of second atom: '))
            except:
                print('Translate using default Atom (first atom in monomer_A and first atom in monomer_B)')
                #print(monomer_A)
                atom_A = monomer_A[0]
                #print ('This is atom_A:',atom_A)
                #atom_B = monomer_B[0]
                print ('This is atom_B:',atom_B)
                dis_new = 10



        config_new = copy.deepcopy(configs[0]) #Everytime need a new list with same vaule need deepcopy.
        #print(monomer_A)
        #print('New atom_A: ',atom_A)
        #print('New atom_B: ',atom_B)

        dis = self.distance(config_new,atom_A,atom_B)
        v_AB = np.array(self.vector(config_new, atom_A, atom_B))
        v_A = np.array(self.vector(config_new, atom_A))
        v_B = np.array(self.vector(config_new, atom_B))
        t = (dis_new-dis)/np.sqrt((np.sum(np.square(v_AB))))

        monomer_B_new = list()
        for atom in monomer_B:
            v_atom = np.array(self.vector(config_new, atom))
            v_atom_new = v_atom + v_AB * t
            config_new[2][atom-1][1] = v_atom_new
        return config_new # One none-nested configuraiton list

    def monomers(self, monomer_A = False, monomer_B = False):
        """To assign monomer groups and return as two lists of number. Use default groups if no assign.

            Example: monomer = [[[3, 4, 5], [6, 1, 2]]]
        """

        monomers=list()

        try:
            if len(self.monomer_AB) is not 0:
                monomers = self.monomer_AB
                 #Meaning monomer has already been assigned. Once assigned forever assigned.
        except:
            print('Monomer is not determined.')
            #if monomer_A is False: #For the convinience of other function's usage
            try:
                a = raw_input()#To try out if can catch input
                self.order('Hit Enter to continue:')
                print('(Enter integers and separate them by whitespace)')
                str1 = raw_input('What atoms are in first monomer: ')
                str2 = raw_input('Waht atoms are in second monomer: ')
            except:
                print('Warming: Using test arguments. Please use terminal to catch input.')
                atom_A= 3# Arguemnt for test
                atom_B= 6# Argument for test
                str1 = '3 4 5'# Argument for test
                str2 = '6 2 1'# Argument for test
            monomer_A = list()
            monomer_B = list()
            for line in str1.split():
                monomer_A.append(int(line.strip()))
            for line in str2.split():
                monomer_B.append(int(line.strip()))
            monomers=[monomer_A, monomer_B]
            self.monomer_AB = monomers
        #print monomers

        #monomer_A = [atom_A]
        #monomer_B = [atom_B]

        return monomers
        #Example: monomers = [[[3, 4, 5], [6, 1, 2]]]

    def dissociation(self,dis_new_lower=False,dis_new_upper=False,atomA=False,atomB=False,config=False):


            config = self.configs_sorted[0]
            self.prt(config)
            self.molden(config)
            print('----Translate along the A-B direction')
            #atom_A = int(raw_input('What is the number of atom A: '))
            #atom_B = int(raw_input('What is the number of atom B: '))
            print('(Input integer with whitespace)')
            #str1 = raw_input('The monomer that contains A also contains: ')
            #str2 = raw_input('The monomer that contains B also contains: ')

            #self.order(configs)
            #for dis in range(0,10,0.5):
        #        self.translate(configs=False,atom_A=3, atom_B = 6, monomer_A = '4 5', monomer_B = '1 2', dis_lower = 2, dis_upper = 10, dis_new_lower=dis, dis_new_upper=dis)

    def add_expand(self, n_configs=False,configs=False,  first_n_configs=False, monomer=False, dis_lower=False, dis_uppder=False, dis_new_lower=False, dis_new_upper=False):
        import numpy as np
        configs_new = list()
        configs = self.configs_check(configs)

        if dis_lower is not False:
            try:
                a = raw_input()
                dis_lower = float(raw_input('The original distance_min (Angstrom) you want is : '))
                dis_upper = float(raw_input('The original distance_max (Angstrom) you want is: '))
                dis_new_lower = float(raw_input('The original distance_new_min (Angstrom) you want is: '))
                dis_new_lower = float(raw_input('The original distance_new_man (Angstrom) you want is: '))
            except:
                print('Using default dis boundaries')
                dis_new_lower = float(6)
                dis_upper = float(5)
                dis_lower = float(2)
                dis_new_upper = float(9)

        monomers = self.monomers(configs)
        print(monomers[0])
        atom_A = monomers[0][0]
        atom_B = monomers[1][0]

        if n_configs is not False:
            n_configs = int(n_configs)
        else:
            n_configs = len(configs)

        configs_count = 0
        configs_ok_count = 0
        rand_pool = len(configs) + 1
        while configs_ok_count < n_configs:
            configs_count += 1
            #print(np.random.randint(1,rand_pool))
            config = configs[np.random.randint(1,rand_pool)-1]
            #self.prt(config)
            dis = self.distance(config,atom_A,atom_B)
            if dis_lower <= dis <= dis_upper:
                configs_ok_count += 1
                dis_new = np.random.uniform(dis_new_lower,dis_new_upper)
                config_new = self.translate(atom_A,atom_B,dis_new,config)
                configs_new.append(config_new)
                #print('{:5.2f} {:5.2f} {:5.2f} {:5.2f}'.format(self.distance(config,atom_A,atom_B),self.distance(config_new,atom_A,atom_B),dis_new,self.distance(config_new,atom_A,atom_B)-self.distance(config,atom_A,atom_B)))
                print('New distance: {:5.2f}'.format(self.distance(config_new,atom_A,atom_B)))
        print('{:d}/{:d} configurations are returned as list.'.format(len(configs_new),configs_ok_count))
        print('*Add point: Expand finished.')
        return configs_new #List of configs that have just been expanded.

    def vector(self,config,atom_A,atom_B=False):
        """For a given configuration and two atoms A and B, give back the vector BA.

        """
        import numpy as np

        if atom_B is False:
            atom_A = int(atom_A) - 1
            return np.array(config[2][atom_A][1])
        else:
            atom_A = int(atom_A) - 1
            atom_B = int(atom_B) - 1
            atom_A = np.array(config[2][atom_A][1])
            atom_B = np.array(config[2][atom_B][1])
            vector21 = atom_B-atom_A
            return vector21

    def molden(self,configs=False):
        configs = self.configs_check(configs)
        self.write('plot.temp',configs)
        #self.cl('/raid/molden4.7/molden plot.temp')#Please use the absolut path of molden becasue it is not considerd as 'installed'
        self.cl('/Users/qingfengwang/Desktop/molden/molden plot.temp')

    def cl(self,command):
        #ip::string, command line as string input
        #op::string, return value is the output of command line
        #Notice, each time when change dire.ctly, cl starts from currect directory.
        #Use three \' if you want to input multiple line
        import subprocess
        import os
        import shlex
        arg = shlex.split(command)
        p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        (output, err) = p.communicate()
        print(output)
        return output

    def plot(self,configs = False,binwidth=False):
        import numpy as np
        import matplotlib.pyplot as plt

        configs = self.configs_check(configs)

        if configs is not False:
            energy_array = list()
            for config in configs:
                energy_array.append(config[1][0][0]* self.hartree_to_cm)
            #print(energy_array)
        else:
            energy_array = self.energy_array_cm
        fig = plt.figure()
        ax = fig.add_subplot(111)

        x = energy_array
        #x = np.random.normal(0,1,1000)
        #print(x)

        if binwidth is False:
            binwidth = 50
        else:
            binwidth = int(binwidth)

        print('\nCurrent binwidth is {:d} cm-1\n'.format(binwidth))


        numBins = (self.energy_highest_cm - self.energy_lowest_cm) // binwidth
        #print(numBins)
        #numBins = 100
        ax.hist(x,numBins,color='green',alpha=0.8)
        ax.set_xlabel("Energy(cm$^{-1}$)")
        ax.set_ylabel("Frequency")
        ax.annotate('Lowest Energy',xy=(energy_array[0],200),xytext=(energy_array[0], 100),arrowprops=dict(facecolor='black',shrink =0.005))
        ax.annotate('Highest Energy',xy=(energy_array[-1],200),xytext=(energy_array[-1], 100),arrowprops=dict(facecolor='black',shrink =0.005))
        fig = plt.gcf()
        plt.show()
        decision = raw_input('''Do you want to save the file? (Enter 'y' to save, enter others to skip)''')
        if decision is 'y':
            filename = raw_input('Please specify .eps (1200 dpi) filename: ').strip()

            fig.savefig(filename, format='eps', dpi=1200)
            print('Plot saved to {}.'.format(filename))

        #print(energy_array[0])

    def v2b(self,configs=False):
        """It returns configs of two monomers"""

        configs = self.configs_check(configs)
        monomers = self.monomers(configs)
        monomer_A = monomers[0]

        configs_A = list()
        molecule_count_A = len(monomer_A)

        for config in configs:
            coordinate=list()
            for atom in monomer_A:
                coordinate.append(config[2][atom-1])
            configs_A.append([[molecule_count_A],[[0]],coordinate])


        monomer_B = monomers[1]
        configs_B = list()
        molecule_count_B = len(monomer_B)
        for config in configs:
            coordinate=list()
            for atom in monomer_B:
                coordinate.append(config[2][atom-1])
            configs_B.append([[molecule_count_B],[[0]],coordinate])

        self.prt(configs_A)
        self.prt(configs_B)


"""Test arguemnts"""





#dis_new_lower = float(6)
#dis_upper = float(5)
#dis_lower = float(2)
#dis_new_upper = float(9)

#train_x = 'testpoint_v2b_co2h2o.dat'

#train_x = 'dimer_47358.abE'

#a = configs(train_x,first_n_configs=2000)
#a = configs(train_x)
#a.plot()
#b = a.list()
#a.prt(b)
#a.dissociation()

#c = a.translate(config= b,dis = 10)
#a.prt(c)
#print(a.distance(c,3,6))
#a.order()
#a.monomers()
#print(a.monomer_AB)
#a.monomers()
#a.add_expand(1000,dis_lower=dis_lower,dis_uppder=dis_upper,dis_new_lower=dis_new_lower,dis_new_upper=dis_new_upper)

#new = a.v2b()
#a.switch(b)
#a.write('configsA',new)
#a.prt(new[0])
#a.plot()
#a.prt(b)

#a.prt(b)
#c = a.add_expand(2000,configs=b)
#a.prt(c)
#a.monomer()
#a.write('6-9A_2000.xyz',a.sort(c))
