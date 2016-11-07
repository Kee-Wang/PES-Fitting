#!/usr/bin/env python
import subprocess
import os
import shlex
import numpy as np
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
train_x = 'pts.dat'
f = open(train_x)
nol = 0
for line in f:
  nol=nol+1
  if nol==1:
    natom=int(line)
nconfig = nol/(natom+2)
f.close()


