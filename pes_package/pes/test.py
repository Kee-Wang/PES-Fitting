def cl(command):
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

cl('which molden')
#cl('/raid/molden4.7/molden plot.temp')
