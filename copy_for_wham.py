#Select files from COLVAR output and feed for WHAM
import numpy
import os
import shutil
import subprocess
import sys
import glob
import os.path
from subprocess import call
import re

correl_time = 1
skipfr  = 100
thresholdfr = 10000
kspr = 20.0
prefix  = 'US_dir2_' # Keep the prefix without number
colfile = 'out.colvars.traj'
colfile_name = 'umbsampling'
maindir = os.getcwd()
scratchdir = '/scratch.global/vaidya/'
rvalarr = []

whamexec = '/home/dorfmank/vsethura/wham/wham/wham'
whamcopy = maindir + '/wham'

if not os.path.isfile(whamcopy):
    shutil.copy2(whamexec,whamcopy)

whammain = 'whaminp.txt'
dir_names = maindir + '/' + prefix + '*'
list_of_dirs = glob.glob(dir_names) 

if list_of_dirs:
    
    wham_dir = maindir + '/wham_ana'
    fmain = open(whammain,'w')
    if not os.path.exists(wham_dir):
        os.mkdir(wham_dir)

    for dirval in range(len(list_of_dirs)):
            
        splice  = list_of_dirs[dirval].split('_')
        rval = splice[len(splice)-1]
        rvalarr.append(float(rval))
        
    for dirval in sorted(rvalarr):

        dirname = maindir + '/' + prefix + str(dirval)
        os.chdir(dirname)
        fylname = dirname + '/' + colfile
        print(dirname)

        if os.path.isfile(fylname):
            
            wham_fyl = wham_dir + '/' + 'whaminp_' + str(dirval)
            fr = open(fylname)
            fw = open(wham_fyl,'w')
            
            line = fr.readline()
            header = ' '.join(line.split()).split(' ')

            for nvar in range(len(header)):

                if header[nvar] == colfile_name:

                    rvar = nvar-1
                    
                elif header[nvar] == 'step':

                    stepvar = nvar-1 

            
            fr.close()
            
            with open(fylname,'r') as fyl:

                for sval in range(skipfr):

                    next(fyl)

                linecntr = 0
                for line in fyl:

                    data = ' '.join(line.split()).split(' ')
                    
                    if not data[rvar].isalpha():

                        fw.write(str(data[stepvar])+'\t'+str(data[rvar])+'\n')
                        linecntr = linecntr + 1

            fw.close()
            
            if linecntr < thresholdfr:
                print("Number of lines do not meet thresold criteria"\
                      ,linecntr,"in",wham_fyl)
            else:
                fmain.write(wham_fyl + '\t' +str(dirval) +'\t'+str(kspr)+'\t'+str(correl_time)+'\n')

            
