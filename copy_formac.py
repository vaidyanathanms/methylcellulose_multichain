#Select latest files of a given type and remove the coordinate from
#file name. For ease in MATLAB postprocessing
import numpy
import os
import shutil
import subprocess
import sys
import glob
import os.path
from subprocess import call

free_chains  = [16,32,48,64,80,100,150]#,200]
free_mons    = 30
graft_chains = 64
graft_mons   = 40
tail_mons    = 10
nsalt        = 510
archarr      = [1,2,3,4]

maindir = os.getcwd()
scratchdir = '/scratch.global/vaidya/'
lmpdir = '/home/dorfmank/vsethura/allfiles/files_pework/src_lmp'
   
for ifree in range(len(free_chains)):
    
    print("Free Number of Chains: ", free_chains[ifree])
    workdir1 = scratchdir + 'pework'
    
    if not os.path.isdir(workdir1):
        sys.exit(workdir1, "not found")

    workdir2 = workdir1 + '/n_' + str(free_chains[ifree])
	
    if not os.path.isdir(workdir2):
        sys.exit(workdir2, "not found")
        
    for iarch in range(len(archarr)):

        free_ntot = free_chains[ifree]*free_mons
        graft_ntot = graft_chains*(graft_mons-tail_mons)

        if archarr[iarch] == 1:
            print("Archval: Block_Block")
            dirstr = 'bl_bl'
            fylstr = 'block_block'
        elif archarr[iarch] == 2:
            print("Archval: Block_Alter")
            dirstr = 'bl_al'
            fylstr = 'block_alter'
        elif archarr[iarch] == 3:
            print("Archval: Alter_Block")
            dirstr = 'al_bl'
            fylstr = 'alter_block'
        elif archarr[iarch] == 4:
            print("Archval: Alter_Alter")
            dirstr = 'al_al'
            fylstr = 'alter_alter'
        else:
            print("Unknown Architecture")
            
        workdir3 = workdir2 + '/' + dirstr
            
        if not os.path.isdir(workdir3):
            sys.exit(workdir3,"not found")
                
        os.chdir(workdir3)
        destdir = os.getcwd()

        print("Copying Files")
            
        # * Heart of the Code
        
        resdir = destdir + '/results_' + str(free_chains[ifree]) + '_' \
            + dirstr

        if not os.path.exists(resdir):
            sys.exit(resdir, "not found")

        superdir = workdir2 + '/results_' + str(free_chains[ifree]) + '_' \
            + dirstr

        if not os.path.exists(resdir):
            os.mkdir(superdir)


        f_name = resdir + '/adsfracv2_*'
        list_of_files = glob.glob(f_name) 

        if list_of_files:

            latest_fyl = max(list_of_files, key=os.path.getctime)
            newf_name = resdir + '/adsfracv2.lammpstrj'

            if os.path.isfile(newf_name):
                print(newf_name)
                os.remove(newf_name)
                shutil.copy(latest_fyl, newf_name)
            else:
                shutil.copy(latest_fyl, newf_name)

            finf_name = superdir + '/adsfracv2.lammpstrj'
            if os.path.isfile(finf_name):
                print(finf_name)
                os.remove(finf_name)
                shutil.copy(latest_fyl, finf_name)
            else:
                shutil.copy(latest_fyl, finf_name)

        f_name = resdir + '/dens_*'
        list_of_files = glob.glob(f_name) 

        if list_of_files:

            latest_fyl = max(list_of_files, key=os.path.getctime)
            newf_name = resdir + '/dens.lammpstrj'
            if os.path.isfile(newf_name):
                print(newf_name)
                os.remove(newf_name)
                shutil.copy(latest_fyl, newf_name)
            else:
                shutil.copy(latest_fyl, newf_name)

            finf_name = superdir + '/dens.lammpstrj'
            if os.path.isfile(finf_name):
                print(finf_name)
                os.remove(finf_name)
                shutil.copy(latest_fyl, finf_name)
            else:
                shutil.copy(latest_fyl, finf_name)

        f_name = resdir + '/grpdens_*'
        list_of_files = glob.glob(f_name) 

        if list_of_files:

            latest_fyl = max(list_of_files, key=os.path.getctime)
            newf_name = resdir + '/grpdens.lammpstrj'
            if os.path.isfile(newf_name):
                print(newf_name)
                os.remove(newf_name)
                shutil.copy(latest_fyl, newf_name)
            else:
                shutil.copy(latest_fyl, newf_name)

            finf_name = superdir + '/grpdens.lammpstrj'
            if os.path.isfile(finf_name):
                print(finf_name)
                os.remove(finf_name)
                shutil.copy(latest_fyl, finf_name)
            else:
                shutil.copy(latest_fyl, finf_name)


