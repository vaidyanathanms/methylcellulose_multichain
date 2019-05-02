import numpy
import os
import shutil
import subprocess
import sys
import glob

from subprocess import call

#------------------Input Directory Details---------------------------

trialnum = 7
temp     = 50
free_chains  = [8]
py_order = 2 #1-Alternate,2-Block,3-Nucleated,0-random,4-jammed
py_dist = ['0.03']#,'0.03','0.07']#,'0.05','0.06','0.09'] 

#-----------------Input variables for analysis-----------------------

nframes     = 4000
skipfr      = 0
nwater      = 0
nmons       = 1000
rgfreq      = 1
rgsysfreq   = 1
eigflag     = 1
eigindflag  = 1
globeigflag = 1
#-----------------Global Directories---------------------------------

maindir = os.getcwd()
scratchdir = '/scratch.global/vaidya/'
jobdir = '/home/dorfmank/vsethura/allfiles/files_cgmc/src_lmp'
lmpdir = '/home/dorfmank/vsethura/mylammps/src'

#-----------------Main Analysis--------------------------------------

for dval in range(len(py_dist)):  

    workdir1 = scratchdir + 'cg_mccalc'
    
    if not os.path.isdir(workdir1):
        os.mkdir(workdir1)

    workdir1 = workdir1 + '/temp_' + str(temp)
    if not os.path.isdir(workdir1):
        print("Error:", workdir1,"not found")
        continue

    workdir1 = workdir1 + '/diffconfig' 
    if not os.path.isdir(workdir1):
        print("Error:", workdir1,"not found")
        continue

    workdird2 = workdir1 + '/n_' + str(nmons)
    if not os.path.isdir(workdird2):
        print("Error:", workdird2,"not found")
        continue

    workdirn2 = workdird2 + '/trial' + str(trialnum)
    if not os.path.isdir(workdirn2):
        print("Error:", workdirn2,"not found")
        continue

    for chainval in range(len(free_chains)):

        workdir2 = workdirn2 + '/n' + str(free_chains[chainval])
        if not os.path.isdir(workdir2):
            print("Error:", workdir2,"not found")
            continue

        if py_order == 2:

            arch = 'block'

        elif py_order == 1:

            arch = 'alternate'

        elif py_order == 0:

            arch = 'random'

        elif py_order == 3:

            arch = 'nucleate'

        elif py_order == 4:

            arch = 'jammed'
    
        print("Analyzing TrialNum/distval/nchains,Order: ",\
              trialnum,py_dist[dval],free_chains[chainval],arch)

        workdir3 = workdir2 + '/' + arch
        if not os.path.isdir(workdir3):
            print(workdir3, "not found")
            continue

        workdir4 = workdir3 + '/dist_' + str(py_dist[dval])
        if not os.path.isdir(workdir4):
            print(workdir4, "not found")
            continue

	
        os.chdir(workdir4)
        destdir = os.getcwd()

        print( "Copying Files")
            
        # Copy Files
        srcfyl = maindir + '/main.f90'
        desfyl = destdir + '/main.f90'
        shutil.copy2(srcfyl, desfyl)
        
        srcfyl = maindir + '/params.f90'
        desfyl = destdir + '/params.f90'
        shutil.copy2(srcfyl, desfyl)

        srcfyl = maindir + '/ran_numbers.f90'
        desfyl = destdir + '/ran_numbers.f90'
        shutil.copy2(srcfyl, desfyl)
        
        srcfyl = maindir + '/anainp_var.txt'
        desfyl = destdir + '/anainp_var.txt'
        shutil.copy2(srcfyl, desfyl)

        srcfyl = lmpdir  + '/lmp_mesabi'
        desfyl = destdir + '/lmp_mesabi'
        shutil.copy2(srcfyl, desfyl)


        srcfyl = jobdir  + '/jobana.sh'
        desfyl = destdir + '/jobana.sh'
        shutil.copy2(srcfyl, desfyl)


        # Manipulate Files
            
        restart_list = glob.glob("restart*")
        if not restart_list:
            print("ERROR: No restart files found in", destdir)
            continue

        latest_rest = max(restart_list, key=os.path.getctime)

        dataname = "replicated_"+str(free_chains[chainval])+"_"\
                   +str(py_dist[dval])+".data"
        subprocess.call(["mpirun","-np","24","./lmp_mesabi","-r",\
                         latest_rest,dataname])

        if os.path.exists(dataname):
            

            print( "Dataname for LAMMPS: ", dataname)
		
        else:
		
            print( "ERROR:", dataname, "not found in", destdir)
            continue


        #First create jobana.sh and keep adding executables
        launch_fyl = 'jobana2.sh'
        fr = open('jobana.sh','r')
        fw = open(launch_fyl,'w')
        fid = fr.read().replace("pydist",py_dist[dval]).\
              replace("pynchains",str(free_chains[chainval])).\
              replace("pyarch",arch)
        fw.write(fid)
        fw.close()
        fr.close()

        fana = open(launch_fyl,'a')
        fana.write("\n")
        # * means all if need specific format then *.csv
        trajnames = destdir + '/dump*.lammpstrj'
        list_of_files = glob.glob(trajnames) 

        if not list_of_files:
            print("Warning, No files found in ", destdir)
            continue

        for trajcnt in range(len(list_of_files)):

#            traj_str = max(list_of_files, key=os.path.getctime)
            traj_str = list_of_files[trajcnt]
            traj_arr = traj_str.split("/")

            print( "trajval", traj_arr[len(traj_arr)-1])
            
            analaunch_fyl = 'anainp' + str(trajcnt+1) + '.txt'
            fr = open('anainp_var.txt','r')
            fw = open(analaunch_fyl,'w')
            fid = fr.read().replace("py_datafyle",dataname).\
                  replace("py_trajfyle",traj_arr[len(traj_arr)-1]).\
                  replace("py_nframes",str(nframes)).\
                  replace("py_skipfr",str(skipfr)).\
                  replace("py_nchains",str(free_chains[chainval])).\
                  replace("py_water",str(nwater)).\
                  replace("py_atperchain",str(nmons)).\
                  replace("py_rgfreq",str(rgfreq)).\
                  replace("py_rgsysfreq",str(rgsysfreq)).\
                  replace("py_eigflag",str(eigflag)).\
                  replace("py_indeigflag",str(eigindflag)).\
                  replace("py_globeigflag",str(globeigflag))            
            fw.write(fid)
            fw.close()
            fr.close()


            
            subprocess.call(["ifort","-r8","-mkl","-qopenmp","ran_numbers.f90",
                             "params.f90","main.f90","-o","ana.o"])

            lyne = "./ana.o " + analaunch_fyl+"\n"
            fana.write(lyne)
            fana.write("wait\n")


        fana.write("\n")
        lynedir = "outresults_"+py_dist[dval]+"_"\
                  +str(free_chains[chainval])+ "_" + arch + "\n"
        fana.write("mkdir " + lynedir)

        lyne = "mv " + "eig* " + lynedir
        fana.write(lyne)
        lyne = "mv " + "avgeig* " + lynedir
        fana.write(lyne)
        lyne = "mv " + "indeig* " + lynedir
        fana.write(lyne)
        lyne = "mv " + "rg* " + lynedir
        fana.write(lyne)
        lyne = "mv " + "glob* " + lynedir
        fana.write(lyne)
        fana.close()

        subprocess.call(["qsub","jobana2.sh"])
        os.chdir(maindir)
