import numpy
import os
import shutil
import subprocess
import sys
import glob

from subprocess import call

restart = 0 # For restarting from given configurations
trialnum = 21
trial_ax = 1 #Different from trial num.. Depends upon how many input configs are available
temp    = 50
nmons    = 1000
free_chains  = [4]#,8]
py_order     = 2 #2-Alternate, 1-Block, 3-Nucleated site, 0-random
py_dist = ['0.03']#,'0.05','0.07']##,'0.07','0.09','0.12','0.14','0.16','0.18','0.20','0.24','0.28'] #0.18','0.21','0.24']#,'0.05'] 
percflex = 0.5

if trial_ax == 1:
    majoraxisval = ['0.615','0.784','0.080']# for trial1
elif trial_ax == 2:
    majoraxisval = ['0.196','0.644','0.739']# for trial2
elif trial_ax == 3:
    majoraxisval = ['0.115','0.079','0.990']# for trial3

datafyl = "inpdata_trial" + str(trial_ax)+ ".data"
maindir = os.getcwd()
scratchdir = '/scratch.global/vaidya/'
lmpdir = '/home/dorfmank/vsethura/allfiles/files_cgmc/src_lmp'
   

inpcheck = maindir + '/' + datafyl
if not os.path.isfile(inpcheck):
    print(datafyl, "not present")
    sys.exit()

for dval in range(len(py_dist)):
    
    workdir1 = scratchdir + 'cg_mccalc'
    
    if not os.path.isdir(workdir1):
        os.mkdir(workdir1)

    workdir1 = workdir1 + '/temp_' + str(temp)
    if not os.path.isdir(workdir1):
        os.mkdir(workdir1)

    workdir1 = workdir1 + '/diffconfig' 
    if not os.path.isdir(workdir1):
        os.mkdir(workdir1)


    workdird2 = workdir1 + '/n_' + str(nmons)
    if not os.path.isdir(workdird2):
        os.mkdir(workdird2)

    workdirn2 = workdird2 + '/trial' + str(trialnum)
    if not os.path.isdir(workdirn2):
        os.mkdir(workdirn2)

    for chainval in range(len(free_chains)):

        workdir2 = workdirn2 + '/n' + str(free_chains[chainval])
        if not os.path.isdir(workdir2):
            os.mkdir(workdir2)

        if py_order == 2:

            arch = 'alternate'

        elif py_order == 1:

            arch = 'block'

        elif py_order == 4:
            
            arch = 'jammed'

        elif py_order == 0:

            arch = 'random'

        elif py_order == 3:

            arch = 'nucleate'

        if py_order == 3:
            py_percent = 1.0
        elif py_order == 4 or percflex == 0.0:
            py_percent = 0.0
            print("In Jammed/Parallel Configuration")
        else:
            py_percent = (free_chains[chainval]*percflex+1)/(free_chains[chainval])

        workdir3 = workdir2 + '/' + arch
        if not os.path.isdir(workdir3):
            os.mkdir(workdir3)

        workdir4 = workdir3 + '/dist_' + str(py_dist[dval])
        if not os.path.isdir(workdir4):
            os.mkdir(workdir4)

        os.chdir(workdir4)
        destdir = os.getcwd()

        if restart == 0:

            print( "Starting from beginning")

            print( "Copying Files")

            # Copy Files
            srcfyl = maindir + '/replicate_diff_geometry.f90'
            desfyl = destdir + '/replicate_diff_geometry.f90'
            shutil.copy2(srcfyl, desfyl)

            srcfyl = maindir + '/ran_numbers.f90'
            desfyl = destdir + '/ran_numbers.f90'
            shutil.copy2(srcfyl, desfyl)

            srcfyl = maindir + '/' + datafyl
            desfyl = destdir + '/' + datafyl
            shutil.copy2(srcfyl, desfyl)
            srcfyl = maindir + '/inprep_var.txt'
            desfyl = destdir + '/inprep_var.txt'
            shutil.copy2(srcfyl, desfyl)

            lmpdataname = "replicated_" + datafyl + "_" + py_dist[dval] 

            srcfyl = lmpdir  + '/in.cgrun1_var'
            desfyl = destdir + '/in.cgrun1_var'
            shutil.copy2(srcfyl, desfyl)
            srcfyl = lmpdir + '/in.cgrun2'
            desfyl = destdir + '/in.cgrun2'
            shutil.copy2(srcfyl, desfyl)
            srcfyl = lmpdir + '/in.cgrun3'
            desfyl = destdir + '/in.cgrun3'
            shutil.copy2(srcfyl, desfyl)
            srcfyl = lmpdir + '/in.cgrun4'
            desfyl = destdir + '/in.cgrun4'
            shutil.copy2(srcfyl, desfyl)
            srcfyl = lmpdir + '/in.cgpair_MC'
            desfyl = destdir + '/in.cgpair_MC'
            shutil.copy2(srcfyl, desfyl)

            if free_chains[chainval] != 2:
                srcfyl = lmpdir + '/jobmain.sh'
                desfyl = destdir + '/jobmain.sh'
                shutil.copy2(srcfyl, desfyl)
            else:
                srcfyl = lmpdir + '/jobmain_2chains.sh'
                desfyl = destdir + '/jobmain.sh'
                shutil.copy2(srcfyl, desfyl)
                    
            srcfyl = '/home/dorfmank/vsethura/mylammps/src/lmp_mesabi'
            desfyl = destdir + '/lmp_mesabi'
            shutil.copy2(srcfyl, desfyl)


            # Manipulate Files

            print( "Copy Successful - Manipulating Input Files")

            launch_fyl = 'in.cgrun1'
            fr = open('in.cgrun1_var','r')
            fw = open(launch_fyl,'w')
            fid = fr.read().replace("lmpdatafyl",str(lmpdataname))
            fw.write(fid)
            fw.close()
            fr.close()


            launch_fyl = 'inprep.txt'
            fr = open('inprep_var.txt','r')
            fw = open(launch_fyl,'w')
            fid = fr.read().replace("py_nchains",str(free_chains[chainval]-1)).\
                  replace("py_xval",str(majoraxisval[0])).\
                  replace("py_yval",str(majoraxisval[1])).\
                  replace("py_zval",str(majoraxisval[2])).\
                  replace("py_dist",str(py_dist[dval])).\
                  replace("py_perc",str(py_percent)).\
                  replace("py_inpdata",str(datafyl)).\
                  replace("py_order",str(py_order))
            fw.write(fid)
            fw.close()
            fr.close()

            # Generate Required Files for LAMMPS

            print( "Running FORTRAN script for generating Datafile")
            subprocess.call(["ifort","-r8","-qopenmp","-check",
                             "-traceback","ran_numbers.f90",
                             "replicate_diff_geometry.f90","-o","lmpinp.o"])

            subprocess.call(["./lmpinp.o","inprep.txt"])


            # Removing unnecessary files


            files = glob.glob(destdir +'/*var*')
            for f in files:
                os.remove(f)
                files = glob.glob(destdir +'/*.mod')
            for f in files:
                os.remove(f)
                files = glob.glob(destdir +'/*.o')
            for f in files:
                os.remove(f)

            print( "Submitting Jobs..")

            subprocess.call(["qsub","jobmain.sh"])

            os.chdir(maindir)

        elif restart == 1:

            print ("Restarting jobs ..")

            # Find if restart file is present

            rest2files = destdir +'/restart3_melt2'
            rest1files = destdir +'/restart3_melt1'

            if os.path.isfile(rest2files):

                srcfyl = lmpdir + '/in.cgrun5'
                desfyl = destdir + '/in.cgrun5'
                shutil.copy2(srcfyl, desfyl)
                srcfyl = lmpdir + '/in.cgpair_MC'
                desfyl = destdir + '/in.cgpair_MC'
                shutil.copy2(srcfyl, desfyl)
                srcfyl = lmpdir + '/jobre.sh'
                desfyl = destdir + '/jobre.sh'
                shutil.copy2(srcfyl, desfyl)

                print( "Submitting Jobs..")

                subprocess.call(["qsub","jobre.sh"])

            elif os.path.isfile(rest1files):

                srcfyl = rest1files
                desfyl = rest2files
                shutil.copy2(srcfyl, desfyl)

                srcfyl = lmpdir + '/in.cgrun5'
                desfyl = destdir + '/in.cgrun5'
                shutil.copy2(srcfyl, desfyl)
                srcfyl = lmpdir + '/in.cgpair_MC'
                desfyl = destdir + '/in.cgpair_MC'
                shutil.copy2(srcfyl, desfyl)
                srcfyl = lmpdir + '/jobre.sh'
                desfyl = destdir + '/jobre.sh'
                shutil.copy2(srcfyl, desfyl)

                print( "Submitting Jobs..")

                subprocess.call(["qsub","jobre.sh"])


            else:

                print("restart3_melt* not found in", destdir)
                print("switching to failsafe mode")

                srcfyl = lmpdir + '/in.cgrun3'
                desfyl = destdir + '/in.cgrun3'
                srcfyl = lmpdir + '/in.cgrun4'
                desfyl = destdir + '/in.cgrun4'
                shutil.copy2(srcfyl, desfyl)
                srcfyl = lmpdir + '/in.cgpair_MC'
                desfyl = destdir + '/in.cgpair_MC'
                shutil.copy2(srcfyl, desfyl)
                srcfyl = lmpdir + '/jobrestart.sh'
                desfyl = destdir + '/jobrestart.sh'
                shutil.copy2(srcfyl, desfyl)

                print( "Submitting Jobs..")

                subprocess.call(["qsub","jobfailsafe.sh"])

