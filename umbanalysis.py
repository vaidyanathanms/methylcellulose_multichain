#-----Automating SMD and US Calculations--------------------------
#-----Version-1: April-23-2018------------------------------------


import numpy
import os
import shutil
import subprocess
import sys
import glob
import os.path
from subprocess import call
import re
import random

#--SMD or UmbrellaSampling---------------------------------------

logic_umb = 1 #0 - SMD for initializing US; 1 - US
smddir    = 3 #2 for smd2 and so on. if 0 it is smd

#-Temp/Trial/config Details--------------------------------------
temp     = 50
trialnum = 61
config   = 'diffconfig'

#--Chain Details-------------------------------------------------
numchains    = [2]
monsperchain = 1000
archarr      = [2]#,2,3] #1-alternate, 2-block, 3 -jammed
initdist     = '0.08'


#--Umbrella Sampling Details-------------------------------------
if smddir != 0:
    USdirnameprefix  = 'US_dir' + str(smddir) + '_'
else:
    USdirnameprefix  = 'US_dir_'

umb_centers = [16.5,19.5,22.5] #[2.0,3.0,4.0,5.0,6.0,7.0,9.0,11.0,13.0,15.0,18.0,21.0,24.0]
force_constant =  15.0

#--US details for FORTRAN CODE-----------------------------------
tolval  = 1.0
axisval = 0
comval  = 1


#--SMD Details---------------------------------------------------
targinit = 1.0
targfin  = 40.0

#--Main Directories----------------------------------------------
maindir = os.getcwd()
scratchdir = '/scratch.global/vaidya/'
lmpdir = '/home/dorfmank/vsethura/allfiles/files_cgmc/src_lmp'

for ifree in range(len(numchains)):

    
    print("Free Number of Chains: ", numchains[ifree])
    workdir1 = scratchdir + 'cg_mccalc/temp_'+str(temp)+'/'+ \
        config+'/n_'+str(monsperchain)+'/trial'+str(trialnum)

    #--Ask user to input final chain number if nchain>2--------
    
    if logic_umb == 0:

        if numchains[ifree] > 2:
            pychinp = input("Input the reference chain")
            ref_init = (int(pychinp)-1)*monsperchain + 1
            pychinp = input("Input the chain that needs to be moved")
            print("moving chain number", pychinp)
            pyinit   = (int(pychinp)-1)*monsperchain + 1
        else:
            ref_init = 1
            pyinit = monsperchain
    
        pyfin  = pyinit + monsperchain - 1
        ref_fin = ref_init + monsperchain - 1
        print("Moving chain init/last mon for SMD",pyinit,pyfin)
        print("Ref chain init/last mon for SMD",ref_init,ref_fin)

    if not os.path.isdir(workdir1):
        print(workdir1,"not found")
        sys.exit()

    workdir2 = workdir1 + '/n' + str(numchains[ifree])
	
    if not os.path.isdir(workdir2):
        print(workdir2,"not found")
        sys.exit()
        
    for iarch in range(len(archarr)):

        ranval = random.random()
        refchain = 1
        free_ntot = numchains[ifree]*monsperchain

        if archarr[iarch] == 1:
            print( "Archval: Alternate")
            dirstr = 'alternate'
            fylstr = 'alter'
        elif archarr[iarch] == 2:
            print( "Archval: Block")
            dirstr = 'block'
            fylstr = 'block'
        elif archarr[iarch] == 3:
            print( "Archval: Jammed")
            dirstr = 'jammed'
            fylstr = 'jam'
        else:
            print( "Unknown Architecture")
          
        coredir = workdir2 + '/' + dirstr + '/dist_' + initdist 

        if not os.path.isdir(coredir):
            print(coredir,"not found")
            sys.exit()

        if smddir == 0:
            smdworkdir = coredir + '/smd'
        else:
            smdworkdir = coredir + '/smd' + str(smddir)

        if not os.path.isdir(smdworkdir):
            if logic_umb == 1:
                print(smdworkdir,"not found")
                sys.exit()
            else:
                os.mkdir(smdworkdir)
                
        os.chdir(smdworkdir)
        destdir = os.getcwd()
          

#--------Prepare for US Simulations--------------------------------------

        if logic_umb == 1:
            
            print("Preparing for US Simulations ...")
            ref_data = smdworkdir + "/datainp.data"
            
            if not os.path.exists(ref_data):
                print("RefData: ", ref_data, "not found")
                sys.exit()

            for ucen in range(len(umb_centers)):
                
                umbdir = coredir + '/' + USdirnameprefix + \
                         str(umb_centers[ucen])
            
                print("US directory under construction", umbdir)
                if not os.path.isdir(umbdir):
                    os.mkdir(umbdir)
                
                os.chdir(umbdir)
                umbworkdir = os.getcwd()
            
                srcfyl = lmpdir  + '/in.umbsample_var'
                desfyl = umbworkdir + '/in.umbsample_var'
                shutil.copy2(srcfyl, desfyl)
                
                srcfyl = lmpdir  + '/umbinp_var'
                desfyl = umbworkdir + '/umbinp_var'
                shutil.copy2(srcfyl, desfyl)
                
                srcfyl = lmpdir  + '/jobumb.sh'
                desfyl = umbworkdir + '/jobumb.sh'
                shutil.copy2(srcfyl, desfyl)
                
                srcfyl = maindir + '/extract_conf.f90'
                desfyl = smdworkdir + '/extract_conf.f90'
                shutil.copy2(srcfyl,desfyl)

                srcfyl = maindir + '/extract_params.f90'
                desfyl = smdworkdir + '/extract_params.f90'
                shutil.copy2(srcfyl,desfyl)

                srcfyl = maindir + '/extract_inp_var.txt'
                desfyl = smdworkdir + '/extract_inp_var.txt'
                shutil.copy2(srcfyl,desfyl)

                srcfyl = lmpdir + '/in.cgpair_MC'
                desfyl = umbworkdir + '/in.cgpair_MC'
                shutil.copy2(srcfyl,desfyl)
                              
                srcfyl = '/home/dorfmank/vsethura/mylammps/src/lmp_mesabi'
                desfyl = umbworkdir + '/lmp_mesabi'
                shutil.copy2(srcfyl, desfyl)


#---Now check whether SMD simulations are okay-----------------------------
#---Call FORTRAN CODE - extract_conf.f90recursively with different
#---centers---

#obtain the reference chain from SMD simulation inputs and feed here

                os.chdir(smdworkdir)
                fyl_smd = smdworkdir + '/smdcolfile.inp'
                
                print("Ref SMD working directory for US: ",smdworkdir)
                with open(fyl_smd) as fp:
                    for lyne in fp:
                        if "group2" in lyne:
                            content = re.findall(r'\d+',lyne)
                            
                pyinit = int(content[1])
                pyfin  = int(content[2])

                with open(fyl_smd) as fp:
                    for lyne in fp:
                        if "group1" in lyne:
                            content = re.findall(r'\d+',lyne)

                ref_init = int(content[1])
                ref_fin  = int(content[2])
                print("ref chain init/last mon for US",ref_init,ref_fin)
                print("moving chain init/last mon for US",pyinit,pyfin)
                
                init_fyl = 'extract_inp.txt'
                fr = open('extract_inp_var.txt','r')
                fw = open(init_fyl,'w')
                fid = fr.read().replace("py_smddata","datainp.data").\
                    replace("py_refinit",str(ref_init)).\
                    replace("py_reffin", str(ref_fin)).\
                    replace("py_init",str(pyinit)).\
                    replace("py_fin",str(pyfin)).\
                    replace("py_savedist",str(umb_centers[ucen])).\
                    replace("py_tol", str(tolval)).\
                    replace("py_axis",str(axisval)).\
                    replace("py_comflag", str(comval))

                fw.write(fid)
                fw.close()
                fr.close()


                subprocess.call(["ifort","-mkl","-qopenmp","-r8",\
                                 "extract_params.f90","extract_conf.f90",\
                                 "-o","extract.o"])
                subprocess.call(["./extract.o","extract_inp.txt"])


                checkinitdata = smdworkdir + "/init_datafile"
                desfyl = umbworkdir + '/init_datafile'

                if os.path.exists(checkinitdata):

                    print(checkinitdata, "prepared for US simulations.")
                    shutil.copy2(checkinitdata,desfyl)
                    
                    #Delete & keep a copy under different name for
                    #future copying if necessary

                    prev_US_fyl = smdworkdir + '/previnit_datafile'
                    shutil.copy2(checkinitdata,prev_US_fyl)
                    os.remove(checkinitdata)
                    
                else: #if US could not find a file near that place -
                    #use SMD initialized file or latest US file
                    
                    print(checkinitdata, "not prepared for", \
                              numchains[ifree],"at",umb_centers[ucen])

                    prev_US_fyl = smdworkdir + '/previnit_datafile'

                    if os.path.exists(prev_US_fyl):

                        print("Copying previous US file for these purposes")
                        shutil.copy2(prev_US_fyl,desfyl)

                    else:

                        print("Copying SMD file for these purposes")
                        shutil.copy2(ref_data,desfyl)
                        
                #---Copy datafiles for backup before analyzing
                
                datadir = maindir + '/USbackupdir'

                if not os.path.isdir(datadir):
                    os.mkdir(datadir)

                datadir2 = datadir + '/n_' + str(numchains[ifree])

                if not os.path.isdir(datadir2):
                    os.mkdir(datadir2)

                datafyle = datadir2 + '/n_' + str(numchains[ifree])\
                    + '_' + dirstr + '_' +  str(umb_centers[ucen])

                shutil.copy2(desfyl,datafyle)
                
                os.chdir(umbworkdir)
                umbmain_fyl = 'umb_colfile.inp'
                fr = open('umbinp_var','r')
                fw = open(umbmain_fyl,'w')
                fid = fr.read().replace("py_init",str(pyinit)).\
                    replace("py_fin",str(pyfin)).\
                    replace("py_refinit",str(ref_init)).\
                    replace("py_reffin",str(ref_fin)).\
                    replace("py_cval",str(umb_centers[ucen])).\
                    replace("py_fcon", str(force_constant))
                fw.write(fid)
                fw.close()
                fr.close()
                
                umblmp_fyl = 'in.umbsample'
                fr = open('in.umbsample_var','r')
                fw = open(umblmp_fyl,'w')
                fid = fr.read().replace("py_init",str(pyinit)).\
                    replace("py_fin",str(pyfin)).\
                    replace("ref_init",str(ref_init)).\
                    replace("ref_fin", str(ref_init)) 
                fw.write(fid)
                fw.close()
                fr.close()
                
                print("Submitting US jobs for",numchains[ifree],"at" \
                          ,umb_centers[ucen])
                subprocess.call(["qsub", "jobumb.sh"])
                
        else:


#----------------Preparing for SMD Simulations-----------------------------

            print("Preparing for SMD simulations ...")
            print("Copying restart files")
            print("Coredir:\t", coredir)

            srcfyl = '/home/dorfmank/vsethura/mylammps/src/lmp_mesabi'
            desfyl = smdworkdir + '/lmp_mesabi'
            shutil.copy2(srcfyl, desfyl)

            srcfyl = maindir + '/archive_restart'+'/restartfiles_n_' + \
                str(numchains[ifree]) + '_' +  dirstr + '.tar.gz'

# Check files are already extracted - this takes the most time

            archfiles = destdir + '/panfs/roc/scratch/vaidya/pework/n_'+ \
                str(numchains[ifree]) + "/" +dirstr +'/restartfiles_n_'+\
                str(numchains[ifree]) + '_' + dirstr +  '/' +  \
                'archival_*.restart'

            corerestfiles = coredir +'/restart*'
            corerest_list = glob.glob(corerestfiles)
            
            list_of_files = glob.glob(archfiles)
            oldtime = 0

            if list_of_files:

                for fyl in list_of_files:

                    splitdata = fyl.split("/")
                    fylname  = splitdata[len(splitdata)-1]
                    timedata = fylname.split("_")
                    timeval  = float(timedata[1].split(".")[0])

                    if timeval > oldtime:

                        oldtime = timeval

                finarchfyl = destdir + '/panfs/roc/scratch/vaidya/pework/n_'+ \
                    str(numchains[ifree]) + "/" +dirstr +'/restartfiles_n_'+\
                    str(numchains[ifree]) + '_' + dirstr +  '/' +  \
                    'archival_' + str(int(oldtime)) + '.restart'

                desfyl = destdir + '/archival_' + str(int(oldtime)) + '.restart'
                shutil.copy2(finarchfyl,desfyl)

                subprocess.call(["mpirun","-np","24", "./lmp_mesabi", "-r",
                                 desfyl, "datainp.data"])


            elif os.path.exists(srcfyl):
		
                print( "Restart file for LAMMPS: ", srcfyl)
                print("Current Dir:", os.getcwd())
                desfyl= destdir+'/restartfiles_n_'+str(numchains[ifree]) \
                    + '_' +  dirstr + '.tar.gz'
            
                shutil.copy2(srcfyl, desfyl)

# Untar restart files and copy it as restart1
            
                subprocess.call(["tar","-xvzf", desfyl])

                archfiles = destdir + '/panfs/roc/scratch/vaidya/pework/n_'+ \
                    str(numchains[ifree]) + "/" +dirstr +'/restartfiles_n_'+\
                    str(numchains[ifree]) + '_' + dirstr +  '/' +  \
                    'archival_*.restart'

                list_of_files = glob.glob(archfiles)
            
                for fyl in list_of_files:

                    splitdata = fyl.split("/")
                    fylname  = splitdata[len(splitdata)-1]
                    timedata = fylname.split("_")
                    timeval  = float(timedata[1].split(".")[0])
                    
                    if timeval > oldtime:

                        oldtime = timeval

                finarchfyl = destdir + '/panfs/roc/scratch/vaidya/pework/n_'+ \
                    str(numchains[ifree]) + "/" +dirstr +'/restartfiles_n_'+\
                    str(numchains[ifree]) + '_' + dirstr +  '/' +  \
                    'archival_' + str(int(oldtime)) + '.restart'

                desfyl = destdir + '/archival_' + str(int(oldtime)) + '.restart'
                shutil.copy2(finarchfyl,desfyl)
                
                #Check the coredir
            elif corerest_list:

                print("Successfully found restart files in the coredir")
                latest_rest = max(corerest_list,key=os.path.getctime)
                desfyl = destdir + '/smdarchival.restart'
                shutil.copy2(latest_rest,desfyl)
                subprocess.call(["mpirun","-np","24", "./lmp_mesabi", "-r",
                                 desfyl, "datainp.data"])

            else:
		
                print(srcfyl, "not found in", maindir)
                os.chdir(maindir)


                            
            srcfyl = lmpdir   + '/in.smd_var'
            desfyl = smdworkdir + '/in.smd_var'
            shutil.copy2(srcfyl, desfyl)

            srcfyl = lmpdir   + '/smdcolfile_var.inp'
            desfyl = smdworkdir + '/smdcolfile_var.inp'
            shutil.copy2(srcfyl, desfyl)

            srcfyl = lmpdir   + '/jobsmd.sh'
            desfyl = smdworkdir + '/jobsmd.sh'
            shutil.copy2(srcfyl, desfyl)

            srcfyl = lmpdir     + '/in.cgpair_MC'
            desfyl = smdworkdir + '/in.cgpair_MC'
            shutil.copy2(srcfyl,desfyl)

            smdlmp_fyl = 'in.smd'
            fr = open('in.smd_var','r')
            fw = open(smdlmp_fyl,'w')
            fid = fr.read().replace("py_init",str(pyinit)).\
                replace("py_fin",str(pyfin)).\
                replace("ref_init",str(ref_init)).\
                replace("ref_fin", str(ref_fin))
            fw.write(fid)
            fw.close()
            fr.close()
                        
            smdcol_fyl = 'smdcolfile.inp'
            fr = open('smdcolfile_var.inp','r')
            fw = open(smdcol_fyl,'w')
            fid = fr.read().replace("py_init",str(pyinit)).\
                replace("py_fin",str(pyfin)).\
                replace("ref_init",str(ref_init)).\
                replace("ref_fin", str(ref_fin)).\
                replace("py_targinit",str(targinit)).\
                replace("py_targfin",str(targfin)).\
                replace("py_fcon", str(force_constant))
            fw.write(fid)
            fw.close()
            fr.close()

            print("Submitting SMD jobs for", numchains[ifree])
            subprocess.call(["qsub", "jobsmd.sh"])
