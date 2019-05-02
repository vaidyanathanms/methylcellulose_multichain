import numpy
import os
import shutil
import subprocess
import sys
import glob

from subprocess import call

restart = 0 # For restarting from given configurations
temp    = '25.0' #Important to give decimal at the end
nchains = 1
nmons   = [800,1000] #[100,200,300,500,600,800,1000]#,1200,1500]

datafyl = "data_N.data"
maindir = os.getcwd()
scratchdir = '/scratch.global/vaidya/'
lmpdir = '/home/dorfmank/vsethura/allfiles/files_cgmc/src_lmp'
denstype = 1 #1-1.5*contour len, 2-(5000/600^3), 3- ntot/600^3   

for nval in range(len(nmons)):
    
    workdir1 = scratchdir + 'cg_larson25'
    
    if not os.path.isdir(workdir1):
        os.mkdir(workdir1)

    workdir1 = workdir1 + '/temp_' + str(temp)
    if not os.path.isdir(workdir1):
        os.mkdir(workdir1)

    workdirdd2 = workdir1 + '/n_' + str(nmons[nval])
    if not os.path.isdir(workdirdd2):
        os.mkdir(workdirdd2)

    workdird2 = workdirdd2 + '/denstyle_' + str(denstype)
    if not os.path.isdir(workdird2):
        os.mkdir(workdird2)
    

    os.chdir(workdird2)

    destdir = os.getcwd()

    if restart == 0:

        print( "Starting from beginning")

        print( "Copying Files")
        
        # Copy All files

        #--- Copy files for interaction paramters
        srcfyl = maindir + '/cgparams.txt'
        desfyl = destdir + '/cgparams.txt'
        shutil.copy2(srcfyl, desfyl)
        srcfyl = maindir + '/create_infile.f90'
        desfyl = destdir + '/create_infile.f90'
        shutil.copy2(srcfyl, desfyl)
        srcfyl = maindir + '/infile_params_var.f90'
        desfyl = destdir + '/infile_params_var.f90'
        shutil.copy2(srcfyl, desfyl)


        #--- Copy files for datafile
        srcfyl = maindir + '/lmp_params_var.f90'
        desfyl = destdir + '/lmp_params_var.f90'
        shutil.copy2(srcfyl, desfyl)
        srcfyl = maindir + '/lammps_inp.f90'
        desfyl = destdir + '/lammps_inp.f90'
        shutil.copy2(srcfyl, desfyl)
        srcfyl = maindir + '/ran_numbers.f90'
        desfyl = destdir + '/ran_numbers.f90'
        shutil.copy2(srcfyl, desfyl)


        #--- Copy LAMMPS files
        fyle = lmpdir + '/in.cgrun1_lar_var'
        if not os.path.isfile(fyle):
            print(fyle, "not found")
            sys.exit()
        
        srcfyl = lmpdir  + '/in.cgrun1_lar_var'
        desfyl = destdir + '/in.cgrun0'
        shutil.copy2(srcfyl, desfyl)
        srcfyl = lmpdir + '/in.cgrun2_lar'
        desfyl = destdir + '/in.cgrun2'
        shutil.copy2(srcfyl, desfyl)
        srcfyl = lmpdir + '/in.cgrun3_lar'
        desfyl = destdir + '/in.cgrun3'
        shutil.copy2(srcfyl, desfyl)
        srcfyl = lmpdir + '/in.cgrun4_lar'
        desfyl = destdir + '/in.cgrun4'
        shutil.copy2(srcfyl, desfyl)
        srcfyl = lmpdir + '/in.cgpair_MC'
        desfyl = destdir + '/in.cgpair_MC'
        shutil.copy2(srcfyl, desfyl)
        srcfyl = lmpdir + '/jobsingle.sh'
        desfyl = destdir + '/jobsingle.sh'
        shutil.copy2(srcfyl, desfyl)

        
        srcfyl = '/home/dorfmank/vsethura/mylammps/src/lmp_mesabi'
        desfyl = destdir + '/lmp_mesabi'
        shutil.copy2(srcfyl, desfyl)


        # Manipulate Files

        print("Copy Successful - Manipulating Input Files")

        # Generate paircoefficient files

        launch_fyl = 'infile_params.f90'
        fr = open('infile_params_var.f90','r')
        fw = open(launch_fyl,'w')
        fid = fr.read().replace("py_nmons",str(nmons[nval]))
        fw.write(fid)
        fw.close()
        fr.close()
        
        print( "Running FORTRAN script for generating Datafile")
        subprocess.call(["ifort","-r8","-qopenmp","infile_params.f90",
                         "create_infile.f90","-o","pairinp.o"])
        
        subprocess.call(["./pairinp.o",str(temp)])
        


        # Generate Datafiles


        if denstype == 1:
            densval = (nchains*nmons[nval]/((1.5*nmons[nval])**(3.0)))
        elif denstype == 2:
            densval = (5000/(600.0**(3.0)))
        else:
            densval = (nchains*nmons[nval]/(600.0**(3.0)))

        launch_fyl = 'lmp_params.f90'
        fr = open('lmp_params_var.f90','r')
        fw = open(launch_fyl,'w')
        fid = fr.read().replace("py_nchains",str(nchains)).\
            replace("py_nmons",str(nmons[nval])).\
            replace("py_dens",str(densval))
        fw.write(fid)
        fw.close()
        fr.close()
        
        # Generate Required Files for LAMMPS

        print( "Running FORTRAN script for generating Datafile")
        subprocess.call(["ifort","-r8","-qopenmp","ran_numbers.f90",
                         "lmp_params.f90","lammps_inp.f90","-o","lmpinp.o"])
        
        subprocess.call(["./lmpinp.o","lj"])
        

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


        lmpdataname = "MCdata_" + str(nmons[nval]) + ".lj"

        if not os.path.isfile(lmpdataname):
            print ("Warning: File not found", lmpdataname)
            sys.exit()
        
        launch_fyl = 'in.cgrun1'
        fr = open('in.cgrun0','r')
        fw = open(launch_fyl,'w')
        fid = fr.read().replace("lmpdatafyl",str(lmpdataname))
        fw.write(fid)
        fw.close()
        fr.close()

        print( "Submitting Jobs..")

        subprocess.call(["qsub","jobsingle.sh"])

        os.chdir(maindir)

    elif restart == 1:

        print ("Restarting jobs ..")
        
        # Find if restart file is present
    
        rest2files = destdir +'restart3_melt2'
        rest1files = destdir +'restart3_melt1'
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


