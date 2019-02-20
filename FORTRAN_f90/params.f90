MODULE PARAMS_SOLVATEDMC

  USE OMP_LIB
  IMPLICIT NONE

  ! Required Input Variables

  INTEGER :: initdist
  INTEGER :: nframes, skipfr
  INTEGER :: nwater, nchains
  INTEGER :: atperchain
  INTEGER :: nproc

  !Structural analysis input details

  INTEGER :: rdffreq,rmaxbin,npairs
  REAL    :: rclus_cut,rvolavg,rdomcut,rbinval,oxycut
  INTEGER :: rgfreq, oxyfreq, oxytype,eigfreq
  INTEGER :: oxyrestype, oxyrescut
  INTEGER :: rg_s_freq
  INTEGER :: densfreq, dens_axis, ndentypes, maxden_bin
  REAL    :: normdens,denbinavg

  !Diffusivity details
  
  INTEGER :: diff_nframes

  ! All flags
  
  INTEGER :: rdfcalc, rgcalc, oxycalc, rescalc
  INTEGER :: globeig,indeig,eigcalc
  INTEGER :: rgall, rgavg, rgsys
  INTEGER :: denscalc
  INTEGER :: diffcalc

  ! File names and unit Numbers
  
  CHARACTER(LEN = 256) :: ana_fname,data_fname,traj_fname,log_fname
  CHARACTER(LEN = 256) :: rdf_fname, dum_fname
  INTEGER, PARAMETER :: anaread = 2,   logout = 3
  INTEGER, PARAMETER :: inpread = 100, rgwrite = 400,rgavgwrite = 300
  INTEGER, PARAMETER :: dumwrite = 200, rgswrite = 250
  INTEGER, PARAMETER :: eigallwrite = 280,eigwrite=275

  !Math Constants

  REAL*8, PARAMETER :: pival  = 3.14159265359
  REAL*8, PARAMETER :: pi2val = 2.0*pival

  !Global analysis variables and arrays

  INTEGER :: atomflag, velflag, bondflag, anglflag, dihdflag,imprflag
  INTEGER :: ntotatoms, ntotbonds, ntotangls,ntotdihds,ntotimprs
  INTEGER :: ntotatomtypes,ntotbondtypes,ntotangltypes,ntotdihdtypes&
       &,ntotimprtypes

  !Lammps trajectory file read details

  REAL :: box_xl,box_yl,box_zl, boxval
  INTEGER*8 :: timestep

  !Structural variables

  INTEGER :: rdfpaircnt
  REAL    :: rvolval
  
  !Structural Average Variables

  REAL :: re2ave, re4ave, rg2ave, rg4ave, b2ave
  
  !Required Arrays - LAMMPS

  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: rxyz_lmp, vel_xyz, charge_lmp&
       &,masses
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: bond_lmp, angl_lmp,&
       & dihd_lmp, impr_lmp,aidvals
  CHARACTER,ALLOCATABLE,DIMENSION(:,:) :: keywords
  REAL,ALLOCATABLE,DIMENSION(:):: boxx_arr, boxy_arr,boxz_arr
  REAL*8,ALLOCATABLE,DIMENSION(:,:):: rx_time, ry_time, rz_time

  !Required Arrays - Structural Quantities

  REAL,ALLOCATABLE,DIMENSION(:,:):: rdfarray, densarray
  REAL,ALLOCATABLE,DIMENSION(:) ::  histoxyarray,tplot_cf
  INTEGER,ALLOCATABLE,DIMENSION(:,:):: pairs_rdf, autocf
  INTEGER, ALLOCATABLE, DIMENSION(:) :: lenarray,oxyarray,oxyresarray&
       &, dentyp_arr
  REAL, ALLOCATABLE, DIMENSION(:,:) :: eigarray

END MODULE PARAMS_SOLVATEDMC
