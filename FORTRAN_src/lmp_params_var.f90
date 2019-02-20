!-----To generate input parameters for LAMMPS CG-MethylCellulose-----
!----- Main File: lammps_inp.f90-------------------------------------
!------Version: Jan-25-2018------------------------------------------
!--------------------------------------------------------------------

MODULE PARAMS

  USE RAN_NUMBERS

  IMPLICIT NONE

! Parameter data for creating the data file

  INTEGER, PARAMETER :: N = py_nchains
  INTEGER, PARAMETER :: M = py_nmons
  INTEGER, PARAMETER :: totpart = N*M
  REAL,    PARAMETER :: DS_MC = 1.8
  REAL,    PARAMETER :: DS_fac = 2.1
  REAL,    PARAMETER :: mean_tol = 0.01
  REAL,    PARAMETER :: basemass = 162.1406
  REAL,    PARAMETER :: graftperMC = 0.0

! Box details

  REAL, PARAMETER :: insidebox = 100.0 !start of all polymers
  REAL, PARAMETER :: density = py_dens !!5000.0/(600.0**3.0)
  REAL :: boxl_x, boxl_y, boxl_z
  REAL :: volbox

! Flags for creating the data file

  INTEGER, PARAMETER :: grafts = 0
  INTEGER, PARAMETER :: stretched = 0
  INTEGER, PARAMETER :: numatomtypes = 8
  INTEGER, PARAMETER :: numbondtypes = 1
  INTEGER, PARAMETER :: numangltypes = 1
  INTEGER, PARAMETER :: numdihdtypes = 1
  INTEGER, PARAMETER :: bondtype = 1
  INTEGER, PARAMETER :: angltype = 1
  INTEGER, PARAMETER :: dihdtype = 1
  INTEGER, PARAMETER :: outfile  = 17

! Global Arrays involved in creating data file
  
  REAL,    DIMENSION(1:N,1:3*M) :: rxyz, uxyz
  INTEGER, DIMENSION(1:N,1:3*M) :: ix
  INTEGER, DIMENSION(1:N,1:M) :: atype
  INTEGER, DIMENSION(1:numatomtypes) :: cntatomtype

! Character Arrays for creating the data file name

  CHARACTER (LEN = 5)  :: ext
  CHARACTER (LEN = 60 ):: datafile

  TYPE (RAN_SAVE) :: X
  INTEGER(4) :: S
  
END MODULE PARAMS
