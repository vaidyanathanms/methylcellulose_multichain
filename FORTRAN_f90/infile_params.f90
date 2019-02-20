!--- To generate pair coeff files for CG MC--------------------------
!--- Input required temperature -------------------------------------
!--- Version: Nov-29-2017--------------------------------------------
!********************************************************************

MODULE INPPARAMS

  IMPLICIT NONE
  
  INTEGER, PARAMETER :: N = 100
  INTEGER, PARAMETER :: natomtypes = 8
  INTEGER, PARAMETER :: grafts = 0

  REAL, DIMENSION(1:natomtypes) :: epsval, sigma, rcut
  REAL, PARAMETER :: bondcoeff = 500
  REAL, PARAMETER :: bondeq = 1.0
  REAL, PARAMETER :: anglcoeff =15
  REAL, PARAMETER :: thetaeq = 165
  REAL, PARAMETER :: dihdcoeff = 2.0
  INTEGER, PARAMETER :: neq = 1
  INTEGER, PARAMETER :: deq = 1

  REAL :: tempval

END MODULE INPPARAMS

!--------------------------------------------------------------------
