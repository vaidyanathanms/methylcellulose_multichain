!-------------To generate replicated system for CG MethylCellulose---
!-------------Newer Version compared to replicated_poly.f90----------
!-------------Has the capability to control the # of replicas--------
!-------------Version: Dec-13-2017-----------------------------------
!********************************************************************

MODULE PARAMS_REP

  USE RAN_NUMBERS
  IMPLICIT NONE

  INTEGER :: axismove !1-x,2-y,3-z
  REAL :: movedist,moveperc !default 5%--20% of box size
  INTEGER:: densflag, moveflag, reprandom
  REAL :: boxxl, boxyl, boxzl, densval
  INTEGER :: alternate
  REAL :: nxvec, nyvec, nzvec

  INTEGER :: ntotatoms, ntotbonds, ntotangls, ntotdihds
  INTEGER :: natomtypes, nbondtypes,nangltypes, ndihdtypes
  INTEGER :: anglflag, dihdflag,atomflag,bondflag
  INTEGER :: natomsperpoly, nbondsperpoly, nanglsperpoly,ndihdsperpoly
  INTEGER :: n_rings, nnewchains, ntotchains, nflexchains
  REAL    :: flexper
  REAL    :: rgring, rxring, ryring, rzring

  CHARACTER(LEN=256) :: data_fname,log_fname,ana_fname,out_fname
  INTEGER, PARAMETER :: anaread = 100, inpread=200, logout = 300&
       &,outfile=500

  REAL, ALLOCATABLE, DIMENSION(:,:) :: angltypearr,dhdtypearr&
       &,bondtypearr, rxyz_lmp, Masses, charge_lmp

  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: bond_lmp,angl_lmp,dihd_lmp&
       &,aidvals, imgflag, ixyz_lmp

  TYPE (RAN_SAVE) :: X
  INTEGER(4) :: S

END MODULE PARAMS_REP

!--------------------------------------------------------------------

PROGRAM REPLICATE

  USE PARAMS_REP
  IMPLICIT NONE

  CALL READINP()
  CALL READDATAFILE()
  CALL SYSTEM_CLOCK(S)
  CALL RAN_INIT(S,X)
  CALL DECIDE_REPLICA_DIST()
  CALL COMPUTE_RG_RING()

  IF(flexper == 0.0)  THEN
     CALL GENERATE_REPLICA_POS()
  ELSE 
     IF(alternate == 0) THEN
        CALL GENERATE_REPLICA_POS()
        CALL GENERATE_FLEXIBLE_POS()
     ELSE
        CALL GENERATE_BLOCK_POS()
     END IF
  END IF

  CALL WRAPALL()
  CALL GENERATE_IMAGE_FLAGS()
  CALL GENERATE_ALL_TOPO()
  CALL WRITELAMMPSDATA()
  CALL DEALLOCATEARRAYS()
  
END PROGRAM REPLICATE

!--------------------------------------------------------------------

SUBROUTINE READINP()

  USE PARAMS_REP
  IMPLICIT NONE
  
  INTEGER :: nargs,ierr,logflag,scrflag
  CHARACTER(256) :: dumchar
  LOGICAL :: fexist

  flexper = 0.0

  nargs = IARGC()
  IF(nargs .NE. 1) STOP "Input incorrect"

  logflag = 0; scrflag = 0; densflag = 0; moveflag = 0
  reprandom = 0; alternate = 0

  CALL GETARG(nargs,ana_fname)

  OPEN(unit = anaread,file=trim(ana_fname),action="read",status="old"&
       &,iostat=ierr)
  
  IF(ierr /= 0) THEN

     PRINT *, trim(ana_fname), "not found"
     STOP

  END IF

  DO

     READ(anaread,*,iostat=ierr) dumchar

     IF(ierr .LT. 0) EXIT

     IF(trim(dumchar) == 'input_datafile') THEN

        READ(anaread,*,iostat=ierr) data_fname
        IF(ierr /= 0) STOP "Datafile not found"
        
     ELSEIF(trim(dumchar) == 'log_file') THEN

        READ(anaread,*,iostat=ierr) log_fname
        logflag = 1

     ELSEIF(trim(dumchar) == 'majoraxis') THEN

        READ(anaread,*,iostat=ierr) axismove
        READ(anaread,*,iostat=ierr) nxvec,nyvec,nzvec

        IF(axismove == 1) THEN
           nxvec = 1; nyvec = 0; nzvec = 0
        ELSEIF (axismove == 2) THEN
           nxvec = 0; nyvec = 1; nzvec = 0
        ELSEIF (axismove == 3) THEN
           nxvec = 0; nyvec = 0; nzvec = 1
        END IF

        PRINT *, nxvec, nyvec, nzvec

     ELSEIF(trim(dumchar) == 'num_newchains') THEN

        READ(anaread,*,iostat=ierr) nnewchains

     ELSEIF(trim(dumchar) == 'movereplica') THEN

        moveflag = 1
        READ(anaread,*,iostat=ierr) movedist
        moveperc = movedist

     ELSEIF(trim(dumchar) == 'replica_random') THEN

        READ(anaread,*,iostat=ierr) reprandom

     ELSEIF(trim(dumchar) == 'ordering') THEN

        READ(anaread,*,iostat=ierr) alternate

        IF(alternate == 0) PRINT *, "MODE:Random initial rings/blocks"
        IF(alternate == 1) PRINT *, "MODE:Chains & Rings in blocks"
        IF(alternate == 2) PRINT *, "MODE:Chains & Rings alternating"
        IF(alternate == 3) PRINT *, "MODE:One chain nucleation"
        IF(alternate == 4) PRINT *, "MODE:Jammed Rings"
        IF(alternate .GT. 4) STOP "Unknown identifier for ordering"

     ELSEIF(trim(dumchar) == 'percentflex') THEN
        
        READ(anaread,*,iostat=ierr) flexper

        IF(flexper .LT. 0 .OR. flexper .GT. 1) STOP "percent flex shou&
             &ld be between 0 and 1"

     ELSEIF(trim(dumchar)=='density') THEN
        
        READ(anaread,*,iostat=ierr) densval
        densflag = 1

     ELSE
        
        PRINT *, "unknown keyword", trim(dumchar)
        STOP

     END IF

  END DO

  IF(reprandom == 1 .AND. alternate /= 0) STOP "Conflicting options fo&
       &r replicate and ordering"

  IF(logflag == 0) log_fname = trim(adjustl("logrep"//data_fname//"&
       &.log"))

  OPEN(unit = logout,file=trim(log_fname),action="write",status="repla&
       &ce",iostat=ierr)

  WRITE(logout,*) "Analysis input file read finished .."
  WRITE(logout,*) "Basic Datafile: ", trim(data_fname)

  ntotchains = nnewchains + 1
  nflexchains = nnewchains*flexper
  n_rings = ntotchains - nflexchains

  WRITE(logout,*) "Total number of chains: ", ntotchains
  WRITE(logout,*) "New number of chains: ", nnewchains
  WRITE(logout,*) "Number of Rings: ", n_rings
  WRITE(logout,*) "Flexible number of chains: ", nflexchains

  PRINT *, "Total number of chains: ", ntotchains
  PRINT *, "New number of chains: ", nnewchains
  PRINT *, "Number of Rings: ", n_rings
  PRINT *, "Flexible number of chains: ", nflexchains

  IF(reprandom == 1) WRITE(logout,*) "Rings at random"
  IF(reprandom == 0) WRITE(logout,*) "Rings in ordered fashion"

  IF(reprandom) THEN
     IF(axismove == 1) PRINT *, "Replicating along X axis"
     IF(axismove == 2) PRINT *, "Replicating along Y axis"
     IF(axismove == 3) PRINT *, "Replicating along Z axis"
     IF(axismove == 4) PRINT *, "Replicating along a given vector"
  END IF

END SUBROUTINE READINP

!--------------------------------------------------------------------

SUBROUTINE READDATAFILE()

  USE PARAMS_REP
  IMPLICIT NONE

  INTEGER :: i,ierr,u,AllocateStatus,j
  INTEGER :: flag, cntr, nwords
  INTEGER :: aid, molid,atype,ix,iy,iz
  REAL    :: charge,rx,ry,rz
  REAL    :: xlo,xhi,ylo,yhi,zlo,zhi
  CHARACTER(256) :: rline,dumchar
  INTEGER :: angtypeval, dhdtypeval, bontypeval,reftype
  REAL :: kval, thetaval, rval
  REAL :: k1, k2, k3, k4

  atomflag = 0; bondflag = 0; anglflag = 0; dihdflag = 0
  OPEN(unit=inpread,file = trim(adjustl(data_fname)),action =&
       & 'read', status='old',iostat=ierr) 

  IF(ierr .NE. 0) THEN
     PRINT *, trim(data_fname) 
     STOP "Data file not found"
  END IF

  WRITE(logout,*) "Datafile used is :", trim(adjustl(data_fname))

  READ(inpread,*)
  READ(inpread,*)

  DO i = 1,8

     READ(inpread,*) u, dumchar
     
     IF(dumchar == "atoms") THEN
        natomsperpoly = u
     ELSEIF(dumchar == "bonds") THEN
        nbondsperpoly = u
     ELSEIF(dumchar == "angles") THEN
        nanglsperpoly = u
     ELSEIF(dumchar == "dihedrals") THEN
        ndihdsperpoly = u
     ELSEIF(dumchar == "atom" .OR. dumchar == "atomtypes") THEN
        natomtypes = u
     ELSEIF(dumchar == "bond" .OR. dumchar == "bondtypes") THEN
        nbondtypes = u
     ELSEIF(dumchar == "angle" .OR. dumchar == "atomtypes") THEN
        nangltypes = u
     ELSEIF(dumchar == "dihedral" .OR. dumchar == "dihedraltypes") THEN
        ndihdtypes = u
     ELSEIF(dumchar == "Masses") THEN
        
        ALLOCATE(masses(natomtypes,1),stat = AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate masses"
        
        DO j = 1,natomtypes
           
           READ(inpread,*) u, masses(u,1)
           
        END DO

     END IF

  END DO

  IF(natomsperpoly == 0) STOP "No atoms found"
  
  READ(inpread,*) 
  READ(inpread,*) xlo, xhi
  READ(inpread,*) ylo, yhi
  READ(inpread,*) zlo, zhi
  
  boxxl = xhi-xlo
  boxyl = yhi-ylo
  boxzl = zhi-zlo
  
  IF(densflag == 0) PRINT *, "Boxsizes", boxxl, boxyl, boxzl

  WRITE(logout,*) "STATISTICS: "
  WRITE(logout,*) "natoms/natomtypes", natomsperpoly,natomtypes
  WRITE(logout,*) "nbonds/nbondtypes", nbondsperpoly,nbondtypes
  WRITE(logout,*) "nangls/nangltypes", nanglsperpoly,nangltypes
  WRITE(logout,*) "ndihds/ndihdtypes", ndihdsperpoly,ndihdtypes

  PRINT*, "STATISTICS: "
  PRINT *, "natoms/natomtypes", natomsperpoly,natomtypes
  PRINT *, "nbonds/nbondtypes", nbondsperpoly,nbondtypes
  PRINT *, "nangls/nangltypes", nanglsperpoly,nangltypes
  PRINT *, "ndihds/ndihdtypes", ndihdsperpoly,ndihdtypes

  flag = 0

  ntotatoms = natomsperpoly*ntotchains
  ntotbonds = nbondsperpoly*ntotchains
  ntotangls = nanglsperpoly*ntotchains
  ntotdihds = ndihdsperpoly*ntotchains

  PRINT*, "FULL STATISTICS: "
  PRINT *, "natoms/natomtypes", ntotatoms,natomtypes
  PRINT *, "nbonds/nbondtypes", ntotbonds,nbondtypes
  PRINT *, "nangls/nangltypes", ntotangls,nangltypes
  PRINT *, "ndihds/ndihdtypes", ntotdihds,ndihdtypes


  IF(densflag == 1) THEN

     boxxl = (REAL(ntotatoms)/densval)**(1.0/3.0)
     boxyl = (REAL(ntotatoms)/densval)**(1.0/3.0)
     boxzl = (REAL(ntotatoms)/densval)**(1.0/3.0)

     PRINT *, "Boxsizes", boxxl, boxyl, boxzl
     PRINT *, "Target Density", densval
     
  END IF

  PRINT *, "density", REAL(ntotatoms)/(boxxl*boxyl*boxzl)

  CALL ALLOCATE_ARRAYS()

  DO 
     
     READ(inpread,*,iostat=ierr) dumchar

     IF(ierr .LT. 0) EXIT


     IF(trim(dumchar) == "Atoms") THEN
             
        atomflag = 1
        print *, "Reading ", trim(dumchar), " info"

        DO j = 1,natomsperpoly

           READ(inpread,*) aid,molid,atype,rx,ry,rz,ix,iy,iz

           aidvals(aid,1)  = aid
           aidvals(aid,2)  = molid
           aidvals(aid,3)  = atype
           rxyz_lmp(aid,1) = rx
           rxyz_lmp(aid,2) = ry
           rxyz_lmp(aid,3) = rz
        
        END DO

     END IF

     IF(trim(dumchar) == "Masses") THEN

        ALLOCATE(masses(natomtypes,1),stat = AllocateStatus)
        IF(AllocateStatus/=0) STOP "did not allocate masses"

        DO j = 1,natomtypes
           
           READ(inpread,*) u, masses(u,1)
           
        END DO

     END IF


     IF(trim(dumchar) == "Bonds") THEN
             
        print *, "Reading ", trim(dumchar), " info"

        DO j = 1,nbondsperpoly

           READ(inpread,*) bond_lmp(j,1),bond_lmp(j,2),bond_lmp(j,3)&
                &,bond_lmp(j,4)

        END DO

     END IF


     IF(trim(dumchar) == "Angles") THEN
             
        print *, "Reading ", trim(dumchar), " info"

        DO j = 1,nanglsperpoly

           READ(inpread,*) angl_lmp(j,1),angl_lmp(j,2),angl_lmp(j,3)&
                &,angl_lmp(j,4),angl_lmp(j,5)

        END DO

     END IF

     IF(trim(dumchar) == "Dihedrals") THEN
             
        print *, "Reading", trim(dumchar), "info"

        DO j = 1,ndihdsperpoly

           READ(inpread,*) dihd_lmp(j,1),dihd_lmp(j,2),dihd_lmp(j,3)&
                &,dihd_lmp(j,4),dihd_lmp(j,5), dihd_lmp(j,6)

        END DO

     END IF

  END DO


  CLOSE(inpread)

END SUBROUTINE READDATAFILE

!--------------------------------------------------------------------

SUBROUTINE DECIDE_REPLICA_DIST()

  USE PARAMS_REP
  IMPLICIT NONE

  IF(reprandom == 0) THEN

     IF(moveflag == 0) THEN
        IF(axismove == 1) movedist = 0.05*boxxl
        IF(axismove == 2) movedist = 0.05*boxyl
        IF(axismove == 3) movedist = 0.05*boxzl
        IF(axismove == 4) movedist = 0.05*(boxxl+boxyl+boxzl)/3.0
        moveperc = 0.05
     ELSE
        IF(axismove == 1) movedist = movedist*boxxl
        IF(axismove == 2) movedist = movedist*boxyl
        IF(axismove == 3) movedist = movedist*boxzl
        IF(axismove == 4) movedist = movedist*(boxxl+boxyl+boxzl)/3.0

     END IF

  ELSE

     IF(moveflag == 0) THEN
        movedist = 0.05*MIN(boxxl,boxyl,boxzl)
     ELSE
        movedist = movedist*MIN(boxxl,boxyl,boxzl)
     END IF

  END IF

  PRINT *, "Distance between nearest rings", movedist

END SUBROUTINE DECIDE_REPLICA_DIST

!--------------------------------------------------------------------

SUBROUTINE GENERATE_REPLICA_POS()

  USE PARAMS_REP
  IMPLICIT NONE

  INTEGER :: i, j, k
  REAL :: rx, ry, rz, boxdir, ranval

  PRINT *, "1"
  DO i = 1,n_rings

     IF(reprandom==0) THEN
     
        DO j = 1, natomsperpoly
           
           k = natomsperpoly*(i-1) + j
           
           IF(axismove == 1) THEN
              
              rxyz_lmp(k,1) = rxyz_lmp(j,1)+movedist*REAL(i-1)
              rxyz_lmp(k,2) = rxyz_lmp(j,2)
              rxyz_lmp(k,3) = rxyz_lmp(j,3)
              

           ELSEIF(axismove == 2) THEN
              
              rxyz_lmp(k,1) = rxyz_lmp(j,1)
              rxyz_lmp(k,2) = rxyz_lmp(j,2)+movedist*REAL(i-1)
              rxyz_lmp(k,3) = rxyz_lmp(j,3)
              
           ELSEIF(axismove == 3) THEN
              
              rxyz_lmp(k,1) = rxyz_lmp(j,1)
              rxyz_lmp(k,2) = rxyz_lmp(j,2)
              rxyz_lmp(k,3) = rxyz_lmp(j,3)+movedist*REAL(i-1)

           ELSEIF(axismove == 4) THEN

              rxyz_lmp(k,1) = rxyz_lmp(j,1) + nxvec*movedist*REAL(i-1)
              rxyz_lmp(k,2) = rxyz_lmp(j,2) + nyvec*movedist*REAL(i-1)
              rxyz_lmp(k,3) = rxyz_lmp(j,3) + nzvec*movedist*REAL(i-1)

           ELSE
           
              PRINT *, "Unknown axis number", axismove
              STOP

           END IF


        END DO
              
     ELSE
           
        ranval = RAN1(X)
        
        DO j = 1, natomsperpoly
           
           k = natomsperpoly*(i-1) + j
           
           aidvals(k,1) = k
           aidvals(k,2) = i
           aidvals(k,3) = aidvals(j,3)
           
           IF(ranval .LE. 1.0/3.0) THEN
              
              rxyz_lmp(k,1) = rxyz_lmp(j,1)+movedist*REAL(i-1)
              rxyz_lmp(k,2) = rxyz_lmp(j,2)
              rxyz_lmp(k,3) = rxyz_lmp(j,3)
              
           ELSEIF(ranval .LE. 2.0/3.0) THEN
              
              rxyz_lmp(k,1) = rxyz_lmp(j,1)
              rxyz_lmp(k,2) = rxyz_lmp(j,2)+movedist*REAL(i-1)
              rxyz_lmp(k,3) = rxyz_lmp(j,3)
              
           ELSE
              
              rxyz_lmp(k,1) = rxyz_lmp(j,1)
              rxyz_lmp(k,2) = rxyz_lmp(j,2)
              rxyz_lmp(k,3) = rxyz_lmp(j,3)+movedist*REAL(i-1)
              
           END IF
           
        END DO
        
     END IF
     
  END DO
  

END SUBROUTINE GENERATE_REPLICA_POS

!--------------------------------------------------------------------

SUBROUTINE GENERATE_FLEXIBLE_POS()

  USE PARAMS_REP
  IMPLICIT NONE

  INTEGER :: i,j,k,u,v,kinit,nearflag
  REAL, PARAMETER :: r0init  = 0.97
  REAL, PARAMETER :: rmaxsq  = r0init*r0init
  REAL, PARAMETER :: math_pi = 3.14159265359
  REAL :: theta, phi
  REAL :: rx2ij,ry2ij,rz2ij,ri2j,rxe2e,rye2e,rze2e,re2e
  REAL :: rxij, ryij, rzij, rij
  PRINT *, "2"        
  i = n_rings+1; u = 0

  DO WHILE(i .LE. ntotchains)
     
     nearflag = 0

     DO WHILE(nearflag == 0)

        k = (i-1)*natomsperpoly + 1; kinit = k
        
        rxyz_lmp(k,1) = RAN1(X)*boxxl
        rxyz_lmp(k,2) = RAN1(X)*boxyl
        rxyz_lmp(k,3) = RAN1(X)*boxzl
        
        ! If distance from replica is too far, recreate flexible
        ! chains .Not even image distance. Just normal distance
        ! should be close

        rxij = rxyz_lmp(k,1) - rxring
        ryij = rxyz_lmp(k,2) - ryring
        rzij = rxyz_lmp(k,3) - rzring
        rij = sqrt(rxij**2 + ryij**2 + rzij**2)

        IF(rij .LT. REAL((2+u)*movedist)) nearflag = 1

     END DO
     
     k = k + 1
     
     theta       = math_pi*RAN1(X)
     phi         = 2*math_pi*RAN1(X)

     rxyz_lmp(k,1) = rxyz_lmp(k-1,1) + r0init*sin(theta)*cos(phi)
     rxyz_lmp(k,2) = rxyz_lmp(k-1,2) + r0init*sin(theta)*sin(phi)
     rxyz_lmp(k,3) = rxyz_lmp(k-1,3) + r0init*cos(theta)
     
     k = k + 1

     DO WHILE (k .LE. i*natomsperpoly)
        
        theta       = math_pi*RAN1(X)
        phi         = 2*math_pi*RAN1(X)

        rxyz_lmp(k,1) = rxyz_lmp(k-1,1) + r0init*sin(theta)*cos(phi)
        rxyz_lmp(k,2) = rxyz_lmp(k-1,2) + r0init*sin(theta)*sin(phi)
        rxyz_lmp(k,3) = rxyz_lmp(k-1,3) + r0init*cos(theta)
        
        rx2ij = rxyz_lmp(k-2,1) - rxyz_lmp(k,1)
        ry2ij = rxyz_lmp(k-2,2) - rxyz_lmp(k,2)
        rz2ij = rxyz_lmp(k-2,3) - rxyz_lmp(k,3)
        ri2j  = rx2ij**2 + ry2ij**2 + rz2ij**2

        IF(ri2j > 1.0404) THEN
           
           k = k + 1
           
        END IF
        
     END DO
     
     k = i*natomsperpoly
     rxe2e = rxyz_lmp(k,1) - rxyz_lmp(kinit,1)
     rye2e = rxyz_lmp(k,2) - rxyz_lmp(kinit,2)
     rze2e = rxyz_lmp(k,3) - rxyz_lmp(kinit,3)
     
     re2e = rxe2e**2 + rye2e**2 + rze2e**2
     
     IF(sqrt(re2e) .LT. 0.5*REAL(natomsperpoly)*r0init) THEN

        i = i + 1
        u = u + 1

     END IF

  END DO
     
END SUBROUTINE GENERATE_FLEXIBLE_POS

!--------------------------------------------------------------------

SUBROUTINE GENERATE_BLOCK_POS()

  USE PARAMS_REP
  IMPLICIT NONE

  INTEGER :: mult,atype
  INTEGER :: i,j,k,u,v,kinit,nearflag,kref1,kmid
  REAL, PARAMETER :: r0init  = 0.97
  REAL, PARAMETER :: rmaxsq  = r0init*r0init
  REAL, PARAMETER :: math_pi = 3.14159265359
  REAL :: theta, phi,totmass
  REAL :: rx2ij,ry2ij,rz2ij,ri2j,rxe2e,rye2e,rze2e,re2e
  REAL :: rxij, ryij, rzij, rij
  REAL, DIMENSION(n_rings) :: ringxcom,ringycom,ringzcom
  REAL :: rxinit, ryinit, rzinit
  REAL :: dispx, dispy, dispz

  IF(alternate == 1 .OR. alternate == 4) mult = 1
  IF(alternate == 2) mult = 2
  PRINT *, "3", mult, alternate
  !Generate Rings first 
  DO i = 1,n_rings

     DO j = 1, natomsperpoly
           
        k = natomsperpoly*(i-1) + j
        
        IF(axismove == 1) THEN
           
           rxyz_lmp(k,1) = rxyz_lmp(j,1)+movedist*REAL(mult*(i-1))
           rxyz_lmp(k,2) = rxyz_lmp(j,2)
           rxyz_lmp(k,3) = rxyz_lmp(j,3)
           
        ELSEIF(axismove == 2) THEN
           
           rxyz_lmp(k,1) = rxyz_lmp(j,1)
           rxyz_lmp(k,2) = rxyz_lmp(j,2)+movedist*REAL(mult*(i-1))
           rxyz_lmp(k,3) = rxyz_lmp(j,3)
           
        ELSEIF(axismove == 3) THEN
           
           rxyz_lmp(k,1) = rxyz_lmp(j,1)
           rxyz_lmp(k,2) = rxyz_lmp(j,2)
           rxyz_lmp(k,3) = rxyz_lmp(j,3)+movedist*REAL(mult*(i-1))
           
        ELSEIF(axismove == 4) THEN

           rxyz_lmp(k,1) = rxyz_lmp(j,1)+nxvec*movedist*REAL(mult*(i&
                &-1))
           rxyz_lmp(k,2) = rxyz_lmp(j,2)+nyvec*movedist*REAL(mult*(i&
                &-1))
           rxyz_lmp(k,3) = rxyz_lmp(j,3)+nzvec*movedist*REAL(mult*(i&
                &-1))

        ELSE
           
           PRINT *, "Unknown axis number", axismove
           STOP
           
        END IF
           
     END DO

  END DO

  IF(axismove == 4) THEN

     ringxcom = 0.0; ringycom = 0.0; ringzcom = 0.0
     totmass = 0.0
     DO i = 1,n_rings
        
        DO j = 1, natomsperpoly
           
           k = natomsperpoly*(i-1) + j
           
           atype = aidvals(j,3)
           ringxcom(i) = ringxcom(i) + rxyz_lmp(k,1)*Masses(atype,1)
           ringycom(i) = ringycom(i) + rxyz_lmp(k,2)*Masses(atype,1)
           ringzcom(i) = ringzcom(i) + rxyz_lmp(k,3)*Masses(atype,1)
           totmass = totmass + Masses(atype,1)
           
        END DO

        ringxcom(i) = ringxcom(i)/totmass
        ringycom(i) = ringycom(i)/totmass
        ringzcom(i) = ringzcom(i)/totmass

        PRINT *, "center of mass"
        PRINT *, i, ringxcom(i), ringycom(i), ringzcom(i)

     END DO
           

  END IF

  IF(alternate == 2) THEN ! Alternating Config
     !Generating Flexible Chains
     i = n_rings + 1; u = 0

     DO WHILE(i .LE. ntotchains)

        k = (i-1)*natomsperpoly + 1; kinit = k
        
        IF(axismove == 1) THEN
           
           rxyz_lmp(k,1) = rxring + movedist*REAL(2*u+1)
           rxyz_lmp(k,2) = ryring + RAN1(X)*rgring
           rxyz_lmp(k,3) = rzring + RAN1(X)*rgring
           
        ELSEIF(axismove == 2) THEN
           
           rxyz_lmp(k,1) = rxring + RAN1(X)*rgring
           rxyz_lmp(k,2) = ryring + movedist*REAL(2*u+1)
           rxyz_lmp(k,3) = rzring + RAN1(X)*rgring
           
        ELSEIF(axismove == 3) THEN
           
           rxyz_lmp(k,1) = rxring + RAN1(X)*rgring
           rxyz_lmp(k,2) = ryring + RAN1(X)*rgring
           rxyz_lmp(k,3) = rzring + movedist*REAL(2*u+1)
           
        ELSEIF(axismove == 4) THEN

           rxyz_lmp(k,1) = rxring + movedist*REAL(2*u+1)*nxvec
           rxyz_lmp(k,2) = ryring + movedist*REAL(2*u+1)*nyvec
           rxyz_lmp(k,3) = rzring + movedist*REAL(2*u+1)*nzvec

        END IF
        
        k = k + 1
        
        theta       = math_pi*RAN1(X)
        phi         = 2*math_pi*RAN1(X)
        
        rxyz_lmp(k,1) = rxyz_lmp(k-1,1) + r0init*sin(theta)*cos(phi)
        rxyz_lmp(k,2) = rxyz_lmp(k-1,2) + r0init*sin(theta)*sin(phi)
        rxyz_lmp(k,3) = rxyz_lmp(k-1,3) + r0init*cos(theta)
        
        k = k + 1
        
        DO WHILE (k .LE. i*natomsperpoly)
           
           theta       = math_pi*RAN1(X)
           phi         = 2*math_pi*RAN1(X)
           
           rxyz_lmp(k,1) = rxyz_lmp(k-1,1) + r0init*sin(theta)*cos(phi)
           rxyz_lmp(k,2) = rxyz_lmp(k-1,2) + r0init*sin(theta)*sin(phi)
           rxyz_lmp(k,3) = rxyz_lmp(k-1,3) + r0init*cos(theta)
           
           rx2ij = rxyz_lmp(k-2,1) - rxyz_lmp(k,1)
           ry2ij = rxyz_lmp(k-2,2) - rxyz_lmp(k,2)
           rz2ij = rxyz_lmp(k-2,3) - rxyz_lmp(k,3)
           ri2j  = rx2ij**2 + ry2ij**2 + rz2ij**2
           
           IF(ri2j > 1.0404) THEN
              
              k = k + 1
              
           END IF
           
        END DO
        
        k = i*natomsperpoly
        rxe2e = rxyz_lmp(k,1) - rxyz_lmp(kinit,1)
        rye2e = rxyz_lmp(k,2) - rxyz_lmp(kinit,2)
        rze2e = rxyz_lmp(k,3) - rxyz_lmp(kinit,3)
        
        re2e = rxe2e**2 + rye2e**2 + rze2e**2
        
        IF(sqrt(re2e) .LT. 0.5*REAL(natomsperpoly)*r0init) THEN
           
           i = i + 1
           u = u + 1
           
        END IF
        
     END DO
        
  ELSE ! Block wise arrangement
     !Generating Flexible Chains     
     i = n_rings + 1; u = 0

     DO WHILE(i .LE. ntotchains)

        k = (i-1)*natomsperpoly + 1; kinit = k
        kref1 = n_rings*natomsperpoly + 1

        IF(axismove == 1) THEN
           
           rxyz_lmp(k,1) = rxring + (n_rings+u)*movedist
           rxyz_lmp(k,2) = ryring + RAN1(X)*rgring
           rxyz_lmp(k,3) = rzring + RAN1(X)*rgring
           
        ELSEIF(axismove == 2) THEN

           rxyz_lmp(k,1) = rxring + RAN1(X)*rgring
           rxyz_lmp(k,2) = ryring + (n_rings+u)*movedist
           rxyz_lmp(k,3) = rzring + RAN1(X)*rgring

        ELSEIF(axismove == 3) THEN
           
           rxyz_lmp(k,1) = rxring + RAN1(X)*rgring
           rxyz_lmp(k,2) = ryring + RAN1(X)*rgring
           rxyz_lmp(k,3) = rzring + (n_rings+u)*movedist
           
        ELSEIF(axismove == 4) THEN

           rxyz_lmp(k,1) = rxring + (n_rings+u)*movedist*nxvec
           rxyz_lmp(k,2) = ryring + (n_rings+u)*movedist*nyvec
           rxyz_lmp(k,3) = rzring + (n_rings+u)*movedist*nzvec
        
        END IF

        k = k + 1

        theta       = math_pi*RAN1(X)
        phi         = 2*math_pi*RAN1(X)
        
        rxyz_lmp(k,1) = rxyz_lmp(k-1,1) + r0init*sin(theta)*cos(phi)
        rxyz_lmp(k,2) = rxyz_lmp(k-1,2) + r0init*sin(theta)*sin(phi)
        rxyz_lmp(k,3) = rxyz_lmp(k-1,3) + r0init*cos(theta)
        
        k = k + 1
        
        DO WHILE (k .LE. i*natomsperpoly)
           
           theta       = math_pi*RAN1(X)
           phi         = 2*math_pi*RAN1(X)
           
           rxyz_lmp(k,1) = rxyz_lmp(k-1,1) + r0init*sin(theta)*cos(phi)
           rxyz_lmp(k,2) = rxyz_lmp(k-1,2) + r0init*sin(theta)*sin(phi)
           rxyz_lmp(k,3) = rxyz_lmp(k-1,3) + r0init*cos(theta)
           
           rx2ij = rxyz_lmp(k-2,1) - rxyz_lmp(k,1)
           ry2ij = rxyz_lmp(k-2,2) - rxyz_lmp(k,2)
           rz2ij = rxyz_lmp(k-2,3) - rxyz_lmp(k,3)
           ri2j  = rx2ij**2 + ry2ij**2 + rz2ij**2
           
           IF(ri2j > 1.0404) THEN
              
              k = k + 1
              
           END IF
           
        END DO
        
        k = i*natomsperpoly
        rxe2e = rxyz_lmp(k,1) - rxyz_lmp(kinit,1)
        rye2e = rxyz_lmp(k,2) - rxyz_lmp(kinit,2)
        rze2e = rxyz_lmp(k,3) - rxyz_lmp(kinit,3)
        
        re2e = rxe2e**2 + rye2e**2 + rze2e**2
        
        IF(sqrt(re2e) .LT. 0.5*REAL(natomsperpoly)*r0init) THEN
           
           i = i + 1
           u = u + 1
           
        END IF

     END DO

     !Now shift such that the middle monomer is along the line joining
     !the two rings - to make it symmetrical

     DO i = n_rings+1,ntotchains
        
        kinit = (i-1)*natomsperpoly + 1
        kmid = kinit+INT(0.5*natomsperpoly)
        dispx = rxyz_lmp(kinit,1) - rxyz_lmp(kmid,1)
        dispy = rxyz_lmp(kinit,2) - rxyz_lmp(kmid,2)
        dispz = rxyz_lmp(kinit,3) - rxyz_lmp(kmid,3)

        DO j = 1,natomsperpoly

           k = kinit + j - 1
           
           rxyz_lmp(k,1) = rxyz_lmp(k,1) + dispx
           rxyz_lmp(k,2) = rxyz_lmp(k,2) + dispy
           rxyz_lmp(k,3) = rxyz_lmp(k,3) + dispz

        END DO

     END DO

  END IF
     
END SUBROUTINE GENERATE_BLOCK_POS

!--------------------------------------------------------------------

SUBROUTINE WRAPALL()

  USE PARAMS_REP
  IMPLICIT NONE

  INTEGER :: i

  DO i = 1,ntotatoms
     
     rxyz_lmp(i,1) = rxyz_lmp(i,1) - boxxl*floor(rxyz_lmp(i,1)/boxxl)
     rxyz_lmp(i,2) = rxyz_lmp(i,2) - boxyl*floor(rxyz_lmp(i,2)/boxyl)
     rxyz_lmp(i,3) = rxyz_lmp(i,3) - boxzl*floor(rxyz_lmp(i,3)/boxzl)
          
  END DO

END SUBROUTINE WRAPALL

!--------------------------------------------------------------------

SUBROUTINE GENERATE_IMAGE_FLAGS()

  USE PARAMS_REP
  IMPLICIT NONE

  INTEGER :: i,j,k
  REAL :: rx, ry, rz
! Make image flags

  DO i = 1,ntotchains

     k = natomsperpoly*(i-1) + 1
     
     ixyz_lmp(k,1) = 0
     ixyz_lmp(k,2) = 0
     ixyz_lmp(k,3) = 0
     
     DO j = 1, natomsperpoly-1
        
        k = natomsperpoly*(i-1) + j

        rx = rxyz_lmp(k,1) - rxyz_lmp(k+1,1)
        ry = rxyz_lmp(k,2) - rxyz_lmp(k+1,2)
        rz = rxyz_lmp(k,3) - rxyz_lmp(k+1,3)

        CALL IMGFLAGS(rx,ixyz_lmp(k,1),boxxl,ixyz_lmp(k+1,1))
        CALL IMGFLAGS(ry,ixyz_lmp(k,2),boxyl,ixyz_lmp(k+1,2))
        CALL IMGFLAGS(rz,ixyz_lmp(k,3),boxzl,ixyz_lmp(k+1,3))
        
     END DO
     
  END DO


END SUBROUTINE GENERATE_IMAGE_FLAGS

!--------------------------------------------------------------------

SUBROUTINE IMGFLAGS(dist,img,boxl,nx)

  USE PARAMS_REP
  IMPLICIT NONE
  
  REAL, INTENT(IN) :: dist,boxl
  INTEGER, INTENT(IN) :: img
  INTEGER, INTENT(OUT) :: nx
  
  IF(dist > boxl/2) THEN
     
     nx = img + 1
     
  ELSEIF(dist < -boxl/2) THEN
     
     nx = img - 1
     
  ELSE

     nx = img
     
  END IF
  
END SUBROUTINE IMGFLAGS

!--------------------------------------------------------------------

SUBROUTINE GENERATE_ALL_TOPO()

  USE PARAMS_REP
  IMPLICIT NONE

  INTEGER :: i, j, k
  INTEGER :: idval, adder

! Make Atom Information

  DO i = 1,ntotchains
     
     DO j = 1, natomsperpoly

        adder = natomsperpoly*(i-1) 
        idval = j + adder
        aidvals(idval,1) = idval
        aidvals(idval,2) = i
        aidvals(idval,3) = aidvals(j,3)

     END DO

  END DO

! Make Bond Topology

  DO i = 1,ntotchains
     
     DO j = 1, nbondsperpoly
        
        adder = natomsperpoly*(i-1) 
        idval = nbondsperpoly*(i-1) + j
        bond_lmp(idval,1) = idval
        bond_lmp(idval,2) = bond_lmp(j,2)
        bond_lmp(idval,3) = bond_lmp(j,3) + adder
        bond_lmp(idval,4) = bond_lmp(j,4) + adder

     END DO

  END DO

! Make Angle Topology

  DO i = 1,ntotchains
     
     DO j = 1, nanglsperpoly
        
        adder = natomsperpoly*(i-1)
        idval = nanglsperpoly*(i-1) + j
        angl_lmp(idval,1) = idval
        angl_lmp(idval,2) = angl_lmp(j,2)
        angl_lmp(idval,3) = angl_lmp(j,3) + adder
        angl_lmp(idval,4) = angl_lmp(j,4) + adder
        angl_lmp(idval,5) = angl_lmp(j,5) + adder

     END DO

  END DO

! Make Dihedral Topology

  DO i = 1,ntotchains
     
     DO j = 1, ndihdsperpoly
        
        adder = natomsperpoly*(i-1)
        idval = ndihdsperpoly*(i-1) + j
        dihd_lmp(idval,1) = idval
        dihd_lmp(idval,2) = dihd_lmp(j,2)
        dihd_lmp(idval,3) = dihd_lmp(j,3) + adder
        dihd_lmp(idval,4) = dihd_lmp(j,4) + adder
        dihd_lmp(idval,5) = dihd_lmp(j,5) + adder
        dihd_lmp(idval,6) = dihd_lmp(j,6) + adder

     END DO

  END DO

END SUBROUTINE GENERATE_ALL_TOPO

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_RG_RING()

  USE PARAMS_REP
  IMPLICIT NONE

  INTEGER :: i,j,molid,atype
  REAL :: rgxx, rgyy, rgzz, rgsq
  REAL :: rxcm, rycm, rzcm, totmass

  rgxx = 0.0; rgyy =0.0; rgzz = 0.0; rgsq = 0.0; totmass = 0.0
  rxcm = 0.0; rycm =0.0; rzcm = 0.0
  
  PRINT *, "Atoms/molecule: ", natomsperpoly
  
  DO i = 1,natomsperpoly

     molid = aidvals(i,2)
     atype = aidvals(i,3)
     totmass = totmass + masses(atype,1)

     rxcm = rxcm + rxyz_lmp(i,1)*masses(atype,1)
     rycm = rycm + rxyz_lmp(i,2)*masses(atype,1)
     rzcm = rzcm + rxyz_lmp(i,3)*masses(atype,1)

  END DO

  rxcm = rxcm/totmass
  rycm = rycm/totmass
  rzcm = rzcm/totmass

  DO i = 1,natomsperpoly

     molid = aidvals(i,2)
     atype = aidvals(i,3)

     rgxx = rgxx + (masses(atype,1)*((rxyz_lmp(i,1)-rxcm))**2)
     rgyy = rgyy + (masses(atype,1)*((rxyz_lmp(i,2)-rycm))**2)
     rgzz = rgzz + (masses(atype,1)*((rxyz_lmp(i,3)-rzcm))**2)

     rgsq = rgsq + masses(atype,1)*((rxyz_lmp(i,1)-rxcm)**2 +&
          & (rxyz_lmp(i,2)-rycm)**2 + (rxyz_lmp(i,3)-rzcm)**2)

  END DO

  rgsq = rgsq/totmass
  rgring = sqrt(rgsq)
  rxring = rxcm; ryring = rycm; rzring = rzcm

  PRINT *, "Radius of gyration of base ring, ", rgring
  PRINT *, "Radius of gyration of base ring (nm), ", 0.515*rgring
     
END SUBROUTINE COMPUTE_RG_RING

!--------------------------------------------------------------------

SUBROUTINE WRITELAMMPSDATA()

  USE PARAMS_REP
  IMPLICIT NONE

  INTEGER :: ierr,i
  INTEGER, PARAMETER :: nullval = 0
  REAL*8 :: csum
  CHARACTER(LEN=256) :: lmpdatafile
  CHARACTER(LEN=5)   :: typchar
  CHARACTER(LEN=4)   :: distchar

  print *, "movedist: ", moveperc
  WRITE(distchar,"(F4.2)") moveperc

  lmpdatafile = "replicated_"//trim(adjustl(data_fname))//"_"&
       &//distchar

  OPEN(unit = outfile,file=trim(lmpdatafile),action="write"&
       &,status="replace",iostat=ierr)

  IF(ierr /= 0) STOP "LAMMPS data file not found"

  WRITE(logout,*) "Generating PDB File: ", trim(lmpdatafile)
2000 FORMAT(2(F15.8,1X),2X,A)

  WRITE(outfile,*) "Input Configuration for Methylcellulose"
  WRITE(outfile,*)

  WRITE(outfile,"(I0,1X,A5)") ntotatoms, "atoms"
  WRITE(outfile,"(I0,1X,A5)") ntotbonds, "bonds"
  WRITE(outfile,"(I0,1X,A6)") ntotangls, "angles"
  WRITE(outfile,"(I0,1X,A9)") ntotdihds, "dihedrals"
  WRITE(outfile,"(I0,1X,A9)") nullval,"impropers"

  WRITE(outfile,"(I0,1X,A10)") natomtypes,"atom types"
  WRITE(outfile,"(I0,1X,A10)") nbondtypes,"bond types"
  WRITE(outfile,"(I0,1X,A11)") nangltypes,"angle types"
  WRITE(outfile,"(I0,1X,A14)") ndihdtypes,"dihedral types"
  WRITE(outfile,"(I0,1X,A14)") nullval,"improper types"

  
  WRITE (outfile,'(2(F15.8,1X),2X,A)') 0.00, boxxl, "xlo xhi"
  WRITE (outfile,'(2(F15.8,1X),2X,A)') 0.00, boxyl, "ylo yhi"
  WRITE (outfile,'(2(F15.8,1X),2X,A)') 0.00, boxzl, "zlo zhi"

  WRITE(outfile,*)

  WRITE(outfile,*) " Masses"

  WRITE(outfile,*)

  
  DO i = 1, natomtypes

     WRITE(typchar,'(I0)') i
     WRITE(outfile,"(I0,1X,F14.9)") i,masses(i,1)

  END DO

  WRITE(outfile,*)

  WRITE(outfile,*) " Atoms"

  WRITE(outfile,*)

  csum = 0.0

  DO i = 1,ntotatoms
     
     WRITE(typchar,'(I0)') aidvals(i,3)
   
     WRITE(outfile,'(3(I0,1X),3(F15.9,1X),3(I0,1X))') aidvals(i,1)&
          &,aidvals(i,2),aidvals(i,3),rxyz_lmp(i,1), rxyz_lmp(i,2),&
          & rxyz_lmp(i,3),ixyz_lmp(i,1),ixyz_lmp(i,2),ixyz_lmp(i,3)
     
  END DO
     

  WRITE(outfile,*)
  WRITE(outfile,*) " Bonds"
  WRITE(outfile,*)

  DO i = 1,ntotbonds
     
     WRITE(outfile,'(4(I0,1X))') bond_lmp(i,1),bond_lmp(i,2)&
          &,bond_lmp(i,3),bond_lmp(i,4)
     
  END DO

  WRITE(outfile,*)
  WRITE(outfile,*) " Angles"
  WRITE(outfile,*)


  DO i = 1,ntotangls
     
     WRITE(outfile,'(5(I0,1X))') angl_lmp(i,1)&
          &,angl_lmp(i,2),angl_lmp(i,3),angl_lmp(i,4),angl_lmp(i,5)
     
  END DO

  WRITE(outfile,*)
  WRITE(outfile,*) " Dihedrals"
  WRITE(outfile,*)


  DO i = 1,ntotdihds
     
     WRITE(outfile,'(6(I0,1X))') dihd_lmp(i,1),dihd_lmp(i,2)&
          &,dihd_lmp(i,3),dihd_lmp(i,4),dihd_lmp(i,5),dihd_lmp(i,6)
     
  END DO


  CLOSE(outfile)


END SUBROUTINE WRITELAMMPSDATA

!--------------------------------------------------------------------

SUBROUTINE ALLOCATE_ARRAYS()

  USE PARAMS_REP
  IMPLICIT NONE

  INTEGER :: AllocateStatus

  ALLOCATE(aidvals(ntotatoms,3),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate angl_lmp"
  ALLOCATE(rxyz_lmp(ntotatoms,3),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate rxyz_lmp"
  ALLOCATE(ixyz_lmp(ntotatoms,3),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate ixyz_lmp"
  ALLOCATE(charge_lmp(ntotatoms,1),stat = AllocateStatus)
  IF(AllocateStatus/=0) STOP "did not allocate charge_lmp"


  IF(ntotbonds /= 0) THEN
     ALLOCATE(bond_lmp(ntotbonds,4),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate bond_lmp"
     ALLOCATE(bondtypearr(nbondtypes,3),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate bondtypearr"
  ELSE
     ALLOCATE(bondtypearr(1,1),stat = AllocateStatus)
     DEALLOCATE(bondtypearr)
     ALLOCATE(bond_lmp(1,1),stat = AllocateStatus)
     DEALLOCATE(angl_lmp)
  END IF
  

  IF(ntotangls /= 0) THEN
     ALLOCATE(angl_lmp(ntotangls,5),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate angl_lmp"
     ALLOCATE(angltypearr(nangltypes,3),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate angltypearr"
  ELSE
     ALLOCATE(angltypearr(1,1),stat = AllocateStatus)
     DEALLOCATE(angltypearr)
     ALLOCATE(angl_lmp(1,1),stat = AllocateStatus)
     DEALLOCATE(angl_lmp)
  END IF
     
  IF(ntotdihds /= 0) THEN
     ALLOCATE(dihd_lmp(ntotdihds,6),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate dihd_lmp"
     ALLOCATE(dhdtypearr(ndihdtypes,5),stat = AllocateStatus)
     IF(AllocateStatus/=0) STOP "did not allocate dhdtypearr"
  ELSE
     ALLOCATE(dhdtypearr(1,1),stat = AllocateStatus)
     DEALLOCATE(dihd_lmp)
     ALLOCATE(dihd_lmp(1,1),stat = AllocateStatus)
     DEALLOCATE(dihd_lmp)
  END IF


END SUBROUTINE ALLOCATE_ARRAYS

!--------------------------------------------------------------------

SUBROUTINE DEALLOCATEARRAYS()

  USE PARAMS_REP
  IMPLICIT NONE

  DEALLOCATE(rxyz_lmp)
  DEALLOCATE(ixyz_lmp)
  DEALLOCATE(aidvals)
  DEALLOCATE(bond_lmp)
  DEALLOCATE(angl_lmp)
  DEALLOCATE(dihd_lmp)

END SUBROUTINE DEALLOCATEARRAYS

!--------------------------------------------------------------------
