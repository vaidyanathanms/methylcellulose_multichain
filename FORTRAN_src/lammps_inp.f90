!-----To generate input parameters for LAMMPS CG-MethylCellulose-----
!----- Parameter File: lmp_params.f90--------------------------------
!------Version: Jan-25-2018------------------------------------------
!--------------------------------------------------------------------

PROGRAM LAMMPSINP

  USE PARAMS

  IMPLICIT NONE
  
  LOGICAL :: input_coor = .false.
  REAL :: bondl, bondlsq
  INTEGER :: ierror,narg

  bondlsq = 0

  CALL SYSTEM_CLOCK(S)

  narg = IARGC()
  
  IF(narg == 1) THEN

     OPEN (unit = outfile, file = "lmp_input.txt", status ="replace",&
          & action="write",iostat=ierror)
     
     IF(ierror /= 0) STOP "Cannot open lmp_input.txt"
        
     CALL CREATEFILE(narg)
     CALL COMPUTE_BOX()
     CALL INPCOR(input_coor)
     IF(grafts == 0) THEN
        CALL CREATE_ATYPE()
     ELSE
        CALL CREATE_ATYPE()
        CALL CREATE_ATYPE_DEF()
     END IF
     CALL BOND_LEN(bondl, bondlsq)
     
     WRITE (outfile,*) "Total particles", totpart
     WRITE (outfile,*) "N = ", N, "M = ", M
        
     CALL LMP_COORD()

     CLOSE (unit = outfile)

  ELSE

     PRINT *, "Please specify the extension"

  END IF

END PROGRAM LAMMPSINP

!--------------------------------------------------------------------

SUBROUTINE LMP_COORD()

  USE PARAMS
  
  IMPLICIT NONE
  
  INTEGER :: i,j, k, ierror
  INTEGER ::  bondid, atomid, molid, anglid, dihdid, DSmon,numatomtot
  REAL :: rx, ry, rz
  REAL :: ranval, massval
  
  i = 1
  
  atomid = 0

20 FORMAT(5X,I0,2X,A)
22 FORMAT(5X,I0,2X,A)
24 FORMAT(5X,I0,2X,F14.6,2X,A)
  
  OPEN (unit=10, file = datafile, status="replace",action=&
       &"write",iostat = ierror)
  
  IF(ierror /= 0) STOP "Failed to open datafile"
     
  WRITE (10,'(A36,1X,A5,1X,A18,1X,F6.3)') "Data for CG-MC simulations &
       &using ",ext, "potential and DS: ", DS_MC 
  WRITE (10,*) 
  WRITE (10,20) totpart, "atoms"

  IF(numbondtypes /= 0) THEN
     WRITE (10,20) N*(M-1), "bonds"
  ELSE
     WRITE (10,20) 0, "bonds"
  END IF

  IF(numangltypes /= 0) THEN
     WRITE (10,20) N*(M-2), "angles"
  ELSE
     WRITE (10,20) 0, "angles"
  END IF

  IF(numdihdtypes /= 0) THEN
     WRITE (10,20) N*(M-3), "dihedrals"
  ELSE
     WRITE (10,20) 0, "dihedrals"
  END IF

  WRITE (10,20) 0, "impropers"
  IF(grafts == 0) THEN
     WRITE (10,20) numatomtypes, "atom types"
     numatomtot = numatomtypes
  ELSE
     WRITE (10,20) 1 + numatomtypes, "atom types"
     numatomtot = 1+numatomtypes
  END IF
  WRITE (10,20) numbondtypes, "bond types"
  WRITE (10,22) numangltypes, "angle types"
  WRITE (10,22) numdihdtypes, "dihedral types"
  WRITE (10,22) 0, "improper types"

  WRITE (10,*)
  WRITE (10,24) 0, boxl_x, "xlo xhi"
  WRITE (10,24) 0, boxl_y, "ylo yhi"
  WRITE (10,24) 0, boxl_z, "zlo zhi"
  WRITE (10,*)
  WRITE (10,*) "Masses"
  WRITE (10,*)
  
!!$Fetching image information
  
  
  ix(1,1) = 0
  ix(1,2) = 0
  ix(1,3) = 0
  
  DO i = 1,N
     
     ix(i,1) = 0
     ix(i,2) = 0
     ix(i,3) = 0
     
     DO j = 1,M-1
        
        k = (j-1)*3
        rx = rxyz(i,k+1) - rxyz(i,k+4)
        ry = rxyz(i,k+2) - rxyz(i,k+5)
        rz = rxyz(i,k+3) - rxyz(i,k+6)
        
        ix(i,k+4) = IMGFLAGS(rx,ix(i,k+1),boxl_x)
        ix(i,k+5) = IMGFLAGS(ry,ix(i,k+2),boxl_y)
        ix(i,k+6) = IMGFLAGS(rz,ix(i,k+3),boxl_z)
        
     END DO
     
  END DO
  

  ! Writing Masses

  DO i = 1,numatomtot

     IF(i == 1) THEN

        DSmon = 0

     ELSEIF(i .LE. 4) THEN

        DSmon = 1

     ELSEIF(i .LE. 7) THEN

        DSmon = 2

     ELSEIF(i == 8) THEN

        DSmon = 3

     ELSEIF(i== 9) THEN

        DSmon = 1.5

     END IF

     massval = 12.0106*(6+DSmon)+15.9994*5+1.008*(10-2.0*DSmon)

     WRITE(10,'(I0,1X,F14.8)') i, massval/massval !massval/basemass

  END DO
  
  ! Writing atomic corrdinates
  
  WRITE (10,*) 
  WRITE (10,*) "Atoms"
  WRITE (10,*)

100 FORMAT(I7,1X,I5,1X,I1,1X,F9.6,1X,F9.6,1X,F9.6,1X,&
         &I2,1X,I2,1X,I2)
  
200 FORMAT(I7,1X,I5,1X,I1,1X,F9.6,1X,F9.6,1X,F9.6,1X,&
         &I2,1X,I2,1X,I2)
  
  DO i = 1,N
     
     DO j = 1,M
        
        k = (j-1)*3
        atomid = atomid + 1

        rx = rxyz(i,k+1) + boxl_x*ix(i,k+1)
        ry = rxyz(i,k+2) + boxl_y*ix(i,k+2)
        rz = rxyz(i,k+3) + boxl_z*ix(i,k+3)

        WRITE(10,'(3(I0,1X),3(F14.6,1X),3(I0,1X))') atomid, i,&
             & atype(i,j), rxyz(i,k+1), rxyz(i,k+2), rxyz(i,k+3),ix(i&
             &,k+1), ix(i,k+2), ix(i,k+3)
        
     END DO
     
  END DO

  IF(numbondtypes /= 0) THEN

     ! Writing Bond Details  
     
     bondid = 0
     atomid = 0
     WRITE (10,*)
     WRITE (10,*) "Bonds"
     WRITE (10,*)
     
     DO i = 1,N
        
        DO j = 1,M-1
           
           bondid = bondid + 1
           atomid = atomid + 1
           
           WRITE(10,'(4(I0,2X))') bondid, bondtype, atomid, atomid+1
           
        END DO
        
        atomid = atomid + 1
        
     END DO

  END IF

  IF(numangltypes /= 0) THEN

     ! Writing Angle Details

     anglid = 0
     atomid = 0
     WRITE (10,*)
     WRITE (10,*) "Angles"
     WRITE (10,*)
     
     DO i = 1,N
        
        DO j = 1,M-2
           
           anglid = anglid + 1
           atomid = atomid + 1
           
           WRITE(10,'(5(I0,2X))') anglid, angltype, atomid, atomid+1,&
                & atomid+2
           
        END DO
        
        atomid = atomid + 2
        
     END DO

  END IF

  IF(numdihdtypes /= 0) THEN

     ! Writing Dihedral Details
     
     dihdid = 0
     atomid = 0
     WRITE (10,*)
     WRITE (10,*) "Dihedrals"
     WRITE (10,*)
     
     DO i = 1,N
        
        DO j = 1,M-3
           
           dihdid = dihdid + 1
           atomid = atomid + 1
           
           WRITE(10,'(6(I0,2X))') dihdid, dihdtype, atomid, atomid+1,&
                & atomid+2,atomid+3
           
        END DO
        
        atomid = atomid + 3
        
     END DO
     
     CLOSE(unit = 10)

  END IF
  
CONTAINS
  
  INTEGER FUNCTION IMGFLAGS(dist, img,boxl)
    
    USE PARAMS
    
    IMPLICIT NONE
    
    REAL, INTENT(IN) :: dist,boxl
    INTEGER, INTENT(IN) :: img
    INTEGER :: nx
    
    IF(dist > boxl/2) THEN
       
       nx = img + 1
       
    ELSEIF(dist < -boxl/2) THEN
       
       nx = img - 1
       
    ELSE
       
       nx = img
       
    END IF
    
    IMGFLAGS = nx
    
  END FUNCTION IMGFLAGS
  
END SUBROUTINE LMP_COORD

!--------------------------------------------------------------------

SUBROUTINE COMPUTE_BOX()

  USE PARAMS

  IMPLICIT NONE

  PRINT *, "Target density: ", density

  volbox = totpart/density
  
  boxl_x = volbox**(1.0/3.0)
  boxl_y = volbox**(1.0/3.0)
  boxl_z = volbox**(1.0/3.0)

END SUBROUTINE COMPUTE_BOX

!--------------------------------------------------------------------

SUBROUTINE INPCOR(input)
  
  USE PARAMS

  IMPLICIT NONE
  
  INTEGER :: i,j,k,u,v,ierror
  LOGICAL,INTENT(IN) :: input 
  REAL, PARAMETER :: r0init  = 0.97
  REAL, PARAMETER :: r0sq3   = r0init/sqrt(3.0)
  REAL, PARAMETER :: rmaxsq  = r0init*r0init
  REAL, PARAMETER :: math_pi = 3.14159265359
  REAL :: theta, phi
  REAL :: rx2ij,ry2ij,rz2ij,ri2j,rxe2e,rye2e,rze2e,re2e
  CALL RAN_INIT(S,X)
  IF(input==.true.) THEN
     
     OPEN (unit=20, file="newdata.txt", status="old", action =&
          &"read",iostat=ierror)
     
     IF (ierror /= 0) THEN
        
        PRINT*, "Failed to open inpcor.txt"
     
     END IF

     DO i = 1,9
        
        READ(20,*) 

     END DO

     DO i = 1,N
        
        DO j = 1,M
           
           k = (j-1)*3
           
           READ(20,*) u, v, rxyz(i,k+1), rxyz(i,k+2), rxyz(i,k+3)
           
        END DO
        
     END DO
        
     CLOSE (unit = 20)
        
  ELSE IF(stretched == 0) THEN
     
     PRINT *, "Random Initial Configuration : NRRW"
     OPEN (unit=21, file="theta_phi.txt", status="replace", action =&
          &"write",iostat=ierror)
    
     OPEN (unit=28, file="re2e.txt", status="replace", action =&
          &"write",iostat=ierror)

     i = 1
     DO WHILE(i <= N)
        
        k = 1
        
        rxyz(i,k)   = RAN1(X)*insidebox
        rxyz(i,k+1) = RAN1(X)*insidebox
        rxyz(i,k+2) = RAN1(X)*insidebox

        k = k + 3

        theta       = math_pi*RAN1(X)
        phi         = 2*math_pi*RAN1(X)
        WRITE(21,*), i,"2",theta,phi
        rxyz(i,k)   = rxyz(i,k-3) + r0init*sin(theta)*cos(phi)
        rxyz(i,k+1) = rxyz(i,k-2) + r0init*sin(theta)*sin(phi)
        rxyz(i,k+2) = rxyz(i,k-1) + r0init*cos(theta)
        
        k = k + 3
        j = 3
        WRITE(21,*), i,j,theta,phi
        DO WHILE (j <= M)
           
           theta       = math_pi*RAN1(X)
           phi         = 2*math_pi*RAN1(X)
           rxyz(i,k)   = rxyz(i,k-3) + r0init*sin(theta)*cos(phi)
           rxyz(i,k+1) = rxyz(i,k-2) + r0init*sin(theta)*sin(phi)
           rxyz(i,k+2) = rxyz(i,k-1) + r0init*cos(theta)
           
           rx2ij = rxyz(i,k) - rxyz(i,k-6)
           ry2ij = rxyz(i,k+1) - rxyz(i,k-5)
           rz2ij = rxyz(i,k+2) - rxyz(i,k-4)
           ri2j  = rx2ij**2 + ry2ij**2 + rz2ij**2
           
           IF(ri2j > 1.0404) THEN
              
              k = k + 3
              j = j + 1
              WRITE(21,*), i,j,theta,phi

           END IF

        END DO

        k = (M-1)*3
        rxe2e = rxyz(i,k+1) - rxyz(i,1)
        rye2e = rxyz(i,k+2) - rxyz(i,2)
        rze2e = rxyz(i,k+3) - rxyz(i,3)
        
        re2e = rxe2e**2 + rye2e**2 + rze2e**2
        
        
        IF(sqrt(re2e) .LT. 0.5*REAL(M)*r0init) THEN
           !Maximum stretch possible is when ree = Mb
           !So we can allow upto 50% of max stretch
           !This condition needs to be checked if no of chains are
           !very less
           WRITE(28,*) i, re2e,M*r0init
           i = i + 1
           

        END IF

     END DO
     
     CLOSE(unit = 28)
     
  ELSE

     PRINT *, "Random Initial Configuration - Stretched"

     DO i = 1, N
        
        
        k = 1
        
        rxyz(i,k)   = RAN1(X)*boxl_x
        rxyz(i,k+1) = RAN1(X)*boxl_y
        rxyz(i,k+2) = RAN1(X)*boxl_z
        
        k = k + 3

        DO j = 2,M

           rxyz(i,k)   = rxyz(i,k-3) + r0sq3
           rxyz(i,k+1) = rxyz(i,k-2) + r0sq3
           rxyz(i,k+2) = rxyz(i,k-1) + r0sq3
           
           k = k+3
        END DO

     END DO


     
  END IF
  ! PBC
  
  DO i = 1,N
     
     DO j = 1,M 
        
        k = (j-1)*3
        
        rxyz(i,k+1) = rxyz(i,k+1) - boxl_x*floor(rxyz(i,k+1)/boxl_x)
        rxyz(i,k+2) = rxyz(i,k+2) - boxl_y*floor(rxyz(i,k+2)/boxl_y)
        rxyz(i,k+3) = rxyz(i,k+3) - boxl_z*floor(rxyz(i,k+3)/boxl_z)
        
     END DO
     
  END DO
  
  CLOSE(unit =21)

END SUBROUTINE INPCOR

!--------------------------------------------------------------------------------------


SUBROUTINE BOND_LEN(bondl, bondlsq)

  USE PARAMS

  IMPLICIT NONE

  INTEGER :: i, j, k
  REAL :: bondl, bondlsq
  REAL :: rx, ry, rz
  REAL :: end2end

  bondlsq = 0
!!$  PRINT *, "entered bondl loop "

  OPEN (unit = 14, file = "bond_lammps.txt", status = "replace", action= &
       &"write")

  OPEN (unit = 15, file = "length_check.txt", status = "replace", action= &
       &"write")


  DO i = 1,N

     bondlsq = 0.0
     DO j = 1,M-1

        k = (j-1)*3

        rx = rxyz(i,k+1) - rxyz(i,k+4)
        ry = rxyz(i,k+2) - rxyz(i,k+5)
        rz = rxyz(i,k+3) - rxyz(i,k+6)

        rx = rx - boxl_x*ANINT(rx/boxl_x)
        ry = ry - boxl_y*ANINT(ry/boxl_y)
        rz = rz - boxl_z*ANINT(rz/boxl_z)

        bondl   = rx**2 + ry**2 + rz**2
        bondlsq = bondlsq + bondl

        WRITE(14,*) i, j, j+1, bondl

     END DO

     rx = rxyz(i,k+1) - rxyz(i,1)
     ry = rxyz(i,k+2) - rxyz(i,2)
     rz = rxyz(i,k+3) - rxyz(i,3)

     rx = rx - boxl_x*ANINT(rx/boxl_x)
     ry = ry - boxl_y*ANINT(ry/boxl_y)
     rz = rz - boxl_z*ANINT(rz/boxl_z)

     end2end = rx**2 + ry**2 + rz**2

     WRITE(15,*) i, bondlsq

  END DO

  bondlsq = bondlsq/(N*(M-1))

!!$  PRINT *, bondlsq

  CLOSE (unit = 14)
  CLOSE (unit = 15)

END SUBROUTINE BOND_LEN
        
!--------------------------------------------------------------------
SUBROUTINE CREATEFILE(narg)

  USE PARAMS
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: narg
  CHARACTER (LEN = 6) :: nmon_char
  CHARACTER (LEN = 10) :: f_char
  CHARACTER (LEN = 7) :: prefix = "MCdata_"
  CALL GETARG(narg,ext)
  
  WRITE(nmon_char,"(I0)") M
  WRITE(f_char,"(F5.2)") DS_MC

  PRINT *, "Preparing data file for ", ext, "simulations with", M,&
       & "monomers and", N, "chains ... for DS: ", f_char

  datafile = prefix//trim(adjustl(nmon_char))//"."//ext
     
  PRINT *, "Data file generated for", ext, " simulation is ",&
       & trim(datafile)
  
END SUBROUTINE CREATEFILE

!--------------------------------------------------------------------

SUBROUTINE CREATE_ATYPE()

  USE PARAMS
  IMPLICIT NONE

  INTEGER :: i,j
  INTEGER, DIMENSION(1:N,1:M) :: dstype
  REAL :: dssum, err, rannum, ranval

  dstype = -1

! Create DS for each monomer

  DO i = 1,N

     err = 9999
     
     PRINT *, "Generating sequence for chain ", i
     
     DO WHILE(err .GT. mean_tol)

        dssum = 0.0
        j = 1
        
        DO WHILE(j .LE. M)

           dstype(i,j) = ANINT(DS_fac*DS_MC*RAN1(X))
           dssum = dssum + REAL(dstype(i,j))
           
           IF(j == M .AND. dstype(i,j) == 4) CYCLE

           IF(dstype(i,j) .GT. 4 .OR. dstype(i,j) .LT. 0) CYCLE

           IF(dstype(i,j) .GT. 4 .OR. dstype(i,j) .LT. 0) THEN

              PRINT *, "Wrong dstype", i,j, dstype(i,j)
              STOP

           END IF

           !IF DS = 4, split into two monomers. Possible combinations
           ! are {(3,1),(1,3),(2,2)}

           IF(dstype(i,j) == 4) THEN
              
              ranval = ran1(X)

              IF(ranval .LE. 0.33333) THEN

                 dstype(i,j) = 2
                 dstype(i,j+1) = 2

              ELSEIF(ranval .LE. 0.66666) THEN
                 
                 dstype(i,j) = 1
                 dstype(i,j+1) = 3

              ELSE

                 dstype(i,j) = 3
                 dstype(i,j+1) = 1

              END IF
              
              j = j + 2

           ELSE
              
              j = j + 1

           END IF

        END DO

        dssum = dssum/REAL(M)
        err = abs(dssum-DS_mc)

     END DO

     WRITE(outfile,*) "DS of Chain ", i, "is", dssum

  END DO

! Create Atype

  cntatomtype = 0

  DO i = 1,N

     DO j = 1,M

        IF(dstype(i,j) == 0) THEN 
           
           atype(i,j) = 1
           cntatomtype(1) = cntatomtype(1)+1

        ELSEIF(dstype(i,j) == 3) THEN

           atype(i,j) = 8
           cntatomtype(8) = cntatomtype(8)+1
           
        ELSEIF(dstype(i,j) == 1) THEN
           
           rannum = RAN1(X)

           IF(rannum .LE. 1.0/3.0) THEN

              atype(i,j) = 2
              cntatomtype(2) = cntatomtype(2)+1

           ELSEIF(rannum .LE. 2.0/3.0) THEN

              atype(i,j) = 3
              cntatomtype(3) = cntatomtype(3)+1

           ELSE

              atype(i,j) = 4
              cntatomtype(4) = cntatomtype(4)+1

           END IF

        ELSEIF(dstype(i,j) == 2) THEN

           rannum = RAN1(X)

           IF(rannum .LE. 1.0/3.0) THEN

              atype(i,j) = 5
              cntatomtype(5) = cntatomtype(5)+1

           ELSEIF(rannum .LE. 2.0/3.0) THEN

              atype(i,j) = 6
              cntatomtype(6) = cntatomtype(6)+1

           ELSE

              atype(i,j) = 7
              cntatomtype(7) = cntatomtype(7)+1

           END IF

        ELSE

           PRINT *, "Unknown DS value", i, j, dstype(i,j)

        END IF

     END DO

  END DO

  PRINT *, "Statistics for different atomtypes"

  DO i = 1,numatomtypes

     PRINT *, "Number of atomtypes of type", i, "is", cntatomtype(i)
     
  END DO

  WRITE(outfile,*) "Statistics"

  DO i = 1,N

     cntatomtype = 0
     WRITE(outfile,*) "Chain: ", i

     DO j = 1,M

        cntatomtype(atype(i,j)) = cntatomtype(atype(i,j)) + 1

     END DO

     DO j = 1,numatomtypes

        WRITE(outfile,*) "Natomtypes of Type ", j, ": ",&
             & cntatomtype(j)

     END DO

  END DO
     
END SUBROUTINE CREATE_ATYPE

!--------------------------------------------------------------------

SUBROUTINE CREATE_ATYPE_DEF()

  USE PARAMS
  IMPLICIT NONE

  INTEGER :: i,j,k,monval
  INTEGER, PARAMETER :: defmons = INT(M*graftperMC)
  REAL, DIMENSION(1:defmons) :: dsnewval
  REAL :: dstot

  PRINT *, "Number of grafted monomers: ", defmons

  DO j = 1,N

     DO i = 1,defmons

        monval = 1+INT(RAN1(X)*M)

        atype(i,monval) = 9

     END DO

  END DO

  DO i = 1,N

     dstot = 0.0

     DO j = 1,M

        k = (i-1)*M + j
        
        IF(atype(i,j) == 1) THEN
           
           dstot = dstot + 0.0
           
        ELSEIF(atype(i,j) .LE. 4) THEN
           
           dstot = dstot + 1.0
           
        ELSEIF(atype(i,j) .LE. 7) THEN
           
           dstot = dstot + 2.0
           
        ELSEIF(atype(i,j) == 8) THEN
           
           dstot = dstot + 3.0
           
        ELSEIF(atype(i,j) == 9) THEN
           
           dstot = dstot + INT(3.0*RAN1(X))
           
        ELSE
           
           PRINT *, "Unknown atype", i,j

        END IF

     END DO

     PRINT *, "Revised DS of chain  " , i, " is: ", dstot/REAL(M)

  END DO

     
END SUBROUTINE CREATE_ATYPE_DEF

!--------------------------------------------------------------------
