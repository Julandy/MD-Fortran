PROGRAM adp

IMPLICIT NONE
! Data dictionary: variables

CHARACTER(len=20) :: inputfile
CHARACTER(len=20) :: outputfile
REAL :: ca, cb, cc, cg, radg
INTEGER :: ierror, frames, l, na, f
INTEGER, DIMENSION(100000) :: id, atype
CHARACTER(len=3), DIMENSION(100000) :: element
REAL, DIMENSION(100000):: x, y, z
INTEGER :: a1, a2, bin, bin2, bin3
REAL :: z0, z00, volume, d
REAL, DIMENSION(10000) :: nt1, nt2, nt3, nt4, nt5 
REAL, DIMENSION(10000) :: nt6, nt7, nt8, nt9, nt10
REAL, DIMENSION(10000) :: nt11, nt12, nt13, nt14
REAL, DIMENSION(10000) :: nt15, nt16, nt17, nt18, nt19

! Defined variables, file names, starting parameters

frames = 2000      !nframes
na = 15645         !natoms
z00 = 82.          !"zero" point
volume = 153573.45 !cell volume

ca = 40.05375
cb = 40.144601
cc = 103.18517
cg = 67.760679
radg = cg*(3.141593/180.)

nt1 = 0.
nt2 = 0.
nt3 = 0.
nt4 = 0.
nt5 = 0.
nt6 = 0.
nt7 = 0.
nt8 = 0.
nt9 = 0.
nt10 = 0.
nt11 = 0.
nt12 = 0.
nt13 = 0.
nt14 = 0.
nt15 = 0.
nt16 = 0.
nt17 = 0.
nt18 = 0.
nt19 = 0.

inputfile = "14U-oh-glu-Na.dump"
!outputfile = "up.adp"
outputfile = "down.adp" !Make changes in the body!

OPEN (9, FILE=inputfile, STATUS='OLD', ACTION='READ', IOSTAT=ierror)

errorcheck: IF (ierror > 0) THEN
   WRITE (*,1010) inputfile
   1010 FORMAT (1X, 'Error: file ', A,' does not exist!')
ELSE
   ! File is opened.

OPEN (10, FILE=outputfile, STATUS='NEW', ACTION='WRITE', IOSTAT=ierror)

DO f=1,frames !Frames loop

  DO l = 1,9  !info lines
    READ(9,*) 
  END DO

  DO a1=1,na   !memorize a frame
    READ(9,*) id(a1), atype(a1), element(a1), x(a1), y(a1), z(a1)   
  END DO

  DO a2=1,na

   z0 = z00              
   IF (z(a2) > z0) THEN  !Changes for up and down here
     DO
     IF (z(a2) < z0) EXIT
     z(a2) = z(a2) - cc
     END DO
   ELSE IF (z(a2) < (z0-cc)) THEN
     DO
     IF (z(a2) > (z0-cc)) EXIT
     z(a2) = z(a2) + cc
     END DO
   END IF 
   
   z0 = z00
  
   DO bin=1,150             !Changes for up and down here
     IF ((z(a2) <= z0) .AND. (z(a2) > (z0 - 0.1))) THEN
       IF (atype(a2) == 1) THEN
          nt1(bin) = nt1(bin) + 1.
       ELSE IF (atype(a2) == 2) THEN
          nt2(bin) = nt2(bin) + 1. 
       ELSE IF (atype(a2) == 3) THEN
          nt3(bin) = nt3(bin) + 1. 
       ELSE IF (atype(a2) == 4) THEN
          nt4(bin) = nt4(bin) + 1. 
       ELSE IF (atype(a2) == 5) THEN
          nt5(bin) = nt5(bin) + 1. 
       ELSE IF (atype(a2) == 6) THEN
          nt6(bin) = nt6(bin) + 1. 
       ELSE IF (atype(a2) == 7) THEN
          nt7(bin) = nt7(bin) + 1. 
       ELSE IF (atype(a2) == 8) THEN
          nt8(bin) = nt8(bin) + 1. 
       ELSE IF (atype(a2) == 9) THEN
          nt9(bin) = nt9(bin) + 1. 
       ELSE IF (atype(a2) == 10) THEN
          nt10(bin) = nt10(bin) + 1. 
       ELSE IF (atype(a2) == 11) THEN
          nt11(bin) = nt11(bin) + 1. 
       ELSE IF (atype(a2) == 12) THEN
          nt12(bin) = nt12(bin) + 1. 
       ELSE IF (atype(a2) == 13) THEN
          nt13(bin) = nt13(bin) + 1.
       ELSE IF (atype(a2) == 14) THEN
          nt14(bin) = nt14(bin) + 1.
       ELSE IF (atype(a2) == 15) THEN
          nt15(bin) = nt15(bin) + 1.
       ELSE IF (atype(a2) == 16) THEN
          nt16(bin) = nt16(bin) + 1.
       ELSE IF (atype(a2) == 17) THEN
          nt17(bin) = nt17(bin) + 1.
       ELSE IF (atype(a2) == 18) THEN
          nt18(bin) = nt18(bin) + 1.
       ELSE IF (atype(a2) == 19) THEN
          nt19(bin) = nt19(bin) + 1.
       END IF
     END IF
       z0 = z0-0.1
   END DO     
  END DO
END DO

DO bin2=1,150
   nt1(bin2) = (nt1(bin2)/(2000.*volume))
   nt2(bin2) = (nt2(bin2)/(2000.*volume))
   nt3(bin2) = (nt3(bin2)/(2000.*volume))
   nt4(bin2) = (nt4(bin2)/(2000.*volume))
   nt5(bin2) = (nt5(bin2)/(2000.*volume))
   nt6(bin2) = (nt6(bin2)/(2000.*volume))
   nt7(bin2) = (nt7(bin2)/(2000.*volume))
   nt8(bin2) = (nt8(bin2)/(2000.*volume))
   nt9(bin2) = (nt9(bin2)/(2000.*volume))
   nt10(bin2) = (nt10(bin2)/(2000.*volume))
   nt11(bin2) = (nt11(bin2)/(2000.*volume))
   nt12(bin2) = (nt12(bin2)/(2000.*volume))
   nt13(bin2) = (nt13(bin2)/(2000.*volume))
   nt14(bin2) = (nt14(bin2)/(2000.*volume))
   nt15(bin2) = (nt15(bin2)/(2000.*volume))
   nt16(bin2) = (nt16(bin2)/(2000.*volume))
   nt17(bin2) = (nt17(bin2)/(2000.*volume))
   nt18(bin2) = (nt18(bin2)/(2000.*volume))
   nt19(bin2) = (nt19(bin2)/(2000.*volume))

END DO 

d=0.1
DO bin3=1,150
  WRITE(10,1020) d, nt1(bin3), nt2(bin3), nt3(bin3), nt4(bin3),  nt5(bin3), &
     nt6(bin3), nt7(bin3), nt8(bin3), nt9(bin3), nt10(bin3), nt11(bin3), &
     nt12(bin3), nt13(bin3), nt14(bin3), nt15(bin3), nt16(bin3), &
     nt17(bin3), nt18(bin3), nt19(bin3)
  1020 FORMAT('', F4.1, 2X, 19E10.3)
  d = d + 0.1
END DO
END IF errorcheck
END PROGRAM adp
