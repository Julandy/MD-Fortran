PROGRAM rdf

!Calculation of rdf or pair of atoms using lammps custom trajectory file

IMPLICIT NONE

!Variables

CHARACTER(len=20) :: filename
CHARACTER(len=20) :: outputfile
REAL :: ca, cb, cc, cg, radg
REAL :: ha, hb, hc, dhalf, trdf
INTEGER :: frames, f, lines, a, b, c, natoms
INTEGER, DIMENSION(100000) :: id, atype
CHARACTER(len=3), DIMENSION(100000) :: element
REAL, DIMENSION(100000):: x, y, z
INTEGER :: ntypes, ta, att, ierror, nti
INTEGER :: bin, bin2, bin3
REAL, DIMENSION(10000) :: nt1, nt2, nt3, nt4, nt5 
REAL, DIMENSION(10000) :: nt6, nt7, nt8, nt9, nt10
REAL, DIMENSION(10000) :: nt11, nt12, nt13, nt14
REAL, DIMENSION(10000) :: nt15, nt16, nt17, nt18, nt19
REAL :: dist, r0, r1, nframes, numt
REAL :: nden1, nden2, nden3, nden4, nden5, nden6
REAL :: nden7, nden8, nden9, nden10, nden11, nden12
REAL :: nden13, nden14, nden15, nden16, nden17, nden18, nden19

!Adjustables variables
frames = 2000
nframes = 2000.
natoms = 15645
ntypes = 19
trdf = 12          !atom type for rdf
numt = 6.         !total number of atoms of this type
nden1 = 0.00234416
nden2 = 0.00093766
nden3 = 0.00586039
nden4 = 0.02932799
nden5 = 0.00187532
nden6 = 0.05885132
nden7 = 0.00046883
nden8 = 0.00046883
nden9 = 0.00044930
nden10 = 0.00019535
nden11 = 0.00003907
nden12 = 0.00003907
nden13 = 0.00007814
nden14 = 0.00019535
nden15 = 0.00003907
nden16 = 0.00019535
nden17 = 0.00007814
nden18 = 0.00023442
nden19 = 0.00019535

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

!Cell parameters

ca = 40.05375
cb = 40.144601
cc = 103.18517
cg = 67.760679
radg = cg*(3.141593/180.)
ha = ca/2.
hb = cb/2.
hc = cc/2.
dhalf = (ha**2. + hb**2. + hc**2.)** (1./2.) 


!Files

filename = "14U-oh-glu-Na.dump"
outputfile = "U.rdf"

OPEN (9, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=ierror)

errorcheck: IF (ierror > 0) THEN
    WRITE (*,1010) filename
    1010 FORMAT (1x, 'Error: file', A, 'does not exist!')
ELSE
    OPEN (10, FILE=outputfile, STATUS='NEW', ACTION='WRITE', IOSTAT=ierror)

!Loop over frames

l0: DO f=1,frames
   
!Loop over empty lines in a frame

l1:   DO lines=1,9
        READ(9,*)
      END DO l1

!Loop in a frame to get arrays for each line-atom

l2:   DO a=1,natoms
        READ(9,*) id(a), atype(a), element(a), x(a), y(a), z(a)   
      END DO l2

!Loop in arrays, processing and counting

l3:   DO b=1,natoms

typea:  IF (atype(b) == trdf) THEN
     
!Loop over all atoms in the frame and calculate distances

l4:     DO c=1,natoms

!Wrap coordinates (half-box criteria)

!PBC along Z-coordinate

       IF (z(c) > (z(b)+hc)) THEN
         DO
         IF (z(c) < (z(b)+hc)) EXIT
         z(c) = z(c) - cc
         END DO
       ELSE IF (z(c) < (z(b)-hc)) THEN
         DO
         IF (z(c) > (z(b)-hc)) EXIT
         z(c) = z(c) + cc
         END DO
       END IF  

!PBC along Y-coordinate
       
       IF (y(c) > (y(b)+hb*SIN(radg))) THEN
         DO
         IF (y(c) < (y(b)+hb*SIN(radg))) EXIT
         y(c) = y(c) - (cb*SIN(radg))
         x(c) = x(c) - (cb*COS(radg))
         END DO 
       ELSE IF (y(c) < (y(b)-hb*SIN(radg))) THEN
         DO
         IF (y(c) > (y(b)-hb*SIN(radg))) EXIT
         y(c) = y(c) + (cb*SIN(radg))
         x(c) = x(c) + (cb*COS(radg))
         END DO
       END IF

!PBC along X-coordinate
      
       IF (x(c) > (x(b)+ha+hb*COS(radg))) THEN
         DO 
         IF (x(c) < (x(b)+ha+hb*COS(radg))) EXIT
         x(c) = x(c) - ca
         END DO
       END IF 
       IF (x(c) < (x(b)-ha-hb*COS(radg))) THEN
         DO
         IF (x(c) > (x(b)-ha-hb*COS(radg))) EXIT
         x(c) = x(c) + ca
         END DO
       END IF
       IF ((x(c) > (x(b)+ha)) .AND. (x(c) < (x(b)+ha+hb*COS(radg))) &
       .AND. ((y(c)/(x(c)-ca)) < TAN(radg))) THEN
        x(c) = x(c) - ca
       ELSE IF ((x(c) < (x(b)+ha)) .AND. (x(c) > (x(b)+ha-hb*COS(radg))) &
       .AND. (y(c) < y(b)) &
       .AND. ((y(c)/(x(c)-ca)) < TAN(radg))) THEN
        x(c) = x(c) - ca
       ELSE IF ((x(c) < (x(b)-ha)) .AND. (x(c) > (x(b)-ha-hb*COS(radg))) &
       .AND. ((y(c)/x(c)) > TAN(radg))) THEN
        x(c) = x(c) + ca
       ELSE IF ((x(c) > (x(b)-ha)) .AND. (x(c) < (x(b)-ha+hb*COS(radg))) &
       .AND. (y(c) > y(b)) &
       .AND. ((y(c)/x(c)) > TAN(radg))) THEN
        x(c) = x(c) + ca
       END IF
       
!Evaluate distances
          
  dist = ((x(c)-x(b))**2+(y(c)-y(b))**2+(z(c)-z(b))**2)**(1./2.) 
  r0 = 0.
  r1 = 0.1
      
!Assign to a bin
       
     DO bin=1,150 
       IF ((dist > r0) .AND. (dist <=r1)) THEN

!Assign to a hystogram

       IF (atype(c) == 1) THEN
          nt1(bin) = nt1(bin) + 1.
       ELSE IF (atype(c) == 2) THEN
          nt2(bin) = nt2(bin) + 1. 
       ELSE IF (atype(c) == 3) THEN
          nt3(bin) = nt3(bin) + 1. 
       ELSE IF (atype(c) == 4) THEN
          nt4(bin) = nt4(bin) + 1. 
       ELSE IF (atype(c) == 5) THEN
          nt5(bin) = nt5(bin) + 1. 
       ELSE IF (atype(c) == 6) THEN
          nt6(bin) = nt6(bin) + 1. 
       ELSE IF (atype(c) == 7) THEN
          nt7(bin) = nt7(bin) + 1. 
       ELSE IF (atype(c) == 8) THEN
          nt8(bin) = nt8(bin) + 1. 
       ELSE IF (atype(c) == 9) THEN
          nt9(bin) = nt9(bin) + 1. 
       ELSE IF (atype(c) == 10) THEN
          nt10(bin) = nt10(bin) + 1. 
       ELSE IF (atype(c) == 11) THEN
          nt11(bin) = nt11(bin) + 1. 
       ELSE IF (atype(c) == 12) THEN
          nt12(bin) = nt12(bin) + 1.
       ELSE IF (atype(c) == 13) THEN
          nt13(bin) = nt13(bin) + 1.
       ELSE IF (atype(c) == 14) THEN
          nt14(bin) = nt14(bin) + 1.
       ELSE IF (atype(c) == 15) THEN
          nt15(bin) = nt15(bin) + 1.
       ELSE IF (atype(c) == 16) THEN
          nt16(bin) = nt16(bin) + 1.
       ELSE IF (atype(c) == 17) THEN
          nt17(bin) = nt17(bin) + 1.
       ELSE IF (atype(c) == 18) THEN
          nt18(bin) = nt18(bin) + 1.
       ELSE IF (atype(c) == 19) THEN
          nt19(bin) = nt19(bin) + 1. 
       END IF
       END IF
       r0 = r0+0.1
       r1 = r1+0.1
    END DO
      
    END DO l4
    END IF typea
    END DO l3
END DO l0

!Calculate rdf (avarage over number, density and number of frames)
   r0 = 0.
   r1 = 0.1  
   DO bin2=1,150
    nt1(bin2) = ((nt1(bin2)/nframes)/numt)/(4.*3.141593*nden1*((r1-r0)**2.)*0.1)
    nt2(bin2) = ((nt2(bin2)/nframes)/numt)/(4.*3.141593*nden2*((r1-r0)**2.)*0.1)
    nt3(bin2) = ((nt3(bin2)/nframes)/numt)/(4.*3.141593*nden3*((r1-r0)**2.)*0.1)
    nt4(bin2) = ((nt4(bin2)/nframes)/numt)/(4.*3.141593*nden4*((r1-r0)**2.)*0.1)
    nt5(bin2) = ((nt5(bin2)/nframes)/numt)/(4.*3.141593*nden5*((r1-r0)**2.)*0.1)
    nt6(bin2) = ((nt6(bin2)/nframes)/numt)/(4.*3.141593*nden6*((r1-r0)**2.)*0.1)
    nt7(bin2) = ((nt7(bin2)/nframes)/numt)/(4.*3.141593*nden7*((r1-r0)**2.)*0.1)
    nt8(bin2) = ((nt8(bin2)/nframes)/numt)/(4.*3.141593*nden8*((r1-r0)**2.)*0.1)
    nt9(bin2) = ((nt9(bin2)/nframes)/numt)/(4.*3.141593*nden9*((r1-r0)**2.)*0.1)
    nt10(bin2) = ((nt10(bin2)/nframes)/numt)/(4.*3.141593*nden10*((r1-r0)**2.)*0.1)
    nt11(bin2) = ((nt11(bin2)/nframes)/numt)/(4.*3.141593*nden11*((r1-r0)**2.)*0.1)
    nt12(bin2) = ((nt12(bin2)/nframes)/numt)/(4.*3.141593*nden12*((r1-r0)**2.)*0.1)
    nt13(bin2) = ((nt13(bin2)/nframes)/numt)/(4.*3.141593*nden13*((r1-r0)**2.)*0.1)
    nt14(bin2) = ((nt14(bin2)/nframes)/numt)/(4.*3.141593*nden14*((r1-r0)**2.)*0.1) 
    nt15(bin2) = ((nt15(bin2)/nframes)/numt)/(4.*3.141593*nden15*((r1-r0)**2.)*0.1)
    nt16(bin2) = ((nt16(bin2)/nframes)/numt)/(4.*3.141593*nden16*((r1-r0)**2.)*0.1)
    nt17(bin2) = ((nt17(bin2)/nframes)/numt)/(4.*3.141593*nden17*((r1-r0)**2.)*0.1)
    nt18(bin2) = ((nt18(bin2)/nframes)/numt)/(4.*3.141593*nden18*((r1-r0)**2.)*0.1)
    nt19(bin2) = ((nt19(bin2)/nframes)/numt)/(4.*3.141593*nden19*((r1-r0)**2.)*0.1)
   r1 = r1 + 0.1
   END DO

!Write results into a file
  
  r1=0.1
  DO bin3=1,150
     WRITE(10,1020) r1, nt1(bin3), nt2(bin3), nt3(bin3), nt4(bin3),  nt5(bin3), &
     nt6(bin3), nt7(bin3), nt8(bin3), nt9(bin3), nt10(bin3), nt11(bin3), &
     nt12(bin3), nt13(bin3), nt14(bin3), nt15(bin3), nt16(bin3), & 
     nt17(bin3), nt18(bin3), nt19(bin3)
     1020 FORMAT('', F4.1, 2X, 19E10.3) 
     r1 = r1 + 0.1 
  END DO

END IF errorcheck
END PROGRAM rdf
