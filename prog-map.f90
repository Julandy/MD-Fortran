PROGRAM map

!To calculate density maps along xy direction in the chosen z slab
!Uses all the trajectory file, reads and analyses it line by line

IMPLICIT NONE

!Data dictionary
CHARACTER(len=30) :: inputfile
CHARACTER(len=30) :: outputfile
CHARACTER(len=4) :: element
INTEGER :: atoms, frames, atom_type, atom
INTEGER :: xbins, ybins
INTEGER :: ierror, id, ff, i, j, k, l, m, n
REAL :: z1, z2
REAL :: x, y, z, xstep, ystep
REAL :: x1, x2, y1, y2, num
INTEGER, DIMENSION (500, 500) :: c
REAL :: cella, cellb, cellc, cellg, radg

!Custom input
inputfile = "14U-oh-glu-Na.dump"
atoms = 15645
frames = 2000
outputfile = "smap-U2-0.3.map"
atom_type = 12
z1 = 104.0
z2 = 102.5
xbins = 133 !80 !40   
ybins = 123 !74 !37   
xstep = 0.3 !0.5 !0.1 
ystep = 0.3 !0.5 !0.1 
cella = 40.05375
cellb = 40.144601
cellc = 103.18517
cellg = 67.760679
radg = cellg*(3.141593/180.)
c = 0.

OPEN (9, FILE= inputfile, STATUS='OLD', ACTION='READ', IOSTAT=ierror)

errorcheck: IF (ierror < 0) THEN
            WRITE (*,1010) inputfile
            1010 FORMAT (1X, 'No ', A)
            ELSE !File is opened

OPEN (10, FILE=outputfile, STATUS='NEW', ACTION='WRITE', IOSTAT=ierror)

L1: DO i=1,frames
    READ (9,*)
    READ (9,*)
    READ (9,*)
    READ (9,*)
    READ (9,*)
    READ (9,*)
    READ (9,*)
    READ (9,*)
    READ (9,*)
L2: DO j=1,atoms
    READ (9,*) id, ff, element, x, y, z

!PBC in z direction for top/up surface

      IF (z < z2 .AND. z < (z2-cellc/2.)) THEN
         z = z + cellc 
      END IF

F1:   IF (ff == atom_type .AND. z>=z2 .AND. z<=z1) THEN

!PBC in y direction
   
     IF (y < 0.) THEN
L3:   DO 
      IF (y > 0.) EXIT
      y = y + (cellb*SIN(radg))
      x = x + (cellb*COS(radg))
      END DO L3
     ELSE IF (y > (cellb*SIN(radg))) THEN
L4:   DO
      IF (y <= (cellb*SIN(radg))) EXIT          
      y = y - (cellb*SIN(radg))
      x = x - (cellb*COS(radg))
      END DO L4
    END IF 

!PBC in x-direction

     IF (x > cella) THEN
     DO
     IF (x < cella) EXIT
     x = x - cella
     END DO
     END IF
     IF (x < 0.) THEN
L6:  DO
     IF (x > 0.) EXIT
     x = x + cella
     END DO L6
     END IF

!Counting loops

     x1 = 0.
     x2 = xstep + x1 
     y1 = 0.
     y2 = ystep + y1
L7:  DO k=1,xbins
     IF (x>=x1 .AND. x<x2) THEN
L8:  DO l = 1,ybins
     IF (y>=y1 .AND. y<y2) THEN
     c(k,l) = c(k,l) + 1
     END IF
     y1 = y1 + ystep
     y2 = y2 + ystep
     END DO L8
     END IF
     x1 = x1 + xstep
     x2 = x2 + xstep
     END DO L7
     x1 = 0.
     x2 = xstep + x1
     y1 = 0.
     y2 = ystep + y1
     END IF F1
     END DO L2
     END DO L1

!Writing loops

    x2 = xstep
    y2 = ystep

END IF errorcheck

L9: DO m=1,xbins
L10: DO n=1,ybins
     WRITE(10,1030) x2, y2, c(m,n)
     1030 FORMAT (1X, F9.4, 1X, F9.4, 1X, I9)
     y2 = y2 + ystep
     END DO L10
     y2 = ystep
     x2 = x2 + xstep
    END DO L9

END PROGRAM map   
     









