!Copyright (c) 2000 Cornell University
!Authors:
!Daniel P. Loucks (dpl3@cornell.edu), Marshall Taylor, Peter French
!This program is free software under the General Public Licence, GPL (>=v2)
!Read the 'GPL License.txt' file distributed with this source code for a full license statement.
!
!	 *************************************************************************************
        INTEGER*4 FUNCTION FINTERP(MODE,XARRAY,YARRAY,NUMPTS,XVAL,YVAL)
C
C       Input:   MODE I*2  0= don't extrapolate function beyond endpts
C                          1= do extrapolate, use slope of 1st (or last)
C                              two points
C                          3= do extrapolate, but force to 0,0 on low end
C       Output:  YVAL corresponding to XVAL's position within functions
C                     YARRAY vs. XARRAY
C       Modifications (yymmdd):
C       930126:PNF - Check for 0 divide.
C       930201:PNF - Check NUMPTS, force YVAL to be YARRAY(1), reset
C                     bounds check from GE  LT to GT  LE.
!	  001515:PNF 
!
!	 *************************************************************************************
      implicit none

      INCLUDE 'IRAS_SYS.INC'     

!  INPUT
        INTEGER*4 MODE,NUMPTS
        REAL*4    XARRAY(IAGMAX), YARRAY(IAGMAX), XVAL 
!  OUTPUT
        REAL*4    YVAL
C
!  CALL: TSTBIT()
C
        LOGICAL*1  TSTBIT
C
        INTEGER*4   I, STATUS
        REAL*4      DX, DYDX, XMIN, XMAX
        LOGICAL*1   EXTRAPOLATE, ZEROEXTR
!-------------------------------------------------------------------------
        STATUS = FAIL
        YVAL = 0.0
        IF (NUMPTS .LT. 1) GO TO 9999
        EXTRAPOLATE = TSTBIT( MODE,0 )
        ZEROEXTR    = TSTBIT( MODE,1 )
        DYDX = 0.0
        YVAL = YARRAY(1)
        STATUS = SUCCES
        XMIN = 99999999.
        XMAX = -9999999.
        DO I = 1, NUMPTS
          XMIN = MIN(XARRAY(I),XMIN)
          XMAX = MAX(XARRAY(I),XMAX)
        ENDDO
C
C       If solution point outside bounds of interpolation arrays
C       then assign YVAL to appropriate endpoint value or extrapolate
C       end segment to solution point.
C
C       Below beginning point
        IF(XVAL.LT.XMIN)THEN
           IF (EXTRAPOLATE) THEN
               IF (ZEROEXTR) THEN
                   DX = XARRAY(1)
                   IF (DX.NE.0.0) DYDX = YARRAY(1)/DX
                   YVAL = DYDX*XVAL
               ELSE
                   DX = XARRAY(2)-XARRAY(1)
                   IF (DX.NE.0.0) DYDX = (YARRAY(2)-YARRAY(1))/DX
                   YVAL = YARRAY(1) + DYDX*(XVAL-XARRAY(1))
               END IF
           ENDIF
C
C       Above ending point
        ELSEIF(XVAL.GT.XMAX) THEN
           IF (.NOT.EXTRAPOLATE) THEN
               YVAL = YARRAY(NUMPTS)
           ELSE
               DX = XARRAY(NUMPTS)-XARRAY(NUMPTS-1)
               IF (DX.NE.0.0) DYDX = (YARRAY(NUMPTS)
     1                                -YARRAY(NUMPTS-1))/DX
               YVAL = YARRAY(NUMPTS) + DYDX*(XVAL-XARRAY(NUMPTS))
          ENDIF
C
C       Within function range
        ELSE
          DO I = 1, NUMPTS-1
             DX  = XARRAY(I+1)-XARRAY(I)
             IF (DX.NE.0.0) DYDX = (YARRAY(I+1)-YARRAY(I))/DX
             IF (XVAL.GE.XARRAY(I) .AND. XVAL.LE.XARRAY(I+1))THEN
                 YVAL = YARRAY(I) + DYDX*(XVAL-XARRAY(I))
                 GO TO 9999
             ENDIF
          ENDDO
        ENDIF
C
9999    FINTERP = STATUS
        RETURN
        END
C
C
C+RPA+RPA+RPA+RPA+RPA+RPA+RPA+RPA+RPA+RPA+RPA+RPA+RPA+RPA+RPA+RPA
C
      LOGICAL*1 FUNCTION TSTBIT(VARIABLE,BIT)
C
      INTEGER*4 VARIABLE, BIT, MASK
      TSTBIT = .FALSE.
      MASK = 2**BIT
      IF(IAND(VARIABLE,MASK).GT.0)TSTBIT=.TRUE.
      RETURN
      END
C

