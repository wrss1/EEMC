!Copyright (c) 2000 Cornell University
!Authors:
!Daniel P. Loucks (dpl3@cornell.edu), Marshall Taylor, Huicheng Zhou,Peter French
!This program is free software under the General Public Licence, GPL (>=v2)
!Read the 'GPL License.txt' file distributed with this source code for a full license statement.
!
!	 *************************************************************************************
      SUBROUTINE SIMLNK( LINK, QIN, QOUT)
C
C      USE:       Computes downstream link flow (link outflow).
C
C      INPUT:     INTEGER*4   LINK  - Link number,
C                 Real*4      QIN - Link inflow,
C                 IMPLICIT INPUTS:
C                  Capacity of the link (CAPL),
C                  Current within-year period (SYSSTAT(T)),
C                  Link flow-function for current period,
C                  If routing:
C                   Sub-link storage volumes at beginning of simulation step.
C                   Routing coefficients (CI, CL, CN)
C                   Total volume link at beginning of simulation step
C                   Fraction of day per simulation step (DayPerTS)
C                 COMMON:
C                   MAX_NRTRES,STEPPRPRD,CAPL(LINK), LNKVOL(LINK)
C
C      OUTPUT:    Real*4   QOUT  - Link outflow,
C                 IMPLICIT OUPUTS:
C                   If not routing:
C                      Final link volume = 0.0
C                      QOUT = QIN - losses
C                   If routing:
C                      Sub-link storages volumes at end of simulation step.
C                      Total volume link at end of simulation step
C
C      NOTES:      Based on previous SIMLNK by GCP and DPL
C                  USES VECTOR FN_RDATA(x,4)
C
C     Modifications:
C       MRT - Dec, 92. Fraction of day simulation and database routines.
C       PNF - 930125:  Change FN_RDATA from ,5 to ,4
C       MRT - 930429:  Revised routing equation and fixed redistribution
C                      of link volume when the number of routing reservoirs
C                      change.
C       MRT - 931112:  Fixed potential routing problem by checking to insure
C                      that either (CI*storage+CL*inflow) >1.0 or CN>1.0
C                      and, if not, modifing equation so as not to raise
C                      term in () to the CNth power.
C       MRT - 940202   Fixed problem with link-loss calculation.
C       HCZ - 000529   Add two methods to compute loss
C                      Width-Based method:
C                        (1)Input data method: loss depth and Flow, time of flow, width
C                        (2)Computed method: loss depth and cross section data
C                                also compute flow width, depth, velocity in channel
!
!	 *************************************************************************************
	IMPLICIT NONE
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
C
!  INPUT
       INTEGER*4 LINK
       REAL*4    QIN
       !COMMON:
       ! LossMethod(LINK): 0-Flow_Based; 1-input data of Width-Based; 2-computed of Width-Based
       ! DayPerTS
       ! LEVAP_PTS(LINK), LINK_FLOW(*,LINK), LINK_EVAP(*,LINK), LINKWidth(*,LINK)
       ! LinkLoss(Link)
       ! L_Method(LINK):Link routing. 1-default, 3-parameter; 2-cascading reservoirs, 4-parameter
       ! L_a(link), L_b(link), L_c(link)
       ! L_NRTRES(LINK),LAST_NRTRES(LINK),LNKVOL(LINK),SUB_LVOL(JJ,LINK)
!  OUTPUT
       REAL*4    QOUT
       !COMMON:
       ! LNKVOL(LINK),SUB_LVOL(JJ,LINK)
C  CALL: functions or subroutines
       INTEGER*4 FINTERP
       !ComputeLinkWDV()
!  Local variables:
       INTEGER*4    JJ, ST
       REAL*4       EVAPLN, TOTFLW, FRACTION_LOST, VOLUM, TRVTM
       REAL*4       CUM_DIST(MAX_NRTRES), CUM_VOLUME(MAX_NRTRES)
       REAL*4       POSITION, VOLUME, Tot_VOLUME, AVG_Q
       REAL*4       INQ, avStorage,aVol, LWidth


!------------------------------------------------------------------------
      
	LastLinkVOL(Link) = LnkVol(Link)
C
C     Compute total inflow and incremental link flow assuming current
C     link storage were emptied during the current within-year period.
!      TOTFLW = QIN + LNKVOL(LINK)
!      ADJFLW = (QIN + (LNKVOL(LINK)/STEPPRPRD)) * INT_TO_USER
!     Now above assumption is changed to current inflow
      TOTFLW = QIN

!     compute link loss
!	Evgenii added if conditoin so only computes loss when loss enabled, 100305
      if (iflinkloss(LINK)==.true. .or. LossMethod(LINK) == 2) then  
		call ComputeLinkLoss(Link, QIn, EvapLn)
      else
		EVAPLN=0
	endif
	
	TLossL(Link) = TLossL(Link)+ EVAPLN
C	
!	If no routing do not perform routing, condition added by Evgenii 100305
	IF (L_Method(LINK)==0) then
		QOut=QIN- EvapLN
		LnkVol(Link)=0.0
		goto 9999
	endif
	
      IF (L_Method(LINK).eq.2) GOTO 222   !routing by method 2: cascading reservoirs
C     Routing by default method(L_Method(LINK)=1): a(volume-detention(c))^b.
      avStorage = MAX(0., LnkVol(link) + QIN  - EvapLN)
      !convert to user units because the coefficients are calibrated by user units
	aVol = MAX(0.,(avStorage -  L_c(link))/LinkUserUnit(UVol))
      if (aVol < 1.0 .and. L_B(link)<1.0) then
        QOUT = MAX(0.,L_a(Link)*aVol)
      else
        QOUT = MAX(0.,L_a(link)*aVol**L_b(link))
      END if
      QOut = QOut*LinkUserUnit(UVol)         !convert back to internal units
      LnkVol(Link) = MAX(0.,avStorage - QOut)
      GOTO 9999

C     Through here only if routing by method 2: cascading reservoirs--(a*inflow+b*Vol)^c
C     Allocate this loss to flow and each res. stor. vol

222   continue

C     If the current and previous number of routing reservoirs are
C     not equal it is necessary to redistribute the link volume into
C     the current routing reservoirs evenly.
      Tot_VOLUME = 0.0
      IF(L_NRTRES(LINK).NE.LAST_NRTRES(LINK))THEN
         DO JJ = 1, LAST_NRTRES(LINK)
           Tot_VOLUME = Tot_VOLUME + SUB_LVOL(JJ,LINK)
         END do
         DO JJ = 1, L_NRTRES(LINK)
           SUB_LVOL(JJ,LINK) = Tot_VOLUME/L_NRTRES(LINK)
         END do
      ENDIF
      LAST_NRTRES(LINK) = L_NRTRES(LINK)
C
C     Compute link outflow and final link storage: (a*inflow+b*Vol)^c
      !L_a-inflow; L_b-volume; L_c-exponent
      INQ = Qin    !for the first reservoir
      DO JJ = 1,L_NRTRES(LINK)
         QOUT = (L_A(link)*InQ + L_b(LINK)*SUB_LVOL(JJ,LINK))
     &          /LinkUserUnit(UVol)
         IF(QOUT>=1.0) QOUT=QOUT**MAX(0.,Min(1.,L_c(Link)))
         QOUT = MAX(0.0, QOUT)*LinkUserUnit(UVol)
         SUB_LVOL(JJ,LINK) = MAX(0.,SUB_LVOL(JJ,LINK) + INQ
     &                      - EvapLn/L_NRTRES(LINK) - QOut)
         inQ = QOUT
      ENDDO
C
      aVOL = 0.
      DO JJ = 1,L_NRTRES(LINK)
        aVOL = aVOL + SUB_LVOL(JJ,LINK)
      ENDDO
      LNKVOL(LINK) = MAX(0.0, aVOL)
C
      !update LastWidth(Link), lastDepth() and LastVelocity() for next time-step
      !only LastWidth(Link) is used to compute loss, and others for output
9999  if (LossMethod(Link)==2) call ComputeLinkWDV(Link,Qin, QOut)
      RETURN
      END


!*************************************************************************
	subroutine ComputeLinkWDV(Link,Qin,Qout)
!  After routing, compute:
!  Width, Depth and Velocity of flow in Link
!  Authors: DPL, HCZ        May 29, 2000
      IMPLICIT NONE
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      INTEGER*4 Link
      REAL*4    QIn, QOut
      !COMMON:
      ! BaseWidth(Link),ChannelDepth(),LSlope(),RSlope(),UpLSlope(),UpRSlope()
      ! LastLinkVol(Link),LinkVol(),LnkLen()
!  OUTPUT
      !COMMON: LastWidth(Link),LastDepth(),LastVelocity()
!  Used Functions:
      INTEGER*4 FINTERP
!  Local
      INTEGER*4 ST
      REAL*4 AvVolume,TimeOfFlow,XArea,Velocity,Width,Depth
      REAL*4 ACmax,CWmax,aa,bb,cc,AXArea,ADepth
!-------------------------------------------------------------------------

      AvVolume = (LastLinkVol(Link)+LnkVol(Link))/2.0
      if (AvVolume.le.0.0) then
        Width = 0.0; Depth = 0.0; Velocity = 0.0
      else
        if (LnkLen(link)>0.and.QOut>0) then
          XArea = avVolume*1000000/LnkLen(Link)     !square meters
          TimeOfFlow = avVolume/QOut                !how many time-steps
          Velocity = LnkLen(Link)/TimeOfFlow        !m/time-step
        else
          XArea = 0; Velocity = 0;
        end if
        ACmax = ChannelDepth(Link)*(BaseWidth(Link) !max area of base channel
     &         +ChannelDepth(Link)*(LSlope(Link)+RSlope(Link))/2.0)
        CWmax = BaseWidth(Link)                     !max width of base channel
     &         +ChannelDepth(Link)*(LSlope(Link)+RSlope(Link))
        if (XArea.le.0.0) then
          Width = 0; Depth = 0
        else
          if (XArea.le.ACMax) then                  !flow within base channel
            aa = 0.5*(LSlope(Link)+RSlope(Link))
            bb = BaseWidth(Link)
            cc = XArea
            Depth =MAX(0.,(bb+(bb*bb+4*aa*cc)**0.5)/(2*aa))
            Width = BaseWidth(Link) + Depth*(LSlope(Link)+RSlope(Link))
          else                                      !flow exceeding base channel
            AXArea = XArea - ACmax
            aa = 0.5*(UpLSlope(Link)+UpRSlope(Link))
            bb = CWmax
            cc = AXArea
            ADepth =MAX(0.,(-bb+(bb*bb+4*aa*cc)**0.5)/(2*aa))
            Width = CWmax + ADepth*(UpLSlope(Link)+UpRSlope(Link))
            Depth = ChannelDepth(Link)+ADepth
          end if
        end if
      end if
      !update Common variables
      LastWidth(Link) = Width
      LastDepth(Link) = Depth
      LastVelocity(Link) = Velocity
      return
      end


!*************************************************************************
      subroutine ComputeLinkLoss(Link, Qin, EvapLn)
!  After routing, compute:
!  Width, Depth and Velocity of flow in Link
!  Authors: DPL, HCZ        May 29, 2000
!           APH            120514   
      IMPLICIT NONE
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      INTEGER*4 Link
      REAL*4    Qin         !inflow
      !COMMON:
!  OUTPUT
      REAL*4    EvapLn
      !COMMON:
!  Used Functions:
      INTEGER*4 FINTERP
!  Local
      INTEGER*4 ST
      REAL*4  LWidth, ADJFLW
      REAL*4  FLOW_ARRAY(IAGMAX), EVAP_ARRAY(IAGMAX)      !Anthony added 120514
      INTEGER*4  n                        !Anthony added 120514
!-------------------------------------------------------------------------
      EVAPLN = 0.
      ADJFLW = QIN / DayPerTS   !10^6 m^3/time-step -> 10^6 m^3/day: for intepolation
      select case (LossMethod(LINK))
        case (0)                          !flow-based
C         Interpolate flow loss function for LINK and compute this
C         time step's link loss
          if (LEVAP_PTS(LINK)>0) then
            IF (LEVAP_PTS(LINK) .EQ. 1) THEN
              EVAPLN = LINK_EVAP(1,LINK)
            ELSE
              FLOW_ARRAY = 0.         !Anthony added 120514
              EVAP_ARRAY = 0.         !Anthony added 120514
              
              do n=1,LEVAP_PTS(LINK)                  !Anthony added 120514
              FLOW_ARRAY(n) = LINK_FLOW(n,LINK)       !Anthony added 120514
              EVAP_ARRAY(n) = LINK_EVAP(n,LINK)       !Anthony added 120514
              enddo
              
              ST = FINTERP( 1, FLOW_ARRAY, EVAP_ARRAY,        !Anthony amended 120514
     1                  LEVAP_PTS(LINK),ADJFLW, EVAPLN )
            END IF
C           Convert flow user unit to internal unit: 10^6 m^3/day -> 10^6 m^3/time-step
            EVAPLN = EVAPLN * DayPerTS
          end if
        case (1)                          !Width-based: input data
          if (LEVAP_PTS(LINK)>0) then
            IF (LEVAP_PTS(LINK) .EQ. 1) THEN
              LWidth = LINKWidth(1,LINK)
            ELSE
              ST = FINTERP( 1, LINK_FLOW(1,LINK), LINKWidth(1,LINK),
     1                  LEVAP_PTS(LINK),ADJFLW, LWidth )
            END IF
            !***Note: LinkLoss unit has been converted to m/day after being read
            EVAPLN = LinkLoss(Link)* LWidth * LnkLen(Link)               !m^3/Day
            EVAPLN = EVAPLN*DayPerTS/1.0e6                               !10^6 m^3/time-step
          end if
        case (2)  !Width-based: computed in terms of cross section data
          !***Note: LinkLoss unit has been converted to m/day after being read
          EVAPLN = LinkLoss(Link)*LastWidth(Link)*LnkLen(Link)           !m^3/Day
          EVAPLN = EVAPLN * DayPerTS / 1.0e6                             !10^6 m^3/time-step
      end select

      END subroutine
      
 