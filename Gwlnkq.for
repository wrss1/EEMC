!Copyright (c) 2009, 2010 by University College London, Cornell University
!Authors:
!Daniel P. Loucks, Huicheng Zhou(dpl3@cornell.edu), Marshall Taylor,
!Evgenii Matrosov (evgenii.matrosov@ucl.ac.uk), Julien Harou (j.harou@ucl.ac.uk)
!This program is free software under the General Public Licence, GPL (>=v2)
!Read the 'GPL License.txt' file distributed with this source code for a full license statement.
!
      subroutine GWLinkNodeSim(DSTO,DINFLW,L_SOLVED)

!     Get node inflow from all inlinks
!     Seperated from FLWSIM.FOR
!     Modifications (yymmdd):
!	PL 000611 Last change before Evgenii

!	 ********************************************************************
      IMPLICIT NONE
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  CALL: NONE
!  INPUT
      REAL*4    DSTO(NODMAX), DINFLW(NODMAX)
      LOGICAL*1 L_SOLVED(LNKMAX)
      !COMMON: TNODES,LINKS,GWNODE(NN), TOTIN(NN), INFLOW(NN)
      !        INLINK(NN,LI),GWLINK(LN), DEQLN(LN), DBQLN(LN)
      !        NIN(LN), NOUT(LN), PowerLink(LN), HPCAP(LN)
!  OUTPUT
      !LOGICAL*1 L_SOLVED(LNKMAX)
      !REAL*4    DSTO(NODMAX), DINFLW(NODMAX), INFLOW(NN), INFLOWTS(NN)
      !COMMON:   DEQLN(LN), DBQLN(LN), TOTREL(NN)
!  Local
      INTEGER*4 LI,LN,LO,NN
!-------------------------------------------------------------------------

C    Compute gw link flows on links connecting two gw nodes.
        DO LN = 1,LINKS
           IF (GWLINK(LN))  THEN
                   
            IF (GWNODE(NIN(LN)).AND.GWNODE(NOUT(LN))) THEN
C              Compute flow based on storage in nodes
               CALL GWLNKQ(LN,DSTO(NIN(LN)),DSTO(NOUT(LN)),DBQLN(LN))
C              Compute hydropower
               IF (.NOT.L_SOLVED(LN).and.
     1              ((PowerLink(LN) .and.HPCAP(LN).GT.0.)
     &              .or.(PumpLink(LN) .and.PConst(LN).GT.0.))) THEN
                  CALL HYDSIM(LN,DSTO(NIN(LN)),DSTO(NOUT(LN)),DBQLN(LN))
               END IF
               DEQLN(LN) = DBQLN(LN)
            END IF
          END IF
        ENDDO
C
C     Next compute storage volumes in groundwater nodes.
        DO NN = 1,TNODES
          IF (GWNODE(NN)) THEN
            INFLOW(NN) = MAX(0.0, INFLOW(NN) + DINFLW(NN))
            INFLOWTS(NN) = MAX(0.0, INFLOWTS(NN) + DINFLW(NN))
            IF (TOTIN(NN) .GT. 0) THEN
                DO LI = 1,TOTIN(NN)
                   LN = INLINK(NN,LI)
                   DSTO(NN) = DSTO(NN)+DEQLN(LN)
                   IF (DEQLN(LN) .GT. 0.) THEN
                     INFLOW(NN) = INFLOW(NN)+DEQLN(LN)
                     INFLOWTS(NN) = INFLOWTS(NN)+DEQLN(LN)
                   ELSE
                     TOTREL(NN) = TOTREL(NN)-DEQLN(LN)
                   END IF
                ENDDO
            END IF
            IF (TOTOUT(NN) .GT. 0) THEN
                DO LO = 1,TOTOUT(NN)
                   LN = OUTLNK(NN,LO)
                   DSTO(NN) = DSTO(NN)-DBQLN(LN)
                   IF (DBQLN(LN) .GT. 0.) THEN
                     TOTREL(NN) = TOTREL(NN)+DBQLN(LN)
                   ELSE
                     INFLOW(NN) = INFLOW(NN)-DBQLN(LN)
                     INFLOWTS(NN) = INFLOWTS(NN)-DBQLN(LN)
                   END IF
                ENDDO
            END IF
          END IF
       ENDDO
      end subroutine


C
!************************************************************************
      SUBROUTINE GWLNKQ(LINK, STOIN, STOUT, GWFLOW)
! Use four methods to compute flow throuth a groundwater link!
! Created by DPL & HCZ        June 2, 2000
      USE vars
      IMPLICIT NONE
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      INTEGER*4  LINK, L
      REAL*4     STOIN, STOUT
      !COMMON:
      ! GWMethod(Link): 0-Policy data; 1-GW-GW Horizontal; 2-GW-GW Vertical
      !                 3-GW-Storage Node; 4-GW-River/stream reaches
!  OUTPUT
      REAL*4     GWFLOW
!  CALL: GWLNKQ0(),GWLNKQ1(),GWLNKQ2(),GWLNKQ3(),GWLNKQ4()
!------------------------------------------------------------------------
      IF(nBufferedLines>0) THEN           !Anthony added 310112 to ensure GWMethod(link) has a non-default value
          DO L=1, nGW
              IF(pGW(L)%LinkID==LinkID(Link))then
                GWMethod(link) = pGW(L)%GWMethod
              END IF
          END DO
      END IF
      
          select case (GWMethod(Link))
          case (0)                                  !policy data
            call GWLNKQ0(LINK, STOIN, STOUT, GWFLOW)
          case (1)                                  !GW-GW Horizontal
            call GWLNKQ1(LINK, STOIN, STOUT, GWFLOW)
          case (2)                                  !GW-GW Vertical
            call GWLNKQ2(LINK, STOIN, STOUT, GWFLOW)
          case (3)                                  !GW-Storage Node
            call GWLNKQ3(LINK, STOIN, STOUT, GWFLOW)
          case (4)                                  !GW-River/stream reaches
            call GWLNKQ4(LINK, STOIN, STOUT, GWFLOW)
      end select
      return
      end

!************************************************************************
      SUBROUTINE GWLNKQ1(LINK, STOIN, STOUT, GWFLOW)
! for GW-GW horizontal
      USE Vars
      IMPLICIT NONE
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      INTEGER*4  LINK
      REAL*4     STOIN, STOUT
      !COMMON: DayPerTS,GWLink(Link),GWNode(Node)
      ! NARVO_PTS(Node),NODE_VOL(*,Node), NODE_ELEV(*,Node)
      ! GWK(Link),GWLength(Link),GWWidth(link),GWElev(Link)
      ! NIn(Link), NOut(Link)
!  OUTPUT
      REAL*4     GWFLOW
!  Functions
      INTEGER*4  FINTERP
!  Local
      REAL*4     Elev1, Elev2
      INTEGER*4  ST, i, NI, NO, L
!------------------------------------------------------------------------
      GWFLOW = 0.0
      !some check
      if (.not.GWLink(Link)) GOTO 9999
      IF(.not.GWNode( Nin(Link)).OR..not.GWNode(NOut(Link)))GOTO 9999
      if (NARVO_PTS(Nin(Link))<=0.or.NARVO_PTS(Nout(Link))<=0) GOTO 9999

      IF(nBufferedLines>0) THEN               !Anthony added 310112 to ensure GWWidth, GWK and GWLength have values 
          DO L=1, nGW
              IF(pGW(L)%LinkID==LinkID(Link))then
                GWWidth(link) = pGW(L)%GWWidth
                GWK(link) = pGW(L)%GWK
                GWLength(link) = pGW(L)%GWLength
              END IF
          END DO
      END IF
      
      !interpolate elevation for end nodes
      NI = NIn(Link); NO = NOut(Link)
      Elev1 = 0.0; Elev2 = 0.0
      ST = FINTERP( 1, NODE_VOL(1,NI), NODE_ELEV(1,NI),
     &         NARVO_PTS(NI), STOIN, Elev1 )
      ST = FINTERP( 1, NODE_VOL(1,NO), NODE_ELEV(1,NO),
     &         NARVO_PTS(NO), STOUT, Elev2 )
      Elev1 = MAX(0., Elev1 - GWElev(Link));
      Elev2 = MAX(0., Elev2 - GWElev(Link));
      if (GWWidth(link)>0) GWFLOW = GWK(Link) * (Elev1 + Elev2)/2.
     &         * GWLength(Link) * (Elev1 - Elev2)/GWWidth(link)
      !unit conversion: m^3/day -> 10^6 m^3/time-step
      GWFLOW = GWFLOW * DayPerTS/1.0E6
      IF (GWFLOW>0) THEN
        GWFLOW = MIN(StoIn, GWFlow)
      else
        GWFlow = MIN(StOut, ABS(GWFlow))*(-1.0)
      END IF
9999  continue
      return
      end

!************************************************************************
      SUBROUTINE GWLNKQ2(LINK, STOIN, STOUT, GWFLOW)
! for GW-GW vertical
      USE Vars
      IMPLICIT NONE
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      INTEGER*4  LINK
      REAL*4     STOIN, STOUT
      !COMMON: DayPerTS,GWLink(Link),GWNode(Node)
      ! NARVO_PTS(Node),NODE_VOL(*,Node), NODE_ELEV(*,Node)
      ! GWK(Link),GWLength(Link),GWWidth(link)
!  OUTPUT
      REAL*4     GWFLOW
!  Functions
      INTEGER*4  FINTERP
!  Local
      REAL*4     Elev1, Elev2
      INTEGER*4  NI, NO, ST, i, L
!------------------------------------------------------------------------
      GWFLOW = 0.0
      !some check
      if (.not.GWLink(Link)) GOTO 9999
      NI = Nin(Link); NO = NOut(Link)
      IF(.not.GWNode(NI).OR..not.GWNode(NO)) GOTO 9999
      if (NARVO_PTS(NI)<=0.or.NARVO_PTS(NO)<=0) GOTO 9999

      IF(nBufferedLines>0) THEN               !Anthony added 310112 to ensure GWWidth, GWK and GWLength have values 
          DO L=1, nGW
              IF(pGW(L)%LinkID==LinkID(Link))then
                GWWidth(link) = pGW(L)%GWWidth
                GWK(link) = pGW(L)%GWK
                GWLength(link) = pGW(L)%GWLength
              END IF
          END DO
      END IF
      
      !interpolate elevation for end nodes
      Elev1 = 0.0; Elev2 = 0.0
      ST = FINTERP( 1, NODE_VOL(1,NI), NODE_ELEV(1,NI),
     &         NARVO_PTS(NI), STOIN, Elev1 )
      ST = FINTERP( 1, NODE_VOL(1,NO), NODE_ELEV(1,NO),
     &         NARVO_PTS(NO), STOUT, Elev2 )
      if (GWWidth(link)>0)
     &   GWFLOW = GWK(Link)*GWLength(Link)*(Elev1-Elev2)/GWWidth(link)
      !Note: Actually GWLength(Link) is an area.
      !unit conversion: m^3/day -> 10^6 m^3/time-step
      GWFLOW = GWFLOW * DayPerTS/1.0E6
      IF (GWFLOW>0) THEN
        GWFLOW = MIN(StoIn, GWFlow)
      else
        GWFlow = MIN(StOut, ABS(GWFlow))*(-1.0)
      END IF
9999  continue
      return
      end

!************************************************************************
      SUBROUTINE GWLNKQ3(LINK, STOIN, STOUT, GWFLOW)
! For GW-Storage Node
      USE Vars
      IMPLICIT NONE
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'

!  INPUT
      INTEGER*4  LINK
      REAL*4     STOIN, STOUT
      !COMMON: DayPerTS,GWLink(Link),GWNode(Node),CapN(Node)
      ! NARVO_PTS(Node),NODE_VOL(*,Node),NODE_ELEV(*,Node), NODE_area(*,Node)
      ! GWK(Link),GWLength(Link),GWWidth(link),GWElev(Link)
!  OUTPUT
      REAL*4     GWFLOW
!  Functions
      INTEGER*4  FINTERP
!  Local
      REAL*4     Elev1, Elev2, PLength, NArea, PI
      INTEGER*4  NI, NO, ST, i, L
      PARAMETER (PI=3.1425926)
!------------------------------------------------------------------------
      GWFLOW = 0.0
      !some check
      if (.not.GWLink(Link)) GOTO 9999
      NI = NIn(Link); NO = NOut(Link)
      if (.not.GWNode(NI).and..not.GWNode(NO)) GOTO 9999
      if (NARVO_PTS(NI)<=0.or.NARVO_PTS(NO)<=0) GOTO 9999

      IF(nBufferedLines>0) THEN               !Anthony added 310112 to ensure GWWidth, GWK, GWLength and GWElev have values 
          DO L=1, nGW
              IF(pGW(L)%LinkID==LinkID(Link))then
                GWWidth(link) = pGW(L)%GWWidth
                GWK(link) = pGW(L)%GWK
                GWLength(link) = pGW(L)%GWLength
                GWElev(link) = pGW(L)%GWElev
              END IF
          END DO
      END IF
      
      !interpolate elevation for end nodes
      Elev1 = 0.0; Elev2 = 0.0
      ST = FINTERP( 1, NODE_VOL(1,NI), NODE_ELEV(1,NI),
     &         NARVO_PTS(NI), STOIN, Elev1 )
      ST = FINTERP( 1, NODE_VOL(1,NO), NODE_ELEV(1,NO),
     &         NARVO_PTS(NO), STOUT, Elev2 )
      Elev1 = MAX(0., Elev1 - GWElev(Link));
      Elev2 = MAX(0., Elev2 - GWElev(Link));
      !node area for storage node
      NArea = 0.0
      if (CapN(NI)>0 .and..not.GWNode(NI))then
        ST = FINTERP( 1, NODE_VOL(1,NI), NODE_area(1,NI),
     &         NARVO_PTS(NI), STOUT, NArea )
      else
        ST = FINTERP( 1, NODE_VOL(1,NO), NODE_area(1,NO),
     &         NARVO_PTS(NO), STOUT, NArea )
      END if
      PLength = MAX(2*GWLength(Link), 2*(PI*NArea)**0.5)
      if (GWWidth(link)>0) GWFLOW = GWK(Link) * (Elev1 + Elev2)/2.
     &         * PLength * (Elev1 - Elev2)/(GWWidth(link)/2.0)
      !unit conversion: m^3/day -> 10^6 m^3/time-step
      GWFLOW = GWFLOW * DayPerTS/1.0E6
      IF (GWFLOW>0) THEN
        GWFLOW = MIN(StoIn, GWFlow)
      else
        GWFlow = MIN(StOut, ABS(GWFlow))*(-1.0)
      END IF
9999  continue

      return
      end

!************************************************************************
      SUBROUTINE GWLNKQ4(LINK, STOIN, STOUT, GWFLOW)
! for GW-River/Stream reaches
!Evgenii- This method requires routing
      USE Vars
      IMPLICIT NONE
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      INTEGER*4  LINK
      REAL*4     STOIN, STOUT
      !COMMON: DayPerTS,GWLink(Link),GWNode(Node)
      ! NARVO_PTS(Node),NODE_VOL(*,Node), NODE_ELEV(*,Node)
      ! GWK(Link),GWLength(Link),GWWidth(link),GWElev(Link)
!  OUTPUT
      REAL*4     GWFLOW
!  Functions
      INTEGER*4  FINTERP
!  Local
      REAL*4     Elev1, Elev2, MaxDepth
      INTEGER*4  Ni, NO, ST, i, L

!------------------------------------------------------------------------
      GWFLOW = 0.0
      !some check
      if (.not.GWLink(Link)) GOTO 9999
      NI = NIn(Link); NO = NOut(Link)
      IF(.not.GWNode(NI).and..not.GWNode(NO)) GOTO 9999

      IF(nBufferedLines>0) THEN               !Anthony added 310112 to ensure GWWidth, GWK, GWLength and GWElev have values 
          DO L=1, nGW
              IF(pGW(L)%LinkID==LinkID(Link))then
                GWWidth(link) = pGW(L)%GWWidth
                GWK(link) = pGW(L)%GWK
                GWLength(link) = pGW(L)%GWLength
                GWElev(link) = pGW(L)%GWElev
              END IF
          END DO
      END IF
      
      !interpolate elevation for end nodes
      Elev1 = 0.0; Elev2 = 0.0
      IF(.not.GWNode(NI))then
        MaxDepth = 0.0
        DO i = 1, LINKS  !Evgenii changed TotIn(i) to LINKS 110714
		!For all the inflow links of the node in question (which are not demand or GW links) 
		IF (NOUT(i)==NI .and. .not.GWLink(i) .and. .not.Dmdlink(i)) !Evgenii added NOUT(i)=NI as check 110714
     &            MaxDepth = MAX(MaxDepth, LastDepth(i)) !Evgenii - lastDepth defined in simlnk.for in ComputeLinkLDW, only if lossmethod2 is activated for link
        END DO
        Elev1 = NElev(NI) + MaxDepth  !Evgenii - NElev is defined in Node Definitions in iras.inp
      else !Evgenii - If GWnode, use rating table 
        if (NARVO_PTS(NI)<=0) GOTO 9999
        ST = FINTERP( 1, NODE_VOL(1,NI), NODE_ELEV(1,NI),
     &         NARVO_PTS(NI), STOIN, Elev1 )
      END if
      IF(.not.GWNode(NO))then
        MaxDepth = 0.0
        DO i = 1, LINKS  !Evgenii changed TotIn(NO) to LINKS 110714
          IF(NOUT(i)==NO .and. .not.GWLink(i).and..not.Dmdlink(i))  !Evgenii added NOUT(i)=NO as check 110714
     &            MaxDepth = MAX(MaxDepth, LastDepth(i))
        END DO
        Elev2 = NElev(NO) + MaxDepth
      else !If GW node
        if (NARVO_PTS(NO)<=0) GOTO 9999
        ST = FINTERP( 1, NODE_VOL(1,NO), NODE_ELEV(1,NO),
     &         NARVO_PTS(NO), STOUT, Elev2 )
      END if
      Elev1 = MAX(0., Elev1 - GWElev(Link));
      Elev2 = MAX(0., Elev2 - GWElev(Link));
      if (GWWidth(link)>0)GWFLOW = GWK(Link) * (Elev1 + Elev2)/2.
     &         * 2*GWLength(Link) * (Elev1 - Elev2)/(GWWidth(link)/2.0)
      !unit conversion: m^3/day -> 10^6 m^3/time-step
      GWFLOW = GWFLOW * DayPerTS/1.0E6
      IF (GWFLOW>0) THEN
        GWFLOW = MIN(StoIn, GWFlow)
      else
        GWFlow = MIN(StOut, ABS(GWFlow))*(-1.0)
      END IF
9999  continue
      return
      end


!************************************************************************
C
      SUBROUTINE GWLNKQ0(LINK, STOIN, STOUT, GWFLOW)
C
C      USE:    Computes flow in groundwater link
C
C      INPUT:  INTEGER*4    LINK - Link number,
C              Real*4       STOIN - Available water at the "in"
C                                   node of the link, accumulated over
C                                   number of DAYS PER TimeStep
C              Real*4       STOUT - Available water at the "out"
C                                   node of the link, accumulated over
C                                   number of DAYS PER TimeStep
C              IMPLICIT INPUT:
C                 Link's groundwater transfer function.
C
C      OUTPUT: Real*4       GWFLOW - Link flow for whatever number of
C                                    days have been simulated
C
C      NOTES:  Groundwater transfers are a function of "available
C              water" at both ends of the groundwater link. "Available
C              water" is assumed to be current storage at a groundwater
C              node and storage plus inflow at a surface water node.
C
C              This routine interpolates the user input table.
C
C      Modifications:
C        PNF   930207    Use FINTERP to interpolate function (don't
C                        extrapolate).
C        MRT   940802    GW links generalized to bi-directional links
C                        and can connect surface water nodes.     
C	   EM    090527    Repaired error so aquifers and wetlands are supported
!------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      INTEGER*4      LINK
      REAL*4         STOIN, STOUT
      !COMMON:
!  OUTPUT
      REAL*4         GWFLOW
C  Functions:
      INTEGER*4      FINTERP
C  Local variables:
      INTEGER*4      NI, NO, I, IN_INDEX, OUT_INDEX, ST
      REAL*4         AVAIL_IN, AVAIL_OUT, FROM(2), TO(2), IN_VAL(2),
     1                OUT_VAL(2), INTRP(2)
      LOGICAL        FIXED_IN, FIXED_OUT
C
C     Note: AVAIL_WATER(1,x) is the FROM node; (2,x) is the TO node
C           GWTRNSFR is (FROM,TO)
      INTEGER*4      GWPOINTS(2)
      REAL*4         AVAIL_WATER(2,IAGMAX), GWTRNSFR(IAGMAX,IAGMAX)
C
!------------------------------------------------------------------------
C
      GWFLOW = 0.0
C
C     Identify incoming and outgoing nodes
      NI = NIN(LINK);  NO = NOUT(LINK)
      IF (NI.LE.0 .OR. NO.LE.0) GO TO 9999

      AVAIL_IN  = STOIN;   AVAIL_OUT = STOUT
C     If either NI or NO is not a storage node
!     then convert AVAIL_* to 10^6 m^3/day for interpolation
      IF (CAPN(NI).LE.0.0 .AND. .NOT.GWNODE(NI)) !Evgenii added  .OR. .NOT.GWNODE(NI) 090527 because GWNODES dont have capacities
     &		 AVAIL_IN = AVAIL_IN / DayPerTS 
      IF (CAPN(NO).LE.0.0 .AND. .NOT.GWNODE(NO))  !Evgenii added  .OR. .NOT.GWNODE(NI) 090527 because GWNODES dont have capacities
     &			 AVAIL_OUT = AVAIL_OUT / DayPerTS 
C
!     get transfer data for the bi-directional link
	!GWPOINTS(2),AVAIL_WATER(2,IAGMAX), GWTRNSFR(IAGMAX,IAGMAX)
      call getTransferData(Link, GWPoints, Avail_Water, GWTRNSFR)
      IF (GWALO_PTS(LINK) .LE. 0) GO TO 9999
C
C     Transfer function retrieved - perform interpolation.
      IN_INDEX = 0
      FIXED_IN = .FALSE.
      DO I = 1, GWPOINTS(1)-1
        IF (AVAIL_WATER(1,I).LE.AVAIL_IN .AND.
     1      AVAIL_IN.LE.AVAIL_WATER(1,I+1)) THEN
              IN_INDEX=I
              GO TO 20
        END IF
      ENDDO
      IF (AVAIL_IN .GT. AVAIL_WATER(1,GWPOINTS(1))) THEN
          IN_INDEX = GWPOINTS(1)
          FIXED_IN = .TRUE.
      END IF
C
20    OUT_INDEX = 0
      FIXED_OUT = .FALSE.
      DO I = 1, GWPOINTS(2)-1
        IF (AVAIL_WATER(2,I).LE.AVAIL_OUT .AND.
     1      AVAIL_OUT.LE.AVAIL_WATER(2,I+1)) THEN
              OUT_INDEX=I
              GO TO 30
        END IF
      ENDDO
      IF (AVAIL_OUT .GT. AVAIL_WATER(2,GWPOINTS(2))) THEN
          OUT_INDEX = GWPOINTS(2)
          FIXED_OUT = .TRUE.
      END IF
C
30    IF (IN_INDEX.GT.0 .AND. OUT_INDEX.GT.0) THEN
          IF (.NOT.FIXED_IN .AND. .NOT.FIXED_OUT) THEN
C             Bilinear interpolate function
              FROM(1) = AVAIL_WATER(1,IN_INDEX)
              FROM(2) = AVAIL_WATER(1,IN_INDEX+1)
              TO(1) = AVAIL_WATER(2,OUT_INDEX)
              TO(2) = AVAIL_WATER(2,OUT_INDEX+1)
              IN_VAL(1)  = GWTRNSFR(IN_INDEX,OUT_INDEX)
              IN_VAL(2)  = GWTRNSFR(IN_INDEX+1,OUT_INDEX)
              OUT_VAL(1)  = GWTRNSFR(IN_INDEX,OUT_INDEX+1)
              OUT_VAL(2)  = GWTRNSFR(IN_INDEX+1,OUT_INDEX+1)
              ST = FINTERP( 0, FROM, IN_VAL, 2, AVAIL_IN, INTRP(1) )
              ST = FINTERP( 0, FROM, OUT_VAL, 2, AVAIL_IN, INTRP(2) )
              ST = FINTERP( 0, TO, INTRP, 2, AVAIL_OUT, GWFLOW )
          ELSE IF (FIXED_IN .AND. FIXED_OUT) THEN
C             Use the max value
              GWFLOW = GWTRNSFR(IN_INDEX,OUT_INDEX)
          ELSE IF (FIXED_IN .AND. .NOT.FIXED_OUT) THEN  !120123 added .NOT.FIXED_OUT
C             Max input, interpolate for output
              TO(1) = AVAIL_WATER(2,OUT_INDEX)
              TO(2) = AVAIL_WATER(2,OUT_INDEX+1)
              INTRP(1) = GWTRNSFR(IN_INDEX,OUT_INDEX)
              INTRP(2) = GWTRNSFR(IN_INDEX,OUT_INDEX+1)
              ST = FINTERP( 0, TO, INTRP, 2, AVAIL_OUT, GWFLOW )
          ELSE
C             Max output, interpolate for input
              FROM(1) = AVAIL_WATER(1,IN_INDEX)
              FROM(2) = AVAIL_WATER(1,IN_INDEX+1)
              INTRP(1) = GWTRNSFR(IN_INDEX,OUT_INDEX)
              INTRP(2) = GWTRNSFR(IN_INDEX+1,OUT_INDEX)
              ST = FINTERP( 0, FROM, INTRP, 2, AVAIL_IN, GWFLOW )
          END IF
C         convert it to internal units: 10^6 m^3/day -> 10^6 m^3/time-step
          GWFLOW = GWFLOW * DAYPERTS
      ENDIF
C
 9999 RETURN
      END

!************************************************************************
    
      Subroutine getTransferData(Link, GWPoints, Avail_Water, GWTRNSFR)
!  Get transfer data for a bi-directional link
      implicit none
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUIT
      INTEGER*4 link
      !COMMON: SysFileName,PolicyGrpID(), LinkPolicy0(), LinkID()
      !        GWALO_PTS(), GWFromVol(),GWToVol(),GWFlowFromTo(IAGMAX*IAGMAX)
!  OUTPUT
      INTEGER*4 GWPOINTS(2)
      REAL*4    AVAIL_WATER(2,IAGMAX), GWTRNSFR(IAGMAX,IAGMAX)
!  Loacal
      INTEGER*4 i,j,ii
      LOGICAL*1 success
!-------------------------------------------------------------------------
      !get data to: GWALO_PTS(), GWFromVol(),GWToVol(),GWFlowFromTo(IAGMAX*IAGMAX)
      !                                   order by GWFromVol(), GWToVol()
      
      if (link==21) then
          
        continue
          
      end if
      
      call ReadTransfer (
     &         PolicyGrpID(SysStat(Sim_Year)),
     &         LinkPolicy0(Transfer0,Link),LinkID(Link),link,success)
      IF (GWALO_PTS(LINK) .LE. 0) GO TO 9999

      !for fromNode: find GWPoints(1),AVAIL_WATER(1,*)
      GWPoints(1) = 1
      AVAIL_WATER(1,1) = GWFromVol(1)
      do i = 2, GWALO_PTS(LINK)
        if (ABS(GWFromVol(i)-AVAIL_WATER(1,GWPoints(1)))>0.000008) then ! Threshold altered by Anthony to allow smaller values to be differentiated between 120115
          GWPoints(1) = GWPoints(1) + 1
          AVAIL_WATER(1,GWPoints(1)) = GWFromVol(i)
          if (GWPoints(1)>=IAGMAX) exit
        end if
      end do

      !for ToNode: find GWPoints(2),AVAIL_WATER(2,*)
      GWPoints(2) = 1
      AVAIL_WATER(2,1) = GWToVol(1)
      do i = 2, GWALO_PTS(LINK)
        !first line is enough in the table to get AVAIL_WATER(2,*)
        if (ABS(GWToVol(i)-GWToVol(i-1))>0.000008) then ! Threshold altered by Anthony to allow smaller values to be differentiated between 120115 ! Anthony replaced AVAIL_WATER(2,GWPoints(2)) with GWToVol(i-1) 310112
          GWPoints(2) = GWPoints(2) + 1
          AVAIL_WATER(2,GWPoints(2)) = GWToVol(i)
         if (GWPoints(2)>=(GWALO_PTS(LINK)/GWPoints(1))) exit ! This division was added as a hack by Anthony - it requires specific ordering of entries in .inp transfer table (see file format info) Anthony 120115
        else
          exit
        end if
      end do

      !find transferring flow: GWFlowFromTo(from,to)
      ii = 0
      do i = 1, GWPoints(1)        !from node
        do j = 1, GWPoints(2)      !to node
          ii = ii + 1
          GWTRNSFR(i,j) = GWFlowFromTo(ii)
        end do
      end do
9999  return
      end
