!Copyright (c) 2009, 2010 by University College London, Cornell University
!Authors:
!G Pegram, Daniel P. Loucks (dpl3@cornell.edu), Marshall Taylor, Peter French, Huicheng Zhou
!Evgenii Matrosov (evgenii.matrosov@ucl.ac.uk), Julien Harou (j.harou@ucl.ac.uk)
!This program is free software under the General Public Licence, GPL (>=v2)
!Read the 'GPL License.txt' file distributed with this source code for a full license statement.
!
!	*************************************************************************************
      SUBROUTINE HYDSIM(LINK, INSTOR, OUTSTOR, LNKFLO)
C
C      USE:       Compute energy and power production at hydroelectric
C                 sites (links) and energy consumption at pumping sites.
C
C      INPUT:     INTEGER*4   LINK      - Link number
C                 Real*4      INSTOR    - Current storage in inflow node
C                 Real*4      OUTSTOR   - Current storage in outflow node.
C                 Real*4      LNKFLO    - Current flow in Link
C
C      CALLED BY  FLWSIM
C
C      OUTPUT:    ENERGY, SYSPWR (Acculates during within-year period)
C
C      NOTES:  LNKFLO is specified on input and pumping (flow) is not
C              stopped if LNKFLO is less than the minimum flow which
C              the user specified in the hydropower coefficient table.
C              Hydroelectric production will not occur if ABS(LNKFLO)
C              is less than the specified minimum flow.
C

C      Source file:      HYDSIM.FOR
C      Rewritten:        June 1993
C      Author:           M.TAYLOR
C      Modifications:
C      MRT  940803   Changed calculation of head on link connecting 
C                    nodes without storage capacity and made this
C                    routine work for bi-directional links.
!      PL   000621   Last change before Evgeni
!	 EM   0908--   Added support for wetlands and aquifers
!	 EM   090818   Changed power coefficient to 2725 to be consistend
!				   with mil m3/day
!	 EM   091007   Fixed efficiency for pumps
!	*************************************************************************************

!Energy = AQH: Q--10^6m^3/time-step; H--meters; A = 0.002725*efficiency
!internal units of power & energy: kilowatt(kw) & kwh
      IMPLICIT NONE
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
C     Arguments:
      INTEGER*4  LINK
      REAL*4     INSTOR, OUTSTOR, LNKFLO
C
C     Functions
      INTEGER*4  FINTERP
C
C     Local variables:
      INTEGER*4  ST, UPNODE, DNNODE
      REAL*4     EMAX, POWER, DENERGY, UPHEAD,DNHEAD, DHEAD
      REAL*4     HEADFLOW, LNKFLO2, TMPR4
!
!-------------------------------------------------------------------------
!     Evgenii-090817, IRAS 2010 power/pump hard-wired to metric and Q to millions m^3 per time step
      IF(HPCAP(LINK).LE.0 .AND. PCONST(LINK).LE.0 .or. LNKFLO==0.0)  !Evgenii 100710 added .or. LNKFLO==0.0
     &	GO TO 9999
C     Get "upstream" and "downstream" nodes, check for validity.
      UPNODE = NIN(LINK)
      DNNODE = NOUT(LINK)
      IF(UPNODE.LT.1 .OR. UPNODE.GT.NODMAX)GO TO 9999
      IF(DNNODE.LT.1 .OR. DNNODE.GT.NODMAX)GO TO 9999
C
C     units: 10^6 m^3/sub time-step
      LNKFLO2 = LNKFLO
      IF ((CAPN(UPNODE).LE.0.0 .AND. CAPN(DNNODE).LE.0.) .and.
     &		(.not. GWNODE(DNNODE) .AND. .not. GWNODE(upnode)))THEN !Evgenii added (.not. GWNODE(DNNODE) .AND. .not. GWNODE(upnode), b/c GWNODES don't have a capacity!
C         Both ends not storage nodes => a constant head link.
          UPHEAD = NELEV(UPNODE)
          DNHEAD = NELEV(DNNODE)
      ELSE
C        Compute head associated with initial storage of node
         ST = FINTERP( 1,NODE_VOL(1,UPNODE),NODE_ELEV(1,UPNODE),
     1                 NARVO_PTS(UPNODE),INSTOR,UPHEAD )
         !Evgenii - if not storage node, no rating table, use defined elevation
	   IF(ST.NE.SUCCES)UPHEAD = NELEV(UPNODE) 
         ST = FINTERP( 1,NODE_VOL(1,DNNODE),NODE_ELEV(1,DNNODE),
     1                 NARVO_PTS(DNNODE),OUTSTOR,DNHEAD)
         !Evgenii - if not storage node, no rating table, use defined elevation
	   IF(ST.NE.SUCCES)DNHEAD = NELEV(DNNODE)   
      END IF
     	
      DHEAD = UPHEAD - DNHEAD   !<=0 for pump; >=0 for power
      HEADFLOW = DHEAD*LNKFLO2
C     Do not compute energy production if ABS(LNKFLO2) is less than the 
C     user specified minimum production flow.
C     If HEADFLOW is positive => hydroelectric production, else
C     if HEADFLOW is negative => pumping. 
      IF(HEADFLOW.GT.0.0 .and. Powerlink(link))THEN !Power !Evgenii added .and. Powerlink(link) 100713
         IF(ABS(LNKFLO2) .LT. HPQMIN(LINK))GO TO 9999  !
         DHEAD = MIN(DHEAD, (UPHEAD-TURBINE_ELEV(LINK)))
         !Energy = AQH: Q--10^6m^3/time-step; H--meters; A = 0.002725*efficiency
         !internal units of power & energy: kilowatt & kilowatt.hours
         !DENERGY = 0.002725*ECONST(LINK)*DHEAD*LNKFLO2 
	   !Evgenii 090818, should 0.002725 not be 2725 since 9.81/3600=.002725 is for flow in m3/s, but we have 10^6 m3/s, so its .00275 *10^6
		DENERGY = ECONST(LINK)*2725*DHEAD*LNKFLO2  !This is kWh, Evgenii changed 0.002725 to 2725 
	!------------------------------------------------------
	ELSE !pump
         if (pumplink(link)) !Evgenii  if (pumplink(link))THEN 100713
     &	   DENERGY = 2725/PCONST(LINK)*DHEAD*LNKFLO2   !Evgenii - Pconst is efficiency so must divide for pumping 091007
	  !DENERGY = 0.002725*PCONST(LINK)*DHEAD*LNKFLO2 !Code pre Evgenii

      ENDIF

!	Evgenii - Energies above are in kWh/time-step
C     Check for hydropower plant capacity and plant factor limitations.
      IF(HEADFLOW.GT.0.0 .AND. DENERGY.GT.0.0)THEN
	   PLANT_FACTOR(LINK) = MIN(1.,PLANT_FACTOR(LINK)) !Evgenii- Fraction of time plant operates
         EMAX = HPCAP(LINK)*PLANT_FACTOR(LINK)*24.*DayPerTS !Evgenii- Cap (kW)*PF*24*DayperTS=kWh/ts
         DENERGY = MIN(EMAX,DENERGY) !Evgenii-  checks if energy is above the max
      ENDIF
C
C     Calculate current power and accumulate totals for period.
      POWER = 0.0
      !Evgenii commented out lines below, not supported in IRAS 2010
	!IF(HEADFLOW.GT.0.0)THEN
      !   TMPR4 = MAX(ABS(DHEAD),ABS(HEAD(LINK))) !Evgenii 090818, Head(LINK) is not defined in Iras 2010
      !   IF(TMPR4.GT.0.0)TMPR4 = ABS(DHEAD/TMPR4)
      !   POWER = HPCAP(LINK)*TMPR4
      !ENDIF
      !SYSPWR = SYSPWR + POWER
      
	ENERGY(LINK) = ENERGY(LINK) + DENERGY
C
 9999 RETURN
      END
C
