!Copyright (c) 2009, 2010 by University College London, Cornell University
!Authors:
!G Pegram, Daniel P. Loucks (dpl3@cornell.edu), Marshall Taylor, Peter French, Huicheng Zhou
!Evgenii Matrosov (evgenii.matrosov@ucl.ac.uk), Julien Harou (j.harou@ucl.ac.uk)
!This program is free software under the General Public Licence, GPL (>=v2)
!Read the 'GPL License.txt' file distributed with this source code for a full license statement.
!
!	 *************************************************************************************

! Subroutines in this file:
! * SetDefaultLoss
! * GetLoss
! * EVAPLOSS
! * SEEP

!	Modifications before Evgenii (yymmdd):
!	PL 000622 Last change before Evgenii

!	 *************************************************************************************
      SUBROUTINE SetDefaultLoss (iDay)
! Compute evap & seep and undate node storsge
! Seperated from flwSim.for
      IMPLICIT NONE
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  CALL: EVAPLOSS(), SEEP()
!  INPUT
      INTEGER*4 iDay
      !COMMON: TNODES, Links, SysEvap(iPID), PolicySysEvap(2,*)
!  OUTPUT
      !COMMON: NODE_EVAP(NN), LinkLoss(Link)
!  Local
      INTEGER*4 i, iP
!-------------------------------------------------------------------------
      !find the period
      iP = 1
      do i = 1,nSysEvap
        if (iDay <= PolicySysEvap(2,i)) ip = i
      end do
      
	IF(iP == iPSysEvap) return
      iPSysEvap = iP
      !set default value
      do i = 1, TNodes
        NODE_EVAP(i) = SysEvap(iP)
      end do
      do i = 1, Links
        LinkLoss(i) = SysEvap(iP)
      end do

      do i = 1, TNodes
        NodePolicyChg(Evaporation0,i) = .true.  !in order to read evap data of nodeself
      end do
      do i = 1, Links
        NodePolicyChg(Routing0,i) = .true.     !in order to read evap data of linkself
      end do

      end subroutine

!	 **********************************************************************
      SUBROUTINE GetLoss (DSTO)
! Compute evap & seep and undate node storsge
! Seperated from flwSim.for
      IMPLICIT NONE
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  CALL: EVAPLOSS(), SEEP()
!  INPUT
      REAL*4 DSTO(NODMAX)
      !COMMON: TNODES,GWNODE(NN)
      !        TEVAPN(NN), TSEEPL(NN)
!  OUTPUT
      !COMMON: TEVAPN(NN), TSEEPL(NN)
!  Local
      INTEGER*4 NN
      REAL*4    DEVAP, DSEEP
!-------------------------------------------------------------------------
!*** Compute evaporation, seepage losses: DEVAP and DSEEP
        DO NN = 1,TNODES
          DEVAP = 0.
          DSEEP = 0.
          IF(DSTO(NN) .GT. 0.) THEN
            !IF (.NOT. GWNODE(NN))THEN
            CALL  EVAPLOSS(NN,DSTO(NN),DEVAP)!Evgenii 120118 took out this if statement so aquifers can have evap. this also makes it so aquifers can have evap, which they shouldnt, but as long as the evap lines are not defined for aquifers and sysevap=0 they will have no evap calculated
               !Calculate initial storage volumes less evap.
               DSTO(NN) = DSTO(NN) - DEVAP
            !ENDIF
            !Calculate seepage loss if updated storage is greater zero.
            IF (DSTO(NN) .GE. 0.) THEN
                !Data needed:  NARVO_PTS(NN),NODE_SEEP(i,NN)
                CALL SEEP(NN,DSTO(NN),DSEEP)
                DSTO(NN) = DSTO(NN) - DSEEP
                IF (DSTO(NN) .LT. 0.) THEN
                  DSEEP = DSEEP + DSTO(NN)
                  DSTO(NN) = 0.
                END IF
            ELSE
                DEVAP = DEVAP + DSTO(NN)
                DSTO(NN) = 0.
                DSEEP = 0.
            END IF
          END IF
C
          !Accumulate total period evaporation and seepage
          TEVAPN(NN) = TEVAPN(NN) + DEVAP
          TSEEPL(NN) = TSEEPL(NN) + DSEEP
        ENDDO
      return
      end


!	 **********************************************************************
C
      SUBROUTINE EVAPLOSS(NN,STO,DEVAP)
C
C
C        USE:  Computes simulation time step evaporation loss from
C              node NN having storage, STO.
C
C      INPUT:  NN  - INTEGER*4   Node number
C              STO - Real*4      Current storage
C   Implicit:  Current within-year period,
C              Evaporation rate for current period,
C              Storage volume - Surface Area function for NN,
C              DayPerTS - Days per simulation time step,
C              CAPN(NN)
C
C     Output:  DEVAP - Real*4    Evaporation volume
C
C      Modifications:
C      GCP 4/3/88   - Modify include files
C      DPL 12/29/91 - Daily simulation
C      MRT 11/15/92 - Fraction of day simulation and
C                     incorporation of database functions.
!	 EM  08/10/09 - Added wetlands to Evap and Seepage calculations
C
!------------------------------------------------------------------------
C
C      DECLARATIONS
      IMPLICIT NONE
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
C     Arguments:
      INTEGER*4      NN
      REAL*4         STO, DEVAP

C     Functions:
      INTEGER*4    FINTERP
C
C     Local variables:
      INTEGER*4      ST
      REAL*4         AREAX
!------------------------------------------------------------------------
C
C     Compute partial evaporation loss from storage at
C     reservoir node for beginning of period.
      DEVAP = 0.
C
C     Check if node has any storage capacity. If not, evap = 0.
      IF(CAPN(NN) .LE. 0 .and. .not. GWNODE(NN)) GO TO 9999 !Evgenii added .not. GWNODE because wetland nodes can have evap but they dont have capacity
      IF (NARVO_PTS(NN) .LE. 0) GO TO 9999
C
C     Interpolate volume-area function about STO.
      ST = FINTERP( 1,NODE_VOL(1,NN),NODE_AREA(1,NN),
     1              NARVO_PTS(NN),STO,AREAX )
C
C     Multiply surface area by within-year period evaporation.
!	Area converted from m^2 to mil m^2      
	AREAX=AREAX /1.E6 
      IF (ST.NE.FAIL) DEVAP = AREAX * NODE_EVAP(NN)*DayPerTS

!	Returns evaporation in mil m3/sub-timestep
 9999 RETURN
      END
C
!	 **********************************************************************
C
      SUBROUTINE SEEP(NN,STO,DSEEP)
C
C        USE:  Computes simulation time step seepage loss from
C              node NN having storage, STO.
C
C      INPUT:  NN  - INTEGER*4   Node number
C              STO - Real*4      Current storage
C   Implicit:  Current within-year period,
C              Evaporation rate for current period,
C              Storage volume - Surface Area function for NN,
C              DayPerTS - Days per simulation time step,
C              CAPN(NN)
C
C     Output:  DSEEP - Real*4    Seepage volume
C
!------------------------------------------------------------------------
C
      IMPLICIT NONE
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'

      INTEGER*4      NN
      REAL*4         STO, DSEEP

C     Functions:
      INTEGER*4    FINTERP
C
C     Local variables:
      INTEGER*4      ST
!------------------------------------------------------------------------
C
      DSEEP = 0.
C
C     Check if node has any storage capacity. If not, seepage = 0.
      IF(CAPN(NN) .LE. 0.0 .and. .not. GWNODE(NN)) GO TO 9999 !Evgenii  added .not. GWNODE(NN) 0906 because GWnodes dont have a capacity
C
C     Interpolate volume-seepage function about STO.
      ST = FINTERP( 1, NODE_VOL(1,NN), NODE_SEEP(1,NN),
     1              NARVO_PTS(NN), STO, DSEEP )
C
      DSEEP = DSEEP * DayPerTS
C
 9999 RETURN
      END
C
C
