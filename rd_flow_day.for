!Copyright (c) 2009, 2010 by University College London, Cornell University
!Authors:
!Evgenii Matrosov (evgenii.matrosov@ucl.ac.uk),Julien Harou (j.harou@ucl.ac.uk),
!Daniel P. Loucks (dpl3@cornell.edu), Marshall Taylor, Peter French, Huicheng Zhou
!This program is free software under the General Public Licence, GPL (>=v2)
!Read the 'GPL License.txt' file distributed with this source code for a full license statement.
!
!	 *************************************************************************************
!	Modifications before Evgenii (yymmdd):
!	 000613: PNF	
!
!	 *************************************************************************************
      SUBROUTINE GetTotalFlow(Success)	!110323 Evgenii took out iMonth bc variable passed in was global 
! Evgenii took out flowfile from subroutine arguments (now global variable) 100120          
! for a day
!   1. For current time-step, read flow data from memory and multiply by flow factors.
!   2. Convert flows to Mm3/day
	USE vars
      IMPLICIT NONE
C
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
C Arguments
   !input:
      INTEGER*4 iMonth
      !COMMOM: TNodes, SysStat(NGAUGE1), NodeID(), GageNF(),
      !        GageUserUnit(UFlow), GageIDObs(GagMax)
   !output:
      logical*1 success
C
C Local Variables:
      INTEGER*4 i, j, kg, nGagesSys,Rec
      real*4 flow(GagMax)
      logical*1 Gage(NodMax)
C -----------------------------------------------------------------------
      
	Success = .false.
	!Set total gauges in system and current time-step
	nGagesSys = SysStat(NGAUGE1)
	Rec=SysStat(ntimestep)
	iMonth=SysStat(Month)
	!If past end of flow record reached, end simulation     
      if (rec>SysStat(NRec)) goto 999
	
	!initializing
	


	do i = 1, nGagesSys; flow(i) = 0.0; end do	
	do i = 1, TNodes;Gage(i)=GageNF(i);	end do
	
	    
	!Read flow data from memory for current time step and multiply by flow factor
	do i=1, SysStat(NGAUGE1)
		flow(i)= FlowData(Rec,i)
     		flow(i)=flow(i)*flowfactor(i,sysstat(run),iMonth)
	end do

	! Convert individual gauge units into Mm3/daym using conversion facters specified in iras.gag file
	do i = 1,TNODES
        if (Gage(i)) then
          do kg = 1, nGagesSys				      
            if (NodeID(i)==GageIDObs(kg)) then           !GageIDObs(i,iType) changed to GageIDObs(i), bcs no added flow 100108
              qinn(i)=flow(kg)*GageUserUnit(kg)       !GageUserUnit(kg,itype) changed to GageUserUnit(kg) bcs no added flow 100108
              exit
            end if
          end do
        end if
      end do
	Success=.true.
	RETURN


999   Success=.false.
	RETURN
	END
!	 *************************************************************************************
      SUBROUTINE GetDemandFlow()	
! Evgenii took out flowfile from subroutine arguments (now global variable) 100120          
! for a day
!   1. For current time-step, read flow data from memory and multiply by flow factors.
!   2. Convert flows to Mm3/day
	USE vars
      IMPLICIT NONE
C
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
C Arguments
   !input:
      INTEGER*4 iMonth
      !COMMOM: TNodes, SysStat(NGAUGE1), NodeID(), GageNF(),
      !        GageUserUnit(UFlow), GageIDObs(GagMax)
   !output:
!      logical*1 success
C
C Local Variables:
      INTEGER*4 i, j,rec
C -----------------------------------------------------------------------

	rec=sysstat(ntimestep)     
	!Read demand data from memory 
	do j = 1,tDemSeries
          do i = 1, tnodes				      
            if (NodeID(i)==DemIDObs(j)) then           
              DMD_TARG(i)=DemData(Rec,j)
			call UnitConversion(2,UFlow,DMD_TARG(i))
			exit
            end if
          end do
      end do
	
	

!	Success=.true.
999	Continue
	RETURN
	END

!     !EVGENII TOOK OUT GAGEFLOW MULTIPLIERS 090820
!!!!!!	****************************COMMENTED OUT********************************
!!     Call ReadObservedGageFlow(iYear,iMonth,iDay,flowObserved,Success,
!    &	SimEnd)					!SysStat(RepeatMethod) Evgenii took out repeat method 091217
	
!
!      Evgenii took out repeat method
!      iRepeat = SysStat(RepeatMethod)
	!Evgenii took out iType (as no added flow in IRAS-2010)      
!	!if (iType ==1) then        !Natural flow
! 
!	Added flow, Evgenii took out 090729
!	!else                     
!      ! nGagesSys = SysStat(NGauge2)
!      !	do i=1,TNodes; Gage(i)=GageAF(i); end do
!      !end if
!	 
!
!
!	2. and 3. below are now commented out in IRAS 2010 by Evgenii
!      --2. for the gage nodes not in observed flow file
!      --3. for other nodes
!
!	!!!2. for the gage nodes not in observed flow file while in SysFilename
!      do i = 1,TNODES
!       kn= NodSeq(i)
!        if (Gage(kn).and.(.not.bComputed(kn))) then
!          flowTot(kn) = 0.0
!          do j = 1, nGagesSys
!            if(GageMultiplier(j,iType,kn)<0.0001)cycle
!          !find No. in sequence
!            do kk = 1, TNodes
!              if (GageID(j,iType,kn) == NodeID(kk)) then
!                 kg = kk ;  exit
!             end if
!            end do
!            if (bComputed(kg)) then
!              flowTot(kn)= flowTot(kn)+flowTot(kg)
!     &                           *GageMultiplier(j,iType,kn)
!            end if
!          end do
!          bComputed(kn) = .true.
!       end if

!      end do
!      !!!3. for other nodes : changed from !!!2.
!      do kn = 1,TNODES
!	  if (.not.Gage(kn)) then
!		flowTot(kn) = 0.0
!		do j = 1, nGagesSys	
!		  if(GageMultiplier(j,iType,kn)<0.0001)cycle		
!          !find No. in sequence
!            do kk = 1, TNodes
!			if (GageID(j,iType,kn) == NodeID(kk)) then
!			   kg = kk ; 
!			   exit
!              end if
!            end do
!            if (bComputed(kg)) then
!			flowTot(kn)= flowTot(kn)+flowTot(kg)
!     &                           *GageMultiplier(j,iType,kn)
!            end if
!          end do
!         if (flowTot(kn)<0) flowTot(kn)=0.0
!          bComputed(kn) = .true.
!        end if
!      end do
!*******************************************************************

	


!*************************************************************************	
!	Subroutine  GetDaysFromYearBegin(nDays,EndMonth, EndDay)
!     get days from 1/1 to EndMonth/EndDay
!	Endday is included
!     It calls none.
!      implicit none
! Arguments
   !INPUT
!      INTEGER*4 EndMonth, EndDay
   !OUTPUT
!      INTEGER*4 nDays
!------------------------------------------------------------------------
!      SELECT CASE (EndMonth)
!          CASE (2)
!            nDays = 31
!          CASE (3)
!            nDays = 31+28
!          CASE (4)
!            nDays = 31+28+31
!          CASE (5)
!            nDays = 31+28+31+30
!          CASE (6)
!            nDays = 31+28+31+30+31
!          CASE (7)
!            nDays = 31+28+31+30+31+30
!          CASE (8)
!            nDays = 31+28+31+30+31+30+31
!          CASE (9)
!            nDays = 31+28+31+30+31+30+31+31
!          CASE (10)
!            nDays = 31+28+31+30+31+30+31+31+30
!          CASE (11)
!            nDays = 31+28+31+30+31+30+31+31+30+31
!          CASE (12)
!            nDays = 31+28+31+30+31+30+31+31+30+31+30
!          CASE DEFAULT
!            nDays = 0
!      END SELECT
!      nDays = nDays + EndDay
!      return
!      end

!********************Evgenii, since 091104 This subroutine not called because there are no gagemultipliers
!
!      SUBROUTINE GetIncrementalFlow( flowInc, Success)
! Calculate incremental flow for all nodes in terms of total flow data
! For a node:
! Incremental flow = Total natural inflow
!                    - inflows of all upstream nodes direcly connecting to it
! It calls none.
!      IMPLICIT NONE
!      INCLUDE 'IRAS_SYS.INC'
!      INCLUDE 'NODE.INC'
!      INCLUDE 'LINK.INC'
! Common
      !INTEGER*4 TNodes
      !INTEGER*4 TotIn(), Inlink(i,j)
C Arguments
   !input:
      !real*4 flowInc(NodMax)  also as output
   !output:
!      real*4 flowInc(NodMax)
!      logical*1 success
C
C Local Variables:
!      INTEGER*4 j, kn
!      real*4 flowTot(NodMax)

C -----------------------------------------------------------------------
!      do j = 1, TNodes
!        flowTot(j) = flowInc(j)
!      end do
!	GO TO 87 !EVGENII COMMENTED OUT GAGE MULTIPLIERS 
	!*****************COMMENTED OUT***************************090820
   !   do kn = 1, TNodes
   !     do j = 1, TotIn(kn)
   !      !only surface water links except for diversion links
   !       IF(GWLink(Inlink(kn,j)).or.LinDiv(Inlink(kn,j))) CYCLE
   !       flowInc(kn) = flowInc(kn)- flowTot(Inlink(kn,j))
   !     end do

    !  end do
!The incremental inflow at a node is the difference between the total uncontrolled flow at that node and the sum of all the
!uncontrolled flows at the inflow nodes of all incoming uni-directional non-diversion
!links. The uncontrolled flows are defined by the gage-flows in the flow file and the userdefined
!gage-site multipliers (that convert gage flows in the flow file to flows at gage and
!non-gage nodes). adobe pg 78 of combined IRAS manual
 !**************************************************************************************
!  87	CONTINUE !EVGENII   
!	Success = .true.
!
!999   RETURN
!      END

!This subroutine no longer used in IRAS-2010
!-------------------------------------------------------------------
!      SUBROUTINE ReadObservedGageFlow(iYear,iMonth,iDay, flow, Success,
!     &	SimEnd)		  !iRepeat, iType, filename Evgenii took out repeat method and itype 091217      
! For a day
! Get flow for day and multiply by flow factor if needed
!
! It calls none.
!      IMPLICIT NONE
!      INCLUDE 'IRAS_SYS.INC'
!      INCLUDE 'NODE.INC'
!      INCLUDE 'LINK.INC'
C Arguments
!      CHARACTER*80 FileName
!   !input:
!      INTEGER*4  iRepeat, iType, iYear,iMonth,iDay
      !COMMON:  SysStat(nGages)
      !         SysStat(BYearGage),SysStat(EYearGage),SysStat(YRStart)
   !output:
!      real*4 flow(GagMax)
!      logical*1 success, SimEnd
C  Local Variables:
!      INTEGER*4 i, j, k, iDatafile, FlowYear, FlowMonth, FlowDay
!      INTEGER*4 EndMonth, EndDay
!      CHARACTER*256 aLine
!      CHARACTER*1 Tab
!      CHARACTER*30 aVar
!      INTEGER*4 iYearRead, iMonthRead, iDayRead, nRead
!      INTEGER*4 nCurrentDays, nDaysWanted, nGagesInFile
!      REAL*4 LastFlow(GagMax)
C
C xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C
!      Success = .false.

      !initializing
!	nGagesInFile = SysStat(NGAUGE1)
!      do j = 1, nGagesInFile; flow(j) = 0.0; end do

!	Commented lines below are Evgenii's used to read iras.flw file for day's gauge flow, but
!	now flow data is stored in memory and this is not needed.
      !FlowYear = iYear; FlowMonth=iMonth; FlowDay = iDay
!      nRead = 0       !no record read
!     !iDatafile = 11
!	flwdata=44 !Evgenii added flwdata id 090911
!	READ(UNIT=flwdata,FMT='(A)',ERR=999, END=999) aLine
!	READ(aLine,*)(flow(j), j=1, nGagesInFile)							

!	if (SysStat(SimRec)==SysStat(NRec)) goto 999
!	SysStat(SimRec)=SysStat(SimRec)+1
!	do j=1, nGagesInFile
!		flow(j)= FlowData(SysStat(SimRec),j)
!		if (flowfactor(j,SysStat(month)) /=1.0) then
!     			flow(j)=flow(j)*flowfactor(j,SysStat(month))
!		end if
!	end do


!      Success = .true.
!	RETURN
!999	WRITE(*,*)'Simulation Ended, no more flow record'
!      Success=.true.
!	SimEnd=.true.
!	RETURN
!	END
C
!*************************************************************************
