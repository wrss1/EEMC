!Copyright (c) 2009, 2010 by University College London, Cornell University
!Authors:
!G Pegram, Daniel P. Loucks (dpl3@cornell.edu), Marshall Taylor, Peter French, Huicheng Zhou
!Evgenii Matrosov (evgenii.matrosov@ucl.ac.uk), Julien Harou (j.harou@ucl.ac.uk)
!This program is free software under the General Public Licence, GPL (>=v2)
!Read the 'GPL License.txt' file distributed with this source code for a full license statement.
!
!	 *************************************************************************************
! Subroutines in this file:
! * ReadSimDef
! * read_network_data()
! * GetUserUnits
! * ReadGageYearUnits
! * readGageUnits
! * InitVariables
! * read_year_policies()
! * InitPolicy
! * SetCurrentPolicy()
! * InitLinkRouting

!************************************************************************
      subroutine ReadSimDef(success)
! Read definitions from IRIS.DEF
! 1. to get simulation years, number of sub-time steps, number of sub-time steps,
! number of days per time step, and number of runs.

! 2. 
      implicit none
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      
	!CHARACTER(len=80) :: filename !Filename (DefFile) is a common variable
!  OUTPUT
      LOGICAL*1 success
      !COMMON:  DefFile, SYSSTAT(YRSTART), SYSSTAT(YREND),nYears
      !          SysStat(NSUBTT), SYSSTAT(NPER),sysstat(nRuns), PolicyGrpID(i)
!  Local
      INTEGER*4 i,j, iLen, iDatafile, nYears,iseq,idemprop,YearChange
      CHARACTER*256 aLine
      CHARACTER*30 aVar
!------------------------------------------------------------------------
      Success = .false.           
	iLen = LEN_TRIM(DefFile) !Character l ength of DefFile name, Evgenii
	
      !--Evgenii commented out below because files names are hardwired
	!IF(INDEX(DefFile,'.') == 0) then   !If no '.'in filename, add .DEF at end
      !  fileName(iLen+1:iLen+4) = '.DEF'
      !END if
      
	iDatafile = 11
		
      !Open formatted flow file
      OPEN(iDatafile,FILE=TRIM(DefFile),STATUS='OLD',ERR=999)
	
      READ(iDatafile, *)aVar
	
      IF(TRIM(aVar).ne.'DEFINITION')then !Checks if first line is DEFINITION
        WRITE(*,*)'The file is not a simulation definition file!'
        CLOSE(iDatafile)
        return
      
	END if
      
	do WHILE (.TRUE.)
        READ(UNIT=iDatafile,FMT='(A)',ERR=999) aLine
        	  IF (aLine(1:1).ne.'!' ) exit
      end do
	
	!Read start year, end year, number of sub-time steps, number of days per time step, and number of runs, Evgenii
	!Evgenii took out SYSTAT(NREP), IRAS 2010 wont have flow record repeat method, but added SYSTAT(NPER), SysStat(NSUBTT)

      READ(aLine,*) SYSSTAT(YRSTART), SYSSTAT(YREND),SYSSTAT(NSUBTT),
     &               SYSSTAT(NPER),SYSSTAT(nRuns),iseq !,idemprop   
	
	!READ(iDatafile, *)aVar !Move down one line
      !Evgenii 091101 Commented out user input file definitions, now they are hardwired in simsys
	!READ(UNIT=iDatafile,FMT=*,ERR=999)aVar, FlowFileName 
	!READ(UNIT=iDatafile,FMT=*,ERR=999)aVar, PolicyFileName
      !READ(UNIT=iDatafile,FMT=*,ERR=999)aVar, SysFileName
	!READ(UNIT=iDatafile,FMT=*,ERR=999)aVar, UnitsFileName
      
      userseq=.false.
      !if user defined node sequence (for priorities)      
      if (iseq==1) userseq=.true.
 
      !if demand source link propagation
      !This makes it so allocation to demand links that are sources do not reduce the
      !step deficit of associated demand nodes allowing you to propagate a demand allocation
      !down a a tree. If you have multiple allocations destined for this demand node, then you can
      !get oversupply when this is enabled
      demprop=.false.
      !if (idemprop==1) demprop=.true.      
      nYears = SYSSTAT(YREND) - SYSSTAT(YRSTART)+ 1
	sysstat(nyear)=nYears
      if (nYears < 1) then
        WRITE(*,*)'Error in simulation years definition!!'
        CLOSE(iDatafile)
        return
      end if

	

      
      !Read Policy Group IDs for each year, Evgenii
      do i = 1, nYears
        READ(UNIT=iDatafile,FMT=*,ERR=999)aVar,PolicyGrpID(i)
      end do
      
      CLOSE(iDatafile)
      Success = .true.
999   continue
	
      end subroutine


!****************************************************************************
      subroutine read_network_data(success)
! read network data from the text file: IRAS.INP
! that is prepared by WIMS
! It calls and UnitConversion() and SelSeq().
!   Get following global variables:
!       TNodes, Links
!       Nodes: --NodeID(i), NElev(i),CapN(i),StoInt(i),iGageNF,RuleSite(i),	
!           iResvNode,iNatLak, iGWNode, iDmdNode, PowerNode(i),
!           PumpNode(i), NName(i)
!		   --TotIn(i), TotOut(i), InLink(i,j), OutLnk(i,j)
!              --Node simulation sequence: NodSeq()
!              **Notation: values of InLink(i,j) and OutLnk(i,j) are in set
!                          {1, 2, ..., Links} other than link ID
!       Links: --LinkID(i),NinID,NOutID,LnkLen(i),CapL(i),StoIlt(i),iLinDiv,  
!                 iDmdLink, iTwoWay, LName(i)

!       GrpResevoirs: nGrpResv,RuleSiteID(),nResvInGrp()
	USE VARS
      implicit none
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  input
      !CHARACTER*80 sysFile
!  output
      logical*1 success
      !COMMON
!  Local variables
      CHARACTER (len = 256):: aLine
      CHARACTER*30 aVar
      INTEGER*4 i, j,idatafile
!     for nodes
      INTEGER*4 iGageNF, iGageAF
      INTEGER*4 iNatLak, iGWNode,iWetLand,iDmdNode,iResvNode
!     For links
      INTEGER*4 NInID, NOutID
      INTEGER*4 iRating, iLinDiv, iDmdLink, iTwoWay, iTrans
      LOGICAL*1 CheckRulesite(NODMax)
!------------------------------------------------------------------------
      success = .false.	
      INPFileID=66 !Evgenii added global INPfile ID, so that INP file can be opened only once. 101115
	idatafile=INPFileID
	OPEN(UNIT=iDatafile, FILE=TRIM(SysFileName), STATUS='old',
     &	 FORM='formatted', ERR= 999)
	
!     determine the number of nodes and the number of links in the network
!     find the 'NODES' line
      do WHILE (.TRUE.)
        READ(UNIT=iDatafile,FMT='(A)',ERR=999) aLine
        if (TRIM(aLine) .eq. 'NODES') EXIT
	
      end do
!     get the number of nodes in the network
      i = 0
      do WHILE (.TRUE.)
	
        READ(UNIT=iDatafile,FMT='(A)',ERR=999) aLine
        if (Index(aLine,'END OF NODES')>=1)exit
        if (aLine(1:1).ne."!") i = i + 1
      end do
	
      TNODES = i
!     get the number of links in the network
      do WHILE (.TRUE.)
        READ(UNIT=iDatafile,FMT='(A)',ERR=999) aLine
        if (TRIM(aLine) .eq. 'LINKS') EXIT
      end do
      i = 0
      do WHILE (.TRUE.)
        READ(UNIT=iDatafile,FMT='(A)',ERR=999) aLine
        if (Index(aLine, 'END OF LINKS') >= 1)exit
        if (aLine(1:1).ne."!") i = i + 1
      end do
      LINKS = i
	
!     read nodes' data and links' data

      rewind (iDatafile, ERR=999)
!     find the 'NODES' line
      do WHILE (.TRUE.)
        READ(UNIT=iDatafile,FMT='(A)',ERR=999) aLine
!	Bring curser down to Nodes line
	  if (TRIM(aLine) .eq. 'NODES') EXIT 
      end do
!     read nodes' data
      do i = 1, TNODES
        GageNF(i) = .FALSE.; !GageAF(i) = .FALSE.;
        ResvNode(i) = .FALSE.; NatLak(i) = .false.
        GWNode(i) = .FALSE.; DmdNode(i) = .false. !WetLand(i) = .FALSE. Evgenii 0904 took out wetland, not used
        do while (.true.)
!         read each node line, saves into Aline exits until no more nodes left: ! symbol
		READ(UNIT=iDatafile,FMT='(A)',ERR=999) aLine
          if (aLine(1:1) .ne. '!') exit 
        END do
        READ(aLine, *) 
     &      NodeID(i), NElev(i),CapN(i),StoInt(i),iGageNF,RuleSite(i),	!Evgenii 0904 took out 2 avar, NodLen(i),iWetLand,NodLen(i),IGageAF,N_Stat(i), because they are not used in simulation
     &       iResvNode,iNatLak, iGWNode, iDmdNode, PowerNode(i),
     &      PumpNode(i), NName(i)
        if (iGageNF > 0)then
		 GageNF(i) = .true.
		 
	  endif
       ! if (iGageAF > 0) GageAF(i) = .true. Evgenii commented this out 0911 because no Added Flow in IRAS 2009
        if (iResvNode > 0) then     
		ResvNode(i) = .true.
		
	  endif	
        if (iNatLak > 0)then
	    NatLak(i) = .true.
		
	  endif
        if (iGWNode > 0) then 
		GWNode(i) = .true. 
		
	  endif

        !if (iWetLand > 0) then !Evgenii commented this out 0904 because wetland not used in simulation (replaced by gwnode),
	  !  WetLand(i) = .true.
	  !endif
		
        if (iDmdNode > 0) DmdNode(i) = .true.
	  !unit conversion
	  CALL UnitConversion(1,ULen, NElev(i) ) 
        CALL UnitConversion(1,ULen, NodLen(i) )
        CALL UnitConversion(1,UVol, CapN(i) )
        CALL UnitConversion(1,UVol, StoInt(i) )
      
	end do
	
!     find the 'LINKS' line
      do WHILE (.TRUE.)
        READ(UNIT=iDatafile,FMT='(A)',ERR=999) aLine
        if (TRIM(aLine) .eq. 'LINKS') EXIT 
      end do
!     read links' data
      do i = 1, LINKS
        LinDiv(i) = .FALSE.; DmdLink(i) = .FALSE.; GWLink(i) = .false.
        TransLink(i)=.FALSE.
        do while (.true.)
          READ(UNIT=iDatafile,FMT='(A)',ERR=999) aLine
          if (aLine(1:1) .ne. '!') exit 
        END do
        READ(aLine, *)
     &      LinkID(i),NinID,NOutID,LnkLen(i),CapL(i),CapLYear(i), !Evgenii 0911 took out aVar,DtnVol(i),iRating, because not in simulation  
     &      StoIlt(i),iLinDiv,iDmdLink, iTwoWay, iTrans, LName(i) !Evgenii 021012 added iTrans for transfer link type
        if (iLinDiv > 0) LinDiv(i) = .true.
        if (iDmdLink > 0) DmdLink(i) = .true.
        if (iTwoWay > 0) GWLink(i) = .true.  !Evgenii - 090818 if GWlink(i)=true, it just means its  bi-directional, not neccessirally GW LINK!
		if (iTrans > 0) then
          TransLink(i) = .true. !Evgenii 021012 added iTrans for transfer link type
          end if
        !unit conversion
        CALL UnitConversion(2,ULen, LnkLen(i) )
        
        if (CapL(i)>0.0) CALL UnitConversion(2,UFlow, CapL(i) ) !flow unit !Evgenii - 1001 This is converted into mil m3/day but should be /time step (unlike all other conversions which need to be in /day)
        if (CapLYear(i)>0.0) CALL UnitConversion(2,UFlow, CapLYear(i) ) !Convert Annual Capactiy from to mil m3/day
!       Evgenii 1001 put line below in to convert CapL() to mil m3/time step from mil m3/day (needed for its later use)
	  CapL(i)=CapL(i)*sysstat(NPER)	  
	  !Convert to volume per year for annual capactiy
	  CapLYear(i)=CapLYear(i)*372.0	      !Right now annual licenses are based on 366 days because of a leap year bug, this needs to be corrected so only leap years are based on 366 days.
	  CALL UnitConversion(2,UVol, StoIlt(i) )
        !CALL UnitConversion(2,UVol, DtnVol(i) ) !DtnVol not used in IRAS2010
        !find Nin(), NOut() by NinID & NOutID
        
	  do j = 1, TNodes
          IF(NodeID(j) == NinID)NIn(i) = j
          IF(NodeID(j) == NOutID)NOut(i) = j
        end do
      end do
      
!      CLOSE(iDatafile) !Evgenii took out CLOSE so file remains open	101115
	 REWIND(iDatafile) !!Evgenii introduced rewind to replace close 101115

	
!     total number of inlinks and outlinks for a node
!     IDs of inlinks and outlinks
      do i = 1, TNodes
        TotIn(i) = 0
        TotOut(i) = 0
        do j = 1, Links
          if (NOut(j) == i)then
            TotIn(i) = TotIn(i) +1;   InLink(i,TotIn(i)) = j
          END if
          if (NIn(j) == i)then
            TotOut(i) = TotOut(i) +1; OutLnk(i,TotOut(i)) = j
          END if
        end do
      end do
	
      if (.not. userseq)then
!       Get simulation sequence for nodes
        CALL SelSeq()
      else
        CALL SetUserSeq()
      end if
!Evgenii moved these lines below from flwsim.for so that this only occurs 
C     Build an array which gives the simulation sequence order.
      DO I = 1, TNODES
        DO J = 1, TNODES
          IF(NODSEQ(J).EQ.I)SEQ_NOD(I) = J
        ENDDO
      ENDDO      
      
!     Print out node sequence
      CALL PrintSeq()
	
	
!     for group of reservoirs, get
      !nGrpResv,nResvInGrp(), RuleSiteID()
      do i = 1, TNodes
        CheckRulesite(i) = .true.
      end do
      nGrpResv = 0;

      do i = 1, TNodes
	  if (NodeID(i) == RuleSite(i)) then     !rule site
          nGrpResv = nGrpResv + 1
          RuleSiteID(nGrpResv) = NodeID(i)
          nResvInGrp(nGrpResv) = 0;
		DO j = 1, TNodes
            IF(CheckRuleSite(j).and.RuleSite(j)==RuleSite(i))then
			nResvInGrp(nGrpResv) = nResvInGrp(nGrpResv) + 1
              CheckRulesite(j) = .false.
            END if
          END DO
        end if
      end do
      !for hydropower and pump, find the links
      do i = 1, TNodes
        if (PowerNode(i)>0) then
          do j = 1, Links
            if (PowerNode(i)==LinkID(j)) then
              PowerLink(j) = .TRUE.;  exit
            end if
          end do
        end if
      end do
      do i = 1, TNodes
        if (PumpNode(i)>0) then
          do j = 1, Links
            if (PumpNode(i)==LinkID(j)) then
              PumpLink(j) = .TRUE.;  exit
            end if
          end do
        end if
      end do

      SysStat(NGAUGE1) = 0  !Number of gauges for natural flow in the system
      SysStat(nDems) = 0    !Number of demand nodes in the system
      do i = 1, TNodes
        if (GAGENF(i))SysStat(NGAUGE1) = SysStat(NGAUGE1) + 1
        if (DmdNode(i))SysStat(nDems) = SysStat(nDems) + 1 !Evgenii 101110 added number of demand nodes in the system
      end do
      !SysStat(NGAUGE2) = 0  !Number of gauges for added flow in the system Evgenii took out added gauge 091026
      !do i = 1, TNodes
      !  if (GAGEAF(i))SysStat(NGAUGE2) = SysStat(NGAUGE2) + 1
      !end do      	
      success = .true.
	return
999   CLOSE(iDatafile)
	return
      end

!************************************************************************
      subroutine GetUserUnits(success)
! Reads unit multipliers from textfile iras.unt
! Get multipliers to convert user units to internal units
      implicit none
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
!      CHARACTER*80 filename
!  OUTPUT
      LOGICAL*1 success
      !COMMON: NodeUserUnit(), LinkUserUnit())
!  Local
      INTEGER*4 i, iDatafile, n, iPos, iComp, iType
      CHARACTER*256 aLine
      CHARACTER*30 aVar
!-------------------------------------------------------------------------
!     The following are read from a system data file
      iDatafile = 11
      Success = .false.
      n = 0
      !Open formatted flow file
      
	OPEN(iDatafile,FILE=TRIM(UnitsFileName),STATUS='OLD',ERR=999) !Evgenii 100205 added END=888
	
      do WHILE (.TRUE.)
	
        READ(UNIT=iDatafile,FMT='(A)',ERR=999,END=888) aLine
	
        IF (INDEX(aLine,'Units:')>0) then
          n = n + 1
	
          iPos = INDEX(aLine, ':')
          read(aLine(iPos+1:),*)iComp, iType
          if (iComp == 1)then
            read(aLine(iPos+1:),*)aVar, aVar,NodeUserUnit(iType)

          else
            read(aLine(iPos+1:),*)aVar, aVar,LinkUserUnit(iType)
          end if
        else
          if (n>0) exit
        end if
      end do

	
      NodeUserUnit(UTime) = 1         !day -> day for flow of time
      LinkUserUnit(UTime) = 1

888      if (n>=12) Success = .true.
999   continue
      CLOSE(iDatafile)
      end subroutine


!************************************************************************
      subroutine readGageUnits(success)
!     This routing reads the following from iras.gag file: 
!    	1. conversion factors from user time-series units to IRAS units (mil m3/day)
!     2. Reads which gage nodes have flow factors.
	USE vars 
      implicit none
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      INTEGER*4 iType
      CHARACTER*80 filename
      !COMMON: TNodes, NName(j), NodeID(j)
!  OUTPUT
      LOGICAL*1 success
      !COMMON: sysstat(BYearGage),SYSSTAT(EYearGage)
      !         SysStat(nGages),GageIDObs(),GageUserUnit()
!  Local
      INTEGER*4 iDatafile, i, j, nGagesInFile,iGageFact(gagmax)
      INTEGER*4 nDemsInFile,iTimeSeries(gagmax),nSysDemands
	CHARACTER*256 aLine
      CHARACTER*30 aVar
      !CHARACTER*20 GageName(GagMax)
      CHARACTER Tab
!------------------------------------------------------------------------
      success = .false.
      iDatafile = 11
      Tab = CHAR(9)
      OPEN(UNIT=iDatafile, FILE=TRIM(FlowFileName), STATUS='old',
     &	 FORM='formatted', ERR= 999)
	continue
	do WHILE (.TRUE.)
        READ(UNIT=iDatafile,FMT='(A)',ERR=999, END=999) aLine
        IF(INDEX(aLine, 'Natural')>0)  EXIT !Evgenii took out "iType ==1" from if conditon bcs no added flow 090722
        !Evgenii took out added flow 090722 so line below commented out.
	  !IF(iType ==2.and.INDEX(aLine, 'Added')>0)  EXIT 
      end do
      !read number of gages
      do WHILE (.TRUE.)
        READ(UNIT=iDatafile,FMT='(A)',ERR=999) aLine
        IF (aLine(1:1).ne.'!' )exit
      end do 
      READ(aline, *) nGagesInFile
	
	!Check to see number of guage nodes defined in .gag file equal number defined in .inp file Evgenii 100108
	if (nGagesInFile /= SysStat(NGAUGE1)) then
	   write(*,*)'Number of gauge nodes in', FlowFileName, nGagesInFile,
     &		   ' do not match with gauges in ',SysFileName,SysStat(NGAUGE1)
	   GOTO 999
      end if 
      !SysStat(nGages1) = nGagesInFile ! Evgenii took out IF(iType==1) 090722
      !IF(iType==2) SysStat(nGages2) = nGagesInFile! Evgenii took out itype 2, its for added flow
      !read gage ID and its Multiplier 
      
	!Evgenii 091212 initialize GagFact variable
	do i=1,nGagesInFile  
		GageFact(i)=.false.
	Enddo
	
	
	do i = 1, nGagesInFile
        !read a line
          do WHILE (.TRUE.)
            READ(UNIT=iDatafile,FMT='(A)',ERR=999, END=999) aLine
            IF(aLine(1:1).ne.'!')exit
          end do
       
	  !check aLine and read in Gage name, conversion factor to IRAS units and see if Flow factors used
		DO j=1, LEN(aLine) 
			 if (aLine(j:j)==tab) aLine(j:j)=' '
		END DO
		
		READ(aline, *)GageName(i), GageUserUnit(i),iGageFact(i) !Evgenii 091212 added iGageFact (tells if gage factors exist for gage), GageUserUnit(i,itype) changed to GageUserUnit(i), bcs no added flow Evgenii 100108
		if (iGageFact(i)>0) GageFact(i)=.true. !Evgenii 091212 tells if gagfile exists for gagenode (if flowfactors for gage node are used)
      end do
	
      !find IDGage by GageName
	do i =1, nGagesInFile
        GageIDObs(i) = 0 !GageIDObs(i) changed to GageIDObs(i), bcs no added flow 100108
  	  do j = 1, TNodes
          !IF(VERIFY(TRIM(GageName(i)),TRIM(NName(j)))==0) then !Makes sure the name of the gage is the node name
          !Evgenii changed above line to the one below so gauge names cannot get confused	  28-3-11
		    IF(TRIM(GageName(i))==TRIM(NName(j))) then
		          GageIDObs(i)=NodeID(j) !!GageIDObs(i,iType) changed to GageIDObs(i), bcs no added flow 100108 Evgenii
		    end if
        end do
        IF(GageIDObs(i)==0)then !GageIDObs(i,iType) changed to GageIDObs(i), bcs no added flow 100108 Evgenii
          WRITE(*,*)'Errors in Gage name: ',TRIM(GageName(i))
          GOTO 999
        end if
      end do

	
	!Allocate demand time series variable
	nSysDemands=sysstat(nDems)
	if (.not.allocated(DemTimeSeries)) Allocate(DemTimeSeries(nSysDemands))
	if (.not.allocated(DemName)) Allocate (DemName(nSysDemands))
	if (.not.allocated(DemIDObs)) Allocate (DemIDObs(nSysDemands))
	do i=1,nSysDemands
		DemTimeSeries(i)=.false.
		DemIDObs(i)=0
	end do
	READ(UNIT=iDatafile,FMT='(A)')aline
	READ(UNIT=iDatafile,FMT='(A)')aline
	READ(aline, *) tDemSeries
	
	!Check to see number of guage nodes defined in .gag file equal number defined in .inp file Evgenii 100108
!	if (nDemsInFile /= SysStat(nDems)) then
!	   write(*,*)'Number of demand nodes in', FlowFileName,
!     &		   ' do not match with gauges in ', SysFileName
!	   GOTO 999
!      end if 


	if (tDemSeries==0) goto  888
	do i = 1, tDemSeries
          !read a line
		READ(UNIT=iDatafile,FMT='(A)',ERR=999, END=999) aLine

       
	  !check aLine and read in demand name if time series used
		!DO j=1, LEN(aLine) 
		!	 if (aLine(j:j)==tab) aLine(j:j)=' '
		!END DO
		
		READ(aline, *)DemName(i)!, iTimeSeries(i) 
		 
      end do
	
      !find IDDem by DemName
	do i =1, nSysDemands
        !DemIDObs(i) = 0 
  	  do j = 1, TNodes
          IF(TRIM(DemName(i))==TRIM(NName(j)))then !Makes sure the name of the demand is the node name
               DemIDObs(i)=NodeID(j)    
		end if
	  end do        
      end do



888   success = .true.
999   CLOSE(iDatafile)
      IF(.not.success) WRITE(*,*)'Error in reading gage or demand data!'
      end subroutine



!************************************************************************
      Subroutine InitVariables
	USE vars
! Initialize some variables
!Run once per system-run (before first time-step)
! Calls: none
      implicit none
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  Local
      INTEGER*4 NN, NODE, JJ, i,j
!------------------------------------------------------------------------
      iPSysEvap = 0   !set current no. of period for system evap
      do i = 1, maxPolicies 
        SysEvap(i) = 0.0
      end do
C     Initialize INFLOW, TOT_TARGET between periods since these are used
C     in FLWSIM to compute inter-period carryover.
      DO NN = 1, TNODES
         INFLOW(NN) = 0.
         STEP_DEFICIT(NN) = 0.
	   DO JJ=1,SYSSTAT(nyear)
			DO i=1,thres_pts_max
				 YearFailEvent(JJ,NN,I)=.false.
		    END DO
	   END DO
      ENDDO
C
C     Read time-independent node function data (vol/area fns) into memory
C     Initialize vectors
      DO NODE = 1,NODMAX
         SUPLY_PTS(NODE) = 0
         LAKEQ_PTS(Node) = 0
         NARVO_PTS(NODE) = 0
	   NODE_EVAP(Node) = 0.0     
	   perf_pts(Node) = 0        !Evgenii 100617
	   Annual_DMD(node)=0
	   Annual_SHRTG(node)=0
	   annualSI(NODE) = 0
	   TSflw_DEFICIT(NODE) = 0
!	   annualSDshrtg(NODE) = 0
 	   TSflw_DEFICIT(NODE) = 0
!	   annualSDdem(NODE) = 0
	   TsSIsum(NODE) = 0
         TS_SHRTG(Node)=0
         TSflw_Surplus(node)=0
	   TS_Surp(NODE) = 0
	   Ts_DMD_Sum(NODE) = 0
	   NODE_EVAP_ON(Node)=.false.
	   nodesourcechange(node)=.false. 
	   DO JJ = 1, IAGMAX
            NODE_VOL(JJ,NODE)  = 0.0
            NODE_AREA(JJ,NODE) = 0.0
            NODE_ELEV(JJ,NODE) = 0.0
            NODE_SEEP(JJ,NODE) = 0.0
         END DO
	   DO JJ = 1, MXSUPLY
		  Source_Type(jj,node)=0
		  
	   ENDDO
         sto_perf_node(node)=.false. !Evgenii 100617
         !flow_perf_node(node)=.false. !Evgenii 120221
	   seep_node(node)=.false.   !Evgenii 100617
	   cons_node(NODE)=.false.   !Evgenii 100617
	   srcpriorities(node)=.false.
	   DO JJ = 1, thres_pts_max  !Evgenii 100617
	      thres_limit(jj,node)=0.0
		  nTime_Steps(jj,node)=0
	      last_T_fail(jj,node)=.false.
		  nTime_Steps_Tot(jj,node)=0
		  max_sto_fail_dur(jj,node)=0
	      nTime_Steps(jj,node)=0
		  nStoreFail(jj,node)=0
		  nTime_Steps_Tot(jj,node)=0
	   End Do

	   




	
	   EnvFlwNode(node)=.false. !added by Evgenii 100601
!	   LAST_STEP_DEFICIT(NODE)=0.0 !Added by Evgenii 100303
      END DO
C        Initialize storage states.  (Sets value of ESTO(NN) )
C        Initialize end storage volume to values read in as data, number
C        of subreservoirs
         DO NN = 1,TNODES
           BSTO(NN) = STOINT(NN);   ESTO(NN) = BSTO(NN) 
         ENDDO

!        ground water
         do i = 1, Links
           !groundwater transfer
           GWMethod(i) = 0
           GWK(i)=0; GWElev(i)=0; GWLength(i)=0; GWWidth(i)=0
		 AnnualEng(i)=0; AnnualCost(i)=0;GlobannCost(i)=0.
		 maxannualcost(i)=0.;ndemnodelink(i)=0;YearQLN(i)=0.
           !link loss for routitng when loss method = 0 and 1
           do j=1,IAGMAX
             LinkWidth(j,i) = 0.0; LINK_EVAP(j,i) = 0.0;
             LINK_FLOW(j,i) = 0.0; LINK_TTV(j,i) = 0.0;
           end do
           !Cross section for routing when lossm ethod = 2
           LossMethod(i) = 0.0; LinkLoss(i)=0.0
           baseWidth(i) = 0.0;  channelDepth(i) = 0.0
           LSlope(i) = 0.0;     RSlope(i) = 0.0;
           UpLSlope(i)=0.0;     UpRSlope(i) = 0.0
		   !Costing and power !Added by Evgenii 100713
		   EngTot(i)= 0.0; EngTotHydSim(i)=0.0; 
	 	   CostTot(i)=0.0; TotPower(i)=0.0;TotPower(i)=0.0;
		   FlowEng(i)=0.0;AnnualEng(i)=0.0;
		   AnnualCost(i)=0.0;
	       TotPower(i)=0.0;AnnCostInc(i)=0.0;
	       maxannualcost(i)=0.0;GlobannCost(i)=0.0;AnnualEng(i)=0
	       EngTotHydSim(i)=0.0;EngTot(i)=0.0;totflow(i)=0.0		   
         end do
      RETURN
      END

!************************************************************************
      Subroutine InitLinkRouting
! Initializing link volume and sub volumes
! Calls: none
      implicit none
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  Local
      INTEGER*4 i,j
!------------------------------------------------------------------------
C     Initialize number of previous routing reservoirs
      DO I = 1,LNKMAX
       !for loss computing
        LastWidth(i) = 0.0;LastDepth(i) = 0.0;LastVelocity(i) = 0.0
        LNKVOL(I) = STOILT(I); LastLinkVOL(i) = LNKVOL(I)
        LAST_NRTRES(I) = 0 !This is the last number of routing reservoirs -Evgenii
      ENDDO

      !initializing for routing. 
      DO i = 1,LINKS
        if (L_NRTRES(i)>0) then
          IF (NIN(i).GT.0 .AND. NOUT(i).GT.0) THEN
            DO j = 1,L_NRTRES(i)
              SUB_LVOL(j,i) =  LNKVOL(i)/L_NRTRES(i)
            ENDDO
          END IF
          LAST_NRTRES(i) = L_NRTRES(i)
        END if
      END do
      return
      end

	
!************************************************************************
      subroutine InitPolicy
      USE vars
	implicit none
      include 'IRAS_SYS.INC'
      include 'NODE.INC'
      include 'LINK.INC'
!  Local
      INTEGER*4 i,nn, jj
!------------------------------------------------------------------------
      !Initializing for Policy
      
      !Kang add for improving performance
      !Initialize the Min and Max ID of current policy and policy type
      MinNodePolicyID = 1
      MaxNodePolicyID = MaxPolicies
      MinLinkPolicyID = 1
      MaxLinkPolicyID = MaxPolicies
	
	do jj = 1,MaxPolicies   !for system evaporation
        PolicySysEvap(1,jj) = 0.0
        PolicySysEvap(2,jj) = 0.0
      END do
      do NN = 1, tnodes !NODMAX ! !          !for node policy (or Period)
        do jj = 1,MaxPolicyTypes
          NodePolicy0(jj,NN) = 1
          NodePolicyChg(jj,NN) = .true.
          do i = 1, MaxPolicies
            NodePolicyBeg(i,jj,NN) = 0.0
          end do
        end do
      end do
      do NN = 1,links !LNKMAX ! !           !for link policy (or period)
        do jj = 1,MaxPolicyTypes
          LinkPolicy0(jj,NN) = 1
          LinkPolicyChg(jj,NN) = .true.
          do i = 1, MaxPolicies
            LinkPolicyBeg(i,jj,NN) = 0.0
          end do
        end do
      end do
	
      return
      end subroutine


!************************************************************************
      subroutine SetCurrentPolicy(iDay)   !Kang modify (rDay)
      use vars
	implicit none
      include 'IRAS_SYS.INC'
      include 'NODE.INC'
      include 'LINK.INC'
!  INPUT
      !Kang modify 
      !The following modification is not reasonable because it results in the difference
      !of parameter type between caller who uses INTEGER*4 and callee who uses INTEGER*4
      !Consequently, rDay is random particularly in release version
      !This kind of problem can make a system unstable
      !Intel Fortran compiler helps developers check this kind of problem by setting 
      !"Diagnostics/Language Usage Warnings/Check Routine Interfaces" to "Yes"
      !By the way, there are many warnings about noalignments which can be disabled by 
      !setting "Diagnostics/Language Usage Warnings/Warning for Unaligned Data" to "No"
      !integer*4 rDay ! 0904 Evgenii changed rDay to integer  (otherwise error in compile)
      INTEGER*4 iDay

      !COMMON: TNodes, Links, NodePolicyBeg(),LinkPolicyBeg()
!  OUTPUT
      !COMMON
	!  NodePolicy0(Policytype,Node), NodePolicyChg()
      !  LinkPolicy0(Policytype,Link), LinkPolicyChg()
!  Local
      INTEGER*4 NN, PT, PP
      INTEGER*4 rDay       !Kang add
      real*4 rrDay
!------------------------------------------------------------------------
      rDay = iDay       !Kang add
      rrDay=real(iDay)
      !Kang modify for improving performance 20100701
      !Here should use TNODES instead of NODMax  
      !do nn = 1, NODMax
      do nn = 1, TNODES
        do PT = 1, MaxPolicyTypes
          !Kang modify for improving performance 20100701
          !do PP = 1, maxPolicies
          do PP = MinNodePolicyID, MaxNodePolicyID
            !Kang modify for improving performance 20100701
            !Difference between the 1st and 2nd condition??
            if (NodePolicyBeg(PP,PT,NN)>1.0.and.
     &          rDay>= NodePolicyBeg(PP,PT,NN).and.
     &          rrDay< NodePolicyEnd(PP,PT,NN).and.
     &          NodePolicy0(PT,NN)<PP) then
              NodePolicyChg(PT,NN) = .true.
              NodePolicy0(PT,NN) = PP
              exit
            end if
          end do
        end do
      end do
      
      !Kang modify for improving performance 20100701
      !Here should use LINKS instead of LNKMax        
      !do nn = 1, LNKMax
      do nn = 1, LINKS
        do PT = 1, MaxPolicyTypes
          !Kang modify for improving performance 20100701
          !do PP = 1, maxPolicies
          do PP = MinLinkPolicyID, MaxLinkPolicyID
            !Kang modify for improving performance 20100701
            !Difference between the 1st and 2nd condition??
            if (LinkPolicyBeg(PP,PT,NN)>1.0.and.
     &          rDay >= LinkPolicyBeg(PP,PT,NN).and.
     &          rrDay< LinkPolicyEnd(PP,PT,NN).and.
     &          LinkPolicy0(PT,NN)<PP) then
              LinkPolicyChg(PT,NN) = .true.
              LinkPolicy0(PT,NN) = PP
              exit
            end if
          end do
        end do
      end do
      return
      end subroutine


!************************************************************************
      SUBROUTINE read_year_policies(iPolicyGrp,success)
! Read policies for nodes and links for a given policyGrp ID from iras.pol
      use VARS
	implicit none
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      !character*80 FileName
      INTEGER*4 iPolicyGrp
!  OUTPUT
      LOGICAL*1 success
      !COMMON: PolicyGrp()
!  CALL:
      LOGICAL ReadPolicyOneLine    !Function   !Kang modify for improving performance
!  Local
      INTEGER*4 i, iDatafile, iPGrp,iPType,iPID,iCompType,iNodeLinkID
      REAL*4    BegDT, EndDT
      CHARACTER*30 aVar
      CHARACTER*256 aLine
      INTEGER  iLine               !Kang modify for improving performance
!------------------------------------------------------------------------
      success = .false.
      
      !Kang modofy for improving performace 
!      iDatafile = 11
!      OPEN(UNIT=iDatafile, FILE=TRIM(PolicyFileName), STATUS='old',
!     &	 FORM='formatted', ERR= 999)
      IF(totalLinePolicyFile>0 .AND. ASSOCIATED(pPolicyData)) THEN
        iLine = 1
      ELSE
        iDatafile = iPolicyFile
        REWIND(UNIT=iDatafile)
	END IF
	
      !Kang add for improving performance
      !Initialize the Min and Max ID of current policy and policy type
      MinNodePolicyID = MaxPolicies
      MaxNodePolicyID = 1
      MinLinkPolicyID = MaxPolicies
      MaxLinkPolicyID = 1
      
      !read time periods for system evaporation, default data for any node and link
	!puts data in to real array PolicySysEvap
      nSysEvap = 0
      do WHILE (.TRUE.)
        !Kang modify for improving performance
        IF(totalLinePolicyFile>0 .AND. ASSOCIATED(pPolicyData)) THEN
            IF(ReadPolicyOneLine(iLine, aLine)) GOTO 999
        ELSE    
        READ(UNIT=iDatafile,FMT='(A)',ERR=999, END=999) aLine
        END IF
        
        if (Index(aLine, 'SysEvaporation:') >= 1) exit
      end do
      nSysEvap = 1; i = 1
      !Read policy ID, Beg day, end day
	read(aLine,*)aVar,iPID,PolicySysEvap(1,iPID),PolicySysEvap(2,iPID)
	
      do WHILE (.TRUE.)
        !Kang modify for improving performance
        IF(totalLinePolicyFile>0 .AND. ASSOCIATED(pPolicyData)) THEN
            IF(ReadPolicyOneLine(iLine, aLine)) GOTO 999
        ELSE    
        READ(UNIT=iDatafile,FMT='(A)',ERR=999, END=999) aLine
        END IF
        
        if (Index(aLine, 'SysEvaporation:') == 0) then
          nSysEvap = i;  exit
        END if
        i = i + 1
        read(aLine,*)aVar,iPID,PolicySysEvap(1,i),PolicySysEvap(2,i)
      end do
	
      !find line: '!POLICIES'
     
	do WHILE (.TRUE.)
        !Kang modify for improving performance
        IF(totalLinePolicyFile>0 .AND. ASSOCIATED(pPolicyData)) THEN
            IF(ReadPolicyOneLine(iLine, aLine)) GOTO 999
        ELSE    
        READ(UNIT=iDatafile,FMT='(A)',ERR=999, END=999) aLine
        END IF
        
	    if (Index(aLine, 'POLICIES') >= 1) exit
	end do
      
	do WHILE (.TRUE.)
        !Kang modify for improving performance
        IF(totalLinePolicyFile>0 .AND. ASSOCIATED(pPolicyData)) THEN
            IF(ReadPolicyOneLine(iLine, aLine)) GOTO 999
        ELSE    
  	  READ(UNIT=iDatafile,FMT='(A)',ERR=999, END=888) aLine
        END IF
        
        if (Index(aLine, 'END OF POLICIES') >= 1) exit

	  IF(aLine(1:1).eq.'!') cycle
        read(aLine,*)aVar, iPGrp
     
	  if (iPGrp.eq.iPolicyGrp) then
          read(aLine,*)aVar,aVar,iPType,iPID,iCompType,iNodeLinkID, 
     &                 BegDT, EndDT
          if (iCompType .eq. 1) then    !for nodes
          !Kang add for improving performance
          !Decrease Min ID and increase Max ID of current policy 
          if(MinNodePolicyID>iPID) MinNodePolicyID = iPID
          if(MaxNodePolicyID<iPID) MaxNodePolicyID = iPID
          
		  do i = 1, TNodes  !find the node no. in network sequence
			if (iNodeLinkID.eq.NodeID(i)) then
                NodePolicyBeg(iPID,iPType,i)= BegDT     !(policy, PolicyType, Node)
                NodePolicyEnd(iPID,iPType,i)= EndDT     !(policy, PolicyType, Node)
                EXIT
              end if
            end do
          end if
          if (iCompType .eq. 2) then    !for links
          !Kang add for improving performance
          !Decrease Min ID and increase Max ID of current policy 
          if(MinLinkPolicyID>iPID) MinLinkPolicyID = iPID
          if(MaxLinkPolicyID<iPID) MaxLinkPolicyID = iPID
          
            do i = 1, Links   !find the link no. in network sequence
              if (iNodeLinkID.eq.LinkID(i)) then
                LinkPolicyBeg(iPID,iPType,i)= BegDT     !(policy, PolicyType, Node)
                LinkPolicyEnd(iPID,iPType,i)= EndDT     !(policy, PolicyType, Node)
                EXIT
              end if
            end do
          end if

        end if
      end do
888   success = .true.
      !Kang add for improving performance
      !Validate Min ID and Max ID of current policy 
      if(MinNodePolicyID<=0) MinNodePolicyID = 1
      if(MaxNodePolicyID>MaxPolicies) MinNodePolicyID = MaxPolicies
      if(MinLinkPolicyID<=0) MinLinkPolicyID = 1
      if(MaxLinkPolicyID>MaxPolicies) MaxLinkPolicyID = MaxPolicies
      
999   CONTINUE !Kang modify for improving performace CLOSE(UNIT=iDataFile)
      if (.not. success) then
        WRITE(*,*)'PolicyGrp data reading failed:',PolicyFileName
      end if
      return
      end
      
!************************************************************************
      subroutine SetSourceLinkNodes()
!Sets the demand nodes for which each source link is responsible
!Written by Evgenii Matrosov (C) 2011 University College London
      implicit none
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  Local
      INTEGER*4 ln, i,j,n,nn
      logical*1 endsearch
      CHARACTER*30 aVar
!------------------------------------------------------------------------
       do j=1,links
            ndemnodelink(j)=0
       end do
      !Count number of nodes for which links are responsible
      do nn=1,tnodes
        if (DMDNODE(nn)) then
            do j=1,SUPLY_Links_pts(nn)            
                ln=SUPL_Link(j,nn)
                ndemnodelink(ln)=ndemnodelink(ln)+1
            end do
        end if         
      end do 
      
      !Now set the dem nodes for which each source link is responsible
      do ln=1,links
        n=0
        endsearch=.false.
        !if demand link is a source for one or more nodes
        if(ndemnodelink(ln)>0 .and. dmdlink(ln)) then
                do i=1,tnodes
                    !Check if there was a policy change, 
                    if (nodesourcechange(i)==.true.) then
                        do j=1,SUPLY_Links_pts(i)				   
			                if(SUPL_Link(j,i)==ln) then
                                n=n+1
                                LinkSourceNode(ln,n)=i !set nodes for each node
                                if (n >= ndemnodelink(ln)) then 
                                    endsearch=.true.
                                    goto 40    
                                end if    !all dem nodes have been found for link                                 
                            end if
                        end do
                    end if
               end do    
        end if
40      continue
      end do
      	
      return
      end     
!************************************************************************ 
      subroutine SetSourceResNodes()
!Sets the demand nodes for which each source link is responsible
!Written by Evgenii Matrosov (C) 2011 University College London
      implicit none
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  Local
      INTEGER*4 sn, i,j,n,nn
      logical*1 endsearch
      CHARACTER*30 aVar
!------------------------------------------------------------------------
       do j=1,tnodes
            ndemnoderes(j)=0
       end do
       !Count number of nodes for which source reservoirs are responsible
      do nn=1,tnodes
         if (DMDNODE(nn)) then
            do j=1,SUPLY_pts(nn)            
                sn=SUPL_NODE(j,nn)
                if (ResvNode(sn))then
                   ndemnoderes(sn)=ndemnoderes(sn)+1
                end if
            end do
         end if
      end do 
          
      do sn=1,tnodes
        n=0
        endsearch=.false.
        if(ndemnoderes(sn)>0) then
                do i=1,tnodes
                    if (nodesourcechange(i)==.true.) then
                        do j=1,SUPLY_pts(i) 				   
			                if(SUPL_NODE(j,i)==sn) then
                                n=n+1
                                ResSourceNode(sn,n)=i !set nodes for each node
                                if (n >= ndemnoderes(sn)) then 
                                    endsearch=.true.
                                    goto 40    
                                end if    !all source nodes have been found for link                                 
                            end if
                        end do
                    end if
               end do    
        end if
40     continue
      end do
      	
      return
      end     
!************************************************************************ 
**************************NO LONGER CALLED in IRAS 2010 -Evgenii******************************
!      subroutine ReadGageYearUnits(success)
!     read data from gage file (.GAG) (text file):
!	1. Historical flow units to IRAS units (mil m3/day) conversions
!	2. Tells which gauges have flow factors enabled
!  CALL:
!   * readGageUnits
!      implicit none
!      INCLUDE 'IRAS_SYS.INC'
!      INCLUDE 'NODE.INC'
!      INCLUDE 'LINK.INC'
!  INPUT
      !CHARACTER*80 filename
      !COMMON: TNodes, NName(j), NodeID(j)
!  OUTPUT
!      LOGICAL*1 success
      !COMMON: 
      !         SysStat(nGauge1),GageIDObs(),GageUserUnit()
!  Local
!      INTEGER*4 iDatafile
!      CHARACTER*256 aLine
!      CHARACTER*30 aVar
!------------------------------------------------------------------------
!      success = .false.
!      iDatafile = 11
!      OPEN(UNIT=iDatafile, FILE=TRIM(FlowFilename), STATUS='old',
!     &	 FORM='formatted', ERR= 999)
	
!      READ(UNIT=iDatafile,FMT=*,ERR=999)aVar  !skip a line 
!      IF(INDEX(aVar, 'FLOWFILE')== 0) GOTO 999
!      do WHILE (.TRUE.)
!	  READ(UNIT=iDatafile,FMT='(A)',ERR=999) aLine
!        IF (aLine(1:1).ne.'!' ) exit
!      end do

!	Evgenii commented out lines below 0911	
!     read the beginning year and end year
!     sysstat(BYearGage),SYSSTAT(EYearGage) - beginning year and end year
!	READ(aLine, *) sysstat(BYearGage),SYSSTAT(EYearGage)
!	CLOSE(iDatafile)
      
!	call readGageUnits(FlowFileName, 1, success)   !for natural flow
!      IF(.not.success) GOTO 999
	
      !call readGageUnits(FlowFileName, 2, success)   !for added flow
      !IF(.not.success) GOTO 999

!      success = .true.
!999   continue
!      IF(.not.success) WRITE(*,*)'Error in reading gage data!'
!      end subroutine


!************************************************************************
!Kang add for improving performance.
      LOGICAL FUNCTION ReadPolicyOneLine(iLine, aLine)
      USE VARS
	implicit none
      INCLUDE 'IRAS_SYS.INC'
      
!  INPUT
      CHARACTER*256, POINTER::pData(:)
      INTEGER iLine      !output too
!  OUTPUT
      CHARACTER*256 aLine 
      LOGICAL bDataEnd
!------------------------------------------------------------------------
      bDataEnd = .FALSE.
      
      IF(iLine > totalLinePolicyFile) THEN
        bDataEnd = .TRUE.
      ELSE
        aLine = pPolicyData(iLine)
        iLine = iLine + 1
      END IF  
      
      ReadPolicyOneLine = bDataEnd

      END FUNCTION
      