!Copyright (c) 2009,2010 by University College London, Cornell University
!Authors:
!G Pegram, Daniel P. Loucks (dpl3@cornell.edu), Marshall Taylor, Peter French, Huicheng Zhou
!Evgenii Matrosov (evgenii.matrosov@ucl.ac.uk), Julien Harou (j.harou@ucl.ac.uk), Anil Dikshit
!This program is free software under the General Public Licence, GPL (>=v2)
!Read the 'GPL License.txt' file distributed with this source code for a full license statement.
!
!	 *************************************************************************************
!	Modifications before Evgenii (yymmdd):
!	 000621: PNF	Last change 
!	               

C -----------------------------------------------------------------------
      SUBROUTINE Read_Simulation_Data(PolicyGrp,success) !Kang modify 20100630 FileName,PolicyGrp,success)
! Read simulation data for nodes and links from (.INP) file
! for current PolicyGrp & policy
! [Notation]: All input data are converted to the following internal units:
!             Length: m; Area: m^2; Vol: 10^6 m^3; Flow: 10^6 m^3/day
! It calls the following subroutines:
! FOR SYSTEM
!   (*) ReadSysEvaporation()
        !PolicySysEvap(Policy), SysEvap()
! FOR LINKS
!   (*) readLinkRouting()
        !L_Method(),L_NRTRES(LINK), L_CI(), L_CL(), L_CN()
!   (*) ReadLinkRating()
        !LEVAP_PTS(LINK), LINK_FLOW(JJ,LINK), LINK_EVAP(JJ,LINK)
        !LINK_TTV(jj, LINK), LinkWidth(JJ,Link)
!   (*) ReadTransfer()
        !read aquifer data  (bi-directional link)
        !for Policy Data:GWFromVol(IAGMAX*IAGMAX), GWToVol(IAGMAX*IAGMAX),GWFlowFromTo()
        !*** Order by GWFromVol, GWToVol
!   (*) ReadGW()                             !Anthony changed from GWData 240112
        !for Darcy Law: GWMethod(),GWK(),GWElev(),GWLength(),GWWidth()
!   (*) readLinkLoss()
        !LinkLoss(Link),LossMethod
!   (*) ReadCrossSection()
        !LossMethod(Link),BaseWidth(),ChannelDepth()
        !LSlope(), RSlope(),UpLSlope(), UpRSlope()
! FOR NODES
!   (*) ReadGage()
        !GageID(GagMax,iType, NodMax), GageMultiplier(GagMax,iType, NodMax)
!   (*) readAllocation ()
        !read allocation data for a node and its outlinks
        !ALLOC_PTS(NODE), NODE_ALLO(JJ,NODE), NODE_CONSUMP(I,NN)
        !LINK_VTOL(JJ,LINK), jj=1,ALLOC_PTS(NODE)
!   (*) readNodeRating()
        !NARVO_PTS(NODE), NODE_ELEV(JJ,NODE), NODE_AREA(JJ,NODE)
        !            NODE_VOL(JJ,NODE),  NODE_SEEP(JJ,NODE)
        !LAKEQ_PTS(NN), NODE_LAKEQ(I,NN)
!   (*) ReadTarget()         for demand node
        !DMD_TARG(NN), DMD_T_CO(NN)
!   (*) ReadTargetSource()
        !SUPLY_PTS(NN), SUPL_NODE(COUNT,NN), SUPL_FRAC(COUNT,NN)
!   (*) ReadPowerPump
        ! HPCAP(LNKMAX), ECONST(LNKMAX), PLANT_FACTOR(LNKMAX)
        ! HPQMIN(), Intake_elev(),TURBINE_ELEV(), Outlet_Elev()
        ! HPQMax(), HPQMin()
!   (*) ReadReleaseRule()
        !RuleSiteID(RNMAX), Res_Rule(8,ZNMAX,RNMAX)
!   (*) ReadGrpBalance()
	  !GrpVOL_PTS(RNMAX),ResvIDInGrp(RNMAX),GrpVol(*,RNMAX),BalanceVol(*,GRMAX,RNMAX)
!   (*) readNodeEvaporation ()
        !read eveporation data to NODE_EVAP(NODE)
!-------------------------------------------------------------------------
      USE vars
	IMPLICIT NONE
	INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      INTEGER*4 PolicyGrp
      !Kang modify 20100630
      !CHARACTER*80 FileName
      INTEGER*4 inFile
!  OUTPUT
      logical*1 success
      !COMMON: see above
!  Local Variables
      INTEGER*4  i, j, iType, iNode
      !Kang add 20100629
      REAL time_begin, time_end
	REAL*4 timings(10)
C _______________________________________________________________________
      success = .false.
	!Kang add 20100629
	infile=INPFileID
      cRead_Sim_Data = cRead_Sim_Data+1;
      CALL CPU_TIME (time_begin)
      
!     READ LINK DATA
	!call cpu_time(timings(1))
	!write(*,*) 'Now starting timestep ', sysstat(ntimestep)
	DO i = 1,LINKS
        DO iType = 1, MaxPolicyTypes
          !first check if reading data needed
		call cpu_time(timings(1))
          if(.not.LinkPolicyChg(iType,i))cycle
          select case (iType)
            case (Routing0)
              IF(GWLINK(i))cycle
              call readLinkRouting(inFile,PolicyGrp,  
     &           LinkPolicy0(iType,i),LinkID(i),i,success)
              LinkPolicyChg(iType,i) = .false.
            case (Transfer0)
              IF(.not.GWLINK(i))cycle 
!              call ReadTransfer (TRIM(FileName),PolicyGrp, !Commented out because called from within GW calculation
!     &             LinkPolicy0(iType,i),LinkID(i),success)
              call ReadGW(inFile,PolicyGrp,                   ! Anthony changed to GW 230112
     &             LinkPolicy0(iType,i),LinkID(i),i,success)
            case (Rate0)
              IF(GWLINK(i))cycle
              ! first read Evaporation (Loss) data and loss compute_method

			call readLinkLoss(inFile,PolicyGrp,				   !Evgenii 110615 		NodePolicy0 should be replaced with linkpolicy
     &             NodePolicy0(iType,i),LinkID(i),i,success)
              
			if (LossMethod(i)==2) then
			
                		call readCrossSection(inFile,PolicyGrp,		!!Evgenii 110615 		NodePolicy0 should be replaced with linkpolicy
     &                  	NodePolicy0(iType,i),LinkID(i),i,success)
             		 else
                		call readLinkRating(inFile,PolicyGrp,
     &          LinkPolicy0(iType,i),LinkID(i),i,LossMethod(i),success)
              		end if
              LinkPolicyChg(iType,i) = .false.
		  case (Cost0)
              IF(GWLINK(i))cycle 
              call ReadCostLink(inFile,PolicyGrp,
     &             LinkPolicy0(iType,i),LinkID(i),i,success)
          end select
		call cpu_time(timings(2))
        END DO
	END do
!	call cpu_time(timings(2))
!     READ NODE DATA
      DO i = 1,TNODES
        do iType = 1, MaxPolicyTypes
          !first check whether reading data needed
          if(.not.NodePolicyChg(iType,i))cycle
          select case (iType)
!            case (Gage0)		Evgenii took out Gage reading 090923
!              call readGage(inFile,PolicyGrp,
!     &             NodePolicy0(iType,i),NodeID(i),i,1,success)
!              call readGage(inFile,PolicyGrp,
!     &             NodePolicy0(iType,i),NodeID(i),i,2,success)
!              NodePolicyChg(iType,i) = .false.
              case (Allocation0)
              call readAllocation(inFile,PolicyGrp,
     &             NodePolicy0(iType,i),NodeID(i),i,success)
              NodePolicyChg(iType,i) = .false.
            case (Rate0)
              IF(CapN(i).le.0.0 .and. .not. GWNODE(i))cycle !Evgenii added  ".and. .not. GWNODE(i)" because GWnodes do not have a capactiy
              call readNodeRating(inFile,PolicyGrp,
     &             NodePolicy0(iType,i),NodeID(i),i,success)
              NodePolicyChg(iType,i) = .false.
            case (Power0)
              IF(PowerNode(i)==0)cycle
              call ReadPowerPump (inFile,PolicyGrp,
     &           NodePolicy0(iType,i),NodeID(i),PowerNode(i),1,success)
              NodePolicyChg(iType,i) = .false.
            case (Pump0)
              IF(PumpNode(i)==0)cycle
              call ReadPowerPump (inFile,PolicyGrp,
     &           NodePolicy0(iType,i),NodeID(i),PumpNode(i),0,success)
              NodePolicyChg(iType,i) = .false.
            case (Target0)
              IF(.not.DMDNODE(i))cycle
              call ReadTarget(inFile,PolicyGrp,
     &             NodePolicy0(iType,i),NodeID(i),i,success)
              NodePolicyChg(iType,i) = .false.
            case (Source0)
              IF(.not.DMDNODE(i))cycle
              call ReadTargetSource(inFile,PolicyGrp,
     &             NodePolicy0(iType,i),NodeID(i),i,success)
              NodePolicyChg(iType,i) = .false.
              nodesourcechange(i)=.true. !Evgenii added 110918 for link sources optimization
              
            case (Evaporation0)
              call readNodeEvaporation(inFile,PolicyGrp,
     &             NodePolicy0(iType,i),NodeID(i),i,success)
              NodePolicyChg(iType,i) = .false.
	!Evgenii added performance measures 100610
	      case (Performance0)
              IF(.not.DMDNODE(i))cycle
			call readPerformance(inFile,PolicyGrp,
     &             NodePolicy0(iType,i),NodeID(i),i,success)
              NodePolicyChg(iType,i) = .false.
          	 !Evgenii added demand reduction measures 100720
	      case (DemRed0)
              IF(.not.DMDNODE(i))cycle
			call readDemRed(inFile,PolicyGrp,
     &             NodePolicy0(iType,i),NodeID(i),i,success)
              NodePolicyChg(iType,i) = .false.
          end select
        end do
      END do
!	call cpu_time(timings(3))
!     Read reservoir release rules for a group of reservoirs
!     Then read Balance data
      do i = 1, nGrpResv
        do j = 1, TNodes
          if (NodeID(j)==RuleSiteID(i))then
            iNode = j; exit
          end if
        end do
        IF(NodePolicyChg(Rule0,iNode)) then
    
		call ReadReleaseRule(inFile,PolicyGrp,
     &         NodePolicy0(Rule0,iNode),i, RuleSiteID(i), success)
          NodePolicyChg(Rule0,iNode) = .false.
        END if
        IF(NodePolicyChg(Balance0,iNode)) then
		call ReadGrpBalance(inFile,PolicyGrp,
     &         NodePolicy0(Balance0,iNode),i, RuleSiteID(i), success)
          NodePolicyChg(Balance0,iNode) = .false.
        END if
	
      end do
      success = .true.
 !     call cpu_time(timings(4))
      !Kang add 20100629
!      CALL CPU_TIME (time_end)
!      tRead_Sim_Data = tRead_Sim_Data +(time_end-time_begin)
	!write(*,fmt=9) timings(2)-timings(1), timings(3)-timings(2),
!   & timings(4)-timings(3)!, timings(5)-timings(4)
9     FORMAT(10E12.3)	  

999   RETURN
      END
C
!******************FOR SYSTEM*********************************************
      subroutine readSysEvaporation(success)
!  read default evaporation for surface nodes from iras.inp
!  these are default seasonal evaporation rates which are repeated for every year
!  these are therefore "policies" and not "policy groups" as defined in IRAS
!  CALL: None
      USE vars
	implicit none
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      CHARACTER*80  filename
!  OUTPUT
      logical*1     success
      !COMMON:      SysEvap(*)
!  Local variables
      INTEGER*4     i, iPID, iDatafile
      CHARACTER*256 aLine
      CHARACTER*20  aVar

!------------------------------------------------------------------------
      !Kang add 20100629
      cSysEvap = cSysEvap +1
      
      success = .true.  !because of initializing
      do i = 1, MaxPolicies
        SysEvap(i) = 0.0
      end do
      iDatafile = INPFileID
	rewind(iDatafile)!Evgenii 10115 took out opening of INP file and replaced with rewind
!      OPEN(UNIT=iDatafile, FILE=TRIM(SysFilename), STATUS='old',
!     &	 FORM='formatted', ERR= 999)
      !find SysEvaporation: line
      i = 0
      do WHILE (.TRUE.)
        READ(UNIT=iDatafile,FMT='(A)',ERR=999, END=999) aLine
        if (Index(aLine,'SysEvaporation:')==1) then
          i = i + 1
          read(aLine,*)aVar, iPID
          if (iPID > MaxPolicies ) exit
          read(aLine,*)aVar, iPID, SysEvap(iPID)
          
		!unit conversion to m
		call UnitConversion(1, ULoss, SysEvap(iPID))
        else
          if (i > 0 ) exit
        END if
      end do
	
!999   CLOSE(UNIT=iDataFile)
999   REWIND(UNIT=iDataFile)!Evgenii 10115 replaced close with rewind

      return
      end
	
!************************** FOR LINKS ***********************************
      Subroutine  ReadLinkRating( iDatafile,PolicyGrp, Policy
     &   ,IDLink, LINK, iLossMethod,success)
! For a link:
! Read data of flow and lost to common variables.
!  CALL: None
      USE vars
	implicit none
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      !Kang modify 20100630
      !CHARACTER*80 FileName
      INTEGER*4 PolicyGrp, Policy, IDLINK, LINK,iLossMethod
!  OUTPUT
      LOGICAL*1 success
      !COMMON:
      ! LEVAP_PTS(LINK),  LINK_FLOW(JJ,LINK),  LINK_EVAP(JJ,LINK)
      ! TTV_PTS(LINK), LINK_TTV(jj,LINK)
!  Local
      INTEGER*4 i, n, GroupRead, IDLinkRead, iDatafile, iPos
      INTEGER*4 CompType, CompTypeRead, PolicyRead
      CHARACTER*30 aVar
      CHARACTER*256 aLine
      INTEGER*4 L
!------------------------------------------------------------------------
      success = .false.
      !Kang add 20100629
      cLinkRating = cLinkRating +1
      CompType = 2  !2-for link. 1-for node
      !Kang modify 20100630
!      iDatafile = 11
!      OPEN(UNIT=iDatafile, FILE=TRIM(Filename), STATUS='old',
!     &	 FORM='formatted', ERR= 999)
      n = 0
      IF(bUseBufferedBinaryData) THEN
        DO L=1, nRating
          IF(pRating(L)%GroupID==PolicyGrp .and.
     &          pRating(L)%Policy==Policy .and.
     &          pRating(L)%CompType==CompType .and.
     &          pRating(L)%ID==IDLink)THEN
            n = n + 1
            iflinkloss(link)=.true. !evgenii added to know if loss enabled on link 090721
            if (iLossMethod == 0) then  
                LINK_EVAP(n,LINK)=pRating(L)%AreaOrEvapor   !aVar,aVar, Evgenii took out 2 avars
                LINK_FLOW(n,LINK)=pRating(L)%VolOrFlow
                CALL UnitConversion(2, UFlow, LINK_EVAP(n,LINK) ) 
                CALL UnitConversion(2, UFlow, LINK_FLOW(n,LINK) ) 
!               CALL UnitConversion(2, UTime, LINK_TTV(n,LINK) ) Evgenii commented out 090702
		    else                        !get flow, time of flow, width    
 			  !Evgenii took out 4 avar to clean up inp 090810
 			  !&  LINK_TTV(n,LINK) Evgenii commented out 090702
              LinkWidth(n,LINK)=pRating(L)%ElevOrWidth
              LINK_FLOW(n,LINK)=pRating(L)%VolOrFlow
												         
              !unit conversion
              CALL UnitConversion(2, ULen, LinkWidth(n,LINK) )   !m
!             CALL UnitConversion(2, UTime, LINK_TTV(n,LINK) )  !day Evgenii commented out 090702
              CALL UnitConversion(2, UFlow, LINK_FLOW(n,LINK))   !10^6 m^3/day  
	      end if
            if (n >= IAGMAX ) exit
          END if
       END DO
      ELSE
      IF(nBufferedLines==0) REWIND(UNIT=iDatafile)
      
      !find line 'Rating:...'
	
      !Kang modify for improving performance 
      L = 1
      IF(Rating_BL>0) L = Rating_BL
      do WHILE (.TRUE.)
        IF(nBufferedL   ines>0) THEN
            IF(L>nBufferedLines) EXIT
            aLine = pFileData(L)
            L = L + 1
        ELSE    
        READ(UNIT=iDatafile,FMT='(A)',ERR=999, END=999) aLine
        END IF

        if (Index(aLine, headRating) >= 1) then
          IF(Rating_BL<0) Rating_BL = L-1  
          iPos = INDEX(aLine, ':')
          read(aLine(iPos+1:),*)GroupRead, PolicyRead,
     &         CompTypeRead, IDLinkRead
          IF(GroupRead==PolicyGrp.and.PolicyRead==Policy.and.
     &      CompTypeRead==CompType.and.IDLinkRead==IDLink)then
            n = n + 1
            iflinkloss(link)=.true. !evgenii added to know if loss enabled on link 090721
            if (iLossMethod == 0) then  
              read(aLine(iPos+1:),*)aVar, aVar, aVar, aVar,
     &          aVar,LINK_EVAP(n,LINK),LINK_FLOW(n,LINK)       !Evgenii took out 4 avar to clean up inp 090810
     &          												 !LINK_TTV(n,LINK) Evgenii took out time, not used
       
              CALL UnitConversion(2, UFlow, LINK_EVAP(n,LINK) ) 
              CALL UnitConversion(2, UFlow, LINK_FLOW(n,LINK) ) 
!             CALL UnitConversion(2, UTime, LINK_TTV(n,LINK) ) Evgenii commented out 090702
            
		  else                        !get flow, time of flow, width
              read(aLine(iPos+1:),*)aVar, aVar, aVar, aVar,     !Evgenii took out 4 avar to clean up inp 090810
     &             LinkWidth(n,LINK),aVar,LINK_FLOW(n,LINK)     !&  LINK_TTV(n,LINK) Evgenii commented out 090702
 															         
              !unit conversion
              CALL UnitConversion(2, ULen, LinkWidth(n,LINK) )   !m
!             CALL UnitConversion(2, UTime, LINK_TTV(n,LINK) )  !day Evgenii commented out 090702
              CALL UnitConversion(2, UFlow, LINK_FLOW(n,LINK))   !10^6 m^3/day  
	   end if
            if (n >= IAGMAX ) exit
          END if
        else
          IF(n > 0) exit
        END if
      end do
      END IF
      
      success = .true.
999   if (iLossMethod==0.or.iLossMethod==1) LEVAP_PTS(LINK) = n
      TTV_PTS(LINK) = n
      !Kang modify 20100630 CLOSE(UNIT=iDataFile)      
      return
      end


!*************************************************************************
      subroutine readLinkLoss(iDatafile,GroupID,PolicyID,
     &     	     IDLink, iReadSeq, success)
!  read gage contribution multipliers for a node
!  CALL: None
      USE vars
	implicit none
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      !Kang modify 20100630
      !CHARACTER*80 FileName
      INTEGER*4     GroupID, PolicyID, iReadSeq, IDLink
!  OUTPUT
      logical*1     success
      !COMMON:      LinkLoss(iReadSeq), LossMethod()
!  Local variables
      INTEGER*4     i, j, iPos, iDatafile, comp_type
      INTEGER*4     GroupIDRead, PolicyIDRead, IDLinkRead
      CHARACTER*256 aLine
      CHARACTER*20  aVar
      INTEGER*4 L

!------------------------------------------------------------------------
      success = .false.
      !Kang add 20100629
      cLinkLoss = cLinkLoss + 1
      !Kang modify 20100630 
!      iDatafile = 11
!      OPEN(UNIT=iDatafile, FILE=TRIM(filename), STATUS='old',
!     &	 FORM='formatted', ERR= 999)
      IF(bUseBufferedBinaryData) THEN      
          DO L=1, nEvapor
              IF(pEvapor(L)%GroupID==GroupID
     &          .and.pEvapor(L)%Policy==PolicyID
     &          .and.pEvapor(L)%CompType==2
     &          .and.pEvapor(L)%ID==IDLink)then
                LinkLoss(iReadSeq) = pEvapor(L)%Evaporation
                LossMethod(iReadSeq) = pEvapor(L)%LossMethod    !Evgenii - Link loss is for loss methods 1 and 2 
		        !unit conversion
                CALL UnitConversion(2, ULoss, LinkLoss(iReadSeq) )   !m/day
                exit
              END IF
          END DO
      ELSE
      IF(nBufferedLines==0) REWIND(UNIT=iDatafile)
      
      !find Evaporation: line
      !Kang modify for improving performance
      L = 1
      IF(Evaporation_BL>0) L = Evaporation_BL      
      do WHILE (.TRUE.)
        IF(nBufferedLines>0) THEN
            IF(L>nBufferedLines) EXIT
            aLine = pFileData(L)
            L = L + 1
        ELSE    
        READ(UNIT=iDatafile,FMT='(A)',ERR=999, END=999) aLine
        END IF

        if (Index(aLine,headEvaporation)==1) then
         
        IF(Evaporation_BL<0) Evaporation_BL = L - 1
              
		iPos = INDEX(aLine, ':')
	
          read(aLine(iPos+1:),*)GroupIDRead,PolicyIDRead,Comp_type,
     &         IDLinkRead
          IF(GroupIDRead == GroupID.and.PolicyIDRead==PolicyID
     &        .and. Comp_type==2 .and.IDLinkRead==IDLink) then
            read(aLine(iPos+1:),*)aVar, aVar, aVar,aVar,
     &                LinkLoss(iReadSeq),LossMethod(iReadSeq) !Evgenii - Link loss is for loss methods 1 and 2 
		  !unit conversion
            CALL UnitConversion(2, ULoss, LinkLoss(iReadSeq) )   !m/day
            exit
          END if
        END if
      end do
      END IF
	
      success = .true.
999   continue !Kang modify 20100630    CLOSE(UNIT=iDataFile)
      return
      end

!*************************************************************************
      subroutine readCrossSection(iDatafile,GroupID,PolicyID,
     &     	     IDLink, iReadSeq, success)
!  read gage contribution multipliers for a node
!  CALL: None
      USE vars
	implicit none
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      CHARACTER*80  filename
      INTEGER*4     GroupID, PolicyID, iReadSeq, IDLink
!  OUTPUT
      logical*1     success
      !COMMON:      BaseWidth(iReadSeq),ChannelDepth(),LSlope(),RSlope()
      !             UpLSlope(),UpRSlope()
!  Local variables
      INTEGER*4     i, j, iPos, iDatafile, comp_type
      INTEGER*4     GroupIDRead, PolicyIDRead, IDLinkRead
      CHARACTER*256 aLine
      CHARACTER*20  aVar
      INTEGER*4 L

!------------------------------------------------------------------------
      success = .false.
      !Kang add 20100629
      cCrossSection = cCrossSection + 1
      !Kang modify 20100630 
!      iDatafile = 11
!      OPEN(UNIT=iDatafile, FILE=TRIM(filename), STATUS='old',
!     &	 FORM='formatted', ERR= 999)
      IF(bUseBufferedBinaryData) THEN      
          DO L=1, nCross
              IF(pCross(L)%GroupID==GroupID
     &          .and.pCross(L)%Policy==PolicyID
     &          .and.pCross(L)%CompType==2
     &          .and.pCross(L)%LinkID==IDLink)then
                BaseWidth(iReadSeq) = pCross(L)%BaseWidth
                ChannelDepth(iReadSeq) = pCross(L)%ChannelDepth
                LSlope(iReadSeq) = pCross(L)%LSlope
                RSlope(iReadSeq) = pCross(L)%RSlope
                UpLSlope(iReadSeq) = pCross(L)%UpLSlope
                UpRSlope(iReadSeq) = pCross(L)%UpRSlope
     
                CALL UnitConversion(2, ULen, BaseWidth(iReadSeq) )   !m
                CALL UnitConversion(2, ULen, ChannelDepth(iReadSeq) )   !m
                exit
              END IF
          END DO
      ELSE
      IF(nBufferedLines==0) REWIND(UNIT=iDatafile)
      
      !find CrossSection: line
      !Kang modify for improving performance
      L = 1
      IF(CrossSection_BL>0) L = CrossSection_BL      
      do WHILE (.TRUE.)
        IF(nBufferedLines>0) THEN
            IF(L>nBufferedLines) EXIT
            aLine = pFileData(L)
            L = L + 1
        ELSE    
        READ(UNIT=iDatafile,FMT='(A)',ERR=999, END=999) aLine
        END IF
        
        if (Index(aLine,headCross)==1) then
          IF(CrossSection_BL<0) CrossSection_BL = L - 1   
          iPos = INDEX(aLine, ':')
          read(aLine(iPos+1:),*)GroupIDRead,PolicyIDRead,Comp_type,
     &         IDLinkRead
		IF(GroupIDRead == GroupID.and.PolicyIDRead==PolicyID
     &        .and. Comp_type==2 .and.IDLinkRead==IDLink) then
            read(aLine(iPos+1:),*)aVar, aVar, aVar,aVar,
     &         BaseWidth(iReadSeq),ChannelDepth(iReadSeq),
     &         LSlope(iReadSeq),RSlope(iReadSeq),
     &         UpLSlope(iReadSeq),UpRSlope(iReadSeq)
            CALL UnitConversion(2, ULen, BaseWidth(iReadSeq) )   !m
            CALL UnitConversion(2, ULen, ChannelDepth(iReadSeq) )   !m
            exit
          END if
        END if
      end do
	END IF
	
      success = .true.
999   continue !Kang modify 20100630   CLOSE(UNIT=iDataFile)
      return
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Subroutine readLinkRouting(iDatafile, PolicyGrp,Policy,
     &           IDLink, iReadSeq, success)
! read link routing number and its parms from IRAS.inp
!  CALL: None
      USE vars
	implicit none
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      character*80 FileName
      INTEGER*4 PolicyGrp,Policy,IDLink, iReadSeq
!  OUTPUT
      LOGICAL*1 success
      !COMMON: L_NRTRES(LINK), L_CI(), L_CL(), L_CN()
!  local
      INTEGER*4 i, iDatafile, iPos, iMethod, NumParms, CompType
      INTEGER*4 GroupRead, PolicyRead, IDLinkRead, CompTypeRead
      CHARACTER*30 aVar
      CHARACTER*256 aLine
      INTEGER*4 L
!------------------------------------------------------------------------
      success = .true.    !always true because of default values
      !Set defaults for method 1
      CompType = 2
      L_method(iReadSeq)=0; L_NRTRES(iReadSeq) = 0 !Evgenii, 100305 changed default L_method(iReadSeq) to 0 instead of 1
      L_a(iReadSeq) = 1; L_b(iReadSeq) = 1; L_C(iReadSeq) = 0 !Evgenii 090703 changed L_C to 0 from 1 so default is no link volume
      !Kang add 20100629
      cLinkRouting = cLinkRouting +1
      
      !Kang modify 20100630 
!      iDatafile = 11
!      OPEN(UNIT=iDatafile, FILE=TRIM(filename), STATUS='old',
!     &	 FORM='formatted', ERR= 999)
      IF(bUseBufferedBinaryData) THEN
        DO L=1, nRouting
          IF(pRouting(L)%GroupID==PolicyGrp .and.
     &          pRouting(L)%Policy==Policy .and.
     &          pRouting(L)%CompType==CompType .and.
     &          pRouting(L)%LinkID==IDLink)THEN
                
              L_method(iReadSeq)=pRouting(L)%iMethod  !Anthony fixed reading in L_method 120116
              IF (pRouting(L)%iMethod == 2) THEN   !cascading reservoirs:(a*inflow+b*Vol)^c
                !Evgenii took out one avar to clean up inp
                !a-flow; b-vol; c-exponent                
                L_NRTRES(iReadSeq)=pRouting(L)%L_NRTRES        !Anthony fixed reading in L_NRTRES 120116
                L_a(iReadSeq) = min(1.,pRouting(L)%L_a)
                L_b(iReadSeq) = min(1.,pRouting(L)%L_b)
                L_c(iReadSeq) = min(1.,pRouting(L)%L_c)        !Evgenii - Volume and flow are not converted, user units used
                exit
                !Evgenii 100305 put iMethod =1 condition below
		      ELSE IF (pRouting(L)%iMethod == 1)THEN                      !default method: a(volume-detention)^b
                !Evgenii took out avar
                !a-liner parm;b-exponent;c-detention storage
                L_a(iReadSeq) = min(1.,pRouting(L)%L_a)
                L_b(iReadSeq) = pRouting(L)%L_b
                L_C(iReadSeq) = pRouting(L)%L_c 
			    CALL UnitConversion(2, UVol, L_C(iReadSeq))
                exit
            end if
          END if
       END DO
      ELSE

      IF(nBufferedLines==0) REWIND(UNIT=iDatafile)
      
      !find line 'Routing:...'
	
      !Kang modify for improving performance
      L = 1
      IF(Routing_BL>0) L = Routing_BL      
      do WHILE (.TRUE.)
        IF(nBufferedLines>0) THEN
            IF(L>nBufferedLines) EXIT
            aLine = pFileData(L)
            L = L + 1
        ELSE    
        READ(UNIT=iDatafile,FMT='(A)',ERR=999, END=999) aLine
        END IF
       
	  if (Index(aLine, headRouting) >= 1) then
	    IF(Routing_BL<0) Routing_BL = L - 1
		iPos = INDEX(aLine, ':')
          read(aLine(iPos+1:),*)GroupRead, PolicyRead,CompTypeRead, !comptype->1 node, 2 link
     &          IDLinkRead, iMethod !,NumParms !NumParms is not used in code, Evgenii
        	linkroute(idlinkread)=imethod 
		IF(GroupRead==PolicyGrp.and.PolicyRead==Policy.and.
     &       CompTypeRead==CompType.and.IDLinkRead==IDLink) then
			
            if (iMethod == 2) then   !cascading reservoirs:(a*inflow+b*Vol)^c
              read(aLine(iPos+1:),*)aVar,aVar,aVar,aVar,
     &          L_method(iReadSeq),L_NRTRES(iReadSeq),     !Evgenii took out one avar to clean up inp
     &          L_a(iReadSeq), L_b(iReadSeq),L_C(iReadSeq) !a-flow; b-vol; c-exponent
              L_a(iReadSeq) = min(1.,L_a(iReadSeq))
              L_b(iReadSeq) = min(1.,L_b(iReadSeq))
              L_c(iReadSeq) = min(1.,L_c(iReadSeq))        !Evgenii - Volume and flow are not converted, user units used
              exit
            !Evgenii 100305 put iMethod =1 condition below
		  else if (iMethod == 1)then                      !default method: a(volume-detention)^b
              
			read(aLine(iPos+1:),*)aVar,aVar,aVar,aVar,
     &          L_method(iReadSeq),				!Evgenii took out avar			
     &          L_a(iReadSeq),L_b(iReadSeq),L_C(iReadSeq) !a-liner parm;b-exponent;c-detention storage
		 	L_a(iReadSeq) = min(1.,L_a(iReadSeq))
			CALL UnitConversion(2, UVol, L_C(iReadSeq))
              exit
            end if
          END if
        END if
      end do
      END IF
	
999   continue !Kang modify 20100630    CLOSE(UNIT=iDataFile)
      return
      end

!******************************************************************
      Subroutine  ReadTransfer (PolicyGrp, Policy,
     &            IDLink, link, success)
! For a bi-directional link:
! Read transfer data
!  CALL: None
      USE vars
	implicit none
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      INTEGER*4 PolicyGrp, Policy, IDLINK, link
!  OUTPUT
      LOGICAL*1 success
      !COMMON:
      !GWFromVol(IAGMAX*IAGMAX), GWToVol(IAGMAX*IAGMAX),GWFlowFromTo()
!  Local
      INTEGER*4 n , iDatafile, iPos
      INTEGER*4 PolicyRead, GroupRead, IDLinkRead
      CHARACTER*30 aVar
      CHARACTER*256 aLine
      INTEGER*4 iFromVol, iToVol
      INTEGER*4 L
!------------------------------------------------------------------------
      success = .false.
      !for unit conversion
      !check whether GWFromVol(n) and GWToVol(n) have flow unit if non-storage
      iFromVol = UVol; iToVol = UVol
      if (CapN(Nin(Link)) > 0 .OR. GWNODE(Nin(Link)) ) then !Evgenii put in .OR. GWNODE(NOut(Link) because GWNodes dont have capacity
        iFromVol = UVol
      else
        iFromVol = UFlow
      end if
      if (CapN(NOut(Link))> 0 .OR. GWNODE(NOut(Link))  ) then !Evgenii put in .OR. GWNODE(NOut(Link) because GWNodes dont have capacity
        iToVol = UVol
      else
        iToVol = UFlow
      end if
      n = 0
      
      IF(nBufferedLines>0) THEN      				!Anthony added 150112 to read transfer data to memory
          DO L=1, nTransfer
              IF(pTransfer(L)%LinkID==IDLink) then  ! Anthony changed the mechanism here as LinkID more likely to vary than PolicyGrp/Policy
                  if (pTransfer(L)%GroupID==PolicyGrp
     &                .and.pTransfer(L)%Policy==Policy)then
                      n = n + 1
                      GWFromVol(n) = pTransfer(L)%GWFromVol
                      GWToVol(n) = pTransfer(L)%GWToVol
			        GWFlowFromTo(n)=pTransfer(L)%GWFlowFromTo
                      IF (GWFromVol(n)>0) THEN
                          CALL UnitConversion(2, iFromVol, GWFromVol(n))
                      END IF
                      IF (GWToVol(n)>0) THEN
                          CALL UnitConversion(2, iToVol, GWToVol(n))
                      END IF
                      IF (GWFlowFromTo(n)>0) THEN
                          CALL UnitConversion(2, UFlow, GWFlowFromTo(n))
                      endif
                  ENDIF
              END IF
          END DO
      ELSE

      IF(nBufferedLines==0) REWIND(UNIT=iDatafile)
      !Kang add 20100629
      cTransfer = cTransfer +1
      
      ! iDatafile = 11
      ! OPEN(UNIT=iDatafile, FILE=TRIM(sysfilename), STATUS='old',
      !&	 FORM='formatted', ERR= 999)
      
      n = 0
      !find line 'Transfer:...'
      DO WHILE (.TRUE.)
        READ(UNIT=iDatafile,FMT='(A)',ERR=999, END=999) aLine
        if (Index(aLine, 'Transfer:') >= 1) then
          iPos = INDEX(aLine, ':')
          read(aLine(iPos+1:),*)GroupRead, PolicyRead,IDLinkRead
          IF(GroupRead==PolicyGrp.and.PolicyRead==Policy.and.
     &       IDLinkRead==IDLink)then
            n = n + 1
            read(aLine(iPos+1:),*)aVar,aVar,aVar,GWFromVol(n),
     &          GWToVol(n),GWFlowFromTo(n)
		  !write(*,*)'ReadTransfer textread, link ',IDLink	
            !unit conversion
            CALL UnitConversion(2, iFromVol, GWFromVol(n))
            CALL UnitConversion(2, iToVol, GWToVol(n))
            CALL UnitConversion(2, UFlow, GWFlowFromTo(n)) 
          END if
        else
          if (n>0) exit
        END if
      END DO
      END IF

      success = .true.
999   GWALO_PTS(Link) = n
      !Kang modify 20100630    CLOSE(UNIT=iDataFile)
      close (iDatafile) 
      return
      end

!******************************************************************
	Subroutine  ReadGW ( iDatafile,PolicyGrp, Policy,
     &            IDLink, link, success)
! For a bi-directional link:
! Read transfer data
!  CALL: None
      USE vars
	implicit none
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      character*80 FileName
      INTEGER*4 PolicyGrp, Policy, IDLINK, link
!  OUTPUT
      LOGICAL*1 success
      !COMMON:
      !GWMethod(Link),GWK(),GWElev(),GWLength(),GWWidth()
!  Local
      INTEGER*4 i, j, iDatafile, iPos
      INTEGER*4 PolicyRead, GroupRead, IDLinkRead
      CHARACTER*30 aVar
      CHARACTER*256 aLine
      INTEGER*4 L
!------------------------------------------------------------------------
      success = .false.
      
      !Kang add 20100629
      cGW = cGW + 1
      
      IF(nBufferedLines>0) THEN      
          DO L=1, nGW
              IF(pGW(L)%GroupID==GroupRead
     &          .and.pGW(L)%Policy==PolicyRead
     &          .and.pGW(L)%LinkID==IDLink)then
                GWMethod(link) = pGW(L)%GWMethod
                GWK(link) = pGW(L)%GWK
			  GWElev(link)=pGW(L)%GWElev
                GWLength(link)=pGW(L)%GWLength
                GWWidth(link)=pGW(L)%GWWidth
                CALL UnitConversion(2, UK, GWK(link) )   !m/day
                CALL UnitConversion(2, ULen, GWElev(link) )   !m
                    if (GWMethod(Link)==2) then   !GW-GW vertical: common area
                      CALL UnitConversion(2, UArea, GWLength(link))   !m^2
                    else
                      CALL UnitConversion(2, ULen, GWLength(link))   !m
                    end if
                CALL UnitConversion(2, ULen, GWWidth(link) )
                exit
              END IF
          END DO
      ELSE
      !Kang modify 20100630 
!      iDatafile = 11
!      OPEN(UNIT=iDatafile, FILE=TRIM(filename), STATUS='old',
!     &	 FORM='formatted', ERR= 999)
      IF(nBufferedLines==0) REWIND(UNIT=iDatafile)
      	
      !find line 'Groundwater:...'
      !Kang modify for improving performance
      L = 1
      IF(GW_BL>0) L = GW_BL      
      do WHILE (.TRUE.)
        IF(nBufferedLines>0) THEN
            IF(L>nBufferedLines) EXIT
            aLine = pFileData(L)
            L = L + 1
        ELSE    
        READ(UNIT=iDatafile,FMT='(A)',ERR=999, END=999) aLine
        END IF
       
        if (Index(aLine, 'Groundwater:') >= 1) then
          IF(GW_BL<0) GW_BL = L - 1      
          iPos = INDEX(aLine, ':')
          read(aLine(iPos+1:),*)GroupRead, PolicyRead,aVar,IDLinkRead
          IF(GroupRead==PolicyGrp.and.PolicyRead==Policy.and.
     &       IDLinkRead==IDLink)then
            read(aLine(iPos+1:),*)aVar,aVar,aVar,aVar,
     &          GWMethod(Link),GWK(Link),GWElev(Link),
     &          GWLength(Link),GWWidth(Link)
            !unit conversion
			!write(*,*)'ReadGW textread, link ',IDLink	
		  CALL UnitConversion(2, UK, GWK(Link) )   !m/day
            CALL UnitConversion(2, ULen, GWElev(Link) )   !m
            if (GWMethod(Link)==2) then   !GW-GW vertical: common area
              CALL UnitConversion(2, UArea, GWLength(Link) )   !m^2
            else
              CALL UnitConversion(2, ULen, GWLength(Link) )   !m
            end if
            CALL UnitConversion(2, ULen, GWWidth(Link) )   !m
            exit
          END if
        END if
      END do
      END IF
	
      success = .true.
999   continue
      !Kang modify 20100630    CLOSE(UNIT=iDataFile)
      return
      end

!*************************************************************************
      subroutine readCostlink(iDatafile,GroupID,PolicyID,
     &     	     IDLink, iReadSeq, success)
!  read cost and energy of flow through link
!  CALL: None
      USE vars
	implicit none
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      CHARACTER*80  filename
      INTEGER*4     GroupID, PolicyID, iReadSeq, IDLink
!  OUTPUT
      logical*1     success
      !COMMON:      
      !             
!  Local variables
      INTEGER*4     i, j, iPos, iDatafile
      INTEGER*4     GroupIDRead, PolicyIDRead, IDLinkRead
      CHARACTER*256 aLine
      CHARACTER*20  aVar
      INTEGER*4 L

!------------------------------------------------------------------------
      success = .false.
      
      !Kang modify for improving performance
!      iDatafile = 11
!      OPEN(UNIT=iDatafile, FILE=TRIM(filename), STATUS='old',
!     &	 FORM='formatted', ERR= 999)
	!write(*,*) 'ReadCostLink nBuff ', nBufferedLines
      IF(nBufferedLines>0) THEN      
          DO L=1, nCost
              IF(pCost(L)%GroupID==GroupID
     &          .and.pCost(L)%Policy==PolicyID
     &          .and.pCost(L)%LinkID==IDLink)then
                FlowCost(iReadSeq) = pCost(L)%FlowCost
                FlowEng(iReadSeq) = pCost(L)%FlowEng
			  AnnCostInc(iReadSeq)=pCost(L)%AnnCostInc
                exit
              END IF
          END DO
      ELSE

	!JRK added this newline

      IF(nBufferedLines==0) REWIND(UNIT=iDatafile)
      
      !find Cost: line
	!write(*,*) 'Cost_BL ', Cost_BL 
      !Kang modify for improving performance
      L = 1
      IF(Cost_BL>0) L = Cost_BL      
      do WHILE (.TRUE.)
        IF(nBufferedLines>0) THEN
            IF(L>nBufferedLines) EXIT
            aLine = pFileData(L)
            L = L + 1
        ELSE    
        READ(UNIT=iDatafile,FMT='(A)',ERR=999, END=999) aLine
        END IF

        if (Index(aLine,headCost)==1) then
          IF(Cost_BL<0) Cost_BL = L - 1
          iPos = INDEX(aLine, ':')
          read(aLine(iPos+1:),*)GroupIDRead,PolicyIDRead,
     &         IDLinkRead
		!write(*,*)'ReadCostLink textread, link ',IDLink
		IF(GroupIDRead == GroupID.and.PolicyIDRead==PolicyID
     &        .and.IDLinkRead==IDLink) then
            read(aLine(iPos+1:),*)aVar, aVar,aVar,
     &         FlowCost(iReadSeq),FlowEng(iReadSeq),AnnCostInc(iReadSeq)
            exit
          END if
        END if
      end do
      END IF
      
      success = .true.
999   continue
      !Kang modify 20100630    CLOSE(UNIT=iDataFile)
      return
      end
    
!************************ FOR NODES *************************************

      Subroutine  readNodeRating ( iDatafile, PolicyGrp, Policy,
     &    IDNode, Node, success)
! For a node:
! read data of stage(elevation), Area, Storage, Seepage (loss) to common variables.
!  CALL: None
      USE vars
	implicit none
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      character*80 FileName
      INTEGER*4 PolicyGrp, Policy, IDNode, Node
      ! PolicyGrp--related to simulation year
      ! Node--changes from 1 to TNodes determined by node reading sequence
!  OUTPUT
      LOGICAL*1 success
      !COMMON:
      ! NARVO_PTS(NODE),NODE_ELEV(JJ,NODE),NODE_AREA(JJ,NODE),NODE_VOL(JJ,NODE),NODE_SEEP(JJ,NODE)
      ! LAKEQ_PTS(NN), NODE_LAKEQ(I,NN)
!  local
      INTEGER*4 i, n , GroupRead, IDNodeRead, iDatafile, iPos
      INTEGER*4 CompType, CompTypeRead, PolicyRead
      CHARACTER*30 aVar
      CHARACTER*256 aLine
      INTEGER*4 L
!------------------------------------------------------------------------
      success = .false.
      CompType = 1  !1-for node; 2-for link
      !Kang add 20100629
      cNodeRating = cNodeRating +1
 
      n = 0
      !Kang modify 20100630 
!      iDatafile = 11
!      OPEN(UNIT=iDatafile, FILE=TRIM(filename), STATUS='old',
!     &	 FORM='formatted', ERR= 999)
      IF(bUseBufferedBinaryData) THEN
        DO L=1, nRating
          IF(pRating(L)%GroupID==PolicyGrp .and.
     &          pRating(L)%Policy==Policy .and.
     &          pRating(L)%CompType==CompType .and.
     &          pRating(L)%ID==IDNode)THEN
            n = n + 1
            NODE_ELEV(n,NODE)=pRating(L)%ElevOrWidth
            NODE_AREA(n,NODE)=pRating(L)%AreaOrEvapor   !aVar,aVar, Evgenii took out 2 avars
            NODE_VOL(n,NODE)=pRating(L)%VolOrFlow
            NODE_SEEP(n,NODE)=pRating(L)%Seep 
            Node_MaxQ(n,NODE)=pRating(L)%MaxQ
            NODE_LAKEQ(n,NODE)=pRating(L)%LakeQ         !aVar,!Evgenii added Node_MaxQ 090804, max release (for reservoirs)
           
            if (NODE_SEEP(n,NODE)>0)  seep_node(node) = .true.  !Evgenii 090720 true if node has seepage (used for node output)
!	      if (Node_MaxQ(n,NODE)>0)  maxflow(node) = .true.    !Evgenii 090720 true if node has maxflow
		    call UnitConversion(1, ULen, NODE_ELEV(n,NODE) )
            call UnitConversion(1, UArea, NODE_AREA(n,NODE) )
            call UnitConversion(1, UVol, NODE_VOL(n,NODE) )
            call UnitConversion(1, UFlow, NODE_SEEP(n,NODE) )
            call UnitConversion(1, UFlow, Node_MaxQ(n,NODE) )   !Evgenii added Node_MaxQ conversion 090804,This is converted into mil m3/day 
!           Evgenii 100603 put line below in to convert CapL() to mil m3/time step from mil m3/day (needed for its later use)
		    !Node_MaxQ(n,NODE)=Node_MaxQ(n,NODE)*sysstat(NPER)	        
		    call UnitConversion(1, UFlow, NODE_LAKEQ(n,NODE) )
            if (n >= IAGMAX ) exit
          END if
       END DO
      ELSE
      IF(nBufferedLines==0) REWIND(UNIT=iDatafile)
      
      !find line 'Rating:...'
      !Kang modify for improving performance
      L = 1
      IF(Rating_BL>0) L = Rating_BL      
      do WHILE (.TRUE.)
        IF(nBufferedLines>0) THEN
            IF(L>nBufferedLines) EXIT
            aLine = pFileData(L)
            L = L + 1
        ELSE    
        READ(UNIT=iDatafile,FMT='(A)',ERR=999, END=999) aLine
        END IF
       
        if (Index(aLine, headRating) >= 1) then
          IF(Rating_BL<0) Rating_BL = L -1      
          iPos = INDEX(aLine, ':')
          read(aLine(iPos+1:),*)GroupRead, PolicyRead,
     &         CompTypeRead, IDNodeRead
          IF(GroupRead==PolicyGrp.and.PolicyRead==Policy.and.
     &      CompTypeRead==CompType.and.IDNodeRead==IDNode)then
            n = n + 1
		  read(aLine(iPos+1:),*)aVar, aVar, aVar, aVar,
     &        NODE_ELEV(n,NODE), NODE_AREA(n,NODE),             !aVar,aVar, Evgenii took out 2 avars
     &        NODE_VOL(n,NODE),NODE_SEEP(n,NODE), 
     &        Node_MaxQ(n,NODE),NODE_LAKEQ(n,NODE)              !aVar,!Evgenii added Node_MaxQ 090804, max release (for reservoirs)
           
            if (NODE_SEEP(n,NODE)>0)  seep_node(node) = .true.  !Evgenii 090720 true if node has seepage (used for node output)
!	      if (Node_MaxQ(n,NODE)>0)  maxflow(node) = .true.    !Evgenii 090720 true if node has maxflow
		  call UnitConversion(1, ULen, NODE_ELEV(n,NODE) )
            call UnitConversion(1, UArea, NODE_AREA(n,NODE) )
            call UnitConversion(1, UVol, NODE_VOL(n,NODE) )
            call UnitConversion(1, UFlow, NODE_SEEP(n,NODE) )
            call UnitConversion(1, UFlow, Node_MaxQ(n,NODE) )   !Evgenii added Node_MaxQ conversion 090804,This is converted into mil m3/day 
!       Evgenii 100603 put line below in to convert CapL() to mil m3/time step from mil m3/day (needed for its later use)
		  !Node_MaxQ(n,NODE)=Node_MaxQ(n,NODE)*sysstat(NPER)	        
		  call UnitConversion(1, UFlow, NODE_LAKEQ(n,NODE) )
            if (n >= IAGMAX ) exit
          END if
        else
          IF(n > 0) exit
        END if
      end do
      END IF
      
      success = .true.
      NARVO_PTS(Node) = n
      LAKEQ_PTS(Node) = n
	
999   continue !Kang modify 20100630    CLOSE(UNIT=iDataFile)
      return
      end


!*************************************************************************
      subroutine readNodeEvaporation(iDatafile,GroupID,PolicyID,
     &     	     IDNode, iReadSeq, success)
!  read gage contribution multipliers for a node
!  CALL: None
      USE vars
	implicit none
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      CHARACTER*80  filename
      INTEGER*4     GroupID, PolicyID, iReadSeq, IDNode
!  OUTPUT
      logical*1     success
      !COMMON:      NODE_EVAP(iReadSeq)
!  Local variables
      INTEGER*4     i, j, iPos, iDatafile, comp_type
      INTEGER*4     GroupIDRead, PolicyIDRead, IDNodeRead
      CHARACTER*256 aLine
      CHARACTER*20  aVar
      INTEGER*4 L

!------------------------------------------------------------------------
      success = .false.
      !Kang add 20100629
      cNodeEvap = cNodeEvap +1
      !Kang modify 20100630 
!      iDatafile = 11
!      OPEN(UNIT=iDatafile, FILE=TRIM(filename), STATUS='old',
!     &	 FORM='formatted', ERR= 999)
      IF(bUseBufferedBinaryData) THEN      
          DO L=1, nEvapor
              IF(pEvapor(L)%GroupID==GroupID
     &          .and.pEvapor(L)%Policy==PolicyID
     &          .and.pEvapor(L)%CompType==1
     &          .and.pEvapor(L)%ID==IDNode)then
                NODE_EVAP(iReadSeq) = pEvapor(L)%Evaporation
                if(NODE_EVAP(iReadSeq)>0.) NODE_EVAP_ON(iReadSeq)=.true.
                call UnitConversion(1, ULoss, NODE_EVAP(iReadSeq) )    ! -> m/day
                exit
              END IF
          END DO
      ELSE
      
      IF(nBufferedLines==0) REWIND(UNIT=iDatafile)
      
      !find Evaporation: line
      !Kang modify for improving performance
      L = 1
      IF(Evaporation_BL>0) L = Evaporation_BL      
      do WHILE (.TRUE.)
        IF(nBufferedLines>0) THEN
            IF(L>nBufferedLines) EXIT
            aLine = pFileData(L)
            L = L + 1
        ELSE    
        READ(UNIT=iDatafile,FMT='(A)',ERR=999, END=999) aLine
        END IF
        
        if (Index(aLine,headEvaporation)==1) then
          IF(Evaporation_BL<0) Evaporation_BL = L - 1              
          iPos = INDEX(aLine, ':')
          read(aLine(iPos+1:),*)GroupIDRead,PolicyIDRead,Comp_type,
     &         IDNodeRead
          IF(GroupIDRead == GroupID.and.PolicyIDRead==PolicyID
     &        .and. Comp_type==1 .and.IDNodeRead==IDNode) then
            read(aLine(iPos+1:),*)aVar, aVar, aVar,aVar,
     &                      NODE_EVAP(iReadSeq)
            call UnitConversion(1, ULoss, NODE_EVAP(iReadSeq) )    ! -> m/day
            exit
          END if
        END if
      end do
      END IF
      success = .true.
999   continue !Kang modify 20100630    CLOSE(UNIT=iDataFile)
      return
      end


!******************************************************************
      Subroutine  readAllocation ( iDatafile, PolicyGrp, Policy,
     &    IDNode, Node, success)
! For a node:
! Read data of allocation and consumption to common variables.
! It calls none.
! 100211 Evgenii rewrote code so that multiple diversions are supported and also consumption on diversion nodes.
!	   This code does not have any checks to make sure that each diversion/consumption table has the same x values 
!	   (allocatable outflow of node). For now this is up to the user.
      USE vars
	implicit none
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      character*80 FileName
      INTEGER*4 PolicyGrp, Policy, Node, IDNode
      !COMMON: TOTOUT(), OUTLNK(i, j)
!  OUTPUT
      LOGICAL*1 success
      !COMMON
      ! ALLOC_PTS(NODE), NODE_ALLO(JJ,NODE), NODE_CONSUMP(I,NN)
      ! For node's outlinks:
      ! LINK_VTOL(JJ,NODE), jj=1,ALLOC_PTS(NODE)
! local
      INTEGER*4 i, n ,m,AlinkIDOld, GroupIDRead, PolicyRead
      INTEGER*4 NodeIDRead, iDatafile, iPos,t_cons_node
      CHARACTER*30 aVar
      CHARACTER*256 aLine
      INTEGER*4 aLinkID
      REAL*4 aLinkAllo, ANodeOutput, ALastNodeOutput
      INTEGER*4 L
!------------------------------------------------------------------------
      success = .false.

      !Kang add 20100629
      cAllocation = cAllocation + 1
      !Kang modify 20100630 
!      iDatafile = 11
!      OPEN(UNIT=iDatafile, FILE=TRIM(filename), STATUS='old',
!     &	 FORM='formatted', ERR= 999)
      n = 0
	m=0
      IF(bUseBufferedBinaryData) THEN
        DO L=1, nAllocation
          IF(pAllocation(L)%GroupID==PolicyGrp .and.
     &          pAllocation(L)%Policy==Policy .and.
     &          pAllocation(L)%NodeID==IDNode)then
              aNodeOutput = pAllocation(L)%NodeOutput
              aLinkID = pAllocation(L)%LinkID
              aLinkAllo = pAllocation(L)%LinkAllo
		      if (aLinkID/=0) then
			    if (aLinkID/=ALinkIDOld) n=0	
			    do i = 1, TOTOUT(NODE)
			      if (aLinkID == LinkID(OUTLNK(NODE,i))) then
			        n=n+1
				    NODE_ALLO(n,node)=aNodeOutput
				    call UnitConversion(1,UFlow, NODE_ALLO(n,node))
				    LINK_VTOL(n, OUTLNK(NODE,i)) = aLinkAllo
			        call UnitConversion(1,UFlow,
     &				    LINK_VTOL(n,OUTLNK(NODE,i)))
				    AlinkIDOld=alinkID
			      end if
                  end do

		      else if (aLinkID==0) then
                m=m+1
			    NODE_ALLO(m,node)=aNodeOutput
			    call UnitConversion(1,UFlow, NODE_ALLO(m,node))
			    NODE_CONSUMP(m,NODE) = aLinkAllo
                call UnitConversion(1,UFlow, NODE_CONSUMP(m,NODE))
			    cons_node(NODE)=.true.
		      END if
          END IF
        END DO
      ELSE
      IF(nBufferedLines==0) REWIND(UNIT=iDatafile)
      
	AlinkIDOld=0
      !find line 'Allocation:...'
      !Kang modify for improving performance
      L = 1
      IF(Allocation_BL>0) L = Allocation_BL
      do WHILE (.TRUE.)
        IF(nBufferedLines>0) THEN
            IF(L>nBufferedLines) EXIT
            aLine = pFileData(L)
            L = L + 1
        ELSE    
        READ(UNIT=iDatafile,FMT='(A)',ERR=999, END=888) aLine
        END IF

	  if (Index(aLine, headAllocation) >= 1) then
	    IF(Allocation_BL<0) Allocation_BL = L - 1
          iPos = INDEX(aLine, ':')
          read(aLine(iPos+1:),*)GroupIDRead, PolicyRead, aVar,NodeIDRead
          IF(GroupIDRead == PolicyGrp.and.PolicyRead == Policy
     &     	     .and. NodeIDRead==IDNode) then
            read(aLine(iPos+1:),*)aVar, aVar, aVar, aVar,
     &        aNodeOutput,aLinkID,aLinkAllo
            
		  if (aLinkID/=0) then
			if (aLinkID/=ALinkIDOld) n=0	
			do i = 1, TOTOUT(NODE)
			  if (aLinkID == LinkID(OUTLNK(NODE,i))) then
			    n=n+1
				NODE_ALLO(n,node)=aNodeOutput
				call UnitConversion(1,UFlow, NODE_ALLO(n,node))
				LINK_VTOL(n, OUTLNK(NODE,i)) = aLinkAllo
			    call UnitConversion(1,UFlow,
     &				 LINK_VTOL(n,OUTLNK(NODE,i)))
				AlinkIDOld=alinkID
			  end if
              end do

		  else if (aLinkID==0) then
              m=m+1
			NODE_ALLO(m,node)=aNodeOutput
			call UnitConversion(1,UFlow, NODE_ALLO(m,node))
			NODE_CONSUMP(m,NODE) = aLinkAllo
              call UnitConversion(1,UFlow, NODE_CONSUMP(m,NODE))
			cons_node(NODE)=.true.
		  END if
		end if
	  end if
	  !Kang remover the following statement 
	  !because it makes later "Allocation" not be read	  
	  !if (Index(aLine, 'END OF ALLOCATION') >= 1) exit
	 
	enddo
	END IF
	
888   ALLOC_PTS(NODE) = n
	if (cons_node(NODE)) ALLOC_PTS(NODE) = m
      success = .true.
999   continue !Kang modify 20100630    CLOSE(UNIT=iDataFile)
      return
      end


!******************************************************************
      Subroutine  ReadTarget ( iDatafile, PolicyGrp, Policy,
     &            IDNode, Node, success)
! For a node:
! read data of demand target and carry-over
!  CALL: None
      USE vars
	implicit none
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      character*80 FileName
      INTEGER*4 PolicyGrp, Policy, IDNode, Node
      ! PolicyGrp--related to simulation year
      ! Node--changes from 1 to TNodes determined by node reading sequence
!  OUTPUT
      LOGICAL*1 success
      !COMMON:
      !DMD_TARG(NN), DMD_T_CO(NN)
!  local
      INTEGER*4 GroupRead, IDNodeRead, iDatafile, iPos
      INTEGER*4 CompType, CompTypeRead, PolicyRead,EnvFlow
      CHARACTER*30 aVar
      CHARACTER*256 aLine
      INTEGER*4 L
!------------------------------------------------------------------------
      success = .false.
      CompType = 1  !1-for node; 2-for link
	EnvFlow=0 !Evgenii added 110513
	!Kang add 20100629 
	cTarget = cTarget +1
      !Kang modify 20100630 
!      iDatafile = 11
!      OPEN(UNIT=iDatafile, FILE=TRIM(filename), STATUS='old',
!     &	 FORM='formatted', ERR= 999)
      IF(bUseBufferedBinaryData) THEN
      if (nodeid(node)==20752.and. policy ==2)then
       continue
       end if
        DO L=1, nTarget
          IF(pTarget(L)%GroupID==PolicyGrp.and.pTarget(L)%Policy==Policy
     &          .and.pTarget(L)%CompType==CompType
     &          .and.pTarget(L)%NodeID==IDNode)then
            DMD_TARG(Node) = pTarget(L)%targ
            DMD_T_CO(Node) = pTarget(L)%t_co
            if (pTarget(L)%srcpriority==0)then
                srcpriorities(node)=.false.
            else if (pTarget(L)%srcpriority==1)  then
                srcpriorities(node)=.true. 
            end if
		    !Evgenii 100708 added refill trigger for non rule site grp resvrs, at which  
		    !storage level (fraction of capacity) of rule site reservoir can this reservior start refilling.
		    if (RESVNODE(node).and.RuleSite(node)/=nodeid(node)) then 
                RefilTrig(Node) = pTarget(L)%RefilTrig
		    end if                 
		  DemInc(node)=pTarget(L)%DemInc
		  if (capn(node)==0.0 .and.DemInc(node)/=0.0 
     &		  .and.sysstat(SIM_YEAR)>1)then
			  if (DemInc(node)>0) then
				DMD_TARG(Node)=DMD_TARG(Node)*(1+DemInc(node))
     &			  **(sysstat(SIM_YEAR)-1)
				else if (DemInc(node)<0)then
					DMD_TARG(Node)=DMD_TARG(Node)/((1+DemInc(node))
     &			  **(sysstat(SIM_YEAR)-1))
				end if

		  end if
			              
            !Evgenii- Aquifer could also be made a possible target node
            if(NatLak(node).or.resvnode(node))then          !Evgenii added following lines so volume targets could be used 090805
     			    call UnitConversion(1,Uvol,DMD_TARG(Node) )
		    else 
			    call UnitConversion(2,UFlow,DMD_TARG(Node) )  
			    if (pTarget(L)%EnvFlow>0)EnvFlwNode(node)=.true.!Evgenii added 110513
		    end if

            exit
          END IF
        END DO
        
      ELSE
        
      IF(nBufferedLines==0) REWIND(UNIT=iDatafile)
      
      !find line 'Target:...'
      !Kang modify for improving performance
      L = 1
      IF(Target_BL>0) L = Target_BL
      do WHILE (.TRUE.)
        IF(nBufferedLines>0) THEN
            IF(L>nBufferedLines) EXIT
            aLine = pFileData(L)
            L = L + 1
        ELSE    
        READ(UNIT=iDatafile,FMT='(A)',ERR=999, END=999) aLine
        END IF
        
        if (Index(aLine, headTarget)>= 1) then
          IF(Target_BL<0) Target_BL = L - 1
          iPos = INDEX(aLine, ':')
          read(aLine(iPos+1:),*)GroupRead, PolicyRead,
     &         CompTypeRead, IDNodeRead
          IF(GroupRead==PolicyGrp.and.PolicyRead==Policy.and.
     &           CompTypeRead==CompType.and.IDNodeRead==IDNode)then
            read(aLine(iPos+1:),*)aVar, aVar, aVar, aVar,
     &        DMD_TARG(Node), DMD_T_CO(Node),EnvFlow ,RefilTrig(node),
     &    			DemInc(Node)
    ! &		     !Evgenii EnvFlow added 110513

		  !Evgenii 100708 added refill trigger for non rule site grp resvrs, at which  
		  !storage level (fraction of capacity) of rule site reservoir can this reservior start refilling.
		  !if (RESVNODE(node).and.RuleSite(node)/=nodeid(node)) then 
		!	  read(aLine(iPos+1:),*)aVar,aVar,aVar,aVar,aVar,aVar, 
     & 	!		  aVar,RefilTrig(node)
		!  end if     
            !Evgenii- Aquifer could also be made a possible target node
		  if(NatLak(node).or.resvnode(node))then          !Evgenii added following lines so volume targets could be used 090805
     			call UnitConversion(1,Uvol,DMD_TARG(Node) )
		  else 
			call UnitConversion(2,UFlow,DMD_TARG(Node) )  
			if (EnvFlow>0)EnvFlwNode(node)=.true.!Evgenii added 110513
		  end if

            exit
          END if
        END if
      end do
      END IF
      
	success = .true.
999   continue !Kang modify 20100630    CLOSE(UNIT=iDataFile)
      return
      end

!************************************************************************
      Subroutine  ReadTargetSource (iDatafile, PolicyGrp, Policy,
     &            IDNode, Node, success)
! For a node:
! read data of demand target suply sources
!  CALL: None
      USE vars
	implicit none
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      character*80 FileName
      INTEGER*4 PolicyGrp, Policy, IDNode, Node
      ! PolicyGrp--related to simulation year
      ! Node--changes from 1 to TNodes determined by node reading sequence
!  OUTPUT
      LOGICAL*1 success
      !COMMON:
      !SUPLY_PTS(NN), SUPL_NODE(COUNT,NN), SUPL_FRAC(COUNT,NN)
!  local
      INTEGER*4 i, j, n, m, GroupRead, IDNodeRead, iDatafile, iPos
      INTEGER*4 CompType, CompTypeRead, PolicyRead,iprop
      CHARACTER*30 aVar
      CHARACTER*256 aLine
      INTEGER*4 L
!------------------------------------------------------------------------
      success = .false.
      CompType = 1  !1-for node; 2-for link
      SUPLY_PTS(Node) = 0
      SUPLY_Links_pts(node) =0
      !Kang add 20100629
      cTargetSource = cTargetSource + 1
      !Kang modify 20100630 
!      iDatafile = 11
!      OPEN(UNIT=iDatafile, FILE=TRIM(filename), STATUS='old',
!     &	 FORM='formatted', ERR= 999)

      n = 0
      m = 0
      IF(bUseBufferedBinaryData) THEN
        DO L=1, nSource
          IF(pSource(L)%GroupID==PolicyGrp.and.pSource(L)%Policy==Policy
     &          .and.pSource(L)%CompType==CompType
     &          .and.pSource(L)%NodeID==IDNode)then
                    
            if(pSource(L)%Source_Type==1) then
                n = n + 1  
                !Source_Type(n,Node)=pSource(L)%Source_type
                SUPL_NODE(n,Node) = pSource(L)%Supl_Node
		        !Source_Type(n,Node) = pSource(L)%Source_Type
                SUPL_FRAC(n,Node) = pSource(L)%Supl_Frac
		        MaxOutTS(n,Node)=  pSource(L)%MaxOutTS	 !110715 Evgenii added MaxOutTS
		        MaxOutYear(n,Node)=  pSource(L)%MaxOutYear	 !280915 Evgenii added MaxOutYear
                !SrcPriority(n,Node)=pSource(L)%SrcPriority
		        !SrcPriority(n,node)= pSource(L)%SrcPriority
		        if(MaxOutTS(n,node)>=0.0) then
     		            call UnitConversion(1,UFlow,MaxOutTS(n,Node)) !Converts to Mm3/day
     		            MaxOutTS(n,node)=MaxOutTS(n,node)*DAYSPRPRD	!Converts to Mm3/TS
     		        end if	    
     		        if(MaxOutYear(n,node)>=0.0) then
     		            call UnitConversion(1,UFlow,MaxOutYear(n,Node)) !Converts to Mm3/day 
     		            MaxOutYear(n,node)=MaxOutYear(n,node)*372.0	!Converts to Mm3/year  !Converts to Mm3/year !Right now annual licenses are based on 371 days because of a leap year bug and also that when running with a weekly time step, so years are 3*7 time steps long, this needs to be corrected so only leap years are based on 366 days.
     		        end if
		    else if(pSource(L)%Source_Type==2) then
		        m = m + 1		        
		        !Source_Type(n,Node)=pSource(L)%Source_type
		        SUPL_Link(m,node)= pSource(L)%Supl_Node !This is the link source ID for node
		        SUPL_FRAC_Link(m,node) = pSource(L)%Supl_Frac
		        LNKPROP(m,node)=.false. !Evgenii added LNKPROP120126 for link propogation links
		        if (pSource(L)%LnkProp>0)then
		            LNKPROP(m,node)=.true.
		        end if	
		        !MaxOut(m,node)=  pSource(L)%MaxOut			        		        		   	    
		    end if
		  END IF
        END DO
      ELSE
      IF(nBufferedLines==0) REWIND(UNIT=iDatafile)
      
      !find line 'Source:...'
      !Kang modify for improving performance
      L = 1
      IF(Source_BL>0) L = Source_BL
      do WHILE (.TRUE.)
        IF(nBufferedLines>0) THEN
            IF(L>nBufferedLines) EXIT
            aLine = pFileData(L)
            L = L + 1
        ELSE    
        READ(UNIT=iDatafile,FMT='(A)',ERR=999, END=999) aLine
        END IF

        if (Index(aLine, headSource) >= 1) then
          IF(Source_BL<0) Source_BL = L - 1
          iPos = INDEX(aLine, ':')
          read(aLine(iPos+1:),*)GroupRead, PolicyRead,
     &         CompTypeRead, IDNodeRead
          IF(GroupRead==PolicyGrp.and.PolicyRead==Policy.and.
     &      CompTypeRead==CompType.and.IDNodeRead==IDNode)then
            n = n + 1
            !Evgenii 110404 added source type to distinguish between node and link sources
		  read(aLine(iPos+1:),*)aVar, aVar, aVar, aVar,
     &        SUPL_NODE(n,Node), Source_Type(n,Node), SUPL_FRAC(n,Node),
     &		MaxOutTS(n,Node),MaxOutYear(n,Node),iprop !, SrcPriority(n,Node)														!Evgenii added Maxout 110714
			call UnitConversion(1,UFlow,MaxOutTS(n,Node)) 		 !110715 Evgenii added MaxOut
			call UnitConversion(1,UFlow,MaxOutYear(n,Node)) 	
		END if
        else
          if (n > 0) exit
        END if
      end do
      END IF
      !

      success = .true.
999   SUPLY_PTS(Node) = n !Number of supply reserviors for node
      SUPLY_Links_pts(node) = m !Number of supply links for node
      
      !Kang modify 20100630 CLOSE(UNIT=iDataFile)
            
      !Set SUPL_NODE() from NodeID to NodeSeq
         
      if (n > 0) then
        do i = 1, n
          do j = 1, TNodes
             if (NodeID(j) == SUPL_NODE(i,Node)) then
               SUPL_NODE(i,Node) = j
               exit
             end if
!            if (NodeID(j) == SUPL_NODE(i,Node)!.and.
!     &		   Source_Type(i,Node)==1) then  !Evgenii added Source_Type(n,Node)=1 so links can be sources too and they dont change their ID here
!			   SUPL_NODE(i,Node) = j
!              exit
!            end if
          end do
        end do
      end if
   
      if (m > 0) then !if (m > 0) then node has one or more supply links     
        do i = 1, m
          do j = 1, Links
             if (LinkID(j) == SUPL_Link(i,node)) then !Change from Link ID in INP (user ID) to internal IRAS ID
                SUPL_Link(i,Node)=j
                !ndemnodelink(j)=ndemnodelink(j)+1 !Increase number of nodes link supplies
               exit
             end if
          end do
        end do
      end if     
      
      return
      end

!******************************************************************
      Subroutine  ReadPowerPump (iDatafile,PolicyGrp, Policy,
     &            IDNode, PowerLink1, iPower, success)
C  for Hydropower/pump
!  CALL: None
      USE vars
	implicit none
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      character*80 FileName
      INTEGER*4 PolicyGrp, Policy, IDNode, PowerLink1, iPower
!  OUTPUT
      LOGICAL*1 success
      !COMMON:
      ! HPCAP(LNKMAX), ECONST(LNKMAX), PLANT_FACTOR(LNKMAX)
      ! HPQMIN(LNKMAX), TURBINE_ELEV(LNKMAX), PCONST(LNKMAX)
!  Local
      INTEGER*4 i, j, iDatafile, iPos, iLinkSeq
      INTEGER*4 PolicyRead, GroupRead, IDNodeRead
      CHARACTER*30 aVar, sFind
      CHARACTER*256 aLine
      INTEGER*4 L
!-------------------------------------------------------------------------
      success = .false.
      !find no. in link sequence
      do i = 1, Links
        IF(LinkID(i) == PowerLink1) then
          iLinkSeq = i;         exit
        END if
      end do
      !Kang add 20100629
      cPowerPump = cPowerPump +1
      !Kang modify 20100630 
!      iDatafile = 11
!      OPEN(UNIT=iDatafile, FILE=TRIM(filename), STATUS='old',
!     &	 FORM='formatted', ERR= 999)
      IF(bUseBufferedBinaryData) THEN
        IF(iPower==1)then
            DO L=1, nPower
              IF(pPower(L)%GroupID==PolicyGrp.and.
     &          pPower(L)%Policy==Policy.and.
     &          pPower(L)%NodeID==IDNode)then
                  PowerLink(iLinkSeq) = .true.
                  HPCAP(iLinkSeq) = pPower(L)%HPCAP
                  PLANT_FACTOR(iLinkSeq) = pPower(L)%PLANT_FACTOR
                  Intake_elev(iLinkSeq) = pPower(L)%Intake_elev
                  Turbine_elev(iLinkSeq) = pPower(L)%Turbine_elev
                  Outlet_elev(iLinkSeq) = pPower(L)%Outlet_elev
                  ECONST(iLinkSeq) = pPower(L)%ECONST
                  HpQMin(iLinkSeq) = pPower(L)%HpQMin					!, HpQMax(iLinkSeq) Evenii took out max because not used
                  exit
              END IF    
            END DO
        ELSE    
            DO L=1, nPump
              IF(pPump(L)%GroupID==PolicyGrp.and.
     &          pPump(L)%Policy==Policy.and.
     &          pPump(L)%NodeID==IDNode)then
                  PumpLink(iLinkSeq) = .true.
                  IF(PowerLink(iLinkSeq)) then
                      PCONST(iLinkSeq) = pPump(L)%ECONST
                  ELSE
                      HPCAP(iLinkSeq) = pPump(L)%HPCAP
                      PLANT_FACTOR(iLinkSeq) = pPump(L)%PLANT_FACTOR
                      Intake_elev(iLinkSeq) = pPump(L)%Intake_elev
                      Turbine_elev(iLinkSeq) = pPump(L)%Turbine_elev
                      Outlet_elev(iLinkSeq) = pPump(L)%Outlet_elev
                      PCONST(iLinkSeq) = pPump(L)%ECONST
                      HpQMin(iLinkSeq) = pPump(L)%HpQMin					!, HpQMax(iLinkSeq) Evenii took out max because not used
                  END IF
                  exit
              END if
            END DO
        END IF
      ELSE
      IF(nBufferedLines==0) REWIND(UNIT=iDatafile)
      
      !Kang modify for improving performance
      L = 1
      
      if (iPower ==1) then
        sFind = headPower
        IF(Power_BL>0) L = Power_BL
      else
        sFind = headPump
        IF(Pump_BL>0) L = Pump_BL
      end if
      !find line 'Power:' or 'Pump:'
      !Kang modify for improving performance
      do WHILE (.TRUE.)
        IF(nBufferedLines>0) THEN
            IF(L>nBufferedLines) EXIT
            aLine = pFileData(L)
            L = L + 1
        ELSE    
        READ(UNIT=iDatafile,FMT='(A)',ERR=999, END=999) aLine
        END IF
        
        if (Index(aLine, TRIM(sFind)) >= 1) then
          IF (iPower ==1) THEN
            IF(Power_BL<0) Power_BL = L - 1
          ELSE
            IF(Pump_BL<0) Pump_BL = L - 1
          END IF
          iPos = INDEX(aLine, ':')
          read(aLine(iPos+1:),*)GroupRead, PolicyRead, aVar,IDNodeRead
          IF(GroupRead==PolicyGrp.and.PolicyRead==Policy
     &                           .and.IDNodeRead==IDNode)then
            IF(iPower==1)then
              PowerLink(iLinkSeq) = .true.
              read(aLine(iPos+1:),*)aVar, aVar, aVar,aVar,
     &           HPCAP(iLinkSeq),PLANT_FACTOR(iLinkSeq),
     &           Intake_elev(iLinkSeq),Turbine_elev(iLinkSeq),
     &           Outlet_elev(iLinkSeq),ECONST(iLinkSeq),
     &           HpQMin(iLinkSeq)					!, HpQMax(iLinkSeq) Evenii took out max because not used
            else
              PumpLink(iLinkSeq) = .true.
              IF(PowerLink(iLinkSeq)) then
                read(aLine(iPos+1:),*)aVar, aVar, aVar,aVar,
     &              aVar,aVar, aVar, aVar,aVar,PCONST(iLinkSeq),
     &              aVar							!, aVar evgenii took out spot for max
              else
                read(aLine(iPos+1:),*)aVar, aVar, aVar,aVar,
     &              HPCAP(iLinkSeq),PLANT_FACTOR(iLinkSeq),
     &              Intake_elev(iLinkSeq),Turbine_elev(iLinkSeq),
     &              Outlet_elev(iLinkSeq),PCONST(iLinkSeq),
     &              HpQMin(iLinkSeq)				!, HpQMax(iLinkSeq) Evenii took out max because not used
              END if
            END if
            exit
          END if
        END if
      end do
      END IF
      
      success = .true.
      call UnitConversion(1,UPower, HPCAP(iLinkSeq) )
      call UnitConversion(1,ULen, Intake_elev(iLinkSeq) )
      call UnitConversion(1,ULen, Turbine_elev(iLinkSeq) )
      call UnitConversion(1,ULen, Outlet_elev(iLinkSeq) )
      call UnitConversion(1,UFlow, HpQMin(iLinkSeq) )
      call UnitConversion(1,UFlow, HpQMax(iLinkSeq) )
999   continue !Kang modify 20100630    CLOSE(UNIT=iDataFile)
      return
      end


!************************************************************************
      Subroutine ReadReleaseRule(iDatafile, PolicyGrp, Policy,
     &           ResvGrp, IDNode, success)
!  Read group release rule for a group of reservoirs
!  CALL: None
      USE vars
	implicit none
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  Input:
      CHARACTER*80 FileName
      INTEGER*4 PolicyGrp,Policy,ResvGrp,IDNode
!  Output:
      LOGICAL*1 success
      !COMMON: Rule_PTS(ResvGrp), res_rule(*,Resv,ResvGrp)
!  Local:
      INTEGER*4  GroupRead, NodeIDRead, PolicyRead
      INTEGER*4  i, j, n , iDatafile, iPos
      CHARACTER*30 aVar
      CHARACTER*256 aLine
      INTEGER*4 L
!------------------------------------------------------------------------
      success = .false.
      !Kang add 20100629
      cReleaseRule = cReleaseRule + 1
      !Kang modify 20100630 
!      iDatafile = 11
!      OPEN(UNIT=iDatafile, FILE=TRIM(filename), STATUS='old',
!     &	 FORM='formatted', ERR= 999)
      n = 0
      IF(bUseBufferedBinaryData) THEN
        DO L=1, nRule
          IF(pRule(L)%GroupID==PolicyGrp.and.pRule(L)%Policy==Policy
     &          .and.pRule(L)%NodeID==IDNode)then
            n = n + 1
            DO I=1, 8 
                res_rule(I,n,ResvGrp) = pRule(L)%res_rule(I)
            END DO
    
		  !Unit conversion
            call unitConversion(1,UVol, res_rule(1,n,ResvGrp) )
            call unitConversion(1,UVol, res_rule(3,n,ResvGrp) )
            call unitConversion(1,UVol, res_rule(5,n,ResvGrp) )
            call unitConversion(1,UVol, res_rule(7,n,ResvGrp) )
            call unitConversion(1,UFlow, res_rule(2,n,ResvGrp) )
            call unitConversion(1,UFlow, res_rule(4,n,ResvGrp) )
            call unitConversion(1,UFlow, res_rule(6,n,ResvGrp) )
		  call unitConversion(1,UFlow, res_rule(8,n,ResvGrp) )

		  if (n >= ZNMAX ) exit
          END IF
        END DO
      ELSE
      IF(nBufferedLines==0) REWIND(UNIT=iDatafile)
      
     
	!find Rule line
      !Kang modify for improving performance
      L = 1
      IF(Rule_BL>0) L = Rule_BL
      do WHILE (.TRUE.)
        IF(nBufferedLines>0) THEN
            IF(L>nBufferedLines) EXIT
            aLine = pFileData(L)
            L = L + 1
        ELSE    
        READ(UNIT=iDatafile,FMT='(A)',ERR=999, END=999) aLine
        END IF

        if (Index(aLine, headRule) >= 1) then
          IF(Rule_BL<0) Rule_BL = L - 1
          iPos = INDEX(aLine, ':')
          read(aLine(iPos+1:),*)GroupRead, PolicyRead, NodeIDRead
	    !write(*,*)GroupRead, PolicyRead, NodeIDRead !evgenii
          IF(GroupRead == PolicyGrp.and.PolicyRead==Policy
     &          .and. NodeIDRead==IDNode) then
		  n = n + 1
            read(aLine(iPos+1:),*)aVar, aVar, aVar, 
     &        res_rule(1,n,ResvGrp),res_rule(2,n,ResvGrp),
     &        res_rule(3,n,ResvGrp),res_rule(4,n,ResvGrp),
     &        res_rule(5,n,ResvGrp),res_rule(6,n,ResvGrp),
     &        res_rule(7,n,ResvGrp),res_rule(8,n,ResvGrp)
    
		  !Unit conversion
            call unitConversion(1,UVol, res_rule(1,n,ResvGrp) )
            call unitConversion(1,UVol, res_rule(3,n,ResvGrp) )
            call unitConversion(1,UVol, res_rule(5,n,ResvGrp) )
            call unitConversion(1,UVol, res_rule(7,n,ResvGrp) )
            call unitConversion(1,UFlow, res_rule(2,n,ResvGrp) )
            call unitConversion(1,UFlow, res_rule(4,n,ResvGrp) )
            call unitConversion(1,UFlow, res_rule(6,n,ResvGrp) )
		  call unitConversion(1,UFlow, res_rule(8,n,ResvGrp) )

		  if (n >= ZNMAX ) exit
          END if
	  else
          IF(n > 0) exit
        END if
      end do
      END IF
	
      success = .true.
999   Rule_PTS(ResvGrp) = n
      !Kang modify 20100630    CLOSE(UNIT=iDataFile)
      return
      end

!******************************************************************
      Subroutine  readGrpBalance (iDatafile, PolicyGrp, Policy,
     &            iGrpSeq, RuleGrp, success)
! Read balancing function for a group of reservoirs
! It calls none.
      USE vars
	implicit none
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      character*80 FileName
      INTEGER*4 PolicyGrp, Policy, RuleGrp, iGrpSeq
!  OUTPUT
      LOGICAL*1 success
      !COMMON: GrpVOL_PTS(RNMAX),ResvIDInGrp(RNMAX),
      !        GrpVol(IAGMAX,RNMAX),BalanceVol(points, nResv, nGrp)
!  local
      INTEGER*4 i, n , GroupIDRead, PolicyRead, bal
      INTEGER*4 RuleGrpRead, iDatafile, iPos
      CHARACTER*30 aVar
      CHARACTER*1000 aLine
      INTEGER*4 L
!------------------------------------------------------------------------
      success = .false.
	n=0
      !Kang add 20100629
      cGrpBalance = cGrpBalance +1
      !Kang modify 20100630 
!      iDatafile = 11
!      OPEN(UNIT=iDatafile, FILE=TRIM(filename), STATUS='old',
!     &	 FORM='formatted', ERR= 999)
      IF(bUseBufferedBinaryData) THEN
        DO L=1, nBalance
          IF(pBalance(L)%GroupID==PolicyGrp.and.pBalance(L)%Policy==
     &      Policy.and.pBalance(L)%RuleID==RuleGrp)then
            n = n + 1
            BalMethod(n,IGrpSeq)=BalanceArray(nBalance)%BalMeth
		    GrpVol(n,IGrpSeq)=BalanceArray(nBalance)%GroupVol	                                							
            !BalCols=CountColumns(pBalance(L)%charBalance)
		    aline=pBalance(L)%charBalance
		    iPos = INDEX(aLine, ':')
		    read(aLine(iPos+1:),*)aVar, aVar, aVar,avar,avar,          
     &         ((ResvIDInGrp(i,iGrpSeq),UseResVol(i,iGrpSeq),
     &		   RuleResVol(n,i,iGrpSeq),BalanceVol(n,i,iGrpSeq)),
     &		  i=1,nResvInGrp(iGrpSeq))  ! !  
 		    call unitConversion(1,UVol, GrpVol(n,iGrpSeq) )
            do i = 1,nResvInGrp(iGrpSeq)
              call unitConversion(1,UVol, BalanceVol(n,i,iGrpSeq) )
              call unitConversion(1,UVol, RuleResVol(n,i,iGrpSeq) ) !Evgenii 121120 added a converstion for RuleResVol(n,i,iGrpSeq)
            end do
		    if (n >= IAGMAX ) then
              exit
              end if
		  end if      
	   end do
	else
	continue
      IF(nBufferedLines==0) REWIND(UNIT=iDatafile)
      

      !find line 'Balance:...'
      !Kang modify for improving performance
      L = 1
      IF(Balance_BL>0) L = Balance_BL
      do WHILE (.TRUE.)
        IF(nBufferedLines>0) THEN
            IF(L>nBufferedLines) EXIT
            aLine = pFileData(L)
            L = L + 1
        ELSE    
		READ(UNIT=iDatafile,FMT='(A)',ERR=999, END=999) aLine
        END IF
       
	  if (Index(aLine, 'Balance:') >= 1) then
	    IF(Balance_BL<0) Balance_BL = L - 1
          iPos = INDEX(aLine, ':')
          read(aLine(iPos+1:),*)GroupIDRead,PolicyRead,RuleGrpRead
          IF(GroupIDRead == PolicyGrp.and.PolicyRead == Policy.and.
     &          RuleGrpRead==RuleGrp) then
           !do i=1,nResvInGrp(iGrpSeq)
		  n = n + 1
          !if (bal==0) then !later add in read bal so when balance is on lead reservior you dont need to put in its own balance
		read(aLine(iPos+1:),*)aVar, aVar, aVar,
     &         BalMethod(n,IGrpSeq),GrpVol(n,IGrpSeq), !Evgenii 1007285 added balance method choice; 0 - Old iras balance method (group storage), 1 - Balance depends only on rule reservoir storage
     &         ((ResvIDInGrp(i,iGrpSeq),RuleResVol(n,i,iGrpSeq),
     &	     BalanceVol(n,i,iGrpSeq)),i=1,nResvInGrp(iGrpSeq))                     		  		  
		  call unitConversion(1,UVol, GrpVol(n,iGrpSeq) )
            do i = 1,nResvInGrp(iGrpSeq)
              call unitConversion(1,UVol, BalanceVol(n,i,iGrpSeq) )
            end do
            if (n >= IAGMAX ) exit
          END if
        else
          if (n > 0) exit
        END if
      end do
	end if
      success = .true.
999   GrpVOL_PTS(iGrpSeq) = n
      !Kang modify 20100630    CLOSE(UNIT=iDataFile)
      return
      end

!******************************************************************
       Subroutine  readPerformance ( iDatafile, PolicyGrp, Policy,
     &    IDNode, Node, success)
! Created by Evgenii 100617 to read in performance thresholds
! For storage and demand nodes:
! read data for performance thresholds
!  CALL: None
      USE vars
	implicit none
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      character*80 FileName
      INTEGER*4 PolicyGrp, Policy, IDNode, Node
      ! PolicyGrp--related to simulation year
      ! Node--changes from 1 to TNodes determined by node reading sequence
!  OUTPUT
      LOGICAL*1 success

!  local
      INTEGER*4 i, n , GroupRead, IDNodeRead, iDatafile, iPos
      INTEGER*4 PolicyRead
      CHARACTER*30 aVar
      CHARACTER*256 aLine
      INTEGER*4 L
!------------------------------------------------------------------------
      success = .false.
      
      !Kang modify 20100630
!      iDatafile = 11
!      OPEN(UNIT=iDatafile, FILE=TRIM(Filename), STATUS='old',
!     &	 FORM='formatted', ERR= 999)
      n = 0
      IF(bUseBufferedBinaryData) THEN
        DO L=1, nPerformance
          IF(pPerformance(L)%GroupID==PolicyGrp .and.
     &          pPerformance(L)%Policy==Policy .and.
     &          pPerformance(L)%NodeID==IDNode)then
            n = n + 1
            thres_limit(n,node) = pPerformance(L)%thres_limit
            !if (capn(node)>0.0) then
            sto_perf_node(node) = .true.
           ! else
             !   flow_perf_node(node) = .true.
            !end if
            if (n >= thres_pts_max) exit
          END IF
        END DO
 !     ELSE
	!IF(nBufferedLines==0) REWIND(UNIT=iDatafile)
	!		
 !     !find line 'Performance:...'
 !     !Kang modify for improving performance
 !     L = 1
 !     IF(Performance_BL>0) L = Performance_BL
 !     do WHILE (.TRUE.)
 !       IF(nBufferedLines>0) THEN
 !           IF(L>nBufferedLines) EXIT
 !           aLine = pFileData(L)
 !           L = L + 1
 !       ELSE    
 !       READ(UNIT=iDatafile,FMT='(A)',ERR=999, END=999) aLine
 !       END IF
 !
 !       if (Index(aLine, headPerformance) >= 1) then
 !         IF(Performance_BL<0) Performance_BL = L - 1
 !         iPos = INDEX(aLine, ':')
 !         read(aLine(iPos+1:),*)GroupRead, PolicyRead,IDNodeRead          
 !         IF(GroupRead==PolicyGrp.and.PolicyRead==Policy.and.
 !    &      IDNodeRead==IDNode)then
 !           n = n + 1
	!	  read(aLine(iPos+1:),*)aVar, aVar, aVar, 
 !    &        thres_limit(n,node)             
 !			sto_perf_node(node) = .true.
 !           if (n >= thres_pts_max) exit
 !         END if
 !       else
 !         IF(n > 0) exit
 !       END if
 !     end do
      END IF
      
      success = .true.
      perf_pts(Node) = n

	
999   continue !Kang modify 20100630   CLOSE(UNIT=iDataFile)
      return
      end

     
!******************************************************************	
       Subroutine  readDemRed ( iDatafile, PolicyGrp, Policy,
     &    IDNode, Node, success)
! Created by Evgenii 100617 to read in performance thresholds
! For storage and demand nodes:
! read data for performance thresholds
!  CALL: None
      USE vars
	implicit none
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      character*80 FileName
      INTEGER*4 PolicyGrp, Policy, IDNode, Node
      ! PolicyGrp--related to simulation year
      ! Node--changes from 1 to TNodes determined by node reading sequence
!  OUTPUT
      LOGICAL*1 success

!  local
      INTEGER*4 i, n , GroupRead, IDNodeRead, iDatafile, iPos
      INTEGER*4 PolicyRead
      CHARACTER*30 aVar
      CHARACTER*256 aLine
      INTEGER*4 L
!------------------------------------------------------------------------
      success = .false.
      
      !Kang modify 20100630
!      iDatafile = 11
!      OPEN(UNIT=iDatafile, FILE=TRIM(Filename), STATUS='old',
!     &	 FORM='formatted', ERR= 999)
      n = 0
      IF(bUseBufferedBinaryData) THEN
        DO L=1, nDemRed
          IF(pDemRed(L)%GroupID==PolicyGrp .and.
     &          pDemRed(L)%Policy==Policy .and.
     &          pDemRed(L)%NodeID==IDNode)then
            n = n + 1
            DemSourceID(node) = pDemRed(L)%DemSourceID
		  Dem_Thres_limit(n,node)  = pDemRed(L)%Dem_Thres_limit
		  DemRedAmt(n,node)  = pDemRed(L)%DemRedAmt
            DemRed_node(node) = .true.
            if (n >= thres_pts_max) exit
		END IF
        END DO
      ELSE
	IF(nBufferedLines==0) REWIND(UNIT=iDatafile)
			
      !find line 'DemRed:...'
   
      L = 1
      IF(DemRed_BL>0) L = DemRed_BL
      do WHILE (.TRUE.)
        IF(nBufferedLines>0) THEN
            IF(L>nBufferedLines) EXIT
            aLine = pFileData(L)
            L = L + 1
        ELSE    
        READ(UNIT=iDatafile,FMT='(A)',ERR=999, END=999) aLine
        END IF

        if (Index(aLine, headDemRed) >= 1) then
          IF(DemRed_BL<0) DemRed_BL = L - 1
          iPos = INDEX(aLine, ':')
          read(aLine(iPos+1:),*)GroupRead, PolicyRead,IDNodeRead          
          IF(GroupRead==PolicyGrp.and.PolicyRead==Policy.and.
     &      IDNodeRead==IDNode)then
            n = n + 1
		  read(aLine(iPos+1:),*)aVar, aVar, aVar, 
     &        DemSourceID(node),Dem_Thres_limit(n,node),
     &   		DemRedAmt(n,node)  
			DemRed_node(node) = .true.
            if (n >= thres_pts_max) exit
          END if
        else
          IF(n > 0) exit
        END if
      end do
      END IF
      
      success = .true.
      DemRedPts(Node) = n

	
999   continue !Kang modify 20100630   CLOSE(UNIT=iDataFile)
      return
      end
!******************************************************************
	SUBROUTINE UnitConversion(iComp,iType, x)
      implicit none
      INCLUDE 'IRAS_SYS.INC'
!  INPUT
      INTEGER*4 iComp, iType         !iComp = 1 for nodes; 2 for links
      REAL*4 x
!  OUTPUT
      !x
	select case (iType)
          case (ULen,UArea,UVol,UFlow,ULoss,UPower,UK,UTime)
            if (iComp == 1) x = x * NodeUserUnit(iType) 
            if (iComp == 2) x = x * LinkUserUnit(iType)
          case default
            x = x     !no change
      end select
	
      END SUBROUTINE
      
      !******************************************************************
	SUBROUTINE UndoUnitConversion(iComp,iType, x)
      implicit none
      INCLUDE 'IRAS_SYS.INC'
!  INPUT
      INTEGER*4 iComp, iType         !iComp = 1 for nodes; 2 for links
      REAL*4 x
!  OUTPUT
      !x
	select case (iType)
          case (ULen,UArea,UVol,UFlow,ULoss,UPower,UK,UTime)
            if (iComp == 1) x = x / NodeUserUnit(iType)
            if (iComp == 2) x = x / LinkUserUnit(iType)
          case default
            x = x     !no change
      end select
	
      END SUBROUTINE

!      subroutine readGage(Filename,GroupID,PolicyID,
!     &     	     IDNode, iReadSeq,  success) !iType, Evgenii took out iType (no added guage anymore) 091026
!  read gage contribution multipliers for a node
!  CALL: None
!      implicit none
!      INCLUDE 'IRAS_SYS.INC'
!      INCLUDE 'NODE.INC'
!      INCLUDE 'LINK.INC'
!  INPUT
!      CHARACTER*80  filename
!      INTEGER*4     GroupID, PolicyID, iReadSeq, IDNode, iType
      !COMMON:      SysStat(P_NGauge)
!  OUTPUT
!      logical*1     success
      !COMMON:      GageID(*,iReadSeq), GageMultiplier(*,iReadSeq)
!  Local variables
!      INTEGER*4     i, j, iPos, iDatafile, nGagesSys
!      INTEGER*4     GroupIDRead, PolicyIDRead, IDNodeRead, nGagesRead
!      CHARACTER*256 aLine
!      CHARACTER*20  aVar, sName

!------------------------------------------------------------------------
!      success = .false.
      
	!if (iType ==1) then !Evgenii commnented this out because no more added gauge 091026
!        sName = 'NaturalFlowGage:'
!        nGagesSys = SysStat(nGauge1)
      !else !Evgenii commnented this out because no more added gauge 091026
      !  sName = 'AddedFlowGage:' 
      !  nGagesSys = SysStat(nGauge2)
!      	end if
!      iDatafile = 11
!      OPEN(UNIT=iDatafile, FILE=TRIM(filename), STATUS='old',
!     &	 FORM='formatted', ERR= 999)
!      nGagesRead = 0
      !find sName line
!	do WHILE (.TRUE.)
        
!	  READ(UNIT=iDatafile,FMT='(A)',ERR=999, END=999) aLine
!	  if (Index(aLine,TRIM(sName))>=1) then
!          iPos = INDEX(aLine, ':')
!          read(aLine(iPos+1:),*)GroupIDRead,PolicyIDRead,IDNodeRead
!		IF(GroupIDRead == GroupID.and.PolicyIDRead==PolicyID
!     &                     .and.IDNodeRead==IDNode) then
!            nGagesRead = nGagesRead + 1
!            read(aLine(iPos+1:),*)aVar, aVar, aVar,
!     &            GageID(nGagesRead,iType,iReadSeq), aVar,
!     &            GageMultiplier(nGagesRead,iType,iReadSeq)
     	
!		  if (nGagesRead == nGagesSys) exit	
!		END if
!        else
!          IF(nGagesRead > 0)exit
!        END if
!      end do
	
!      success = .true.
!999   CLOSE(UNIT=iDataFile)
!      return
!      end
!******************************************************************