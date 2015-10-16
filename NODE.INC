!Copyright (c) 2008, 2009 by University College London, Cornell University
!Authors:
!Daniel P. Loucks (dpl3@cornell.edu), Marshall Taylor, Peter French, 
!Evgenii Matrosov evgenii.matrosov@ucl.ac.uk), Julien Harou (j.harou@ucl.ac.uk), HC Zhou
!This program is free software under the General Public Licence, GPL (>=v2)
!Read the 'GPL License.txt' file distributed with this source code for a full license statement.
C+IRIS+IRIS+IRIS+IRIS+IRIS+IRIS+IRIS+IRIS+IRIS+IRIS+IRIS+IRIS+IRIS+IRIS
C
C      File: NODE.INC
C      (Prerequisites are FND_VARS.INC and SIMDIM.INC)
C
C+IRIS+IRIS+IRIS+IRIS+IRIS+IRIS+IRIS+IRIS+IRIS+IRIS+IRIS+IRIS+IRIS+IRIS
C  ************************* SYSTEM DATA *****************************
C     Total number of nodes defined in network, number of nodes
C     defined in simulation sequence, total number of links in netwrk,
C     node number associated with sequence index NX.
      INTEGER*4      TNODES, NODES, LINKS
      logical cons_node(NODMAX), seep_node(NODMAX),maxflow(nodmax) !Evgenii 090717
	  !Total gage and storage nodes Evgenii 090717
C    NODSEQ = vector of node sequencing ( use RD_FND(SYSGRP,NNSEQ, ) )
      INTEGER*4  NODSEQ(NODMAX),SEQ_NOD(NODMAX)
      INTEGER*4 FFFlag(NODMAX) ! added by FF August 2012
      INTEGER*4 FFSN(NODMAX) 
C
C
C  ************************* HYDROPOWER  *****************************
C     Total system wide hydropower production.
      REAL*4         SYSPWR
C
C  ************************* HYDROPOWER  *****************************

C     This file declares the parameters and variables which are used to
C     define a node.
C
C
C                        NODE DEFINITION
C     Total number of inflow links, total number of outflow links,
C     identifier for links flowing in to, and out of a node,
C     identifier for node at the inflow, and outflow end of a link.
C     Screen coordinates of nodes (X and Y) for schematic representation.

C     Node name
!     Evgenii 2009 21 4 changed character length to 95 format  
	  CHARACTER(len=100)    NNAME(NODMAX)
      
	  INTEGER*4 NodeID(NODMAX)         !Zhou added
C     Node Position
!     Evgenii 2009 21 took out *8 after integer and replaced with  ::       
	  INTEGER ::      IX(NODMAX), IY(NODMAX)

C     Number of inflow links and outflow links and the link numbers
      INTEGER*4      TOTIN(NODMAX), TOTOUT(NODMAX),
     1               INLINK(NODMAX,ILKMAX), OUTLNK(NODMAX,OLKMAX)
C
C     Index of the reservoir release rule curve which this node references.
CCC      INTEGER*4      RRN(NODMAX)
C     Release rule site data for one node...
C
      REAL*4      CAPN(NODMAX)         !Reservoir or Lake capacity
      REAL*4      STOINT(NODMAX)       !Initial Storage
      REAL*4      NELEV(NODMAX)        !Node elevation
      REAL*4      NODLEN(NODMAX)       !Node length (res & lakes)
C
      INTEGER*4   MASTER_NODE(NODMAX)  !Subsystem group number
      INTEGER*4   N_STAT(NODMAX)       !Current display status
      INTEGER*4   SN_STAT(NODMAX)      !Simulation display status
      INTEGER*4   RULESITE(NODMAX)     !=NodeID for RuleSite
                                       !reservoirs in a same group have the same RuleSite
      INTEGER*4   POWERNODE(NODMAX)    !=NodeID for hydropower
      INTEGER*4   PUMPNODE(NODMAX)     !=NodeID for pump?

      LOGICAL*1   GAGENF(NODMAX)         !Is node a gage?
!      LOGICAL*1   GAGEAF(NODMAX)         !Is node a gage?

      LOGICAL*1   NATLAK(NODMAX)       !Is node a natural lake?
      LOGICAL*1   GWNODE(NODMAX)       !Is node a groundwater node?
      LOGICAL*1   DMDNODE(NODMAX)      !Is node a demand node?
      LOGICAL*1   RESVNODE(NODMAX)     !Is node a reservoir?
      LOGICAL*1   WETLAND(NODMAX)      !Is node a wetland?
      LOGICAL*1   EnvFlwNode(NODMAX)  !Is it an Environmental Flow Demand node added by Evgenii 110513
C
C     Node threshold data (lower and upper)
!      REAL*4     N_THRESHOLD(2,MAX_THRESHOLD)
C
C
C     Number of gage sites in system,
C     pointer identifying flow file location or order of the ith gage
C     as defined in system definition and simulation sequence.
C     number of records in replicate,
C     number of replicates in flow file,
C     number of gages in flow file,
C     number of records in replicate in unformatted flow file,
C     array to read flow file.
      INTEGER*4      NGAUGE,tDemSeries
      INTEGER*4      GAGEND(GAGMAX)
C
C ---------------------------------------------------------------
      COMMON /NODE_1/IX,IY,TOTIN,TOTOUT,INLINK,OUTLNK,
     1               CAPN, NODLEN, STOINT, MASTER_NODE, N_STAT,
     2               GAGENF,  NATLAK, PowerNode, PumpNode,
     3               GWNODE, DMDNODE, RESVNODE, WETLAND,
     4               NELEV, SN_STAT, RULESITE,EnvFlwNode !GageAF,
C
      COMMON /NODE_2/NNAME, NodeID    !Zhou+: NodeID
C
      COMMON /GAGE_1/NGAUGE, GAGEND
C
      COMMON/SYSTEM_1/TNODES,NODES,LINKS,SYSPWR,NodSeq,FFFlag,cons_node, 
     1				seep_node,tDemSeries,SEQ_NOD,FFSN

!
C+IRAS/RPA+IRAS/RPA+IRAS/RPA+IRAS/RPA+IRAS/RPA+IRAS/RPA+IRAS/RPA+IRAS/RPA
C
C      Include file FLOSIM.INC
C
C      File contains declarations of variables necessary to run
C      FLOw SIMulation routines in an Interactive River Simulation
C      Program.
C  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C
C            - beginning and ending node storage
      REAL*4   BSTO(NODMAX),                      ! beginning node storage
     1         ESTO(NODMAX),                      ! ending node storage
C            - ending link storage volumes
     1         TOT_LVOL(LNKMAX),
     1         LNKVOL(LNKMAX),
     1         LastLinkVOL(LNKMAX),
C            - total evaporation losses,
     3         TEVAPN(NODMAX),                    ! total evaporation losses
C            - flow at the beginning and end of link,
     5         BQLN(LNKMAX),
     5         EQLN(LNKMAX),
C            - total reservoir release,
     6         TOTREL(NODMAX),
     7         CONSUMPTION(NODMAX),
C            - node total inflow.
     8         INFLOW(NODMAX),
               !inflow within a time-step
     &         InflowTS(NODMax),
C            - node total seepage.
     9         TSEEPL(NODMAX),
     &         ReleaseTS(nodmax),TOT_SUPL_R1(nodmax) !Evgenii 100608 added ReleaseTS (subtime step outflow of node),
      Real*4 YearQLN(LNKMAX) !Annual flow through link
	  Real*4	 TotSrcRelease(MXSUPLY,nodmax)	 
	  Real*4 SrcResRelease(nodmax,nodmax)
		 
C
C      Unregulated natural flow at node
C      Unregulated gage flow at gage site,
C      incremental flow at each node - down from all upstream nodes,
C      cummulative local flow between a site and its upstream reservoirs.
      REAL*4 QINN(NODMAX),NQINN(NODMAX),AddQINN(NODMAX)
C
C     Beginning and Ending Link Flow
      REAL*4   DBQLN(LNKMAX) , DEQLN(LNKMAX)
C
C
      
!Zhou+
!     for gage data
!      INTEGER*4 GageID(GagMax, 2, NodMax)          !for every node

!Evgenii 2009 21 took out *8 after integer and replaced with  ::      
!	  real ::    GageMultiplier(GagMax, 2, NodMax)

!     Policy for system evaporation
      INTEGER*4 nSysEvap, iPSysEvap
      real*4 PolicySysEvap(2,MaxPolicies)     !begin and end date: 0.000-364.999
      real*4 SysEvap(MaxPolicies)
!     Policy group for each simulation year
      INTEGER*4 PolicyGrpID(YRMAX)
!     policy data for a given policy group
      real*4    NodePolicyBeg(MaxPolicies,MaxPolicyTypes,NODMAX)
      real*4    NodePolicyEnd(MaxPolicies,MaxPolicyTypes,NODMAX)
      real*4    LinkPolicyBeg(MaxPolicies,MaxPolicyTypes,LNKMAX)
      real*4    LinkPolicyEnd(MaxPolicies,MaxPolicyTypes,LNKMAX)
      INTEGER*4 NodePolicy0(MaxPolicyTypes,NODMAX) !current node policy for a policyGrp
      INTEGER*4 LinkPolicy0(MaxPolicyTypes,LNKMAX) !current link policy for a policyGrp
      LOGICAL*1 NodePolicyChg(MaxPolicyTypes,NODMAX) !current node policy changed?
      LOGICAL*1 LinkPolicyChg(MaxPolicyTypes,LNKMAX) !current link policy changed?

!     Group of reservoirs
      INTEGER*4 nGrpResv                    !number of groups
      INTEGER*4 nResvInGrp(RNMAX)           !number of reservoirs in a group
      INTEGER*4 RuleSiteID(RNMAX)           !NodeID for rule site
      INTEGER*4 Rule_PTS(RNMAX)             !Number of release rule volumes
      real*4  Res_Rule(8,ZNMAX,RNMAX)     !release rule !evgenii changed to real
                ! 8-beg and end: Min_Stor, Min_rel, max_stor, and Max_rel
      real*4  lastRelease(nodmax)  !release of reservoir group evgenii
      real*4  RelOverMax(nodmax)   !Release over max release from reservoir Evgeni
	real::  RefilTrig(nodmax)    !Evgenii 100708 added RefilTrig for non rule site group reserviors dictating at what
								 !rule reservior storage they can be start to be refilled (in percentage)
!for balance function
      	
      INTEGER*4 ResvIDInGrp(GRMAX,RNMAX)
	  INTEGER*4 UseResVol(GRMAX,RNMAX)
      INTEGER*4 GrpVOL_PTS(RNMAX),BalMethod(GRMAX,RNMAX)          !Number of group volume storages
      real*4    GrpVol(IAGMAX,RNMAX)        !group volumes
      Real*4    RuleResVol(IAGMAX,GRMAX,RNMAX) !Volume of rule reservoir for 2nd balance method Evgenii added 100727
	real*4    BalanceVol(IAGMAX,GRMAX,RNMAX) !(points,nResv,nGrp), each reservoir balancing
	
      
!!Evgenii---------------------------
 !     integer*4 flwdatamax !,flwfactmax
!      parameter (flwdatamax = 73000) !365 days x 200 years
!      parameter (flwfactmax = 10000*gagmax)
 !     real, allocatable, dimension(:,:):: flowdata	  !Evgenii 101014 made flowdata allocatable
!      real  flowfactor(gagmax,12)
 
!Storage node performance variables
	logical:: sto_perf_node(nodmax),last_T_fail(thres_pts_max,nodmax)
	!logical*1 flow_perf_node(nodmax)
	real:: thres_limit(thres_pts_max,nodmax),min_store_reach(nodmax)  
	INTEGER*4 perf_pts(nodmax),nStoreFail(thres_pts_max,nodmax),
     &		  nTime_Steps_Tot(thres_pts_max,nodmax),
     &		  max_sto_fail_dur(thres_pts_max,nodmax),
     &		  nTime_Steps(thres_pts_max,nodmax)
!Flow demand node performance variables		
	real:: Annual_DMD(nodmax),Annual_SHRTG(nodmax),annualSI(nodmax),
     &	   TSflw_DEFICIT(nodmax), !annualSDshrtg(nodmax)
     &      TsSIsum(nodmax),TS_SHRTG(nodmax),TS_surp(nodmax),    !annualSDdem(nodmax)
     &	   Ts_DMD_Sum(nodmax),minYearSD(nodmax),TSflw_Surplus(nodmax) 

! Demand Reduction variables		  !Evgenii added 100719
	INTEGER*4 DemSourceID(nodmax),DemRedPts(nodmax)              
	real::    Dem_Thres_limit(thres_pts_max,nodmax),
     &    	DemRedAmt(thres_pts_max,nodmax)
	logical:: DemRed_node(nodmax)
!  Common 
	COMMON /Performance/ sto_perf_node,last_T_fail,thres_limit,
     &	min_store_reach,perf_pts,nStoreFail,nTime_Steps_Tot,
     &	max_sto_fail_dur,nTime_Steps,minYearSD !,flow_perf_node

	COMMON /DemRed/ DemSourceID, Dem_Thres_limit, DemRedAmt,
     &				DemRed_node,DemRedPts	 
	COMMON /DMD_Performance/ Annual_DMD,Annual_SHRTG,annualSI,
     &	   TSflw_DEFICIT,TS_SHRTG,TS_Surp,			!annualSDshrtg,annualSDdem,
     &	   Ts_DMD_Sum,TsSIsum,TSflw_Surplus
C
      COMMON / BL12 / BSTO,ESTO, TSEEPL,TEVAPN, BQLN,EQLN,TOTREL,
     1                INFLOW,InflowTS, LNKVOL,LastLinkVOL,
     2                DBQLN,DEQLN, TOT_LVOL, CONSUMPTION,
     3                QINN, NQINN, AddQinn,ReleaseTS,TOT_SUPL_R1,
     4				  YearQLN
      COMMON /Policy_1/ PolicyGrpID
      COMMON /SysEvap_1/ PolicySysEvap, SysEvap, nSysEvap,iPSysEvap
      COMMON /PolicyNode_1/ NodePolicyChg, NodePolicy0
      COMMON /PolicyNode_2/ NodePolicyBeg, NodePolicyEnd
      COMMON /PolicyLink_1/ LinkPolicyChg, LinkPolicy0
      COMMON /PolicyLink_2/ LinkPolicyBeg, LinkPolicyEnd

      COMMON /GrpRes_1/nGrpResv,nResvInGrp,RuleSiteID,Rule_PTS,Res_Rule,
     1                 lastRelease,RelOverMax,RefilTrig         !Evgenii 100708 added RefilTrig	 ,flowdata

      COMMON /balance_1/GrpVOL_PTS,ResvIDInGrp,UseResVol,GrpVol,
     & 	  BalanceVol
      COMMON /Gage_1/  RuleResVol,BalMethod !flowdata !GageMultiplier,!GageID,
    	COMMON /source/ TotSrcRelease, SrcResRelease !Evgenii 110715

C
C-IRAS/RPA-IRAS/RPA-IRAS/RPA-IRAS/RPA-IRAS/RPA-IRAS/RPA-IRAS/RPA-IRAS/RPA
C
