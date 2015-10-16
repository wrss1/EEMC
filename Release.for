!	 Interactive River-Aquifer Simulation 2010 (IRAS-2010)
!	 Copyright (c) 2009, 2010 by University College London, Cornell University
!	 Authors:
!	 Daniel P. Loucks (dpl3@cornell.edu), Marshall Taylor, Peter French, Huicheng Zhou,
!	 Evgenii Matrosov (evgenii.matrosov@ucl.ac.uk), Julien Harou (j.harou@ucl.ac.uk),
!	 This program is free software under the General Public Licence, GPL (>=v2)
!	 Read the 'GPL License.txt' file distributed with this source code for a full license statement.
!
!	 **********************************************************************************************
!
C
      SUBROUTINE RELEASE(STORAGES, DEMANDS, TARGETS)
C
C      USE: Obtain calculation of target releases from surface water
C           storage nodes based on beginning of day storages.
C           This routine: 1) sets TARGETS equal to the demand (supplemental -Evgenii) release
C           from a reservoir and adjusts a temporary arrary containing
C           storages accounting for that release, 2) uses the new
C           storage to calculate the minimum allowable release from
C           each reservoir, 3) solves the group release and balancing
C           functions for user defined group release rule sites, and
C           4) adjusts TARGETS to account for supply-driven releases
C           and the minimum release requirements.
C
C           NOTE:  This routine is called once per simulation time step.
C
C      INPUT:
C              Real*4     STORAGES() - Current storage at each node
C              Real*4     DEMANDS()  - Today's "demanded" supplemental
C                                      flows for each source node.
C
C      OUTPUT:     Real*4   TARGETS() - Target releases at all nodes
C                           (Targets are by definition zero at non-
C                            storage nodes and at groundwater nodes)
C                            Units are internal volumes.
C
C      Modifications (yymmdd):
C	 940227: PNF	Check for NODE=0 in loops 60,100
!      000610: PL	    Last change before Evgenii:  PL 
!      090805: EM     Added maxflow for reservoirs  
!
!------------------------------------------------------------------------
C
      IMPLICIT NONE
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  Inpput:
       REAL*4      STORAGES(NODMAX),DEMANDS(NODMAX)
!  Output:
	 real*4      TARGETS(NODMAX)
C
C      Functions:
       INTEGER*4   FINTERP, GT_GRP_REL, BALANCE_GRP
C
C      Local variables:
       INTEGER*4    NUMPTS, ST, I, NODE
       REAL*4       TMP_STOR(NODMAX),MIN_REL(NODMAX),MAX_REL(nodmax)
       REAL*4       OUTFLOW, GRP_RELEASE, TRIAL,max_Q,AVAIL_REL(nodmax)
	 REAL*4       MAX_RELTS(nodmax)
!------------------------------------------------------------------------

C     Copy storages into local working array
      DO NODE = 1, TNODES
         TMP_STOR(NODE) = STORAGES(NODE)
         TARGETS(NODE)  = 0.0
         MIN_REL(NODE)  = 0.0
	   MAX_REL(NODE)  = 0.0  !Evgenii 091104 added MAX_REL
      ENDDO
C
C     Determine the demand and minimum target releases for eac
      DO NODE = 1, TNODES
C
!***    Demand Releases for all Source Reservoirs.
        IF(RESVNODE(NODE) .AND. TMP_STOR(NODE).GT.0.0
     1     .AND. DEMANDS(NODE).GT.0.0)THEN                    
           TARGETS(NODE) = MIN(TMP_STOR(NODE),DEMANDS(NODE))             !Evgenii - Targets becomes either the rest of the storage or supplimental release 
           PastSTSDemsRelease(node)=TARGETS(NODE)
C          Adjust TMP_STOR for later use in supply-driven release eqns.
           TMP_STOR(NODE) = TMP_STOR(NODE) - TARGETS(NODE)               !Evgenii - Temp storage adjusted for targets 
        ENDIF

C       Lakes, Independent Reservoirs and Group Reservoir Minimum Releases
        IF((RESVNODE(NODE).OR.NATLAK(NODE)).AND.TMP_STOR(NODE).GT.0)THEN
           IF (NARVO_PTS(NODE).LT.2) GO TO 25
           IF (LAKEQ_PTS(NODE).LT.2) GO TO 25
            NUMPTS = NARVO_PTS(NODE)
            IF (NUMPTS .NE. LAKEQ_PTS(NODE)) GO TO 25
            ST = FINTERP(3,NODE_VOL(1,NODE),NODE_LAKEQ(1,NODE),
     1                   NUMPTS,TMP_STOR(NODE),OUTFLOW)
            IF(ST.NE.SUCCES)THEN
                WRITE(*,*)'FAILED RELEASE FUNC INTERP, NODE: ',NODE
                GO TO 25
            ENDIF
            MIN_REL(NODE)=MAX(0.0, OUTFLOW*DayPerTS)
        ENDIF

!	  Independent Reservoirs and Group Reservoir Max Releases added by Evgenii 091020 
!	  This max release code only works when the reservior releases the same amount each sub-time step
!	  if the iterative step_deficit calculation is used (the non EnvFlow one, see releasefromdeficit subroutine)
!       this will not release the full maximum release.     	  
	  IF(RESVNODE(NODE).AND.TMP_STOR(NODE).GT.0)THEN 
           IF (NARVO_PTS(NODE).LT.2) GO TO 25
           IF (LAKEQ_PTS(NODE).LT.2) GO TO 25
            NUMPTS = NARVO_PTS(NODE)
            IF (NUMPTS .NE. LAKEQ_PTS(NODE)) GO TO 25
            ST = FINTERP(3,NODE_VOL(1,NODE),Node_MaxQ(1,NODE),
     1                   NUMPTS,TMP_STOR(NODE),max_Q)
		  IF(ST.NE.SUCCES)THEN
               WRITE(*,*)'FAILED RELEASE FUNC INTERP, NODE: ',NODE
                GO TO 25
            ENDIF
		  MAX_RELTS(NODE)=MAX(0.0,max_Q*sysstat(NPER))
		  MAX_REL(NODE)=MAX(0.0, max_Q*DayPerTS)		 
	  ENDIF
 25     CONTINUE
      ENDDO !Evgenii - end do for all nodes
      
!89	Continue
C
C     Now solve for TRIAL releases from reservoirs in one or more
C     reservoir groups. Note that group release is based upon adjusted
C     storages (after the demand releases) and that the solution is
C     performed in simulation sequence order (giving "priority" to
C     rulesites earlier in the simulation sequence
	DO 60 I = 1, TNODES
         NODE = NODSEQ(I)
         IF (NODE .EQ. 0) GO TO 60
         !Now for rules sites	   
	   IF(RULESITE(NODE).eq.NodeID(Node))THEN    

            ST = GT_GRP_REL(NODE,TMP_STOR,GRP_RELEASE)			!Evgenii - calls reservoir rules release, with TMP_STOR adjusted for TARGETS 
           IF(ST.NE.SUCCES) GO TO 45 !Evgenii changed 50 to 45, so if no rule table, balance still performed 100728!Grp_release in volume/sub-time step units evgenii
C           If rulesite is identified as a demand site then add
C           the required supplimental release
		  
		  !Evgenii commented out line below so that DEMANDS are 
		  !not double counted for rule sites (TMP_STOR takes DEMANDS
		  !into account in the form of TARGETS. Double counting
		  !becomes obvious in BALANCE_GRP and further down in RELEASE
!		  GRP_RELEASE = GRP_RELEASE  + DEMANDS(NODE)			
            
		  !Reservoirs should be balanced even if no rule release
	      !(there can still  be DEMANDS release) so Evgenii changed
		  !GRP_RELEASE<=0.0 to GRP_RELEASE<0.0 in line below 100127

45		  IF(GRP_RELEASE<0.0)GO TO 50				
		  ST = BALANCE_GRP(NODE,TMP_STOR,GRP_RELEASE,STORAGES)
            !Evgenii - total group Storage now decreased by grp-release

		  IF(ST.NE.SUCCES)GO TO 50
         ENDIF
50       CONTINUE
60     CONTINUE
C
C      Now combine the TARGET and TRIAL releases and reconcile
C      with minimum releases...
       DO 100 I = 1, TNODES
          NODE = NODSEQ(I)
          IF (NODE .EQ. 0) GO TO 100
          TRIAL   = STORAGES(NODE) - TMP_STOR(NODE) 
          if(nodeid(node)==201902 .and. sysstat(19)==4.and.
     &      sysstat(20)==167) then
          continue
          end if
          TRIAL   = MAX(0.0, TRIAL) 
          TRIAL   = MAX(0.0, TRIAL - TARGETS(NODE)) 
          TARGETS(NODE) = TARGETS(NODE) + TRIAL 
		!Evgenii put in max release below
          if(RESVNODE(NODE) ) then
     			if (sysstat(sstep)<sysstat(NSUBTT)) then
			   TARGETS(NODE) = min(TARGETS(NODE),MAX_REL(NODE))
			else if (sysstat(sstep)==sysstat(NSUBTT)) then
			   AVAIL_REL(NODE) = MAX(0.0,MAX_RELTS(NODE)-TOTREL(NODE)) !TOTREL is cumulative release upto this subtime step,Max_Q is i mil mil m3/time step
			   AVAIL_REL(NODE) = MAX(0.0, AVAIL_REL(NODE))
			   TARGETS(NODE) = min(TARGETS(NODE),AVAIL_REL(NODE))
			endif
		endif
	
			TARGETS(NODE) = MAX(TARGETS(NODE),MIN_REL(NODE)) 
              TARGETS(NODE) = MAX(TARGETS(NODE),0.0) 
!	if(RESVNODE(NODE) )then 
!		step_deficit(node)=step_deficit(node)+TARGETS(NODE)
!	endif
100    CONTINUE
C
 9999 RETURN
      END

C     
!********************************************************************************
C
      INTEGER*4 FUNCTION GT_GRP_REL(RULENODE,STORAGES,RELEASE) !,
!     & REALSTO)
C
C       Code pulled out of RELMRT.FOR and restructured. 940121
C
C     Inputs:  RULENODE - INTEGER*4 - The rulesite at which the group
C                                     release policy and balancing function
C                                     is defined.
C              STORAGES - REAL*4() -  A vector of storages at each node
C                                     in the schematic.
C              RELEASE - REAL*4    -  The trial release to be made from
C                                     the reservoir group. Note that RELEASE
C                                     is returned in FLWSIMs internal flow
C                                     units (volume/simulation time step).
C      Note:   If current group storage is not within the lower and upper
C              storage bounds of a user-defined release zone then the
C              function returns with status SUCCES and RELEASE=0.
C	 Last change before Evgenii:  PL   31 May 2000    4:10 pm
!-------------------------------------------------------------------------------------

C     Arguments:
      implicit none
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  Input:
      INTEGER*4 RULENODE           !={1, TNodes}: seqence no.
      REAL*4    STORAGES(nodmax)!,REALSTO(NODMAX) !Evgenii 090729 changed * to nodemax
      !COMMON: CurrentDay
      !        NodePolicy0(Rule0,RuleNode),NodePolicyBeg(PP,Rule0,RuleNode),NodePolicyEnd()
      !        RuleSiteID(iGrp), Rule_PTS(iGrpResv), Res_Rule(8,ZNMax,RNMAX)
      !        RuleSite(Node), NodeID(Node)
!  Output:
      REAL*4    RELEASE
C
C     Local Variables:
      INTEGER*4 iGrpResv, ST, I, iNode, iZone, iItem, pp,nres
      REAL*4    GRP_STORAGE, FRACT_P, FRACT_S, DELTA_STOR
      REAL*4    MIN_REL_I, MIN_STOR_I, MAX_REL_I, MAX_STOR_I,
     1          MIN_REL_F, MIN_STOR_F, MAX_REL_F, MAX_STOR_F
      REAL*4    MIN_STOR, MAX_STOR, MIN_REL, MAX_REL
      REAL*4    TPolicyDays, CurrentPDay

!-------------------------------------------------------------------------------------
C     Initialization
      GT_GRP_REL = FAIL
      RELEASE = 0.0
	!Evgenii added line below  090825 to check rulesite
	IF(RULESITE(RULENODE)==0)GO TO 9999 
C     Return with FAIL status if specified site not a rulesite.
      IF(RULESITE(RULENODE).ne.NodeID(RuleNode))GO TO 9999
      !find serial no. for the Group
      do i = 1, nGrpResv
        IF(RuleSiteID(i)==RULESITE(RULENODE))then
          iGrpResv = i; exit
        END if
      end do
!***  Compute the current group storage.
       GRP_STORAGE = 0.0
       !DO I = 1, nResvInGrp(iGrpResv) !Evgenii 090728 took out this erroneus do loop, 
	 !do loop below suffices for GRP_STORAGE calculation
        nres=0
        do iNode = 1, TNodes
           IF(RuleSite(RuleNode).eq.RuleSite(iNode))THEN
             GRP_STORAGE = GRP_STORAGE + STORAGES(iNode)
             nres=nres+1
             if (nres >= nResvInGrp(iGrpResv))exit !Evgenii added to shorten loop run time 120215
           ENDIF
        end do
       !ENDDO

C     Completed fraction of total days for current policy
      !get total days for current policy: TPolicyDays
      pp = NodePolicy0(Rule0,RuleNode)
      TPolicyDays = NodePolicyEnd(pp,Rule0,RuleNode)
     &             -NodePolicyBeg(PP,Rule0,RuleNode)
      
	!get current policy day
      CurrentPDay=MAX(0.,SysStat(Time)-NodePolicyBeg(PP,Rule0,RuleNode))
     
      !fraction
      FRACT_P = MIN(1.,CurrentPDay/TPolicyDays)
      !Interpolate release
      DO iZone = 1, Rule_PTS(iGrpResv)
            MIN_STOR_I = Res_Rule(1,iZone,iGrpResv)
            MIN_REL_I  = Res_Rule(2,iZone,iGrpResv)
            MAX_STOR_I = Res_Rule(3,iZone,iGrpResv)
            MAX_REL_I  = Res_Rule(4,iZone,iGrpResv)
            MIN_STOR_F = Res_Rule(5,iZone,iGrpResv)
            MIN_REL_F  = Res_Rule(6,iZone,iGrpResv)
            MAX_STOR_F = Res_Rule(7,iZone,iGrpResv)
            MAX_REL_F  = Res_Rule(8,iZone,iGrpResv)
            MIN_STOR = MIN_STOR_I*(1.0-FRACT_P) + MIN_STOR_F*(FRACT_P)
            MAX_STOR = MAX_STOR_I*(1.0-FRACT_P) + MAX_STOR_F*(FRACT_P)
            MIN_REL  = MIN_REL_I*(1.0-FRACT_P) + MIN_REL_F*(FRACT_P)
            MAX_REL  = MAX_REL_I*(1.0-FRACT_P) + MAX_REL_F*(FRACT_P)
            IF(GRP_STORAGE.GE.MIN_STOR.AND.GRP_STORAGE.LE.MAX_STOR)THEN
                 DELTA_STOR = MAX_STOR - MIN_STOR
                 FRACT_S = 0.
                 IF(DELTA_STOR.GT.0.0)THEN
                   FRACT_S = (GRP_STORAGE-MIN_STOR)/DELTA_STOR
                 ENDIF
                 RELEASE = MIN_REL*(1-FRACT_S)+MAX_REL*FRACT_S
C                Convert user release to internal units
                 RELEASE = RELEASE * DayPerTS 
                 GT_GRP_REL = SUCCES
                 GO TO 9999
            ENDIF
      ENDDO
C
 9999   RETURN
        END

!********************************************************************************
        INTEGER*4 FUNCTION BALANCE_GRP(RULENODE,STORAGES,TOTAL_RELEASE, 
     & REALSTO)

C       Use:     This routine uses the reservoir balancing table specified 
C                for RULENODE, the current reservoir group STORAGES, and 
C                the specified TOTAL_RELEASE to be made from the group to 
C                compute new target storages for each reservoir in the group.
C                Be sure to release only enough water to meet its allocation of
C                TOTAL_RELEASE even if a larger release would be required
C                to reach the balancing table's target storage. 
C
C       Inputs:  RULENODE - INTEGER*4 - The rulesite at which the group
C                                       release policy and balancing function
C                                       is defined.
C                STORAGES - REAL*4()    - The current storage at each node.
C                TOTAL_RELEASE - REAL*4 - The release to be made from the
C                                         reservoir group. This value is
C                                         in volume units.
C
C       Outputs:  STORAGES - The final storages which (if feasible) meet
C                            the group's release requirement and takes water
C                            from the reservoirs whose storage exceeds their
C                            current target storages.
C
C       Function value: The routine will return FAIL if: 
C                        o RULENODE has not been designated as a rulesite,
C                        o nGrpResv.LE.1 or TOTAL_RELEASE.le.0.001
C                       Otherwise, the routine returns SUCCES.
C       
C	  Last change before Evgenii:  PL   17 May 2000   11:25 am
C
!-------------------------------------------------------------------------------------
      implicit none
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      INTEGER*4  RULENODE
      REAL*4     STORAGES(NODMAX), TOTAL_RELEASE,REALSTO !Evgenii changed storages(*) to storages(nodmax)
      !COMMON: ResvNode(),RuleSiteID(i),RULESITE(RULENODE)
      ! nGrpResv,nResvInGrp(RNMAX),ResvIDInGrp(RNMAX),
      ! GrpVOL_PTS(RNMAX),GrpVol(IAGMAX,RNMAX),BalanceVol(points,nResv,nGrp)
!  OUTPUT
      !REAL*4     STORAGES(*)
C  Functions:
      INTEGER*4   FINTERP
C  Locals
      INTEGER*4  ST, NUMPTS
      INTEGER*4  I, J, GRPRESV(NRMAX), TMP_NODE
      INTEGER*4  iGrpResv, RESV_COUNT
      REAL*4   GRPSTOR(IAGMAX),RSVSTOR(IAGMAX,NRMAX)
	Real*4   RuleResStor(IAGMAX,RNMAX)
      REAL*4   Surplus(NODMAX),Deficit(NODMAX),TARGET_STORAGE(NRMAX)
      REAL*4   NEW_STORAGE, TOTAL_STORAGE, TotSurplus, TotDeficit
      REAL*4   NetSurplus, targets(nodmax),RealTRuleVol, partstorage
!-------------------------------------------------------------------------------------
C     Initialization...
      BALANCE_GRP = FAIL
      partstorage=0.0
C     If RULENODE not a rulecurve site then exit routine with FAIL
      IF(RULESITE(RULENODE).ne.NodeID(RuleNode))GO TO 9999
	!Evgenii  replaced <=0.000 with <0.0 below so even if no
	!rule release, balancing function continues
	if (TOTAL_RELEASE<0.0) GOTO 9999

	!Evgenii commented if below 090804, because it didnt allow balancing function to continue-
!
!      if (ResvNode(RuleNode)) then
!        STORAGES(RuleNode) = STORAGES(RuleNode) - TOTAL_RELEASE
!        STORAGES(RuleNode) = MAX(0., STORAGES(RuleNode))
!        BALANCE_GRP = SUCCES
!     end if

C     If no reservoir group, returns fail
      IF(nGrpResv<1)GO TO 9999 !Evgenii made this < instead of <= that was before 040809

      !find serial no. for the Group
      do i = 1, nGrpResv
        IF(RuleSiteID(i)==RULESITE(RULENODE))then	    
          iGrpResv = i; exit
        END if
      end do
C
   
	RESV_COUNT = nResvInGrp(iGrpResv)
      IF(RESV_COUNT.LT.1)GO TO 9999

      !GRPRESV(I)={1,...,TNodes}: serial node No. for reservoir node in the group
      TOTAL_STORAGE = 0.0 
      !Compute group's total storage
	
	!Evgenii added if below so that single reserviors do not need balance tables. 12.04.11
	IF(RESV_COUNT == 1)THEN
		TOTAL_STORAGE=STORAGES(RULENODE)
	else if (RESV_COUNT > 1) then
		do i = 1, RESV_COUNT
		!in this loop put the partial group volume calculation
            do j = 1, TNodes	
	        if (NodeID(j) == ResvIDInGrp(i, iGrpResv)) then
	           GRPRESV(I) = j !Reservoir IDs in the group in the order they are in in the balance table in INP
			   TOTAL_STORAGE = TOTAL_STORAGE + STORAGES(J) !Evgenii changed STORAGES(GRPRESV(I)) to J, the IRASID is needed here (J) not the Userset INP id (GRPRESV(I)) 201112
			   if (UseResVol(i, iGrpResv)==1)then
				 Partstorage=Partstorage+STORAGES(J)
			   end if
                 if (J==RULENODE) then   !Evgenii 100728 put this in for data 2nd balancing method, 121120 Evgenii changed if condition to J==RULENODE to fix Node ID bug
				   RealTRuleVol=storages(j)
			   end if
	           exit
	        end if
		 end do
	   end do
	end if
C
C     If only one reservoir in the group
      IF(RESV_COUNT.EQ.1)THEN
        TMP_NODE = RULENODE !GRPRESV(1)  !Evgenii replaced GRPRESV(1) with RULENODE out so balance tables are not needed for single reserviors 12.04.11
        STORAGES(TMP_NODE) = STORAGES(TMP_NODE) - TOTAL_RELEASE 
	  STORAGES(TMP_NODE) = MAX(0.0, STORAGES(TMP_NODE)) 
        BALANCE_GRP = SUCCES
        GO TO 9999
      ENDIF
C
C     Set balancing functions for interpolation routine.
      !Set NUMPTS, GRPSTOR(i), RSVSTOR(i,j)
      NUMPTS = GrpVOL_PTS ( iGrpResv )
      DO I = 1, NUMPTS !I is the balance volume
         GRPSTOR(i) = GrpVol(I,iGrpResv)
         DO J = 1, nResvInGrp(iGrpResv) !J is the reservoir number
           RSVSTOR(i,j) = BalanceVol(i,j,iGrpResv) 
	     RuleResStor(i,j)=RuleResVol(I,j,iGrpResv)  !Added by evgenii for 2nd balance method 100728
         ENDDO
      ENDDO

      NEW_STORAGE = TOTAL_STORAGE - TOTAL_RELEASE 
      RealTRuleVol=RealTRuleVol - TOTAL_RELEASE
      Partstorage=Partstorage - TOTAL_RELEASE
	TotSurplus = 0.0;  TotDeficit = 0.0
      DO I = 1, RESV_COUNT
         IF(GRPRESV(I).GT.0 .AND. GRPRESV(I).LE.NODMAX)THEN
           TARGET_STORAGE(I) = STORAGES(GRPRESV(I))
           if (BalMethod(1,iGrpResv)==0)then		!Evgenii added if statement to choose balance method 100728
			ST = FINTERP(3,GRPSTOR,RSVSTOR(1,I),NUMPTS,NEW_STORAGE,	   !NEW_STORAGE is the total group storage-rule release for balance method 0
     1                   TARGET_STORAGE(I))
			!GROSTOR - Grpstorage balance table, RESVSTOR reservoir storage balance table,NEW_STORAGE - the total group storage-rule release

		 else if (BalMethod(1,iGrpResv)==1)	then   !Evgenii added if statement to choose balance method 100728
     			ST =FINTERP(3,RuleResStor(1,I),RSVSTOR(1,I),NUMPTS,	   !RealTRuleVol is the group storage-rule releasefor balance method 1
     1                   RealTRuleVol,TARGET_STORAGE(I))
 		 else if (BalMethod(1,iGrpResv)==2)	then   !Evgenii added if statement to choose balance method 130731
     			ST =FINTERP(3,GRPSTOR,RSVSTOR(1,I),NUMPTS,	   !partstorage is the partial group storage for balancing method 2
     1                   partstorage,TARGET_STORAGE(I))                                  
		else 
			goto 9999
		 end if
           Surplus(I) = !Evgenii, I is the order of the Reservoirs in the balance table. 
     &           MAX(0.0, STORAGES(GRPRESV(I))-TARGET_STORAGE(I)) !121120 Evgenii changed Surplus(I) to Surplus(GRPRESV(I)), bug fix to make Surplas ID consistant with the rest of the formula (INP order ID)
           !TotSurplus = TotSurplus + Surplus(i) Evgenii commented TotSurplus and TotDeficit out 100314, no longer used
           Deficit(I) = 
     &           MAX(0.0,TARGET_STORAGE(I)-STORAGES(GRPRESV(I)))  !121120 Evgenii changed deficit(I) to deficit(GRPRESV(I)), bug fix to make Surplas ID consistant with the rest of the formula (INP order ID)
           !TotDeficit = TotDeficit + Deficit(i)
         ENDIF
      ENDDO
      !NetSurplus = MAX(0.0, TotSurplus - TotDeficit) !Evgenii commented TotSurplus and TotDeficit out 100314, no longer used
C
C     compute the final storages
      DO I = 1, RESV_COUNT
         TMP_NODE = GRPRESV(I)
         IF(TMP_NODE.GT.0 .AND. TMP_NODE.LE.NODMAX)THEN
           !Evgenii 091010 made it so rulesite only releases group release 
		 !Rulesite should only release Rule release, not balancing release
		  IF(TMP_NODE.EQ.RULENODE)THEN								
               STORAGES(TMP_NODE) = STORAGES(TMP_NODE)-TOTAL_RELEASE
		  else
		  STORAGES(TMP_NODE) =STORAGES(TMP_NODE)-
     &						  Surplus(I) !*NetSurplus/Total_RELEASE
            !100127 Evgenii commented out '*NetSurplus/Total_RELEASE' because 
	      !other wise balancing would not function without rule release
		  !(Total_RELEASE). Now total surplus is attempted to be released
		  !in addition to any DEMANDS release calculted in RELEASE but 
		  !later constrained by MIN and MAX releases
		  endif
		  STORAGES(TMP_NODE)=MAX(0.0,STORAGES(TMP_NODE))
         ENDIF
      ENDDO
      BALANCE_GRP = SUCCES
 9999 RETURN
      END



!*************************************************************************
       subroutine ReleaseFromDeficit(SUPL_RELEAS,TOT_SUPL_R,DSTO,
     &   PotSrcRel,potQINN) 
    
! Compute this step's deficit, additional releases from sources for deficits at demand nodes
! Seperated from FLWSIM.FOR
      IMPLICIT NONE
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  CALL: NONE
!  INPUT
      REAL*4 DSTO(NODMAX),lastIN(nodmax),PotSrcRel(NODMAX)
      !COMMON: TNODES,CAPN(NN), QInn(nn)
      !        DMD_TARG(NN),DMD_T_co(NN),STEP_DEFICIT(NN)
      !        SUPLY_PTS(NN), SUPL_FRAC(I,NN),
!  OUTPUT
      REAL*4 SUPL_RELEAS(NODMAX), TOT_SUPL_R(NODMAX)
      !COMMON: TOT_SUPL_R(NN)
!  Local
!      Functions:
       INTEGER*4   FINTERP, GT_GRP_REL, BALANCE_GRP

      INTEGER*4 i, j, NN, S_NODE, STEP,STEPS_LEFT,GrpRruleID,ls,check
     &  ,glink,LL,LN
	Real*4 EXPECT_INFL(NODMAX),TMP_TARGET,TOT_Target(NODMAX),AvailTSRel
	Real*4 TOT_DEFICIT(NODMAX),RuleStoFrac,Reduction,SrcRelease
	Real*4 AvailYearRel,gagesrcflow,potQINN(NODMAX),gagesrcflowTS(lnkmax)
	Real*4 Step_Def_diff,OrigTemp_STEP_DEFICIT,Temp_STEP_DEFICIT(nodmax),
     &	AVAIL_CAP,AVAIL_CAP_year
	logical*1 success
!-------------------------------------------------------------------------
      step=sysstat(sstep)
      STEPS_LEFT=STEPPRPRD-STEP+1 
	!Evgenii added line above from 1995 code 090807, +1 becaues this is the beginning of the step.
      do ln=1,links !120320 put in do loop (gagesrcflowTS wasnt resetting)
        gagesrcflowTS(ln)=0.0  
      end do
C     Initialize supplemental releases for this step
      DO NN = 1,TNODES
          SUPL_RELEAS(NN) = 0.0
      END DO  
	DO NN = 1,TNODES
          IF (DMDNODE(NN)) THEN
		        !compute deficits
				!For storage nodes				 
				IF (CAPN(NN) .GT. 0.0) THEN     !for storage node
				!Node with storage: Deficit is diff. between target and storage + release from previous subtime step in the current time-step (+releaseTS(nn) added by Evgenii)
       			  if (step > 1) then
					TOT_DEFICIT(NN) = DMD_TARG(NN)-DSTO(NN) 
 !    &					+releaseTS(nn)           !+ release from previous subtime step in the current time-step (+releaseTS(nn)) added by Evgenii 100608 to account for lag
					
					!Evgenii 100708 added check to see if non-rule site resrvrs in group are above filling trigger 
					if (RESVNODE(nn).and.RuleSite(nn)/=NodeID(nn)
     &					.and.TOT_DEFICIT(NN)>0.0)then
						 !calculate storage frac of rule site reservor
						 do i=1,tnodes
							if(NodeID(i)==RuleSite(nn)) GrpRruleID=i
						 end do		
						 RuleStoFrac=DSTO(GrpRruleID)/CAPN(GrpRruleID)
		                 !If rulesite storage frac less than refill trigger set deficit to 0
						 if(RuleStoFrac<RefilTrig(nn)) then
							TOT_DEFICIT(NN)=0.0
						 end if
					end if	 	
					STEP_DEFICIT(NN) = MAX( 0.0, TOT_DEFICIT(NN)) 
! 			        STEP_DEFICIT(NN) = TOT_DEFICIT(NN)/steps_left	
				  end if
				!For non-storage nodes
				else
					!Determin demand reduction due to low supply, if any
					Reduction=0.0
					if(DemRed_node(nn)) then
						 do i=1,tnodes
							if(NodeID(i)==DemSourceID(nn))then
								do j=1,DemRedPts(nn)
     &								!must add check so that source node here was a capacity 
									if(DSTO(i)<(Dem_Thres_limit(j,nn)
     &										 *capn(i)))then
										Reduction=max(reduction,
     &										 DemRedAmt(j,nn))
									end if
								end do 			
							end if								
						  end do
					end if
					if(EnvFlwNode(NN)) then !Evgenii 100513 Added environmental flow step deficit
				!Code below is IRAS-2000 code, Evgenii 100301
				!IRAS 2000 CODE below, doesnt take into accout passive water to reduce deficit	
					  STEP_DEFICIT(NN)=(DMD_TARG(NN)-
     & 					  DMD_TARG(NN)*Reduction-QINN(NN))*DAYPERTS
					  STEP_DEFICIT(NN) = MAX(0.,STEP_DEFICIT(NN))    
					end if
                    !Now the algorithm that does take into account passive water
					if (step > 1.and. .not. EnvFlwNode(NN)) then 
C              Have to wait until after the 1st step before we know what
C              deficit at site NN was.
C              Expected "natural" remaining inflow is average over
C              previous steps of (cumulative inflow - cumulative
C              supplemental release)*steps left.
		      !Natural inflow neglects any attempted extra releases resulting from step-decicits from previous sub-time steps (TOT_SUPL_R), Evgenii 100303
			  !The expected inflow calculation does not yet take into account extra diversion from source demand links 110517 evgenii
						EXPECT_INFL(nn) = STEPS_LEFT*(INFLOW(NN)-
     &                           TOT_SUPL_R(NN))/(STEP-1)
						!Evgenii put in line below to prevent expected inflow from going negative 110926
						EXPECT_INFL(nn)=max(0.0,EXPECT_INFL(nn))
						!DMD_TARG() is  mil m3 per day, converted to mil m3 per time-step by DAYSPRPRD, Evgenii 100303
						TOT_DEFICIT(NN)=(DMD_TARG(NN)-		  
     &						DMD_TARG(NN)*Reduction)*DAYSPRPRD
     &						-INFLOW(NN)- EXPECT_INFL(NN)			
						TOT_DEFICIT(NN) = MAX( 0.0, TOT_DEFICIT(NN)) 
				        STEP_DEFICIT(NN) = TOT_DEFICIT(NN)/steps_left                 
					END IF 					 
				endif 		  		
	         !Step_Deficit calculation is ignoring excess water Julien and Evgenii 090810
C              If a demand node, compute required releases for upstream reservoirs.
C              If a deficit, compute required supplemental releases
				    
				IF (STEP_DEFICIT(NN).GT.0.0)THEN 				 		 
				 Temp_STEP_DEFICIT(NN)=STEP_DEFICIT(NN)				
				 DO I = 1,SUPLY_PTS(NN)
					 SrcRelease=0
					 !TotSrcRelease restes to 0 at at the first subtime step in flwsim.for.
				  IF(SYSSTAT(sstep)==0) TotSrcRelease(I,NN)=0.0	!Evgenii 150711
				    S_NODE = SUPL_NODE(I,NN)    !SUPL_NODE = {1,...TNodes}, not NodeID                
				    IF (S_NODE.GT.0) THEN	
				      !If source node is a non reservoir or lake guage node, reduce the Temp_STEP_DEFICIT by the ammount the gauge can provide, added by Evgenii 120123
				      if (GageNF(s_node) .and. .not.ResvNode(s_node).and..not.
     &				    NatLak(s_node))then								          		         				
		                check=0
		                glink=0
		                !If the demand node is directly downstream of supply node
		                do LL = 1,TOTOUT(s_node)
                           LN = OUTLNK(s_node,LL)
                           if (Nout(ln)==nn) then
                              check=check+1
                              glink=LN
                              exit                                
                           end if
                        end do    
                       !If the demand node is not directly downstream, but link coming out of source gauge is a source link for the demand node
                       if (check==0) then                        
		                 do ll = 1,SUPLY_Links_pts(NN)
		                   if (Nin(SUPL_Link(ll,nn))==s_node)then
                             glink=SUPL_Link(ll,nn)
                             check=check+1
                             exit
                           end if
                         end do
                       end if
                       if (check==1) then
		               !Make sure the demand link responsible for taking the flow from the gauge node is limited by its capacity
		               !Capacities for gauge sources are set by the Link capacity, not in the source lines!!!
		               AVAIL_CAP = MAX(0.0, CapL(glink)- BQLN(glink)
     &		                     -gagesrcflowTS(glink))
		               Avail_Cap_Year=MAX(0.0, CapLYear(glink)- 
     &		                     YearQLN(glink)-gagesrcflowTS(glink))
		               !Update potential contribution from this guage node taking into account annual and TS capacities
		               potqinn(s_node)=MIN(potqinn(s_node)*
     &		                   SUPL_FRAC(I,NN),AVAIL_CAP,Avail_Cap_Year)
				       !Assign the potential guage flow
				       gagesrcflow=min(Temp_STEP_DEFICIT(NN),potqinn(s_node))     
                         gagesrcflow=max(gagesrcflow,0.)  !Added by Evgenii 120214 to prevent from going negative
		               !Reduce temp_step_deficit to reflect potential guage contribution
		               Temp_STEP_DEFICIT(NN)=Temp_STEP_DEFICIT(NN)
     & 	                   -gagesrcflow !potqinn(s_node)
		               Temp_STEP_DEFICIT(NN)=max(Temp_STEP_DEFICIT(NN),0.) !Added by Evgenii 120214 to prevent from going negative
		               potqinn(s_node)=potqinn(s_node)-gagesrcflow
		               potqinn(s_node)=max(potqinn(s_node),0.0)	
		               !Accumulate the total expected inflow into the link resulting from gauge source contribution (to make sure capacity is obeyed)
		               gagesrcflowTS(glink)=gagesrcflowTS(glink)
     &  			               +gagesrcflow   	       
		               end if				              
		               cycle
		              end if !End (GageNF(s_node) .and. .not.ResvNode(s_node)etc..
					    
				     !For reservoir sources			   
			         if (.not. GageNF(s_node)) then		!Evgenii 210312, so g			    
                           !Evgenii 150711 added SrcRelease 
					   SrcRelease=0.0
					   SrcRelease=SUPL_FRAC(I,NN)*Temp_STEP_DEFICIT(NN)
			         			     				  	     
		               !Evgenii 150711 added if below 
		               !if Trail_Release brings the source over its cumulative maximum time-step contribution
		               !reduce Trail_Release so that it brings the contribution to the maximum
	                 if(MaxOutTS(I,NN)>=0.0) then
		                 AvailTSRel=max((MaxOutTS(I,NN)-TotSrcRelease(I,NN))
     &		                 ,0.0)
		                 SrcRelease=min(SrcRelease,AvailTSRel)
	                 end if 
	                 if(MaxOutYear(I,NN)>=0.0) then
			            AvailYearRel=max((MaxOutYear(I,NN)-AnnualOut(I,NN)),0.0)
      			        SrcRelease=min(SrcRelease,AvailYearRel)
				       end if                    
                      !IF (SrcRelease+TotSrcRelease(I,NN)>
!     &			       (MaxOut(I,NN)*DAYSPRPRD))THEN     					
! 			           SrcRelease=MaxOut(I,NN)*DAYSPRPRD-TotSrcRelease(I,NN)
!	      	      END IF					     					         					       					     					   
		             !Now reduce deficit by what you think the reservior can release which is its maximum release rate at its current volume
		             !or its current volume
	                SrcRelease=max(0.0,SrcRelease)							    		
	                SUPL_RELEAS(S_NODE)=SUPL_RELEAS(S_NODE)+SrcRelease
                      !Evgenii 150711 changed SUPL_FRAC(I,NN)*STEP_DEFICIT(NN) to SrcRelease 
	                TotSrcRelease(i,nn)=TotSrcRelease(i,nn)+SrcRelease !SUPL_FRAC(I,NN)*STEP_DEFICIT(NN)
                      AnnualOut(i,nn)=AnnualOut(i,nn)+SrcRelease 
                      !Set OrigTemp_STEP_DEFICIT which will be used later after Temp_STEP_DEFICIT is modified                 				 
		              !Evgenii- Need to write check to make sure sources that use priorities are never also used with contribution factors	        			 
         		      !In the case of source priorities
         		      
         		      !Set how much res is release for this demand node Evgenii added  121017
         		      SrcResRelease(s_node,nn)=SrcRelease
         		      if (srcpriorities(nn)==.true.) then
        		       OrigTemp_STEP_DEFICIT=Temp_STEP_DEFICIT(NN)  
         		       !Reduce Temp_STEP_DEFICIT by the amount that the reservior can release 
       			       Temp_STEP_DEFICIT(NN)=Temp_STEP_DEFICIT(NN)
     &     				 -PotSrcRel(S_node) !*SUPL_FRAC(I,NN)
                       !Make sure it doesnt go below 0
                       Temp_STEP_DEFICIT(NN)=
     &                 max(0.0,Temp_STEP_DEFICIT(NN))
                       !Step_Def_diff is how much of last subtime step's release went to node nn
                       Step_Def_diff=OrigTemp_STEP_DEFICIT
     &                    -Temp_STEP_DEFICIT(NN)         
     				       !How much of last subtime steps do we have left for lower priority nodes
     				       PotSrcRel(S_node)=PotSrcRel(S_node)-Step_Def_diff	  		
				       !Here add the code to limit SUPL_RELEAS                        
                      end if
                     end if	!end if not gauge node                                                
                    END IF !end if S_NODE>0 
				 ENDDO                                     

                   TOT_SUPL_R(NN) = TOT_SUPL_R(NN) + STEP_DEFICIT(NN)
                   NonPropSTEP_DEFICIT(NN)=STEP_DEFICIT(NN)
				END IF !End if Step Deficit
          END IF !End if demand node	
      ENDDO        
	end subroutine
!*************************************************************************
       subroutine PotentialStoSourceRel(Storages,PotSrcRel,potQINN,
     &       success) 
    
! Computes the potential releases from reserviors for this sub time step
! Also computes the potential timestep release for gauge nodes to demand nodes
      IMPLICIT NONE
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  CALL: NONE
!  INPUT

      Real*4 Storages(nodmax)
!  OUTPUT
      Real*4 PotSrcRel(nodmax),potQINN(nodmax)
      logical*1 success
!  Function
      integer*4 FINTERP  
!  Local
      INTEGER*4 NUMPTS, ST, nn
      Real*4 MAX_RELTS, MAX_REL,max_Q

!-------------------------------------------------------------------------
      MAX_RELTS=0.;MAX_REL=0.;max_Q=0.
      success=.false.
      do nn=1,tnodes
         if(storages(nn)>0.0 .and. capn(nn)>0.)then
            if(NARVO_PTS(nn)<2) then 
                go to 25
            end if
            if(LAKEQ_PTS(nn)<2) then 
                go to 25
            endif
            NUMPTS = NARVO_PTS(nn)
            if (NUMPTS /= LAKEQ_PTS(nn)) then
                go to 25
            end if
            ST = FINTERP(3,NODE_VOL(1,nn),Node_MaxQ(1,nn),
     1               NUMPTS,Storages(nn),max_Q)
	      if(ST/=SUCCES)then
               write(*,*)'FAILED RELEASE FUNC INTERP, NODE: ',nn
               go to 35
            end if		  
		      !Total max release over entire time step
		      MAX_RELTS=max(0.0,max_Q*sysstat(NPER))
		      !Max release evenly divided into subtime step
		      MAX_REL=max(0.0, max_Q*DayPerTS)		 
	  
             !if not on the last subtime step
             if (sysstat(sstep)<sysstat(NSUBTT)) then
               PotSrcRel(nn) = MAX_REL
              !If on the last subtime step
             else if (sysstat(sstep)==sysstat(NSUBTT)) then
               PotSrcRel(nn) = MAX(0.0,MAX_RELTS-TOTREL(nn)) !TOTREL is cumulative release upto this subtime step,Max_Q is i mil mil m3/time step
            endif
             PotSrcRel(nn) = min(Storages(nn),PotSrcRel(nn))
	       PotSrcRel(nn) = max(PotSrcRel(nn),0.0)
	  else if (GageNF(nn)) then
	       potqinn(nn)=qinn(nn)* DayPerTS  !Potential time step release to a demand node if gauge is a source, Evgenii 120124
	       potqinn(nn)=max(potqinn(nn),0.0)
	  end if	  	  	
	end do	
25	success=.true.	   
35    continue
      return
      end subroutine
!*************************************************************************
       subroutine UpdateDTRELfromSrcRes(NN,DTREL) 
    

      IMPLICIT NONE
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  CALL: NONE
!  INPUT
      integer*4 NN
      real*4    DTREL(NODMAX)
!  OUTPUT


!  Function

!  Local
      integer*4 DN, i, S_NODE
      real*4    AdjSrcRel

!-------------------------------------------------------------------------
      !Later this needs to be optimized, there should be a var that stores all the demand nodes for which a res is responsible for
      !similar to the var used for links. 
      do DN=1,tnodes
        if (DMDNODE(DN))then
            do i = 1,SUPLY_PTS(DN)
               S_NODE = SUPL_NODE(I,DN)    
               if (NN==S_NODE)then
                 !now we found a demand for which the reservoir NN is responsible for
                 !Is the step_deficit going to be already fully or partiall satisfied from other flow?
                 if (SrcResRelease(NN,DN)>NonPropSTEP_DEFICIT(DN)
     &                 *SUPL_FRAC(I,DN))then
                    !if it is, remove previously added release for this dem node   
                    DTREL(NN)=DTREL(NN)-SrcResRelease(NN,DN)
                    !add the current step deficit instead
                    DTREL(NN)=DTREL(NN)+
     &                   max(NonPropSTEP_DEFICIT(DN)*SUPL_FRAC(I,DN),0.)
                    DTREL(NN)=max(DTREL(NN),0.0)
                 end if                         
               end if
            end do
        end if
      end do            
35    continue
      return
      end subroutine

!       subroutine AllocateToDeficit(DemLink_Allocation) 
!    
!! Compute this step's deficit, additional releases from sources for deficits at demand nodes
!! Seperated from FLWSIM.FOR
!      IMPLICIT NONE
!      INCLUDE 'IRAS_SYS.INC'
!      INCLUDE 'NODE.INC'
!      INCLUDE 'LINK.INC'
!!  CALL: NONE
!!  INPUT
!
!      !COMMON: TNODES,CAPN(NN), QInn(nn)
!      !        DMD_TARG(NN),DMD_T_co(NN),STEP_DEFICIT(NN)
!      !        SUPLY_PTS(NN), SUPL_FRAC(I,NN),
!!  OUTPUT
!
!      !COMMON: TOT_SUPL_R(NN)
!!  Local
!      INTEGER*4 i, j, LN, NN, S_Link 
!
!	Real*4 TOT_DEFICIT(NODMAX), RuleStoFrac,Def_Allocation(Lnkmax),LinkAllo
!	Real*4 Step_Def_diff,OrigTemp_STEP_DEFICIT,Temp_STEP_DEFICIT(nodmax)
!!-------------------------------------------------------------------------
!
!!     Initialize supplemental releases for this step
!      DO LN = 1,Links
!          Def_Allocation(LN) = 0.0
!      END DO  
!	DO NN = 1,TNODES
!!       If a demand node has a deficit, compute required demand link allocations.
!		IF (DMDNODE(NN) .and. STEP_DEFICIT(NN)>0.0)THEN 
!		 !For link priority set Temp_STEP_DEFICIT(NN)    
!		 Temp_STEP_DEFICIT(NN)=STEP_DEFICIT(NN)
!		 DO I = 1,SUPLY_Links_pts(NN)
!			 S_Link = SUPL_Link(I,NN)      !SUPL_NODE = {1,...TNodes}, not NodeID
!		     IF (S_Link.GT.0) THEN		
!			    LinkAllo=SUPL_FRAC_Link(I,NN)*Temp_STEP_DEFICIT(NN)
!		     END IF
!		     !Limit to time-step limit
!		     IF (LinkAllo+TotLinkAllo(I,NN)>(MaxAlloTS(I,NN)*DAYSPRPRD))THEN  	       					
! 	           LinkAllo=MaxAlloTS(I,NN)*DAYSPRPRD-TotLinkAllo(I,NN)
!		     END IF		     
!			 LinkAllo=max(0.0,LinkAllo)		 													 
!		     !Add to total allocation to demand link
!		     TotLink_Allocation(S_Link)=TotLink_Allocation(S_Link)+LinkAllo
!             !Set OrigTemp_STEP_DEFICIT which will be used later after Temp_STEP_DEFICIT is modified
!             OrigTemp_STEP_DEFICIT=Temp_STEP_DEFICIT(NN) 				 
! 			 !Reduce Temp_STEP_DEFICIT by the amount that the reservior released last subtimestep
! 			 Temp_STEP_DEFICIT(NN)=Temp_STEP_DEFICIT(NN)-LinkAllo !LinkAllo is not right here, we need something that predicts how much is actually available
!             !Make sure it doesnt go below 0
!             Temp_STEP_DEFICIT(NN)=max(0.0,Temp_STEP_DEFICIT(NN))
!                    
!             !Step_Def_diff is how much of last subtime step's release went to node nn
!              Step_Def_diff=OrigTemp_STEP_DEFICIT -Temp_STEP_DEFICIT(NN)
!                      
! 		     !How much of last subtime steps do we have left for lower priority nodes
! 		      PastSTSLinkAllo(S_Link)=PastSTSLinkAllo(S_Link)
!     &         				 -Step_Def_diff	 
!		 ENDDO		 
!		END IF !End if Step Deficit
!      ENDDO
!
!	end subroutine