      subroutine ChangeNonPolicyData()
!Annual changes for variables that are between the NODES LINKS and END OF NODES and END OF LINKS
!lines. This offers another method to incorporate annual changes than the policy group algorithm.
!Written by Evgenii Matrosov (C) 2011 University College London
      use vars
      implicit none
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  Local
      INTEGER*4 ln,i,j,nn
      CHARACTER*30 aVar
!------------------------------------------------------------------------
      do i = 1, nNonPolChange
       !Does this change happen this year?
       if (pNonPolChange(i)%year==sysstat(year)) then
        !For nodes
        if (pNonPolChange(i)%CompType==1) then
          do nn=1,tnodes
             if (pNonPolChange(i)%CompID==nodeid(nn)) then
                select case(trim(pNonPolChange(i)%Var_Name))
                    !For lake and reservior capacity
                    case ('CAPN')
                        capn(nn)=pNonPolChange(i)%Var_Value
                end select  
             end if                               
          end do !all nodes        
        !For links
        else if (pNonPolChange(i)%CompType==2) then
          do ln=1,links
             if (pNonPolChange(i)%CompID==linkid(ln)) then
                select case(trim(pNonPolChange(i)%Var_Name))
                    !For lake and reservior capacity
                    case ('CapL')
                        CapL(nn)=pNonPolChange(i)%Var_Value
                    case ('CapLYear')
                        CapLYear(nn)=pNonPolChange(i)%Var_Value                        
                end select  
              end if                               
          end do !all links             
        end if !Node or link
       end if !Year
      end do  !All nNonPolChange    	
      return
      end        
  
!************************************************************************
      subroutine ChangePolicyData()
!Annual changes for variables that are between the NODES LINKS and END OF NODES and END OF LINKS
!lines. This offers another method to incorporate annual changes than the policy group algorithm.
!Written by Evgenii Matrosov (C) 2011 University College London
      use vars
      use CHANGES
      implicit none
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  Local
      INTEGER*4 CompType,ipol,lc,lp,L,i,nn,iType,RealYear,iyear
      INTEGER*4 tablerowp(thres_pts_max),rowp,rowc
      INTEGER*4 tablerowc(thres_pts_max)
      INTEGER*4 GroupID,n,np,nc,ln
!------------------------------------------------------------------------
      RealYear=sysstat(year)
      iyear=sysstat(sim_year)
      !If year has a policy change

      if(YearPolChg(iyear,1)) then
        CompType = 1
        GroupID=1
        !Cycle through nodes to check which node has change
        do nn=1,tnodes
          !Cycle through policy types to check which policy type has change
          do itype=1,MaxPolicyTypes
              !So if this year, this policy type and this node has a change begin change algorithm
              if(YearNodePolicyChg(CompType,iyear,iType,nn))then
                  !The policy type
                  select case(iType)  
                  !Demand reduction tables
                  case (DemRed0)
                      !Check to make sure demand reduction array and change demand reduction array are allocated
                      if(.not. allocated(DemRedArrayC).and.
     &                      .not. allocated(DemRedArray)) goto 999
                      !Count the number of lines 
                      do ipol=1,MaxPolicies       !This should be changed to 'total policis for node and parameter type', but need to build a counter for this first                           
                        nc=0
                        do Lc=1, nDemRedC
                          if(pDemRedC(Lc)%GroupID==GroupID
     &                     .and.  pDemRedC(Lc)%Policy==ipol
     &                     .and.  pDemRedC(Lc)%NodeID==nodeid(nn)) then
                              nc = nc + 1
                              tablerowc(nc)=Lc
                          end if
                        end do
                       
                       
                        np=0
                        do Lp=1, nDemRed
                          if(pDemRed(Lp)%GroupID==GroupID
     &                     .and.  pDemRed(Lp)%Policy==ipol
     &                     .and.  pDemRed(Lp)%NodeID==nodeid(nn)) then
                              np = np + 1
                              tablerowp(np)=Lp
                          end if
                        end do
                       
                        if (np/=nc) goto 999
                       !index of first row of table for node, policy etc.   
                        Lc=tablerowc(1)
                        lp=tablerowp(1)
                       !now find corresponding row in DemRed array
                        do n=0,np-1
                           pDemRed(lp+n)%DemSourceID=
     &                     pDemRedC(Lc+n)%DemSourceID 
                           pDemRed(lp+n)%Dem_Thres_limit=
     &                     pDemRedC(Lc+n)%Dem_Thres_limit           
                           pDemRed(lp+n)%DemRedAmt=
     &                     pDemRedC(Lc+n)%DemRedAmt 
                        end do 
                      end do  !end max policies  
                  case (Allocation0)
                      !Check to make sure demand reduction array and change demand reduction array are allocated
                      if(.not. allocated(AlloArrayC).and.
     &                      .not. allocated(AlloArray)) goto 999
                      !Count the number of lines 
                      do ipol=1,MaxPolicies       !This should be changed to 'total policis for node and parameter type', but need to build a counter for this first                           
                        nc=0
                        do Lc=1, nAllocationC
                          if(pAllocationC(Lc)%GroupID==GroupID
     &                     .and.  pAllocationC(Lc)%Policy==ipol
     &                     .and.pAllocationC(Lc)%NodeID==nodeid(nn))then
                              nc = nc + 1
                              tablerowc(nc)=Lc
                          end if
                        end do
                       
                       
                        np=0
                        do Lp=1, nAllocation
                          if(pAllocation(Lp)%GroupID==GroupID
     &                     .and.  pAllocation(Lp)%Policy==ipol
     &                     .and.pAllocation(Lp)%NodeID==nodeid(nn)) then
                              np = np + 1
                              tablerowp(np)=Lp
                          end if
                        end do
                       
                        if (np/=nc) goto 999
                       !index of first row of table for node, policy etc.   
                        Lc=tablerowc(1)
                        lp=tablerowp(1)
                       !now find corresponding row in DemRed array
                        do n=0,np-1
                           pAllocation(lp+n)%NodeOutput=
     &                       pAllocationC(Lc+n)%NodeOutput 
                           pAllocation(lp+n)%LinkID=
     &                       pAllocationC(Lc+n)%LinkID           
                           pAllocation(lp+n)%LinkAllo=
     &                       pAllocationC(Lc+n)%LinkAllo 
                        end do 
                      end do  !end max policies                        
                  case (Rate0)    
                      !Check to make sure demand reduction array and change demand reduction array are allocated
                      if(.not. allocated(RatingArrayC).and.
     &                      .not. allocated(RatingArray)) goto 999
                      !Count the number of lines 
                      do ipol=1,MaxPolicies       !This should be changed to 'total policis for node and parameter type', but need to build a counter for this first                           
                        nc=0
                        do Lc=1, nRatingC
                          if(pRatingC(Lc)%GroupID==GroupID
     &                     .and.  pRatingC(Lc)%Policy==ipol
     &                     .and.pRatingC(Lc)%ID==nodeid(nn)
     &                     .and. pRatingC(Lc)%CompType==CompType)then   
                              nc = nc + 1
                              tablerowc(nc)=Lc
                          end if
                        end do
                       
                       
                        np=0
                        do Lp=1, nRating
                          if(pRating(Lp)%GroupID==GroupID
     &                     .and.  pRating(Lp)%Policy==ipol
     &                     .and. pRating(Lp)%ID==nodeid(nn)
     &                     .and. pRating(Lp)%CompType==CompType)then
                              np = np + 1
                              tablerowp(np)=Lp
                          end if
                        end do
                       
                        if (np/=nc) goto 999
                       !index of first row of table for node, policy etc.   
                        Lc=tablerowc(1)
                        lp=tablerowp(1)
                       !now find corresponding row in DemRed array
                        do n=0,np-1
                           pRating(lp+n)%ElevOrWidth=
     &                       pRatingC(Lc+n)%ElevOrWidth 
                           pRating(lp+n)%AreaOrEvapor=
     &                       pRatingC(Lc+n)%AreaOrEvapor           
                           pRating(lp+n)%VolOrFlow=
     &                       pRatingC(Lc+n)%VolOrFlow 
                           pRating(lp+n)%Seep=
     &                       pRatingC(Lc+n)%Seep 
                           pRating(lp+n)%MaxQ=
     &                       pRatingC(Lc+n)%MaxQ           
                           pRating(lp+n)%LakeQ=
     &                       pRatingC(Lc+n)%LakeQ      
                        end do 
                      end do  !end max policies                                          
                  case (Power0)
                      !Check to make sure demand reduction array and change demand reduction array are allocated
                      if(.not. allocated(PowerArrayC).and.
     &                      .not. allocated(PowerArray)) goto 999
                      !Count the number of lines 

                                             
                      do ipol=1,MaxPolicies                      
                        !Find index of Target value for policy and node in question in change Target array
                        rowc=0
                        rowp=0
                        do Lc=1, nPowerC
                          if(pPowerC(Lc)%GroupID==GroupID
     &                     .and.  pPowerC(Lc)%Policy==ipol
     &                     .and.  pPowerC(Lc)%NodeID==nodeid(nn)) then
                              rowc=Lc
                              exit
                          end if
                        end do
                       
                       !Find index of Target value for policy and node in question in Target array               
                       do Lp=1, nPower
                          if(pPower(Lp)%GroupID==GroupID
     &                     .and.  pPower(Lp)%Policy==ipol
     &                     .and.  pPower(Lp)%NodeID==nodeid(nn)) then
                            rowp=Lp
                            exit
                          end if
                       end do
                       

                       !now make the change
                        if (rowc >0 .and. rowp>0) then
                          pPower(rowp)%HPCAP=
     &                     pPowerC(rowc)%HPCAP 
                          pPower(rowp)%PLANT_FACTOR=
     &                     pPowerC(rowc)%PLANT_FACTOR           
                          pPower(rowp)%Intake_elev=
     &                     pPowerC(rowc)%Intake_elev 
                          pPower(rowp)%Turbine_elev=
     &                     pPowerC(rowc)%Turbine_elev 
                          pPump(rowp)%Outlet_elev=
     &                     pPowerC(rowc)%Outlet_elev           
                          pPower(rowp)%ECONST=
     &                     pPowerC(rowc)%ECONST  
                          pPower(rowp)%HpQMin=
     &                     pPowerC(rowc)%HpQMin       
                        end if                     
                      end do  !end max policies                   
                  case (Pump0)
                      !Check to make sure demand reduction array and change demand reduction array are allocated
                      if(.not. allocated(PumpArrayC).and.
     &                      .not. allocated(PumpArray)) goto 999
                      !Count the number of lines 

                                             
                      do ipol=1,MaxPolicies                      
                        !Find index of Target value for policy and node in question in change Target array
                        rowc=0
                        rowp=0
                        do Lc=1, nPumpC
                          if(pPumpC(Lc)%GroupID==GroupID
     &                     .and.  pPumpC(Lc)%Policy==ipol
     &                     .and.  pPumpC(Lc)%NodeID==nodeid(nn)) then
                              rowc=Lc
                              exit
                          end if
                        end do
                       
                       !Find index of Target value for policy and node in question in Target array               
                       do Lp=1, nPump
                          if(pPump(Lp)%GroupID==GroupID
     &                     .and.  pPump(Lp)%Policy==ipol
     &                     .and.  pPump(Lp)%NodeID==nodeid(nn)) then
                            rowp=Lp
                            exit
                          end if
                       end do
                       

                       !now make the change
                        if (rowc >0 .and. rowp>0) then
                          pPump(rowp)%HPCAP=
     &                     pPumpC(rowc)%HPCAP 
                          pPump(rowp)%PLANT_FACTOR=
     &                     pPumpC(rowc)%PLANT_FACTOR           
                          pPump(rowp)%Intake_elev=
     &                     pPumpC(rowc)%Intake_elev 
                          pPump(rowp)%Turbine_elev=
     &                     pPumpC(rowc)%Turbine_elev 
                          pPump(rowp)%Outlet_elev=
     &                     pPumpC(rowc)%Outlet_elev           
                          pPump(rowp)%ECONST=
     &                     pPumpC(rowc)%ECONST  
                          pPump(rowp)%HpQMin=
     &                     pPumpC(rowc)%HpQMin       
                        end if                     
                      end do  !end max policies                       
                  case (Target0)                   
                      !Check to make sure demand reduction array and change demand reduction array are allocated
                      if(.not. allocated(TargetArrayC).and.
     &                      .not. allocated(TargetArray)) goto 999
                      !Count the number of lines 

                                             
                      do ipol=1,MaxPolicies                      
                        !Find index of Target value for policy and node in question in change Target array
                        rowc=0
                        rowp=0
                        do Lc=1, nTargetC
                          if(pTargetC(Lc)%GroupID==GroupID
     &                     .and.  pTargetC(Lc)%Policy==ipol
     &                     .and.  pTargetC(Lc)%NodeID==nodeid(nn)) then
                              rowc=Lc
                              exit
                          end if
                        end do
                       
                       !Find index of Target value for policy and node in question in Target array               
                       do Lp=1, nTarget
                          if(pTarget(Lp)%GroupID==GroupID
     &                     .and.  pTarget(Lp)%Policy==ipol
     &                     .and.  pTarget(Lp)%NodeID==nodeid(nn)) then
                            rowp=Lp
                            exit
                          end if
                       end do
                       

                       !now make the change
                        if (rowc >0 .and. rowp>0) then
                          pTarget(rowp)%targ=
     &                     pTargetC(rowc)%targ 
                          pTarget(rowp)%t_co=
     &                     pTargetC(rowc)%t_co           
                          pTarget(rowp)%EnvFlow=
     &                     pTargetC(rowc)%EnvFlow 
                          pTarget(rowp)%RefilTrig=
     &                     pTargetC(rowc)%RefilTrig 
                          pTarget(rowp)%DemInc=
     &                     pTargetC(rowc)%DemInc           
                          pTarget(rowp)%srcpriority=
     &                     pTargetC(rowc)%srcpriority  
                        end if                     
                      end do  !end max policies                        
                  case (Source0)
                      !Check to make sure demand reduction array and change demand reduction array are allocated
                      if(.not. allocated(SourceArrayC).and.
     &                      .not. allocated(SourceArray)) goto 999
                      !Count the number of lines 
                      do ipol=1,MaxPolicies       !This should be changed to 'total policis for node and parameter type', but need to build a counter for this first                           
                        nc=0
                        do Lc=1, nSourceC
                          if(pSourceC(Lc)%GroupID==GroupID
     &                     .and.  pSourceC(Lc)%Policy==ipol
     &                     .and.pSourceC(Lc)%NodeID==nodeid(nn)
     &                     .and. pSourceC(Lc)%CompType==CompType)then   
                              nc = nc + 1
                              tablerowc(nc)=Lc
                          end if
                        end do
                       
                       
                        np=0
                        do Lp=1, nSource
                          if(pRating(Lp)%GroupID==GroupID
     &                     .and.  pSource(Lp)%Policy==ipol
     &                     .and. pSource(Lp)%NodeID==nodeid(nn)
     &                     .and. pSource(Lp)%CompType==CompType)then
                              np = np + 1
                              tablerowp(np)=Lp
                          end if
                        end do
                       
                        if (np/=nc) goto 999
                       !index of first row of table for node, policy etc.   
                        Lc=tablerowc(1)
                        lp=tablerowp(1)
                       !now find corresponding row in source array
                        do n=0,np-1
                           pSource(lp+n)%Supl_Frac=
     &                       pSourceC(Lc+n)%Supl_Frac 
                           pSource(lp+n)%MaxOutTS=
     &                       pSourceC(Lc+n)%MaxOutTS           
                           pSource(lp+n)%MaxOutYear=
     &                       pSourceC(Lc+n)%MaxOutYear 
                           pSource(lp+n)%LnkProp= !Evgenii added LNKPROP120126 for link propogation links
     &                       pSourceC(Lc+n)%LnkProp           
                        end do 
                      end do  !end max policies                        
                  case (Rule0)
                      !Check to make sure demand reduction array and change demand reduction array are allocated
                      if(.not. allocated(RuleArrayC).and.
     &                      .not. allocated(RuleArray)) goto 999
                      !Count the number of lines 
                      do ipol=1,MaxPolicies       !This should be changed to 'total policis for node and parameter type', but need to build a counter for this first                           
                        nc=0
                        do Lc=1, nRuleC
                          if(pRuleC(Lc)%GroupID==GroupID
     &                     .and.  pRuleC(Lc)%Policy==ipol
     &                     .and.pRuleC(Lc)%NodeID==nodeid(nn))then
                              nc = nc + 1
                              tablerowc(nc)=Lc
                          end if
                        end do
                       
                       
                        np=0
                        do Lp=1, nRule
                          if(pSource(Lp)%GroupID==GroupID
     &                     .and.  pRule(Lp)%Policy==ipol
     &                     .and. pRule(Lp)%NodeID==nodeid(nn))then
                              np = np + 1
                              tablerowp(np)=Lp
                          end if
                        end do
                       
                        if (np/=nc) goto 999
                       !index of first row of table for node, policy etc.   
                        Lc=tablerowc(1)
                        lp=tablerowp(1)
                       !now find corresponding row in DemRed array
                        do n=0,np-1
                           pRule(lp+n)%res_rule(1)=
     &                       pRuleC(Lc+n)%res_rule(1) 
                           pRule(lp+n)%res_rule(2)=
     &                       pRuleC(Lc+n)%res_rule(2)           
                           pRule(lp+n)%res_rule(3)=
     &                       pRuleC(Lc+n)%res_rule(3)   
                           pRule(lp+n)%res_rule(4)=
     &                       pRuleC(Lc+n)%res_rule(4) 
                           pRule(lp+n)%res_rule(5)=
     &                       pRuleC(Lc+n)%res_rule(5)           
                           pRule(lp+n)%res_rule(6)=
     &                       pRuleC(Lc+n)%res_rule(6)      
                           pRule(lp+n)%res_rule(7)=
     &                       pRuleC(Lc+n)%res_rule(7) 
                           pRule(lp+n)%res_rule(8)=
     &                       pRuleC(Lc+n)%res_rule(8)                    
                        end do 
                      end do  !end max policies                        
                          
                  case (Balance0)
                      !Check to make sure demand reduction array and change demand reduction array are allocated
                      if(.not. allocated(BalanceArrayC).and.
     &                      .not. allocated(BalanceArray)) goto 999
                      !Count the number of lines 
                      do ipol=1,MaxPolicies       !This should be changed to 'total policis for node and parameter type', but need to build a counter for this first                           
                        nc=0
                        do Lc=1, nBalanceC
                          if(pBalanceC(Lc)%GroupID==GroupID
     &                     .and.  pBalanceC(Lc)%Policy==ipol
     &                     .and.pBalanceC(Lc)%RuleID==nodeid(nn))then  
                              nc = nc + 1
                              tablerowc(nc)=Lc
                          end if
                        end do
                       
                       
                        np=0
                        do Lp=1, nBalance
                          if(pBalance(Lp)%GroupID==GroupID
     &                     .and.  pBalance(Lp)%Policy==ipol
     &                     .and. pBalance(Lp)%RuleID==nodeid(nn))then
                              np = np + 1
                              tablerowp(np)=Lp
                          end if
                        end do
                       
                        if (np/=nc) goto 999
                       !index of first row of table for node, policy etc.   
                        Lc=tablerowc(1)
                        lp=tablerowp(1)
                       !now find corresponding row in DemRed array
                        do n=0,np-1
                           pBalance(lp+n)%BalMeth=
     &                       pBalanceC(Lc+n)%BalMeth 
                           pBalance(lp+n)%GroupVol=
     &                       pBalanceC(Lc+n)%GroupVol           
                           pBalance(lp+n)%charBalance=
     &                       pBalanceC(Lc+n)%charBalance      
                        end do 
                      end do  !end max policies                              
                  case (Evaporation0)
                      !Check to make sure demand reduction array and change demand reduction array are allocated
                      if(.not. allocated(EvapArrayC).and.
     &                      .not. allocated(EvapArray)) goto 999
                      !Count the number of lines 

                                             
                      do ipol=1,MaxPolicies                      
                        !Find index of Target value for policy and node in question in change Target array
                        rowc=0
                        rowp=0
                        do Lc=1, nEvaporC
                          if(pEvaporC(Lc)%GroupID==GroupID
     &                     .and.  pEvaporC(Lc)%Policy==ipol
     &                     .and.  pEvaporC(Lc)%ID==nodeid(nn)
     &                     .and. pEvaporC(Lc)%CompType==CompType)then
                              rowc=Lc
                              exit
                          end if
                        end do
                       
                       !Find index of Target value for policy and node in question in Target array               
                       do Lp=1, nEvapor
                          if(pEvapor(Lp)%GroupID==GroupID
     &                     .and.  pEvapor(Lp)%Policy==ipol
     &                     .and.  pEvapor(Lp)%ID==nodeid(nn) 
     &                     .and. pEvapor(Lp)%CompType==CompType)then   
                            rowp=Lp
                            exit
                          end if
                       end do
                       

                       !now make the change
                        if (rowc >0 .and. rowp>0) then
                          pEvapor(rowp)%Evaporation=
     &                     pEvaporC(rowc)%Evaporation              
                        end if                     
                      end do  !end max policies                        
                  case (Performance0)
                      !Check to make sure demand reduction array and change demand reduction array are allocated
                      if(.not. allocated(PerfArrayC).and.
     &                      .not. allocated(PerfArray)) goto 999
                      !Count the number of lines 
                      do ipol=1,MaxPolicies       !This should be changed to 'total policis for node and parameter type', but need to build a counter for this first                           
                        nc=0
                        do Lc=1, nPerformanceC
                          if(pPerformanceC(Lc)%GroupID==GroupID
     &                     .and.  pPerformanceC(Lc)%Policy==ipol
     &                    .and.pPerformanceC(Lc)%NODEID==nodeid(nn))then  
                              nc = nc + 1
                              tablerowc(nc)=Lc
                          end if
                        end do
                       
                       
                        np=0
                        do Lp=1, nPerformance
                          if(pPerformance(Lp)%GroupID==GroupID
     &                     .and.  pPerformance(Lp)%Policy==ipol
     &                     .and.pPerformance(Lp)%NODEID==nodeid(nn))then
                              np = np + 1
                              tablerowp(np)=Lp
                          end if
                        end do
                       
                        if (np/=nc) goto 999
                       !index of first row of table for node, policy etc.   
                        Lc=tablerowc(1)
                        lp=tablerowp(1)
                       !now find corresponding row in DemRed array
                        do n=0,np-1
                           pPerformance(lp+n)%thres_limit=
     &                       pPerformanceC(Lc+n)%thres_limit     
                        end do               
                      end do  !end max policies
                  end select  
               end if !YearNodePolicyChg(iyear,iType,nn)
              end do   
              
          end do !nn =1,tnodes
      end if !if YearPolChg
      !links
      if(YearPolChg(iyear,2)) then
        CompType = 2
        GroupID=1
        !Cycle through nodes to check which node has change
        do ln=1,links
          !Cycle through policy types to check which policy type has change
          do itype=1,MaxPolicyTypes
              !So if this year, this policy type and this node has a change begin change algorithm
              if(YearNodePolicyChg(CompType,iyear,iType,ln))then
                  !The policy type
                  select case(iType)
                  case (Routing0)
                  !Routing doesnt need midrun changes    
                  case (Transfer0)
                      !Check to make sure demand reduction array and change demand reduction array are allocated
                      if(.not. allocated(TransferArrayC).and.
     &                      .not. allocated(TransferArray)) goto 999
                      !Count the number of lines                       
                      do ipol=1,MaxPolicies                      
                        !Find index of Target value for policy and node in question in change Target array
                        rowc=0
                        rowp=0
                        do Lc=1, nCostC
                          if(pTransferC(Lc)%GroupID==GroupID
     &                     .and.  pTransferC(Lc)%Policy==ipol
     &                     .and.  pTransferC(Lc)%linkID==linkid(ln))then
                              rowc=Lc
                              exit
                          end if
                        end do
                       
                       !Find index of Target value for policy and node in question in Target array               
                       do Lp=1, nCost
                          if(pTransfer(Lp)%GroupID==GroupID
     &                     .and.  pTransfer(Lp)%Policy==ipol
     &                     .and.  pTransfer(Lp)%linkID==linkid(ln))then   
                            rowp=Lp
                            exit
                          end if
                       end do
                       

                       !now make the change
                        if (rowc >0 .and. rowp>0) then
                          pTransfer(rowp)%GWFromVol=
     &                     pTransferC(rowc)%GWFromVol              
                          pTransfer(rowp)%GWToVol=
     &                     pTransferC(rowc)%GWToVol
                          pTransfer(rowp)%GWFlowFromTo=
     &                     pTransferC(rowc)%GWFlowFromTo                 
                        end if                     
                      end do  !end max policies        
				  
                  case (Rate0)    
                     !Check to make sure demand reduction array and change demand reduction array are allocated
                      if(.not. allocated(RatingArrayC).and.
     &                      .not. allocated(RatingArray)) goto 999
                      !Count the number of lines 
                      do ipol=1,MaxPolicies       !This should be changed to 'total policis for node and parameter type', but need to build a counter for this first                           
                        nc=0
                        do Lc=1, nRatingC
                          if(pRatingC(Lc)%GroupID==GroupID
     &                     .and.  pRatingC(Lc)%Policy==ipol
     &                     .and.pRatingC(Lc)%ID==linkid(ln)
     &                     .and. pRatingC(Lc)%CompType==CompType)then
                              nc = nc + 1
                              tablerowc(nc)=Lc
                          end if
                        end do
                       
                       
                        np=0
                        do Lp=1, nRating
                          if(pRating(Lp)%GroupID==GroupID
     &                     .and.  pRating(Lp)%Policy==ipol
     &                     .and. pRating(Lp)%ID==linkid(ln)
     &                     .and. pRating(Lp)%CompType==CompType)then
                              np = np + 1
                              tablerowp(np)=Lp
                          end if
                        end do
                       
                        if (np/=nc) goto 999
                       !index of first row of table for node, policy etc.   
                        Lc=tablerowc(1)
                        lp=tablerowp(1)
                       !now find corresponding row in DemRed array
                        do n=0,np-1
                           pRating(lp+n)%ElevOrWidth=
     &                       pRatingC(Lc+n)%ElevOrWidth 
                           pRating(lp+n)%AreaOrEvapor=
     &                       pRatingC(Lc+n)%AreaOrEvapor           
                           pRating(lp+n)%VolOrFlow=
     &                       pRatingC(Lc+n)%VolOrFlow      
                        end do 
                      end do  !end max policies                                          
                  
                  case (Cost0)
                      !Check to make sure demand reduction array and change demand reduction array are allocated
                      if(.not. allocated(CostArrayC).and.
     &                      .not. allocated(CostArray)) goto 999
                      !Count the number of lines                       
                      do ipol=1,MaxPolicies                      
                        !Find index of Target value for policy and node in question in change Target array
                        rowc=0
                        rowp=0
                        do Lc=1, nCostC
                          if(pCostC(Lc)%GroupID==GroupID
     &                     .and.  pCostC(Lc)%Policy==ipol
     &                     .and.  pCostC(Lc)%linkID==linkid(ln))then
                              rowc=Lc
                              exit
                          end if
                        end do
                       
                       !Find index of Target value for policy and node in question in Target array               
                       do Lp=1, nCost
                          if(pCost(Lp)%GroupID==GroupID
     &                     .and.  pCost(Lp)%Policy==ipol
     &                     .and.  pCost(Lp)%linkID==linkid(ln))then   
                            rowp=Lp
                            exit
                          end if
                       end do
                       

                       !now make the change
                        if (rowc >0 .and. rowp>0) then
                          pCost(rowp)%FlowCost=
     &                     pCostC(rowc)%FlowCost              
                          pCost(rowp)%FlowEng=
     &                     pCostC(rowc)%FlowEng
                          pCost(rowp)%AnnCostInc=
     &                     pCostC(rowc)%AnnCostInc                 
                        end if                     
                      end do  !end max policies                        
                  
                  case (5)
                     !No policy defined for links 
                  case (6)
                     !No policy defined for links
                  case (7)
                     !No policy defined for links 
                  case (8)
                     !No policy defined for links 
                  case (9)
                     !No policy defined for links 
                  case (Evaporation0)
                      !Check to make sure demand reduction array and change demand reduction array are allocated
                      if(.not. allocated(EvapArrayC).and.
     &                      .not. allocated(EvapArray)) goto 999
                      !Count the number of lines 

                                             
                      do ipol=1,MaxPolicies                      
                        !Find index of Target value for policy and node in question in change Target array
                        rowc=0
                        rowp=0
                        do Lc=1, nEvaporC
                          if(pEvaporC(Lc)%GroupID==GroupID
     &                     .and.  pEvaporC(Lc)%Policy==ipol
     &                     .and.  pEvaporC(Lc)%ID==linkid(ln)
     &                     .and. pEvaporC(Lc)%CompType==CompType)then 
                              rowc=Lc
                              exit
                          end if
                        end do
                       
                       !Find index of Target value for policy and node in question in Target array               
                       do Lp=1, nEvapor
                          if(pEvapor(Lp)%GroupID==GroupID
     &                     .and.  pEvapor(Lp)%Policy==ipol
     &                     .and.  pEvapor(Lp)%ID==linkid(ln)
     &                     .and.  pEvapor(Lp)%CompType==CompType)then  
                            rowp=Lp
                            exit
                          end if
                       end do
                       

                       !now make the change
                        if (rowc >0 .and. rowp>0) then
                          pEvapor(rowp)%Evaporation=
     &                     pEvaporC(rowc)%Evaporation           
                          pEvapor(rowp)%LossMethod=
     &                     pEvaporC(rowc)%LossMethod       
                        end if                     
                      end do  !end max policies                         
                  case (11)
                  end select  
              end if  !YearNodePolicyChg(CompType,iyear,iType,ln)
            end do 
          end do !links
      end if !YearPolChg(iyear,2)
      
999   continue
      return
      end       