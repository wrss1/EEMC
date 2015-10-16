	Module CHANGES
	!Declares global allocatable variables
	IMPLICIT NONE
	SAVE
        TYPE ELEMENT_TARGET_CHANGE
		  SEQUENCE 
            INTEGER*4   :: year, GroupID, Policy, CompType, NodeID
            INTEGER*4   :: srcpriority
            REAL*4      :: Targ, T_CO, EnvFlow, RefilTrig, DemInc
        END TYPE ELEMENT_TARGET_CHANGE
        
        
        
        TYPE ELEMENT_EVAPORATION_CHANGE
            SEQUENCE 
		  INTEGER*4   :: year,GroupID, Policy, CompType, ID        !for both Node and Link
            REAL*4      :: Evaporation
            INTEGER*4   :: LossMethod                           !only for Link
        END TYPE ELEMENT_EVAPORATION_CHANGE        
        !structure defintion of link routing
        TYPE ELEMENT_LROUTING_CHANGE
            SEQUENCE 
		    INTEGER*4   :: year,GroupID, Policy, CompType, LinkID
		    INTEGER*4   :: iMethod
            REAL*4      :: L_NRTRES,L_a, L_b,L_c
        END TYPE ELEMENT_LROUTING_CHANGE
        !structure defintion of Allocation
        TYPE ELEMENT_ALLOCATION_CHANGE
            SEQUENCE 
		    INTEGER*4   :: year, GroupID, Policy, CompType, NodeID
            REAL*4      :: NodeOutput
            INTEGER*4   :: LinkID
            REAL*4      :: LinkAllo
        END TYPE ELEMENT_ALLOCATION_CHANGE
        !structure defintion of Source
        TYPE ELEMENT_SOURCE_CHANGE
            SEQUENCE 
		    INTEGER*4   :: year,GroupID, Policy, CompType, NodeID, Supl_Node, 
     &		  Source_Type, SrcPriority
            REAL*4      :: Supl_Frac,MaxOutTS, MaxOutYear
            INTEGER*4   :: LnkProp !Evgenii added LNKPROP120126 for link propogation links            
        END TYPE ELEMENT_SOURCE_CHANGE
        !structure defintion of rating
        !CompType=1 for node; CompType=2 for link
        !It's better that Node and Link use different key words
        !(e.g. "NodeRating" and "LinkRating" because the data is different).  
        TYPE ELEMENT_RATING_CHANGE
            SEQUENCE 
		    INTEGER*4   :: year,GroupID, Policy, CompType, ID        !for both Node and Link
            REAL*4      :: ElevOrWidth, AreaOrEvapor, VolOrFlow !for both Node and Link
            REAL*4      :: Seep, MaxQ, LakeQ                    !only for Node            
        END TYPE ELEMENT_RATING_CHANGE
        !structure defintion of CrossSection
        TYPE ELEMENT_CROSS_CHANGE
            SEQUENCE 
		    INTEGER*4   :: year, GroupID, Policy, CompType, LinkID
            REAL*4      :: BaseWidth, ChannelDepth, LSlope, RSlope
            REAL*4      :: UpLSlope, UpRSlope
        END TYPE ELEMENT_CROSS_CHANGE
        !structure defintion of Performance
        TYPE ELEMENT_PERFORMANCE_CHANGE
            SEQUENCE 
		    INTEGER*4   :: year, GroupID, Policy, NodeID
            REAL*4      :: thres_limit
        END TYPE ELEMENT_PERFORMANCE_CHANGE
        !structure defintion of Cost
        TYPE ELEMENT_COST_CHANGE
            SEQUENCE 
		    INTEGER*4   :: year, GroupID, Policy, LinkID
            REAL*4      :: FlowCost, FlowEng,AnnCostInc
        END TYPE ELEMENT_COST_CHANGE
        !structure defintion of Power and Pump
        !Both Power and Pump have the same data structure
        TYPE ELEMENT_POWER_PUMP_CHANGE
            SEQUENCE 
		    INTEGER*4 :: year, GroupID, Policy, NodeID 
            REAL*4    :: HPCAP,PLANT_FACTOR,Intake_elev,Turbine_elev
            REAL*4    :: Outlet_elev,ECONST,HpQMin            
        END TYPE ELEMENT_POWER_PUMP_CHANGE
        !structure defintion of Rule
        TYPE ELEMENT_RULE_CHANGE
            SEQUENCE 
		    INTEGER*4 :: year, GroupID, Policy, NodeID
            REAL*4    :: res_rule(8)
        END TYPE ELEMENT_RULE_CHANGE

	  TYPE ELEMENT_DEMRED_CHANGE	   !Added by Evgenii 100719
            SEQUENCE 
		    INTEGER*4   :: year,GroupID, Policy, NodeID
            INTEGER*4   :: DemSourceID
		    REAL*4      :: Dem_Thres_limit,DemRedAmt
        END TYPE ELEMENT_DEMRED_CHANGE

	  TYPE ELEMENT_BALANCE_CHANGE	   !Evgenii 110504
          SEQUENCE 
		  INTEGER*4   :: year, GroupID, Policy, RuleID, BalMeth
		  character*256 :: charBalance
 		  REAL*4 :: GroupVol
        END TYPE ELEMENT_BALANCE_CHANGE 
        
        TYPE ELEMENT_GW_CHANGE        !Added by Anthony 230112
		  SEQUENCE 
            INTEGER*4   :: year, GroupID, Policy, CompType, LinkID
            INTEGER*4   :: GWMethod
            REAL*4      :: GWK, GWElev, GWLength, GWWidth
         END TYPE ELEMENT_GW_CHANGE
         
         !   Defining structure of transfer table
         TYPE ELEMENT_TRANSFER_CHANGE        !Added by Anthony 230112
		  SEQUENCE 
            INTEGER*4   :: year, GroupID, Policy, LinkID
            REAL*4      :: GWFromVol, GWToVol, GWFlowFromTo
        END TYPE ELEMENT_TRANSFER_CHANGE

        INTEGER*4::nTargetC
        TYPE (ELEMENT_TARGET_CHANGE), POINTER::pTargetC(:)
        INTEGER*4::nEvaporC
        TYPE (ELEMENT_EVAPORATION_CHANGE), POINTER::pEvaporC(:)
        INTEGER*4::nRoutingC
        TYPE (ELEMENT_LROUTING_CHANGE), POINTER::pRoutingC(:)
        INTEGER*4::nAllocationC
        TYPE (ELEMENT_ALLOCATION_CHANGE), POINTER::pAllocationC(:)
        INTEGER*4::nSourceC
        TYPE (ELEMENT_SOURCE_CHANGE), POINTER::pSourceC(:)
        INTEGER*4::nRatingC
        TYPE (ELEMENT_RATING_CHANGE), POINTER::pRatingC(:)
        INTEGER*4::nCrossC
        TYPE (ELEMENT_CROSS_CHANGE), POINTER::pCrossC(:)
        INTEGER*4::nPerformanceC
        TYPE (ELEMENT_PERFORMANCE_CHANGE), POINTER::pPerformanceC(:)
        INTEGER*4::nCostC
        TYPE (ELEMENT_COST_CHANGE), POINTER::pCostC(:)
        INTEGER*4::nPowerC
        TYPE (ELEMENT_POWER_PUMP_CHANGE), POINTER::pPowerC(:)
        INTEGER*4::nPumpC
        TYPE (ELEMENT_POWER_PUMP_CHANGE), POINTER::pPumpC(:)
        INTEGER*4::nRuleC
        TYPE (ELEMENT_RULE_CHANGE), POINTER::pRuleC(:)
	    INTEGER*4::nDemRedC
	    TYPE (ELEMENT_DEMRED_CHANGE), POINTER::pDemRedC(:)
	    INTEGER*4::nBalanceC
	    TYPE (ELEMENT_BALANCE_CHANGE), POINTER::pBALANCEC(:)	  !Evgenii 110504  
        INTEGER*4::nGWC
        TYPE (ELEMENT_GW_CHANGE), POINTER::pGWC(:) !Added by Anthony 230112
        INTEGER*4::nTransferC
        TYPE (ELEMENT_TRANSFER_CHANGE), POINTER::pTransferC(:) !Added by Anthony 230112
              
      TYPE (ELEMENT_TARGET_CHANGE),DIMENSION(:)
     &      ,ALLOCATABLE, TARGET::TargetArrayC
	  TYPE (ELEMENT_EVAPORATION_CHANGE),DIMENSION(:)
     &	, ALLOCATABLE, TARGET::EvapArrayC
	  TYPE (ELEMENT_LROUTING_CHANGE),DIMENSION(:)
     &	, ALLOCATABLE, TARGET::RoutingArrayC
	  TYPE (ELEMENT_ALLOCATION_CHANGE),DIMENSION(:)
     &	, ALLOCATABLE, TARGET::AlloArrayC
	  TYPE (ELEMENT_SOURCE_CHANGE),DIMENSION(:)
     &	, ALLOCATABLE, TARGET::SourceArrayC
	  TYPE (ELEMENT_RATING_CHANGE),DIMENSION(:)
     &	, ALLOCATABLE, TARGET::RatingArrayC
  	  TYPE (ELEMENT_CROSS_CHANGE),DIMENSION(:)
     &	, ALLOCATABLE, TARGET::CrossArrayC
      TYPE (ELEMENT_PERFORMANCE_CHANGE),DIMENSION(:)
     &	, ALLOCATABLE, TARGET::PerfArrayC
	  TYPE (ELEMENT_COST_CHANGE),DIMENSION(:)
     &	, ALLOCATABLE, TARGET::CostArrayC
	  TYPE (ELEMENT_POWER_PUMP_CHANGE),DIMENSION(:)
     &	, ALLOCATABLE, TARGET::PowerArrayC
	  TYPE (ELEMENT_POWER_PUMP_CHANGE),DIMENSION(:)
     &	, ALLOCATABLE, TARGET::PumpArrayC	
	  TYPE (ELEMENT_RULE_CHANGE),DIMENSION(:)
     &	, ALLOCATABLE, TARGET::RuleArrayC
	  TYPE (ELEMENT_DEMRED_CHANGE),DIMENSION(:)
     &	, ALLOCATABLE,TARGET::DemRedArrayC		!Evgenii 100719
 	  TYPE (ELEMENT_BALANCE_CHANGE),DIMENSION(:)
     & 	, ALLOCATABLE,TARGET::BalanceArrayC	!Evgenii 110504
      TYPE (ELEMENT_GW_CHANGE),DIMENSION(:)  !Added by Anthony 230112
     &      ,ALLOCATABLE, TARGET::GWArrayC     !Added by Anthony 230112
      TYPE (ELEMENT_TRANSFER_CHANGE),DIMENSION(:)     !Added by Anthony 230112
     &      ,ALLOCATABLE, TARGET::TransferArrayC     !Added by Anthony 230112
      
      END MODULE CHANGES
              