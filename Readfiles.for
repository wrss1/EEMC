	subroutine readFlowFile(success)
	USE vars
	IMPLICIT NONE	
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
	!Local
	INTEGER*4 iFlowFile
	INTEGER*4 ierr,j
	logical*1 success
	character*10000 aLine
	INTEGER*4 LineCounter, nRead
	!Call
	LOGICAL CountFileLines
!------------------------------------------------------------
	!Read flow data
	success=.false.
	iFlowFile = 242
      OPEN(UNIT=iFlowFile, FILE=TRIM(flwdataFileName), STATUS='old',
     &	 FORM='formatted', ERR= 999,ACTION='READ')
    		
      !Count the total number of .flw file
      LineCounter = 0
      IF(.NOT. CountFileLines(iFlowFile, LineCounter)) GOTO 999 
	
	!Allocate flowdata variable to number of lines of code, number of gauges nodes
	if(.not.allocated(flowdata) )
     &  Allocate (flowdata(LineCounter,sysstat(NGAUGE1)))
	!Read inflows into memory from iras.flw file	
	nRead = 0       !no record read
	rewind (unit=iFlowFile)
	do WHILE (.TRUE.)
  	  READ(UNIT=iFlowFile,FMT='(A)',iostat=ierr) aLine
	  IF(ierr ==0) THEN 
		nRead=nRead+1
		READ(aLine,*)(FlowData(nRead,j), j=1, SysStat(NGAUGE1))	
        ELSE IF (ierr<0) THEN 
            EXIT
        else if (ierr>0) then
		continue	
	  END IF   
	end do
	sysstat(nRec)=nread

888   success = .true.
999   CLOSE(UNIT=iFlowFile)
      if (.not. success) then
        WRITE(*,*)'Flow data reading failed'
      end if
      return
	end subroutine

!**********************************************************
	
	Subroutine ReadFlwFactors(success)
!	Evgenii created this subroutine on 091104, it reads iras.flw
!	gauge flow factors for run into memory
      use vars
	implicit none
	INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT

      
!  OUTPUT
      LOGICAL*1 success
      !COMMON: flowdata
!  Local
      INTEGER*4 iFactfile,j,i,nr
	CHARACTER*10000 aLine
      CHARACTER*50 aVar
	character(len=100) GageFactName
!--------------------------------------------------------------
	success = .false.
	!Open flow factor files for gauge nodes with flow factors
!	do i=1,SysStat(NGAUGE1)
!		if (GageFact(i)) then
!		  GageFactName(i)=trim(GageName(i))//'.fac'
!		  outunit(i)=100+i
!		  OPEN(UNIT = outunit(i),FILE=GageFactName(i),STATUS='old',
!    &		  FORM='formatted')
!		endif
!	enddo
      if (.not. allocated(flowfactor))
     &	Allocate(flowfactor(sysstat(NGAUGE1),sysstat(nRuns),12))
	
	!EVGENII - Must add check to make sure number of flowfactor records match number of runs
	!Initialize flow factos to 1
!	do j=1,gagmax
!		do k=1,12
!			flowfactor(j,k)=1.0
!		end do
!	end do
	iFactfile=100
	!Read flow factors into memory
  	do i=1,SysStat(NGAUGE1)
		!If Gauge node has flowfactors
		if (GageFact(i)) then
			!Name file
			GageFactName=trim(GageName(i))//'.fac'
		    !Read flow factors for each run from file into memory
			OPEN(UNIT = iFactfile,FILE=trim(GageFactName),
     &		      STATUS='old',FORM='formatted')

			do nr=1,sysstat(nRuns)		  			
			   READ(UNIT=iFactfile,FMT='(A)',ERR=999, END=999)aLine
     			   read(aLine,*)avar,(flowfactor(i,nr,j),j=1,12)
			end do 			
	    	CLOSE(UNIT=iFactfile)
		!If Gauge node doesnt have flowfactors
		else 
			do nr=1,sysstat(nRuns)		  
				do j=1,12; flowfactor(i,nr,j)=1 ;end do
			end do	
		  endif
	enddo


      success = .true.
 999  if (.not. success) then
        WRITE(*,*)'Flow factor reading failed ',GageFactName
      end if
      return
      end subroutine

!**********************************************************
	
	Subroutine ReadChangesFile(success)
!	Evgenii created this subroutine on 091104, it reads iras.flw
!	gauge flow factors for run into memory
      use vars
      use CHANGES
	implicit none
	INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
!  OUTPUT
      LOGICAL*1 success
!  Local
      INTEGER*4 iChangefile,i,ierr,iPos,y,p,nn,tyears,j,begyear,realyear
      INTEGER*4 ln
	CHARACTER*256 aLine
	CHARACTER*30 aVar
!--------------------------------------------------------------
	success = .false.
	iChangefile = 262
      OPEN(UNIT=iChangefile, FILE=TRIM('input/iras.chg'), STATUS='old',
     &	 FORM='formatted', ERR= 999,ACTION='READ')
	
      NULLIFY(pNonPolChange)
      nNonPolChange = 0
!     Determine the number of non-policy varialbe changes in the network
      do WHILE (.TRUE.)
        READ(UNIT=iChangefile,FMT='(A)',ERR=999) aLine
        if (TRIM(aLine) == 'NON-POLICY CHANGES') EXIT
      end do
      i = 0
      do WHILE (.TRUE.)
        READ(UNIT=iChangefile,FMT='(A)',ERR=999) aLine
        if (Index(aLine,'END OF NON-POLICY CHANGES')>=1)exit
        if (aLine(1:1).ne."!") i = i + 1
      end do
      nNonPolChange=i
     
      if(.not. allocated(NonPolChangeArray))
     &  ALLOCATE (NonPolChangeArray(nNonPolChange))         
	rewind (iChangefile, ERR=999)
!     find the 'NON-POLICY CHANGES' line
      do WHILE (.TRUE.)
        READ(UNIT=iChangefile,FMT='(A)',ERR=999) aLine
!	Bring curser down to NON-POLICY CHANGES line
	  if (TRIM(aLine) .eq. 'NON-POLICY CHANGES') EXIT 
      end do
!     Store non-policy change data
      DO i = 1, nNonPolChange
	  READ(UNIT=iChangefile,FMT='(A)',iostat=ierr) aLine 
	  IF(ierr ==0) THEN 
               READ(aLine,*)NonPolChangeArray(i)%CompID,  
     &            NonPolChangeArray(i)%CompType,
     &            NonPolChangeArray(i)%Year,
     &            NonPolChangeArray(i)%Var_Name,
     &            NonPolChangeArray(i)%Var_Value  
        ELSE IF(ierr<0) then 
            EXIT           
        END IF
      END DO

      if (ALLOCATED(NonPolChangeArray)) pNonPolChange => 
     &    NonPolChangeArray(:)      


      !Now read in policy variable changes
      nTargetC  = 0
      NULLIFY(pTargetC)
      nEvaporC = 0
      NULLIFY(pEvaporC)
      nRoutingC = 0
      NULLIFY(pRoutingC)
      nAllocationC = 0
      NULLIFY(pAllocationC)
      nSourceC= 0
      NULLIFY(pSourceC)
      nRatingC = 0
      NULLIFY(pRatingC)
      nCrossC = 0
      NULLIFY(pCrossC)
      nPerformanceC = 0
      NULLIFY(pPerformanceC)
      nCostC = 0
      NULLIFY(pCostC)
      nPowerC = 0
      NULLIFY(pPowerC)
      nPumpC = 0
      NULLIFY(pPumpC)
      nRuleC = 0
      NULLIFY(pRuleC)
      nDemRedC=0			 
      NULLIFY(pDemRedC)		
	  nBalanceC=0			
      NULLIFY(pBalanceC)
      nGWC  = 0					! Added by Anthony 230112
      NULLIFY(pGWC)				! Added by Anthony 230112
      nTransferC  = 0			! Added by Anthony 230112
      NULLIFY(pTransferC)		! Added by Anthony 230112

      !Count the number of each type of change
      DO WHILE (.TRUE.)
        READ(UNIT=iChangefile,FMT='(A)',iostat=ierr) aLine !Evgenii took out ERR=999
        IF(ierr ==0) THEN 
            IF(INDEX(aLine,headTarget)>0) THEN
                nTargetC = nTargetC+1
            ELSE IF(INDEX(aLine,headEvaporation)==1) THEN   
                nEvaporC = nEvaporC+1
            ELSE IF(INDEX(aLine,headRouting)==1) THEN   
                nRoutingC = nRoutingC+1
            ELSE IF(INDEX(aLine,headAllocation)==1) THEN   
                nAllocationC = nAllocationC+1
            ELSE IF(INDEX(aLine,headSource)==1) THEN   
                nSourceC = nSourceC+1
            ELSE IF(INDEX(aLine,headRating)==1) THEN   
                nRatingC = nRatingC+1
            ELSE IF(INDEX(aLine,headCross)==1) THEN   
                nCrossC = nCrossC+1
            ELSE IF(INDEX(aLine,headPerformance)==1) THEN   
                nPerformanceC = nPerformanceC+1
            ELSE IF(INDEX(aLine,headCost)==1) THEN   
                nCostC = nCostC+1
            ELSE IF(INDEX(aLine,headPower)==1) THEN   
                nPowerC = nPowerC+1
            ELSE IF(INDEX(aLine,headPump)==1) THEN   
                nPumpC = nPumpC+1
            ELSE IF(INDEX(aLine,headRule)==1) THEN   
                nRuleC = nRuleC+1
		    ELSE IF(INDEX(aLine,headDemRed)==1) THEN  !Added by Evgenii 101907  
                nDemRedC= nDemRedC+1
		    ELSE IF(INDEX(aLine,headBalance)==1) THEN  !Added by Evgenii 101907  
                nBalanceC = nBalanceC+1
            ELSE IF(INDEX(aLine,headGW)==1) THEN       !Added by Anthony 230112
                nGWC = nGWC+1           				!Added by Anthony 230112
            ELSE IF(INDEX(aLine,headTransfer)==1) THEN     !Added by Anthony 230112
                nTransferC = nTransferC+1                 !Added by Anthony 230112
            END IF    
        ELSE IF (ierr<0) THEN 
            EXIT
        END IF    
      END DO      

        !Allocate memory of some parts of inp file
        IF(nTargetC>0) ALLOCATE (TargetArrayC(nTargetC))
        IF(nEvaporC>0) ALLOCATE (EvapArrayC(nEvaporC))
        IF(nRoutingC>0) ALLOCATE (RoutingArrayC(nRoutingC))
        IF(nAllocationC>0) ALLOCATE (AlloArrayC(nAllocationC))
        IF(nSourceC>0) ALLOCATE (SourceArrayC(nSourceC))
        IF(nRatingC>0) ALLOCATE (RatingArrayC(nRatingC))
        IF(nCrossC>0) ALLOCATE (CrossArrayC(nCrossC))
        IF(nPerformanceC>0) ALLOCATE (PerfArrayC(nPerformanceC))
        IF(nCostC>0) ALLOCATE (CostArrayC(nCostC))
        IF(nPowerC>0) ALLOCATE (PowerArrayC(nPowerC))
        IF(nPumpC>0) ALLOCATE (PumpArrayC(nPumpC))
        IF(nRuleC>0) ALLOCATE (RuleArrayC(nRuleC))
        IF(nDemRedC>0) ALLOCATE (DemRedArrayC(nDemredC))  !Added by Evgenii 101907 
        IF(nBalanceC>0) ALLOCATE (BalanceArrayC(nBalanceC))!Added by Evgenii 110504 
        IF(nGWC>0) ALLOCATE (GWArrayC(nGWC))   								!Added by Anthony 230112
        IF(nTransferC>0) ALLOCATE (TransferArrayC(nTransferC))            !Added by Anthony 230112
        !set the file position to the beginning of file
        REWIND(UNIT=iChangefile)
        
        nTargetc = 0
        nEvaporC = 0
        nRoutingC = 0
        nAllocationC = 0
        nSourceC = 0
        nRatingC = 0
        nCrossC = 0
        nCrossC = 0
        nPerformanceC = 0
        nCostC = 0
        nPowerC = 0
        nPumpC = 0
        nRuleC = 0
        nDemredC=0	
	    nBalanceC=0	
        nGWc = 0        	!Added by Anthony 230112
        nTransferc = 0        !Added by Anthony 230112
	  
	!Now read changes into memory
        IF(nTargetC>0) ALLOCATE (TargetArrayC(nTargetC))
        IF(nEvaporC>0) ALLOCATE (EvapArrayC(nEvaporC))
        IF(nRoutingC>0) ALLOCATE (RoutingArrayC(nRoutingC))
        IF(nAllocationC>0) ALLOCATE (AlloArrayC(nAllocationC))
        IF(nSourceC>0) ALLOCATE (SourceArrayC(nSourceC))
        IF(nRatingC>0) ALLOCATE (RatingArrayC(nRatingC))
        IF(nCrossC>0) ALLOCATE (CrossArrayC(nCrossC))
        IF(nPerformanceC>0) ALLOCATE (PerfArrayC(nPerformanceC))
        IF(nCostC>0) ALLOCATE (CostArrayC(nCostC))
        IF(nPowerC>0) ALLOCATE (PowerArrayC(nPowerC))
        IF(nPumpC>0) ALLOCATE (PumpArrayC(nPumpC))
        IF(nRuleC>0) ALLOCATE (RuleArrayC(nRuleC))
        IF(nDemRedC>0) ALLOCATE (DemRedArrayC(nDemRedC))  !Added by Evgenii 101907 
        IF(nBalanceC>0) ALLOCATE (BalanceArrayC(nBalanceC))!Added by Evgenii 110504 
        IF(nGWC>0) ALLOCATE (GWArrayC(nGWC))        !Added by Anthony 230112
        IF(nTransferC>0) ALLOCATE (TransferArrayC(nTransferC))        !Added by Anthony 230112
        !set the file position to the beginning of file
        REWIND(UNIT=iChangefile)
        

        nTargetC = 0
        nEvaporC = 0
        nRoutingC = 0
        nAllocationC = 0
        nSourceC = 0
        nRatingC = 0
        nCrossC = 0
        nCrossC = 0
        nPerformanceC = 0
        nCostC = 0
        nPowerC = 0
        nPumpC = 0
        nRuleC = 0
        nDemRedC=0		 !Added by Evgenii 100719 
	    nBalanceC=0	 !Added by Evgenii 110504
        nGWC = 0  		!Added by Anthony 230112
        nTransferC = 0  ! Added by Anthony 230112
		
        tyears=sysstat(nyear)
       do i=1,2 
        do y=1,tyears
          do p=1,MaxPolicyTypes
              do nn=1,tnodes
                  YearNodePolicyChg(i,y,p,nn)=.false.
              end do
          end do
        end do
       end do 
      do i = 1, tyears
          do j=1,2
              YearPolChg(i,j)=.false.
          end do
      end do
      begyear=sysstat(YRSTART)
        DO WHILE (.TRUE.)
          READ(UNIT=iChangefile,FMT='(A)',iostat=ierr) aLine !Evgenii took out ERR=999
		IF(ierr ==0) THEN 
                iPos = INDEX(aLine,headTarget)
                IF(iPos>=1 .AND. ALLOCATED(TargetArrayC)) THEN

                    nTargetC = nTargetC + 1
                    iPos = iPos + LEN(headTarget)
                    strTemp = aLine(iPos:)
                    read(strTemp,*)TargetArrayC(nTargetC)%Year,
     &                    TargetArrayC(nTargetC)%GroupID,
     &                                   TargetArrayC(nTargetC)%Policy,
     &                                  TargetArrayC(nTargetC)%CompType,
     &                                    TargetArrayC(nTargetC)%NodeID, 
     &                                    TargetArrayC(nTargetC)%targ,
     &                                    TargetArrayC(nTargetC)%t_co, 
     &                                   TargetArrayC(nTargetC)%EnvFlow,
     &								      TargetArrayC(nTargetC)%RefilTrig,
     &   								  TargetArrayC(nTargetC)%DemInc,
     &                                TargetArrayC(nTargetC)%srcpriority
                    
                    do nn=1,tnodes
				      if (nodeid(nn)==TargetArrayC(nTargetC)%NodeID)exit
                    end do
                    realyear=TargetArrayC(nTargetC)%Year-begyear+1
                    YearNodePolicyChg(1,realyear,Target0,nn)=.true.
                    YearPolChg(realyear,1)=.true.
                CYCLE
                END IF
                
                iPos = INDEX(aLine,headEvaporation)
                IF(iPos==1 .AND. ALLOCATED(EvapArrayC)) THEN

                    nEvaporC = nEvaporC + 1
                    iPos = iPos + LEN(headEvaporation)
                    read(aLine(iPos:),*)EvapArrayC(nEvaporC)%Year,
     &                                    EvapArrayC(nEvaporC)%GroupID, 
     &                                    EvapArrayC(nEvaporC)%Policy,
     &                                    EvapArrayC(nEvaporC)%CompType,
     &                                    EvapArrayC(nEvaporC)%ID, 
     &                                  EvapArrayC(nEvaporC)%Evaporation
                    IF(EvapArrayC(nEvaporC)%CompType == 2) THEN
                       read(aLine(iPos:),*)aVar,aVar,aVar,aVar,aVar,
     &                                   EvapArrayC(nEvaporC)%LossMethod
     
                        do ln=1,links
				            if (linkid(ln)==EvapArrayC(nEvaporC)%ID)exit
                        end do
                        realyear=EvapArrayC(nEvaporC)%Year-begyear+1
                        YearNodePolicyChg(2,realyear,Evaporation0,ln)
     &                        =.true.
                        YearPolChg(realyear,2)=.true.
                    ELSE
                        EvapArrayC(nEvaporC)%LossMethod = 0
 !                       read(aLine(iPos:),*)aVar,aVar,aVar,aVar,aVar,
 !    &                                   EvapArrayC(nEvaporC)%LossMethod
     
                        do nn=1,tnodes
				            if (nodeid(nn)==EvapArrayC(nEvaporC)%ID)exit
                        end do
                        realyear=EvapArrayC(nEvaporC)%Year-begyear+1
                        YearNodePolicyChg(1,realyear,Evaporation0,nn)
     &                        =.true.
                        YearPolChg(realyear,1)=.true.                        
                    END IF
                CYCLE
                END IF
                
                iPos = INDEX(aLine,headRouting)
                IF(iPos>=1 .AND. ALLOCATED(RoutingArrayC)) THEN

                    nRoutingC = nRoutingC + 1
                    iPos = iPos + LEN(headRouting)
                    read(aLine(iPos:),*)RoutingArrayC(nRoutingC)%Year,
     &                    RoutingArrayC(nRoutingC)%GroupID, 
     &                                RoutingArrayC(nRoutingC)%Policy,
     &                                RoutingArrayC(nRoutingC)%CompType,
     &                                RoutingArrayC(nRoutingC)%LinkID, 
     &                                RoutingArrayC(nRoutingC)%iMethod,
     &                                RoutingArrayC(nRoutingC)%L_NRTRES,
     &                                RoutingArrayC(nRoutingC)%L_a,
     &                                RoutingArrayC(nRoutingC)%L_b,
     &                                RoutingArrayC(nRoutingC)%L_c
                CYCLE
                END IF

                iPos = INDEX(aLine,headAllocation)
                IF(iPos>=1 .AND. ALLOCATED(AlloArrayC)) THEN

                    nAllocationC = nAllocationC + 1
                    iPos = iPos + LEN(headAllocation)
                   read(aLine(iPos:),*)AlloArrayC(nAllocationC)%year,
     &                  AlloArrayC(nAllocationC)%GroupID, 
     &                                AlloArrayC(nAllocationC)%Policy,
     &                                AlloArrayC(nAllocationC)%CompType,
     &                                AlloArrayC(nAllocationC)%NodeID, 
     &                              AlloArrayC(nAllocationC)%NodeOutput,
     &                                AlloArrayC(nAllocationC)%LinkID,
     &                                AlloArrayC(nAllocationC)%LinkAllo
     
                   do nn=1,tnodes
				    if (nodeid(nn)==AlloArrayC(nAllocationC)%NodeID)exit
                   end do                   
                   realyear=AlloArrayC(nAllocationC)%Year-begyear+1
                   YearNodePolicyChg(1,realyear,Allocation0,nn)=.true.
                   YearPolChg(realyear,1)=.true.
                CYCLE
                END IF
                iPos = INDEX(aLine,headSource)
                IF(iPos>=1 .AND. ALLOCATED(SourceArrayC)) THEN
                    nSourceC = nSourceC + 1
                    iPos = iPos + LEN(headSource)
                    read(aLine(iPos:),*)SourceArrayC(nSourceC)%Year, 
     &                                  SourceArrayC(nSourceC)%GroupID, 
     &                                SourceArrayC(nSourceC)%Policy,
     &                                SourceArrayC(nSourceC)%CompType,
     &                                SourceArrayC(nSourceC)%NodeID, 
     &                                SourceArrayC(nSourceC)%Supl_Node,
     & 								  SourceArrayC(nSourceC)%Source_Type,	
     &                                SourceArrayC(nSourceC)%Supl_Frac,
     &								  SourceArrayC(nSourceC)%MaxOutTS,			
     &                                SourceArrayC(nSourceC)%MaxOutYear,
     &                                SourceArrayC(nSourceC)%LnkProp !Evgenii added LNKPROP120126 for link propogation links
                do nn=1,tnodes
				      if (nodeid(nn)==SourceArrayC(nSourceC)%NodeID)exit
                end do
                realyear=SourceArrayC(nSourceC)%Year-begyear+1
                YearNodePolicyChg(1,realyear,Source0,nn)=.true.
                YearPolChg(realyear,1)=.true.   
                CYCLE             
                END IF
                
                iPos = INDEX(aLine,headRating)
                IF(iPos>=1 .AND. ALLOCATED(RatingArrayC)) THEN
                    nRatingC = nRatingC + 1
                    iPos = iPos + LEN(headRating)
                    read(aLine(iPos:),*)RatingArrayC(nRatingC)%Year,
     &                    RatingArrayC(nRatingC)%GroupID,
     &                                RatingArrayC(nRatingC)%Policy,
     &                                RatingArrayC(nRatingC)%CompType,
     &                                RatingArrayC(nRatingC)%ID, 
     &                               RatingArrayC(nRatingC)%ElevOrWidth,
     &                              RatingArrayC(nRatingC)%AreaOrEvapor,
     &                                RatingArrayC(nRatingC)%VolOrFlow
                        
                    IF(RatingArrayC(nRatingC)%CompType == 1) THEN
                       read(aLine(iPos:),*)avar,
     &                       aVar,aVar,aVar,aVar,aVar,
     &                                aVar,aVar,
     &                                RatingArrayC(nRatingC)%Seep,
     &                                RatingArrayC(nRatingC)%MaxQ,
     &                                RatingArrayC(nRatingC)%LakeQ
                        do nn=1,tnodes
				              if (nodeid(nn)==RatingArrayC(nRatingC)%ID)exit
                        end do
                        realyear=RatingArrayC(nRatingC)%Year-begyear+1
                        YearNodePolicyChg(1,realyear,Rate0,nn)=.true.
                        YearPolChg(realyear,1)=.true.       
                    ELSE
                       RatingArrayC(nRatingC)%Seep = 0.0
                       RatingArrayC(nRatingC)%MaxQ = 0.0
                       RatingArrayC(nRatingC)%LakeQ = 0.0
                       do ln=1,links
				              if (linkid(ln)==RatingArrayC(nRatingC)%ID)exit
                       end do
                       realyear=RatingArrayC(nRatingC)%Year-begyear+1
                       YearNodePolicyChg(2,realyear,Rate0,ln)=.true.
                       YearPolChg(realyear,2)=.true.                           
                    END IF
                CYCLE
                END IF
  
                iPos = INDEX(aLine,headCross)
                IF(iPos==1 .AND. ALLOCATED(CrossArrayC)) THEN
                    nCrossC = nCrossC + 1
                    iPos = iPos + LEN(headCross)
                    read(aLine(iPos:),*)CrossArrayC(nCrossC)%Year,
     &                    CrossArrayC(nCrossC)%GroupID, 
     &                                CrossArrayC(nCrossC)%Policy,
     &                                CrossArrayC(nCrossC)%CompType,
     &                                CrossArrayC(nCrossC)%LinkID, 
     &                                CrossArrayC(nCrossC)%BaseWidth,
     &                                CrossArrayC(nCrossC)%ChannelDepth, 
     &                                CrossArrayC(nCrossC)%LSlope,
     &                                CrossArrayC(nCrossC)%RSlope,
     &                                CrossArrayC(nCrossC)%UpLSlope,
     &                                CrossArrayC(nCrossC)%UpRSlope
                CYCLE
                END IF
                
               iPos = INDEX(aLine,headPerformance)
                IF(iPos>=1 .AND. ALLOCATED(PerfArrayC)) THEN
                    nPerformanceC = nPerformanceC + 1
                    iPos = iPos + LEN(headPerformance)
                    read(aLine(iPos:),*)PerfArrayC(nPerformanceC)%Year,
     &                    PerfArrayC(nPerformanceC)%GroupID,
     &                              PerfArrayC(nPerformanceC)%Policy,
     &                              PerfArrayC(nPerformanceC)%NodeID, 
     &                             PerfArrayC(nPerformanceC)%thres_limit
                    do nn=1,tnodes
			              if (nodeid(nn)==PerfArrayC(nPerformanceC)%nodeid)exit
                    end do
                    realyear=PerfArrayC(nPerformanceC)%Year-begyear+1
                    YearNodePolicyChg(1,realyear,Performance0,nn)=.true.
                    YearPolChg(realyear,1)=.true.                                   
                CYCLE
                END IF
                
                iPos = INDEX(aLine,headCost)
                IF(iPos==1 .AND. ALLOCATED(CostArrayC)) THEN
                    nCostC = nCostC + 1
                    iPos = iPos + LEN(headCost)
                    read(aLine(iPos:),*)CostArrayC(nCostC)%Year,
     &                    CostArrayC(nCostC)%GroupID,
     &                              CostArrayC(nCostC)%Policy,
     &                              CostArrayC(nCostC)%LinkID, 
     &                              CostArrayC(nCostC)%FlowCost,
     &                              CostArrayC(nCostC)%FlowEng,
     &							  CostArrayC(nCostC)%AnnCostInc
                    do ln=1,links
			              if (linkid(ln)==CostArrayC(nCostC)%linkid)exit
                    end do
                    realyear=CostArrayC(nCostC)%Year-begyear+1
                    YearNodePolicyChg(2,realyear,Cost0,ln)=.true.
                    YearPolChg(realyear,2)=.true.               
                CYCLE
                END IF

               iPos = INDEX(aLine,headPower)
                IF(iPos>=1 .AND. ALLOCATED(PowerArrayC)) THEN
                    nPowerC = nPowerC + 1
                    iPos = iPos + LEN(headPower)
                    read(aLine(iPos:),*)PowerArrayC(nPowerC)%Year,
     &                   PowerArrayC(nPowerC)%GroupID,
     &                              PowerArrayC(nPowerC)%Policy,aVar,
     &                              PowerArrayC(nPowerC)%NodeID, 
     &                              PowerArrayC(nPowerC)%HPCAP,
     &                              PowerArrayC(nPowerC)%PLANT_FACTOR,
     &                              PowerArrayC(nPowerC)%Intake_elev,
     &                              PowerArrayC(nPowerC)%Turbine_elev,
     &                              PowerArrayC(nPowerC)%Outlet_elev,
     &                              PowerArrayC(nPowerC)%ECONST,
     &                              PowerArrayC(nPowerC)%HpQMin
                
                        do nn=1,tnodes
				              if (nodeid(nn)==PowerArrayC(nPowerC)%nodeid)exit
                        end do
                        realyear=PowerArrayC(nPowerC)%Year-begyear+1
                        YearNodePolicyChg(1,realyear,Power0,nn)=.true.
                        YearPolChg(realyear,1)=.true.  
                       CYCLE
                END IF

                iPos = INDEX(aLine,headPump)
                IF(iPos>=1 .AND. ALLOCATED(PumpArrayC)) THEN
                   nPumpC = nPumpC + 1
                    iPos = iPos + LEN(headPump)
                    read(aLine(iPos:),*)PumpArrayC(nPumpC)%Year,
     &                    PumpArrayC(nPumpC)%GroupID,
     &                              PumpArrayC(nPumpC)%Policy,aVar,
     &                              PumpArrayC(nPumpC)%NodeID, 
     &                              PumpArrayC(nPumpC)%HPCAP,
     &                              PumpArrayC(nPumpC)%PLANT_FACTOR,
     &                              PumpArrayC(nPumpC)%Intake_elev,
     &                              PumpArrayC(nPumpC)%Turbine_elev,
     &                              PumpArrayC(nPumpC)%Outlet_elev,
     &                              PumpArrayC(nPumpC)%ECONST,
     &                              PumpArrayC(nPumpC)%HpQMin
                        do nn=1,tnodes
				              if (nodeid(nn)==PumpArrayC(nPumpC)%nodeid)exit
                        end do
                        realyear=PumpArrayC(nPumpC)%Year-begyear+1
                        YearNodePolicyChg(1,realyear,Pump0,nn)=.true.
                        YearPolChg(realyear,1)=.true.                  
                CYCLE
                END IF

                iPos = INDEX(aLine,headRule)
                IF(iPos>=1 .AND. ALLOCATED(RuleArrayC)) THEN
                    nRuleC = nRuleC + 1
                    iPos = iPos + LEN(headRule)
				  read(aLine(iPos:),*)RuleArrayC(nRuleC)%Year,
     &				  RuleArrayC(nRuleC)%GroupID,
     &                              RuleArrayC(nRuleC)%Policy,
     &                              RuleArrayC(nRuleC)%NodeID, 
     &                              RuleArrayC(nRuleC)%res_rule(1),
     &                              RuleArrayC(nRuleC)%res_rule(2),
     &                              RuleArrayC(nRuleC)%res_rule(3),
     &                              RuleArrayC(nRuleC)%res_rule(4),
     &                              RuleArrayC(nRuleC)%res_rule(5),
     &                              RuleArrayC(nRuleC)%res_rule(6),
     &                              RuleArrayC(nRuleC)%res_rule(7),
     &                              RuleArrayC(nRuleC)%res_rule(8)                   				
                   do nn=1,tnodes
				    if (nodeid(nn)==RuleArrayC(nRuleC)%NodeID)exit
                   end do
                    realyear=RuleArrayC(nRuleC)%Year-begyear+1
                    YearNodePolicyChg(1,realyear,Rule0,nn)=.true.
                    YearPolChg(realyear,1)=.true.                
                CYCLE
                END IF
                
			  iPos = INDEX(aLine,headDemRed)		   !Added by Evgenii 100719 
                IF(iPos>=1 .AND. ALLOCATED(DemRedArrayC)) THEN
                    nDemRedC = nDemRedC + 1
                    iPos = iPos + LEN(headDemRed)
                    read(aLine(iPos:),*)DemRedArrayC(nDemRedC)%Year,
     &                    DemRedArrayC(nDemRedC)%GroupID,
     &                            DemRedArrayC(nDemRedC)%Policy,
     &                            DemRedArrayC(nDemRedC)%NodeID, 
     &                            DemRedArrayC(nDemRedC)%DemSourceID,
     &							DemRedArrayC(nDemRedC)%Dem_Thres_limit,
     &                            DemRedArrayC(nDemRedC)%DemRedAmt
                   do nn=1,tnodes
				    if (nodeid(nn)==DemRedArrayC(nDemRedC)%NodeID)exit
                   end do
                    realyear=DemRedArrayC(nDemRedC)%Year-begyear+1
                    YearNodePolicyChg(1,realyear,DemRed0,nn)=.true.
                    YearPolChg(realyear,1)=.true.
     							


                CYCLE
                END IF
			  
			  iPos = INDEX(aLine,headBalance)		   !Added by Evgenii 110504 
                IF(iPos>=1 .AND. ALLOCATED(BalanceArrayC)) THEN
                    nBalanceC = nBalanceC + 1
                    iPos = iPos + LEN(headBalance)
				    !count number of characters in avar
				    BalanceArrayC(nBalanceC)%charBalance=trim(aLine)
                    read(aLine(iPos:),*)BalanceArrayC(nBalanceC)%Year,
     &                            BalanceArrayC(nBalanceC)%GroupID,
     &                            BalanceArrayC(nBalanceC)%Policy,
     &                            BalanceArrayC(nBalanceC)%RuleID, 
     &                            BalanceArrayC(nBalanceC)%BalMeth,
     &							  BalanceArrayC(nBalanceC)%GroupVol
                    do nn=1,tnodes
				    if (nodeid(nn)==BalanceArrayC(nBalanceC)%RuleID)exit
                   end do
                    realyear=BalanceArrayC(nBalanceC)%Year-begyear+1
                    YearNodePolicyChg(1,realyear,Balance0,nn)=.true.
                    YearPolChg(realyear,1)=.true.                
                CYCLE
                END IF
				
			iPos = INDEX(aLine,headGW) !Added by Anthony 230112
          IF(iPos>=1 .AND. ALLOCATED(GWArrayC)) THEN       !Anthony - possible that iPos>=1 needs to be changed to iPos==1
              nGWC = nGWC + 1
              iPos = iPos + LEN(headGW)
              read(aLine(iPos:),*)GWArrayC(nGWC)%Year,
     &                    GWArrayC(nGWC)%GroupID,
     &                        GWArrayC(nGWC)%Policy,
     &                        GWArrayC(nGWC)%CompType,
     &                        GWArrayC(nGWC)%LinkID,
     &                        GWArrayC(nGWC)%GWMethod,
     &                        GWArrayC(nGWC)%GWK, 
     &                        GWArrayC(nGWC)%GWElev,
     &                        GWArrayC(nGWC)%GWLength,
     &                        GWArrayC(nGWC)%GWWidth
                    do ln=1,links
			              if (linkid(ln)==GWArrayC(nGWC)%linkid)exit
                    end do
                    realyear=GWArrayC(nGWC)%Year-begyear+1
                    YearNodePolicyChg(2,realyear,GW0,ln)=.true.
                    YearPolChg(realyear,2)=.true. 
                    	                      
                CYCLE
                END IF
                
                iPos = INDEX(aLine,headTransfer)      !Added by Anthony 230112
                IF(iPos>=1 .AND. ALLOCATED(TransferArrayC)) THEN   !Anthony - possible that iPos>=1 needs to be changed to iPos==1
                    nTransferC = nTransferC + 1
                    iPos = iPos + LEN(headTransfer)
                    read(aLine(iPos:),*)TransferArrayC(nTransferC)%Year,
     &                    TransferArrayC(nTransferC)%GroupID,
     &                    TransferArrayC(nTransferC)%Policy,
     &                    TransferArrayC(nTransferC)%LinkID, 
     &                    TransferArrayC(nTransferC)%GWFromVol,
     &                    TransferArrayC(nTransferC)%GWToVol, 
     &                    TransferArrayC(nTransferC)%GWFlowFromTo
                    do ln=1,links
			              if (linkid(ln)==TransferArrayC(nTransferC)%linkid)exit
                    end do
                    realyear=TransferArrayC(nTransferC)%Year-begyear+1
                    YearNodePolicyChg(2,realyear,Transfer0,ln)=.true.
                    YearPolChg(realyear,2)=.true.
					
                CYCLE
                END IF	
				
            ELSE IF(ierr<0) THEN 
                EXIT
            END IF
          END DO
          IF(ALLOCATED(TargetArrayC)) pTargetC => TargetArrayC(:)
          IF(ALLOCATED(EvapArrayC)) pEvaporC => EvapArrayC(:)
          IF(ALLOCATED(RoutingArrayC)) pRoutingC => RoutingArrayC(:)
          IF(ALLOCATED(AlloArrayC)) pAllocationC => AlloArrayC(:)
          IF(ALLOCATED(SourceArrayC)) pSourceC => SourceArrayC(:)
          IF(ALLOCATED(RatingArrayC)) pRatingC => RatingArrayC(:)
          IF(ALLOCATED(CrossArrayC)) pCrossC => CrossArrayC(:)
          IF(ALLOCATED(PerfArrayC)) pPerformanceC => PerfArrayC(:)
          IF(ALLOCATED(CostArrayC)) pCostC => CostArrayC(:)
          IF(ALLOCATED(PowerArrayC)) pPowerC => PowerArrayC(:)
          IF(ALLOCATED(PumpArrayC)) pPumpC => PumpArrayC(:)
          IF(ALLOCATED(RuleArrayC)) pRuleC => RuleArrayC(:)
	      IF(ALLOCATED(DemRedArrayC)) pdemRedC => DemRedArrayC(:)
		  IF(ALLOCATED(BalanceArrayC)) pBalanceC=> BalanceArrayC(:)		
          IF(ALLOCATED(GWArrayC)) pGWC => GWArrayC(:)          					!Added by Anthony 230112
          IF(ALLOCATED(TransferArrayC)) pTransferC => TransferArrayC(:)           !Added by Anthony 230112
	close (iChangefile)
	success=.true.
999	continue
	return
	end subroutine
      


!**********************************************************


	Subroutine readDemFile(success)
!	Evgenii created this subroutine on 091104, it reads iras.flw
!	  gauge flow data into memory
      USE vars

	implicit none

	INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT

      
!  OUTPUT
      LOGICAL*1 success
      !COMMON: flowdata
!  Local
	INTEGER*4 j, iDemFile
      INTEGER*4 nGagesInFile,ierr
      INTEGER*4 LineCounter,nRead
	CHARACTER*10000 aLine
      CHARACTER*30 aVar
!  Call
	LOGICAL CountFileLines
!------------------------------------------------------------------------
      	!Read flow data
	success=.false.
	iDemFile = 33
      OPEN(UNIT=iDemFile, FILE='iras.dem', STATUS='old',
     &	 FORM='formatted', ERR= 999,ACTION='READ')
    		
      !Count the total number of .dem file
      LineCounter = 0
      IF(.NOT. CountFileLines(iDemFile, LineCounter)) GOTO 999 
	
	!Allocate flowdata variable to number of lines of code, number of gauges nodes
	Allocate (DemData(LineCounter,tDemSeries))
	!Read inflows into memory from iras.dem file	
	nRead = 0       !no record read
	rewind (unit=iDemFile)
	do WHILE (.TRUE.)
  	  READ(UNIT=iDemFile,FMT='(A)',iostat=ierr) aLine
	  IF(ierr ==0) THEN 
		nRead=nRead+1
		READ(aLine,*)(DemData(nRead,j), j=1, tDemSeries)	
        ELSE IF (ierr<0) THEN 
            EXIT
        else if (ierr>0) then
		continue	
	  END IF   
	end do

	if (sysstat(nRec)/=nread) then 
		write(*,*)'Demand time series not equal in length to flow time
     &	     time-series'
		goto 999
	end if
	
888   success = .true.
999   CLOSE(UNIT=iDemFile)
      if (.not. success) then
        WRITE(*,*)'Demand data reading failed'
      end if
      return
      end subroutine

	!Open .pol file. NOTE: the value of iDatafile should be different from others
 !************************************************************************     
	Subroutine readPolFile(success)
!	Evgenii created this subroutine on 091104, it reads iras.flw
!	  gauge flow data into memory
      USE vars

	implicit none
	
	INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT

      
!  OUTPUT
      LOGICAL*1 success
      !COMMON: flowdata
!  Local
	INTEGER*4 LineCounter
      INTEGER*4 nRead,ierr,ny,nn
	INTEGER*4 iPos, nNeedBufferedLinesInText

	CHARACTER*10000 aLine
      CHARACTER*30 aVar
!  Call
	LOGICAL CountFileLines

!Kang modify for improving performance
	!Initialize some gobal variables
!----------------------------------------------------------------
	
	iPolicyFile = 88
      OPEN(UNIT=iPolicyFile, FILE=TRIM(PolicyFileName), STATUS='old',
     &	 FORM='formatted', ERR= 999,ACTION='READ')
    		
      !Count the total number of .pol file
      LineCounter = 0
      IF(.NOT. CountFileLines(iPolicyFile, LineCounter)) GOTO 999
    		
      !If the total lines are less than the maximum of lines allowed for speed
      !allocate an array and then read all data into an array
      IF(LineCounter <= MaxLineForSpeed) THEN
        totalLinePolicyFile = LineCounter        
        IF(totalLinePolicyFile>0) THEN
            if (.not.allocated(PolicyFileContent))
     &        ALLOCATE (PolicyFileContent(totalLinePolicyFile))
        END IF
        
        !set the file position to the beginning of file
        REWIND(UNIT=iPolicyFile)
        
        LineCounter = 0
        
        DO WHILE (ALLOCATED(PolicyFileContent))
            READ(UNIT=iPolicyFile,FMT='(A)',iostat=ierr) aLine !Evgenii took out ERR=999
            IF(ierr ==0) THEN                 
                LineCounter = LineCounter+1
                PolicyFileContent(LineCounter) = TRIM(aLine)
            ELSE IF(ierr<0) THEN 
                EXIT
            END IF
          END DO
          
          pPolicyData => PolicyFileContent(:)
    	END IF
	success=.true.
999	continue
	End Subroutine
!************************************************************************
	Subroutine readINPFile(success)
      !This subroutine reads iras.inp	data into memory
      USE vars
	implicit none	
	INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT

      
!  OUTPUT
      LOGICAL*1 success
      !COMMON: flowdata
!  Local
	INTEGER*4 j,i,id,gagid,LineCounter, iDatafile
      INTEGER*4 nRead,ierr,ny,nn,ncols
	INTEGER*4 iPos, nNeedBufferedLinesInText

	CHARACTER*10000 aLine
      CHARACTER*30 aVar
!  Call
	LOGICAL CountFileLines
	INTEGER*4 CountColumns
!Kang modify for improving performance
	!Initialize some gobal variables
!----------------------------------------------------------------
      !Kang add 20100629

	!write(*,*) 'Welcome to readinpdata. Enjoy your stay'
            
      cRead_Sim_Data = 0
      cSysEvap = 0
      cLinkRating = 0
      cLinkLoss = 0
      cCrossSection = 0
      cLinkRouting = 0
      cTransfer = 0
      cGW = 0
      cNodeRating = 0
      cNodeEvap = 0
      cAllocation = 0
      cTarget = 0
      cTargetSource = 0
      cPowerPump = 0
      cReleaseRule = 0
      cGrpBalance = 0
      tRead_Sim_Data = 0



	NULLIFY(pFileData)
      nBufferedLines = 0
      nNeedBufferedLinesInText = 0
      Rating_BL = -1 
      Evaporation_BL = -1
      CrossSection_BL = -1
      Routing_BL = -1
      GW_BL = -1
      Allocation_BL = -1
      Target_BL = -1
      Source_BL = -1
      Power_BL = -1
      Pump_BL = -1
      Rule_BL = -1
      Balance_BL = -1
      Performance_BL = -1
	  Cost_BL = -1      
      Transfer_BL = -1

      
      totalLinePolicyFile = 0
      NULLIFY(pPolicyData)
            
      bUseBufferedBinaryData = .FALSE.      
      nTarget  = 0
      NULLIFY(pTarget)
      nEvapor = 0
      NULLIFY(pEvapor)
      nRouting = 0
      NULLIFY(pRouting)
      nAllocation = 0
      NULLIFY(pAllocation)
      nSource = 0
      NULLIFY(pSource)
      nRating = 0
      NULLIFY(pRating)
      nCross = 0
      NULLIFY(pCross)
      nPerformance = 0
      NULLIFY(pPerformance)
      nCost = 0
      NULLIFY(pCost)
      nPower = 0
      NULLIFY(pPower)
      nPump = 0
      NULLIFY(pPump)
      nRule = 0
      NULLIFY(pRule)
      nDemRed=0				 !Added by Evgenii 101907 
      NULLIFY(pDemRed)		 !Added by Evgenii 101907 
	  nBalance=0				 !Added by Evgenii 110504
      NULLIFY(pBalance)
      nGW  = 0          !Added by Anthony 230112
      NULLIFY(pGW)      !Added by Anthony 230112
      nTransfer  = 0          !Added by Anthony 230112
      NULLIFY(pTransfer)      !Added by Anthony 230112
	!Open .inp file. NOTE: the value of iDatafile should be different from others
	iDatafile = INPFileID



	!write(*,*) 'We are trying to open this file.'
	!write(*,*) 'Systemfilename is ', SysFilename
	!write(*,*) 'iDatafile is ', iDatafile
      !OPEN(UNIT=iDatafile, FILE=TRIM(SysFilename), STATUS='old',
    !&	 FORM='formatted', ERR= 999,ACTION='READ')
	REWIND (iDatafile) !JRK changed, 3 December 2010
	!write(*,*) 'We have opened this file.'
      
      !Count the total number of .inp file
      !Count the number of some parts in .inp file
      LineCounter = 0
      DO WHILE (.TRUE.)
        READ(UNIT=iDatafile,FMT='(A)',iostat=ierr) aLine !Evgenii took out ERR=999
        IF(ierr ==0) THEN 
            LineCounter = LineCounter+1
            IF(INDEX(aLine,headTarget)>0) THEN
                nTarget = nTarget+1
            ELSE IF(INDEX(aLine,headEvaporation)==1) THEN   
                nEvapor = nEvapor+1
            ELSE IF(INDEX(aLine,headRouting)==1) THEN   
                nRouting = nRouting+1
            ELSE IF(INDEX(aLine,headAllocation)==1) THEN   
                nAllocation = nAllocation+1
            ELSE IF(INDEX(aLine,headSource)==1) THEN   
                nSource = nSource+1
            ELSE IF(INDEX(aLine,headRating)==1) THEN   
                nRating = nRating+1
            ELSE IF(INDEX(aLine,headCross)==1) THEN   
                nCross = nCross+1
            ELSE IF(INDEX(aLine,headPerformance)==1) THEN   
                nPerformance = nPerformance+1
            ELSE IF(INDEX(aLine,headCost)==1) THEN   
                nCost = nCost+1
            ELSE IF(INDEX(aLine,headPower)==1) THEN   
                nPower = nPower+1
            ELSE IF(INDEX(aLine,headPump)==1) THEN   
                nPump = nPump+1
            ELSE IF(INDEX(aLine,headRule)==1) THEN   
                nRule = nRule+1
		    ELSE IF(INDEX(aLine,headDemRed)==1) THEN  !Added by Evgenii 101907  
                nDemRed = nDemRed+1
		    ELSE IF(INDEX(aLine,headBalance)==1) THEN  !Added by Evgenii 101907  
                nBalance = nBalance+1
            ELSE IF(INDEX(aLine,headGW)==1) THEN          		!Added by Anthony 230112
                nGW = nGW+1                     				!Added by Anthony 230112
            ELSE IF(INDEX(aLine,headTransfer)==1) THEN             !Added by Anthony 230112
                nTransfer = nTransfer+1                           !Added by Anthony 230112
            ELSE
                nNeedBufferedLinesInText = nNeedBufferedLinesInText + 1
            END IF    
            ELSE IF (ierr<0) THEN 
                 EXIT
            END IF    
      END DO
      
	!write(*,*) 'LineCounter before if ', LineCounter
	!write(*,*) 'MaxLineForSpeed ',MaxLineForSpeed
      !If the total lines are less than the maximum of lines allowed for speed
      !allocate an array and then read all data into an array
      IF(LineCounter <= MaxLineForSpeed) THEN
	  !write(*,*) 'nbufferedLines, before=', nBufferedLines
	  nBufferedLines = nNeedBufferedLinesInText
	  !write(*,*) 'nbufferedLines, after=', nBufferedLines
        bUseBufferedBinaryData = .TRUE.
        
        !Allocate memory of some parts of inp file
        !SysFileContent stores data in text format
        !TargetArray and EvapArray store data of "Target" in binary format
        IF(nBufferedLines>0 .and..not.allocated(SysFileContent)) 
     &   ALLOCATE (SysFileContent(nBufferedLines))
        IF(nTarget>0.and..not.allocated(TargetArray)) 
     &   ALLOCATE (TargetArray(nTarget))
        IF(nEvapor>0.and..not.allocated(EvapArray)) 
     &   ALLOCATE (EvapArray(nEvapor))
        IF(nRouting>0.and..not.allocated(RoutingArray))
     &   ALLOCATE (RoutingArray(nRouting))
        IF(nAllocation>0.and..not.allocated(AlloArray))
     &   ALLOCATE (AlloArray(nAllocation))
        IF(nSource>0.and..not.allocated(SourceArray))
     &   ALLOCATE (SourceArray(nSource))
        IF(nRating>0.and..not.allocated(RatingArray))
     &   ALLOCATE (RatingArray(nRating))
        IF(nCross>0.and..not.allocated(CrossArray)) 
     &   ALLOCATE (CrossArray(nCross))
        IF(nPerformance>0.and..not.allocated(PerfArray)) 
     &   ALLOCATE (PerfArray(nPerformance))
        IF(nCost>0.and..not.allocated(CostArray)) 
     &   ALLOCATE (CostArray(nCost))
        IF(nPower>0.and..not.allocated(PowerArray))
     &   ALLOCATE (PowerArray(nPower))
        IF(nPump>0.and..not.allocated(PumpArray))
     &   ALLOCATE (PumpArray(nPump))
        IF(nRule>0.and..not.allocated(RuleArray))
     &   ALLOCATE (RuleArray(nRule))
        IF(nDemRed>0.and..not.allocated(DemRedArray)) 
     &   ALLOCATE (DemRedArray(nDemred))  !Added by Evgenii 101907 
        IF(nBalance>0.and..not.allocated(BalanceArray)) 
     &   ALLOCATE (BalanceArray(nBalance))!Added by Evgenii 110504 
        IF(nGW>0.and..not.allocated(GWArray)) 
     &   ALLOCATE (GWArray(nGW))             							!Added by Anthony 230112
        IF(nTransfer>0.and..not.allocated(TransferArray))
     &   ALLOCATE (TransferArray(nTransfer))                       !Added by Anthony 230112
        !set the file position to the beginning of file
         REWIND(UNIT=iDatafile)
        
        LineCounter = 0
        nTarget = 0
        nEvapor = 0
        nRouting = 0
        nAllocation = 0
        nSource = 0
        nRating = 0
        nCross = 0
        nCross = 0
        nPerformance = 0
        nCost = 0
        nPower = 0
        nPump = 0
        nRule = 0
        nDemred=0		 !Added by Evgenii 100719 
	    nBalance=0	 !Added by Evgenii 110504
        nGW = 0         	!Added by Anthony 230112
        nTransfer = 0         !Added by Anthony 230112

        DO WHILE (.TRUE.)
          READ(UNIT=iDatafile,FMT='(A)',iostat=ierr) aLine !Evgenii took out ERR=999
		IF(ierr ==0) THEN 
                iPos = INDEX(aLine,headTarget)
                IF(iPos>=1 .AND. ALLOCATED(TargetArray)) THEN

                    nTarget = nTarget + 1
                    iPos = iPos + LEN(headTarget)
                    strTemp = aLine(iPos:)
                    read(strTemp,*)TargetArray(nTarget)%GroupID,  
     &                                    TargetArray(nTarget)%Policy,
     &                                    TargetArray(nTarget)%CompType,
     &                                    TargetArray(nTarget)%NodeID, 
     &                                    TargetArray(nTarget)%targ,
     &                                    TargetArray(nTarget)%t_co, 
     &                                    TargetArray(nTarget)%EnvFlow,
     &								      TargetArray(nTarget)%RefilTrig,
     &   								  TargetArray(nTarget)%DemInc,
     &                                TargetArray(nTarget)%srcpriority
                    
                    !IF(CountColumns(strTemp)==8) THEN
                    !  read(strTemp,*)aVar,aVar,aVar,aVar,aVar,aVar,aVar,
!     &              !                      TargetArray(nTarget)%RefilTrig
                    !ELSE
                    !    TargetArray(nTarget)%RefilTrig = 0.0
                    !END IF
				  
				                      
                CYCLE
                END IF
                
                iPos = INDEX(aLine,headEvaporation)
                IF(iPos==1 .AND. ALLOCATED(EvapArray)) THEN

                    nEvapor = nEvapor + 1
                    iPos = iPos + LEN(headEvaporation)
                    read(aLine(iPos:),*)EvapArray(nEvapor)%GroupID, 
     &                                    EvapArray(nEvapor)%Policy,
     &                                    EvapArray(nEvapor)%CompType,
     &                                    EvapArray(nEvapor)%ID, 
     &                                    EvapArray(nEvapor)%Evaporation
                    IF(EvapArray(nEvapor)%CompType == 2) THEN
                        read(aLine(iPos:),*)aVar,aVar,aVar,aVar,aVar,
     &                                    EvapArray(nEvapor)%LossMethod
                    ELSE
                        EvapArray(nEvapor)%LossMethod = 0
                    END IF
                CYCLE
                END IF
                
                iPos = INDEX(aLine,headRouting)
                IF(iPos>=1 .AND. ALLOCATED(RoutingArray)) THEN

                    nRouting = nRouting + 1
                    iPos = iPos + LEN(headRouting)
                    read(aLine(iPos:),*)RoutingArray(nRouting)%GroupID, 
     &                                RoutingArray(nRouting)%Policy,
     &                                RoutingArray(nRouting)%CompType,
     &                                RoutingArray(nRouting)%LinkID, 
     &                                RoutingArray(nRouting)%iMethod,
     &                                RoutingArray(nRouting)%L_NRTRES,
     &                                RoutingArray(nRouting)%L_a,
     &                                RoutingArray(nRouting)%L_b,
     &                                RoutingArray(nRouting)%L_c
                CYCLE
                END IF

                iPos = INDEX(aLine,headAllocation)
                IF(iPos>=1 .AND. ALLOCATED(AlloArray)) THEN

                    nAllocation = nAllocation + 1
                    iPos = iPos + LEN(headAllocation)
                    read(aLine(iPos:),*)AlloArray(nAllocation)%GroupID, 
     &                                AlloArray(nAllocation)%Policy,
     &                                AlloArray(nAllocation)%CompType,
     &                                AlloArray(nAllocation)%NodeID, 
     &                                AlloArray(nAllocation)%NodeOutput,
     &                                AlloArray(nAllocation)%LinkID,
     &                                AlloArray(nAllocation)%LinkAllo
                CYCLE
                END IF
                iPos = INDEX(aLine,headSource)
                IF(iPos>=1 .AND. ALLOCATED(SourceArray)) THEN
                    nSource = nSource + 1
                    iPos = iPos + LEN(headSource)
                    read(aLine(iPos:),*)SourceArray(nSource)%GroupID, 
     &                                SourceArray(nSource)%Policy,
     &                                SourceArray(nSource)%CompType,
     &                                SourceArray(nSource)%NodeID, 
     &                                SourceArray(nSource)%Supl_Node,
     & 								  SourceArray(nSource)%Source_Type,	!Evgenii 110404 added source_type
     &                                SourceArray(nSource)%Supl_Frac,
     &								  SourceArray(nSource)%MaxOutTS,			
     &                                SourceArray(nSource)%MaxOutYear,
     &                                SourceArray(nSource)%LnkProp
                CYCLE
                END IF
                
                iPos = INDEX(aLine,headRating)
                IF(iPos>=1 .AND. ALLOCATED(RatingArray)) THEN
                    nRating = nRating + 1
                    iPos = iPos + LEN(headRating)
                    read(aLine(iPos:),*)RatingArray(nRating)%GroupID,
     &                                RatingArray(nRating)%Policy,
     &                                RatingArray(nRating)%CompType,
     &                                RatingArray(nRating)%ID, 
     &                                RatingArray(nRating)%ElevOrWidth,
     &                                RatingArray(nRating)%AreaOrEvapor,
     &                                RatingArray(nRating)%VolOrFlow
                        
                    IF(RatingArray(nRating)%CompType == 1) THEN
                       read(aLine(iPos:),*)aVar,aVar,aVar,aVar,aVar,
     &                                aVar,aVar,
     &                                RatingArray(nRating)%Seep,
     &                                RatingArray(nRating)%MaxQ,
     &                                RatingArray(nRating)%LakeQ
                    ELSE
                       RatingArray(nRating)%Seep = 0.0
                       RatingArray(nRating)%MaxQ = 0.0
                       RatingArray(nRating)%LakeQ = 0.0
                    END IF
                CYCLE
                END IF
  
                iPos = INDEX(aLine,headCross)
                IF(iPos==1 .AND. ALLOCATED(CrossArray)) THEN
                    nCross = nCross + 1
                    iPos = iPos + LEN(headCross)
                    read(aLine(iPos:),*)CrossArray(nCross)%GroupID, 
     &                                CrossArray(nCross)%Policy,
     &                                CrossArray(nCross)%CompType,
     &                                CrossArray(nCross)%LinkID, 
     &                                CrossArray(nCross)%BaseWidth,
     &                                CrossArray(nCross)%ChannelDepth, 
     &                                CrossArray(nCross)%LSlope,
     &                                CrossArray(nCross)%RSlope,
     &                                CrossArray(nCross)%UpLSlope,
     &                                CrossArray(nCross)%UpRSlope
                CYCLE
                END IF
                
               iPos = INDEX(aLine,headPerformance)
                IF(iPos>=1 .AND. ALLOCATED(PerfArray)) THEN
                    nPerformance = nPerformance + 1
                    iPos = iPos + LEN(headPerformance)
                    read(aLine(iPos:),*)PerfArray(nPerformance)%GroupID,
     &                              PerfArray(nPerformance)%Policy,
     &                              PerfArray(nPerformance)%NodeID, 
     &                              PerfArray(nPerformance)%thres_limit
                CYCLE
                END IF
                
                iPos = INDEX(aLine,headCost)
                IF(iPos==1 .AND. ALLOCATED(CostArray)) THEN
                    nCost = nCost + 1
                    iPos = iPos + LEN(headCost)
                    read(aLine(iPos:),*)CostArray(nCost)%GroupID,
     &                              CostArray(nCost)%Policy,
     &                              CostArray(nCost)%LinkID, 
     &                              CostArray(nCost)%FlowCost,
     &                              CostArray(nCost)%FlowEng,
     &							  CostArray(nCost)%AnnCostInc
                CYCLE
                END IF

               iPos = INDEX(aLine,headPower)
                IF(iPos>=1 .AND. ALLOCATED(PowerArray)) THEN
                    nPower = nPower + 1
                    iPos = iPos + LEN(headPower)
                    read(aLine(iPos:),*)PowerArray(nPower)%GroupID,
     &                              PowerArray(nPower)%Policy,aVar,
     &                              PowerArray(nPower)%NodeID, 
     &                              PowerArray(nPower)%HPCAP,
     &                              PowerArray(nPower)%PLANT_FACTOR,
     &                              PowerArray(nPower)%Intake_elev,
     &                              PowerArray(nPower)%Turbine_elev,
     &                              PowerArray(nPower)%Outlet_elev,
     &                              PowerArray(nPower)%ECONST,
     &                              PowerArray(nPower)%HpQMin
                CYCLE
                END IF

                iPos = INDEX(aLine,headPump)
                IF(iPos>=1 .AND. ALLOCATED(PumpArray)) THEN
                   nPump = nPump + 1
                    iPos = iPos + LEN(headPump)
                    read(aLine(iPos:),*)PumpArray(nPump)%GroupID,
     &                              PumpArray(nPump)%Policy,aVar,
     &                              PumpArray(nPump)%NodeID, 
     &                              PumpArray(nPump)%HPCAP,
     &                              PumpArray(nPump)%PLANT_FACTOR,
     &                              PumpArray(nPump)%Intake_elev,
     &                              PumpArray(nPump)%Turbine_elev,
     &                              PumpArray(nPump)%Outlet_elev,
     &                              PumpArray(nPump)%ECONST,
     &                              PumpArray(nPump)%HpQMin
                CYCLE
                END IF

                iPos = INDEX(aLine,headRule)
                IF(iPos>=1 .AND. ALLOCATED(RuleArray)) THEN
                    nRule = nRule + 1
                    iPos = iPos + LEN(headRule)
				  read(aLine(iPos:),*)RuleArray(nRule)%GroupID,
     &                              RuleArray(nRule)%Policy,
     &                              RuleArray(nRule)%NodeID, 
     &                              RuleArray(nRule)%res_rule(1),
     &                              RuleArray(nRule)%res_rule(2),
     &                              RuleArray(nRule)%res_rule(3),
     &                              RuleArray(nRule)%res_rule(4),
     &                              RuleArray(nRule)%res_rule(5),
     &                              RuleArray(nRule)%res_rule(6),
     &                              RuleArray(nRule)%res_rule(7),
     &                              RuleArray(nRule)%res_rule(8)                   				
                CYCLE
                END IF
                
			  iPos = INDEX(aLine,headDemRed)		   !Added by Evgenii 100719 
                IF(iPos>=1 .AND. ALLOCATED(DemRedArray)) THEN
                    nDemRed = nDemRed + 1
                    iPos = iPos + LEN(headDemRed)
                    read(aLine(iPos:),*)DemRedArray(nDemRed)%GroupID,
     &                            DemRedArray(nDemRed)%Policy,
     &                            DemRedArray(nDemRed)%NodeID, 
     &                            DemRedArray(nDemRed)%DemSourceID,
     &							  DemRedArray(nDemRed)%Dem_Thres_limit,
     &                            DemRedArray(nDemRed)%DemRedAmt
     							


                CYCLE
                END IF
			  
			  iPos = INDEX(aLine,headBalance)		   !Added by Evgenii 110504 
                IF(iPos>=1 .AND. ALLOCATED(BalanceArray)) THEN
                    nBalance = nBalance + 1
                    iPos = iPos + LEN(headBalance)
				  !count number of characters in avar
				  BalanceArray(nBalance)%charBalance=trim(aLine)
                    read(aLine(iPos:),*)BalanceArray(nBalance)%GroupID,
     &                            BalanceArray(nBalance)%Policy,
     &                            BalanceArray(nBalance)%RuleID, 
     &                            BalanceArray(nBalance)%BalMeth,
     &							  BalanceArray(nBalance)%GroupVol
                CYCLE
                END IF
                iPos = INDEX(aLine,headGW)                          !Added by Anthony 230112
                IF(iPos>=1 .AND. ALLOCATED(GWArray)) THEN

                    nGW = nGW + 1
                    iPos = iPos + LEN(headGW)
                    read(aLine(iPos:),*)GWArray(nGW)%GroupID,
     &                             GWArray(nGW)%Policy,
     &                             GWArray(nGW)%CompType,
     &                             GWArray(nGW)%LinkID,
     &                             GWArray(nGW)%GWMethod,
     &                             GWArray(nGW)%GWK,
     &                             GWArray(nGW)%GWElev,
     &                             GWArray(nGW)%GWLength,
     &                             GWArray(nGW)%GWWidth
                              
      			                      
                CYCLE
                END IF
                
                
              iPos = INDEX(aLine,headTransfer)                          !Added by Anthony 230112
              IF(iPos>=1 .AND. ALLOCATED(TransferArray)) THEN

                  nTransfer = nTransfer + 1
                  iPos = iPos + LEN(headTransfer)
                  read(aLine(iPos:),*)TransferArray(nTransfer)%GroupID,
     &                             TransferArray(nTransfer)%Policy,
     &                             TransferArray(nTransfer)%LinkID, 
     &                             TransferArray(nTransfer)%GWFromVol,
     &                             TransferArray(nTransfer)%GWToVol, 
     &                             TransferArray(nTransfer)%GWFlowFromTo
                			                      
                CYCLE
                END IF

                LineCounter = LineCounter+1
                SysFileContent(LineCounter) = TRIM(aLine)
            ELSE IF(ierr<0) THEN 
                EXIT
            END IF
          END DO
          IF(ALLOCATED(SysFileContent)) pFileData => SysFileContent(:)
          IF(ALLOCATED(TargetArray)) pTarget => TargetArray(:)
          IF(ALLOCATED(EvapArray)) pEvapor => EvapArray(:)
          IF(ALLOCATED(RoutingArray)) pRouting => RoutingArray(:)
          IF(ALLOCATED(AlloArray)) pAllocation => AlloArray(:)
          IF(ALLOCATED(SourceArray)) pSource => SourceArray(:)
          IF(ALLOCATED(RatingArray)) pRating => RatingArray(:)
          IF(ALLOCATED(CrossArray)) pCross => CrossArray(:)
          IF(ALLOCATED(PerfArray)) pPerformance => PerfArray(:)
          IF(ALLOCATED(CostArray)) pCost => CostArray(:)
          IF(ALLOCATED(PowerArray)) pPower => PowerArray(:)
          IF(ALLOCATED(PumpArray)) pPump => PumpArray(:)
          IF(ALLOCATED(RuleArray)) pRule => RuleArray(:)
		  IF(ALLOCATED(DemRedArray)) pdemRed => DemRedArray(:)
		  IF(ALLOCATED(BalanceArray)) pBalance => BalanceArray(:)		
					 !Added by Evgenii 100719 
          IF(ALLOCATED(GWArray)) pGW => GWArray(:)      					!Added by Anthony 230112
          IF(ALLOCATED(TransferArray)) pTransfer => TransferArray(:)      !Added by Anthony 230112
	end if
	rewind (iDatafile)
	success=.true.
999	continue
	return
	end subroutine

!************************************************************************
!Kang add for improving performance.
      LOGICAL FUNCTION CountFileLines(iFile, iLines)
!  INPUT
      INTEGER*4 iFile
!  OUTPUT
      INTEGER*4 iLines
      LOGICAL bSuccess
!  LOCAL
      CHARACTER*256 aLine
!------------------------------------------------------------------------
      bSuccess = .FALSE.
      iLines = 0
      DO WHILE (.TRUE.)
        READ(UNIT=iFile,FMT='(A)',iostat=ierr) aLine !Evgenii took out ERR=999
        IF(ierr ==0) THEN 
            iLines = iLines+1
        ELSE IF (ierr<0) THEN 
            EXIT
        END IF    
      END DO
      
      bSuccess = .TRUE.
      
999   CountFileLines = bSuccess

      END FUNCTION
      
!************************************************************************
!Kang add for improving performance.
      INTEGER FUNCTION CountColumns(aLine)
!  INPUT
      CHARACTER*256 aLine
!  OUTPUT
      INTEGER nColumns
!  LOCAL
      INTEGER I, nLen
      LOGICAL bStart
      CHARACTER*256 strTemp
      CHARACTER aChar
!------------------------------------------------------------------------
      nColumns = 0
      strTemp = aLine
      nLen = LEN(strTemp)
      bStart = .FALSE.
      DO I=1, nLen
        aChar = strTemp(I:I)
        IF(ICHAR(aChar)==32 .OR. ICHAR(aChar)==9) THEN
            IF(bStart) THEN
                nColumns = nColumns + 1
                bStart = .FALSE.
            END IF
        ELSE 
            IF(ICHAR(aChar)==10 .OR. ICHAR(aChar)==13 .OR. 
     &         aChar=='!') THEN
                IF(bStart) THEN
                    nColumns = nColumns + 1
                    bStart = .FALSE.
                    EXIT
                END IF
            ELSE
                bStart = .TRUE.
            END IF     
        END IF         
      END DO
      
      CountColumns = nColumns

      END FUNCTION      