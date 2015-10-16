!Copyright (c) 2009, 2010 by University College London, Cornell University
!Authors:
!Evgenii Matrosov (evgenii.matrosov@ucl.ac.uk), Julien Harou (j.harou@ucl.ac.uk), 
!This program is free software under the General Public Licence, GPL (>=v2)
!Read the 'GPL License.txt' file distributed with this source code for a full license statement.
!**************************************************************************************************
	subroutine StoragePerformance()
	!Created by Evgenii 1006
	!This subroutine calculates the time-step performance measures for storage nodes
	USE vars
	IMPLICIT NONE
	INCLUDE 'IRAS_SYS.INC'
	INCLUDE 'NODE.INC'
	INCLUDE 'LINK.INC'
!  Common 
!storage_limit(i,nn),perf_pts(nn),last_T_fail(i,nn),sto_perf_node(nn),nStoreFail(i,nn),
!min_store_reach(nn),nTime_Steps_Tot(i,nn),max_sto_fail_dur(i,nn)
!must make min-store_reach initialized to initial storage of storage node at beginning of simulation

!  INPUT
!	real::    ESTO(nodmax)
!  Local
	INTEGER*4 nn,i,j	
	real::    DMD_TS(nodmax)
	logical:: failure(thres_pts_max,nodmax), infailure=.false.
!------------------------------------------------------------------------

	Do i=1,tnodes
		do j=1,thres_pts_max
			failure(j,i)=.false.
		end do
		!for flow demands get target time step demand flow (mil m3/time-step)
		if (dmdnode(i) .and. capn(i)==0.0)
     &		DMD_TS(i)=DMD_TARG(i)*DAYSPRPRD
	end do
		!for first time step  set minimum of performance variable at  targets
	if (sysstat(nTimeStep)==1)then
	  do i=1,tnodes
	  	 !for storage
		 if (CAPN(i)>0.0)then
			 min_store_reach(i)=max(0.0,STOINT(i)) !100617 added by Evgenii for performance measures
		 else
			 min_store_reach(i)=DMD_TS(i)
		 end if
	  end do
	end if
	DO NN = 1, TNODES
	 IF(sto_perf_node(nn)) THEN  !must initialize sto_perf_node(nn) as false
		   !For storage nodes- If current storage less than record minimum (in percentage)
		 if (capn(nn)>0.0) then 
			  min_store_reach(nn)=
     &		    min(min_store_reach(nn),esto(nn)/CAPN(NN))
		   	   
			  min_store_reach(nn)=
     &			max(min_store_reach(nn),0.0)

		   !For flow nodes (in percentage)
		 else 
			  min_store_reach(nn)=
     &			min(min_store_reach(nn),inflow(nn)/DMD_TS(nn))
			  min_store_reach(nn)=
     &			max(min_store_reach(nn),0.0)
		  
		 end if	 
		 do i=1,perf_pts(nn)			
		   !if below threshold for both storage and flow target nodes set infailure to true		   
		   infailure=.false.
		   if(capn(nn)>0.0 .and. 
     & 		   ESTO(NN)<=thres_limit(i,nn)*CAPN(NN))then
     			   infailure=.true.
			end if

     		   if (capn(nn)==0.0.and.inflow(nn)<=thres_limit(i,nn)
     &	       *DMD_TS(nn))then
				infailure=.true.
		   end if
		   !if node in failure for current threshold
		   if(infailure) then
			 YearFailEvent(sysstat(SIM_YEAR),nn,i)=.true.
			!if new failure (previous sub-time step was not failure)
		 	 if (.not. last_T_fail(i,nn)) then
				failure(i,nn)=.true.
				!increase count of total reservoir failure (for the whole run)
				nTime_Steps_Tot(i,nn)=nTime_Steps_Tot(i,nn)+1
				!start count for duration of specific failure
				nTime_Steps(i,nn)=nTime_Steps(i,nn)+1
				!Increase total failures in run by one
				nStoreFail(i,nn)=nStoreFail(i,nn)+1
				!Increase failures for year by one
				!Must initialize variable below to 0
					
				
			!if failure and previous sub-time step was also failure
			 else
				failure(i,nn)=.true.
				!increase count of total reservoir failure (for the whole run)
				nTime_Steps_Tot(i,nn)=nTime_Steps_Tot(i,nn)+1
				!increase count for duration of specific failure
				nTime_Steps(i,nn)=nTime_Steps(i,nn)+1
			 end if
	  !if above threshold
		   else if (.not. infailure) then  !if above threshold
			 failure(i,nn)=.false.
			 if (last_T_fail(i,nn)) then 
				!replace max failure duration if this failure longer than any other failure in the simulation
				max_sto_fail_dur(i,nn)=max(max_sto_fail_dur(i,nn),
     &				nTime_Steps(i,nn))
				!Reset counter for specific fail
				nTime_Steps(i,nn)=0
			 end if
	       end if !if below thres_limit(i,nn)
	     !after failure analysis, update failure variable for future sum-time step	  
		 last_T_fail(i,nn)=failure(i,nn)
		
		end do !end do for all pts
	  end if !sto_perf_node(nn)
	end do
	DO NN = 1, TNODES
	!now calculate addition performance criteria for flow target nodes
	    if(.not.ResvNode(NN).and..not.natlak(nn) .and. DMDNODE(nn)) THEN
			!Annual demand target
			Annual_DMD(nn)=Annual_DMD(nn)+DMD_TARG(NN)*DAYSPRPRD
			!Annual shortage
			Annual_SHRTG(nn)=Annual_SHRTG(nn)+TSflw_DEFICIT(NN)
			!Time step SI variable
			if(DMD_TARG(NN)/=0.0)then
				TsSIsum(nn)=TsSIsum(nn)+
     &			(TSflw_DEFICIT(NN)/(DMD_TARG(NN)*DAYSPRPRD))**2
		    end if
			!Time-step shortage
			Ts_SHRTG(nn)=TS_SHRTG(nn)+TSflw_DEFICIT(NN)
              !Time-step surplus
              TS_Surp(nn)=TS_Surp(nn)+TSflw_Surplus(NN)
			!Time-step demand
			Ts_DMD_Sum(nn)=Ts_DMD_Sum(nn)+DMD_TARG(NN)*DAYSPRPRD
	  end if
	end do !End for all nodes

   
	END
	!************************************************************************
	subroutine annualPerformance()
	!Created by Evgenii 1006
	!
	IMPLICIT NONE
	INCLUDE 'IRAS_SYS.INC'
	INCLUDE 'NODE.INC'
	INCLUDE 'LINK.INC'
	!  INPUT

	!  OUTPUT

	!  LOCAL
	real:: YearSD
	INTEGER*4 LN, NN
!------------------------------------------------------------------------
	!Annual cost of energy
	DO LN = 1, LINKS			
		!Incorporate energy price increase	
		AnnualCost(ln)=(AnnualEng(ln)*flowcost(ln)
     &		  *(1+AnnCostInc(ln))**(sysstat(SIM_YEAR)-1))
		!Discount energy costs
		!AnnualCost(ln)=AnnualCost(ln)/(1.035**(sysstat(SIM_YEAR)-1))


		GlobannCost(ln)=GlobannCost(ln)+AnnualCost(ln)
		!maxannualcost(ln)=max(AnnualCost(ln),maxannualcost(ln))
		
		
		!Reset accumulation variable for new year
		AnnualCost(ln)=0.0
		AnnualEng(ln)=0.0
	END DO

	DO NN = 1, TNODES
	  if (.not.ResvNode(NN).and..not.natlak(nn) .and. DMDNODE(nn)) then
!   		 annualSI(NN)=annualSI(NN)+(Annual_SHRTG(nn)
!     &		 /Annual_DMD(nn))**2
	     !annualSDshrtg(nn)=annualSDshrtg(nn)+Annual_SHRTG(nn)
		 !annualSDdem(nn)=annualSDdem(nn)+Annual_DMD(nn)
		 !Reset flow demand annual vairiables to 0 for new year
		 !Calculate annual SD for year that just occured
		 if(sysstat(SIM_YEAR)==1) minYearSD(nn)=100 !sysstat(SIM_YEAR)==1 changed from sysstat(SIM_YEAR)==2, thanks to Joe K. for finding this bug 290712
		 YearSD=0.0
		 if (Annual_DMD(nn)/=0.0) then
			YearSD=(1-(Annual_SHRTG(NN)/Annual_DMD(nn)))*100
		 else	
			YearSD=0.0
		 end if
		 
		 minYearSD(nn)=min(minYearSD(nn),YearSD)
		 Annual_SHRTG(NN)=0.0
		 Annual_DMD(nn)=0.0
	  end if
	end do
	end
!************************************************************************
	subroutine LinkCost_and_Power()
		!Created by Evgenii 1006
	!
	IMPLICIT NONE
	INCLUDE 'IRAS_SYS.INC'
	INCLUDE 'NODE.INC'
	INCLUDE 'LINK.INC'
	!  INPUT

	!  OUTPUT

	!  LOCAL

	INTEGER*4 LN
!------------------------------------------------------------------------
	DO LN = 1, LINKS
!	write(*,*)ln
	!Power expended (only .not.powerlink(ln), so as to not count energy from pumping/power calculations)
		if (FlowEng(ln)>0.0.and. .not. powerlink(ln).and.
     &	    .not.GWLINK(ln).and.BQLN(ln)>0.0) then
			  EngTot(ln)=EngTot(ln)+FlowEng(ln)*BQLN(ln)*1E6
			  AnnualEng(ln)=AnnualEng(ln)+FlowEng(ln)*BQLN(ln)*1E6
			  
					
                    
		!Accumulate pumping energy from power/pump calc in hydsim.for
		else if(pumplink(ln) .and. ENERGY(ln)< 0.0) then
		      EngTotHydSim(ln)=EngTotHydSim(ln)+ENERGY(ln)
			  AnnualEng(ln)=AnnualEng(ln)+abs(ENERGY(ln))
			  		
		end if
		totflow(ln)=totflow(ln)+ BQLN(ln)
		

	!Costs
	!	if(FlowCost(ln)>0.0.and..not.GWLINK(ln).and.BQLN(ln)>0.0)then
	!	  CostTot(ln)=CostTot(ln)+FlowCost(ln)*BQLN(ln)*1E6
	!	end if	
		!Power obtained 
		if(ENERGY(ln)>0.0)then
			TotPower(ln)=TotPower(ln)+ENERGY(ln)
		end if


	END DO
	End
!************************************************************************
	subroutine PerformanceOutput(csustain, cflow, cdemand)
	!Created by Evgenii 1006
	!This subroutine prints out the time step results to a file found in the /out directory.
	USE vars
	IMPLICIT NONE
	INCLUDE 'IRAS_SYS.INC'
	INCLUDE 'NODE.INC'
	INCLUDE 'LINK.INC'
	!  Common 
	!storage_limit(i,nn),perf_pts(nn),last_T_fail(i,nn),perf_node(nn),nStoreFail(i,nn),
	!min_store_reach(nn),nTime_Steps_Tot(i,nn),max_sto_fail_dur(i,nn)
	!must make min-store_reach initialized to initial storage of storage node at beginning of simulation

	!  INPUT

	!  Local
	real*4 ave_sto_fail(thres_pts_max,nodmax),Ntime, Nfail,SD(nodmax)
	real*4 SI(nodmax),SIts(nodmax),SDts(nodmax),totfineng(lnkmax)
	real*4 totgain(lnkmax),totexp(lnkmax),globPower,globcost
	real*4 globeng,globrev,reliability(thres_pts_max,nodmax)
	real*4 nAnnualFailreal,years
	INTEGER*4 nn,ln,ioutperf,nperfNodes,s,f,i,j,nSimYears,ny
	INTEGER*4 TotSimYears, nAnnualFail(thres_pts_max,nodmax),flowrecs
	logical:: calcFail(lnkmax)
	character(len=30)xrun
	character(len=30)performance_filename
      character*1 csustain, cdemand
      character*2 cflow 

	real*4 minDem,PowerCost,londem,reliabilityTS(thres_pts_max,nodmax)
      real*4 SumSD,SumSurp
      integer*4 ENVTotalFailures,ENVmaxTimeFail,ENVTotalTime
      integer*4 ioutperfMOEA
      integer*4 LOS1TotalFailures,LOS1maxTimeFail,LOS1TotalTime
      integer*4 LOS3TotalFailures,LOS3maxTimeFail,LOS3TotalTime      
      integer*4 LOS3TotalFailuresAbberton,LOS3maxTimeFailAbberton
      integer*4 LOS3TotalTimeAbberton
      integer*4 LOS3TotalFailuresHanning,LOS3maxTimeFailHanning
      integer*4 LOS3TotalTimeHanning      
      real*4 minEssex, minAll, userflow
      !real*4 avgglobcost         
      !real*4 avgSumSD            
      !real*4 avgENVTotalFailures 
      !real*4 avgENVmaxTimeFail   
      !real*4 avgENVTotalTime     
      !real*4 avgLOS1TotalFailures
      !real*4 avgLOS1maxTimeFail  
      !real*4 avgLOS1TotalTime    
      !real*4 avgLOS3TotalFailures
      !real*4 avgLOS3maxTimeFail  
      !real*4 avgLOS3TotalTime   

      

!------------------------------------------------------------------------
	ioutperf=30
	globeng=0.0;globrev=0.0 
	globPower=0.0;globcost=0.0
	do ln=1,links
	  calcFail(ln)=.false.
	  totgain(ln)=0.0
	end do
	
	!First make sure that maximum failure duration is true
	!This needs to be done because performance nodes are below 
	!failure thresholds during the last sub-time step, 	max_sto_fail_dur
	!was not updated in the  StoragePerformance
	
	do ny=1,sysstat(nyear)
	  DO NN = 1, TNODES
		  do i=1,perf_pts(nn)
				nAnnualFail(i,nn)=0
		  end do
	  end do
	end do
	do ny=1,sysstat(nyear)
	  DO NN = 1, TNODES
	    IF(sto_perf_node(nn)) THEN  
		  do i=1,perf_pts(nn)
			 !YearFailEvent is how many failures for each zone occured in each year for each performance node
			 if(YearFailEvent(ny,nn,i)==.true.) then
				nAnnualFail(i,nn)=nAnnualFail(i,nn)+1
			 end if
		  end do
	    end if
	  end do
	end do
	  
		
	
	
	DO NN = 1, TNODES
	 IF(sto_perf_node(nn)) THEN  !must initialize sto_perf_node(nn) as false
		 do i=1,perf_pts(nn)
			!replace max failure duration if this failure longer than any other failure in the simulation
				max_sto_fail_dur(i,nn)=max(max_sto_fail_dur(i,nn),
     &				nTime_Steps(i,nn)) 		
		end do !end do for all pts
	 end if !sto_perf_node(nn)
	end do

	
	xrun=' '
	write(xrun,*)sysstat(run)
	xrun=adjustl(xrun)
!	Add run number to output file name (gaugename.out)
		performance_filename='performance'//'-'//cdemand//'-'
     &     //csustain//'-'//trim(cflow)//'.out'
	!performance_filename='performance.out'
	OPEN(UNIT = ioutperf, FILE ='output/'//trim(performance_filename), 
     &	STATUS='replace') 
	!calculate total time in failures for each failure zone
	!For storage nodes
	write(UNIT=ioutperf,FMT=*)'Performance Nodes'
	DO NN = 1, TNODES
	  IF(sto_perf_node(nn)) THEN 
!		nsto_perfNodes=nsto_perfNodes+1
		do i=1,perf_pts(nn)
			nTime=0
			nFail=0
			!Must use dummy real variable to convert integers into reals before doind division
			nAnnualFailreal=real(nAnnualFail(i,nn))
		    years= sysstat(nyear)
			reliability(i,nn)=1- nAnnualFailreal/years 
            !Time step reliability 
            reliabilityTS(i,nn)=1 - real(nTime_Steps_Tot(i,nn))
     &         /real(sysstat(nrec))
			if (nStoreFail(i,nn)>0) then
				!Must use dummy real variable to convert integers into reals before doind division
				nTime=nTime_Steps_Tot(i,nn)
				Nfail=nStoreFail(i,nn)
				ave_sto_fail(i,nn)=nTime/Nfail
				!If system was always below fail threshold, max_sto_fail_dur=ave_sto_fail so that max_fail_dur wouldnt be 0
				if(max_sto_fail_dur(i,nn)==0.0) then
					max_sto_fail_dur(i,nn)=ave_sto_fail(i,nn)
				end if
			else
				ave_sto_fail(i,nn)=nStoreFail(i,nn)			
			end if
			 
		end do
	  end if
      END DO
      ENVTotalFailures=0.0
      ENVmaxTimeFail=0.0
      ENVTotalTime=0.0

      SumSD=0.0
	  MinEssex=0.0
	  minAll = 0.0
      LOS1TotalFailures=0;LOS1maxTimeFail=0;LOS1TotalTime=0
      LOS3TotalFailures=0;LOS3maxTimeFail=0;LOS3TotalTime=0    
      LOS3TotalFailuresAbberton=0;LOS3maxTimeFailAbberton=0
      LOS3TotalTimeAbberton=0
      LOS3TotalFailuresHanning=0;LOS3maxTimeFailHanning=0
      LOS3TotalTimeHanning=0      
	!Write storage performance measures
	minEssex=0.0
      sumsurp=0.0
	do i=1,tnodes 
		NN = NodSEQ(I)
		IF(sto_perf_node(nn)) THEN
            write(UNIT=ioutperf,FMT=9)nname(nn),
     &        'Minimum_Dem_Reached (percent):', min_store_reach(nn)*100
            write(UNIT=ioutperf,FMT=10)'Threshold',
     &            'Total Failures for Threshold',
     &            'Ave. Time Steps in Failures',
     &            'Max Time Spent in Failure',
     &          'Total Time Spent in Failure',
     &            'Annual Reliability',
     &          'Time Step Reliability'
             
              if (trim(nname(nn))=='Abberton_Res') then
                  minEssex=minEssex+min_store_reach(nn)*100
                  minAll=minAll+min_store_reach(nn)*100
              end if
              !if (trim(nname(nn))=='Hanningfield_Res') then
                !  minEssex=minEssex+min_store_reach(nn)*100
               !   minAll=minAll+min_store_reach(nn)*100
              !end if
                   if(trim(nname(nn))=='Alton_Res')then
                      minAll=minAll+min_store_reach(nn)*100
                   end if
                  if(trim(nname(nn))=='Ardliegh_Reservior')then
                      minAll=minAll+min_store_reach(nn)*100
                   end if             
                   if(trim(nname(nn))=='Ravens_Hollow')then
                      minAll=minAll+min_store_reach(nn)*100
                   end if
                   if(trim(nname(nn))=='Pitsford')then
                      minAll=minAll+min_store_reach(nn)*100
                   end if
                   if(trim(nname(nn))=='Rutland')then
                      minAll=minAll+min_store_reach(nn)*100
                   end if             
                   if(trim(nname(nn))=='Grafham_Res')then
                     minAll=minAll+min_store_reach(nn)*100     
                   end if      
                   if(trim(nname(nn))=='Covenham_Res')then
                      minAll=minAll+min_store_reach(nn)*100
                   end if
	 
			do j=1,perf_pts(nn)
				 write(UNIT=ioutperf,FMT=11)
     &			 j,nStoreFail(j,nn),ave_sto_fail(j,nn)				
     &             ,max_sto_fail_dur(j,nn),nTime_Steps_Tot(j,nn),
     &			 reliability(j,nn),reliabilityTS(j,nn)
                   
                   !Environmental flow performance
                   if(j==4.and.trim(nname(nn))=='Pitsford')then
                   ENVTotalFailures=ENVTotalFailures+nStoreFail(j,nn)
        ENVmaxTimeFail=ENVmaxTimeFail+max_sto_fail_dur(j,nn)
                   ENVTotalTime=ENVTotalTime+nTime_Steps_Tot(j,nn)
                   end if
                   if(j==4.and.trim(nname(nn))=='Rutland')then
                   ENVTotalFailures=ENVTotalFailures+nStoreFail(j,nn)
        ENVmaxTimeFail=ENVmaxTimeFail+max_sto_fail_dur(j,nn)
                   ENVTotalTime=ENVTotalTime+nTime_Steps_Tot(j,nn)                   
                   end if              
                   if(j==4.and.trim(nname(nn))=='Grafham_Res')then
                   ENVTotalFailures=ENVTotalFailures+nStoreFail(j,nn)
        ENVmaxTimeFail=ENVmaxTimeFail+max_sto_fail_dur(j,nn)
                   ENVTotalTime=ENVTotalTime+nTime_Steps_Tot(j,nn)                   
                   end if     
                   
                   
                   !LOS 1, j=3, hosepipe ban
                   
                   if(j==3.and.trim(nname(nn))=='Alton_Res')then
                   LOS1TotalFailures=LOS1TotalFailures+nStoreFail(j,nn)
        LOS1maxTimeFail=LOS1maxTimeFail+max_sto_fail_dur(j,nn)
                   LOS1TotalTime=LOS1TotalTime+nTime_Steps_Tot(j,nn)
                   end if
                  if(j==3.and.trim(nname(nn))=='Ardliegh_Reservior')then
                   LOS1TotalFailures=LOS1TotalFailures+nStoreFail(j,nn)
        LOS1maxTimeFail=LOS1maxTimeFail+max_sto_fail_dur(j,nn)
                   LOS1TotalTime=LOS1TotalTime+nTime_Steps_Tot(j,nn)                  
                   end if              
                   if(j==3.and.trim(nname(nn))=='Ravens_Hollow')then
                   LOS1TotalFailures=LOS1TotalFailures+nStoreFail(j,nn)
        LOS1maxTimeFail=LOS1maxTimeFail+max_sto_fail_dur(j,nn)
                   LOS1TotalTime=LOS1TotalTime+nTime_Steps_Tot(j,nn)  
                   end if
                   if(j==3.and.trim(nname(nn))=='Pitsford')then
                   LOS1TotalFailures=LOS1TotalFailures+nStoreFail(j,nn)
        LOS1maxTimeFail=LOS1maxTimeFail+max_sto_fail_dur(j,nn)
                   LOS1TotalTime=LOS1TotalTime+nTime_Steps_Tot(j,nn)
                   end if
                   if(j==3.and.trim(nname(nn))=='Rutland')then
                   LOS1TotalFailures=LOS1TotalFailures+nStoreFail(j,nn)
        LOS1maxTimeFail=LOS1maxTimeFail+max_sto_fail_dur(j,nn)
                   LOS1TotalTime=LOS1TotalTime+nTime_Steps_Tot(j,nn)                 
                   end if              
                   if(j==3.and.trim(nname(nn))=='Grafham_Res')then
                   LOS1TotalFailures=LOS1TotalFailures+nStoreFail(j,nn)
        LOS1maxTimeFail=LOS1maxTimeFail+max_sto_fail_dur(j,nn)
                   LOS1TotalTime=LOS1TotalTime+nTime_Steps_Tot(j,nn)                           
                   end if       
                   if(j==3.and.trim(nname(nn))=='Covenham_Res')then
                   LOS1TotalFailures=LOS1TotalFailures+nStoreFail(j,nn)
        LOS1maxTimeFail=LOS1maxTimeFail+max_sto_fail_dur(j,nn)
                   LOS1TotalTime=LOS1TotalTime+nTime_Steps_Tot(j,nn)
                   end if
                   !Level 2 for Essex
                   if(j==2.and.trim(nname(nn))=='Abberton_Res')then
                   LOS1TotalFailures=LOS1TotalFailures+nStoreFail(j,nn)
        LOS1maxTimeFail=LOS1maxTimeFail+max_sto_fail_dur(j,nn)
                   LOS1TotalTime=LOS1TotalTime+nTime_Steps_Tot(j,nn)                            
                   end if       
          !         if(j==2.and.trim(nname(nn))=='Hanningfield_Res')then
          !         LOS1TotalFailures=LOS1TotalFailures+nStoreFail(j,nn)
        !LOS1maxTimeFail=max(LOS1maxTimeFail,max_sto_fail_dur(j,nn))
         !          LOS1TotalTime=LOS1TotalTime+nTime_Steps_Tot(j,nn)
         !          end if       
                   
                   
                       
                   !LOS 3, j=5, rationing
                   if(j==5.and.trim(nname(nn))=='Alton_Res')then
                   LOS3TotalFailures=LOS3TotalFailures+nStoreFail(j,nn)
        LOS3maxTimeFail=LOS3maxTimeFail+max_sto_fail_dur(j,nn)
                   LOS3TotalTime=LOS3TotalTime+nTime_Steps_Tot(j,nn)
                   end if
                  if(j==5.and.trim(nname(nn))=='Ardliegh_Reservior')then
                   LOS3TotalFailures=LOS3TotalFailures+nStoreFail(j,nn)
        LOS3maxTimeFail=LOS3maxTimeFail+max_sto_fail_dur(j,nn)
                   LOS3TotalTime=LOS3TotalTime+nTime_Steps_Tot(j,nn)                   
                   end if              
                   if(j==5.and.trim(nname(nn))=='Ravens_Hollow')then
                   LOS3TotalFailures=LOS3TotalFailures+nStoreFail(j,nn)
        LOS3maxTimeFail=LOS3maxTimeFail+max_sto_fail_dur(j,nn)
                   LOS3TotalTime=LOS3TotalTime+nTime_Steps_Tot(j,nn) 
                   end if
                   if(j==5.and.trim(nname(nn))=='Pitsford')then
                   LOS3TotalFailures=LOS3TotalFailures+nStoreFail(j,nn)
        LOS3maxTimeFail=LOS3maxTimeFail+max_sto_fail_dur(j,nn)
                   LOS3TotalTime=LOS3TotalTime+nTime_Steps_Tot(j,nn)
                   end if
                   if(j==5.and.trim(nname(nn))=='Rutland')then
                    LOS3TotalFailures=LOS3TotalFailures+nStoreFail(j,nn)
        LOS3maxTimeFail=LOS3maxTimeFail+max_sto_fail_dur(j,nn)
                   LOS3TotalTime=LOS3TotalTime+nTime_Steps_Tot(j,nn)                  
                   end if              
                   if(j==5.and.trim(nname(nn))=='Grafham_Res')then
                   LOS3TotalFailures=LOS3TotalFailures+nStoreFail(j,nn)
        LOS3maxTimeFail=LOS3maxTimeFail+max_sto_fail_dur(j,nn)
                   LOS3TotalTime=LOS3TotalTime+nTime_Steps_Tot(j,nn)                           
                   end if       
                   if(j==5.and.trim(nname(nn))=='Covenham_Res')then
                   LOS3TotalFailures=LOS3TotalFailures+nStoreFail(j,nn)
        LOS3maxTimeFail=LOS3maxTimeFail+max_sto_fail_dur(j,nn)
                   LOS3TotalTime=LOS3TotalTime+nTime_Steps_Tot(j,nn)
                   end if
                   !Level 2 for Essex
                   if(j==4.and.trim(nname(nn))=='Abberton_Res')then
                   LOS3TotalFailures=LOS3TotalFailures+nStoreFail(j,nn)
        LOS3maxTimeFail=LOS3maxTimeFail+max_sto_fail_dur(j,nn)
                   LOS3TotalTime=LOS3TotalTime+nTime_Steps_Tot(j,nn)                           
                  end if       
          !         if(j==4.and.trim(nname(nn))=='Hanningfield_Res')then
          !         LOS3TotalFailuresHanning=nStoreFail(j,nn)
          !         LOS3maxTimeFail=max_sto_fail_dur(j,nn)
          !         LOS3TotalTime=nTime_Steps_Tot(j,nn)
          !         end if                    
              end do
                   
              write(UNIT=ioutperf,fmt=*)
		end if
      end do
      !Abberton and Hanningfield performance is aggregate
      !if(LOS3TotalFailuresAbberton>1.and.LOS3TotalFailuresHanning>1)then
      !    LOS3TotalFailures=LOS3TotalFailures+1
      !end if
	WRITE(UNIT=ioutperf,FMT=*)
	write(UNIT=ioutperf,FMT=*)'SI and SD for flow demands'
	write(UNIT=ioutperf,FMT=12)'Node Name',       !,'SD' 'SI Annual',
     &	'SI Time-Step','SD Time-Step','Min Annual SD','Surplus'
	flowrecs=sysstat(nrec)
	DO NN = 1, TNODES
	  IF(DMDNODE(nn).and. .not.ResvNode(nn) .and.  
     &	  .not. gwnode(nn) .and. .not.  natlak(nn))THEN 	  ! IF(CAPN(NN) == 0.0 .and. DMDNODE(nn))
		!If simulations start in the middle of the year, these calculations are slightly off
		!nSimYears=SYSSTAT(YREND)-SYSSTAT(YRSTART)+1 
		!SI(nn)=100/nSimYears*annualSI(NN)
		!SD(nn)=(1-(annualSDshrtg(nn)/annualSDdem(nn)))*100
		if 	(TsSIsum(nn)>0.0) then
			SIts(nn)=100/real(flowrecs)*TsSIsum(nn)
		else
			 SIts(nn)=0.0
		end if
		if 	(Ts_DMD_Sum(nn)>0.0) then
			!SDts(nn)=(1-(Ts_SHRTG(nn)/Ts_DMD_Sum(nn)))*100
			SDts(nn)=Ts_SHRTG(nn)

		else	
			SDts(nn)=0.0
		end if
		write(UNIT=ioutperf,FMT=13)nname(nn),SIts(nn) 
     &		,SDts(nn),minYearSD(nn),TS_Surp(nn)
       end if
      if (trim(nname(nn))=='Dem_Fenland')then
          SumSD=SumSD+SDts(nn)
          SumSurp=SumSurp+TS_Surp(nn)
      end if
      if(trim(nname(nn))=='Dem_NNorfolk_Coast') then 
          SumSD=SumSD+SDts(nn)
          SumSurp=SumSurp+TS_Surp(nn)
      end if    
      if (trim(nname(nn))=='NorthfolkRuralDem') then 
          SumSD=SumSD+SDts(nn) 
          SumSurp=SumSurp+TS_Surp(nn)
      end if    
      if(trim(nname(nn))=='Dem_Norwich_Broads') then 
          SumSD=SumSD+SDts(nn)
          SumSurp=SumSurp+TS_Surp(nn)
      end if     
      if (trim(nname(nn))=='Dem_NorthernCentral') then 
          SumSD=SumSD+SDts(nn)
          SumSurp=SumSurp+TS_Surp(nn)
      end if    
      if (trim(nname(nn))=='Dem_Hartismere') then 
          SumSD=SumSD+SDts(nn)
          SumSurp=SumSurp+TS_Surp(nn)
      end if    
      if (trim(nname(nn))=='Dem_Blyth') then 
          SumSD=SumSD+SDts(nn)
          SumSurp=SumSurp+TS_Surp(nn)
      end if    
      if (trim(nname(nn))=='Dem_Lincoln_Inland_SE') then 
          SumSD=SumSD+SDts(nn)
          SumSurp=SumSurp+TS_Surp(nn)
      end if    
      if(trim(nname(nn))=='Dem_C_Lincoln_North') then 
          SumSD=SumSD+SDts(nn)
          SumSurp=SumSurp+TS_Surp(nn)
      end if    
      if(trim(nname(nn))=='Dem_Camb_and_WSuffolk') then 
          SumSD=SumSD+SDts(nn)
          SumSurp=SumSurp+TS_Surp(nn)
      end if    
      if (trim(nname(nn))=='Dem_Affinity1-4') then 
          SumSD=SumSD+SDts(nn)
          SumSurp=SumSurp+TS_Surp(nn)
      end if    
      if (trim(nname(nn))=='Dem_Affinity5') then 
          SumSD=SumSD+SDts(nn)
          SumSurp=SumSurp+TS_Surp(nn)
      end if    
      if (trim(nname(nn))=='Dem_West_Lincolnshire') then
           SumSD=SumSD+SDts(nn)
           SumSurp=SumSurp+TS_Surp(nn)
      end if
      !if(trim(nname(nn))=='Dem_Essex_AW')SumSD=SumSD+SDts(nn)
      if(trim(nname(nn))=='Dem_Cambridge') then 
          SumSD=SumSD+SDts(nn)
          SumSurp=SumSurp+TS_Surp(nn)
      end if    


      if(trim(nname(nn))=='Rutland_Dem') then 
          SumSurp=SumSurp+TS_Surp(nn)
      end if     
      if(trim(nname(nn))=='Dem_SouthRuthamford') then 
          SumSurp=SumSurp+TS_Surp(nn)
      end if     
      if(trim(nname(nn))=='Dem_ELincoln_Coast') then 
          SumSurp=SumSurp+TS_Surp(nn)
      end if     
      if(trim(nname(nn))=='Dem_Essex') then 
          SumSurp=SumSurp+TS_Surp(nn)
      end if     
      if(trim(nname(nn))=='Dem_East_Suffolk') then 
          SumSurp=SumSurp+TS_Surp(nn)
      end if 

      if(trim(nname(nn))=='Dem_AffinityEast') then 
          SumSurp=SumSurp+TS_Surp(nn)
      end if 
      

      
      if(trim(nname(nn))=='Dem_Affinity6') then 
          SumSurp=SumSurp+TS_Surp(nn)
      end if       
      

      end do

	!Now costs and power for links
	WRITE(UNIT=ioutperf,FMT=*)
	write(UNIT=ioutperf,FMT=*)'Energy Gain and Loss'
	write(UNIT=ioutperf,FMT=14)'Link Name','Total Flow',
     &	'Total Power Generated (MWh)',       !,'SD'
     &	'Total Power Consumed (MWh)','TotPower Revenue(1000s)',
     &	'Tot Power Cost(M)'
	DO ln = 1, Links
		if(EngTotHydSim(ln)<0.0 .and. EngTot(ln)>0.0)then
			calcFail(ln)=.true.
		end if
		
		!Combine energy arrays into TotFinEng (total KWh for run)
		!EngTotHydSim is pumping energy calculated by pumping algorithm	(negative value since its energy consummed)
		!EngTot is calculated by (total flow on link x FLOWENG(KWh/m3) ) 
		TotFinEng(ln)=-EngTotHydSim(ln)+EngTot(ln)
		if(TotFinEng(ln)>0.0)then
		 !totexp(ln)=TotFinEng(ln)*flowcost(ln)
		 !globeng=globeng+TotFinEng(ln)
		 !globcost=globcost+GlobannCost(ln)
		end if
          !Power links
		if(powerlink(ln))then
			totgain(ln)=flowcost(ln)*TotPower(ln)
			globPower=globPower+TotPower(ln)
			globrev=globrev+totgain(ln)
		endif
		if (powerlink(ln).or.pumplink(ln).or.FlowEng(ln)/=0.0.or.
     &		flowcost(ln)/=0.0) then	
             userflow=Totflow(ln)
             call UndoUnitConversion(2,uFlow,userflow)              
             globeng=globeng+userflow*FlowEng(ln)
             globcost=globcost+userflow*flowcost(ln)     
			write(UNIT=ioutperf,FMT=15)lname(ln),userflow,
     &	    totpower(ln)/1000000,userflow*FlowEng(ln)/1000,
     &      totgain(ln)/1000,userflow*flowcost(ln)/1000000
		end if	
	end do
	TotSimYears=(SYSSTAT(YREND)-SYSSTAT(YRSTART))+1 !+1 added to this equation 2907112, thanks to Joe K. for spotting this

	WRITE(UNIT=ioutperf,FMT=*)
	write(UNIT=ioutperf,FMT=18)'Tot Energy (GWh) Loss', 'Gain'
	write(UNIT=ioutperf,FMT=19)globeng/1000000,
     &	globPower/1000000
	write(UNIT=ioutperf,FMT=18)'Annual Ave. Loss (GWh)','Gain '
	write(UNIT=ioutperf,FMT=19)globeng/1000000/TotSimYears,
     &	  globPower/1000/TotSimYears
     &	 	
	WRITE(UNIT=ioutperf,FMT=*)
	write(UNIT=ioutperf,FMT=20)'Total Cost (1000s)', 
     &	'Ave. Annual Cost (1000s)','Total Revenue (1000s)', 
     &	'Ave. Annual Revenue (1000s)'
	write(UNIT=ioutperf,FMT=21)globcost/1000000,
     & globcost/TotSimYears/1000000,globrev/1000,
     & globrev/TotSimYears/1000
	close(unit=ioutperf)
	 






      
                   !Write MOEA performance file
              ioutperfMOEA=132
         OPEN(UNIT = ioutperfMOEA, FILE ='output/MEOA_Performance.out', 
     &      STATUS='replace')        
       !    Write MOEA performance file
            write(UNIT = ioutperfMOEA,FMT=22)'Operating_costs_M£',
     &      'Non_Storage_Reliabilty','ENVTotalFailures','ENVmaxTimeFail'
     &     ,'ENVTotalTime','LOS1TotalFailures','LOS1maxTimeFail',
     &     'LOS1TotalTime','LOS3TotalFailures','LOS3maxTimeFail',
     &      'LOS3TotalTime'
            write(UNIT =ioutperfMOEA,FMT=23)globcost/1000000,SumSD,
     &      ENVTotalFailures,ENVmaxTimeFail,ENVTotalTime,
     &      LOS1TotalFailures,LOS1maxTimeFail,
     &      LOS1TotalTime,LOS3TotalFailures,LOS3maxTimeFail,
     &      LOS3TotalTime
          write(UNIT = ioutperfMOEA,FMT=22)'minEssex','Min. storage all'
     &          ,'Surplus'
          write(UNIT =ioutperfMOEA,FMT=24)minEssex, minAll,SumSurp
          close(unit=ioutperfMOEA)           
         
	 















      


9     FORMAT(A30,A30,F30.2)
10    FORMAT(A30,45A30)
11	FORMAT(I30,I30,F30.2,2I30,2F30.3)
12	FORMAT(A30,5A30)
13    FORMAT(A30,5F30.2)
14	FORMAT(A30,6A30)
15    FORMAT(A30,6F30.2)

16    FORMAT(1A30,1F30.2)	
17    FORMAT(A30,I30,F30.2)	
18    FORMAT(2A30)
19    FORMAT(2F30.4)	
20    FORMAT(4A30)
21    FORMAT(4F30.4)		
22	FORMAT(11A23)
23	FORMAT(F23.4,F23.2,9I23)       
24	FORMAT(F23.4,F23.2,F23.2)     
	end