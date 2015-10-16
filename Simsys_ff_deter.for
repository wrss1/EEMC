!Copyright (c) 2009,2010 by University College London, Cornell University
!Authors:
!G Pegram, Daniel P. Loucks (dpl3@cornell.edu), Marshall Taylor, Evgenii Matrosov (evgenii.matrosov@ucl.ac.uk),
!Julien Harou (j.harou@ucl.ac.uk), Peter French, Huicheng Zhou
!This program is free software under the General Public Licence, GPL (>=v2)
!Read the 'GPL License.txt' file distributed with this source code for a full license statement.

! modified by Fai Fung 08/2012 - added lines to write to netcdf
! to compile on Oxford systems:
!     ifort modules.for Allocations.for changes_mod.for Changes.for Finterp.for Flwsim.for Gwlnkq.for Hydsim.for InitSys.for Loss.for iras_variables.for output_iras_netcdf.f90 Output_ff.for Performance.for rd_flow_day.for Read_sim_data.for Readfiles.for Release.for Selseq_rules.for source_priorities.f90 Selseq_ff.for Simlnk.for sort.f90 Simsys_ff.for -o iras2010 -lnetcdf -lnetcdff
!
!	 *************************************************************************************

!This file includes:
! SIMULATION
! SIMSYS
! GetMonthDay

      PROGRAM SIMULATION
! CALL:
!  * ReadSimDef()
!  * GetUserUnits()
!  * read_network_data()
!  * ReadGageYearUnits()
!  * SIMSYS()
      !Kang remove 20100629 because the following statement is not supported on linux
    !  USE DFLIB
	USE VARS
	USE CHANGES
	IMPLICIT NONE
      
	INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
 
!  Local

	LOGICAL*1 success
	real*4 time_end, time_begin
	real*4 time_sim_beg,time_sim_end
	INTEGER*4 sim,r
      character*1 csustain,  cdemand
      character*2 cflow
      integer*4 flow, demand, sustain,n,nn,ln,s
	
!------------------------------------------------------------------------     
      !Begin running time calculation, Evgenii
	CALL CPU_TIME (time_begin) 
	!moved from Initsys() for WREA study
      DefFile		   = 'input/iras.def'				    !System defintions, Evgenii
	FlowFileName   = 'input/iras.gag'					!Gauge file now just gauge conversion factor defintions, Evgenii
      PolicyFileName = 'input/iras.pol'					!Seasonal and annual parameter changes defined in policy file, Evgenii
      !SysFileName    = 'iras.inp'					!Network definitions, Evgenii
	UnitsFileName  = 'input/iras.unt'					!Evgenii added units file
	!flwdatafilename = 'iras.flw'				!Evgenii added flw file, these are time-series flows at gauge nodes that used to be in Gage file.      
      r=0
      do demand=3,3          
          !Do a case on each demand, multiplying by 0.95, 1,and 1.05
      do sustain=3,3          
          ! change INP name to 1,2,3
          select case (sustain)
          case (1)
              SysFileName    = 'input/iras1.inp' !Confirmed
          case (2)
              SysFileName    = 'input/iras2.inp' !Likely
          case (3)
              SysFileName    = 'input/iras3.inp' !Unknown
          end select
      do flow=1,1
          ! change flow file to 1-20, 2-20. 3-20
          write(csustain,'(i1)')sustain
          write(cflow,'(i2)')flow
          if(flow<10)write(cflow,'(i1)')flow
          if(flow>9)write(cflow,'(i2)')flow
          write(cdemand,'(i1)')demand
          flwdatafilename='input/iras'//trim(csustain)//'-'
     &         //trim(cflow)//'.flw'
          
	!System initialization
	call INITSYS(time_begin)
      select case (demand)
      case (1)
          do n=1,ntarget
		!Fenland
		if(pTarget(n)%NodeID==60410)then
			pTarget(n)%targ=pTarget(n)%targ*0.95	
		end if
        !North_Norfolk_Coast
		if(pTarget(n)%NodeID==60510)then
			pTarget(n)%targ=pTarget(n)%targ*0.95	
		end if		
		!Norwich_and_the_Broads
		if(pTarget(n)%NodeID==60610)then
			pTarget(n)%targ=pTarget(n)%targ*0.95	
          end if											    						
		!Norfolk_Rural
		if(pTarget(n)%NodeID==60710)then
			pTarget(n)%targ=pTarget(n)%targ*0.95	
          end if
          		!Cambridgeshire_and_West_Suffolk
		if(pTarget(n)%NodeID==60810)then
			pTarget(n)%targ=pTarget(n)%targ*0.95	
		end if
        !East_Suffolk
		if(pTarget(n)%NodeID==60910)then
			pTarget(n)%targ=pTarget(n)%targ*0.95	
		end if		
		!Essex_AW
		if(pTarget(n)%NodeID==601010)then
			pTarget(n)%targ=pTarget(n)%targ*0.95	
          end if						
		!CambridgeWater
		if(pTarget(n)%NodeID==601110)then
			pTarget(n)%targ=pTarget(n)%targ*0.95	
		end if
        !AffinityEast
		if(pTarget(n)%NodeID==601210)then
			pTarget(n)%targ=pTarget(n)%targ*0.95	
		end if		
		!Affinity5
		if(pTarget(n)%NodeID==601510)then
			pTarget(n)%targ=pTarget(n)%targ*0.95	
          end if				          
 		!NorthernCentral
		if(pTarget(n)%NodeID==601610)then
			pTarget(n)%targ=pTarget(n)%targ*0.95	
          end if
        !Hartismere
		if(pTarget(n)%NodeID==601710)then
			pTarget(n)%targ=pTarget(n)%targ*0.95	
		end if	          
        !Blyth
		if(pTarget(n)%NodeID==601810)then
			pTarget(n)%targ=pTarget(n)%targ*0.95	
		end if		
		!Ravens_Hollow_Dem
		if(pTarget(n)%NodeID==60759)then
			pTarget(n)%targ=pTarget(n)%targ*0.95	
          end if				         
  		!Pitsford
		if(pTarget(n)%NodeID==60756)then
			pTarget(n)%targ=pTarget(n)%targ*0.95	
          end if
        !Grafham_Dem
		if(pTarget(n)%NodeID==782)then
			pTarget(n)%targ=pTarget(n)%targ*0.95	
		end if	          
        !Rutland
		if(pTarget(n)%NodeID==60752)then
			pTarget(n)%targ=pTarget(n)%targ*0.95	
		end if		
		!Dem_Lincoln_Inland_SE
		if(pTarget(n)%NodeID==60010)then
			pTarget(n)%targ=pTarget(n)%targ*0.95	
          end if		         
 		!East_Lincoln_Coast
		if(pTarget(n)%NodeID==60110)then
			pTarget(n)%targ=pTarget(n)%targ*0.95	
          end if				         
  		!West_Lincolnshire
		if(pTarget(n)%NodeID==602510)then
			pTarget(n)%targ=pTarget(n)%targ*0.95	
          end if
        !Central_Lincolnshire_North
		if(pTarget(n)%NodeID==60310)then
			pTarget(n)%targ=pTarget(n)%targ*0.95	
		end if	          
        !Affinity6
		if(pTarget(n)%NodeID==601310)then
			pTarget(n)%targ=pTarget(n)%targ*0.95	
		end if		
		!Affinity1-4
		if(pTarget(n)%NodeID==601410)then
			pTarget(n)%targ=pTarget(n)%targ*0.95	
          end if		         
		!Essex
		if(pTarget(n)%NodeID==601910)then
			pTarget(n)%targ=pTarget(n)%targ*0.95	
          end if	          
          end do
      case (2)
          continue
      case (3)
          do n=1,ntarget
		!Fenland
		if(pTarget(n)%NodeID==60410)then
			pTarget(n)%targ=pTarget(n)%targ*1.05	
		end if
        !North_Norfolk_Coast
		if(pTarget(n)%NodeID==60510)then
			pTarget(n)%targ=pTarget(n)%targ*1.05	
		end if		
		!Norwich_and_the_Broads
		if(pTarget(n)%NodeID==60610)then
			pTarget(n)%targ=pTarget(n)%targ*1.05	
          end if											    						
		!Norfolk_Rural
		if(pTarget(n)%NodeID==60710)then
			pTarget(n)%targ=pTarget(n)%targ*1.05	
          end if
          		!Cambridgeshire_and_West_Suffolk
		if(pTarget(n)%NodeID==60810)then
			pTarget(n)%targ=pTarget(n)%targ*1.05	
		end if
        !East_Suffolk
		if(pTarget(n)%NodeID==60910)then
			pTarget(n)%targ=pTarget(n)%targ*1.05	
		end if		
		!Essex_AW
		if(pTarget(n)%NodeID==601010)then
			pTarget(n)%targ=pTarget(n)%targ*1.05	
          end if						
		!CambridgeWater
		if(pTarget(n)%NodeID==601110)then
			pTarget(n)%targ=pTarget(n)%targ*1.05	
		end if
        !AffinityEast
		if(pTarget(n)%NodeID==601210)then
			pTarget(n)%targ=pTarget(n)%targ*1.05	
		end if		
		!Affinity5
		if(pTarget(n)%NodeID==601510)then
			pTarget(n)%targ=pTarget(n)%targ*1.05	
          end if				          
 		!NorthernCentral
		if(pTarget(n)%NodeID==601610)then
			pTarget(n)%targ=pTarget(n)%targ*1.05	
          end if
        !Hartismere
		if(pTarget(n)%NodeID==601710)then
			pTarget(n)%targ=pTarget(n)%targ*1.05	
		end if	          
        !Blyth
		if(pTarget(n)%NodeID==601810)then
			pTarget(n)%targ=pTarget(n)%targ*1.05	
		end if		
		!Ravens_Hollow_Dem
		if(pTarget(n)%NodeID==60759)then
			pTarget(n)%targ=pTarget(n)%targ*1.05	
          end if				         
  		!Pitsford
		if(pTarget(n)%NodeID==60756)then
			pTarget(n)%targ=pTarget(n)%targ*1.05	
          end if
        !Grafham_Dem
		if(pTarget(n)%NodeID==782)then
			pTarget(n)%targ=pTarget(n)%targ*1.05	
		end if	          
        !Rutland
		if(pTarget(n)%NodeID==60752)then
			pTarget(n)%targ=pTarget(n)%targ*1.05	
		end if		
		!Dem_Lincoln_Inland_SE
		if(pTarget(n)%NodeID==60010)then
			pTarget(n)%targ=pTarget(n)%targ*1.05	
          end if		         
 		!East_Lincoln_Coast
		if(pTarget(n)%NodeID==60110)then
			pTarget(n)%targ=pTarget(n)%targ*1.05	
          end if				         
  		!West_Lincolnshire
		if(pTarget(n)%NodeID==602510)then
			pTarget(n)%targ=pTarget(n)%targ*1.05	
          end if
        !Central_Lincolnshire_North
		if(pTarget(n)%NodeID==60310)then
			pTarget(n)%targ=pTarget(n)%targ*1.05	
		end if	          
        !Affinity6
		if(pTarget(n)%NodeID==601310)then
			pTarget(n)%targ=pTarget(n)%targ*1.05	
		end if		
		!Affinity1-4
		if(pTarget(n)%NodeID==601410)then
			pTarget(n)%targ=pTarget(n)%targ*1.05	
          end if		         
		!Essex
		if(pTarget(n)%NodeID==601910)then
			pTarget(n)%targ=pTarget(n)%targ*1.05	
          end if	          
          end do          
      end select
      do ln=1,links
	    YearQLN(ln)=0.0
	end do
	do nn=1,tnodes
	do s=1,MXSUPLY
		AnnualOut(s,nn)=0.0
      end do
      end do
	!r=0
	sysstat(run)=0 
	do sim=1,sysstat(nRuns) 
		 call cpu_time(time_sim_beg)
		 r=r+1
		 CALL  SIMSYS (success, csustain, cflow, cdemand)
		 !Print performance output file
		 call PerformanceOutput(csustain, cflow, cdemand)				 
		 CALL CPU_TIME ( time_sim_end )
		 write(*,*)'Run number ',r,' ',time_sim_end-
     &			time_sim_beg,'From start: ',time_sim_end-
     &			time_begin
           
           
           CLOSE(UNIT=INPFileID);CLOSE(UNIT=iPolicyFile)
           CLOSE(UNIT=iOutFile);CLOSE(UNIT=ioutNodes)
           CLOSE(UNIT=ioutlinks)
           call de_allocate_vars()
	end do !End of flow factor runs loop	

      end do !demands
	end do !sustain
	end do !flow	
	
	
	if (.not.success) then
		WRITE(*,*)''
		WRITE(*,*)'Some errors occured in simulation!'
	end if
	!Close INP and POL files
	!CLOSE(UNIT=INPFileID);CLOSE(UNIT=iPolicyFile)
	call de_allocate_vars()
888	CALL CPU_TIME ( time_end )
	WRITE(*,*)'Time of operation was ',time_end - time_begin,'seconds' 
	 	 
      STOP

      END PROGRAM

!	 *************************************************************************************
      SUBROUTINE INITSYS(time_begin)
 	!This subroutine initializes the system. It is run once, even in the case of multiple
	!runs. It initializes system variables, reads the input files into memory and allocates
	!allocatable arrays.

	USE vars
      implicit none	
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
	!Local
	logical*1 success
	real*4 time_begin,time_end_readfile
	INTEGER*4 i,ly
 !------------------------------------------------------------------------ 
	!DefFile		   = 'iras.def'				    !System defintions, Evgenii
	!FlowFileName   = 'iras.gag'					!Gauge file now just gauge conversion factor defintions, Evgenii
 !     PolicyFileName = 'iras.pol'					!Seasonal and annual parameter changes defined in policy file, Evgenii
 !     SysFileName    = 'iras.inp'					!Network definitions, Evgenii
	!UnitsFileName  = 'iras.unt'					!Evgenii added units file
	!flwdatafilename = 'iras.flw'				!Evgenii added flw file, these are time-series flows at gauge nodes that used to be in Gage file.

	!read file names and year policyGrp
      CALL ReadSimDef(success)			
      if (.not. success) then 
		write(*,*)'Error in reading iras.def file'
		stop
	end if
	

	!Get units conversion multipliers which convert user units to internal unit just after having been read	(iras.unt)
      call GetUserUnits(success)
	IF(.not.success)then
        WRITE(*,*)'Error in reading units file:', UnitsFileName 
      STOP
      END if


 !     read network data from text file: iras.inp (INP file is opened and closed once here, performance loss- Evgenii 101115)	
	CALL read_network_data (success) 
      if (.not.success) then
        WRITE(*,*)'Network data reading failed from file:',sysFilename
        STOP
      END if
	!Read INP file
	call readINPFile(success)
	if ( .not.success) then																			
			WRITE(*,*)'Error in reading INP file'
			STOP
      END if
	!Read Pol file
	call readPolFile(success)
	if ( .not.success) then																			
			WRITE(*,*)'Error in reading POL file'
			STOP
      END if
	
	!Read changes file
	call readChangesFile(success)
	if ( .not.success) then																			
			WRITE(*,*)'Error in reading CHG file'
			STOP
      END if
	CALL allocate_vars (success) 
      if (.not.success) then
        WRITE(*,*)'Variable allocation failed'
        STOP
      END if

	!following needs network data
      
	!Evgenii commented out call to ReadGageYearUnits, no longer needed in IRAS 2010
	!Now readGageUnits called directly (instead of from ReadGageYearUnits)
	!call ReadGageYearUnits(success)
	
	!read IRAS.gag gauge conversion factors and determine which gauges have flowfactors if any (iras.gag)
	call readGageUnits(success)   
      if (.not.success) then
        WRITE(*,*)'Error in reading gage file:',FlowFilename
        STOP
      END if

	!Read the iras.flw file into memory
	call readFlowFile(success)
	if ( .not.success) then																			
			WRITE(*,*)'Error in reading flow file'
			STOP
      END if

	!Read demand file if there demand time-series
	if (tDemSeries>0) then
		call readDemFile(success)
		if ( .not.success) then																			
			WRITE(*,*)'Error in reading demand file'
			STOP
		END if
	end if
	!Make directory Out
	!Kang remove 20100629 because the following statement is not supported on linux
	!makeDR=MAKEDIRQQ('Out') 


	
	!Read flow factors for run from guage node specific factor files (gaugename.fac)
	call ReadFlwFactors(success)
	if ( .not.success) then																			
		WRITE(*,*)'Error in reading factor file'
		STOP
      END if

	
	!Kang add 20100629
	CALL CPU_TIME (time_end_readfile)
	WRITE(*,*)'Time of reading param files was ',
     &         time_end_readfile - time_begin,'seconds' 
	
      success = .false.
	sysstat(nleap)=0      
	!Leap years added 101801 by Evgenii
	!Initialize total leap years in simulation system variable
	!find which years are leap yeras
	ly=sysstat(yrstart)

	!Do for all years
	do i=1,(SYSSTAT(YREND)-SYSSTAT(YRSTART)+1)
		!Identify leap years
		LeapYear(i)=.false.
		IF(MOD(ly,100)/=0.AND.MOD(ly,4)==0)then  
     			LeapYear(i)=.true.
			sysstat(nleap)=sysstat(nleap)+1	
		else if	(MOD(ly,400)==0)then
			LeapYear(i)=.true.
			sysstat(nleap)=sysstat(nleap)+1	
		endif
		!advance to next year
		ly=ly+1
	enddo 
		!Set total days in timestep, Evgenii changed this to NPER, in previous IRAS hardwired to one 1
	DAYSPRPRD =  SYSSTAT(NPER) 
		
	!STEPPRPRD: number of sub-time steps in time-step
	!if it is less than nMinSteps, then set to nMinSteps
 	STEPPRPRD = max(SYSSTAT(NSUBTT),nMinSteps) !Evgenii replaced MAX(nMinSteps,10.0) 091001 by current code
		
	!Number of days per sub-time step
	DAYPERTS = DAYSPRPRD/STEPPRPRD	!Number of days per sub-time step
	
	!Find total number of days in simulation period (including leap years)
	SYSSTAT(NUMDAYS)=(SYSSTAT(YREND)-SYSSTAT(YRSTART)+1)*365
     &	+SYSSTAT(nLeap) !Total number of days in simulation period
	
	!read default evaporation rates for surface nodes from iras.inp (subroutine found in Read_sim_Data.for)
	call readSysEvaporation(success)
	
	return
	end
!	 *************************************************************************************
	SUBROUTINE SIMSYS (success, csustain, cflow, cdemand)
! revised from old SIMSYS in order to read flow data from a text file
!       flow data file: IRAS.GAG
! It calls the following subroutines:
!  * ReadFlwFactors
!  * InitPolicy
!  * readSysEvaporation
!  * read_year_policies
!  * InitVariables
!  * Read_Simulation_Data
!  * InitLinkRouting
!  * GetMonthDay
!  * FLWSIM()
!  * DayOutputText()
!  *DayOutputNode
!  *DayOutputLink

 	USE vars
!        use output_iras_netcdf, only: nodes_netcdf_filename,
!     &       links_netcdf_filename
      implicit none
	
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
C
C    Arguments:
!  Input:
      !COMMON: CHARACTER*80 PolicyFileName,FlowFileName, SysFileName
!  Output:
      LOGICAL*1     success
C    Local variables:
      !Kang add 20100629
	REAL:: time_begin, time_end, time_temp
	INTEGER*4 iDatafile, iFlowfile
	CHARACTER*256 aLine
	CHARACTER*20  aVar
      INTEGER*4 sim,i,nn,s
      INTEGER*4 iYear, iDay,idayY,iMonth,iDayM,ileap,ln    
	character(len=30) NodesFileName,LinksFileName
	character(len=30)xrun
       character(len=30) DemandsFileName
       integer*4 ioutDemands
      character*1 csustain, cdemand
      character*2 cflow 
C
!  CALL:
      LOGICAL CountFileLines   !Function !Kang modify for improving performance
      INTEGER CountColumns     !Function !Kang modify for improving performance
!------------------------------------------------------------------------ 

      iOutFile = 20
	ioutNodes = 243 !Evgenii added ioutNodes for nodes output 090721 ! Anthony changed to 243 from 22 due to apparent clash meaning node output blank
	iOutLinks = 24 !Evgenii added iOutlinks for links output 090721

	!Debuging output file
	!open (UNIT=iOutFile, FILE='IRASTMP.OUT', STATUS='Replace')
      !WRITE(iOutFile,*)"IRAS OUTPUT"

	!Loop for total number of runs (flow factor runs), Evgenii 0911

		sysstat(run)=sysstat(run)+1
		SysStat(SimRec)=0
		xrun=' '
		write(xrun,*)sysstat(run)
		xrun=adjustl(xrun)
		
		!Add run number to output file name (gaugename.out)
		LinksFileName='output/links_run'//'-'//cdemand//'-'
     &     //csustain//'-'//trim(cflow)//'.out'
		NodesFileName='output/nodes_run'//'-'//cdemand//'-'
     &     //csustain//'-'//trim(cflow)//'.out'
		
!                nodes_netcdf_filename='nodes_run'//trim(xrun)//'.nc' ! - name of output netcdf file FF 07/2012               
!                links_netcdf_filename='links_run'//trim(xrun)//'.nc' ! 

		!Open links and nodes out file for new run, Evgenii 090721 
		!OPEN(UNIT=iOutLinks,FILE=trim(LinksFileName),STATUS='replace') 
		!OPEN(UNIT=iOutNodes,FILE=trim(NodesFileName),STATUS='replace')
                
		!EVGENII - Activate four lines below (and calls for DayOutputNodeTS and DayOutputLinkTS
!	towards end of FLWSIM in flwsim.for) for diagnostic sub-time step OUTPUT
	 ! ioutnodests=32
	!  ioutlinksts=33                                                
	!  OPEN(UNIT = iOutnodesTS, FILE ="NodesTS.out", STATUS='replace') 
	 ! OPEN(UNIT = iOutlinksTS, FILE ='LinksTS.out', STATUS='replace') 
      
!      ioutDemands=9899
!      DemandsFileName='demands_run'//trim(xrun)//'.out'  ! added by FF 07/2012
!      open(unit=ioutdemands,file=trim(DemandsFileName),
!     &      form="unformatted", status='replace') 
		
!		!For date debugging
		!OPEN(UNIT=200,FILE=trim('debug.txt'),STATUS='replace')
		
	  !Kang add 20100629
	  !CALL CPU_TIME (time_begin)
	  	
		!Initialize policy variables (subroutine found in initsys.for)
		call InitPolicy
                
		!read the first year policies (subroutine found in initsys.for)
      	call read_year_policies(PolicyGrpID(1),success) 

      	if (.not. success)GOTO 9999

		!PolicygroupID was specified in iras.def, it describes the year each policy is applicable 

		call InitVariables
	
		!read all simulation data for policyGrp(1) because of InitLinkRouting() (subroutine found in read_sim_data.for)
        !Kang modify 20100630
 		!call Read_Simulation_Data(SysFileName,PolicyGrpID(1),success) !read_sim_data.for
        call Read_Simulation_Data(PolicyGrpID(1),success) !read_sim_data.for
       
 
        call changeseq

	call SetSourceLinkNodes() !Evgenii added 110915 to search link sources
 		!Kang add 20100629
 		!CALL CPU_TIME ( time_end )
 		!time_total = time_end - time_begin
	
		call InitLinkRouting  !related to initial volume & number of reservoirs for routing

!***		Loop over number of runs and simulation YEARS ...   
!------------Code below heavily changed by Evgenii for variable time step 090706----------------	
                
	!   Initialize total leap years passed to 0
		ileap=0
		iday=0
	!	Initialize day in year (between 1 and 365, 366 for leap)
		sysstat(time)=0
		sysstat(sstep)=0
	
		!Initialize simulation with first year
		iyear=SYSSTAT(YRSTART)  
		SYSSTAT(YEAR) = iYear
		
		!Total years simulated
		SYSSTAT(SIM_YEAR)=1 

	    sysstat(ntimestep)=0
        
		!Do from day 1 to total number of days in simulation, in increments of time period
		do iDay = 1,SYSSTAT(numdays),SYSSTAT(NPER) 
        
			sysstat(nday)=iday
		    sysstat(ntimestep)=sysstat(ntimestep)+1
			nTS=sysstat(ntimestep)
			if(sysstat(ntimestep)==(sysstat(nrec)+1))then
				WRITE(*,*)'Simulation Ended, no more flow record'
				goto 301
			end if
			!change year when timestep brings a new year
			!For normal years
			IF (.not.leapyear(SYSSTAT(SIM_YEAR)).and.  
     &             (SysStat(Time)+SYSSTAT(NPER))>365)then
     				 call annualPerformance()
				 iyear=iyear+1  
				 SYSSTAT(YEAR)=iYear
				 !update simulation year
				 SYSSTAT(SIM_YEAR) =SYSSTAT(YEAR)-SYSSTAT(YRSTART) + 1 
				 
				 call InitPolicy 			     
			     call read_year_policies(  !Read policies for new year
     &                       PolicyGrpID(SysStat(Sim_Year)),success)
				 
				 
			     call ChangeNonPolicyData		
			     call ChangePolicyData		 
				 !Reset annual licenses
				 do ln=1,links
				    YearQLN(ln)=0.0
				 end do
				 do nn=1,tnodes
				    do s=1,MXSUPLY
				        AnnualOut(s,nn)=0.0
				    end do
				 end do						 
			!leap year
			else if (leapyear(SYSSTAT(SIM_YEAR)).and.  
     &             (SysStat(Time)+SYSSTAT(NPER))>366)then	
				 call annualPerformance()
				 iyear=iyear+1  
				 SYSSTAT(YEAR)=iYear 
				 !add to total years simulated
				 ileap=ileap+1
				 !update simulation year
				 SYSSTAT(SIM_YEAR) =SYSSTAT(YEAR)-SYSSTAT(YRSTART) + 1  
				 call InitPolicy 
			     call read_year_policies(  !Read policies for new year
     &                       PolicyGrpID(SysStat(Sim_Year)),success)
			     call ChangeNonPolicyData
			     call ChangePolicyData	
				 do ln=1,links
				    YearQLN(ln)=0.0
				 end do
				 !Reset annual licenses
				 do nn=1,tnodes
				    do s=1,MXSUPLY
				        AnnualOut(s,nn)=0.0
				    end do
				 end do				 
			Endif			
			
			!Set day in year (between 1 and 365)
			SysStat(Time) = iDay-(SYSSTAT(SIM_YEAR)-1)*365
			!Adjust for leap years simulated
			SysStat(Time)=SysStat(time)-ileap
			idayY=SysStat(Time)   

			!Get date			
			 !For leap years  
			if (leapyear(SYSSTAT(SIM_YEAR)))then
				call GetMonthDayleap(SysStat(Time), SysStat(Month), 
     &				SysStat(Day))  
			 !For normal years
			else
      			call GetMonthDay(SysStat(Time), SysStat(Month),
     &			  SysStat(Day))
			end if
			iMonth = SysStat(Month); iDayM = SysStat(Day)
			!for date debugging
			!WRITE(200,*)iYEAR,iMonth,iDayM,sysstat(time),iday,ileap	

C		  FLWSIM performs the fMonthlow simulation over the
C		  specified number of simulation time steps: one step = n days
		  !Simulation data are related to simulation day.
		  !So, data such as allocation, rating are read if needed in FLWSIM()
		  
		  !Advance by one time step
		  SysStat(SimRec)=SysStat(SimRec)+1
		  CALL FLWSIM(PolicyGrpID(SysStat(Sim_Year)),
     &	           iYear,idayY,success)         

			 IF(.not.success) GOTO 9999 
		     !Kang remove because I am afraid it maybe have side effect on performance
		     !WRITE(*,6010)iYEAR,iMonth,iDayM

!            **Write simulation results for a within-year period
!		   Evgenii- If statement below is only for IRAS diagnostics, not needed for  simulatins 
!		   if (SYSSTAT(SIM_YEAR)==1)    !only output results in tmp out file for first year 
!    &		Call DayOutputText(iOutFile,iYEAR,iMonth,iDayM, NodSEQ)
         
			!Evgenii - Write time-step results to Nodes outfile 090721
             ! call DayOutputNode(ioutNodes,iDay) 
			!Evgenii - Write time-step results to New Links outfile 090721
			!call DayOutputLink(ioutLinks,iDay)
            
            ! write time-step results to demands outfile - FF 120718
            !call output_ff_writetobinary_dmd(ioutDemands) 
200		 CONTINUE
C		  End of within-year loop
C
300		ENDDO

	!close (iOutFile); !close (ioutNodes); close (ioutLinks)			
301   success = .true.
9999  continue      
      return
C
6010  FORMAT('+','Year: ',I5,', Month: ',I2, ', Day: ',I2)

6011  FORMAT('+','Day: ',I5,',Year: ',I5)
999	WRITE(*,*)'Error in reading gage flow!'
!       write (*,*)'Press any key to exit... '
!       read (*,*)
      END
C
C


!************************************************************************
      subroutine GetMonthDay(nDays, iMonth, iDay)
!  INPUT
      INTEGER*4 nDays
!  OUTPUT
      INTEGER*4 iMonth, iDay
!------------------------------------------------------------------------
      iMonth = 1
      iDay = 1
     
!
	
	SELECT CASE (nDays)
          CASE (1:31)
            iMonth = 1; iDay = nDays
          CASE (32:31+28)
            iMonth = 2; iDay = nDays - 31
          CASE (60:59+31)
            iMonth = 3; iDay = nDays - 59
          CASE (91:90+30)
            iMonth = 4; iDay = nDays - 90
          CASE (121:120+31)
            iMonth = 5; iDay = nDays - 120
          CASE (152:151+30)
            iMonth = 6; iDay = nDays - 151
          CASE (182:181+31)
            iMonth = 7; iDay = nDays - 181
          CASE (213:212+31)
            iMonth = 8; iDay = nDays - 212
          CASE (244:243+30)
            iMonth = 9; iDay = nDays - 243
          CASE (274:273+31)
            iMonth = 10; iDay = nDays - 273
          CASE (305:304+30)
            iMonth = 11; iDay = nDays - 304
          CASE (335:365)
            iMonth = 12; iDay = nDays - 334
      END SELECT
      end subroutine
		 

!************************************************************************
      subroutine GetMonthDayLeap(nDays, iMonth, iDay)
!Evgenii modified GetMonthDay for leap years 101809
!For Leap Years
!  INPUT
      INTEGER*4 nDays
!  OUTPUT
      INTEGER*4 iMonth, iDay
!------------------------------------------------------------------------
      iMonth = 1
      iDay = 1
     
!
	
	SELECT CASE (nDays)
          CASE (1:31)
            iMonth = 1; iDay = nDays
          CASE (32:60) !31+29)
            iMonth = 2; iDay = nDays - 31
          CASE (61:91) !59+31)
            iMonth = 3; iDay = nDays - 60
          CASE (92:121) !90+30)
            iMonth = 4; iDay = nDays - 91
          CASE (122:152) !120+31)
            iMonth = 5; iDay = nDays - 121
          CASE (153:182) !151+30)
            iMonth = 6; iDay = nDays - 152
          CASE (183:213) !181+31)
            iMonth = 7; iDay = nDays - 182
          CASE (214:244) !212+31)
            iMonth = 8; iDay = nDays - 213
          CASE (245:274) !243+30)
            iMonth = 9; iDay = nDays - 244
          CASE (275:305) !273+31)
            iMonth = 10; iDay = nDays - 274
          CASE (306:335) !304+30)
            iMonth = 11; iDay = nDays - 304
          CASE (336:366)
            iMonth = 12; iDay = nDays - 335
      END SELECT
      end subroutine
		 
		 
