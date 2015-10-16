!System dimension
!Copyright (c) 2008, 2009 by University College London, Cornell University
!Authors:
!Daniel P. Loucks (dpl3@cornell.edu), Marshall Taylor, Peter French, 
!Evgenii Matrosov evgenii.matrosov@ucl.ac.uk), Julien Harou (j.harou@ucl.ac.uk), HC Zhou
!This program is free software under the General Public Licence, GPL (>=v2)
!Read the 'GPL License.txt' file distributed with this source code for a full license statement.
			
	  CHARACTER*80 PolicyFileName,FlowFileName,FlwDataFilename
	  Integer*4 iOutFile,ioutNodes, ioutlinks
        CHARACTER*80 SysFileName,UnitsFileName,DefFile !Evgenii added UnitsFileName,flwdatafilename
	  COMMON /FileName_1/PolicyFileName,FlowFileName,DefFile,
     &                     SysFileName,UnitsFileName,FlwDataFilename

	  			 
!       Max number of policies for a type of data for a given policy group
        INTEGER*4  MaxPolicies, MaxPolicyTypes
 !      PARAMETER (MaxPolicies = 5)      !for debugging
        PARAMETER (MaxPolicies = 366)      !for final
        PARAMETER (MaxPolicyTypes = 11)  !for gage, allocation, rating, etc.
!zhou
C
C     LNKMAX - Maximum number of links
C     NODMAX -       "           nodes
C     TMAX   -       "           time periods per year
C     YRMAX  -       "           years for which dynamic functions can be defined
!	nMinSteps - Number of sub-time steps added by Evgenii
      INTEGER*4      LNKMAX, NODMAX, TMAX, YRMAX, nMinSteps
!      for easily debugging example 'Sys_a'
!      PARAMETER     (LNKMAX = 5, NODMAX = 5,
!     2               YRMAX  = 2, nMinSteps = 10)
      PARAMETER     (LNKMAX = 1500, NODMAX = 1500,
     2               YRMAX  = 300, nMinSteps = 10)
C
C     GAGMAX - Maximum number of gauge sites,
      INTEGER*4      GAGMAX
      PARAMETER     (GAGMAX  = 1000)
!     PARAMETER     (GAGMAX  = 3)               !for debugging
C
C     Maximum number of flow-routing reservoirs in a link.    
      INTEGER*4  MAX_NRTRES
      PARAMETER (MAX_NRTRES=10)

C     NRMAX  - reservoirs and release sites in a group,
C     ILKMAX - maximum number of incoming
C     OLKMAX - maximum number of outgoing links
C     IAGMAX - maximum number of arguments on each function needing interpolation,
C              and maximum number of reference nodes where balancing functions
C              could be defined.
      INTEGER*4      NRMAX, ILKMAX, OLKMAX, IAGMAX
      PARAMETER     (NRMAX  =  9,
     2               ILKMAX =  100, OLKMAX =  100, IAGMAX =  25 )
	 Integer*4 UserIDMAX			!Maximum user set ID Evgenii 121120
	 Parameter (UserIDMAX=9999999)
C
C
C     RNMAX - Max # reservoir rule curve sites
C     GRMAX - Max # reservoirs in a group                               
C     ZNMAX - Max # reservoir rule curve zones
      INTEGER*4  RNMAX,  GRMAX, ZNMAX
      PARAMETER (RNMAX=60, GRMAX=15, ZNMAX=8) !Evgenii 091118 changed RNMAX to 5 instead of 2
C
C
C     Maximum number of sites which can supply water to a demand node
      INTEGER*4  MXSUPLY
      PARAMETER (MXSUPLY = 100 )
C
	logical*1 ::    LeapYear(YRMAX)
C+IRIS+IRIS+IRIS+IRIS+IRIS+IRIS+IRIS+IRIS+IRIS+IRIS+IRIS+IRIS+IRIS+IRIS+IRIS
C
C      File: SYS_STAT.INC
C
!Zhou moved from "CAP_PARM.INC"
        INTEGER*4   SUCCES, FAIL
        PARAMETER  (SUCCES=0, FAIL=1)
!Zhou
        !for PolicyType
        INTEGER*4 DemRed0, Allocation0, Rate0, Power0, Pump0      !for Node
        INTEGER*4 Target0,Source0,Rule0, Balance0, Evaporation0 !for node
        INTEGER*4 Routing0, Transfer0,Performance0,Cost0,GW0    !for link also Rate0: same as Node !Evgenii added Performance and Cost 100617
        !for nodes
        PARAMETER(DemRed0       = 1)	!Evgenii changed Gage0  to DemRed0	100720
        PARAMETER(Allocation0   = 2)
        PARAMETER(Rate0         = 3)
        PARAMETER(Power0        = 4)
        PARAMETER(Pump0         = 5)
        PARAMETER(Target0       = 6)
        PARAMETER(Source0       = 7)
        PARAMETER(Rule0         = 8)
        PARAMETER(Balance0      = 9)
        PARAMETER(Evaporation0  = 10)
	    PARAMETER(Performance0  = 11)
        !for links
        PARAMETER(Routing0      = 1)
        PARAMETER(Transfer0     = 2)
	    PARAMETER(Cost0         = 4)

  
      
        INTEGER*4 SIM_YEAR,YEAR,TIME, Month, Day, BYearGage, EYearGage
        INTEGER*4 YRSTART,YREND,nGages1,nGages2, RepeatMethod,nRuns
        INTEGER*4  NGAUGE1, NSUBTT, NPER, NREP,NUMDAYS,NRec,SimRec,run
	  INTEGER*4 nleap,nday,sstep, ntimestep, nyear,nDems
        PARAMETER (YEAR=1)      !"real" current year (e.g. 1980,...)
        PARAMETER (SIM_YEAR=2)  !# of the year being simulated (1,2,...)
        PARAMETER (Month = 3)
        PARAMETER (Day = 4)	  !Day in month
        PARAMETER (TIME=5)      !number of days from year beginning
        PARAMETER (YRSTART=6)  !Actual start (year) of flow simulation
        PARAMETER (YREND=7)    !Last year of simulation
        PARAMETER (NGAUGE1=8) !# gauges in the system
        PARAMETER (NUMDAYS=9)  !Evgenii 090706 added this, number of days in simulation
        PARAMETER (NSUBTT=10)  !# of siumulation time steps Evgenii added 091026
        PARAMETER (NRec=11)     !number of flow records
        PARAMETER (SimRec=12)     !current flow record
        PARAMETER (NPER=13)     !time-perid length in days, evgenii
	  PARAMETER (NREP=14)     !Number of replicates
	  PARAMETER (nRuns=15)    !Number of model runs (flow factor sets) added by Evgenii 091115
	  PARAMETER (Run=16)    !Current model run  added by Evgenii 091115
	  PARAMETER (nleap=17) !number of leap years in simulation
	  PARAMETER (nday =18)     !day of simulation
	  PARAMETER (sstep =19)    !Sub-time step
	  PARAMETER (ntimestep=20)    !current timestep
	  PARAMETER (nyear=21)	   !Total years in simulation
	  PARAMETER (nDems =22)  !# gauges in gage file, Evgeni took out nGages2 (No added flow in IRAS 2010) 091216
       !PARAMETER (RepeatMethod = 9) !1-use last record; 2-return to beginning Evgenii took out RepeatMethod 091216
!	  PARAMETER (NGAUGE2=15) !# gauges of added flow in the system Evgenii took out added gauges
!       PARAMETER (BYearGage=12)!beginning year in gage, file
!       PARAMETER (EYearGage=13)!end year in gage, file
!       PARAMETER (nGages1 =7)  !# gauges in gage file

	  
	  INTEGER*4   MAX_STATS
        PARAMETER   (MAX_STATS=22)
        INTEGER*4   SYSSTAT(MAX_STATS)
	!Integer   step !Evgenii added this 090807	
	
        !for converting user units to internal units
        !internal units: length: m; area: s.m; volume: 10^6 cubic meters; flow: 10^6 CMS/day
        INTEGER*4 ULen, UArea, UVol, Uflow, ULoss, UK, UPower
        INTEGER*4 UTime
        PARAMETER (ULen  = 1)
        PARAMETER (UArea = 2)
        PARAMETER (UVol  = 3)
        PARAMETER (UFlow = 4)
        PARAMETER (ULoss = 5)       !to m/day for node andlink evap and for some link loss methods
        PARAMETER (UPower= 6)     
        PARAMETER (UK    = 7)       !hydraulic conductivity
        PARAMETER (Utime = 8) 
        
! Evgenii 220409 changed UK parameter from 5 to 7, eliminated repeat of value 5
        INTEGER*4 maxUnits
        PARAMETER (maxUnits = 15)
        real*4    NodeUserUnit(maxUnits)
        real*4    LinkUserUnit(maxUnits)
        real*4    GageUserUnit(GAGMAX) !Evgenii took out 2nd dimension of GageUserUnit because added flow taken out, used to be GageUserUnit(GAGMAX,2) 100108
        INTEGER*4 GageIDObs(GagMax)               !flow gage file, Evgenii took out 2nd dimension of GageUserUnit because added flow taken out, used to be GageIDObs(GagMax,2) 100108
        logical :: GageFact(gagmax) !Evgenii 091212 
        CHARACTER*100 GageName(GagMax) !Evgenii 091212 moved GageName to global variables from local variable in subroutine readGageUnits
	  INTEGER*4 outunit(gagmax)
C     Factors to convert user flow units to internal (volume/simulation
C     time step) units and back. The values are computed in FLWSIM.

!	For performance measures
	INTEGER*4 thres_pts_max    
	parameter (thres_pts_max=10)


        REAL*4      DAYSPRPRD   ! Days per period
        REAL*4      STEPPRPRD   ! Steps per period
        REAL*4      DAYPERTS    ! Days per simulation time step
	  Real::      nTS		  !	Current TS Evgenii
	  integer NINV, flwdata !Number of within year periods, flowfile id	 Evgenii
		!integer NREPMX !Number of replicates evgenii
	    integer*4 ioutnodests, ioutlinksts
        logical*1 userseq, demprop
      
	!For changes
	logical*1 YearPolChg(YRMAX,2)
	logical*1 YearNodePolicyChg(2,YRMAX,MaxPolicyTypes,nodmax)    

        COMMON /SYSSTAT_1/SYSSTAT,DAYPERTS,DAYSPRPRD,STEPPRPRD,nTS
     &   ,userseq,demprop
        COMMON /Unit_1/ NodeUserUnit,LinkUserUnit,LeapYear
        COMMON /Unit_2/ GageIDObs,GageUserunit,GageName,GageFact,outunit
     &    ,ioutnodests,ioutlinksts,iOutFile,ioutNodes, ioutlinks
	    COMMON /Change/ YearPolChg,YearNodePolicyChg

