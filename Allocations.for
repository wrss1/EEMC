!Copyright (c) 2009, 2010 by University College London, Cornell University
!Authors:
!Evgenii Matrosov (evgenii.matrosov@ucl.ac.uk), Julien Harou (j.harou@ucl.ac.uk), 
!This program is free software under the General Public Licence, GPL (>=v2)
!Read the 'GPL License.txt' file distributed with this source code for a full license statement.
	subroutine allocate_vars(success)
	USE vars
	!Created by Evgenii 1006
	!
	
	IMPLICIT NONE
	INCLUDE 'IRAS_SYS.INC'
	INCLUDE 'NODE.INC'
	INCLUDE 'LINK.INC'
	!  INPUT

	!  OUTPUT

	!  LOCAL
	logical:: success
	real:: YearSD
	INTEGER*4 LN, NN,ny,i
 !	-----------------------------------------------------------
	if (.not. allocated(YearFailEvent) )
     &  Allocate (YearFailEvent(sysstat(nyear),tnodes,thres_pts_max))
	do ny=1,sysstat(nyear)
	  do nn = 1, tnodes 
		  do i=1,thres_pts_max
			 YearFailEvent(ny,nn,i)=.false. !Made this a logical variable Evgenii
		  end do
	  end do
	end do

	end
!******************************************************************
	subroutine de_allocate_vars()
	USE vars
	!Created by Evgenii 101115
	!	
	IMPLICIT NONE
	INCLUDE 'IRAS_SYS.INC'
	INCLUDE 'NODE.INC'
	INCLUDE 'LINK.INC'

		
      !Deallocate memory of some dynamic arrays
      IF(ALLOCATED(SysFileContent)) DEALLOCATE(SysFileContent)
      IF(ALLOCATED(TargetArray))  DEALLOCATE(TargetArray)
      IF(ALLOCATED(EvapArray))  DEALLOCATE(EvapArray)
      IF(ALLOCATED(RoutingArray)) DEALLOCATE(RoutingArray)
      IF(ALLOCATED(AlloArray))  DEALLOCATE(AlloArray)
      IF(ALLOCATED(SourceArray))  DEALLOCATE(SourceArray)
      IF(ALLOCATED(RatingArray))  DEALLOCATE(RatingArray)
      IF(ALLOCATED(CrossArray))  DEALLOCATE(CrossArray)
      IF(ALLOCATED(PerfArray))  DEALLOCATE(PerfArray)
      IF(ALLOCATED(CostArray))  DEALLOCATE(CostArray)
      IF(ALLOCATED(PowerArray))  DEALLOCATE(PowerArray)
      IF(ALLOCATED(PumpArray))  DEALLOCATE(PumpArray)
      IF(ALLOCATED(RuleArray))  DEALLOCATE(RuleArray)
	IF(ALLOCATED(DemRedArray))  DEALLOCATE(DemRedArray)	 !Evgenii 100719
      IF(ALLOCATED(PolicyFileContent)) DEALLOCATE(PolicyFileContent)
	!Evgenii's variables
	IF(ALLOCATED(flowdata)) DEALLOCATE(flowdata)
	IF(ALLOCATED(flwfactmax)) DEALLOCATE(flwfactmax)
	IF(ALLOCATED(YearFailEvent)) DEALLOCATE(YearFailEvent)
	IF(ALLOCATED(DemName)) DEALLOCATE(DemName)
	IF(ALLOCATED(DemTimeSeries)) DEALLOCATE(DemTimeSeries)
	IF(ALLOCATED(DemIDObs)) DEALLOCATE(DemIDObs)
	IF(ALLOCATED(DemData)) DEALLOCATE(DemData)
	IF(ALLOCATED(flowfactor)) DEALLOCATE(flowfactor)
	IF(ALLOCATED(BalanceArray)) DEALLOCATE(BalanceArray)
      
	return
	end