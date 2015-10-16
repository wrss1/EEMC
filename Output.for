!!************************************************************************
	subroutine DayOutputNode(iOutfileBin,iDay)
	!Created by Evgenii 0907
	!This subroutine prints out the time step results to a file found in the /out directory.
!    modified by FF 08/2012 - added call to output_node_netcdf and output_link_netcdf
!	use output_iras_netcdf, only: output_node_netcdf
      IMPLICIT NONE
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      INTEGER*4:: iOutFileBin,iDay,NN_SEQ(NODMAX)
!  Local
      INTEGER*4 i,j,NN
      INTEGER*4 j0, j1,k,k1,k2,k3,k4,k5
	INTEGER*4 tevap_nodes, t_storage, t_gage,tseep
	INTEGER*4 storage_id(nodmax),gage_id(nodmax),t_cons,ln
	INTEGER*4 cons_id(NODMAX),evap_id(NODMAX),seep_id(NODMAX)
	real:: storage_beg(NODMAX), storage_end(NODMAX),gage(nodmax)
	real:: TOT_ALLOC,convbqln
	real:: consump(nodmax),tevap(NODMAX),seepage(NODMAX)
	real :: unallocated(NODMAX),EndNodeSto(NODMAX),BegNodeSto(NODMAX)
	character(45) storage_Name(NODMAX),gage_name(nodmax),unallC(nodmax)  !Name of storage
	character(45) cons_name(NODMAX),evap_name(NODMAX)
	character(45) flowin(NODMAX),outflow(nodmax),gagein(NODMAX),
     &	begsto(NODMAX),endsto(NODMAX),cons(NODMAX),evap(NODMAX),
     &	seep_name(NODMAX),seep(NODMAX)
!------------------------------------------------------------------------
	t_cons=0; tevap_nodes=0; t_storage=0; t_gage=0; tseep=0
	
	j=0; k=0; k1=0; k2=0; k3=0; k4=0; k5=0	
	
	do i=1,tnodes
	  !Convert units from Mm3/time-step to what guage input units were
        
	  call UndoUnitConversion(1,uFlow,INFLOW(i))
        INFLOW(i)=INFLOW(i)/SYSSTAT(NPER)
        call UndoUnitConversion(1,uFlow,TEVAPN(i))
        TEVAPN(i)=TEVAPN(i)/SYSSTAT(NPER)
        call UndoUnitConversion(1,uFlow,TSEEPL(i))
        TSEEPL(i)=TSEEPL(i)/SYSSTAT(NPER)
        call UndoUnitConversion(1,uFlow,TOTREL(i))
        TOTREL(i)=TOTREL(i)/SYSSTAT(NPER)
        call UndoUnitConversion(1,uFlow,CONSUMPTION(i))
        CONSUMPTION(i)=CONSUMPTION(i)/SYSSTAT(NPER)
        BegNodeSto(i)=BSTO(i)
        EndNodeSto(i)=ESTO(i)
        call UndoUnitConversion(1,uVol,BegNodeSto(i))
        call UndoUnitConversion(1,uVol,EndNodeSto(i))
      end do
      do i=1,tnodes !make arrays for the output for each category
		NN = NodSEQ(I)
		IF (NN .EQ. 0) GO TO 100
		!Makes arrays for output of beginning storage (only for storage nodes)
		if (ResvNode(i).or. NatLak(i) .or. GWNODE(i)) then
			t_storage=t_storage+1
			k=k+1
			storage_beg(k)=BegNodeSto(i)  !Beg. storage array
			storage_Name(k)=TRIM(NNAME(i)) !Name for begstorage nodes
			storage_id(k)=nn   !node id
			storage_end(k)=EndNodeSto(i)
		endif

		if (GageNF(i) == .true.) then
			t_gage=t_gage+1
			k2=k2+1
			Gage_Name(k2)=TRIM(NNAME(i)) 
			Gage_id(k2)=nn
			gage(k2)=NQinn(i) 
		endif
		
		if (cons_node(i)== .true.)then
			t_cons=t_cons+1
			k3=K3+1
			cons_name(k3)=TRIM(NNAME(i)) 
			cons_id(k3)=nn
			consump(k3)=CONSUMPTION(i) 
		endif

		if (ResvNode(i).or. NatLak(i) .or. GWNODE(i))then
			tevap_nodes=tevap_nodes+1
			k4=K4+1
			evap_name(k4)=TRIM(NNAME(i)) 
			evap_id(k4)=nn
			tevap(k4)=TEvapn(i) 
		endif
	
		if (seep_node(i)== .true.)then !seep_node defined read_sim_data line 723
			tseep=tseep+1
			k5=k5+1
			seep_name(k5)=TRIM(NNAME(i)) 
			seep_id(k5)=nn
			seepage(k5)=TSeepl(i) 
		endif
		TOT_ALLOC = 0.0
          DO J = 1,TOTOUT(NN) !Find non allocated flow
			 LN = OUTLNK(NN,J)
			 convBQLN=BQLN(ln)
			 call UndoUnitConversion(1,uFlow,convbqln)
               convbqln=convbqln/SYSSTAT(NPER)
             TOT_ALLOC = TOT_ALLOC + convBQLN
          ENDDO
          
          TOT_ALLOC = TOT_ALLOC - TOTREL(NN)
          unallocated(nn)=ABS(TOT_ALLOC)
          !call UndoUnitConversion(1,uFlow,unallocated(nn))
	enddo	
!      call output_node_netcdf(BegNodeSto, EndNodeSto, unallocated) ! added by FF 07/2012
!      return
    
	if (sysstat(nday)==1) then 
		do i=1,tnodes
		flowin(i)='Inflow';outflow(i)='Outflow';unallC(i)='unAll_flow'
		enddo

		!do i=1,t_gage !t_gage defined in initsys.for
		!gagein(i)='Flow_Input'
		!enddo

		do i=1,t_storage !t_storage defined in initsys.for
		begsto(i)='Beg_Vol'; endsto(i)='End_Vol'
		enddo

		
		do i=1,t_cons
		cons(i)='Cons'
		enddo
		
		do i=1,tevap_nodes
		evap(i)='Evap'
		enddo

		do i=1,tseep
		seep(i)='Seepage'
		enddo
	
		write(ioutfilebin,*)'Run number = ', sysstat(run)
		write(ioutfilebin,FMT=11,advance='no')     !Write the output category                 
     &	(begsto(i),i=1,t_storage),(endsto(i),i=1,t_storage),
     &	(cons(i),i=1,t_cons)
     &	,(evap(i),i=1,tevap_nodes),(seep(i),i=1,tseep),  
     &	(flowin(i),outflow(i),i=1,tnodes)                             !(endsto(i),i=1,t_storage)
     & !   (unallC(i),i=1,tnodes) !(gagein(i),i=1,t_gage),
     			
		
		write(ioutfilebin,*)
		
		write(ioutfilebin,FMT=10,advance='no')'Day', !Write the name of node
     &   (storage_Name(i),i=1,t_storage),(storage_Name(i),i=1,t_storage)
     &	,(cons_name(i),i=1,t_cons)	
     &	,(evap_name(i),i=1,tevap_nodes),(seep_name(i),i=1,tseep),
     &	(NName(i),NName(i),i=1,tnodes)					!(storage_Name(i),i=1,t_storage)
     &  !  (NName(i),i=1,tnodes)!(Gage_Name(i),i=1,t_gage)	
     			
	write(ioutfilebin,*)
	endif !ends if loop for day=1
	

	write(ioutfilebin,FMT=9,advance='no')iday, !Write the data for each day 
     &	(storage_beg(i),i=1,t_storage),(storage_end(i),i=1,t_storage)
     &	,(consump(i),i=1,t_cons)
     &	,(tevap(i),i=1,tevap_nodes),(seepage(i),i=1,tseep),
     &	(inflow(i),TOTREL(i),i=1,tnodes)
 !    &  (unallocated(i),i=1,tnodes) !(gage(i),i=1,t_gage)                      !(storage_end(i),i=1,t_storage),
     	
	write(ioutfilebin,*)
	
9     FORMAT(I6,10000F51.3)
10    FORMAT(A5,5X,10000A51)
11	FORMAT(8X, 10000A51)


 100   Continue
	 end subroutine
!***************************************************************************************

	subroutine DayOutputLink(iOutfileBin,iDay)
!
	!Created by Evgenii 0907
	
      IMPLICIT NONE
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      INTEGER*4:: iOutFileBin,iDay,NN_SEQ(Nodmax) 
!  Local
      INTEGER*4 i,j,NN
      INTEGER*4 j0, j1,k,k1,k2,k3,k4
	integer :: tloss,tend,tenergy,tgeo,troute
	real:: Loss(LNKMAX),endln(lnkmax),energyln(LNKMAX),volume(LNKMAX),  
     &		depth(LNKMAX),width(LNKMAX),velocity(LNKMAX),linkvol(LNKMAX) 	
	character(20) begflow(LNKMAX),endflow(LNKMAX),cEnd(LNKMAX)  !Name of storage
	character(20) Loss_Name(LNKMAX),CLoss(LNKMAX),end_name(LNKMAX),
     &		energy_name(LNKMAX), cenergy(LNKMAX),cdepth(LNKMAX),
     &        cwidth(LNKMAX),cvelocity(LNKMAX), cvolume(LNKMAX),
     &		GeoName(LNKMAX),Volume_Name(LNKMAX)
		
	
!------------------------------------------------------------------------
	tloss=0; tend=0; tenergy=0; tgeo=0; troute=0 !Resets the count troute=0 
	j=0; k=0; k1=0; k2=0; k3=0; k4=0

      DO i = 1,LINKS
        !Even if loss conversions used Uloss to convert in read_sim_data before, loss is in flow units so Uflow is used now)
        call UndoUnitConversion(2,uPower,ENERGY(i))
        call UndoUnitConversion(2,uFlow,BQLN(i))
        BQLN(i)=BQLN(i)/SYSSTAT(NPER)
        call UndoUnitConversion(2,uFlow,EQLN(i))
        EQLN(i)=EQLN(i)/SYSSTAT(NPER)
        !volumes need to be kept in iras units for future timestep
        linkvol(i)=TOT_LVOL(i)
        call UndoUnitConversion(2,uVol,linkvol(i))
        call UndoUnitConversion(2,uFlow,TLossL(i))
        TLossL(i)=TLossL(i)/SYSSTAT(NPER)
       ENDDO
	do i=1,links !make arrays for the output for each category
		
		!Makes arrays for output of loss output (only for links with loss)
		!Only if there is linkloss or if a rating table is present
		if (iflinkloss(i)==.true. .or. LossMethod(i) == 2) then   
			tloss=tloss+1
			k=k+1
			Loss(k)=TLossL(i)  !Beg. storage array
			Loss_Name(k)=TRIM(LName(i)) !Name for begstorage nodes
		endif
		
	!Makes arrays for output of endflow (only if endflow != begflow)
		if (iflinkloss(i)==.true. .or. LossMethod(i) == 2.or.
     &		L_Method(i)>= 1) then   !Only goes into loop if there is linkloss or routing
			tend=tend+1
			k1=k1+1
			Endln(k1)=EQLN(i)  !Beg. storage array
			End_Name(k1)=TRIM(LName(i)) !Name for begstorage nodes
		endif
		if (powerlink(i)==.true. .or. pumplink(i)== .true.) then !for energy
			tenergy=tenergy+1
			k2=k2+1
			energyln(k2)=ENERGY(i)  !Beg. storage array
			Energy_Name(k2)=TRIM(LName(i)) !Name for begstorage nodes
		endif
		!for routing links
		if (L_Method(i)>= 1) then
			troute=troute+1
			k3=k3+1
			volume(k3)=linkvol(i)
			Volume_Name(k3)=TRIM(LName(i))
		endif
		!for loss method 2
		if (LossMethod(i)== 2) then
				tgeo=tgeo+1
				k4=k4+1
				GeoName(k4)=TRIM(LName(i))
				depth(k4)=LastDepth(i)
				width(k4)=LastWidth(i)
				velocity(k4)=LastVelocity(i)
		endif	
	enddo	
      
!      call output_link_netcdf(  ! added by FF 07/2012
!     &      tloss,tend,tenergy,tgeo,troute,
!     &      Loss,endln,energyln,volume,depth,width,velocity,linkvol,
!     &      Loss_Name,end_name,energy_name,GeoName,Volume_Name)

!      return
      

	if (sysstat(nday)==1) then  !Initialize column names	
		do i=1,links
		begflow(i)='Beg_Flow'
		enddo

		do i=1,tLoss
		CLoss(i)='Loss'
		enddo

		do i=1,tend
		Cend(i)='End_Flow'
		enddo

		do i=1,tenergy
		Cenergy(i)='Energy'
		enddo
		
		do i=1,troute
		Cvolume(i)='Volume'
		enddo 

		do i=1,tgeo
		cdepth(i)='Depth'
		cwidth(i)='Width'
		cvelocity(i)='Velocity'
		enddo

		write(ioutfilebin,*)'Run number = ', SysStat(run)
		
		write(ioutfilebin,FMT=11,advance='no')     !Write the output category                 
     &	(begflow(i),i=1,links),(Cend(i),i=1,tend),
     &	(CLoss(i),i=1,tLoss),(cenergy(i),i=1,tenergy),
     &	(Cvolume(i),i=1,troute),(Cwidth(i),i=1,tgeo),
     &	(Cdepth(i),i=1,tgeo),(Cvelocity(i),i=1,tgeo)	
		write(ioutfilebin,*)
		
		write(ioutfilebin,FMT=10,advance='no')'Day', !Write the name of node
     &	(LName(i),i=1,links),(end_Name(i),i=1,tend),
     &	(Loss_Name(i),i=1,tLoss),(energy_name(i),i=1,tenergy),
     &	(Volume_Name(i),i=1,troute),(GeoName(i),i=1,tgeo),
     &	(GeoName(i),i=1,tgeo),(GeoName(i),i=1,tgeo)
	write(ioutfilebin,*)
	endif !ends if loop for day=1
	

	write(ioutfilebin,FMT=9,advance='no')iday,(BQLN(i),i=1,links), !Write the data for each day 
     &	(endln(i),i=1,tend),(loss(i),i=1,tloss), 
     &    (energyln(i),i=1,tenergy),(volume(i),i=1,troute),
     &	(width(i),i=1,tgeo),(depth(i),i=1,tgeo),(velocity(i),i=1,tgeo)
    
	
	write(ioutfilebin,*)


	
9     FORMAT(I6,10000F51.3)
10    FORMAT(A5,47X,10000A51)
11	FORMAT(8X, 10000A51)


100   Continue
	end subroutine


	subroutine DayOutputLinkTS(iDay,step)
!------------------------------------------------------------------------------------
	!Created by Evgenii 0907
	
      IMPLICIT NONE
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      INTEGER*4 iOutFileBin,iDay,NN_SEQ(Nodmax) 
!  Local
      INTEGER*4 i,j,NN,step
      INTEGER*4 j0, j1,k,k1,k2,k3,k4
	INTEGER*4 tloss,tend,tenergy,tgeo,troute
	real*4 Loss(LNKMAX),endln(lnkmax),energyln(LNKMAX),volume(LNKMAX),  
     &		depth(LNKMAX),width(LNKMAX),velocity(LNKMAX),EnergyConv(LNKMAX),
     &      BQLNconv(LNKMAX),TOT_LVOLconv(LNKMAX),TLossLconv(LNKMAX),
     &      EQLNconv(LNKMAX),DBQLNconv(LNKMAX) 	
	character(10) begflow(LNKMAX),endflow(LNKMAX),cEnd(LNKMAX)  !Name of storage
	character(10) Loss_Name(LNKMAX),CLoss(LNKMAX),end_name(LNKMAX),
     &		energy_name(LNKMAX), cenergy(LNKMAX),cdepth(LNKMAX),
     &        cwidth(LNKMAX),cvelocity(LNKMAX), cvolume(LNKMAX),
     &		GeoName(LNKMAX),Volume_Name(LNKMAX),begflowTS(LNKMAX)
		
	
!------------------------------------------------------------------------
	
	tloss=0; tend=0; tenergy=0; tgeo=0; troute=0 !Resets the count troute=0 


		
	j=0; k=0; k1=0; k2=0; k3=0; k4=0

	 DO i = 1,LINKS
        !Even if loss conversions used Uloss to convert in read_sim_data before, loss is in flow units so Uflow is used now)
        
        EnergyConv(i)=ENERGY(i)
        call UndoUnitConversion(2,uPower,EnergyConv(i))
        BQLNconv(i)=BQLN(i)
        call UndoUnitConversion(2,uFlow,BQLNconv(i))
        BQLNconv(i)=BQLNconv(i)/SYSSTAT(NPER)
        EQLNconv(i)=EQLN(i)
        call UndoUnitConversion(2,uFlow,EQLNconv(i))
        EQLNconv(i)=EQLNconv(i)/SYSSTAT(NPER)
        !volumes need to be kept in iras units for future timestep
        TOT_LVOLconv(i)=TOT_LVOL(i)
        call UndoUnitConversion(2,uVol,TOT_LVOLconv(i))
        TLossLconv=(tLossL(i))
        call UndoUnitConversion(2,uFlow,TLossLconv(i))
        TLossLconv(i)=TLossLconv(i)/SYSSTAT(NPER)
        DBQLNconv(i)=DBQLN(i)
        call UndoUnitConversion(2,uFlow,DBQLNconv(i))
        DBQLNconv(i)=DBQLNconv(i)/SYSSTAT(NPER)
        

       ENDDO

	
	
	do i=1,links !make arrays for the output for each category
		
		!Makes arrays for output of loss output (only for links with loss)
		!Only if there is linkloss or if a rating table is present
		if (iflinkloss(i)==.true. .or. LossMethod(i) == 2) then   
			tloss=tloss+1
			k=k+1
			Loss(k)=TLossLconv(i)  !Beg. storage array
			Loss_Name(k)=TRIM(LName(i)) !Name for begstorage nodes
		endif
		
	!Makes arrays for output of endflow (only if endflow != begflow)
		if (iflinkloss(i)==.true. .or. LossMethod(i) == 2.or.
     &		L_Method(i)>= 1) then   !Only goes into loop if there is linkloss or routing
			tend=tend+1
			k1=k1+1
			Endln(k1)=EQLNconv(i)  !Beg. storage array
			End_Name(k1)=TRIM(LName(i)) !Name for begstorage nodes
		endif
		if (powerlink(i)==.true. .or. pumplink(i)== .true.) then !for energy
			tenergy=tenergy+1
			k2=k2+1
			energyln(k2)=ENERGYconv(i)  !Beg. storage array
			Energy_Name(k2)=TRIM(LName(i)) !Name for begstorage nodes
		endif
		!for routing links
		if (L_Method(i)>= 1) then
			troute=troute+1
			k3=k3+1
			volume(k3)=TOT_LVOLconv(i)
			Volume_Name(k3)=TRIM(LName(i))
		endif
		!for loss method 2
		if (LossMethod(i)== 2) then
				tgeo=tgeo+1
				k4=k4+1
				GeoName(k4)=TRIM(LName(i))
				depth(k4)=LastDepth(i)
				width(k4)=LastWidth(i)
				velocity(k4)=LastVelocity(i)
		endif	
	enddo	
	if (sysstat(nday)==1 .and. step==1) then  !Initialize column names	
		do i=1,links
		begflow(i)='Beg_Flow'
		begflowTS(i)='Beg_FlowTS'
		enddo

		do i=1,tLoss
		CLoss(i)='Loss'
		enddo

		do i=1,tend
		Cend(i)='End_Flow'
		enddo

		do i=1,tenergy
		Cenergy(i)='Energy'
		enddo
		
		do i=1,troute
		Cvolume(i)='Volume'
		enddo 

		do i=1,tgeo
		cdepth(i)='Depth'
		cwidth(i)='Width'
		cvelocity(i)='Velocity'
		enddo

		write(ioutlinksts,*)'Run number = ', SysStat(run)
		
		write(ioutlinksts,FMT=11,advance='no')     !Write the parameter name                 
     &	(begflow(i),begflowTS(i),i=1,links),(Cend(i),i=1,tend),
     &	(CLoss(i),i=1,tLoss),(cenergy(i),i=1,tenergy),
     &	(Cvolume(i),i=1,troute),(Cwidth(i),i=1,tgeo),
     &	(Cdepth(i),i=1,tgeo),(Cvelocity(i),i=1,tgeo)	
		write(ioutlinksts,*)
		
		write(ioutlinksts,FMT=10,advance='no')'Day','Step', !Write the name of node
     &	(LName(i),LName(i),i=1,links),(end_Name(i),i=1,tend),
     &	(Loss_Name(i),i=1,tLoss),(energy_name(i),i=1,tenergy),
     &	(Volume_Name(i),i=1,troute),(GeoName(i),i=1,tgeo),
     &	(GeoName(i),i=1,tgeo),(GeoName(i),i=1,tgeo)
	write(ioutlinksts,*)
	endif !ends if loop for day=1
	

	write(ioutlinksts,FMT=9,advance='no')sysstat(nday),step, !Write the data for each day 
     &	(BQLNconv(i),DBQLNconv(i),i=1,links),(endln(i),i=1,tend), 
     &  (loss(i),i=1,tloss),(energyln(i),i=1,tenergy),
     &  (volume(i),i=1,troute),
     &	(width(i),i=1,tgeo),(depth(i),i=1,tgeo),(velocity(i),i=1,tgeo)
    
	
	write(ioutlinksts,*)


	
9     FORMAT(2I6,10000F50.4)
10    FORMAT(2A5,9X,10000A50)
11	FORMAT(8X, 10000A50)
100   Continue
	end subroutine
!************************************************************************
     	subroutine DayOutputNodeTS(iDay,step,DSTO,DTREL)
	!Created by Evgenii 0907
	! USED TO SEE OUTPUT OF SUB_TIME_STEPS
      IMPLICIT NONE
      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
!  INPUT
      INTEGER*4 iDay,NN_SEQ(NODMAX)      
!  Local
      INTEGER*4 i,j,NN,step
      INTEGER*4 j0, j1,k,k1,k2,k3,k4,k5
	INTEGER*4 tevap_nodes, t_storage, t_gage,tseep
	INTEGER*4 storage_id(nodmax),gage_id(nodmax),t_cons
	INTEGER*4 cons_id(NODMAX),evap_id(NODMAX),seep_id(NODMAX)
	real:: storage_beg(NODMAX), storage_end(NODMAX),gage(nodmax)
     &  ,DTREL(NODMAX)	
	real:: DSTO(nodmax),inflowtsconv(nodmax),INFLOWconv(nodmax)
	real:: TEVAPNconv(nodmax),TSEEPLconv(nodmax),TOTRELconv(nodmax)
	real:: CONSUMPTIONconv(nodmax),BSTOconv(nodmax),DSTOconv(nodmax)
	real:: consump(nodmax),tevap(NODMAX),seepage(NODMAX),DTRELconv(NODMAX) !Begining and end storage of storage nodes
	character(10) storage_Name(NODMAX),gage_name(nodmax)  !Name of storage
	character(10) cons_name(NODMAX),evap_name(NODMAX)
	character(10) flowin(NODMAX),outflow(nodmax),gagein(NODMAX),
     &	begsto(NODMAX),endsto(NODMAX),cons(NODMAX),evap(NODMAX),
     &	seep_name(NODMAX),seep(NODMAX),inflowst(nodmax),outflowTS(nodmax)
!------------------------------------------------------------------------
	t_cons=0; tevap_nodes=0; t_storage=0; t_gage=0; tseep=0
	
	j=0; k=0; k1=0; k2=0; k3=0; k4=0; k5=0	
	do i=1,tnodes
	  !Convert units from Mm3/time-step to what guage input units were
        inflowtsconv(i)=inflowts(i)
        call UndoUnitConversion(1,uFlow,inflowtsconv(i))
        inflowtsconv(i)=inflowtsconv(i)/SYSSTAT(NPER)
        INFLOWconv(i)=INFLOW(i)
        call UndoUnitConversion(1,uFlow,INFLOWconv(i))
        INFLOWconv(i)=INFLOWconv(i)/SYSSTAT(NPER)
        TEVAPNconv(i)=TEVAPN(i)
        call UndoUnitConversion(1,uFlow,TEVAPNconv(i))
        TEVAPNconv(i)=TEVAPNconv(i)/SYSSTAT(NPER)
        TSEEPLconv(i)=TSEEPL(i)
        call UndoUnitConversion(1,uFlow,TSEEPLconv(i))
        TSEEPLconv(i)=TSEEPLconv(i)/SYSSTAT(NPER)
        TOTRELconv(i)=TOTREL(i)
        call UndoUnitConversion(1,uFlow,TOTRELconv(i))        
        TOTRELconv(i)=TOTRELconv(i)/SYSSTAT(NPER)
        DTRELconv(i)=DTREL(i)
        call UndoUnitConversion(1,uFlow,DTRELconv(i))
        DTRELconv(i)=DTRELconv(i)/SYSSTAT(NPER)
        CONSUMPTIONconv(i)=CONSUMPTION(i)
        call UndoUnitConversion(1,uFlow,CONSUMPTIONconv(i))  
        CONSUMPTIONconv(i)=CONSUMPTIONconv(i)/SYSSTAT(NPER)
        BSTOconv(i)=BSTO(i)
        DSTOconv(i)=DSTO(i)
        call UndoUnitConversion(1,uVol,BSTOconv(i))
        call UndoUnitConversion(1,uVol,DSTOconv(i))
        
      end do
 
	
	
	
	do i=1,tnodes !make arrays for the output for each category
		NN = NODSEQ(I)
		IF (NN .EQ. 0) GO TO 100
		!Makes arrays for output of beginning storage (only for storage nodes)
		if (ResvNode(i) .or. NatLak(i) .or. GWNODE(i)) then
			t_storage=t_storage+1
			k=k+1
			storage_beg(k)=BSTOconv(i)  !Beg. storage array
			storage_Name(k)=TRIM(NNAME(i)) !Name for begstorage nodes
			storage_id(k)=nn   !node id
			storage_end(k)=DSTOconv(i)
		endif

		
	
		if (GageNF(i) == .true.) then
			t_gage=t_gage+1
			k2=k2+1
			Gage_Name(k2)=TRIM(NNAME(i)) 
			Gage_id(k2)=nn
			gage(k2)=NQinn(i)
			call UndoUnitConversion(1,uFlow,gage(k2)) 
		endif
		
		if (cons_node(i)== .true.)then
			t_cons=t_cons+1
			k3=K3+1
			cons_name(k3)=TRIM(NNAME(i)) 
			cons_id(k3)=nn
			consump(k3)=CONSUMPTIONconv(i) 
		endif

		if (ResvNode(i).or. NatLak(i) .or. GWNODE(i))then
			tevap_nodes=tevap_nodes+1
			k4=K4+1
			evap_name(k4)=TRIM(NNAME(i)) 
			evap_id(k4)=nn
			tevap(k4)=TEVAPNconv(i) 
		endif
	
		if (seep_node(i)== .true.)then !seep_node defined read_sim_data line 723
			tseep=tseep+1
			k5=k5+1
			seep_name(k5)=TRIM(NNAME(i)) 
			seep_id(k5)=nn
			seepage(k5)=TSeeplconv(i) 
		endif
	enddo	

	if (sysstat(nday)==1 .and. step==1) then 
		do i=1,tnodes
		    flowin(i)='Inflow';outflow(i)='Outflow';inflowst(i)='Inflow_TS'
		    outflowTS(i)='OutflowTS'
		enddo

		do i=1,t_gage !t_gage defined in initsys.for
		gagein(i)='Flow_Input'
		enddo

		do i=1,t_storage !t_storage defined in initsys.for
		begsto(i)='Beg_Vol'; endsto(i)='End_Vol'
		enddo

		
		do i=1,t_cons
		cons(i)='Cons'
		enddo
		
		do i=1,tevap_nodes
		evap(i)='Evap'
		enddo

		do i=1,tseep
		seep(i)='Seepage'
		enddo
		write(ioutnodests,FMT=11,advance='no')     !Write the output category                 
     &	(begsto(i),i=1,t_storage),(endsto(i),i=1,t_storage), 
     &	(flowin(i),inflowst(i),outflow(i),outflowTS(i),i=1,tnodes),
     &	(gagein(i),i=1,t_gage),  
     &	(cons(i),i=1,t_cons),(evap(i),i=1,tevap_nodes),		
     &	(seep(i),i=1,tseep)
		write(ioutnodests,*)
		
		write(ioutnodests,FMT=10,advance='no')'Day ','step', !Write the name of node
     &	(storage_Name(i),i=1,t_storage),(storage_Name(i),i=1,t_storage),
     &	(NName(i),NName(i),NName(i),NName(i),i=1,tnodes),
     &  (Gage_Name(i),i=1,t_gage)    
     &	,(cons_name(i),i=1,t_cons),	
     &	(evap_name(i),i=1,tevap_nodes),(seep_name(i),i=1,tseep)
		
		write(ioutnodests,*)
	endif !ends if loop for day=1
	

	write(ioutnodests,FMT=9,advance='no')sysstat(nday),step, !Write the data for each day 
     &	(storage_beg(i),i=1,t_storage),(storage_end(i),i=1,t_storage),
     &	 (inflowconv(i),inflowtsconv(i),TOTRELconv(i),
     &  DTRELconv(i),i=1,tnodes),(gage(i),i=1,t_gage),
     &	(consump(i),i=1,t_cons),
     &	(tevap(i),i=1,tevap_nodes),(seepage(i),i=1,tseep)
	
	write(ioutnodests,*)

	
9     FORMAT(2I8,10000F50.4)
10    FORMAT(2A8,8X,10000A50)
11	FORMAT(12X, 10000A50)


 100   Continue
	 end subroutine
!***************************************************************************************
!     Subroutine below not called in IRAS-2010
!	 subroutine DayOutputBinary(iOutfileBin,bWriteHeader,iDay)
! Read gage data from gage file
! Compute incremental flows for all nodes
! Save incremental flows to a temporary binary file
! Binary file format:
! Year1, Day1, Inc natural flow, inc added flow of node 1(# in simulation sequence)
! Year1, Day1, Inc natural flow, inc added flow of node 2

!      IMPLICIT NONE
!      INCLUDE 'IRAS_SYS.INC'
!      INCLUDE 'NODE.INC'
!      INCLUDE 'LINK.INC'
!  INPUT
!      INTEGER*4 iOutFileBin,iDay              !iDay=1,365
!      LOGICAL*1 bWriteHeader
!  Local
!      INTEGER*4 i,j,iRecNo,TNRecNL1,TNRecNodeAll,nHeader
!      INTEGER*4 NodeLines, LinkLines, j0, j1
!------------------------------------------------------------------------
      !record: Year, Day, incremental natural flow, incremental added flow
      !        i*2, i*2, r*4, r*4
      !RecLen = BIT_SIZE(i)*2 + BIT_SIZE(i)*2*2
!      nHeader = 1
!      NodeLines = FLOOR((-1)*REAL(TNodes)/20.0)*(-1)
!	LinkLines = FLOOR((-1)*REAL(Links)/20.0)*(-1)
      
!	if (bWriteHeader) then
	!first line of header
!        WRITE(iOutfileBin,rec=1) SYSSTAT(YRSTART), SYSSTAT(YREND)
!     &        ,TNodes, Links
        !next for nodeID(): 20 data/line
!        do i = 1, NodeLines
!          j0 = (i-1)*20 + 1
!          j1 = MIN(TNodes, j0+20-1)
!          WRITE(iOutfileBin,rec=nHeader+i)(NodeID(j),j=j0,j1)
!        end do
!        nHeader = nHeader + NodeLines
        !next for linkID(): 20 data/line
!        do i = 1, LinkLines
!          j0 = (i-1)*20 + 1
!          j1 = MIN(Links, j0+20-1)
!          WRITE(iOutfileBin,rec=nHeader+i)(LinkID(j),j=j0,j1)
!        end do
!        return
!      end if

!     nHeader = 1 + NodeLines + LinkLines

      !total num of records for a node/link in whole simulation
!      TNRecNL1 =  (SYSSTAT(YREND) - SYSSTAT(YRSTART) + 1)*365
      !total num of records for all nodes
!      TNRecNodeAll = TNodes*TNRecNL1
      !write nodes' results
!      do i = 1, TNodes
        !compute record No
!        iRecNo=nHeader+(i-1)*TNRecNL1+365*(SysStat(Sim_Year)-1) + iDay
!        WRITE(iOutFileBin, rec=iRecNo) NodeID(i),NName(i)
!     &       ,SysStat(Year),SysStat(Month),SysStat(Day)
!     &       ,inflow(i),NQinn(i),AddQinn(i)
!     &       ,TOTREL(i),CONSUMPTION(i),TEvapn(i)+TSeepl(i)
!     &       ,BSTO(i),ESTO(i)
!      end do
      !write links' results
!      do i = 1, Links
        !compute record No
!        iRecNo = nHeader + TNRecNodeAll
!     &         + (i-1)*TNRecNL1 + 365*(SysStat(Sim_Year)-1) + iDay
!        WRITE(iOutFileBin, rec=iRecNo) LinkID(i),LName(i)
!     &       ,SysStat(Year),SysStat(Month),SysStat(Day)
!     &       ,BQLN(i),EQLN(i),TLossL(i),TOT_LVOL(i)
!     &       ,LastDepth(i), LastWidth(i), LastVelocity(i)
!      end do
!      end subroutine
!Copyright (c) 2009, 2010 by University College London, Cornell University
!Authors:
!Evgenii Matrosov (evgenii.matrosov@ucl.ac.uk), Julien Harou (j.harou@ucl.ac.uk), 
!Daniel P. Loucks (dpl3@cornell.edu), Marshall Taylor, Peter French, 
!This program is free software under the General Public Licence, GPL (>=v2)
!Read the 'GPL License.txt' file distributed with this source code for a full license statement.
!
!	 *************************************************************************************
!      SUBROUTINE DayOutputText(iOutFile,iYear,iMonth,iDay,NN_SEQ)
!C
!C     940109:PNF  Adjust flow variables for flow factor
!C     950116:PNF  TXCLSC only if PR_LUN = 5
!C
!!	 *************************************************************************************
!C
!      IMPLICIT NONE
!      INCLUDE 'IRAS_SYS.INC'
!      INCLUDE 'NODE.INC'
!      INCLUDE 'LINK.INC'
!
!C     Arguments
!      INTEGER*4    iOutFile, iYear,iMonth,iDay, NN_SEQ(nodmax)  !Evgenii changed (8) to (nodemax) 090804
!      LOGICAL*1    Append
!C
!
!C     Local variables
!      REAL*4       TOT_ALLOC
!      INTEGER*4    I, J, JJ, NN, LN, ST, MM
!      LOGICAL*1    STOR
!C
!!-------------------------------------------------------------------------
!C
!	
!	WRITE(UNIT=iOutFile,FMT=*)'' !Evgenii 092204 changed * to FMT=*
!	WRITE(UNIT=iOutFile,FMT=10) iYear, iMonth, iDay
!      DO 100 I = 1,NODES
!         NN = NN_SEQ(I)
!         IF (NN .EQ. 0) GO TO 100
!         STOR = .FALSE.
!         IF (CAPN(NN).GT.0. .AND. .NOT.GWNODE(NN)) STOR = .TRUE.
!         WRITE(UNIT=iOutFile,FMT=*)TRIM(NNAME(NN)) !!Evgenii 090422 changed * to FMT=*
!         WRITE(UNIT=iOutFile,FMT=11)
!         WRITE(UNIT=iOutFile,FMT=15)
!     1               NN, QINN(NN), INFLOW(NN),
!     1               BSTO(NN),TEvapn(NN)+TSeepl(nn),CONSUMPTION(NN),
!     1               TOTREL(NN), ESTO(NN)
!C
!C
!C        Do the links
!         TOT_ALLOC = 0.0
!         DO J = 1,TOTOUT(NN)
!            LN = OUTLNK(NN,J)
!C
!            IF(j==1)WRITE(UNIT=iOutFile,FMT=40)
!C
!            WRITE(UNIT=iOutFile,FMT=45) LN, LNAME(LN), BQLN(LN),
!     1                    EQLN(LN), TOT_LVOL(LN), ENERGY(LN)
!            TOT_ALLOC = TOT_ALLOC + BQLN(LN)
!C
!         ENDDO
!         TOT_ALLOC = TOT_ALLOC - TOTREL(NN)
!         IF(TOT_ALLOC.NE.0.0)THEN
!           WRITE(UNIT=iOutFile,FMT=6010)NN,ABS(TOT_ALLOC)
!         ENDIF
!	   write(UNIT=iOutFile,fmt=*)' '
!C
!100   CONTINUE
!C
!9999  RETURN
!C
!12	FORMAT(' ',I2) !Evgenii 090422
!13	FORMAT(' ',I60) !Evgenii
!10    FORMAT('Year: ',I5,', Month: ',I2,', Day: ',I2)
!11    FORMAT('NN   I Inflo   T Inflo    Beg Stor   T.Loss',
!     &       '     Consumpt.  Outflow    End Stor')
!15    FORMAT(I2,1X,7G11.5)
!40    FORMAT(5X,'LN   Name    Beg Flow   End Flow   Volume',
!     1      '     Energy')
!45    FORMAT(5X,I2,1X,A8,1X,4G11.3)
!6010  FORMAT(' Unallocated water at node ',I2,' = ',G12.3)
!C
!      END
!C
