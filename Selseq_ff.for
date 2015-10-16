!Copyright (c) 2000 Cornell University
!Authors:
!Daniel P. Loucks (dpl3@cornell.edu), Marshall Taylor, Peter French
!This program is free software under the General Public Licence, GPL (>=v2)
!Read the 'GPL License.txt' file distributed with this source code for a full license statement.
!
!	 *************************************************************************************

C
      SUBROUTINE SELSEQ
C
C      USE:      Used to define the sequence of nodes to be
C                simulated in IRIS
C
C      INPUT:    rules   nodes that determine the processing order:
C                        Where rule(:,1) must be processed before
C                        rule(:,2) FF 10/2012
C                nrules  number of rules FF 10/2012
C
C      OUTPUT:     NNSEQ (from database read), NODES
C
C      NOTES:      Based on algorithm developed by J. Andreu
C
C      Modifications:
C      MRT  9201         Major mods to incorporate database manager.
C      MRT  921022       Updated for looping networks as per DPL's most
C                        recent revisions.
!	 PL   000518	   Last channge


!	 *************************************************************************************

C
!     Input or known global variables:
!       TNodes, Links
!       TotIn(), TotOut()
!       InLink(i,j), OutLnk(i,j)
!     Output global variable:
!       NodSeq()-Sequence of nodes to be simulated in IRAS

!	 *************************************************************************************
      implicit none

      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
C
C      Local Variables:
      INTEGER*4   ST, NN, LL, LN, SL(LNKMAX)
      INTEGER*4   FLAG, CHFLAG, LOOPND, SN(NODMAX)
      INTEGER*4   I, TYPE
      INTEGER*4   NSEQ(NODMAX) 
      INTEGER*4   streamorder(NODMAX) ! added by FF 10/2012
      INTEGER*4   streamorderlink(LNKMAX) ! added by FF 09/2012
      LOGICAL:: FFFLAGDONE  ! added by FF 08/2012
      INTEGER:: KK ! added by FF 08/2012
      integer::pass
C
!------------------------------------------------------------------------
C
      pass = 0
      FFFLag = 0
      TYPE   = 0
      LOOPND = 0

C
C     Check system's status.
      IF (TNODES .LE. 0 ) RETURN
C
C
C     Initialize variables. No sequence defined yet.
  1   CONTINUE
      DO I = 1,TNODES
         NSEQ(I) = 0
         SN(I) = 0
         FFFLAG(I) = -1 ! added by FF 08/2012
         streamorder(i) = 0 ! added by FF 09/2012
         streamorderlink(i) = 0 ! added by FF 09/2012
      END DO
C
      DO LN = 1,LINKS
        SL(LN) = 0
      END DO
      NODES = 0
      CHFLAG = 0      
C
C     Define sequence.
  10  CONTINUE

      IF (CHFLAG .EQ. 2) CHFLAG = 3   !  no change last two passes
      IF (CHFLAG .EQ. 0) CHFLAG = 2   !  no change last pass
      IF (CHFLAG .EQ. 1) CHFLAG = 0   !  at least one node added last pass
      DO 20 NN = 1,TNODES
C        Check if already in sequence.
         IF (SN(NN) .EQ. 1) GO TO 20
C        Check if initial node with no incoming links.
         IF(TOTIN(NN) .EQ. 0) THEN
            NODES = NODES + 1
            NSEQ(NODES) = NN   
            SN(NN) = 1
            CHFLAG = 1
            FFFLAG(NN) = 0 ! added by FF 08/2012
            streamorder(nn) = 1 ! added by FF 09/2012
C
C           Flag all outgoing links from this node
            DO LL = 1,TOTOUT(NN)
               LN = OUTLNK(NN,LL)
               SL(LN) = 1
               streamorderlink(LN) = 1
            END DO
            GO TO 20
         END IF
C
C        Not a start node. Check if all incoming links flaged.
         FLAG = 0
         DO LL = 1,TOTIN(NN)
            LN = INLINK(NN,LL)
            IF(SL(LN) .NE. 1) FLAG = 1
         END DO
C
C        If flag is still 0 this means all incoming links come
C        from nodes already in the sequence.
C        If flag is 1 and CHFLAG is 3, indicates a end-of-loop node
C        Highlight this node and leave highlighted.
         IF(FLAG.EQ.0 .OR. (FLAG.EQ.1 .AND. CHFLAG.EQ.3)) THEN

            NODES = NODES + 1

            NSEQ(NODES) = NN  
            SN(NN) = 1 
            
            DO LL = 1,TOTIN(NN)
               LN = INLINK(NN,LL)
               IF (streamorder(NN) < streamorderlink(LN)) THEN
                  streamorder(NN) = streamorderlink(LN)
               END IF
            END DO              ! added by FF
            streamorder(NN) = streamorder(NN) + 1

            FFFLAG(NN) = 0  ! added by FF 08/2012
            IF (FLAG == 1 .and. chflag == 3) THEN 

               IF (.NOT.FFFLAGDONE) THEN
                  DO KK = 1,TNODES
                     IF (FFFLAG(KK) == -1) THEN
                        FFFLAG(KK) = 1
                     END IF
                  END DO 
                  FFFLAGDONE = .TRUE.
               END IF
               FFFLAG(NN) = 2
            END IF
 	    
            IF(CHFLAG .EQ. 3) LOOPND = LOOPND + 1
            CHFLAG = 1
C
C           Flag all outgoing links
            DO LL = 1,TOTOUT(NN)
               LN = OUTLNK(NN,LL)
               SL(LN) = 1
               streamorderlink(ln) = streamorder(nn)
            END DO
            GO TO 20
         END IF
C
  20  CONTINUE
C
      IF (CHFLAG .EQ. 3 )GO TO 9999
      IF (NODES .LT. TNODES) GO TO 10
C

 9999 CONTINUE
!      ST = WR_FND(SYSGRP,NNSEQ,0,0,NSEQ,TNODES, 1, TYPE)
      do i = 1, TNodes
        NODSEQ(i) = Nseq(i)
      end do
	
      FFSN = SN

      !CALL changeseq(streamorder) this has now been moved to after source priorities are read

      open(8989,file="iras.str")
      do i = 1, tnodes
         write(8989,*) Nodeid(i), streamorder(i)
      end do
      close(8989)
      
      RETURN
      END

!	 *************************************************************************************
      SUBROUTINE PrintSeq()
C
C      USE:     Prints the nodes sequence to file seq.out
C               
C      INPUT:    None
C
C      OUTPUT:     none
C
C      NOTES:      Evgenii Matrosov 261011
C
!	 *************************************************************************************
      implicit none

      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
C
C      Local Variables:
       INTEGER*4   ioutSeq,i

!------------------------------------------------------------------------
      ioutSeq=35
      OPEN(UNIT=ioutSeq, FILE='seq.out',STATUS='replace') 
      do i=1,tnodes
        write(ioutSeq, fmt=10)NodeID(NodSeq(i)), nname(NodSeq(i)), 
     & 	    FFFlag(NodSeq(i))
      end do

      close (ioutSeq)      

10    FORMAT(I20,A50,i2)    

      return  
      end
!	 *************************************************************************************      
      SUBROUTINE SetUserSeq()
C
C      USE:     Prints the nodes sequence to file seq.out
C               
C      INPUT:    None
C
C      OUTPUT:     none
C
C      NOTES:      Evgenii Matrosov 261011
C
!	 *************************************************************************************
      implicit none

      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
C
C      Local Variables:
       INTEGER*4   userID(nodmax),iINSeq,i,j,LineCounter
       logical*1 success,CountFileLines
!------------------------------------------------------------------------
      iINSeq=54
      OPEN(UNIT=iINSeq, FILE='input/iras.seq',STATUS='old') 
      !IF(.NOT. CountFileLines(ioutSeq, LineCounter)) GOTO 999 
      do i=1,tnodes
        read(iINSeq, *)userID(i)
      end do
      do i=1,tnodes
        do j=1,tnodes
            if(userID(i)==nodeid(j)) then           
                NODSEQ(i)=j
            end if
        end do
      end do      
      success=.true.
999   close (iINSeq)      
      if (.not.success) then
        write(*,*)'The number of nodes in the iras.seq file do not 
     &    match the number of nodes in the iras.inp file'
      end if

      return  
      end 

!	 *************************************************************************************

C
      SUBROUTINE SELSEQ_transfers(trans_start, trans_source, 
     &     trans_nsource, 
     &     ntrans, nsource)       ! include transfers priorities FF 10/2012
C
C      USE:      Used to define the sequence of nodes to be
C                simulated in IRIS
C
C      INPUT:    rules   nodes that determine the processing order:
C                        Where rule(:,1) must be processed before
C                        rule(:,2) FF 10/2012
C                nrules  number of rules FF 10/2012
C
C      OUTPUT:     NNSEQ (from database read), NODES
C
C      NOTES:      Based on algorithm developed by J. Andreu
C
C      Modifications:
C      MRT  9201         Major mods to incorporate database manager.
C      MRT  921022       Updated for looping networks as per DPL's most
C                        recent revisions.
!	 PL   000518	   Last channge

!      FF   08/2012      Modifications to add flags where process order fails
!                        - FFFLAG = 0, no problem
!                        - FFFLAG = 1, problem upstream of this node
!                        - FFFLAG = 2, problem located near this node

!      FF   10/2012     Modification to allow changing of sequence order
!                       according to transfers.  Where trans_sources must be
!                       processed first before trans_start

!	 *************************************************************************************

C
!     Input or known global variables:
!       TNodes, Links
!       TotIn(), TotOut()
!       InLink(i,j), OutLnk(i,j)
!     Output global variable:
!       NodSeq()-Sequence of nodes to be simulated in IRAS

!	 *************************************************************************************
      implicit none

      INCLUDE 'IRAS_SYS.INC'
      INCLUDE 'NODE.INC'
      INCLUDE 'LINK.INC'
C
C      Local Variables:
      INTEGER*4   ST, NN, LL, LN, SL(LNKMAX)
      INTEGER*4   FLAG, CHFLAG, LOOPND, SN(NODMAX)
      INTEGER*4   I, TYPE
      INTEGER*4   NSEQ(NODMAX) 

!     Added by FF 10/2012
      integer::ntrans                ! number of transfer links
      integer::nsource               ! max no. of sources per transfer link
      integer::trans_start(ntrans)   ! start node of transfer link
      integer::trans_source(ntrans,nsource) ! source nodes for end node of transfer link
      integer::trans_nsource(ntrans) ! no. of source nodes for each end node of transfer link

      integer::streamorder(NODMAX)     ! stream order
      integer::streamorderlink(LNKMAX) ! stream order link
      logical::ffflagdone
      integer::KK 
      integer::pass
      integer::itrans
      integer::isource
C
!------------------------------------------------------------------------
C
      pass = 0
      FFFLag = 0
      TYPE   = 0
      LOOPND = 0

C
C     Check system's status.
      IF (TNODES .LE. 0 ) RETURN
C
C
C     Initialize variables. No sequence defined yet.
  1   CONTINUE
      DO I = 1,TNODES
         NSEQ(I) = 0
         SN(I) = 0
         FFFLAG(I) = -1 ! added by FF 08/2012
         streamorder(i) = 0 ! added by FF 09/2012
         streamorderlink(i) = 0 ! added by FF 09/2012
      END DO
C
      DO LN = 1,LINKS
        SL(LN) = 0
      END DO
      NODES = 0
      CHFLAG = 0      
C
C     Define sequence.
  10  CONTINUE

      IF (CHFLAG .EQ. 2) CHFLAG = 3   !  no change last two passes
      IF (CHFLAG .EQ. 0) CHFLAG = 2   !  no change last pass
      IF (CHFLAG .EQ. 1) CHFLAG = 0   !  at least one node added last pass
      DO 20 NN = 1,TNODES
C        Check if already in sequence.
         IF (SN(NN) .EQ. 1) GO TO 20
C        Check if initial node with no incoming links.
         IF(TOTIN(NN) .EQ. 0) THEN
            NODES = NODES + 1
            NSEQ(NODES) = NN   
            SN(NN) = 1
            CHFLAG = 1
            FFFLAG(NN) = 0 ! added by FF 08/2012
            streamorder(nn) = 1 ! added by FF 09/2012
C
C           Flag all outgoing links from this node
            DO LL = 1,TOTOUT(NN)
               LN = OUTLNK(NN,LL)
               SL(LN) = 1
               streamorderlink(LN) = 1
            END DO
            GO TO 20
         END IF
C
C        Not a start node. Check if all incoming links flaged.
         FLAG = 0
         DO LL = 1,TOTIN(NN)
            LN = INLINK(NN,LL)
            IF(SL(LN) .NE. 1) FLAG = 1
         END DO

C     only set flag to zero if transfer priorities are statisfied - FF 10/2012
         if (ntrans > 0) then
            do itrans = 1, ntrans
               if (nn == trans_start(itrans)) then ! check whether current node is affected by transfers

                  do isource = 1, trans_nsource(itrans)
                     if (sn(trans_source(itrans,isource)) == 0) then ! check whether all other sources have been calculated yet
                        flag = 1
                     end if
                  end do
               end if
            end do
!            do itrans = 1, ntrans
!               if (trans_start(itrans) == 744) then ! check whether current node is affected by transfers
!                  print *, nodeid(nn),flag,chflag,trans_nsource(itrans)
!                  do isource = 1, trans_nsource(itrans)
!                     print *, nodeid(trans_source(itrans,isource)),
!     &                    sn(trans_source(itrans,isource))
!                  end do
!                  print *,"----------------------------"
!               end if
!            end do
         end if
C
C        If flag is still 0 this means all incoming links come
C        from nodes already in the sequence.
C        If flag is 1 and CHFLAG is 3, indicates a end-of-loop node
C        Highlight this node and leave highlighted.
         IF(FLAG.EQ.0 .OR. (FLAG.EQ.1 .AND. CHFLAG.EQ.3)) THEN

            NODES = NODES + 1

            NSEQ(NODES) = NN  
            SN(NN) = 1 
            
C     store the stream order of each node - FF 10/2012
            DO LL = 1,TOTIN(NN)
               LN = INLINK(NN,LL)
               IF (streamorder(NN) < streamorderlink(LN)) THEN
                  streamorder(NN) = streamorderlink(LN)
               END IF
            END DO 
            streamorder(NN) = streamorder(NN) + 1

            FFFLAG(NN) = 0  ! added by FF 08/2012
            IF (FLAG == 1 .and. chflag == 3) THEN 

               IF (.NOT.FFFLAGDONE) THEN
                  DO KK = 1,TNODES
                     IF (FFFLAG(KK) == -1) THEN
                        FFFLAG(KK) = 1
                     END IF
                  END DO 
                  FFFLAGDONE = .TRUE.
               END IF
               FFFLAG(NN) = 2               
            END IF
 	    
            IF(CHFLAG .EQ. 3) LOOPND = LOOPND + 1
            CHFLAG = 1
C
C           Flag all outgoing links
            DO LL = 1,TOTOUT(NN)
               LN = OUTLNK(NN,LL)
               SL(LN) = 1
               streamorderlink(ln) = streamorder(nn) ! added by FF 10/2012
            END DO
            GO TO 20
         END IF
C
  20  CONTINUE
C
      IF (CHFLAG .EQ. 3 )GO TO 9999
      IF (NODES .LT. TNODES) GO TO 10
C

 9999 CONTINUE
!      ST = WR_FND(SYSGRP,NNSEQ,0,0,NSEQ,TNODES, 1, TYPE)
      do i = 1, TNodes
        NODSEQ(i) = Nseq(i)
      end do
	
      FFSN = SN


      open(8989,file="iras.str1")
      do i = 1, tnodes
         write(8989,*) Nodeid(i), streamorder(i)
      end do
      close(8989)

      OPEN(UNIT=9000, FILE='seq1.out',STATUS='replace') 
      do i=1,tnodes
        write(9000, fmt=109)NodeID(NodSeq(i)),nname(NodSeq(i)),
     & 	    FFFlag(NodSeq(i))
      end do

      close (9000)

 109  FORMAT(I20,A50,i2)    
      
      RETURN
      END

