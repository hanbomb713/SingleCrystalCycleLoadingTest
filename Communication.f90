!****************************************************************************
! This file contains all the subroutines which do the communication jobs.
!****************************************************************************
Module CommunicationMD
	USE ZQMPI
	use TreeLeaf
	CONTAINS
!****************************************************************************
!   bcastSomeVariables:
!                 Broadcase some variables read in dataReader by master.
!         Called by BcastAllGVs (broadcast global variables).
!****************************************************************************
subroutine bcastSomeVariables
        use ZQMPI
        use vectors
        use MPIModule,only:iNTotalDPs,iNMAXS
        use variables,only:MAx_Quad,N_Times,N_TimesCopy,iLoop_Time,checkNeiBur, &
                        iBeginLoop,Mu,NU,LATTICE,Delta_Sig, MAX_NODE,&
                        A_Cube,Mobility,DTime,StrainRate,ElasticCons,dcrit,iNPoint,&
                        MAX_PLANE,IntegrateMethod,APPLIED_SIG,RMIN,totalStrain,totalTime, &
                        logAnnihilation, logArrangeLoop,CycFreq,dpCycF0,dpCycStrain, MAX_NEIGHBOR,  &
                        Load,iArr,iObstacle,iInertial,totalStress,Loading_Dir,iMatType,PEIERLS,MOB_R,PEI_R, &
			ATRatio,TRatio,NGRatio
        implicit none

        integer,dimension(19) :: tmpControl
        double precision,dimension(27) ::tmpPara

        if(iMyid==iMaster)then
                iNMaxS = iNPoint/(iNumprocs-1)+1

                tmpControl(1)=MAX_QUAD
                tmpControl(2)=N_TIMES
                tmpControl(3)=N_TIMESCOPY
                tmpControl(4)=ILOOP_TIME
                tmpControl(5)=CheckNeiBur
                tmpControl(6)=iNMaxS
                tmpControl(7)=iBeginLoop
                tmpControl(8)=iNTotalDPs
                tmpCOntrol(9)=MAX_NODE
                tmpControl(10)=MAX_PLANE
                tmpControl(11)=IntegrateMethod
                tmpControl(12)=logAnnihilation
                tmpControl(13)=logArrangeLoop
                tmpControl(14)=MAX_NEIGHBOR
                tmpControl(15)=Load
                tmpControl(16)=iArr
                tmpControl(17)=iObstacle
                tmpControl(18)=iInertial
		tmpControl(19)=iMatType

                tmpPara(1)=MU
                tmpPara(2)=NU
                tmpPara(3)=LATTICE
                tmpPara(4)=totalStress
                tmpPara(5)=DELTA_SIG
                tmpPara(6)=A_CUBE
                tmpPara(7)=MOBILITY
                tmpPara(8)=DTIME
                tmpPara(9)=StrainRate
                tmpPara(10)=ElasticCons
                tmppara(11)=dcrit
                tmpPara(12)=APPLIED_SIG
                tmpPara(13)=RMIN
                tmpPara(14)=totalStrain
                tmpPara(15)=totalTime
                tmpPara(16)=CycFreq
                tmpPara(17)=dpCycF0
                tmpPara(18)=dpCycStrain
                tmpPara(19)=Loading_Dir%v(1)
                tmpPara(20)=Loading_Dir%v(2)
                tmpPara(21)=Loading_Dir%v(3)
		tmpPara(22)=PEIERLS
		tmpPara(23)=MOB_R
		tmpPara(24)=PEI_R
		tmpPara(25)=ATRatio
		tmpPara(26)=TRatio
		tmpPara(27)=NGRatio
        endif


        call MPI_BCAST(tmpControl(1),19,MPI_INTEGER,iMaster,MPI_COMM_WORLD,iIerr)

        call MPI_BCAST(tmpPara(1),27,MPI_DOUBLE_PRECISION,iMaster,MPI_COMM_WORLD,iIerr)

        if(iMyid/=iMaster)then

                MAX_QUAD         = tmpControl(1)
                N_TIMES          = tmpControl(2)
                N_TIMESCOPY          = tmpControl(3)
                ILOOP_TIME       = tmpControl(4)
                CheckNeiBur      = tmpControl(5)
                iNMaxS           = tmpControl(6)
                iBeginLoop       = tmpControl(7)
                iNTotalDPs               = tmpControl(8)
               MAX_NODE                 = tmpCOntrol(9)
                MAX_PLANE                = tmpControl(10)
                IntegrateMethod          = tmpControl(11)
                logAnnihilation         =       tmpControl(12)
                logArrangeLoop          =       tmpControl(13)
                MAX_NEIGHBOR            =       tmpControl(14)
                Load                            =       tmpControl(15)
                iArr                            =       tmpControl(16)
                iObstacle                       =       tmpControl(17)
                iInertial                       =       tmpControl(18)
		iMatType			=	tmpControl(19)

                MU                  = tmpPara(1)
                NU                  = tmpPara(2)
                LATTICE             = tmpPara(3)
                totalStress         = tmpPara(4)
                DELTA_SIG           = tmpPara(5)
                A_CUBE              = tmpPara(6)
                MOBILITY            = tmpPara(7)
                DTIME               = tmpPara(8)
                StrainRate          = tmpPara(9)
                ElasticCons         = tmpPara(10)
                dcrit               = tmppara(11)
                APPLIED_SIG                 = tmpPara(12)
                RMIN                        = tmpPara(13)
                totalStrain                     = tmpPara(14)
                totalTime                       = tmpPara(15)
                CycFreq                         = tmpPara(16)
                dpCycF0                         = tmpPara(17)
                dpCycStrain                     = tmpPara(18)
                Loading_Dir%v(1)                = tmpPara(19)
                Loading_Dir%v(2)                = tmpPara(20)
                Loading_Dir%v(3)                = tmpPara(21)
		PEIERLS				= tmpPara(22)
		MOB_R				= tmpPara(23)
		PEI_R				= tmpPara(24)
		ATRatio				= tmpPara(25)
		TRatio				= tmpPara(26)
		NGRatio				= tmpPara(27)
        endif

end subroutine bcastSomeVariables
!**************************************************************
!	BcastStepVariables: 
!		Broadcast some variables that change every timestep.
!**************************************************************
subroutine BcastStepVariables
	use ZQMPI
	use MPIModule,only:iNmaxS
	implicit none

	call MPI_BCAST(iNMaxS,1,MPI_INTEGER,iMaster,MPI_COMM_WORLD,iIerr)

end subroutine BcastStepVariables
!***********************************************************************
!    BcastGlobalTree: 
!		Broadcast the global tree to slaves.
!       Called by DistributeGLobalTree.
!***********************************************************************
subroutine BcastGlobalTree
	use ZQMPI
	use vectors
	use MPIModule,only:tTreeRoot,iNumTNode
	use treeTransfer,only:tmpTTLogNL,tmpTTInts,tmpTTDps
	use TreeMD, only : PackUnPackTree
	implicit none

	integer :: status

	call MPI_Bcast(iNumTNode,1,MPI_INTEGER,iMaster,MPI_COMM_WORLD,iIerr)

	allocate(tmpTTLogNL(iNumTNode))
	allocate(tmpTTInts(5*iNumTNode))
	allocate(tmpTTDps(6*iNumTNode))

	if(iMyid==iMaster)then    ! master: pack the data and send
		call PackUnPackTree(tTreeRoot) 
	endif

	!//////////////////////////////////////////////////////////////////////////
	! Broadcast info
	call MPI_Bcast(tmpTTLogNL(1),iNumTNode,MPI_LOGICAL,iMaster,MPI_COMM_WORLD,iIerr)

	call MPI_Bcast(tmpTTInts,iNumTNode*5,MPI_INTEGER,iMaster,MPI_COMM_WORLD,iIerr)

	call MPI_Bcast(tmpTTDps,iNumTNode*6,MPI_DOUBLE_PRECISION,iMaster,MPI_COMM_WORLD,iIerr)

	if(iMyid/=iMaster)then
		allocate(tTreeRoot)
		tTreeRoot%iIndex=1
		NULLIFY(tTreeRoot%tpChildL)
		NULLIFY(tTreeRoot%tpChildR)
		call PackUnPackTree(tTreeRoot)
	endif

	deallocate(tmpTTLogNL,tmpTTInts,tmpTTDps)
	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

end subroutine BcastGlobalTree
!***************************************************************************
!	BcastGPlanes:
!		Distribute general glide planes to all the slaves.
!		If update, only do this for the newly generated glide planes.
!
!		iNTotalGPs is initially zero. Each time the GPs are distributed, iNTotalGPs
!		will be assigned to the value of iNPlane. If INPlane changes later, the
!		difference betwee iNTotalGPs and INPlanes, i.e., the newly generated planes,
!		will be updated only.
!***************************************************************************  
subroutine BcastGPlanes
	use ZQMPI
	use vectors
	use variables,only:GPlanes,iNPlane
	use MPIModule,only:iNTotalGPs
	use GPTransferMD,only:tmpGPDPs,tmpGPInts
	implicit none

	integer::i,j

	!////////////////////////////////////////////////////
	! Broadcast the number of new planes
	call MPI_BCAST(iNPlane,1,MPI_INTEGER,iMaster,MPI_COMM_WORLD,iIerr)

	if(iNPlane>0)then !update if there is any new plane.
		i=iNPlane
		j=i*3
		allocate(tmpGPDps(j))      !origin
		allocate(tmpGPInts(i*7))   !miller, inumPequiv, and 5 lists

		if(iMyid==iMaster)then
			call PackUnpackGPs(j)  !j changed
		endif

		!//////////////////////////////////////////////////////
		call MPI_Bcast(j,1,MPI_INTEGER,iMaster,MPI_COMM_WORLD,iIerr)
		call MPI_BCAST(tmpGPDPs,i*3,MPI_DOUBLE_PRECISION,iMaster,MPI_COMM_WORLD,iIerr)
		call MPI_BCAST(tmpGPInts,j,MPI_INTEGER,iMaster,MPI_COMM_WORLD,iIerr)
		
		if(iMyid/=iMaster)then
			call PackUnpackGPs(j)
		endif
		deallocate(tmpGPDps,tmpGPInts)

	endif

	iNTotalGPs=iNPlane
	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

end subroutine BcastGPlanes
!****************************************************************************
!	PackUnpackGPs:
!		Pack glide planes for distribution if on master;
!		Unpack glide planes if on slaves.
!		On master:
!****************************************************************************
subroutine PackUnpackGPs(num1)
	use ZQMPI
	use vectors
	use variables,only:iNPlane,GPlanes
	use MPIModule,only:iNTotalGPs
	use GPTransferMD,only:tmpGPDps,tmpGPInts
	implicit none

	integer::i,j,k,num1

	i=iNPlane
	num1=0

	!/////////////////////////////////////////
	!Pack planes together on master.
	if(iMyid==iMaster)then 
		do j=1,i
			tmpGPDPs(j*3-2)	=	GPlanes(j)%Origin%v(1)
			tmpGPDPs(j*3-1)	=	GPlanes(j)%Origin%v(2)
			tmpGPDPs(j*3)		=	GPlanes(j)%Origin%v(3)
			tmpGPInts(num1+1)		=	GPlanes(j)%iMiller
			tmpGPInts(num1+2)		=	GPlanes(j)%iNumPEquiv
			num1=num1+2
			if(GPlanes(j)%iNumPEquiv>0)then
				do k=1,GPlanes(j)%iNumPEquiv
					tmpGPInts(num1+k)=GPlanes(j)%iIDPEquiv(k)
				enddo
				num1=num1+GPlanes(j)%iNumPEquiv
			endif

		enddo
	else !Unpack planes on slaves
		do j=1,i
			GPlanes(j)%Origin%v(1)	=	tmpGPDPs(j*3-2)
			GPlanes(j)%Origin%v(2)	=	tmpGPDPs(j*3-1)
			GPlanes(j)%Origin%v(3)	=	tmpGPDPs(j*3)
			GPlanes(j)%iMiller		=	tmpGPInts(num1+1)
			GPlanes(j)%iNumPEquiv	=	tmpGPInts(num1+2)
			num1=num1+2
			if(GPlanes(j)%iNumPEquiv>0)then
				do k=1,GPlanes(j)%iNumPEquiv
					GPlanes(j)%iIDPEquiv(k)=tmpGPInts(num1+k)
				enddo
				num1=num1+GPlanes(j)%iNumPEquiv
			endif
		enddo
	endif

end subroutine PackUnpackGPs
!*********************************************************************
!	BcastGIndexArray:
!		Broadcast global index array to each processor.
!	
!		Broadcast total number of dislocation points
!*********************************************************************
subroutine BcastGIndexArray
	use ZQMPI
	use vectors
	use MPIModule,only:iNTotalDPs,iaGLIndex
	implicit none

	integer::i,j,status
	integer,allocatable,dimension(:)::tmpIndex

	!///////////////////////////////////////////
	!Broadcast total number of dps
	call MPI_BCAST(iNTotalDPs,1,MPI_INTEGER,iMaster,MPI_COMM_WORLD,iIerr)
	
	i=iNTotalDPs*4
	allocate(tmpIndex(i))

	!construct temp array
	if(iMyid==iMaster)then
		do j=1,iNTotalDPs
			tmpIndex((j-1)*4+1)	=	iaGLIndex(j,1)
			tmpIndex((j-1)*4+2)	=	iaGLIndex(j,2)
			tmpIndex((j-1)*4+3)	=	iaGLIndex(j,3)
			tmpIndex((j-1)*4+4)	=	iaGLIndex(j,4)
		enddo
	endif

	!///////////////////////////////////////////
	!Broadcast the index array
	call MPI_BCAST(tmpIndex,i,MPI_INTEGER,iMaster,MPI_COMM_WORLD,iIerr)

	!construct global index array from tmp array
	if(iMyid/=iMaster)then
		do j=1,iNTotalDPs
			iaGLIndex(j,1)=tmpIndex((j-1)*4+1)
			iaGLIndex(j,2)=tmpIndex((j-1)*4+2)
			iaGLIndex(j,3)=tmpIndex((j-1)*4+3)
			iaGLIndex(j,4)=tmpIndex((j-1)*4+4)
		enddo
	endif

	deallocate(tmpIndex,STAT=status)
	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

end subroutine BCastGIndexArray
!************************************************************************
!	DistributeDPs:
!		Master distribute Dps to local processors.
!
!		Arrays of millers, glide planes, burgers vectors have already
!		been constructed on local processors.
!		
!		Components needed to be transfered:
!			Integers:7 
!					iID,iLoopID,itype,iPlane,iBurgers,iBeginP,iEndP
!
!			Vectors:  5 ->15 double precisions
!					tvPL,tvTL,tvAcce,tvPreV,tvPreVT
!------------------------------------------------------------------------
! Modified: ZQ Wang, 01/26/06
!	Added the tvacce,tvPreV, and dpmass into the list of transfer. Now, we have:
!	Integers: 7
!	Vectors: 5--->15 double precisions  --|
!					                      |===>totally 16
!	Double precision: 1                ---|
!	Also, added two variables, idpCount,iIntCount to remember the numbers.
!
!	02/14/06:
!	Modified the output frequency.
!	06/02/06:
!	5 vectors to transfer
!	
!	9/25/06:
!		lToCal only communicates once here. No need to transfer back
!	after calculation is finished.
!************************************************************************
subroutine DistributeDPs
	use ZQMPI
	use vectors
	use MPIModule,only:tpBoxroots
	use variables,only:DPoints,iloop
	implicit none

	integer :: i,j,k,itmp,idpCount,iIntCount

	integer,dimension(:),pointer :: tmpNumDPs            !number of DPs
	double precision,dimension(:),pointer :: tmpDoubles  !for dp & vector
	integer,dimension(:),pointer :: tmpInts			  !for integer

	! master does
	if(iMyid == iMaster)then
		Do i=1, iNumprocs-1
			! num of subs
			allocate(tmpNumDPs(1))
			tmpNumDPs(1) = tpBoxRoots(i)%p%iNumDPs
			itmp=tmpNumDPs(1)
			if(mod(iloop,5)==0)write(*,*)iMyid,": Number of Dps for Proc ",i," is:",itmp  !2/14/06

			idpCount=16     	! 5 vectors + 1 dp, 06/02/06
			iIntCount=8		! 7 integers, 05/9/06

			allocate(tmpDoubles(int(itmp*idpCount)))  
			allocate(tmpInts(int(itmp*iIntCount)))  

			! assign DP values to tmp arrays
			do j=1,itmp
				!////////////////////
				! Integers
				tmpInts((j-1)*iIntCount+1)	=	tpBoxRoots(i)%p%tapDPs(j)%p%iID
				tmpInts((j-1)*iIntCount+2)	=	tpBoxRoots(i)%p%tapDPs(j)%p%iLoopID
				tmpInts((j-1)*iIntCount+3)	=	tpBoxRoots(i)%p%tapDPs(j)%p%itype
				tmpInts((j-1)*iIntCount+4)	=	tpBoxRoots(i)%p%tapDPs(j)%p%iPlane
				tmpInts((j-1)*iIntCount+5)	=	tpBoxRoots(i)%p%tapDPs(j)%p%iBurgers
				tmpInts((j-1)*iIntCount+6)	=	tpBoxRoots(i)%p%tapDPs(j)%p%iBeginP
				tmpInts((j-1)*iIntCount+7)	=	tpBoxRoots(i)%p%tapDPs(j)%p%iEndP
				if(tpBoxRoots(i)%p%tapDPs(j)%p%lToCal)then
					tmpInts((j-1)*iIntCount+8)	=	1
				else
					tmpInts((j-1)*iIntCount+8)	=	0
				endif

				!//////////////////////////
				! Doubles
				do k=1,3
					tmpDoubles((j-1)*idpCount+k)=tpBoxRoots(i)%p%tapDPs(j)%p%tvPL%v(k)
				enddo
				do k=1,3
					tmpDoubles((j-1)*idpCount+k+3)=tpBoxRoots(i)%p%tapDPs(j)%p%tvTL%v(k)
				enddo
				do k=1,3
					tmpDoubles((j-1)*idpCount+k+6)=tpBoxRoots(i)%p%tapDPs(j)%p%tvAcce%v(k)
				enddo
				do k=1,3
					tmpDoubles((j-1)*idpCount+k+9)=tpBoxRoots(i)%p%tapDPs(j)%p%tvPreV%v(k)
				enddo
				do k=1,3
					tmpDoubles((j-1)*idpCount+k+12)=tpBoxRoots(i)%p%tapDPs(j)%p%tvPreVT%v(k)
				enddo

				tmpDoubles((j-1)*idpCount+16)=tpBoxRoots(i)%p%tapDPs(j)%p%dpMass
								
			enddo
			
			!////////////////////////////////////////
			! Distribute tmp arrays to processors

			! number of dps
			call ZQMPI_Isend(tmpNumDPs,1,i,0)

			! Integers
			call ZQMPI_Isend(tmpInts,itmp*iIntCount,i,0)

			! Doubles
			call ZQMPI_Isend(tmpDoubles,itmp*idpCount,i,0)

			!/////////////////////////////////////////////
			! Nullify variables after distribution for each processor.
			NULLIFY(tmpNumDPs,tmpInts,tmpDoubles) ! nullify these and get ready for next one

		enddo
	endif

end subroutine DistributeDPs

!****************************************************************************
!   ReceiveLocalDps: 
!		Done by slaves.
!		These DPs are real computational DPs. The number is assigned to iNPoint.
!		Later, there will be other Dps, they are connection and neighbor DPs,
!		they are not real computational Dps, just for calculation reference.
!
!		INPoint is refreshed each time step.
!----------------------------------------------------------------------------
! Modified: ZQ Wang, 01/26/06
!	Added tvacce, tvPrev and dpmass to the list of transfer
!
!****************************************************************************
subroutine ReceiveLocalDPs
	use ZQMPI
	use vectors
	use variables,only:iNPoint,DPoints
	use MPIModule,only:iNTotalDPs
	implicit none

	double precision,dimension(:),pointer :: tmpDoubles  !for dp & vector
	integer,dimension(:),pointer :: tmpInts			  !for integer
	integer :: i,j,k,status,iIntCount,iDPCount

	if(iMyid/=iMaster)then
		!////////////////////////////////////////////////////
		! Receive number of DPs
		call MPI_Recv(iNPoint,1,MPI_INTEGER,iMaster,MPI_ANY_TAG,MPI_COMM_WORLD,iStat,iIerr)

		!//////////////////////////////
		! allocate variables
		iIntCount=8    !5/9/06
		iDPcount=16   !05/9/06

		allocate(tmpDoubles(iNPoint*iDPCount),tmpInts(iNPoint*iIntCount))

		!//////////////////////////////
		! receives DPs in tmp arrays
		call MPI_Recv(tmpInts,iNPoint*iIntCount,MPI_INTEGER,iMaster,MPI_ANY_TAG,MPI_COMM_WORLD,iStat,iIerr)
		call MPI_Recv(tmpDoubles,iNpoint*iDPCount,MPI_DOUBLE_PRECISION,iMaster,MPI_ANY_TAG,MPI_COMM_WORLD,iStat,iIerr)

		!/////////////////////////////////
		! assign values to real local DPs
		do j=1,iNPoint
			!////////////////////
			! Integers
			DPoints(j)%iID			=	tmpInts((j-1)*iIntCount+1)
			DPoints(j)%iLoopID		=	tmpInts((j-1)*iIntCount+2)
			DPoints(j)%itype		=	tmpInts((j-1)*iIntCount+3)
			DPoints(j)%iPlane		=	tmpInts((j-1)*iIntCount+4)
			DPoints(j)%iBurgers		=	tmpInts((j-1)*iIntCount+5)
			DPoints(j)%iBeginP		=	tmpInts((j-1)*iIntCount+6)
			DPoints(j)%iEndP		=	tmpInts((j-1)*iIntCount+7)
			DPoints(j)%iNumNei		=	0
			if(tmpInts((j-1)*iIntCount+8)==1)then
				DPoints(j)%lToCal=.true.
			else
				DPoints(j)%lToCal=.false.
			endif
			!//////////////////////////
			! Doubles
			do k=1,3
				DPoints(j)%tvPL%v(k)=tmpDoubles((j-1)*iDPCount+k)
			enddo
			do k=1,3
				DPoints(j)%tvTL%v(k)=tmpDoubles((j-1)*iDPCount+k+3)
			enddo
			do k=1,3
				DPoints(j)%tvAcce%v(k)=tmpDoubles((j-1)*iDPCount+k+6)
			enddo
			do k=1,3
				DPoints(j)%tvPreV%v(k)=tmpDoubles((j-1)*iDPCount+k+9)
			enddo
			do k=1,3
				DPoints(j)%tvPreVT%v(k)=tmpDoubles((j-1)*iDPCount+k+12)
			enddo
			DPoints(j)%dpMass=tmpDoubles((j-1)*iDPCount+16)

			call setGlobal(j)
		enddo

		deallocate(tmpInts,tmpDoubles,STAT=status)
		iNTotalDPs=iNPoint
	endif

end subroutine ReceiveLocalDPs
!******************************************************************************
!    TransGhostTree: 
!		 Build ghost tree and then finish local essential tree,done by slaves.
!        Including information exchange between slaves.
!        Attention: here, tNode should be included. This is different
!        from distribution of global tree. Thus, while building
!        local tree, tNode should be assigned the value of tapsubs
!
!future: add multipole into exchanged info.
!
!			Determination of close box:
!				The distance between the centers of two boxes should be larger
!				than the sum of the following variables: 
!				1.dcrit,
!				2.maxlength of rearrangement of dislocaion line, 
!				3.the sizes of two boxes.
!*******************************************************************************
subroutine TransGhostTree
	use ZQMPI
	use vectors
	use MPIModule,only:tLocalRoot,iNumCloseBox,iIndexOfCloseBox,tpBoxRoots,iNumTNode
	use variables,only:dcrit,iNPoint,iloop
	use TreeMD,only : DetermineGhost, BuildGhost
	use LoopRearrangeMD,only:dpMaxAveLength
	use GhostTransMD,only:iNumOfNodeSend,tmpSendLog,tmpSendInt,tmpSendDps, &
							iNumOfNodeRece,tmpReceLog,tmpReceInt,tmpReceDps
	implicit none

	integer :: status

	integer :: i,j,k,ii,jj,kk,itemp,jtemp,ktemp
	type(vector) :: tempV1,tempV2,tempV3
	double precision :: tmpS1,tmpS2,tmpS3
	double precision :: GetDist

	if(iMyid/=iMaster)then

		! check distance between boxes and establish a list of close boxes.
		tmpS1=MAG(tLocalRoot%tSizeV)
		iNumCloseBox=0
		do i=1,iNumprocs-1
			if(i/=iMyid)then
				tmpS2=MAG(tpBoxRoots(i)%p%tSizev)
				tmpS3=GetDist(tLocalRoot%tCenterv,tpBoxRoots(i)%p%tCenterv)
				!///////////////////////////////////////////////
				! Determine the distance between two boxes.
				if(tmpS3<=tmpS1+tmpS2+dcrit+dpMAXAveLength)then  ! add to list
					iNumCloseBox=iNumCloseBox+1
					iIndexOfCloseBox(iNumCloseBox)=i
				endif
			end if
		end do
	endif

	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

	if(iMyid/=iMaster)then
		!//////////////////////////////////////////
		! send tree to others
		do i=1,iNumCloseBox
			j = iIndexOfCloseBox(i)

			allocate(iNumOfNodeSend(1))
			allocate(tmpSendLog(iNPoint*2+100))
			allocate(tmpSendInt(5*iNPoint*2+100))     !2 means the total number of tree nodes is 2n-1
			allocate(tmpSendDps(6*iNPoint*2+100))

			iNumOfNodeSend(1) = 0
			tmpS1 = MAG(tpBoxRoots(j)%p%tSizev)

			call DetermineGhost(tLocalRoot,j,tmpS1)

			k=iNumOfNodeSend(1)
			call ZQMPI_Isend(iNumOfNodeSend,1,j,0)

			if(k>0)then
				!///////////////////////////////////////////////////////////
				! Communications
				call ZQMPI_Isend(tmpSendLog,k,j,0)

				call ZQMPI_Isend(tmpSendInt,k*5,j,0)

				call ZQMPI_Isend(tmpSendDps,k*6,j,0)

				! Nullify pointers and prepare for next one
				NULLIFY(tmpSendLog,tmpSendInt,tmpSendDps)
			else
				deallocate(tmpSendLog,tmpSendInt,tmpSendDps)
			endif

			Nullify(iNumOfNodeSend)

		enddo

		!/////////////////////////////////////////
		! receive trees from others 
		do i=1,iNumCloseBox

			! call receive....iNumOfNodeRece....length : 1
			call MPI_Recv(iNumOfNodeRece,1,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,iStat,iIerr)
			j=iStat(MPI_SOURCE)

			if(iNumOfNodeRece>0)then
				allocate(tmpReceLog(iNumOfNodeRece))
				allocate(tmpReceInt(5*iNumOfNodeRece))
				allocate(tmpReceDps(6*iNumOfNodeRece))

				! call receive....tmpReceLog....length : iNumofNodeRece
				call MPI_Recv(tmpReceLog,iNumOfNodeRece,MPI_LOGICAL,j,MPI_ANY_TAG,MPI_COMM_WORLD,iStat,iIerr)

				! call receive....tmpReceInt....length : iNumOfNodeRece*8
				call MPI_Recv(tmpReceInt,iNumOfNodeRece*5,MPI_INTEGER,j,MPI_ANY_TAG,MPI_COMM_WORLD,iStat,iIerr)

				! call receive....tmpReceDps....length : iNumOfNodeRece*9
				call MPI_Recv(tmpReceDps,iNumOfNodeRece*6,MPI_DOUBLE_PRECISION,j,MPI_ANY_TAG,MPI_COMM_WORLD,iStat,iIerr)

				!////////////////////////////////////////////
				! Attach ess tree to the appropriate root boxes
				tpBoxRoots(j)%p%iChild(1) = iNumTNode+1
				tpBoxRoots(j)%p%iChild(2) = iNumTNode+2
				tpBoxRoots(j)%p%lNodeOrLeaf = .true.
				
				call BuildGhost(tpBoxRoots(j)%p)
				
				iNumTNode = iNumTNode + iNumOfNodeRece
				
				deallocate(tmpReceLog,tmpReceInt,tmpReceDps,STAT=status)
			endif
		enddo
	endif
	
	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

end subroutine TransGhostTree
!******************************************************************************
!  UpdateStrainStress: called by UpdateMechanics. Done by everybody
!******************************************************************************
subroutine UpdateStrainStress1
        use ZQMPI
        use vectors
        use MPIModule
        use variables,only:DeltaStrain,SIG_APP,ElasticCons,Dtime,N_TIMES,STRAINRate,MU,zero,SR1,SR2, &
                                        StrainRate_P,StrainCurrent,totalStrain,totalTime,Delta_SIG,Load, &
                                        loading_ES,totalStress,ILOOP
        use FunctionMD,only:DpFCyclicForce
        implicit none

        double precision :: tmpDp
        double precision::hertmp
        double precision::sigtmp
!       double precision::dpFCyclicForce
        double precision::dpTmp(2)
        type(matrix)::tmpSig,trans_ES
        integer::i,j,l,m

        ! update global strain on master: tmpDp contains total plastic strain
        call MPI_ALLReduce(DeltaStrain,tmpDp,1,MPI_DOUBLE_PRECISION,MPI_SUM,  &
                                                        MPI_COMM_WORLD,iIerr)

        if(iMyid==iMaster)then
                StrainRate_p=tmpDp/(N_times*DTIME)
                strainCurrent = strainCurrent+tmpDp
                TotalStrain=TotalStrain+strainRate*N_Times*Dtime
                dpTmp(1)=totalStrain
                dpTmp(2)=StrainRate_P
        endif

        call MPI_BCAST(dpTmp(1),2,MPI_DOUBLE_PRECISION,iMaster,MPI_COMM_WORLD,iIerr)
        if(iMyid/=iMaster)then
                totalStrain=dpTmp(1)
                StrainRate_P=dpTmp(2)
        endif

        SR1=totalStress

        if(Load==1)then ! update stress ---from strain rate
               if(ILOOP<=400 .or. ILOOP>1200)then ! cyclic loading positive strain rate
                StrainRate=1D2
               elseif(ILOOP>400 .and. ILOOP<=1200)then ! cyclic loading negative strainrate
                StrainRate=-1D2
               endif
write(*,*) "Strainrate= ",StrainRate, "Iloop_time= ",ILOOP               
totalStress = totalStress     &
                                                +ElasticCons*Dtime*N_TIMES*(StrainRate-StrainRate_p)/Mu
        elseif(Load==2)then ! update stress----from delta_sig
                totalStress = totalStress+Delta_Sig/MU

        elseif(Load==3)then ! update stress----from real value
                totalStress = 0
        elseif(Load==4)then ! update stress-----Cyclic
                totalStress = dpFCyclicForce()
        endif

        SR2=totalStress

        !Calculate the stress status in the global coordinate system
        trans_es=Trans(Loading_ES)
	tmpSig=zero
        tmpSig%v(1)%v(1)=totalStress
        SIG_APP=zero
        do i=1,3;do j=1,3
             do l=1,3;do m=1,3
                 Sig_APP%v(i)%v(j)=Sig_APP%v(i)%v(j)+Trans_ES%v(i)%v(l)*Trans_ES%v(j)%v(m)*tmpSig%v(l)%v(m)
             enddo;enddo
        enddo;enddo

end subroutine UpdateStrainStress1
!******************************************************************************
!  UpdateStrainStress: called by UpdateMechanics. Done by everybody
!	02/06/2009: Temp to correct stress calculation.
!			totalstress, totalstrain variables don't change.
!	02/17/2009: Temp to calculate Uniaxial stress loading:
!			Plastic strain from dislocations;
!			Elastic strain from loading stresses, only sigma_11 not zero;
!			Total strain is from both plastic and elastic strain and then stress-strain curve known.
!	04/13/2009: Updated to calculate uniaxial stress from strain rate loading
!******************************************************************************
subroutine UpdateStrainStress(tmpStrain)
        use ZQMPI
        use vectors
        use MPIModule
        use variables,only:DeltaStrain,SIG_APP,ElasticCons,Dtime,N_TIMES,STRAINRate,MU,zero,SR1,SR2, &
                                        StrainRate_P,StrainCurrent,totalStrain,totalTime,Delta_SIG,Load, &
                                        loading_ES,totalStress,NU
        use FunctionMD,only:DpFCyclicForce
        implicit none

        double precision :: tmpDp
        double precision::hertmp
        double precision::sigtmp
!       double precision::dpFCyclicForce
        type(matrix)::tmpSig,trans_ES
	type(matrix)::tmpStrain          !strain matrix in loading coordinate system
	type(matrix)::localStrainRate_p,localStrain_p,localStrainRate,localStrain_e,localStress
        integer::i,j,l,m
	double precision::lamda,E_O_M, E1,E2 !E/Mu, E1 and E2 are used to calculate strain from stress
	double precision::epsilon_kk

        ! update global strain on master: tmpDp contains total plastic strain
	!Elastic modulus
	lamda=nu/(1+nu)/(1-2*nu)  ! it is the coefficient to calculate lamda from E, lamda=E*lamda, EOM will be normalized by MU,
                                  ! as below.
	E_O_M=ElasticCons/MU
	E1=0.5/(1+nu)      !simplified (lamda+mu)/(3*lamda+2*mu), these two are used to calculate strain from stress
	E2=-0.5*nu/(1+nu)  !simplified -1/2*lamda/2/(3*lamda+2*mu)

	!/////////////////////////////////////
	!plastic strain is always needed for load 1 or 2
	!----total plastic strain-----
	do i=1,3
		do j=1,3
			if(iMyid==iMaster)then
				DeltaStrain=0.0
			else
				DeltaStrain=tmpStrain%v(i)%v(j)
			endif
			call MPI_ALLReduce(DeltaStrain,tmpDp,1,MPI_DOUBLE_PRECISION,MPI_SUM,  &
                                                        MPI_COMM_WORLD,iIerr)
			localStrain_p%v(i)%v(j)=tmpDP
		enddo
	enddo
	!-----total plastic strain rate-------
        StrainRate_p=localStrain_p%v(1)%v(1)/(N_times*DTIME)
        strainCurrent = strainCurrent+localStrain_p%v(1)%v(1)
	localStrainRate_p=(1.0/(N_Times*DTime))*LocalStrain_P
	!/////////////////////////////////////////

        SR1=totalStress

        if(Load==1)then ! update stress ---from strain rate!  Total strain is from stran rate
		!----Total applied strain, local and global forms------
		TotalStrain=TotalStrain+strainRate*N_Times*Dtime
        	localStrainRate=zero
        	localStrainRate%v(1)%v(1)=strainRate

		!Elastic strain changes
		localStrain_e=Dtime*N_Times*(localStrainRate-localStrainRate_p)

		!stress from hook's law: local stress change due to elastic strain change
		localStress=zero
		!sigma_11
                localStress%v(1)%v(1) = E_O_M*localStrain_e%v(1)%v(1)

		!total stress_11
		totalStress=totalStress+localStress%v(1)%v(1)

        elseif(Load==2)then ! update stress----from delta_sig
		!Stress from step loading
		localStress=zero
!		if(totalStress>=250d6/MU)then
!			localStress%v(1)%v(1)=0.0
!		else
			localStress%v(1)%v(1)=Delta_Sig/MU
!		endif
                totalStress = totalStress+localStress%v(1)%v(1)

		!First calculate elastic strain based on current stress changes
		localStrain_e=zero
		localStrain_e%v(1)%v(1)=E1*localStress%v(1)%v(1)
		localStrain_e%v(2)%v(2)=E2*localStress%v(1)%v(1)
		localStrain_e%v(3)%v(3)=E2*localStress%v(1)%v(1)

		!Then, total strain
		totalStrain=totalStrain+localStrain_e%v(1)%v(1)+localStrain_p%v(1)%v(1)

        elseif(Load==3)then ! update stress----from real value
                totalStress = 0
        elseif(Load==4)then ! update stress-----Cyclic
                totalStress = dpFCyclicForce()
        endif

        SR2=totalStress

        !Calculate the stress status in the global coordinate system
        trans_es=Trans(Loading_ES)
        tmpSig=zero
	do i=1,3;do j=1,3
             do l=1,3;do m=1,3
                 tmpSig%v(i)%v(j)=tmpSig%v(i)%v(j)+Trans_ES%v(i)%v(l)*Trans_ES%v(j)%v(m)*localStress%v(l)%v(m)
             enddo;enddo
        enddo;enddo

	Sig_APP=Sig_APP+tmpSIG

end subroutine UpdateStrainStress
!**************************************************************************
! AddDensityTogether: 
!	Combine local densities together.
!**************************************************************************
subroutine AddDensityTogether
	use ZQMPI
	use MPIModule,only:dTmpDensi,dDensity,dTmpSDensi,dTmpEDensi,dScrewDensity,dEdgeMixDensity,iIerr
	implicit none

	call MPI_AllReduce(dTmpDensi,dDensity,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
									MPI_COMM_WORLD,iIerr)

	call MPI_AllReduce(dTmpSDensi,dScrewDensity,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
									MPI_COMM_WORLD,iIerr)
	call MPI_AllReduce(dTmpEDensi,dEdgeMixDensity,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
									MPI_COMM_WORLD,iIerr)


end subroutine AddDensityTogether
!******************************************************************************
!	CollectDPsFrom:
!		Used in update, master collect Dps from slaves.
!	On Master:
!		iNtotalPlanes=iNPlane after distribution of glide planes. iNPlane
!	may change in this update process. iNTotalPlanes will not change during 
!	the running of the program.
!
!	On Slaves:
!		iNTotalPlanes=iNPlane at the initial receiving of the planes. During
!	the running of the program, iNPlane will change. But iNTotalPlanes will not
!	change until the distribution of the iNPlane.
!
!	Process to update DPs:
!		First, to establish a connection between the master and one slave point.
!		1. Master receive new planes, update iNPlane and Gplanes.
!		2. Master send back current iNPlane such that the slave knows the index of
!	the new planes.
!		3. Slaves begin to group DPs and update glide plane information.
!		4. Slaves begin to send pointns back to master.56wexcvbnmz/. 
!--------------------------------------------------------------------------
! Modified: ZQ Wang, 01/26/06
!    Added tvAcce, tvPreV and dpMass.
!		Added idpCount to reprent the number of DPs in transfer.
! 
!******************************************************************************
subroutine CollectDPsFrom
	use ZQMPI
	use variables,only:iNPoint,iNPlane,GPlanes,DPoints,MAX_NEIGHBOR
	use MPIMOdule,only:iNTotalGPs
	implicit none

	integer::i,j,k,l,m,n,iSrc,idpCount
	integer::iNumNewPl,itmpNPoint,itmpNum,iNumTotalNei
	integer,dimension(:),pointer::iTmpRecv  !
	integer,dimension(:),pointer::iTmpRecvNei
	double precision,dimension(:),pointer::dpTmpRecv !

	do i=1,iNumProcs-1
		!////////////////////////////////////////////////
		! Establish a connection between master and slave
		call MPI_Recv(itmpNum,1,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,iStat,iIerr)
		iTmpNPoint=itmpNum ! number of point from slave processors
		iSrc=iStat(MPI_SOURCE)

		if(iTmpNPoint > 0)then
			!////////////////////////////////////////////////////////
			! Master receive DPs from slave processors
			iDPCount=25         !04/20/06          !31      !01/26/06
			allocate(iTmpRecv(iTmpNPoint*3+1),iTmpRecvNei(iTmpNPoint*MAX_NEIGHBOR), &
				dpTmpRecv(iTmpNPoint*iDPCount))

			call MPI_Recv(iTmpRecv,iTmpNPoint*3+1,MPI_INTEGER,iSrc,MPI_ANY_TAG,MPI_COMM_WORLD,iStat,iIerr)
			iNumTotalNei=iTmpRecv(iTmpNPoint*3+1)

			call MPI_Recv(dpTmpRecv,iTmpNPoint*iDPCount,MPI_DOUBLE_PRECISION,iSrc,MPI_ANY_TAG,MPI_COMM_WORLD,iStat,iIerr)

			if(iNumTotalNei>0)then
				call MPI_Recv(iTmpRecvNei,iNumTotalNei,MPI_INTEGER,iSrc,MPI_ANY_TAG,MPI_COMM_WORLD,iStat,iIerr)
			endif

			!///////////////////////////////////////////////
			! Begin to update points
			iNumTotalNei=0
			do j=1,iTmpNPOint
				k					=	iTmpRecv((j-1)*3+1)
				DPoints(k)%iPlane	=	iTmpRecv((j-1)*3+2)
				DPoints(k)%iNumNei	=	iTmpRecv((j-1)*3+3)

				do l=1,3
					DPoints(k)%tvPL%v(l)	=	dpTmpRecv((j-1)*iDPCount+l)
				enddo
				do l=1,3
					DPoints(k)%tvTL%v(l)	=	dpTmpRecv((j-1)*iDPCount+l+3)
				enddo
				do l=1,3
					DPoints(k)%tvPG%v(l)	=	dpTmpRecv((j-1)*iDPCount+l+6)
				enddo
				do l=1,3
					DPoints(k)%tvTG%v(l)	=	dpTmpRecv((j-1)*iDPCount+l+9)
				enddo
				do l=1,3
					DPoints(k)%force%v(l)	=	dpTmpRecv((j-1)*iDPCount+l+12)
				enddo
				do l=1,3
					DPoints(k)%tvAcce%v(l)	=	dpTmpRecv((j-1)*iDPCount+l+15)
				enddo
				do l=1,3
					DPoints(k)%tvPreV%v(l)	=	dpTmpRecv((j-1)*iDPCount+l+18)
				enddo
				do l=1,3
                                        DPoints(k)%tvPreVT%v(l)  =       dpTmpRecv((j-1)*iDPCount+l+21)
                                enddo
				
				DPoints(k)%dpMass	=	dpTmpRecv((j-1)*iDPCount+iDPCount)
				!///////////////////////////////
				! Neighbor list
				do l=1,DPoints(k)%iNumNei
					DPoints(k)%iaIDOfNei(l)=iTmpRecvNei(iNumTotalNei+l)   !already in global
				enddo
				iNumTotalNei=iNumTotalNei+DPoints(k)%iNumNei
			enddo
			deallocate(iTmpRecv,iTmpRecvNei,dpTmpRecv)
		endif
	enddo

end subroutine CollectDPsFrom
!*****************************************************************************
!	SendDpsToMaster:
!		Used in update, Slaves send DPs to master.
!-----------------------------------------------------------------------------
! Modified: ZQ Wang,01/26/06
!	Added dpmass, tvPreV and tvacce.
!*****************************************************************************
subroutine SendDPsToMaster
	use ZQMPI
	use variables,only:iNPoint,iNPlane,GPlanes,DPoints,MAX_NEIGHBOR
	use MPIMOdule,only:iNTotalGPs,iNTotalDPs
	implicit none

	integer::i,j,k,iDPCount
	integer::iNumTotalNei
	integer,dimension(:),pointer::itmpSendNei
	integer,dimension(:),pointer::iTmpSend
	double precision,dimension(:),pointer::dpTmpSend

	!/////////////////////////////////////////////////////
	! Send a number of new planes to the master.
	allocate(iTmpSend(1))
	itmpSend(1)=iNPoint
	call ZQMPI_Isend(itmpSend,1,iMaster,0)
	nullify(iTmpSend)

	!////////////////////////////////////////////////////
	!Update planes to master if any new.
	if(iNPoint>0)then
		!/////////////////////////////////////////////////////////////////
		!	 Begin to group DPs and update glide plane number.
		!	 Connection index does not need to be updated here. They only
		!	 change in short range inter.
		!	Integer::iID,iPlane,iNumNei; iBurgers,iLoopID and itype do not change here, actually iBeginP and iEndP for a frank-read source do not change either.
		!	Vector:: tvPL,tvPL2,tvTL,tvTL2,tvPG,tvTG,tvTG2,force
		!	Integer::Neighbor List
		!   
		!	01/26/06:
		!	Vector:: tvAcce
		!	Vector:: tvPreV
		!	DP:		dpMass
		!   Total vector: 10===>30 dp+1dp===>31 dp
		!	
		iDPCount=25   !04/20/06
		allocate(iTmpSend(iNPoint*3+1),dpTmpSend(iNPoint*iDPCount),iTmpSendNei(iNPoint*MAX_NEIGHBOR))
		iNumTotalNei=0
		do i=1,iNPoint
			iTmpSend((i-1)*3+1)=DPoints(i)%iID
			iTmpSend((i-1)*3+2)=DPoints(i)%iPlane
			iTmpSend((i-1)*3+3)=DPoints(i)%iNumNei
		
			!//////////////////////////////////////////
			! Neighbor list
			do j=1,DPoints(i)%iNumNei
				iTmpSendNei(iNumTotalNei+j)=DPoints(DPoints(i)%iaIDOfNei(j))%iID !global
			enddo
			iNumTotalNei=iNumTotalNei+DPoints(i)%iNumNei

			!//////////////////////////////
			! vectors
			do j=1,3
				dpTmpSend((i-1)*iDPCount+j)=DPoints(i)%tvPL%v(j)
			enddo
			do j=1,3
				dpTmpSend((i-1)*iDPCount+j+3)=DPoints(i)%tvTL%v(j)
			enddo
			do j=1,3
				dpTmpSend((i-1)*iDPCount+j+6)=DPoints(i)%tvPG%v(j)
			enddo
			do j=1,3
				dpTmpSend((i-1)*iDPCount+j+9)=DPoints(i)%tvTG%v(j)
			enddo
			do j=1,3
				dpTmpSend((i-1)*iDPCount+j+12)=DPoints(i)%force%v(j)
			enddo
			do j=1,3
				dpTmpSend((i-1)*iDPCount+j+15)=DPoints(i)%tvAcce%v(j)
			enddo
			do j=1,3
				dpTmpSend((i-1)*iDPCount+j+18)=DPoints(i)%tvPreV%v(j)
			enddo
			do j=1,3
                                dpTmpSend((i-1)*iDPCount+j+21)=DPoints(i)%tvPreVT%v(j)
                        enddo
			dpTmpSend((i-1)*iDPCount+iDPCount)=DPoints(i)%dpMass
		enddo
		!////////////////////////////////////
		! Send info to master
		iTmpSend(iNPoint*3+1)=iNumTotalNei
		call ZQMPI_Isend(iTmpSend,iNPoint*3+1,iMaster,0)
		call ZQMPI_ISend(dpTmpSend,iNPoint*iDPCount,iMaster,0)
		if(iNumTotalNei>0)then
			call ZQMPI_ISend(iTmpSendNei,iNumTotalNei,iMaster,0)
		endif
		nullify(iTmpSend,dpTmpSend,iTmpSendNei)
	endif
end subroutine SendDPsToMaster
!*****************************************************************
! UpdateGlidePlanes:
!		Called by dynamics. Used to update glide plane.
!		Master receives from slaves a flag indicator to see
!		if new glide planes generated.
!
!		Plane indices of the DPs are also updated.
!		INTotalGPs on both master and slave are both assigned to iNPlane
!		after the broadcast of the new planes.
!
!		Update list of equivalent planes:
!		1. if new plane, add this plane and get its index, at the same time,
!			construct its list of equivalent planes.
!		2. Update plane index and list on each local processors.
!		3. Combine lists from each processor for old planes, new plane's list
!			has been constructed in step 1.
!*****************************************************************
subroutine UpdateGlidePlanes
	use ZQMPI
	use vectors
	use MPIModule,only:iNTotalGPs
	use Variables,only:iNPlane,iNPoint,DPoints,GPlanes,Miller
	implicit none

	integer::i,j,k,itmpSr,ii,jj,kk
	integer::iTmpNumNewPL,iTmpNumIndex
	integer,dimension(:),pointer::iNumNewPlane
	integer,dimension(:),pointer::iTmpSend,iTmpRecv
	integer,dimension(:),pointer::iTmpIndex
	double precision,dimension(:),pointer::dpTmpSend,dpTmpRecv
	type(vector)::tmpVector
	logical::lFlag

	logical::lInEquivPList

	!////////////////////////////////////
	if(iMyid==iMaster)then
		do i=1,iNumProcs-1
			!////////////////////////////////////////////////////////////////////////////////
			! Receive new plane numbers from slaves
			call MPI_Recv(iTmpNumNewPL,1,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,iStat,iIerr)
			itmpSr=iStat(MPI_SOURCE)

			if(iTmpNumNEwPL>0)then
				!//////////////////////////////////////////////////////////////////////
				! Receive plane info from slaves.
				! iTmpSend here is the new index for new planes.

				allocate(iTmpRecv(iTmpNumNewPL),dpTmpRecv(iTmpNumNewPL*3))
				allocate(iTMpSend(itmpNumNewPL),iTmpIndex(itmpNumNewPL))
				call MPI_Recv(iTmpRecv,iTmpNumNEwPL,MPI_INTEGER,itmpSr,MPI_ANY_TAG,MPI_COMM_WORLD,iStat,iIerr)
				call MPI_Recv(dpTmpRecv,iTmpNumNEwPL*3,MPI_DOUBLE_PRECISION,itmpSr,MPI_ANY_TAG,MPI_COMM_WORLD,iStat,iIerr)
				do j=1,iTmpNumNewPL
					tmpVector%v(1)=dpTmpRecv((j-1)*3+1)
					tmpVector%v(2)=dpTmpRecv((j-1)*3+2)
					tmpVector%v(3)=dpTmpRecv((j-1)*3+3)
					lFlag=.false.
					do k=iNTotalGPs+1,iNPlane
						!if the plane is already generated from other processors
						if(GPlanes(k)%iMiller==iTmpRecv(j) .and.  &
					dabs((tmpVector-GPlanes(k)%Origin)*UV(Miller(iTmpRecv(j))%v))<1.d0)then
							iTmpSend(j)=k
							lFlag=.true.
							exit
						endif
					enddo
					!the plane is not in the list yet
					if(.not. lFlag)then
						iNPlane=iNPlane+1
						GPlanes(iNPlane)%Origin=tmpVector
						GPlanes(iNPlane)%iMiller=iTmpRecv(j)
						iTmpSend(j)=iNPlane
					endif
					iTmpIndex(j)=iTmpSend(j)
				enddo
				deallocate(iTmpRecv,dpTmpRecv)


				!////////////////////////////////////////////////////
				! Send update info to slave, i.e. new index of plane
				call ZQMPI_Isend(iTmpSend,iTmpNumNewPL,itmpSr,0)
				nullify(iTmpSend)

				!////////////////////////////////////////////////////////
				! Receive list of equivalent planes for these new planes
				call MPI_Recv(iTmpNumIndex,1,MPI_INTEGER,itmpSr,MPI_ANY_TAG,MPI_COMM_WORLD,iStat,iIerr)
				allocate(itmpRecv(iTmpNumIndex))
				call MPI_Recv(iTmpRecv,iTmpNumIndex,MPI_INTEGER,itmpSr,MPI_ANY_TAG,MPI_COMM_WORLD,iStat,iIerr)

				! update list for new planes
				k=0
				do j=1,iTmpNumNewPL
					k=k+1
					jj=iTmpIndex(j)
					if(iTmpRecv(k)>0)then
						do ii=1,iTmpRecv(k)
							kk=iTmpRecv(k+ii)
							if(.not.lInEquivPList(jj,kk))call AddToEquivPList(jj,kk)
							if(.not.lInEquivPList(kk,jj))call AddToEquivPList(kk,jj)
						enddo
					endif
					k=k+iTmpRecv(k)
				enddo

				deallocate(itmpRecv,itmpIndex)
			endif

			!/////////////////////////////////////////////////
			! Update list of equivalent planes for old planes
			allocate(itmpRecv(iNTotalGPs*6))
			call MPI_Recv(ii,1,MPI_INTEGER,itmpSr,MPI_ANY_TAG,MPI_COMM_WORLD,iStat,iIerr)
			call MPI_Recv(iTmpRecv,ii,MPI_INTEGER,itmpSr,MPI_ANY_TAG,MPI_COMM_WORLD,iStat,iIerr)
			ii=0
			do j=1,iNTotalGPs
				ii=ii+1
				if(iTmpRecv(ii)>0)then
					do k=1,iTmpRecv(ii)
						jj=iTmpRecv(ii+k)
						if(.not.lInEquivPList(j,jj))call AddToEquivPList(j,jj)
					enddo
					ii=ii+iTmpRecv(ii)
				endif
			enddo
			deallocate(itmpRecv)
		enddo
	!///////////////////
	!Slave:
	else
		!/////////////////////////////
		! Send new plane number
		allocate(iNumNewPlane(1))
		iNumNewPlane(1)=iNPlane-iNTotalGPs
		iTmpNumNEwPL=iNumNewPlane(1)
		call ZQMPI_Isend(iNumNewPlane,1,iMaster,0)
		nullify(iNumNewPlane)
		
		!/////////////////////////////////
		! Send new plane info to master
		if(iTmpNumNewPL>0)then
			allocate(iTmpSend(iTmpNumNEwPL),dpTmpSend(iTmpNumNewPL*3))
			allocate(iTmpRecv(iTmpNumNewPL),iTmpIndex(6*iTmpNumNewPL))
			do i=1,iTmpNumNEwPL
				do j=1,3
					dpTmpSend((i-1)*3+j)=GPlanes(iNTotalGPs+i)%Origin%v(j)
				enddo
				iTmpSend(i)=GPlanes(iNTotalGPs+i)%iMiller
			enddo
			call ZQMPI_Isend(itmpSend,iTmpNumNewPL,iMaster,0)
			call ZQMPI_Isend(dpTmpSend,iTmpNumNewPL*3,iMaster,0)
			nullify(itmpSend,dpTmpSend)

			!////////////////////////////////////////
			! Try to update plane indices of points
			call MPI_Recv(iTmpRecv,iTmpNumNewPL,MPI_INTEGER,iMaster,MPI_ANY_TAG,MPI_COMM_WORLD,iStat,iIerr)
			do j=1,iNPoint
				if(Dpoints(j)%iPlane>iNTotalGPs)then
					DPoints(j)%iPlane=iTmpRecv(DPoints(j)%iPlane-iNTotalGPs)
				endif
			enddo
			do j=1,iNPlane
				do k=1,GPlanes(j)%iNumPEquiv
					if(GPlanes(j)%iIDPEquiv(k)>iNTotalGPs)then
						GPlanes(j)%iIDPEquiv(k)=iTmpRecv(GPlanes(j)%iIDPEquiv(k)-iNTotalGPs)
					endif
				enddo
			enddo
			deallocate(iTmpRecv)

			!//////////////////////////////////////////////
			! Arrange list of equivalent planes to send
			allocate(itmpSend(1))
			ii=0
			do j=1,iTmpNumNewPL
				ii=ii+1
				jj=iNTotalGPs+j
				itmpIndex(ii)=GPlanes(jj)%iNumPEquiv
				do i=1,GPlanes(jj)%iNumPEquiv
					ii=ii+1
					itmpIndex(ii)=GPlanes(jj)%iIDPEquiv(i)
				enddo
			enddo
			itmpSend(1)=ii
			call ZQMPI_ISend(iTmpSend,1,iMaster,0)
			call ZQMPI_Isend(iTmpIndex,ii,iMaster,0)
			nullify(iTmpSend,iTmpIndex)
		endif

		!///////////////////////////////////////////////////////////
		! Arrange list of equivalent planes for old planes to send
		allocate(itmpIndex(iNTotalGPs*6),itmpSend(1))
		ii=0
		do i=1,iNTotalGPs
			ii=ii+1
			itmpIndex(ii)=GPlanes(i)%iNumPEquiv
			if(itmpIndex(ii)>0)then
				do j=1,itmpIndex(ii)
					itmpIndex(ii+j)=GPlanes(i)%iIDPEquiv(j)
				enddo
				ii=ii+GPlanes(i)%iNumPEquiv
			endif
		enddo
		itmpSend(1)=ii
		call ZQMPI_ISend(itmpSend,1,iMaster,0)
		call ZQMPI_ISend(itmpIndex,ii,iMaster,0)
		nullify(iTmpIndex,itmpSend)
	endif

	call MPI_Barrier(MPI_COMM_WORLD,iIerr)	
	!///////////////////////////////////////////////////
	! broadcast new planes from master to all slaves.
	allocate(iTmpSend(1))
	if(iMyid==iMaster)then
		iTmpSend(1)=iNPlane-iNTotalGPs
	endif
	call MPI_BCast(iTmpSend,1,MPI_INTEGER,iMaster,MPI_COMM_WORLD,iIerr)
	iTmpNumNewPL=iTmpSend(1)
	deallocate(iTmpSend)
	
	if(iTmpNumNewPL>0)then
		allocate(iTmpsend(iTmpNumNewPL),dpTmpSend(iTmpNumNewPL*3))
		if(iMyid==iMaster)then
			do j=1,iTmpNumNewPL
				iTmpSend(j)=GPlanes(iNTotalGPs+j)%iMiller
				do k=1,3
					dpTmpSend((j-1)*3+k)=GPlanes(iNTotalGPs+j)%Origin%v(k)
				enddo
			enddo
		endif

		call MPI_BCast(iTmpSend,iTmpNumNewPL,MPI_INTEGER,iMaster,MPI_COMM_WORLD,iIerr)
		call MPI_BCast(dpTmpSend,iTmpNumNewPL*3,MPI_DOUBLE_PRECISION,iMaster,MPI_COMM_WORLD,iIerr)

		if(iMyid/=iMaster)then
			do j=1,iTmpNumNewPL
				GPlanes(iNTotalGPs+j)%iMiller=iTmpSend(j)
				do k=1,3
					GPlanes(iNTotalGPs+j)%Origin%V(K)=dpTmpSend((j-1)*3+k)
				enddo
			enddo
			iNPlane=iNTotalGPs+iTmpNumNewPL
		endif
		deallocate(iTmpSend,dpTmpSend)
	endif
	
	!/////////////////////////////////////////
	! broadcast lists of equivalent planes
	allocate(iTmpSend(iNPlane*6))
	ii=0
	if(iMyid==iMaster)then
		do i=1,iNPlane
			ii=ii+1
			iTmpSend(ii)=GPlanes(i)%iNumPEquiv
			if(GPlanes(i)%iNumPEquiv>0)then
				do j=1,GPlanes(i)%iNumPEquiv
					ii=ii+1
					iTmpSend(ii)=GPlanes(i)%iIDPEquiv(j)
				enddo
			endif
		enddo
	endif

	call MPI_BCast(ii,1,MPI_INTEGER,iMaster,MPI_COMM_WORLD,iIerr)

	call MPI_BCast(iTmpSend,ii,MPI_INTEGER,iMaster,MPI_COMM_WORLD,iIerr)

	if(iMyid/=iMaster)then
		ii=0
		do i=1,iNPlane
			ii=ii+1
			GPlanes(i)%iNumPEquiv=iTmpSend(ii)
			if(GPlanes(i)%iNumPEquiv>0)then
				do j=1,GPlanes(i)%iNumPEquiv
					ii=ii+1
					GPlanes(i)%iIDPEquiv(j)=iTmpSend(ii)
				enddo
			endif
		enddo
	endif
	deallocate(iTmpSend)

	!//////////////////////////////
	! for both master and slaves.
	iNTotalGPs=iNPlane
	
end subroutine UpdateGlidePlanes
!**************************************************************
!	CommConnPoints:
!		Only communicate the connection points.
!**************************************************************
subroutine CommConnPoints
	use ZQMPI
	use vectors
	use variables,only:iNPoint,DPoints
	use MPIModule,only:iaGLIndex,iNumTransP,iaTransList,iNTotalDPs,iNVLPoints
	implicit none

	integer::i,j,k,i1,iToSend,iNumToRecv,iSrc,iDPSCount,iIntSCount
	integer::iNumBoxChange
	integer,dimension(:),pointer::iNumToSend     			!number of connection points needed.
	integer,dimension(:),pointer::iIntToSend			!index of connections needed and sent to other source processors
	double precision,dimension(:),pointer::dpToSend
	integer,dimension(:),allocatable::iTmp				!
	double precision,dimension(:),allocatable::dpTmp
	integer,dimension(:,:),allocatable::iNumIndex				!number of connections need to be sent to other destination processors
	integer,dimension(:,:),allocatable::iIndex			!index of connection points need to be sent, at most all points need to be sent

	!////////////////////////////////////
	!arrange indices
	if(iMyid/=iMaster)then

		allocate(iTmp(iNPoint*2))          !index of the connection points needed to be transferred,at most iNumTransP, which is at most 2 for
							!each point. Could just allocate iNumTransP, but later we will use the same array to send points
							!to others, this number may be large than iNumTransP, so we just allocate the largest number for 
							!both case.
		iNumBoxChange=0

		!Send list
		do i=1,iNumProcs-1

			if(i/=iMyid)then
				!/////////////////////////////////
				!determine the points to request
				iToSend=0
				do j=1,iNumTransP
					if(iaGLIndex(iaTransList(j),1)==i)then
						iToSend=iToSend+1
						iTmp(iToSend)=iaTransList(j)
					endif
				enddo
				
				!....send number and points.....
				if(iToSend>0)then
					iNumBoxChange=iNumBoxChange+1

					!number
					allocate(iNumToSend(1))
					iNumToSend(1)=iToSend
					call ZQMPI_ISend(iNumToSend,1,i,0)
					nullify(iNumToSend)

					!points
					allocate(iIntToSend(iToSend))
					do j=1,iToSend
						iIntToSend(j)=iTmp(j)
					enddo
					call ZQMPI_ISend(iIntToSend,iToSend,i,0)
					nullify(iIntToSend)

				endif

			endif

		enddo

		!....receive index.....index of points that this processor needs to send
		if(iNumBoxChange>0)then

			allocate(iNumIndex(iNumBoxChange,2),iIndex(iNumBoxChange,iNPoint))  !iIndex, at most iNpoint need to be sent		
			do i=1,iNumBoxChange
				call MPI_Recv(iToSend,1,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,iStat,iIerr)
				iSrc=iStat(MPI_SOURCE)
				iNumIndex(i,1)=iSrc
				iNumIndex(i,2)=iToSend
				
				if(iToSend>0)then
					call MPI_Recv(iTmp,iToSend,MPI_INTEGER,iSrc,MPI_ANY_TAG,MPI_COMM_WORLD,iStat,iIerr)
					do j=1,iToSend
						iIndex(i,j)=iTmp(j)
					enddo
				endif
			enddo

		endif
	endif

	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

	!//////////////////////////////////////////////
	!Indices arranged, begin to send points
	if(iMyid/=iMaster)then
		!.....Send To Others.....
		iDPSCount=16
		iIntSCount=7
		if(iNumBoxChange>0)then
			do i=1,iNumBoxChange

				iSrc=iNumIndex(i,1)
				iToSend=iNumIndex(i,2)

				if(iToSend>0)then
			
					allocate(dpToSend(iToSend*iDPSCount),iIntToSend(iToSend*iIntSCount),iNumToSend(1))
					iNumToSend(1)=iToSend

					do j=1,iToSend

						i1=iaGLIndex(iIndex(i,j),2)
						iIntToSend((j-1)*7+1)=DPoints(i1)%iPlane
						iIntToSend((j-1)*7+2)=DPoints(i1)%iID
						iIntToSend((j-1)*7+3)=DPoints(i1)%iLoopID
						iIntToSend((j-1)*7+4)=DPoints(i1)%iType
						iIntToSend((j-1)*7+5)=DPoints(i1)%iBurgers
						iIntToSend((j-1)*7+6)=DPoints(i1)%iBeginP
						iIntToSend((j-1)*7+7)=DPoints(i1)%iEndP
					
						do k=1,3
							dpToSend((j-1)*iDPSCount+k)=DPoints(i1)%tvPL%v(k)
						enddo
						do k=1,3
							dpToSend((j-1)*iDPSCount+3+k)=DPoints(i1)%tvTL%v(k)
						enddo
						do k=1,3
							dpToSend((j-1)*iDPSCount+6+k)=DPoints(i1)%tvAcce%v(k)
						enddo
						do k=1,3
							dpToSend((j-1)*iDPSCount+9+k)=DPoints(i1)%tvPreV%v(k)
						enddo
						do k=1,3
							dpToSend((j-1)*iDPSCount+12+k)=DPoints(i1)%tvPreVT%v(k)
						enddo
						dpToSend((j-1)*iDPSCount+iDPSCount)=DPoints(i1)%dpMass
					enddo
					call ZQMPI_ISend(iNumToSend,1,iSrc,0)
					call ZQMPI_ISend(iIntToSend,iToSend*iIntSCount,iSrc,0)
					call ZQMPI_ISend(dpToSend,iToSend*iDPSCount,iSrc,0)
					nullify(iIntToSend,dpToSend,iNumToSend)

				endif
			enddo
			deallocate(iNumIndex,iIndex)

			!Receive
			do i=1,iNumBoxChange
				call MPI_Recv(iToSend,1,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,iStat,iIerr)
				iSrc=iStat(MPI_SOURCE)
				allocate(dpToSend(iToSend*iDPSCount),iIntToSend(iToSend*iIntSCount))
	
				call MPI_Recv(iIntToSend,iToSend*iIntSCount,MPI_INTEGER,iSrc,MPI_ANY_TAG,MPI_COMM_WORLD,iStat,iIerr)
	
				call MPI_Recv(dpToSend,iToSend*iDPSCount,MPI_DOUBLE_PRECISION,iSrc,MPI_ANY_TAG,MPI_COMM_WORLD,iStat,iIerr)
	
				do j=1,iToSend
					i1=iNTotalDPs+j
					DPoints(i1)%iPlane	=iIntToSend((j-1)*iIntSCount+1)
					DPoints(i1)%iID		=iIntToSend((j-1)*iIntSCount+2)
					DPoints(i1)%iLoopID	=iIntToSend((j-1)*iIntSCount+3)
					DPoints(i1)%iType	=iIntToSend((j-1)*iIntSCount+4)
					DPoints(i1)%iBurgers	=iIntToSend((j-1)*iIntSCount+5)
					DPoints(i1)%iBeginP	=iIntToSend((j-1)*iIntSCount+6)
					DPoints(i1)%iEndP	=iIntToSend((j-1)*iIntSCount+7)
					DPoints(i1)%iNumNei	=0
	
					do k=1,3
						DPoints(i1)%tvPL%v(k)		=dpToSend((j-1)*iDPSCount+k)
					enddo
					do k=1,3
						DPoints(i1)%tvTL%v(k)		=dpToSend((j-1)*iDPSCount+3+k)
					enddo
					do k=1,3
						DPoints(i1)%tvAcce%v(k)		=dpToSend((j-1)*iDPSCount+6+k)
					enddo
					do k=1,3
						DPoints(i1)%tvPreV%v(k)		=dpToSend((j-1)*iDPSCount+9+k)
					enddo
					do k=1,3
						DPoints(i1)%tvPreVT%v(k)	=dpToSend((j-1)*iDPSCount+12+k)
					enddo
					DPoints(i1)%dpMass=dpToSend((j-1)*iDPSCount+iDPSCount)

					call setGlobal(i1)
				enddo
	
				iNTotalDPs=iNTotalDPs+iToSend
				deallocate(dpToSend,iIntToSend)
	
			enddo
		endif

		deallocate(iTmp)
		iNVLPoints=iNTotalDPs
	endif

end subroutine CommConnPoints
!**************************************************************
!	CommGhostPoints:
!		Only communicate the ghost points.
!**************************************************************
subroutine CommGhostPoints
	use ZQMPI
	use variables,only:iNPoint,DPoints
	use MPIModule,only:iaGLIndex,iNumTransP,iaTransList,iNTotalDPS
	implicit none

	integer::i,j,k,i1,i2,i3,iToSend,iNumToRecv,iSrc
	integer::iNumBoxChange
	integer,dimension(:),pointer::iNumToSend
	integer,dimension(:),pointer::iIntToSend
	double precision,dimension(:),pointer::dpToSend
	integer,dimension(:),allocatable::iTmp
	double precision,dimension(:),allocatable::dpTmp
	integer,dimension(iNumProcs)::iNumIndex
	integer,dimension(iNumProcs,5000)::iIndex

	if(iMyid/=iMaster)then

		allocate(iTmp(5000))    !??????
		iNumBoxChange=0

		!Send list
		do i=1,iNumProcs-1

			if(i/=iMyid)then
				!....determine the points to request.....
				iToSend=0
				do j=1,iNumTransP
					if(iaGLIndex(iaTransList(j),1)==i)then
						iToSend=iToSend+1
						iTmp(iToSend)=iaTransList(j)
					endif
				enddo
				
				!....number....
				allocate(iNumToSend(1))
				iNumToSend(1)=iToSend
				call ZQMPI_ISend(iNumToSend,1,i,0)
				nullify(iNumToSend)

				!....points....
				if(iToSend>0)then
					iNumBoxChange=iNumBoxChange+1

					allocate(iIntToSend(iToSend))

					do j=1,iToSend
						iIntToSend(j)=iTmp(j)
					enddo

					call ZQMPI_ISend(iIntToSend,iToSend,i,0)

					nullify(iIntToSend)
				endif

			endif

		enddo
		! receive index
		do i=1,iNumProcs-1
			if(i/=iMyid)then
				call MPI_Recv(iToSend,1,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,iStat,iIerr)
				iSrc=iStat(MPI_SOURCE)
				iNumIndex(iSrc)=iToSend
				
				if(iToSend>0)then
					call MPI_Recv(iTmp,iToSend,MPI_INTEGER,iSrc,MPI_ANY_TAG,MPI_COMM_WORLD,iStat,iIerr)
					do j=1,iToSend
						iIndex(iSrc,j)=iTmp(j)
					enddo
				endif
			endif

		enddo
	endif

	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

	if(iMyid/=iMaster)then

		!Send To Others
		do i=1,iNumProcs-1
			if(i/=iMyid)then
				iToSend=iNumIndex(i)
				iSrc=i
				if(iToSend>0)then
					
					allocate(dpToSend(iToSend*12),iIntToSend(7*iToSend),iNumToSend(1))
					iNumToSend(1)=iToSend

					do j=1,iToSend

						i1=iaGLIndex(iIndex(i,j),2)
						iIntToSend((j-1)*7+1)=DPoints(i1)%iPlane
						iIntToSend((j-1)*7+2)=DPoints(i1)%iID
						iIntToSend((j-1)*7+3)=DPoints(i1)%iLoopID
						iIntToSend((j-1)*7+4)=DPoints(i1)%iType
						iIntToSend((j-1)*7+5)=DPoints(i1)%iBurgers
						if(DPoints(i1)%iBeginP/=-1)then
							iIntToSend((j-1)*7+6)=DPoints(DPoints(i1)%iBeginP)%iID
						else
							iIntToSend((j-1)*7+6)=DPoints(i1)%iBeginP
						endif
						if(DPoints(i1)%iEndP/=-1)then	
							iIntToSend((j-1)*7+7)=DPoints(DPoints(i1)%iEndP)%iID
						else
							iIntToSend((j-1)*7+7)=DPoints(i1)%iEndP
						endif
							
						!actually, ghost neighbors only need position information to calculate interaction,
						!no dynamics related information is needed, such as mass and acceleration.
						do k=1,3
							dpToSend((j-1)*12+k)=DPoints(i1)%tvPL%v(k)
						enddo
						
						do k=1,3
							dpToSend((j-1)*12+3+k)=DPoints(i1)%tvTL%v(k)
						enddo

						do k=1,3
							dpToSend((j-1)*12+6+k)=DPoints(i1)%tvPG%v(k)
						enddo

						do k=1,3
							dpToSend((j-1)*12+9+k)=DPoints(i1)%tvTG%v(k)
						enddo

					enddo
					call ZQMPI_ISend(iNumToSend,1,iSrc,0)
					call ZQMPI_ISend(iIntToSend,iToSend*7,iSrc,0)
					call ZQMPI_ISend(dpToSend,iToSend*12,iSrc,0)
					nullify(iIntToSend,dpToSend,iNumToSend)

				endif

			endif

		enddo

		!Receive

		do i=1,iNumBoxChange
			call MPI_Recv(iToSend,1,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,iStat,iIerr)
			iSrc=iStat(MPI_SOURCE)
			allocate(dpToSend(iToSend*12),iIntToSend(iToSend*7))

			call MPI_Recv(iIntToSend,iToSend*7,MPI_INTEGER,iSrc,MPI_ANY_TAG,MPI_COMM_WORLD,iStat,iIerr)

			call MPI_Recv(dpToSend,iToSend*12,MPI_DOUBLE_PRECISION,iSrc,MPI_ANY_TAG,MPI_COMM_WORLD,iStat,iIerr)

			do j=1,iToSend
				DPoints(iNTotalDPs+j)%iPlane	=iIntToSend((j-1)*7+1)
				DPoints(iNTotalDPs+j)%iID	=iIntToSend((j-1)*7+2)
				DPoints(iNTotalDPs+j)%iLoopID	=iIntToSend((j-1)*7+3)
				DPoints(iNTotalDPs+j)%iType	=iIntToSend((j-1)*7+4)
				DPoints(iNTotalDPs+j)%iBurgers	=iIntToSend((j-1)*7+5)
				DPoints(iNTotalDPs+j)%iBeginP	=iIntToSend((j-1)*7+6)
				DPoints(iNTotalDPs+j)%iEndP	=iIntToSend((j-1)*7+7)
				DPoints(iNTotalDPs+j)%iNumNei	=0

				do k=1,3
					DPoints(iNTotalDPs+j)%tvPL%v(k)		=dpToSend((j-1)*12+k)
				enddo
				do k=1,3
					DPoints(iNTotalDPs+j)%tvTL%v(k)		=dpToSend((j-1)*12+3+k)
				enddo
				do k=1,3
					DPoints(iNTotalDPs+j)%tvPG%v(k)		=dpToSend((j-1)*12+6+k)
				enddo
				do k=1,3
					DPoints(iNTotalDPs+j)%tvTG%v(k)		=dpToSend((j-1)*12+9+k)
				enddo

			enddo

			iNTotalDPs=iNTotalDPs+iToSend
			nullify(dpToSend,iIntToSend)

		enddo
		deallocate(iTmp)
	endif

	call MPI_Barrier(MPI_COMM_WORLD,iIerr)
end subroutine CommGhostPoints
!********************************************************
! GhostNeighbors:
!	New way to determine close boxes (tree nodes).
!
!??????????????????????????????????????????????????????????
! Updated:ZQ Wang,
!	 5/10/06
!	Check the corresponding variables in send-recive pairs
!********************************************************
subroutine GhostNeighbors
	use ZQMPI
	use vectors
	use MPIModule,only:tLocalRoot,iNumCloseBox,iIndexOfCloseBox,tpBoxRoots,iNumTNode, &
						iaGLIndex,iNTotalDPs,iNMaxS
	use variables,only:dcrit,iNPoint,iloop,DPoints,ConnectVec,ConnectPlane
	use LoopRearrangeMD,only:dpMaxAveLength
	use GhostTransMD,only:iNumOfNodeSend,tmpSendLog,tmpSendInt,tmpSendDps, &
							iNumOfNodeRece,tmpReceLog,tmpReceInt,tmpReceDps
	use FunctionMD,only:PointPBC
	implicit none

	integer :: status,iSrc

	integer :: i,j,k,ii,jj,kk,itemp,jtemp,ktemp
	integer::i1,i2,i3,i4,iNext
	type(vector) :: tempV1,tempV2,tempV3
	type(vector)::tmpVec1,tmpVec2
	double precision :: tmpS1,tmpS2,tmpS3,dpTmp1
	double precision :: GetDist

	integer,dimension(:),pointer::tmpIntP1,tmpIntP2
	double precision,dimension(:),pointer::tmpDPP

	integer::tmpNumClose
	integer,dimension(iNumProcs)::tmpIDClose
	type(vector)::tvP11,tvP12,tvT11,tvT12,tvP21,tvP22,tvT21,tvT22

	double precision::u,shap(4,3)

	logical::InMyNeighbor !,lSegBoxDis
	logical::lTmp,lVirtualLocalPoint

	!new Find neighbors
	integer,allocatable,dimension(:,:)::iNumTP
	integer,allocatable,dimension(:,:)::iTPId   ! iTPId(:)
	type(vector),allocatable,dimension(:,:)::tTPPos   !tTPPos(:)


	!////////////////////////////
	!Determine the close boxes.
	if(iMyid/=iMaster)then

		!....check distance between boxes and establish a list of close boxes....
		tmpS1=MAG(tLocalRoot%tSizeV)
		iNumCloseBox=0
		do i=1,iNumprocs-1
			if(i/=iMyid)then
				tmpS2=MAG(tpBoxRoots(i)%p%tSizev)
				tmpS3=GetDist(tLocalRoot%tCenterv,tpBoxRoots(i)%p%tCenterv)
				!///////////////////////////////////////////////
				! Determine the distance between two boxes.
				! Add to list, according to distance between boxes or distance between    
				! segments in local box and distant boxes.
				if(tmpS3<=tmpS1+tmpS2+dcrit+dpMAXAveLength)then  
					iNumCloseBox=iNumCloseBox+1
					iIndexOfCloseBox(iNumCloseBox)=i
				endif
			end if
		end do
	endif

	do i=1,iNumProcs-1

		tmpNumClose=iNumCloseBox

		call MPI_Bcast(tmpNumClose,1,MPI_INTEGER,i,MPI_COMM_WORLD,iIerr)
		if(tmpNumClose>0)then
			!tmp close box
			if(i==iMyid)then
				do j=1,iNumCloseBox
					tmpIDClose(j)=iIndexOfCloseBox(j)
				enddo
			endif
			call MPI_Bcast(tmpIDClose,tmpNumClose,MPI_INTEGER,i,MPI_COMM_WORLD,iIerr)

			if(i/=iMyid.and.iMyid/=iMaster)then
				do j=1,tmpNumClose   !throught list to see if it is in the list of close boxes.
					if(tmpIDClose(j)==iMyid)then ! if i am in, check if that box also in my list
						lTmp=.false.
						do k=1,iNumCloseBox
							if(iIndexOfCloseBox(k)==i)then
								lTmp=.true.
								exit
							endif
						enddo
						if(.not.lTmp)then
							iNumCloseBox=iNumCloseBox+1
							iIndexOfCloseBox(iNumCloseBox)=i
						endif
						exit	
					endif
				enddo
			endif

		endif
	enddo

	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

	!Transfer close box points
	if(iMyid/=iMaster)then
		do i=1,iNumCloseBox
			iSrc=iIndexOfCloseBox(i)
			allocate(tmpIntP1(2))
			tmpIntP1(1)=iNTotalDPs
			tmpIntP1(2)=iNPoint
			call ZQMPI_ISend(tmpIntP1,2,iSrc,0)
			nullify(tmpIntP1)
			if(iNtotalDPs>0)then			
				allocate(tmpIntP2(iNTotalDPs*5),tmpDPP(iNTotalDPs*6))
				
				do j=1,iNTotalDPs
					tmpIntP2((j-1)*5+1)=DPoints(j)%iID
					tmpIntP2((j-1)*5+2)=DPoints(j)%iType
					tmpIntP2((j-1)*5+3)=DPoints(j)%iPlane
					tmpIntP2((j-1)*5+4)=DPoints(j)%iBeginP
					tmpIntP2((j-1)*5+5)=DPoints(j)%iEndP

					tmpDPP((j-1)*6+1)=DPoints(j)%tvPG%v(1)
					tmpDPP((j-1)*6+2)=DPoints(j)%tvPG%v(2)
					tmpDPP((j-1)*6+3)=DPoints(j)%tvPG%v(3)
					tmpDPP((j-1)*6+4)=DPoints(j)%tvTG%v(1)
					tmpDPP((j-1)*6+5)=DPoints(j)%tvTG%v(2)
					tmpDPP((j-1)*6+6)=DPoints(j)%tvTG%v(3)
				enddo

				call ZQMPI_ISend(tmpIntP2,iNTotalDPs*5,iSrc,0)

				call ZQMPI_ISend(tmpDPP,iNTotalDPs*6,iSrc,0)
				nullify(tmpIntP2,tmpDPP)
			endif
		enddo
	endif

	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

	!--------receive points from neighbor tree nodes (processors).
	if(iMyid/=iMaster)then
		allocate(tmpIntP1(2))

		!Maximum number of DPs is iNMaxS for new split method. It is defined.
		allocate(iNumTP(iNumCloseBox,2),iTPID(iNumCloseBox,3*(iNmaxS+5)*5),tTPPos(iNumCloseBox,3*(iNmaxS+5)*2))


		do i=1,iNumCloseBox
			call MPI_Recv(tmpIntP1,2,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,iStat,iIerr)
			iSrc=iStat(MPI_SOURCE)
			iNumTP(i,1)=0
			iNumTP(i,2)=0

			if(tmpIntP1(1)>0)then
				allocate(tmpIntP2(tmpIntP1(1)*5),tmpDPP(tmpIntP1(1)*6))
				call MPI_Recv(tmpIntP2,tmpIntP1(1)*5,MPI_INTEGER,iSrc,MPI_ANY_TAG,MPI_COMM_WORLD,iStat,iIerr)

				call MPI_Recv(tmpDPP,tmpIntP1(1)*6,MPI_DOUBLE_PRECISION,iSrc,MPI_ANY_TAG,MPI_COMM_WORLD,iStat,iIerr)

				iNumTP(i,1)=tmpIntP1(1)   !total number of DPs received, including connection points, etc.
				iNumTP(i,2)=tmpIntP1(2)   !number of real segment DPs received.

				do j=1,iNumTP(i,1)
					iTPID(i,(j-1)*5+1)=tmpIntP2((j-1)*5+1)
					iTPID(i,(j-1)*5+2)=tmpIntP2((j-1)*5+2)
					iTPID(i,(j-1)*5+3)=tmpIntP2((j-1)*5+3)
					iTPID(i,(j-1)*5+4)=tmpIntP2((j-1)*5+4)
					iTPID(i,(j-1)*5+5)=tmpIntP2((j-1)*5+5)

					tTPPos(i,(j-1)*2+1)%v(1)=tmpDPP((j-1)*6+1)
					tTPPos(i,(j-1)*2+1)%v(2)=tmpDPP((j-1)*6+2)
					tTPPos(i,(j-1)*2+1)%v(3)=tmpDPP((j-1)*6+3)
					tTPPos(i,(j-1)*2+2)%v(1)=tmpDPP((j-1)*6+4)
					tTPPos(i,(j-1)*2+2)%v(2)=tmpDPP((j-1)*6+5)
					tTPPos(i,(j-1)*2+2)%v(3)=tmpDPP((j-1)*6+6)
				enddo
				deallocate(tmpIntP2,tmpDPP)
			endif
		enddo
		deallocate(tmpIntP1)
	endif

	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

	!===============================================================================
	!-----looking neighbors from neighbor tree nodes (processors)
	if(iMyid/=iMaster)then
		do i=1,iNumCloseBox
			!close non-local neighbors, only consider those segments which begins with ghost
			!local points.
			do ii=1,iNTotalDPs !here should be all points, including connection points
				if(DPoints(ii)%iEndP/=-1)then
					tvP11=DPoints(ii)%tvPG
					tvT11=DPoints(ii)%tvTG
					call TransNextGlobal(ii,tvP12,tvT12)
	
					do j=1,iNumTP(i,2)
						k=iTPId(i,(j-1)*5+1) !%iID
						if(.not.lVirtualLocalPoint(k))then
						if(.not. InMyNeighbor(ii,k).and.iTPID(i,(j-1)*5+5)/=-1)then
							!try to get the next global of the close points
							tvP21=tTPPos(i,(j-1)*2+1)
							tvT21=tTPPos(i,(j-1)*2+2)
	
							iNext=iTPID(i,(j-1)*5+5)
							tvP22=tTPPos(i,(iNext-1)*2+1)
							tvT22=tTPPos(i,(iNext-1)*2+2)
	
							do i1=1,21
								u=0.05*(i1-1)
								call getshape(shap,u)
								tmpVec1=shap(1,1)*tvP11+shap(2,1)*tvT11+shap(3,1)*tvP12+shap(4,1)*tvT12
								tmpVec1=PointPBC(tmpVec1)
								do i2=1,21
									u=0.05*(i2-1)
									call getshape(shap,u)
									tmpVec2=shap(1,1)*tvP21+shap(2,1)*tvT21+shap(3,1)*tvP22+shap(4,1)*tvT22
									tmpVec2=PointPBC(tmpVec2)
	
									if( GetDist(tmpVec1,tmpVec2)<dcrit+dpMaxAveLength)then
										call AddToMyNeighbor(ii,k)
										goto 1  !for next ghost segment
									endif
								enddo
							enddo
						endif
						endif
	
	1					continue
	
					enddo !j
				endif
			enddo !end ii	
		enddo
		deallocate(iTPID,tTPPos)
	endif

end subroutine GhostNeighbors
!***************************************************************************
!	lsegBoxDis:
!		Determine if the distance between segments in local box and other boxes
!		satisfies the criteria for close-box
!***************************************************************************
logical function lsegBoxDis(i)
	use ZQMPI
	use variables,only:iNPoint,DPoints,dCrit
	use MPIModule,only:tLocalRoot,iNumCloseBox,iIndexOfCloseBox,tpBoxRoots,iNTotalDPs
	use LoopRearrangeMD,only:dpMaxAveLength
	implicit none

	double precision::tmpS2,tmpS3
	integer::i,j,k,iEnd
	type(vector)::tvTmp1,tvTmp2
	type(vector)::tvP1,tvP2,tvT1,tvT2
	double precision::u,shap(4,3)

	double precision::getDist

	!////////////////////////////////////////
	lSegBoxDis=.false.
	tmpS2=MAG(tpBoxRoots(i)%p%tSizev)
	
	do i=1,iNTotalDPs

		if(DPoints(i)%iEndP/=-1)then !for each dislocation segments
			!....Ending points....
			tvP1=DPoints(i)%tvPG
			tvT1=DPoints(i)%tvTG
			call TransNextGlobal(i,tvP2,tvT2)

			!....Check distance....
			do j=1,11
				u=0.1*(j-1)
				call getshape(shap,u)
				tvTmp1=shap(1,1)*tvP1+shap(2,1)*tvT1+shap(3,1)*tvP2+shap(4,1)*tvT2
				tmpS3=GetDist(tvTmp1,tpBoxRoots(i)%p%tCenterv)
				if(tmpS3<tmpS2+dcrit+dpMAXAveLength)then
					lSegBoxDis=.true.
					exit
				endif
			enddo
		endif

		if(lSegBoxDis)exit
	enddo


end function lsegBoxDis
!*************************************************************************
!	CombineVelocityNumber:
!		Combine the counted number for different velocities.
!-------------------------------------------------------------------------
! Modified: ZQ Wang, 01/28/06
!	Change the way to count number and calculate length for different
!	velocity distributions.
!	02/26/06:
!	Increased dimension for velocity distribution.
!*************************************************************************
Subroutine CombineVelocityNumber
	use ZQMPI
	use Variables,only:dpVCountL,dpVCountN, iVCountN
	implicit none

	integer::i,j,k,iSum
	double precision::dpTmp,dpSum

	iSum=0
	dpSum=0.d0
	do i=1,100
		call MPI_AllReduce(iVCountN(i),k,1,MPI_INTEGER,MPI_SUM, &
									MPI_COMM_WORLD,iIerr)
		iVCountN(i)=k
		iSum=iSum+k
		call MPI_AllReduce(dpVCountL(i),dpTmp,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                                                                        MPI_COMM_WORLD,iIerr)
		dpVCountL(i)=dpTmp
		dpSum=dpSum+dpTmp
	enddo

	do i=1,100
		dpVCountL(i)=dpVCountL(i)/dpSum*100   !%
		dpVCountN(i)=(iVCountN(i)*1.d0)/(iSum*1.d0)*100  !%
	enddo
end subroutine CombineVelocityNumber
!=============================================================
END module CommunicationMD
