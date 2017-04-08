!**************************************************************************
!    Build_Global_Tree: 
!                    1. Build global part of the tree structure,
!                    2. Determine which dislocation goes to which process.
! 
!                      This operation is not needed to be done for every
!                      step. By rebuilding the tree, work load is balanced.
!**************************************************************************
Module SplitbuildMD
  use ZQMPI
  use TreeLeaf
	CONTAINS

!************************************************************************
!    Build_Global_Tree:
!************************************************************************
subroutine Build_Global_Tree
	use ZQMPI
	use TreeMD, only : NewBuildGlobalTree
	implicit none

	!MASTER:       
	if(iMyid==iMaster)then
		! Build tree, establish process-dislcoation list.
		call NewBuildGlobalTree
	endif
	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

	!distribute global tree, dislocation points, etc
	call DistributeGlobalTree
	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

end subroutine Build_Global_Tree
!***********************************************************************
!     DistributeGlobalTree: Distribute global tree, dislocation points,
!							glide planes, global index array to slaves,
!                           done by master.
!           what we have: tree structures, determine which dislocation goes
!                         where and maintain information.
!***********************************************************************
subroutine DistributeGlobalTree
	use ZQMPI
	use CommunicationMD,only : BcastGlobalTree,BCastGPlanes
	implicit none

	!////Distribute glide plane information.
	call BCastGPlanes
	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

	!////Distribute global index array for dislocation points
	call DistributeGIndexArray
	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

	    !////Distribute global tree; master brcast, slaves receives
	call BcastGlobalTree
	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

	!////Distribute subsegments; master sends, slaves receive
	call DistributeDPoints
	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

end subroutine DistributeGlobalTree
!***********************************************************************
!	DistributeGIndexArray:
!		1. Construct the index array;
!		2. Broadcast the index array to each processor.
!***********************************************************************
subroutine DistributeGIndexArray 
	use ZQMPI
	use vectors
	use CommunicationMD,only:BcastGIndexArray
	use TreeMD,only:CreateGIndexArray
	implicit none


	!/////////////////////////////////////////////
	!Create index array.
	if(iMyid==iMaster)then
		call CreateGIndexArray
	endif

	!/////////////////////////////////////////////
	!Broadcast index array.
	call BcastGIndexArray

end subroutine DistributeGIndexArray
!***********************************************************************
!     DistributeDPoints:   Pass dislocation points to slaves, "i" is the 
!                          ID of process.
!		                   what do we have:
!			              Groups of dislocation points from building of the
!						  global tree. Just send them
!***********************************************************************
subroutine DistributeDPoints
	use ZQMPI
	use CommunicationMD,only : DistributeDPs,ReceiveLocalDPs
	implicit none

	!////master does
	if(iMyid==iMaster)then
		call DistributeDPs       
	else  !////slaves do
		call ReceiveLocalDPs
	endif
	
	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

end subroutine DistributeDPoints
!***********************************************************************
!     Build_Local_Tree:
!                  1. After receive the data, slave processes build their
!                     own local tree and combine it with the global one.
!                  2. Calculate the multipole information
!                  3. Find neighbors and distribute the neighbors (dislocation
!                     position and multipole not in the same box) to others.
!         
!          if update, two things need to be done:
!                   1. calculate the multipole
!                   2. send multipole and dislocation position
!          this subroutine is executed by slaves and information is exchanged
!          between two slaves
!************************************************************************
subroutine Build_Local_Tree
	use ZQMPI
	use TreeMD,only : BuildLocalTree
	use CommunicationMD,only : TransGhostTree
	use variables,only:iloop
	implicit none

	!//////////////////////////////////////////////////////////////
	! Build local tree; slaves do it.
	call BuildLocalTree   
	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

!==========================================================================
	!///////////////////////////////////////////////////////////////
	! Get connection points
	call TransferConnectionPoints
	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

	!///////////////////////////////////////////////////////////////
	! Determine close Box
	call NewFindNeighbors2
	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

	!///////////////////////////////////////////////////////////////
	! Transfer ghost neighbors
	call TransGhostNeighbors
	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

!===========================================================================
	!///////////////////////////////////////////////////////////////
	!Build ghost tree here.
	!call TransGhostTree
	!call MPI_Barrier(MPI_Comm_WORlD,iIerr)

	!///////////////////////////////////////////////////////////////
	! Trees are now ready to be used.
	! Find tree neighbors if not update
	!call FindTreeNeighbors
	!call NewFindNeighbors
	!call MPI_Barrier(MPI_COMM_WORLD,iIerr)

	! Calculate the multipole information
	!if(iMyid/=iMaster)then
	! call CalculateMultipole
	!endif
	!call MPI_Barrier(MPI_COMM_WORLD,iIerr)
!===========================================================================
end subroutine Build_Local_Tree
!************************************************************************
!	TransferConnectionPoints
!************************************************************************
subroutine TransferConnectionPoints
	use ZQMPI
	use variables,only:iNPoint,DPoints
	use MPIModule,only:iaGLIndex,iNumTransP,iNTotalDPs
	use CommunicationMD,only:CommConnPoints
	implicit none

	integer::i,j,k,l,m,iID
	logical::lInTransList
	logical::lTmp


	if(iMyid/=iMaster)then
		iNumTransP=0

		do i=1,iNPoint
			!////////////////////////////////////////////////////
			! Check connection point of local points, at this point iBeginP and iEndP are still in global.
			do j=1,2
				if(j==1)then
					k=DPoints(i)%iBeginP
				else	
					k=DPoints(i)%iEndP
				endif
				
				call addToTransList(k)

			enddo

		enddo
	endif

	call MPI_Barrier(MPI_COMM_WORLD,iIErr)
	call CommConnPoints
	
	call MPI_Barrier(MPI_COMM_WORLD,iIErr)
	!transfer connection to local index
	! first 1 to iNpoint, then iNpoint+1 to iNTotalDPS
	if(iMyid/=iMaster)then
		do i=1,iNPoint
			if(DPoints(i)%iBeginP/=-1)then
				iID=DPoints(i)%iBeginP
				if(iaGLIndex(iID,1)==iMyid)then
					DPoints(i)%iBeginP=iaGLIndex(iID,2)   !here only change itself because the
									  !connected points will also be surfed.
				endif
			endif

			if(DPoints(i)%iEndP/=-1)then
				iID=DPoints(i)%iEndP
				if(iaGLIndex(iID,1)==iMyid)then
					DPoints(i)%iEndP=iaGLIndex(iID,2)
				endif
			endif
		enddo
		do i=iNPoint+1,iNTotalDPs
			if(DPoints(i)%iBeginP/=-1)then
				iID=DPoints(i)%iBeginP
				lTmp=.false.
				if(iaGLIndex(iID,1)==iMyid)then
					DPoints(i)%iBeginP=iaGLIndex(iID,2)
					DPoints(iaGLIndex(iID,2))%iEndP=i   !also change the connected local points
					lTmp=.true.
				else
					do j=iNPoint+1,iNTotalDPs
						if(DPoints(j)%iID==iID)then
							DPoints(i)%iBeginP=j
							lTmp=.true.
							exit
						endif
					enddo
				endif
				if(.not.lTmp)then
					DPoints(i)%iBeginP=-1
				endif
			endif
			if(DPoints(i)%iEndP/=-1)then
				iID=DPoints(i)%iEndP
				lTmp=.false.
				if(iaGLIndex(iID,1)==iMyid)then
					DPoints(i)%iEndP=iaGLIndex(iID,2)
					DPoints(iaGLIndex(iID,2))%iBeginP=i   !also change the connected local points
					lTmp=.true.
				else
					do j=iNPoint+1,iNTotalDPs
						if(DPoints(j)%iID==iID)then
							DPoints(i)%iEndP=j
							lTmp=.true.
							exit
						endif
					enddo
				endif
				if(.not.lTmp)then
					DPoints(i)%iEndP=-1
				endif
			endif
		enddo
	endif
	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

end subroutine TransferConnectionPoints
!************************************************************************
!	NewFindNeighbors2: 
!		Transfer all points in close boxes to current box and find the 
!	neighbors
!************************************************************************
subroutine NewFindNeighbors2
	use ZQMPI
	use vectors
	use variables,only:iNPoint,DPoints,dcrit
	use communicationMD,only:ghostneighbors
	use LoopRearrangeMD,only:dpMaxAveLength
	use MPIModule,only:iaGLIndex,iNTotalDPs
	use FunctionMD,only:PointPBC
	implicit none
	
	integer::i,j,k,ii
	integer::i1,i2,i3,i4
	type(vector)::tmpVec1,tmpVec2
	double precision::getDist
	logical::InMyNeighbor

	type(vector)::tvP11,tvP12,tvT11,tvT12,tvP21,tvP22,tvT21,tvT22
	double precision::u,shap(4,3)


	!.....find local neighbors....
	if(iMyid/=iMaster)then

		do i=1,iNTotalDPs
			DPoints(i)%iNumNei=0
		enddo

		do i=1,iNTotalDPs !here should be all points, including connection points
		if(DPoints(i)%iEndP/=-1)then
			tvP11=DPoints(i)%tvPG
			tvT11=DPoints(i)%tvTG
			call TransNextGlobal(i,tvP12,tvT12)

			do j=1,iNTotalDPs !here should be all points, including connection points
			if(DPoints(j)%iEndP/=-1)then
				tvP21=DPoints(j)%tvPG
				tvT21=DPoints(j)%tvTG
				call TransNextGlobal(j,tvP22,tvT22)
				if(i<j .and. .not. InMyNeighbor(i,DPoints(j)%iID))then
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

							if(GetDist(tmpVec1,tmpVec2)<dcrit+dpMaxAveLength)then
								call AddToMyNeighbor(i,DPoints(j)%iID)
								call AddToMyNeighbor(j,DPoints(i)%iID)
								goto 1   !next
							endif
						enddo
					enddo
				endif
1				continue
			endif
			enddo !j
		endif
		enddo !i
	endif

	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

	!....find ghost neighbors....
	call GhostNeighbors
	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

end subroutine NewFindNeighbors2
!************************************************************************
! TransGhostNeighbors:
!************************************************************************
subroutine TransGhostNeighbors
	use ZQMPI
	use Vectors
	use variables,only:DPoints
	use MPIModule,only:iNumTransP,iNTotalDPs,iaGLIndex,iNVLPoints
	use CommunicationMD,only:CommGhostPoints
	implicit none

	integer::i,j,k,l,m
	logical::lVirtualLocalPoint

	if(iMyid/=iMaster)then
		iNumTransP=0
		!//////////////////////////////////////////////////
		! Check neighbor points
		do i=1,iNVLPoints         !at this time, iNTotalDPs includes connection points
			do j=1,DPoints(i)%iNumNei
				k=DPoints(i)%iaIDOfNei(j)
				m=k
				if(.not. lVirtualLocalPoint(k))then
					call addToTransList(k)
					!//////////////////////////////////////
					! Check connection points of neighbors
					do l=2,2   !only next point
						!////////////////////////////////////
						! For neighbors' connection points, because neighbor points are not
						! on local now, its connection points are determined from iaGLIndex
						k=iaGlIndex(m,2+l)  
						if(.not.lVirtualLocalPoint(k))then 
							call addToTransList(k)
						endif
					enddo
				endif
			enddo
		enddo
	endif

	call MPI_Barrier(MPI_COMM_WORLD,iIerr)
	call CommGhostPoints

	call MPI_BARRIER(MPI_COMM_WORLD,iIerr)

	if(iMyid/=iMaster)then
		call ChangeToLocalIndex
	endif

	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

end subroutine TransGhostNeighbors
!*****************************************************************
!	ChangeToLocalIndex:
!		Here do two things:
!			1. change the neighbor list to local index
!			2. change the ghost neighbor connection to local index,
!				and before here, the local connection index is 
!				already local.
!
!		Attention:
!			Local points may apprear in the ghost neighbor connections,
!		so they also need to be reuqired.
!		
!		local and connected points don't need to change connection here.
!	Their relationship has been established!!!
!*****************************************************************
subroutine ChangeToLocalIndex
	use ZQMPI
	!use vectors
	use variables,only:iNPoint,DPoints
	use MPIModule,only:iNtotalDPs,iaGLIndex,iNVLPoints
	implicit none

	integer::i,j,k,iID
	logical::lStat
	logical::lVirtualLocalPoint

	!.1.
	do i=1,iNVLPoints     ! here should be iNPoint and the connection points
		do j=1,DPoints(i)%iNumNei
			iID=DPoints(i)%iaIDOfNei(j)
			if(iaGLIndex(iID,1)==iMyid)then
				DPoints(i)%iaIDOfNei(j)=iaGLIndex(iID,2)
			elseif(lVirtualLocalPoint(iID))then
				do k=iNPoint+1,iNVLPoints
					if(DPoints(k)%iID==iID)then
						DPoints(i)%iaIDOfNei(j)=k
						exit
					endif
				enddo
			else
				do k=iNVLPoints+1,iNTotalDPs
					if(DPoints(k)%iID==iID)then
						DPoints(i)%iaIDOfNei(j)=k
						exit
					endif
				enddo
			endif
		enddo
	enddo
		
	!.2. Ghost relations, endP should check iNVLPoints
	do i=iNVLPoints+1,iNTotalDPs
		do j=1,2
			if(j==1)then
				iID=DPoints(i)%iBeginP
			else
				iID=DPoints(i)%iEndP
			endif
			lStat=.false.
			do k=iNVLPoints+1,iNTotalDPs
				if(iID==DPoints(k)%iID)then
					if(j==1)then
						DPoints(i)%iBeginP=k
					else	
						DPoints(i)%iEndP=k
					endif
					lStat=.true.
					exit
				endif
			enddo
			!if not in ghost
			if(.not.lStat)then
				if(j==2)then
					do k=iNPoint+1,iNVLPoints
						if(iID==DPoints(i)%iID)then
							DPoints(i)%iEndP=k
							lStat=.true.
							exit
						endif
					enddo
				endif
			endif
			if(.not.lStat)then
				if(j==1)then
					DPoints(i)%iBeginP=-1
				else
					DPoints(i)%iEndP=-1
				endif
			endif
		enddo
	enddo

end subroutine ChangeToLocalIndex
!*****************************************************************
!	lInTransList:
!		Check if the index input is already in the transfer list.
!*****************************************************************
logical function lInTransList(iIndex)
	use ZQMPI
	use vectors
	use variables,only:
	use MPIModule,only:iaTransList,iNumTransP
	implicit none

	integer,intent(in)::iIndex
	integer::i

	lInTransList=.false.
	do i=1,iNumTransP
		if(iaTransList(i)==iIndex)then
			lInTransList=.true.
			exit
		endif

	enddo
	return

end function lInTransList
!*****************************************************************
!	AddToTransList:
!		Add a index to the transfer list of points that should be
!	transferred from other processors.
!		iIndex is the global index in the array iaGLIndex.
!*****************************************************************
subroutine AddToTransList(iIndex)
	use ZQMPI
	use MPIModule,only:iaTransList,iNumTransP,iaGLIndex
	implicit none

	integer,intent(in)::iIndex
!	logical::lInTransList

	if(iIndex/=-1)then !has connection points
		if(iaGLIndex(iIndex,1)/=iMyid)then !not on local
			if(.not. lInTransList(iIndex))then !not in the trans list
				!////////////////////////
				!add k to translist
				iNumTransP=iNumTransP+1
				iaTransList(iNumTransP)=iIndex
			endif
		endif
	endif

end subroutine AddToTransList
!==================================================================
end Module SplitbuildMD
