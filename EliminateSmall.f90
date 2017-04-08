!*********************************************************
!	EliminateSmall:
!		If a close loop is very small. Just eliminate it. Because
!	It is difficult to perform dynamics.
!		if a point is pointing to itself, also eliminate it.
!		
!		lToEliminate = .true.: to eliminate
!---------------------------------------------------------------
! Modified: ZQ Wang, 01/26/06
!	Added checksingle to eliminate isolated single point, this is 
!	especially for the rearrange case with less point than previous 
!	configuration.
!	
!	5//06: Change conditions for open loops
!*********************************************************
subroutine Elimination
	use ZQMPI
	use vectors
	use variables,only:iNPoint,DPoints,iDPTypeFixed,MAX_NODE,MAX_PLANE
	implicit none

	integer::i,j,k
	logical,allocatable,dimension(:)::lToEliminate

	allocate(lToEliminate(MAX_NODE*MAX_PLANE))
	do i=1,iNPoint
		DPoints(i)%lStat=.false.
		lToEliminate(i)=.false.
	enddo

	
	!open loop, beginning point must have a next point
	do i=1,iNPoint
		if(.not. DPoints(i)%lStat .and. DPoints(i)%iBeginP==-1 .and. &
									DPoints(i)%iEndP/=-1)then!beginning point of a open loop
			call setOpenToTrue(i,lToEliminate)
		endif
	enddo

	!check close
	do i=1,iNPoint
		if(.not.DPoints(i)%lStat .and. DPoints(i)%iBeginP/=-1)then
			call checkClose(i,lToEliminate)
		endif
	enddo

	!check single and all other
	do i=1,iNPoint
		if(.not.DPoints(i)%lStat)then
			call checkSingle(i,lToEliminate)
		endif
	enddo
	call ToEliminate(lToEliminate)

	deallocate(lToEliminate)
end subroutine Elimination
!***************************************************************
!checkClose:
!	Check if the close loop should be eliminated. Set those to 
!	true in lToE.
!	
!	1. If the close loop is very small;
!	2. if the close loop has only one point.
!	
!	so, length should be calculated.
!
!  Modified: ZQ Wang, 
!	5/17/06: changed the length to delete
!***************************************************************
subroutine checkClose(PointID,lToE)
	use ZQMPI
	use vectors
	use variables,only:DPoints,iNPoint,iDPTypeFixed
	use FunctionMD,only:DeterConnVecQ
	implicit none

	integer,intent(in)::POintID
	logical,dimension(iNPoint)::lToE

	integer::i,j,k
	integer::iNumIndex
	integer,dimension(200)::tmpIndex   !if max_node <200, then 200 is enough
	type(vector)::PG1,PG2,tmpVec

	integer::iENterID,iCur,iNext
	logical::lContinue

	double precision::tmpLength

	iEnterID=PointID
	iCUr=PointID

	lContinue=.true.
	iNumIndex=0
	
	tmpLength=0.d0

	do while(lContinue)
		iNumIndex=iNumIndex+1
		tmpIndex(iNumIndex)=iCur
		
		DPoints(iCur)%lStat=.true.
		DPoints(iCur)%iType=iDPTypeFixed              !set close to fixed.
		iNext=DPoints(icur)%iEndP
		iCur=iNext
		
		if(iNumIndex>1)then
			PG1=DPoints(tmpIndex(iNumIndex-1))%tvPG
			PG2=DPoints(tmpIndex(iNumIndex))%tvPG
			tmpVec=DeterConnVecQ(PG2,PG1)
			PG2=PG2-tmpVec

			tmpLength=tmpLength+MAG(PG2-PG1)
		endif
		
		if(iNext==iEnterID)then
			lContinue=.false.
		endif

	enddo

	if(iNumIndex<4 .or. tmpLength<60.d0)then ! if num<4, then eliminate this loop  !5/17/06: <60.d0; r=10a
		do i=1,iNumIndex
			lToE(tmpIndex(i))=.true.
		enddo
	endif

end subroutine checkClose
!***************************************************************
!	setiOpenToTrue:
!		set all connected point to true, for open loop
!***************************************************************
subroutine setOpenToTrue(PointID,lToE)
	use ZQMPI
	use vectors
	use variables,only:DPoints,iDPTypeFree,iDPTypeFixed,iNPoint
	implicit none

	integer,intent(in)::PointID
	integer::I,j,k
	logical,dimension(iNPoint)::lToE

	logical::lContinue

	integer::iEnterID,iCur,iNext

	lContinue=.true.
	iCur=PointID
	iEnterID=PointID

	do while(lContinue)
		DPoints(iCur)%lStat=.true.

		!--------------------------------------------------------------
		! Set type
		if(DPoints(iCur)%iBeginP==-1 .or. DPoints(iCur)%iEndP==-1)then
			DPoints(iCur)%iType=iDPTypeFixed
		else
			DPoints(iCur)%iType=iDPTypeFree
		endif


		iNext=DPoints(iCur)%iEndP
		if(iNext==-1)then
			lContinue=.false.
		endif
		iCur=iNext
	enddo

end subroutine setOpenToTrue
!****************************************************************
!	checkSingle:
!		check if a point is a single isolated point. If yes,
!	eliminate this point.
!---------------------------------------------------------------
! Modified: ZQ Wang, 01/26/06
!	New
!	ZQ Wang, 09/27/06
!	Comment out the connection condition. Remove all non-open non-close loops.
!****************************************************************
subroutine checkSingle(PointID,lToE)
	use ZQMPI
	use vectors
	use variables,only:DPoints,iNPoint
	implicit none

	integer,intent(in)::PointID
	logical,dimension(iNPoint)::lToE

	!/////////////////////////////////
	! 
	DPoints(PointID)%lStat=.true.

!	if(DPoints(PointID)%iBeginP==-1.and.DPoints(PointID)%iEndP==-1)then
		lToE(PointID)=.true.
!	endif

end subroutine checkSingle
!***************************************************************
!	ToEliminate:
!	
!	With the flags have been set, now to perform the action to 
!	delete dislocation points.	
!	
!	The procedure is as following:
!	1. Do from end to beginning
!	2. while move a point forward, only those points connected to 
!		it was reset.
!---------------------------------------------------------------
! Modified: ZQ Wang, 01/26/06
!	Added the tvAcce,tvPreV,dpMass
!***************************************************************
subroutine toEliminate(lToE)
	use ZQMPI
	use vectors
	use variables,only:iNPoint,DPoints
	implicit none

	logical,intent(in),dimension(iNPoint)::lToE

	integer::i,j,k
	integer::ibegin,iEnd,iCur

	j=iNPoint

	do i=j,1,-1
		if(lToE(i))then
			do k=i,iNPoint-1
!				DPoints(k)=DPoints(k+1)
				!new assign
				DPoints(k)%iloopID	=	DPoints(k+1)%iloopID
				DPoints(k)%iType	=	DPoints(k+1)%iType
				DPoints(k)%tvPL		=	DPoints(k+1)%tvPL
				DPoints(k)%tvTL		=	DPoints(k+1)%tvTL
				DPoints(k)%tvPG		=	DPoints(k+1)%tvPG
				DPoints(k)%tvTG		=	DPoints(k+1)%tvTG
				DPoints(k)%iPlane	=	DPoints(k+1)%iPlane
				DPoints(k)%iBurgers	=	DPoints(k+1)%iBurgers
				DPoints(k)%lToCal	=	DPoints(k+1)%lToCal
				!inertial
				DPoints(k)%tvAcce	=	DPoints(k+1)%tvAcce
				DPoints(k)%tvPreV	=	DPoints(k+1)%tvPreV
				DPoints(k)%dpMass	=	DPoints(k+1)%dpMass
				DPoints(k)%tvPreVT	=	DPoints(k+1)%tvPreVT
				!

				DPoints(k)%iBeginP	=	DPoints(k+1)%iBeginP
				DPoints(k)%iEndP	=	DPoints(k+1)%iEndP

				DPoints(k)%iID=k

				!reset connection
				iBegin=DPoints(k)%iBeginP
				iEnd=DPoints(k)%iEndP

				if(iBegin/=-1)DPoints(iBegin)%iEndP=k
				if(iEnd/=-1)DPoints(iEnd)%iBeginP=k
			enddo
			iNPoint=iNPoint-1
		endif
	enddo

end subroutine toEliminate

