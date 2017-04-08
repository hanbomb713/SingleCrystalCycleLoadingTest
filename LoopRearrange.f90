!*****************************************************************
! LoopRearrange:
!	This subroutine is used to re-distribute the nodes on dislocation
!	loops. After nodes are re-distributed, tangents at points are 
!	calculated by using tangent_vec methods to asure the C2 continuity.
!*****************************************************************
subroutine LoopReArrange
	use ZQMPI
	use vectors
	use variables,only:iNPoint,DPoints
	implicit none

	integer::I
	external LocalRearrange
	
	call SurfLoop(LocalRearrange)

	do i=1,iNPoint
		call SetGlobal(i)
	enddo

	!
	call Elimination

	!
	call RegroupDP

end subroutine LoopReArrange
!****************************************************************
!	LocalRearrange:
!		Loop id has been determined. Now begin to rearrange the nodes
!	distribution for this particular loop.
!--------------------------------------------------------------------
! Modified: ZQ WANG, 01/26/06
!	Added tvTmpAcce,tvTmpPreV,dpTmpMass
!****************************************************************
subroutine LocalRearrange(PointID)
	use ZQMPI
	use vectors
	use variables,only:iNPoint, DPoints, MAx_Node,iDPTypeFree,&
						ConnectPlane,ConnectVec,zero
	use LoopRearrangeMD,only:itmpNumP,iTmpNewNumP,tvTmpPL,tvTmpTL,iTmpIndex, lTmpToCal,&
				iTmpLooptype,itmpPlane,itmpBurgers,itmpType1,itmpType2,tvTmpAcce,tvTmpPreV,dpTmpMass,tvTmpPreVT
	use FunctionMD,only:DeterConnVec
	implicit none

	integer,intent(in)::PointID
	integer::iEnterID,iCur,iBegin,iNext
	integer::i,j,k,looptype,state
	integer::id1,id2
	logical::lContinue
	logical::lSegCriteria
	logical::StopSurfLoop

	lSegCriteria	=	.false.
	lContinue		=	.true.
	iEnterID		=	PointID
	iCur			=	PointID
	iBegin			=	DPoints(PointID)%iBeginP
	lTmpToCal=.false.

	if(DPoints(iEnterID)%iBeginP/=-1)then
		looptype=0
	else
		looptype=1
	endif

	itmpLooptype=looptype
	itmpBurgers=DPoints(PointID)%iBurgers
	itmpPlane=DPoints(PointID)%iPlane

	iTmpNumP		=	0
	
	do while(lCOntinue)
		iNext=DPoints(iCur)%iEndP

		lContinue=.not.StopSurfLoop(iCur,iEnterID,lSegCriteria)
		!....remember indexes and move on....
		iTmpNumP=iTmpNumP+1
		if(DPoints(iCur)%lToCal)lTmpToCal=.true.
		iTmpIndex(iTmpNumP)=iCur

		iCur=iNext
	end do

	!....get tmpTL,tmpPL, and arrange
	if(.not. lContinue .and. looptype==1 .and. iTmpNumP>1)then   !here specify the close loop condition.
		
		iTmpNewNumP=iTmpNumP
		itmpType1=DPoints(iTmpIndex(1))%iType
		itmpType2=DPoints(iTmpIndex(iTmpNump))%iType
		if(itmpLoopType==1)then
			call copyEndPoint(iTmpIndex(1),iTmpIndex(iTmpNump))
		endif
		! assign values
		ConnectPlane=DPoints(PointID)%iPlane
		ConnectVec=zero%v(1)
		do i=1,iTmpNumP
			id2=itmpIndex(i)
			if(i>1)then
				j=i-1
				id1=itmpIndex(j)
				ConnectVec=ConnectVec+DeterConnVec(id2,id1) !i in front of j
			endif
			k=1

			call GetRALocal(id2,tvTmpPL(i),tvTmpTL(i),k)
			tvTmpAcce(i)=DPoints(id2)%tvAcce
			tvTmpPreV(i)=DPoints(id2)%tvPreV
			tvTmpPreVT(i)=DPoints(id2)%tvPreVT
			dpTmpMass(i)= DPoints(id2)%dpMass
		enddo

		!first need to calculate the TL
		call getTangent(itmpNewNumP,tvTmpPL,tvTmpTL,itmpLooptype)   !for itmpNewNum>3

		call LoopReArrangementT
	
		call DPReposition

	endif

end subroutine LocalRearrange
!*******************************************************************
!	LoopRearrangementT
!		Loop TL, PL have been determined. move nodes to different
!	positions and compute the tangents.
!*******************************************************************
subroutine LoopReArrangementT
	use ZQMPI
	implicit none


	!....determine the number of node and  average length....
	call checkAveLengthNodes

	!....move nodes to their positions....
	call moveNodes

	!....compute the tangent vectors from position vectors....
	call computeTangent

end subroutine LoopReArrangementT
!******************************************************************
!    checkAveLength:
!          determine the number of nodes and average length for rearrangement
!
!	Modified: ZQ Wang
!	06/06/06:
!	While calculate the intermediate points, normalize the tangent vectors
!	with the segment length.
!******************************************************************
subroutine checkAveLengthNodes
	use ZQMPI
	use vectors
	use variables,only:NPoint_I,Max_Node
	use LoopReArrangeMD,only:dpAveLength,iNoOfNodes,dpMaxAveLength,                       &
							itmpIndex,itmpNumP,itmpNewNumP,tvTmpPL,tvTmpTL,&
							itmpLooptype
	use ANNIHILATIONMD,only:iAnniCount
	use CrossslipMD,only:iCSCount
	implicit none

	integer::i,j,k
	double precision::shap(4,3),u0
	type(vector)::p1,p2
	double precision::dpsegL

	!///////////////////////////////////////////////

	!Calculate the total length
	p1=tvTmpPL(1)
	dpAveLength=0.d0
	do i=1,iTmpNumP-1
		dpSegL=MAG(tvTmpPL(i)-tvTmpPL(i+1))
		do j=1,20
			u0=j*0.05
			call getshape(shap,u0)
!--------------------------------------------------------------------------------------------
		       if(iAnniCount>0 .or. iCSCount>0)then
			p2=shap(1,1)*tvTmpPL(i)+shap(2,1)*dpSegL*UV(tvTmpTL(i))+ &
			   shap(3,1)*tvTmpPL(i+1)+shap(4,1)*dpSegL*UV(tvTmpTL(i+1))
			else
                       p2=shap(1,1)*tvTmpPL(i)+shap(2,1)*tvTmpTL(i)+ &
                           shap(3,1)*tvTmpPL(i+1)+shap(4,1)*tvTmpTL(i+1)
			endif
!-------------------------------------------------------------------------------------------
			dpAveLength=dpAveLength+MAG(p2-p1)
			p1=p2
		enddo
	enddo
	if(itmpLooptype==0)then
		dpSegL=MAG(tvTmpPL(iTmpNumP)-tvTmpPL(1))
		do j=1,20
			u0=j*0.05
			call getshape(shap,u0)

!---------------------------------------------------------------------------------------------------
			if(iAnniCount>0 .or. iCSCount>0)then
			p2=shap(1,1)*tvTmpPL(iTmpNumP)+shap(2,1)*dpSegL*UV(tvTmpTL(iTmpNumP))+ &
				shap(3,1)*tvTmpPL(1)+shap(4,1)*dpSegL*UV(tvTmpTL(1))
			else
			p2=shap(1,1)*tvTmpPL(iTmpNumP)+shap(2,1)*tvTmpTL(iTmpNumP)+ &
                                shap(3,1)*tvTmpPL(1)+shap(4,1)*tvTmpTL(1)
			endif
!---------------------------------------------------------------------------------------------------

			dpAveLength=dpAveLength+MAG(p2-p1)
			p1=p2
		enddo
	endif

	!Determine the number of nodes
	iNoOfNodes=iTmpNumP

	if(iTmpNumP<MAX_Node)then
		i=int(dpAvelength/dpMaxAveLength)

		!Add points if nodes needed is more than real nodes
		!05/27/06
		if(i>iTmpNumP)then
			iNoOfNodes=iTmpNumP+1  !    abs(i-iTmpNumP)
		!Remove points if nodes needed is less than real nodes
		elseif(i<iTmpNumP)then
			iNoOfNodes=iTMpNumP-1  !    abs(i-iTMpNumP)
		endif
	endif
	if(iNoOfNodes>Max_Node)iNoOfNodes=Max_Node
	if(iNoOfNodes<NPoint_I)iNoOfNodes=NPoint_I

	!Calculate new
	dpAveLength=dpAveLength/(iNoOfNodes-itmplooptype)

end subroutine checkAveLengthNodes
!*******************************************************************
!	moveNodes:
!			move nodes to their positions. Compute the position vector
!	Here, we generate temporary nodes tvaTmpP. We arrange them to fit the 
!	average length and find a group of selected points from tvaTmpP. These 
!	points are assigned to the final list of new points.
!	Alone with them are the mass, acceleration, and previous velocity.
!---------------------------------------------------------------------------
! Modified: ZQ Wang, 01/26/06
!    Added the calculation of mass, acceleration and velocity
!
!	06/06/06
!	Normalized the boundary tangent vectors with the segment length.
!*******************************************************************
subroutine moveNodes
	use ZQMPI
	use vectors
	use variables,only:max_node
	use LoopReArrangeMD,only:dpAveLength,iNoOfNodes,&
				itmpNumP,itmpNewNumP,itmpLooptype,tvTmpPL,tvTmpTL,tvTmpAcce,tvTmpPreV,dpTmpMass,tvTmpPreVT
	use ANNIHILATIONMD,only:iAnniCount
	use CrossslipMD,only:iCSCount
	implicit none

	type(vector),dimension(1000)::tvaTmpP,tvaTmpAcce,tvaTmpPreV,tvaTmpPreVT     !allocate 1000, should be ok if # of seg less than 200
	double precision,dimension(1000)::dpaTmpMass
	type(vector),dimension(MAX_Node)::tvaTmpPL,tvaTmpAcceL,tvaTmpPreVL,tvaTmpPreVTL
	double precision,dimension(MAX_Node)::dpaTmpMassL
	type(vector)::tvTmpP
	integer::iNumPSeg
	integer::i,j,k,i1,i2
	double precision::u0,shap(4,3)
	double precision::dpTmpSegLength,dpAveLength1,dpAveLength2
	double precision::dpSegL					! 06/06/06, to normalize the tangent vector

	!//////////////////////////////////////////////////////////////////
	!....generate tmp points.....
	iNumPSeg=500/(itmpNumP-itmpLooptype)	
	if(iNumPSeg<5)iNumPSeg=5

	i1=0  ! total number of subpoints
	do i=1,itmpNumP-1
		dpSegL=MAG(tvTmpPL(i)-tvTmpPL(i+1))
		do j=1,iNumPSeg
			u0=(j-1)*1.d0/(iNumPSeg-1)
			call getshape(shap,u0)

!---------------------------------------------------------------------------------------------------
			if(iAnniCount>0 .or. iCSCount>0)then
			tvTmpP=shap(1,1)*tvTmpPL(I)+shap(2,1)*dpSegL*UV(tvTmpTL(I))+&
					shap(3,1)*tvTmpPL(I+1)+shap(4,1)*dpSegL*UV(tvTmpTL(I+1))
			else
			tvTmpP=shap(1,1)*tvTmpPL(I)+shap(2,1)*tvTmpTL(I)+&
                                        shap(3,1)*tvTmpPL(I+1)+shap(4,1)*tvTmpTL(I+1)
			endif
!---------------------------------------------------------------------------------------------------
			i1=i1+1
			tvaTmpP(i1)=tvTmpP
			!
			tvaTmpAcce(i1)=(1-u0)*tvTmpAcce(I)+u0*tvTmpAcce(I+1)
			tvaTmpPreV(i1)=(1-u0)*tvTmpPreV(I)+u0*tvTmpPreV(I+1)
			tvaTmpPreVT(i1)=(1-u0)*tvTmpPreVT(I)+u0*tvTmpPreVT(I+1)
			dpaTmpMass(i1)=(1-u0)*dpTmpMass(I)+u0*dpTmpMass(I+1)
			!
		enddo
	enddo

	if(itmpLooptype==0)then
		dpSegL=MAG(tvTmpPL(itmpNumP)-tvTmpPL(1))
		do j=1,iNumPSeg
			u0=(j-1)*1.d0/(iNumPSeg-1)
			call getshape(shap,u0)

!----------------------------------------------------------------------------------------------------------
			if(iAnniCount>0 .or. iCSCount>0)then
			tvTmpP=shap(1,1)*tvTmpPL(itmpNumP)+shap(2,1)*dpSegL*UV(tvTmpTL(itmpNumP))+&
					shap(3,1)*tvTmpPL(1)+shap(4,1)*dpSegL*UV(tvTmpTL(1))
			else
			tvTmpP=shap(1,1)*tvTmpPL(itmpNumP)+shap(2,1)*tvTmpTL(itmpNumP)+&
                                        shap(3,1)*tvTmpPL(1)+shap(4,1)*tvTmpTL(1)
			endif
!----------------------------------------------------------------------------------------------------------
			i1=i1+1
			tvaTmpP(i1)=tvTmpP
			!
			tvaTmpAcce(i1)=(1-u0)*tvTmpAcce(itmpNumP)+u0*tvTmpAcce(1)
			tvaTmpPreV(i1)=(1-u0)*tvTmpPreV(itmpNumP)+u0*tvTmpPreV(1)
			tvaTmpPreVT(i1)=(1-u0)*tvTmpPreVT(itmpNumP)+u0*tvTmpPreVT(1)
			dpaTmpMass(i1)=(1-u0)*dpTmpMass(itmpNumP)+u0*dpTmpMass(1)
			!
		enddo
	endif
			
	!///////////////////////////////////////////////////////////////////
	!....pick and position points....
	!///////
	!i2 is the number of points that have been applied to the line
	i2=1
	tvaTmpPL(1)=tvaTmpP(1)
	!
	tvaTmpAcceL(1)=tvaTmpAcce(1)
	tvaTmpPreVL(1)=tvaTmpPreV(1)
	tvaTmpPreVTL(1)=tvaTmpPreVT(1)
	dpaTmpMassL(1)=dpaTmpMass(1)
	!
	dpTmpSegLength=0.d0
	do i=2,i1-1
		!/////////////////////////////////////
		! count for the last second point
		dpTmpSegLength=dpTmpSegLength+MAG(tvaTmpP(i)-tvaTmpP(i-1))

		!///////////////////////////
		! if the dpTmpSegLength larger than avelength, take the mid
		! points.
		!////////////////////
		if(dpTmpSegLength>dpAveLength)then
			i2=i2+1
			tvaTmpPL(i2)=0.5d0*(tvaTmpP(i-1)+tvaTmpP(i))
			!
			tvaTmpAcceL(i2)=0.5d0*(tvaTmpAcce(i-1)+tvaTmpAcce(i))
			tvaTmpPreVL(i2)=0.5d0*(tvaTmpPreV(i-1)+tvaTmpPreV(i))
			tvaTmpPreVTL(i2)=0.5d0*(tvaTmpPreVT(i-1)+tvaTmpPreVT(i))
			dpaTmpMassL(i2)=0.5d0*(dpaTmpMass(i-1)+dpaTmpMass(i))
			!
			!....record the new beginning seg length....
			if(itmplooptype==1 .and. i2==2)then
				dpAveLength1=dpTmpSegLength
			endif
			! reset the dpTmpSegLength
			dpTmpSegLength=MAG(tvaTmpP(i)-tvaTmpPL(i2)) 

			if((itmpLooptype==1.and.i2>=iNoOfNodes-1 .and.i2<Max_node) .or. &
				(itmpLooptype==0 .and. i2>=iNoOfNodes.and.i2<=Max_node))then
				exit
			endif
		endif
	enddo
	!.... last point....
	if(itmpLooptype==1)then
		i2=i2+1
		tvaTmpPL(i2)=tvaTmpP(i1)
		!
		tvaTmpAcceL(i2)=tvaTmpAcce(i1)
		tvaTmpPreVL(i2)=tvaTmpPreV(i1)
		tvaTmpPreVTL(i2)=tvaTmpPreVT(i1)
		dpaTmpMassL(i2)=dpaTmpMass(i1)
		!
	endif
	if(itmplooptype==1)then
		dpAveLength2=dpTmpSegLength
	endif

	!/////////////////////////////////////////////////////////////////////////////////
	!.....all nodes have been reArranged......
	itmpNewNumP=i2
	do i=1,i2
		tvTmpPL(I)=tvaTmpPL(i)
		!
		tvTmpAcce(I)=tvaTmpAcceL(i)
		tvTmpPreV(I)=tvaTmpPreVL(i)
		tvTmpPreVT(I)=tvaTmpPreVTL(i)
		dpTmpMass(I)=dpaTmpMassL(i)
		!
	enddo

	!Assign boundary conditions for tangent vectors
	!Tangent vectors are normalized with the segment length.       06/06/06.
	if(itmpLooptype==1)then
		!first point
		dpSegL=MAG(tvTmpPL(1)-tvTmpPL(2))

!-------------------------------------------------------
		if(iAnniCount>0 .or. iCSCount>0)then
		tvTmpTL(1)=dpSegL*UV(tvTmpTL(1))
		else
		tvTmpTL(1)=tvTmpTL(1)
		endif
!-------------------------------------------------------

		!last point
		dpSegL=MAG(tvTmpPL(itmpNewNumP-1)-tvTmpPL(itmpNewNumP))

!------------------------------------------------------------------
		if(iAnniCount>0 .or. iCSCount>0)then
		tvTmpTL(itmpNewNumP)=dpSegL*UV(tvTmpTL(itmpNumP))
		else
		tvTmpTL(itmpNewNumP)=tvTmpTL(itmpNumP)
		endif

	endif

end subroutine moveNodes
!*******************************************************************
!	computeTangent:
!			compute the tangent vector from the position vector
!          9/9/03: currently for only open loops
!*******************************************************************
subroutine ComputeTangent
	use ZQMPI
	use vectors
	use LoopRearrangeMD,only:itmpNewNumP,tvtmpPL,tvTmpTL,itmpLooptype
	
	call getTangent(itmpNewNumP,tvTmpPL,tvTmpTL,itmpLooptype)

end subroutine ComputeTangent
!*****************************************************************
!	DPRepostion:
!		To put dislocation points into array, and finalize loop
!	annihilation
!		Pay attention to the end points of open loop
!
!------------------------------------------------------------------
! Modified: ZQ Wang, 01/26/06
!	Added the tvTmpAcce,tvTmpPreV,dpTmpMass.
!*****************************************************************
subroutine DPReposition
	use ZQMPI
	use vectors
	use variables,only:iNPoint,DPoints,iDPTypeFree
	use LoopRearrangeMD,only:itmpNumP,itmpNewNumP,tvTmpPL,tvTmpTL,itmpIndex,&
					itmpLooptype,itmpPlane,itmpBurgers,itmpType1,itmpType2, &
					tmpDP,tvTmpAcce,tvTmpPreV,dpTmpMass,tvTmpPreVT,lTmpToCal
	implicit none
	integer::i,j,k
	integer::iBegin,iEnd
	integer::ii,jj,kk


	!.....addd more indexes.....	
	do i=1,itmpNewNumP-itmpNumP
		itmpIndex(itmpNump+i)=iNPoint+i
	enddo
	!....assign new values to points....
	do i=1,itmpNewNumP
		j=itmpIndex(i)
		DPoints(j)%iID=j
		DPoints(j)%iBurgers=itmpBurgers
		DPoints(j)%lStat=.true.
		DPoints(j)%lToCal=lTmpToCal
		DPoints(j)%iNumNei=0
		DPoints(j)%iType=iDPtypeFree

		DPoints(j)%tvPL=tvTmpPL(i)
		DPoints(j)%tvTL=tvTmpTL(i)
		DPoints(j)%iPlane=itmpPlane
		
		!
		DPoints(j)%tvAcce=tvTmpAcce(i)
		DPoints(j)%tvPreV=tvTmpPreV(i)
		DPoints(j)%tvPreVT=tvTmpPreVT(i)
		DPoints(j)%dpMass=dpTmpMass(i)
		!
		!points connection
		if(i>1)then
			DPoints(j)%iBeginP=itmpIndex(i-1)
		elseif(itmpLooptype/=1)then
			DPoints(j)%iBeginP=itmpIndex(itmpNewNumP)
		endif

		if(i<itmpNewNumP)then
			DPoints(j)%iEndP=itmpIndex(i+1)
		elseif(itmpLooptype/=1)then
			DPoints(j)%iEndP=itmpIndex(1)
		endif

		if(itmpLooptype==1 .and.(i==1 .or. i==itmpNewNumP))then
			if(i==1)then
				ii=1
				DPoints(j)%iBeginP=tmpDP(ii)%iBeginP
				if(DPoints(j)%iBeginP/=-1)DPoints(DPoints(j)%iBeginP)%iEndP=j
			else
				ii=2
				DPoints(j)%iEndP=tmpDP(ii)%iEndP
				if(Dpoints(j)%iEndP/=-1)DPoints(DPoints(j)%iEndP)%iBeginP=j
			endif
			DPoints(j)%iType=tmpDP(ii)%iType

		endif
		
	enddo

	if(itmpNewNumP>itmpNumP)then
		iNPoint=iNPoint+itmpNewNumP-itmpNumP
	else
		do i=1,itmpNumP-itmpNewNumP
			DPoints(itmpIndex(itmpNewNumP+i))%iBeginP=-1
			DPoints(itmpIndex(itmpNewNumP+i))%iEndP=-1
			DPoints(itmpIndex(itmpNewNumP+i))%iType=iDPTypeFree
			DPoints(itmpIndex(itmpNewNumP+i))%lStat=.true.
		enddo
	endif


end subroutine
!***************************************************************
! copyEndPoint
!***************************************************************
subroutine copyEndPoint(i1,i2)
	use ZQMPI
	use vectors
	use variables,only:DPoints,MAX_Neighbor
	use LoopRearrangeMD,only:tmpDP
	implicit none

	integer,intent(in)::i1,i2

	integer::ii(2)
	integer::i,j,k

	ii(1)=i1
	ii(2)=i2

	do i=1,2
		j=ii(i)
		tmpDP(i)=DPoints(j)
	enddo

end subroutine
