!**************************************************************
!	NewAnni:
!		Annihilation implementation.
!**************************************************************
subroutine NewAnni
	use ZQMPI
	use vectors
	use variables,only:iNPoint,DPoints
	use ANNIHILATIONMD,only:iAnniCount,dpAnniLength
	implicit none

	integer::i,j,k
	logical::IsCloseLoopPoint

	iAnniCount=0
	dpAnniLength=0
	! initialize state
	j=iNPoint
	do i=1,iNPoint
		DPoints(i)%lStat=.false.
	enddo

	! go each point segment to check annihilation
	! Dpoints lStat is only set to true after annihilation implemented.
	do i=1,j
		if(IsCloseLoopPoint(i))then  !close loop stay theres
			DPoints(i)%lStat=.true.
		endif
		if(DPoints(i)%iEndP/=-1 .and. .not. DPoints(i)%lStat)then 
			call ImplementAnni(i)
		endif
	enddo

	!
	call Elimination

	!
	call RegroupDP

	!	
	do i=1,iNPoint
		call setGlobal(i)
	enddo
end subroutine NewAnni
!******************************************************
!	ImplementAnni
!		Check for the input segment, if annihilation
!	conditions are satisfied, perform implementation.
!*******************************************************
subroutine ImplementAnni(PointID)
	use ZQMPI
	use vectors
	use variables,only:DPoints,iloop
	use ANNIHILATIONMD,only:iAnniCount,lLogAnni,dpAnniLength,dpAnniDis
	implicit none

	integer,intent(in)::PointID
	integer::i,j,k
	double precision::u(2)
	logical::lAnni
	integer::iAnniIndex

	logical::lHaveSameBur		! if have same burgers vector
	logical::lSameLoopSeg	! if on smae loop
	logical::lConnectedSeg	! if connected loop
	double precision::dpShortDis	!check short distance

	lAnni=.false.
	dpShortDis=1000.d0
	iAnniIndex=-1

	do i=1,DPoints(PointID)%iNumNei
		j=DPoints(PointID)%iaIDOfNei(i)
		if(.not. DPoints(j)%lStat)then ! have not been surfed yet.
				call CheckSegDisAngle(PointID,j,lAnni,dpShortDis,iAnniIndex,U)
		endif
	enddo

	if(lAnni)then
		iAnniCount=iAnniCount+1
		if(.not. lLogAnni .and. lAnni)lLogAnni=lAnni
		call ImplementAnniL(PointID,iAnniIndex,U)
	endif
end subroutine ImplementAnni
!**********************************************************
!	CheckSegDisAngle:
!		check distance between segments and remember the dpshortDis,
!	if shorter than critical distance, then check angles. 
!	lAnni=true only if distance and angle both satisfy the condition.
!
!	4 return cases:
!		not equivalent planes, not same burgers vectors, connected 
!	segments and end-points.
!
!	!--------------------------
!	!4/15/08:
!	If two are within a distance that with current velocity, they will annihilate next step, annihilation will
!	also be implemented this step.
!**********************************************************
subroutine CheckSegDisAngle(id11,id21,lAnni,dpShortDis,iAnniIndex,dpAU)
	use ZQMPI
	use vectors
	use variables,only:DPoints,N_Times, DTIME,GPlanes, Miller,iloop
	use AnnihilationMD,only:dpAnniDis,dpAnniAngle
	implicit none

	! arguments
	integer,intent(in)::id11,id21
	logical::lAnni
	double precision::dpShortDis
	double precision,dimension(2)::dpAU
	integer::iAnniIndex

	double precision::VectorAngle
	! local variables
	integer::i,j
	integer::id12,id22
	integer::i_p1,i_p2,i_l1,i_l2,i_m1,i_m2
	double precision::shap(4,3),u01,u02
	type(vector)::PL(2,2),TL(2,2)
	type(vector)::P1,t1,p2,t2,pn1,pn2
	double precision::dpAngle,dpDis,dpDis2
	double precision::getDist
	logical:: lInEquivPList
	logical::lHaveSameBur		! if have same burgers vector
	logical::lSameLoopSeg	! if on smae loop
	logical::lConnectedSeg	! if connected loop
	double precision::dpFComputeLength

	!4/15/08
	type(vector)::dpApproxV1,dpApproxV2
	type(vector)::PV(2,2),TV(2,2)
	integer::itmpGaussP1,itmpGaussP2  !here to determine how many gauss points the code should use to check the distance
					  !between two segments. The segments should be divided into sub-segments at the length
					  !of dpAnniDis (the annihilation distance)
	!-----------------------------------------------------------
	
	id12=DPoints(id11)%iEndP
	id22=DPoints(id21)%iEndP

	i_p1=DPoints(id11)%iPlane
	i_p2=DPoints(id21)%iPlane

	!4/15/08
	i_m1=GPlanes(i_p1)%iMiller
	i_m2=GPlanes(i_p2)%iMiller

	!//////////////////////////////////////////////
	! if not equivalent plane, return
	if(.not. lInEquivPList(i_p1,i_p2))return

	!//////////////////////////////////////////////
	! here to avoid repeated computing, set id11 <=id21 to perform.
	if(id12==-1.or. id22==-1 .or. id11>id21)return

	if(lConnectedSeg(id11,id21))return

	if(.not.lHaveSameBur(id11,id21))return

!	if(.not. lSameLoopSeg(id11,id21))return


	!/////////////////////////////////////////////////////
	!Initial Position values.
!	DPoints(id21)%lStat=.true.
	!
	PL(1,1)=DPoints(id11)%tvPG
	TL(1,1)=DPoints(id11)%tvTG
	call TransNextGlobal(id11,PL(1,2),TL(1,2))
	!4/15/08
	PV(1,1)=DPoints(id11)%tvPreV;	PV(1,1)%v(3)=0.0
	TV(1,1)=DPoints(id11)%tvPreVT;	TV(1,1)%v(3)=0.0
	PV(1,2)=DPoints(id12)%tvPreV;	PV(1,2)%v(3)=0.0
	TV(1,2)=DPoints(id12)%tvPreVT;	TV(1,2)%v(3)=0.0

	PL(2,1)=DPoints(id21)%tvPG
	TL(2,1)=DPoints(id21)%tvTG
	call TransNextGlobal(id21,PL(2,2),TL(2,2))
	!4/15/08
	PV(2,1)=DPoints(id21)%tvPreV;	PV(2,1)%v(3)=0.0
	TV(2,1)=DPoints(id21)%tvPreVT;	TV(2,1)%v(3)=0.0
	PV(2,2)=DPoints(id22)%tvPreV;	PV(2,2)%v(3)=0.0
	TV(2,2)=DPoints(id22)%tvPreVT;	TV(2,2)%v(3)=0.0

	!////////////////////////////////////////////////
	!First determine tmpGaussP
	itmpGaussP1=int(dpFComputeLength(id11)/dpAnniDis)
	itmpGaussP2=int(dpFComputeLength(id21)/dpAnniDis)

	!///////////////////////////////////////////////
	!Check distance and angle
	do i=1,itmpGaussP1
		u01=(i-1)*(1.0/(itmpGaussP1-1))
		call getshape(shap,u01)
		p1=shap(1,1)*PL(1,1)+shap(2,1)*TL(1,1)+&
			shap(3,1)*PL(1,2)+shap(4,1)*TL(1,2)
		t1=shap(1,2)*PL(1,1)+shap(2,2)*TL(1,1)+&
			shap(3,2)*PL(1,2)+shap(4,2)*TL(1,2)

		!4/15/08, calculate velocity vector and determine next possible position
		dpApproxV1=shap(1,1)*PV(1,1)+shap(2,1)*TV(1,1)+&
			   shap(3,1)*PV(1,2)+shap(4,1)*TV(1,2)
		pn1=p1+N_TIMES*DTIME*Trans(Miller(i_m1)%ES)*dpApproxV1

		do j=1,itmpGaussP2
			u02=(j-1)*(1.0/(itmpGaussP2-1))
			call getshape(shap,u02)
			p2=shap(1,1)*PL(2,1)+shap(2,1)*TL(2,1)+&
				shap(3,1)*PL(2,2)+shap(4,1)*TL(2,2)
			t2=shap(1,2)*PL(2,1)+shap(2,2)*TL(2,1)+&
				shap(3,2)*PL(2,2)+shap(4,2)*TL(2,2)

			!4/15/08
			dpApproxV2=shap(1,1)*PV(2,1)+shap(2,1)*TV(2,1)+&
			   shap(3,1)*PV(2,2)+shap(4,1)*TV(2,2)
			pn2=p2+N_TIMES*DTIME*Trans(Miller(i_m2)%ES)*dpApproxV2
			
				
			!//////////////////////////
			! check distance
			dpDis=getDist(p1,p2)
			dpDis2=getDist(pn1,pn2)
			if(dpDis<dpAnniDis .or. dpDis2<dpAnniDis)then
				dpAngle=VectorAngle(t1,t2)
				!Condition satisfied
				if((dpDis<dpShortDis .or. dpDis2<dpShortDis).and. dpAngle>dpAnniAngle)then
				!if(dpDis<dpShortDis)then
					lAnni=.true.
					dpShortDis=dpDis
					if(dpDis2<dpDis)dpShortDis=dpDis2
					iAnniIndex=id21
					dpAU(1)=u01
					dpAU(2)=u02
				endif
			endif
		enddo
	enddo

end subroutine CheckSegDisAngle
!***************************************************************
!	ImplementAnniL
!		Real implementation of annihilation.
!		1. If two same plane, connect different points
!		2. If two parallel plane, move them onto one.
!		3. If two intersecting planes, added one intersection point
!****************************************************************
subroutine ImplementAnniL(id11,id21,U)
	use ZQMPI
	use vectors
	use variables,only:DPoints,iNPoint,GPlanes,MAX_Plane,MAX_Node,iDPTypeFixed
	use AnnihilationMD,only:dpAnniLength
	implicit none

	integer,intent(in)::id11,id21
	double precision,dimension(2)::U

	integer::iP1,iP2,im1,im2
	integer::id12,id22

	integer::i,j,k
	double precision::shap(4,3)
	double precision::u1,u2
	type(vector)::p1,t1,p2,t2,p3,t3,p4,t4
	type(vector)::PL(2),TL(2)
	type(vector)::PLL(2,2),TLL(2,2)
	double precision::getdist

	iP1=DPoints(id11)%iPlane
	ip2=DPoints(id21)%iPlane
	im1=GPlanes(ip1)%iMiller
	im2=GPlanes(ip2)%iMiller
	id12=DPoints(id11)%iEndP
	id22=DPoints(id21)%iEndP
	if(id12/=-1)DPoints(id12)%lStat=.true.
	if(id22/=-1)DPoints(id22)%lStat=.true.

	if(im1==im2)then !same miller plane, just change connection
		DPoints(id11)%iEndP=id22
		DPoints(id22)%iBeginP=id11
		DPoints(id21)%iEndP=id12
		DPoints(id12)%iBeginP=id21
		DPoints(id11)%lStat=.true.
		DPoints(id21)%lStat=.true.
		dpAnniLength=dpAnniLength+getDist(DPoints(id11)%tvPG,DPoints(id12)%tvPG)+  &
					  getDist(DPoints(id21)%tvPG,DPoints(id22)%tvPG)
	else ! intersecting plane
		return
		if(iNPoint<Max_Node*MAX_Plane)then
			!
			call GetInterPoint(id11,id21,U,PL)
			PLL(1,1)=PL(1)
			PLL(1,2)=PL(2)

			u1=U(1)
			U(1)=U(2)
			U(2)=u1

			call GetInterPoint(id21,id11,U,PL)
			PLL(2,1)=PL(1)
			PLL(2,2)=PL(2)

			! create new points
			iNPoint=iNPoint+2
			i=iNPoint-1
			DPoints(i)%iID=i
			DPoints(i)%iType=iDPTYpeFixed
			DPoints(i)%iPlane=ip2
			DPoints(i)%iBeginP=id11
			DPoints(i)%iEndP=id22
			DPoints(i)%lStat=.true.
			DPoints(i)%iNumNei=0
			DPoints(i)%iBurgers=DPoints(id11)%iBurgers
			DPoints(i)%tvPL=PLL(1,1)

			DPoints(id11)%iEndP=i
			DPoints(id22)%iBeginP=i
			DPoints(id11)%lStat=.true.

			!
			i=iNpoint
			DPoints(i)%iID=i
			DPoints(i)%iType=iDPTYpeFixed
			DPoints(i)%iPlane=ip1
			DPoints(i)%iBeginP=id21
			DPoints(i)%iEndP=id12
			DPoints(i)%lStat=.true.
			DPoints(i)%iNumNei=0
			DPoints(i)%iBurgers=DPoints(id21)%iBurgers
			DPoints(i)%tvPL=PLL(2,1)

			DPoints(id21)%iEndP=i
			DPoints(id12)%iBeginP=i
			DPoints(id21)%lStat=.true.
		endif
	endif
end subroutine ImplementAnniL

!********************************************************************
! GetInterPoint:
!	Get the position and tangent vectors for intersect point.
!	The global position is relative to point1.
!********************************************************************
subroutine GetInterPoint(id11,id21,U,PLI)
	use ZQMPI
	use vectors
	use variables,only:DPoints,GPlanes,Miller,zero
	use FunctionMD,only:DeterConnVecQ
	implicit none

	integer,intent(in)::id11,id21
	type(vector)::PLI(2)
	double precision,dimension(2)::U

	integer::id12,id22
	integer::i,j,k
	double precision::shap(4,3),u1,u2
	type(vector)::p1,p2,p3,p4,t1,t2,t3,t4
	integer::ip1,ip2,im1,im2,ip,im
	type(vector)::PL(2,2),TL(2,2),tmpVec,connectVec

	ip1=DPoints(id11)%iPlane
	ip2=DPoints(id21)%iPlane
	im1=GPlanes(ip1)%iMiller
	im2=GPlanes(ip2)%iMiller
	id12=DPoints(id11)%iEndP
	id22=DPoints(id21)%iEndP

	u1=U(1) !0.9*U(1)
	call getshape(shap,u1)
	PL(1,1)=DPoints(id11)%tvPG
	TL(1,1)=DPoints(id11)%tvTG
	call TransNextGlobal(id11,PL(1,2),TL(1,2))
	p1=shap(1,1)*PL(1,1)+shap(2,1)*TL(1,1)+shap(3,1)*PL(1,2)+shap(4,1)*TL(1,2)
	t1=shap(1,2)*PL(1,1)+shap(2,2)*TL(1,1)+shap(3,2)*PL(1,2)+shap(4,2)*TL(1,2)

	u1=U(2)+0.9*(1-U(2))
	call getshape(shap,u1)
	PL(2,1)=DPoints(id21)%tvPG
	TL(2,1)=DPoints(id21)%TVtG
	call TransNextGlobal(id21,PL(2,2),TL(2,2))
	p2=shap(1,1)*PL(2,1)+shap(2,1)*TL(2,1)+shap(3,1)*PL(2,2)+shap(4,1)*TL(2,2)
	t2=shap(1,2)*PL(2,1)+shap(2,2)*TL(2,1)+shap(3,2)*PL(2,2)+shap(4,2)*TL(2,2)

	!get global position, relative to point2
	connectVec=DeterConnVecQ(p1,p2)
	p1=p1-connectVec
	p3=0.5d0*(p1+p2)

	! loop to find the proper intersection point
	do i=1,10
		if(mod(i,2)==0)then
			ip=ip1
			im=im1
			tmpVec=connectVec
		else
			ip=ip2
			im=im2
			tmpVec=zero%v(1)
		endif
		!to local,adjust v(3), then back to global, if relative to p3, 
		! a transfer vector 'tmpvec' should be added into the formula.
		p4=Miller(im)%ES*(p3-GPlanes(ip)%Origin+tmpVec)
		p4%v(3)=0.d0
		p3=TRANS(Miller(im)%ES)*p4+GPlanes(ip)%Origin-tmpVec
	enddo

	! pay attention here: 2,1 are different.
	PLI(2)=Miller(im1)%ES*(p3-GPlanes(ip1)%Origin+connectVec)
	PLI(1)=Miller(im2)%ES*(p3-GPlanes(ip2)%Origin)

end subroutine
