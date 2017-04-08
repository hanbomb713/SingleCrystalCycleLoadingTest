!****************************************************************
!Start to work on March 30, 2006
!****************************************************************
!
! CrossSlip:
!	In these file contains subroutines dealing with cross-
!	slip.
!	
! Updated: ZQ Wang 
!	04/04/06: ok
!**********************************************************
subroutine CrossSlip
	use ZQMPI
	use vectors
	use variables,only:iNPoint
	use CrossslipMD,only:iCSCount,dpCSLength,iNumTotalScrew,iNumTotalScrew1
	implicit none
	
	external CrossSlipL
	integer::i
	!//////////////////////////////
	iNumTotalScrew=0
	iNumTotalScrew1=0
	iCSCount=0
	dpCSLength=0.0

	call SurfLoop(CrossSlipL)

	call Elimination

	do i=1,iNPoint
		call setGlobal(i)
	enddo

end subroutine CrossSlip
!**************************************************************
!	CrossSlipL:
!		Find the index of each loop and try to find the screw 
!	segment of this loop.
!
! Updated: ZQ Wang
!       04/04/06: ok. No need to change global system coordinates to 
!		local if they are used to look for the screw segments.
!
!	04/04/06:
!		Modified the way to calculate the length of the loop.
!	04/05/06: ok
!**************************************************************
subroutine CrossSlipL(PointID)
	use ZQMPI
	use vectors
	use variables,only:DPoints,iDPTypeFree,iloop
	use crossSlipMD,only:itmpIndex,itmpNumP,itmpNewNumP,tvTmpPL,&
						itmpLooptype,itmpBurgers,itmpPlane,lScrewSeg,&
						lScrewForce,lCrossSlip,lToImplement
	implicit none

	integer,intent(in)::PointID

	integer::i,j,k,looptype
	integer::iBegin,iEnd,iNext,iCur,iEnterID
	logical::lContinue,lSegCriteria
	logical::StopSurfLoop

	double precision::daLength
	type(vector),dimension(2)::PL,TL

	iEnterID=PointID
	iCur=iEnterID
	lCOntinue=.true.
	lSegCriteria=.false.
	itmpNumP=0
	lScrewSeg=.false.
	lScrewForce=.false.
	lCrossSlip=.false.
	lToImplement=.false.
	daLength=0.d0

	if(DPoints(iEnterID)%iType==iDPTypeFree)then
		looptype=0
	else
		looptype=1
	endif

	itmpLooptype=looptype
	itmpBurgers=DPoints(PointID)%iBurgers
	itmpPlane=DPoints(PointID)%iPlane

	!/////////////////////////////////////////////////////
	! Surf the loop and assign the index of point to itmpIndex
	! Here, we know all the points. 
	do while(lContinue)
		itmpNumP=itmpNumP+1
		itmpIndex(itmpNumP)=iCur
		lContinue=.not.StopSurfLoop(iCur,iEnterID,lSegCriteria)

		iNext=DPoints(iCur)%iEndP
		iCur=iNext
	enddo

	iTmpNewNumP=iTmpNumP
	! assign values	of global 
	do i=1,iTmpNumP
		tvTmpPL(i)=DPoints(iTmpIndex(i))

		!//////////////////////////////////
		! 4/4/06
		PL(1)=DPoints(iTmpIndex(i))%tvPL
		TL(1)=DPoints(iTmpIndex(i))%tvTL
		if(DPoints(iTmpIndex(i))%iEndP/=-1)then
			call TransNextLocal(iTmpIndex(i),PL(2),TL(2))
			daLength=daLength+MAG(PL(2)-PL(1))
		endif
	enddo

	!------------------------------------
	if(daLength>=1000.d0)then
		call getCrossSlipMiller
		call LookingForScrew
		call calculateProbability
		call implementCrossSlip
	endif

end subroutine CrossSlipL
!******************************************************************
!   getCrossSlipMiller                                            *
!         Get the cross-slip miller index                         *
!
!	!//////////////////////////////////////////////////////////
!	!FCC Crystals
!         Four kind of miller indexes. On each glide plane for    *
!         each burgers vector, there is only one cross-slip plane.*
!         Each glide plane has 3 possible Burgers vectors. Thus,  *
!         Each glide plane may have 3 possible cross-slip planes  *
!         According to the Burgers vectors.                       *
!
!	After this, we can calculate the cross-slip force.
! Updated: ZQ Wang
!       04/04/06: ok
!
!	!//////////////////////////////////////////////////////////
!	!BCC Crystals 
!	!For each BCC (110) slip plane, there are two possible cross slip
!	!planes for each possible Burgers vector, and their information should be recorded until the possibilities
!	! of cross-slip have been calculated and one specific plane has been selected.
!	!The information includes: millers, resolved forces.
!******************************************************************
subroutine getCrossSlipMiller
	use ZQMPI
	use vectors
	use Variables,only:Miller,GPlanes,Burgers,iNumMiller,iMatType,iBCC,iFCC
	use CrossSlipMD,only:itmpPlane,itmpBurgers,iCSMiller,inumCSM
	implicit none

	type(vector)::tmpMiller,CSMillerTMP
	type(vector)::tmpBurgers
	integer::i

	tmpMiller=Miller(GPlanes(itmpPlane)%iMiller)%v
	tmpBurgers=Burgers(itmpBurgers)

	!//////////////////////////////////////
	!FCC materials
	if(iMatType==iFCC)then
		iNumCSM=1
		if(tmpMiller%v(1)*tmpMiller%v(3)>0 .and. tmpMiller%v(2)*tmpMiller%v(3)>0)then
		! this is miller index of (111)
			if(dabs(tmpBurgers%v(2))<1.d-5)then
				CSMillerTMP%v(1)=1.d0;CSMillerTMP%v(2)=-1.d0;CSMillerTMP%v(3)=1.d0
			else if(dabs(tmpBurgers%v(3))<1.d-5)then
				CSMillerTMP%v(1)=-1.d0;CSMillerTMP%v(2)=-1.d0;CSMillerTMP%v(3)=1.d0
			else if(dabs(tmpBurgers%v(1))<1.d-5)then
				CSMillerTMP%v(1)=-1.d0;CSMillerTMP%v(2)=1.d0;CSMillerTMP%v(3)=1.d0
			endif
		elseif(tmpMiller%v(1)*tmpMiller%v(3)>0 .and. tmpMiller%v(2)*tmpMiller%v(3)<0)then
		! this is miller index of (1-11)
			if(dabs(tmpBurgers%v(2))<1.d-5)then
				CSMillerTMP%v(1)=1.d0;CSMillerTMP%v(2)=1.d0;CSMillerTMP%v(3)=1.d0
			else if(dabs(tmpBurgers%v(1))<1.d-5)then
				CSMillerTMP%v(1)=-1.d0;CSMillerTMP%v(2)=-1.d0;CSMillerTMP%v(3)=1.d0
			else if(dabs(tmpBurgers%v(3))<1.d-5)then
				CSMillerTMP%v(1)=-1.d0;CSMillerTMP%v(2)=1.d0;CSMillerTMP%v(3)=1.d0
			endif
		elseif(tmpMiller%v(1)*tmpMiller%v(3)<0 .and. tmpMiller%v(2)*tmpMiller%v(3)<0)then
		! this is miller index of (-1-11)
			if(dabs(tmpBurgers%v(3))<1.d-5)then
				CSMillerTMP%v(1)=1.d0;CSMillerTMP%v(2)=1.d0;CSMillerTMP%v(3)=1.d0
			else if(dabs(tmpBurgers%v(1))<1.d-5)then
				CSMillerTMP%v(1)=1.d0;CSMillerTMP%v(2)=-1.d0;CSMillerTMP%v(3)=1.d0
			else if(dabs(tmpBurgers%v(2))<1.d-5)then
				CSMillerTMP%v(1)=-1.d0;CSMillerTMP%v(2)=1.d0;CSMillerTMP%v(3)=1.d0
			endif
		elseif(tmpMiller%v(1)*tmpMiller%v(3)<0 .and. tmpMiller%v(2)*tmpMiller%v(3)>0)then
		! this is miller index of (-111)
			if(dabs(tmpBurgers%v(1))<1.d-5)then
				CSMillerTMP%v(1)=1.d0;CSMillerTMP%v(2)=1.d0;CSMillerTMP%v(3)=1.d0
			else if(dabs(tmpBurgers%v(3))<1.d-5)then
				CSMillerTMP%v(1)=1.d0;CSMillerTMP%v(2)=-1.d0;CSMillerTMP%v(3)=1.d0
			else if(dabs(tmpBurgers%v(2))<1.d-5)then
				CSMillerTMP%v(1)=-1.d0;CSMillerTMP%v(2)=-1.d0;CSMillerTMP%v(3)=1.d0
			endif
		endif
		!Check the miller index
		do i=1,iNumMiller
			if(CSMillerTmp%v(1)==Miller(i)%v%v(1).and.CSMillerTmp%v(2)==Miller(i)%v%v(2)  &
								.and.CSMillerTmp%v(3)==Miller(i)%v%v(3))then
				iCSMiller(1)=i
				exit
			endif
		enddo	
	!///////////////////////////////////////
	! BCC materials: 
	! Current only for {110} system
	elseif(iMatType==iBCC)then
		iNumCSM=2
		if(tmpMiller%v(1)*tmpMiller%v(2)>0)then
		! this is miller index of (110)
			if(tmpBurgers%v(2)*tmpBurgers%v(3)>0)then
				!b=[1-1-1]
				iCSMiller(1)=3; iCSMiller(2)=6
			else
				!b=[1-11]
				iCSMiller(1)=4; iCSMiller(2)=5
			endif
		elseif(tmpMiller%v(1)*tmpMiller%v(2)<0)then
		! this is miller index of (1-10)
			if(tmpBurgers%v(2)*tmpBurgers%v(3)>0)then
				!b=[111]
				iCSMiller(1)=4; iCSMiller(2)=6
			else
				!b=[11-1]
				iCSMiller(1)=3; iCSMiller(2)=5
			endif
		elseif(tmpMiller%v(1)*tmpMiller%v(3)>0)then
		! this is miller index of (101)
			if(tmpBurgers%v(2)*tmpBurgers%v(3)>0)then
				!b=[1-1-1]
				iCSMiller(1)=1; iCSMiller(2)=6
			else
				!b=[11-1]
				iCSMiller(1)=2; iCSMiller(2)=5
			endif
		elseif(tmpMiller%v(1)*tmpMiller%v(3)<0)then
		! this is miller index of (10-1)
			if(tmpBurgers%v(2)*tmpBurgers%v(3)>0)then
				!b=[111]
				iCSMiller(1)=2; iCSMiller(2)=6
			else
				!b=[1-11]
				iCSMiller(1)=1; iCSMiller(2)=5
			endif
		elseif(tmpMiller%v(2)*tmpMiller%v(3)>0)then
		! this is miller index of (011)
			if(tmpBurgers%v(1)*tmpBurgers%v(3)>0)then
				!b=[1-11]
				iCSMiller(1)=1; iCSMiller(2)=4
			else
				!b=[11-1]
				iCSMiller(1)=2; iCSMiller(2)=3
			endif
		elseif(tmpMiller%v(2)*tmpMiller%v(3)<0)then
		! this is miller index of (01-1)
			if(tmpBurgers%v(1)*tmpBurgers%v(3)>0)then
				!b=[111]
				iCSMiller(1)=2; iCSMiller(2)=4
			else
				!b=[1-1-1]
				iCSMiller(1)=1; iCSMiller(2)=3
			endif
		endif
	
	endif

end subroutine getCrossSlipMiller
!******************************************************************
!   LookingForScrew                                               *
!	If the angle between the Burgers vector and the line vector
!	is less than 15 degree, we say that the segment is with 
!	screw character.
!
! Updated: ZQ Wang
!       04/04/06: Modified the way to set the starting point to look
!		  for screw segment on close loops.
!	04/05/06: ok
!
!******************************************************************
subroutine LookingForScrew
	use ZQMPI
	use vectors
	use variables,only:PI,DPoints,Burgers
	use CrossSlipMD,only:NumScrew,itmpLooptype,&
						itmpNumP,itmpIndex,itmpBurgers
	implicit none

	integer::i,j,k,iSt
	double precision::u0,shap(4,3),u1,getdist
	type(vector)::R
	type(vector),dimension(10)::P
	type(vector),dimension(2)::tmpP,tmpT
	type(vector)::tmpBurgers
	logical::beforeIsScrew
	integer::screwCounter,segN1,segN0
	
	NumScrew=0
	tmpBurgers=Burgers(itmpBurgers)

	beforeIsScrew=.false.
	screwCounter=0
	if(itmpLooptype==1)then  ! open loop from beginning
		do i=1,itmpNumP-itmpLooptype
			if(screwCounter==0 .and.i==1)then
				P(1)=DPoints(itmpIndex(i))%tvPG
				u0=0.0
				segN0=i
			endif

			tmpP(1)=DPoints(itmpIndex(i))%tvPG
			tmpT(1)=DPoints(itmpIndex(i))%tvTG
			call TransNextGlobal(itmpIndex(i),tmpP(2),tmpT(2))

			do j=1,5
				u1=0.2*j
				call getshape(shap,u1)
				R=shap(1,1)*tmpP(1)+shap(2,1)*tmpT(1)+&
				  shap(3,1)*tmpP(2)+shap(4,1)*tmpT(2)
				P(2)=R
				u1=u1
				segN1=i
				!05/27/06......if too close to ending point, just neglect it because of numerical instability
				if(getDist(R,DPoints(itmpIndex(1))%tvPG)<=300.d0 .or. getDist(R,DPoints(itmpIndex(itmpNumP))%tvPG)<=300.d0)then
					beforeISScrew=.false.
				else
					call DetermineScrew(P(1),P(2),beforeIsScrew,screwCounter,u0,u1,segN0,segN1)
				endif
				u0=u1
			        segN0=segN1
			        P(1)=P(2)
			enddo
		enddo
	else   ! close loop should be carefully considered to find correct screw segment
		! if first point is screw

		!/////////////////////////////////////////////////
		! To see if the beginning point belongs to a screw. If yes, move the starting point from somewhere else.
		u0=0.2
		tmpP(1)=DPoints(itmpIndex(1))%tvPG
		tmpT(1)=DPoints(itmpIndex(1))%tvTG
		call TransNextGlobal(itmpIndex(1),tmpP(2),tmpT(2))
		call getshape(shap,u0)
		P(2)=shap(1,1)*tmpP(1)+shap(2,1)*tmpT(1)+shap(3,1)*tmpP(2)+shap(4,1)*tmpT(2)
		P(1)=tmpP(1)

		!.....set the start point....
		if(dabs(UV(P(2)-P(1))*UV(tmpBurgers))>dcos(0.08334*PI))then ! is screw
			iSt=itmpNumP*3/4+1
		else
			iSt=1
		endif
		!///////////////////////////////////////////////////

		if(screwCounter==0)then   ! set the first segment to begin searching
			k=iSt
			P(1)=DPoints(itmpIndex(k))%tvPG
			u0=0.0
			segN0=k
		endif

		do i=1,itmpNumP
			k=iSt+i-1
			if(k>itmpNumP)then  
				k=mod(k,itmpNumP)
			endif

			tmpP(1)=DPoints(itmpIndex(k))%tvPG
			tmpT(1)=DPoints(itmpIndex(k))%tvTG
			tmpP(2)=DPoints(itmpIndex(mod(k,itmpNumP)+1))%tvPG
			tmpT(2)=DPoints(itmpIndex(mod(k,itmpNumP)+1))%tvTG

			do j=1,5
				u1=0.2*j
				call getshape(shap,u1)
				R=shap(1,1)*tmpP(1)+shap(2,1)*tmpT(1)+&
				  shap(3,1)*tmpP(2)+shap(4,1)*tmpT(2)
				P(2)=R
				u1=u1
				segN1=k
				call DetermineScrew(P(1),P(2),beforeIsScrew,screwCounter,u0,u1,segN0,segN1)
				u0=u1
			        segN0=segN1
			        P(1)=P(2)
			enddo
		enddo
	endif

end subroutine LookingForScrew

!******************************************************************
!  DetermineScrew:                                                *
!         Determine if a segment is screw.                        *
!         Record the segment no. and u0                           *
!        
!         If the line sense vector is aligned within 15 degree 
!	with Burgers vector, it is a screw segment.
!
! Updated: ZQ Wang
!       04/04/06: ok
!
!******************************************************************
subroutine DetermineScrew(P1,P2,beforeIsScrew,screwCounter,u0,u1,segN0,segN1)
	use ZQMPI
	use vectors
	use variables,only:PI,Burgers
	use CrossSlipMD,only:itmpBurgers,UScrew,SegScrew,NumScrew,tvLastScrewPoint
	implicit none

	type(vector)::P1,P2,R1,R2
	logical::beforeIsScrew
	integer::screwCounter,segN0,segN1
	double precision::u0,u1,ang,getDist

	R1=UV(P2-P1)
	R2=UV(Burgers(itmpBurgers))

	if(NumScrew<=3)then  !only if NumScrew less than 3
		ang=R1*R2
		if(dabs(ang)>dcos(0.08334*PI))then  ! is screw segment
		
			if(beforeIsScrew)then ! screw continue
				UScrew(NumScrew,2)=u1          !ending of the screw
				SegScrew(NumScrew,2)=segN1
			elseif(NumScrew>0 .and. getDist(P2,tvLastScrewPoint)<300.d0)then    !screw segment, to see if two screw segments are too close. 
			!6/6/06								    !If yes, connect them, otherwise create new.
				UScrew(NumScrew,2)=u1
				SegScrew(NumScrew,2)=segN1
			else  ! create new screw and new nonscrew
				NumScrew=NumScrew+1
				UScrew(NumScrew,1)=u0
				UScrew(NumScrew,2)=u1
				SegScrew(NumScrew,1)=segN0
				SegScrew(NumScrew,2)=segN1
			endif
			beforeIsScrew=.true.
			tvLastScrewPoint=P2   !6/6/06
		else
			beforeIsScrew=.false.
		endif
	
	endif
	screwCounter=screwCounter+1
end subroutine DetermineScrew
!******************************************************************
!   CalculateProbability                                          *
!         Calculate probability of screw segment to cross-slip.   *
!         Rearrange the segments.                                 *
!
!
! Updated: ZQ Wang
!       04/04/06: When calculate the force, compare the force on 
!	cross-slip plane and original-slip plane. If cross-slip
!	force is smaller than the original force, there is no need
!	to calculate the probability and no need to implement the 
!	cross-slip. 
!
!	04/05/06: working on this subroutine now.
!	04/08/06:
!	Here, CSYN is set to false if: 1). force on cross-slip
!	plane is smaller than that on the original plane; 2). the probability
!	calculated does not allow dislocation to cross-slip.
!******************************************************************
subroutine CalculateProbability
	use ZQMPI
	use vectors
	use CrossSlipMD,only:FinalScrewN,NumScrew
	implicit none

	integer::i,j,k
	double precision::TauCS
	double precision::LengthCS,LengthCSNew

	FinalScrewN=0
!	print *,"iNumScrew:",NumScrew
	do i=1,NumScrew
		call CSForceLength(i,TauCS,LengthCS,LengthCSNew)

		call ComputeProbability(i,TauCS,LengthCS,LengthCSNew)
		
		call ScrewRearrangement(i)
	enddo
end subroutine CalculateProbability
!******************************************************************
!    CSForceLength:                                               *
!       Get the force on cross-slip plane and length of the screw *
!       segments.  Also, the miller index of the cross-slip plane *
!       is got.                                                   *
!!
! Updated: ZQ Wang
!       04/04/06: Compare forces before return. if cross-slip force
!	is small, no need for further calculation.
!	4/12/06:
!	tmpForce is the force on such a length, average force should
!	be divided by the length.
!
!******************************************************************
subroutine CSForceLength(ScrewID,TauCS,LengthCS,LengthCSNew)
	use ZQMPI
	use vectors
	use variables,only:DPoints,Miller,GPLanes,zero,Burgers,MU
	use CrossSlipMD,only:segScrew,itmpIndex,itmpBurgers,&
	                      UScrew,itmpNumP,itmpLooptype,CSYN,iCSMiller,iNumCSM
	implicit none
	
	integer::ScrewID
	double precision::TauCS,LengthCS, LengthCSNew,iDealFR
	double precision::dpMAGCF,dpMAGOF,dpTmpForceCom
	integer::i,j,k,iNext,iTmpForceCom
	double precision::shap(4,3)
	type(vector)::P1,P2
	type(vector),dimension(2)::tmpP,tmpT,tvScrewEnd
	type(vector)::tmpForce,tmpOForce,tmpCForce(3)   !force on original plane (tmpOForce) and on cross-slip plane (tmpCForce)
	double precision::GetDist

	CSYN=.false.

	!//////////////////////////////////////////////////////
	! the length
	LengthCS=0.d0
	LengthCSNew=0.d0

	if(segScrew(ScrewID,1)<segScrew(ScrewID,2))then
		!1.
		call getshape(shap,UScrew(ScrewID,1))
		i=segScrew(sCrewID,1)
		tmpP(1)=DPoints(itmpIndex(i))%tvPL
		tmpT(1)=DPoints(itmpIndex(i))%tvTL
		call TransNextLocal(itmpIndex(i),tmpP(2),tmpT(2))
		LengthCSNew=LengthCSNew+MAG(tmpP(1)-tmpP(2))
		
		tmpP(1)=DPoints(itmpIndex(i))%tvPG
		tmpT(1)=DPoints(itmpIndex(i))%tvTG
		call TransNextGlobal(itmpIndex(i),tmpP(2),tmpT(2))
		P1=shap(1,1)*tmpP(1)+shap(2,1)*tmpT(1)+shap(3,1)*tmpP(2)+shap(4,1)*tmpT(2)
		tvScrewEnd(1)=P1
		!2.
		do i=segScrew(sCrewID,1)+1,segScrew(sCrewID,2)-1
			tmpP(1)=DPoints(itmpIndex(i))%tvPL
			tmpT(1)=DPoints(itmpIndex(i))%tvTL
			call TransNextLocal(itmpIndex(i),tmpP(2),tmpT(2))
			LengthCSNew=LengthCSNew+MAG(tmpP(1)-tmpP(2))
		enddo
	
		!3.
		call getshape(shap,UScrew(ScrewID,2))
		i=segScrew(sCrewID,2)
		tmpP(1)=DPoints(itmpIndex(i))%tvPL
		tmpT(1)=DPoints(itmpIndex(i))%tvTL
		call TransNextLocal(itmpIndex(i),tmpP(2),tmpT(2))	
		LengthCSNew=LengthCSNew+MAG(tmpP(1)-tmpP(2))

		tmpP(1)=DPoints(itmpIndex(i))%tvPG
		tmpT(1)=DPoints(itmpIndex(i))%tvTG
		call TransNextGlobal(itmpIndex(i),tmpP(2),tmpT(2))
		P2=shap(1,1)*tmpP(1)+shap(2,1)*tmpT(1)+shap(3,1)*tmpP(2)+shap(4,1)*tmpT(2)
		tvScrewEnd(2)=P2
	elseif(segScrew(ScrewID,1)==segScrew(ScrewID,2))then

		i=segScrew(sCrewID,1)
		tmpP(1)=DPoints(itmpIndex(i))%tvPL
		tmpT(1)=DPoints(itmpIndex(i))%tvTL
		call TransNextLocal(itmpIndex(i),tmpP(2),tmpT(2))

		!1
		call getshape(shap,UScrew(ScrewID,1))
		P1=shap(1,1)*tmpP(1)+shap(2,1)*tmpT(1)+shap(3,1)*tmpP(2)+shap(4,1)*tmpT(2)
		!2
		call getshape(shap,UScrew(ScrewID,2))
		P2=shap(1,1)*tmpP(1)+shap(2,1)*tmpT(1)+shap(3,1)*tmpP(2)+shap(4,1)*tmpT(2)
		
		LengthCSNew=LengthCSNew+MAG(P1-P2)

		!----1
		i=segScrew(sCrewID,1)
		tmpP(1)=DPoints(itmpIndex(i))%tvPG
		tmpT(1)=DPoints(itmpIndex(i))%tvTG
		call TransNextGlobal(itmpIndex(i),tmpP(2),tmpT(2))

		call getshape(shap,UScrew(ScrewID,1))
		P1=shap(1,1)*tmpP(1)+shap(2,1)*tmpT(1)+shap(3,1)*tmpP(2)+shap(4,1)*tmpT(2)
		call getshape(shap,UScrew(ScrewID,2))
		P2=shap(1,1)*tmpP(1)+shap(2,1)*tmpT(1)+shap(3,1)*tmpP(2)+shap(4,1)*tmpT(2)
		
		tvScrewEnd(1)=P1
		tvScrewEnd(2)=P2
	else
		!1
		call getshape(shap,UScrew(ScrewID,1))
		i=segScrew(sCrewID,1)
		tmpP(1)=DPoints(itmpIndex(i))%tvPL
		tmpT(1)=DPoints(itmpIndex(i))%tvTL
		call TransNextLocal(itmpIndex(i),tmpP(2),tmpT(2))
		LengthCSNew=LengthCSNew+MAG(tmpP(1)-tmpP(2))

		tmpP(1)=DPoints(itmpIndex(i))%tvPG
		tmpT(1)=DPoints(itmpIndex(i))%tvTG
		call TransNextGlobal(itmpIndex(i),tmpP(2),tmpT(2))
		P1=shap(1,1)*tmpP(1)+shap(2,1)*tmpT(1)+shap(3,1)*tmpP(2)+shap(4,1)*tmpT(2)
		tvScrewEnd(1)=P1

		!2
		do i=segScrew(ScrewID,1)+1,itmpNumP-itmpLooptype
			tmpP(1)=DPoints(itmpIndex(i))%tvPL
			tmpT(1)=DPoints(itmpIndex(i))%tvTL
			call TransNextLocal(itmpIndex(i),tmpP(2),tmpT(2))
			LengthCSNew=LengthCSNew+MAG(tmpP(1)-tmpP(2))
		enddo
		do i=1,segScrew(ScrewID,2)-1
			tmpP(1)=DPoints(itmpIndex(i))%tvPL
			tmpT(1)=DPoints(itmpIndex(i))%tvTL
			call TransNextLocal(itmpIndex(i),tmpP(2),tmpT(2))
			LengthCSNew=LengthCSNew+MAG(tmpP(1)-tmpP(2))
		enddo
		!3
		call getshape(shap,UScrew(ScrewID,2))
		i=segScrew(sCrewID,2)
		tmpP(1)=DPoints(itmpIndex(i))%tvPL
		tmpT(1)=DPoints(itmpIndex(i))%tvTL
		call TransNextLocal(itmpIndex(i),tmpP(2),tmpT(2))
		LengthCSNew=LengthCSNew+MAG(tmpP(1)-tmpP(2))

		tmpP(1)=DPoints(itmpIndex(i))%tvPG
		tmpT(1)=DPoints(itmpIndex(i))%tvTG
		call TransNextGlobal(itmpIndex(i),tmpP(2),tmpT(2))
		P2=shap(1,1)*tmpP(1)+shap(2,1)*tmpT(1)+shap(3,1)*tmpP(2)+shap(4,1)*tmpT(2)
		tvScrewEnd(2)=P2
	endif
	LengthCS=GetDist(tvScrewEnd(1),tvScrewEnd(2))
!	print *,"LengthCS:",LengthCS,",NewL:",LengthCSNew,LengthCS/LengthCSNew

	!////////////////////////////////////////////////////////////////////////
	! the force. taForce is global one. we should change it to the cross-slip 
	! plane and the original plane
	! 04/08/06:
	! Here we use an approximation: the average force is simply the total force divided by
	! the number of nodes. The actual force should be the nodal forces time the length and 
	! divided by the total length.
	k=0

	tmpForce=zero%v(1)
	do j=1,iNumCSM
		tmpCForce(j)=zero%v(1)
	enddo
	tmpOForce=zero%v(1)
	if(SegScrew(ScrewID,1)>SegScrew(ScrewID,2))then
		do j=SegScrew(ScrewID,1),itmpNumP-itmpLooptype
			k=k+1
			tmpForce=tmpForce+DPoints(itmpIndex(j))%force
		enddo
		do j=1,SegScrew(ScrewID,2)
			k=k+1
			tmpForce=tmpForce+DPoints(itmpIndex(j))%force
		enddo
	else
		do j=SegScrew(ScrewID,1),SegScrew(ScrewID,2)
			k=k+1
			tmpForce=tmpForce+DPoints(itmpIndex(j))%force
		enddo
	endif

!--------!////////////////////////////////////////////////////////////////
	!Force should be averaged by dividing it by the length
	tmpForce=1./LengthCSNew*tmpForce   !Global force

	!.....Transfer to in-plane local force
	tmpOForce=Trans(Miller(GPlanes(DPoints(itmpIndex(2))%iPlane)%iMiller)%ES)*tmpForce
	tmpOForce%v(3)=0
	
	dpTmpForceCom=0.0
	iTmpForceCom=1
	do j=1,iNumCSM
		tmpCForce(j)=TRANS(Miller(iCSMiller(j))%ES)*TmpForce
		tmpCForce(j)%v(3)=0
		if(MAG(tmpCForce(j))>dpTmpForceCom)then    !Find the largest force, choose that miller
			dpTmpForceCom=MAG(tmpCForce(j))
			iTmpForceCom=j
		endif
	enddo

	!....Force perpendecular to the dislocation
	!New: with BCC: dpMAGCF should be the largest one of iNumCSM cross-slip forces: fcc: 1, bcc: 2-3
	dpMAGCF=dsqrt(MAG(tmpCForce(iTmpForceCom))**2-(tmpCForce(iTmpForceCom)*UV(Burgers(itmpBurgers)))**2)
	dpMAGOF=dsqrt(MAG(tmpOForce)**2-(tmpOForce*UV(Burgers(itmpBurgers)))**2)

!-------!////////////////////////////////////////////////////////////////////////////




!	print *,"CF/OF:",dpMAGCF*1d-6,dpMAGOF*1d-6,dpMAGCF/dpMAGOF
	!.....
	! If the force on cross-slip plane is larger than that on original plane, set CSYN to true.
	if(dpMAGCF*1d-6>dpMAGOF*1d-6+1 .and. LengthCS>300.d0)then
		
		TauCS=dpMAGCF*1.414*1d-6    !unit of MPa. P-K force*1.414 --->corresponding stress

		!ideal strength for FR
		iDealFR=MU*dsqrt(2.d0)/2/LengthCS*1d-6   !mu*b/L, unit of MPa
	
!		print *,"TauCS:",TauCS,",iDeal:",iDealFR,TAUCS/iDealFR

		!if real force is larger than critical ideal strength of FR, check force distribution
		if(TauCS>iDealFR*0.6)then  ! take alpha=0.5
		
			j=0
			k=0
			if(SegScrew(ScrewID,1)>SegScrew(ScrewID,2))then
				do i=SegScrew(ScrewID,1),itmpNumP-itmpLooptype
					tmpForce=DPoints(itmpIndex(i))%force
					tmpP(1)=DPoints(itmpIndex(i))%tvPL
					tmpT(1)=DPoints(itmpIndex(i))%tvTL
					call TransNextLocal(itmpIndex(i),tmpP(2),tmpT(2))
					tmpForce=TRANS(Miller(iCSMiller(iTmpForceCom))%ES)*(1./(MAG(tmpP(1)-tmpP(2)))*tmpForce)
					tmpForce%v(3)=0
					if(dsqrt(MAG(tmpForce)**2-(tmpForce*UV(Burgers(itmpBurgers)))**2)>TauCS*1d6/1.414)then
						j=j+1
					else
						k=k+1
					endif		
				enddo
				do i=1,SegScrew(ScrewID,2)
					tmpForce=DPoints(itmpIndex(i))%force
					tmpP(1)=DPoints(itmpIndex(i))%tvPL
					tmpT(1)=DPoints(itmpIndex(i))%tvTL
					call TransNextLocal(itmpIndex(i),tmpP(2),tmpT(2))
					tmpForce=TRANS(Miller(iCSMiller(iTmpForceCom))%ES)*(1./(MAG(tmpP(1)-tmpP(2)))*tmpForce)
					tmpForce%v(3)=0
					if(dsqrt(MAG(tmpForce)**2-(tmpForce*UV(Burgers(itmpBurgers)))**2)>TauCS*1d6/1.414)then
						j=j+1
					else
						k=k+1
					endif		
				enddo
			else
				do i=SegScrew(ScrewID,1),SegScrew(ScrewID,2)
					tmpForce=DPoints(itmpIndex(i))%force
					tmpP(1)=DPoints(itmpIndex(i))%tvPL
					tmpT(1)=DPoints(itmpIndex(i))%tvTL
					call TransNextLocal(itmpIndex(i),tmpP(2),tmpT(2))
					tmpForce=TRANS(Miller(iCSMiller(itmpForceCom))%ES)*(1./(MAG(tmpP(1)-tmpP(2)))*tmpForce)
					tmpForce%v(3)=0
					if(dsqrt(MAG(tmpForce)**2-(tmpForce*UV(Burgers(itmpBurgers)))**2)>TauCS*1d6/1.414)then
						j=j+1
					else
						k=k+1
					endif		
				enddo
			endif
!			print *,"j,k:",j,k
			if(j>=k)then
				CSYN=.true.   !now criteria is satisfied: TauCS > Critical Stress, CF>OF, Distribution OK.
			endif
		endif
	endif

	!record the miller with the largest force as the cross-slip miller, stored in the first element of the array
	iCSMiller(1)=iCSMiller(itmpForceCom)
end subroutine CSForceLength
!******************************************************************
!   ComputeProbability                                            *
!
!	Updated: ZQ Wang 
!	04/04/06: Need to know the probablity.
!
!	04/08/06:
!	1). If CSYN from CSForceLength is not true, that means that the force
!	on cross-slip plane is not large enough.
!	2). If Yes, calculate probability
!	3). At the end, check if CSYN is true. If any CSYN for each segment is 
!	true, global lToImplement should be set as true.
!******************************************************************
subroutine ComputeProbability(ScrewID,TauCS,LengthCS)
	use ZQMPI
	use vectors
	use variables,only:iloop,MU,iMatType,iFCC,iBCC,DTime,N_Times
	use CrossSlipMD,only:Length0,TAUIII,CSYN,lToImplement,iNumTotalScrew,iNumTotalScrew1
	use functionMD,only:funcBCCCSL
	implicit none

	integer,intent(in)::ScrewID
	double precision,intent(in)::TauCS,LengthCS
	double precision::Pcs,dpRandomPcs
	double precision::dpTmpSH

	if(CSYN)then
		if(iMatType==iBCC)then
			Length0=500.0
			TauIII=funcBCCCSL(length0)
		endif
		Pcs=5*N_Times*DTime*LengthCS/Length0*dexp(1.21*(TauCS-TauIII))

		!/////////////////////////////////////
		! Set a criterion to activate cross-slip
		call RandomGenerator(dpRandomPcs)
		if(Pcs<dpRandomPcs)then
			CSYN=.false.   !not crossslip
		endif
	endif

	!/////////////////////////////////////
	! set global lToImplement to be true if any CSYN is true
	if(CSYN)then
		iNumTotalScrew=iNumTotalScrew+1
		if(Pcs>=1)iNumTotalScrew1=iNumTotalScrew1+1
		lToImplement=.true.
!		print *,"lToImplement=",lToImplement
	endif
	
end subroutine ComputeProbability
!******************************************************************
!  RandomGenerator:
!	there will be different versions.
!
!	Updated: ZQ Wang
!	4/12/06:
!	Here to generate a random number that will be used to select
!	the cross-slip probability.
!******************************************************************
subroutine RandomGenerator(dpRandomPcs)
	use ZQMPI
	implicit none

	double precision::dpRandomPcs

	call random_number(dpRandomPcs)


end subroutine RandomGenerator


!******************************************************************
!   ScrewRearrangement:                                           *
!         Here we will get final value of non-cross-slip(NCS)     *
!         number. This number will be used to assign the last NCS *
!         to old dislocation.                                     *
!
! Updated: ZQ Wang
! 	4/4/06: need more work
!
!	4/8/06:
!	1). If this screw segment is going to cross-slip, add it to the
!	final list of segments that are going to cross-slip;
!	2). If this screw segment is not going to cross-slip, remove it from 
!	the list and change the list of non-cross-slip segment.
!	3). Still working on.
!******************************************************************
subroutine ScrewRearrangement(i)
	use ZQMPI
	use vectors
	use CrossSlipMD,only:FinalScrewN,FinalScrewSeg,SegScrew,FinalScrewU,&
						UScrew,FinalScrewMiller,iCSMiller,itmpLooptype,&
						itmpNumP,CSYN
	implicit none

	integer,intent(in)::i

	if(CSYN .and. FinalScrewN<3)then ! this segment will be implemented.
		
		FinalScrewN=FinalScrewN+1
		FinalScrewSeg(FinalScrewN,1)=SegScrew(i,1)
		FinalScrewSeg(FinalScrewN,2)=SegScrew(i,2)
		FinalScrewU(FinalScrewN,1)=UScrew(i,1)
		FinalScrewU(FinalScrewN,2)=UScrew(i,2)
		FinalScrewMiller(FinalScrewN)=iCSMiller(1)
	endif

end subroutine ScrewRearrangement
!******************************************************************
!   ImplementCrossSlip:                                           
!	After implementation, the number of neighbor is set to zero.  
!	
!	Updated:ZQ Wang
!	4/15/06:
!	After determineInter, use tvTmpPL instead of DPoints.
!******************************************************************
subroutine ImplementCrossSlip
	use ZQMPI
	use vectors
	use variables,only:GPlanes,DPoints
	use CrossSlipMD,only:ScrewN,FinalScrewN,FinalScrewMiller,FinalScrewSeg,&
						FinalScrewU,itmpBurgers,itmpIndex,itmpLooptype,&
						lToImplement,iCSCount
	implicit none

	integer::I,J,K

	integer::i_p1,i_p2,i_m1,i_m2
	type(vector),dimension(ScrewN*3)::tvTmpInterP
	integer::iNIndexCp,iNInterP

	!//////////////////////////////////////////////////////////////
	iNInterP=0 ! number of intersection points
	i_p1=DPoints(itmpIndex(1))%iPlane
	i_m1=GPlanes(i_p1)%iMiller
	!/////////////////////////////////////////////////////////////
	! If no segment to crossslip, return
	! If there is at least one segment to cross-slip, cross-slip counter adds one.
	if(.not.lToImplement .or. FinalScrewN==0)return
	iCSCount=iCSCount+FinalScrewN

	do i=1,FinalScrewN
		! determine the intersection point and glide planes, millers
		call DetermineInter(i,i_p1,i_m1,i_p2,i_m2,iNInterP,tvTmpInterP)
		!////////////////////////////////////////////////////////////
		! Non-cross-slip part
		if(itmpLooptype==1 .or. (itmpLooptype==0 .and. i>1))then
			call GenerateNonCSSeg(i,i_p1,i_m1,i_p2,i_m2,iNInterP,tvTmpInterP)
		endif
		!///////////////////////////////////////////////////////////
		! generate cross-slip part.
		call GenerateCSSeg(i,i_p1,i_m1,i_p2,i_m2,iNInterP,tvTmpInterP)
		!///////////////////////////////////////////////////////////
		if(i==FinalScrewN)then
			j=i+1
			call GenerateNonCSSeg(j,i_p1,i_m1,i_p2,i_m2,iNInterP,tvTmpInterP)
		endif
	enddo

		
end subroutine ImplementCrossSlip
!******************************************************************
!	DetermineInter:
!		Determine the cross-slip plane(miller,origin), determine the intersection
!	point.
!
!	Updated: ZQ Wang
!	4/4/06: ok
!	4/8/06: made modification
!******************************************************************
subroutine DetermineInter(iCurscrewN,i_p1,i_m1,i_p2,i_m2,iNInter,tvTmpInterP)
	use ZQMPI
	use vectors
	use variables,only:iNPlane,GPlanes,Miller,A_CUBE,MAX_Plane, &
					Dpoints,iloop,iMatType,iBCC,iFCC
	use CrossSlipMD,only:FinalScrewN,FinalScrewSeg,FinalScrewU,ScrewN,&
						FinalScrewMiller,tvTmpPL,itmpNumP,  &
						itmpIndex
	implicit none

	integer::i,j,k
	integer::i1,i2,i3,i4
	integer,intent(in)::iCurscrewN,i_p1,i_m1
	integer::i_p2,i_m2
	integer::iNInter
	type(vector),dimension(ScrewN*3)::tvTmpInterP
	type(vector)::tmpOrigin,tmpMiller,tmpP1,tmpP2,tmpP3,tmpP4,tmpP5
	type(vector),dimension(2)::tmpP,tmpT
	double precision::shap(4,3),u

	logical::lFlag

	!....tmp origin....	
	u=FinalScrewU(iCurScrewN,1)
	i1=FinalScrewSeg(iCurScrewN,1)
	tmpP(1)=DPoints(itmpIndex(i1))%tvPG
	tmpT(1)=DPoints(itmpIndex(i1))%tvTG
	call TransNextGlobal(itmpIndex(i1),tmpP(2),tmpT(2))

	call getshape(shap,u)
	tmpP1=shap(1,1)*tmpP(1)+shap(2,1)*tmpT(1)+ &
			shap(3,1)*tmpP(2)+shap(4,1)*tmpT(2)

	!...............
	u=FinalScrewU(iCurScrewN,2)
	i1=FinalScrewSeg(iCurScrewN,2)
	tmpP(1)=DPoints(itmpIndex(i1))%tvPG
        tmpT(1)=DPoints(itmpIndex(i1))%tvTG
	call TransNextGlobal(itmpIndex(i1),tmpP(2),tmpT(2))

	call getshape(shap,u)
	tmpP2=shap(1,1)*tmpP(1)+shap(2,1)*tmpT(1)+ &
                        shap(3,1)*tmpP(2)+shap(4,1)*tmpT(2)

	!...............
	tmpP3=0.5d0*(tmpP1+tmpP2) !tmpOrigin
	
	!...find miller index.....
	i_m2=FinalScrewMiller(iCurScrewN)

	!....begin to compute origin......
	!....Map tmpP5 to the local system, set pl%v(3)=0, recompute pg....
	!....is just the origin.....
	if(iMatType==iFCC)then
		if(i_m2==1)then
			tmpP4%v(1)=0.d0;tmpP4%v(2)=0.d0;tmpP4%v(3)=0.d0
		elseif(i_m2==2)then
			tmpP4%v(1)=0.d0;tmpP4%v(2)=A_Cube;tmpP4%v(3)=0.d0
		elseif(i_m2==3)then
			tmpP4%v(1)=A_Cube;tmpP4%v(2)=A_Cube;tmpP4%v(3)=0.d0
		else
			tmpP4%v(1)=A_Cube;tmpP4%v(2)=0.d0;tmpP4%v(3)=0.d0
		endif
	elseif(iMatType==iBCC)then
		if(i_m2==1)then
			tmpP4%v(1)=0.0;tmpP4%v(2)=0.0;tmpP4%v(3)=0.0
		elseif(i_m2==2)then
			tmpP4%v(1)=0.0;tmpP4%v(2)=A_CUBE;tmpP4%v(3)=0.0
		elseif(i_m2==3)then
			tmpP4%v(1)=0.0;tmpP4%v(2)=A_Cube;tmpP4%v(3)=0.0
		elseif(i_m2==4)then
			tmpP4%v(1)=A_Cube;tmpP4%v(2)=A_Cube;tmpP4%v(3)=0.0
		elseif(i_m2==5)then
			tmpP4%v(1)=A_Cube;tmpP4%v(2)=A_Cube;tmpP4%v(3)=0.0
		elseif(i_m2==6)then
			tmpP4%v(1)=A_Cube;tmpP4%v(2)=0.0;tmpP4%v(3)=0.0
		endif
	endif

	tmpP5=Miller(i_m2)%ES*(tmpP4-tmpP3)
	tmpP5%v(3)=0.d0
	tmpOrigin=TRANS(Miller(i_m2)%ES)*tmpP5+tmpP3

	!....determine plane id....
	lFlag=.true.
	do i=1,iNPlane
		if(GPlanes(i)%iMiller==i_m2 .and. MAG(tmpOrigin-GPlanes(i)%Origin)==0)then
			i_p2=i
			lFlag=.false.
			exit
		endif
	enddo

	if(lFlag)then ! add new plane
		if(iNPlane<MAX_Plane)then
			iNPlane=iNplane+1
			GPlanes(iNPlane)%Origin=tmpOrigin
			GPlanes(iNPLane)%iMiller=i_m2
			i_p2=iNPlane
		else
			print *,"Plane number exceeds maximum number...."
			stop
		endif
	endif

	!.....FInd inter from tmpP1,tmpP2.....
	do i=1,10
		! map to plane 1
		tmpP1=Miller(i_m1)%ES*(tmpP1-GPlanes(i_p1)%Origin)
		tmpP2=Miller(i_m1)%ES*(tmpP2-GPlanes(i_p1)%Origin)
		tmpP1%v(3)=0.d0
		tmpP2%v(3)=0.d0
		tmpP1=Trans(Miller(i_m1)%ES)*tmpP1+GPlanes(i_p1)%Origin
		tmpP2=Trans(Miller(i_m1)%ES)*tmpP2+GPlanes(i_p1)%Origin
		! map to plane 2
		tmpP1=Miller(i_m2)%ES*(tmpP1-GPlanes(i_p2)%Origin)
		tmpP2=Miller(i_m2)%ES*(tmpP2-GPlanes(i_p2)%Origin)
		tmpP1%v(3)=0.d0
		tmpP2%v(3)=0.d0
		tmpP1=Trans(Miller(i_m2)%ES)*tmpP1+GPlanes(i_p2)%Origin
		tmpP2=Trans(Miller(i_m2)%ES)*tmpP2+GPlanes(i_p2)%Origin
	enddo

	iNInter=iNInter+2
	tvTmpInterP(iNInter-1)=tmpP1
	tvTmpInterP(iNInter)=tmpP2
	
end subroutine DetermineInter
!******************************************************************
!	GenerateCSSeg:
!		Interpoints and plane number, miller number have been
!	generated.
!
!	Be careful about u:
!		1). In the beginning, if u==0.d0, 
!			the point should be included in;
!		2). In the end, if u==1.d0, 
!			the next point should be included in.
!	This is reverse for non-screw segment. 
!	
!	Old index for CS part are within the bound of the end points of CS part.
!
!	Updated: ZQ Wang
!	4/4/06: Position between local and global should be differrentiated
!******************************************************************
subroutine GenerateCSSeg(iCurScrewSeg,i_p1,i_m1,i_p2,i_m2,iNInterP,tvTmpInterP)
	use ZQMPI
	use vectors
	use variables,only:NPOint_I,DPoints,iNPoint,MAX_NODE,Miller,GPlanes,&
					iDPtypeFree,iDPTypeFixed,zero,dpEMass0
	use CrossSlipMD,only:FinalScrewSeg,FinalScrewU,itmpNumP,ScrewN, &
						itmpIndex,FinalScrewN,itmpLooptype,itmpBurgers,dpCSLength
	implicit none

	integer,intent(in)::iCurScrewSeg,i_p1,i_m1,i_p2,i_m2,iNInterP
	type(vector),intent(in),dimension(ScrewN*3)::tvTmpInterP

	integer::i,j,k
	integer::iNumP,iNumNewP
	logical::lFlag
	type(vector),allocatable,dimension(:)::PL,TL

	!
	integer,dimension(MAX_NODE*10)::iLocalIndex
	integer::iLIbegin,iLIEnd
	!....check how many points are on the screw segment....
	lFlag=.true.
	iNumP=0

	if(FinalScrewSeg(iCurScrewSeg,1)>FinalScrewSeg(iCurScrewSeg,2))then
		iNumP=itmpNumP-FinalScrewSeg(iCurScrewSeg,1)+ &
				FinalScrewSeg(iCurScrewSeg,2)
	else
		iNumP=FinalScrewSeg(iCurScrewSeg,2)-FinalScrewSeg(iCurScrewSeg,1)
	endif
	iLIBegin=mod(FinalScrewSeg(iCurScrewSeg,1)+1,itmpnumP)   !index on old loop, these point should be kept when
								! generate new loops
	if(iLIBegin==0)iLIBegin=itmpNumP
	iLIEnd=FinalScrewSeg(iCurScrewSeg,2)       !

	if(FinalScrewU(iCurScrewSeg,1)==0.d0)then
		iNumP=iNumP+1
		iLIBegin=FinalScrewSeg(iCurScrewSeg,1)
	endif
	if(FinalScrewU(iCurScrewSeg,2)==1.d0)then
		iNumP=iNumP+1
		iLIEnd=mod(FinalScrewSeg(iCurScrewSeg,2)+1,itmpNumP)
		if(iLIEnd==0)iLIEnd=itmpNUmP
	endif
	
!	if(iNumP>Max_Node)iNumP=Max_Node
	if(iNumP<NPOint_I)then
		iNumNewP=NPOINT_I
		lFlag=.false.
	else
		iNumNewP=iNumP
	endif

	!initialize PL,TL
	allocate(PL(iNumNewP),TL(iNumNewP))
	PL(1)=Miller(i_m2)%ES*(tvTmpInterP(2*iCurScrewSeg-1)-GPlanes(i_p2)%Origin)
	PL(2)=Miller(i_m2)%ES*(tvTmpInterP(2*iCurScrewSeg)-GPlanes(i_p2)%Origin)
	dpCSLength=dpCSLength+MAG(PL(1)-PL(2))

	call setNewInitDis(iNumNewP,PL,TL)

	!/////////////////////////////////////////////////////////
	!assign value to DPOints

	!---------------------------------------------
	!index of DPs
	! Put new points at the end of the index array
	do i=1,iNumNewP
		if(i<=iNumP)then
			j=mod(iLIBegin+i-1,itmpNumP)
			if(j==0)j=itmpNumP
			iLocalIndex(i)=itmpIndex(j)
		else
			iLocalIndex(i)=iNPoint+i-iNumP
		endif
	enddo
	iNPoint=iNPoint+iNUmNewP-iNumP
	!------------------------------------------------

	!assign values to DPs
	do i=1,iNumNewP
		j=iLocalIndex(i)

		DPoints(j)%iBurgers=itmpBurgers
		DPoints(j)%iNumNei=0
		DPoints(j)%tvAcce=zero%v(1)
		DPoints(j)%tvPreV=zero%v(1)
		DPoints(j)%tvPreVT=zero%v(1)
		DPoints(j)%dpMass=dpEMass0

		!lstat
		DPoints(j)%lStat=.true.

		!....Point type....
		DPoints(j)%iType=iDPTypeFree

		!fixed first or last point 
		if(i==1.or.i==iNumNewP)then
			DPoints(j)%iType=iDPtypeFixed
		endif

		!....Position....
			DPoints(j)%tvPL=PL(i)
			DPoints(j)%tvTL=TL(i)
			DPoints(j)%iPlane=i_p2
		!....Connection....
		if(i>1)then ! beginning, no worry about if it is inter.
			DPoints(j)%iBeginP=iLocalIndex(i-1)
		else
			DPoints(j)%iBeginP=-1
		endif

		if(i<iNumNewP)then ! ending, no worry about if it is inter.
			DPoints(j)%iEndP=iLocalIndex(i+1)
		else
			DPoints(j)%iEndP=-1
		endif
		DPoints(j)%iNumNei=0
	enddo

	deallocate(PL,TL)

end subroutine GenerateCSSeg
!******************************************************************
!	GenerateNonCSSeg:
!		1).Generate non-screw cross-slip dislocations.
!		2).Also, connect to screw dislocations.
!
!	For close loop: called after first screwsegment. If last screw,
!	should consider the last non-screw.
!
!	For open, consider last screw specially too.
!    
!	iNCSBegin: old index that the second point on the non-screw part, the 
!		  first point is the last of the screw or the first point on the loop.
!		  It is used to assign index
!		  for iNumP points should be re-filled when calculating new 
!		  positions.
!	iNumP: old index of points on the non-screw part, does not include the
!	       boundary points, the beginning or end, so that is why iNumAddPoints
!		is set
!	iNumAddPoints: to have the beginning or end points for the non-screw parts.
!		       These points are shared between screw and non-screw parts, they
!			are not counted in iNumP
!
!	Update: ZQ Wang
!	4/8/06:
!	Working on this.
!******************************************************************
subroutine GenerateNonCSSeg(iCurScrewSeg,i_p1,i_m1,i_p2,i_m2,iNInterP,tvTmpInterP)
	use ZQMPI
	use vectors
	use variables,only:MAX_NODE,NPOINT_I,iNPoint,DPoints,iDPTypeFree,&
						iDPTypeFixed,MIller,GPlanes,dpEMass0,zero
	use CrossSlipMD,only:ScrewN,itmpLooptype,FinalScrewN,FinalScrewSeg,&
				FinalScrewU,itmpIndex,itmpNumP, &
				itmpBUrgers
	implicit none

	integer,intent(in)::iCurScrewSeg,i_p1,i_m1,i_p2,i_m2,iNInterP
	type(vector),intent(in),dimension(ScrewN*3)::tvTmpInterP

	integer::i,j,k
	integer::iNumP,iNumNewP,iNumAddPoint   !iNumP: number of old index of nodes that belongs to a non-screw part
	logical::lFlag
	type(vector),allocatable,dimension(:)::PL,TL

	!
	integer,dimension(MAX_NODE*10)::iLocalIndex
	integer::iLIbegin,iLIEnd   !segments from which the non-screw parts start and end
	integer::iNCSBegin   !beginning local old index for non-screw parts.
	integer::iCSBegin
	double precision::dpLIBeginU,dpLIEndU

	iNumP=0
	iNumNewP=0
	!....find the beginning index for current screw segment.
	if(iCurScrewSeg<=FinalScrewN)then
		iCSBegin=mod(FinalScrewSeg(iCurScrewSeg,1)+1,itmpNumP)
		if(iCSBegin==0)iCSBegin=itmpNumP
		if(FinalScrewU(iCurScrewSeg,1)==0.d0)then
			iCSBegin=FinalScrewSeg(iCurScrewSeg,1)
		endif
	elseif(itmpLooptype==1)then !open loop, last one
		iCSBegin=itmpNumP
	else !close loop, point to first seg.
		iCSBegin=mod(FinalScrewSeg(1,1)+1,itmpNumP)
		if(iCSBegin==0)iCSBegin=itmpNumP
		if(FinalScrewU(1,1)==0.d0)then
			iCSBegin=FinalScrewSeg(1,1)
		endif
	endif
		

	!....if open loop first seg is screwsegment, then exit,....
	!....this case no before-connected non-screw segment.......
	if(itmpLooptype==1)then
		if(iCurScrewSeg==1 .and. &
			FinalScrewSeg(iCurScrewSeg,1)==1 .and. &
			FinalScrewU(iCurScrewSeg,1)==0.d0)then
			return
		endif
		if(iCurScrewSeg==FinalScrewN+1 .and.  &
			FinalScrewSeg(iCurScrewSeg-1,2)==itmpNumP-1 .and. &
			FinalScrewU(iCurScrewSeg-1,2)==1.d0)then
			return
		endif
	endif

	!....determine how many points are on this segment....
	!....open loop begin and end segments are different......
	!....others are same.....................................
	!....should add some special consideration for ending....

	if(itmpLooptype==1 .and. iCurScrewSeg==1)then  !open loop the first edge may be possibly from the first point.
		iNumP=FinalScrewSeg(iCurScrewSeg,1)
		iLIBegin=1
		dpLIBeginU=0.d0
		iLIEnd=FinalScrewSeg(iCurScrewSeg,1)
		dpLIEndU=FinalScrewU(iCurScrewSeg,1)
		iNCSBegin=1
		if(finalScrewU(iCurScrewSeg,1)==0.d0)then
			iNumP=FinalScrewSeg(iCurScrewSeg,1)-1
			iLIEnd=FinalScrewSeg(iCurScrewSeg,1)-1
			dpLIEndU=1.d0
		endif
	elseif(itmpLooptype==1 .and. iCurScrewSeg==FinalScrewN+1)then   !open loop last non-screw possibly ends at last point.
		iNumP=itmpNumP-FinalScrewSeg(iCurScrewSeg-1,2)
		iLIBegin=FinalScrewSeg(iCurScrewSeg-1,2)
		dpLIBeginU=FinalScrewU(iCurScrewSeg-1,2)
		iLIEnd=iTmpNumP-1
		dpLIEndU=1.0
		iNCSBegin=FinalScrewSeg(iCurScrewSeg-1,2)+1
		if(finalScrewU(iCurScrewSeg-1,2)==1.d0)then
			iNumP=iNumP-1
			iLIBegin=iLIBegin+1
			dpLIBeginU=0.d0
			iNCSBegin=FinalScrewSeg(iCurScrewSeg-1,2)+2
		endif
	elseif(itmplooptype==1)then
		iNumP=FinalScrewSeg(iCurScrewSeg,1)-FinalScrewSeg(iCurScrewSeg-1,2)
		iLIBegin=FinalScrewSeg(iCurScrewSeg-1,2)
		dpLIBeginU=FinalScrewU(iCurScrewSeg-1,2)
		iLIEnd=FinalScrewSeg(iCurScrewSeg,1)
		dpLIEndU=FinalScrewU(iCurScrewSeg,1)
		iNCSBegin=FinalScrewSeg(iCurScrewSeg-1,2)+1
		if(finalScrewU(iCurScrewSeg-1,2)==1.d0)then
			iNumP=iNumP-1
			iLIBegin=iLIBegin+1
			dpLIBeginU=0.d0
			iNCSBegin=FinalScrewSeg(iCurScrewSeg-1,2)+2
		endif
		if(finalScrewU(iCurScrewSeg,1)==0.d0)then
			iNumP=iNUmP-1
			iLIEnd=FinalScrewSeg(iCurScrewSeg,1)-1
			dpLIEndU=1.d0
		endif
	elseif(itmplooptype==0)then
		j=iCurScrewSeg
		k=iCurScrewSeg-1
		if(iCurScrewSeg==1)then
			k=FinalScrewN
		endif

		if(FinalScrewSeg(j,1)<FinalScrewSeg(k,2))then
			iNumP=itmpNumP-FinalScrewSeg(k,2)+FinalScrewSeg(j,1)
		else
			iNumP=FinalScrewSeg(j,1)-FinalScrewSeg(k,2)
		endif

		iLIBegin=FinalScrewSeg(k,2)
		dpLIBeginU=FinalScrewU(k,2)
		iLIEnd=FinalScrewSeg(j,1)
		dpLIEndU=FinalScrewU(j,1)
		iNCSBegin=mod(FinalScrewSeg(k,2)+1,itmpNumP)
		if(iNCSBegin==0)iNCSBegin=itmpNumP
		if(finalScrewU(k,2)==1.d0)then
			iNumP=iNumP-1
			iLIBegin=mod(iLIBegin+1,itmpNumP)
			if(iLIBegin==0)iLIBegin=itmpNumP
			dpLIBeginU=0.d0
			iNCSBegin=mod(FinalScrewSeg(k,2)+2,itmpNumP)
			if(iNCSBegin==0)iNCSBegin=itmpNumP
		endif
		if(finalScrewU(j,1)==0.d0)then
			iNumP=iNUmP-1
			iLIEnd=iLIEnd-1
			if(iLIEnd<=0)iLIEnd=itmpNumP
			dpLIEndU=1.d0
		endif

			
	endif

	iNumNewP=iNumP
			
	if(iNumNewP<NPOINT_I)then
		iNumNewP=NPoint_I
	endif

	do i=1,iNumNewP
		if(i<=iNumP)then
			j=mod(iNCSBegin+i-1,itmpNumP)
			if(j==0)j=itmpNumP
			iLocalIndex(i)=itmpIndex(j)   !iNCSBegin takes of different situtations, here don't worry, just use the 
                                                                  !the same expression.
		else
			iLocalIndex(i)=iNPoint+i-iNumP
		endif
	enddo

	iNPoint=iNPoint+iNumNewP-iNumP

	!....begin to select points....		
	allocate(PL(iNumNewP),TL(iNumNewP))
	call SelectPoints(iNumNewP,iLIBegin,iLIEnd,dpLIBeginU,dpLIEndU,PL,TL)

	call getTangent(iNumNewP,PL,TL,itmpLooptype)
	!....assign value of points....
	do i=1,iNumNewP
		j=iLocalIndex(i)
		
		DPoints(j)%iNumNei=0
		DPoints(j)%iBurgers=itmpBurgers
		DPoints(j)%tvAcce=zero%v(1)
                DPoints(j)%tvPreV=zero%v(1)
		DPoints(j)%tvPreVT=zero%v(1)
                DPoints(j)%dpMass=dpEMass0
		
		!lstat
		DPoints(j)%lStat=.true.

		!point type
		DPOints(j)%iType=iDPTypeFree
		if(i==1.or.i==iNumNewP)then
			DPoints(j)%itype=iDPTypeFixed
		endif
		
		!connection
		if(i<iNumNewP)then
			DPoints(j)%iEndP=iLocalIndex(i+1)
		else
			DPoints(j)%iEndP=-1
		endif
		if(i>1)then
			DPoints(j)%iBeginP=iLocalIndex(i-1)
		else
			DPoints(j)%iBeginP=-1
		endif

		!assign PL,TL value
		DPoints(j)%tvPL=PL(i)
		DPoints(j)%tvTL=TL(i)
		DPoints(j)%iPlane=i_p1
	
	enddo

	deallocate(PL,TL)
end subroutine GenerateNonCSSeg

!*******************************************************************
! SeclectPoints:
!	select points from the non-screw segment
!
!	Updated: ZQ Wang
!	4/9/06:
!	Most of the code is ok.
!	Need more work on get the local coordinates for the loop, not the segment.
!
!*******************************************************************
subroutine SelectPoints(iNumNewP,iLIBegin,iLIEnd,dpLIBeginU,dpLIEndU,PL,TL)
	use ZQMPI
	use vectors
	use variables,only:MAX_NODE,DPOints
	use CrossSlipMD,only:itmpIndex,itmpNumP,tvTmpPL
	implicit none

	integer,intent(in)::iNumNewP,iLIBegin,iLIEnd
	double precision,intent(in)::dpLIBeginU,dpLIEndU
	type(vector),dimension(MAX_NODE)::PL,TL

	type(vector),dimension(Max_Node*30)::tmpV,tmpT   !i1=20, iNumI at most Max_node, so Max_node*30 is enough
	type(vector)::tvPL1,tvTL1,tvPL2,tvTL2
	double precision::bu,eu
	double precision::shap(4,3),u
	integer::i,j,k,i1,i2,ii
	integer::iCur,iNext
	integer::iNumI
	integer::iID1
	type(DisPoint)::tmpDP

	!/////////////////////////////////////////////////////////////////

	!///////////////////////////////////////////
	do i=1,itmpNumP
		tmpDP=DPoints(itmpIndex(i))
		DPoints(itmpIndex(i))=tvTmpPL(i)
		tvTmpPL(i)=tmpDP
	enddo

	!//////////////////////////////////////////////
	i1=20

	if(iLIEnd<iLIBegin)then
		iNumI=itmpNumP-iLIBegin+1+iLIEnd
	else
		iNumI=iLIEnd-iLIBegin+1
	endif
	!....first point.....
	j=1
	u=dpLIBeginU
	iCur=itmpIndex(iLiBegin)
	iNext=DPoints(iCur)%iEndP

	tvPL1=DPoints(iCur)%tvPL;tvTL1=DPoints(iCur)%tvTL
!	call TransNextLocal(iCur,tvPL2,tvTL2)
	tvPL2=DPoints(iNext)%tvPL;tvTL2=DPoints(iNext)%tvTL

	call getshape(shap,u)
	tmpV(j)=shap(1,1)*tvPL1+shap(2,1)*tvTL1+shap(3,1)*tvPL2+shap(4,1)*tvTL2
	tmpT(j)=shap(1,2)*tvPL1+shap(2,2)*tvTL1+shap(3,2)*tvPL2+shap(4,2)*tvTL2

	do i=1,iNumI
		iID1=iLIBegin+i-1
		if(iID1>itmpNumP)then
			iID1=mod(iID1,itmpNumP)
		endif

		bu=0.d0
		eu=1.d0
		if(iID1==iLIBegin)then
			bu=dpLIBeginU
		endif
		if(iID1==iLIEnd)then
			eu=dpLIEndU
		endif

		iCur=itmpIndex(iID1)
		iNext=DPoints(iCur)%iEndP
		tvPL1=DPoints(iCur)%tvPL;tvTL1=DPoints(iCur)%tvTL
!		call TransNextLocal(iCur,TvPL2,TvTL2)
		tvPL2=DPoints(iNext)%tvPL;tvTL2=DPoints(iNext)%tvTL

		do ii=1,i1
			j=j+1
			u=bu+(eu-bu)/(1.0*i1)*ii
			call getshape(shap,u)
			tmpV(j)=shap(1,1)*tvPL1+shap(2,1)*tvTL1+ &
					shap(3,1)*tvPL2+shap(4,1)*tvTL2
			tmpT(j)=shap(1,2)*tvPL1+shap(2,2)*tvTL1+ &
					shap(3,2)*tvPL2+shap(4,2)*tvTL2
		enddo
	enddo

	k=(j-1)/(iNumNewP-1)
	PL(1)=tmpV(1)
	TL(1)=tmpT(1)
	PL(iNumNewP)=tmpV(j)
	TL(iNumNewP)=tmpT(j)
	i2=0
	i1=1
	do i=2,j-1
		i2=i2+1
		if(mod(i2,k)==0 .and. i1<=iNumNewP-2)then
			i1=i1+1
			PL(i1)=tmpV(i)
			TL(i1)=tmpT(i)
		endif
	enddo

 	!Exchange points back
	do i=1,itmpNumP
                tmpDP=DPoints(itmpIndex(i))
                DPoints(itmpIndex(i))=tvTmpPL(i)
                tvTmpPL(i)=tmpDP
        enddo

end subroutine SelectPoints
!*******************************************************************

!1) lStat should be carefully considered. check first beginning point 
!	and last ending point.
