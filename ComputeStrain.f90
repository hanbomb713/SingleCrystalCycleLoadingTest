!********************************************************************
!	This file contains subroutines that are used to calculated the 
!	plastic strain caused by dislocation slip.
!********************************************************************
!********************************************************************
!	ComputeStrain:
!		Do calculation for plastic strain.
!-------------------------------------------------------------------
!Modified: ZQ Wang, 2/6/06
!	Modified the strain calculation with inertial effect.
!
!********************************************************************
subroutine ComputeStrain
	use ZQMPI
	use vectors
	use variables,only:iNPoint,DPoints
	implicit none

	integer::i,iEnd,IOUT
	logical::IsCloseLoopPoint,lSign

	do i=1,iNPoint
		iEnd=DPoints(i)%iEndP

		if(iEnd/=-1 .and. DPoints(i)%iPlane==DPoints(iEnd)%iPlane)then   !1/18/04
			lSign=IsCloseLoopPoint(i)
			if(.not. lSign)then  !this criteria eliminates the close-loop strain because they are fixed.
				call ComputeStrainL(i)
			endif
		endif
	enddo

end subroutine ComputeStrain
!********************************************************************************
! ComputeStrainL:
!	Compute the strain of a point (leading a segment).
!********************************************************************************
subroutine ComputeStrainL(PointID)
	use ZQMPI
	use vectors
	use variables,only:iNPoint,DPoints,GPlanes,Miller,Burgers,DTime,MAX_QUAD,POS, &
						wt,A_cube,strainField,zero
	use DynamicsMD,only:DMatrix
	implicit none

	integer,intent(in)::PointID
	
	type(vector)::PG1,TG1,PG2,TG2,tvUnitVelocity,tvPV1,tvPV2,tvTV1,tvTv2,tvPV3,tvTv3,RU,UnitN
	double precision::dpVelocity,DS,shap(4,3),dpVolume,u
	integer::i,j,k,i_m,i_p,i_b,iNext
	type(matrix)::dpTmpS1
	
	iNext=DPoints(PointID)%iEndP
	i_p=DPoints(PointID)%iPlane
	i_b=DPoints(PointID)%iBurgers
	i_m=GPlanes(i_p)%iMiller
	PG1=DPoints(PointID)%tvPG
	TG1=DPoints(PointID)%tvTG
	call TransNextGlobal(PointID,PG2,TG2)
	dpVolume=A_CUBE**3

	!/////////////////////////////////////////////////////////
	! Average velocity during the time step of ending points.
       tvPV1=TRANS(Miller(i_m)%ES)*(DPoints(PointID)%tvPreV-0.5d0*DTIME*DPoints(PointID)%tvAcce)
       tvTV1=TRANS(Miller(i_m)%ES)*DPoints(PointID)%tvPreVT
       tvPV2=TRANS(Miller(i_m)%ES)*(DPoints(iNext)%tvPreV-0.5d0*DTime*DPoints(iNext)%tvAcce)
       tvTV2=TRANS(Miller(i_m)%ES)*DPoints(iNext)%tvPreVT

	!///////////////////////////////////////////////////////
	! Integration of strain over the segment
	do i=1,MAX_QUAD
		u=(Pos(i)+1)*0.5
		call getshape(shap,u)
		!//////////////////////////////////////////
		! DS
		RU=shap(1,2)*PG1+shap(2,2)*TG1+shap(3,2)*PG2+shap(4,2)*TG2
		DS=DSQRT(RU*RU)

		!//////////////////////////////////////
                ! Velocity and unit vector of the velocity
                tvPV3=shap(1,1)*tvPV1+shap(2,1)*tvTV1+shap(3,1)*tvPV2+shap(4,1)*tvTV2
                dpVelocity=MAG(tvPV3)
		if(dpVelocity==0)then
			dpTmpS1=zero
			goto 1
		endif
		!//////////////////////////////
		! Unit N
		tvUnitVelocity=UV(tvPV3)
		UnitN=UV(tvUnitVelocity**RU)

		!//////////////////////////////
		! Strain
		dpTmpS1=(-0.25/dpVolume*dpVelocity*DS*DTime*wt(i))*(UnitN/Burgers(i_b)+Burgers(i_b)/UnitN)

1		continue
		Strainfield=StrainField+dpTmpS1

	enddo
end subroutine ComputeStrainL
