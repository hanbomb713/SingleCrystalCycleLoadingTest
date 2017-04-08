!***************************************************************************
! This file contains functions with return value as a defined special types.
!***************************************************************************
module FunctionMD
	use ZQMPI
	use vectors
	implicit none

	contains
!*********************************************************
! DetermineConnVec:
!	To determine the connection vector.
!	Input: two point ids
!	Output: the connection vector
!
!	The loop-global position of point 2 is computed as 
!	lgp=gp-cv
!
!	Input:
!	id1: the point whose loop-global position are going to
!		computed.
!	id2: the point which is the reference point with respect
!		the point of id1 is going to be mapped.
!*********************************************************
type(vector) function DeterConnVec(id1,id2)
	use ZQMPI
	use vectors
	use variables,only:zero,DPoints,A_Cube
	implicit none

	integer,intent(in)::id1,id2
	integer::i,j,k
	double precision::dpTmp,dpTmp1

	dpTmp=0.5*A_Cube

	DeterConnVec=Zero%v(1)

	do i=1,3
		dpTmp1=Dpoints(id1)%tvPG%v(i)-DPoints(id2)%tvPG%v(i)
		if(dabs(dpTmp1)>dpTmp)then
			DeterConnVec%v(i)=A_cube*dpTmp1/dabs(dpTmp1)
		endif
	enddo

end function DeterConnVec
!*************************************************************
!	DeterConnVecQ
!*************************************************************
type(vector) function DeterConnVecQ(q1,q2)
	use ZQMPI
	use vectors
	use variables,only:zero,A_Cube
	implicit none

	type(vector),intent(in)::q1,q2
	integer::i,j,k
	double precision::dpTmp,dpTmp1

	dpTmp=0.5*A_Cube

	DeterConnVecQ=Zero%v(1)

	do i=1,3
		dpTmp1=q1%v(i)-q2%v(i)
		if(dabs(dpTmp1)>dpTmp)then
			DeterConnVecQ%v(i)=A_cube*dpTmp1/dabs(dpTmp1)
		endif
	enddo

end function DeterConnVecQ
!****************************************************
! PointPBC:
!	Apply periodic boundary condition for the point
!	input. If the point is inside the box, it is ok. 
!	otherwise, map it back to the simulation box.
!****************************************************
type(vector) function PointPBC(Q)
	use ZQMPI
	use vectors
	use variables,only:A_Cube
	implicit none

	type(vector),intent(in)::Q
	integer::i
	
	PointPBC=Q
	do i=1,3
		if(Q%v(i)<0.d0)then
			PointPBC%v(i)=dmod(Q%v(i),A_Cube)+A_Cube
		elseif(Q%v(i)>A_cube)then
			PointPBC%v(i)=dmod(Q%v(i),A_Cube)
		endif
	enddo

end function PointPBC
!**************************************************************
!dpFCyclicForce:
!   To provide the force in the cyclic loading.
!   Here the strategy is that a constant strain rate is applied,
!   but the direction of the loading is changing back and forth 
!   after a while.
!   That should be translated to a force.
!
!   Total strain is known. Once it reaches a value, the loading is
!   reversed.
!***************************************************************
double precision function dpFCyclicForce()
	use ZQMPI
	use vectors
	use variables,only:totalStress,ElasticCons,Dtime,N_Times,  &
				StrainRate,strainRate_p,MU,dpCycStrain,totalStrain
	implicit none

	!Update stress first based on previous strain rate
	dpFCyclicForce = totalStress     &
                                                +ElasticCons*Dtime*N_TIMES*(StrainRate-StrainRate_p)/Mu
	!Check if the strain loading needs to be reversed
	if(dabs(totalStrain)>=dpCycStrain)then
		StrainRate=-StrainRate
	endif

end function dpFCyclicForce
!******************************************************************
!******************************************************************
logical function funcIsScrew(PointID,T_Q)
	use ZQMPI
	use vectors
	use variables,only:BURGERS,DPOINTS,PI
	implicit none

	type(vector)::T_Q
	integer::PointID

	funcIsScrew=.true.

	if(dabs(UV(T_Q)*UV(Burgers(DPoints(PointID)%iBurgers)))<dcos(45.0/180.0*PI))funcIsScrew=.false. !11/29/08

end function funcIsScrew
!******************************************************************
double precision function funcBCCCSL(length0)
	use ZQMPI
	use vectors
	use variables,only:MU
	implicit none

	double precision::length0


	if(length0==500.0)then
		funcBCCCSL=0.001*MU*1d-6  !in MPa
	elseif(length0==50.0)then
		funcBCCCSL=0.01*MU*1d-6   !in MPa
	endif

end function funcBCCCSL
!******************************************************************
!
!To see if force on the anti-twinning direction: 
!Here lists twinning direction
!******************************************************************
logical function funcIsAntitwin(PointID,tmpF)
	use ZQMPI
	use vectors
	use variables,only:DPoints,Miller, GPlanes, Burgers
	implicit none

	integer,intent(in)::PointID
	type(vector),intent(in)::tmpF
	integer::tmpMiller,tmpBurgers
	type(vector)::tmpLocalF
	type(vector)::tmpLocalB   !local burgers vector
	double precision::VectorAngle,dpAngle

	funcIsAntitwin=.false.

	!Local Force
	tmpMiller=GPlanes(DPoints(PointID)%iPlane)%iMiller
	tmpLocalF=tmpF
	tmpLocalF%v(3)=0.0

	! share same twinning direction: burgers 1
	if(tmpMiller==17 .or. tmpMiller==14 .or. tmpMiller==8)then
		tmpBurgers=1
	elseif(tmpMiller==15 .or. tmpMiller==12 .or. tmpMiller==10)then
		tmpBurgers=7
	elseif(tmpMiller==18 .or. tmpMiller==11 .or. tmpMiller==9)then
		tmpBurgers=4
	elseif(tmpMiller==16 .or. tmpMiller==13 .or. tmpMiller==7)then
		tmpBurgers=6
	endif

	!Local Twinning direction
	tmpLocalB=Miller(tmpMiller)%ES*Burgers(tmpBurgers)

	tmpLocalB=UV(tmpLocalB)

	!check angle
	dpAngle=VectorAngle(tmpLocalF,tmpLocalB)

	!
	if(dpAngle>90.0)then
		funcIsAntitwin=.true.
	else
		funcIsAntitwin=.false.
	endif
	
end function funcIsAntitwin
!******************************************************************
end module
