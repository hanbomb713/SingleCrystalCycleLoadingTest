!**********************************************************************
! Vectors:
!	Definition of new types of variables, new calculation symbols and rules
!	for new defined types.
!
!************************************************************************
MODULE vectors 
	use ZQMPI
	implicit none

	!....type vector....
	TYPE vector
		DOUBLE PRECISION,DIMENSION(3)::v
	END TYPE vector

	!....type matrix.....
	TYPE matrix
		TYPE(vector),DIMENSION(3)::v
	END TYPE matrix

	!....type tensor....
	TYPE tensor
		TYPE(matrix),DIMENSION(3)::v   
	END TYPE tensor
	
	!...dislocation points....
	type DisPoint
		integer::iID             !global index
		integer::iLoopID		 !loop id
		integer::itype
		
		type(vector)::tvPL
		type(vector)::tvTL
		type(vector)::tvPG
		type(vector)::tvTG
		type(vector)::tvPreV         !Previous step velocity, for acceleration
		type(vector)::tvPreVT
		type(vector)::tvAcce
		type(vector)::force

		double precision::dpMass

		integer::iPlane
		integer::iBurgers
		integer::iBeginP
		integer::iEndP
		integer::iNumNei
		integer,allocatable,dimension(:)::iaIDOfNei

		logical::lStat
		logical::lToCal
	end type DisPoint
	
	type DPPointer
		type(DisPoint),pointer::p
	end type

	!....type glide plane.....
	type GlidePlane
		integer::iMiller
		type(vector)::Origin
		integer::iNumPEquiv
		integer,dimension(5)::iIDPEquiv
	end type GlidePlane

	!....type miller index....
	type MillerIndex
		type(vector)::V
		type(matrix)::ES
	end type MillerIndex

	!.....interfaces........

	INTERFACE OPERATOR(+)
		MODULE PROCEDURE v1_plus_v2
		MODULE PROCEDURE m1_plus_m2
		MODULE PROCEDURE t1_plus_t2
	END INTERFACE
	INTERFACE OPERATOR(-)
		MODULE PROCEDURE v1_minus_v2
		MODULE PROCEDURE m1_minus_m2
		MODULE PROCEDURE t1_minus_t2
	END INTERFACE
	INTERFACE OPERATOR(*)
		MODULE PROCEDURE v1_dot_v2;MODULE PROCEDURE real_times_v1
		MODULE PROCEDURE m1_mul_v2;MODULE PROCEDURE m1_mul_m2
		MODULE PROCEDURE t1_mul_v2;MODULE PROCEDURE real_times_t1
		MODULE PROCEDURE real_times_m1
	END INTERFACE
	INTERFACE OPERATOR(/)
		MODULE PROCEDURE v1_div_real;MODULE PROCEDURE m1_div_real
		MODULE PROCEDURE t1_div_real
		MODULE PROCEDURE v1_outerproduct_v2
		MODULE PROCEDURE m1_outerproduct_v2
		MODULE PROCEDURE v1_outerproduct_m2
	END INTERFACE
	INTERFACE OPERATOR(//)
		MODULE PROCEDURE m1_outerproduct1_v2
	END INTERFACE
	INTERFACE OPERATOR(**)
		MODULE PROCEDURE v1_cross_v2
		MODULE PROCEDURE theta_rotate_v1
	END INTERFACE
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	CONTAINS
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

	TYPE(vector) FUNCTION v1_minus_v2(a,b)
	TYPE(vector), INTENT(IN)::a,b
	integer::i
	DO i=1,3;v1_minus_v2%v(i)=a%v(i)-b%v(i);END DO
	END FUNCTION v1_minus_v2

	TYPE(matrix) FUNCTION m1_minus_m2(a,b)
	TYPE(matrix), INTENT(IN)::a,b
	integer::i,j
	DO i=1,3;DO j=1,3
	m1_minus_m2%v(i)%v(j)=a%v(i)%v(j)-b%v(i)%v(j)
	END DO;END DO
	END FUNCTION m1_minus_m2

	TYPE(tensor) FUNCTION t1_minus_t2(a,b)
	TYPE(tensor), INTENT(IN)::a,b
	integer::i,j,k
	DO i=1,3;DO j=1,3;DO k=1,3
	t1_minus_t2%v(i)%v(j)%v(k)=a%v(i)%v(j)%v(k)-b%v(i)%v(j)%v(k)
	END DO;END DO;END DO
	END FUNCTION t1_minus_t2

	TYPE(vector) FUNCTION v1_plus_v2(a,b)
	TYPE(vector), INTENT(IN)::a,b
	integer::i
	DO i=1,3;v1_plus_v2%v(i)=a%v(i)+b%v(i);END DO
	END FUNCTION v1_plus_v2

	TYPE(matrix) FUNCTION m1_plus_m2(a,b)
	TYPE(matrix), INTENT(IN)::a,b
	integer::i,j
	DO i=1,3;DO j=1,3
	m1_plus_m2%v(i)%v(j)=a%v(i)%v(j)+b%v(i)%v(j)
	END DO;END DO
	END FUNCTION m1_plus_m2

	TYPE(tensor) FUNCTION t1_plus_t2(a,b)
	TYPE(tensor), INTENT(IN)::a,b
	integer::i,j,k
	DO i=1,3;DO j=1,3;DO k=1,3
	t1_plus_t2%v(i)%v(j)%v(k)=a%v(i)%v(j)%v(k)+b%v(i)%v(j)%v(k)
	END DO;END DO;END DO
	END FUNCTION t1_plus_t2

	TYPE(vector) FUNCTION real_times_v1(a,b)
	DOUBLE PRECISION,INTENT(IN)::a
	TYPE(vector), INTENT(IN)::b
	integer::i
	DO i=1,3;real_times_v1%v(i)=a*b%v(i);END DO
	END FUNCTION real_times_v1

	TYPE(matrix) FUNCTION real_times_m1(a,b)
	DOUBLE PRECISION,INTENT(IN)::a
	TYPE(matrix), INTENT(IN)::b
	integer::i,j
	DO i=1,3;DO j=1,3
	real_times_m1%v(i)%v(j)=a*b%v(i)%v(j)
	END DO;END DO
	END FUNCTION real_times_m1

	TYPE(tensor) FUNCTION real_times_t1(a,b)
	DOUBLE PRECISION,INTENT(IN)::a
	TYPE(tensor), INTENT(IN)::b
	integer::i,j,k
	DO i=1,3;DO j=1,3;DO k=1,3
	real_times_t1%v(i)%v(j)%v(k)=a*b%v(i)%v(j)%v(k)
	END DO;END DO;END DO
	END FUNCTION real_times_t1

	TYPE(vector) FUNCTION v1_div_real(a,b)
	DOUBLE PRECISION,INTENT(IN)::b
	TYPE(vector), INTENT(IN)::a
	integer::i
	DO i=1,3;v1_div_real%v(i)=a%v(i)/b;END DO
	END FUNCTION v1_div_real

	TYPE(matrix) FUNCTION m1_div_real(a,b)
	DOUBLE PRECISION,INTENT(IN)::b
	TYPE(matrix), INTENT(IN)::a
	integer::i,j
	DO i=1,3;DO j=1,3
	m1_div_real%v(i)%v(j)=a%v(i)%v(j)/b
	END DO;END DO
	END FUNCTION m1_div_real

	TYPE(tensor) FUNCTION t1_div_real(a,b)
	DOUBLE PRECISION,INTENT(IN)::b
	TYPE(tensor), INTENT(IN)::a
	integer::i,j,k
	DO i=1,3;DO j=1,3;DO k=1,3
	t1_div_real%v(i)%v(j)%v(k)=a%v(i)%v(j)%v(k)/b
	END DO;END DO;END DO
	END FUNCTION t1_div_real

	TYPE(matrix) FUNCTION v1_outerproduct_v2(a,b)
	TYPE(vector), INTENT(IN)::a,b
	integer::i,j
	DO i=1,3;DO j=1,3
	v1_outerproduct_v2%v(i)%v(j)=a%v(i)*b%v(j)
	END DO;END DO
	END FUNCTION v1_outerproduct_v2

	TYPE(tensor) FUNCTION m1_outerproduct_v2(a,b)
	TYPE(matrix), INTENT(IN)::a
	TYPE(vector), INTENT(IN)::b
	integer::i,j,k
	DO i=1,3;DO j=1,3;DO k=1,3
	m1_outerproduct_v2%v(i)%v(j)%v(k)=a%v(i)%v(j)*b%v(k);
	END DO;END DO;END DO
	END FUNCTION m1_outerproduct_v2

	TYPE(tensor) FUNCTION m1_outerproduct1_v2(a,b)
	TYPE(matrix), INTENT(IN)::a
	TYPE(vector), INTENT(IN)::b
	integer::i,j,k
	DO i=1,3;DO j=1,3;DO k=1,3
	m1_outerproduct1_v2%v(i)%v(j)%v(k)=a%v(k)%v(i)*b%v(j);
	END DO;END DO;END DO
	END FUNCTION m1_outerproduct1_v2

	TYPE(tensor) FUNCTION v1_outerproduct_m2(a,b)
	TYPE(vector), INTENT(IN)::a
	TYPE(matrix), INTENT(IN)::b
	integer::i,j,k
	DO i=1,3;DO j=1,3;DO k=1,3
	v1_outerproduct_m2%v(i)%v(j)%v(k)=b%v(j)%v(k)*a%v(i);
	END DO;END DO;END DO
	END FUNCTION v1_outerproduct_m2

	TYPE(vector) FUNCTION VECT(a,b)
	TYPE(tensor), INTENT(IN)::a
	INTEGER::b,i
	DO i=1,3
	SELECT CASE(b)
	CASE(1)
	VECT%v(i)=a%v(i)%v(1)%v(1)+a%v(i)%v(2)%v(2)+a%v(i)%v(3)%v(3)
	CASE(2)
	VECT%v(i)=a%v(1)%v(i)%v(1)+a%v(2)%v(i)%v(2)+a%v(3)%v(i)%v(3)
	CASE(3)
	VECT%v(i)=a%v(1)%v(1)%v(i)+a%v(2)%v(2)%v(i)+a%v(3)%v(3)%v(i)
	END SELECT
	END DO
	END FUNCTION VECT

	DOUBLE PRECISION FUNCTION v1_dot_v2(a,b)
	TYPE(vector), INTENT(IN)::a,b
	v1_dot_v2=a%v(1)*b%v(1)+a%v(2)*b%v(2)+a%v(3)*b%v(3)
	END FUNCTION v1_dot_v2

	TYPE(vector) FUNCTION v1_cross_v2(a,b)
	TYPE(vector), INTENT(IN)::a,b
	v1_cross_v2%v(1)=a%v(2)*b%v(3)-b%v(2)*a%v(3)
	v1_cross_v2%v(2)=a%v(3)*b%v(1)-b%v(3)*a%v(1)
	v1_cross_v2%v(3)=a%v(1)*b%v(2)-b%v(1)*a%v(2)
	END FUNCTION v1_cross_v2

	TYPE(vector) FUNCTION theta_rotate_v1(a,b)
	DOUBLE PRECISION, INTENT(IN)::a
	TYPE(vector), INTENT(IN)::b
	TYPE(vector)::c
	c%v(1)=0;c%v(2)=0;c%v(3)=1
	theta_rotate_v1=dcos(a)*b+dsin(a)*c**b
	END FUNCTION theta_rotate_v1

	TYPE(vector) FUNCTION m1_mul_v2(a,b)
	TYPE(matrix), INTENT(IN)::a
	TYPE(vector), INTENT(IN)::b
	integer::i
	DO i=1,3;m1_mul_v2%v(i)=a%v(i)*b;END DO
	END FUNCTION m1_mul_v2

	TYPE(matrix) FUNCTION TRANS(a)
	TYPE(matrix),INTENT(IN)::a
	integer::i,j
	DO i=1,3;DO j=1,3;TRANS%v(i)%v(j)=a%v(j)%v(i);END DO;END DO
	END FUNCTION TRANS

	TYPE(matrix) FUNCTION m1_mul_m2(a,b)
	TYPE(matrix), INTENT(IN)::a,b
	integer::i
	DO i=1,3;m1_mul_m2%v(i)=TRANS(b)*a%v(i);END DO
	END FUNCTION m1_mul_m2

	TYPE(matrix) FUNCTION t1_mul_v2(a,b)
	TYPE(tensor), INTENT(IN)::a
	TYPE(vector), INTENT(IN)::b
	integer::i
	DO i=1,3;t1_mul_v2%v(i)=a%v(i)*b;END DO
	END FUNCTION t1_mul_v2

	DOUBLE PRECISION FUNCTION MAG(A)
	TYPE(vector),INTENT(IN)::A
	MAG=DSQRT(A%v(1)**2+A%v(2)**2+A%v(3)**2)
	END FUNCTION MAG

	TYPE(vector) FUNCTION UV(A)
		TYPE(vector),INTENT(IN)::A
		UV=A/MAG(A)
	END FUNCTION UV

	type(vector) function GetUV(a,b,c)

		type(vector):: a,b,tmp1,tmp2

		double precision::tempACube,tmpDP,c
		integer::i

		tempACube=c*0.5d0

		!first, we calculate the traslation vector
		do i=1,3
			tmp1%v(i)=0.0
			tmpDp=a%v(i)-b%v(i)
			if(dabs(tmpDp)>=tempACube)then
				tmp1%v(i)=tmpDp/dabs(tmpDp)*c
			endif
		enddo
		!now, translate the former vector to a new vector

		tmp2=tmp1+b

		!Then, calculate the unit vector
		GetUV=UV(a-tmp2)
	end function GetUV

END MODULE vectors

!********************************************************************
! Variables:
!	Define global veriables used in the code.
!
!*********************************************************************
MODULE VARIABLES 
	use ZQMPI
	USE VECTORS
	IMPLICIT NONE

	!//////////////////////////////////////////////
	! Following are control related variables

	DOUBLE PRECISION::DTIME  !timestep

	!Double precision,parameter::DMaxSpeed=9.85D12   !  3560 m/s, sound speed in copper
	Double precision,parameter::DMaxSpeed=1.1d13   !  3989 m/s, sound speed in copper

	!MAX_NODE: MAXIMUM NUMBER OF NODES FOR A LOOP
	!N_NODE: Number of nodes for loops without obstacles 

	DOUBLE PRECISION::APPLIED_SIG,MU,NU,LATTICE,A_CUBE,MOBILITY, MOB_R,PEIERLS,PEI_R,    &
	               PI,RMIN,U,INTERVAL,INTERVAL2,DELTA_SIG,DENSITY_DIS !zq Aug. 2
	Double precision:: dpLSoundV,dpTSoundV,dpEMass0   !01/24/06, longitudinal sound speed, transverse sound speed,  
								    ! and the effective mass0 

	integer::MAX_QUAD,MAX_PLANE,MAX_NODE,MAX_NEIGHBOR

	INTEGER::N_TIMES,N_TIMESCOPY,ILOOP_TIME,CheckNeiBur,IntegrateMethod,iBeginLoop,iloop
	INTEGER::logAnnihilation, logArrangeLoop,iObstacle,iInertial,iCrossslip   !iInertial=1,high speed, inertial effect should be included.
                                                                                  !icrossslip=1, turn on cross-slip.
	INTEGER::I_P,IL,NPOINT_I,iLoopNum, Load,iArr

	integer,parameter::iDPTypeFree=1
	integer,parameter::iDPTypeFixed=2
	integer,parameter::iDPTypeInter=3

	TYPE(MATRIX)::D,ZERO,SIG_APP,SIG
	!Gauss Qudrature Points
	DOUBLE PRECISION, DIMENSION(300)::POS,WT

	!system control
	! Following are mechanics related variables
	double precision::StrainRate,ElasticCons,StrainRate_p, &
			  DeltaStrain,StrainCurrent

	type(vector)::Loading_Dir                    !loading direction
        type(matrix)::Loading_ES
	double precision::totalStrain,totalStress,totalTime
	TYPE(MATRIX)::StrainField !plastic strain

	double precision::SR1,SR2
	
	!///////////////////////////////////////////
	! Following are dislocation related variables, including slip system for specific materials

	! slip system variables
	integer::iMatType !material tyep, BCC or FCC
	integer,parameter::iBCC=1
	integer,parameter::iFCC=2
	integer,parameter::iHCP=3

	!....burgers vector, miller index, dimension(:)......
	integer::iNumBurgers,iNumMiller    !They are variables depending on iMatType
	TYPE(VECTOR),dimension(:),allocatable::BURGERS
	type(millerIndex),dimension(:),allocatable::Miller

	!----------------------------
	!This is for each simulation
	!....glide planes
	integer::iNPLANE
	type(GlidePlane),allocatable,dimension(:)::GPlanes
	
	!...Dislocation points........
	integer::iNPoint
	type(DisPoint),allocatable,dimension(:),target::DPoints
	
	!///////////////////////////////////////////////
	!	connectVec: defined as vector from first point to new point
	!	only used to get local points
	integer::connectPlane
	type(vector)::connectVec

	!//////////////////////////////////////////////
	double precision::dcrit

	!/////////////////////////////////////////////////
	! Cyclic deformation
	double precision::CycFreq,dpCycF0
	double precision::dpCycStrain

	!/////////////////////////////////////////////////
	!	For velocity statistics  
	!      Modified: 01/28/06
	!	02/26/06
	double precision::dpVSeg(100)    ! 100 speed intervals
	double precision::dpVCountL(100) ! length count, percentage
	integer::iVCountN(100)  ! number of nodes count
	double precision::dpVCountN(100)  ! number of nodes, percentage


	!///////////////////////////////////
	!111908: ATRatio in Dynamics.f90
	double precision::ATRatio,TRatio    !Anti-twin ratio and twin ratio
	double precision::NGRatio

END MODULE VARIABLES
!***************************************************************
!    ANNIHILATIONMD:
!		Used for annihilation subroutines
!***************************************************************
MODULE ANNIHILATIONMD
	use ZQMPI
	implicit none
	
	integer::iAnniCount
	double precision::dpAnniLength
	double precision,parameter::dpAnniDis=100.d0
	double precision,parameter::dpAnniAngle=90.d0

	logical::lLogAnni

end MODULE ANNIHILATIONMD
!****************************************************************
!  LoopReArrangeMD:
!		Used for loopRearrangement subroutines.
!****************************************************************
MODULE LoopRearrangeMD
	use ZQMPI
	use vectors
	implicit none
	
	integer,parameter::iMaxNodeDim=200

	! record variables
	integer::iTmpLooptype,itmpBurgers,itmpPlane
	integer,dimension(iMaxNodeDim)::iTmpIndex ! to record indexes of points in DPoints
	integer::iTmpNumP,iTmpNewNumP ! tmp number of points on dislocations
	type(vector),dimension(iMaxNodeDim)::tvTmpPL ! tmp Pl for dislocation.
	type(vector),dimension(iMaxNodeDim)::tvTmpTL ! tmp Tl for dislocation.
	type(vector),dimension(iMaxNodeDim)::tvTmpAcce  !tmp acceleration for dislocation
	type(vector),dimension(iMaxNodeDim)::tvTmpPrev  !tmp previous velocity
	type(vector),dimension(iMaxNodeDim)::tvTmpPreVT  !tmp previous velocity of tangent
	double precision,dimension(iMaxNodeDim):: dpTmpMass   !tmp mass
	logical::lTmpToCal
	integer::itmpType1,itmpType2

	!control variables
	double precision,parameter::dpMaxAveLength=500.d0
	integer::iNoOfNodes
	double precision::dpAveLength

	type(DisPoint)::tmpDP(2)

end module LoopRearrangeMD
!************************************************
! DynamicsMD:
!	Variables used in solving the EOM.
!************************************************
module DynamicsMD
	use ZQMPI
	use vectors
	implicit none
	
	type MatrixPointer
		double precision,allocatable,dimension(:)::FG
		double precision,allocatable,dimension(:,:)::KG
		double precision,allocatable,dimension(:,:)::QG
		type(vector)::tvPLRate
		type(vector)::tvTLRate
	end type

	type(MatrixPointer),allocatable,dimension(:)::DMatrix

	integer::NEQ

	!local matrix does not change
	double precision,dimension(2,8,8)::KL
	double precision,dimension(2,8)::FL

	integer::iLastLoop

end module DynamicsMD
!****************************************************************************
!OutputMD
!Modified: ZQ Wang
!	02/16/06: 
!	Added two parameters to control the output frequency for general
!	output and past calculation record.
!		
!
!****************************************************************************
module OutputMD
	use ZQMPI
	use vectors
	implicit none

	type(vector)::tyPrePoint
	integer::iDisOutFreq                ! Determine how often dislocations are written into files.
	integer::iGeneralOutFreq            ! How often general outputs are written into new files.
	integer::iPastFreq                  ! How often past intermediate data should be written into files.
end module OutputMD

!****************************************************************
!  CrossSlipMD:
!               Used for cross-slip subroutines.
!       Updated: ZQ Wang
!       4/8/06:
!       Added a variable lToImplement to determine if cross-slip
!       should be implemented according to the force or probability.
!****************************************************************
MODULE CrossSlipMD
	use ZQMPI
	use vectors
	implicit none

	integer,parameter::iMaxNodeDim=200
	integer,dimension(iMaxNodeDim)::iTmpIndex ! to record indexes of points in DPoints
	integer::itmpNumP,itmpNewNumP
	integer::iTmpLooptype,itmpBurgers,itmpPlane
	integer::iCSCount,iNumTotalScrew,iNumTotalScrew1
	double precision::dpCSLength

	type(DisPoint),dimension(iMaxNodeDim)::tvTmpPL ! tmp dislocation position
	logical::lScrewSeg,lScrewForce,lCrossSlip,lToImplement,CSYN  !lToImplement is for entire loop, CSYN for each screw segment.
	type(vector)::tvLastScrewPoint

        !!!!!
	integer,parameter::ScrewN=6
	double precision,parameter::Beta=1
	double precision::TauIII=32.d0      !in MPa, 32.0MPa for fcc coper, for bcc, calculated in function
	double precision::Length0=2766d0    !for fcc copper, for BCC, change in code
	integer,dimension(3)::iCSMiller
	integer::iNumCSM   !number of available cross slip miller, different for bcc and fcc


	double precision::UScrew(ScrewN,2)
	integer::SegScrew(ScrewN,2)
	type(vector)::PScrew(ScrewN,2)
	integer::NumScrew

	integer::FinalScrewN
	integer::FinalScrewSeg(ScrewN,2)
	double precision::FinalScrewU(ScrewN,2)
	integer::FinalScrewMiller(ScrewN)

end module CrossSlipMD
!**************************************************************************
Module SurfMD
	use ZQMPI
	implicit none

	integer::iSurfCounter

end Module SurfMD
!****************************************************************************
