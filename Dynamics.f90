!************************************************************
!	DynamicsMD:
!		This module contains subroutines that are used to solve
!		dislocation dynamics, i.e., to solve equations of motion, 
!		to do numerical integration.
!************************************************************
Module DynamicsSolverMD
!	use ZQMPI
!	implicit none

	contains
!****************************************************************
!	DynamicsMain:
!		Main subroutine called by simulation to solve dislocation
!	dynamics.
!****************************************************************
SUBROUTINE DynamicsMain
	use ZQMPI
	USE VECTORS
	USE VARIABLES,only:N_times
	IMPLICIT NONE

	INTEGER::IOUT

	!=====Load is applied after this line===========
	!//////////////////////
	!   Solve dynamics
	DO IOUT = 1, N_TIMES
		if(iMyid/=iMaster)then
			!....solve dynamics....
			call ExplicitInt

			!....Compute strain....
			CALL ComputeStrain 

			!.....Check equilibrium conditions.......
			if(lMetaEquilibrium)exit

		end if
	end do

	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

END SUBROUTINE DynamicsMain
!***************************************************************
!  ExplicitInt:   
!	Explicit integration method.
!***************************************************************
subroutine ExplicitInt
	use ZQMPI
	use vectors
	use variables,only:DPoints,zero,iloop
	use DynamicsMD,only:DMatrix
	use MPIMOdule,only:iNTotalDPs,iNVLPoints
	implicit none
	
	integer::i,j,k
	integer::iDim

	!...initialization.....
	do i=1,iNVLPoints
		iDim=4	
		do j=1,iDim
			DMatrix(i)%FG(j)=0.d0
			do k=1,iDim
				DMatrix(i)%KG(j,k)=0.d0
			enddo
			do k=1,8
				DMatrix(i)%QG(j,k)=0.d0
			enddo
		enddo
		DMatrix(i)%tvPLRate=zero%v(1)
		DMatrix(i)%tvTLRate=zero%v(1)
	enddo

	!....Construct matrix for each point....
	call ConstructMatrix

	!....solve matrix....
	call SolveMatrix

	!....numerical integration.....
	call NumerIntegration

end subroutine ExplicitInt
!*****************************************************
! ConstructMatrix:
!	Construct stiffness, force matrix for each dislocation
!	point.
!*****************************************************
subroutine ConstructMatrix
	use ZQMPI
	use vectors
	use variables,only:iNPoint,DPoints
	implicit none

	integer::i
	logical::IsCloseLoopPoint

	do i=1,iNpoint
		if(.not. IsCloseLoopPoint(i) .and. DPoints(i)%lToCal)call ConstructLocal(i)
	enddo

end subroutine ConstructMatrix
!********************************************************
!  SolveMatrix:
!	Stiffness, force matrixes have been constructed and
!	iteration method will be applied to solve the matrix.
!
!	Each point has matrix, its connected points have change
!	rates. 
!********************************************************
subroutine SolveMatrix
	use ZQMPI
	use vectors
	use variables,only:iNPoint,DPoints,zero
	use DynamicsMD,only:DMatrix
	implicit none
	
	integer::I,j,iSolverCount
	logical::IsCloseLoopPoint

	iSolverCount=10
	do i=1,iSolverCount
		do j=1,iNPoint

			if(.not. IsCloseLoopPoint(j) .and. DPoints(j)%lToCal)then
				call SolveLocal(j)
			else
				DMatrix(i)%tvPLRate=zero%v(1)
				DMatrix(i)%tvTLRate=zero%v(1)
				DPoints(i)%tvAcce=zero%v(1)
				DPoints(i)%tvPreV=zero%v(1)
				DPoints(i)%tvPreVT=zero%v(1)
			endif

		enddo
	enddo

end subroutine SolveMatrix
!***************************************************************
!	NumerIntegration:
!		Numerical integration is performed on the local positions.
!	Global positions are then calculated based on local positions.
!
!------------------------------------------------------------------
! Modified: ZQ Wang, 01/25/06
!	Added the subroutine to calculate the effective mass and acceleration
!	of dislocations.
!	01/30/06:
!	The verlet method, first we solve velocity, sencond we solve acceleration,
!	finally, we solve position.
!
!***************************************************************
subroutine NumerIntegration
	use ZQMPI
	use vectors
	use variables,only:iNPoint,DPoints,DTime,DMaxSpeed,iInertial,dpTSoundV,   &
						Miller,GPlanes,iNPlane,iloop,zero
	use DynamicsMD,only:DMatrix
	implicit none

	integer::I,j
	integer::ip1,ip2,im1,im2
	type(vector)::tvTAcce

	!....explicit integration.... 
	do I=1,iNPoint

		tvTAcce=zero%v(1)
		!....control maximum values....

		if(iInertial==0)then   !no inertial effect, consider velocity control
			if(MAG(DMatrix(i)%tvPLRate)>DMaxSpeed)then
				DMatrix(i)%tvPLRate=DMaxSpeed*UV(DMatrix(i)%tvPLRate)
			endif
			if(MAG(DMatrix(i)%tvTLRate)>DMaxSpeed)then
				DMatrix(i)%tvTLRate=DMaxSpeed*UV(DMatrix(i)%tvTLRate)
			endif
		else            !with inertial effect, avoid abrupt changes.
			if( MAG(DPoints(i)%tvAcce)>0 .and. (MAG(DMatrix(i)%tvPLRate)>10d0*MAG(DPoints(i)%tvAcce) .or.     &
                              MAG(DPoints(i)%tvPreV+DTime*DMatrix(i)%tvPLRate)>dpTSoundV))then
				DMatrix(i)%tvPLRate=zero%v(1)    !DPoints(i)%tvAcce
			endif
			if(MAG(DPoints(i)%tvPreVT)>0 .and. MAG(DPoints(i)%tvPreVT+DTime*DMatrix(i)%tvTLRate)>10*MAG(DPoints(i)%tvPreVT) .or. &
                          MAG(DPoints(i)%tvPreVT+DTime*DMatrix(i)%tvTLRate)>dpTSoundV)then
				DMatrix(i)%tvTLRate=zero%v(1)
			endif
		endif

		!....compute the effective mass and acceleration.....01/25/06
                if(iInertial==1)then
                        call computeMassAcce(i,tvTAcce)
                endif

		!....numerical integration....
		do j=1,2
			DPoints(i)%tvPL%v(j)=DMatrix(i)%tvPLRate%v(j)*DTIME+DPoints(i)%tvPL%v(j)-0.5*DPoints(i)%tvAcce%v(j)*DTime*DTime
		enddo
		do j=1,2
			DPoints(i)%tvTL%v(j)=DMatrix(i)%tvTLRate%v(j)*DTIME+DPoints(i)%tvTL%v(j)-0.5*tvTAcce%v(j)*DTime*Dtime
		enddo

		DPoints(i)%tvPL%v(3)=0.d0
		DPoints(i)%tvTL%v(3)=0.d0

		!....compute the global positions
		call setGlobal(i)

		DPoints(i)%tvPreV=DMatrix(i)%tvPLRate
		DPoints(i)%tvPreVT=DMatrix(i)%tvTLRate

	enddo   !iPoint

end subroutine NumerIntegration
!*********************************************************
!	ConstructLocalPre:
!		Construct local matrix for one specific point.
!*********************************************************
subroutine ConstructLocal(iCur)
	use ZQMPI
	use vectors
	use variables,only:DPoints
	use DynamicsMD,only:KL,FL,DMatrix
	implicit none

	integer,intent(in)::iCur
	
	integer::iNext,iPre
	integer::i,j,k,N

	double precision,dimension(8)::FFL
	double precision,dimension(8,8)::KKL

	N=8

	do i=1,8
		FFL(i)=0.d0
		do j=1,8
			KKL(i,j)=0.d0
		enddo
	enddo
		
	do i=1,2
		do j=1,8
			FL(i,j)=0.d0
			do k=1,8
				KL(i,j,k)=0.d0
			enddo
		enddo
	enddo

	iPre=DPoints(iCur)%iBeginP
	iNext=DPoints(iCur)%iEndP
	
	!....beginning segment matrix, local....
	if(iPre/=-1)then
		if(DPoints(iPre)%iEndP/=-1)then
			! previous point has not been calculated.
			! first point. Previous matrix needed to be evaluated.
			i=iPre
			call Force(N,FFL,i)

			call Stiffness(N,KKL,i)

			do i=1,8
				FL(1,i)=FFL(i)
				do j=1,8
					KL(1,i,j)=KKL(i,j)
				enddo
			enddo
		endif
	endif

	!....following segment matrix, local......special for open loop last point....
	if(iNext/=-1)then !

		call Force(N,FFL,iCur)
	
		call Stiffness(N,KKL,iCur)
	
		do i=1,8
			FL(2,i)=FFL(i)
			do j=1,8
				KL(2,i,j)=KKL(i,j)
			enddo
		enddo
	endif
		
	!....Construct the matrix for dislocation point.....
	!....not intersection point, dimension is 4.....
	do i=1,4
		do j=1,4
			DMatrix(iCur)%KG(i,j)=KL(1,i+4,j+4)+KL(2,i,j)
			DMatrix(iCur)%QG(i,j)=KL(1,i+4,j)
			DMatrix(iCur)%QG(i,j+4)=KL(2,i,j+4)
		enddo
		DMatrix(iCur)%FG(i)=FL(1,i+4)+FL(2,i)
	enddo

end subroutine ConstructLocal
!*********************************************************
!	Force:
!		To compute the force matrix for a segment. It is
!	the local one in previous simulations.
!*********************************************************
SUBROUTINE FORCE(NN,FE,PointID)
	use ZQMPI
	USE VECTORS
	USE VARIABLES,only:Burgers,DPoints,GPlanes,Miller,POS,WT,MAX_QUAD, &
						MU,ZERO,CheckNeiBur
	use FunctionMD,only:PointPBC
	IMPLICIT NONE
	
	INTEGER,intent(in):: NN,PointID
	DOUBLE PRECISION, DIMENSION(NN)::FE

	type(vector)::P1,TT1,P2,TT2
	TYPE(VECTOR)::R,RU,RUU,RR
	TYPE(VECTOR)::Q,T_Q,B_Q,MILLER_Q,TOTAL_F_L,Pforce
	TYPE(MATRIX)::SIG_Q,SIG_QQ

	DOUBLE PRECISION,DIMENSION(4,3)::SHAP
	DOUBLE PRECISION::DS,KAPPA,U
	INTEGER::I,II,I_P,I_B,I_M

	I_P=DPoints(PointID)%iPlane
	I_B=DPoints(PointID)%iBURGERS
	I_M=GPlanes(i_p)%iMiller

	B_Q=Burgers(I_B)
	MILLER_Q=Miller(I_M)%V

	DO I=1,NN
		FE(I)=0.0
	ENDDO

	!....assign values....special for intersection points as the second point.....
	! For being as the first point, it is all the same.

	P1=DPoints(PointID)%tvPL
	TT1=DPoints(PointID)%tvTL
	call TransNextLocal(PointID,P2,TT2)
	
	DPoints(PointID)%force=zero%v(1)
		
	DO I=1,MAX_QUAD
		!....geometry at the point....
		U=(POS(I)+1)/2.0D0

		CALL GETSHAPE(SHAP,U)
		R=SHAP(1,1)*P1+SHAP(2, 1)*TT1+SHAP(3,1)*P2+SHAP(4,1)*TT2
		RU=SHAP(1,2)*P1+SHAP(2,2)*TT1+SHAP(3,2)*P2+SHAP(4,2)*TT2
		RUU=SHAP(1,3)*P1+SHAP(2,3)*TT1+SHAP(3,3)*P2+SHAP(4,3)*TT2
		DS=DSQRT(RU*RU)
		RR=RU**RUU
		KAPPA=DSQRT(RR*RR)/(DSQRT(RU*RU))**3
		IF(RR%V(3) .GT. 0) THEN
			KAPPA=-KAPPA
		ENDIF

		!....change tangent vector to global one.....
		T_Q=TRANS(Miller(I_M)%ES)*RU
		T_Q=UV(T_Q)
		Q=TRANS(Miller(I_M)%ES)*R+GPlanes(I_P)%ORIGIN
		Q=PointPBC(Q)

		!....stress matrix....
		SIG_Q=ZERO;SIG_QQ=ZERO

		!....other f-r source, interaction....
		if(CheckNeiBur==1)then
			IF(DPoints(PointID)%iNumNei .NE. 0)THEN
				CALL ComputeInteraction(Q,SIG_QQ,PointID)
				SIG_Q=SIG_Q+SIG_QQ
			ENDIF
		endif		

		CALL COMPUTE_FORCE1(I_M,I,PointID,Q,T_Q,B_Q,MILLER_Q,KAPPA,SIG_Q,TOTAL_F_L,Pforce)
		TOTAL_F_L=MU*TOTAL_F_L

		!...record force......the force for the total length of this segment
		Pforce=MU*Pforce
		DPoints(PointID)%force=DPoints(PointID)%force+WT(i)*0.5*DS*Pforce

		!....local force matrix....
		DO II=1,4 
			FE(II*2-1)=FE(II*2-1)      &
						   +WT(I)*0.5*DS*SHAP(II,1)*(TOTAL_F_L%V(1))
			FE(II*2)=FE(II*2)          &
						   +WT(I)*0.5*DS*SHAP(II,1)*(TOTAL_F_L%V(2))
		ENDDO   
		
	ENDDO
	
END SUBROUTINE FORCE
!**********************************************************
!	Stiffness:
!		To compute the stiffness matrix for a segment. It is
!	the local one in previous simulations.	
!**********************************************************
SUBROUTINE STIFFNESS(NN,KE,PointID)
	use ZQMPI
	USE VECTORS
	USE VARIABLES,only:DPoints
	IMPLICIT NONE

	INTEGER,intent(in):: NN,PointID
	TYPE(VECTOR):: P1,TT1,P2,TT2
	DOUBLE PRECISION,DIMENSION(NN,NN)::KE
	INTEGER ::I,J

	DO I=1,NN
		Do J=1,NN
			KE(I,J)=0.0D0
		ENDDO
	ENDDO

	!....assign values....
	P1=DPoints(PointID)%tvPL
	TT1=DPoints(PointID)%tvTL
	call TransNextLocal(PointID,P2,TT2)

	CALL ELEMENTSTIFF(PointID,KE,P1,TT1,P2,TT2)

END SUBROUTINE STIFFNESS
!*************************************************************
!	ElementStiff:
!		Real calculation of element stiffnees matrix
!------------------------------------------------------------
! Modified: ZQ Wang, 02/07/06
!	Changed stiffness matrix calculation for inertial effect.
! 	Basically, an effective mobility is defined.
!
!*************************************************************
SUBROUTINE ELEMENTSTIFF(PointID,KE,P1,TT1,P2,TT2)
	use ZQMPI
	USE VECTORS
	USE VARIABLES,only:MAX_QUAD,POS,WT,MOBILITY,DPoints,iInertial,DTime,iMATType,iBCC,MOB_R,Miller,GPlanes
	USE functionMD,only:funcIsScrew
	IMPLICIT NONE

	integer,intent(in)::PointID
	TYPE(VECTOR),intent(in):: P1,P2,TT1,TT2
	type(vector)::R,RU,T_Q
	
	DOUBLE PRECISION, DIMENSION(8,8)::KE
	DOUBLE PRECISION, DIMENSION(4,3)::SHAP
	DOUBLE PRECISION :: DS,U,dpMass,dpEMobility,dpTmpR  !dpMass is effective mass, dpFMass is a function
							! dpEMobility is effective mobility

	INTEGER :: I,J,II,JJ
	


	DO I=1,8;DO J=1,8
		KE(I,J)=0.0D0
	ENDDO; ENDDO

	DO I=1,MAX_QUAD

		!
		U=(POS(I)+1)*0.5d0
		CALL GETSHAPE(SHAP,U)
		R=SHAP(1,1)*P1+SHAP(2,1)*TT1+SHAP(3,1)*P2+SHAP(4,1)*TT2
		RU=SHAP(1,2)*P1+SHAP(2,2)*TT1+SHAP(3,2)*P2+SHAP(4,2)*TT2
		T_Q=TRANS(Miller(GPlanes(DPoints(PointID)%iPlane)%iMiller)%ES)*RU
		DS=DSQRT(RU*RU)

		!
		
		if(iMatType==iBCC .and. .not. funcIsScrew(PointID,T_Q))then !for mobility for BCC case
			dpTmpR=MOB_R
		else
			dpTmpR=1.0   !for mobility (edge and screw) for FCC case
		endif
		!
		if(iInertial==1)then
			! New Method, April 20, 2006
			dpMass=dpFMass(PointID,I)
			dpEMobility=dpMass
		else
			dpEMobility=MOBILITY*dpTmpR
		endif

		!
		DO II=1,4      !  1,4
			DO JJ=1,4
				KE(II*2-1,JJ*2-1)=KE(II*2-1,JJ*2-1)   &
			          +WT(I)*0.5*dpEMobility*SHAP(II,1)*SHAP(JJ,1)*DS

				KE(ii*2,jj*2)=KE(ii*2,jj*2) &
					  +WT(I)*0.5*dpEMobility*SHAP(II,1)*SHAP(JJ,1)*DS 
			ENDDO
		ENDDO   
	ENDDO

END SUBROUTINE ELEMENTSTIFF
!*************************************************************
!	SolveLocal:
!		solve local matrix, applying boundary conditions
!*************************************************************
subroutine SolveLocal(PointID)
	use ZQMPI
	use vectors
	use variables,only:DPoints,iNPoint,iDPTypeFixed,iDPTypeFree,GPlanes,Miller
	use DynamicsMD,only:DMatrix
	implicit none

	integer,intent(in)::PointID
	double precision,dimension(8)::dpTmpQDot
	double precision,allocatable,dimension(:,:)::KLL
	double precision,allocatable,dimension(:)::FLL
	integer::i,j,k,N,i_p1,i_p2,i_m1,i_m2
	integer::iBegin,iEnd,iType

	double precision::dpRatio
	type(vector)::tmpVector

	!initialization
	iBegin	=	DPoints(pointID)%iBeginP
	iEnd	=	DPoints(pointID)%iEndP
	iType	=	DPoints(pointID)%iType
	
	N=4
	
	allocate(KLL(N,N),FLL(N))

	do i=1,N
		FLL(i)=0.d0
		do j=1,N
			KLL(i,j)=0.d0
		enddo
	enddo
	do i=1,8
		dpTmpQDot(i)=0.d0
	enddo

	!
	! here apply connection boundary conditions 1
	if(iBegin==-1)then
		dpTmpQDot(1)=0.d0
		dpTmpQDot(2)=0.d0
		dpTmpQDot(3)=0.d0
		dpTmpQDot(4)=0.d0
	else
		dpTmpQDot(1)=DMatrix(iBegin)%tvPLRate%v(1)
		dpTmpQDot(2)=DMatrix(iBegin)%tvPLRate%v(2)
		dpTmpQDot(3)=DMatrix(iBegin)%tvTLRate%v(1)
		dpTmpQDot(4)=DMatrix(iBegin)%tvTLRate%v(2)
	endif

	if(iEnd==-1)then
		dpTmpQDot(5)=0.d0
		dpTmpQDot(6)=0.d0
		dpTmpQDot(7)=0.d0
		dpTmpQDot(8)=0.d0
	else
		dpTmpQDot(5)=DMatrix(iEnd)%tvPLRate%v(1)
		dpTmpQDot(6)=DMatrix(iEnd)%tvPLRate%v(2)
		dpTmpQDot(7)=DMatrix(iEnd)%tvTLRate%v(1)
		dpTmpQDot(8)=DMatrix(iEnd)%tvTLRate%v(2)
	endif
		
	! copy matrix and compute FL
	do i=1,N
		do j=1,N
			KLL(i,j)=DMatrix(PointID)%KG(i,j)
		enddo

		FLL(i)=DMatrix(PointID)%FG(i)
		do j=1,8
			FLL(i)=FLL(i)-DMatrix(PointID)%QG(i,j)*dpTmpQDot(j)
		enddo
	enddo

	! here apply boundary condition 2, fixed boundary conditions
	if(iType==iDPTypeFixed)then
		do i=1,2
			do j=1,4
				KLL(i,j)=0.d0
			enddo
			FLL(i)=0.d0
			KLL(i,i)=1.d0
		enddo
	endif

	!
	call Gauss1(N,KLL,FLL)

	!update rates
	if(iType==iDPTypeFixed)then
		DMatrix(PointID)%tvPLRate%v(1)=0.d0
		DMatrix(PointID)%tvPLRate%v(2)=0.d0
	else
		DMatrix(PointID)%tvPLRate%v(1)=FLL(1)
		DMatrix(PointID)%tvPLRate%v(2)=FLL(2)
	endif
	DMatrix(PointID)%tvPLRate%v(3)=0.d0
	
	DMatrix(PointID)%tvTLRate%v(1)=FLL(3)
	DMatrix(PointID)%tvTLRate%v(2)=FLL(4)
	DMatrix(PointID)%tvTLRate%v(3)=0.d0

	deallocate(FLL,KLL)

end subroutine SolveLocal
!**********************************************************************
!	Compute_force1:
!		Compute all the forces acting on the dislocation segments.
!		Inertial_f is in the local system.
!------------------------------------------------------------------------
! Modified: ZQ WANG, 01/25/06
!	Added the inertial force calculation for high speed dislocations.
!	Added iQuad, iInertial, INERTIAL_F, and the subroutine
!**********************************************************************
SUBROUTINE COMPUTE_FORCE1(I_M,iQuad,PointID,Q,T_Q,B_Q,MILLER_Q,KAPPA,SIG_Q,TOTAL_F_L,Pforce)
	use ZQMPI
	USE VECTORS
	USE VARIABLES,only:NU,SIG_APP,PI,Miller,iObstacle,MU,zero,iInertial,iMatType,IBCC,PEIERLS,PEI_R,DPoints, &
				GPlanes,ATRatio, TRatio
 	use functionMD,only:funcIsScrew, funcIsAntitwin
	implicit none

	integer,intent(in)::I_M,iQuad,PointID
	TYPE(VECTOR),intent(in)::Q,T_Q,B_Q,MILLER_Q
	TYPE(VECTOR)::APPLIED_F,SELF_F,PK_F,PEIERLS_F,INERTIAL_F,TOTAL_F_L,Pforce
	TYPE(VECTOR)::TMPVEC
	TYPE(MATRIX)::SIG_Q
	
	DOUBLE PRECISION,intent(in)::KAPPA
	DOUBLE PRECISION::CA,C2A,E_ALPHA,E2_ALPHA,C1,C2,C3
	DOUBLE PRECISION::DPObstacle
	DOUBLE PRECISION::PEIERLS_L,tmpF
	double precision::Per_TwinRate    !ratio for twinning and anti-twinning Peierls stresses.
	!///////////////////////////////////////////////////////////////////////
	!Peach-Kohler, interaction force from other dislocations and obstacles.
        PK_F=(SIG_Q*B_Q)**T_Q

	!////////////////////////////////
	!Self-force
	IF (DABS(KAPPA) .NE. 0) THEN
		CA=T_Q*UV(B_Q)
		C2A=2*CA*CA-1      
		E_ALPHA=B_Q*B_Q/(4*PI*(1-NU))*(1-NU*CA*CA)
		E2_ALPHA=B_Q*B_Q*NU/(2*PI*(1-NU))*C2A
		C1=(E_ALPHA+E2_ALPHA)*DLOG(2*8/MAG(B_Q)/DABS(KAPPA))
		C2=B_Q*B_Q*(21+CA*CA)/(64*PI)
		C3=B_Q*B_Q*C2A/(2*PI)
		SELF_F=KAPPA*(-C1+C2+C3)*(UV(MILLER_Q)**T_Q)
	ELSE
		SELF_F%V(1)=0.0D0  
		SELF_F%V(2)=0.0D0  
		SELF_F%V(3)=0.0D0  
	ENDIF  

	!////////////////////////////////
	!Applied force
	APPLIED_F=(SIG_APP*B_Q)**T_Q

	!////////////////////////////////
	! 12/17/08
	!Non-glide force
!	APPLIED_F=(1.0+NonGlide_F(PointID,APPLIED_F,B_Q,T_Q,I_M))*APPLIED_F

	!////////////////////////////////
	!Peierls force
	PEIERLS_F=zero%v(1)

	!/////////////////////////////////////
	!Inertial force, return in local form.
	INERTIAL_F=zero%v(1)
	if(iInertial==1)then
		call ComputeInertialForce(iQuad,PointID,T_Q,INERTIAL_F)
	endif

	!/////////////////////////////////////
	!total force
	TOTAL_F_L=Miller(i_m)%ES*(APPLIED_F+SELF_F +PK_F)+INERTIAL_F !change it to local
	Pforce=Applied_F+PK_F

	!/////////////////////////////////////////
	!Here to determine if the force in the direction of applied force is larger than the Peierls barrier.
	!If not, it should be say that the total force is zero, otherwise the total force is the applied force
	!substracting the Peierls barrier.
	!
	!For <112> planes, twinning and anti-twinning directions should have different Peierls stresses.
	if(iMatType==iBCC)then
		tmpF=MAG(Total_F_L)

		!If not <112> plane, Per_TwinRate=1. Peierls stress is as original
		!If <112> twinning direction, Per_TwinRate=0.75: easy direction;
		!If <112> anti-twinning direction, Per_TwinRate=1.25: hard direction;
		
		if(GPlanes(DPoints(PointID)%iPlane)%iMiller<=6)then   !<110>
			Per_TwinRate=1.0
		elseif(funcIsAntitwin(PointID,Total_F_L))then   !<112> hard
			! Per_TwinRate= ATRatio  !5/14/09
			Per_TwinRate=fun112ATP(SIG_Q,SIG_APP,B_Q,T_Q,i_m)
		else                                      !<112> easy
			Per_TwinRate= TRatio   !0.25
		endif


		if(funcIsScrew(PointID,T_Q))then   !screw dislocations
			Peierls_L=PEIERLS*MAG(B_Q)*Per_TwinRate
		else					!edge dislocations, Peierls force can be changed
			Peierls_L=PEIERLS*MAG(B_Q)*PEI_R !*Per_TwinRate
		endif

		if(tmpF<=PEierls_L)then   !applied force less than barrier, don't move, force is zero
			Total_F_L=zero%v(1)
		else				!applied force larger than barrier, move, with reduced force
			Total_F_L=(tmpF-Peierls_L)*UV(Total_F_L)
		endif
	endif

	!Distributed obstacle force, updated by ZQ Wang, 11/11/2005
	if(iObstacle==1)then
		TMPVEC=TOTAL_F_L
		TMPVEC%v(3)=0.d0
		DPObstacle=40D6/MU*MAG(B_Q)
		if(MAG(TMPVEC)>DPOBstacle)then
			TOTAL_F_L=TOTAL_F_L+(-1.d0)*DPObstacle*UV(TMPVEC)
		else
			TOTAL_F_L=zero%v(1)
		endif
	endif

END SUBROUTINE COMPUTE_FORCE1 
!*******************************************************************************
!	NonGlide_F
!		To account for the non-glide force effect on screw dislocations
!	return (tau_2+tau_3)*0.5/tau_1. 0.5 could be other values.
!*******************************************************************************
double precision function  NonGlide_F(PointID,APPLIED_F,B_Q,T_Q,i_m)
	use ZQMPI
	use vectors
	use variables,only:Miller,DPoints,StrainRate,NGRatio
	implicit none


	type(vector)::B_Q,T_Q,APPLIED_F
	integer,intent(in)::i_m,PointID

	integer::i,j,I_B
	integer::i_mm(2)
	double precision::sign_mm(2)
	type(vector)::APP_f(3)   !on the other 2 planes


	NonGlide_F=0.0

	!Check planes
	I_B=DPoints(PointID)%iBurgers

        if(I_B==3 .or. I_B==7)then
                if(i_m==1)then
                        i_mm(1)=3;i_mm(2)=6
                        sign_mm(1)=1;sign_mm(2)=-1
                elseif(i_m==3)then
                        i_mm(1)=1;i_mm(2)=6
                        sign_mm(1)=1;sign_mm(2)=-1
                elseif(i_m==6)then
                        i_mm(1)=1;i_mm(2)=3
                        sign_mm(1)=-1;sign_mm(2)=-1
                endif
        elseif(i_b==4 .or. i_b==8)then
                if(i_m==1)then
                        i_mm(1)=4;i_mm(2)=5
                        sign_mm(1)=1;sign_mm(2)=-1
                elseif(i_m==4)then
                        i_mm(1)=1;i_mm(2)=5
                        sign_mm(1)=1;sign_mm(2)=-1
                elseif(i_m==5)then
                        i_mm(1)=1;i_mm(2)=4
                        sign_mm(1)=-1;sign_mm(2)=-1
                endif
        elseif(i_b==1 .or. i_b==5)then
                if(i_m==2)then
                        i_mm(1)=4;i_mm(2)=6
                        sign_mm(1)=-1;sign_mm(2)=-1
                elseif(i_m==4)then
                        i_mm(1)=2;i_mm(2)=6
                        sign_mm(1)=-1;sign_mm(2)=1
                elseif(i_m==6)then
                        i_mm(1)=2;i_mm(2)=4
                        sign_mm(1)=-1;sign_mm(2)=1
                endif
        elseif(i_b==2 .or. I_b==6)then
                if(i_m==2)then
                        i_mm(1)=3;i_mm(2)=5
                        sign_mm(1)=-1;sign_mm(2)=-1
                elseif(i_m==3)then
                        i_mm(1)=2;i_mm(2)=5
                        sign_mm(1)=-1;sign_mm(2)=1
                elseif(i_m==5)then
                        i_mm(1)=2;i_mm(2)=3
                        sign_mm(1)=-1;sign_mm(2)=1
                endif
	endif
		

	!Local force
	do i=1,2
		APP_f(i)=Miller(i_mm(i))%ES*APPLIED_F
		APP_f(i)%v(3)=0.0
	enddo
	APP_f(3)=APPLIED_F
	APP_f(3)%v(3)=0.0

	!Ratio: should include sign
	NonGlide_F=(MAG(APP_f(1))+MAG(APP_f(2)))/MAG(APP_f(3))*NGRatio
	if(StrainRate<0)then
		NonGlide_f=-1.0*NonGlide_f
	endif

end function NonGlide_F
!*************************************************************
!	ComputeInteraction
!		Compute the interaction between dislocation neighbors.
!*************************************************************
subroutine ComputeInteraction(Q,SIG_Q,PointID)
	use ZQMPI
	use vectors
	use variables,only:DPoints,GPlanes,Burgers,Miller,interval,pos,wt,RMin,PI,NU,MAX_QUAD,&
						A_Cube,D
	use FunctionMD,only:PointPBC
	implicit none

	type(vector),intent(in)::Q
	type(matrix)::SIG_Q
	integer,intent(in)::PointID

	integer::i,j,k,I_B,I_M,iNext

	DOUBLE PRECISION::R,S,TM,FACTOR0, FACTOR1
	double precision,dimension(4,3)::shap
	double precision::getDist,u

	TYPE(VECTOR)::E,P,T,T_G,BUR
	type(vector)::PG1,PG2,TG1,TG2
	TYPE(MATRIX)::TBE,BET

	!....loop for neighbors....
	do i=1,DPOints(PointID)%iNumNei
		
		j=DPoints(PointID)%iaIDOfNei(i)

		if(DPoints(PointID)%iLoopID/=DPoints(j)%iLoopID .and. DPoints(j)%iEndP/=-1)then

			I_B		=	DPoints(j)%iBurgers
			I_M		=	GPlanes(DPoints(j)%iPlane)%iMiller
			iNext		=	DPoints(j)%iEndP
			BUR		=	BURGERS(I_B)
			PG1		=	DPoints(j)%tvPG
			TG1		=	DPoints(j)%tvTG
			call TransNextGlobal(j,PG2,TG2)

			!----calculate the force-----
			DO K=1,MAX_QUAD
				U=(POS(K)+1)*0.5d0
				call getshape(shap,u)
				FACTOR0=INTERVAL*WT(K)
				
				P = shap(1,1)*PG1+shap(2,1)*TG1     &
				   +shap(3,1)*PG2+shap(4,1)*TG2
				P = PointPBC(P)

				T_G=shap(1,2)*PG1+shap(2,2)*TG1   &
				   +shap(3,2)*PG2+shap(4,2)*TG2

				R=GetDist(Q,P)
				if(R<1.d-6)goto 1
				E=GetUV(Q,P,A_Cube)

				TM=MAG(T_G)
				if(tm==0.0d0)goto 1
				T=UV(T_G)

				IF (R<RMIN) THEN
					R=RMIN
				END IF

				TBE=T**BUR/E
				TBE=TBE+TRANS(TBE)
				BET=BUR**E/T
				BET=BET+TRANS(BET)
				S=T**BUR*E

				FACTOR1=FACTOR0*TM/(4*PI*R*R)

				SIG_Q=SIG_Q+FACTOR1*(BET+(TBE-S*D-3*S*(E/E))/(1-NU))

1				continue

			ENDDO !MAX_QUAD
		endif
	enddo

end subroutine ComputeInteraction
!**********************************************************
! ComputeInertialForce:
!	Compute the inertial force acting on dislocations for
!	high speed deformation.
!
!	For each integral point, their mass, acceleration are 
!	calculated first and then the inertial force are computed
!	from them.
!-----------------------------------------
!	Q: Position of the point
!	iQuad: ID of the quadrature point
!	PointID: ID of the segment.
!	INTERTIAL_F: inertial force 
!-----------------------------------------
! Modified: ZQ Wang, 1/25/06
!	New.
!    Inertial force is now -BV+MV/DT
!**********************************************************
subroutine ComputeInertialForce(iQuad,PointID,T_Q,INERTIAL)
	use ZQMPI
	use vectors
	use variables,only:DPoints,pos,DTIme,MU,MOBILITY,MOB_R,iMatType,IBCC
	USE functionMD,only:funcIsScrew
	implicit none

	type(vector)::INERTIAL,T_Q
	integer::iQuad,PointID

	double precision::dpTMass
	type(vector)::tvPreV

	double precision::Upos,dpTmpR
	integer::iNext
	double precision::shap(4,3)

	!//////////////////////////////////////////////////
	! Mass and acceleration are calculated first.
	Upos=(pos(iQuad)+1)/2
	iNext=DPoints(PointID)%iEndP
	dpTMass=DPoints(PointID)%dpMass*(1-Upos)+DPoints(iNext)%dpMass*UPos
	call getshape(shap,Upos)
	tvPreV=shap(1,1)*DPoints(PointID)%tvPreV+shap(2,1)*DPoints(PointID)%tvPreVT+      &
		shap(3,1)*DPoints(iNext)%tvPreV+shap(4,1)*DPoints(iNext)%tvPreVT

	!////////////////////////////////////////////////////////////
	! April 20, 2006. New Method
	if(iMatType==iBCC .and. .not. funcIsScrew(PointID,T_Q))then !for mobility for BCC case
		dpTmpR=MOB_R
	else
		dpTmpR=1.0   !for mobility (edge and screw) for FCC case
	endif

	INERTIAL=(-1.d0*MOBILITY*dpTmpR)*tvPreV/MU

end subroutine ComputeInertialForce
!**********************************************************
! computeMassAcce:
!	Compute the effective mass and acceleration of a dislocation
!	based on its velocity. The mass is coming from both edge 
!	and screw character.
!-----------------------------------------------------------
! Modified: ZQ Wang, 01/25/06
!	New
!    02/16/06:
!	added ratio control. The reason here is to avoid numerical
!	instabilities. It has been found that when the velocity is 
!	relatively low, less than 1d-3*SoundV, the mass doesn't 
!	change much from the mass at 1d-3*SoundV.
!**********************************************************
subroutine computeMassAcce(pointID,tvTAcce)
	use ZQMPI
	use vectors
	use variables,only:DPoints,dpLSoundV,dpTSoundV,dpEMass0,MU,BURGERS,  &
						DTime,PI,zero
	use DynamicsMD,only:DMatrix
	implicit none

	integer,intent(in)::pointID
	double precision::dpMagV,dpGamma,dpGammaL
	double precision::dpFEdge,dpFScrew
	double precision::vectorAngle, angle
	double precision::dpMassE,dpMassS
	double precision::ratio,ratioL
	type(vector)::tvTmp1,tvTmp2
	type(vector)::tvTAcce


!////////////////////////////////////////////////////////
	!Transfer tacce and calculate tvPlrate
	DPoints(PointID)%tvAcce=DMatrix(PointID)%tvPLRate
	
!	if(DPoints(PointID)%iBeginP==-1 .or. DPoints(PointID)%iEndP==-1)then
!		tvTAcce=0.01d0*DMatrix(PointID)%tvTLRate
!	else
		tvTAcce=DMatrix(PointID)%tvTLRate
!	endif
	
	DMatrix(PointID)%tvPLRate=DPoints(PointID)%tvPreV+DTime*DPoints(PointID)%tvAcce
	DMatrix(PointID)%tvTLRate=DPoints(PointID)%tvPreVT+DTime*tvTAcce
	
! April 20, 2006, above, new method
!-----------------------------------------------------------------------------------------------	
	!////////////////////////////////////////////////////////
	! Determine the percentage from edge and screw
	tvTmp1=BURGERS(DPoints(PointID)%iBurgers)   !direction of burgers vector
	tvTmp2=DPoints(PointID)%tvTG               !direction of tangent vector
	angle=vectorAngle(tvTmp1,tvTmp2)/180.0*PI
	dpFScrew=(dcos(angle))**2
	dpFEdge=(dsin(angle))**2

	!///////////////////////////////////////
	!	Velocity
	dpMagV=MAG(DMatrix(PointID)%tvPLRate)

	if(dpMagV>0.99*dpTSoundV)then   ! 02/17/06
		DMatrix(PointID)%tvPLRate=DPoints(PointID)%tvPreV
		dpMagV=MAG(DMatrix(PointID)%tvPLRate)
	endif

	ratio=dpMagV/dpTSoundV
	ratioL=dpMagV/dpLSoundV
	if(ratio<1d-3)ratio=1d-3
	if(ratioL<1d-3)ratioL=1d-3
	dpGammaL=0.d0
	dpGamma=0.d0

	!///////////////////////////////////////
	!
	if(dpMagV/=0.d0)then
		dpGammaL=dsqrt(dabs(1-ratioL**2))
		dpGamma=dsqrt(dabs(1-ratio**2))
	
		dpMassS=dpEMass0/ratio**2*(-1/dpGamma+1/dpGamma**3) 
		dpMassE=dpEMass0/ratio**4                            &
			*(-8*dpGammaL-20/dpGammaL+4/dpGammaL**3+7*dpGamma+25/dpGamma-11/dpGamma**3+3/dpGamma**5)

		!///////////////////////////////////////////////////////
		!mass
		DPoints(PointID)%dpMass=dpMassS*dpFScrew+dpMassE*dpFEdge

		!/////////////////////////////////////////////////////////
		! acceleration, the acceleration is in the local system
		DPoints(PointID)%tvAcce=1.d0/DTime*(DMatrix(PointID)%tvPLRate-DPoints(PointID)%tvPreV)
	else
		DPoints(PointID)%dpMass=dpEMass0
		DPoints(PointID)%tvAcce=zero%v(1)
	endif

end subroutine computeMassAcce
!*************************************************************
! dpFMass:
! Dislocation mass at a point on dislocation loop.
!----------------------------------
! Modified:ZQ Wang,1/30/06
!	New
!*************************************************************
double precision function dpFMass(PointID,iQuad)
	use ZQMPI
	use vectors
	use variables,only:DPoints,pos
	implicit none

	integer,intent(in)::PointID,iQuad
	double precision::Upos
	integer::iNext

	Upos=(pos(iQuad)+1.0)/2.0
        iNext=DPoints(PointID)%iEndP
        dpFMass=DPoints(PointID)%dpMass*(1-Upos)+DPoints(iNext)%dpMass*UPos

end function dpFMass
!***************************************************************************
! lMetaEquilibrium
!   This file is to check if all the points have reached an equilibrium state
!   so that the computation does not have to continue.
!****************************************************************************
logical function lMetaEquilibrium()
	use ZQMPI
	use vectors
	use variables,only:iNPoint,DPoints,DTime
	implicit none

	integer::i
	double precision::dpE,dpE2

	dpE=1d9
	dpE2=2*dpE*dpE
	lMetaEquilibrium=.true.

	do i=1,iNPoint
		if(MAG(DPoints(i)%tvPreV)>dpE .or. MAG(DPoints(i)%tvPreVT)>dpE .or. MAG(DPoints(i)%tvAcce)>dpE2)then
		lMetaEquilibrium=.false.
			exit
		endif
	enddo

end function lMetaEquilibrium
!**************************************************************************
! Normal vector of stress on screw to determine the Peierls stress
!**************************************************************************
double precision function fun112ATP(SIG_Q,SIG_APP,B_Q,T_Q,i_m)
	use ZQMPI
	use vectors
	use variables,only:Peierls,Miller,Loading_Dir
	implicit none

	type(matrix)::SIG_Q,SIG_APP
	type(vector)::B_Q,T_Q
	type(matrix)::SIG_Tmp
	type(vector)::Vec_tmp
	type(vector)::UvN,UvB,UvNor
	double precision:: dp1,dp2,dpRate
	integer::i_m
	integer::i_sign

	SIG_tmp=SIG_Q+SIG_APP
	vec_tmp=SIG_Tmp*B_Q-((SIG_Tmp*B_Q)**T_Q)  !this is normal
	dp1=MAG(vec_tmp)
	dp2=MAG(Peierls*B_Q)

	UvN=UV(Miller(i_m)%v)
	UvB=UV(B_Q)
	UvNor=UV(vec_tmp)

	if((UvN**UvB)*UvNor>0)then
		i_sign=1
	else
		i_sign=-1
	endif

	if(Loading_Dir%v(1)/=0 .and. Loading_Dir%v(2)==0.0 .and. Loading_Dir%v(3)==0.0)i_sign=-1*i_sign

	dp1=i_sign*dp1

	dpRate=dp1/dp2

	if(dpRate<-0.055)then
		fun112ATP=5000
	elseif(dpRate>=-0.055 .and. dpRate<=0.46)then
		fun112ATP=1.12*exp(0.0567/(dpRate-0.5))
	elseif(dpRate>0.46)then
		fun112ATP=0.275
	endif


end function fun112ATP
!**************************************************************************
end module DynamicsSolverMD
