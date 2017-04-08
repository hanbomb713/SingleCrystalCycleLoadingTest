!***********************************************************
!   This file contains subroutines that deal with output.    *
!***********************************************************
!***********************************************************
!   OutputDis: Dislocations output                         *
!***********************************************************
Subroutine OutputDis
	use vectors
	implicit none

	external OutputSegment

	call SurfLoop(OutputSegment)

End subroutine OutputDis
!**************************************************************
!	OutputSegment
!		Output the segment related to current point.
!**************************************************************
subroutine OutPutSegment(PointID)
	use ZQMPI
	use vectors
	use variables,only:GPlanes,Burgers,DPoints,max_node
	use OUTPUTMD,only:tyPrePoint
	implicit none

	integer,intent(in)::PointID
	
	integer::I,j,mapstop
	integer::iEnterID,iCur,iBegin,iEnd
	logical::lContinue
	logical::lSegCriteria
	logical::StopSurfLoop

	integer::iInter
	double precision::dpInter

	if(max_node>10)then
		iInter=6
		dpInter=0.2
	elseif(max_node==10)then
		iInter=11
		dpInter=0.1
	else
		iInter=6
		dpInter=0.2
	endif

	lSegCriteria=.true.
	!.....Output of initial dislocation microstructure....

	tyPrePoint%v(1)=1d10
	tyPrePoint%v(2)=1d10
	tyPrePoint%v(3)=1d10

	iEnterID	=	PointID
	iCur		=	PointID
	lContinue	=	.true.

	do while(lContinue)
		iEnd	=	DPoints(iCur)%iEndP
		DPoints(iCur)%lStat=.true.		

		!stop
		lContinue=.not.StopSurfLoop(iCur,iEnterID,lSegCriteria)

		call Map_Segment(iCur,iInter,dpInter)
		iCur=iEnd

	enddo

end subroutine OutputSegment

!******************************************************************
!	Map_Segment
!		segment has been determined in Outputsegment. Now, it
!	is going to be mapped
!******************************************************************
subroutine Map_segment(PointID,iInter,dpInter)
	use ZQMPI
	use vectors
	use variables,only:DPoints,GPlanes, &
						Miller,ConnectVec
	implicit none

	integer,intent(in)::PointID
	integer,intent(in)::iInter
	double precision,intent(in)::dpInter

	type(vector)::pointL1,pointL2,pointL3,pointL4,TangL3,TangL4,RR

	double precision::du1,du2,u0,shap(4,3)
	integer::i_p,i_m
	integer::i,j,k,mapstop
	integer::iNext

	iNext=DPoints(PointID)%iEndP
	i_p=DPoints(PointID)%iPlane
	i_m=GPlanes(i_p)%iMiller

	!....global output....
	PointL3=DPoints(pointID)%tvPG
	TangL3=DPoints(pointID)%tvTG
	call TransNextGlobal(pointID,PointL4,Tangl4)

	DO i=1,iInter
		u0=dpInter*(i-1)
		call getshape(shap,u0)
		RR=SHAP(1,1)*PointL3+SHAP(2,1)*TangL3 &
			+SHAP(3,1)*PointL4+SHAP(4,1)*TangL4

		if(I==1)then
			PointL1=RR;du1=u0
		else if(I.ne.1)then
			PointL2=RR;du2=u0 !z
			mapstop=0
			call map_point(PointL1,PointL2,PointL3,PointL4,TangL3,   & !z  
							TangL4,du1,du2,mapstop) !Z
			PointL1=PointL2;du1=du2 !z
		end if

	ENDDO
end subroutine Map_Segment
!***********************************************************************
!   Map_point: output with periodic boundary conditions.               *
!***********************************************************************
recursive subroutine map_point(Point1,Point2,Point3,Point4,Tang3,Tang4,   &
							    U1,U2,mapstop)
	use ZQMPI
	use vectors
	implicit none

	Type(vector),intent(in)::Point1,Point2,Point3,Point4,Tang3,Tang4
	integer::mapstop
	DOUBLE PRECISION,intent(in)::u1,u2

	double precision::shap(4,3)
	double precision::U0,UU0
	integer::ii
	integer::i,j,k,MNUM
	TYPE(VECTOR)::CRITICAL_VAL,RR,RR1
	TYPE(vector),dimension(20)::MID_CONNECT_P
	double precision,dimension(20)::UUU

	mapstop=mapstop+1

	call IS_IN_SAME(Point1,Point2,point3,Point4,Tang3,Tang4,U1,U2, &
					CRITICAL_VAL,ii)

	if(mapstop>=10)then
		call map_piece(point1,point2)
	else IF(ii==0)THEN !NOT IN THE SAME CUBE,case 3
		CALL FIND_MIDPOINTS(POINT1,POINT2,POINT3,POINT4,TANG3,TANG4,U1,U2,    &
		CRITICAL_VAL,MID_CONNECT_P,UUU,MNUM)

		RR1=point1
		UU0=U1
		do i=1,MNUM
			call Map_point(RR1,MID_CONNECT_P(i),Point3,Point4,   &
							Tang3,Tang4,UU0,UUU(I),mapstop)
			RR1=MID_CONNECT_P(i)
			UU0=UUU(I)
		enddo

		call MAP_POINT(RR1,POINT2,Point3,Point4,Tang3,Tang4,UU0,U2,mapstop)
	ELSE !IN THE SAME CUBE,case 1&2, 2 can be 3
		call Map_piece(Point1,Point2)
	ENDIF
	return
END SUBROUTINE map_point
!***********************************************************************
!  IS_In_SAME:                                                         *
!***********************************************************************
SUBROUTINE IS_IN_SAME(POINT1,POINT2,Point3,Point4,Tang3,tang4,U1,U2, &
						CRITICAL_VAL,I)
	use ZQMPI
	USE VECTORS
	USE VARIABLES,only:A_CUBE
	IMPLICIT NONE

	TYPE(VECTOR),intent(in)::POINT1,POINT2,point3,point4,Tang3,Tang4
	double precision,intent(in)::U1,U2
	INTEGER,intent(out)::I !I=0 MEANS NOT IN THE SAME CUBE
	TYPE(VECTOR),intent(out)::CRITICAL_VAL !=1 MEANS HAS BOUNDARY VALUES AND NOT IN THE SAME CUBE
	INTEGER::J,K,L
	INTEGER::CONTROL_A,CONTROL_B
	TYPE(VECTOR),DIMENSION(2)::TMP_P
	INTEGER,DIMENSION(2,3,3)::CUB_IN

	CRITICAL_VAL%V(1)=1.D0
	CRITICAL_VAL%V(2)=1.D0
	CRITICAL_VAL%V(3)=1.D0

	! NOT IN THE SAME CUBE,HAS BOUNDARY VALUES
	TMP_P(1)=POINT1
	TMP_P(2)=POINT2

	I=1 ! DEFAULT,IN THE SAME CUBE

	! DETERMINE WHICH CUBE THE POINT IS IN
	DO J=1,2 ! point
		DO K=1,3 ! coordinate
			IF(TMP_P(J)%V(K)>=0)THEN
				CONTROL_A=INT(TMP_P(J)%V(K)/A_CUBE)
				CONTROL_B=CONTROL_A+1
				CUB_IN(J,K,1)=CONTROL_A
				CUB_IN(J,K,2)=CONTROL_A
				CUB_IN(J,K,3)=CONTROL_A

				IF(DABS(TMP_P(J)%V(K)-CONTROL_A*A_CUBE)<2)THEN
					CUB_IN(J,K,2)=CONTROL_A-1
				ENDIF
				IF(DABS(TMP_P(J)%V(K)-CONTROL_B*A_CUBE)<2)THEN
					CUB_IN(J,K,3)=CONTROL_B
				ENDIF
			ELSE
				CONTROL_A=INT(TMP_P(J)%V(K)/A_CUBE)-1
				CONTROL_B=CONTROL_A-1
				CUB_IN(J,K,1)=CONTROL_A
				CUB_IN(J,K,2)=CONTROL_A
				CUB_IN(J,K,3)=CONTROL_A

				IF(DABS(TMP_P(J)%V(K)-CONTROL_A*A_CUBE)<2)THEN
					CUB_IN(J,K,2)=CONTROL_B
				ENDIF
				IF(DABS(TMP_P(J)%V(K)-(CONTROL_A+1)*A_CUBE)<2)THEN
					CUB_IN(J,K,3)=CONTROL_A+1
				ENDIF
			ENDIF
		ENDDO
	ENDDO

	! LOOK TO SEE IF THE TWO POINTS ARE IN THE SAME CUBE
	DO J=1,3
		DO K=1,3
			DO L=1,3
				IF(CUB_IN(1,J,K)==CUB_IN(2,J,L))THEN
					CRITICAL_VAL%V(J)=0.D0 !IN THE SAME CUBE
					GOTO 4
				ENDIF
			ENDDO
		ENDDO
4		continue
	ENDDO

	DO J=1,3
	IF(CRITICAL_VAL%V(J)==1.D0)THEN
	I=0
	GOTO 5    ! not in the same cube, exit
	ENDIF
	ENDDO

5	continue
END SUBROUTINE IS_IN_SAME
!************************************************************************
!     FIND_MIDPOINTS: only used for case 3                              *
!************************************************************************
SUBROUTINE FIND_MIDPOINTS(POINT1,POINT2,POINT3,POINT4,TANG3,TANG4,   &
						     U1,U2,CRITICAL_VAL,MID_CONNECT_P,UUU,MNUM)
	use ZQMPI
	USE VECTORS
	USE VARIABLES,only:A_CUBE
	IMPLICIT NONE

	TYPE(VECTOR),intent(in)::POINT1,POINT2,POINT3,POINT4,TANG3,TANG4,CRITICAL_VAL
	DOUBLE PRECISION,intent(in)::U1,U2
	DOUBLE PRECISION,DIMENSION(20)::UUU
	TYPE(VECTOR),DIMENSION(20)::MID_CONNECT_P
	INTEGER::MNUM

	TYPE(VECTOR)::TMP_POINT
	INTEGER::I,J,K
	DOUBLE PRECISION::TMP_U1
	DOUBLE PRECISION,DIMENSION(4,3)::SHAP
	INTEGER::CONTROL_A,CONTROL_B
	integer::returnCheck

	MNUM=0

	DO I=1,3 ! FOR THREE COORDINATES
		IF(CRITICAL_VAL%V(I).EQ.1.D0)THEN ! HAS BOUNDARY VALUES
			J=INT(point1%v(I)/A_CUBE);K=INT(point2%v(i)/A_CUBE)
			if(point1%v(i)<0.d0)then
				J=j-1
			endif
			if(point2%v(i)<0.d0)then
				K=k-1
			endif

			CONTROL_A=MIN(J,K)
			CONTROL_B=MAX(J,K)

			DO J=CONTROL_A+1,CONTROL_B
				! FIND THE BOUNDARY VALUES SEPARATING THE TWO POINTS
				TMP_POINT%V(I)=J*A_CUBE

				! FROM THE BOUNDARY VALUE TO FIND THE U FOR THIS POINT
				if(MNUM>=20)exit
				CALL GET_U1(POINT3%V(I),POINT4%V(I),TANG3%V(I),TANG4%V(I),&
							U1,U2,TMP_POINT%V(I),UUU(MNUM+1),returncheck)

				if(returncheck==0)continue
				MNUM=MNUM+1
				CALL GETSHAPE(SHAP,UUU(MNUM))

				! CONSTRUCT THE BOUNDARY VALUE
				MID_CONNECT_P(MNUM)=SHAP(1,1)*POINT3+SHAP(2,1)*TANG3  &
									+SHAP(3,1)*POINT4+SHAP(4,1)*TANG4

			ENDDO
		ENDIF
	ENDDO

	!....rearrage the order of the points....
	DO I=1,MNUM
		TMP_U1=UUU(I)
		K=I

		DO J=I+1,MNUM
			!COMPARE Us AND PUT THE SMALLEST ONE TO THE FIRST POSITION
			IF(TMP_U1>UUU(J))THEN
				UUU(I)=UUU(J)
				UUU(J)=TMP_U1
				TMP_U1=UUU(I)

				TMP_POINT=MID_CONNECT_P(I)
				MID_CONNECT_P(I)=MID_CONNECT_P(J)
				MID_CONNECT_P(J)=TMP_POINT
			ENDIF
		ENDDO
	ENDDO
    
END SUBROUTINE FIND_MIDPOINTS
!***********************************************************************
!    map_piece:                                                        *
!***********************************************************************
subroutine map_piece(point1,point2)
	use ZQMPI
	use vectors
	use variables,only:A_CUBE
	use OUTPUTMD,only:tyPrePoint
	use FunctionMD,only:PointPBC
	implicit none

	integer::i,j
	type(vector),dimension(2)::piece_points
	type(vector)::point1,point2
	type(vector)::tmp_p
	integer,dimension(3)::cub_current


	!////////////////////////////////////////////////////////////
	!	This piece is already be determined to be in the same cube
	!	so we don't care about this problem.
	!	We only care how to map a whole piece in the same cube
	!	solve problem step by step     
	tmp_p=(point1+point2)/2.d0
	piece_points(1)=PointPBC(point1)
	piece_points(2)=PointPBC(point2)

	if(tyPrePoint%v(1)==1d10)then  !only for the beginning point of a loop
		write(11,*)"zone"
		write(11,*)piece_points(1)
		tyPrePoint=piece_points(1)
	endif

	if(MAG(tyPrePoint-piece_points(2))>0.5*A_Cube)then
		write(11,*)"zone"
	endif

	write(11,*)piece_points(2)
	tyPrePoint=piece_points(2)

end subroutine map_piece
!=======================================================================
subroutine get_U1(Point3,Point4,Tang3,Tang4,U01,U02,critical_val,U1,&
					returncheck)
	use ZQMPI
	use vectors
	implicit none

	double precision,intent(in)::Point3,Point4,Tang3,Tang4,critical_val !U1 IS OUTPUT VALUE
	double precision::fa,fb,fc,fd,fpp,fq
	double precision::U0,u1
	double precision,intent(in)::U01,u02 !U01 IS INPUT VALUE
	
	INTEGER::loopcon,returncheck

	loopcon=0
	returncheck=1

	U0=U01+0.5d0*dabs(U01-u02)

	fa=2.*Point3+Tang3-2.*Point4+Tang4
	fb=-3.*Point3-2.*Tang3+3.*Point4-Tang4
	fc=Tang3
	fd=Point3-critical_val

	!Newton method to get the solution of a higher order function.

	do 
		fpp=fa*U0**3+fb*U0**2+fc*U0+fd
		fq=3.*fa*U0**2+2.*fb*U0+fc

		if(dabs(fq)<1.0d-8)then
			returncheck=0
			exit
		endif
		U1=U0-fpp/fq
		loopcon=loopcon+1

		if(loopcon>=500)then
			returncheck=0
			exit
		endif

		if(dabs(U1-U0)>1.d-4)then
			U0=U1
		else
			returncheck=1
			exit
		end if
	enddo  

end subroutine get_U1
   
