!*******************************************************************
!   This file contains subroutines that initialize the simulation,
!   including data read, memory allocation, global variable 
!	distribution etc.
!
!	After initialization, all processors have global
!	variables. Master processor has dislocations and next step is 
!	to decompose the compuational domain and distribute computational
!	load to individual slave processors.
!************************************************************
module InitialMD
	use zqmpi

	contains
!************************************************************
!  SimulationInit:                                          
!************************************************************
subroutine SimulationInit
	use ZQMPI
	implicit none

	!....default values....
	call Defaults
	
	!....Master read data....
	if(iMyid==iMaster)then
		call Open_IO_files
		
		call DATA_Reader
	endif

	call MPI_Barrier(MPI_COMM_WORLD,iIerr)
	!....Global variables.....
	call BCastAllGVs

	!....Gauss Points....
	call quadrature

	!....Simulation variables initialization.....
	call SimulationInitial

	call MPI_Barrier(MPI_COMM_WORLD,iIerr)
	
end subroutine SimulationInit
!************************************************************
!   DEFAULTS:  
!		Set default values of variables.                 
!		Both master and slaves.
!---------------------------------------------
! Modified: ZQ Wang, 1/25/06
!           set iInertial and iObstacle.
!	01/28/06:
!	Added velocity count initialization
!	02/16/06:
!	Added output control parameters.
!	02/26/06:
!	Increased the dimension for velocity statistics.
!	04/08/06:
!	Added variable iCrossslip
!************************************************************
SUBROUTINE DEFAULTS
	use ZQMPI
	USE VARIABLES
	USE VECTORS
	use MPIModule,only:iNTotalGPs
	use OutputMD,only:iDisOutFreq,iGeneralOutFreq,iPastFreq
	IMPLICIT NONE

	integer::i,j
	double precision::dpTmp1
	!////////////////////////////
	do i=1,3
		do j=1,3
			D%v(i)%v(j)=int(i/j)*int(j/i)
			ZERO%V(I)%V(J)=0.0D0
		end do
	end do

	PI	          =  2.*DASIN(1.D0)
	INTERVAL        =  0.5D0
	INTERVAL2       =  INTERVAL*INTERVAL
	RMIN            =  DSQRT(2.0D0)/2.d0

	!material
	iMatType		=iFCC
	MU		      =  50D9
	NU		      =  0.3
	LATTICE 	  =  3.615D-10
	CycFreq		  =	 50
	dpCycF0		  =  50d6
	dpCycStrain	  =   1.2D-3
	dpLSoundV	 =	1.1D13
	dpTSoundV	 =	6.82D12
	dpEMass0	=	MU*0.5/dpTSoundV**2   !0.5 is the square of the burgers vector
	PEIERLS		=	0
	MOB_R		=	1
	PEI_R		=	1
	ATRatio		=	1.25
	TRatio		=	0.25
	NGRatio		=	0.5

	!speed intervals,1/28/06
	dpTmp1=dpTSoundV/100.0
	do i=1,100
		dpVSeg(i)=i*dpTmp1
		dpVCountL(i)=0.d0
		iVCountN(i)=0
		dpVCountN(i)=0.d0
	enddo

	!simulation 1
	A_CUBE  	  =  10000.0D0
	SIG_APP       =  ZERO
	StrainField	=zero
	StrainCurrent	=0	!the plastic strain
	iNTotalGPs	=0


	!simulation 2
	MAX_QUAD	  =  16
	MAX_NODE	  =  20
	MAX_PLANE	  =  100
	MAX_NEIGHBOR  =  50
	DELTA_SIG     =  10D6
	DTIME         =   1.0D-14
	iArr		  =	5
	iDisOutFreq	=50	
	iGeneralOutFreq	= 50
	iPastFreq	= 50
	iObstacle	=	0
	iInertial	=	0
	iCrossslip	=	0
	!simulation 3
	IntegrateMethod =  1
	N_TIMES			=  100
	NPOINT_I		=  10
	dcrit           =  100.d0
	logArrangeLoop	=	1
	logAnnihilation	=	1
	Load			=	1
	Loading_Dir%v(1)	=1
	Loading_Dir%v(2)	=0
	Loading_Dir%v(3)	=0

END SUBROUTINE DEFAULTS

!****************************************************
!	Initialize Slip Systems
!	Slip systems for a specific material defined by iMatType
!	Put it together with Allocate_dimension
!****************************************************
subroutine InitialSS	
	USE ZQMPI
	use vectors
	use variables,only:iMatType,Burgers,iNumBurgers,iNumMiller,Miller,D,iFCC,iBCC
	implicit none
	
	integer::i,j
	
	!/////////////////////////////////////
	!Burgers
	if(iMatType==iFCC)then
		!1/2{110} direction
		iNumBurgers=12
		allocate(Burgers(iNumBurgers))
		Burgers(1)%v(1)=0.d0
		Burgers(1)%v(2)=0.5d0
		Burgers(1)%v(3)=0.5d0
	
		Burgers(2)%v(1)=0.5d0
		Burgers(2)%v(2)=0.d0
		Burgers(2)%v(3)=0.5d0
	
		Burgers(3)%v(1)=0.5d0
		Burgers(3)%v(2)=0.5d0
		Burgers(3)%v(3)=0.d0
	
		Burgers(4)%v(1)=0.d0
		Burgers(4)%v(2)=-0.5d0
		Burgers(4)%v(3)=0.5d0
	
		Burgers(5)%v(1)=-0.5d0
		Burgers(5)%v(2)=0.d0
		Burgers(5)%v(3)=0.5d0
	
		Burgers(6)%v(1)=-0.5d0
		Burgers(6)%v(2)=0.5d0
		Burgers(6)%v(3)=0.d0

		do i=7,12
			j=i-6
			Burgers(i)%v(1)=-1d0*Burgers(j)%v(1)
			Burgers(i)%v(2)=-1d0*Burgers(j)%v(2)
			Burgers(i)%v(3)=-1d0*Burgers(j)%v(3)
		enddo
	elseif(iMatType==iBCC)then
		!1/2{111} direction
		iNumBurgers=8
		allocate(Burgers(iNumBurgers))
		Burgers(1)%v(1)=0.5d0
		Burgers(1)%v(2)=0.5d0
		Burgers(1)%v(3)=0.5d0
	
		Burgers(2)%v(1)=0.5d0
		Burgers(2)%v(2)=0.5d0
		Burgers(2)%v(3)=-0.5d0
	
		Burgers(3)%v(1)=0.5d0
		Burgers(3)%v(2)=-0.5d0
		Burgers(3)%v(3)=-0.5d0
	
		Burgers(4)%v(1)=0.5d0
		Burgers(4)%v(2)=-0.5d0
		Burgers(4)%v(3)=0.5d0

		do i=5,8
			j=i-4
			Burgers(i)%v(1)=-1d0*Burgers(j)%v(1)
			Burgers(i)%v(2)=-1d0*Burgers(j)%v(2)
			Burgers(i)%v(3)=-1d0*Burgers(j)%v(3)
		enddo
	endif

	!////////////////////////////////
	!Miller
	if(iMatType==iFCC)then

		!{111} planes
		iNumMiller=4
		allocate(Miller(iNumMiller))
		Miller(1)%v%v(1)=1d0
		Miller(1)%v%v(2)=1d0
		Miller(1)%v%v(3)=1d0
	
		Miller(2)%v%v(1)=1d0
		Miller(2)%v%v(2)=-1d0
		Miller(2)%v%v(3)=1d0
	
		Miller(3)%v%v(1)=-1d0
		Miller(3)%v%v(2)=-1d0
		Miller(3)%v%v(3)=1d0
	
		Miller(4)%v%v(1)=-1d0
		Miller(4)%v%v(2)=1d0
		Miller(4)%v%v(3)=1d0

	elseif(iMatType==iBCC)then
		
		!{110} & {112} planes
		iNumMiller=18
		allocate(Miller(iNumMiller))
		Miller(1)%v%v(1)=1d0
		Miller(1)%v%v(2)=1d0
		Miller(1)%v%v(3)=0d0
	
		Miller(2)%v%v(1)=1d0
		Miller(2)%v%v(2)=-1d0
		Miller(2)%v%v(3)=0d0
	
		Miller(3)%v%v(1)=1d0
		Miller(3)%v%v(2)=0d0
		Miller(3)%v%v(3)=1d0
	
		Miller(4)%v%v(1)=1d0
		Miller(4)%v%v(2)=0d0
		Miller(4)%v%v(3)=-1d0

		Miller(5)%v%v(1)=0d0
		Miller(5)%v%v(2)=1d0
		Miller(5)%v%v(3)=1d0

		Miller(6)%v%v(1)=0d0
		Miller(6)%v%v(2)=1d0
		Miller(6)%v%v(3)=-1d0
	
		!<112>
		Miller(7)%v%v(1)=-1d0
		Miller(7)%v%v(2)=-1d0
		Miller(7)%v%v(3)=-2d0
	
		Miller(8)%v%v(1)=1d0
		Miller(8)%v%v(2)=1d0
		Miller(8)%v%v(3)=-2d0
	
		Miller(9)%v%v(1)=1d0
		Miller(9)%v%v(2)=-1d0
		Miller(9)%v%v(3)=-2d0

		Miller(10)%v%v(1)=-1d0
		Miller(10)%v%v(2)=1d0
		Miller(10)%v%v(3)=-2d0
	
		Miller(11)%v%v(1)=1d0
		Miller(11)%v%v(2)=2d0
		Miller(11)%v%v(3)=1d0
	
		Miller(12)%v%v(1)=-1d0
		Miller(12)%v%v(2)=-2d0
		Miller(12)%v%v(3)=1d0
	
		Miller(13)%v%v(1)=-1d0
		Miller(13)%v%v(2)=2d0
		Miller(13)%v%v(3)=1d0

		Miller(14)%v%v(1)=1d0
		Miller(14)%v%v(2)=-2d0
		Miller(14)%v%v(3)=1d0
	
		Miller(15)%v%v(1)=2d0
		Miller(15)%v%v(2)=1d0
		Miller(15)%v%v(3)=1d0
	
		Miller(16)%v%v(1)=2d0
		Miller(16)%v%v(2)=-1d0
		Miller(16)%v%v(3)=1d0
	
		Miller(17)%v%v(1)=-2d0
		Miller(17)%v%v(2)=1d0
		Miller(17)%v%v(3)=1d0

		Miller(18)%v%v(1)=-2d0
		Miller(18)%v%v(2)=-1d0
		Miller(18)%v%v(3)=1d0

		!

	endif
	
	!....compute ES for each Miller direction....
	do i=1,iNumMiller
		Miller(i)%ES%v(3)=UV(Miller(i)%v)
		Miller(i)%ES%v(1)=UV(d%v(3)**Miller(i)%v)
		Miller(i)%ES%v(2)=Miller(i)%ES%v(3)**Miller(i)%ES%v(1)
	enddo

END SUBROUTINE InitialSS
!***************************************************************
!	BCastAllGVs: 
!		Both master and slaves
!		Broadcast and alloate global variables.
!	
!		iaGLIndex will not be deallocated for each simulation step. 
!		But, after each time the tree is rebuilt, its content will 
!		be refreshed. It is same for iaTransList.
!***************************************************************
subroutine BCastAllGvs
	use ZQMPI
	use CommunicationMD,only:bcastSomeVariables
	implicit none

	call bcastSomeVariables     ! in Comuncation file

	! Allocate dimensions
	if(iMyid/=iMaster)THEN
		CALL InitialSS
		CALL ALLOCATE_DIMENSIONS
	endif

end subroutine BCAstAllGVs
!*************************************************************************
!    SimulationInitial:   
!		Initialize some variables needed by simulation.
!		Both master and slaves.
!*************************************************************************
subroutine SimulationInitial
	use ZQMPI
	implicit none

	call TransLoadingDir

end subroutine SimulationInitial

!***************************************************************
!   OPEN_IO_FILES:                                             *
!		Only on master.
!	Initial Input: single digit;   (material, geom_input, iconornew)
!	Geometry output: from 10 to 19;
!	Other output: from 20 to 39;
!	Past intermediate files: 111, 112, 113, 114, etc;
!
!---------------------------------------------------------------
! Modified: ZQ Wang, 01/28/06
!	Added new file 28.
!	02/14/06:
!	removed old 25 and 27, change 26 and 28 to 25 and 27, defined
!	new rules above.
!***************************************************************
SUBROUTINE OPEN_IO_FILES
	use ZQMPI
	implicit none

	!Input files
	OPEN(UNIT=1,FILE="material_input.txt",action="read")
	OPEN(UNIT=2,FILE="interaction_geom_input.txt",action="read")
	open(unit=3,file="iconornew.txt",action="read")

	!Output files
	OPEN(UNIT=11,FILE="fr_source_g.dat")

	open(unit=20,file="strain_stress.dat")
	open(unit=21,file="time_stress.dat")
	open(unit=22,file="time_strain.dat")
	open(unit=23,file="strain_afreq.dat")
	open(unit=24,file="strain_density.dat")
	open(unit=25,file="VeoL.dat")
	open(unit=26,file="VeoN.dat")
	open(unit=27,file="stress_density.dat")
	open(unit=28,file="strain_cs.dat")
	open(unit=30,file="strain_csl.dat")
	open(unit=31,file="strain_al.dat")
	open(unit=29,file="strain_loopnum.dat")
	open(unit=32,file="strain_csn.dat")
	open(unit=33,file="strain_screwd.dat")
	open(unit=34,file="strain_edgemixd.dat")
END SUBROUTINE OPEN_IO_FILES
!******************************************************************
!   DATA_READER:  
!		Read material constants, dislocation positions, and simulation
!		control data.                               
!		Only on master.
!-----------------------------------------------------------
!  Modified: ZQ Wang
!	02/16/06:
!	Rearranged the file input sequence.
!	04/08/06:
!	added variable iCrossslip
!******************************************************************
SUBROUTINE DATA_READER
	use ZQMPI
	USE VECTORS
	USE VARIABLES
	use MPIModule,only:iNTotalDPs
	use OutputMD,only:iDisOutFreq,iGeneralOutFreq,iPastFreq
	IMPLICIT NONE

	integer::I,iconOrNew
	double precision::ld1,ld2,ld3

	integer::iEnd,iBegin
	integer,dimension(5000)::iCount
	logical::lOpen

	NAMELIST /MATERIAL/iMatType,MU,NU,LATTICE,A_CUBE,APPLIED_SIG,RMIN,DELTA_SIG,MOBILITY,DTIME,PEI_R,	&
			   ElasticCons,StrainRate,dcrit,CycFreq,dpCycF0,dpCycStrain,ld1,ld2,ld3,PEIERLS,MOB_R,ATRatio,TRatio, &
				NGRatio

	NAMELIST /DIMENSIONS/MAX_QUAD,MAX_NODE, MAX_PLANE,MAX_NEIGHBOR,           &
	                     N_TIMES,ILOOP_TIME,NPOint_I,iArr,iObstacle,iInertial,iCrossslip,iDisOutFreq, &
				IntegrateMethod,logAnnihilation, logArrangeLoop,Load,CheckNeiBur,iPastFreq,iGeneralOutFreq

	!........Initial loading direction.......
	ld1=0
	ld2=0
	ld3=0

	!//////////////////////////////////////////////////
	!Read data in

	!....Read Simulation control data....
	READ(1,MATERIAL) 
	READ(1,DIMENSIONS)
        read(3,*)iconOrnew

	call InitialSS
	call allocate_dimensions

	!....import dislocation data....
	iNPlane=0
	iNPoint=0

	if(iconOrnew==1)then       ! new simulation
		
		CALL DislocationReader
		iBeginLoop=0
		totalStrain=0.d0
		totalTime=0.d0
	else                         ! continued simulation
	
		call Read_Past
	
	endif

	!///////////////////////////////////////////
	!Some preliminary calculations
    	!....compute PG,TG for each dislocation point....
    	Do i=1,iNPoint
        	call SetGlobal(i)
    	enddo

	!....Total number of npoints........	
	iNTotalDPs    =iNPoint               ! different with iNPoint on local processors, same on master
	N_TIMESCOPY   =N_TIMES

	!--------Loading conditions---------------
	Loading_Dir%v(1)=ld1;Loading_Dir%v(2)=ld2;Loading_Dir%v(3)=ld3
	if(Mag(Loading_Dir)==0)Loading_Dir%v(1)=1
	Loading_Dir=UV(Loading_Dir)


	!.....
	PEIERLS=PEIERLS/MU

	close(1)
	close(2)
	close(3)

END SUBROUTINE DATA_READER
!****************************************************************
!   Allocate_dimensions:  
!		Allocate memories for arrays.         
!****************************************************************
SUBROUTINE ALLOCATE_DIMENSIONS
	use ZQMPI
	USE VECTORS
	USE VARIABLES,only:GPLANES,DPoints,MAx_Plane,MAX_Node, &
						MAX_Neighbor,zero,dpEMass0
	use MPIModule,only:iaGLIndex,iaTransList,tpBoxRoots,iIndexOfCloseBox
	use DynamicsMD,only:Dmatrix
	use LoopRearrangeMD,only:tmpDP
	implicit none

	integer::i

	! for each glide plane
	Allocate(GPlanes(int(3*MAX_Plane))) ! one point may have 3 glide planes if it is mapped using PBC.
	do i=1,3*MAX_PLANE
		GPlanes(i)%iNumPEquiv=0
	enddo

	! for each DPs
	ALLOCATE(DPoints(int(MAX_PLANE*MAX_NODE)))
	do i=1,Max_Plane*Max_Node
		DPoints(i)%iNumNei=0
		DPoints(i)%tvAcce=zero%v(1)
		DPoints(i)%tvPreV=zero%v(1)
		DPoints(i)%tvPreVT=zero%v(1)
		DPoints(i)%dpMass=dpEMass0
		DPoints(i)%lToCal=.true.
		allocate(DPoints(i)%iaIDOfNei(MAX_NEIGHBOR))
	enddo

	allocate(iaGLIndex(int(MAX_PLANE*MAX_NODE),4))

	allocate(tpBoxRoots(iNumProcs),iIndexOfCloseBox(iNumProcs))

	!iaTranslist
	if(iMyid/=iMaster)then
		allocate(iaTransList(int(Max_Plane*Max_Node)))
		allocate(Dmatrix(MAX_Plane*MAX_NODE))
		do i=1,MAX_PLANE*MAX_NODE
			allocate(Dmatrix(i)%FG(4))
			allocate(Dmatrix(i)%KG(4,4))
			allocate(Dmatrix(i)%QG(4,8))
		enddo
	endif

	! for tmpdp
	if(iMyid==iMaster)then
		allocate(tmpDP(1)%iaIDOfNei(MAX_Neighbor),tmpDP(2)%iaIDOfNei(MAX_Neighbor))
	endif

END SUBROUTINE ALLOCATE_DIMENSIONS
!****************************************************************
!   DislocationReader:  
!	Read dislocation positions from files.
!****************************************************************
SUBROUTINE DislocationReader
	use ZQMPI
	use vectors
	implicit none

	CHARACTER ::DEFECT_KIND*20, ENDDEFECT*9
	
	integer::i,j,k
	INTEGER :: PLANEID, BurID,NODES,Plane_Number,number_frs
	double precision :: circ_radius
	Type(vector):: BURGERL,plane_miller,ORIGINL,circle_cntr,P_Local(2)

	number_frs=0
	DO WHILE(DEFECT_KIND .NE. 'DEFECT' )
		READ(2,*) DEFECT_KIND
		IF(DEFECT_KIND .EQ. 'DEFECT')THEN
			DO WHILE (DEFECT_KIND .NE. 'ENDDEFECT')

				READ(2,*) DEFECT_KIND
				IF (DEFECT_KIND .EQ. 'FRS')     GOTO 1
				IF (DEFECT_KIND .EQ. 'ENDDEFECT') GOTO 2

1				continue

			    IF(DEFECT_KIND .EQ. 'frs' ) THEN
					READ(2,*)DEFECT_KIND, PLANE_NUMBER, NODES,  &
			                 plane_miller, BURGERL, ORIGINL,  &
							 circ_radius,circle_cntr%v(1), circle_cntr%v(2)
					call DetermineGPlanes(OriginL,plane_miller,PlaneID)
					
					call DetermineBurgers(BurgerL,BurID)
					
					DO k = 1, 2
						READ (2,*) J, P_local(k)%v(1),P_local(k)%v(2)
					END DO
					number_frs=number_frs+1
					call AddNewDislocation(PlaneID,BurID,number_frs,P_local)

				END IF
			enddo 
		endif
	enddo

2	CONTINUE

END SUBROUTINE DislocationReader
!**************************************************************
! DetermineBurgers:
!	Determine the index of the burgers vectors, return index.
!**************************************************************
Subroutine DetermineBurgers(Bur,iIndex)
	use ZQMPI
	use vectors
	use variables,only:Burgers,iNumBurgers
	implicit none

	type(vector),intent(in)::Bur
	integer::iIndex
	integer::i,j
	
	iIndex=1

	do i=1,iNumBurgers
		if(MAG(Burgers(i)-Bur)==0.d0)then
			iIndex=i
			exit
		endif
	enddo

end subroutine DetermineBurgers
!********************************************************************
!	AddNewDislocation:
!		Add new dislocation point to the arrays.
!********************************************************************
subroutine AddNewDislocation(PlaneID,BurID,number_frs,PLocal)
	use ZQMPI
	use vectors
	use variables,only:DPoints,iNPoint,NPoint_I
	implicit none

	integer,intent(in)::PlaneID,BurID,number_frs
	type(vector),intent(in)::PLocal(2)

	integer::i,j,k,nPoint
	type(vector),dimension(NPOint_I)::PL,TL
	type(vector)::tmpP,tmpT

	nPoint=NPoint_I
	PL(1)=PLocal(1)
	PL(2)=PLocal(2)
	call setNewInitDis(NPoint,PL,TL)

	!....add new points to the system....
	call AddNewDisPoint(PlaneID,BurID,number_frs,NPOint,PL,TL)

end subroutine AddNewDislocation
!*********************************************************************
!	AddNewDisPoint:
!		Add a new dislocation point to the array.
!*********************************************************************
subroutine AddNewDisPoint(PlaneID,BurID,number_frs,NPoint,PL,TL)
	use ZQMPI
	use vectors
	use variables,only:DPoints,iNPoint,iDPTypeFree,iDPTypeFixed
	implicit none

	integer,intent(in)::PlaneID,BurID,number_frs,NPOINT
	type(vector),intent(in),dimension(NPOINT)::PL,TL

	integer::i,j,k

	do i=1,NPOINT
		iNPOint=iNPoint+1
		j=iNPoint
		
		!...assign values....
		DPoints(j)%iID=iNPoint
		DPoints(j)%iLoopID=number_frs
		if(i==1 .or. i==NPOINT)then
			DPoints(j)%iType=iDPTypeFixed
		else
			DPoints(j)%iType=iDPtypeFree
		endif

		DPoints(j)%tvPL=PL(i)
		DPoints(j)%tvTL=TL(i)

		DPoints(j)%iPlane=PlaneID
		DPoints(j)%iBurgers=BurID

		DPoints(j)%iBeginP=iNPoint-1
		DPoints(j)%iEndP=iNPoint+1
		if(i==1)DPoints(j)%iBeginP=-1     ! beginning is fixed
		if(i==NPoint)DPoints(j)%iEndP=-1  ! ending is fixed
	enddo
	
end subroutine AddNewDisPoint
!************************************************************************
! TransLoadingDir(ld1,ld2,ld3)
!
!    If loading is not in the default direction, then transfer it to the 
!global coordinate system.
!************************************************************************
subroutine TransLoadingDir
	use ZQMPI
	use vectors
	use variables,only:loading_dir,Loading_ES,APPLIED_SIG,MU,D,SIG_APP,totalStress,zero
	implicit none

	integer::i,j,l,m
	double precision::ld1,ld2,ld3
	type(vector)::tmpzV
	type(matrix)::tmpSig,trans_es

	totalStress=Applied_Sig/MU

	! Calcualte the transforamtion matrix
	call ComputeES(loading_ES,loading_Dir)

	!Calculate stress state in global coordinate system
	tmpSig%v(1)%v(1)=totalStress
	trans_es=Trans(Loading_ES)
	SIG_APP=zero

	do i=1,3;do j=1,3
		do l=1,3;do m=1,3
			Sig_APP%v(i)%v(j)=SIG_APP%v(i)%v(j)+Trans_ES%v(i)%v(l)*Trans_ES%v(j)%v(m)*tmpSig%v(l)%v(m)
		enddo;enddo
	enddo;enddo

end subroutine TransLoadingDir
!********************************************************************
! ComputeES(a,b)
!********************************************************************
subroutine ComputeES(a,b)
	use ZQMPI
	use vectors
	use variables,only:zero,D,PI
	implicit none

	type(matrix)::a
	type(vector),intent(in)::b
	type(vector)::c,tmpZV
	type(vector)::v1,v2
	double precision::dpAngle

	c=UV(b)

	if(c%v(3)==0.0)then
		v1%v(1)=1.0;v1%v(2)=0.d0;v1%v(3)=0.0
		v2%v(1)=0.0;v2%v(2)=1.d0;v2%v(3)=0.0
		dpAngle=dacos(c*v1)
		if(v2*c<0)dpAngle=2*PI-dpAngle

		a%v(1)%v(1)=dcos(dpAngle);a%v(1)%v(2)=dsin(dpAngle);a%v(1)%v(3)=0.d0
		a%v(2)%v(1)=-dsin(dpAngle);a%v(2)%v(2)=dcos(dpAngle);a%v(2)%v(3)=0.d0
		a%v(3)%v(1)=0.d0;a%v(3)%v(2)=0.d0;a%v(3)%v(3)=1.d0
	else

		! the local y direction is the real z direction
		a%v(3)=UV(c)
		a%v(1)=UV(d%v(3)**c)
		a%v(2)=a%v(3)**a%v(1)

		!
		tmpZV%v(1)=0;tmpZV%v(2)=1;tmpZV%v(3)=0     !local y
		c=Trans(a)*tmpZV              !global z

		!
		a%v(3)=UV(c)
		a%v(1)=UV(d%v(3)**c)
		a%v(2)=a%v(3)**a%v(1)
	endif

end subroutine ComputeES
!****************************************************************
end module InitialMD
