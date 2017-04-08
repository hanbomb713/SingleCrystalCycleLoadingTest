!*************************************************************************
!   Simulation: 
!		main calculation subroutine. 
!*************************************************************************
MODULE SimulationMD
	use ZQMPI
CONTAINS

!****************************************************
! Simulation main subroutine
!----------------------------------------------------
! Modified: ZQ Wang
!	02/16/06:
!	Removed some outputs
!****************************************************
subroutine Simulation
	use ZQMPI
	use variables,only:iBeginLoop,iLoop_time,iloop,iNPoint,iNplane,GPlanes
	use PerformanceMD
	implicit none

	double precision::dT1,dt2,dt3,dt4

	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

	dTotalTime = 0.d0
	dDynamicsTime = 0.d0
	dUpdateTime = 0.d0
	dBuildTreeTime = 0.d0
	
	! stress or strain update loop
	do iLoop = iBeginLoop+1, iBeginLoop+iLoop_time
				
		!....step output and step initialization.......
		if(iMyid==iMaster)then
			write(*,*)"----------------------------------------"
			write(*,*)iMyid,":ILOOP=",iloop,";Total DPs=",iNPoint,",N_Plane=",iNPlane
		endif
		call LoopStepInitial
		call MPI_Barrier(MPI_COMM_WORLD,iIerr)
		
		!....Test performance, beginning time.....
		dBeginTime = MPI_Wtime()

		!////////////////////////////////////////////////////////////////////
		!....split space and build global tree, distribute dislocations....
		call Split_Space_Build_Tree
		call MPI_Barrier(MPI_COMM_WORLD,iIerr)

		!------------------------------------------------
		dEndTime    = MPI_Wtime()
		dt1=dEndTime-dBeginTime
		dBuildTreeTime  = dBuildTreeTime+dt1
		dTotalTime  = dTotalTime+dt1
		dBeginTime = dEndTime

		Call MPI_Barrier(MPI_COMM_WORLD,iIerr)

		!------------------------------------------------
		dEndTime= MPI_WTime()
		dt2=dEndTime-dBeginTime
		dTotalTime = dTotalTime + dt2
		dBeginTime = dEndTime
					
		!///////////////////////////////////////////////////////////////////////
		!....every process has got enough information and begin calculation.....
		call Dynamics_Solver
		Call MPI_Barrier(MPI_COMM_WORLD,iIerr)
		!------------------------------------------------
		dEndTime = MPI_Wtime()
		dt3 = dEndTime-dBeginTime
		dDynamicsTime = dDynamicsTime + dt3
		dTotalTime = dTotalTime + dt3
		dBeginTime = dEndTime
		
		!//////////////////////////////////////////////////////////////////////
		!....update and output information....
		call Update_Output_Statistics
		call MPI_Barrier(MPI_COMM_WORLD,iIerr)
		!------------------------------------------------
		dEndTime = MPI_Wtime()
		dt4 = dEndTime - dBeginTime
		dUpdateTime = dUpdateTime + dt4
		dTotalTime = dTotalTime + dt4
		dBeginTime = dEndTime

		!/////////////////////////////////////////////////
		!....Output of time, dynamics and Update time....
		if(iMyid==iMaster)then
			write(*,*)iMyid,":Build Tree Time:",dt1,dBuildTreeTime/iloop
			write(*,*)iMyid,":Dynamics S Time:",dt3,dDynamicsTime/iLoop
			write(*,*)iMyid,":Update Pro Time:",dt4,dUpdateTime/iLoop
			write(*,*)iMyid,":Averaege l Time:",dTotalTime+pastTotalTime,&
			(dTotalTime+pastTotalTime)/iLoop
		end if		
		call MPI_Barrier(MPI_COMM_WORLD,iIerr)

	enddo

end subroutine Simulation
!*************************************************************************
!    Split_Space_BUild_Tree: 
!*************************************************************************
subroutine Split_Space_Build_Tree
	use ZQMPI
	use SplitbuildMD,only:Build_Global_Tree,Build_Local_Tree
	use variables,only:iloop
	implicit none

	!....Build global tree and distribute it. No need to build every step....
	call Build_Global_Tree

	!....Build loacal tree, calculate multipole information, if updated, only calculate new....
	!....multipole information....
	call Build_Local_Tree
	
end subroutine Split_Space_Build_Tree
!******************************************************************************
!   LoopStepInitial: 
!			Each loop step after the space is resplit, the local variables, such
!           as strain, should be reset to be zero. Stress keeps the same.
!           Master keeps the total strain.
!
!			iNTotalDPs is always the same as total number of DPs on master.
!******************************************************************************
subroutine LoopStepInitial
	use ZQMPI
	use Variables,only: zero,StrainCurrent,StrainField,iNPoint,iNPlane
	use MPIMOdule,only: iNTotalDPs,iNMaxS,iNumTransP
	use CommunicationMD,only:BcastStepVariables
	use AnnihilationMD,only:iAnniCount,dpAnniLength
	use CrossslipMD,only:iCSCount,dpCSLength,iNumTotalScrew,iNumTotalScrew1
	implicit none

	if(iMyid/=iMaster)then
		StrainCurrent	=	0.d0
		StrainField		=	Zero
		iNumTransP		=	0
		iNPOint			=	0
		iNPlane			=	0
	else
		iAnniCount=0
		dpAnniLength=0
		iCSCount=0
		iNumTotalScrew=0
		iNumTotalScrew1=0
		dpCSLength=0
		iNTotalDPs	= iNPoint
		iNMaxS		= iNPoint/(iNumProcs-1)
	endif

	call BcastStepVariables

end subroutine LoopStepInitial
!*************************************************************************
!   Dynamics_Solver: 
!		Do loop for each dislocation on current processor.
!       done by slaves
!*************************************************************************
subroutine Dynamics_Solver
	use ZQMPI
	use DynamicsSolverMD,only:DynamicsMain
	implicit none

	call DynamicsMain

end subroutine Dynamics_Solver

!*************************************************************************
!   Update_Output_Statistics:
!*************************************************************************
subroutine Update_Output_Statistics
	use ZQMPI
	use UpdateMD
	implicit none

	!....Dislocation in the process are automatically update due to the calculation.
	!....Update the subsegment in each process, including information exchange.
	!....No multipole calculation is performed in this subroutine but leave it to the 
	!....subroutine of Split_Space_Build_Tree. 
	
	!....update strain, stress
	call UpdateMechanics
	call MPI_Barrier(MPI_COMM_WORLD,iIerr)
	
	! Statistics, perform statistics in each process and then gather them to the master
	! and output information to files.
	call Statistics
	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

	! Update dislocations in master
	call UpdateMasterDis
	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

	! Remove allocate variables
	call RemoveVariables
	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

	!/////////////////////////////////////
	!Short range interactions
	call ShortRangeInter
	call MPI_BArrier(MPI_COMM_WORLD,iIerr)
	
	!////////////////////////////////////
	! Ouput, output is done by master.
	call Output
	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

	!////////////////////////////////////
	! Adjust time step, N_TIMES
!	call AdjustTimeStep
!	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

end subroutine Update_Output_Statistics

!*************************************************************************

END MODULE SimulationMD
  
