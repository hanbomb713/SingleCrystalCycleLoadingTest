!********************************************************************
!	Microplasticity: 
!		Parallel dislocation dynamics simulation code for study of 
!	crystal plasticity. Written at LANL and UCLA by ZQ Wang.
!********************************************************************
PROGRAM MICROPLASTICITY
	use ZQMPI
	use InitialMD,only:SimulationInit
	use SimulationMD,only:Simulation
	IMPLICIT NONE

	!/////////////////////////////////////
	! Begin program
	!////////////////////////////////////

	!////////// MPI Initialization ////////
	call ZQMPI_Init
	call MPI_Barrier(MPI_COMM_WORLD,iIerr)
	
	!///////// Simulation Initialization /////////
	call SimulationInit
	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

	!///////// Main Simulation part //////////////
	call Simulation
	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

	!///////// Program is Ended. /////////////////
	if(iMyid==iMaster)then
	 write(*,*)"Program is ended successfully."
	endif

	!///////// MPI Finalization /////////////
	call ZQMPI_Finalization
	
		
END PROGRAM MICROPLASTICITY
