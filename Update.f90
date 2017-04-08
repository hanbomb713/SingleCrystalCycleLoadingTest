!Update and Output and Statistics
!***********************************************************************
!   Update: update strain, stress;
!           update subsegments, update local dislocations. Only positions
!**************************************************************************
Module UpdateMD
	use ZQMPI

CONTAINS
!**************************************************************************
!   UpdateMechanics:   update strain, stress. Master gather these data and
!                      brcast new values to each process.
!**************************************************************************
subroutine UpdateMechanics
	use ZQMPI
	use vectors
	use variables,only:StrainField,StrainCurrent,DeltaStrain,&
						StrainRate_p,DTime,N_Times,totalTime,loading_ES,zero
	use CommunicationMD,only : UpdateStrainStress1
	implicit none

	integer::i,j,l,m
	type(matrix)::tmpStrain


	totalTime=totalTime+N_TIMES*DTIME
	tmpStrain=zero
	if(iMyid/=iMaster)then    ! slaves do

		!Strain for the applied direction, calculated from global coordinate system
		do i=1,3;do j=1,3;
			do l=1,3;do m=1,3
				tmpStrain%v(i)%v(j)=tmpStrain%v(i)%v(j)+loading_ES%v(i)%v(l)*loading_ES%v(j)%v(m)*strainField%v(l)%v(m)
			enddo;enddo
		enddo;enddo

		StrainCurrent = tmpStrain%v(1)%v(1)        !plastic strain
		DeltaStrain = StrainCurrent                !plastic delta strain
		StrainRate_p = DeltaStrain/(DTIME*N_TIMES) !plastic strain rate,for each iout, so no iout needed here.
	else
		StrainRate_p = 0.d0
		DeltaStrain =0.d0
	endif

	call UpdateStrainStress1

end subroutine UpdateMechanics
!**************************************************************************
!   UpdateMechanics:   update strain, stress. Master gather these data and
!                      brcast new values to each process.
!	02/06/2009: to correct stress calculation
!**************************************************************************
subroutine UpdateMechanics1
        use ZQMPI
        use vectors
        use variables,only:StrainField,StrainCurrent,DeltaStrain,&
                                                StrainRate_p,DTime,N_Times,totalTime,loading_ES,zero
        use CommunicationMD,only : UpdateStrainStress
        implicit none

        integer::i,j,l,m
        type(matrix)::tmpStrain


        totalTime=totalTime+N_TIMES*DTIME
        tmpStrain=zero
        if(iMyid/=iMaster)then    ! slaves do

                !Strain for the applied direction, calculated from global coordinate system
                do i=1,3;do j=1,3;
                        do l=1,3;do m=1,3
                          tmpStrain%v(i)%v(j)=tmpStrain%v(i)%v(j)+loading_ES%v(i)%v(l)*loading_ES%v(j)%v(m)*strainField%v(l)%v(m)
                        enddo;enddo
                enddo;enddo
        endif

        call UpdateStrainStress(tmpStrain)

end subroutine UpdateMechanics1
!**************************************************************************
!   Statistics: first do it locally, then master will gather them together.
!					1. Density, average distance
!					2. Time
!**************************************************************************
subroutine Statistics
	use ZQMPI
	use vectors
	use MPIModule,only:dTotalLength,dTmpDensi,dDensity,dScrewDensity,dEdgeMixDensity,dTmpSDensi,dTmpEDensi
	use variables,only:A_cube,Lattice,iNPoint,DPoints,Miller,GPlanes
	use CommunicationMD,only:AddDensityTogether
	use FunctionMD,only:funcIsScrew
	implicit none

	integer::I,j,iCurrentD,iNext
	type(vector)::s1,s2,RR,PL(2),TL(2),tt
	double precision::dUU0
	double precision,dimension(4,3)::shap

	!/////////////////////////////
	! desity:
	dTotalLength=0.d0
	dTmpDensi=0.d0
	dDensity=0.d0
	dScrewDensity=0.d0
	dTmpSDensi=0.d0
	dEdgeMixDensity=0.d0
	dTmpEDensi=0.d0
	if(iMyid/=iMaster)then
		DO iCurrentD=1,iNPoint
			iNext=DPoints(iCurrentD)%iEndP
			if(iNext/=-1 .and. DPoints(iCurrentD)%iPlane==DPoints(iNext)%iPlane)then
				PL(1)=DPoints(iCurrentD)%tvPL
				TL(1)=DPoints(iCurrentD)%tvTL
				call TransNextLocal(iCurrentD,PL(2),TL(2))
				DO I=1,11
					dUU0=0.1*(I-1)
					CALL GETSHAPE(SHAP,dUU0)
					RR=SHAP(1,1)*PL(1)+shap(2,1)*TL(1)+ &
						shap(3,1)*PL(2)+shap(4,1)*TL(2)
					tt=SHAP(1,2)*PL(1)+SHAP(2,2)*TL(1)+ &
						shap(3,2)*PL(2)+shap(4,2)*TL(2)
					tt=Miller(Gplanes(DPoints(iCurrentD)%iPlane)%iMiller)%ES*tt
					if(I==1)then
						s1=RR
					else
						s2=RR
						dTotalLength=dTotalLength+Mag(S2-S1)
						if(funcIsScrew(iCurrentD,tt))then
							dTmpSDensi=dTmpSDensi+MAG(S2-S1)
						else
							dTmpEDensi=dTmpEDensi+MAG(S2-S1)
						endif
						s1=s2
					endif
				enddo
			endif
		enddo

		dTmpDensi=dTotalLength/(A_Cube**3)/((Lattice*100)**2)   ! in unit of cm
		dTmpSDensi=dTmpSDenSi/(A_Cube**3)/((Lattice*100)**2)   ! in unit of cm
		dTmpEDensi=dTmpEDenSi/(A_Cube**3)/((Lattice*100)**2)   ! in unit of cm

	endif

	call AddDensityTogether

	call VelocityStat

end subroutine Statistics
!***************************************************************************
!	VelocityStat:
!		Count the number of nodes and calculate length for different
!	velocity distributions.
!-----------------------------------------------------------------------------
! Modified: ZQ Wang, 01/28/06
!	New way to count.
!	02/26/06
!	Increased dimension for velocity distribution.
!***************************************************************************
subroutine VelocityStat
	use ZQMPI
	use Vectors
	use Variables,only:dpVseg,dpVCountL,iVCountN,iNPoint,DPoints
	use DynamicsMD,only:DMatrix
	use CommunicationMD,only:CombineVelocityNumber
	implicit none

	integer::i,j,k
	double precision::dpTmp1,dpTmp2
	double precision,dimension(4,3)::shap
	double precision::dpFComputeLength
	
	!//////////////////////////////////////////////
	do i=1,100
		dpVCountL(i)=0.d0
		iVCountN(i)=0
	enddo
	dpTmp1=0.d0
	dpTmp2=0.d0


	!////////////////////////////////////////
	!	Count the number of nodes at different speeds
	!	Calculate the length of dislocations at different speeds.
	if(iMyid/=iMaster)then
		do i=1,iNPoint
			dpTmp1=MAG(DMatrix(i)%tvPLRate)
			!////////////////////////////////////
			! number of nodes
			j=dpTmp1/dpVSeg(1)+1      ! check which speed slot should the point be put in.

			if(j>100)then
				iVCountN(100)=iVCountN(100)+1
			elseif(j<1)then
				iVCountN(1)=iVCountN(1)+1
			else
				iVCountN(j)=iVCountN(j)+1
			endif

			!////////////////////////////////////
			! length of dislocations
			if(DPoints(i)%iEndP/=-1)then
				
				dpTmp2=MAG(DMatrix(DPoints(i)%iEndP)%tvPLRate)
				dpTmp1=0.5*(dpTmp1+dpTmp2)
				j=dpTmp1/dpVseg(1)+1

				if(j>100)then
					dpVCountL(100)=dpVCountL(100)+dpFComputeLength(i)
				elseif(j<1)then
					dpVCountL(1)=dpVCountL(1)+dpFComputeLength(i)
				else
					dpVCountL(j)=dpVCountL(j)+dpFComputeLength(i)
				endif
			endif
		enddo
	endif

	!///////////////////////////////////////
	!	Combine the number of
	call CombineVelocityNumber

end subroutine VelocityStat
!***********************************************************************
!  UpdateMasterDis: Send Dislocation back to master. This is only done
!                   for each stress step.
!	
!	On Master:
!		iNtotalPlanes=iNPlane after distribution of glide planes. iNPlane
!	may change in this update process. iNTotalPlanes will not change during 
!	the running of the program.
!
!	On Slaves:
!		iNTotalPlanes=iNPlane at the initial receiving of the planes. During
!	the running of the program, iNPlane will change. But iNTotalPlanes will not
!	change until the distribution of the iNPlane.
!
!	Process to update DPs:
!		First, to establish a connection between the master and one slave point.
!		1. Master receive new planes, update iNPlane and Gplanes.
!		2. Master send back current iNPlane such that the slave knows the index of
!	the new planes.
!		3. Slaves begin to group DPs and update glide plane information.
!		4. Slaves begin to send pointns back to master.56wexcvbnmz/. 
!***********************************************************************
subroutine UpdateMasterDis
	use ZQMPI
	use CommunicationMD,only: CollectDPsFrom, SendDPsToMaster, UpdateGlidePlanes
	implicit none
	
	call MPI_BARRIER(MPI_COMM_WORLD,iIerr)
	call UpdateGlidePlanes

	call MPI_Barrier(MPI_COMM_WORLD,iIerr)
	if(iMyid==iMaster)then
		! Gather data, if the first time, data has been stored in the allocated
		! array. The first time to build tree, data from file. Not the first 
		! time to build tree, data from different process. Each time build the tree
		call CollectDPsFrom
	else 	! Slave
		call SendDPsToMaster
	endif
	
	call MPI_Barrier(MPI_COMM_WORLD,iIerr)

end subroutine UpdateMasterDis
!*****************************************************************
!	ShortRangeInter:
!		Crossslip, annihilation etc.
!	Updated: ZQ Wang
!	04/08/06:
!	Added variable iCrossslip
!*****************************************************************
subroutine ShortRangeInter
	use ZQMPI
	use vectors
	use variables,only:iloop,logAnnihilation, logArrangeLoop,iArr,iCrossslip
	use AnnihilationMD,only:lLogAnni,iAnniCount
	use CrossSlipMD,only:iCSCount
	implicit none

	if(iMyid==iMaster)then

		lLogAnni=.false.

		!....cross slip....
		if(iCrossSlip==1 .and. mod(iloop,5)==0)then
			call CrossSlip
			if(iCSCount>0)lLogAnni=.true.
		endif
		
		!....Annihilation:self or interact....
		if(mod(iloop,1)==0)then
			if(logAnnihilation==1)then
				call NewAnni
				if(iAnniCount>0)lLogAnni=.true.
			endif
		endif

		!....rearrange points on dislocations...
		if(mod(iloop,iArr)==0 .or. lLogAnni)then
			if(logArrangeLoop==1)then
				call LoopReArrange
			endif
		endif
	endif

end subroutine ShortRangeInter

!***********************************************************
!   Output: 
!			Done on master.
!		Dislocations output.                         
!		Screen output.
!----------------------------------------------------------
! Modified: ZQ Wang, 01/28/06
!    New velocity distribution output
!***********************************************************
Subroutine Output
	use ZQMPI
	use vectors
	use MPIModule,only:dDensity,dScrewDensity,dEdgeMixDensity
	use variables,only:iloop,MU,StrainCurrent,totalStress,StrainRate_P,totalStrain,totalTime, &
						N_Times,DTime,dpVCountL,dpVCountN,iloopNum
	use ANNIHILATIONMD,only:iAnniCount,dpAnniLength
	use OutputMD,only:iDisOutFreq,iGeneralOutFreq
	use crossslipMD,only:iCSCount,dpCSLength,iNumTotalScrew,iNumTotalScrew1
	implicit none
	
	integer::i

	if(iMyid==iMaster)then
		
		!//////////////////////////////////////////////
		!Screen Ouput
		write(*,*)iMyid,"Strain&Stress:",TotalStrain,totalStress*MU*1d-6,"MPa"
		write(*,*)iMyid,"Plastic Strain:",StrainCurrent,StrainRate_p
		write(*,*)iMyid,"Density:",dDensity/1d10,"cm/cm3"

		!///////////////////////////////////////////////
		!Write_past
		call write_Past

		!////////////////////////////////
		!open new files
		call newOutputFiles
		
		!///////////////////////////////////////////////
		!OUTPUT STRAIN, STRAIN RATE, AND STRESS
		write(20,*)totalStrain,totalStress*MU    !stress-strain
		write(21,*)totalTime,totalStress*MU      !time-stress
		write(22,*)totalTime,strainCurrent       !total plastic strain
		write(23,*)totalStrain,iAnniCount        
		write(31,*)totalStrain,dpAnniLength
		write(24,*)totalStrain,dDensity/10d9     ! in unit of 10^10 cm/cm^3
		write(27,*)totalStress,dDensity/10d9
		write(33,*)totalStrain,dScrewDensity/10d9
		write(34,*)totalStrain,dEdgeMixDensity/10d9
		write(28,*)totalStrain,iCSCount
		write(32,*)totalStrain,StrainCurrent,iCSCount,iNumTotalScrew,iNumTotalScrew1
		write(30,*)totalStrain,dpCSLength
		write(29,*)totalStrain,iLoopNum
		
		!...Velocity count output....
		if(mod(iloop,iGeneralOutFreq)==0)then
			do i=1,100
				if(i==1)then
					write(25,*)0,0
					write(26,*)0,0
				endif
				write(25,*)(i-1)*1.d0/100.d0,dpVCountL(i)
				write(26,*)(i-1)*1.d0/100.d0,dpVCountN(i)
				
				write(25,*)i*1.d0/100.d0,dpVCountL(i)
				write(26,*)i*1.d0/100.d0,dpVCountN(i)

				if(i==100)then
					write(25,*)1.0,0
					write(26,*)1.0,0
				endif
			enddo
		endif
		
!		call OutputVelocityCenter

		!/////////////////////////////////////////////
		!Dislocation Dynamics
		if(mod(iloop,iDisOutFreq)==0)then
			call outputdis
		endif

	endif

End subroutine Output
!**************************************************************
!   AdjustTimeStep:                                           *
!		Done on master and slaves.							  *
!**************************************************************
subroutine AdjustTimeStep
	use ZQMPI
	use vectors
	use variables,only:N_TIMES,N_TIMESCOPY,SR1,SR2,DTIME
	implicit none

	if(dabs(SR2)<dabs(SR1))then   ! stress decreases, decrease the time step
		N_TIMES=0.618*N_TIMES
		if(N_TIMES<=200)N_TIMES=200
	else
		N_TIMES=N_TIMES/0.618
		if(N_TIMES>=N_TIMESCOPY)N_TIMES=N_TIMESCOPY
	endif
	SR1=SR2

end subroutine AdjustTimeStep

!***********************************************************************
!   RemoveVariables: If next step the tree is to be rebuilt, currently 
!                    allocate variables should be deallocate.
!          Done by slaves.
!          New: this step should be done before rebuild the tree because
!               master needs to gather the dislcoations back.
!
!
!	Variables allocated in the program:
!			(variable)	->	(Deallocation)
!		MPIMODULE:
!			iaGLIndex	->	No
!			iaTransList	->	No
!			tpBoxRoots	->	No
!			iIndexOfCloseBox	->	No
!			tTreeRoot	->	Yes
!		Variables:
!			GPlanes		->	No
!			DPoints		->	No
!			DMatrix 	->	No
!***********************************************************************
subroutine RemoveVariables
	use ZQMPI
	use MPIModule,only:tTreeRoot,tLocalRoot,tpBoxRoots
	implicit none

	integer :: status,i

	call DeallocateTree(tTreeRoot)

	if(iMyid/=iMaster)then 
		Nullify(tLocalRoot)
	endif
	do i=1,iNumProcs
		Nullify(tpBoxRoots(i)%p)
	enddo

end subroutine RemoveVariables
!***********************************************************************
!    DeallocateTree: deallocate trees on processes
!***********************************************************************
recursive subroutine DeallocateTree(root)
	use ZQMPI
	use TreeLeaf
	implicit none

	integer :: status
	type(tTreeLeaf),pointer::root

	if(root%lNodeOrLeaf)then                 ! if it is node, deallocate its children first.
		call DeallocateTree(root%tpChildL)
		call DeallocateTree(root%tpChildR)
	end if                                     ! if it is leaf, deallocate it now.
	if(root%iNumDPs/=0)deallocate(root%tapDPs, STAT=status)
	deallocate(root,STAT=status)

end subroutine DeallocateTree
!***************************************************************
!	NewOutputFiles:
!		Close old files and open new files for output.
!---------------------------------------------------------------
! Modified: ZQ Wang, 01/28/06
!	Added new file 28.
!	02/16/06:
!	Modified some file numbers.
!***************************************************************
subroutine NewOutputFiles
	use ZQMPI
	use vectors
	use variables,only:iloop
	use OutputMD,only:iDisOutFreq,iGeneralOutFreq
	implicit none

	integer::i1,i2,i3,i4,i5
	character(LEN=40)::filename

	i1=iloop/10000
	i2=(iloop-i1*10000)/1000
	i3=(iloop-i1*10000-i2*1000)/100
	i4=(iloop-i1*10000-i2*1000-i3*100)/10
	i5=mod(iloop,10)
	!....geometry files.....
	if(mod(iloop,iDisOutFreq)==0)then
		close(11)  

		FILENAME="G"// char(i1+ichar('0'))//     &
			  char(i2+ichar('0'))//char(i3+ichar('0'))//char(i4+ichar('0'))//char(i5+ichar('0'))//".dat"
		open(unit=11,file=FILENAME)

	endif
	
	!....mechanics files.....
	if(mod(iloop,iGeneralOutFreq)==0.and.iloop>1)then
		close(20)
		FILENAME="strain_stress"// char(i1+ichar('0'))//     &
		          char(i2+ichar('0'))//char(i3+ichar('0'))//char(i4+ichar('0'))//char(i5+ichar('0'))//".dat"
		OPEN(unit=20, FILE=FILENAME) 

		close(21)
		FILENAME="time_stress"// char(i1+ichar('0'))//     &
		          char(i2+ichar('0'))//char(i3+ichar('0'))//char(i4+ichar('0'))//char(i5+ichar('0'))//".dat"
		OPEN(unit=21, FILE=FILENAME) 

		close(22)
		FILENAME="time_strain"// char(i1+ichar('0'))//     &
		          char(i2+ichar('0'))//char(i3+ichar('0'))//char(i4+ichar('0'))//char(i5+ichar('0'))//".dat"
		OPEN(unit=22, FILE=FILENAME) 

		close(23)
		FILENAME="strain_afreq"// char(i1+ichar('0'))//     &
		          char(i2+ichar('0'))//char(i3+ichar('0'))//char(i4+ichar('0'))//char(i5+ichar('0'))//".dat"
		OPEN(unit=23, FILE=FILENAME) 

		close(24)
		FILENAME="strain_density"// char(i1+ichar('0'))//     &
		          char(i2+ichar('0'))//char(i3+ichar('0'))//char(i4+ichar('0'))//char(i5+ichar('0'))//".dat"
		OPEN(unit=24, FILE=FILENAME) 

		close(25)
		FILENAME="VeoL"// char(i1+ichar('0'))//     &
		          char(i2+ichar('0'))//char(i3+ichar('0'))//char(i4+ichar('0'))//char(i5+ichar('0'))//".dat"
		OPEN(unit=25, FILE=FILENAME) 

		close(26)
                FILENAME="VeoN"// char(i1+ichar('0'))//     &
                          char(i2+ichar('0'))//char(i3+ichar('0'))//char(i4+ichar('0'))//char(i5+ichar('0'))//".dat"
                OPEN(unit=26, FILE=FILENAME)

                close(27)
                FILENAME="stress_density"// char(i1+ichar('0'))//     &
                          char(i2+ichar('0'))//char(i3+ichar('0'))//char(i4+ichar('0'))//char(i5+ichar('0'))//".dat"
                OPEN(unit=27, FILE=FILENAME)

                close(28)
                FILENAME="strain_cs"// char(i1+ichar('0'))//     &
                          char(i2+ichar('0'))//char(i3+ichar('0'))//char(i4+ichar('0'))//char(i5+ichar('0'))//".dat"
                OPEN(unit=28, FILE=FILENAME)

		close(30)
                FILENAME="strain_csl"// char(i1+ichar('0'))//     &
                          char(i2+ichar('0'))//char(i3+ichar('0'))//char(i4+ichar('0'))//char(i5+ichar('0'))//".dat"
                OPEN(unit=30, FILE=FILENAME)

                close(29)
                FILENAME="strain_loopnum"// char(i1+ichar('0'))//     &
                          char(i2+ichar('0'))//char(i3+ichar('0'))//char(i4+ichar('0'))//char(i5+ichar('0'))//".dat"
                OPEN(unit=29, FILE=FILENAME)

		close(31)
                FILENAME="strain_al"// char(i1+ichar('0'))//     &
                          char(i2+ichar('0'))//char(i3+ichar('0'))//char(i4+ichar('0'))//char(i5+ichar('0'))//".dat"
                OPEN(unit=31, FILE=FILENAME)
		
		close(32)
		FILENAME="strain_csn"// char(i1+ichar('0'))//     &
                          char(i2+ichar('0'))//char(i3+ichar('0'))//char(i4+ichar('0'))//char(i5+ichar('0'))//".dat"
                OPEN(unit=32, FILE=FILENAME)

		close(33)
                FILENAME="strain_screwd"// char(i1+ichar('0'))//     &
                          char(i2+ichar('0'))//char(i3+ichar('0'))//char(i4+ichar('0'))//char(i5+ichar('0'))//".dat"
                OPEN(unit=33, FILE=FILENAME)
		
		close(34)
		FILENAME="strain_edgemixd"// char(i1+ichar('0'))//     &
                          char(i2+ichar('0'))//char(i3+ichar('0'))//char(i4+ichar('0'))//char(i5+ichar('0'))//".dat"
                OPEN(unit=34, FILE=FILENAME)

	endif

end subroutine NewOutputFiles
!*********************************************************************
!OutputvelocityCenter
!
!*********************************************************************
subroutine OutputVelocityCenter
        use ZQMPI
        use vectors
        use variables,only:iloop,N_Times,iNPoint,DPoints,DTime,Lattice
        implicit none

        integer::i,j,k

        k=1
!       k=2   !1: single case; 2 dipole case
        if(k==1)then
                if(mod(iNpoint,2)==0)then   !even point, odd segment
                        i=iNPoint/2
                        j=i+1
                        write(999,*)int(iloop*N_Times*DTime/1d-12),0.5*(MAG(DPoints(i)%tvPreV)+MAG(DPoints(j)%tvPreV))*Lattice
                else
                        i=iNPoint/2+1
                        write(999,*)int(iLoop*N_Times*DTime/1d-12),MAG(DPoints(i)%tvPreV)*Lattice
                endif
        else
                i=iNPoint/2/2
                j=i+1
                write(999,*)int(iloop*N_Times*DTime/1d-12),0.5*(MAG(DPoints(i)%tvPreV)+MAG(DPoints(j)%tvPreV))*Lattice
                i=iNPoint/2+i
                j=i+1
                write(998,*)int(iloop*N_Times*DTime/1d-12),0.5*(MAG(DPoints(i)%tvPreV)+MAG(DPoints(j)%tvPreV))*Lattice
        endif

end subroutine OutputVelocityCenter

!********************************************************************* 
END Module UpdateMD
