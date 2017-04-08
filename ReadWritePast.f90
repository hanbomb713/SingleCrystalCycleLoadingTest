!**************************************************************************
!  Write_Past:                                                            
!		Write current simulation information into files for continued 
!		calculation.
!-----------------------------------------------------------------------
! Modified: ZQ Wang, 01/24/06
!		Added strainRate to the list. This is for control of cyclic deformation
!	to adjust the direction of the applied load.
!
!	---------
!	01/26/06:
!		Added the dpMass and tvAcce for dislocation points.
!**************************************************************************
subroutine Write_Past
	use ZQMPI
	use vectors
	use variables,only:Sig_APP,MU,TotalStrain,totalStress,TotalTime, StrainRate,iloop,&
					iNPoint,iNPlane,DPoints,GPlanes,strainCurrent
	use OutputMD,only:iPastFreq
	implicit none
	
	integer::i,j,k,i1,i2,i3,i4,i5
	character(len=40)::FILENAME
	logical::lPPast

	lPPast=.false.

	!....open file....
	i=mod(iloop,2)
	if(i==1)then
		FILENAME="past_material0.txt"
		open(unit=111,file=filename)
		FILENAME="past_geometry0.txt"
		open(unit=112,file=filename)
	else
		FILENAME="past_material1.txt"
		open(unit=111,file=filename)
		FILENAME="past_geometry1.txt"
		open(unit=112,file=filename)
	endif
	!keep permanent files.
	if(mod(iloop,iPastFreq)==0)then 		!Every iPastFreq steps keep permanent past files.
		lPPast=.true.
		i1=iloop/10000
		i2=(iloop-i1*10000)/1000
		i3=(iloop-i1*10000-i2*1000)/100
		i4=(iLoop-i1*10000-i2*1000-i3*100)/10
		i5=iLoop-i1*10000-i2*1000-i3*100-i4*10
		FILENAME="past_material"//char(i1+ichar('0'))//             &
				char(i2+ichar('0'))//char(i3+ichar('0'))//         &
				char(i4+ichar('0'))//char(i5+ichar('0'))//".txt"
		open(unit=113,file=filename)
		FILENAME="past_geometry"//char(i1+ichar('0'))//             &
				char(i2+ichar('0'))//char(i3+ichar('0'))//         &
				char(i4+ichar('0'))//char(i5+ichar('0'))//".txt"
		open(unit=114,file=filename)
	endif
	
	!....simulation control....
	!Begin loop, sig_app(total stres),total strain, strain rate, total simulation time,
	!no. of planes,no. of nodes
	write(111,*)iloop,totalStress*MU,totalStrain,StrainCurrent,StrainRate,TotalTime,iNPlane,iNPoint
	if(lPPast)write(113,*)iloop,totalStress*MU,totalStrain,StrainCurrent,StrainRate,TotalTime,iNPlane,iNPoint

	do i=1,iNPlane
		write(112,*)GPlanes(i)%iMiller,GPlanes(i)%Origin,GPlanes(i)%iNumPEquiv
		j=GPlanes(i)%iNumPEquiv
		if(j>0)write(112,*)(GPlanes(i)%iIDPEquiv(k),k=1,j)		

		if(lPPast)then
			write(114,*)GPlanes(i)%iMiller,GPlanes(i)%Origin,GPlanes(i)%iNumPEquiv
			j=GPlanes(i)%iNumPEquiv
			if(j>0)write(114,*)(GPlanes(i)%iIDPEquiv(k),k=1,j)		
		endif
	enddo

	do i=1,iNPoint
		write(112,*)DPoints(i)%iID,DPoints(i)%iLoopID,DPoints(i)%itype,DPoints(i)%tvPL,DPoints(i)%tvTL, &
					DPoints(i)%iPlane,DPoints(i)%iBurgers, &
					DPoints(i)%iBeginP,DPoints(i)%iEndP,DPoints(i)%tvAcce,DPoints(i)%tvPreV,DPoints(i)%dpMass,   &
					DPoints(i)%tvPreVT,DPoints(i)%lToCal
		if(lPPast)write(114,*)DPoints(i)%iID,DPoints(i)%iLoopID,DPoints(i)%itype,DPoints(i)%tvPL,DPoints(i)%tvTL, &
					DPoints(i)%iPlane,DPoints(i)%iBurgers, &
					DPoints(i)%iBeginP,DPoints(i)%iEndP,DPoints(i)%tvAcce,DPoints(i)%tvPreV,DPoints(i)%dpMass,  &
					DPoints(i)%tvPreVT,DPoints(i)%lToCal
	enddo

	close(111)
	close(112)
	if(lPPast)then
		close(113)
		close(114)
	endif

end subroutine Write_Past
!**********************************************************************
!	Read_Past:
!		Read previous information from files to continue simulation.
!**********************************************************************
subroutine Read_Past
	use ZQMPI
	use vectors
	use variables,only:APPLIED_SIG,iNPoint,iNPlane,GPlanes, StrainRate,  &
						DPoints,totalStrain,iBeginLoop,totalTime,StrainCurrent
	implicit none

	integer::i,j,k
	character(len=40)::FILENAME

	!....open file.....
	read(3,*)i
	if(i==1)then
		open(unit=111,file="past_material1.txt")
		open(unit=112,file="past_geometry1.txt")
	else
		open(unit=111,file="past_material0.txt")
		open(unit=112,file="past_geometry0.txt")
	endif

	!Begin loop, sig_app(total stres),total strain, total simulation time,
	!no. of nodes, no. of planes
	read(111,*)iBeginLoop, APPLIED_SIG,totalStrain,StrainCurrent,StrainRate,totalTime,iNPlane,iNPoint

	do i=1,iNPlane
		read(112,*)GPlanes(i)%iMiller,GPlanes(i)%Origin,GPlanes(i)%iNumPEquiv
		if(GPlanes(i)%iNumPEquiv>0)then
			read(112,*)(GPlanes(i)%iIDPEquiv(j),j=1,GPlanes(i)%iNumPEquiv)
		endif
	enddo

	do i=1,iNPoint
		Read(112,*)DPoints(i)%iID,DPoints(i)%iLoopID,DPoints(i)%itype,DPoints(i)%tvPL,DPoints(i)%tvTL, &
					DPoints(i)%iPlane,DPoints(i)%iBurgers, &
					DPoints(i)%iBeginP,DPoints(i)%iEndP,DPoints(i)%tvAcce,DPoints(i)%tvPreV,DPoints(i)%dpMass, &
					DPoints(i)%tvPreVT,DPoints(i)%lToCal
	enddo

	close(111)
	close(112)

end subroutine Read_Past
