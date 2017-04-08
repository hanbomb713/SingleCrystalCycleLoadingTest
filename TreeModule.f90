!***************************************************************************************************
!   TreeMD:
!         1. Defines the tree structure
!         2. Defines some basic tree subroutines
!***************************************************************************************************
Module TreeMD
  use TreeLeaf
CONTAINS

!***********************************************************************
!     NewBuildGlobalTree:  Build global tree with all the dislocation info
!                       has already been put on master.
!               Done on master now.
!***********************************************************************
subroutine NewBuildGlobalTree
        use ZQMPI
        use vectors
        use MPIModule,only:tTreeRoot
        implicit none

        if(iMyid==iMaster)then

                call NewMakeNewRoot(tTreeRoot)            ! Build new tree,i=1 means it is root
                call NewSplitDomain(tTreeRoot)

        endif
end subroutine NewBuildGlobalTree
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!	NewMakeNewRoot
!	For each level of the tree, a pointer array is allocated to point to 
!	the DPs within this node. It will be deallocated at next level.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine NewMakeNewRoot(root)
        use ZQMPI
        use vectors
        use MPIModule,only:iCurrentNodeIndex,iNumTNode,tpBoxRoots,iNMaxS
        use variables,only:A_Cube,iNPoint,DPoints
        implicit none

        type(tTreeLeaf),pointer::root
        integer ::i

        allocate(root)
        root%lNodeOrLeaf=.false.          ! first assign false value as leaf, later modify it

        root%iID = 1
        root%iIndex=1
	root%iLevel=iNumProcs-1            !---iLevel is indicating how many leaves finally below this level.!
        iCurrentNodeIndex = 1
        tpBoxRoots(root%iID)%p=>root

        root%tCenterv%v(1)=0.5d0 * A_Cube
        root%tCenterv%v(2)=0.5d0 * A_Cube
        root%tCenterv%v(3)=0.5d0 * A_Cube
        root%tSizev%v(1)=0.5d0 * A_Cube
        root%tSizev%v(2)=0.5d0 * A_Cube
        root%tSizev%v(3)=0.5d0 * A_Cube
        NULLIFY(root%tpChildL)
        NULLIFY(root%tpChildR)

	!Allocate pointers for dislocation points within this root
        root%iNumDPs=iNPoint
        allocate(root%tapDPs(iNPoint))
	do i=1,iNPoint
		root%tapDPs(i)%p=>DPoints(i)
	enddo
        iNumTNode=1

end subroutine NewMakeNewRoot
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!  New Spliting of domains
!
!	0). Already known the DPs in this level.
!
!	1). Start from tTreeRoot, at each level of the tree, the split ratio is 
!	calculated. Then the domain is split according to this ratio. The procedure 
!	is repeated until the lowest level is reached.
!
!	2). iLevel will be used to indicate how many final leaves will be under this level.
!	So, the leaves will have iLevel=1.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
recursive subroutine NewSplitDomain(root)
	use ZQMPI
	use TreeLeaf 
	implicit none

	type(tTreeLeaf),pointer::root

	!------------------------------
	!if already lowest level, return
	if(root%iLevel==1)return
	!-----------------------------
	!otherwise, split

	call AllocateNewChildren(root)

	call CalculateRatioSplit(root)

	call NewSplitDomain(root%tpChildL)
	call NewSplitDomain(root%tpChildR)

end subroutine NewSplitDomain
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!      AllocateNewChildren
!--------------------------------------------------------------------------------
subroutine AllocateNewChildren(root)
	use ZQMPI
	use MPIModule,only:iNumTNode,tpBoxRoots,iCurrentNodeIndex,iNMaxS
	implicit none

	type(tTreeLeaf),pointer::root
	integer::iRight,iLeft,itemp

!	!Allocate children
	allocate(root%tpChildL)
	allocate(root%tpChildR)

	root%tpChildL%tParent=>root
	root%tpChildR%tParent=>root
        !--------modify former leaf to node
        root%lNodeOrLeaf=.TRUE.
        root%tpChildL%lNodeOrLeaf=.false.
        root%tpChildR%lNodeOrLeaf=.false.

        !--------set iID
        root%tpChildL%iID = root%iID
        root%iID = 0
        iCurrentNodeIndex = iCurrentNodeIndex +1
        root%tpChildR%iID = iCurrentNodeIndex

        !--------set pointer to the processor nodes
        if(iMyid==iMaster)then
                tpBoxRoots(root%tpChildL%iID)%p=>root%tpChildL   !?
                tpBoxRoots(root%tpChildR%iID)%p=>root%tpChildR   !?
        end if

	!--------set total number of nodes
	iNumTNode=iNumTNode+2
        root%tpChildL%iIndex=iNumTNode-1
        root%tpChildR%iIndex=iNumTNode

        root%iChild(1)=iNumTNode-1
        root%iChild(2)=iNumTNode

        root%tpChildL%iNumDPs=0
        root%tpChildR%iNumDPs=0

        NULLIFY(root%tpChildL%tpChildL)
        NULLIFY(root%tpChildL%tpChildR)
        NULLIFY(root%tpChildR%tpChildL)
        NULLIFY(root%tpChildR%tpChildR)

	!--------set level
	iLeft=root%iLevel/2
	iRight=root%iLevel-iLeft
	root%tpChildL%iLevel=iLeft
	root%tpChildR%iLevel=iRight

end subroutine AllocateNewChildren
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! CalculateRatioSplit
!	Only deal with DPs in the nodes.    
!
!	Already have the level number, calculate the split ration 
!	and divide the space
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
subroutine CalculateRatioSplit(root)
	use ZQMPI
	use vectors
	use TreeLeaf
	implicit none

	type(tTreeLeaf),pointer::root
	type(vector)::tvSplitP          !the point where the space is going to be split
	double precision::dpRatio
	integer::dpSD          !split direction
	integer::i
	integer::status
	integer::iMaxIndex,iMinIndex                        !pointing to the index with max and min distances
	integer::iCIndex                                    !the index to compare to.
	integer::iNumAss                                    !the number of DPs assigned to the left
	double precision::dpTmp

	!--------------------------------
	!determine the split direction, always split along the largest size
	dpSD=1
	do i=2,3
		if(dabs(root%tsizev%v(i))>dabs(root%tsizev%v(dpSD)))then
			dpSD=i
		endif
	enddo

	!--------------------------------
	!calculate ratio
	dpRatio=(root%tpChildL%iLevel*1.d0)/(root%iLevel*1.d0)
	root%tpChildL%iNumDPs = root%iNumDPs*dpRatio
	root%tpChildR%iNumDPs = root%iNumDPs-root%tpChildL%iNumDPs
	allocate(root%tpChildL%tapDPS(root%tpChildL%iNumDPS+20))
	allocate(root%tpChildR%tapDPS(root%tpCHildR%iNumDPS+20))

	!----------------------------------
	!arrange DPs so to divide according to numbers of iNumDPS
	allocate(itIndex(root%iNumDPS,2))
	do i=1,root%iNumDPS
		itIndex(i,1)=-1;itIndex(i,2)=-1
	enddo

	iMaxIndex=1
	iMinIndex=1
	!-------arrange
	do i=2,root%iNumDPS    ! compare distance and put them in right order
		if(root%tapDPS(i)%p%tvPG%v(dpSD)>=root%tapDPS(iMaxIndex)%p%tvPG%v(dpSD))then     !largest
			itIndex(iMaxIndex,2)=i
			iMaxIndex=i
		elseif(root%tapDPS(i)%p%tvPG%v(dpSD)<root%tapDPS(iMinIndex)%p%tvPG%v(dpSD))then    !smallest
			itIndex(iMinIndex,1)=i
			iMinIndex=i
		else                                                                        !in the middle
			iCIndex=1       ! 1 is alway the root, iCIndex control the tree index
			do 
				if(root%tapDPS(i)%p%tvPG%v(dpSD)>=root%tapDPS(iCIndex)%p%tvPG%v(dpSD))then  !larger
					if(itIndex(iCIndex,2)==-1)then    !bottom, add large
						itIndex(iCIndex,2)=i
						exit
					else                               !continue going left
						iCIndex=itIndex(iCindex,2)
					endif
				else            !smaller
					if(itIndex(iCIndex,1)==-1)then    !bottom, add small
						itIndex(iCIndex,1)=i
						exit
					else
						iCIndex=itIndex(iCindex,1)
					endif
				endif
			enddo
		endif
	enddo

	!--------count and calculate the position to split
	iNumAss=0
	iCIndex=1
	call AssignTreeDPs(root,iNumAss,iCindex)
	deallocate(itIndex,STAT=status)

	!--------size of the two children
	root%tpChildL%tsizeV=root%tsizeV
	root%tpChildR%tsizeV=root%tsizeV
	root%tpChildL%tCenterV=root%tCenterV
	root%tpChildR%tCenterV=root%tCenterV

	dpTmp=0.5*(root%tpChildL%tapDPs(root%tpChildL%iNumDPs)%p%tvPG%v(dpSD)+root%tpChildR%tapDPs(1)%p%tvPG%v(dpSD))
	root%tpChildL%tsizeV%v(dpSD)=0.5*dabs(root%tCenterV%v(dpSD)-root%tSizeV%v(dpSD)-dpTmp)
	root%tpChildR%tsizeV%v(dpSD)=root%tsizev%v(dpSD)-root%tpChildL%tsizev%v(dpSD)

	root%tpChildL%tcenterV%v(dpSD)=dpTmp-root%tpChildL%tsizeV%v(dpSD)
	root%tpChildR%tcenterV%v(dpSD)=dpTmp+root%tpChildR%tsizeV%v(dpSD)

	!----------------------------------
	!reset parent
	root%iNumDPs=0
	deallocate(root%tapDPs,STAT=status)

end subroutine CalculateRatioSplit
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! AssignTreeDPs(root,iNumAss)
!	
!	recursively go through the itIndexArray
!
!	iNumAss enters as zero
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
recursive subroutine AssignTreeDPs(root,iNumAss,iCIndex)
	use ZQMPI
	use TreeLeaf
	use vectors
	implicit none

	type(tTreeLeaf),pointer::root
	integer::iNumAss
	integer,intent(in)::iCIndex

	integer::iRight,iLeft

	iRight=itIndex(iCIndex,2);iLeft=itIndex(iCIndex,1)

	!surf left
	if(iLeft/=-1)then
		call AssignTreeDPs(root,iNumAss,iLeft)
	endif

	!surf right
	if(iRight/=-1)then
		call AssignTreeDPs(root,iNumAss,iRight)
	endif

	!Assign this DP.
	iNumAss=iNumAss+1
	if(iNumAss<=root%tpChildL%iNumDPs)then   !assign to left
		root%tpChildL%tapDPs(iNumAss)%p=>root%tapDPs(iCIndex)%p
	else                                    !
		root%tpChildR%tapDPs(iNumAss-root%tpChildL%iNumDPs)%p=>root%tapDPs(iCindex)%p
	endif

end subroutine AssignTreeDPs
!******************************************************************************
!    PutSubInTree:   put dislocation points in tree represented by argument
!					 "root".
!           arguments:
!                    i--1.Global;2.local. So that each node contains a certain
!                       number of subs associated with global or local choices.
!					iID: global index of the dislocation points.
!					root: the tree node where the root should be first put in.
!******************************************************************************
recursive subroutine PutSubInTree(iID, root, i)
	use ZQMPI
	use vectors
	use TreeLeaf
	use MPIModule,only:iNMaxS,iCurrentNodeIndex
	use variables,only:iNPoint,DPoints
	implicit none

	integer,intent(in) :: iID,i
	type(tTreeLeaf),Pointer :: root

	integer :: itemp,jtemp
	logical :: lgAdd

	if(i==1)then
		itemp=iNMaxS         ! Number of points in the leaf
		jtemp=iNumprocs-1    ! Number of leaves
	else
		itemp=1              ! number of points in the leaf
		jtemp=iNPoint   ! number of leaves
	endif

	lgAdd=.false.

	if(.not. root%lNodeOrLeaf)then
		if(i==2)then
			!/////////////////////////////////////////
			!if the number of points is less than 1 or not overlap points
			if(root%iNumDPs<1)then
				lgAdd=.true.
			elseif(.not.(MAG(root%tapDPs(1)%p%tvPG-DPoints(iID)%tvPG)<1.d-6))then
				lgAdd=.true.
			endif
		else if(i==1)then
			lgAdd=.true.
		endif

		if(lgAdd)then	
			root%iNumDPs=root%iNumDPs+1
			root%tapDPs(root%iNumDPs)%p=>DPoints(iID)                        ! ?

			if(root%iNumDPs>itemp .and. iCurrentNodeIndex<jtemp)call SplitRoot(root,i)
		endif
	
	elseif(IsInBox(DPoints(iID)%tvPG,root%tpChildL%tCenterv,root%tpChildL%tSizev))then

		call PutSubInTree(iID,root%tpChildL,i)

	else

		call PutSubInTree(iID,root%tpChildR,i)

	endif

end subroutine PutSubInTree
!***********************************************************************
!    SplitRoot: 
!		Divide the root and split the points into two groups.
!    
!     iID is set to zero for those non-leaf nodes.
!***********************************************************************
subroutine SplitRoot(root,i)
	use ZQMPI
	use vectors
	use TreeLeaf
	use MPIModule,only:iNMaxS,iNumTNode,iCurrentNodeIndex,tpBoxRoots
	use variables,only:zero
	implicit none

	integer,intent(in):: i
	integer :: status,itemp
	type(tTreeLeaf),pointer::root

	if(i==1)then ! global tree
		itemp = iNMaxS
	else
		itemp=1
	endif

	! split once, add two node to the
	iNumTNode = iNumTNode+2

	! allocate new children nodes
	allocate(root%tpChildL)
	allocate(root%tpChildR)

	!modify former leaf to node
	root%lNodeOrLeaf=.TRUE.
    root%tpChildL%lNodeOrLeaf=.false.
    root%tpChildR%lNodeOrLeaf=.false.

	!set iID
	root%tpChildL%iID = root%iID
	root%iID = 0
	iCurrentNodeIndex = iCurrentNodeIndex +1
	root%tpChildR%iID = iCurrentNodeIndex

	!///////////////////////////////////////////////
	!set pointer to the processor nodes
	if(iMyid==iMaster)then
		tpBoxRoots(root%tpChildL%iID)%p=>root%tpChildL   !?
		tpBoxRoots(root%tpChildR%iID)%p=>root%tpChildR   !?
	end if

	root%tpChildL%iIndex=iNumTNode-1
	root%tpChildR%iIndex=iNumTNode

	root%iChild(1)=iNumTNode-1
	root%iChild(2)=iNumTNode

	root%tpChildL%iNumDPs=0
	root%tpChildR%iNumDPs=0

	NULLIFY(root%tpChildL%tpChildL)
	NULLIFY(root%tpChildL%tpChildR)
	NULLIFY(root%tpChildR%tpChildL)
	NULLIFY(root%tpChildR%tpChildR)

	allocate(root%tpChildL%tapDPs(int(5*itemp)))             !?  later. the number 5 here is also related to allocation in SplitBuild ---iTPID, etc.
	allocate(root%tpCHildR%tapDPs(int(5*itemp)))             !?  later

	! calculate the center of dividing space, move subs.
	call MoveSegmentsToLR(root,i)      

	! after moving DPs to children, set number to zero and deallocate
	! array.
	root%iNumDPs=0
	deallocate(root%tapDPs,STAT=status)

end subroutine SplitRoot
!***********************************************************************
!    MovSegmentsToLR: move segments, set children sizes and centers
!***********************************************************************
subroutine MoveSegmentsToLR(root,iGLcase)
	use ZQMPI
	use vectors
	use TreeLeaf
	use FunctionMD,only:PointPBC
	use MPIModule,only:iNMaxS
	use variables,only:zero
	implicit none
    
	type(tTreeLeaf),pointer::root
	integer,intent(in)::iGLcase
	integer:: i,j,status,NumDPsInNode
	type(vector):: dptemp,dptemp2,dptemp3
	type(vector):: tMassCenter,tMassCenter2
	double precision :: itemp
	integer::ii,jj(3),kk

	!///////////////////////////////////////////////
	if(iGLcase==1)then
		NumDPsInNode=iNMaxS         ! Number of DPs should be in the leaf
	else
		NumDPsInNode=1              ! number of DPs should be in the leaf
	endif

	tMassCenter=zero%v(1)
	
	!---calculate tMassCenter
	do i=1, root%iNumDPs
		tMassCenter=tMassCenter+PointPBC(root%tapDPs(i)%p%tvPG)
	enddo
	tMassCenter=1.d0/(root%iNumDPs*1.d0)*tMassCenter

	dptemp2=PointPBC(root%tapDPs(1)%p%tvPG)
	dptemp3=PointPBC(root%tapDPs(2)%p%tvPG)
	if(iMyid/=iMaster)then
		tMassCenter=0.5d0*(dptemp2+dptemp3)
	endif
	tMassCenter2=tMassCenter

	!		
	! split method 1: 1, 2, 3 in turn

	!---split method 2: randomly choose one largest deviation.

	ii=0
	do i=1,3
		if(tMassCenter%v(i)/=dptemp2%v(i).and.tMassCenter%v(i)/=dptemp3%v(i))then    ! Here can be played.
			ii=ii+1
			jj(ii)=i-1
		endif
	end do

	if(ii/=0)then
		j=jj(1)
		itemp=dabs(root%tSizev%v(j+1))
		do i=2,ii
			if(itemp<=dabs(root%tSizev%v(jj(i)+1)))then
				itemp=dabs(root%tSizev%v(jj(i)+1))
				j=jj(i)
			endif
		enddo
	else
		do i=1,3
			if(dptemp2%v(i)/=dptemp3%v(i))then
				j=i-1
				exit
			endif
		enddo
	endif
	!---end of split------------------
		
	!---reset tMassCenter, only keep information of j+1
	kk=j
	do i=1,3
		if(i/=j+1)then
			tMassCenter%v(i) = root%tCenterv%v(i)
		end if
	end do
	!---end tMassCenter

	! new center and size of children
	dptemp=tMassCenter

	dptemp%v(j+1)=root%tCenterv%v(j+1)-root%tSizev%v(j+1)
	root%tpChildR%tCenterv = 0.5d0 * (dptemp+tMassCenter)
	root%tpChildR%tSizev=root%tSizev
	root%tpChildR%tSizev%v(j+1)=          &
		dabs(0.5d0 *(dptemp%v(j+1)-tMassCenter%v(j+1)))


	dptemp%v(j+1)=root%tCenterv%v(j+1)+root%tSizev%v(j+1)
	root%tpChildL%tCenterv=0.5d0 * (dptemp+tMassCenter)
	root%tpChildL%tSizev=root%tSizev
	root%tpChildL%tSizev%v(j+1)=          &
		dabs(0.5d0 *(dptemp%v(j+1)-tMassCenter%v(j+1)))
   
	! move subsegments
	do i=1, root%iNumDPs
		if(IsInBox(root%tapDPs(i)%p%tvPG,root%tpChildL%tCenterv,root%tpChildL%tSizev)  .and. &
		root%tpChildL%iNumDPs<NumDPsInNode)then
			root%tpChildL%iNumDPs		=	root%tpChildL%iNumDPs+1
			j							=	root%tpChildL%iNumDPs
			root%tpChildL%tapDPs(j)%p	=>	root%tapDPs(i)%p
		else if(root%tpChildR%iNumDPs<NumDPsInNode)then
			root%tpChildR%iNumDPs		=	root%tpChildR%iNumDPs+1
			j							=	root%tpChildR%iNumDPs
			root%tpChildR%tapDPs(j)%p	=>	root%tapDPs(i)%p
		endif
		NULLIFY(root%tapDPs(i)%p)
	enddo

end subroutine MoveSegmentsToLR

!***********************************************************************
!    IsInBox: check if DP is in the node box
!		rootC: center of the tree node;
!		rootS: size of the tree node.
!***********************************************************************
logical function IsInBox(R,rootC,rootS)
	use ZQMPI
	use vectors
	use TreeLeaf
	use variables,only:A_Cube
	implicit none

	integer :: i,j
	type(vector),intent(in)::R,rootC,rootS
	type(vector) :: tmpVector

	IsInBox=.true.

	do i=1,3
		if(R%v(i)<0.d0)then
			tmpVector%v(i)=dmod(R%v(i),A_CUBE)+A_CUBE
		else
			tmpVector%v(i)=dmod(R%v(i),A_Cube)
		end if
	end do

	do i=1,3
		if(dabs(tmpVector%v(i)-rootC%v(i))>rootS%v(i))then
			IsInBox=.false.
			exit
		end if
	enddo

	return

end function IsInBox
!***********************************************************************
!     BuildLocalTree: Build local tree structure and combine with the
!                     global one. The combination is done in transfering
!                     global tree to the process by pointing the local 
!                     root to the global tree.
!
!              what we have: a set of DPs and global tree.
!							tLocalroot is pointing to the node on global
!							tree which is representing the current 
!							processor.               
!***********************************************************************
subroutine BuildLocalTree
	use ZQMPI
	use MPIModule,only:tLocalRoot,iCurrentNodeIndex
	use variables,only:iNPoint,iloop
	implicit none

	integer :: i
	integer :: itemp
	
	if(iMyid/=iMaster)then
		i=2
		tLocalRoot%iNumDPs = 0
		if(iNPoint>0)then
			tLocalRoot%iID=1
			allocate(tLocalRoot%tapDPs(2))
		endif

		iCurrentNodeIndex = 1

		do itemp=1,iNPoint
			i=2
			call PutSubInTree(itemp,tLocalRoot,i)
		enddo

	endif

end subroutine BuildLocalTree

!***********************************************************************
!    WalkTree: walk through the tree to find neighbors
!         Arguments:
!            root: current node to be inspected;
!            iID: id of current dislocation point being considered
!            lmp: flag indicating whether the node should be added as 
!                 an neighbor. =0: multipole; =1: direct neighbor
!
!		The DP neighbors are now in global index. These points will be transferred 
!		to local processor later.
!***********************************************************************
recursive subroutine WalkTree(root, iID,lmp)
	use ZQMPI
	use vectors
	use TreeLeaf
	use variables,only:DPoints,dcrit
	use LoopRearrangeMD,only:dpMaxAveLength
	implicit none

	type(tTreeLeaf),pointer :: root
	integer,intent(in) :: iID
	integer :: lmp	! if both return false, this node will be considered as multipole
					! if this node is considered as multipole, it will return true.
					! lmp means add or not. true means add.
	integer :: llmp1, llmp2
	integer :: i,j
	double precision :: GetDist

	lmp= 0      ! not add, as multipole

	if(root%lNodeOrLeaf)then                  ! if it is a node.  
		if(OverLapNode(root, iID))then       ! overlaped, check child
			call WalkTree(root%tpChildL,iID,llmp1)
			call WalkTree(root%tpChildR,iID,llmp2)

			if(llmp1 ==0 .and. llmp2 ==0)then   ! if both are multipole, this one should be multipole
							! wait until upper level decides if add it along or 
							! combine with others. return lmp =0 as indicator
				lmp = 0
				if(root%iIndex==1)then     ! root
					call AddNeighborToDis(root,iID,lmp)
				end if
				return
			else if(llmp1==0)then       ! following cases are both having one direct neighbors, lmp =1
				call AddNeighborToDis(root%tpChildL,iID,llmp1)
			else if(llmp2==0)then
				call AddNeighborToDis(root%tpChildR,iID,llmp2)
			endif

			lmp = 1
			return
		else                          ! not overlap, add it as multipole
			lmp= 0
			return
		end if
	else    ! leaf, two cases: 1. only DP; 2. Multipole. check distance
		if(root%iNumDPs==0 .or. root%iID == 0) then   ! this is a multipole leaf
			lmp = 0
		else if(GetDist(root%tapDPs(root%iNumDPs)%p%tvPG,DPoints(iID)%tvPG) <= dcrit+dpMaxAveLength)then  !???
			lmp = 1  
			call AddNeighborToDis(root,iID,lmp)
		else
			lmp = 0
		endif

		return
	endif

end subroutine WalkTree

!***********************************************************************
!    OverLapNode: check if subsegment and a node in the tree is overlap 
!                 or close even not overlap
!                    
!		iID: index of DP
!***********************************************************************
logical function OverLapNode(root, iID)
	use ZQMPI
	use vectors
	use TreeLeaf
	use variables,only:DPoints,dcrit
	use LoopRearrangeMD,only:dpMaxAveLength
	implicit none

	type(tTreeLeaf),pointer :: root
	integer,intent(in) :: iID
	double precision :: dtmp,GetDist
	integer :: i

	OverLapNode = .false.
	dtmp = 0.d0

	dtmp=MAG(root%tSizev)+dcrit+dpMaxAveLength

	if(GetDist(DPoints(iID)%tvPG,root%tCenterv) <= dtmp)then
		OverLapNode = .true.
	endif

	return 

end function OverLapNode
!****************************************************************************
!    AddNeighborToDis: 
!			Add neighbors to the neighbor lists of dislocation points.
!           lmp=1: should add the neighbor as a dislocation point;
!			lmp=2: should add the neighbor as a multipole expansion.
!****************************************************************************
subroutine AddNeighborToDis(root,iID,lmp)
	use ZQMPI
	use vectors
	use TreeLeaf
	use variables,only:DPoints,Max_Neighbor
	implicit none

	type(tTreeLeaf),pointer :: root
	integer, intent(in) :: iID
	integer, intent(in) :: lmp

	integer :: i,ktemp

	!///////////////////////////////////////
	i = 1
	if( .not. InNeighborList(root, iID,lmp))then
		
		!....this means should be added as single....
		if(lmp == 1)then     
			if(DPoints(iID)%iNumNei<Max_Neighbor)then
				DPoints(iID)%iNumNei=DPoints(iID)%iNumNei+1
				ktemp=DPoints(iID)%iNumNei
				!//////////////////////////////////////////////////
				!Add the global index of the neighbor to the neighbor list.
				!The index will be changed to local index after all the 
				!NC points are transferred to local.
				DPoints(iID)%iaIDOfNei(ktemp)=root%tapDPs(root%iNumDPs)%p%iID
			end if

		!....this means should be added as multipole....
		else                 
			 ! Multipole part should be included here.
		endif
	endif

end subroutine AddNeighborToDis
!********************************************************************************
!    InNeighborList: chech if a DP is already in the neighbor lists
!                    of dislocations.
!********************************************************************************
logical function InNeighborList(root,iID,lmp)
	use ZQMPI
	use vectors
	use TreeLeaf
	use variables,only:DPoints
	implicit none

	type(tTreeLeaf),pointer:: root
	integer,intent(in)::iID
	integer,intent(in)::lmp

	integer :: i,j,itemp,ktemp,iloopID

	InNeighborList = .false.

	iloopID=DPoints(iID)%iloopID
	j=root%iNumDPs

	if(lmp==1)then
		if(root%tapDPs(j)%p%iLoopID==iloopID)then         ! the same loop
			InNeighborList = .true.
			return
		else    ! if not same loop, check list
			do itemp=1,DPoints(iID)%iNumNei
				ktemp=DPoints(iID)%iaIDOfNei(itemp)
				if(root%tapDPs(j)%p%iID==ktemp)then  ! if already in the list
					InNeighborList = .true.
					return
				endif
			enddo
		end if
	else  ! add as multipole
		InNeighborList = .true.
		return
	endif

end function InNeighborList

!***********************************************************************
!    PackUnPackTree: master packs the tree, slaves unpack the tree,called
!                    by DistributeGlobalTree
!    ok:6-2-06
!***********************************************************************
recursive subroutine PackUnPackTree(root)
	use ZQMPI
	use vectors
	use TreeLeaf
	use variables,only:zero
	use MPIModule,only:tpBoxRoots,tLocalRoot
	use treeTransfer,only:tmpTTInts,tmpTTDps,tmpTTLogNL
	implicit none
	integer :: i

	type(tTreeLeaf),pointer :: root

	!///////////////////////////////////////////////
	if(iMyid==iMaster)then     ! pack
		i=root%iIndex

		tmpTTLogNL(i)=root%lNodeOrLeaf
		
		tmpTTInts(1+(i-1)*5)=root%iID
		tmpTTInts(2+(i-1)*5)=root%iIndex
		tmpTTInts(3+(i-1)*5)=root%iNumDPs
		tmpTTInts(4+(i-1)*5)=root%iChild(1)
		tmpTTInts(5+(i-1)*5)=root%iChild(2)

		tmpTTDps(1+(i-1)*6)=root%tCenterv%v(1)
		tmpTTDps(2+(i-1)*6)=root%tCenterv%v(2)
		tmpTTDps(3+(i-1)*6)=root%tCenterv%v(3)
		tmpTTDps(4+(i-1)*6)=root%tSizev%v(1)
		tmpTTDps(5+(i-1)*6)=root%tSizev%v(2)
		tmpTTDps(6+(i-1)*6)=root%tSizev%v(3)

		if(root%lNodeOrLeaf)then			
			call PackUnPackTree(root%tpChildL)

			call PackUnPackTree(root%tpChildR)
		endif
	else                       ! unpack

		i=root%iIndex           ! preassigned

		root%lNodeOrLeaf  = tmpTTLogNL(i)

		root%iID          = tmpTTInts(1+(i-1)*5)
		root%iNumDPs	  = tmpTTInts(3+(i-1)*5) 
		root%iChild(1)    = tmpTTInts(4+(i-1)*5)
		root%iCHild(2)    = tmpTTInts(5+(i-1)*5)

		root%tCenterv%v(1)  = tmpTTDps(1+(i-1)*6)
		root%tCenterv%v(2)  = tmpTTDps(2+(i-1)*6)
		root%tCenterv%v(3)  = tmpTTDps(3+(i-1)*6)
		root%tSizev%v(1)    = tmpTTDps(4+(i-1)*6)
		root%tSizev%v(2)    = tmpTTDps(5+(i-1)*6)
		root%tSizev%v(3)    = tmpTTDps(6+(i-1)*6)

		NULLIFY(root%tpChildL)
		NULLIFY(root%tpChildR)

		if(root%lNodeOrLeaf)then     ! if it is node, then construct the children
			allocate(root%tpChildL)
			allocate(root%tpChildR)
			root%tpChildL%iIndex = root%iChild(1)
			root%tpChildR%iIndex = root%iChild(2)

			call PackUnPackTree(root%tpChildL)
			call PackUnPackTree(root%tpChildR)
		else
			tpBoxRoots(root%iID)%p => root                ! assign local pointers

			if(iMyid==root%iID)then
				! The attribute of lNodeOrLeaf will be changed later.
				! only if it is local root or in the future it is a close neighbor
				tLocalRoot => root
				! we can assign iNumSubsHere here, then we don't need to send it in the SendSubSegments and ReceiveSubSegments.
				! iNumSubsHere = root%iNumSubs
			endif
		endif
                root%iID=0 !global leaves set to zero. none zero only for lowest level node containing dislocations.
                root%iNumDPs=0
	endif

end subroutine PackUnPackTree
!******************************************************************************
!	CreateGIndexArray:
!		On master.
!		Create the global index array.
!		The processor nodes of the global tree on each local processor
!		does not have the links to DPs.
!		
!		iaGlIndex(,1)=processor
!		iaGLindex(,2)=local index
!		iaGLindex(,3)=previous point
!		iaGLindex(,4)=next point
!******************************************************************************    
subroutine CreateGIndexArray
	use ZQMPI
	use vectors
	use MPIModule,only:iaGLIndex,tpBoxRoots
	implicit none

	integer::i,j,k
	
	!////////////////////////////////////
	!For each processor
	do i=1,iNumprocs-1
		!///////////////////////////////////////////
		!For each dislocation point on the processor
		do j=1,tpBoxRoots(i)%p%iNumDPs
			k=tpBoxRoots(i)%p%tapDPs(j)%p%iID
			iaGLIndex(k,1)=i
			iaGLIndex(k,2)=j
			iaGLIndex(k,3)=tpBoxRoots(i)%p%tapDPs(j)%p%iBeginP
			iaGLIndex(k,4)=tpBoxRoots(i)%p%tapDPs(j)%p%iEndP
		enddo
	enddo
	
end subroutine CreateGIndexArray

!**************************************************************************
!   DetermineGhost: determine the ghost ndoes that need to be sent to close
!                   boxes. 
!         Once enter this code, first add info to the list, then determine 
!         if it is time to stop.
!         called by TransGhostTree
!**************************************************************************
recursive subroutine DetermineGhost(root,j,tmpS1)
	use ZQMPI
	use TreeLeaf
	use vectors
	use LoopRearrangeMD,only:dpMaxAveLength
	use variables,only:dcrit
	use MPIModule,only:tpBoxRoots
	use GhostTransMD,only:iNumOfNodeSend,tmpSendLog,tmpSendInt,tmpSendDps
	implicit none

	type(tTreeLeaf),pointer :: root
	integer,intent(in) :: j      ! which box

	double precision :: tmpS1,tmpS2,tmpS3
	intent(in)::tmpS1
	integer :: tmpIndex(2)
	integer :: ii,jj,kk,itemp,jtemp
	double precision :: GetDist

	! add some info to the list, root must be node because if it is leaf
	! the code has been terminated in above level.
	do ii=1,2
		iNumOfNodeSend(1) = iNumOfNodeSend(1) +1
		tmpIndex(ii) = iNumOfNodeSend(1)

		if(ii==1)then
			tmpSendLog(iNumOfNodeSend(1)) = root%tpChildL%lNodeOrLeaf

			tmpSendInt(1+(iNumOfNodeSend(1)-1)*5) = root%tpChildL%iID
			tmpSendInt(2+(iNumOfNodeSend(1)-1)*5) = root%tpChildL%iIndex
			tmpSendInt(3+(iNumOfNodeSend(1)-1)*5) = root%tpChildL%iNumDPs
			tmpSendInt(4+(iNumOfNodeSend(1)-1)*5) = 0
			tmpSendInt(5+(iNumOfNodeSend(1)-1)*5) = 0

			tmpSendDps(1+(iNumOfNodeSend(1)-1)*6) = root%tpChildL%tCenterv%v(1)
			tmpSendDps(2+(iNumOfNodeSend(1)-1)*6) = root%tpChildL%tCenterv%v(2)
			tmpSendDps(3+(iNumOfNodeSend(1)-1)*6) = root%tpChildL%tCenterv%v(3)
			tmpSendDps(4+(iNumOfNodeSend(1)-1)*6) = root%tpChildL%tSizev%v(1)
			tmpSendDps(5+(iNumOfNodeSend(1)-1)*6) = root%tpChildL%tSizev%v(2)
			tmpSendDps(6+(iNumOfNodeSend(1)-1)*6) = root%tpChildL%tSizev%v(3)

		else

			tmpSendLog(iNumOfNodeSend(1)) = root%tpChildR%lNodeOrLeaf

			tmpSendInt(1+(iNumOfNodeSend(1)-1)*5) = root%tpChildR%iID
			tmpSendInt(2+(iNumOfNodeSend(1)-1)*5) = root%tpChildR%iIndex
			tmpSendInt(3+(iNumOfNodeSend(1)-1)*5) = root%tpChildR%iNumDPs
			tmpSendInt(4+(iNumOfNodeSend(1)-1)*5) = 0
			tmpSendInt(5+(iNumOfNodeSend(1)-1)*5) = 0

			tmpSendDps(1+(iNumOfNodeSend(1)-1)*6) = root%tpChildR%tCenterv%v(1)
			tmpSendDps(2+(iNumOfNodeSend(1)-1)*6) = root%tpChildR%tCenterv%v(2)
			tmpSendDps(3+(iNumOfNodeSend(1)-1)*6) = root%tpChildR%tCenterv%v(3)
			tmpSendDps(4+(iNumOfNodeSend(1)-1)*6) = root%tpChildR%tSizev%v(1)
			tmpSendDps(5+(iNumOfNodeSend(1)-1)*6) = root%tpChildR%tSizev%v(2)
			tmpSendDps(6+(iNumOfNodeSend(1)-1)*6) = root%tpChildR%tSizev%v(3)

		endif
	end do

	! check to see if stop should be made.
	do ii=1,2
		if(ii==1)then
			! if it is node, check distance to see if stop
			if(root%tpChildL%lNodeOrLeaf)then
				tmpS2 = GetDist(root%tpChildL%tCenterv,tpBoxRoots(j)%p%tCenterv)
				tmpS3 = MAG(root%tpChildL%tSizev)
				if(tmpS2<=tmpS3+tmpS1+dcrit+dpMAXAveLength)then ! close, check children, set children index
					tmpSendInt(4+(tmpIndex(ii)-1)*5) = iNumOfNodeSend(1) +1
					tmpSendInt(5+(tmpIndex(ii)-1)*5) = iNumOfNodeSend(1) +2

					call DetermineGhost(root%tpChildL,j,tmpS1)
				else          ! not close, end code, modify the node to leaf
					tmpSendInt(4+(tmpIndex(ii)-1)*5) = 0
					tmpSendInt(5+(tmpIndex(ii)-1)*5) = 0
					tmpSendLog(tmpIndex(ii)) = .false.
				end if
			end if     ! if not node, just end here, tmpSendLog is right, index is right.
		else
			if(root%tpChildR%lNodeOrLeaf)then  ! if it is node, check distance to see if stop
				tmpS2 = GetDist(root%tpChildR%tCenterv,tpBoxRoots(j)%p%tCenterv)
				tmpS3 = MAG(root%tpChildR%tSizev)
				if(tmpS2<=tmpS3+tmpS1+dcrit+dpMAXAveLength)then         ! close, check children, set children index
					tmpSendInt(4+(tmpIndex(ii)-1)*5) = iNumOfNodeSend(1) +1
					tmpSendInt(5+(tmpIndex(ii)-1)*5) = iNumOfNodeSend(1) +2
					call DetermineGhost(root%tpChildR,j,tmpS1)
				else                                     ! not close, end code, modify the node to leaf
					tmpSendInt(4+(tmpIndex(ii)-1)*5) = 0
					tmpSendInt(5+(tmpIndex(ii)-1)*5) = 0
					tmpSendLog(tmpIndex(ii)) = .false.
				end if
			end if                                      ! if not node, just end here, tmpSendLog is right, index is right.
		end if
	enddo

end subroutine DetermineGhost

!**************************************************************************
!   BuildGhost: build the ghost tree with nodes from close boxes.
!**************************************************************************
recursive subroutine BuildGhost(root)
	use ZQMPI
	use vectors
	use TreeLeaf
	use MPIModule,only:iNumTNode
	use GhostTransMD,only:tmpReceLog,tmpReceInt,tmpReceDps
	implicit none

	type(tTreeLeaf),pointer::root

	integer :: i,j,k,ii,kk

	allocate(root%tpChildL)
	allocate(root%tpChildR)
	nullify(root%tpChildL%tpChildR)
	nullify(root%tpChildL%tpChildL)
	nullify(root%tpChildR%tpChildL)
	nullify(root%tpChildR%tpChildR)

	! construct children
	i = root%iChild(1)-iNumTNode   ! left index
	root%tpChildL%lNodeOrLeaf = tmpReceLog(i)
	root%tpCHildL%iID         = tmpReceInt(1+(i-1)*5)
	root%tpChildL%iIndex      = tmpReceInt(2+(i-1)*5)
	root%tpChildL%iNumDPs    = tmpReceInt(3+(i-1)*5)
	root%tpChildL%iChild(1)   = tmpReceInt(4+(i-1)*5)+iNumTNode
	root%tpChildL%iChild(2)   = tmpReceInt(5+(i-1)*5)+iNumTNode

	root%tpChildL%tCenterv%v(1) = tmpReceDps(1+(i-1)*6)
	root%tpChildL%tCenterv%v(2) = tmpReceDps(2+(i-1)*6)
	root%tpChildL%tCenterv%v(3) = tmpReceDps(3+(i-1)*6)
	root%tpChildL%tSizev%v(1)   = tmpReceDps(4+(i-1)*6)
	root%tpChildL%tSizev%v(2)   = tmpReceDps(5+(i-1)*6)
	root%tpChildL%tSizev%v(3)   = tmpReceDps(6+(i-1)*6)

	i = root%iChild(2)-iNumTNode   ! right index
	root%tpChildR%lNodeOrLeaf = tmpReceLog(i)
	root%tpCHildR%iID         = tmpReceInt(1+(i-1)*5)
	root%tpChildR%iIndex      = tmpReceInt(2+(i-1)*5)
	root%tpChildR%iNumDPs    = tmpReceInt(3+(i-1)*5)
	root%tpChildR%iChild(1)   = tmpReceInt(4+(i-1)*5)+iNumTNode
	root%tpChildR%iChild(2)   = tmpReceInt(5+(i-1)*5)+iNumTNode

	root%tpChildR%tCenterv%v(1) = tmpReceDps(1+(i-1)*6)
	root%tpChildR%tCenterv%v(2) = tmpReceDps(2+(i-1)*6)
	root%tpChildR%tCenterv%v(3) = tmpReceDps(3+(i-1)*6)
	root%tpChildR%tSizev%v(1)   = tmpReceDps(4+(i-1)*6)
	root%tpChildR%tSizev%v(2)   = tmpReceDps(5+(i-1)*6)
	root%tpChildR%tSizev%v(3)   = tmpReceDps(6+(i-1)*6)

	! check if need to stop : if node, construct child, if leaf, stop.
	if(root%tpChildL%lNodeOrLeaf)call BuildGhost(root%tpChildL)
	if(root%tpChildR%lNodeOrLeaf)call BuildGhost(root%tpChildR)

end subroutine BuildGhost
!***************************************
END MODULE TreeMD
