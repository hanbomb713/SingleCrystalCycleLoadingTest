!*******************************************************************
! Variables and structures specially defined for parallel computing.
!*******************************************************************
!*****************************************************************
!	TreeLeaf:
!		Definition of tree structures
!
!	iID: it corresponds to which processor it is on if it is a leaf. Otherwise
!		it is zero. 
!	iIndex: it is the index of every node. not zero for everyone.
!*****************************************************************
MODULE TreeLeaf
	use ZQMPI
	use vectors
	implicit none

	!6/14/06
	integer,dimension(:,:),allocatable::itIndex         !tmp index for arranging, (iNumDPs,2),1:index of smaller; 2: index of larger.

	Type tTreePointer
		type(tTreeLeaf),pointer :: p
	end Type tTreePointer

	TYPE tTreeLeaf
		!/////Tree info///////
		logical:: lNodeOrLeaf           ! true for node, false for leaf
		integer:: iID                   ! id only for leaves. node:zero. global leaves(which corresponds to processors):zero
		integer:: iIndex                 !index of all nodes(including leaves)
		integer:: iLevel
		TYPE(vector):: tCenterv        ! center of the box
		TYPE(vector):: tSizev          ! size of the box,half x,y,z
		integer,dimension(2) :: iChild
		TYPE(tTreeLeaf),pointer::tpChildL,tpChildR
		TYPE(tTreeLeaf),pointer::tParent

		!////container info///////
		integer:: iNumDPs              ! how many nodes in the leaf
		TYPE(DPPointer),pointer,dimension(:) :: tapDPs
		
		! Multipole information here
		! pointers to nodes here

	end TYPE tTreeLeaf

End MODULE TreeLeaf
!****************************************************************
! TreeTransfer:
!	Used in the process to transfer tree structures.
!****************************************************************
MODULE treeTransfer
	use ZQMPI
	implicit none

	logical,allocatable,dimension(:) :: tmpTTLogNL
	integer,allocatable,dimension(:) :: tmpTTInts  
	double precision, allocatable, dimension(:) :: tmpTTDps  ! (6,:)

END MODULE treeTransfer
!***************************************************************
!	GPTransferMD:
!		Variables used in the distribution of glide planes.
!***************************************************************
Module GPTransferMD
	use ZQMPI
	implicit none

	double precision,dimension(:),allocatable::tmpGPDps
	integer,dimension(:),allocatable::tmpGPInts

end module GPTransferMD
!***************************************************************
!	GhostTransMD:
!		variables used in the distribution of ghost tree.
!***************************************************************
Module GhostTransMD
	use ZQMPI
	use vectors
	implicit none

	integer::iNumOfNodeRece
	logical,allocatable,dimension(:)::tmpReceLog
	integer,allocatable,dimension(:)::tmpReceInt
	double precision,allocatable,dimension(:)::tmpReceDPs

	integer,dimension(:),pointer::iNumOfNodeSend
	logical,dimension(:),pointer::tmpSendLog
	integer,dimension(:),pointer::tmpSendInt
	double precision,dimension(:),pointer::tmpSendDps

End module ghostTransMD
!***************************************************************
! MPIModule:
!	variables used in parallel version, such as index arrays,
!tree structures.
!***************************************************************
module MPIMOdule
	use ZQMPI
	use TreeLeaf
	implicit none

	! 1. 
	integer,allocatable,dimension(:,:)::iaGLIndex !iaGLIndex(iNG,4): iNG->GLobal index; 1->processor;
                                                      !2->local index. 
                                                      !global index array: To store the local processor
                                                      !and local index of global indexes.
	integer::iNTotalDPs ! total number of DPs, different with INPoint, which means local points
	integer::iNTotalGPs ! total number of glide planes. If there is any change of iNPlane, the later
                            ! process will only update the new planes.
	integer::iNVLPoints !iNPoints + connection points, no ghost points

	! 2. transfer points
	integer::iNumTransP                             !number of transfer points
	integer,allocatable,dimension(:)::iaTransList   !transfer list of points from other processors

	! 3. tree variables
	integer :: iNMaxS                                             !max subsegment number in a node(leaf).
	type(tTreeLeaf),POINTER :: tTreeRoot                          !pointer to the root of global tree.
	type(tTreeLeaf),pointer :: tLocalRoot                         !pointer to the root of local tree.
	type(tTreePointer), allocatable,dimension(:):: tpBoxRoots     ! pointers to the other roots
	integer :: iCurrentNodeIndex                           ! used in master to assign index to each box.
	integer :: iNumTNode                                   ! number of total node(leaves) in the global tree
	integer,allocatable,dimension(:):: iIndexOfCloseBox
	integer::iNumCLoseBox

	!4. statistics variables
	double precision::dTotalLength, dTmpDensi,dDensity,dTmpEDensi,dTmpSDensi,dScrewDensity,dEdgeMixDensity

end module MPIModule
!************************************************************************
!	PerformanceMD:
!		Variables used in test of the code performance.
!************************************************************************
Module PerformanceMD
  use ZQMPI
  implicit none

  double precision:: dBeginTime,dENdTime,dPeriodTime,dTotalTime,dRatio1,dRatio2
  double precision:: dBuildTreeTime,dDynamicsTime, dUpdateTime
  double precision::pastTotalTime

end Module PerformanceMD
