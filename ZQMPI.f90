!********************************************************************
!	ZQMPI:
!		MPI variables, subroutines used in the code. General interfaces
!	provided.
!********************************************************************
Module ZQMPI
  include "mpif.h"

  ! MPI communication data  
  type ZQMPIDATA
     integer :: iTypeTag
     integer :: iDesSrc
     integer :: iNumData
     integer :: iMpiTag
     logical,dimension(:),pointer :: laLogData
     integer,dimension(:),pointer :: iaIntData
     double precision,dimension(:),pointer :: dpaDPData
     integer,pointer :: req
  end type ZQMPIDATA
  
  ! MPI processes variables
  integer :: iMyid, iIerr, iNumprocs, iStat(MPI_STATUS_SIZE) ! general var
  integer :: iMaster            ! master ID
 
  type(ZQMPIDATA),allocatable,dimension(:) :: ZQMPI_SndRcv
  integer :: iNumZQSR

  integer,parameter :: iMaxNumZQSR = 1000
  integer,parameter :: iIntegerType = 1
  integer,parameter :: iLogicalType = 2
  integer,parameter :: iDoubleType  = 3

  interface ZQMPI_Isend
     module procedure ZQMPI_Isendl
     module procedure ZQMPI_Isendi
     module procedure ZQMPI_Isendd
  end interface

  interface ZQMPI_Recv
     module procedure ZQMPI_Recvl
     module procedure ZQMPI_Recvi
     module procedure ZQMPI_Recvd
  end interface

  interface ZQMPI_Bcast
     module procedure ZQMPI_Bcastl
     module procedure ZQMPI_Bcasti
     module procedure ZQMPI_Bcastd
  end interface

CONTAINS
  !************************************************************************
  !	ZQMPI_Init:
  !		Initialize MPI.
  !************************************************************************
  subroutine ZQMPI_Init
    implicit none

    call MPI_Init(iIerr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,iMyid,iIerr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,iNumprocs,iIerr)

    allocate(ZQMPI_SndRcv(iMaxNumZQSR))
    iMaster = 0
    iNumZQSR = 0

    if(iMyid==iMaster)write(*,*)iMyid,":iNumprocs:",iNumprocs
 
  end subroutine ZQMPI_Init

  !************************************************************************
  !	ZQMPI_Finalization:
  !		Finalize MPI.
  !************************************************************************
  subroutine ZQMPI_Finalization
	implicit none
	
	call MPI_Finalize(iIerr)

   end subroutine ZQMPI_Finalization
  !************************************************************************
  !	ZQMPI_ISendl:
  !		Isend logical type.
  !************************************************************************
  subroutine ZQMPI_Isendl(tmpData,tmpNum,tmpDesSrc,tmpTag)
    implicit none
    integer :: tmpNum,tmpDesSrc,tmpTag
    logical,dimension(:),pointer :: tmpData
    integer :: i,status

    if(iNumZQSR/=0) call ZQMPI_Clean

    if(iNumZQSR<iMaxNumZQSR)then
       iNumZQSR = iNumZQSR + 1
       ZQMPI_SndRcv(iNumZQSR)%iTypeTag = iLogicalType
       ZQMPI_SndRcv(iNumZQSR)%iDesSrc  = tmpDesSrc
       ZQMPI_SndRcv(iNumZQSR)%iNumData = tmpNum
       ZQMPI_SndRcv(iNumZQSR)%iMpiTag  = tmpTag
       ZQMPI_SndRcv(iNumZQSR)%laLogData => tmpData
       allocate(ZQMPI_SndRcv(iNumZQSR)%req)

       call MPI_Isend(ZQMPI_SndRcv(iNumZQSR)%laLogData,tmpNum,MPI_LOGICAL,tmpDesSrc,tmpTAG,   &
            MPI_COMM_WORLD,ZQMPI_SndRcv(iNumZQSR)%req,iIerr)

    endif

  end subroutine ZQMPI_ISENDL

  !************************************************************************
  !	ZQMPI_ISendI:
  !		ISend integer type.
  !************************************************************************
  subroutine ZQMPI_Isendi(tmpData,tmpNum,tmpDesSrc,tmpTag)
    implicit none
    integer :: tmpNum,tmpDesSrc,tmpTag
    integer,dimension(:),pointer :: tmpData
    integer :: i,status

    if(iNumZQSR/=0)   call ZQMPI_Clean
    if(iNumZQSR<iMaxNumZQSR)then
       iNumZQSR = iNumZQSR + 1
       ZQMPI_SndRcv(iNumZQSR)%iTypeTag = iIntegerType
       ZQMPI_SndRcv(iNumZQSR)%iDesSrc  = tmpDesSrc
       ZQMPI_SndRcv(iNumZQSR)%iNumData = tmpNum
       ZQMPI_SndRcv(iNumZQSR)%iMpiTag  = tmpTag
       ZQMPI_SndRcv(iNumZQSR)%iaIntData => tmpData
       allocate(ZQMPI_SndRcv(iNumZQSR)%req)

       call MPI_Isend(ZQMPI_SndRcv(iNumZQSR)%iaIntData,tmpNum,MPI_INTEGER,tmpDesSrc,tmpTAG,   &
            MPI_COMM_WORLD,ZQMPI_SndRcv(iNumZQSR)%req,iIerr)

    endif

  end subroutine ZQMPI_ISENDI

  !************************************************************************
  !	ZQMPI_ISendD:
  !		ISend double precision type.
  !************************************************************************
  subroutine ZQMPI_Isendd(tmpData,tmpNum,tmpDesSrc,tmpTag)
    implicit none
    integer :: tmpNum,tmpDesSrc,tmpTag
    double precision,dimension(:),pointer :: tmpData
    integer :: i,status

    if(iNumZQSR/=0) call ZQMPI_Clean

    if(iNumZQSR<iMaxNumZQSR)then
       iNumZQSR = iNumZQSR + 1
       ZQMPI_SndRcv(iNumZQSR)%iTypeTag = iDoubleType
       ZQMPI_SndRcv(iNumZQSR)%iDesSrc  = tmpDesSrc
       ZQMPI_SndRcv(iNumZQSR)%iNumData = tmpNum
       ZQMPI_SndRcv(iNumZQSR)%iMpiTag  = tmpTag
       ZQMPI_SndRcv(iNumZQSR)%dpaDPData => tmpData
       allocate(ZQMPI_SndRcv(iNumZQSR)%req)

       call MPI_Isend(ZQMPI_SndRcv(iNumZQSR)%dpaDPData,tmpNum,MPI_DOUBLE_PRECISION,tmpDesSrc,tmpTAG,   &
            MPI_COMM_WORLD,ZQMPI_SndRcv(iNumZQSR)%req,iIerr)

    endif

  end subroutine ZQMPI_ISENDD

  !************************************************************************
  subroutine ZQMPI_Recvl(tmpData,tmpNum,tmpSr,tmpTag)
    implicit none
    integer :: tmpSr,tmpNum,tmpTag
    logical,dimension(tmpNum) :: tmpData

    call MPI_Recv(tmpData,tmpNum,MPI_LOGICAL,tmpSr,tmpTag,MPI_COMM_WORLD,iStat,iIerr)
    if(tmpSr==MPI_ANY_SOURCE)then
       tmpSr = iStat(MPI_SOURCE)
    endif

  end subroutine ZQMPI_RECVL

  !************************************************************************
  subroutine ZQMPI_Recvi(tmpData,tmpNum,tmpSr,tmpTag)
    implicit none
    integer :: tmpSr,tmpNum,tmpTag
    integer,dimension(tmpNum) :: tmpData

    call MPI_Recv(tmpData,tmpNum,MPI_INTEGER,tmpSr,tmpTag,MPI_COMM_WORLD,iStat,iIerr)
    if(tmpSr==MPI_ANY_SOURCE)then
       tmpSr = iStat(MPI_SOURCE)
    endif

  end subroutine ZQMPI_RECVI

  !************************************************************************
  subroutine ZQMPI_Recvd(tmpData,tmpNum,tmpSr,tmpTag)
    implicit none
    integer :: tmpSr,tmpNum,tmpTag
    double precision,dimension(tmpNum) :: tmpData

    call MPI_Recv(tmpData,tmpNum,MPI_DOUBLE_PRECISION,tmpSr,tmpTag,MPI_COMM_WORLD,iStat,iIerr)
    if(tmpSr==MPI_ANY_SOURCE)then
       tmpSr = iStat(MPI_SOURCE)
    endif

  end subroutine ZQMPI_RECVD

  !************************************************************************
  subroutine ZQMPI_Bcastl(tmpData,tmpNum,tmpSr)
    implicit none
    integer :: tmpNum,tmpSr
    logical,dimension(tmpNum) :: tmpData
    
    call MPI_Bcast(tmpData,tmpNum,MPI_LOGICAL,tmpSr,MPI_COMM_WORLD,iIerr)

  end subroutine ZQMPI_BCASTL

  !************************************************************************
  subroutine ZQMPI_Bcasti(tmpData,tmpNum,tmpSr)
    implicit none
    integer :: tmpNum,tmpSr
    integer,dimension(tmpNum) :: tmpData
    
    call MPI_Bcast(tmpData,tmpNum,MPI_INTEGER,tmpSr,MPI_COMM_WORLD,iIerr)

  end subroutine ZQMPI_BCASTI

  !************************************************************************
  subroutine ZQMPI_Bcastd(tmpData,tmpNum,tmpSr)
    implicit none
    integer :: tmpNum,tmpSr
    double precision,dimension(tmpNum) :: tmpData
    
    call MPI_Bcast(tmpData,tmpNum,MPI_DOUBLE_PRECISION,tmpSr,MPI_COMM_WORLD,iIerr)

  end subroutine ZQMPI_BCASTD

  !************************************************************************
  !	ZQMPI_Clean:
  !		Clean the Isend operation before to save memory. 
  !************************************************************************
  subroutine ZQMPI_Clean
    implicit none
    logical :: lFlag
    integer :: status
    integer :: i,j,k
    
    k = iNumZQSR

    do i=k,1,-1

       j=ZQMPI_SndRcv(i)%req
       call MPI_Test(j,lFlag,iStat,iIerr)

       if(lFlag)then         ! Finished, free memory, if not, go on to next
          if(ZQMPI_SndRcv(i)%iTypeTag==iLogicalType)then
             deallocate(ZQMPI_SndRcv(i)%laLogData,STAT=status)
          else if(ZQMPI_SndRcv(i)%iTypeTag==iIntegerType)then
             deallocate(ZQMPI_SndRcv(i)%iaIntData,STAT=status)
          else
             deallocate(ZQMPI_SndRcv(i)%dpaDPData,STAT=status)
          endif
          
          deallocate(ZQMPI_SndRcv(i)%req,STAT=status)
          
          do j=i,iNumZQSR-1    ! move those after this.

             ZQMPI_SndRcv(j)%iTypeTag = ZQMPI_SndRcv(j+1)%iTypeTag
             ZQMPI_SndRcv(j)%iDesSrc  = ZQMPI_SndRcv(j+1)%iDesSrc
             ZQMPI_SndRcv(j)%iNumData = ZQMPI_SndRcv(j+1)%iNumData
             ZQMPI_SndRcv(j)%iMpiTag  = ZQMPI_SndRcv(j+1)%iMpiTag

             if(ZQMPI_SndRcv(j+1)%iTypeTag==iLogicalType)then

                ZQMPI_SndRcv(j)%laLogData => ZQMPI_SndRcv(j+1)%laLogData
                NULLIFY(ZQMPI_SndRcv(j+1)%laLogData)

             else if(ZQMPI_SndRcv(j+1)%iTypeTag==iIntegerType)then

                ZQMPI_SndRcv(j)%iaIntData => ZQMPI_SndRcv(j+1)%iaIntData
                NULLIFY(ZQMPI_SndRcv(j+1)%iaIntData)

             else

                ZQMPI_SndRcv(j)%dpaDPData => ZQMPI_SndRcv(j+1)%dpaDPData
                NULLIFY(ZQMPI_SndRcv(j+1)%dpaDPData)

             endif
             
             ZQMPI_SndRcv(j)%req => ZQMPI_SndRcv(j+1)%req
             NULLIFY(ZQMPI_SndRcv(j+1)%req)
             !  ZQMPI_SndRcv(j)%req = ZQMPI_SndRcv(j+1)%req
          end do
          iNumZQSR = iNumZQSR - 1
       end if
    end do

    if(iNumZQSR>=iMaxNumZQSR)then
       write(*,*)"Please increase the iMaxNumZQSR!"
       stop
    endif

  end subroutine ZQMPI_Clean

end Module ZQMPI
