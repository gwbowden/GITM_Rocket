!----------------------------------------------------------------------
! Author: George Bowden, UNSW Canberra
! 
! Comments: Variables used to implement hard acoustic source at node on
!           GITM grid
!----------------------------------------------------------------------

module ModAcoustic

  use ModInputs, only : iCharLen_
  use ModIoUnit, only : UnitTmp_
  use ModSatellites, only: nMaxSatInputLines, nMaxSatPos

  implicit none

  integer, parameter :: dblprec = selected_real_kind(14,200)
  real, allocatable                :: AcousticPos(:,:)
  real, allocatable                :: AcousticValues(:,:)
  real (kind=dblprec), allocatable :: AcousticTime(:)
  real, dimension(3)               :: CurrentAcousticPos
  real, dimension(12)              :: CurrentAcousticValues

  integer   :: nAcousticLines

contains

  subroutine init_mod_acoustic

    if(allocated(AcousticTime)) return
    allocate( &
         AcousticTime(nMaxSatInputLines), &
         AcousticPos(3, nMaxSatInputLines), &
         AcousticValues(12, nMaxSatInputLines))

  end subroutine init_mod_acoustic

end module ModAcoustic
