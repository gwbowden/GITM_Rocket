!----------------------------------------------------------------------
! Author: George Bowden, UNSW Canberra
! 
! Comments: Variables used to track rocket exhaust plumes through the 
!           GITM grid
!----------------------------------------------------------------------

module ModRocket

  use ModInputs, only : iCharLen_
  use ModIoUnit, only : UnitTmp_
  use ModSatellites, only: nMaxSatInputLines, nMaxSatPos

  implicit none

  integer, parameter :: dblprec = selected_real_kind(14,200)
  real, allocatable                :: RocketPos(:,:)
  real, allocatable                :: RocketMF(:,:)
  real (kind=dblprec), allocatable :: RocketTime(:)
  real                             :: CurrentRocketMF(5)
  real                             :: CurrentRocketPos(3)
  
  integer   :: nRocketLines
!!$  integer   :: CurrentRocketiLon, CurrentRocketiLat, CurrentRocketiAlt
!!$  integer   :: CurrentRocketiBlock
!!$  integer   :: r, rlon, rlat, ralt
  integer   :: r

contains

  subroutine init_mod_rocket

    if(allocated(RocketTime)) return
    allocate( &
         RocketTime(nMaxSatInputLines), &
         RocketPos(3, nMaxSatInputLines), &
         RocketMF(5, nMaxSatInputLines))

  end subroutine init_mod_rocket

end module ModRocket
