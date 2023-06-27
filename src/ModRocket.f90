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
  real, dimension(5)               :: CurrentRocketMF
  real, dimension(3)               :: CurrentRocketPos
  real, dimension(4)               :: CRocket1, CRocket2, CRocket3, CRocket4
!!$  real, dimension(4)               :: DtRocket
  
  integer   :: nRocketLines
!!$  integer   :: CurrentRocketiLon, CurrentRocketiLat, CurrentRocketiAlt
!!$  integer   :: CurrentRocketiBlock
!!$  integer   :: r, rlon, rlat, ralt
  real      :: r
  real      :: Ha

contains

  subroutine init_mod_rocket

    if(allocated(RocketTime)) return
    allocate( &
         RocketTime(nMaxSatInputLines), &
         RocketPos(3, nMaxSatInputLines), &
         RocketMF(5, nMaxSatInputLines))

  end subroutine init_mod_rocket

end module ModRocket
