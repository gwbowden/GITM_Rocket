!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!-----------------------------------------------------------------------------
! $Id: main.f90,v 1.10 2013/10/24 18:36:35 agburr Exp $
!
! Author: Aaron Ridley
!
! Comments: GITM main
!
! AGB 3/31/13: Added call for RCMR data assimilation
! AGB 10/23/13: Adapted RCMR call to new format
!-----------------------------------------------------------------------------

program GITM

  use ModInputs
  use ModTime
  use ModGITM
  use ModMpi
  use ModRCMR
  use ModSatellites, only: SatCurrentDat, SatAltDat, nRCMRSat
  use ModEUV, only: sza

  implicit none

  integer :: iBlock
  integer :: iii

  ! ------------------------------------------------------------------------
  ! initialize stuff
  ! ------------------------------------------------------------------------
  
  call init_mpi
  call start_timing("GITM")
  call delete_stop

  call init_planet
  call set_defaults

  call read_inputs(cInputFile)
  call set_inputs

  call initialize_gitm(CurrentTime)

  call write_output

  call report("Starting Main Time Loop",0)

  ! ------------------------------------------------------------------------
  ! Run for a few iterations
  ! ------------------------------------------------------------------------

  do while (CurrentTime < EndTime)

     if (UseRocketExhaust) then
        do iBlock = 1, nBlocks
           call move_rocket(iBlock)
        end do
     endif

     if (UseAcousticSource) then
        do iBlock = 1, nBlocks
           call acoustic_source(iBlock)
        end do
     endif

     call calc_pressure

     !!! We may have to split cMax and Dt calculation!!!
     if(RCMRFlag) then
        Dt = 2
     else
        Dt = FixedDt
     end if

     call calc_timestep_vertical
     if (.not. Is1D) call calc_timestep_horizontal

     if(RCMRFlag) then
        call run_RCMR
     endif

     call advance

     if (.not.IsFramework) then
        call check_stop
     endif
     
     iStep = iStep + 1

     call write_output
  end do

  ! ------------------------------------------------------------------------
  ! Finish run
  ! ------------------------------------------------------------------------

  call finalize_gitm

end program GITM

!============================================================================
! The following subroutines are here so that we can use SWMF library routines
! Also some features available in SWMF mode only require empty subroutines
! for compilation of the stand alone code.
!============================================================================

subroutine CON_stop(StringError)

  implicit none
  character (len=*), intent(in) :: StringError
  call stop_gitm(StringError)

end subroutine CON_stop

subroutine CON_set_do_test(String,DoTest,DoTestMe)
  implicit none
  character (len=*), intent(in)  :: String
  logical          , intent(out) :: DoTest, DoTestMe

  DoTest = .false.; DoTestMe = .false.

end subroutine CON_set_do_test

subroutine CON_io_unit_new(iUnit)

  implicit none
  integer, intent(in) :: iUnit

  return

end subroutine CON_io_unit_new

!---------------------------------------------------------------------------

