!----------------------------------------------------------------------
! Author: George Bowden, UNSW Canberra
! 
! Comments: Routine to read in rocket exhaust plume data in GITM
! 
! Includes: read_rocket: reads space, time and mass flow rate data from
! a rocket exhaust plume file specified in the UAM.in file.
! move_rocket: determine the location of the rocket and calculate 
! chemical source terms corresponding to the rocket at that location
! Sections of the following shamelessly taken from elsewhere in GITM
!----------------------------------------------------------------------

subroutine read_rocket(iError)

  use ModRocket
  use ModInputs, only: iDebugLevel, iCharLen_, cRocketFile
  use ModIoUnit, only : UnitTmp_
  use ModSatellites, only: nMaxSatInputLines, nMaxSatPos
  use ModGITM, only: iEast_, iNorth_, iUp_
  use ModConstants, only: pi

  implicit none

  integer, intent(out)      :: iError
  character (len=iCharLen_) :: cLine
  integer                   :: itime(7)
  real                      :: pos(3)
  real                      :: massflow(5)
  real (kind = dblprec)     :: OldTime, NewTime
  logical                   :: IsStartFound, IsFine

  call init_mod_rocket

  iError = 0
  nRocketLines = 0

  if (iDebugLevel > 2) &
          write(*,*) "Reading Rocket Exhaust Plume File : ",cRocketFile

  open(unit=UnitTmp_,file=cRocketFile,status="old",iostat = iError)
  ! Check that the unit used is correct

  if (iError /= 0) then
     write(*,*) "Error opening rocket exhaust plume file : ",cRocketFile
     return
  endif

  IsStartFound = .false.

  do while (.not. IsStartFound)
     cLine = ""
     read(UnitTmp_, *, iostat = iError) cLine
     if (iError /= 0) IsStartFound = .true.
     if (index(cline,"#START") > 0) IsStartFound = .true.
  enddo

  if (iError /= 0) then
     write(*,*) "Error finding #START in rocket exhaust plume file : ", cRocketFile
     close(UnitTmp_)
     return
  endif

  OldTime = 0.0
  IsFine = .false.
  
  do while (iError == 0)
     
     read(UnitTmp_,*,iostat=iError) itime, pos, massflow
     if (iError == 0) then
        
        IsFine = .true.
        
        ! Convert position data from degrees and km to radians and m
!!$        if (Pos(iEast_) < 0.0) Pos(iEast_) = Pos(iEast_) + 360.0
        Pos(iEast_)  = Pos(iEast_)*pi/180.0
        Pos(iNorth_) = Pos(iNorth_)*pi/180.0
        Pos(iUp_)    = Pos(iUp_)*1000.0
        call time_int_to_real(iTime, NewTime)
        
        if (NewTime /= OldTime) then
           ! Ignore repeated times
           
           nRocketLines = nRocketLines + 1

           if (nRocketLines > nMaxSatInputLines) then
              write(*,*) "Too Many Lines in rocket exhaust plume file : ",cRocketFile
              iError = 1
           else
              RocketPos(1, nRocketLines) = Pos(1)
              RocketPos(2, nRocketLines) = Pos(2)
              RocketPos(3, nRocketLines) = Pos(3)
              RocketMF(1, nRocketLines) = massflow(1)
              RocketMF(2, nRocketLines) = massflow(2)
              RocketMF(3, nRocketLines) = massflow(3)
              RocketMF(4, nRocketLines) = massflow(4)
              RocketMF(5, nRocketLines) = massflow(5)
              RocketTime(nRocketLines) = NewTime
              if (iDebugLevel > 3) then
                 write(*,*) "====> Rocket line, time, position, flow quantities: ", &
                      nRocketLines, NewTime, Pos, massflow
              end if
           endif
           
        endif

        OldTime = NewTime
        
     endif

  enddo

  close(UnitTmp_)
  if (IsFine) iError = 0

end subroutine read_rocket

subroutine move_rocket

  use ModTime, only: CurrentTime, tSimulation
  use ModRocket
  use ModGITM
  use ModConstants, only: pi
  use ModPlanet, only: R_Earth
  use ModInputs, only: iDebugLevel, sigmaLat, sigmaLon, sigmaAlt

  implicit none
  
  integer             :: iLine
!  real                :: r, rlon, rlat, ralt
  integer             :: iBlock
  
  if ( CurrentTime >= RocketTime(1) .and. &
       CurrentTime <= RocketTime(nRocketLines)) then

     if (iDebugLevel > 3) then
        write(*,*) "====>  Moving rocket"
     end if

     iLine = 1

     if (CurrentTime == RocketTime(nRocketLines)) then
        iLine = nRocketLines
     else
        do while (RocketTime(iLine+1) < CurrentTime)
           iLine = iLine + 1
        enddo
     endif

     if ((iLine > 0) .and. (iLine < nRocketLines)) then
        r = 1-(CurrentTime - RocketTime(iLine)) / &
             (RocketTime(iLine+1) - RocketTime(iLine))
        CurrentRocketPos = &
             r*RocketPos(:,iLine) &
             +(1-r)*RocketPos(:,iLine+1)
        CurrentRocketMF = &
             r*RocketMF(:,iLine) &
             +(1-r)*RocketMF(:,iLine+1)
     else
        r = 0.0
        CurrentRocketPos = RocketPos(:,iLine)
        CurrentRocketMF = RocketMF(:,iLine)
     endif

     ! Divide mass flow rate by normalising constant (for Gaussian source)
     CurrentRocketMF = CurrentRocketMF / ((pi**1.5) &
          *(sigmaLat*sigmaLon*sigmaAlt &
          *cos(CurrentRocketPos(2))*(R_Earth+CurrentRocketPos(3))**2))

!!$     if (iDebugLevel > 3) then
!!$        write(*,*) "====>  nBlocks:", nBlocks
!!$        do iBlock=1,nBlocks
!!$           write(*,*) "====>  Minimum Lon, Lat, Alt:", &
!!$                (Longitude(0,iBlock)+Longitude(1,iBlock))/2,&
!!$                (Latitude(0,iBlock)+Latitude(1,iBlock))/2,&
!!$                Altitude_GB(1, 1, 0, iBlock),&
!!$                ", for block:", iBlock
!!$           write(*,*) "====>  Maximum Lon, Lat, Alt:", &
!!$                (Longitude(nLons,iBlock)+Longitude(nLons+1,iBlock))/2,&
!!$                (Latitude(nLats,iBlock)+Latitude(nLats+1,iBlock))/2,&
!!$                Altitude_GB(1, 1, nAlts, iBlock),&
!!$                ", for block:", iBlock
!!$        enddo
!!$             
!!$     end if

!!$     ! Determine current rocket cell
!!$     call LocationIndex(CurrentRocketPos(1),CurrentRocketPos(2), &
!!$          CurrentRocketiBlock,CurrentRocketiLon,CurrentRocketiLat,rLon,rLat)
!!$     call BlockAltIndex(CurrentRocketPos(3),CurrentRocketiBlock, &
!!$          CurrentRocketiLon,CurrentRocketiLat,CurrentRocketiAlt,rAlt)
!!$
!!$     if (iDebugLevel > 3) then
!!$        write(*,*) "====>  Rocket iLon, iLat, iAlt, iBlock", &
!!$             CurrentRocketiLon, CurrentRocketiLat, CurrentRocketiAlt, CurrentRocketiBlock
!!$     end if
!!$
!!$     if ((CurrentRocketiAlt .eq. -1) .and. (rAlt .lt. -0.999)) then
!!$
!!$        CurrentRocketiLon = -1e6
!!$        CurrentRocketiLat = -1e6
!!$        CurrentRocketiAlt = -1e6
!!$        CurrentRocketiBlock = -1e6
!!$        CurrentRocketMF = (/0.0, 0.0, 0.0, 0.0, 0.0/)
!!$
!!$     else
!!$        
!!$        ! Divide mass flow rate by volume of cell (for cell based source)
!!$        CurrentRocketMF = CurrentRocketMF / &
!!$             CellVolume(CurrentRocketiLon, CurrentRocketiLat, &
!!$             CurrentRocketiAlt, CurrentRocketiBlock)
!!$        ! Alternatively divide by
!!$             (dLonDist_GB(CurrentRocketiLon, CurrentRocketiLat, &
!!$             CurrentRocketiAlt, CurrentRocketiBlock) &
!!$             * dLatDist_GB(CurrentRocketiLon, CurrentRocketiLat, &
!!$             CurrentRocketiAlt, CurrentRocketiBlock) &
!!$             * dAlt_GB(CurrentRocketiLon, CurrentRocketiLat, &
!!$             CurrentRocketiAlt, CurrentRocketiBlock))
!!$        
!!$
!!$     endif
!!$
!!$     if (iDebugLevel > 3) then
!!$        write(*,*) "====> Rocket time, position, flow quantities: ", &
!!$             CurrentTime, CurrentRocketPos, CurrentRocketMF
!!$     end if
     
  else

!!$     CurrentRocketiLon = -1e6
!!$     CurrentRocketiLat = -1e6
!!$     CurrentRocketiAlt = -1e6
!!$     CurrentRocketiBlock = -1e6
     CurrentRocketPos = RocketPos(:,1)
     CurrentRocketMF = (/0.0, 0.0, 0.0, 0.0, 0.0/)
     
  endif

end subroutine move_rocket
