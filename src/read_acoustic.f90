!----------------------------------------------------------------------
! Author: George Bowden, UNSW Canberra
! 
! Comments: Routine to read in acoustic data to GITM
! 
! Includes: read_acoustic: reads data from acoustic source input file
! specified in the UAM.in file.
! acoustic_source: modify source node temperature and species number
! density values.
! Sections of the following shamelessly taken from elsewhere in GITM
!----------------------------------------------------------------------

subroutine read_acoustic(iError)

  use ModAcoustic
  use ModInputs, only: iDebugLevel, iCharLen_, cAcousticFile
  use ModIoUnit, only : UnitTmp_
  use ModSatellites, only: nMaxSatInputLines, nMaxSatPos
  use ModGITM, only: iEast_, iNorth_, iUp_
  use ModConstants, only: pi

  implicit none

  integer, intent(out)      :: iError
  character (len=iCharLen_) :: cLine
  integer                   :: itime(7)
  real                      :: pos(3)
  real                      :: nodevals(12)
  real (kind = dblprec)     :: OldTime, NewTime
  logical                   :: IsStartFound, IsFine

  call init_mod_acoustic

  iError = 0
  nAcousticLines = 0

  if (iDebugLevel > 2) &
       write(*,*) "Reading acoustic source file : ",cAcousticFile

  open(unit=UnitTmp_,file=cAcousticFile,status="old",iostat = iError)
  ! Check that the unit used is correct

  if (iError /= 0) then
     write(*,*) "Error opening acoustic source file : ",cAcousticFile
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
     write(*,*) "Error finding #START in acoustic source file : ", cAcousticFile
     close(UnitTmp_)
     return
  endif

  OldTime = 0.0
  IsFine = .false.

  do while (iError == 0)

     read(UnitTmp_,*,iostat=iError) itime, pos, nodevals
     if (iError == 0) then

        IsFine = .true.

        ! Convert position data from degrees and km to radians and m
        Pos(iEast_)  = Pos(iEast_)*pi/180.0
        Pos(iNorth_) = Pos(iNorth_)*pi/180.0
        Pos(iUp_)    = Pos(iUp_)*1000.0
        call time_int_to_real(itime, NewTime)

        if (NewTime /= OldTime) then
           ! Ignore repeated times

           nAcousticLines =  nAcousticLines + 1

           if (nAcousticLines > nMaxSatInputLines) then
              write(*,*) "Too Many Lines in acoustic source file :", cAcousticFile
           else
              AcousticPos(1, nAcousticLines) = pos(1)
              AcousticPos(2, nAcousticLines) = pos(2)
              AcousticPos(3, nAcousticLines) = pos(3)
              AcousticValues(1, nAcousticLines) = nodevals(1)
              AcousticValues(2, nAcousticLines) = nodevals(2)
              AcousticValues(3, nAcousticLines) = nodevals(3)
              AcousticValues(4, nAcousticLines) = nodevals(4)
              AcousticValues(5, nAcousticLines) = nodevals(5)
              AcousticValues(6, nAcousticLines) = nodevals(6)
              AcousticValues(7, nAcousticLines) = nodevals(7)
              AcousticValues(8, nAcousticLines) = nodevals(8)
              AcousticValues(9, nAcousticLines) = nodevals(9)
              AcousticValues(10, nAcousticLines) = nodevals(10)
              AcousticValues(11, nAcousticLines) = nodevals(11)
              AcousticValues(12, nAcousticLines) = nodevals(12)
              AcousticTime(nAcousticLines) = NewTime
              if (iDebugLevel > 3) then
                 write(*,*) "====> Acoustic line, time, position, node values: ", &
                      nAcousticLines, NewTime, pos, nodevals
              end if
              
           endif
           
        endif

        OldTime = NewTime
        
     endif

  enddo

  close(UnitTmp_)
  if (IsFine) iError = 0

end subroutine read_acoustic



subroutine acoustic_source(iBlock)

  use ModSizeGitm
  use ModTime, only: CurrentTime, tSimulation
  use ModAcoustic
  use ModGITM
  use ModPlanet, only: RBody
  use ModInputs, only: iDebugLevel, UseSoftSource
  use ModMpi
  use ModConstants, only: pi, Univ_Gas_Constant

  implicit none

  integer, intent(in)   :: iBlock

  integer             :: iLine
  integer             :: iilon, iilat, iialt

  integer   :: CurrentAcousticiLon, CurrentAcousticiLat, CurrentAcousticiAlt
  integer   :: CurrentAcousticiBlock
  real      :: r
  real      :: rlon, rlat, ralt
  real      :: volel

  integer   :: iError

  if ( CurrentTime >= AcousticTime(1) .and. &
       CurrentTime <= AcousticTime(nAcousticLines)) then

     iLine = 1

     if (CurrentTime == AcousticTime(nAcousticLines)) then
        iLine = nAcousticLines
     else
        do while (AcousticTime(iLine+1) .lt. CurrentTime)
           iLine = iLine + 1
        enddo
     endif

     if ((iLine .gt. 0) .and. (iLine .lt. nAcousticLines)) then
        r = 1-(CurrentTime - AcousticTime(iLine)) / &
             (AcousticTime(iLine+1) - AcousticTime(iLine))
        CurrentAcousticPos = &
             r*AcousticPos(:,iLine) &
             +(1-r)*AcousticPos(:,iLine+1)
        CurrentAcousticValues = &
             r*AcousticValues(:,iLine) &
             +(1-r)*AcousticValues(:,iLine+1)
     else
        r = 0.0
        CurrentAcousticPos = AcousticPos(:,iLine)
        CurrentAcousticValues = AcousticValues(:,iLine)
     endif

     ! Determine current rocket cell
     call BlockLocationIndex(CurrentAcousticPos(1),CurrentAcousticPos(2), &
          iBlock,CurrentAcousticiLon,CurrentAcousticiLat,rLon,rLat)
     call BlockAltIndex(CurrentAcousticPos(3),iBlock, &
          CurrentAcousticiLon,CurrentAcousticiLat,CurrentAcousticiAlt,rAlt)

     if (iDebugLevel > 3) then
        write(*,*) "====>  Acoustic iLon, iLat, iAlt, rLon, rLat, rAlt, iBlock:", &
             CurrentAcousticiLon, CurrentAcousticiLat, CurrentAcousticiAlt, &
             rLon, rLat, rAlt, iBlock
        write(*,*) "====>  LonStart, LonEnd", Longitude(-1,iBlock), Longitude(nLons+2,iBlock)
        write(*,*) "====>  LatStart, LatEnd", Latitude(-1,iBlock), Latitude(nLats+2,iBlock)
        write(*,*) "====>  AltStart, AltEnd", &
             Altitude_GB(CurrentAcousticiLon, CurrentAcousticiLat, -1,iBlock), &
             Altitude_GB(CurrentAcousticiLon, CurrentAcousticiLat, nAlts+2,iBlock)
     endif

     if (UseSoftSource) then
        ! Soft source specifies rate of addition of mass and thermal energy at a
        ! particular grid point
        if ((rLon > -0.5) .and. (rLat > -0.5)) then
           iilon = CurrentAcousticiLon+nint(rLon)
           iilat = CurrentAcousticiLat+nint(rLat)
           iialt = CurrentAcousticiAlt+nint(rAlt)
           volel = CellVolume(iilon,iilat,iialt,iBlock)
           Temperature(iilon,iilat,iialt,iBlock) = &
                Temperature(iilon,iilat,iialt,iBlock) + &
                (Gamma(iilon,iilat,iialt,iBlock) - 1)* &
                CurrentAcousticValues(1)*dt/volel &
                / NDensity(iilon,iilat,iialt,iBlock) &
                / TempUnit(iilon,iilat,iialt) / Univ_Gas_Constant
           NDensityS(iilon,iilat,iialt,iCO2_,iBlock) = &
                NDensityS(iilon,iilat,iialt,iCO2_,iBlock) + &
                CurrentAcousticValues(2)*dt/volel
           NDensityS(iilon,iilat,iialt,iH_,iBlock) = &
                NDensityS(iilon,iilat,iialt,iH_,iBlock) + &
                CurrentAcousticValues(3)*dt/volel
           NDensityS(iilon,iilat,iialt,iHe_,iBlock) = &
                NDensityS(iilon,iilat,iialt,iHe_,iBlock) + &
                CurrentAcousticValues(4)*dt/volel
           NDensityS(iilon,iilat,iialt,iN2_,iBlock) = &
                NDensityS(iilon,iilat,iialt,iN2_,iBlock) + &
                CurrentAcousticValues(5)*dt/volel
           NDensityS(iilon,iilat,iialt,iN_2D_,iBlock) = &
                NDensityS(iilon,iilat,iialt,iN_2D_,iBlock) + &
                CurrentAcousticValues(6)*dt/volel
           NDensityS(iilon,iilat,iialt,iN_2P_,iBlock) = &
                NDensityS(iilon,iilat,iialt,iN_2P_,iBlock) + &
                CurrentAcousticValues(7)*dt/volel
           NDensityS(iilon,iilat,iialt,iN_4S_,iBlock) = &
                NDensityS(iilon,iilat,iialt,iN_4S_,iBlock) + &
                CurrentAcousticValues(8)*dt/volel
           NDensityS(iilon,iilat,iialt,iNO_,iBlock) = &
                NDensityS(iilon,iilat,iialt,iNO_,iBlock) + &
                CurrentAcousticValues(9)*dt/volel
           NDensityS(iilon,iilat,iialt,iO2_,iBlock) = &
                NDensityS(iilon,iilat,iialt,iO2_,iBlock) + &
                CurrentAcousticValues(10)*dt/volel
           NDensityS(iilon,iilat,iialt,iO_1D_,iBlock) = &
                NDensityS(iilon,iilat,iialt,iO_1D_,iBlock) + &
                CurrentAcousticValues(11)*dt/volel
           NDensityS(iilon,iilat,iialt,iO_3P_,iBlock) = &
                NDensityS(iilon,iilat,iialt,iO_3P_,iBlock) + &
                CurrentAcousticValues(12)*dt/volel
        endif
     else
        ! Hard source specifies density and temperature at particular grid point
        if ((rLon > -0.5) .and. (rLat > -0.5)) then
           iilon = CurrentAcousticiLon+nint(rLon)
           iilat = CurrentAcousticiLat+nint(rLat)
           iialt = CurrentAcousticiAlt+nint(rAlt)
           Temperature(iilon,iilat,iialt,iBlock) = CurrentAcousticValues(1)
           NDensityS(iilon,iilat,iialt,iCO2_,iBlock) = CurrentAcousticValues(2)
           NDensityS(iilon,iilat,iialt,iH_,iBlock) = CurrentAcousticValues(3)
           NDensityS(iilon,iilat,iialt,iHe_,iBlock) = CurrentAcousticValues(4)
           NDensityS(iilon,iilat,iialt,iN2_,iBlock) = CurrentAcousticValues(5)
           NDensityS(iilon,iilat,iialt,iN_2D_,iBlock) = CurrentAcousticValues(6)
           NDensityS(iilon,iilat,iialt,iN_2P_,iBlock) = CurrentAcousticValues(7)
           NDensityS(iilon,iilat,iialt,iN_4S_,iBlock) = CurrentAcousticValues(8)
           NDensityS(iilon,iilat,iialt,iNO_,iBlock) = CurrentAcousticValues(9)
           NDensityS(iilon,iilat,iialt,iO2_,iBlock) = CurrentAcousticValues(10)
           NDensityS(iilon,iilat,iialt,iO_1D_,iBlock) = CurrentAcousticValues(11)
           NDensityS(iilon,iilat,iialt,iO_3P_,iBlock) = CurrentAcousticValues(12)
        endif
     endif
     
  endif

end subroutine acoustic_source
