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

subroutine move_rocket(iBlock)

  use ModSizeGitm
  use ModTime, only: CurrentTime, tSimulation
  use ModRocket
  use ModGITM
  use ModConstants, only: pi
  use ModPlanet, only: R_Earth, &
       nSpecies, nSpeciesTotal, iH2O_, iCO2_, iOH_, iH2_, Diff0, DiffExp
!!$  use ModInputs, only: iDebugLevel, sigmaLat, sigmaLon, sigmaAlt
  use ModInputs, only: iDebugLevel, tRocket
  use ModMpi

  implicit none
  
  integer, intent(in)   :: iBlock

  integer             :: iLine
!  real                :: r, rlon, rlat, ralt
  integer             :: iSpec, jSpec, kSpec
  real                :: DRocket, DtRocket
  real                :: Dij
  real                :: N0 = 1.28E+19
  real                :: NDRocket, NDSRocket
  real                :: Hi
  real                :: TempRocket, gRocket, meanmRocket

  integer   :: CurrentRocketiLon, CurrentRocketiLat, CurrentRocketiAlt
  integer   :: CurrentRocketiBlock
  real      :: rlon, rlat, ralt

  ! CO numerical diffusion coeffecients 
  real, parameter, dimension(nSpecies) :: Diff0CO = 1.0e17 * &
       !------------------------------------------------------------------------+
       !  0       02      N2      N       NO      He      H2O     H2      CO2
       !------------------------------------------------------------------------+
       (/0.9691, 0.8458, 0.7035, 1.2134, 0.5775, 2.3485, 0.8877, 2.4966, 0.5537/)
  real, parameter, dimension(nSpecies) :: DiffExpCO = &
       !--------------------------------------------------------------------------+
       !  0      02     N2     N      NO     He     H2O    H2     CO2
       !--------------------------------------------------------------------------+
       (/0.75,  0.75,  0.75,  0.75,  0.75,  0.75,  0.75,  0.75,  0.75/)

  integer   :: iError
  real      :: LocalVar, GlobalVar
  
  if ( CurrentTime >= RocketTime(1) .and. &
       CurrentTime <= RocketTime(nRocketLines)) then

     if (iDebugLevel > 3) then
        write(*,*) "====>  Moving rocket"
     end if

     iLine = 1

     if (CurrentTime == RocketTime(nRocketLines)) then
        iLine = nRocketLines
     else
        do while (RocketTime(iLine+1) .lt. CurrentTime)
           iLine = iLine + 1
        enddo
     endif

     if ((iLine .gt. 0) .and. (iLine .lt. nRocketLines)) then
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
!!$     CurrentRocketMF = CurrentRocketMF / ((pi**1.5) &
!!$          *(sigmaLat*sigmaLon*sigmaAlt &
!!$          *cos(CurrentRocketPos(2))*(R_Earth+CurrentRocketPos(3))**2))

     ! Determine current rocket cell
     call BlockLocationIndex(CurrentRocketPos(1),CurrentRocketPos(2), &
          iBlock,CurrentRocketiLon,CurrentRocketiLat,rLon,rLat)
     call BlockAltIndex(CurrentRocketPos(3),iBlock, &
          CurrentRocketiLon,CurrentRocketiLat,CurrentRocketiAlt,rAlt)
     
     if (iDebugLevel > 3) then
        write(*,*) "====>  Rocket iLon, iLat, iAlt, rLon, rLat, rAlt, iBlock:", &
             CurrentRocketiLon, CurrentRocketiLat, CurrentRocketiAlt, &
             rLon, rLat, rAlt, iBlock
        write(*,*) "====>  LonStart, LonEnd", Longitude(-1,iBlock), Longitude(nLons+2,iBlock)
        write(*,*) "====>  LatStart, LatEnd", Latitude(-1,iBlock), Latitude(nLats+2,iBlock)
        write(*,*) "====>  AltStart, AltEnd", &
             Altitude_GB(CurrentRocketiLon, CurrentRocketiLat, -1,iBlock), &
             Altitude_GB(CurrentRocketiLon, CurrentRocketiLat, nAlts+2,iBlock)
     endif

     ! Determine constants in approximate analytic expression for gas concentration
     ! See Bernhardt 1979

!!$     if (CurrentRocketiAlt .lt. 0) then
     if ((rLon .lt. 0.0) .or. (rLat .lt. 0.0) .or. (rAlt .lt. 0.0)) then
!!$        DtRocket = 1.0e12
        Ha = 1.0e12
        Hi = 1.0e12
        CRocket1 = (/0.0, 0.0, 0.0, 0.0/)
        CRocket2 = (/0.0, 0.0, 0.0, 0.0/)
        CRocket3 = (/0.0, 0.0, 0.0, 0.0/)
        CRocket4 = (/0.0, 0.0, 0.0, 0.0/)
        CurrentRocketMF = (/0.0, 0.0, 0.0, 0.0, 0.0/)
     else
        ! Interpolate to find density at rocket location (admittedly this should be done in a subroutine)
        NDRocket = (  rLon)*(  rLat)*(  rAlt)* &
             NDensity(CurrentRocketiLon,CurrentRocketiLat,CurrentRocketiAlt,iBlock) + &
             (1-rLon)*(  rLat)*(  rAlt)* &
             NDensity(CurrentRocketiLon+1,CurrentRocketiLat  ,CurrentRocketiAlt  ,iBlock) + &
             (  rLon)*(1-rLat)*(  rAlt)* &
             NDensity(CurrentRocketiLon  ,CurrentRocketiLat+1,CurrentRocketiAlt  ,iBlock) + &
             (1-rLon)*(1-rLat)*(  rAlt)* &
             NDensity(CurrentRocketiLon+1,CurrentRocketiLat+1,CurrentRocketiAlt  ,iBlock) + &
             (  rLon)*(  rLat)*(1-rAlt)* &
             NDensity(CurrentRocketiLon  ,CurrentRocketiLat  ,CurrentRocketiAlt+1,iBlock) + &
             (1-rLon)*(  rLat)*(1-rAlt)* &
             NDensity(CurrentRocketiLon+1,CurrentRocketiLat  ,CurrentRocketiAlt+1,iBlock) + &
             (  rLon)*(1-rLat)*(1-rAlt)* &
             NDensity(CurrentRocketiLon  ,CurrentRocketiLat+1,CurrentRocketiAlt+1,iBlock) + &
             (1-rLon)*(1-rLat)*(1-rAlt)* &
             NDensity(CurrentRocketiLon+1,CurrentRocketiLat+1,CurrentRocketiAlt+1,iBlock)
        TempRocket = (  rLon)*(  rLat)*(  rAlt)* &
             Temperature(CurrentRocketiLon,CurrentRocketiLat,CurrentRocketiAlt,iBlock)* &
             TempUnit(CurrentRocketiLon  ,CurrentRocketiLat  ,CurrentRocketiAlt  ) + &
             (1-rLon)*(  rLat)*(  rAlt)* &
             Temperature(CurrentRocketiLon+1,CurrentRocketiLat  ,CurrentRocketiAlt  ,iBlock)* &
             TempUnit(CurrentRocketiLon+1,CurrentRocketiLat  ,CurrentRocketiAlt  ) + &
             (  rLon)*(1-rLat)*(  rAlt)* &
             Temperature(CurrentRocketiLon  ,CurrentRocketiLat+1,CurrentRocketiAlt  ,iBlock)* &
             TempUnit(CurrentRocketiLon  ,CurrentRocketiLat+1,CurrentRocketiAlt  ) + &
             (1-rLon)*(1-rLat)*(  rAlt)* &
             Temperature(CurrentRocketiLon+1,CurrentRocketiLat+1,CurrentRocketiAlt  ,iBlock)* &
             TempUnit(CurrentRocketiLon+1,CurrentRocketiLat+1,CurrentRocketiAlt  ) + &
             (  rLon)*(  rLat)*(1-rAlt)* &
             Temperature(CurrentRocketiLon  ,CurrentRocketiLat  ,CurrentRocketiAlt+1,iBlock)* &
             TempUnit(CurrentRocketiLon  ,CurrentRocketiLat  ,CurrentRocketiAlt+1) + &
             (1-rLon)*(  rLat)*(1-rAlt)* &
             Temperature(CurrentRocketiLon+1,CurrentRocketiLat  ,CurrentRocketiAlt+1,iBlock)* &
             TempUnit(CurrentRocketiLon+1,CurrentRocketiLat  ,CurrentRocketiAlt+1) + &
             (  rLon)*(1-rLat)*(1-rAlt)* &
             Temperature(CurrentRocketiLon  ,CurrentRocketiLat+1,CurrentRocketiAlt+1,iBlock)* &
             TempUnit(CurrentRocketiLon  ,CurrentRocketiLat+1,CurrentRocketiAlt+1) + &
             (1-rLon)*(1-rLat)*(1-rAlt)* &
             Temperature(CurrentRocketiLon+1,CurrentRocketiLat+1,CurrentRocketiAlt+1,iBlock)* &
             TempUnit(CurrentRocketiLon+1,CurrentRocketiLat+1,CurrentRocketiAlt+1)
        gRocket = (  rLon)*(  rLat)*(  rAlt)* &
             Gravity_GB(CurrentRocketiLon,CurrentRocketiLat,CurrentRocketiAlt,iBlock) + &
             (1-rLon)*(  rLat)*(  rAlt)* &
             Gravity_GB(CurrentRocketiLon+1,CurrentRocketiLat  ,CurrentRocketiAlt  ,iBlock) + &
             (  rLon)*(1-rLat)*(  rAlt)* &
             Gravity_GB(CurrentRocketiLon  ,CurrentRocketiLat+1,CurrentRocketiAlt  ,iBlock) + &
             (1-rLon)*(1-rLat)*(  rAlt)* &
             Gravity_GB(CurrentRocketiLon+1,CurrentRocketiLat+1,CurrentRocketiAlt  ,iBlock) + &
             (  rLon)*(  rLat)*(1-rAlt)* &
             Gravity_GB(CurrentRocketiLon  ,CurrentRocketiLat  ,CurrentRocketiAlt+1,iBlock) + &
             (1-rLon)*(  rLat)*(1-rAlt)* &
             Gravity_GB(CurrentRocketiLon+1,CurrentRocketiLat  ,CurrentRocketiAlt+1,iBlock) + &
             (  rLon)*(1-rLat)*(1-rAlt)* &
             Gravity_GB(CurrentRocketiLon  ,CurrentRocketiLat+1,CurrentRocketiAlt+1,iBlock) + &
             (1-rLon)*(1-rLat)*(1-rAlt)* &
             Gravity_GB(CurrentRocketiLon+1,CurrentRocketiLat+1,CurrentRocketiAlt+1,iBlock)
        meanmRocket = (  rLon)*(  rLat)*(  rAlt)* &
             MeanMajorMass(CurrentRocketiLon,CurrentRocketiLat,CurrentRocketiAlt) + &
             (1-rLon)*(  rLat)*(  rAlt)* &
             MeanMajorMass(CurrentRocketiLon+1,CurrentRocketiLat  ,CurrentRocketiAlt  ) + &
             (  rLon)*(1-rLat)*(  rAlt)* &
             MeanMajorMass(CurrentRocketiLon  ,CurrentRocketiLat+1,CurrentRocketiAlt  ) + &
             (1-rLon)*(1-rLat)*(  rAlt)* &
             MeanMajorMass(CurrentRocketiLon+1,CurrentRocketiLat+1,CurrentRocketiAlt  ) + &
             (  rLon)*(  rLat)*(1-rAlt)* &
             MeanMajorMass(CurrentRocketiLon  ,CurrentRocketiLat  ,CurrentRocketiAlt+1) + &
             (1-rLon)*(  rLat)*(1-rAlt)* &
             MeanMajorMass(CurrentRocketiLon+1,CurrentRocketiLat  ,CurrentRocketiAlt+1) + &
             (  rLon)*(1-rLat)*(1-rAlt)* &
             MeanMajorMass(CurrentRocketiLon  ,CurrentRocketiLat+1,CurrentRocketiAlt+1) + &
             (1-rLon)*(1-rLat)*(1-rAlt)* &
             MeanMajorMass(CurrentRocketiLon+1,CurrentRocketiLat+1,CurrentRocketiAlt+1)
        Ha = -Boltzmanns_Constant*TRocket/(gRocket*meanmRocket)

        if ((CurrentRocketMF(1)+CurrentRocketMF(2)+CurrentRocketMF(3)+CurrentRocketMF(4)) < 1.0e-9) then
           CurrentRocketMF(5) = 0.0
        else
           CurrentRocketMF(5) = CurrentRocketMF(5)/(Mass(iH2O_)*CurrentRocketMF(1)+Mass(iH2_)*CurrentRocketMF(2)+ &
                Mass(iCO2_)*CurrentRocketMF(3)+Mass(iCO_)*CurrentRocketMF(4))
        endif

        if (iDebugLevel > 3) then
           write(*,*) "====>  NDRocket:", NDRocket
           write(*,*) "====>  TempRocket:", TempRocket
           write(*,*) "====>  gRocket:", gRocket
           write(*,*) "====>  meanmRocket:", meanmRocket
           write(*,*) "====>  CurrentRocketMF(5):", CurrentRocketMF(5)
        endif
        
        do iSpec = 1, 4
           select case (iSpec)
           case (1)
              kSpec = iH2O_
           case (2)
              kSpec = iH2_
           case (3)
              kSpec = iCO2_
           case (4)
              kSpec = iCO_
           end select
           DRocket = 0.0
           do jSpec = 1, nSpecies
              if (kSpec .ne. jSpec) then
                 NDSRocket = (  rLon)*(  rLat)*(  rAlt)* &
                      NDensityS(CurrentRocketiLon  ,CurrentRocketiLat  ,CurrentRocketiAlt  ,jSpec,iBlock) + &
                      (1-rLon)*(  rLat)*(  rAlt)* &
                      NDensityS(CurrentRocketiLon+1,CurrentRocketiLat  ,CurrentRocketiAlt  ,jSpec,iBlock) + &
                      (  rLon)*(1-rLat)*(  rAlt)* &
                      NDensityS(CurrentRocketiLon  ,CurrentRocketiLat+1,CurrentRocketiAlt  ,jSpec,iBlock) + &
                      (1-rLon)*(1-rLat)*(  rAlt)* &
                      NDensityS(CurrentRocketiLon+1,CurrentRocketiLat+1,CurrentRocketiAlt  ,jSpec,iBlock) + &
                      (  rLon)*(  rLat)*(1-rAlt)* &
                      NDensityS(CurrentRocketiLon  ,CurrentRocketiLat  ,CurrentRocketiAlt+1,jSpec,iBlock) + &
                      (1-rLon)*(  rLat)*(1-rAlt)* &
                      NDensityS(CurrentRocketiLon+1,CurrentRocketiLat  ,CurrentRocketiAlt+1,jSpec,iBlock) + &
                      (  rLon)*(1-rLat)*(1-rAlt)* &
                      NDensityS(CurrentRocketiLon  ,CurrentRocketiLat+1,CurrentRocketiAlt+1,jSpec,iBlock) + &
                      (1-rLon)*(1-rLat)*(1-rAlt)* &
                      NDensityS(CurrentRocketiLon+1,CurrentRocketiLat+1,CurrentRocketiAlt+1,jSpec,iBlock)
                 if (iSpec .eq. 4) then
                    Dij = (1.0e-04)*&              ! Scales the Dij from cm^2/s -> m^2/s
                         (Diff0CO(jSpec)*(TempRocket**DiffExpCO(jSpec))) / &
                         (NDRocket*(1.0e-06))     ! Converts to #/cm^-3
                 else
                    Dij = (1.0e-04)*&              ! Scales the Dij from cm^2/s -> m^2/s
                         (Diff0(kSpec,jSpec)*(TempRocket**DiffExp(kSpec,jSpec))) / &
                         (NDRocket*(1.0e-06))     ! Converts to #/cm^-3
                 endif
                 if (iDebugLevel > 3) then
                    write(*,*) "====>  D(", jSpec, ",", kSpec, "):", Dij
                    write(*,*) "====>  NDSRocket:", NDSRocket
                 endif
                 DRocket = DRocket+(NDSRocket/Dij)
              endif
           enddo
           DRocket = NDRocket/DRocket
           NDSRocket = (  rLon)*(  rLat)*(  rAlt)* &
                NDensityS(CurrentRocketiLon  ,CurrentRocketiLat  ,CurrentRocketiAlt  ,kSpec,iBlock) + &
                (1-rLon)*(  rLat)*(  rAlt)* &
                NDensityS(CurrentRocketiLon+1,CurrentRocketiLat  ,CurrentRocketiAlt  ,kSpec,iBlock) + &
                (  rLon)*(1-rLat)*(  rAlt)* &
                NDensityS(CurrentRocketiLon  ,CurrentRocketiLat+1,CurrentRocketiAlt  ,kSpec,iBlock) + &
                (1-rLon)*(1-rLat)*(  rAlt)* &
                NDensityS(CurrentRocketiLon+1,CurrentRocketiLat+1,CurrentRocketiAlt  ,kSpec,iBlock) + &
                (  rLon)*(  rLat)*(1-rAlt)* &
                NDensityS(CurrentRocketiLon  ,CurrentRocketiLat  ,CurrentRocketiAlt+1,kSpec,iBlock) + &
                (1-rLon)*(  rLat)*(1-rAlt)* &
                NDensityS(CurrentRocketiLon+1,CurrentRocketiLat  ,CurrentRocketiAlt+1,kSpec,iBlock) + &
                (  rLon)*(1-rLat)*(1-rAlt)* &
                NDensityS(CurrentRocketiLon  ,CurrentRocketiLat+1,CurrentRocketiAlt+1,kSpec,iBlock) + &
                (1-rLon)*(1-rLat)*(1-rAlt)* &
                NDensityS(CurrentRocketiLon+1,CurrentRocketiLat+1,CurrentRocketiAlt+1,kSpec,iBlock)
           DRocket = DRocket*(1-(NDSRocket/NDRocket))
           DtRocket = tRocket*DRocket
           CurrentRocketMF(iSpec) = CurrentRocketMF(iSpec) / ((4*pi*DtRocket)**1.5)
           Hi = -Boltzmanns_Constant*TRocket/(gRocket*Mass(kSpec))
           CRocket1(iSpec) = (0.75/Ha)+(0.5/Hi)
           CRocket2(iSpec) = Ha**2/DtRocket
           CRocket3(iSpec) = 1.0/(4.0*DtRocket)
           CRocket4(iSpec) = (1.0/Ha-1.0/Hi)**2*DtRocket/4.0
        enddo
        
     endif
     
     LocalVar = Ha
     call MPI_ALLREDUCE(LocalVar, GlobalVar, 1, MPI_REAL, MPI_MIN, &
          iCommGITM, iError)
     Ha = GlobalVar
     ! Find correct value of constants in approximate analytic expression
     do iSpec = 1, 4
        LocalVar = CRocket1(iSpec)
        call MPI_ALLREDUCE(LocalVar, GlobalVar, 1, MPI_REAL, MPI_MAX, &
          iCommGITM, iError)
        CRocket1(iSpec) = GlobalVar
        LocalVar = CRocket2(iSpec)
        call MPI_ALLREDUCE(LocalVar, GlobalVar, 1, MPI_REAL, MPI_MAX, &
          iCommGITM, iError)
        CRocket2(iSpec) = GlobalVar
        LocalVar = CRocket3(iSpec)
        call MPI_ALLREDUCE(LocalVar, GlobalVar, 1, MPI_REAL, MPI_MAX, &
          iCommGITM, iError)
        CRocket3(iSpec) = GlobalVar
        LocalVar = CRocket4(iSpec)
        call MPI_ALLREDUCE(LocalVar, GlobalVar, 1, MPI_REAL, MPI_MAX, &
          iCommGITM, iError)
        CRocket4(iSpec) = GlobalVar
     enddo
     
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
     
!!$  elseif (CurrentTime < RocketTime(1))
!!$  elseif (CurrentTime > RocketTime(nRocketLines))
  else

     if (iDebugLevel > 3) then
        write(*,*) "====>  CurrentTime outside rocket file times"
     end if
     CurrentRocketPos = RocketPos(:,1)
     CurrentRocketMF = (/0.0, 0.0, 0.0, 0.0, 0.0/)
     CRocket1 = (/0.0, 0.0, 0.0, 0.0/)
     CRocket2 = (/0.0, 0.0, 0.0, 0.0/)
     CRocket3 = (/0.0, 0.0, 0.0, 0.0/)
     CRocket4 = (/0.0, 0.0, 0.0, 0.0/)
     Ha = 1.0e12
     Hi = 1.0e12
     
  endif

  if (iDebugLevel > 3) then
     write(*,*) "====>  CurrentRocketPos:", CurrentRocketPos(1), CurrentRocketPos(2), CurrentRocketPos(3)
  end if
  do iSpec = 1, 4
     if (iDebugLevel > 3) then
        write(*,*) "====>  iSpec:", iSpec
        write(*,*) "====>  CRocket:", CRocket1(iSpec), CRocket2(iSpec), CRocket3(iSpec), CRocket4(iSpec)
        write(*,*) "====>  Hi, Ha:", Hi, Ha
        write(*,*) "====>  CurrentRocketMF:", CurrentRocketMF(iSpec)
        write(*,*) "====>  DRocket:", DRocket
     end if
  enddo

end subroutine move_rocket
