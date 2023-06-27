!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModPlanet

  use ModConstants
  use ModOrbital
  use ModSizeGITM, only: nAlts

  implicit none

  integer, parameter :: iO_3P_   = 1
  integer, parameter :: iO2_     = 2
  integer, parameter :: iN2_     = 3
  integer, parameter :: iN_4S_   = 4
  integer, parameter :: iNO_     = 5
  integer, parameter :: iHe_     = 6
  integer, parameter :: iH2O_    = 7
  integer, parameter :: iH2_     = 8
  integer, parameter :: iCO2_    = 9
  
  integer, parameter :: nSpecies = 9

  integer, parameter :: iN_2D_ = 10
  integer, parameter :: iN_2P_ = 11
  integer, parameter :: iH_    = 12
!  integer, parameter :: iAr_  = 10
!  integer, parameter :: iCO2_  = 13
!  integer, parameter :: iO_1D_ = 14
!  integer, parameter :: iCO_   = 15
!  integer, parameter :: iOH_   = 16
  integer, parameter :: iO_1D_ = 13
  integer, parameter :: iCO_   = 14
  integer, parameter :: iOH_   = 15
  integer, parameter :: nSpeciesTotal = 15

  integer, parameter  :: iO_4SP_ = 1
  integer, parameter  :: iO2P_   = 2
  integer, parameter  :: iN2P_   = 3
  integer, parameter  :: iNP_    = 4
  integer, parameter  :: iNOP_   = 5
  integer, parameter  :: iO_2DP_ = 6
  integer, parameter  :: iO_2PP_ = 7
  integer, parameter  :: iHP_    = 8
  integer, parameter  :: iHeP_   = 9
  integer, parameter  :: ie_     = 10
  integer, parameter  :: nIons   = ie_
  integer, parameter  :: nIonsAdvect = 1
  integer, parameter  :: nSpeciesAll = nSpeciesTotal + nIons - 1

  character (len=20) :: cSpecies(nSpeciesTotal)
  character (len=20) :: cIons(nIons)

  real :: Mass(nSpeciesTotal), MassI(nIons)

  real :: Vibration(nSpeciesTotal)

  integer, parameter :: iE2470_ = 1
  integer, parameter :: iE7320_ = 2
  integer, parameter :: iE3726_ = 3
  integer, parameter :: iE5200_ = 4
  integer, parameter :: iE10400_ = 5
  integer, parameter :: iE6300_ = 6
  integer, parameter :: iE6364_ = 7

  integer, parameter :: nEmissions = 10

  integer, parameter :: i3371_ = 1
  integer, parameter :: i4278_ = 2
  integer, parameter :: i5200_ = 3
  integer, parameter :: i5577_ = 4
  integer, parameter :: i6300_ = 5
  integer, parameter :: i7320_ = 6
  integer, parameter :: i10400_ = 7
  integer, parameter :: i3466_ = 8
  integer, parameter :: i7774_ = 9
  integer, parameter :: i8446_ = 10
  integer, parameter :: i3726_ = 11

  real, parameter :: GC_Earth               = 9.8                    ! m/s^2
  real, parameter :: RP_Earth               = 24.0*3600.0            ! seconds
  real, parameter :: R_Earth                = 6372.0*1000.0          ! meters
  real, parameter :: DP_Earth               = -31100.0e-9            ! nT

  real, parameter :: Gravitational_Constant = GC_Earth
  real, parameter :: Rotation_Period        = RP_Earth
  real, parameter :: RBody                  = R_Earth
  real, parameter :: DipoleStrength         = DP_Earth

  real, parameter :: OMEGABody              = 2.00*pi/Rotation_Period  ! rad/s

  real, parameter :: HoursPerDay = Rotation_Period / 3600.0
  real, parameter :: Tilt = 23.5

  real, parameter :: DaysPerYear = 365.25
  real, parameter :: SecondsPerYear = DaysPerYear * Rotation_Period

  integer, parameter :: iVernalYear   = 1999
  integer, parameter :: iVernalMonth  =    3
  integer, parameter :: iVernalDay    =   21
  integer, parameter :: iVernalHour   =    0
  integer, parameter :: iVernalMinute =    0
  integer, parameter :: iVernalSecond =    0

  ! Old orbital parameters
 !real, parameter :: SunOrbit_A = 1.000110
 !real, parameter :: SunOrbit_B = 0.034221
 !real, parameter :: SunOrbit_C = 0.001280
 !real, parameter :: SunOrbit_D = 0.000719
 !real, parameter :: SunOrbit_E = 0.000077

  !New Orbital Parameters
  !A: semi-major axis in AU
  !B: eccentricity
  !C: Longitude of perihelion
  !D: Mean Longitude
  !E: For calulating actual Longitude
 real, parameter :: SunOrbit_A = 1.0000001124
 real, parameter :: SunOrbit_B = 0.0167
 real, parameter :: SunOrbit_C = 102.94719
 real, parameter :: SunOrbit_D = 100.46435
 real, parameter :: SunOrbit_E = 129597740.63

 real :: semimajoraxis_0 = semimajor_Earth
 real :: eccentricity_0 = eccentricity_Earth
 real :: inclination_0 = inclination_Earth
 real :: longitudePerihelion_0 = longitudePerihelion_Earth
 real :: longitudeNode_0 = longitudeNode_Earth
 real :: meanLongitude_0 = meanLongitude_Earth
 real :: semimajordot = semimajordot_Earth
 real :: eccentricitydot = eccentricitydot_Earth
 real :: inclinationdot = inclinationdot_Earth
 real :: longitudePeriheliondot = longitudePeriheliondot_Earth
 real :: longitudeNodedot = longitudeNodedot_Earth
 real :: meanLongitudedot = meanLongitudedot_Earth


  !Used as a damping term in Vertical solver.
  real :: VertTau(nAlts)

  logical :: IsEarth = .true.
  logical :: IsMars = .false.
  logical :: IsTitan = .false.
  logical :: NonMagnetic = .false.
  real, parameter :: PlanetNum = 0.03

  character (len=10) :: cPlanet = "Earth"

  integer, parameter :: nEmissionWavelengths = 20
  integer, parameter :: nPhotoBins = 190


  ! These are for the neutral friction routine...

  ! These are the numerical coefficients in Table 1 in m^2 instead of cm^2
  ! JMB: Updated the Diff0 and DiffExp to be dynamic use (nspecies, nspecies)
  ! in the shape of the array

!  ! 5-Species version (Omitting Helium)
!  real, parameter, dimension(nSpecies, nSpecies) :: Diff0 = 1.0e17 * reshape( (/ &
!       !   0       02      N2       N      NO
!       !--------------------------------------+
!       0.000,   0.969,  0.969,  0.969,  0.715,&              ! O
!       0.969,   0.000,  0.715,  0.969,  0.715,&              ! O2
!       0.969,   0.715,  0.000,  0.969,  0.527,&              ! N2
!       0.969,   0.969,  0.969,  0.000,  0.969, &             ! N
!       0.715,   0.715,  0.527,  0.969,  0.000 /), &
!       (/nSpecies,nSpecies/) )  ! NO

!  ! These are the exponents
!  real, parameter, dimension(nSpecies, nSpecies) :: DiffExp = reshape( (/ &
!       !   0      02      N2      N     NO
!       !---------------------------------+
!       0.000,  0.774,  0.774, 0.774,  0.750, &              ! O
!       0.774,  0.000,  0.750, 0.774,  0.750, &              ! O2
!       0.774,  0.750,  0.000, 0.774,  0.810, &              ! N2
!       0.774,  0.774,  0.774, 0.000,  0.774, &              ! N
!       0.750,  0.750,  0.810, 0.774,  0.000  /), &
!       (/nSpecies,nSpecies/) )  ! NO

! Unknowns!!!!
!  real, parameter, dimension(5, 5) :: Diff0 = 1.0e17 * reshape( (/ &
!       !   0       02      N2       N      NO
!       !--------------------------------------+
!       0.000,   0.969,  0.969,  0.???,  0.???,&            ! O
!       0.969,   0.000,  0.715,  0.???,  0.???,&            ! O2
!       0.969,   0.715,  0.000,  0.969,  0.527,&            ! N2
!       0.???,   0.???,  0.969,  0.000,  0.???, & ! N
!       0.???,   0.???,  0.527,  0.???,  0.000 /), (/5,5/) )  ! NO
!
!  ! These are the exponents
!  real, parameter, dimension(5, 5) :: DiffExp = reshape( (/ &
!       !   0      02      N2      N     NO
!       !---------------------------------+
!       0.000,  0.774,  0.774, 0.???,  0.???, &              ! O
!       0.774,  0.000,  0.750, 0.???,  0.???, &              ! O2
!       0.774,  0.750,  0.000, 0.774,  0.810, &              ! N2
!       0.???,  0.???,  0.774, 0.000,  0.???, &              ! N
!       0.???,  0.???,  0.810, 0.???,  0.000  /), (/5,5/) )  ! NO

! Colegrove:
!
! real, parameter, dimension(5, 5) :: Diff0 = 1.0e4 * reshape( (/ &
!       !   0      02      N2       N      NO
!       !---------------------------------+
!       0.000,  0.260,  0.260,  0.300,  0.181, &             !  O
!       0.260,  0.000,  0.181,  0.220,  0.181, &             ! O2
!       0.260,  0.181,  0.000,  0.220,  0.181, &             ! N2
!       0.300,  0.220,  0.220,  0.000,  0.181, &             !  N
!       0.181,  0.181,  0.181,  0.181,  0.000 /), (/5,5/) )  ! NO
!
!  ! These are the exponents
!  real, parameter, dimension(5, 5) :: DiffExp = reshape( (/ &
!       ! 0      02     N2
!       !---------------------------------+
!       0.00,  0.75,  0.75, 0.75,  0.75, &             ! O
!       0.75,  0.00,  0.75, 0.75,  0.75, &             ! O2
!       0.75,  0.75,  0.00, 0.75,  0.75, &             ! N2
!       0.75,  0.75,  0.75, 0.00,  0.75, &            !N
!       0.75,  0.75,  0.75, 0.75,  0.0  /), (/5,5/) )  ! NO

!  real, parameter, dimension(4, 4) :: Diff0 = 1.0e4 * reshape( (/ &
!       ! 0      02     N2      N     NO
!       !---------------------------------+
!       0.00,  0.260, 0.260, 0.300, &            ! O
!       0.26,  0.000, 0.181, 0.220, &            ! O2
!       0.26,  0.181, 0.000, 0.220, &            ! N2
!       0.30,  0.220, 0.220, 0.000 /), (/4,4/) )  ! N
!
!  ! These are the exponents
!  real, parameter, dimension(4, 4) :: DiffExp = reshape( (/ &
!       ! 0      02     N2
!       !---------------------------------+
!       0.00,  0.75,  0.75, 0.75, &             ! O
!       0.75,  0.00,  0.75, 0.75, &             ! O2
!       0.75,  0.75,  0.00, 0.75, &             ! N2
!       0.75,  0.75,  0.75, 0.00 /), (/4,4/) )  ! N

!!!!!!!  real,  AltMinIono=100.0 ! in km

  ! JMB:  Updated 06/24/2016
  ! 5-Species version (Omitting Helium)
  ! Updated the N2-O based upon Massman [1998] recommended values
  ! Assumed N-He and  NO-He were the same as N2-He (just a guess)
!!$  real, parameter, dimension(nSpecies, nSpecies) :: Diff0 = 1.0e17 * reshape( (/ &
!!$       !-------------------------------------------+
!!$       !   0       02      N2    N      NO     He
!!$       !-------------------------------------------+
!!$       0.000,   0.969,  0.969,  0.969,  0.715, 3.440, & ! O
!!$       0.969,   0.000,  0.715,  0.969,  0.715, 3.208, & ! O2
!!$       0.969,   0.715,  0.000,  0.969,  0.527, 2.939, & ! N2
!!$       0.969,   0.969,  0.969,  0.000,  0.969, 2.939, & ! N
!!$       0.715,   0.715,  0.527,  0.969,  0.000, 2.939, & ! NO
!!$       3.440,   3.208,  2.939,  2.939,  2.939, 0.000  /), & !He
!!$       (/nSpecies,nSpecies/) )
!!$
!!$
!!$  ! These are the exponents
!!$  real, parameter, dimension(nSpecies, nSpecies) :: DiffExp = reshape( (/ &
!!$       !------------------------------------------+
!!$       !   0      02   N2     N       NO      He
!!$       !------------------------------------------+
!!$       0.000,  0.774,  0.774, 0.774,  0.750, 0.749, &      ! O
!!$       0.774,  0.000,  0.750, 0.774,  0.750, 0.710, &      ! O2
!!$       0.774,  0.750,  0.000, 0.774,  0.810, 0.718, &      ! N2
!!$       0.774,  0.774,  0.774, 0.000,  0.774, 0.718, &      ! N
!!$       0.750,  0.750,  0.810, 0.774,  0.000, 0.718, &      ! NO
!!$       0.749,  0.710,  0.718, 0.718,  0.718, 0.000  /), &  ! He
!!$       (/nSpecies,nSpecies/) )

  real, parameter, dimension(nSpecies, nSpecies) :: Diff0 = 1.0e17 * reshape( (/ &
       !------------------------------------------------------------------+
       !   0       02      N2    N      NO     He     H2O     H2     CO2
       !------------------------------------------------------------------+
       0.000,   0.969,  0.969,  0.969,  0.715, 3.440, 1.2884, 3.8504, 0.7527, & ! O
       0.969,   0.000,  0.715,  0.969,  0.715, 3.208, 1.0569, 2.6921, 0.6963, & ! O2
       0.969,   0.715,  0.000,  0.969,  0.527, 2.939, 0.9049, 2.5495, 0.5632, & ! N2
       0.969,   0.969,  0.969,  0.000,  0.969, 2.939, 1.5414, 3.9523, 0.9847, & ! N
       0.715,   0.715,  0.527,  0.969,  0.000, 2.939, 0.7371, 2.0647, 0.4653, & ! NO
       3.440,   3.208,  2.939,  2.939,  2.939, 0.000, 2.8750, 5.5936, 1.9622, & ! He
       1.2884, 1.0569, 0.9049, 1.5414, 0.7371, 2.8750, 0.0000, 2.9556, 0.7227, & ! H2O
       3.8504, 2.6921, 2.5495, 3.9523, 2.0647, 5.5936, 2.9556, 0.0000, 2.1424, & ! H2
       0.7527, 0.6963, 0.5632, 0.9847, 0.4653, 1.9622, 0.7227, 2.1424, 0.0000/), & ! CO2
       (/nSpecies,nSpecies/) )


  ! These are the exponents
  real, parameter, dimension(nSpecies, nSpecies) :: DiffExp = reshape( (/ &
       !-----------------------------------------------------------------+
       !   0      02   N2     N       NO      He     H2O     H2     CO2
       !-----------------------------------------------------------------+
       0.000,  0.774,  0.774, 0.774,  0.750, 0.749, 0.75, 0.75, 0.75, &      ! O
       0.774,  0.000,  0.750, 0.774,  0.750, 0.710, 0.75, 0.75, 0.75, &      ! O2
       0.774,  0.750,  0.000, 0.774,  0.810, 0.718, 0.75, 0.75, 0.75, &      ! N2
       0.774,  0.774,  0.774, 0.000,  0.774, 0.718, 0.75, 0.75, 0.75, &      ! N
       0.750,  0.750,  0.810, 0.774,  0.000, 0.718, 0.75, 0.75, 0.75, &      ! NO
       0.749,  0.710,  0.718, 0.718,  0.718, 0.000, 0.75, 0.75, 0.75, &      ! He
       0.75,   0.75,   0.75,  0.75,   0.75,  0.75,  0.00, 0.75, 0.75, &      ! H2O
       0.75,   0.75,   0.75,  0.75,   0.75,  0.75,  0.75, 0.00, 0.75, &      ! H2
       0.75,   0.75,   0.75,  0.75,   0.75,  0.75,  0.75, 0.75, 0.00/), &    ! CO2
       (/nSpecies,nSpecies/) )

contains

  subroutine init_planet

    use ModTime

    integer :: itime(7)

    Mass(iH_)    = 1.0 * AMU
    Mass(iHe_)   = 4.0 * AMU
    Mass(iN_4S_) = 14.0 * AMU
    Mass(iO_3P_) = 16.0 * AMU
    Mass(iN_2D_) = Mass(iN_4S_)
    Mass(iN_2P_) = Mass(iN_4S_)
    Mass(iN2_)   = 2*Mass(iN_4S_)
    Mass(iO2_)   = 2*Mass(iO_3P_)
    Mass(iNO_)   = Mass(iN_4S_) + Mass(iO_3P_)
    Mass(iCO2_)  = 12.0*AMU + 2*Mass(iO_3P_)
    Mass(iH2O_)  = 2*Mass(iH_) + Mass(iO_3P_)
    Mass(iOH_)   = Mass(iO_3P_) + Mass(iH_)
    Mass(iCO_)   = Mass(iCO2_) - Mass(iO_3P_)
    Mass(iH2_)   = 2*Mass(iH_)

    cSpecies(iH_)    = "H"
    cSpecies(iHe_)   = "He"
    cSpecies(iO_3P_) = "O(!U3!NP)"
    cSpecies(iO2_)   = "O!D2!N"
    cSpecies(iN2_)   = "N!D2!N"
    cSpecies(iN_4S_) = "N(!U4!NS)"
    cSpecies(iN_2D_) = "N(!U2!ND)"
    cSpecies(iN_2P_) = "N(!U2!NP)"
    cSpecies(iNO_)   = "NO"
    cSpecies(iO_1D_) = "O(!U1!ND)"
    cSpecies(iCO2_)  = "CO!D2!N"
    cSpecies(iH2O_)  = "H2O"
    cSpecies(iOH_)   = "OH"
    cSpecies(iCO_)   = "CO"
    cSpecies(iH2_)   = "H2"
!    cSpecies(iAr_)   = "Ar"

    cIons(iO_4SP_) = "O_4SP_!U+!N"
    cIons(iO2P_)   = "O!D2!U+!N"
    cIons(iN2P_)   = "N!D2!U+!N"
    cIons(iNP_)    = "N!U+!N"
    cIons(iNOP_)   = "NO!U+!N"
    cIons(iO_2DP_) = "O(!U2!ND)!U+!N"
    cIons(iO_2PP_) = "O(!U2!NP)!U+!N"
    cIons(iHP_)    = "H!U+!N"
    cIons(iHeP_)   = "He!U+!N"
    cIons(ie_)     = "e-"

    Vibration(iO_3P_)    = 5.0
    Vibration(iO2_)   = 7.0
    Vibration(iN2_)   = 7.0
    if (nSpecies > 3) Vibration(iN_4S_) = 5.0
    if (nSpecies > 4) Vibration(iNO_)   = 7.0
    if (nSpecies > 5) Vibration(iHe_)   = 5.0

    MassI(iO_4SP_) = Mass(iO_3P_)
    MassI(iO_2DP_) = Mass(iO_3P_)
    MassI(iO_2PP_) = Mass(iO_3P_)
    MassI(iO2P_) = Mass(iO2_)
    MassI(iNP_) = Mass(iN_2D_)
    MassI(iN2P_) = Mass(iN2_)
    MassI(iHP_) = Mass(iH_)
    MassI(iHeP_) = Mass(iHe_)
    MassI(iNOP_) = Mass(iN_4S_) + Mass(iO_3P_)
    MassI(ie_) = Mass_Electron

    VertTau = 1.0e9

    itime = 0
    itime(1) = iVernalYear
    itime(2) = iVernalMonth
    itime(3) = iVernalDay
    itime(4) = iVernalHour
    itime(5) = iVernalMinute
    itime(6) = iVernalSecond
    call time_int_to_real(itime, VernalTime)



  end subroutine init_planet

!! Placeholder subroutines (for Titan specific Phyisics)

  subroutine init_radcooling
  return
  end subroutine init_radcooling

  subroutine init_magheat
  return
  end subroutine init_magheat

  subroutine init_isochem
  return
  end subroutine init_isochem

  subroutine init_aerosol
  return
  end subroutine init_aerosol

  subroutine init_topography
    return
  end subroutine init_topography

end module ModPlanet
