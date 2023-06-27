!\
! ------------------------------------------------------------
! Calculate O+ flux at upper boundary
! ------------------------------------------------------------
!/

! Implementation of TIEGCM O+ flux model

subroutine oplus_flux(lambdam,psi,FOplus)

  use ModInputs, only: PhiDO, PhiNO
  use ModConstants, only: pi
  implicit none

  real              :: Afactor
  real, intent(in)  :: lambdam, psi 
  real, intent(out) :: FOplus

  ! lambdam and psi are input in degrees

  if (abs(lambdam) < 7.5) then
     Afactor = 0.5*(1.0 + sin(abs(lambdam*pi/7.5) - 0.5*pi))
  else 
     Afactor = 1.0
  endif

  if (psi < 80.0) then
     FOplus = Afactor*PhiDO
  else if (psi > 100.0) then
     FOplus = Afactor*PhiNO
  else
     FOplus = Afactor*0.5*((PhiDO + PhiNO) &
          + (PhiDO - PhiNO)*cos(pi*(psi-80.0)/20.0))
  endif

end subroutine oplus_flux
