!Some cosmological parameters used for inflation. This module sets
!their fiducial values because observable predictions have only a weak
!dependence in this (may change with increased data accuracy).
module cosmopar
  use infprec, only : kp
  implicit none
  
  public

  real(kp), parameter :: HubbleSquareRootOf3OmegaRad = 7.4585d-63
  real(kp), parameter :: HubbleSquareRootOf2OmegaRad = sqrt(2._kp/3._kp)*HubbleSquareRootOf3OmegaRad

  real(kp), parameter :: lnMpcToKappa = 130.282_kp

!1MeV
!!  real(kp), parameter :: lnRhoNuc = -196.97
!10MeV
  real(kp), parameter :: lnRhoNuc = -187.747
!100MeV
!!  real(kp), parameter :: lnRhoNuc = -178.55

!only used for reheating using slow-roll (libslowroll)
!COBE quadrupole moment
  real(kp), parameter :: QrmsOverT = 6e-6
!Best scalar amp for slow-roll (update with new constraints)
  real(kp), parameter :: powerAmpScalar = 2.165e-9

end module cosmopar
