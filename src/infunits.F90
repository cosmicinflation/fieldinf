!Some cosmological parameter values used for inflation. This module
!sets their standard values because observable predictions have only a
!weak dependence in this (may change with increased data accuracy). If
!aspic is found, use the one defined in libaspic

module infunits
#if !defined (NOASPIC)
  use cosmopar
#else
  use infprec, only : kp
  implicit none
  
  public

  real(kp), parameter :: HubbleSquareRootOf3OmegaRad = 7.4585d-63
  real(kp), parameter :: HubbleSquareRootOf2OmegaRad = sqrt(2._kp/3._kp)*HubbleSquareRootOf3OmegaRad

  real(kp), parameter :: lnMpcToKappa = 130.282_kp

  real(kp), parameter :: lnMpinGeV=42.334_kp

!1MeV
!!  real(kp), parameter :: lnRhoNuc = -196.97
!10MeV
  real(kp), parameter :: lnRhoNuc = -187.747
!100MeV
!!  real(kp), parameter :: lnRhoNuc = -178.55

#endif

end module infunits