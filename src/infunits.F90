!   This file is part of FieldInf
!
!   Copyright (C) 2005-2021 C. Ringeval
!   
!   FieldInf is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   FieldInf is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with FieldInf.  If not, see <https://www.gnu.org/licenses/>.






!Some cosmological parameter values used for inflation. This module
!sets their standard values because observable predictions have only a
!weak dependence in this (may change with increased data accuracy). If
!aspic is found, use the one defined in libaspic

module infunits
#if !defined (NOASPIC)
  use cosmopar
#else
  use fieldprec, only : kp
  implicit none
  
  public
 
  real(kp), parameter :: HubbleSquareRootOf3OmegaRad = 7.5437d-63
  real(kp), parameter :: HubbleSquareRootOf2OmegaRad = sqrt(2._kp/3._kp)*HubbleSquareRootOf3OmegaRad

  real(kp), parameter :: lnMpcToKappa = 130.282_kp

  real(kp), parameter :: lnMpinGeV=42.334_kp

!1MeV
!!  real(kp), parameter :: lnRhoNuc = -196.98
!10MeV
  real(kp), parameter :: lnRhoNuc = -187.77
!100MeV
!!  real(kp), parameter :: lnRhoNuc = -178.56

#endif

end module infunits
