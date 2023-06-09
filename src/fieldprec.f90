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





module fieldprec
  implicit none

  public
  
!quad precision
!  integer, parameter :: kp = kind(1.0_16)

!double precision
  integer, parameter :: kp = kind(1.0_8)
  
!home made precision: p number of digit
! integer, parameter :: kp = selected_real_kind(p=32)

!default integration accuracy
  real(kp), parameter :: tolkp = 1.d-11

!default length for short strings
  integer, parameter :: lenshort = 6

!pi
  real(kp), parameter :: pi = 3.141592653589793238_kp

!workaround for passing argument to old f77 functions. Only pointer
!can be deferred shape in derived data type.
!Allows to stop integration from conditions coming from called
!functions (find the end of inflation)
  type transfert
     logical :: yesno1,yesno2, yesno3
     logical :: mat, hub, ismax
     integer :: vsign
     integer :: int1, int2, int3
     real(kp) :: real1, real2, real3, real4, real5
     real(kp), dimension(:), pointer :: ptrvector1 => null()
     real(kp), dimension(:), pointer :: ptrvector2 => null()
!reserved
     logical :: check,update
     real(kp) :: xend
  end type transfert


end module fieldprec
