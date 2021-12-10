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





  interface
     subroutine fcn(n, x, y, yprime, otherdata)
       use fieldprec, only : kp,transfert         
       implicit none          
       integer :: n
       real(kp) :: x
       real(kp), dimension(n) :: y, yprime
       type(transfert), optional, intent(inout) :: otherdata
     end subroutine fcn
  end interface

