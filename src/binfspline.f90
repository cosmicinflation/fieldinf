!   This file is part of FieldInf
!
!   Copyright (C) 2005-2021 C. Ringeval (minor changes)
!
!   Copyright (C) 2000 Wolfgang Schadow (see below for the original
!   licensing)
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






! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
!   VERSION 2.2
!
!   f90 VERSION
!
!   This library contains routines for B-spline interpolation in
!   one, two, and three dimensions. Part of the routines are based
!   on the book by Carl de Boor: A practical guide to Splines (Springer,
!   New-York 1978) and have the same calling sequence and names as
!   the corresponding routines from the IMSL library. For documen-
!   tation see the additional files. NOTE: The results in the demo
!   routines may vary slightly on different architectures.
!
!   by W. Schadow 12/04/99
!   last changed by W. Schadow 07/28/2000
!
!
!   Wolfgang Schadow
!   TRIUMF
!   4004 Wesbrook Mall
!   Vancouver, B.C. V6T 2A3
!   Canada
!
!   email: schadow@triumf.ca  or  schadow@physik.uni-bonn.de
!
!   www  : http://www.triumf.ca/people/schadow
!
!
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
!   Copyright (C) 2000 Wolfgang Schadow
!
!   This library is free software; you can redistribute it and/or
!   modify it under the terms of the GNU Library General Public
!   License as published by the Free Software Foundation; either
!   version 2 of the License, or (at your option) any later version.
!
!   This library is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!   Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public
!   License along with this library; if not, write to the
!   Free Software Foundation, Inc., 59 Temple Place - Suite 330,
!   Boston, MA  02111-1307, USA.
!
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!infmod
!remove unused routines
!module numeric

!  integer, parameter :: sgl = kind(1.0)
!  integer, parameter :: dbl = kind(1.0d0)

!end module numeric

module numeric
  use fieldprec
  integer, parameter :: sgl = kind(1.0)
  integer, parameter :: dbl = kp
end module numeric
!end infmod

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


module bspline

!
!  ------------------------------------------------------------------
!
!
!   The following routines are included:
!
!            dbsnak
!
!            dbsint
!            dbsval
!            dbsder
!            dbs1gd
!
!            dbs2in
!            dbs2dr
!            dbs2vl
!            dbs2gd
!
!            dbs3in
!            dbs3vl
!            dbs3dr
!            dbs3gd
!
!  ------------------------------------------------------------------
!

  private

  public dbsnak
  public dbsint, dbsval
  public dbsder


contains


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine dbsnak(nx,xvec,kxord,xknot)

!
!  Compute the `not-a-knot' spline knot sequence.
!  (see de Boor p. 167)
!
!   nx     - number of data points.  (input)
!   xvec   - array of length ndata containing the location of the
!            data points.  (input)
!   kxord  - order of the spline.  (input)
!   xknot  - array of length ndata+korder containing the knot
!            sequence.  (output)
!

    use numeric

    implicit none

    integer, intent(in) :: nx, kxord

    real(kind=dbl), dimension(nx), intent(in)        :: xvec
    real(kind=dbl), dimension(nx+kxord), intent(out) :: xknot

    real(kind=dbl) :: eps
    integer        :: ix

!infmod
!    logical        :: first = .true.
!    save first,eps


!    if (first) then
!       first=.false.
       eps = epsilon(1.0_dbl)
!       write(*,*) "subroutine dbsnak: "
!       write(*,*) "eps = ",eps
!    endif
!end infmod

    if((kxord .lt. 0) .or. (kxord .gt. nx)) then
       write(*,*) "subroutine dbsnak: error"
       write(*,*) "0 <= kxord <= nx is required."
       write(*,*) "kxord = ", kxord, " and nx = ", nx,  " is given."
       stop
    endif

    do ix = 1, kxord
       xknot(ix) = xvec(1)
    end do

    if(mod(kxord,2) .eq. 0) then
       do ix = kxord+1, nx
          xknot(ix) = xvec(ix-kxord/2)
       end do
    else
       do ix = kxord+1, nx
          xknot(ix) = 0.5_dbl * (xvec(ix-kxord/2) + xvec(ix-kxord/2-1))
       end do
    endif

    do ix = nx+1, nx+kxord
       xknot(ix) = xvec(nx) * (1.0_dbl + eps)
    end do

  end subroutine dbsnak


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine dbsint(nx,xvec,xdata,kx,xknot,bcoef)

!
!  Computes the spline interpolant, returning the B-spline coefficients.
!  (see de Boor p. 204)
!
!   nx     - number of data points.  (input)
!   xvec   - array of length nx containing the data point
!            abscissas.  (input)
!   xdata  - array of length ndata containing the data point
!            ordinates.  (input)
!   kx     - order of the spline.  (input)
!            korder must be less than or equal to ndata.
!   xknot  - array of length nx+kx containing the knot
!            sequence.  (input)
!            xknot must be nondecreasing.
!   bscoef - array of length ndata containing the B-spline
!            coefficients.  (output)
!

    use numeric

    implicit none

    integer, intent(in)                          :: nx, kx
    real(kind=dbl), dimension(nx), intent(in)    :: xdata, xvec
    real(kind=dbl), dimension(nx+kx), intent(in) :: xknot
    real(kind=dbl), dimension(nx), intent(out)   :: bcoef

    integer                                :: nxp1, kxm1, kpkm2, leftx, lenq
    integer                                :: ix, ik,ilp1mx, jj, iflag
    real(kind=dbl)                         :: xveci
    real(kind=dbl), dimension((2*kx-1)*nx) :: work


    nxp1  = nx + 1
    kxm1  = kx - 1
    kpkm2 = 2 * kxm1
    leftx = kx
    lenq  = nx * (kx + kxm1)

    do ix = 1, lenq
       work(ix) = 0.0_dbl
    end do

    do  ix = 1, nx
       xveci  = xvec(ix)
       ilp1mx = min0(ix+kx,nxp1)
       leftx   = max0(leftx,ix)
       if (xveci .lt. xknot(leftx)) goto 998
30     if (xveci .lt. xknot(leftx+1)) go to 40
       leftx = leftx + 1
       if (leftx .lt. ilp1mx) go to 30
       leftx = leftx - 1
       if (xveci .gt. xknot(leftx+1)) goto 998
40     call bsplvb (xknot,nx+kx,kx,1,xveci,leftx,bcoef)
       jj = ix - leftx + 1 + (leftx - kx) * (kx + kxm1)
       do ik = 1, kx
          jj       = jj + kpkm2
          work(jj) = bcoef(ik)
       end do
    end do

    call banfac(work,kx+kxm1,nx,kxm1,kxm1,iflag)

    if (iflag .ne. 1) then
       write(*,*) "subroutine dbsint: error"
       write(*,*) "no solution of linear equation system !!!"
       stop
    end if

    do ix = 1, nx
       bcoef(ix) = xdata(ix)
    end do

    call banslv(work,kx+kxm1,nx,kxm1,kxm1,bcoef)

    return

998 write(*,*) "subroutine dbsint:"
    write(*,*) "xknot(ix) <= xknot(ix+1) required."
    write(*,*) ix,xknot(ix),xknot(ix+1)

    stop

  end subroutine dbsint


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  function dbsval(x,kx,xknot,nx,bcoef)

!
!  Evaluates a spline, given its B-spline representation.
!
!   x      - point at which the spline is to be evaluated.  (input)
!   kx     - order of the spline.  (input)
!   xknot  - array of length nx+kx containing the knot
!            sequence.  (input)
!            xknot must be nondecreasing.
!   nx     - number of B-spline coefficients.  (input)
!   bcoef  - array of length nx containing the B-spline
!            coefficients.  (input)
!   dbsval - value of the spline at x.  (output)
!

    use numeric

    implicit none

    integer, intent(in)                          :: nx, kx
    real(kind=dbl)                               :: dbsval
    real(kind=dbl)                               :: x
    real(kind=dbl), dimension(nx+kx), intent(in) :: xknot
    real(kind=dbl), dimension(nx), intent(in)    :: bcoef

    integer                       :: il, ik, ix, leftx
    real(kind=dbl)                :: save1, save2
    real(kind=dbl), dimension(kx) :: work, dl, dr

!
!     check if xknot(i) <= xknot(i+1) and calculation of i so that
!     xknot(i) <= x < xknot(i+1)
!

    leftx = 0

    do ix = 1,nx+kx-1
       if (xknot(ix) .gt. xknot(ix+1)) then
          write(*,*) "subroutine dbsval:"
          write(*,*) "xknot(ix) <= xknot(ix+1) required."
          write(*,*) ix,xknot(ix),xknot(ix+1)
          stop
       endif
       if((xknot(ix) .le. x) .and. (x .lt. xknot(ix+1))) leftx = ix
    end do

    if(leftx .eq. 0) then
       write(*,*) "subroutine dbsval:"
       write(*,*) "ix with xknot(ix) <= x < xknot(ix+1) required."
       write(*,*) "x = ", x
       stop
    endif

    do ik = 1, kx-1
       work(ik) = bcoef(leftx+ik-kx)
       dl(ik)   = x - xknot(leftx+ik-kx)
       dr(ik)   = xknot(leftx+ik) - x
    end do

    work(kx)  = bcoef(leftx)
    dl(kx)    = x - xknot(leftx)

    do ik = 1, kx-1
       save2 = work(ik)
       do il = ik+1, kx
          save1 = work(il)
          work(il) = (dl(il) * work(il) + dr(il-ik) * save2)                  &
               &           / (dl(il) + dr(il - ik))
          save2 = save1
       end do
    end do

    dbsval = work(kx)

  end function dbsval



  function dbsder(iderx,x,kx,xknot,nx,bcoef)

!
!  Evaluates the derivative of a spline, given its B-spline representation.
!
!
!   iderx  - order of the derivative to be evaluated.  (input)
!            in particular, iderx = 0 returns the value of the
!            spline.
!   x      - point at which the spline is to be evaluated.  (input)
!   kx     - order of the spline.  (input)
!   xknot  - array of length nx+kx containing the knot
!            sequence.  (input)
!            xknot must be nondecreasing.
!   nx     - number of B-spline coefficients.  (input)
!   bcoef  - array of length nx containing the B-spline
!            coefficients.  (input)
!   dbsder - value of the iderx-th derivative of the spline at x.
!            (output)
!

    use numeric

    implicit none

    integer, intent(in)                          :: iderx, kx, nx
    real(kind=dbl)                               :: dbsder
    real(kind=dbl), intent(in)                   :: x
    real(kind=dbl), dimension(nx+kx), intent(in) :: xknot
    real(kind=dbl), dimension(nx), intent(in)    :: bcoef

    integer                       :: ix, ik, il, leftx
    real(kind=dbl)                :: save, save1, save2, y, sum, dik
    real(kind=dbl), dimension(kx) :: work, dl, dr,bsp

!
!     check if xknot(i) <= xknot(i+1) and calculation of i so that
!     xknot(i) <= x < xknot(i+1)
!

    leftx = 0
    do ix = 1,nx+kx-1
       if (xknot(ix) .gt. xknot(ix+1)) then
          write(6,*) "subroutine dbsder:"
          write(6,*) "xknot(ix) <= xknot(ix+1) required."
          stop
       endif
       if ((xknot(ix) .le. x) .and. (x .lt. xknot(ix+1))) leftx = ix
    end do

    if (leftx .eq. 0) then
       write(6,*) "subroutine dbsder:"
       write(6,*) "ix with xknot(ix) <= x < xknot(ix+1) required."
       write(6,*) "xknot(1)     = ", xknot(1)
       write(6,*) "xknot(nx+kx) = ", xknot(nx+kx)
       write(6,*) "         x   = ", x
       stop
    endif

    if (iderx .eq. 0) then

       do ik = 1,kx-1
          work(ik) = bcoef(leftx+ik-kx)
          dl(ik)   = x - xknot(leftx+ik-kx)
          dr(ik)   = xknot(leftx+ik) - x
       end do

       work(kx)  = bcoef(leftx)
       dl(kx)    = x - xknot(leftx)

       do ik = 1,kx-1
          save2 = work(ik)
          do il = ik+1,kx
             save1 = work(il)
             work(il) = (dl(il) * work(il) + dr(il-ik) * save2)               &
                  &              / (dl(il) + dr(il - ik))
             save2 = save1
          end do
       end do

       dbsder = work(kx)

    elseif ((iderx .ge. 1) .and. (iderx .lt. kx)) then

       bsp(1) = 1.0_dbl
       do ik = 1,kx-iderx-1
          dr(ik) = xknot(leftx+ik) - x
          dl(ik) = x - xknot(leftx+1-ik)
          save   = bsp(1)
          bsp(1) = 0.0_dbl
          do il = 1, ik
             y         = save / (dr(il) + dl(ik+1-il))
             bsp(il)   = bsp(il) + dr(il) * y
             save      = bsp(il+1)
             bsp(il+1) = dl(ik+1-il) * y
          end do
       end do

       do ik = 1, kx
          work(ik) = bcoef(leftx+ik-kx)
          dr(ik)   = xknot(leftx+ik) - x
          dl(ik)   = x - xknot(leftx+ik-kx)
       end do

       do ik = 1, iderx
          dik   = dble(kx - ik)
          save2 = work(ik)
          do il = ik+1, kx
             save1    = work(il)
             work(il) = dik * (work(il) - save2) /(dl(il) + dr(il-ik))
             save2    = save1
          end do
       end do

       sum = 0.0_dbl

       do ix = 1, kx-iderx
          sum = sum + bsp(ix) * work(iderx+ix)
       end do

       dbsder = sum

    else
       dbsder = 0.0_dbl
    endif

  end function dbsder





  subroutine bsplvb(t,n,jhigh,index,x,left,biatx)

    use numeric

    implicit none

    integer, intent(in) :: n, jhigh, index, left

    real(kind=dbl), intent(in)                    :: x
    real(kind=dbl), dimension(n), intent(in)      :: t
    real(kind=dbl), dimension(jhigh), intent(out) :: biatx

    integer                          :: j = 1
    integer                          :: i, jp1
    real(kind=dbl)                   :: saved, term
    real(kind=dbl), dimension(jhigh) :: dl, dr


    if (index .eq. 1) then
       j = 1
       biatx(1) = 1.0_dbl
       if (j .ge. jhigh) return
    end if

20  jp1 = j + 1

    dr(j) = t(left+j) - x
    dl(j) = x - t(left+1-j)
    saved = 0._dbl

    do i = 1, j
       term     = biatx(i) / (dr(i) + dl(jp1-i))
       biatx(i) = saved + dr(i) * term
       saved    = dl(jp1-i) * term
    end do

    biatx(jp1) = saved
    j          = jp1

    if (j .lt. jhigh) go to 20

  end subroutine bsplvb


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine banfac(w,nroww,nrow,nbandl,nbandu,iflag)

    use numeric

    implicit none

    integer, intent(in)                                  :: nroww,nrow
    integer, intent(in)                                  :: nbandl,nbandu
    integer, intent(out)                                 :: iflag
    real(kind=dbl), dimension(nroww,nrow), intent(inout) :: w

    real(kind=dbl) :: pivot, factor
    integer        :: middle, nrowm1, jmax, kmax, ipk, midmk, i, j, k


    iflag  = 1
    middle = nbandu + 1
    nrowm1 = nrow - 1

    if (nrowm1 .lt. 0) goto 999
    if (nrowm1 .eq. 0) goto 900
    if (nrowm1 .gt. 0) goto 10

10  if (nbandl .gt. 0) go to 30

    do i = 1, nrowm1
       if (w(middle,i) .eq. 0._dbl) go to 999
    end do

    go to 900

30  if (nbandu .gt. 0) go to 60

    do i = 1, nrowm1
       pivot = w(middle,i)
       if(pivot .eq. 0._dbl) go to 999
       jmax = min0(nbandl, nrow - i)
       do j = 1, jmax
          w(middle+j,i) = w(middle+j,i) / pivot
       end do
    end do

    return

60  do i = 1, nrowm1
       pivot = w(middle,i)
       if (pivot .eq. 0._dbl) go to 999
       jmax = min0(nbandl,nrow - i)
       do j = 1,jmax
          w(middle+j,i) = w(middle+j,i) / pivot
       end do

       kmax = min0(nbandu,nrow - i)

       do k = 1, kmax
          ipk    = i + k
          midmk  = middle - k
          factor = w(midmk,ipk)
          do j = 1, jmax
             w(midmk+j,ipk) = w(midmk+j,ipk) - w(middle+j,i)                  &
                  &              * factor
          end do
       end do
    end do

900 if (w(middle,nrow) .ne. 0._dbl) return
999 iflag = 2

  end subroutine banfac


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  subroutine banslv(w,nroww,nrow,nbandl,nbandu,b)

    use numeric

    implicit none

    integer, intent(in)                               :: nroww,nrow
    integer, intent(in)                               :: nbandl,nbandu
    real(kind=dbl), dimension(nroww,nrow), intent(in) :: w
    real(kind=dbl), dimension(nrow), intent(inout)    :: b

    integer :: middle, nrowm1, jmax, i, j

    middle = nbandu + 1
    if (nrow .eq. 1) goto 99
    nrowm1 = nrow - 1
    if (nbandl .eq. 0) goto 30

    do i = 1, nrowm1
       jmax = min0(nbandl, nrow - i)
       do j = 1, jmax
          b(i+j) = b(i+j) - b(i) * w(middle+j,i)
       end do
    end do

30  if (nbandu .gt. 0)  goto 50

    do i = 1, nrow
       b(i) = b(i) / w(1,i)
    end do

    return

50  do i = nrow, 2, -1
       b(i) = b(i)/w(middle,i)
       jmax = min0(nbandu,i-1)
       do j = 1, jmax
          b(i-j) = b(i-j) - b(i) * w(middle-j,i)
       end do
    end do

99  b(1) = b(1) / w(middle,1)

  end subroutine banslv



end module bspline
