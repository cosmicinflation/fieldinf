!   This file is part of FieldInf
!
!   Copyright (C) 2005-2021 C. Ringeval
!
!   Copyright (C) 1985 Kenneth R. Jackson for dverk (public domain)
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





module infsolvers
  use fieldprec, only : kp,transfert  
  implicit none

  private

  public easydverk, tunedverk
  public zbrent



contains
 
 

  subroutine easydverk(n,fcn,x,y,xend,tol,extradata)
    implicit none
     
    integer :: n,ind
    real(kp) :: x,y(n),xend,tol,c(24),w(n,9)    
    type(transfert), optional,intent(inout) :: extradata    

!    external :: fcn

    include 'infsolvers.h'

   
!    ind=1    
    ind=2
    c = 0._kp
    c(3) = epsilon(1._kp)    
 
    call dverk(n,fcn,x,y,xend,tol,ind,c,n,w,extradata)  
    
    if (ind.ne.3) then
       write(*,*) 'easydverk: stop ind = ',ind
       write(*,*) 'try to reduce tolerance'
       stop
    endif    
  end subroutine easydverk
 
 

  subroutine tunedverk(n,fcn,x,y,xend,tol,extradata)
    implicit none
     
    integer :: n,ind
    real(kp) :: x,y(n),xend,tol,c(24),w(n,9)    
    type(transfert), optional,intent(inout) :: extradata    

!    external :: fcn

    include 'infsolvers.h'

   

    ind=2
    c = 0._kp
    c(3) = epsilon(1._kp)
    c(4) = epsilon(1._kp)
 
    call dverk(n,fcn,x,y,xend,tol,ind,c,n,w,extradata)  
    
    if (ind.ne.3) then
       write(*,*) 'unsanedverk: stop ind = ',ind
       write(*,*) 'desesperate accuracy unreachable...'
       write(*,*) 'try tuning c(4) and c(6) in infsolvers:unsanedverk!'
       stop
    endif    
  end subroutine tunedverk
 
  




  subroutine diydverk(n,fcn,x,y,xend,tol,ind,c,extradata)
    implicit none
     
    integer :: n,ind
    real(kp) :: x,y(n),xend,tol,c(24),w(n,9)    
    type(transfert), optional, intent(inout) :: extradata    

!    external :: fcn

    include 'infsolvers.h'

    ind=2
    c = 0._kp

    c(3) = epsilon(1._kp)    
    c(4) = 0.001_kp

    call dverk(n,fcn,x,y,xend,tol,ind,c,n,w,extradata)  
    
    if ((ind.ne.3).and.(ind.ne.7)) then
       write(*,*) 'diydverk: stop ind = ',ind
       stop
    endif    
  end subroutine diydverk




!  function easyfixpnf(xstart,f,fjac,extradata)    
!    use hompack, only : rhodum, rhojacdum,fixpnf    
!    implicit none
   
!    real(kp) :: xstart, easyfixpnf
!    type(transfert), optional :: extradata


!    real(kp), parameter :: tolZero = 1e-5

!    integer, parameter :: n = 1
!    integer :: iflag, trace, nfe
!    real(kp) :: arcre,arcae
!    real(kp) :: ansre,ansae
!    real(kp) :: arclen

!    real(kp), dimension(n) :: a
!    real(kp), dimension(n+1) :: y,yp,yold,ypold

!    real(kp), dimension(1:N,1:N+2) :: qr
!    real(kp), dimension(n) :: alpha
!    real(kp), dimension(n+1) :: tz,w,wp,z0,z1

!    real(kp), dimension(8) :: sspar
   
!    integer, dimension(n+1) :: pivot

!    include 'hominfpack.h'

!    iflag = -1
!    trace = -1
!    sspar = -1

!    arcre = -1.
!    arcae = -1.
!    ansre = tolZero
!    ansae = tolZero

    
!    y(n+1) = xstart

!    call fixpnf(N,Y,IFLAG,ARCRE,ARCAE,ANSRE,ANSAE,TRACE,A,NFE, &
!         ARCLEN,YP,YOLD,YPOLD,QR,ALPHA,TZ,PIVOT,W,WP,Z0,Z1,SSPAR, &
!         f,fjac,rhodum,rhojacdum,extradata)

!    if (iflag.ne.1) then
!       write(*,*)'easyzero: iflag = ',iflag
!       stop
!    endif

!    easyfixpnf = y(n+1)

!  end function easyfixpnf







  subroutine dverk(n,fcn,x,y,xend,tol,ind,c,nw,w,extradata)
    implicit none
    integer :: n, ind, nw, k
    real(kp) :: x, y(n), xend, tol, c(*), w(nw,9), temp
    type(transfert), optional, intent(inout) :: extradata
    

!
!***********************************************************************
!                                                                      *
! note added 11/14/85.                                                 *
!                                                                      *
! if you discover any errors in this subroutine, please contact        *
!                                                                      *
!        kenneth r. jackson                                            *
!        department of computer science                                *
!        university of toronto                                         *
!        toronto, ontario,                                             *
!        canada   m5s 1a4                                              *
!                                                                      *
!        phone: 416-978-7075                                           *
!                                                                      *
!        electroni! mail:                                              *
!        uucp:   {cornell,decvax,ihnp4,linus,uw-beaver}!utcsri!krj     *
!        csnet:  krj@toronto                                           *
!        arpa:   krj.toronto@csnet-relay                               *
!        bitnet: krj%toronto@csnet-relay.arpa                          *
!                                                                      *
! dverk is written in fortran 66.                                      *
!                                                                      *
! the constants dwarf and rreb -- c(10) and c(11), respectively -- are *
! set for a  vax  in  double  precision.  they  should  be  reset,  as *
! described below, if this program is run on another machine.          *
!                                                                      *
! the c array is declared in this subroutine to have one element only, *
! although  more  elements  are  referenced  in this subroutine.  this *
! causes some compilers to issue warning messages.  there is,  though, *
! no  error  provided  c is declared sufficiently large in the calling *
! program, as described below.                                         *
!                                                                      *
! the following external statement  for  fcn  was  added  to  avoid  a *
! warning  message  from  the  unix  f77 compiler.  the original dverk *
! comments and code follow it.                                         *
!                                                                      *
!***********************************************************************
!

!EXTRADATA

!might be dangerous (xlf90) with optional argument, rather use explicit
!interface
!      external fcn
    include 'infsolvers.h'

    if (present(extradata)) extradata%update = .false.
!
!***********************************************************************
!                                                                      *
!     purpose - this is a runge-kutta  subroutine  based  on  verner's *
! fifth and sixth order pair of formulas for finding approximations to *
! the solution of  a  system  of  first  order  ordinary  differential *
! equations  with  initial  conditions. it attempts to keep the global *
! error proportional to  a  tolerance  specified  by  the  user.  (the *
! proportionality  depends  on the kind of error control that is used, *
! as well as the differential equation and the range of integration.)  *
!                                                                      *
!     various options are available to the user,  including  different *
! kinds  of  error control, restrictions on step sizes, and interrupts *
! which permit the user to examine the state of the  calculation  (and *
! perhaps make modifications) during intermediate stages.              *
!                                                                      *
!     the program is efficient for non-stiff systems.  however, a good *
! variable-order-adams  method  will probably be more efficient if the *
! function evaluations are very costly.  such a method would  also  be *
! more suitable if one wanted to obtain a large number of intermediate *
! solution values by interpolation, as might be the case  for  example *
! with graphical output.                                               *
!                                                                      *
!                                    hull-enright-jackson   1/10/76    *
!                                                                      *
!***********************************************************************
!                                                                      *
!     use - the user must specify each of the following                *
!                                                                      *
!     n  number of equations                                           *
!                                                                      *
!   fcn  name of subroutine for evaluating functions - the  subroutine *
!           itself must also be provided by the user - it should be of *
!           the following form                                         *
!              subroutine fcn(n, x, y, yprime)                         *
!              integer n                                               *
!              double precision x, y(n), yprime(n)                     *
!                      *** etc ***                                     *
!           and it should evaluate yprime, given n, x and y            *
!                                                                      *
!     x  independent variable - initial value supplied by user         *
!                                                                      *
!     y  dependent variable - initial values of components y(1), y(2), *
!           ..., y(n) supplied by user                                 *
!                                                                      *
!  xend  value of x to which integration is to be carried out - it may *
!           be less than the initial value of x                        *
!                                                                      *
!   tol  tolerance - the subroutine attempts to control a norm of  the *
!           local  error  in  such  a  way  that  the  global error is *
!           proportional to tol. in some problems there will be enough *
!           damping  of  errors, as well as some cancellation, so that *
!           the global error will be less than tol. alternatively, the *
!           control   can   be  viewed  as  attempting  to  provide  a *
!           calculated value of y at xend which is the exact  solution *
!           to  the  problem y' = f(x,y) + e(x) where the norm of e(x) *
!           is proportional to tol.  (the norm  is  a  max  norm  with *
!           weights  that  depend on the error control strategy chosen *
!           by the user.  the default weight for the k-th component is *
!           1/max(1,abs(y(k))),  which therefore provides a mixture of *
!           absolute and relative error control.)                      *
!                                                                      *
!   ind  indicator - on initial entry ind must be set equal to  either *
!           1  or  2. if the user does not wish to use any options, he *
!           should set ind to 1 - all that remains for the user to  do *
!           then  is  to  declare c and w, and to specify nw. the user *
!           may also  select  various  options  on  initial  entry  by *
!           setting ind = 2 and initializing the first 9 components of *
!           c as described in the next section.  he may also  re-enter *
!           the  subroutine  with ind = 3 as mentioned again below. in *
!           any event, the subroutine returns with ind equal to        *
!              3 after a normal return                                 *
!              4, 5, or 6 after an interrupt (see options c(8), c(9))  *
!              -1, -2, or -3 after an error condition (see below)      *
!                                                                      *
!     c  communications vector - the dimension must be greater than or *
!           equal to 24, unless option c(1) = 4 or 5 is used, in which *
!           case the dimension must be greater than or equal to n+30   *
!                                                                      *
!    nw  first dimension of workspace w -  must  be  greater  than  or *
!           equal to n                                                 *
!                                                                      *
!     w  workspace matrix - first dimension must be nw and second must *
!           be greater than or equal to 9                              *
!                                                                      *
!     the subroutine  will  normally  return  with  ind  =  3,  having *
! replaced the initial values of x and y with, respectively, the value *
! of xend and an approximation to y at xend.  the  subroutine  can  be *
! called  repeatedly  with new values of xend without having to change *
! any other argument.  however, changes in tol, or any of the  options *
! described below, may also be made on such a re-entry if desired.     *
!                                                                      *
!     three error returns are also possible, in which  case  x  and  y *
! will be the most recently accepted values -                          *
!     with ind = -3 the subroutine was unable  to  satisfy  the  error *
!        requirement  with a particular step-size that is less than or *
!        equal to hmin, which may mean that tol is too small           *
!     with ind = -2 the value of hmin  is  greater  than  hmax,  which *
!        probably  means  that the requested tol (which is used in the *
!        calculation of hmin) is too small                             *
!     with ind = -1 the allowed maximum number of fcn evaluations  has *
!        been  exceeded,  but  this  can only occur if option c(7), as *
!        described in the next section, has been used                  *
!                                                                      *
!     there are several circumstances that will cause the calculations *
! to  be  terminated,  along with output of information that will help *
! the user determine the cause of  the  trouble.  these  circumstances *
! involve  entry with illegal or inconsistent values of the arguments, *
! such as attempting a normal  re-entry  without  first  changing  the *
! value of xend, or attempting to re-enter with ind less than zero.    *
!                                                                      *
!***********************************************************************
!                                                                      *
!     options - if the subroutine is entered with ind = 1, the first 9 *
! components of the communications vector are initialized to zero, and *
! the subroutine uses only default values  for  each  option.  if  the *
! subroutine  is  entered  with ind = 2, the user must specify each of *
! these 9 components - normally he would first set them all  to  zero, *
! and  then  make  non-zero  those  that  correspond to the particular *
! options he wishes to select. in any event, options may be changed on *
! re-entry  to  the  subroutine  -  but if the user changes any of the *
! options, or tol, in the course of a calculation he should be careful *
! about  how  such changes affect the subroutine - it may be better to *
! restart with ind = 1 or 2. (components 10 to 24 of c are used by the *
! program  -  the information is available to the user, but should not *
! normally be changed by him.)                                         *
!                                                                      *
!  c(1)  error control indicator - the norm of the local error is  the *
!           max  norm  of  the  weighted  error  estimate  vector, the *
!           weights being determined according to the value of c(1) -  *
!              if c(1)=1 the weights are 1 (absolute error control)    *
!              if c(1)=2 the weights are 1/abs(y(k))  (relative  error *
!                 control)                                             *
!              if c(1)=3 the  weights  are  1/max(abs(c(2)),abs(y(k))) *
!                 (relative  error  control,  unless abs(y(k)) is less *
!                 than the floor value, abs(c(2)) )                    *
!              if c(1)=4 the weights are 1/max(abs(c(k+30)),abs(y(k))) *
!                 (here individual floor values are used)              *
!              if c(1)=5 the weights are 1/abs(c(k+30))                *
!              for all other values of c(1), including  c(1) = 0,  the *
!                 default  values  of  the  weights  are  taken  to be *
!                 1/max(1,abs(y(k))), as mentioned earlier             *
!           (in the two cases c(1) = 4 or 5 the user must declare  the *
!           dimension of c to be at least n+30 and must initialize the *
!           components c(31), c(32), ..., c(n+30).)                    *
!                                                                      *
!  c(2)  floor value - used when the indicator c(1) has the value 3    *
!                                                                      *
!  c(3)  hmin specification - if not zero, the subroutine chooses hmin *
!           to be abs(c(3)) - otherwise it uses the default value      *
!              10*max(dwarf,rreb*max(weighted norm y/tol,abs(x))),     *
!           where dwarf is a very small positive  machine  number  and *
!           rreb is the relative roundoff error bound                  *
!                                                                      *
!  c(4)  hstart specification - if not zero, the subroutine  will  use *
!           an  initial  hmag equal to abs(c(4)), except of course for *
!           the restrictions imposed by hmin and hmax  -  otherwise it *
!           uses the default value of hmax*(tol)**(1/6)                *
!                                                                      *
!  c(5)  scale specification - this is intended to be a measure of the *
!           scale of the problem - larger values of scale tend to make *
!           the method more reliable, first  by  possibly  restricting *
!           hmax  (as  described  below) and second, by tightening the *
!           acceptance requirement - if c(5) is zero, a default  value *
!           of  1  is  used.  for  linear  homogeneous  problems  with *
!           constant coefficients, an appropriate value for scale is a *
!           norm  of  the  associated  matrix.  for other problems, an *
!           approximation to  an  average  value  of  a  norm  of  the *
!           jacobian along the trajectory may be appropriate           *
!                                                                      *
!  c(6)  hmax specification - four cases are possible                  *
!           if c(6).ne.0 and c(5).ne.0, hmax is taken to be            *
!              min(abs(c(6)),2/abs(c(5)))                              *
!           if c(6).ne.0 and c(5).eq.0, hmax is taken to be  abs(c(6)) *
!           if c(6).eq.0 and c(5).ne.0, hmax is taken to be            *
!              2/abs(c(5))                                             *
!           if c(6).eq.0 and c(5).eq.0, hmax is given a default  value *
!              of 2                                                    *
!                                                                      *
!  c(7)  maximum number of function evaluations  -  if  not  zero,  an *
!           error  return with ind = -1 will be caused when the number *
!           of function evaluations exceeds abs(c(7))                  *
!                                                                      *
!  c(8)  interrupt number  1  -  if  not  zero,  the  subroutine  will *
!           interrupt   the  calculations  after  it  has  chosen  its *
!           preliminary value of hmag, and just before choosing htrial *
!           and  xtrial  in  preparation for taking a step (htrial may *
!           differ from hmag in sign, and may  require  adjustment  if *
!           xend  is  near) - the subroutine returns with ind = 4, and *
!           will resume calculation at the point  of  interruption  if *
!           re-entered with ind = 4                                    *
!                                                                      *
!  c(9)  interrupt number  2  -  if  not  zero,  the  subroutine  will *
!           interrupt   the  calculations  immediately  after  it  has *
!           decided whether or not to accept the result  of  the  most *
!           recent  trial step, with ind = 5 if it plans to accept, or *
!           ind = 6 if it plans to reject -  y(*)  is  the  previously *
!           accepted  result, while w(*,9) is the newly computed trial *
!           value, and w(*,2) is the unweighted error estimate vector. *
!           the  subroutine  will  resume calculations at the point of *
!           interruption on re-entry with ind = 5 or 6. (the user  may *
!           change ind in this case if he wishes, for example to force *
!           acceptance of a step that would otherwise be rejected,  or *
!           vice versa. he can also restart with ind = 1 or 2.)        *
!                                                                      *
!***********************************************************************
!                                                                      *
!  summary of the components of the communications vector              *
!                                                                      *
!     prescribed at the option       determined by the program         *
!           of the user                                                *
!                                                                      *
!                                    c(10) rreb(rel roundoff err bnd)  *
!     c(1) error control indicator   c(11) dwarf (very small mach no)  *
!     c(2) floor value               c(12) weighted norm y             *
!     c(3) hmin specification        c(13) hmin                        *
!     c(4) hstart specification      c(14) hmag                        *
!     c(5) scale specification       c(15) scale                       *
!     c(6) hmax specification        c(16) hmax                        *
!     c(7) max no of fcn evals       c(17) xtrial                      *
!     c(8) interrupt no 1            c(18) htrial                      *
!     c(9) interrupt no 2            c(19) est                         *
!                                    c(20) previous xend               *
!                                    c(21) flag for xend               *
!                                    c(22) no of successful steps      *
!                                    c(23) no of successive failures   *
!                                    c(24) no of fcn evals             *
!                                                                      *
!  if c(1) = 4 or 5, c(31), c(32), ... c(n+30) are floor values        *
!                                                                      *
!***********************************************************************
!                                                                      *
!  an overview of the program                                          *
!                                                                      *
!     begin initialization, parameter checking, interrupt re-entries   *
!  ......abort if ind out of range 1 to 6                              *
!  .     cases - initial entry, normal re-entry, interrupt re-entries  *
!  .     case 1 - initial entry (ind .eq. 1 or 2)                      *
!  v........abort if n.gt.nw or tol.le.0                               *
!  .        if initial entry without options (ind .eq. 1)              *
!  .           set c(1) to c(9) equal to zero                          *
!  .        else initial entry with options (ind .eq. 2)               *
!  .           make c(1) to c(9) non-negative                          *
!  .           make floor values non-negative if they are to be used   *
!  .        end if                                                     *
!  .        initialize rreb, dwarf, prev xend, flag, counts            *
!  .     case 2 - normal re-entry (ind .eq. 3)                         *
!  .........abort if xend reached, and either x changed or xend not    *
!  .        re-initialize flag                                         *
!  .     case 3 - re-entry following an interrupt (ind .eq. 4 to 6)    *
!  v        transfer control to the appropriate re-entry point.......  *
!  .     end cases                                                  .  *
!  .  end initialization, etc.                                      .  *
!  .                                                                v  *
!  .  loop through the following 4 stages, once for each trial step .  *
!  .     stage 1 - prepare                                          .  *
!***********error return (with ind=-1) if no of fcn evals too great .  *
!  .        calc slope (adding 1 to no of fcn evals) if ind .ne. 6  .  *
!  .        calc hmin, scale, hmax                                  .  *
!***********error return (with ind=-2) if hmin .gt. hmax            .  *
!  .        calc preliminary hmag                                   .  *
!***********interrupt no 1 (with ind=4) if requested.......re-entry.v  *
!  .        calc hmag, xtrial and htrial                            .  *
!  .     end stage 1                                                .  *
!  v     stage 2 - calc ytrial (adding 7 to no of fcn evals)        .  *
!  .     stage 3 - calc the error estimate                          .  *
!  .     stage 4 - make decisions                                   .  *
!  .        set ind=5 if step acceptable, else set ind=6            .  *
!***********interrupt no 2 if requested....................re-entry.v  *
!  .        if step accepted (ind .eq. 5)                              *
!  .           update x, y from xtrial, ytrial                         *
!  .           add 1 to no of successful steps                         *
!  .           set no of successive failures to zero                   *
!**************return(with ind=3, xend saved, flag set) if x .eq. xend *
!  .        else step not accepted (ind .eq. 6)                        *
!  .           add 1 to no of successive failures                      *
!**************error return (with ind=-3) if hmag .le. hmin            *
!  .        end if                                                     *
!  .     end stage 4                                                   *
!  .  end loop                                                         *
!  .                                                                   *
!  begin abort action                                                  *
!     output appropriate  message  about  stopping  the  calculations, *
!        along with values of ind, n, nw, tol, hmin,  hmax,  x,  xend, *
!        previous xend,  no of  successful  steps,  no  of  successive *
!        failures, no of fcn evals, and the components of y            *
!     stop                                                             *
!  end abort action                                                    *
!                                                                      *
!***********************************************************************
!
!     ******************************************************************
!     * begin initialization, parameter checking, interrupt re-entries *
!     ******************************************************************
!
!  ......abort if ind out of range 1 to 6
         if (ind.lt.1 .or. ind.gt.6) go to 500
!
!        cases - initial entry, normal re-entry, interrupt re-entries
         go to (5, 5, 45, 1111, 2222, 2222), ind
!        case 1 - initial entry (ind .eq. 1 or 2)
!  .........abort if n.gt.nw or tol.le.0
    5       if (n.gt.nw .or. tol.le.0._kp) go to 500
            if (ind.eq. 2) go to 15
!              initial entry without options (ind .eq. 1)
!              set c(1) to c(9) equal to 0
               do 10 k = 1, 9
                  c(k) = 0._kp
   10          continue
               go to 35
   15       continue
!              initial entry with options (ind .eq. 2)
!              make c(1) to c(9) non-negative
               do 20 k = 1, 9
                  c(k) = abs(c(k))
   20          continue
!              make floor values non-negative if they are to be used
               if (c(1).ne.4._kp .and. c(1).ne.5._kp) go to 30
                  do 25 k = 1, n
                     c(k+30) = abs(c(k+30))
   25             continue
   30          continue
   35       continue
!           initialize rreb, dwarf, prev xend, flag, counts
            c(10) = 2._kp**(-56)
            c(11) = 1.d-35
!           set previous xend initially to initial value of x
            c(20) = x
            do 40 k = 21, 24
               c(k) = 0._kp
   40       continue
            go to 50
!        case 2 - normal re-entry (ind .eq. 3)
!  .........abort if xend reached, and either x changed or xend not
   45       if (c(21).ne.0._kp .and. &
                 (x.ne.c(20) .or. xend.eq.c(20))) go to 500
!           re-initialize flag
            c(21) = 0._kp
            go to 50
!        case 3 - re-entry following an interrupt (ind .eq. 4 to 6)
!           transfer control to the appropriate re-entry point..........
!           this has already been handled by the computed go to        .
!        end cases                                                     v
   50    continue
!
!     end initialization, etc.
!
!     ******************************************************************
!     * loop through the following 4 stages, once for each trial  step *
!     * until the occurrence of one of the following                   *
!     *    (a) the normal return (with ind .eq. 3) on reaching xend in *
!     *        stage 4                                                 *
!     *    (b) an error return (with ind .lt. 0) in stage 1 or stage 4 *
!     *    (c) an interrupt return (with ind  .eq.  4,  5  or  6),  if *
!     *        requested, in stage 1 or stage 4                        *
!     ******************************************************************
!
99999 continue
!
!        ***************************************************************
!        * stage 1 - prepare - do calculations of  hmin,  hmax,  etc., *
!        * and some parameter  checking,  and  end  up  with  suitable *
!        * values of hmag, xtrial and htrial in preparation for taking *
!        * an integration step.                                        *
!        ***************************************************************
!
!***********error return (with ind=-1) if no of fcn evals too great
            if (c(7).eq.0._kp .or. c(24).lt.c(7)) go to 100
               ind = -1
               return
  100       continue
!
!           calculate slope (adding 1 to no of fcn evals) if ind .ne. 6
            if (ind .eq. 6) go to 105
!EXTRADATA call fcn(n, x, y, w(1,1))
               if (present(extradata)) extradata%check=.true.               
               call fcn(n, x, y, w(1,1),extradata)             
               if (present(extradata)) then
                  if (extradata%update) xend = extradata%xend
                  extradata%check=.false.                  
               endif
!END EXTRADATA
               c(24) = c(24) + 1._kp
  105       continue
!
!           calculate hmin - use default unless value prescribed
            c(13) = c(3)
            if (c(3) .ne. 0._kp) go to 165
!              calculate default value of hmin
!              first calculate weighted norm y - c(12) - as specified
!              by the error control indicator c(1)
               temp = 0._kp
               if (c(1) .ne. 1._kp) go to 115
!                 absolute error control - weights are 1
                  do 110 k = 1, n
                     temp = max(temp, abs(y(k)))
  110             continue
                  c(12) = temp
                  go to 160
  115          if (c(1) .ne. 2._kp) go to 120
!                 relative error control - weights are 1/abs(y(k)) so
!                 weighted norm y is 1
                  c(12) = 1._kp
                  go to 160
  120          if (c(1) .ne. 3._kp) go to 130
!                 weights are 1/max(c(2),abs(y(k)))
                  do 125 k = 1, n
                     temp = max(temp, abs(y(k))/c(2))
  125             continue
                  c(12) = min(temp, 1._kp)
                  go to 160
  130          if (c(1) .ne. 4._kp) go to 140
!                 weights are 1/max(c(k+30),abs(y(k)))
                  do 135 k = 1, n
                     temp = max(temp, abs(y(k))/c(k+30))
  135             continue
                  c(12) = min(temp, 1._kp)
                  go to 160
  140          if (c(1) .ne. 5._kp) go to 150
!                 weights are 1/c(k+30)
                  do 145 k = 1, n
                     temp = max(temp, abs(y(k))/c(k+30))
  145             continue
                  c(12) = temp
                  go to 160
  150          continue
!                 default case - weights are 1/max(1,abs(y(k)))
                  do 155 k = 1, n
                     temp = max(temp, abs(y(k)))
  155             continue
                  c(12) = min(temp, 1._kp)
  160          continue
               c(13) = 10._kp*max(c(11),c(10)*max(c(12)/tol,abs(x)))
  165       continue
!
!           calculate scale - use default unless value prescribed
            c(15) = c(5)
            if (c(5) .eq. 0._kp) c(15) = 1._kp
!
!           calculate hmax - consider 4 cases
!           case 1 both hmax and scale prescribed
               if (c(6).ne.0._kp .and. c(5).ne.0._kp) &
                    c(16) = min(c(6), 2._kp/c(5))
!           case 2 - hmax prescribed, but scale not
               if (c(6).ne.0._kp .and. c(5).eq.0._kp) c(16) = c(6)
!           case 3 - hmax not prescribed, but scale is
               if (c(6).eq.0._kp .and. c(5).ne.0._kp) c(16) = 2._kp/c(5)
!           case 4 - neither hmax nor scale is provided
               if (c(6).eq.0._kp .and. c(5).eq.0._kp) c(16) = 2._kp
!
!***********error return (with ind=-2) if hmin .gt. hmax
            if (c(13) .le. c(16)) go to 170
               ind = -2
               return
  170       continue
!
!           calculate preliminary hmag - consider 3 cases
            if (ind .gt. 2) go to 175
!           case 1 - initial entry - use prescribed value of hstart, if
!              any, else default
               c(14) = c(4)
               if (c(4) .eq. 0._kp) c(14) = c(16)*tol**(1./6.)
               go to 185
  175       if (c(23) .gt. 1._kp) go to 180
!           case 2 - after a successful step, or at most  one  failure,
!              use min(2, .9*(tol/est)**(1/6))*hmag, but avoid possible
!              overflow. then avoid reduction by more than half.
               temp = 2._kp*c(14)
               if (tol .lt. (2._kp/.9_kp)**6*c(19)) &
                    temp = .9_kp*(tol/c(19))**(1./6.)*c(14)
               c(14) = max(temp, .5_kp*c(14))
               go to 185
  180       continue
!           case 3 - after two or more successive failures
               c(14) = .5_kp*c(14)
  185       continue
!
!           check against hmax
            c(14) = min(c(14), c(16))
!
!           check against hmin
            c(14) = max(c(14), c(13))
!
!***********interrupt no 1 (with ind=4) if requested
            if (c(8) .eq. 0._kp) go to 1111
               ind = 4
               return
!           resume here on re-entry with ind .eq. 4   ........re-entry..
 1111       continue
!
!           calculate hmag, xtrial - depending on preliminary hmag, xend
            if (c(14) .ge. abs(xend - x)) go to 190
!              do not step more than half way to xend
               c(14) = min(c(14), .5_kp*abs(xend - x))
               c(17) = x + sign(c(14), xend - x)
               go to 195
  190       continue
!              hit xend exactly
               c(14) = abs(xend - x)
               c(17) = xend
  195       continue
!
!           calculate htrial
            c(18) = c(17) - x
!
!        end stage 1
!
!        ***************************************************************
!        * stage 2 - calculate ytrial (adding 7 to no of  fcn  evals). *
!        * w(*,2), ... w(*,8)  hold  intermediate  results  needed  in *
!        * stage 3. w(*,9) is temporary storage until finally it holds *
!        * ytrial.                                                     *
!        ***************************************************************
!
            temp = c(18)/1398169080000._kp
!
            do 200 k = 1, n
               w(k,9) = y(k) + temp*w(k,1)*233028180000._kp
  200       continue
            call fcn(n, x + c(18)/6._kp, w(1,9), w(1,2),extradata)
!
            do 205 k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*74569017600._kp &
                    + w(k,2)*298276070400._kp  )
  205       continue
            call fcn(n, x + c(18)*(4._kp/15._kp), w(1,9), w(1,3),extradata)
!
            do 210 k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*1165140900000._kp &
                    - w(k,2)*3728450880000._kp &
                    + w(k,3)*3495422700000._kp )
  210       continue
            call fcn(n, x + c(18)*(2._kp/3._kp), w(1,9), w(1,4),extradata)
!
            do 215 k = 1, n
               w(k,9) = y(k) + temp*( - w(k,1)*3604654659375._kp &
                    + w(k,2)*12816549900000._kp &
                    - w(k,3)*9284716546875._kp &
                    + w(k,4)*1237962206250._kp )
  215       continue
            call fcn(n, x + c(18)*(5._kp/6._kp), w(1,9), w(1,5),extradata)
!
            do 220 k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*3355605792000._kp &
                    - w(k,2)*11185352640000._kp &
                    + w(k,3)*9172628850000._kp &
                    - w(k,4)*427218330000._kp &
                    + w(k,5)*482505408000._kp  )
  220       continue
            call fcn(n, x + c(18), w(1,9), w(1,6),extradata)
!
            do 225 k = 1, n
               w(k,9) = y(k) + temp*( - w(k,1)*770204740536._kp &
                    + w(k,2)*2311639545600._kp &
                    - w(k,3)*1322092233000._kp &
                    - w(k,4)*453006781920._kp &
                    + w(k,5)*326875481856._kp  )
  225       continue
            call fcn(n, x + c(18)/15._kp, w(1,9), w(1,7),extradata)
!
            do 230 k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*2845924389000._kp &
                    - w(k,2)*9754668000000._kp &
                    + w(k,3)*7897110375000._kp &
                    - w(k,4)*192082660000._kp &
                    + w(k,5)*400298976000._kp &
                    + w(k,7)*201586000000._kp  )
  230       continue
            call fcn(n, x + c(18), w(1,9), w(1,8),extradata)
!
!           calculate ytrial, the extrapolated approximation and store
!              in w(*,9)
            do 235 k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*104862681000._kp &
                    + w(k,3)*545186250000._kp &
                    + w(k,4)*446637345000._kp &
                    + w(k,5)*188806464000._kp &
                    + w(k,7)*15076875000._kp &
                    + w(k,8)*97599465000._kp   )
  235       continue
!
!           add 7 to the no of fcn evals
            c(24) = c(24) + 7._kp
!
!        end stage 2
!
!        ***************************************************************
!        * stage 3 - calculate the error estimate est. first calculate *
!        * the  unweighted  absolute  error  estimate vector (per unit *
!        * step) for the unextrapolated approximation and store it  in *
!        * w(*,2).  then  calculate the weighted max norm of w(*,2) as *
!        * specified by the error  control  indicator  c(1).  finally, *
!        * modify  this result to produce est, the error estimate (per *
!        * unit step) for the extrapolated approximation ytrial.       *
!        ***************************************************************
!
!           calculate the unweighted absolute error estimate vector
            do 300 k = 1, n
               w(k,2) = (   w(k,1)*8738556750._kp &
                    + w(k,3)*9735468750._kp &
                    - w(k,4)*9709507500._kp &
                    + w(k,5)*8582112000._kp &
                    + w(k,6)*95329710000._kp &
                    - w(k,7)*15076875000._kp &
                    - w(k,8)*97599465000._kp)/1398169080000._kp
  300       continue
!
!           calculate the weighted max norm of w(*,2) as specified by
!           the error control indicator c(1)
            temp = 0._kp
            if (c(1) .ne. 1._kp) go to 310
!              absolute error control
               do 305 k = 1, n
                  temp = max(temp,abs(w(k,2)))
  305          continue
               go to 360
  310       if (c(1) .ne. 2._kp) go to 320
!              relative error control
               do 315 k = 1, n
                  temp = max(temp, abs(w(k,2)/y(k)))
  315          continue
               go to 360
  320       if (c(1) .ne. 3._kp) go to 330
!              weights are 1/max(c(2),abs(y(k)))
               do 325 k = 1, n
                  temp = max(temp, abs(w(k,2)) &
                       / max(c(2), abs(y(k))) )
  325          continue
               go to 360
  330       if (c(1) .ne. 4._kp) go to 340
!              weights are 1/max(c(k+30),abs(y(k)))
               do 335 k = 1, n
                  temp = max(temp, abs(w(k,2)) &
                       / max(c(k+30), abs(y(k))) )
  335          continue
               go to 360
  340       if (c(1) .ne. 5._kp) go to 350
!              weights are 1/c(k+30)
               do 345 k = 1, n
                  temp = max(temp, abs(w(k,2)/c(k+30)))
  345          continue
               go to 360
  350       continue
!              default case - weights are 1/max(1,abs(y(k)))
               do 355 k = 1, n
                  temp = max(temp, abs(w(k,2)) &
                       / max(1._kp, abs(y(k))) )
  355          continue
  360       continue
!
!           calculate est - (the weighted max norm of w(*,2))*hmag*scale
!              - est is intended to be a measure of the error  per  unit
!              step in ytrial
            c(19) = temp*c(14)*c(15)
!
!        end stage 3
!
!        ***************************************************************
!        * stage 4 - make decisions.                                   *
!        ***************************************************************
!
!           set ind=5 if step acceptable, else set ind=6
            ind = 5
            if (c(19) .gt. tol) ind = 6
!
!***********interrupt no 2 if requested
            if (c(9) .eq. 0._kp) go to 2222
               return
!           resume here on re-entry with ind .eq. 5 or 6   ...re-entry..
 2222       continue
!
            if (ind .eq. 6) go to 410
!              step accepted (ind .eq. 5), so update x, y from xtrial,
!                 ytrial, add 1 to the no of successful steps, and set
!                 the no of successive failures to zero
               x = c(17)
               do 400 k = 1, n
                  y(k) = w(k,9)
  400          continue
               c(22) = c(22) + 1._kp
               c(23) = 0._kp
!**************return(with ind=3, xend saved, flag set) if x .eq. xend
               if (x .ne. xend) go to 405
                  ind = 3
                  c(20) = xend
                  c(21) = 1._kp
                  return
  405          continue
               go to 420
  410       continue
!              step not accepted (ind .eq. 6), so add 1 to the no of
!                 successive failures
               c(23) = c(23) + 1._kp
!**************error return (with ind=-3) if hmag .le. hmin
               if (c(14) .gt. c(13)) go to 415
                  ind = -3
                  return
  415          continue
  420       continue
!
!        end stage 4
!
      go to 99999
!     end loop
!
!  begin abort action
  500 continue
      write(6,*)'Computation stopped in dverk with'
      write(6,*)'ind= tol= ',ind,tol
      write(6,*)'x= n= ',x,n
      write(6,*)'c(13)= xend= ',c(13),xend
      write(6,*)'nw= c(16)= c(20)= ',nw, c(16),c(20)
      write(6,*)'c(22)= c(23)= c(24)= ',c(22),c(23),c(24)
      write(6,*)'y(:)= ',y
!      write(6,*) ind, tol, x, n, c(13), xend, nw, c(16), c(20), &
!           c(22), c(23), c(24), (y(k), k = 1, n)
!  505 format( /// 1h0, 58hcomputation stopped in dverk with the following values - / 1h0, 5hind =, i4, 5x, 6htol  =, 1pd13.6, 5x, 11hx         =,&
!           1pd22.15&
!           / 1h , 5hn   =, i4, 5x, 6hhmin =, 1pd13.6, 5x, 11hxend      =,&
!           1pd22.15&
!           / 1h , 5hnw  =, i4, 5x, 6hhmax =, 1pd13.6, 5x, 11hprev xend =,&
!           1pd22.15&
!           / 1h0, 14x, 27hno of successful steps    =, 0pf8.0&
!           / 1h , 14x, 27hno of successive failures =, 0pf8.0&
!           / 1h , 14x, 27hno of function evals      =, 0pf8.0&
!           / 1h0, 23hthe components of y are&
!           // (1h , 1p5d24.15)                                           )
!
      stop
!
!  end abort action
!
    end subroutine dverk


    



    function zbrent(func,x1,x2,tol,extradata)
      INTEGER ITMAX
      real(kp) zbrent,tol,x1,x2,EPS
      type(transfert) :: extradata 
      
      !real(kp) :: func
!      EXTERNAL func
      PARAMETER (ITMAX=1000000,EPS=1d-15)
      INTEGER iter
      real(kp) a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm      
      
      integer, parameter :: iexmax = 1000
      real(kp), parameter :: tolExpand = 1e-2
      integer :: iex
      logical :: notbracketed
      
      interface
         function func(x,otherdata)
           use fieldprec, only : kp,transfert         
           implicit none                     
           real(kp) :: func
           real(kp), intent(in) :: x
           type(transfert), optional, intent(inout) :: otherdata
         end function func
      end interface

      if (x1.gt.x2) then
         stop 'zbrent: min > max on input!'
      endif
      

      notbracketed = .true.
      iex = 0

      a=x1
      b=x2      

      do while (notbracketed.and.(iex.lt.iexmax))
         fa=func(a,extradata)
         fb=func(b,extradata)
         if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.)) then
            write(*,*)'x1=',a,'f(x1)=',fa
            write(*,*)'x2=',b,'f(x2)=',fb            
            write(*,*)'zbrent: expanding interval!'
            a = a - abs(a)*tolExpand
            b = b + abs(b)*tolExpand
            iex = iex + 1
            notbracketed = .true.
         else
            notbracketed = .false.
         endif
      enddo

      if (notbracketed) stop 'root must be bracketed for zbrent'

      c=b
      fc=fb
      do iter=1,ITMAX
        if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
          c=a
          fc=fa
          d=b-a
          e=d
        endif
        if(abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2.*EPS*abs(b)+0.5*tol
        xm=.5*(c-b)
        if(abs(xm).le.tol1 .or. fb.eq.0.)then
          zbrent=b
          return
        endif
        if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa
          if(a.eq.c) then
            p=2.*xm*s
            q=1.-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          endif
          if(p.gt.0.) q=-q
          p=abs(p)
          if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
          else
            d=xm
            e=d
          endif
        else
          d=xm
          e=d
        endif
        a=b
        fa=fb
        if(abs(d) .gt. tol1) then
          b=b+d
        else
          b=b+sign(tol1,xm)
        endif
        fb=func(b,extradata)
     enddo
     !stop 'zbrent exceeding maximum iterations'
     write(*,*) 'zbrent exceeding maximum iterations'
     zbrent=b
     return
   end function zbrent







 end module infsolvers


