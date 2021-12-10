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




module infbgspline
!splinning and interpolation routines for the background. Apart for
!test they are only used to determine rapidly the bfold of quantum
!mode creation

  use fieldprec, only : kp
  implicit none

  private

!for debugging
  logical, parameter :: display = .false.
  logical, parameter :: dump_file = .false.


!bad idea, force the compiler to make copy of pointer arrays at each
!call of bspline (involves explicit shape array)
!  type splineinfbg
!     integer :: order
!     integer :: bcoefNum
!     integer :: fieldNum
!     real(kp), dimension(:), pointer :: bfoldKnot, epsilon1JFBcoef
!     real(kp), dimension(:), pointer :: hubbleBcoef, epsilon1Bcoef
!     real(kp), dimension(:,:), pointer :: fieldBcoef, fieldDotBcoef
!  end type splineinfbg
  

!to be initialized explicitely from bg computation
  
!  type(splineinfbg) :: splineBg

  integer, save :: order
  integer, save :: bcoefNum
  integer, save :: splineNum
  real(kp), dimension(:), allocatable, save :: bfoldKnot, epsilon1JFBcoef
  real(kp), dimension(:), allocatable, save :: hubbleBcoef, epsilon1Bcoef, epsilon2Bcoef
  real(kp), dimension(:,:), allocatable, save :: fieldBcoef, fieldDotBcoef

  

  public free_infbg_spline,set_infbg_spline, check_infbg_spline
  public splineval_hubble_parameter,splineval_epsilon1,splineval_epsilon1JF
  public splineval_epsilon2,splineval_field,splineval_fielddot

contains



  function check_infbg_spline()
    implicit none
    logical :: check_infbg_spline

    check_infbg_spline = (allocated(bfoldKnot) &
         .or. allocated(epsilon1JFBcoef) &
         .or. allocated(hubbleBcoef) &
         .or. allocated(epsilon1Bcoef) &
         .or. allocated(epsilon2Bcoef) &
         .or. allocated(fieldBcoef) &
         .or. allocated(fieldDotBcoef))
    
  end function check_infbg_spline




  subroutine free_infbg_spline()    
    implicit none    
    if (check_infbg_spline()) then       
       deallocate(bfoldKnot)
       deallocate(hubbleBcoef)
       deallocate(epsilon1Bcoef)
       deallocate(epsilon2Bcoef)
       deallocate(epsilon1JFBcoef)
       deallocate(fieldBcoef)
       deallocate(fieldDotBcoef)       
       if (display) write(*,*)'free_infbg_spline: infbg spline freed'
    else
       stop 'free_infbg_spline: not infbg spline data allocated'
    endif      
   
  end subroutine free_infbg_spline





  subroutine set_infbg_spline(bgEnd,ptrFirstBgdata)
!initialize the spline and set efoldIni, efoldEnd
    use infbgmodel, only : fieldNum
    use infbg, only : infbgdata, infbgphys
    use infbg, only : count_infbg_data
    use bspline
    implicit none 
    type(infbgphys), intent(in) :: bgEnd
    type(infbgdata), pointer :: ptrFirstBgdata
       
    integer :: bgdataNum

    integer :: i
    type(infbgdata), pointer :: ptrRun

!the background computation need to go further than epsilon1 > 1 (spline).
    real(kp), allocatable, dimension(:) :: bcoefTemp, dataTemp
    real(kp), allocatable, dimension(:) :: bfold, hubble, epsilon1, epsilon2, epsilon1JF
    real(kp), allocatable, dimension(:,:) :: field, fieldDot

    i = 0
    ptrRun => null()

    if (check_infbg_spline()) then
       stop 'set_infbg_spline: splineBg already associated'
    endif
    
    if (.not.associated(ptrFirstBgData)) then
       stop 'set_infbg_spline: no data found!'
    endif

    bgdataNum = count_infbg_data(ptrFirstBgdata)
    if (display) then
       write(*,*)
       write(*,*)'-------------------set_infbg_spline-------------------'
       write(*,*)
       write(*,*)'bgdataNum = ',bgdataNum
    endif
    order = 3
    bcoefNum = bgdataNum
    splineNum = fieldNum

    allocate(bfold(bgdataNum))
    allocate(hubble(bgdataNum))
    allocate(epsilon1(bgdataNum))
    allocate(epsilon2(bgdataNum))
    allocate(epsilon1JF(bgdataNum))
    allocate(field(bgdataNum,fieldNum))
    allocate(fieldDot(bgdataNum,fieldNum))

    allocate(bfoldKnot(bgdataNum + order))
    allocate(hubbleBcoef(bcoefNum))
    allocate(epsilon1Bcoef(bcoefNum))
    allocate(epsilon2Bcoef(bcoefNum))
    allocate(epsilon1JFBcoef(bcoefNum))
    allocate(fieldBcoef(bcoefNum,splineNum))
    allocate(fieldDotBcoef(bcoefNum,splineNum))
    
    if (associated(ptrFirstBgdata)) then
       ptrRun => ptrFirstBgdata
       i = 0       
       do while (associated(ptrRun))      
          i = i + 1
          bfold(i) = ptrRun%bg%efold - bgEnd%efold
          if (bfold(i).eq.0.) bfold(i) = epsilon(1.0_kp)
          hubble(i) = ptrRun%bg%hubble
          epsilon1(i) = ptrRun%bg%epsilon1
          epsilon2(i) = ptrRun%bg%epsilon2
          epsilon1JF(i) = ptrRun%bg%epsilon1JF
          field(i,:) = ptrRun%bg%field(:)
          fieldDot(i,:) = ptrRun%bg%fieldDot(:)
          ptrRun => ptrRun%ptr               
       enddo      
         
       if (display) then
          write(*,*)
          write(*,*)'assoNum = ',i
          write(*,*)'bgdataNum = ',bgdataNum
          write(*,*)'splinning range:'
          write(*,*)'bfoldIni = bfoldEndData = ',bfold(1),bfold(i)
          write(*,*)
          write(*,*)'------------------------------------------------------'
          write(*,*)
       endif
       ptrRun => null()
    else
       stop 'no bgdata found'
    endif
 
    call dbsnak(bgdataNum,bfold,order,bfoldKnot)
    call dbsint(bgdataNum,bfold,hubble,order,bfoldKnot &
         ,hubbleBcoef)
    call dbsint(bgdataNum,bfold,epsilon1,order,bfoldKnot &
         ,epsilon1Bcoef)
    call dbsint(bgdataNum,bfold,epsilon2,order,bfoldKnot &
         ,epsilon2Bcoef)
    call dbsint(bgdataNum,bfold,epsilon1JF,order,bfoldKnot &
         ,epsilon1JFBcoef)

    allocate(bcoefTemp(bgDataNum))
    allocate(dataTemp(bgDataNum))

    do i=1,splineNum
       dataTemp(:) = field(:,i)
       call dbsint(bgdataNum,bfold,dataTemp,order,bfoldKnot,bcoefTemp)
       fieldBcoef(:,i) = bcoefTemp

       dataTemp(:) = fieldDot(:,i)
       call dbsint(bgdataNum,bfold,dataTemp,order,bfoldKnot,bcoefTemp)
       fieldDotBcoef(:,i) = bcoefTemp
    enddo

    deallocate(bcoefTemp)
    deallocate(dataTemp)

    deallocate(bfold,hubble,epsilon1,epsilon2,epsilon1JF,field,fieldDot)

  end subroutine set_infbg_spline







  function splineval_hubble_parameter(bfold)   
    use bspline
    implicit none
       
    real(kp) :: splineval_hubble_parameter    
    real(kp), intent(in) :: bfold
        
    splineval_hubble_parameter = dbsval(bfold,order,bfoldKnot &
         ,bcoefNum,hubbleBcoef)
   
  end function splineval_hubble_parameter



   function splineval_epsilon1(bfold)   
    use bspline
    implicit none
       
    real(kp) :: splineval_epsilon1    
    real(kp), intent(in) :: bfold
        
    splineval_epsilon1 = dbsval(bfold,order,bfoldKnot &
         ,bcoefNum,epsilon1Bcoef)
   
  end function splineval_epsilon1

  
  function splineval_epsilon2(bfold)   
    use bspline
    implicit none
       
    real(kp) :: splineval_epsilon2   
    real(kp), intent(in) :: bfold
        
    splineval_epsilon2 = dbsval(bfold,order,bfoldKnot &
         ,bcoefNum,epsilon2Bcoef)
   
  end function splineval_epsilon2



  function splineval_epsilon1JF(bfold)   
    use bspline
    implicit none
          
    real(kp) :: splineval_epsilon1JF    
    real(kp), intent(in) :: bfold
        
    splineval_epsilon1JF = dbsval(bfold,order,bfoldKnot &
         ,bcoefNum,epsilon1JFBcoef)
   
  end function splineval_epsilon1JF




  function splineval_field(bfold)   
    use bspline
    implicit none
         
    real(kp), dimension(splineNum) :: splineval_field
    real(kp), intent(in) :: bfold
    real(kp), dimension(:), allocatable :: bcoefTemp
    integer :: i,bgDataNum
    bgDataNum = bcoefNum

    allocate(bcoefTemp(bgDataNum))

    i = 0
    do i=1,splineNum
       bcoefTemp(:) = fieldBcoef(:,i)
       splineval_field(i) = dbsval(bfold,order,bfoldKnot &
            ,bcoefNum,bcoefTemp)
    enddo
   
    deallocate(bcoefTemp)

  end function splineval_field




  function splineval_fielddot(bfold)
!reminder: this is the efold(bfold) fieldDot = Dfield/De(b)fold
    use bspline
    implicit none
       
    real(kp), dimension(splineNum) :: splineval_fielddot    
    real(kp), intent(in) :: bfold
    real(kp), dimension(:), allocatable :: bcoefTemp
    integer :: i,bgDataNum
    bgDataNum = bcoefNum

    allocate(bcoefTemp(bgDataNum))

    i = 0
    do i=1,splineNum
       bcoefTemp(:) = fieldDotBcoef(:,i)
       splineval_fielddot(i) = dbsval(bfold,order,bfoldKnot &
            ,bcoefNum,bcoefTemp)
    enddo
   
    deallocate(bcoefTemp)

  end function splineval_fielddot



end module infbgspline
