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





module infpowspline
!compute the primordial power spectra for few values of k and set a
!spline at intermediate k

  use fieldprec, only : kp

  implicit none

  private

  logical, parameter :: display = .false.

 
  integer, save :: scalOrder
  integer, save :: scalBcoefNum
  integer, save :: scalModeNum
  real(kp), dimension(:), allocatable, save :: scalLnkmpcKnot
  real(kp), dimension(:,:,:), allocatable, save :: scalPowerBcoef
!$omp threadprivate(scalOrder,scalBcoefNum,scalModeNum)
!$omp threadprivate(scalLnkmpcKnot,scalPowerBcoef)
  
  integer, save :: tensOrder
  integer, save :: tensBcoefNum
  integer, save :: tensModeNum
  real(kp), dimension(:), allocatable, save :: tensLnkmpcKnot
  real(kp), dimension(:), allocatable, save :: tensPowerBcoef  
!$omp threadprivate(tensOrder,tensBcoefNum,tensModeNum)
!$omp threadprivate(tensLnkmpcKnot,tensPowerBcoef)

  public check_power_scal_spline, check_power_tens_spline
  public free_power_scal_spline, free_power_tens_spline
  public set_power_scal_spline, set_power_tens_spline
  public splineval_power_scal, splineval_power_tens
  public splineval_ns_scal, splineval_nt_tens
  public splineval_alphas_scal, splineval_betas_scal

contains
  

  function check_power_scal_spline()
    implicit none
    logical :: check_power_scal_spline

    check_power_scal_spline = (allocated(scalLnkmpcKnot) &
         .or. allocated(scalPowerBcoef))

  end function check_power_scal_spline
      



  subroutine free_power_scal_spline()
    implicit none

    if (check_power_scal_spline()) then
       deallocate(scalLnkmpcKnot)
       deallocate(scalPowerBcoef)
       if (display) write(*,*) 'free_power_scal_spline: powerscalar spline freed'
    else
       write(*,*) 'free_power_scal_spline: not powerscalar spline data allocated'
    endif

  end subroutine free_power_scal_spline




  function check_power_tens_spline()
    implicit none
    logical :: check_power_tens_spline

    check_power_tens_spline = (allocated(tensLnkmpcKnot) &
         .or. allocated(tensPowerBcoef))

  end function check_power_tens_spline




  subroutine free_power_tens_spline()
    implicit none

    if (check_power_tens_spline()) then
       deallocate(tensLnkmpcKnot)
       deallocate(tensPowerBcoef)
       if (display) write(*,*) 'free_power_tens_spline: powertensor spline freed'
    else
       write(*,*) 'free_power_tens_spline: not powertensor spline data allocated'
    endif

  end subroutine free_power_tens_spline




  subroutine set_power_scal_spline(infCosmo,lnkmpcVec)
    use bspline, only : dbsnak, dbsint
    use inftorad, only : inftoradcosmo
    use infpert, only : scalNum, power_spectrum_scal
    implicit none

    type(inftoradcosmo), intent(in) :: infCosmo
    real(kp), dimension(:), intent(in) :: lnkmpcVec
    
    real(kp), dimension(:), allocatable :: lnPowerScalTemp
    real(kp), dimension(:), allocatable :: scalPowerBcoefTemp
    real(kp), dimension(:,:,:), allocatable :: lnPowerScal
    real(kp), dimension(scalNum,scalNum) :: powerScal
    real(kp) :: kmpc
    
    integer :: i,j,k,lnkmpcNum

    lnkmpcNum = size(lnkmpcVec,1)

    scalOrder = 3
    scalBcoefNum = lnkmpcNum
    scalModeNum = scalNum
    
!spline allocation
    if (check_power_scal_spline()) then
       write(*,*) 'set_power_scal_spline: spline already allocated'
       stop
    endif

    allocate(scalLnkmpcKnot(lnkmpcNum + scalOrder))
    allocate(scalPowerBcoef(scalBcoefNum,scalNum,scalNum))

!local
    allocate(lnPowerScal(lnkmpcNum,scalNum,scalNum))

    if (display) then
       write(*,*)'set_power_scal_spline: computing knots'
    endif

!$omp parallel do &
!$omp default(shared) &
!$omp private(i,kmpc,powerScal,j,k) &
!$omp schedule(dynamic)    
    do i=1,lnkmpcNum       
       kmpc = exp(lnkmpcVec(i))
       powerScal = power_spectrum_scal(infCosmo,kmpc)
       do k=1,scalNum
          do j=1,scalNum
             if (powerScal(j,k).eq.0._kp) then
                powerScal(j,k) = tiny(1.0_kp)             
             endif
             lnPowerScal(i,j,k) = log(powerScal(j,k))   
          enddo
       enddo
    enddo
!$omp end parallel do    

    if (display) then
       write(*,*)'set_power_scal_spline: end computing knots'
    endif

    call dbsnak(lnkmpcNum,lnkmpcVec,scalOrder,scalLnkmpcKnot)

!to avoid race conditions    
!    allocate(scalPowerBcoefTemp(scalBcoefNum))
!    allocate(lnPowerScalTemp(lnkmpcNum))

!    do j=1,scalNum
!       do i=1,scalNum
!          lnPowerScalTemp(:) = lnPowerScal(:,i,j)
!          call dbsint(lnkmpcNum,lnkmpcVec,lnPowerScalTemp,scalOrder,scalLnkmpcKnot &
!               ,scalPowerBcoefTemp)
!          scalPowerBcoef(:,i,j) = scalPowerBcoefTemp(:)
!       enddo
!    enddo

    do j=1,scalNum
       do i=1,scalNum
          call dbsint(lnkmpcNum,lnkmpcVec,lnPowerScal(:,i,j),scalOrder,scalLnkmpcKnot &
               ,scalPowerBcoef(:,i,j))
       enddo
    enddo
    
!    deallocate(scalPowerBcoefTemp)
!    deallocate(lnPowerScalTemp)
    deallocate(lnPowerScal)
    

  end subroutine set_power_scal_spline





  subroutine set_power_tens_spline(infCosmo,lnkmpcVec)
    use bspline, only : dbsnak, dbsint
    use infbgmodel, only : fieldNum
    use inftorad, only : inftoradcosmo
    use infpert, only : scalNum, power_spectrum_tens   
    implicit none

    type(inftoradcosmo), intent(in) :: infCosmo
    real(kp), dimension(:), intent(in) :: lnkmpcVec
    
    real(kp), dimension(:), allocatable :: lnPowerTens
    real(kp) :: kmpc, powerTens
    integer :: i,lnkmpcNum

  
    lnkmpcNum = size(lnkmpcVec,1)
    tensOrder = 3
    tensBcoefNum = lnkmpcNum
    tensModeNum = 1
    
!spline allocation
    if (check_power_tens_spline()) then       
       write(*,*) 'set_power_tens_spline: spline already allocated'
       stop
    endif

    allocate(tensLnkmpcKnot(lnkmpcNum + tensOrder))
    allocate(tensPowerBcoef(tensBcoefNum))
   
!local
    allocate(lnPowerTens(lnkmpcNum))

    if (display) then
       write(*,*)'set_power_tens_spline: computing knots'
    endif


!$omp parallel do &
!$omp default(shared) &
!$omp private(i,kmpc,powerTens) &
!$omp schedule(dynamic)
    do i=1,lnkmpcNum       
       kmpc = exp(lnkmpcVec(i))
       powerTens = power_spectrum_tens(infCosmo,kmpc)
       if (powerTens.eq.0._kp) then
          powerTens = tiny(1.0_kp)
       endif
       lnPowerTens(i) = log(powerTens)     
    enddo
!$omp end parallel do    

    if (display) then
       write(*,*)'set_power_scal_spline: end computing knots'
    endif

    call dbsnak(lnkmpcNum,lnkmpcVec,tensOrder,tensLnkmpcKnot)
    
    call dbsint(lnkmpcNum,lnkmpcVec,lnPowerTens,tensOrder,tensLnkmpcKnot &
         ,tensPowerBcoef)

    deallocate(lnPowerTens)
    
  end subroutine  set_power_tens_spline





  function splineval_power_scal(kmpc)
    use bspline, only : dbsval
    implicit none
    real(kp), intent(in) :: kmpc
    real(kp), dimension(scalModeNum,scalModeNum) :: splineval_power_scal    
    real(kp) :: lnkmpc
    real(kp), dimension(:), allocatable :: scalPowerBcoefTemp
    integer :: i,j

    lnkmpc = log(kmpc)  
    
!avoid race condition
!    allocate(scalPowerBcoefTemp(scalBcoefNum))

!    do j=1,scalModeNum
!       do i=1,scalModeNum
!          scalPowerBcoefTemp = scalPowerBcoef(:,i,j)
!          splineval_power_scal(i,j) = exp(dbsval(lnkmpc,scalOrder,scalLnkmpcKnot &
!               ,scalBcoefNum,scalPowerBcoefTemp))
!       enddo
!    enddo

    do j=1,scalModeNum
       do i=1,scalModeNum
          splineval_power_scal(i,j) = exp(dbsval(lnkmpc,scalOrder,scalLnkmpcKnot &
               ,scalBcoefNum,scalPowerBcoef(:,i,j)))
       enddo
    enddo

    
!    deallocate(scalPowerBcoefTemp)

  end function splineval_power_scal




  function splineval_power_tens(kmpc)
    use bspline, only : dbsval
    implicit none
    real(kp), intent(in) :: kmpc
    real(kp) :: splineval_power_tens   
    real(kp) :: lnkmpc

    lnkmpc = log(kmpc)

    splineval_power_tens = exp(dbsval(lnkmpc,tensOrder,tensLnkmpcKnot &
         ,tensBcoefNum,tensPowerBcoef))
    
  end function splineval_power_tens

!the derivative of the spline at kmpc (+1)
  function splineval_ns_scal(kmpc)
    use bspline, only : dbsder
    implicit none
    real(kp), intent(in) :: kmpc
    real(kp), dimension(scalModeNum,scalModeNum) :: splineval_ns_scal    
    real(kp) :: lnkmpc
    integer :: i,j

!first derivative: spectral index    
    integer, parameter :: ider = 1
    
    lnkmpc = log(kmpc)  


    do j=1,scalModeNum
       do i=1,scalModeNum
          splineval_ns_scal(i,j) = dbsder(ider,lnkmpc,scalOrder,scalLnkmpcKnot &
               ,scalBcoefNum,scalPowerBcoef(:,i,j))
       enddo
    enddo

!we have computed d ln(Ps)/d ln(k) which is ns-1    
    splineval_ns_scal = splineval_ns_scal + 1._kp

  end function splineval_ns_scal


  function splineval_nt_tens(kmpc)
    use bspline, only : dbsder
    implicit none
    real(kp), intent(in) :: kmpc
    real(kp) :: splineval_nt_tens   
    real(kp) :: lnkmpc

!first derivative: spectral index    
    integer, parameter :: ider = 1
    
    lnkmpc = log(kmpc)

!d ln(Ph) / d ln(k) which is nt    
    splineval_nt_tens = dbsder(ider,lnkmpc,tensOrder,tensLnkmpcKnot &
         ,tensBcoefNum,tensPowerBcoef)

  end function splineval_nt_tens



  function splineval_alphas_scal(kmpc)
    use bspline, only : dbsder
    implicit none
    real(kp), intent(in) :: kmpc
    real(kp), dimension(scalModeNum,scalModeNum) :: splineval_alphas_scal    
    real(kp) :: lnkmpc
    integer :: i,j

!second derivative: running
    integer, parameter :: ider = 2
    
    lnkmpc = log(kmpc)  


    do j=1,scalModeNum
       do i=1,scalModeNum
          splineval_alphas_scal(i,j) = dbsder(ider,lnkmpc,scalOrder,scalLnkmpcKnot &
               ,scalBcoefNum,scalPowerBcoef(:,i,j))
       enddo
    enddo


  end function splineval_alphas_scal


  function splineval_betas_scal(kmpc)
    use bspline, only : dbsder
    implicit none
    real(kp), intent(in) :: kmpc
    real(kp), dimension(scalModeNum,scalModeNum) :: splineval_betas_scal    
    real(kp) :: lnkmpc
    integer :: i,j

!third derivative: running of running
    integer, parameter :: ider = 3
    
    lnkmpc = log(kmpc)  


    do j=1,scalModeNum
       do i=1,scalModeNum
          splineval_betas_scal(i,j) = dbsder(ider,lnkmpc,scalOrder,scalLnkmpcKnot &
               ,scalBcoefNum,scalPowerBcoef(:,i,j))
       enddo
    enddo

    
  end function splineval_betas_scal
  


end module infpowspline
