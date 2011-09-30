!provides the slow-roll initial field values such that N e-folds of
!inflation occurs.

module infsric
  use infprec, only : kp
  use infbgmodel, only : matterNum
  implicit none

  private

  logical, parameter :: display = .true.
  
  integer, parameter :: efoldBound = 110._kp

  
  public slowroll_initial_matter

contains

  

  function slowroll_initial_matter(infParam,efoldWanted)   
    use infbgmodel, only : infbgparam
    implicit none
    type(infbgparam), intent(in) :: infParam
    real(kp), optional, intent(in) :: efoldWanted
    real(kp), dimension(matterNum) :: slowroll_initial_matter

    real(kp) :: efold
    real(kp), parameter :: efoldDefault = efoldBound
    real(kp), dimension(matterNum) :: matterEnd, matterIni
    

    if (present(efoldWanted)) then
       efold = efoldWanted
    else
       efold = efoldDefault
    endif
    
    select case (infParam%name)

    case ('largef')
       matterIni = 2._kp

    end select
    

    if (display) then
       write(*,*)
       write(*,*)'slowroll_initial_matter: (*kappa)'
       write(*,*)'matterEnd = ',matterEnd
       write(*,*)'matterIni = ',matterIni
       write(*,*)
    endif

    slowroll_initial_matter = matterIni

  end function slowroll_initial_matter


end module infsric
