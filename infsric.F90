!provides the slow-roll initial field values such that N e-folds of
!inflation occurs.

module infsric
  use infprec, only : kp
  use infbgmodel, only : matterNum, infbgparam
  use infbounds, only : field_thbound, field_stopinf
  implicit none

  private

  logical, parameter :: display = .true.
  
  integer, parameter :: efoldBound = 110._kp

  
  public slowroll_initial_matter

contains

  

  function slowroll_initial_matter(infParam,efoldWanted)   
    implicit none
    type(infbgparam), intent(in) :: infParam
    real(kp), optional, intent(in) :: efoldWanted
    real(kp), dimension(matterNum) :: slowroll_initial_matter

    real(kp) :: efold
    real(kp), parameter :: efoldDefault = efoldBound
    real(kp), dimension(matterNum) :: matterEnd, matterIni
    
    slowroll_initial_matter = 0._kp
    
    if (present(efoldWanted)) then
       efold = efoldWanted
    else
       efold = efoldDefault
    endif
    
    select case (infParam%name)

    case ('largef')
       slowroll_initial_matter = lf_initial_field(infParam,efold)

    case ('smallf')
       slowroll_initial_matter = sf_initial_field(infParam,efold)

    case ('kksf')
       slowroll_initial_matter = kksf_initial_field(infParam,efold)

    case ('kklt')
       slowroll_initial_matter = kklt_initial_field(infParam,efold)

    case ('mixlf')
       slowroll_initial_matter = mixlf_initial_field(infParam,efold)

    case default
       stop 'slowroll_initial_matter: model not implemented!'

    end select
    

    if (display) then
       write(*,*)
       write(*,*)'slowroll_initial_matter: (*kappa)'
       write(*,*)'efold= ',efold
       write(*,*)'matterIni = ',slowroll_initial_matter
       write(*,*)
    endif

  end function slowroll_initial_matter



  function lf_initial_field(infParam,efold)
    use lfsrevol, only : lf_x_endinf,lf_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: lf_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: p, xEnd, xIni, bfold

    bfold = -efold
    p = infParam%consts(2)

    xEnd = lf_x_endinf(p)

    if (display) write(*,*)'lf_initial_field: xend= ',xEnd

    xIni = lf_x_trajectory(bfold,xEnd,p)

    lf_initial_field(:) = xIni

  end function lf_initial_field




  function sf_initial_field(infParam,efold)
    use sfsrevol, only : sf_x_endinf,sf_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: sf_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: p, mu, xEnd, xIni, bfold

    bfold = -efold

    p = infParam%consts(2)
    mu = infParam%consts(3)

    if (mu.le.0._kp) stop 'sf_initial_field: improper mu<0!'
    if (p.lt.2._kp) stop 'sf_initial_field: improper p<2!'

    xEnd = sf_x_endinf(p,mu)

    if (display) write(*,*)'sf_initial_field: xend= ',xEnd

    xIni = sf_x_trajectory(bfold,xEnd,p,mu)

    sf_initial_field(:) = xIni * mu

  end function sf_initial_field



  function kksf_initial_field(infParam,efold)
    use kksfsrevol, only : kksf_x_endinf,kksf_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: kksf_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: p, mu, xEnd, xIni, bfold

    bfold = -efold

    p = infParam%consts(2)
    mu = infParam%consts(3)

    if (mu.le.0._kp) stop 'kksf_initial_field: improper mu<0!'
    if (p.lt.0._kp) stop 'kksf_initial_field: improper p<0!'

    xEnd = kksf_x_endinf(p,mu)

    if (display) write(*,*)'kksf_initial_field: xend= ',xEnd

    xIni = kksf_x_trajectory(bfold,xEnd,p,mu)

    kksf_initial_field(:) = xIni * mu

  end function kksf_initial_field
  


  
  function kklt_initial_field(infParam,efold)
    use kkltsrevol, only : kklt_x_endinf,kklt_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: kklt_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: p, mu, xEnd, xEps, xIni, xUv, xStrg, bfold
    real(kp), dimension(2) :: fieldUv, fieldStop

    bfold = -efold

    p = infParam%consts(2)
    mu = infParam%consts(3)

    if (mu.le.0._kp) stop 'kklt_initial_field: improper mu<0!'
    if (p.lt.2._kp) stop 'kklt_initial_field: improper p<0!'
    if (infParam%consts(5).ne.-1._kp) then
       stop 'kklt_initial_field: improper consts(5)'
    endif

    fieldUv = field_thbound(infParam)
    fieldStop = field_stopinf(infParam)

    xUv = fieldUv(1)/mu
    xStrg = fieldStop(1)/mu

    xEps = kklt_x_endinf(p,mu)

!if phistrg occurs before slow-roll violation
    xEnd = max(xStrg,xEps)

    if (display) write(*,*)'kklt_initial_field: xend= ',xEnd

    if (xEnd.gt.xUv) then
       write(*,*)'kklt_initial_field: xEnd > XUv!'
    endif

    xIni = kklt_x_trajectory(bfold,xEnd,p,mu)

    kklt_initial_field(:) = xIni * mu

  end function kklt_initial_field



  function mixlf_initial_field(infParam,efold)
    use mixlfsrevol, only : mixlf_x_endinf,mixlf_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: mixlf_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: p, q, alpha, xEnd, xIni, bfold

    bfold = -efold

#ifndef PP5
    p = infParam%consts(2)
    q = infParam%consts(12)
    alpha = infParam%consts(6)
#else
    stop 'mixlf_initial_field: -DPP12 not defined!'
#endif
    

    xEnd = mixlf_x_endinf(p,q,alpha)

    if (display) write(*,*)'mixlf_initial_field: xend= ',xEnd

    xIni = mixlf_x_trajectory(bfold,xEnd,p,q,alpha)

    mixlf_initial_field(:) = xIni

  end function mixlf_initial_field
  


end module infsric
