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


#ifndef NOSRMODELS
   
    
    select case (infParam%name)

    case ('largef')
       slowroll_initial_matter = lfi_initial_field(infParam,efold)

    case ('mixlf')
       slowroll_initial_matter = mlfi_initial_field(infParam,efold)

    case ('rcmass')
       slowroll_initial_matter = rcmi_initial_field(infParam,efold)
       
    case ('rcquad')
       slowroll_initial_matter = rcqi_initial_field(infParam,efold)

    case ('natinf')
       if (infParam%consts(2).eq.1._kp) then
          slowroll_initial_matter = pni_initial_field(infParam,efold)
       elseif (infParam%consts(2).eq.-1._kp) then
          slowroll_initial_matter = mni_initial_field(infParam,efold)
       else
          stop 'slowroll_initial_matter: natural inflation not found!'
       endif

    case ('exsusy')
       slowroll_initial_matter = esi_initial_field(infParam,efold)

    case ('powlaw')
       slowroll_initial_matter = pli_initial_field(infParam,efold)

    case ('kahmod')
       if (infParam%consts(3).eq.1._kp) then
          slowroll_initial_matter = kmii_initial_field(infParam,efold)
       elseif (infParam%consts(3).eq.4._kp/3._kp) then
          stop 'khamo2 soon!'
       else
          stop 'slowroll_initial_matter: khaler moduli not found!'
       endif

    case ('hfline')
       slowroll_initial_matter = hf1i_initial_field(infParam,efold)

!    case ('smallf')
!       slowroll_initial_matter = sfi_initial_field(infParam,efold)

!    case ('kksf')
!       slowroll_initial_matter = sfbi_initial_field(infParam,efold)

!    case ('kklt')
!       slowroll_initial_matter = bi_initial_field(infParam,efold)



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

#else

    write(*,*) 'slowroll_initial_matter: FieldInf not compiled against libsrmodels!'
    stop 'You have to set non zero initial field values explicitly!'

#endif

  end function slowroll_initial_matter


#ifndef NOSRMODELS


  function lfi_initial_field(infParam,efold)
    use lfisr, only : lfi_x_endinf,lfi_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: lfi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: p, xEnd, xIni, bfold

    bfold = -efold
    p = infParam%consts(2)

    xEnd = lfi_x_endinf(p)

    if (display) write(*,*)'lf_initial_field: xend= ',xEnd

    xIni = lfi_x_trajectory(bfold,xEnd,p)

    lfi_initial_field(:) = xIni

  end function lfi_initial_field



  function mlfi_initial_field(infParam,efold)
    use mlfisr, only : mlfi_x_endinf,mlfi_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: mlfi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: p, q, alpha, xEnd, xIni, bfold

    bfold = -efold


    p = infParam%consts(2)
    q = infParam%consts(4)
    alpha = infParam%consts(3)
    
    xEnd = mlfi_x_endinf(p,q,alpha)

    if (display) write(*,*)'mlfi_initial_field: xend= ',xEnd

    xIni = mlfi_x_trajectory(bfold,xEnd,p,q,alpha)

    mlfi_initial_field(:) = xIni

  end function mlfi_initial_field




  function rcmi_initial_field(infParam,efold)
    use rcmisr, only : rcmi_x_endinf,rcmi_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: rcmi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alpha, xEnd, xIni, bfold

    bfold = -efold
    
    alpha = 0.5_kp * infParam%consts(2)
    
    xEnd = rcmi_x_endinf(alpha)

    if (display) write(*,*)'rcmi_initial_field: xend= ',xEnd

    xIni = rcmi_x_trajectory(bfold,xEnd,alpha)

    rcmi_initial_field(:) = xIni

  end function rcmi_initial_field




  function rcqi_initial_field(infParam,efold)
    use rcqisr, only : rcqi_x_endinf,rcqi_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: rcqi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alpha, xEnd, xIni, bfold

    bfold = -efold
    
    alpha = infParam%consts(2)
    
    xEnd = rcqi_x_endinf(alpha)

    if (display) write(*,*)'rcqi_initial_field: xend= ',xEnd

    xIni = rcqi_x_trajectory(bfold,xEnd,alpha)

    rcqi_initial_field(:) = xIni

  end function rcqi_initial_field



  
  function pni_initial_field(infParam,efold)
    use pnisr, only : pni_x_endinf,pni_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: pni_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: mu, xEnd, xIni, bfold

    bfold = -efold
    
    mu = infParam%consts(3)
    
    xEnd = pni_x_endinf(mu)

    if (display) write(*,*)'pni_initial_field: xend= ',xEnd

    xIni = pni_x_trajectory(bfold,xEnd,mu)

    pni_initial_field(:) = xIni

  end function pni_initial_field



  
  function mni_initial_field(infParam,efold)
    use mnisr, only : mni_x_endinf,mni_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: mni_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: mu, xEnd, xIni, bfold

    bfold = -efold
    
    mu = infParam%consts(3)
    
    xEnd = mni_x_endinf(mu)

    if (display) write(*,*)'mni_initial_field: xend= ',xEnd

    xIni = mni_x_trajectory(bfold,xEnd,mu)

    mni_initial_field(:) = xIni

  end function mni_initial_field


  
  function esi_initial_field(infParam,efold)
    use esisr, only : esi_x_endinf,esi_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: esi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: q, xEnd, xIni, bfold

    bfold = -efold
    
    q = infParam%consts(2)
    
    xEnd = esi_x_endinf(q)

    if (display) write(*,*)'esi_initial_field: xend= ',xEnd

    xIni = esi_x_trajectory(bfold,xEnd,q)

    esi_initial_field(:) = xIni

  end function esi_initial_field



  
  function pli_initial_field(infParam,efold)
    use plisr, only : pli_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: pli_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alpha, XIni, xEnd, bfold
    real(kp), dimension(2) :: fieldStop

    bfold = -efold

    alpha = infParam%consts(2)

    fieldStop = field_stopinf(infParam)
       
    xEnd = fieldStop(1)

    if (display) write(*,*)'pli_initial_field: xend= ',xEnd
   
    xIni = pli_x_trajectory(bfold,xEnd,alpha)

    pli_initial_field(:) = xIni 

  end function pli_initial_field


  
  function kmii_initial_field(infParam,efold)
    use kmiisr, only : kmii_x_endinf,kmii_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: kmii_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alpha, xEnd, xIni, bfold

    bfold = -efold
    
    alpha = infParam%consts(2)
    
    xEnd = kmii_x_endinf(alpha)

    if (display) write(*,*)'kmii_initial_field: xend= ',xEnd

    xIni = kmii_x_trajectory(bfold,xEnd,alpha)

    kmii_initial_field(:) = xIni

  end function kmii_initial_field




  function hf1i_initial_field(infParam,efold)
    use hf1isr, only : hf1i_x_endinf,hf1i_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: hf1i_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alpha, xEnd, xIni, bfold

    bfold = -efold
    
    alpha = infParam%consts(2)
    
    xEnd = hf1i_x_endinf(alpha)

    if (display) write(*,*)'hf1i_initial_field: xend= ',xEnd

    xIni = hf1i_x_trajectory(bfold,xEnd,alpha)

    hf1i_initial_field(:) = xIni

  end function hf1i_initial_field



#ifdef NOYET
  function sfi_initial_field(infParam,efold)
    use sfisr, only : sfi_x_endinf,sfi_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: sfi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: p, mu, xEnd, xIni, bfold

    bfold = -efold

    p = infParam%consts(2)
    mu = infParam%consts(3)

    if (mu.le.0._kp) stop 'sf_initial_field: improper mu<0!'
    if (p.lt.2._kp) stop 'sf_initial_field: improper p<2!'

    xEnd = sfi_x_endinf(p,mu)

    if (display) write(*,*)'sf_initial_field: xend= ',xEnd

    xIni = sfi_x_trajectory(bfold,xEnd,p,mu)

    sfi_initial_field(:) = xIni * mu

  end function sfi_initial_field



  function sfbi_initial_field(infParam,efold)
    use sfbisr, only : sfbi_x_endinf,sfbi_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: sfbi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: p, mu, xEnd, xIni, bfold

    bfold = -efold

    p = infParam%consts(2)
    mu = infParam%consts(3)

    if (mu.le.0._kp) stop 'kksf_initial_field: improper mu<0!'
    if (p.lt.0._kp) stop 'kksf_initial_field: improper p<0!'

    xEnd = sfbi_x_endinf(p,mu)

    if (display) write(*,*)'kksf_initial_field: xend= ',xEnd

    xIni = kksf_x_trajectory(bfold,xEnd,p,mu)

    sfbi_initial_field(:) = xIni * mu

  end function sfbi_initial_field
  


  
  function bi_initial_field(infParam,efold)
    use bisr, only : bi_x_endinf,bi_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: bi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: p, mu, xEnd, xEps, xIni, xUv, xStrg, bfold
    real(kp), dimension(2) :: fieldUv, fieldStop

    bfold = -efold

    p = infParam%consts(2)
    mu = infParam%consts(3)

    if (mu.le.0._kp) stop 'bi_initial_field: improper mu<0!'
    if (p.lt.2._kp) stop 'bi_initial_field: improper p<0!'
    if (infParam%consts(5).ne.-1._kp) then
       stop 'bi_initial_field: improper consts(5)'
    endif

    fieldUv = field_thbound(infParam)
    fieldStop = field_stopinf(infParam)

    xUv = fieldUv(1)/mu
    xStrg = fieldStop(1)/mu

    xEps = bi_x_endinf(p,mu)

!if phistrg occurs before slow-roll violation
    xEnd = max(xStrg,xEps)

    if (display) write(*,*)'bi_initial_field: xend= ',xEnd

    if (xEnd.gt.xUv) then
       write(*,*)'bi_initial_field: xEnd > XUv!'
    endif

    xIni = bi_x_trajectory(bfold,xEnd,p,mu)

    bi_initial_field(:) = xIni * mu

  end function bi_initial_field

  
#endif  


#endif

end module infsric
