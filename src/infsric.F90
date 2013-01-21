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


#ifndef NOASPIC
   
    
    select case (infParam%name)

    case ('largef')
       slowroll_initial_matter = lfi_initial_field(infParam,efold)    

    case ('gmixlf')
       slowroll_initial_matter = gmlfi_initial_field(infParam,efold)

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
       elseif ((infParam%consts(3).eq.4._kp/3._kp) &
            .and.(infParam%consts(5).eq.4._kp/3._kp)) then
          slowroll_initial_matter = kmiii_initial_field(infParam,efold)
       else
          stop 'slowroll_initial_matter: kahler moduli not found!'
       endif

    case ('hfline')
       slowroll_initial_matter = hf1i_initial_field(infParam,efold)

    case ('smallf')
       slowroll_initial_matter = sfi_initial_field(infParam,efold)

    case ('gswli')
       slowroll_initial_matter = li_initial_field(infParam,efold)

    case ('interm')
       slowroll_initial_matter = ii_initial_field(infParam,efold)

    case ('colwei')
       slowroll_initial_matter = cwi_initial_field(infParam,efold)

    case ('higgsi')
       slowroll_initial_matter = hi_initial_field(infParam,efold)

    case ('twisti')
       slowroll_initial_matter = twi_initial_field(infParam,efold)

    case ('tdwell')
       slowroll_initial_matter = dwi_initial_field(infParam,efold)

    case ('mhitop')
       slowroll_initial_matter = mhi_initial_field(infParam,efold)

!    case ('kksf')
!       slowroll_initial_matter = sfbi_initial_field(infParam,efold)

!    case ('kklt')
!       slowroll_initial_matter = bi_initial_field(infParam,efold)

    case ('logmd1')
       slowroll_initial_matter = lmi1_initial_field(infParam,efold)

    case ('logmd2')
       slowroll_initial_matter = lmi2_initial_field(infParam,efold)

    case ('ricci1')
       slowroll_initial_matter = rpi1_initial_field(infParam,efold)

    case ('ricci2')
       slowroll_initial_matter = rpi2_initial_field(infParam,efold)

    case ('betexp')
       slowroll_initial_matter = bei_initial_field(infParam,efold)

    case ('radiag')
       slowroll_initial_matter = rgi_initial_field(infParam,efold)

    case ('nmssmi')
       slowroll_initial_matter = mssmi_initial_field(infParam,efold)

    case ('rinfpt')
       slowroll_initial_matter = ripi_initial_field(infParam,efold)

    case ('gmssmi')
       slowroll_initial_matter = gmssmi_initial_field(infParam,efold)

    case ('bsusyb')
       slowroll_initial_matter = bsusybi_initial_field(infParam,efold)

    case ('tipinf')
       slowroll_initial_matter = ti_initial_field(infParam,efold)

    case ('psenat')
       slowroll_initial_matter = psni_initial_field(infParam,efold)

    case ('nckahi')
       slowroll_initial_matter = ncki_initial_field(infParam,efold)

    case ('hybrid')
       slowroll_initial_matter = vhi_initial_field(infParam,efold)

    case ('dysusy')
       slowroll_initial_matter = dsi_initial_field(infParam,efold)

    case ('arctan')
       slowroll_initial_matter = ai_initial_field(infParam,efold)

    case ('fixnsa')
       slowroll_initial_matter = cnai_initial_field(infParam,efold)

    case ('fixnsb')
       slowroll_initial_matter = cnbi_initial_field(infParam,efold)

    case ('fixnsc')
       slowroll_initial_matter = cnci_initial_field(infParam,efold)

    case ('fixnsd')
       slowroll_initial_matter = cndi_initial_field(infParam,efold)

    case ('nszero')
       slowroll_initial_matter = csi_initial_field(infParam,efold)

    case('oifold')
       slowroll_initial_matter = oi_initial_field(infParam,efold)

    case ('sugrab')
       slowroll_initial_matter = sbi_initial_field(infParam,efold)

    case ('sneus1')
       slowroll_initial_matter = ssi1_initial_field(infParam,efold)

    case ('sneus2')
       slowroll_initial_matter = ssi2_initial_field(infParam,efold)

    case ('sneus3')
       slowroll_initial_matter = ssi3_initial_field(infParam,efold)

    case ('sneus4')
       slowroll_initial_matter = ssi4_initial_field(infParam,efold)

    case ('sneus5')
       slowroll_initial_matter = ssi5_initial_field(infParam,efold)

    case ('sneus6')
       slowroll_initial_matter = ssi6_initial_field(infParam,efold)

    case ('runma1')
       slowroll_initial_matter = rmi1_initial_field(infParam,efold)

    case ('runma2')
       slowroll_initial_matter = rmi2_initial_field(infParam,efold)

    case ('runma3')
       slowroll_initial_matter = rmi3_initial_field(infParam,efold)

    case ('runma4')
       slowroll_initial_matter = rmi4_initial_field(infParam,efold)

    case ('logpo1')
       slowroll_initial_matter = lpi1_initial_field(infParam,efold)

    case ('logpo2')
       slowroll_initial_matter = lpi2_initial_field(infParam,efold)

    case ('logpo3')
       slowroll_initial_matter = lpi3_initial_field(infParam,efold)

!    case ('f-term')
!       slowroll_initial_matter = fterm_initial_field(infParam,efold)

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
    write(*,*)
    write(*,*)'slowroll_initial_matter: '
    write(*,*)'FieldInf not built with libaspic support!'
    write(*,*)'initial conditions functions are provided for the largef models only'
    if (infParam%name.ne.'largef') then
       stop 'You have to set non zero initial field values explicitly!'
    endif

    slowroll_initial_matter = largef_initial_field(infParam,efold)
    
#endif

  end function slowroll_initial_matter



  function largef_initial_field(infParam,efold)   
    use infbgmodel, only : infbgparam
    implicit none
    real(kp), dimension(matterNum) :: largef_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), optional, intent(in) :: efold
    
        
    real(kp) :: xEnd, xIni
    real(kp) :: p,bfold

    bfold = -efold
    p = infParam%consts(2)
   
    xEnd = p/sqrt(2._kp)
    xIni = sqrt(p**2/2._kp + 2._kp*p*efold)

    if (display) write(*,*)'largef_initial_field: xend= ',XEnd

    largef_initial_field(:) = xIni
  end function largef_initial_field



#ifndef NOASPIC


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



  function gmlfi_initial_field(infParam,efold)
    use gmlfisr, only : gmlfi_x_endinf,gmlfi_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: gmlfi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: p, q, alpha, xEnd, xIni, bfold

    bfold = -efold


    p = infParam%consts(2)
    q = infParam%consts(4)
    alpha = infParam%consts(3)
    
    xEnd = gmlfi_x_endinf(p,q,alpha)

    if (display) write(*,*)'gmlfi_initial_field: xend= ',xEnd

    xIni = gmlfi_x_trajectory(bfold,xEnd,p,q,alpha)

    gmlfi_initial_field(:) = xIni

  end function gmlfi_initial_field




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



  function li_initial_field(infParam,efold)
    use lisr, only : li_x_endinf,li_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: li_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alpha, xEnd, xIni, bfold

    bfold = -efold
    
    alpha = infParam%consts(2)
    
    xEnd = li_x_endinf(alpha)

    if (display) write(*,*)'li_initial_field: xend= ',xEnd

    xIni = li_x_trajectory(bfold,xEnd,alpha)

    li_initial_field(:) = xIni

  end function li_initial_field



  function ii_initial_field(infParam,efold)
    use iisr, only : ii_x_trajectory, ii_prior_xendmin
    implicit none
    real(kp), dimension(matterNum) :: ii_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: beta, XIni, xEnd, xEndMin, bfold
    real(kp), dimension(2) :: fieldStop

    bfold = -efold

    beta = infParam%consts(2)

    fieldStop = field_stopinf(infParam)
       
    xEnd = fieldStop(1)

    if (display) write(*,*)'ii_initial_field: xend= ',xEnd

    xEndMin = ii_prior_xendmin(beta,bfold)

    if (xEnd.lt.xEndMin) then
       write(*,*)'xEndMin= bfold= ',xEndMin,bfold
       stop 'ii_initial_field: fieldStop too small!'
    endif

    xIni = ii_x_trajectory(bfold,xEnd,beta)

    ii_initial_field(:) = xIni 

  end function ii_initial_field



  function kmiii_initial_field(infParam,efold)
    use kmiiisr, only : kmiii_x_endinf,kmiii_x_trajectory,kmiii_alphamin
    implicit none
    real(kp), dimension(matterNum) :: kmiii_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alphaMin
    real(kp) :: alpha, beta, xEnd, xIni, bfold

    bfold = -efold
    
    alpha = infParam%consts(2)
    beta = infParam%consts(4)

!model valid for alpha < beta e    

    if (alpha.ge.beta*exp(1._kp)) stop 'kmiii_initial_field: alpha > beta e!'

    alphaMin = kmiii_alphamin(beta)

    if (alpha.lt.alphaMin) then
       write(*,*)'beta= ',beta
       write(*,*)'alpha= alphamin= ',alpha,alphamin
       stop 'kmiii_initial_field: inflation never stops'
    endif

    xEnd = kmiii_x_endinf(alpha,beta)

    if (display) write(*,*)'kmiii_initial_field: xend= ',xEnd

    xIni = kmiii_x_trajectory(bfold,xEnd,alpha,beta)

    kmiii_initial_field(:) = xIni

  end function kmiii_initial_field



  function cwi_initial_field(infParam,efold)
    use cwisr, only : cwi_x_endinf,cwi_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: cwi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alpha, mu, xEnd, xIni, bfold

    bfold = -efold

    alpha = infParam%consts(2)
    mu = infParam%consts(3)

    if (mu.ne.(4._kp*exp(1._kp)/infParam%consts(2))**0.25_kp) then
       stop 'cwi_initial_field: improper Q value!'
    endif

    xEnd = cwi_x_endinf(alpha,mu)

    if (display) write(*,*)'cwi_initial_field: xend= ',xEnd

    xIni = cwi_x_trajectory(bfold,xEnd,alpha,mu)

    cwi_initial_field(:) = xIni

  end function cwi_initial_field
  
  

  function hi_initial_field(infParam,efold)
    use hisr, only : hi_x_endinf,hi_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: hi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: xEnd, xIni, bfold

    bfold = -efold
   
    xEnd = hi_x_endinf()

    if (display) write(*,*)'hi_initial_field: xend= ',xEnd

    xIni = hi_x_trajectory(bfold,xEnd)

    hi_initial_field(:) = xIni
    
  end function hi_initial_field



  function twi_initial_field(infParam,efold)
    use twisr, only : twi_x_trajectory, twi_x_endsr,phi0eps1
    implicit none
    real(kp), dimension(matterNum) :: twi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), dimension(2) :: fieldStop
    real(kp), intent(in) :: efold

    real(kp) :: mu, xEnd, xIni, xEps, bfold

    bfold = -efold
   
    mu = infParam%consts(3)
    
    fieldStop = field_stopinf(infParam)

    xEnd = fieldStop(1)
    if (mu.lt.phi0eps1) then
       xEps = twi_x_endsr(mu)
       if (xEnd.lt.xEps) then
          write(*,*) 'twi_initial_field: xEnd is in slow-roll violation region'
       endif
    endif

    if (display) write(*,*)'twi_initial_field: xend= ',xEnd,xEps

    xIni = twi_x_trajectory(bfold,xEnd,mu)

    twi_initial_field(:) = xIni                   

  end function twi_initial_field



 function dwi_initial_field(infParam,efold)
    use dwisr, only : dwi_x_endinf,dwi_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: dwi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: mu, xEnd, xIni, bfold

    bfold = -efold
   
    mu = infParam%consts(2)
    
    xEnd = dwi_x_endinf(mu)   
  
    if (display) write(*,*)'twi_initial_field: xend= ',xEnd

    xIni = dwi_x_trajectory(bfold,xEnd,mu)

    dwi_initial_field(:) = xIni*mu          

  end function dwi_initial_field



  function mhi_initial_field(infParam,efold)
    use mhisr, only : mhi_x_endinf,mhi_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: mhi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: mu, xEnd, xIni, bfold

    bfold = -efold
   
    mu = infParam%consts(2)
    
    xEnd = mhi_x_endinf(mu)   
  
    if (display) write(*,*)'mhi_initial_field: xend= ',xEnd

    xIni = mhi_x_trajectory(bfold,xEnd,mu)

    mhi_initial_field(:) = xIni*mu          

  end function mhi_initial_field



  function lmi1_initial_field(infParam,efold)
    use lmi1sr, only : lmi1_x_endinf,lmi1_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: lmi1_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alpha, beta, gamma
    real(kp) :: xEnd, xIni, bfold

    bfold = -efold
    
    alpha = infParam%consts(2)
    beta = infParam%consts(3)
    gamma = infParam%consts(4)

!model valid for alpha = 4(1-gamma)
    if (alpha.ne.4._kp*(1._kp - gamma)) stop 'lmi1_initial_field: alpha >< 4(1-gamma)'
   
    xEnd = lmi1_x_endinf(gamma,beta)

    if (display) write(*,*)'lmi1_initial_field: xend= ',xEnd

    xIni = lmi1_x_trajectory(bfold,xEnd,gamma,beta)

    lmi1_initial_field(:) = xIni

  end function lmi1_initial_field



  function lmi2_initial_field(infParam,efold)
    use lmi2sr, only : lmi2_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: lmi2_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold
    real(kp), dimension(2) :: fieldStop

    real(kp) :: alpha, beta, gamma
    real(kp) :: xEnd, xIni, bfold

    bfold = -efold
    
    alpha = infParam%consts(2)
    beta = infParam%consts(3)
    gamma = infParam%consts(4)

    fieldStop = field_stopinf(infParam)

    xEnd = fieldStop(1)
   
!model valid for alpha=4(1-gamma)
    if (alpha.ne.4._kp*(1._kp - gamma)) stop 'lmi2_initial_field: alpha >< 4(1-gamma)'
       
    if (display) write(*,*)'lmi2_initial_field: xend= ',xEnd

    xIni = lmi2_x_trajectory(bfold,xEnd,gamma,beta)

    lmi2_initial_field(:) = xIni

  end function lmi2_initial_field


  function rpi1_initial_field(infParam,efold)
    use rpi1sr, only : rpi1_x_endinf,rpi1_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: rpi1_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: p
    real(kp) :: xEnd, xIni, bfold

    bfold = -efold
    
    p = infParam%consts(2)
    
    if (p.lt.1._kp) stop 'rpi1_initial_field: p<1!'
   
    xEnd = rpi1_x_endinf(p)

    if (display) write(*,*)'rpi1_initial_field: xend= ',xEnd

    xIni = rpi1_x_trajectory(bfold,xEnd,p)

    rpi1_initial_field(:) = xIni * sqrt(3._kp/2._kp)

  end function rpi1_initial_field



  function rpi2_initial_field(infParam,efold)
    use rpi2sr, only : rpi2_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: rpi2_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: p
    real(kp) :: xEnd, xIni, bfold
    real(kp), dimension(2) :: fieldStop

    bfold = -efold
    
    p = infParam%consts(2)

    fieldStop = field_stopinf(infParam)
    xEnd = fieldStop(1)*sqrt(2._kp/3._kp)

    if (p.lt.1._kp) stop 'rpi2_initial_field: p<1!'
       
    if (display) write(*,*)'rpi2_initial_field: xend= ',xEnd

    xIni = rpi2_x_trajectory(bfold,xEnd,p)

    rpi2_initial_field(:) = xIni * sqrt(3._kp/2._kp)

  end function rpi2_initial_field

   

  function bei_initial_field(infParam, efold)
    use beisr, only : bei_x_trajectory, bei_x_endinf
    implicit none
    real(kp), dimension(matterNum) :: bei_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: bfold
    real(kp) :: xEnd, xIni
    real(kp) :: beta, lambda

    bfold = -efold
    beta = 1._kp/infParam%consts(3)
    lambda = infParam%consts(2)*infParam%consts(3)

    xEnd = bei_x_endinf(lambda,beta)

    if (display) write(*,*)'bei_initial_field: xend= ',xEnd

    xIni = bei_x_trajectory(bfold,xEnd,lambda,beta)

    bei_initial_field(:) = xIni

  end function bei_initial_field



  function rgi_initial_field(infParam, efold)
    use rgisr, only : rgi_x_trajectory, rgi_x_endinf
    implicit none
    real(kp), dimension(matterNum) :: rgi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: bfold
    real(kp) :: xEnd, xIni
    real(kp) :: alpha

    bfold = -efold
    alpha = infParam%consts(2)
   

    xEnd = rgi_x_endinf(alpha)

    if (display) write(*,*)'rgi_initial_field: xend= ',xEnd

    xIni = rgi_x_trajectory(bfold,xEnd,alpha)

    rgi_initial_field(:) = xIni

  end function rgi_initial_field


  function mssmi_initial_field(infParam, efold)
    use mssmisr, only : mssmi_x_trajectory, mssmi_x_endinf
    implicit none
    real(kp), dimension(matterNum) :: mssmi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: bfold
    real(kp) :: xEnd, xIni
    real(kp) :: alpha

    bfold = -efold
    alpha = infParam%consts(3)
   

    xEnd = mssmi_x_endinf(alpha)

    if (display) write(*,*)'mssmi_initial_field: xend= ',xEnd

    xIni = mssmi_x_trajectory(bfold,xEnd,alpha)

    mssmi_initial_field(:) = xIni

  end function mssmi_initial_field
  


  function ripi_initial_field(infParam, efold)
    use ripisr, only : ripi_x_trajectory, ripi_x_endinf
    implicit none
    real(kp), dimension(matterNum) :: ripi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: bfold
    real(kp) :: xEnd, xIni
    real(kp) :: alpha

    bfold = -efold
    alpha = infParam%consts(3)
   
    xEnd = ripi_x_endinf(alpha)

    if (display) write(*,*)'ripi_initial_field: xend= ',xEnd

    xIni = ripi_x_trajectory(bfold,xEnd,alpha)

    ripi_initial_field(:) = xIni

  end function ripi_initial_field



  function gmssmi_initial_field(infParam, efold)
    use gmssmisr, only : gmssmi_x_trajectory, gmssmi_x_endinf
    implicit none
    real(kp), dimension(matterNum) :: gmssmi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: bfold
    real(kp) :: xEnd, xIni
    real(kp) :: alpha, beta

    bfold = -efold
    alpha = infParam%consts(3)
    beta = infParam%consts(5)

    xEnd = gmssmi_x_endinf(alpha,beta)

    if (display) write(*,*)'gmssmi_initial_field: xend= ',xEnd

    xIni = gmssmi_x_trajectory(bfold,xEnd,alpha,beta)

    gmssmi_initial_field(:) = xIni

  end function gmssmi_initial_field



  function bsusybi_initial_field(infParam, efold)
    use bsusybisr, only : bsusybi_x_trajectory, bsusybi_xendmax
    implicit none
    real(kp), dimension(matterNum) :: bsusybi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: bfold
    real(kp) :: xEnd, xIni, xEndMax
    real(kp) :: gam
    real(kp), dimension(2) :: fieldStop


    bfold = -efold
    fieldStop = field_stopinf(infParam)
       
    gam = infParam%consts(4)/sqrt(6._kp)

    xEnd = fieldStop(1)
    xEndMax = bsusybi_xendmax(efold,gam)

    if (xEnd.gt.xEndMax) then
       write(*,*)'xEnd= xEndMax= ',xEnd,xEndMax
       stop 'bsusybi_initial_field: XEnd too large'
    endif

    if (display) write(*,*)'bsusybi_initial_field: xend= ',xEnd

    xIni = bsusybi_x_trajectory(bfold,xEnd,gam)

    bsusybi_initial_field(:) = xIni

  end function bsusybi_initial_field



  function ti_initial_field(infParam, efold)
    use tisr, only : ti_x_trajectory, ti_x_endinf
    implicit none
    real(kp), dimension(matterNum) :: ti_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: bfold
    real(kp) :: xEnd, xIni
    real(kp) :: alpha, mu

    bfold = -efold
    alpha = infParam%consts(2)
    mu = infParam%consts(3)

    xEnd = ti_x_endinf(alpha,mu)

    if (display) write(*,*)'ti_initial_field: xend= ',xEnd

    xIni = ti_x_trajectory(bfold,xEnd,alpha,mu)

    ti_initial_field(:) = xIni*mu

  end function ti_initial_field



  function psni_initial_field(infParam, efold)
    use psnisr, only : psni_x_trajectory, psni_x_endinf
    implicit none
    real(kp), dimension(matterNum) :: psni_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: bfold
    real(kp) :: xEnd, xIni
    real(kp) :: alpha, mu

    bfold = -efold
    alpha = infParam%consts(2)
    mu = infParam%consts(3)

    xEnd = psni_x_endinf(alpha,mu)

    if (display) write(*,*)'psni_initial_field: xend= ',xEnd

    xIni = psni_x_trajectory(bfold,xEnd,alpha,mu)

    psni_initial_field(:) = xIni*mu

  end function psni_initial_field



  function ncki_initial_field(infParam, efold)
    use nckisr, only : ncki_x_trajectory, ncki_x_endinf
    use nckisr, only : ncki_x_potmax
    implicit none
    real(kp), dimension(matterNum) :: ncki_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: bfold
    real(kp) :: xEnd, xIni, xPotMax
    real(kp) :: alpha, beta

    bfold = -efold
    alpha = infParam%consts(2)
    beta = infParam%consts(3)

    if (beta.le.0._kp) then
       xpotMax = ncki_x_potmax(alpha,beta)
       write(*,*)
       write(*,*)'ncki_initial_field: '
       write(*,*)'slow-roll violated, initial condition innaccurate!'
       write(*,*)'inflation occurs at the top: xpotMax= ',xpotMax
       write(*,*)
    endif

    xEnd = ncki_x_endinf(alpha,beta)

    if (display) write(*,*)'ncki_initial_field: xend= ',xEnd

    xIni = ncki_x_trajectory(bfold,xEnd,alpha,beta)

    ncki_initial_field(:) = xIni

  end function ncki_initial_field
  


  function vhi_initial_field(infParam, efold)
    use vhisr, only : vhi_x_trajectory, vhi_xendmax, vhi_xendmin
    implicit none
    real(kp), dimension(matterNum) :: vhi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: bfold
    real(kp) :: xEnd, xIni, xEndMax, xEndMin
    real(kp) :: p,mu
    real(kp), dimension(2) :: fieldStop


    bfold = -efold
    fieldStop = field_stopinf(infParam)
       
    p = infParam%consts(2)
    mu = infParam%consts(3)

    xEnd = fieldStop(1)/mu
    xEndMax = vhi_xendmax(efold,p,mu)
    xEndMin = vhi_xendmin(p,mu)

    if ((xEnd.gt.xEndMax).or.(xEnd.lt.xEndmin)) then
       write(*,*)'xEnd= xEndMin= xEndMax ',xEnd,xEndMin,xEndMax
       stop 'vhi_initial_field: XEnd out of bounds'
    endif

    if (display) write(*,*)'vhi_initial_field: xend= ',xEnd

    xIni = vhi_x_trajectory(bfold,xEnd,p,mu)

    vhi_initial_field(:) = xIni*mu

  end function vhi_initial_field

  

  function dsi_initial_field(infParam, efold)
    use dsisr, only : dsi_x_trajectory, dsi_xendmin
    implicit none
    real(kp), dimension(matterNum) :: dsi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: bfold
    real(kp) :: xEnd, xIni, xEndMin
    real(kp) :: p,mu
    real(kp), dimension(2) :: fieldStop


    bfold = -efold
    fieldStop = field_stopinf(infParam)
       
    p = infParam%consts(2)
    mu = infParam%consts(3)

    xEnd = fieldStop(1)/mu   
    xEndMin = dsi_xendmin(efold,p,mu)

    if ((xEnd.lt.xEndmin)) then
       write(*,*)'xEnd= xEndMin=  ',xEnd,xEndMin
       stop 'dsi_initial_field: xEnd out of bounds'
    endif

    if (display) write(*,*)'dsi_initial_field: xend= ',xEnd

    xIni = dsi_x_trajectory(bfold,xEnd,p,mu)

    dsi_initial_field(:) = xIni*mu

  end function dsi_initial_field


  function ai_initial_field(infParam, efold)
    use aisr, only : ai_x_trajectory, ai_x_endinf
    implicit none
    real(kp), dimension(matterNum) :: ai_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: bfold
    real(kp) :: xEnd, xIni
    real(kp) :: mu

    bfold = -efold
    mu = infParam%consts(2)
   
    xEnd = ai_x_endinf(mu)

    if (display) write(*,*)'ai_initial_field: xend= ',xEnd

    xIni = ai_x_trajectory(bfold,xEnd,mu)

    ai_initial_field(:) = xIni*mu

  end function ai_initial_field



  function cnai_initial_field(infParam, efold)
    use cnaisr, only : cnai_x_trajectory, cnai_x_endinf
    use cnaisr, only : cnai_x_potzero
    implicit none
    real(kp), dimension(matterNum) :: cnai_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: bfold
    real(kp) :: xEnd, xIni, xMax
    real(kp) :: alpha

    bfold = -efold
    alpha = infParam%consts(2)
   
    xMax = cnai_x_potzero(alpha)

    xEnd = cnai_x_endinf(alpha)

    if (display) write(*,*)'cnai_initial_field: xend= xmax=',xEnd,xMax

    xIni = cnai_x_trajectory(bfold,xEnd,alpha)

    cnai_initial_field(:) = xIni

  end function cnai_initial_field



  function cnbi_initial_field(infParam, efold)
    use cnbisr, only : cnbi_x_trajectory, cnbi_x_endinf
    implicit none
    real(kp), dimension(matterNum) :: cnbi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: bfold
    real(kp) :: xEnd, xIni
    real(kp) :: alpha

    bfold = -efold
    alpha = infParam%consts(2)
   
    xEnd = cnbi_x_endinf(alpha)

    if (display) write(*,*)'cnbi_initial_field: xend= ',xEnd

    xIni = cnbi_x_trajectory(bfold,xEnd,alpha)

    cnbi_initial_field(:) = xIni

  end function cnbi_initial_field



  function cnci_initial_field(infParam, efold)
    use cncisr, only : cnci_x_trajectory, cnci_xendmin
    implicit none
    real(kp), dimension(matterNum) :: cnci_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: bfold
    real(kp) :: xEnd, xIni, xEndMin
    real(kp) :: alpha
    real(kp), dimension(2) :: fieldStop

    bfold = -efold
    alpha = infParam%consts(2)
    fieldStop = field_stopinf(infParam)
    
    xEnd = fieldStop(1)
    xEndMin = cnci_xendmin(efold, alpha)

    if (xEnd.lt.xEndMin) then
       write(*,*)'xEnd= xEndMin= ',xEnd,xEndMin
       stop 'cnci_initial_field: XEnd too small'
    endif

    if (display) write(*,*)'cnci_initial_field: xend= ',xEnd
    
    xIni = cnci_x_trajectory(bfold,xEnd,alpha)

    cnci_initial_field(:) = xIni

  end function cnci_initial_field



  function cndi_initial_field(infParam, efold)
    use cndisr, only : cndi_x_trajectory, cndi_xendmax
    implicit none
    real(kp), dimension(matterNum) :: cndi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: bfold
    real(kp) :: xEnd, xIni, xEndMax
    real(kp) :: alpha, beta
    real(kp), dimension(2) :: fieldStop

    bfold = -efold
    alpha = infParam%consts(2)
    beta = infParam%consts(3)
    fieldStop = field_stopinf(infParam)
    
    xEnd = fieldStop(1)
    xEndMax = cndi_xendmax(efold, alpha, beta)

    if (xEnd.gt.xEndMax) then
       write(*,*)'xEnd= xEndMax= ',xEnd,xEndMax
       stop 'cndi_initial_field: XEnd too large'
    endif

    if (display) write(*,*)'cndi_initial_field: xend= ',xEnd
    
    xIni = cndi_x_trajectory(bfold,xEnd,alpha,beta)

    cndi_initial_field(:) = xIni

  end function cndi_initial_field



  function csi_initial_field(infParam,efold)
    use csisr, only : csi_x_trajectory, csi_xendmax
    implicit none
    real(kp), dimension(matterNum) :: csi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: bfold
    real(kp) :: xEnd, xIni, xEndMax
    real(kp) :: alpha
    real(kp), dimension(2) :: fieldStop

    bfold = -efold
    fieldStop = field_stopinf(infParam)
       
    alpha = infParam%consts(2)

    xEnd = fieldStop(1)
    xEndMax = csi_xendmax(efold,alpha)

    if ((xEnd.gt.xEndMax)) then
       write(*,*)'xEnd= xEndMax=  ',xEnd,xEndMax
       stop 'dsi_initial_field: xEnd out of bounds'
    endif

    if (display) write(*,*)'csi_initial_field: xend= ',xEnd

    xIni = csi_x_trajectory(bfold,xEnd,alpha)

    csi_initial_field(:) = xIni

  end function csi_initial_field



  function oi_initial_field(infParam,efold)
    use oisr, only : oi_x_trajectory, oi_x_endinf
    implicit none
    real(kp), dimension(matterNum) :: oi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: bfold
    real(kp) :: xEnd, xIni
    real(kp) :: alpha, mu

    bfold = -efold
           
    alpha = infParam%consts(2)
    mu = infParam%consts(3)

    xEnd = oi_x_endinf(alpha,mu)

    if (display) write(*,*)'oi_initial_field: xend= ',xEnd

    xIni = oi_x_trajectory(bfold,xEnd,alpha,mu)

    oi_initial_field(:) = xIni

  end function oi_initial_field


  function sbi_initial_field(infParam,efold)
    use sbisr, only : sbi_x_trajectory, sbi_x_endinf
    implicit none
    real(kp), dimension(matterNum) :: sbi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: bfold
    real(kp) :: xEnd, xIni
    real(kp) :: alpha, mu

    bfold = -efold
           
    alpha = infParam%consts(2)
    mu = infParam%consts(3)

    xEnd = sbi_x_endinf(alpha,mu)

    if (display) write(*,*)'sbi_initial_field: xend= ',xEnd

    xIni = sbi_x_trajectory(bfold,xEnd,alpha,mu)

    sbi_initial_field(:) = xIni

  end function sbi_initial_field

  

  function ssi1_initial_field(infParam,efold)
    use ssi1sr, only : ssi1_x_endinf,ssi1_x_trajectory
    use ssi1sr, only : ssi1_alphamin
    implicit none
    real(kp), dimension(matterNum) :: ssi1_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alpha, beta, alphamin
    real(kp) :: xEnd, xIni, bfold

    bfold = -efold
    
    alpha = infParam%consts(2)
    beta = infParam%consts(3)

    alphamin = ssi1_alphamin(beta)

    if (alpha.lt.alphamin) then
       write(*,*)'ssi1_initial_field: alpha= alphamin= ',alpha,alphamin
       stop
    endif
       
    xEnd = ssi1_x_endinf(alpha,beta)

    if (display) write(*,*)'ssi1_initial_field: xend= ',xEnd

    xIni = ssi1_x_trajectory(bfold,xEnd,alpha,beta)

    ssi1_initial_field(:) = xIni

  end function ssi1_initial_field



  function ssi2_initial_field(infParam,efold)
    use ssi2sr, only : ssi2_x_endinf,ssi2_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: ssi2_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alpha, beta
    real(kp) :: xEnd, xIni, bfold

    bfold = -efold
    
    alpha = infParam%consts(2)
    beta = infParam%consts(3)
  
    xEnd = ssi2_x_endinf(alpha,beta)

    if (display) write(*,*)'ssi2_initial_field: xend= ',xEnd

    xIni = ssi2_x_trajectory(bfold,xEnd,alpha,beta)

    ssi2_initial_field(:) = xIni

  end function ssi2_initial_field



  function ssi3_initial_field(infParam,efold)
    use ssi3sr, only : ssi3_x_endinf,ssi3_x_trajectory
    use ssi3sr, only : ssi3_alphamin
    implicit none
    real(kp), dimension(matterNum) :: ssi3_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alpha, beta, alphamin
    real(kp) :: xEnd, xIni, bfold

    bfold = -efold
    
    alpha = infParam%consts(2)
    beta = infParam%consts(3)

    alphamin = ssi3_alphamin(beta)

    if (alpha.lt.alphamin) then
       write(*,*)'ssi3_initial_field: alpha= alphamin= ',alpha,alphamin
       stop
    endif
       
    xEnd = ssi3_x_endinf(alpha,beta)

    if (display) write(*,*)'ssi3_initial_field: xend= ',xEnd

    xIni = ssi3_x_trajectory(bfold,xEnd,alpha,beta)

    ssi3_initial_field(:) = xIni

  end function ssi3_initial_field



  function ssi4_initial_field(infParam,efold)
    use ssi4sr, only : ssi4_x_endinf,ssi4_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: ssi4_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alpha, beta
    real(kp) :: xEnd, xIni, bfold

    bfold = -efold
    
    alpha = infParam%consts(2)
    beta = infParam%consts(3)
  
    xEnd = ssi4_x_endinf(alpha,beta)

    if (display) write(*,*)'ssi4_initial_field: xend= ',xEnd

    xIni = ssi4_x_trajectory(bfold,xEnd,alpha,beta)

    ssi4_initial_field(:) = xIni

  end function ssi4_initial_field



  function ssi5_initial_field(infParam,efold)
    use ssi5sr, only : ssi5_x_endinf,ssi5_x_trajectory
    use ssi5sr, only : ssi5_alphamax
    implicit none
    real(kp), dimension(matterNum) :: ssi5_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alpha, beta, alphamax
    real(kp) :: xEnd, xIni, bfold

    bfold = -efold
    
    alpha = infParam%consts(2)
    beta = infParam%consts(3)

    alphamax = ssi5_alphamax(beta)

    if (alpha.gt.alphamax) then
       write(*,*)'ssi5_initial_field: alpha= alphamax= ',alpha,alphamax
       stop
    endif
       
    xEnd = ssi5_x_endinf(alpha,beta)

    if (display) write(*,*)'ssi5_initial_field: xend= ',xEnd

    xIni = ssi5_x_trajectory(bfold,xEnd,alpha,beta)

    ssi5_initial_field(:) = xIni

  end function ssi5_initial_field



  function ssi6_initial_field(infParam,efold)
    use ssi6sr, only : ssi6_x_endinf,ssi6_x_trajectory
    use ssi6sr, only : ssi6_alphamax
    implicit none
    real(kp), dimension(matterNum) :: ssi6_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alpha, beta, alphamax
    real(kp) :: xEnd, xIni, bfold

    bfold = -efold
    
    alpha = infParam%consts(2)
    beta = infParam%consts(3)

    alphamax = ssi6_alphamax(beta)

    if (alpha.gt.alphamax) then
       write(*,*)'ssi6_initial_field: alpha= alphamax= ',alpha,alphamax
       stop
    endif
       
    xEnd = ssi6_x_endinf(alpha,beta)

    if (display) write(*,*)'ssi6_initial_field: xend= ',xEnd

    xIni = ssi6_x_trajectory(bfold,xEnd,alpha,beta)

    ssi6_initial_field(:) = xIni

  end function ssi6_initial_field



 function rmi1_initial_field(infParam,efold)
    use rmi1sr, only : rmi1_x_trajectory, rmi1_numacc_xendmax
    implicit none
    real(kp), dimension(matterNum) :: rmi1_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: bfold
    real(kp) :: xEnd, xIni, xEndMaxNum
    real(kp) :: c, mu
    real(kp), dimension(2) :: fieldStop

    bfold = -efold
    fieldStop = field_stopinf(infParam)
       
    mu = infParam%consts(3)
    c = 2._kp*infParam%consts(4)

    xEnd = fieldStop(1)/mu
    xEndMaxNum = rmi1_numacc_xendmax(efold,c,mu)

    if (xEnd.gt.xEndMaxNum) then
       write(*,*)'xEnd= xEndMaxNum=  ',xEnd,xEndMaxNum
       write(*,*)'rmi1_initial_field: xEnd is not reachable at that numerical accuracy!'
    endif

    if (display) write(*,*)'rmi1_initial_field: xend= ',xEnd

    xIni = rmi1_x_trajectory(bfold,xEnd,c,mu)

    rmi1_initial_field(:) = xIni * mu

  end function rmi1_initial_field



  function rmi2_initial_field(infParam,efold)
    use rmi2sr, only : rmi2_x_trajectory, rmi2_numacc_xendmin
    implicit none
    real(kp), dimension(matterNum) :: rmi2_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: bfold
    real(kp) :: xEnd, xIni, xEndMinNum
    real(kp) :: c, mu
    real(kp), dimension(2) :: fieldStop

    bfold = -efold
    fieldStop = field_stopinf(infParam)
       
    mu = infParam%consts(3)
    c = 2._kp*infParam%consts(4)

    xEnd = fieldStop(1)/mu
    xEndMinNum = rmi2_numacc_xendmin(efold,c,mu)

    if (xEnd.lt.xEndMinNum) then
       write(*,*)'xEnd= xEndMinNum=  ',xEnd,xEndMinNum
       write(*,*)'rmi2_initial_field: xEnd is not reachable at that numercal accuracy!'
    endif

    if (display) write(*,*)'rmi2_initial_field: xend= ',xEnd

    xIni = rmi2_x_trajectory(bfold,xEnd,c,mu)

    rmi2_initial_field(:) = xIni * mu

  end function rmi2_initial_field



  function rmi3_initial_field(infParam,efold)
    use rmi3sr, only : rmi3_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: rmi3_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: bfold
    real(kp) :: xEnd, xIni
    real(kp) :: c, mu
    real(kp), dimension(2) :: fieldStop

    bfold = -efold
    fieldStop = field_stopinf(infParam)
       
    mu = infParam%consts(3)
    c = 2._kp*infParam%consts(4)

    xEnd = fieldStop(1)/mu
  
    if (display) write(*,*)'rmi3_initial_field: xend= ',xEnd

    xIni = rmi3_x_trajectory(bfold,xEnd,c,mu)

    rmi3_initial_field(:) = xIni * mu

  end function rmi3_initial_field



  function rmi4_initial_field(infParam,efold)
    use rmi4sr, only : rmi4_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: rmi4_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: bfold
    real(kp) :: xEnd, xIni
    real(kp) :: c, mu
    real(kp), dimension(2) :: fieldStop

    bfold = -efold
    fieldStop = field_stopinf(infParam)
       
    mu = infParam%consts(3)
    c = 2._kp*infParam%consts(4)

    xEnd = fieldStop(1)/mu
  
    if (display) write(*,*)'rmi4_initial_field: xend= ',xEnd

    xIni = rmi4_x_trajectory(bfold,xEnd,c,mu)

    rmi4_initial_field(:) = xIni * mu

  end function rmi4_initial_field
  


  function lpi1_initial_field(infParam,efold)
    use lpi1sr, only : lpi1_x_endinf,lpi1_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: lpi1_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: p,q,mu
    real(kp) :: xEnd, xIni, bfold

    bfold = -efold
    
    p = infParam%consts(2)
    mu = infParam%consts(3)
    q = infParam%consts(4)

    xEnd = lpi1_x_endinf(p,q,mu)

    if (display) write(*,*)'lpi1_initial_field: xend= ',xEnd

    xIni = lpi1_x_trajectory(bfold,xEnd,p,q,mu)

    lpi1_initial_field(:) = xIni*mu

  end function lpi1_initial_field



  function lpi2_initial_field(infParam,efold)
    use lpi2sr, only : lpi2_x_endinf,lpi2_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: lpi2_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: p,q,mu
    real(kp) :: xEnd, xIni, bfold

    bfold = -efold
    
    p = infParam%consts(2)
    mu = infParam%consts(3)
    q = infParam%consts(4)

    xEnd = lpi2_x_endinf(p,q,mu)

    if (display) write(*,*)'lpi2_initial_field: xend= ',xEnd

    xIni = lpi2_x_trajectory(bfold,xEnd,p,q,mu)

    lpi2_initial_field(:) = xIni*mu

  end function lpi2_initial_field



  function lpi3_initial_field(infParam,efold)
    use lpi3sr, only : lpi3_x_endinf,lpi3_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: lpi3_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: p,q,mu
    real(kp) :: xEnd, xIni, bfold

    bfold = -efold
    
    p = infParam%consts(2)
    mu = infParam%consts(3)
    q = infParam%consts(4)

    xEnd = lpi3_x_endinf(p,q,mu)

    if (display) write(*,*)'lpi3_initial_field: xend= ',xEnd

    xIni = lpi3_x_trajectory(bfold,xEnd,p,q,mu)

    lpi3_initial_field(:) = xIni*mu

  end function lpi3_initial_field


 

#ifdef NOYET
  
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

#ifdef TWOFIELDS
  function fterm_initial_field(infParam,efold)    
    use infprec, only : pi
    implicit none
    real(kp), dimension(matterNum) :: fterm_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: kappa, mu

    real(kp) :: lambda, phic, M

    kappa = infParam%consts(2)
    mu = infParam%consts(3)

    lambda = kappa**2  * mu**4
    phic = sqrt(2._kp) * mu
    M = 2._kp * mu
    

    fterm_initial_field(1) = phic
    fterm_initial_field(2) = (lambda**(3._kp/2._kp) / (48._kp * pi**2._kp &
         * sqrt(pi * kappa**4._kp * log(2._kp)/(16._kp * pi**(2._kp)))))**(1._kp/2._kp)


  end function fterm_initial_field
#endif


end module infsric
