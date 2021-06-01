!provides the slow-roll initial field values such that N e-folds of
!inflation occurs.

module infsric
  use fieldprec, only : kp
  use infbgmodel, only : matterNum, infbgparam
  use infbounds, only : field_thbound, field_stopinf
  implicit none

  private

  logical, parameter :: display = .true.
  
  integer, parameter :: efoldBound = 120._kp

  
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

    case ('largef','lfi')
       slowroll_initial_matter = lfi_initial_field(infParam,efold)    

    case ('gmixlf','gmlfi')
       slowroll_initial_matter = gmlfi_initial_field(infParam,efold)

    case ('rcmass','rcmi')
       slowroll_initial_matter = rcmi_initial_field(infParam,efold)
       
    case ('rcquad','rcqi')
       slowroll_initial_matter = rcqi_initial_field(infParam,efold)

    case ('natinf','ni')
       slowroll_initial_matter = ni_initial_field(infParam,efold)

    case ('hybnat','hni1')
       slowroll_initial_matter = hni1_initial_field(infParam,efold)
      
    case ('exsusy','esi')
       slowroll_initial_matter = esi_initial_field(infParam,efold)

    case ('powlaw','pli')
       slowroll_initial_matter = pli_initial_field(infParam,efold)

    case ('kahmod','kmi')
       if (infParam%consts(3).eq.1._kp) then
          slowroll_initial_matter = kmii_initial_field(infParam,efold)
       elseif ((infParam%consts(3).eq.4._kp/3._kp) &
            .and.(infParam%consts(5).eq.4._kp/3._kp)) then
          slowroll_initial_matter = kmiii_initial_field(infParam,efold)
       else
          stop 'slowroll_initial_matter: kahler moduli not found!'
       endif

    case ('hfline','hf1i')
       slowroll_initial_matter = hf1i_initial_field(infParam,efold)

    case ('smallf','sfi')
       slowroll_initial_matter = sfi_initial_field(infParam,efold)

    case ('gswli','li')
       slowroll_initial_matter = li_initial_field(infParam,efold)

    case ('interm','ii')
       slowroll_initial_matter = ii_initial_field(infParam,efold)

    case ('colwei','cwi')
       slowroll_initial_matter = cwi_initial_field(infParam,efold)

    case ('staroi','si')
       slowroll_initial_matter = si_initial_field(infParam,efold)

    case ('twisti','twi')
       slowroll_initial_matter = twi_initial_field(infParam,efold)

    case ('tdwell','dwi')
       slowroll_initial_matter = dwi_initial_field(infParam,efold)

    case ('gdwell','gdwi')
       slowroll_initial_matter = gdwi_initial_field(infParam,efold)       

    case ('mhitop','mhi')
       slowroll_initial_matter = mhi_initial_field(infParam,efold)

!    case ('kksf')
!       slowroll_initial_matter = sfbi_initial_field(infParam,efold)

!    case ('kklt')
!       slowroll_initial_matter = bi_initial_field(infParam,efold)

    case ('logmd1','lmi1')
       slowroll_initial_matter = lmi1_initial_field(infParam,efold)

    case ('logmd2','lmi2')
       slowroll_initial_matter = lmi2_initial_field(infParam,efold)

    case ('ricci1','rpi1')
       slowroll_initial_matter = rpi1_initial_field(infParam,efold)

    case ('ricci2','rpi2')
       slowroll_initial_matter = rpi2_initial_field(infParam,efold)

    case ('corsi1','ccsi1')
       slowroll_initial_matter = ccsi1_initial_field(infParam,efold)

    case ('corsi2','ccsi2')
       slowroll_initial_matter = ccsi2_initial_field(infParam,efold)

    case ('corsi3','ccsi3')
       slowroll_initial_matter = ccsi3_initial_field(infParam,efold)

    case ('betexp','bei')
       slowroll_initial_matter = bei_initial_field(infParam,efold)

    case ('radiag','rgi')
       slowroll_initial_matter = rgi_initial_field(infParam,efold)

    case ('nmssmi','mssmi')
       slowroll_initial_matter = mssmi_initial_field(infParam,efold)

    case ('orifpt','oripi')
       slowroll_initial_matter = oripi_initial_field(infParam,efold)

    case ('gmssmi')
       slowroll_initial_matter = gmssmi_initial_field(infParam,efold)

    case ('bsusyb','bsusybi')
       slowroll_initial_matter = bsusybi_initial_field(infParam,efold)

    case ('tipinf','ti')
       slowroll_initial_matter = ti_initial_field(infParam,efold)

    case ('psenat','psni')
       slowroll_initial_matter = psni_initial_field(infParam,efold)

    case ('nckahi','ncki')
       slowroll_initial_matter = ncki_initial_field(infParam,efold)

    case ('hybrid','vhi')
       slowroll_initial_matter = vhi_initial_field(infParam,efold)

    case ('dysusy','dsi')
       slowroll_initial_matter = dsi_initial_field(infParam,efold)

    case ('arctan','ai')
       slowroll_initial_matter = ai_initial_field(infParam,efold)

    case ('pureai','pai')
       slowroll_initial_matter = pai_initial_field(infParam,efold)
       
    case ('fixnsa','cnai')
       slowroll_initial_matter = cnai_initial_field(infParam,efold)

    case ('fixnsb','cnbi')
       slowroll_initial_matter = cnbi_initial_field(infParam,efold)

    case ('fixnsc','cnci')
       slowroll_initial_matter = cnci_initial_field(infParam,efold)

    case ('fixnsd','cndi')
       slowroll_initial_matter = cndi_initial_field(infParam,efold)

    case ('nszero','csi')
       slowroll_initial_matter = csi_initial_field(infParam,efold)

    case('oifold','oi')
       slowroll_initial_matter = oi_initial_field(infParam,efold)

    case ('sugrab','sbi')
       slowroll_initial_matter = sbi_initial_field(infParam,efold)

    case ('spsyb1','ssbi1')
       slowroll_initial_matter = ssbi1_initial_field(infParam,efold)

    case ('spsyb2','ssbi2')
       slowroll_initial_matter = ssbi2_initial_field(infParam,efold)

    case ('spsyb3','ssbi3')
       slowroll_initial_matter = ssbi3_initial_field(infParam,efold)

    case ('spsyb4','ssbi4')
       slowroll_initial_matter = ssbi4_initial_field(infParam,efold)

    case ('spsyb5','ssbi5')
       slowroll_initial_matter = ssbi5_initial_field(infParam,efold)

    case ('spsyb6','ssbi6')
       slowroll_initial_matter = ssbi6_initial_field(infParam,efold)

    case ('runma1','rmi1')
       slowroll_initial_matter = rmi1_initial_field(infParam,efold)

    case ('runma2','rmi2')
       slowroll_initial_matter = rmi2_initial_field(infParam,efold)

    case ('runma3','rmi3')
       slowroll_initial_matter = rmi3_initial_field(infParam,efold)

    case ('runma4','rmi4')
       slowroll_initial_matter = rmi4_initial_field(infParam,efold)

    case ('logpo1','lpi1')
       slowroll_initial_matter = lpi1_initial_field(infParam,efold)

    case ('logpo2','lpi2')
       slowroll_initial_matter = lpi2_initial_field(infParam,efold)

    case ('logpo3','lpi3')
       slowroll_initial_matter = lpi3_initial_field(infParam,efold)

    case ('ostach','osti')
       slowroll_initial_matter = osti_initial_field(infParam,efold)

    case ('witorh','wri')
       slowroll_initial_matter = wri_initial_field(infParam,efold)

    case ('invmon','imi')
       slowroll_initial_matter = imi_initial_field(infParam,efold)

    case ('nrifpt','ripi')
       slowroll_initial_matter = ripi_initial_field(infParam,efold)

    case ('grifpt','gripi')
       slowroll_initial_matter = gripi_initial_field(infParam,efold)

    case ('branei','bi')
       slowroll_initial_matter = bi_initial_field(infParam,efold)

    case ('kklmmt','kklti')
       slowroll_initial_matter = kklti_initial_field(infParam,efold)

    case ('nform1','nfi1')
       slowroll_initial_matter = nfi1_initial_field(infParam,efold)

    case ('nform2','nfi2')
       slowroll_initial_matter = nfi2_initial_field(infParam,efold)

    case ('nform3','nfi3')
       slowroll_initial_matter = nfi3_initial_field(infParam,efold)

    case ('nform4','nfi4')
       slowroll_initial_matter = nfi4_initial_field(infParam,efold)

    case ('dualsb','di')
       slowroll_initial_matter = di_initial_field(infParam,efold)

    case ('nrcoli','ncli')
       slowroll_initial_matter = ncli_initial_field(infParam,efold)

    case ('mukhai','vfmi')
       slowroll_initial_matter = vfmi_initial_field(infParam,efold)

    case ('axhtop','ahi')
       slowroll_initial_matter = ahi_initial_field(infParam,efold)

    case ('sbkahi','sbki')
       slowroll_initial_matter = sbki_initial_field(infParam,efold)

    case ('fibrei','fi')
       slowroll_initial_matter = fi_initial_field(infParam,efold)

    case ('sduali','sdi')
       slowroll_initial_matter = sdi_initial_field(infParam,efold)

    case ('scaaai','saai')
       slowroll_initial_matter = saai_initial_field(infParam,efold)

    case ('scaabi','sabi')
       slowroll_initial_matter = sabi_initial_field(infParam,efold)

    case ('scaaci','saci')
       slowroll_initial_matter = saci_initial_field(infParam,efold)

    case ('hyperb','hbi')
       slowroll_initial_matter = hbi_initial_field(infParam,efold)

    case ('smearh','shi')
       slowroll_initial_matter = shi_initial_field(infParam,efold)

    case ('rcinfp','rcipi')
       slowroll_initial_matter = rcipi_initial_field(infParam,efold)

    case ('saii1')
       slowroll_initial_matter = saii1_initial_field(infParam,efold)

    case ('saxone','saii2')
       slowroll_initial_matter = saii2_initial_field(infParam,efold)

    case ('saiii1')
       slowroll_initial_matter = saiii1_initial_field(infParam,efold)

    case ('saiii2')
       slowroll_initial_matter = saiii2_initial_field(infParam,efold)

    case ('saxtwo','saiii3')
       slowroll_initial_matter = saiii3_initial_field(infParam,efold)
       
    case ('dblexp','dei')
       slowroll_initial_matter = dei_initial_field(infParam,efold)
       
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


  
  function ni_initial_field(infParam,efold)
    use nisr, only : ni_x_endinf,ni_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: ni_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: mu, xEnd, xIni, bfold

    bfold = -efold
    
    mu = infParam%consts(2)
    
    xEnd = ni_x_endinf(mu)

    if (display) write(*,*)'ni_initial_field: xend= ',xEnd

    xIni = ni_x_trajectory(bfold,xEnd,mu)

    ni_initial_field(:) = xIni

  end function ni_initial_field



  function hni1_initial_field(infParam,efold)
    use hni1sr, only : hni1_x_endinf, hni1_x_trajectory, hni1_alphamin
    implicit none
    real(kp), dimension(matterNum) :: hni1_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alphaMin
    real(kp) :: alpha, mu, xEnd, xIni, bfold

    bfold = -efold

    mu = infParam%consts(2)
    alpha = infParam%consts(3)

    alphaMin = hni1_alphamin(mu)
         
    if (alpha.lt.alphamin) then
       write(*,*)'alphamin= alpha= ',alphaMin, alpha
       stop 'alpha too small to end inflation!'
    endif
   
    xEnd = hni1_x_endinf(alpha,mu)

    if (display) write(*,*)'hni1_initial_field: xend= ',xEnd

    xIni = hni1_x_trajectory(bfold,xEnd,alpha,mu)

    hni1_initial_field(:) = xIni*mu

  end function hni1_initial_field

  
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
    use iisr, only : ii_x_trajectory, ii_xendmin
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

    xEndMin = ii_xendmin(efold,beta)

    if (xEnd.lt.xEndMin) then
       write(*,*)'xEndMin= bfold= ',xEndMin,bfold
       stop 'ii_initial_field: fieldStop too small!'
    endif

    xIni = ii_x_trajectory(bfold,xEnd,beta)

    ii_initial_field(:) = xIni 

  end function ii_initial_field



  function kmiii_initial_field(infParam,efold)
    use kmiiisr, only : kmiii_x_endinf,kmiii_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: kmiii_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alpha, beta, xEnd, xIni, bfold

    bfold = -efold
    
    alpha = infParam%consts(2)
    beta = infParam%consts(4)

!model valid for alpha < beta e    

    if (alpha.ge.beta*exp(1._kp)) stop 'kmiii_initial_field: alpha > beta e!'

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

    if (alpha.ne.(4._kp*exp(1._kp))) then
       stop 'cwi_initial_field: improper alpha value!'
    endif

    xEnd = cwi_x_endinf(alpha,mu)

    if (display) write(*,*)'cwi_initial_field: xend= ',xEnd

    xIni = cwi_x_trajectory(bfold,xEnd,alpha,mu)

    cwi_initial_field(:) = xIni*mu

  end function cwi_initial_field
  
  

  function si_initial_field(infParam,efold)
    use sisr, only : si_x_endinf,si_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: si_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: xEnd, xIni, bfold

    bfold = -efold
   
    xEnd = si_x_endinf()

    if (display) write(*,*)'si_initial_field: xend= ',xEnd

    xIni = si_x_trajectory(bfold,xEnd)

    si_initial_field(:) = xIni
    
  end function si_initial_field



  function twi_initial_field(infParam,efold)
    use twisr, only : twi_x_trajectory, twi_x_epstwounity,phi0eps1
    implicit none
    real(kp), dimension(matterNum) :: twi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), dimension(2) :: fieldStop
    real(kp), intent(in) :: efold

    real(kp) :: mu, xEnd, xIni, xEpsTwo, bfold

    bfold = -efold
   
    mu = infParam%consts(3)
    
    fieldStop = field_stopinf(infParam)

    xEnd = fieldStop(1)
    if (mu.lt.phi0eps1) then
       xEpsTwo = twi_x_epstwounity(mu)
       if (xEnd.lt.xEpsTwo) then
          write(*,*) 'twi_initial_field: xEnd is in slow-roll violation region'
       endif
    endif

    if (display) write(*,*)'twi_initial_field: xend= ',xEnd,xEpsTwo

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

    if (mu.lt.sqrt(8._kp)) then
       write(*,*)'dwi_initial_field: slow-roll violated at x=0!'
    endif
    
    if (display) write(*,*)'dwi_initial_field: xend= ',xEnd

    xIni = dwi_x_trajectory(bfold,xEnd,mu)
    
    dwi_initial_field(:) = xIni*mu          

  end function dwi_initial_field



  function gdwi_initial_field(infParam,efold)
    use gdwisr, only : gdwi_x_endinf,gdwi_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: gdwi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: p, mu, xEnd, xIni, bfold

    bfold = -efold

    p = infParam%consts(3)
    mu = infParam%consts(2)
    
    xEnd = gdwi_x_endinf(p,mu)   
    
    if (display) write(*,*)'gdwi_initial_field: xend= ',xEnd

    xIni = gdwi_x_trajectory(bfold,xEnd,p,mu)

    gdwi_initial_field(:) = xIni*mu          

  end function gdwi_initial_field

  

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



  function ccsi1_initial_field(infParam,efold)
    use ccsi1sr, only : ccsi1_x_endinf,ccsi1_x_trajectory, ccsi1_numacc_xinimax
    implicit none
    real(kp), dimension(matterNum) :: ccsi1_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alpha
    real(kp) :: xEnd, xIni, bfold, xIniNumAcc

    bfold = -efold
    
    alpha = infParam%consts(2)
    
    if (alpha.lt.0._kp) stop 'ccsi1_initial_field: alpha<0!'
   
    xEnd = ccsi1_x_endinf(alpha)

    if (display) write(*,*)'ccsi1_initial_field: xend= ',xEnd

    xIni = ccsi1_x_trajectory(bfold,xEnd,alpha)
    xIniNumAcc = ccsi1_numacc_xinimax(alpha)

    if (xIni.gt.XiniNumAcc) then
       write(*,*)'ccsi1_initial_field: numerical accuracy reached for xIni'
    endif

    ccsi1_initial_field(:) = xIni * sqrt(3._kp/2._kp)

  end function ccsi1_initial_field



  function ccsi2_initial_field(infParam,efold)
    use ccsi2sr, only : ccsi2_x_trajectory, ccsi2_numacc_xendmin
    implicit none
    real(kp), dimension(matterNum) :: ccsi2_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alpha
    real(kp) :: xEnd, xEndMinNum, xIni, bfold
    real(kp), dimension(2) :: fieldStop

    bfold = -efold
    
    alpha = infParam%consts(2)

    fieldStop = field_stopinf(infParam)
    xEnd = fieldStop(1)*sqrt(2._kp/3._kp)

    if (alpha.lt.0._kp) stop 'ccsi2_initial_field: alpha<0!'

    xEndMinNum = ccsi2_numacc_xendmin(efold,alpha)

    if (xEnd.lt.xEndMinNum) then
       write(*,*)'xEnd= xEndMinNum=  ',xEnd,xEndMinNum
       write(*,*)'ccsi2_initial_field: xEnd is not reachable at that numerical accuracy!'
    endif
           
    if (display) write(*,*)'ccsi2_initial_field: xend= ',xEnd

    xIni = ccsi2_x_trajectory(bfold,xEnd,alpha)

    ccsi2_initial_field(:) = xIni * sqrt(3._kp/2._kp)

  end function ccsi2_initial_field


   
  function ccsi3_initial_field(infParam,efold)
    use ccsi3sr, only : ccsi3_x_endinf,ccsi3_x_trajectory
    use ccsi3sr, only : ccsi3_alphamin
    implicit none
    real(kp), dimension(matterNum) :: ccsi3_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alpha, alphaMin
    real(kp) :: xEnd, xIni, bfold

    bfold = -efold
    
    alpha = infParam%consts(2)
    alphaMin = ccsi3_alphamin(efold)
 
    if (alpha.gt.0._kp) stop 'ccsi3_initial_field: alpha>0!'
    
    if (alpha.lt.alphamin) then
       write(*,*)'efold= alphamin(efold)= alpha= ',efold,alphaMin, alpha
       stop 'alpha too small to get enough inflation!'
    endif
   
    xEnd = ccsi3_x_endinf(alpha)

    if (display) write(*,*)'ccsi3_initial_field: xend= ',xEnd

    xIni = ccsi3_x_trajectory(bfold,xEnd,alpha)

    ccsi3_initial_field(:) = xIni * sqrt(3._kp/2._kp)

  end function ccsi3_initial_field



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
    real(kp) :: alpha,mu

    bfold = -efold
    alpha = 1.5_kp*infParam%consts(3)
    mu = infParam%consts(6)

    xEnd = mssmi_x_endinf(mu)

    if (alpha.ne.1._kp) stop 'mssmi_initial_field: alpha><1!'

    if (display) write(*,*)'mssmi_initial_field: xend= ',xEnd

    xIni = mssmi_x_trajectory(bfold,xEnd,mu)

    mssmi_initial_field(:) = xIni*mu

  end function mssmi_initial_field
  


  function oripi_initial_field(infParam, efold)
    use oripisr, only : oripi_x_trajectory, oripi_x_endinf
    implicit none
    real(kp), dimension(matterNum) :: oripi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: bfold
    real(kp) :: xEnd, xIni
    real(kp) :: mu

    bfold = -efold
    mu = infParam%consts(6)
   
    xEnd = oripi_x_endinf(mu)

    if (display) write(*,*)'oripi_initial_field: xend= ',xEnd

    xIni = oripi_x_trajectory(bfold,xEnd,mu)

    oripi_initial_field(:) = xIni*mu

  end function oripi_initial_field



  function gmssmi_initial_field(infParam, efold)
    use gmssmisr, only : gmssmi_x_trajectory, gmssmi_x_endinf
    implicit none
    real(kp), dimension(matterNum) :: gmssmi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: bfold
    real(kp) :: xEnd, xIni
    real(kp) :: alpha, beta, mu

    bfold = -efold
    alpha = 1.5_kp*infParam%consts(3)
    beta = 5._kp*infParam%consts(4)
    mu = infParam%consts(6)

    xEnd = gmssmi_x_endinf(alpha,mu)

    if (alpha.ne.beta) stop 'gmssmi_initial_field: params incorrect!'

    if (display) write(*,*)'gmssmi_initial_field: xend= ',xEnd

    xIni = gmssmi_x_trajectory(bfold,xEnd,alpha,mu)

    gmssmi_initial_field(:) = xIni*mu

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



  function pai_initial_field(infParam, efold)
    use paisr, only : pai_x_trajectory, pai_x_endinf
    implicit none
    real(kp), dimension(matterNum) :: pai_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: bfold
    real(kp) :: xEnd, xIni
    real(kp) :: mu

    bfold = -efold
    mu = infParam%consts(2)
   
    xEnd = pai_x_endinf(mu)

    if (display) write(*,*)'pai_initial_field: xend= ',xEnd

    xIni = pai_x_trajectory(bfold,xEnd,mu)

    pai_initial_field(:) = xIni*mu

  end function pai_initial_field
  


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

    oi_initial_field(:) = xIni * mu

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

  

  function ssbi1_initial_field(infParam,efold)
    use ssbi1sr, only : ssbi1_x_endinf,ssbi1_x_trajectory
    use ssbi1sr, only : ssbi1_alphamin
    implicit none
    real(kp), dimension(matterNum) :: ssbi1_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alpha, beta, alphamin
    real(kp) :: xEnd, xIni, bfold

    bfold = -efold
    
    alpha = infParam%consts(2)
    beta = infParam%consts(3)

    alphamin = ssbi1_alphamin(beta)

    if (alpha.lt.alphamin) then
       write(*,*)'ssbi1_initial_field: alpha= alphamin= ',alpha,alphamin
       stop
    endif
       
    xEnd = ssbi1_x_endinf(alpha,beta)

    if (display) write(*,*)'ssbi1_initial_field: xend= ',xEnd

    xIni = ssbi1_x_trajectory(bfold,xEnd,alpha,beta)

    ssbi1_initial_field(:) = xIni

  end function ssbi1_initial_field



  function ssbi2_initial_field(infParam,efold)
    use ssbi2sr, only : ssbi2_x_endinf,ssbi2_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: ssbi2_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alpha, beta
    real(kp) :: xEnd, xIni, bfold

    bfold = -efold
    
    alpha = infParam%consts(2)
    beta = infParam%consts(3)
  
    xEnd = ssbi2_x_endinf(alpha,beta)

    if (display) write(*,*)'ssbi2_initial_field: xend= ',xEnd

    xIni = ssbi2_x_trajectory(bfold,xEnd,alpha,beta)

    ssbi2_initial_field(:) = xIni

  end function ssbi2_initial_field



  function ssbi3_initial_field(infParam,efold)
    use ssbi3sr, only : ssbi3_x_endinf,ssbi3_x_trajectory
    use ssbi3sr, only : ssbi3_alphamin
    implicit none
    real(kp), dimension(matterNum) :: ssbi3_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alpha, beta, alphamin
    real(kp) :: xEnd, xIni, bfold

    bfold = -efold
    
    alpha = infParam%consts(2)
    beta = infParam%consts(3)

    alphamin = ssbi3_alphamin(beta)

    if (alpha.lt.alphamin) then
       write(*,*)'ssbi3_initial_field: alpha= alphamin= ',alpha,alphamin
       stop
    endif
       
    xEnd = ssbi3_x_endinf(alpha,beta)

    if (display) write(*,*)'ssbi3_initial_field: xend= ',xEnd

    xIni = ssbi3_x_trajectory(bfold,xEnd,alpha,beta)

    ssbi3_initial_field(:) = xIni

  end function ssbi3_initial_field



  function ssbi4_initial_field(infParam,efold)
    use ssbi4sr, only : ssbi4_x_endinf,ssbi4_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: ssbi4_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alpha, beta
    real(kp) :: xEnd, xIni, bfold

    bfold = -efold
    
    alpha = infParam%consts(2)
    beta = infParam%consts(3)
  
    xEnd = ssbi4_x_endinf(alpha,beta)

    if (display) write(*,*)'ssbi4_initial_field: xend= ',xEnd

    xIni = ssbi4_x_trajectory(bfold,xEnd,alpha,beta)

    ssbi4_initial_field(:) = xIni

  end function ssbi4_initial_field



  function ssbi5_initial_field(infParam,efold)
    use ssbi5sr, only : ssbi5_x_endinf,ssbi5_x_trajectory
    use ssbi5sr, only : ssbi5_alphamax
    implicit none
    real(kp), dimension(matterNum) :: ssbi5_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alpha, beta, alphamax
    real(kp) :: xEnd, xIni, bfold

    bfold = -efold
    
    alpha = infParam%consts(2)
    beta = infParam%consts(3)

    alphamax = ssbi5_alphamax(beta)

    if (alpha.gt.alphamax) then
       write(*,*)'ssbi5_initial_field: alpha= alphamax= ',alpha,alphamax
       stop
    endif
       
    xEnd = ssbi5_x_endinf(alpha,beta)

    if (display) write(*,*)'ssbi5_initial_field: xend= ',xEnd

    xIni = ssbi5_x_trajectory(bfold,xEnd,alpha,beta)

    ssbi5_initial_field(:) = xIni

  end function ssbi5_initial_field



  function ssbi6_initial_field(infParam,efold)
    use ssbi6sr, only : ssbi6_x_endinf,ssbi6_x_trajectory
    use ssbi6sr, only : ssbi6_alphamax
    implicit none
    real(kp), dimension(matterNum) :: ssbi6_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alpha, beta, alphamax
    real(kp) :: xEnd, xIni, bfold

    bfold = -efold
    
    alpha = infParam%consts(2)
    beta = infParam%consts(3)

    alphamax = ssbi6_alphamax(beta)

    if (alpha.gt.alphamax) then
       write(*,*)'ssbi6_initial_field: alpha= alphamax= ',alpha,alphamax
       stop
    endif
       
    xEnd = ssbi6_x_endinf(alpha,beta)

    if (display) write(*,*)'ssbi6_initial_field: xend= ',xEnd

    xIni = ssbi6_x_trajectory(bfold,xEnd,alpha,beta)

    ssbi6_initial_field(:) = xIni

  end function ssbi6_initial_field



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



  function osti_initial_field(infParam,efold)
    use ostisr, only : osti_x_endinf,osti_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: osti_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: mu, p
    real(kp) :: xEnd, xIni, bfold

    bfold = -efold
    
    p = infParam%consts(2)
    mu = infParam%consts(3)
   
    if (p.ne.2._kp) then
       stop 'osti_initial_field: p >< 2!'
    endif

    xEnd = osti_x_endinf(mu)

    if (display) write(*,*)'osti_initial_field: xend= ',xEnd

    xIni = osti_x_trajectory(bfold,xEnd,mu)

    osti_initial_field(:) = xIni*mu

  end function osti_initial_field


 
  function wri_initial_field(infParam,efold)
    use wrisr, only : wri_x_endinf,wri_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: wri_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: mu, p, q
    real(kp) :: xEnd, xIni, bfold

    bfold = -efold
    
    p = infParam%consts(2)
    q = infParam%consts(4)
    mu = infParam%consts(3)
   
    if ((p.ne.0._kp).and.(q.ne.2._kp)) then
       write(*,*)'p= q= ',p,q
       stop 'wri_initial_field: wrong parameters value'
    endif

    xEnd = wri_x_endinf(mu)

    if (display) write(*,*)'wri_initial_field: xend= ',xEnd

    xIni = wri_x_trajectory(bfold,xEnd,mu)

    wri_initial_field(:) = xIni*mu

  end function wri_initial_field



  function imi_initial_field(infParam,efold)
    use imisr, only : imi_x_trajectory, imi_xendmin
    implicit none
    real(kp), dimension(matterNum) :: imi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) ::p, XIni, xEnd, bfold, xEndMin
    real(kp), dimension(2) :: fieldStop

    bfold = -efold

    p = infParam%consts(2)

    fieldStop = field_stopinf(infParam)
       
    xEnd = fieldStop(1)

    if (display) write(*,*)'imi_initial_field: xend= ',xEnd
   
    xEndMin = imi_xendmin(efold,p)
    if (xEnd.lt.xEndMin) then
       write(*,*)'xend= xendmin= ',xEnd,xEndMin
       stop 'imi_initial_field: xend too small!'
    endif

    xIni = imi_x_trajectory(bfold,xEnd,p)

    imi_initial_field(:) = xIni 

  end function imi_initial_field


  function ripi_initial_field(infParam, efold)
    use ripisr, only : ripi_x_trajectory, ripi_x_endinf
    implicit none
    real(kp), dimension(matterNum) :: ripi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: bfold
    real(kp) :: xEnd, xIni
    real(kp) :: alpha, mu

    bfold = -efold
    alpha = 2._kp*infParam%consts(4)
    mu = infParam%consts(6)
   
    if (alpha.ne.1._kp) then
       stop 'ripi_initial_field: improper params!'
    endif

    xEnd = ripi_x_endinf(mu)

    if (display) write(*,*)'ripi_initial_field: xend= ',xEnd

    xIni = ripi_x_trajectory(bfold,xEnd,mu)

    ripi_initial_field(:) = xIni*mu

  end function ripi_initial_field


 function gripi_initial_field(infParam, efold)
    use gripisr, only : gripi_x_trajectory, gripi_x_endinf
    implicit none
    real(kp), dimension(matterNum) :: gripi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: bfold
    real(kp) :: xEnd, xIni
    real(kp) :: alpha, mu

    bfold = -efold
    alpha = 2._kp*infParam%consts(4)
    mu = infParam%consts(6)
   
    if (3._kp/4._kp*infParam%consts(3).ne.alpha) then
       stop 'gripi_initial_field: improper params!'
    endif

    xEnd = gripi_x_endinf(alpha,mu)

    if (display) write(*,*)'gripi_initial_field: xend= ',xEnd

    xIni = gripi_x_trajectory(bfold,xEnd,alpha,mu)

    gripi_initial_field(:) = xIni*mu

  end function gripi_initial_field



  function bi_initial_field(infParam,efold)
    use bisr, only : bi_x_epsoneunity, bi_x_trajectory
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
    if (p.lt.0._kp) stop 'bi_initial_field: improper p<0!'

    fieldUv = field_thbound(infParam)
    fieldStop = field_stopinf(infParam)

    xUv = fieldUv(1)/mu
    xStrg = fieldStop(1)/mu

    xEps = bi_x_epsoneunity(p,mu)

!if phistrg occurs before slow-roll violation
    xEnd = max(xStrg,xEps)

    if (display) write(*,*)'bi_initial_field: xend= ',xEnd

    if (xEnd.gt.xUv) then
       write(*,*)'bi_initial_field: xEnd > XUv!'
       write(*,*)'xEnd= xUV= ',xEnd,xUV
       stop
    endif

    xIni = bi_x_trajectory(bfold,xEnd,p,mu)

    if (xIni.gt.xUv) then
       write(*,*)'bi_initial_field: xIni > XUv!'
       write(*,*)'xIni= xUV= ',xIni,xUV
       stop
    endif

    bi_initial_field(:) = xIni * mu
    
  end function bi_initial_field
  
   
  
  function kklti_initial_field(infParam,efold)
    use kkltisr, only : kklti_x_epsoneunity, kklti_x_trajectory
    implicit none
    real(kp), dimension(matterNum) :: kklti_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: p, mu, xEnd, xEps, xIni, xUv, xStrg, bfold
    real(kp), dimension(2) :: fieldUv, fieldStop

    bfold = -efold

    p = infParam%consts(2)
    mu = infParam%consts(3)

    if (mu.le.0._kp) stop 'kklti_initial_field: improper mu<0!'
    if (p.lt.2._kp) stop 'kklti_initial_field: improper p<0!'
    

    fieldUv = field_thbound(infParam)
    fieldStop = field_stopinf(infParam)

    xUv = fieldUv(1)/mu
    xStrg = fieldStop(1)/mu

    xEps = kklti_x_epsoneunity(p,mu)

!if phistrg occurs before slow-roll violation
    xEnd = max(xStrg,xEps)

    if (display) write(*,*)'kklti_initial_field: xend= ',xEnd

    if (xEnd.gt.xUv) then
       write(*,*)'kklti_initial_field: xEnd > XUv!'
       write(*,*)'xEnd= xUV= ',xEnd,xUV
       stop
    endif

    xIni = kklti_x_trajectory(bfold,xEnd,p,mu)

    if (xIni.gt.xUv) then
       write(*,*)'kklti_initial_field: xIni > XUv!'
       write(*,*)'xIni= xUV= ',xIni,xUV
       stop
    endif

    kklti_initial_field(:) = xIni * mu

  end function kklti_initial_field

  

  function nfi1_initial_field(infParam, efold)
    use nfi1sr, only : nfi1_x_trajectory, nfi1_x_endinf
    use nfi1sr, only : nfi1_check_params, nfi1_numacc_amin
    use nfi1sr, only : nfi1_amax
    implicit none
    real(kp), dimension(matterNum) :: nfi1_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: bfold
    real(kp) :: xEnd, xIni
    real(kp) :: a,b, amin,amax

    bfold = -efold
    a = infParam%consts(2)
    b = infParam%consts(3)

    if (.not.nfi1_check_params(a,b)) then
       stop 'nfi1_initial_field: nfi1 requires a>0, b>1'
    endif

    amin = nfi1_numacc_amin(b)
   
    if (a.lt.amin) then
       write(*,*)'nfi1_initial_field: a<amin= ',amin
       write(*,*)'potential values larger than huge'
    endif

    if (b.lt.2._kp) then
       amax = nfi1_amax(efold,b)
    
       if (a.gt.amax) then
          write(*,*)'nfi1_initial_field: a>amax= ',amax
          write(*,*)'not enough efolds for efoldWanted= ',efold
       endif
    endif


    xEnd = nfi1_x_endinf(a,b)

    if (display) write(*,*)'nfi1_initial_field: xend= ',xEnd

    xIni = nfi1_x_trajectory(bfold,xEnd,a,b)

    nfi1_initial_field(:) = xIni

  end function nfi1_initial_field



  function nfi2_initial_field(infParam, efold)
    use nfi2sr, only : nfi2_x_trajectory
    use nfi2sr, only : nfi2_check_params, nfi2_xendmax
    use nfi2sr, only : nfi2_numacc_xendmax
    implicit none
    real(kp), dimension(matterNum) :: nfi2_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold
    
    real(kp), dimension(2) :: fieldStop
    real(kp) :: bfold
    real(kp) :: xEnd, xIni, xEndMax
    real(kp) :: a,b

    bfold = -efold
    a = infParam%consts(2)
    b = infParam%consts(3)

    if (.not.nfi2_check_params(a,b)) then
       stop 'nfi2_initial_field: nfi2 requires a<0, b>1' 
    endif

    fieldStop = field_stopinf(infParam)

    xEnd = fieldStop(1)

    if (display) write(*,*)'nfi2_initial_field: xend= ',xEnd

    xEndMax = nfi2_xendmax(efold,a,b)

    if (xEnd.gt.xEndMax) then
       write(*,*)'xend= xendmax= efold= ',xEnd,xEndMax,efold
       stop 'nfi2_initial_field: xend too large'
    endif

    if (xEnd.gt.nfi2_numacc_xendmax(efold,a,b)) then
       write(*,*)'nfi2_initial_field: potential values larger than huge'
    endif

    xIni = nfi2_x_trajectory(bfold,xEnd,a,b)

    nfi2_initial_field(:) = xIni

  end function nfi2_initial_field



  function nfi3_initial_field(infParam, efold)
    use nfi3sr, only : nfi3_x_trajectory, nfi3_x_endinf
    use nfi3sr, only : nfi3_check_params, nfi3_numacc_absamax
    implicit none
    real(kp), dimension(matterNum) :: nfi3_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: bfold
    real(kp) :: xEnd, xIni
    real(kp) :: a,b, absamax

    bfold = -efold
    a = infParam%consts(2)
    b = infParam%consts(3)

    if (.not.nfi3_check_params(a,b)) then
       stop 'nfi3_initial_field: nfi3 requires a<0, 0<b<1 or a>0, b<0'
    endif

    absamax = nfi3_numacc_absamax(b)

   
    if (abs(a).gt.absamax) then
       write(*,*)'nfi3_initial_field: |a|>amax= ',absamax
       write(*,*)'potential values larger than huge'
    endif

    xEnd = nfi3_x_endinf(a,b)

    if (display) write(*,*)'nfi3_initial_field: xend= ',xEnd

    xIni = nfi3_x_trajectory(bfold,xEnd,a,b)

    nfi3_initial_field(:) = xIni

  end function nfi3_initial_field


  
  function nfi4_initial_field(infParam, efold)
    use nfi4sr, only : nfi4_x_trajectory
    use nfi4sr, only : nfi4_check_params, nfi4_xendmin
    use nfi4sr, only : nfi4_numacc_xendmax
    implicit none
    real(kp), dimension(matterNum) :: nfi4_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp), dimension(2) :: fieldStop

    real(kp) :: bfold
    real(kp) :: xEnd, xIni, xEndMin
    real(kp) :: a,b

    bfold = -efold
    a = infParam%consts(2)
    b = infParam%consts(3)

    if (.not.nfi4_check_params(a,b)) then
       stop 'nfi4_initial_field: nfi4 requires  a>0, 0<b<1 or a<0, b<0'
    endif

    fieldStop = field_stopinf(infParam)

    xEnd = fieldStop(1)

    if (display) write(*,*)'nfi4_initial_field: xend= ',xEnd

    xEndMin = nfi4_xendmin(efold,a,b)

    if (xEnd.lt.xEndMin) then
       write(*,*)'xend= xendmin= efold= ',xEnd,xEndMin,efold
       stop 'nfi4_initial_field: xend too small'
    endif

    if (xEnd.gt.nfi4_numacc_xendmax(a,b)) then
       write(*,*)'nfi4_initial_field: xend too large for numerical accuracy'
    endif

    xIni = nfi4_x_trajectory(bfold,xEnd,a,b)

    nfi4_initial_field(:) = xIni

  end function nfi4_initial_field


  function di_initial_field(infParam, efold)
    use disr, only : di_k2_trajectory, di_k2_epsoneunity, di_x
    implicit none
    real(kp), dimension(matterNum) :: di_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: k2End, k2Ini, f, lambda
    real(kp) :: bfold, xend

    bfold = -efold
    f = infParam%consts(1)
    lambda = infParam%consts(2)

    k2end = di_k2_epsoneunity(f,lambda)
    
    if (display) then
       xend = di_x(k2end)
       write(*,*)'di_initial_field: k2end= matterEnd= ',k2end,xend*lambda
    endif

    k2ini = di_k2_trajectory(bfold,k2end,f,lambda)

    di_initial_field(:) = di_x(k2ini)*lambda

  end function di_initial_field


  function ncli_initial_field(infParam, efold)
    use nclisr, only : ncli_x_trajectory, ncli_x_endinf, ncli_efold_primitive
    use nclisr, only : ncli_xinimax, ncli_phizeromin, ncli_x_epsoneunity
    
    real(kp), dimension(matterNum) :: ncli_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alpha, p, mu, muMin, efoldMax
    real(kp) :: xEnd, xIni, xIniMax, bfold

    bfold = -efold

    alpha = infParam%consts(2)
    mu = infParam%consts(3)
    p = infParam%consts(4)
    
    if (mu.lt.ncli_phizeromin(alpha,p)) then
       write(*,*)'mu= muMin= ',mu,muMin
       stop 'ncli_initial_field: mu too small'
    endif

    xEnd = ncli_x_endinf(alpha,mu,p)

    if (display) write(*,*)'ncli_initial_field: xend= ',xEnd
   
    xIniMax = ncli_xinimax(alpha,mu,p)

    efoldMAx = -ncli_efold_primitive(xend,alpha,mu,p) &
         + ncli_efold_primitive(xiniMax,alpha,mu,p)

    if (efold.gt.efoldMax) then
       write(*,*)'ncli_initial_field: efold > efoldMax in plateau'
       write(*,*)'efold= efoldMax= ',efold,efoldMax
    end if
       
    xIni = ncli_x_trajectory(bfold,xEnd,alpha,mu,p)

    ncli_initial_field(:) = xIni

  end function ncli_initial_field



  function vfmi_initial_field(infParam, efold)
    use vfmisr, only : vfmi_x_endinf, vfmi_x_trajectory

    
    real(kp), dimension(matterNum) :: vfmi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alpha, beta, bfold
    real(kp) :: xEnd, xIni

    bfold = -efold

    alpha = infParam%consts(2)
    beta = infParam%consts(3)
    
    xEnd = vfmi_x_endinf(alpha,beta)
    if (display) write(*,*)'vfmi_initial_field: xend= ',xEnd
            
    xIni = vfmi_x_trajectory(bfold,xEnd,alpha,beta)

    vfmi_initial_field(:) = xIni

  end function vfmi_initial_field



  function ahi_initial_field(infParam, efold)
    use ahisr, only : ahi_x_endinf, ahi_x_trajectory
   
    real(kp), dimension(matterNum) :: ahi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: mu, bfold
    real(kp) :: xEnd, xIni

    bfold = -efold

    mu = infParam%consts(2)
    
    xEnd = ahi_x_endinf(mu)
    if (display) write(*,*)'ahi_initial_field: xend= ', xEnd
            
    xIni = ahi_x_trajectory(bfold,xEnd,mu)

    ahi_initial_field(:) = xIni*mu

  end function ahi_initial_field



  function sbki_initial_field(infParam, efold)
    use sbkisr, only : sbki_x_trajectory, sbki_x_endinf, sbki_efold_primitive
    use sbkisr, only : sbki_xinimax, sbki_alphamin, sbki_alphamax
    
    real(kp), dimension(matterNum) :: sbki_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alpha, alphamin, alphamax, efoldMax
    real(kp) :: xEnd, xIni, xIniMax, bfold

    bfold = -efold

    alpha = infParam%consts(2)
    alphamin = sbki_alphamin()
    alphamax = sbki_alphamax()

    if ((alpha.lt.alphamin).or.(alpha.gt.alphamax)) then
       write(*,*)'alpha= ',alpha
       write(*,*)'alpha should be in the range: ',alphamin, alphamax
       stop 'sbki_initial_field: alpha out of range'
    endif

    xEnd = sbki_x_endinf(alpha)

    if (display) write(*,*)'sbki_initial_field: xend= ',xEnd
   
    xIniMax = sbki_xinimax(alpha)

    efoldMAx = -sbki_efold_primitive(xend,alpha) &
         + sbki_efold_primitive(xiniMax,alpha)

    if (efold.gt.efoldMax) then
       write(*,*)'sbki_initial_field: efold > efoldMax'
       write(*,*)'efold= efoldMax= ',efold,efoldMax
    end if
       
    xIni = sbki_x_trajectory(bfold,xEnd,alpha)

    sbki_initial_field(:) = xIni

  end function sbki_initial_field



  function fi_initial_field(infParam, efold)
    use fisr, only : fi_x_trajectory, fi_x_endinf, fi_efold_primitive
    use fisr, only : fi_efoldmax
    real(kp), dimension(matterNum) :: fi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: delta, n
    real(kp) :: xEnd, xIni, bfold, efoldMax

    bfold = -efold

    delta = infParam%consts(2)
    n = infParam%consts(3) - 1._kp
    
    xEnd = fi_x_endinf(delta,n)

    if (display) write(*,*)'fi_initial_field: xend= ',xEnd
      
    efoldMax = fi_efoldmax(delta,n)

    if (efold.gt.efoldMax) then
       write(*,*)'fi_initial_field: efold > efoldMax'
       write(*,*)'efold= efoldMax= ',efold,efoldMax
    end if
       
    xIni = fi_x_trajectory(bfold,xEnd,delta,n)

    fi_initial_field(:) = xIni

  end function fi_initial_field


  function sdi_initial_field(infParam, efold)
    use sdisr, only : sdi_x_trajectory
    use sdisr, only : sdi_phizeromin

    real(kp), dimension(matterNum) :: sdi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: mu, bfold, mumin
    real(kp) :: xEnd, xIni
    real(kp), dimension(2) :: fieldStop

    bfold = -efold

    mumin = sdi_phizeromin()
    mu = infParam%consts(2)

    if (mu.lt.mumin) then
       write(*,*)'mu= mumin= ',mu, mumin
       stop 'sbki_initial_field: mu too small!'
    endif

    fieldstop = field_stopinf(infParam)

    xEnd = fieldStop(1)/mu
    if (display) write(*,*)'sdi_initial_field: xend= ', xEnd

    xIni = sdi_x_trajectory(bfold,xEnd,mu)

    sdi_initial_field(:) = xIni*mu

  end function sdi_initial_field



  function saai_initial_field(infParam, efold)
    use saaisr, only : saai_x_endinf, saai_x_trajectory

    real(kp), dimension(matterNum) :: saai_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alpha, bfold
    real(kp) :: xEnd, xIni

    bfold = -efold

    alpha = infParam%consts(2)

    xEnd = saai_x_endinf(alpha)
    if (display) write(*,*)'saai_initial_field: xend= ', xEnd

    xIni = saai_x_trajectory(bfold,xend,alpha)

    saai_initial_field(:) = xIni

  end function saai_initial_field



  function sabi_initial_field(infParam, efold)
    use sabisr, only : sabi_x_endinf, sabi_x_trajectory

    real(kp), dimension(matterNum) :: sabi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alpha, p, bfold
    real(kp) :: xEnd, xIni

    bfold = -efold

    p = infParam%consts(2)/2._kp
    alpha = infParam%consts(3)**2/6._kp

    xEnd = sabi_x_endinf(alpha,p)
    if (display) write(*,*)'sabi_initial_field: xend= ', xEnd

    xIni = sabi_x_trajectory(bfold,xend,alpha,p)

    sabi_initial_field(:) = xIni

  end function sabi_initial_field



  function saci_initial_field(infParam, efold)
    use sacisr, only : saci_x_endinf, saci_x_trajectory

    real(kp), dimension(matterNum) :: saci_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alpha, p, bfold
    real(kp) :: xEnd, xIni

    bfold = -efold

    p = infParam%consts(2)/2._kp
    alpha = infParam%consts(3)**2/6._kp

    xEnd = saci_x_endinf(alpha,p)
    if (display) write(*,*)'saci_initial_field: xend= ', xEnd

    xIni = saci_x_trajectory(bfold,xend,alpha,p)

    saci_initial_field(:) = xIni

  end function saci_initial_field



  function hbi_initial_field(infParam, efold)
    use hbisr, only : hbi_x_endinf, hbi_x_trajectory

    real(kp), dimension(matterNum) :: hbi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: mu, p, bfold
    real(kp) :: xEnd, xIni

    bfold = -efold

    p = infParam%consts(2)
    mu = infParam%consts(3)

    xEnd = hbi_x_endinf(p,mu)

    if (display) write(*,*)'hbi_initial_field: xend= ', xEnd

    xIni = hbi_x_trajectory(bfold,xend,p,mu)

    hbi_initial_field(:) = xIni*mu

  end function hbi_initial_field


  function shi_initial_field(infParam, efold)
    use shisr, only : shi_x_endinf, shi_x_trajectory

    real(kp), dimension(matterNum) :: shi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: mu, alpha, bfold
    real(kp) :: xEnd, xIni

    bfold = -efold

    alpha = infParam%consts(2)
    mu = infParam%consts(3)

    xEnd = shi_x_endinf(alpha,mu)

    if (display) write(*,*)'shi_initial_field: xend= ', xEnd

    xIni = shi_x_trajectory(bfold,xend,alpha,mu)

    shi_initial_field(:) = xIni*mu

  end function shi_initial_field



  function rcipi_initial_field(infParam, efold)
    use rcipisr, only : rcipi_x_endinf, rcipi_x_trajectory
    use rcipisr, only : rcipi_efoldmax

    real(kp), dimension(matterNum) :: rcipi_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: alpha, beta, p, bfold
    real(kp) :: xEnd, xIni, efoldMax

    bfold = -efold

    p = infParam%consts(2)
    alpha = infParam%consts(3)
    beta = infParam%consts(4)

    xEnd = rcipi_x_endinf(p,alpha,beta)

    if (display) write(*,*)'rcipi_initial_field: xend= ', xEnd

    efoldMax = rcipi_efoldmax(p,alpha,beta)

    if (efold.gt.efoldMax) then
       write(*,*)'rcipi_initial_field: efold > efoldMax!'
       write(*,*)'efold= efoldMax= ',efold,efoldMax
    end if

    xIni = rcipi_x_trajectory(bfold,xend,p,alpha,beta)

    rcipi_initial_field(:) = xIni

  end function rcipi_initial_field


 function saii1_initial_field(infParam, efold)
    use saii1sr, only : saii1_x_endinf, saii1_x_trajectory
    use saii1sr, only : saii1_numacc_mumin
    
    real(kp), dimension(matterNum) :: saii1_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: mu, alpha, bfold
    real(kp) :: xEnd, xIni

    bfold = -efold

    alpha = infParam%consts(2)
    mu = infParam%consts(3)

    if (mu.lt.saii1_numacc_mumin(efold,alpha)) then
       write(*,*)'saii1_initial_field:'
       write(*,*)'mu is numerically too small!'
    endif
    
    xEnd = saii1_x_endinf(alpha,mu)

    if (display) write(*,*)'saii1_initial_field: xend= ', xEnd

    xIni = saii1_x_trajectory(bfold,xend,alpha,mu)

    saii1_initial_field(:) = xIni*mu

  end function saii1_initial_field

  
  function saii2_initial_field(infParam, efold)
    use saii2sr, only : saii2_x_endinf, saii2_x_trajectory
    use saii2sr, only : saii2_numacc_mumin
    real(kp), dimension(matterNum) :: saii2_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: mu, alpha, bfold
    real(kp) :: xEnd, xIni

    bfold = -efold

    alpha = infParam%consts(2)
    mu = infParam%consts(3)

    if (mu.lt.saii2_numacc_mumin(efold,alpha)) then
       write(*,*)'saii2_initial_field:'
       write(*,*)'mu is numerically too small!'
    endif
    
    xEnd = saii2_x_endinf(alpha,mu)

    if (display) write(*,*)'saii2_initial_field: xend= ', xEnd

    xIni = saii2_x_trajectory(bfold,xend,alpha,mu)

    saii2_initial_field(:) = xIni*mu

  end function saii2_initial_field


 function saiii1_initial_field(infParam, efold)
   use saiii1sr, only : saiii1_x_endinf, saiii1_x_trajectory
   use saiii1sr, only : saiii1_numacc_mumin

    real(kp), dimension(matterNum) :: saiii1_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: mu, alpha, beta, bfold
    real(kp) :: xEnd, xIni

    bfold = -efold

    alpha = infParam%consts(2)
    beta = infParam%consts(3)
    mu = infParam%consts(4)

    if (mu.lt.saiii1_numacc_mumin(efold,alpha,beta)) then
       write(*,*)'saiii1_initial_field:'
       write(*,*)'mu is numerically too small!'
    endif
    
    xEnd = saiii1_x_endinf(alpha,beta,mu)

    if (display) write(*,*)'saiii1_initial_field: xend= ', xEnd

    xIni = saiii1_x_trajectory(bfold,xend,alpha,beta,mu)

    saiii1_initial_field(:) = xIni*mu

  end function saiii1_initial_field



  function saiii2_initial_field(infParam, efold)
   use saiii2sr, only : saiii2_x_endinf, saiii2_x_trajectory
   use saiii2sr, only : saiii2_numacc_mumin

    real(kp), dimension(matterNum) :: saiii2_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: mu, alpha, beta, bfold
    real(kp) :: xEnd, xIni

    bfold = -efold

    alpha = infParam%consts(2)
    beta = infParam%consts(3)
    mu = infParam%consts(4)

    if (mu.lt.saiii2_numacc_mumin(efold,alpha,beta)) then
       write(*,*)'saiii2_initial_field:'
       write(*,*)'mu is numerically too small!'
    endif
    
    xEnd = saiii2_x_endinf(alpha,beta,mu)

    if (display) write(*,*)'saiii2_initial_field: xend= ', xEnd

    xIni = saiii2_x_trajectory(bfold,xend,alpha,beta,mu)

    saiii2_initial_field(:) = xIni*mu

  end function saiii2_initial_field


  function saiii3_initial_field(infParam, efold)
    use saiii3sr, only : saiii3_x_endinf, saiii3_x_trajectory

    real(kp), dimension(matterNum) :: saiii3_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: mu, alpha, beta, bfold
    real(kp) :: xEnd, xIni

    bfold = -efold

    alpha = infParam%consts(2)
    beta = infParam%consts(3)
    mu = infParam%consts(4)
    
    xEnd = saiii3_x_endinf(alpha,beta,mu)

    if (display) write(*,*)'saiii3_initial_field: xend= ', xEnd

    xIni = saiii3_x_trajectory(bfold,xend,alpha,beta,mu)

    saiii3_initial_field(:) = xIni*mu

  end function saiii3_initial_field


  function  dei_initial_field(infParam, efold)
    use deisr, only : dei_x_endinf, dei_x_trajectory

    real(kp), dimension(matterNum) :: dei_initial_field
    type(infbgparam), intent(in) :: infParam
    real(kp), intent(in) :: efold

    real(kp) :: mu , beta, bfold
    real(kp) :: xEnd, xIni

    bfold = -efold

    beta = infParam%consts(2)
    mu = infParam%consts(3)
    
    xEnd = dei_x_endinf(beta,mu)

    if (display) write(*,*)'dei_initial_field: xend= ', xEnd

    xIni = dei_x_trajectory(bfold,xend,beta,mu)

    dei_initial_field(:) = xIni*mu

  end function dei_initial_field
  
  
#endif

#ifdef TWOFIELDS
  function fterm_initial_field(infParam,efold)    
    use fieldprec, only : pi
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
