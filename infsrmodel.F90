!provides the slow-roll initial field values such that N e-folds of
!inflation occurs.

module infsrmodel
  use infprec, only : kp
  use infbgmodel, only : matterNum
  implicit none

  private

  logical, parameter :: display = .true.
  
  integer, parameter :: efoldBound = 110._kp


  public field_stopinf, field_thbound
  public slowroll_initial_matter_lf, slowroll_initial_matter_sf
  public slowroll_initial_matter_hy, slowroll_initial_matter_rm
  public slowroll_initial_matter_kksf, slowroll_initial_matter_kklt
  public slowroll_initial_matter_mix

  public sr_efold_sf, sr_endinf_sf, sr_iniinf_sf
  public sr_efold_kklt, sr_endinf_kklt, sr_iniinf_kklt
  public sr_efold_rm, sr_iniinf_rm
  public sr_efold_mix, sr_endinf_mix, sr_iniinf_mix



 



  
contains

  

  function field_stopinf(infParam)
!return the field stop value in (1) and if it is a maximum (+1) or
!minimum (-1) in (2)
    use infbgmodel, only : infbgparam, matterParamNum
    implicit none
    type(infbgparam), intent(in) :: infParam    
    real(kp), dimension(2) :: field_stopinf

    
    select case (infParam%name)

    case ('largef')
       field_stopinf(1) = infParam%consts(matterParamNum)
       field_stopinf(2) = -1._kp
   
    case ('smallf')
       field_stopinf(1) = infParam%consts(matterParamNum)*infParam%consts(3)
       field_stopinf(2) = 1._kp

    case ('hybrid')
       field_stopinf(1) = infParam%consts(matterParamNum) &
            * slowroll_stopmax_matter_hy(infParam)
       field_stopinf(2) = -1._kp
       
    case ('runmas')

       field_stopinf(1) = infParam%consts(matterParamNum)
!for nu>0
       if (infParam%consts(matterParamNum).lt.infParam%consts(3)) then
          field_stopinf(2) = -1._kp
       else
          field_stopinf(2) = +1._kp
       endif
!reverse if nu<0
       if (infParam%consts(4).lt.0._kp) then
          field_stopinf(2) = -field_stopinf(2)
       endif

    case ('kklmmt')
       field_stopinf(1) = infParam%consts(matterParamNum)*infParam%consts(3)
       field_stopinf(2) = -1._kp
       print *, 'fieldStop = infParam(matterParamNum)'

    case ('mixinf')
       field_stopinf(1) = infParam%consts(matterParamNum)
       field_stopinf(2) = -1._kp


    end select

    

  end function field_stopinf



  function field_thbound(infParam)
!returns a theoretical bound on the allowed field values if any. May be
!from the stochastic regime or uv limit in brane setup
    use infbgmodel, only : infbgparam, matterParamNum
    implicit none
    type(infbgparam), intent(in) :: infParam    
    real(kp), dimension(2) :: field_thbound

    
    select case (infParam%name)

    case ('kklmmt')
       field_thbound(1) = infParam%consts(matterParamNum-1)
       field_thbound(2) = +1._kp
       print *, 'fieldUv = infParam%consts(matterParamNum-1)'

    case default
       stop 'no theoretical field bound implemented for this model!'

    end select

  end function field_thbound




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!large field models
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  function slowroll_initial_matter_lf(infParam,efoldWanted)   
    use infbgmodel, only : infbgparam
    implicit none
    type(infbgparam), intent(in) :: infParam
    real(kp), optional, intent(in) :: efoldWanted
    real(kp), dimension(matterNum) :: slowroll_initial_matter_lf
    
    real(kp), parameter :: efoldDefault = efoldBound
    real(kp) :: matterEnd, matterIni
    real(kp) :: p,mu,efold


    if (present(efoldWanted)) then
       efold = efoldWanted
    else
       efold = efoldDefault
    endif
  
    mu = infParam%consts(3)
    if (mu.ne.0._kp) stop 'slowroll_initial_matter_lf: improper parameters'
       

    p = infParam%consts(2)
   
    matterEnd = p/sqrt(2._kp)
   
    matterIni = sqrt(p**2/2._kp + 2._kp*p*efold)

    if (display) then
       write(*,*)
       write(*,*)'slowroll_initial_matter_lf: (*kappa)'
       write(*,*)'matterEnd = ',matterEnd
       write(*,*)'matterIni = ',matterIni
       write(*,*)
    endif

    slowroll_initial_matter_lf = matterIni

  end function slowroll_initial_matter_lf




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!small field models
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  function slowroll_initial_matter_sf(infParam,efoldWanted)
    use infprec, only : transfert,tolkp
    use inftools, only : zbrent    
    use infbgmodel, only : infbgparam
    implicit none
    type(infbgparam), intent(in) :: infParam
    real(kp), optional, intent(in) :: efoldWanted
    real(kp), dimension(matterNum) :: slowroll_initial_matter_sf

    type(transfert) :: sfData
    real(kp), parameter :: efoldDefault = efoldBound
    real(kp), parameter :: tolFind = tolkp
    
    
    real(kp) :: matterOverMuEnd, matterOverMuIni
    real(kp) :: mini,maxi
    real(kp) :: p, mu, efold

  
    p = infParam%consts(2)    
    mu = infParam%consts(3)

    if (infParam%consts(3).le.0._kp) then
       stop 'slowroll_initial_matter_sf: improper parameters'
    endif
   

    if (p.lt.2._kp) then
       write(*,*) 'slowroll_initial_matter_sf: p = ',p
       stop
    endif

     if (present(efoldWanted)) then
       efold = efoldWanted
    else
       efold = efoldDefault
    endif

!find the end of inflation
      
    mini = epsilon(1._kp)
    maxi = 1._kp + epsilon(1._kp)

    sfData%real1 = p
    sfData%real2 = mu

    matterOverMuEnd = zbrent(sr_endinf_sf,mini,maxi,tolFind,sfData)
   


!find the initial field values efolds before

    mini = epsilon(1._kp)
    maxi = matterOverMuEnd

    sfData%real1 = p
    sfData%real2 = 2._kp*p*efold/mu**2 + sr_efold_sf(matterOverMuEnd,p)
    

    matterOverMuIni = zbrent(sr_iniinf_sf,mini,maxi,tolFind,sfData)
    
    slowroll_initial_matter_sf = matterOverMuIni*mu

    if (display) then
       write(*,*)
       write(*,*)'slowroll_initial_matter_sf: (*kappa) (/mu)'
       write(*,*)'matterEnd = ',matterOverMuEnd*mu,matterOverMuEnd
       write(*,*)'matterIni = ',matterOverMuIni*mu,matterOverMuIni
       write(*,*)
    endif


  end function slowroll_initial_matter_sf




  function sr_endinf_sf(x,sfData)
    use infprec, only : transfert    
    implicit none
    real(kp), intent(in) :: x    
    type(transfert), optional, intent(inout) :: sfData
    real(kp) :: sr_endinf_sf
    real(kp) :: p,mu

    p=sfData%real1
    mu=sfData%real2

    sr_endinf_sf = x**(p-1._kp) + sqrt(2._kp)*mu/abs(p) * (x**p - 1._kp)

  end function sr_endinf_sf
 



  function sr_iniinf_sf(x,sfData)
    use infprec, only : transfert
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: sfData
    real(kp) :: sr_iniinf_sf
    real(kp) :: p

    p=sfData%real1

    sr_iniinf_sf = sr_efold_sf(x,p) - sfData%real2
   
  end function sr_iniinf_sf




  function sr_efold_sf(x,p)
    implicit none
    real(kp), intent(in) :: x,p
    real(kp) :: sr_efold_sf
    
    if (p == 2._kp) then
       sr_efold_sf = x**2 - 2._kp * log(x)
    else
       sr_efold_sf = x**2 + 2._kp/(p-2._kp) * x**(2._kp-p)
    endif
    
  end function sr_efold_sf



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!kklmmt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!some function call the small field ones since kklmmt potential looks
!like the small field ones with p->-p


  function slowroll_initial_matter_kksf(infParam,efoldWanted)
   use infprec, only : transfert,tolkp
    use inftools, only : zbrent    
    use infbgmodel, only : infbgparam
    implicit none
    type(infbgparam), intent(in) :: infParam
    real(kp), optional, intent(in) :: efoldWanted
    real(kp), dimension(matterNum) :: slowroll_initial_matter_kksf

    type(transfert) :: sfData
    real(kp), parameter :: efoldDefault = efoldBound
    real(kp), parameter :: tolFind = tolkp
    
    
    real(kp) :: matterOverMuEnd, matterOverMuIni
    real(kp) :: mini,maxi
    real(kp) :: p, mu, efold

  
    p = infParam%consts(2)    
    mu = infParam%consts(3)

    if (infParam%consts(3).le.0._kp) then
       stop 'slowroll_initial_matter_kksf: improper parameters'
    endif
   

    if (p.lt.0._kp) then
       write(*,*) 'slowroll_initial_matter_kksf: p = ',p
       stop
    endif

     if (present(efoldWanted)) then
       efold = efoldWanted
    else
       efold = efoldDefault
    endif

!find the end of inflation
      
    mini = 1._kp + epsilon(1._kp)
    maxi = 1._kp/epsilon(1._kp)

    sfData%real1 = -p
    sfData%real2 = mu

    matterOverMuEnd = zbrent(sr_endinf_sf,mini,maxi,tolFind,sfData)
   


!find the initial field values efolds before

    mini = matterOverMuEnd
    maxi = 1._kp/epsilon(1._kp)

    sfData%real1 = -p
    sfData%real2 = 2._kp*(-p)*efold/mu**2 + sr_efold_sf(matterOverMuEnd,-p)
    

    matterOverMuIni = zbrent(sr_iniinf_sf,mini,maxi,tolFind,sfData)
    
    slowroll_initial_matter_kksf = matterOverMuIni*mu

    if (display) then
       write(*,*)
       write(*,*)'slowroll_initial_matter_kksf: (*kappa) (/mu)'
       write(*,*)'matterEnd = ',matterOverMuEnd*mu,matterOverMuEnd
       write(*,*)'matterIni = ',matterOverMuIni*mu,matterOverMuIni
       write(*,*)
    endif


  end function slowroll_initial_matter_kksf






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!hybrid models
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




  function slowroll_initial_matter_hy(infParam,efoldWanted)
    use infprec, only : transfert, tolkp
    use inftools, only : zbrent    
    use infbgmodel, only : infbgparam, matterParamNum
    implicit none
    type(infbgparam), intent(in) :: infParam
    real(kp), optional, intent(in) :: efoldWanted
    real(kp), dimension(matterNum) :: slowroll_initial_matter_hy

    type(transfert) :: hyData
    real(kp), parameter :: efoldDefault = efoldBound
    real(kp), parameter :: tolFind = tolkp
    
    real(kp), dimension(2) :: fieldStop
    real(kp) :: matterOverMuMax, matterOverMuIni, matterOverMuStop
    real(kp) :: mini,maxi
    real(kp) :: p, mu
    real(kp) :: muConnex, efold

  
    p = infParam%consts(2)
    mu = infParam%consts(3)

!sanity checks
    if (infParam%consts(3).le.0._kp) then
       stop 'slowroll_initial_matter_hy: improper consts(3)'
    endif

    if (p.lt.2._kp) then
       write(*,*) 'slowroll_initial_matter_hy: p = ',p
       stop
    endif

    if (infParam%consts(matterParamNum).gt.1._kp) then
       stop 'slowroll_initial_matter_hy: improper consts(5)'
    endif


    if (present(efoldWanted)) then
       efold = efoldWanted
    else
       efold = efoldDefault
    endif


!determines the maximum allowed field value to stop hybrid inflation
!(see the called function)
!    matterOverMuStopMax = slowroll_stopmax_matter_hy(infParam) / mu

!for hybrid inflation, it assumes that consts(matparamnum) is in unit of
!mattertopMax
    fieldStop = field_stopinf(infParam)
    matterOverMuStop = fieldStop(1) / mu
  

!upper field value bound
    matterOverMuMax = slowroll_startmax_matter_hy(infParam)


!find matterIni otherwise to get the right number of efolds

    mini = matterOverMuStop
    maxi = matterOverMuMax

    hyData%real1 = p
    hyData%real2 = 2._kp*p*(efold)/mu**2 + sr_efold_hy(matterOverMuStop,p)

    matterOverMuIni =  zbrent(sr_iniinf_hy,mini,maxi,tolFind,hyData)

    slowroll_initial_matter_hy = matterOverMuIni * mu

    if (display) then
       write(*,*)
       write(*,*)'slowroll_initial_matter_hy: (*kappa) (/mu)'       
       write(*,*)'matterStop = ',matterOverMuStop*mu, matterOverMuStop  
       write(*,*)'matterIni  = ',matterOverMuIni*mu, matterOverMuIni
       write(*,*)
    endif


  end function slowroll_initial_matter_hy





  function slowroll_stopmax_matter_hy(infParam)
    use infprec, only : transfert,tolkp
    use inftools, only : zbrent    
    use infbgmodel, only : infbgparam
    implicit none
    type(infbgparam), intent(in) :: infParam
    real(kp) :: slowroll_stopmax_matter_hy
    
    real(kp), parameter :: efoldHybrid = efoldBound

    real(kp), parameter :: tolFind = tolkp
    type(transfert) :: hyData
    real(kp) :: p,mu, mini,maxi
    real(kp) :: matterOverMuMax, matterOverMuStopMax

    p = infParam%consts(2)
    mu = infParam%consts(3)

!sanity checks
    if (infParam%consts(3).le.0._kp) then
       stop 'slowroll_stopmax_matter_hy: improper consts(3)'
    endif

    if (p.lt.2._kp) then
       write(*,*) 'slowroll_stopmax_matter_hy: p = ',p
       stop
    endif

    matterOverMuMax = slowroll_startmax_matter_hy(infParam)/mu

!MatterStopMax is the zero of the sr evolution equation with N --> -N
!and set to efoldHybrid.

    mini = epsilon(1._kp)
    maxi = matterOverMuMax

    hyData%real1 = p
    hyData%real2 = -2._kp*p*(efoldHybrid)/mu**2 + sr_efold_hy(matterOverMuMax,p) 
   
    matterOverMuStopMax = zbrent(sr_iniinf_hy,mini,maxi,tolFind,hyData)

!return the field value, in unit of kappa
    slowroll_stopmax_matter_hy = matterOverMuStopMax * mu

    if (display) then
       write(*,*)
       write(*,*)'slowroll_matter_stopmax_hy: (*kappa) (/mu)'
       write(*,*)'efold definition = ',efoldHybrid
       write(*,*)'matterStopMax = ',matterOverMuStopMax*mu, matterOverMuStopMax      
       write(*,*)
    endif        

  end function slowroll_stopmax_matter_hy






  function slowroll_startmax_matter_hy(infParam)
    use infprec, only : transfert,tolkp
    use inftools, only : zbrent
    use infbgmodel, only : infbgparam
    implicit none
    type(infbgparam), intent(in) :: infParam
    real(kp) :: slowroll_startmax_matter_hy

    real(kp), parameter :: tolFind = tolkp
    type(transfert) :: hyData
    real(kp) :: p,mu,mini,maxi, muConnex
    real(kp) :: matterOverMuTrans, matterOverMuMax

    p = infParam%consts(2)
    mu = infParam%consts(3)

!sanity checks
    if (infParam%consts(3).le.0._kp) then
       stop 'slowroll_startmax_matter_hy: improper consts(3)'
    endif

    if (p.lt.2._kp) then
       write(*,*) 'slowroll_startmax_matter_hy: p = ',p
       stop
    endif


!this is the mu above which eps1 = 1 has no solution and for which
!inflationary domains in field space are simply connected
    muConnex = p/sqrt(8._kp)

!this is the matter/mu for which eps1 is maximum. eps1(matter) is a
!increasing function wrt the field under this value, and decreasing
!above. So hybrid like behaviour only appears for matter <
!matterTrans.  Otherwise, this is a mixture between large field like
!and hybrid inflation during which the sign of eps2 changes
    matterOverMuTrans = (p - 1._kp)**(1._kp/p)

    
!We are looking to the maximum allowed value of the field to stop
!inflation and to get at least "efoldHybrid" efolds of inflation. Here
!we discard the case where inflation may start for matterIni >
!matterTrans since it would be large field like inflation. So the
!ultimate upper limit for matterIni, defined as matterMax, is the
!min(matterTrans,matterOne) where eps1(matterOne) = 1. MatterOne is
!seeked in [0,matterTrans]: this selects only the lower root, the
!other corresponding to the end of a large field like inflation)
    mini = epsilon(1._kp)
    maxi = matterOverMuTrans

    hyData%real1 = p
    hyData%real2 = mu

    if (mu.le.muConnex) then
       matterOverMuMax = zbrent(sr_endinf_hy,mini,maxi,tolFind,hyData)
    else
       matterOverMuMax = matterOverMuTrans
    endif

    slowroll_startmax_matter_hy = matterOverMuMax * mu
    
    if (display) then
       write(*,*)'slowroll_startmax_matter_hy: (*kappa) (/mu)'
       write(*,*)'matterMax = ',matterOverMuMax*mu,matterOverMuMax
       if (mu.lt.muConnex) then
          write(*,*) '<--- due to epsilon1 > 1 above'
       else
          write(*,*) '<--- due to epsilon2 > 0 above'
       endif
    endif

  end function slowroll_startmax_matter_hy



  function sr_endinf_hy(x,hyData)
    use infprec, only : transfert
    implicit none
    real(kp), intent(in) :: x    
    type(transfert), optional, intent(inout) :: hyData
    real(kp) :: sr_endinf_hy
    real(kp) :: p,mu

    p=hyData%real1
    mu=hyData%real2

    sr_endinf_hy = x**(p-1._kp) - sqrt(2._kp)*mu/p * (x**p + 1._kp)

  end function sr_endinf_hy




  function sr_iniinf_hy(x,hyData)
    use infprec, only : transfert
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: hyData
    real(kp) :: sr_iniinf_hy
    real(kp) :: p

    p=hyData%real1
  
    sr_iniinf_hy = sr_efold_hy(x,p) - hyData%real2
   
  end function sr_iniinf_hy




  function sr_efold_hy(x,p)
    implicit none
    real(kp), intent(in) :: x,p
    real(kp) :: sr_efold_hy

    if (p == 2._kp) then
       sr_efold_hy = x**2 + 2._kp * log(x)
    else
       sr_efold_hy = x**2 - 2._kp/(p-2._kp) * x**(2._kp-p)
    endif


  end function sr_efold_hy



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!running mass models
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  function slowroll_initial_matter_rm(infParam,efoldWanted)
    use infprec, only : transfert, tolkp
    use inftools, only : zbrent    
    use infbgmodel, only : infbgparam
    implicit none
    type(infbgparam), intent(in) :: infParam
    real(kp), optional, intent(in) :: efoldWanted
    real(kp), dimension(matterNum) :: slowroll_initial_matter_rm

    type(transfert) :: rmData
    real(kp), parameter :: efoldDefault = efoldBound
    real(kp), parameter :: tolFind = tolkp

    real(kp) :: efold
    real(kp) :: mini,maxi,randnum
    real(kp) :: p,mu,nu,lambda
    real(kp) :: matterOverMuStop,matterStop,matterEnd,matterIni
    real(kp) :: matterZero
    real(kp), dimension(2) :: fieldStop

!should match with the definition of consts!
    p = infParam%consts(2)
    mu = infParam%consts(3)
    nu = infParam%consts(4)

!sanity checks
    if (infParam%consts(3).le.0._kp) then
       stop 'slowroll_initial_matter_rm: improper consts(3)'
    endif

    if (p.lt.2._kp) then
       write(*,*) 'slowroll_initial_matter_rm: p = ',p
       stop
    endif

    if ((abs(infParam%consts(4)).gt.0.5_kp).or.(infParam%consts(4).eq.0._kp)) then
       write(*,*)'slowroll_initial_matter_rm: improper consts(4) = ',infParam%consts(4)
       read(*,*)
    endif
    
    lambda = nu * mu**p

    if (present(efoldWanted)) then
       efold = efoldWanted
    else
       efold = efoldDefault
    endif



!the potential is
! U = M^4 { 1 + nu*[1/p - ln(matter/mu)]*matter^p}


!matterstop required
    fieldStop = field_stopinf(infParam)
    matterStop = fieldStop(1)


    
!for nu>0 there is matterEnd such as eps1(matterEnd) = 1. So the end
!of inflation in that case is either given by matterEnd or by
!matterStop. Note however that eps2 is usually big in that case which
!makes this calculation useless DAMNED!!. The eps1=1 eq has
!one solution only in [mu,matterZero] where U(matterZero) = 0
    if ((nu.gt.0._kp).and.(matterStop.gt.mu)) then
       mini = mu + epsilon(1._kp)
       maxi = 1._kp/epsilon(1._kp)

       rmData%real1 = p
       rmData%real2 = mu
       rmData%real3 = nu              

       matterZero = zbrent(rm_potential,mini,maxi,tolFind,rmData)

       mini = mu
       maxi = matterZero

       rmData%real1 = p
       rmData%real2 = mu
       rmData%real3 = nu

       matterEnd = zbrent(sr_endinf_rm,mini,maxi,tolFind,rmData)
    
       if (display) then
          if (matterEnd.lt.matterStop) then
             write(*,*)'slowroll_initial_matter_rm: epsilon1 stops inflation in SR approx'
          endif
       endif

!useless       matterStop = min(matterStop,matterEnd)
    endif


!find matterIni according to matterStop to get the right number of
!efolds. The choice between the 4 models is done according to the
!value of matterStop. matterIni>1 has been allowed for RM4, with still
!matterStop < 1
   
    if (matterStop.eq.mu) then       
       write(*,*)'slowroll_initial_matter_rm: matterStop/mu = ',matterStop/mu
       stop
    endif

    if (matterStop.lt.mu) then
       mini = epsilon(1._kp)
       maxi = mu - epsilon(1._kp)
    elseif (matterStop.gt.mu) then
       mini = mu + epsilon(1._kp)
       if (nu.gt.0._kp) then
          maxi = 1._kp
       else
          maxi = 1._kp/epsilon(1._kp)
       endif
    else
       stop 'slowroll_initial_matter_rm: error'
    endif

     matterOverMuStop = matterStop/mu
    
!    print *,'mini maxi',mini,maxi

    rmData%real1 = p
    rmData%real2 = mu
    rmData%real3 = nu
    rmData%xend = sr_efold_rm(matterOverMuStop,p,lambda) + 2._kp*p*efold/mu**2

    matterIni =  zbrent(sr_iniinf_rm,mini,maxi,tolFind,rmData)

    slowroll_initial_matter_rm = matterIni

    if (display) then
       write(*,*)
       write(*,*)'slowroll_initial_matter_rm: (*kappa)'       
       write(*,*)'matterZero = ',matterZero
       write(*,*)'matterEnd  = ',matterEnd
       write(*,*)'matterStop = ',matterStop  
       write(*,*)'matterIni  = ',matterIni
       write(*,*)
    endif


  end function slowroll_initial_matter_rm




  function sr_efold_rm(x,p,l)
    use specialinf, only : dp, dei
    implicit none
    real(kp), intent(in) :: x,p,l
    real(kp) :: sr_efold_rm
    real(dp) :: argei1,argei2

!l=nu*mu^p
    
    if (p == 2._kp) then
       argei1 = 2._dp*log(x)
       sr_efold_rm = x**2 - (2._kp/l)*log(abs(log(x))) - dei(argei1)
    else
       argei1 = (2._dp-p)*log(x)
       argei2 = 2._dp*log(x)
       sr_efold_rm = x**2 - (2._kp/l)*dei(argei1) - (2._kp/p)*dei(argei2)       
    endif
    
!    print *,'arg ei',argei1,dei(argei1)

  end function sr_efold_rm




  function sr_endinf_rm(matter,rmData)
!vanishes for eps1=1 in the sr approx
    use infprec, only : transfert    
    implicit none
    real(kp), intent(in) :: matter    
    type(transfert), optional, intent(inout) :: rmData
    real(kp) :: sr_endinf_rm
    real(kp) :: p,mu,nu,x

    p=rmData%real1
    mu=rmData%real2
    nu=rmData%real3

    x = matter/mu
   
    sr_endinf_rm = 1._kp + nu*(1._kp/p - log(x))*matter**p &
         - (nu*p/sqrt(2._kp))*log(x)*matter**(p-1._kp)

  end function sr_endinf_rm
 




  function sr_iniinf_rm(matter,rmData)
    use infprec, only : transfert    
    implicit none
    real(kp), intent(in) :: matter    
    type(transfert), optional, intent(inout) :: rmData
    real(kp) :: sr_iniinf_rm
    real(kp) :: p,mu,nu,x,l

    p=rmData%real1
    mu=rmData%real2
    nu=rmData%real3

    x = matter/mu
    l = nu * mu**p

    sr_iniinf_rm = sr_efold_rm(x,p,l) - rmData%xend

  end function sr_iniinf_rm
 


  function rm_potential(matter,rmData)
    use infprec, only : transfert    
    implicit none
    real(kp), intent(in) :: matter    
    type(transfert), optional, intent(inout) :: rmData
    real(kp) :: rm_potential
    real(kp) :: p,mu,nu,x

    p=rmData%real1
    mu=rmData%real2
    nu=rmData%real3

    x = matter/mu
    
    rm_potential = 1._kp + nu*(1._kp/p - log(x))*matter**p

  end function rm_potential



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!kklmmt models with m2=0: V = M^4 / [1 + (mu/phi)^p]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  function slowroll_initial_matter_kklt(infParam,efoldWanted)
    use infprec, only : transfert, tolkp
    use inftools, only : zbrent    
    use infbgmodel, only : infbgparam, matterParamNum
    implicit none
    type(infbgparam), intent(in) :: infParam
    real(kp), optional, intent(in) :: efoldWanted
    real(kp), dimension(matterNum) :: slowroll_initial_matter_kklt

    type(transfert) :: kkltData
    real(kp), parameter :: efoldDefault = efoldBound
    real(kp), parameter :: tolFind = tolkp
    
    real(kp), dimension(2) :: fieldStop, fieldUv
    real(kp) :: matterOverMuUv, matterOverMuString
    real(kp) :: matterOverMuEps, matterOverMuIni, matterOverMuEnd
    real(kp) :: mini,maxi
    real(kp) :: p, mu
    real(kp) :: efold

  
    p = infParam%consts(2)
    mu = infParam%consts(3)

!sanity checks
    if (infParam%consts(3).le.0._kp) then
       stop 'slowroll_initial_matter_kklt: improper consts(3)'
    endif

    if (p.lt.2._kp) then
       write(*,*) 'slowroll_initial_matter_kklt: p = ',p
       stop
    endif

    if (infParam%consts(5).ne.-1._kp) then
       stop 'slowroll_initial_matter_kklt: improper consts(5)'
    endif


    if (present(efoldWanted)) then
       efold = efoldWanted
    else
       efold = efoldDefault
    endif


!upper field value bound
    fieldUv = field_thbound(infParam)
    matterOverMuUv = fieldUv(1)/mu


!inflation stops at matterOverMuEps1 (its determination is no accurate
!since eps2>1 in that region, but the number of efold in between is
!small. This condition could be replaced by matterOverMuEps2 as well

    mini = 0._kp
    maxi = 1._kp/epsilon(1._kp)

    kkltData%real1 = p
    kkltData%real2 = mu

    matterOverMuEps = zbrent(sr_endinf_kklt,mini,maxi,tolFind,kkltData)


!the branes collide at matterString
    fieldStop = field_stopinf(infParam)
    matterOverMuString = fieldStop(1)/mu
  
    if (matterOverMuString.lt.1.) then
       write(*,*)'slowroll_initial_matter_kklt: matterOverMuString < mu'
    endif


    matterOverMuEnd = max(matterOverMuString,MatterOverMuEps)


    if (matterOverMuEnd.gt.matterOverMuUv) then
       write(*,*)'slowroll_initial_matter_kklt: matterOverMuEnd > matterOverMuUv'
    endif

!find matterIni that gives the wanted number of efolds

    mini = 0._kp
    maxi = 1._kp/epsilon(1._kp)

    kkltData%real1 = p
    kkltData%real2 = p*(efold)/mu**2 + sr_efold_kklt(matterOverMuEnd,p)

    matterOverMuIni =  zbrent(sr_iniinf_kklt,mini,maxi,tolFind,kkltData)

    slowroll_initial_matter_kklt = matterOverMuIni * mu

    if (display) then
       write(*,*)
       write(*,*)'slowroll_initial_matter_kklt: (*kappa) (/mu)'       
       write(*,*)'matterEnd = ',matterOverMuEnd*mu, matterOverMuEnd
       write(*,*)'matterString  = ',matterOverMuString*mu, matterOverMuString
       write(*,*)'matterEps  = ',matterOverMuEps*mu, matterOverMuEps
       write(*,*)'matterIni  = ',matterOverMuIni*mu, matterOverMuIni
       write(*,*)'matterUv  = ',matterOverMuUv*mu, matterOverMuUv
       write(*,*)
    endif


  end function slowroll_initial_matter_kklt




  function sr_endinf_kklt(x,kkltData)
    use infprec, only : transfert    
    implicit none
    real(kp), intent(in) :: x    
    type(transfert), optional, intent(inout) :: kkltData
    logical, parameter :: endinfIsEpsOne=.true.
    real(kp) :: sr_endinf_kklt
    real(kp) :: p,mu

    p=kkltData%real1
    mu=kkltData%real2

    if (endinfIsEpsOne) then
!epsilon1=1
       sr_endinf_kklt = x**(p+1._kp) + x - p/(mu*sqrt(2._kp))
    else
!epsilon2=1
       sr_endinf_kklt = (x**(p+1._kp) + x)**2 - (2._kp*p/mu/mu) &
            *( (p+1._kp)*x**p + 1._kp)
    endif

  end function sr_endinf_kklt
 



  function sr_iniinf_kklt(x,kkltData)
    use infprec, only : transfert
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: kkltData
    real(kp) :: sr_iniinf_kklt
    real(kp) :: p

    p=kkltData%real1

    sr_iniinf_kklt =  sr_efold_kklt(x,p) - kkltData%real2
   
  end function sr_iniinf_kklt




  function sr_efold_kklt(x,p)
    implicit none
    real(kp), intent(in) :: x,p
    real(kp) :: sr_efold_kklt
    
    sr_efold_kklt = 0.5_kp*x**2 + x**(p+2._kp)/(p+2._kp)

    
  end function sr_efold_kklt




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!mixed inflation M^4( phi^p + alpha phi^q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  function slowroll_initial_matter_mix(infParam,efoldWanted)
    use infprec, only : transfert,tolkp
    use inftools, only : zbrent    
    use infbgmodel, only : infbgparam
    implicit none
    type(infbgparam), intent(in) :: infParam
    real(kp), optional, intent(in) :: efoldWanted
    real(kp), dimension(matterNum) :: slowroll_initial_matter_mix

    type(transfert) :: mixData
    real(kp), parameter :: efoldDefault = efoldBound
    real(kp), parameter :: tolFind = tolkp
    
    
    real(kp) :: matterEnd, matterIni
    real(kp) :: mini,maxi
    real(kp) :: alpha, efold, p, q
        
#ifndef PP5
    p = infParam%consts(2)    
    q = infParam%consts(12)
    alpha = infParam%consts(6)
#else
    stop 'slowroll_initial_matter_mix: check potParamNum value!'
#endif
   
          
    if (alpha.lt.0._kp) then
       stop 'mixed inflation: infParam%consts(6) < 0!'
    endif

    if (present(efoldWanted)) then
       efold = efoldWanted
    else
       efold = efoldDefault
    endif

!find the end of inflation
      
    mini = epsilon(1._kp)
    maxi = 1._kp/epsilon(1._kp)

    mixData%real1 = p
    mixData%real2 = q
    mixData%real3 = alpha

    matterEnd = zbrent(sr_endinf_mix,mini,maxi,tolFind,mixData)
   


!find the initial field values efolds before

    mini = matterEnd - epsilon(1._kp)
    maxi = 1._kp/epsilon(1._kp)

    mixData%real1 = p
    mixData%real2 = q
    mixData%real3 = alpha
    mixData%real4 = efold + sr_efold_mix(matterEnd,p,q,alpha)
    
    matterIni = zbrent(sr_iniinf_mix,mini,maxi,tolFind,mixData)
    
    slowroll_initial_matter_mix = matterIni

    if (display) then
       write(*,*)
       write(*,*)'slowroll_initial_matter_mix: (*kappa)'
       write(*,*)'matterEnd = ',matterEnd
       write(*,*)'matterIni = ',matterIni
       write(*,*)
    endif


  end function slowroll_initial_matter_mix




  function sr_endinf_mix(x,mixData)
    use infprec, only : transfert    
    implicit none
    real(kp), intent(in) :: x    
    type(transfert), optional, intent(inout) :: mixData
    real(kp) :: sr_endinf_mix
    real(kp) :: p,q,alpha

    p=mixData%real1
    q=mixData%real2
    alpha=mixData%real3

    sr_endinf_mix = alpha*sqrt(2._kp)*x**(q-p+1._kp) - alpha*q*x**(q-p) &
         + sqrt(2._kp)*x - p

  end function sr_endinf_mix
 



  function sr_iniinf_mix(x,mixData)
    use infprec, only : transfert
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: mixData
    real(kp) :: sr_iniinf_mix
    real(kp) :: p,q,alpha

    p=mixData%real1
    q=mixData%real2
    alpha=mixData%real3

    sr_iniinf_mix = sr_efold_mix(x,p,q,alpha) - mixData%real4
   
  end function sr_iniinf_mix


  function sr_efold_mix(x,p,q,alpha)
    use specialinf, only : hypergeom_2F1
    implicit none
    real(kp), intent(in) :: x,p,q,alpha
    real(kp) :: sr_efold_mix
    
    sr_efold_mix = 0.5_kp*x**2/(p*q)*(p + (q-p) &
         *hypergeom_2F1(1._kp,2._kp/(q-p),1._kp + 2._kp/(q-p),-alpha*q*x**(q-p)/p))
    
  end function sr_efold_mix


end module infsrmodel
