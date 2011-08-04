!solves the slow-roll evolution numerically

module srevol
  use infprec, only : kp,tolkp

  implicit none

  public


contains



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!large field models
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function sr_matter_endinf_lf(p)
    implicit none
    real(kp), intent(in) :: p
    real(kp) :: sr_matter_endinf_lf

    sr_matter_endinf_lf = p/sqrt(2._kp)

  end function sr_matter_endinf_lf


  function sr_matter_lf(bfold,p)    
    implicit none
    real(kp), intent(in) :: bfold, p
    real(kp) :: sr_matter_lf,matterEnd
       
    matterEnd = sr_matter_endinf_lf(p)

    sr_matter_lf = sqrt(-2._kp*p*bfold + matterEnd**2)
    
  end function sr_matter_lf



  function sr_epsilon1_lf(matter,p)
    use infprec, only : kp
    implicit none
    real(kp) :: sr_epsilon1_lf
    real(kp), intent(in) :: matter,p
    
    sr_epsilon1_lf = 0.5_kp*(p/matter)**2        

  end function sr_epsilon1_lf



  function sr_epsilon2_lf(matter,p)
    use infprec, only : kp
    implicit none
    real(kp) :: sr_epsilon2_lf
    real(kp), intent(in) :: matter,p
    
    sr_epsilon2_lf = 2._kp*p/matter**2
    
  end function sr_epsilon2_lf


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!small field models
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  function sr_matterovermu_endinf_sf(p,mu)
    use infprec, only : kp,transfert
    use inftools, only : zbrent
    use infsrmodel, only : sr_endinf_sf
    implicit none
    real(kp), intent(in) :: p,mu
    real(kp) :: sr_matterovermu_endinf_sf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: sfData

    mini = epsilon(1._kp)
    maxi = 1._kp + epsilon(1._kp)

    sfData%real1 = p
    sfData%real2 = mu

    sr_matterovermu_endinf_sf = zbrent(sr_endinf_sf,mini,maxi,tolFind,sfData)
   
  end function sr_matterovermu_endinf_sf



    
  function sr_matterovermu_sf(bfold,matterOverMuEnd,p,mu)
    use infprec, only : kp, transfert
    use inftools, only : zbrent
    use infsrmodel, only : sr_efold_sf,sr_iniinf_sf
    implicit none
    real(kp), intent(in) :: bfold, p, mu, matterOverMuEnd
    real(kp) :: sr_matterovermu_sf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: sfData

    mini = epsilon(1._kp)
    maxi = matterOverMuEnd

    sfData%real1 = p
    sfData%real2 = -2._kp*p*bfold/mu**2 + sr_efold_sf(matterOverMuEnd,p)
    
    sr_matterovermu_sf = zbrent(sr_iniinf_sf,mini,maxi,tolFind,sfData)
    
   
  end function sr_matterovermu_sf
 



  function sr_lnreheat_matterovermu_sf(kmpc,matterOverMuStar,matterOverMuEnd,p,mu)
    use infprec, only : kp
    use infsrmodel, only : sr_efold_sf
    implicit none
    real(kp), intent(in) :: kmpc,matterOverMuStar,matterOverMuEnd,p,mu
    real(kp) :: sr_lnreheat_matterovermu_sf

    real(kp) :: Nfold ,lnC, potEnd, potStar

! ln[1Mpc/sqrt(8pi)/lPl] = 130.282
  real(kp), parameter :: lnMpcToKappa = 130.282_kp
  real(kp), parameter :: lnHubbleSquareRootOf2OmegaRad = -142.983_kp

    potStar = 1._kp - matterOverMuStar**p
    potEnd = 1._kp - matterOverMuEnd**p

    Nfold = (sr_efold_sf(matterOverMuStar,p) &
         - sr_efold_sf(matterOverMuEnd,p))*mu**2/(2.*p)

!    print *,'Nfold=',Nfold

    lnC = lnMpcToKappa + 0.5*lnHubbleSquareRootOf2OmegaRad &
         - 0.25*log(6._kp)

!    print *,'lnC',lnC

    sr_lnreheat_matterovermu_sf = log(kmpc) + Nfold - 0.5*log(potStar/potEnd)&
         - lnC

  end function sr_lnreheat_matterovermu_sf


  
  
  function sr_matterovermu_epstwo_sf(p,mu,epsilon2)
    use infprec, only : kp,transfert
    use inftools, only : zbrent
    use infsrmodel, only : sr_efold_sf,sr_endinf_sf
    implicit none
    real(kp), intent(in) :: p,mu,epsilon2
    real(kp) :: sr_matterovermu_epstwo_sf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: sfData

    mini = epsilon(1._kp)
    maxi = 1._kp + epsilon(1._kp)

    sfData%real1 = p
    sfData%real2 = mu
    sfData%real3 = epsilon2

    sr_matterovermu_epstwo_sf = zbrent(sr_findeps2_sf,mini,maxi,tolFind,sfData)
   
  end function sr_matterovermu_epstwo_sf




  function sr_epsilon1_sf(matterOverMu,p,mu)
    use infprec, only : kp
    implicit none
    real(kp) :: sr_epsilon1_sf
    real(kp), intent(in) :: matterOverMu,p,mu
    
    sr_epsilon1_sf = 0.5_kp*(p/mu * (matterOverMu**(p-1._kp))/(1._kp-matterOverMu**p))**2
    
  end function sr_epsilon1_sf




  function sr_epsilon2_sf(matterOverMu,p,mu)
    use infprec, only : kp
    implicit none
    real(kp) :: sr_epsilon2_sf
    real(kp), intent(in) :: matterOverMu,p,mu
    
    sr_epsilon2_sf = 2._kp*(p/mu**2)*matterOverMu**(p-2._kp) &
         * (p-1._kp + matterOverMu**p)/(1._kp-matterOverMu**p)**2
    
  end function sr_epsilon2_sf



  function sr_findeps2_sf(x,sfData)
    use infprec, only : kp, transfert
    implicit none
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: sfData
    real(kp) :: sr_findeps2_sf

    real(kp) :: p,mu,epsilon2

    p = sfData%real1
    mu = sfData%real2
    epsilon2 = sfData%real3

    sr_findeps2_sf = sr_epsilon2_sf(x,p,mu) - epsilon2

  end function sr_findeps2_sf



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!kklmmt in the small field like approximation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!use is made of the small field functions since kklmmt = sf with p -> -p

 function sr_matterovermu_endinf_kksf(p,mu)
    use infprec, only : kp, transfert       
    use inftools, only : zbrent
    use infsrmodel, only : sr_endinf_sf
    implicit none
    real(kp), intent(in) :: p,mu
    real(kp) :: sr_matterovermu_endinf_kksf   
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: sfData

    mini = 1._kp + epsilon(1._kp)
    maxi = 1._kp/epsilon(1._kp)

    sfData%real1 = -p
    sfData%real2 = mu

    sr_matterovermu_endinf_kksf = zbrent(sr_endinf_sf,mini,maxi,tolFind,sfData)
   
  end function sr_matterovermu_endinf_kksf



    
  function sr_matterovermu_kksf(bfold,matterOverMuEnd,p,mu)
    use infprec, only : kp, transfert    
    use inftools, only : zbrent
    use infsrmodel, only : sr_efold_sf,sr_iniinf_sf
    implicit none
    real(kp) :: sr_matterovermu_kksf
    real(kp), intent(in) :: bfold, p, mu, matterOverMuEnd           
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: sfData

    mini = matterOverMuEnd
    maxi = 1._kp/epsilon(1._kp)

    sfData%real1 = -p
    sfData%real2 = -2._kp*(-p)*bfold/mu**2 + sr_efold_sf(matterOverMuEnd,-p)
    
    sr_matterovermu_kksf = zbrent(sr_iniinf_sf,mini,maxi,tolFind,sfData)
    
  end function sr_matterovermu_kksf
 



  function sr_lnreheat_matterovermu_kksf(kmpc,matterOverMuStar,matterOverMuEnd,p,mu)
    use infprec, only : kp   
    implicit none
    real(kp), intent(in) :: kmpc,matterOverMuStar,matterOverMuEnd,p,mu
    real(kp) :: sr_lnreheat_matterovermu_kksf
 
    sr_lnreheat_matterovermu_kksf & 
         = sr_lnreheat_matterovermu_sf(kmpc,matterOverMuStar,matterOverMuEnd,-p,mu)

  end function sr_lnreheat_matterovermu_kksf


  
  
  function sr_matterovermu_epstwo_kksf(p,mu,epsilon2)
    use infprec, only : kp,transfert
    use inftools, only : zbrent
    use infsrmodel, only : sr_efold_sf,sr_endinf_sf
    implicit none
    real(kp), intent(in) :: p,mu,epsilon2
    real(kp) :: sr_matterovermu_epstwo_kksf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: sfData

    mini = 1._kp + epsilon(1._kp)
    maxi = 1._kp/epsilon(1._kp)

    sfData%real1 = -p
    sfData%real2 = mu
    sfData%real3 = epsilon2

    sr_matterovermu_epstwo_kksf = zbrent(sr_findeps2_sf,mini,maxi,tolFind,sfData)
   
  end function sr_matterovermu_epstwo_kksf




   function sr_epsilon1_kksf(matterOverMu,p,mu)
    use infprec, only : kp
    implicit none
    real(kp) :: sr_epsilon1_kksf
    real(kp), intent(in) :: matterOverMu,p,mu
    
    sr_epsilon1_kksf = sr_epsilon1_sf(matterOverMu,-p,mu)
    
  end function sr_epsilon1_kksf




  function sr_epsilon2_kksf(matterOverMu,p,mu)
    use infprec, only : kp
    implicit none
    real(kp) :: sr_epsilon2_kksf
    real(kp), intent(in) :: matterOverMu,p,mu
    
    sr_epsilon2_kksf = sr_epsilon2_sf(matterOverMu,-p,mu)
    
  end function sr_epsilon2_kksf



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!running mass models
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    
  function sr_matter_rm(bfold,matterStop,p,mu,nu)
    use infprec, only : kp, transfert
    use inftools, only : zbrent
    use infsrmodel, only : sr_efold_rm,sr_iniinf_rm
    implicit none
    real(kp), intent(in) :: bfold, p, mu, nu,matterStop
    real(kp) :: sr_matter_rm
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,matterOverMuStop,lambda
    type(transfert) :: rmData


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
       stop 'sr_matter_rm: error'
    endif

     matterOverMuStop = matterStop/mu
    
!    print *,'mini maxi',mini,maxi

    rmData%real1 = p
    rmData%real2 = mu
    rmData%real3 = nu

    lambda = nu * mu**p

    rmData%xend = sr_efold_rm(matterOverMuStop,p,lambda) &
         + 2._kp*p*(-bfold)/mu**2

    sr_matter_rm =  zbrent(sr_iniinf_rm,mini,maxi,tolFind,rmData)

  end function sr_matter_rm




   function sr_epsilon1_rm(matter,p,mu,nu)
    use infprec, only : kp
    implicit none
    real(kp) :: sr_epsilon1_rm
    real(kp), intent(in) :: matter,p,mu,nu
    
    real(kp) :: x
    
    x = matter

    sr_epsilon1_rm = (nu**2*p**2 * x**(-2 + 2*p) * Log(x/mu)**2) &
         /(2.*(1 + nu*x**p * (1/p - Log(x/mu)))**2)
    
  end function sr_epsilon1_rm



  function sr_epsilon2_rm(matter,p,mu,nu)
    use infprec, only : kp
    implicit none
    real(kp) :: sr_epsilon2_rm
    real(kp), intent(in) :: matter,p,mu,nu

    real(kp) :: x

    x = matter

    sr_epsilon2_rm = (2*nu*p**2*x**(-2 + p) &
         * (p + nu*x**p + Log(x/mu) &
         * ((-1 + p)*p + nu*x**p* (-1 + p*Log(x/mu))))) &
         /(p - nu*x**p * (-1 + p*Log(x/mu)))**2
    
  end function sr_epsilon2_rm



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!kklmmt in the AdS throat approximation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 function sr_matterovermu_endinf_kklt(p,mu)
    use infprec, only : kp, transfert       
    use inftools, only : zbrent
    use infsrmodel, only : sr_endinf_kklt
    implicit none
    real(kp), intent(in) :: p,mu
    real(kp) :: sr_matterovermu_endinf_kklt   
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: kkltData

    mini = 1._kp + epsilon(1._kp)
    maxi = 1._kp/epsilon(1._kp)

    kkltData%real1 = p
    kkltData%real2 = mu

    sr_matterovermu_endinf_kklt = zbrent(sr_endinf_kklt,mini,maxi,tolFind,kkltData)
   
  end function sr_matterovermu_endinf_kklt



    
  function sr_matterovermu_kklt(bfold,matterOverMuEnd,p,mu)
    use infprec, only : kp, transfert    
    use inftools, only : zbrent
    use infsrmodel, only : sr_efold_kklt,sr_iniinf_kklt
    implicit none
    real(kp) :: sr_matterovermu_kklt
    real(kp), intent(in) :: bfold, p, mu, matterOverMuEnd           
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: kkltData


    mini = matterOverMuEnd
    maxi = 1._kp/epsilon(1._kp)

    kkltData%real1 = p
    kkltData%real2 = p*(-bfold)/mu**2 + sr_efold_kklt(matterOverMuEnd,p)

    sr_matterovermu_kklt =  zbrent(sr_iniinf_kklt,mini,maxi,tolFind,kkltData)
    
  end function sr_matterovermu_kklt
 



  function sr_lnreheat_matterovermu_kklt(kmpc,matterOverMuStar,matterOverMuEnd,p,mu)
    use infprec, only : kp
    use infsrmodel, only : sr_efold_kklt
    implicit none
    real(kp), intent(in) :: kmpc,matterOverMuStar,matterOverMuEnd,p,mu
    real(kp) :: sr_lnreheat_matterovermu_kklt

    real(kp) :: Nfold ,lnC, potEnd, potStar

! ln[1Mpc/sqrt(8pi)/lPl] = 130.282
    real(kp), parameter :: lnMpcToKappa = 130.282_kp
    real(kp), parameter :: lnHubbleSquareRootOf2OmegaRad = -142.983_kp

    potStar = 1._kp/(1._kp + matterOverMuStar**p)
    potEnd = 1._kp/(1._kp + matterOverMuEnd**p)

    Nfold = (sr_efold_kklt(matterOverMuStar,p) &
         - sr_efold_kklt(matterOverMuEnd,p))*mu**2/p
    
!    print *,'Nfold=',Nfold

    lnC = lnMpcToKappa + 0.5*lnHubbleSquareRootOf2OmegaRad &
         - 0.25*log(6._kp)

!    print *,'lnC',lnC

    sr_lnreheat_matterovermu_kklt = log(kmpc) + Nfold - 0.5*log(potStar/potEnd)&
         - lnC

  end function sr_lnreheat_matterovermu_kklt





!  function sr_lnreheat_matterovermu_kklt(kmpc,matterOverMuStar,matterOverMuEnd,p,mu)
!    use infprec, only : kp   
!    implicit none
!    real(kp), intent(in) :: kmpc,matterOverMuStar,matterOverMuEnd,p,mu
!    real(kp) :: sr_lnreheat_matterovermu_kklt

!    write(*,*)'warning: using the small field formulae with p->-p'
 
!    sr_lnreheat_matterovermu_kklt & 
!         = sr_lnreheat_matterovermu_sf(kmpc,matterOverMuStar,matterOverMuEnd,-p,mu)

!  end function sr_lnreheat_matterovermu_kklt


  
  
  function sr_matterovermu_epstwo_kklt(p,mu,epsilon2)
    use infprec, only : kp,transfert
    use inftools, only : zbrent
    use infsrmodel, only : sr_efold_sf,sr_endinf_sf
    implicit none
    real(kp), intent(in) :: p,mu,epsilon2
    real(kp) :: sr_matterovermu_epstwo_kklt
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: kkltData

    mini = 1._kp + epsilon(1._kp)
    maxi = 1._kp/epsilon(1._kp)

    kkltData%real1 = p
    kkltData%real2 = mu
    kkltData%real3 = epsilon2

    sr_matterovermu_epstwo_kklt = zbrent(sr_findeps2_kklt,mini,maxi,tolFind,kkltData)
   
  end function sr_matterovermu_epstwo_kklt


  

  function sr_epsilon1_kklt(matterOverMu,p,mu)
    use infprec, only : kp
    implicit none
    real(kp) :: sr_epsilon1_kklt
    real(kp), intent(in) :: matterOverMu,p,mu
    
    sr_epsilon1_kklt = 0.5*(p/mu /(matterOverMu**(p+1._kp) + matterOverMu))**2
    
  end function sr_epsilon1_kklt



  function sr_epsilon2_kklt(matterOverMu,p,mu)
    use infprec, only : kp
    implicit none
    real(kp) :: sr_epsilon2_kklt
    real(kp), intent(in) :: matterOverMu,p,mu
    
    sr_epsilon2_kklt = 2.*(p/mu**2)*matterOverMu**(-p-2._kp) &
         * (p+1._kp + matterOverMu**(-p)) &
         / (1._kp + matterOverMu**(-p))**2
    
  end function sr_epsilon2_kklt




  function sr_findeps2_kklt(x,kkltData)
    use infprec, only : kp, transfert
    implicit none
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: kkltData
    real(kp) :: sr_findeps2_kklt

    real(kp) :: p,mu,epsilon2

    p = kkltData%real1
    mu = kkltData%real2
    epsilon2 = kkltData%real3

    sr_findeps2_kklt = sr_epsilon2_kklt(x,p,mu) - epsilon2

  end function sr_findeps2_kklt



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!mixed inflation model
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  function sr_matter_endinf_mix(p,q,alpha)
    use infprec, only : kp,transfert
    use inftools, only : zbrent
    use infsrmodel, only : sr_endinf_mix
    implicit none
    real(kp), intent(in) :: p,q,alpha
    real(kp) :: sr_matter_endinf_mix
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: mixData

    mini = epsilon(1._kp)
    maxi = 1._kp/epsilon(1._kp)

    mixData%real1 = p
    mixData%real2 = q
    mixData%real3 = alpha

    sr_matter_endinf_mix = zbrent(sr_endinf_mix,mini,maxi,tolFind,mixData)
   
  end function sr_matter_endinf_mix



    
  function sr_matter_mix(bfold,matterEnd,p,q,alpha)
    use infprec, only : kp, transfert
    use inftools, only : zbrent
    use infsrmodel, only : sr_efold_mix,sr_iniinf_mix
    implicit none
    real(kp), intent(in) :: bfold, p, q, alpha, matterEnd
    real(kp) :: sr_matter_mix
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: mixData

    mini = matterEnd - epsilon(1._kp)
    maxi = 1._kp/epsilon(1._kp)


    mixData%real1 = p
    mixData%real2 = q
    mixData%real3 = alpha
    mixData%real4 = -bfold + sr_efold_mix(matterEnd,p,q,alpha)
    
    sr_matter_mix = zbrent(sr_iniinf_mix,mini,maxi,tolFind,mixData)
       
  end function sr_matter_mix
 

  

  function sr_epsilon1_mix(matter,p,q,alpha)
    use infprec, only : kp
    implicit none
    real(kp) :: sr_epsilon1_mix
    real(kp), intent(in) :: matter,p,q,alpha
    
    sr_epsilon1_mix = 0.5_kp/matter**2 &
         *((p+alpha*q*matter**(q-p))/(1+alpha*matter**(q-p)))**2
    
  end function sr_epsilon1_mix




  function sr_epsilon2_mix(matter,p,q,alpha)
    use infprec, only : kp
    implicit none
    real(kp) :: sr_epsilon2_mix
    real(kp), intent(in) :: matter,p,q,alpha
    
    sr_epsilon2_mix = 2_kp/matter**2 &
         * (p-alpha*(p*p+q*(q-1)-p*(2*q+1))*matter**(q-p) &
         + alpha*alpha*q*(matter*matter)**(q-p)) &
         / (1+alpha*matter**(q-p))**2
    
  end function sr_epsilon2_mix




end module srevol
