module inftorad
  use infunits
  use infprec, only : kp
  use infbgmodel, only : fieldNum
  use infbg, only : infbgphys
  implicit none

  private


!we use this everywhere
  type inftoradcosmo
     real(kp) :: bfoldIni
     real(kp) :: bfoldEnd 
     real(kp) :: efoldEndToToday
     real(kp) :: lnEnergyEnd
     type(infbgphys) :: bgIni
     type(infbgphys) :: bgEnd
  end type inftoradcosmo


!physical quantities at hubble exit
  type infhubblexit
     real(kp) :: kmpc
     real(kp) :: bfold
     real(kp) :: hubble
     real(kp) :: epsilon1
     real(kp) :: epsilon1JF
  end type infhubblexit


!for debugging
  logical, parameter :: display = .false.
  logical, parameter :: dump_file = .false.


!some physical constants
  real(kp), parameter :: scaleFactorToday = 1._kp  

! ln[1Mpc/sqrt(8pi)/lPl] = 130.282
!  real(kp), parameter :: lnMpcToKappa = 130.282_kp

!enforce constant equation of state during reheating for large field models
  logical, parameter :: LargeFieldWreh = .false.

  public inftoradcosmo, print_inftoradcosmo
  public scaleFactorToday, lnMpcToKappa
 
  public set_inftorad_cosmo, bfold_hubble_fraction

  
  public infhubblexit, hubble_splinexit
  public print_infhubblexit

contains


  subroutine print_inftoradcosmo(radvar,inname)
    use infbg, only : print_infbgphys
    implicit none
    type(inftoradcosmo), intent(in) :: radvar
    character(len=*), intent(in), optional :: inname
    
    if (present(inname)) then
       write(*,*)'type inftoradcosmo: ',inname
    endif
    write(*,*)'bfoldIni=          ',radvar%bfoldIni
    write(*,*)'bfoldEnd=          ',radvar%bfoldEnd
    write(*,*)'efoldEndToToday=   ',radvar%efoldEndToToday
    write(*,*)'lnRhoEnd=          ',radvar%lnEnergyEnd
    call print_infbgphys(radvar%bgIni,'bgIni=')
    call print_infbgphys(radvar%bgEnd,'bgEnd=')

  end subroutine print_inftoradcosmo



  subroutine print_infhubblexit(hexvar,inname)    
    implicit none
    type(infhubblexit), intent(in) :: hexvar
    character(len=*), intent(in), optional :: inname
    
    if (present(inname)) then
       write(*,*)'type infhubblexit: ',inname
    endif
    write(*,*)'kmpc=          ',hexvar%kmpc
    write(*,*)'bfold=         ',hexvar%bfold
    write(*,*)'hubble=        ',hexvar%hubble
    write(*,*)'epsilon1=      ',hexvar%epsilon1
    write(*,*)'epsilon1(JF)=  ',hexvar%epsilon1JF
!    write(*,*)'epsilon2=      ',hexvar%epsilon2

  end subroutine print_infhubblexit



  function set_inftorad_cosmo(bgParam,bgIni,bgEnd,lnReheat,inferror)
    use infbgmodel, only : infbgparam
    use infbg, only : infbgdata, infbgphys
    use infbgfunc, only :  matter_energy_density,  matter_energy_density_JF
    use infbgspline, only : set_infbg_spline  
    implicit none

    type(inftoradcosmo) :: set_inftorad_cosmo

    type(infbgparam), intent(in) :: bgParam
    type(infbgphys), intent(in) :: bgIni
    type(infbgphys), intent(in) :: bgEnd
    
!ln(aend/areh) + (1/4)ln(rhoend/rhoreh) + (1/4)ln(rhoend). This is
!ln(R), the rescaled reheating parameter, not ln(Rrad)
    real(kp), intent(in) :: lnReheat

!useful for setting hard prior
    integer, optional :: inferror

!rho at the end of inflation   
    real(kp) :: lnEnergyEndInf

!bound on lnReheat
    real(kp) :: lnReheatMin, lnReheatMax

!deviation from rad-like reheating: lnRrad = lnR - 1/4 ln(rhoend)
    real(kp) :: lnRrad

    real(kp) :: wreh, p, thirdminusw
    real(kp) :: ThirdMinusWlnEnergyEndReh, lnEnergyEndReh

!kappaeff^4 x rhonuc with rhonuc~10MeV
!energyNuc = 2.9d-82
!    real(kp), parameter :: lnEnergyNuc = -187.747
    real(kp), parameter :: lnEnergyNuc = lnRhoNuc

    real(kp), dimension(fieldNum) :: velocityEnd, fieldEnd
   
    fieldEnd = bgEnd%field
    velocityEnd = bgEnd%fieldDot*bgEnd%hubble

!if epsilon1=1 at the end of inflation, it should be Hend^2/3 (for one
!field models)

    lnEnergyEndInf = log(matter_energy_density(fieldEnd,velocityEnd))

!comes from -1/3 <wreheat < 1 and rhorheat > rhonuc

    lnReheatMax = (1._kp/12._kp)* (-lnEnergyNuc) &
         - (1._kp/3._kp) * (-lnEnergyEndInf)

    lnReheatMin = - (1._kp/4._kp)* (-lnEnergyNuc)
    
    
    set_inftorad_cosmo%lnEnergyEnd = lnEnergyEndInf  
    set_inftorad_cosmo%bfoldIni = bgIni%efold - bgEnd%efold
    set_inftorad_cosmo%bfoldEnd = 0._kp

    set_inftorad_cosmo%efoldEndToToday = - ln_scale_factor_scale_inv(lnEnergyEndInf) &
         - lnReheat

!TODELETE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    set_inftorad_cosmo%efoldEndToToday = 60._kp
!    write(*,*)'efoldEndToToday fixed to: ',set_inftorad_cosmo%efoldEndToToday
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    set_inftorad_cosmo%bgIni = bgIni
    set_inftorad_cosmo%bgEnd = bgEnd


!hard priors: inflation should occur at energy scale higher than bbn,
!as well as reheating
    if (lnEnergyEndInf.le.lnEnergyNuc) then
       write(*,*)'set_inftorad_cosmo: inflation till bbn (111)'
       if (display) then
          write(*,*)'lnRhoEndInf = ',lnEnergyEndInf
          write(*,*)'lnRhoNuc = ',lnEnergyNuc
       endif
       if (present(inferror)) then          
          inferror = 111
          return
       else          
          stop
       endif
    endif

    if (lnReheat.gt.lnReheatMax) then
       write(*,*)'set_inftorad_cosmo: higher reheating limit reached (w=1 rho=bbn) (112)'
       write(*,*)'lnReheat = ',lnReheat
       write(*,*)'lnReheatMax= ',lnReheatMax       
       if (present(inferror)) then
          inferror = 112
          return
       else
          stop
       endif
    endif

    if (lnReheat.lt.lnReheatMin) then
       write(*,*)'set_inftorad_cosmo: lower reheating limit reached (w=-1/3 rho=bbn) (113)'
       write(*,*)'lnReheat = ',lnReheat
       write(*,*)'lnReheatMin = ',lnReheatMin       
       if (present(inferror)) then
          inferror = 113
          return
       else
          stop
       endif
    endif

    if (LargeFieldWreh) then
       
       if (bgParam%name.ne.'largef') then
          stop 'set_inftorad_cosmo: incompatible models in inftorad!'
       endif

       p = bgParam%consts(2)
       if (p.gt.1) then
          wreh = (p-2._kp)/(p+2._kp)
          thirdminusw = 1._kp/3._kp - wreh                 
       else
          stop 'set_inftorad_cosmo: reheating cannot be parametric when p<1'
       endif

       write(*,*)'set_inftorad_cosmo: enforcing wreh = (p-2)/(p+2) =',wreh

!reheating is radiation
       if (thirdminusw.eq.0._kp) then
          write(*,*)'set_inftorad_cosmo: reheating is radiation: no constraint'
          return
       endif

       lnRrad = lnReheat - 0.25_kp*lnEnergyEndInf

       thirdMinusWlnEnergyEndReh = thirdminusw_ln_energy_endreh(lnRrad,lnEnergyEndInf,wreh)
       lnEnergyEndReh = thirdMinusWlnEnergyEndReh/thirdminusw

       if (lnEnergyEndReh.gt.lnEnergyEndInf) then
          write(*,*)'set_inftorad_cosmo: higher cte reheating limit reached (rho=rhoinf) (114)'
          write(*,*)'lnRhoEndReh = ',lnEnergyEndReh
          write(*,*)'lnRhoEndInf = ',lnEnergyEndInf
          if (present(inferror)) then
             inferror = 114
             return
          else
             stop
          endif
       endif

       if (lnEnergyEndReh.lt.lnEnergyNuc) then
          write(*,*)'set_inftorad_cosmo: lower cte reheating limit reached (rho=bbn) (115)'
          write(*,*)'lnRhoEndReh = ',lnEnergyEndReh
          write(*,*)'lnRhoNuc = ',lnEnergyNuc
          if (present(inferror)) then
             inferror = 115
             return
          else
             stop
          endif
       endif
       
    endif


         
     
    if (display) then
       write(*,*)'--------------- set_inftorad_cosmo ---------------'
       write(*,*)
       write(*,*)'a0 = ',scaleFactorToday
       write(*,*)'bfoldIni = ',set_inftorad_cosmo%bfoldIni
       write(*,*)'bfoldEnd = ',set_inftorad_cosmo%bfoldEnd
       write(*,*)'efoldIni = ',set_inftorad_cosmo%bgIni%efold
       write(*,*)'efoldEnd = ',set_inftorad_cosmo%bgEnd%efold      
       write(*,*)'lnReheatMin = ',lnReheatMin
       write(*,*)'lnReheat = ',lnReheat
       write(*,*)'lnReheatMax = ',lnReheatMax
       write(*,*)       
       if (largeFieldWreh) then
          write(*,*)'wreh = ',wreh
          write(*,*) 'ln(EnergyEndReh)',lnEnergyEndReh/4
       endif
       write(*,*)'ln(EnergyEndInf) = ',lnEnergyEndInf/4
       write(*,*)'ln(EnergyNuc) = ',lnEnergyNuc/4
       write(*,*)
       write(*,*)'efoldEndToToday = ',set_inftorad_cosmo%efoldEndToToday
       write(*,*)'ln[1Mpc/kappaeff] = ',lnMpcToKappa
       write(*,*)       
       write(*,*)'fieldEnd = ',set_inftorad_cosmo%bgEnd%field
       write(*,*)
       write(*,*)'energyEnd^1/2 = ',exp(0.5*set_inftorad_cosmo%lnEnergyEnd)
       write(*,*)'hubbleEnd.3^1/2 = ',sqrt(3.)*set_inftorad_cosmo%bgEnd%hubble
       write(*,*)
       write(*,*)'-------------------------------------------------'
    endif


  end function set_inftorad_cosmo




  function ln_scale_factor_eq()
!all in unit of kappa = sqrt(8pi)/Mpl. H0 et OmegaRad are fixed since
!the incertainties on their measured values is irrelevant compared to
!inflationary params. They can be added in the chains but this would
!render all prim parameters as "slow" in cosmocmc.

    implicit none
        
    real(kp) :: ln_scale_factor_eq
    
!    real(kp), parameter :: HubbleSquareRootOf2OmegaRad = 8.00d-63 !x kappaendinf(#kappatoday)


!this is: 1/4 ln(rhoeq) + ln(aeq/a0)
!             = (1/4) ln(rhoeq) - ln(1+zeq)
!             = -(1/4) ln(2/3) + (1/2) ln[Heq/(1+zeq)^2]
! where used has been made of: rhoeq = (3/2) Heq^2 (in kappa unit)
!
! Heq/(1+zeq)^2 = sqrt(2 Omega0) H0 / sqrt(1+zeq)
! 1+zeq = Omega0/OmegaRadiation0 <---- photons + neutrinos
! H0 = 100h km/s/Mpc = 1.74717d-61 Mpl
! There is no dependency in h and Omega0 in Heq/(1+zeq)^2 = sqrt(2 OmegaRad0) H0 (x kappatoday)

    ln_scale_factor_eq = log(scaleFactorToday) &
         + 0.5_kp*log(HubbleSquareRootOf2OmegaRad) &
         - 0.25_kp*log(2._kp/3._kp)

  end function ln_scale_factor_eq





  function ln_scale_factor_scale_inv(lnEnergyEndInf)
!this is ln(aend/a0) when R=1
!in that case we have:
!    ln(aend/a0) =
! ln(aend/areh) - 1/4 ln(rhoreh) + 1/2 ln(rhoend) <--lnReheat
! + ln(aeq/a0) + 1/4 ln(rhoeq) - 1/2 ln(rhoend) <-- ln_scale_inv

    implicit none
    real(kp), intent(in) :: lnEnergyEndInf
    real(kp) :: ln_scale_factor_scale_inv
    
   
    ln_scale_factor_scale_inv = ln_scale_factor_eq() - 0.5_kp*lnEnergyEndInf

  end function ln_scale_factor_scale_inv
  





  function ln_scale_factor_radreheat(lnEnergyEndInf)
!This is ln(aend/a0) for Rrad=1. Rough approximation of the physical
!scale factor at the end of inflation. Assumes instantenous reheating
!and instantenous transitions, all in unit of kappa =
!sqrt(8pi)/Mpl. The input is supposed to be the energy transfered
!after inflation to the "matter sector", in our case this is the
!energy density of the matter field at the end of inflation. This
!assumes pure GR after the end of inflation.

    implicit none
    
    real(kp), intent(in) :: lnEnergyEndInf
    real(kp) :: ln_scale_factor_radreheat
       
! ln(aend/a0) = ln(aend/areh) + ln(areh/aeq) + ln(aeq/a0)
!             ~ (1/4) ln(rhoeq/rhoend) - ln(1+zeq)
!             = -(1/4) ln(2rhoend/3) + (1/2) ln[Heq/(1+zeq)^2]
! where used has been made of: rhoeq = (3/2) Heq^2 (in kappa unit)
!
! Heq/(1+zeq)^2 = sqrt(2 Omega0) H0 / sqrt(1+zeq)
! 1+zeq = Omega0/OmegaRadiation0 <---- photons + neutrinos
! H0 = 100h km/s/Mpc = 1.74717d-61 Mpl
! There is no dependency in h and Omega0 in Heq/(1+zeq)^2 = sqrt(2 OmegaRad0) H0 (x kappatoday)

    ln_scale_factor_radreheat = ln_scale_factor_eq() - 0.25_kp*lnEnergyEndInf

  end function ln_scale_factor_radreheat



  function ln_scale_factor_corr_radreheat(lnEnergyEndInf,lnEnergyEndReh,wreh)
!This is ln(Rrad): the correction from a rad-like reheating period
!when during reheating one has P=<wre> x rho

!          ln(aend/areh) + (1/4) ln(rhoend/rhoreheat) = (1/4) x
!          (w-1/3)/(w+1) x ln(rhoend/rhoreheat)

    implicit none
    
    real(kp), intent(in) :: lnEnergyEndInf, lnEnergyEndReh, wreh
    real(kp) :: ln_scale_factor_corr_radreheat
    
    ln_scale_factor_corr_radreheat = 0.25_kp * (wreh - 1._kp/3._kp)/(wreh + 1._kp) &
         * (lnEnergyEndInf - lnEnergyEndReh)

  end function ln_scale_factor_corr_radreheat
  




  function thirdminusw_ln_energy_endreh(lnRrad,lnEnergyEndInf,wreh)
!This is the inverse function. It gives (1/3 -w) ln(rhoreh) as a function of ln(Rrad) and ln(rhoend)
    implicit none
    
    real(kp), intent(in) :: lnRrad,lnEnergyEndinf,wreh
    real(kp) :: thirdminusw_ln_energy_endreh

    thirdminusw_ln_energy_endreh = 4._kp*(1._kp+wreh)*lnRrad + (1._kp/3._kp-wreh)*lnEnergyEndInf

  end function thirdminusw_ln_energy_endreh






  function bfold_plus_ln_hubble(bfold,cosmoData)
!returns N + ln[H(N)] - cosmoData%real1, required by zeros finder subroutine:
!zbrent
    use infprec, only : transfert
    use infbgspline, only : check_infbg_spline
    use infbgspline, only : splineval_hubble_parameter
    implicit none

    type(transfert), optional, intent(inout) :: cosmoData
    real(kp), intent(in) :: bfold
    real(kp) :: bfold_plus_ln_hubble

    real(kp) :: hubble    

    if (.not.check_infbg_spline()) then
       stop 'bfold_plus_ln_hubble: bgdata not splined!'
    endif

    hubble = splineval_hubble_parameter(bfold)

    bfold_plus_ln_hubble = bfold + log(hubble) - cosmoData%real1       
    
!    print *,'bfold',bfold,bfold_plus_ln_hubble,cosmoData%real1
!    print *,'hubble^2',hubble**2

  end function bfold_plus_ln_hubble






  function bfold_hubble_fraction(kmpc,infCosmo,kphysOverHubble,inferror)
!return the bfold at which k/aH = kphysOverHubble
    use infsolvers, only : zbrent
    use infprec, only : transfert
    implicit none
    real(kp), intent(in) :: kmpc
    real(kp), intent(inout) :: kphysOverHubble
    type(inftoradcosmo), intent(in) :: infCosmo
    real(kp) :: bfold_hubble_fraction

!hard prior if quantum ic cannot be set
    integer, optional :: inferror

    type(transfert) :: cosmoData    
    real(kp), parameter :: tolBfold = 1e-6
    real(kp) :: bfoldSign


    cosmoData%real1 = log(kmpc) + infCosmo%efoldEndToToday - lnMpcToKappa &
         - log(kphysOverHubble)

    bfoldSign = bfold_plus_ln_hubble(infCosmo%bfoldIni,cosmoData)
 
    if (bfoldSign.ge.0d0) then
       write(*,*)'no crossing found (444) for mode k = ',kmpc, 'at k/aH = ' &
            ,kphysOverHubble       
       write(*,*)'bfoldIni = ',infCosmo%bfoldIni
       if (present(inferror)) then
          inferror = 222
          return
       else
          if (bfoldSign.lt.log(kphysOverHubble)) then
             write(*,*)'setting initial condition at bfoldIni'
             write(*,*)'ln(k/aH) = ', -bfoldSign + log(kphysOverHubble)
             kphysOverHubble = exp(-bfoldSign) * kphysOverHubble            
             bfold_hubble_fraction = infCosmo%bfoldIni
             return
          else
             stop
          endif
       endif
    endif

    bfold_hubble_fraction = zbrent(bfold_plus_ln_hubble,infCosmo%bfoldIni &
         , infCosmo%bfoldEnd ,tolBfold, cosmoData)


  end function bfold_hubble_fraction




  function hubble_splinexit(infCosmo,kmpc)
!spline evaluation of physical quantities at hubble exit (k/aH=1) for the mode kmpc    
    use infio
    use infbgspline, only : splineval_hubble_parameter
    use infbgspline, only : splineval_epsilon1
    use infbgspline, only : splineval_epsilon1JF
    implicit none

    real(kp), intent(in) :: kmpc
    type(infhubblexit) :: hubble_splinexit
    type(inftoradcosmo) :: infCosmo
    real(kp) :: bfoldExit, hubbleExit
    real(kp) :: epsilon1Exit, epsilon1JFExit
    real(kp) :: kphysOverHubbleExit
    real(kp) :: pi
    
    pi = acos(-1._kp)
    kphysOverHubbleExit = 1._kp

    bfoldExit = bfold_hubble_fraction(kmpc,infCosmo,kphysOverHubbleExit)
    hubbleExit = splineval_hubble_parameter(bfoldExit)
    epsilon1Exit = splineval_epsilon1(bfoldExit)
    epsilon1JFExit = splineval_epsilon1JF(bfoldExit)

    if (display) then
!this is the tensor spectrum, up to a negligeable relaxation after Hubble exit
       write(*,*)
       write(*,*)'kmpc = ',kmpc
       write(*,*)'hubble crossing at bfold = ',bfoldExit
       write(*,*)'2(Hk/pi)^2 = ',(2._kp/pi/pi)*hubbleExit**2 
       write(*,*)'epsilon1 = ',epsilon1Exit
       write(*,*)'epsilon1JF = ',epsilon1JFExit
       write(*,*)
    endif

    if (dump_file) then
       call livewrite('hubblexit.dat',kmpc,(2._kp/pi/pi)*hubbleExit**2)
       call livewrite('epsilonexit.dat',kmpc,epsilon1Exit,epsilon1JFExit)
       call livewrite('bfoldexit.dat',kmpc,bfoldExit)
    endif

    hubble_splinexit%kmpc = kmpc
    hubble_splinexit%bfold = bfoldExit
    hubble_splinexit%hubble = hubbleExit
    hubble_splinexit%epsilon1 = epsilon1Exit
    hubble_splinexit%epsilon1JF = epsilon1JFExit

  end function hubble_splinexit


end module inftorad
