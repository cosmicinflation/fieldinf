!This module provides the initial power spectra for CAMB, computed
!from the inflationary module

module InitialPower   
  use precision, only : dl
  use infprec, only : kp 
  use infbgmodel, only : infbgparam
  use infbg, only : infbgphys, infbgdata
  use inftorad, only : inftoradcosmo

  implicit none   

  private
 
  logical, parameter :: display=.true.
   
  character(LEN=*), parameter :: Power_Name = 'power_inf'
  integer, parameter :: nnmax= 1 

!mind that for mcmc
  logical, parameter :: checkStopDefault = .false.
  logical, parameter :: checkBoundDefault = .false.
  logical, parameter :: useSplineDefault = .true.


  type InitialPowerParams
     integer :: nn    
     integer :: bgParamNum
     logical :: checkStop     
     logical :: checkBound
     logical :: useSpline
     integer :: lnkmpcNum
     real(kp) :: lnkmpcMax, lnkmpcMin
     real(kp) :: lnReheat
     real(kp) :: kstar
     type(infbgparam) :: infParam    
  end type InitialPowerParams


  type InitialPowerData
     type(initialpowerparams) :: initP
     type(infbgphys) :: infIni
     type(infbgphys) :: infObs     
     type(infbgphys) :: infEnd     
     type(inftoradcosmo) :: infCosmo
     type(infbgdata), pointer :: ptrBgdata => null()
  end type InitialPowerData


  type ExportInfProp
     real(kp) :: lnEnergyEnd
     real(kp) :: efoldEndToToday
  end type ExportInfProp


  real(dl) :: curv  !Curvature contant, set in InitializePowers     



  type(InitialPowerData), save :: powerD

  
  

  interface operator (==)     
     module procedure inipowerparams_equal
  end interface

  interface operator (/=)
     module procedure inipowerparams_unequal
  end interface

  private operator(==),operator(/=)

  public nnmax
  public SetInfBg, SetInfBgSpline, SetInfCosmo, SetInfScalPow
  public InitializePowers, FreePowers, ScalarPower, TensorPower
  public InitialPowerParams,Power_Descript, Power_Name, SetDefPowerParams

  public exportinfprop, UpdateInfProp

  public InitialPower_ReadParams

contains
       




  function inipowerparams_equal(PinA, PinB)
    use infbgmodel, only : operator(==) 
    implicit none
    type(initialpowerparams), intent(in) :: PinA, PinB
    logical :: inipowerparams_equal

    inipowerparams_equal = ((PinA%infParam == PinB%infParam) &
         .and. (PinA%lnReheat == PinB%lnReheat) &
         .and. (PinA%nn == PinB%nn) &
         .and. (PinA%bgParamNum == PinB%bgParamNum) &
         .and. (PinA%checkStop .eqv. PinB%checkStop) &
         .and. (PinA%checkBound .eqv. PinB%checkBound) &
         .and. (PinA%useSpline .eqv. PinB%useSpline) &
         .and. (PinA%lnkmpcNum == PinB%lnkmpcNum) &
         .and. (PinA%lnkmpcMax == PinB%lnkmpcMax) &
         .and. (PinA%lnkmpcMin == PinB%lnkmpcMin))

  end function inipowerparams_equal




  function inipowerparams_unequal(PinA, PinB)
    implicit none
    type(initialpowerparams), intent(in) :: PinA, PinB
    logical :: inipowerparams_unequal

    inipowerparams_unequal = .not.(inipowerparams_equal(PinA,PinB))
    
  end function inipowerparams_unequal




  subroutine UpdateInfProp(export)
    implicit none
    type(exportinfprop), intent(out) :: export

    export%lnEnergyEnd = powerD%infCosmo%lnEnergyEnd
    export%efoldEndToToday = powerD%infCosmo%efoldEndToToday
   
  end subroutine UpdateInfProp

  



  subroutine SetDefPowerParams(Pin)
    use infbgmodel, only : fieldNum
    use infdilaton, only : dilatonNum
    use infbgmodel, only : infParamNum
    implicit none
    type (InitialPowerParams), intent(out) :: Pin

    Pin%nn = 1
   
    Pin%bgParamNum = infParamNum

!stop inflation according to field values
    Pin%checkStop = checkStopDefault

!impose hard prior from theoretical bounds on field values
    Pin%checkBound = checkBoundDefault
   
!range and sampling of the power spectra spline
    Pin%useSpline = useSplineDefault 
    Pin%lnkmpcMin = -14.
    Pin%lnkmpcMax = 0.
    Pin%lnkmpcNum = 12
      
    Pin%lnReheat = 0.
    Pin%kstar = 0.05

!value for the parameters (see infbg.f90)
    Pin%infParam%name = 'largef'
    Pin%infParam%consts(1:infParamNum) = 0.

    if (dilatonNum.ne.0) then
       stop 'SetDefPowerParams: power_inf.f90 needs edition!'
!       Pin%infParam%conforms(1:dilatonNum) = 1.
    endif

    Pin%infParam%consts(1) = 1e-5
    Pin%infParam%consts(2) = 2.   
    Pin%infParam%consts(3:4) = 0.
    Pin%infParam%consts(5) = 1.

    if (infParamNum.gt.5) then
       Pin%infParam%consts(6:infParamNum) = 0.
    endif

  end subroutine SetDefPowerParams




  subroutine InitializePowers(Pin,acurv,wantScalars,wantTensors)
    implicit none

    type (initialpowerparams) :: Pin
         !Called before computing final Cls in cmbmain.f90
         !Could read spectra from disk here, do other processing, etc.
    real(dl) :: acurv
    logical, optional, intent(in) :: wantScalars, wantTensors
    integer :: inferror
    logical, parameter :: usePstar=.false.
    real(kp), parameter :: Pstar = 1._kp

    inferror = 0

    if (Pin%nn > nnmax) then
       stop 'can only used one initial power spectrum'
    end if
    curv=acurv         

    if (curv.ne.0d0) stop 'flat universe only'

!one background for all k (that's why nnmax=1)
!    print *,'we are in initialize power',(Pin%infParam == powerD%initP%infParam) &
!         ,(Pin%infParam /= powerD%initP%infParam)   
  
    call SetInfBg(Pin,inferror)
    if (inferror.ne.0) stop 'InitializePowers: unproper infbg'

    if (usePstar) then
!this is only for playing with normalised power spectra to Pstar. 
       write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       write(*,*)'InitializePower: renormalising P(k*) to Pstar=',Pstar
       write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       read(*,*)
       call SetInfScalPow(Pin,Pstar)
    endif

    call SetInfBgSpline(Pin)

    call SetInfCosmo(Pin,inferror)
    if (inferror.ne.0) stop 'InitializePowers: unproper inftorad'
      
  end subroutine InitializePowers




  subroutine SetInfBg(Pin,inferror)
    use infbgmodel, only : operator(==)
    use infbgmodel, only : set_infbg_param, matterNum
    use infbounds, only : field_stopinf
    use infbg, only : operator(==)
    use infbg, only : set_infbg_ini, bg_field_evol    
    implicit none

    type (initialpowerparams), intent(in) :: Pin    
    integer, optional :: inferror

    logical :: areTryParamsOk, areInfIniOk, isThBoundOk
    logical :: areParamsSet
    integer, parameter :: infbgPoints = 1000
    type(initialpowerparams) :: Pstack
    real(kp), dimension(2) :: fieldStop
    logical :: stopAtMax
   
!    print *,'Pin',Pin
    if (associated(powerD%ptrBgdata)) then
       if (Pin%infParam == powerD%initP%infParam) then
          if (display) then
             write(*,*)'SetInfBg: same infbg params'
             write(*,*)'Pin%infParam= ',Pin%infParam
          endif
          return
       else
          if (display) then
             write(*,*)'SetInfBg: new infbgparams'
             write(*,*)'Pin%infParam= ',Pin%infParam
          endif
          Pstack = powerD%initP
          call FreePowers(Pstack)
       endif
    endif

!hard prior tests on the parameters
    areTryParamsOk = HardPriorAcceptParam(Pin%infParam,inferror)
    if (.not.areTryParamsOk) then
       if (present(inferror).and.display) write(*,*)'SetInfBg: inferror ',inferror
       return
    endif

       
!update powerD%initP according to Pin (may be modified in)    
    areParamsSet = set_infbg_param(Pin%infParam)
    powerD%initP = Pin
    if (.not.(areParamsSet)) then
       stop 'SetInfBg: params initialisation failed!'
    endif

    powerD%infIni = set_infbg_ini(Pin%infParam)



!hard prior tests on the initial field values
    areInfIniOk = HardPriorAcceptInfIni(powerD%infIni,inferror)
    if (.not.areInfIniOk) then
       if (present(inferror).and.display) write(*,*)'SetInfBg: inferror ',inferror
       powerD%infEnd = powerD%infIni
       return
    endif



!hard prior tests for theoretical bound on initial field values 
    if (powerD%initP%checkBound) then

       isThBoundOk = HardPriorInThBound(powerD%initP%infParam &
            ,powerD%infIni,inferror)

       if (.not.isThBoundOk) then          
          powerD%infEnd = powerD%infIni
          if (present(inferror).and.display) write(*,*)'SetInfBg: inferror ',inferror
          return
       endif

    endif


    if (powerD%initP%checkStop) then
       fieldStop = field_stopinf(powerD%initP%infParam)
       if (fieldStop(2).gt.0._kp) then
          stopAtMax = .true.
       else
          stopAtMax = .false.
       endif

       powerD%infEnd = bg_field_evol(powerD%infIni,infbgPoints &
            ,powerD%infObs,powerD%ptrBgdata,fieldStop(1),stopAtMax)
    else

       powerD%infEnd = bg_field_evol(powerD%infIni,infbgPoints &
            ,powerD%infObs,powerD%ptrBgdata)

    endif

!last hard prior: the number of efolds is not enough
    if (powerD%infEnd == powerD%infIni) then
       if (display) write(*,*)'SetInfBg: infEnd = infIni (4444)'
       if (present(inferror)) inferror = 4444
       return
    endif
   

  end subroutine SetInfBg



  function HardPriorAcceptParam(infParam,inferror)
    use infbgmodel, only : matterParamNum
    implicit none
    logical :: HardPriorAcceptParam
    type(infbgparam), intent(in) :: infParam
    integer, optional :: inferror

    HardPriorAcceptParam = .true.

!if during mcmc, consts(4)=0 is unprobably tested (singular for this value)
    if (infParam%name=='runmas') then
       if (infParam%consts(4).eq.0._kp) then
          if (display) then
             write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
             write(*,*)'HardPriorAcceptParam: infParam%consts(4)=0 (1111)'
             write(*,*)'for running mass model!'
             write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          endif
          if (present(inferror)) inferror = 1111
          HardPriorAcceptParam=.false.
       endif
    endif

    if (infParam%name=='kklmmt') then
       if (infParam%consts(matterParamNum).lt.infParam%consts(3)) then
          if (display) write(*,*)'HardPriorAcceptParam: MatterString < mu (1112)'
          if (present(inferror)) inferror = 1112
          HardPriorAcceptParam=.false.
       endif
    endif


  end function HardPriorAcceptParam



  function HardPriorAcceptInfIni(infIni,inferror)
    implicit none
    logical :: HardPriorAcceptInfIni
    type(infbgphys), intent(in) :: infIni
    integer, optional :: inferror

    HardPriorAcceptInfIni=.true.

    if (infIni%epsilon1.gt.1._kp) then            
       if (display) write(*,*) 'HardPriorAcceptInfIni: epsilon1 > 1 (2222) ' 
       if (present(inferror)) inferror = 2222
       HardPriorAcceptInfIni = .false.
    endif

  end function HardPriorAcceptInfIni




  function HardPriorInThBound(infParam,infIni,inferror)
    use infbgmodel, only : matterNum
    use infbounds, only : field_thbound
    implicit none
    logical :: HardPriorInThBound
    type(infbgparam), intent(in) :: infParam
    type(infbgphys), intent(in) :: infIni
    integer, optional :: inferror

    real(kp), dimension(2) :: fieldBound

    HardPriorInThBound = .true.

    fieldBound = field_thbound(infParam)
    if (fieldBound(2).gt.0._kp) then
       HardPriorInThBound = (maxval(infIni%field(1:matterNum)) &
            .lt.fieldBound(1))
    else
       HardPriorInThBound = (minval(infIni%field(1:matterNum)) &
            .lt.fieldBound(1))
    endif
    
    if (.not.HardPriorInthBound) then
       if (display) write(*,*) 'HardPriorInThBound: bound violated (3333)'
       if (present(inferror)) inferror = 3333
    endif

  end function HardPriorInThBound
  



  subroutine SetInfBgSpline(Pin)
    use infbgmodel, only : operator(==)
    use infbgspline, only : set_infbg_spline,check_infbg_spline
    implicit none

    type (initialpowerparams), intent(in) :: Pin    
   
    if (check_infbg_spline()) return

!sanity checks
    if (Pin%infParam == powerD%initP%infParam) then
       call set_infbg_spline(powerD%infEnd,powerD%ptrBgdata)
    else
       stop 'SetInfBgSpline: Pin%params >< Data params'
    endif

  end subroutine SetInfBgSpline





  subroutine FreeInfBgSpline(Pin)
    use infbgmodel, only : operator(==)
    use infbgspline, only : free_infbg_spline, check_infbg_spline
    implicit none

    type (initialpowerparams), intent(in) :: Pin    
   
    if (.not. check_infbg_spline()) return

!sanity checks
    if (Pin%infParam == powerD%initP%infParam) then
       call free_infbg_spline()
    else
       stop 'FreeInfBgSpline: Pin%params >< Data params'
    endif

  end subroutine FreeInfBgSpline





  subroutine FreeInfBgData(Pin)
    use infbgmodel, only : operator(==)
    use infbg, only : free_infbg_data
    implicit none

    type (initialpowerparams), intent(in) :: Pin    
   
    if (.not.associated(powerD%ptrBgdata)) stop 'FreeInfBgData: no bgdata found!'

!sanity checks
    if (Pin%infParam == powerD%initP%infParam) then
       call free_infbg_data(powerD%ptrBgdata)      
    else
       stop 'FreeInfBgData: Pin%params >< Data params'
    endif

  end subroutine FreeInfBgData
  




  subroutine SetInfCosmo(Pin,inferror)
    use infbgmodel, only : operator(/=)   
    use inftorad, only : set_inftorad_cosmo
    use infpowspline, only : check_power_scal_spline, check_power_tens_spline
    
    implicit none

    type (initialpowerparams), intent(in) :: Pin        
    integer, optional :: inferror

!sanity checks
    if (Pin%infParam /= powerD%initP%infParam) then
       stop 'SetInfCosmo: Pin%params >< Data params'
    endif

    
    if (Pin /= powerD%initP) then 
       if (Pin%useSpline) call FreePowSpline(Pin)
       if (display) then
          write(*,*)'SetInfCosmo: new inftorad params'
          write(*,*)'lnReheat = ',Pin%lnReheat
       endif
       powerD%initP = Pin
    else
       if (display) then
          write(*,*)'SetInfCosmo: same inftorad params'
          write(*,*)'lnReheat = ',Pin%lnReheat
       endif
    endif

    powerD%infCosmo = set_inftorad_cosmo(powerD%initP%infParam &
         ,powerD%infObs,powerD%infEnd,powerD%InitP%lnReheat,inferror)

    if (present(inferror)) then
       if (inferror.ne.0) return
    endif

    if (Pin%useSpline) then
       call SetScalPowSpline(Pin,inferror)
       call SetTensPowSpline(Pin,inferror)
    endif

  end subroutine SetInfCosmo





  subroutine SetTensPowSpline(Pin,inferror)    
    use inftorad, only : bfold_hubble_fraction
    use infpowspline, only : set_power_tens_spline
    use infpowspline, only : check_power_tens_spline
    implicit none
    type(initialpowerparams), intent(in) :: Pin   
    integer, optional :: inferror

    real(kp), dimension(Pin%lnkmpcNum) :: lnkmpcKnots
    
    if (check_power_tens_spline()) return

    call IniPowSpline(Pin,lnkmpcKnots,inferror)

    if (present(inferror)) then
       if (inferror.ne.0) return
    endif

    call set_power_tens_spline(powerD%infCosmo,lnkmpcKnots)
    
!    if (present(inferror)) inferror = 0

  end subroutine SetTensPowSpline





  subroutine SetScalPowSpline(Pin,inferror)        
    use infpowspline, only : set_power_scal_spline
    use infpowspline, only : check_power_scal_spline
    implicit none
    type(initialpowerparams), intent(in) :: Pin   
    integer, optional :: inferror

    real(kp), dimension(Pin%lnkmpcNum) :: lnkmpcKnots

    if (check_power_scal_spline()) return
    
    call IniPowSpline(Pin,lnkmpcKnots,inferror)
       
    if (present(inferror)) then
       if (inferror.ne.0) return
    endif
    
    call set_power_scal_spline(powerD%infCosmo,lnkmpcKnots)       

  end subroutine SetScalPowSpline

  



  
  subroutine IniPowSpline(Pin,lnkmpcVec,inferror)
    use inftorad, only : bfold_hubble_fraction
    implicit none
    type(initialpowerparams), intent(in) :: Pin   
    real(kp), dimension(Pin%lnkmpcNum) :: lnkmpcVec
    integer, optional :: inferror

    integer :: i
    real(kp) :: kmpcMin, bfoldSmallest
    real(kp) :: kphysOverHubbleInit

    kphysOverHubbleInit = 100._kp

 !sanity checks
    if (.not.associated(powerD%ptrBgdata)) then
       stop 'SetInfPowSpline: no InfBg found, call SetInfBg before!'
    else
       if (Pin /= powerD%initP) then
          write(*,*)'SetInfPowSpline: Pin >< Data'
          stop      
       endif
    endif
    if (.not.Pin%useSpline) stop 'SetInfPowSpline: useSpline is F'
         
    do i=1,powerD%initP%lnkmpcNum
       lnkmpcVec(i) = powerD%initP%lnkmpcMin + real(i-1)*(powerD%initP%lnkmpcMax &
            - powerD%initP%lnkmpcMin)/real(powerD%initP%lnkmpcNum-1)
    enddo

!this test that initial conditions at kphysOverHubbleInit can be set,
!or in other words that there are enough efold to do it
    if (present(inferror)) then
       kmpcMin = exp(powerD%InitP%lnkmpcMin)      
       bfoldSmallest = bfold_hubble_fraction(kmpcMin,powerD%infCosmo &
            ,kphysOverHubbleInit,inferror)
    endif

  end subroutine IniPowSpline





  subroutine FreePowSpline(Pin)
    use infpowspline, only : check_power_scal_spline
    use infpowspline, only : check_power_tens_spline
    use infpowspline, only : free_power_scal_spline
    use infpowspline, only : free_power_tens_spline
    implicit none
    type(initialpowerparams), intent(in) :: Pin

    if (Pin%useSPline) then
       if (check_power_scal_spline()) call free_power_scal_spline()
       if (check_power_tens_spline()) call free_power_tens_spline()
    else
       stop 'FreePowSpline: useSpline is false!'
    endif
        
  end subroutine FreePowSpline
  




  subroutine SetInfScalPow(Pin,Pwanted)
!to normalise the scalar spectrum to Pwanted at kpivot
    use infbgmodel, only : operator(/=)
    use infbg, only : rescale_potential   
    use inftorad, only : set_inftorad_cosmo
    use infpert, only : scalNum
    use infpert, only : power_spectrum_scal
    implicit none
    
    type (initialpowerparams), intent(inout) :: Pin 
    real(kp), intent(in) :: Pwanted

    real(kp) :: scale
    real(kp) :: Pscal(scalNum,scalNum)
   
    integer :: ierror

!sanity checks
    if (Pin%infParam /= powerD%initP%infParam) then
       stop 'SetInfScalPow: Pin%params >< Data params'
    endif

    powerD%infCosmo = set_inftorad_cosmo(powerD%initP%infParam, &
         powerD%infObs,powerD%infEnd,powerD%InitP%lnReheat,ierror)
    
!the bg spline should be set to compute the power spectrum (if not
!already) TODO: not optimal at all
    call SetInfBgSpline(Pin)

    Pscal = power_spectrum_scal(powerD%infCosmo,Pin%kstar)
    
    if (Pscal(scalNum,scalNum).ne.0.) then
       scale = Pwanted/Pscal(scalNum,scalNum)
    else
       stop 'SetInfScalPow: Pscal = 0'
    endif

!rescale all the relevant quantities. After this call, everything
!should be as if the potential has been normalised according to Pwanted
    call rescale_potential(scale,powerD%initP%infParam,powerD%infIni &
         ,powerD%infEnd,powerD%infObs,powerD%ptrBgdata)

!update the external infParam
    Pin%infParam = powerD%initP%infParam

!but the background spline is not longer up to date
    call FreeInfBgSpline(Pin)
    
    
!for debugging only
    if (display) then
       write(*,*)
       write(*,*)'SetInfScalPow:'
       write(*,*)'Pwanted = ',Pwanted
       write(*,*)'Pnow = ',Pscal(scalNum,scalNum)
       call SetInfBgSpline(Pin)
       powerD%infCosmo = set_inftorad_cosmo(powerD%InitP%infParam &
            ,powerD%infObs,powerD%infEnd,powerD%InitP%lnReheat,ierror)
       Pscal = power_spectrum_scal(powerD%infCosmo,Pin%kstar)
       call FreeInfBgSpline(Pin)
       write(*,*)'Pafter = ',Pscal(scalNum,scalNum)
       write(*,*)
    endif
     

  end subroutine SetInfScalPow





  subroutine FreePowers(Pin)
    use infbgmodel, only : operator(==)  
    implicit none
    type(initialpowerparams), intent(in) :: Pin

    if (.not.associated(powerD%ptrBgdata)) then
       write(*,*) 'FreePowers: powerD%ptrBgdata not associated!'
       return
    endif

    if (Pin%infParam == powerD%initP%infParam) then       
       if (display) write(*,*) 'FreePowers: freeing infbg spline'
       call FreeInfBgSpline(Pin)

       if (display) write(*,*) 'FreePowers: freeing infbg data'
       call FreeInfBgData(Pin)

       if (powerD%initP%useSpline) then
          if (display) write(*,*) 'FreePowers: freeing powspline data'
          call FreePowSpline(Pin)          
       endif
    else
       write(*,*)
       write(*,*) 'FreePowers: Pin <> Data'
       write(*,*) 'FreePowers: Pin%Param = ',Pin%infParam
       write(*,*) 'FreePowers: Data%Param = ',powerD%initP%infParam
       write(*,*)
       stop       
    endif
    
  end subroutine FreePowers


       



  function ScalarPower(k,in)    
!"in" gives the index of the power to return for this k
!ScalarPower = const for scale invariant spectrum
!The normalization is defined so that for adiabatic perturbations the gradient of the 3-Ricci 
!scalar on co-moving hypersurfaces receives power
! < |D_a R^{(3)}|^2 > = int dk/k 16 k^6/S^6 (1-3K/k^2)^2 ScalarPower(k) 
!In other words ScalarPower is the power spectrum of the conserved curvature perturbation given by
!-chi = \Phi + 2/3*\Omega^{-1} \frac{H^{-1}\Phi' - \Psi}{1+w}
!(w=p/\rho), so < |\chi(x)|^2 > = \int dk/k ScalarPower(k).
!Near the end of inflation chi is equal to 3/2 Psi.
!Here nu^2 = (k^2 + curv)/|curv| 

!This power spectrum is also used for isocurvature modes where 
!< |\Delta(x)|^2 > = \int dk/k ScalarPower(k)
!For the isocurvture velocity mode ScalarPower is the power in the neutrino heat flux.
    use infpert, only : power_spectrum_scal
    use infpert, only : scalNum
    use infpowspline, only : splineval_power_scal
    implicit none

    real(kp), dimension(scalNum,scalNum) :: powerSpectrumScal
    real(dl) :: ScalarPower,k
    integer :: in
   
    if (in.ne.1) stop 'ScalarPower: only 1 Ps allowed'

    if (powerD%initP%useSpline) then
       powerSpectrumScal = splineval_power_scal(k*1._kp)
    else
       powerSpectrumScal = power_spectrum_scal(powerD%infCosmo,k*1._kp)
    endif

    ScalarPower = powerSpectrumScal(scalNum,scalNum)
    
  end function ScalarPower




      
  function TensorPower(k,in)
!TensorPower= const for scale invariant spectrum
!The normalization is defined so that
! < h_{ij}(x) h^{ij}(x) > = \sum_nu nu /(nu^2-1) (nu^2-4)/nu^2 TensorPower(k)
!for a closed model
! < h_{ij}(x) h^{ij}(x) > = int d nu /(nu^2+1) (nu^2+4)/nu^2 TensorPower(k)
!for an open model
!"in" gives the index of the power spectrum to return 
!Here nu^2 = (k^2 + 3*curv)/|curv| 
    use infpert, only : power_spectrum_tens
    use infpowspline, only : splineval_power_tens
    implicit none

    real(dl) :: TensorPower,k   
    real(kp), parameter :: PiByTwo=3.14159265d0/2._kp
   
    integer :: in

    if (in.ne.1) stop 'TensorPower: only 1 Pt allowed'

     if (powerD%initP%useSpline) then
       TensorPower = splineval_power_tens(k*1._kp)
    else
       TensorPower = power_spectrum_tens(powerD%infCosmo,k*1._kp)
    endif
  
    if (curv < 0) TensorPower=TensorPower*tanh(PiByTwo*sqrt(-k**2/curv-3)) 

  end function TensorPower






  function Power_Descript(in, Scal, Tens, Keys, Vals)
!Get parameters describing parameterisation (for FITS file)
    character(LEN=8), intent(out) :: Keys(*)
    real(kp), intent(out) :: Vals(*)
    integer, intent(IN) :: in
    logical, intent(IN) :: Scal, Tens
    integer num, Power_Descript
    num=0
    if ((Scal).or.(Tens)) then       
       num=num+1    
       Keys(num) = 'C1'
       Vals(num) = powerD%initP%infParam%consts(1)
       num=num+1
       Keys(num) = 'C2'
       Vals(num) = powerD%initP%infParam%consts(2)
       num=num+1
       Keys(num) = 'C3'
       Vals(num) = powerD%initP%infParam%consts(3)
       num=num+1
       Keys(num) = 'C4'
       Vals(num) = powerD%initP%infParam%consts(4)
       num=num+1
       Keys(num) = 'C5'
       Vals(num) = powerD%initP%infParam%consts(4)
       num=num+1
!       Keys(num) = 'Conf'
!       Vals(num) = powerD%initP%infParam%conforms(1)
!       num=num+1             
       Keys(num) = 'Field'
       Vals(num) = powerD%initP%infParam%matters(1)    
       num=num+1
       Keys(num) = 'lnReheat'
       Vals(num) = powerD%initP%lnReheat
       num=num+1       
       Keys(num) = 'HubbleEnd'
       Vals(num) = powerD%infCosmo%bgEnd%hubble
       num=num+1
       Keys(num) = 'RhoEndInf'
       Vals(num) = powerD%infCosmo%lnEnergyEnd
       num=num+1
       Keys(num) = 'EfoldEndToToday'
       Vals(num) = powerD%infCosmo%efoldEndToToday
    end if
    Power_Descript = num
    
  end function Power_Descript


  subroutine InitialPower_ReadParams(InitPower, Ini, WantTensors)
    use IniFile
!fields
    use infmatter, only : matterNum
!    use infdilaton, only : dilatonNum
    real :: mu,nu
!end fields

    Type(InitialPowerParams) :: InitPower
    Type(TIniFile) :: Ini
    logical, intent(in) :: WantTensors
    integer i
              
!fields

    InitPower%nn = Ini_Read_Int('initial_power_num',1)
    if (InitPower%nn>nnmax) stop 'Too many initial power spectra - increase nnmax in InitialPower'
    InitPower%infParam%name = Ini_Read_String('inf_model_name')
    InitPower%bgParamNum = Ini_Read_Int('inf_param_number',3)
    
    do i=1, InitPower%bgParamNum
       InitPower%infParam%consts(i) &
            = Ini_Read_Double_Array_File(Ini,'inf_param',i,0._dl)    		  
    end do
    
!    do i=1,dilatonNum
!       InitPower%infParam%conforms(i) &
!            = Ini_Read_Double_Array_File(Ini,'inf_conform',i,1d0)
!    enddo

    do i=1,matterNum
       InitPower%infParam%matters(i) &
            = Ini_Read_Double_Array_File(Ini,'inf_matter',i,1d1)
    enddo
    
    InitPower%checkBound = Ini_Read_Logical('inf_check_bound',.false.)
    InitPower%checkStop = Ini_Read_Logical('inf_check_stop',.false.)
    
    InitPower%lnReheat = Ini_Read_Double('inf_ln_reheat',0d0)
    
    InitPower%useSpline = Ini_Read_Logical('power_spectra_spline',.false.)
    if (InitPower%useSpline) then
       InitPower%lnkmpcMin = Ini_Read_Double('spline_lnkmpc_min',-14d0)
       InitPower%lnkmpcMax = Ini_Read_Double('spline_lnkmpc_max',0d0)
       InitPower%lnkmpcNum = Ini_Read_Int('spline_lnkmpc_num',10)
    endif
    

!convenient rescalings    
!    InitPower%infParam%consts(1) = 10.**InitPower%infParam%consts(1)
!    InitPower%infParam%consts(4) = 10.**InitPower%infParam%consts(4)
    mu = InitPower%infParam%consts(3)
    nu = InitPower%infParam%consts(4)
    if ((mu.ne.0.).and.(nu.eq.0.)) then
       InitPower%infParam%matters = InitPower%infParam%matters &
            *abs(mu)
    endif
    
!end fields


  end  subroutine InitialPower_ReadParams



end module InitialPower
