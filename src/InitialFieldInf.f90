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

!This module provides the initial power spectra for CAMB, computed
!from the inflationary module

module InitialPower
  use Precision
  use MpiUtils, only : MpiStop
  use classes
  use precision, only : dl
  use fieldprec, only : kp 
  use infbgmodel, only : infbgparam
  use infbg, only : infbgphys, infbgdata
  use inftorad, only : inftoradcosmo

  implicit none

  private

  type InitialPowerParams
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
     real(kp) :: lnM
     real(kp) :: lnRrad
     real(kp) :: lnRhoEnd
     real(kp) :: efoldEndToToday
  end type ExportInfProp

  
  logical, parameter :: display=.false.

  logical, parameter, public :: checkStopDefault = .false.
  logical, parameter, public :: checkBoundDefault = .false.
  logical, parameter, public :: useSplineDefault = .true.


  
  Type, extends(TInitialPower) :: TInitialFieldInf

     type(initialpowerparams) :: ipp
         
     real(dl), private :: curv = 0._dl !curvature parameter
     

   contains
         
     procedure :: Init => TInitialFieldInf_Init
     procedure, nopass :: PythonClass => TInitialFieldInf_PythonClass
     procedure, nopass :: SelfPointer => TInitialFieldInf_SelfPointer
     procedure :: ScalarPower => TInitialFieldInf_ScalarPower
     procedure :: TensorPower => TInitialFieldInf_TensorPower
     procedure :: ReadParams => TInitialFieldInf_ReadParams
     procedure :: SetDefPowerParams => TInitialFieldInf_SetDefPowerParams
     procedure :: SetInfBg => TInitialFieldInf_SetInfBg
     procedure :: SetInfCosmo => TInitialFieldInf_SetInfCosmo
     procedure :: SetInfBgSpline => TInitialFieldInf_SetInfBgSpline
     procedure :: FreeInfBgSpline => TInitialFieldInf_FreeInfBgSpline
     procedure :: FreeInfBgData => TInitialFieldInf_FreeInfBgData
     procedure :: SetInfScalPow => TInitialFieldInf_SetInfScalPow
     procedure :: FreePowers => TInitialFieldInf_FreePowers
     procedure :: Effective_ns => TInitalPowerLaw_Effective_ns
     procedure :: Effective_as => TInitalPowerLaw_Effective_as
     procedure :: Effective_bs => TInitalPowerLaw_Effective_bs
  end Type TInitialFieldInf


  Type, extends(TInitialPower) :: TSplinedInitialPower
     real(dl) :: effective_ns_for_nonlinear = -1._dl !used for halofit
     real(dl) :: kmin_scalar, kmax_scalar
     real(dl) :: kmin_tensor, kmax_tensor
     class(TSpline1D), allocatable :: Pscalar, Ptensor
   contains
     procedure :: SetScalarTable => TSplinedInitialPower_SetScalarTable
     procedure :: SetTensorTable => TSplinedInitialPower_SetTensorTable
     procedure :: SetScalarLogRegular => TSplinedInitialPower_SetScalarLogRegular
     procedure :: SetTensorLogRegular => TSplinedInitialPower_SetTensorLogRegular
     procedure :: ScalarPower => TSplinedInitialPower_ScalarPower
     procedure :: TensorPower => TSplinedInitialPower_TensorPower
     procedure :: HasTensors => TSplinedInitialPower_HasTensors
     procedure :: Effective_ns => TSplinedInitialPower_Effective_ns
     procedure, nopass :: PythonClass => TSplinedInitialPower_PythonClass
     procedure, nopass :: SelfPointer => TSplinedInitialPower_SelfPointer
  end Type TSplinedInitialPower

  
!a local static allocated ressource
  type(InitialPowerData), save :: powerD


  
  interface operator (==)     
     module procedure inipowerparams_equal
  end interface operator (==)

  interface operator (/=)
     module procedure inipowerparams_unequal
  end interface operator (/=)

  private operator(==),operator(/=)

  
  public TInitialFieldInf, ExportInfProp


  public UpdateInfProp
  

contains
 
  
  function TInitialFieldInf_PythonClass()
    character(LEN=:), allocatable :: TInitialFieldInf_PythonClass

    TInitialFieldInf_PythonClass = 'InitialPower'

  end function TInitialFieldInf_PythonClass


  subroutine TInitialFieldInf_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TInitialFieldInf), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

  end subroutine TInitialFieldInf_SelfPointer



  function inipowerparams_equal(PinA, PinB)
    use infbgmodel, only : operator(==) 
    implicit none
    type(initialpowerparams), intent(in) :: PinA, PinB
    logical :: inipowerparams_equal

    inipowerparams_equal = ((PinA%infParam == PinB%infParam) &
         .and. (PinA%lnReheat == PinB%lnReheat) &
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

    export%lnM = log(powerD%initP%infParam%consts(1))
    export%lnRrad = powerD%initP%lnReheat &
         - 0.25*powerD%infCosmo%lnEnergyEnd
    export%lnRhoEnd = powerD%infCosmo%lnEnergyEnd
    export%efoldEndToToday = powerD%infCosmo%efoldEndToToday
   
  end subroutine UpdateInfProp

  
  
  subroutine TInitialFieldInf_SetDefPowerParams(this)
    implicit none
    class(TInitialFieldInf) :: this

    call SetDefPowerParams(this%ipp)

  end subroutine TInitialFieldInf_SetDefPowerParams
    

  subroutine SetDefPowerParams(Pin)
    use infbgmodel, only : fieldNum
    use infdilaton, only : dilatonNum
    use infbgmodel, only : infParamNum
    implicit none
    type (InitialPowerParams), intent(out) :: Pin
      
    Pin%bgParamNum = infParamNum

!stop inflation according to field values
    Pin%checkStop = checkStopDefault

!impose hard prior from theoretical bounds on field values
    Pin%checkBound = checkBoundDefault
   
!range and sampling of the power spectra spline
    Pin%useSpline = useSplineDefault 
    Pin%lnkmpcMin = -14._kp
    Pin%lnkmpcMax = 1._kp
    Pin%lnkmpcNum = 13
      
    Pin%lnReheat = 0._kp
    Pin%kstar = 0.05_kp

!value for the parameters (see infbg.f90)
    Pin%infParam%name = 'largef'
    Pin%infParam%consts(1:infParamNum) = 0._kp

    if (dilatonNum.ne.0) then
       stop 'SetDefPowerParams: InitialFieldInf.f90 needs edition!'
!       Pin%infParam%conforms(1:dilatonNum) = 1.
    endif

    Pin%infParam%consts(1) = 1d-5
    Pin%infParam%consts(2) = 2._kp
    Pin%infParam%consts(3:4) = 0._kp
    Pin%infParam%consts(5) = 1._kp

    if (infParamNum.gt.5) then
       Pin%infParam%consts(6:infParamNum) = 0._kp
    endif

  end subroutine SetDefPowerParams




  
  subroutine TInitialFieldInf_Init(this, Params)
    use results
    use constants, only : c
    class(TInitialFieldInf) :: this
    class(TCAMBParameters), intent(in) :: Params

    integer :: inferror
    logical, parameter :: usePstar=.false.
    real(kp), parameter :: Pstar = 2e-9
    
    
    select type(Params)
    class is (CAMBParams)
       !Curvature parameter if non-flat
       this%curv = -Params%Omk/((c/1000)/Params%H0)**2
       
    class default
       call Mpistop( 'TInitialFieldInf_Init: class type not found!' )
       
    end select


    inferror = 0

    if (this%curv.ne.0._dl) call mpistop ('TInitialFieldInf_Init: inflation in flat universe only!')

!one background for all k (that's why nnmax=1)
!    print *,'we are in initialize power',(Pin%infParam == powerD%initP%infParam) &
!         ,(Pin%infParam /= powerD%initP%infParam)   
  
    call this%SetInfBg(inferror)
    
    if (inferror.ne.0) call mpistop( 'TInitialFieldInf_Init: unproper infbg')

    if (usePstar) then
!this is only for playing with normalised power spectra to Pstar. 
       write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       write(*,*)'TInitialFieldInf_Init: renormalising P(k*) to Pstar=',Pstar
       write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       read(*,*)
       call this%SetInfScalPow(Pstar)
    endif

    call this%SetInfBgSpline()

    call this%SetInfCosmo(inferror)
    if (inferror.ne.0) call mpistop ('TInitialFieldInf_Init: unproper inftorad')

  end subroutine TInitialFieldInf_Init


  subroutine TInitialFieldInf_SetInfBg(this, inferror)
    implicit none
    class(TInitialFieldInf) :: this
    integer, optional :: inferror

    call SetInfBg(this%ipp,inferror)

  end subroutine TInitialFieldInf_SetInfBg



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
       call mpistop('SetInfBg: params initialisation failed!')
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
  


  subroutine TInitialFieldInf_SetInfBgSpline(this)
    implicit none
    class(TInitialFieldInf) :: this

    call SetInfBgSpline(this%ipp)

  end subroutine TInitialFieldInf_SetInfBgSpline

  

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
       call mpistop('SetInfBgSpline: Pin%params >< Data params')
    endif

  end subroutine SetInfBgSpline


  subroutine TInitialFieldInf_FreeInfBgSpline(this)
    implicit none
    class(TInitialFieldInf) :: this

    call FreeInfBgSpline(this%ipp)

  end subroutine TInitialFieldInf_FreeInfBgSpline


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
       call mpistop('FreeInfBgSpline: Pin%params >< Data params')
    endif

  end subroutine FreeInfBgSpline


  subroutine TInitialFieldInf_FreeInfBgData(this)
    implicit none
    class(TInitialFieldInf) :: this

    call FreeInfBgData(this%ipp)

  end subroutine TInitialFieldInf_FreeInfBgData



  subroutine FreeInfBgData(Pin)
    use infbgmodel, only : operator(==)
    use infbg, only : free_infbg_data
    implicit none

    type (initialpowerparams), intent(in) :: Pin    
   
    if (.not.associated(powerD%ptrBgdata)) call mpistop('FreeInfBgData: no bgdata found!')

!sanity checks
    if (Pin%infParam == powerD%initP%infParam) then
       call free_infbg_data(powerD%ptrBgdata)      
    else
       call mpistop('FreeInfBgData: Pin%params >< Data params')
    endif

  end subroutine FreeInfBgData
  

  subroutine TInitialFieldInf_SetInfCosmo(this, inferror)
    implicit none
    class(TInitialFieldInf) :: this
    integer, optional :: inferror

    call SetInfCosmo(this%ipp,inferror)

  end subroutine TInitialFieldInf_SetInfCosmo


  subroutine SetInfCosmo(Pin,inferror)
    use infbgmodel, only : operator(/=)   
    use inftorad, only : set_inftorad_cosmo
    use infpowspline, only : check_power_scal_spline, check_power_tens_spline
    
    implicit none

    type (initialpowerparams), intent(in) :: Pin        
    integer, optional :: inferror

!sanity checks
    if (Pin%infParam /= powerD%initP%infParam) then
       call mpistop('SetInfCosmo: Pin%params >< Data params')
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
       call mpistop('SetInfPowSpline: no InfBg found, call SetInfBg before!')
    else
       if (Pin /= powerD%initP) then
          call mpistop('SetInfPowSpline: Pin >< Data')
       endif
    endif
    if (.not.Pin%useSpline) stop 'SetInfPowSpline: useSpline is F'
         
    do i=1,powerD%initP%lnkmpcNum
       lnkmpcVec(i) = powerD%initP%lnkmpcMin + real(i-1,kp)*(powerD%initP%lnkmpcMax &
            - powerD%initP%lnkmpcMin)/real(powerD%initP%lnkmpcNum-1,kp)
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
  


  subroutine TInitialFieldInf_SetInfScalPow(this, Pwanted)
    implicit none
    class(TInitialFieldInf) :: this
    real(kp), intent(in) :: Pwanted
    
    call SetInfScalPow(this%ipp,Pwanted)

  end subroutine TInitialFieldInf_SetInfScalPow


  

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


  subroutine TInitialFieldInf_FreePowers(this)
    implicit none
    class(TInitialFieldInf) :: this

    call FreePowers(this%ipp)

  end subroutine TInitialFieldInf_FreePowers

  

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



  


  
  function TInitialFieldInf_ScalarPower(this, k)
 
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

    class(TInitialFieldInf) :: this
    real(dl), intent(in) :: k
    real(dl) TInitialFieldInf_ScalarPower

    real(kp), dimension(scalNum,scalNum) :: powerSpectrumScal


    if (powerD%initP%useSpline) then
       powerSpectrumScal = splineval_power_scal(real(k,kp))
    else
       powerSpectrumScal = power_spectrum_scal(powerD%infCosmo,real(k,kp))
    endif

    TInitialFieldInf_ScalarPower = powerSpectrumScal(scalNum,scalNum)

  end function TInitialFieldInf_ScalarPower

  

  

  function TInitialFieldInf_TensorPower(this,k)

    !TensorPower= const for scale invariant spectrum
    !The normalization is defined so that
    ! < h_{ij}(x) h^{ij}(x) > = \sum_nu nu /(nu^2-1) (nu^2-4)/nu^2 TensorPower(k)
    !for a closed model
    ! < h_{ij}(x) h^{ij}(x) > = int d nu /(nu^2+1) (nu^2+4)/nu^2 TensorPower(k)
    !for an open model
    !Here nu^2 = (k^2 + 3*curv)/|curv|

    use infpert, only : power_spectrum_tens
    use infpowspline, only : splineval_power_tens
    use constants
    implicit none

    class(TInitialFieldInf) :: this
    real(dl), intent(in) :: k
    real(dl) TInitialFieldInf_TensorPower

    
    if (powerD%initP%useSpline) then
       TInitialFieldInf_TensorPower = splineval_power_tens(real(k,kp))
    else
       TInitialFieldInf_TensorPower = power_spectrum_tens(powerD%infCosmo,real(k,kp))
    endif

  end function TInitialFieldInf_TensorPower




  
  function CompatKey(Ini, name)
    class(TIniFile), intent(in) :: Ini
    character(LEN=*), intent(in) :: name
    character(LEN=:), allocatable :: CompatKey
    !Allow backwards compatibility with old .ini files where initial power parameters were arrays

    if (Ini%HasKey(name//'(1)')) then
       CompatKey = name//'(1)'
       if (Ini%HasKey(name)) call MpiStop('Must have one of '//trim(name)//' or '//trim(name)//'(1)')
    else
       CompatKey = name
    end if
  end function CompatKey



  
  subroutine TInitialFieldInf_ReadParams(this, Ini)
    use IniObjects
    use infmatter, only : matterNum

    class(TInitialFieldInf) :: this
    class(TIniFile), intent(in) :: Ini
    logical :: WantTensors
    
    real :: mu,nu
    integer :: i

    WantTensors = Ini%Read_Logical('get_tensor_cls', .true.)

    this%ipp%infParam%name = Ini%Read_String('inf_model_name')
    this%ipp%bgParamNum = Ini%Read_Int('inf_param_number',3)

    do i=1,this%ipp%bgParamNum
       this%ipp%infParam%consts(i) = Ini%Read_Double_Array('inf_param',i,0._dl)
    end do

!    do i=1,dilatonNum
!       this%ipp%infParam%conforms(i) &
!            = Ini%Read_Double_Array('inf_conform',i,1._dl)
!    enddo
    
    do i=1,matterNum
       this%ipp%infParam%matters(i) &
            = Ini%Read_Double_Array('inf_matter',i,1d1)
    enddo


    this%ipp%checkBound = Ini%Read_Logical('inf_check_bound',.false.)
    this%ipp%checkStop = Ini%Read_Logical('inf_check_stop',.false.)
    
    this%ipp%lnReheat = Ini%Read_Double('inf_ln_reheat',0._dl)
    
    this%ipp%useSpline = Ini%Read_Logical('power_spectra_spline',.false.)
    if ( this%ipp%useSpline) then
       this%ipp%lnkmpcMin = Ini%Read_Double('spline_lnkmpc_min',-14._dl)
       this%ipp%lnkmpcMax = Ini%Read_Double('spline_lnkmpc_max',0._dl)
       this%ipp%lnkmpcNum = Ini%Read_Int('spline_lnkmpc_num',10)
    endif

!convenient rescalings    
!    this%ipp%infParam%consts(1) = 10.**this%ipp%infParam%consts(1)
!    this%ipp%infParam%consts(4) = 10.**this%ipp%infParam%consts(4)
    mu = this%ipp%infParam%consts(3)
    nu = this%ipp%infParam%consts(4)
    if ((mu.ne.0.).and.(nu.eq.0.)) then
       this%ipp%infParam%matters = this%ipp%infParam%matters &
            *abs(mu)
    endif    


  end subroutine TInitialFieldInf_ReadParams



  function TInitalFieldInf_Effective_ns(this)
    use infpowspline, only : splineval_ns_scal
    class(TInitialFieldInf) :: this
    real(dl) :: TInitalPowerLaw_Effective_ns

!only possible if spline is on (safest way to compute numerical derivatives)    
    if (.not.(this%ipp%useSpline)) then
       write(*,*)'useSpline= ',this%ipp%useSpline
       call Mpistop( 'TInitalFieldInf_Effective_ns: enable spline to compute ns!' )
    end if
       
    TInitialFieldInf_Effective_ns = splineval_ns_scal(real(this%ipp%kstar,kp))   

  end function TInitalPowerLaw_Effective_ns

!same as before for the running and running of the running
  function TInitalFieldInf_Effective_as(this)
    use infpowspline, only : splineval_alphas_scal
    class(TInitialFieldInf) :: this
    real(dl) :: TInitalPowerLaw_Effective_as

    if (.not.(this%ipp%useSpline)) then
       write(*,*)'useSpline= ',this%ipp%useSpline
       call Mpistop( 'TInitalFieldInf_Effective_ns: enable spline to compute as!' )
    end if
       
    TInitialFieldInf_Effective_as = splineval_alphas_scal(real(this%ipp%kstar,kp))   

  end function TInitalFieldInf_Effective_as

  function TInitalFieldInf_Effective_bs(this)
    use infpowspline, only : splineval_betas_scal
    class(TInitialFieldInf) :: this
    real(dl) :: TInitalPowerLaw_Effective_bs

    if (.not.(this%ipp%useSpline)) then
       write(*,*)'useSpline= ',this%ipp%useSpline
       call Mpistop( 'TInitalFieldInf_Effective_as: enable spline to compute bs!' )
    end if
       
    TInitialFieldInf_Effective_bs = splineval_betas_scal(real(this%ipp%kstar,kp))   

  end function TInitalFieldInf_Effective_bs
  

  
  subroutine TSplinedInitialPower_SelfPointer(cptr, P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TSplinedInitialPower), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

  end subroutine TSplinedInitialPower_SelfPointer

  logical function TSplinedInitialPower_HasTensors(this)
    class(TSplinedInitialPower) :: this

    TSplinedInitialPower_HasTensors = allocated(this%Ptensor)

  end function TSplinedInitialPower_HasTensors

  function TSplinedInitialPower_ScalarPower(this, k)
    class(TSplinedInitialPower) :: this
    real(dl), intent(in) ::k
    real(dl) TSplinedInitialPower_ScalarPower

    if (k <= this%kmin_scalar) then
       TSplinedInitialPower_ScalarPower = this%Pscalar%F(1)
    elseif (k >= this%kmax_scalar) then
       TSplinedInitialPower_ScalarPower = this%Pscalar%F(this%Pscalar%n)
    else
       TSplinedInitialPower_ScalarPower = this%Pscalar%Value(k)
    end if

  end function TSplinedInitialPower_ScalarPower

  function TSplinedInitialPower_TensorPower(this, k)
    class(TSplinedInitialPower) :: this
    real(dl), intent(in) ::k
    real(dl) TSplinedInitialPower_TensorPower

    if (k <= this%kmin_tensor) then
       TSplinedInitialPower_TensorPower = this%Ptensor%F(1)
    elseif (k >= this%kmax_tensor) then
       TSplinedInitialPower_TensorPower = this%Ptensor%F(this%Ptensor%n)
    else
       TSplinedInitialPower_TensorPower = this%Ptensor%Value(k)
    end if

  end function TSplinedInitialPower_TensorPower

  function TSplinedInitialPower_PythonClass()
    character(LEN=:), allocatable :: TSplinedInitialPower_PythonClass

    TSplinedInitialPower_PythonClass = 'SplinedInitialPower'

  end function TSplinedInitialPower_PythonClass

  subroutine TSplinedInitialPower_SetScalarTable(this, n, k, PK)
    class(TSplinedInitialPower) :: this
    integer, intent(in) :: n
    real(dl), intent(in) :: k(n), PK(n)

    if (allocated(this%Pscalar)) deallocate(this%Pscalar)
    if (n>0) then
       allocate(TCubicSpline::this%Pscalar)
       select type (Sp => this%Pscalar)
       class is (TCubicSpline)
          call Sp%Init(k,PK)
       end select
       this%kmin_scalar = k(1)
       this%kmax_scalar = k(n)
    end if

  end subroutine TSplinedInitialPower_SetScalarTable


  subroutine TSplinedInitialPower_SetTensorTable(this, n, k, PK)
    class(TSplinedInitialPower) :: this
    integer, intent(in) :: n
    real(dl), intent(in) :: k(n), PK(n)

    if (allocated(this%PTensor)) deallocate(this%PTensor)
    if (n>0) then
       allocate(TCubicSpline::this%PTensor)
       select type (Sp => this%PTensor)
       class is (TCubicSpline)
          call Sp%Init(k,PK)
       end select
       this%kmin_tensor = k(1)
       this%kmax_tensor = k(n)
    end if

  end subroutine TSplinedInitialPower_SetTensorTable

  subroutine TSplinedInitialPower_SetScalarLogRegular(this, kmin, kmax, n, PK)
    class(TSplinedInitialPower) :: this
    integer, intent(in) :: n
    real(dl), intent(in) ::kmin, kmax, PK(n)

    if (allocated(this%Pscalar)) deallocate(this%Pscalar)
    if (n>0) then
       allocate(TLogRegularCubicSpline::this%Pscalar)
       select type (Sp => this%Pscalar)
       class is (TLogRegularCubicSpline)
          call Sp%Init(kmin, kmax, n, PK)
       end select
       this%kmin_scalar = kmin
       this%kmax_scalar = kmax
    end if

  end subroutine TSplinedInitialPower_SetScalarLogRegular


  subroutine TSplinedInitialPower_SetTensorLogRegular(this, kmin, kmax, n, PK)
    class(TSplinedInitialPower) :: this
    integer, intent(in) :: n
    real(dl), intent(in) ::kmin, kmax, PK(n)

    if (allocated(this%Ptensor)) deallocate(this%Ptensor)
    if (n>0) then
       allocate(TLogRegularCubicSpline::this%Ptensor)
       select type (Sp => this%Ptensor)
       class is (TLogRegularCubicSpline)
          call Sp%Init(kmin, kmax, n, PK)
       end select
       this%kmin_tensor = kmin
       this%kmax_tensor = kmax
    end if

  end subroutine TSplinedInitialPower_SetTensorLogRegular

  function TSplinedInitialPower_Effective_ns(this)
    use config
    class(TSplinedInitialPower) :: this
    real(dl) :: TSplinedInitialPower_Effective_ns

    if (this%effective_ns_for_nonlinear==-1._dl) then
       call GlobalError('TSplinedInitialPower: effective_ns_for_nonlinear not set',error_inital_power)
    else
       TSplinedInitialPower_Effective_ns = this%effective_ns_for_nonlinear
    end if
  end function TSplinedInitialPower_Effective_ns



    
end module InitialPower
