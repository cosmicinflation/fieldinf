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






module infbg
  use fieldprec, only : kp, tolkp
  use infbgmodel, only : matterNum, dilatonNum, fieldNum
  use infbgfunc, only : field_second_derivative, velocity_derivative
  use infbgfunc, only : hubble_parameter_square, slowroll_first_parameter_JF
  use infbgfunc, only : slowroll_first_parameter, slowroll_second_parameter
  use infbgfunc, only : slowroll_third_parameter

!background evolution in the Einstein FLRW Frame for multifields:
!scalar gravity + matter fields.

  implicit none

  private

  abstract interface
     function sr_param(x,dx,switch)
       use fieldprec, only : kp
       use infbgmodel, only : fieldNum
       real(kp) :: sr_param
       real(kp), dimension(fieldNum), intent(in) :: x, dx
       logical, intent(in) :: switch
     end function sr_param
  end interface

  
   
!for debugging
  logical, parameter :: display = .false.
  logical, parameter :: dump_file = .false.

!default value for integration settings, can be accessed with some
!dedicated public function below
  logical, save :: BgEvolCheckForMatterStop = .false.
  real(kp), save :: BgEvolMatterStopValue = tolkp
  logical, save :: BgEvolStopAtMax = .false.
  
  logical, save :: BgEvolCheckForHubbleStop = .false.
  real(kp), save :: BgEvolHubbleStopValue = 0._kp
  integer, save :: BgEvolStopVelocitySign = 0

!use other slow-roll parameter than eps1
  logical, save :: BgEvolUseOtherEpsilon = .false.
  real(kp), save :: BgEvolEpsilonStopValue = 1._kp
  procedure(sr_param), pointer :: ptr_alternate_stop_parameter => slowroll_second_parameter

!default maximum number of efolds  
  real(kp), save :: BgEvolEfoldMaxiStop = 1000._kp

!default number of e-folds to store after stopping criterion
  real(kp), save :: BgEvolEfoldExploreOsc = 0._kp
  

  
!to store snapshot (ini or end, or more)
  type infbgphys    
     sequence
     real(kp) :: efold, hubble
     real(kp) :: epsilon1, epsilon1JF
     real(kp) :: epsilon2, epsilon3
     real(kp), dimension(fieldNum) :: field
     real(kp), dimension(fieldNum) :: fieldDot
  end type infbgphys


!to store the bg integration as chained list
  type infbgdata
     type(infbgphys) :: bg
     type(infbgdata), pointer :: ptr => null()
  end type infbgdata
  

  interface operator (==)     
     module procedure infbgphys_equal
  end interface


  interface operator (/=)
     module procedure infbgphys_unequal
  end interface



  public infbgdata, infbgphys
  public operator(==),operator(/=)
  public free_infbg_data, count_infbg_data
  public print_infbgphys

  public set_infbg_ini
  public rescale_potential
  public bg_field_evol, bg_field_dot_coupled

  public set_bgfieldevol_matterstop, set_bgfieldevol_efoldexploreosc
  public set_bgfieldevol_hubblestop, set_bgfieldevol_epsilonstop
  public set_bgfieldevol_useotherepsilon, set_bgfieldevol_efoldmaxistop
  

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!recursivity need enough stacksize for big lists, otherwise it
!segfaults. Only needed to store and free the data in memory.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  function infbgphys_equal(infbgphysA, infbgphysB)
    implicit none
    type(infbgphys), intent(in) :: infbgphysA, infbgphysB
    logical :: infbgphys_equal

    infbgphys_equal = ((infbgphysA%efold == infbgphysB%efold) &
         .and. (infbgphysA%hubble == infbgphysB%hubble) &
         .and. (infbgphysA%epsilon1 == infbgphysB%epsilon1) &
         .and. (infbgphysA%epsilon1JF == infbgphysB%epsilon1JF) &
         .and. (infbgphysA%epsilon2 == infbgphysB%epsilon2) &
!         .and. (infbgphysA%epsilon3 == infbgphysB%epsilon3) &
         .and. all(infbgphysA%field == infbgphysB%field) &
         .and. all(infbgphysA%fieldDot == infbgphysB%fieldDot))

  end function infbgphys_equal



  function infbgphys_unequal(infbgphysA, infbgphysB)
    implicit none
    type(infbgphys), intent(in) :: infbgphysA, infbgphysB
    logical :: infbgphys_unequal

    infbgphys_unequal = ((infbgphysA%efold /= infbgphysB%efold) &
         .or. (infbgphysA%hubble /= infbgphysB%hubble) &
         .or. (infbgphysA%epsilon1 /= infbgphysB%epsilon1) &
         .or. (infbgphysA%epsilon1JF /= infbgphysB%epsilon1JF) &
         .or. (infbgphysA%epsilon2 /= infbgphysB%epsilon2) &
!         .or. (infbgphysA%epsilon3 /= infbgphysB%epsilon3) &
         .or. any(infbgphysA%field /= infbgphysB%field) &
         .or. any(infbgphysA%fieldDot /= infbgphysB%fieldDot))

  end function infbgphys_unequal



  recursive subroutine free_infbg_data(ptrFirst)
    implicit none
    type(infbgdata), pointer :: ptrFirst

    if (associated(ptrFirst%ptr)) call free_infbg_data(ptrFirst%ptr)
    deallocate(ptrFirst)

  end subroutine free_infbg_data



  recursive function count_infbg_data(ptrFirst) result(bgdataCount)
    implicit none    
    integer :: bgdataCount
    type(infbgdata), pointer :: ptrFirst

    bgdataCount = 1
    if (associated(ptrFirst%ptr)) then
       bgdataCount =  count_infbg_data(ptrFirst%ptr) + 1       
    endif

  end function count_infbg_data


  subroutine print_infbgphys(infvar,inname)
    implicit none
    type(infbgphys), intent(in) :: infvar
    character(len=*), intent(in), optional :: inname

    if (present(inname)) then
       write(*,*)'infbgphys type: ',inname
    endif
    write(*,*)'efold=          ',infvar%efold
    write(*,*)'hubble=         ',infvar%hubble
    write(*,*)'epsilon1=       ',infvar%epsilon1
    write(*,*)'epsilon1(JF)=   ',infvar%epsilon1JF
    write(*,*)'epsilon2=       ',infvar%epsilon2
    write(*,*)'epsilon3=       ',infvar%epsilon3
    write(*,*)'field value=    ',infvar%field
    write(*,*)'Dfield/Defold=  ',infvar%fieldDot


  end subroutine print_infbgphys

  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!inflation settings: initial conditions, rescaling, normalisation...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  function set_infbg_ini(infParam,fieldDot)
!to start on the attractor. This is epsilon2 = 0 for epsilon1<<1
!initialy.
    use infbgmodel, only : infbgparam
    use infsigma, only : metric_inverse
    use infpotential, only : deriv_ln_potential
    use infsric, only : slowroll_initial_matter

    implicit none
    type(infbgparam), intent(in) :: infParam
    real(kp), dimension(fieldNum), intent(in), optional :: fieldDot
    type(infbgphys) :: set_infbg_ini
    
    real(kp) :: hubbleSquareIni
    type(infbgphys) :: infIni


    infIni%field(1:matterNum) = infParam%matters
    infIni%field(matterNum+1:fieldNum) = infParam%conforms


!if the matter fields are set to 0, use slow-roll guesses
    if (all(infParam%matters == 0._kp)) then
       infIni%field(1:matterNum) = slowroll_initial_matter(infParam)
    endif

!if on input "fieldDot" is given then fix the initial velocities
!accordingly, otherwise, starts with initial velocities on the
!slow-roll attractor
    if (present(fieldDot)) then
       infIni%fieldDot = fieldDot
    else
       infIni%fieldDot &
            = - matmul(metric_inverse(infIni%field),deriv_ln_potential(infIni%field))
    endif

    hubbleSquareIni = hubble_parameter_square(infIni%field,infIni%fieldDot,.false.)

    if (hubbleSquareIni.ge.0._kp) then
       infIni%hubble = sqrt(hubbleSquareIni)       
    else     
       write(*,*)'Field= FieldDot= ',infIni%field,infIni%fieldDot
       write(*,*)'H^2 = ',hubbleSquareIni
       stop 'H^2 < 0, check initial condition'
    endif

    infIni%epsilon1 = slowroll_first_parameter(infIni%field,infIni%fieldDot, .false.)
    infIni%epsilon2 = slowroll_second_parameter(infIni%field,infIni%fieldDot, .false.)
    infIni%epsilon3 = slowroll_third_parameter(infIni%field,infIni%fieldDot, .false.)
    infIni%epsilon1JF = slowroll_first_parameter_JF(infIni%field, infIni%fieldDot,.false.)
    infIni%efold = 0._kp

    set_infbg_ini = infIni
   
    if (display) then
       if (infIni%epsilon1.lt.epsilon(1._kp)) then
          write(*,*)
          write(*,*)'set_infbg_ini: epsilon1 < accuracy',infIni%epsilon1
          write(*,*)
       endif
    endif

  end function  set_infbg_ini





  subroutine rescale_potential(scale,infParam,infIni,infEnd,infObs,ptrBgdata)
!update all relevant data such as Unew = scale*Uold
    use infbgmodel, only : infbgparam, set_infbg_param
    implicit none
    type(infbgparam), intent(inout) :: infParam
    type(infbgphys), intent(inout) :: infIni,infEnd,infObs
    type(infbgdata), optional, pointer :: ptrBgdata
    real(kp), intent(in) :: scale

    type(infbgdata), pointer :: ptrRun    
    logical :: updateBgParams

    ptrRun => null()

!see infbgmodel U propto M^4
    infParam%consts(1) = infParam%consts(1)*scale**0.25_kp

    updateBgParams = set_infbg_param(infParam)
    if (.not.updateBgParams) stop 'rescale_potential: updating params failed!'
    
    infIni%hubble = infIni%hubble * sqrt(scale)
    infEnd%hubble = infEnd%hubble * sqrt(scale)
    infObs%hubble = infObs%hubble * sqrt(scale)

    if (present(ptrBgData)) then
       if (associated(ptrBgdata)) then
          ptrRun => ptrBgdata
          do while (associated(ptrRun))
             ptrRun%bg%hubble = ptrRun%bg%hubble * sqrt(scale)
             ptrRun => ptrRun%ptr
          enddo
          ptrRun => null()
       else
          stop 'rescale_potential_by: data not found'
       endif
    endif

  end subroutine rescale_potential


 
  

!convenience functions to allow user alteration of various integrator
!criteria

  
  subroutine set_bgfieldevol_matterstop(switch, valuestop, stopatmax)
    implicit none
    logical, intent(in) :: switch
    real(kp), intent(in) :: valuestop
    logical, intent(in) :: stopatmax

    BgEvolCheckForMatterStop = switch
    BgEvolMatterStopValue = valuestop
    BgEvolStopAtMax = stopatmax
    
    write(*,*)'infbg: default checkMatterStop sets to: ',switch
    write(*,*)'infbg: default stopMatterValue sets to: ',valuestop
    write(*,*)'infbg: default stopAtMax sets to:       ',stopatmax
       
  end subroutine set_bgfieldevol_matterstop


  
  subroutine set_bgfieldevol_hubblestop(switch, valuestop, stopsign)
    implicit none
    logical, intent(in) :: switch
    real(kp), intent(in) :: valuestop
    integer, intent(in), optional :: stopsign

    BgEvolCheckForHubbleStop = switch
    BgEvolHubbleStopValue = valuestop
    if (present(stopsign)) BgEvolStopVelocitySign = stopsign
    
    write(*,*)'infbg: checkHubbleStop sets to: ',switch
    write(*,*)'infbg: HubbleStopValue sets to: ',valuestop
    if (present(stopsign)) write(*,*)'infbg: stopVelocitySign sets to: ',stopsign
       
  end subroutine set_bgfieldevol_hubblestop
 

  

  subroutine set_bgfieldevol_useotherepsilon(switch,epsname)
    implicit none
    logical, intent(in) :: switch
    character(len=*), intent(in), optional :: epsname
    
    BgEvolUseOtherEpsilon = switch

    write(*,*)'infbg: useOtherEpsilon sets to: ',switch

    if (present(epsname)) then
       select case(epsname)

       case('epsilon1JF')
          ptr_alternate_stop_parameter => slowroll_first_parameter_JF
       case('epsilon2')
          ptr_alternate_stop_parameter => slowroll_second_parameter
       case('epsilon3')
          ptr_alternate_stop_parameter => slowroll_third_parameter
       case default
          stop 'set_bgfieldevol_useotherepsilon: name not found!'

       end select
    end if
    
  end subroutine set_bgfieldevol_useotherepsilon



  
  subroutine set_bgfieldevol_epsilonstop(eps)
    implicit none
    real(kp), intent(in) :: eps

    BgEvolEpsilonStopValue = eps

    write(*,*)'infbg: setting EpsilonStop= ',eps

  end subroutine set_bgfieldevol_epsilonstop




  subroutine set_bgfieldevol_efoldmaxistop(efold)
    implicit none
    real(kp), intent(in) :: efold

    BgEvolEfoldMaxiStop = efold

    write(*,*)'infbg: setting EfoldMaxiStop= ',efold

  end subroutine set_bgfieldevol_efoldmaxistop



  subroutine set_bgfieldevol_efoldexploreosc(efold)
    implicit none
    real(kp), intent(in) :: efold

    BgEvolEfoldExploreOsc = efold

    write(*,*)'infbg: setting EfoldExploreOsc= ',efold

  end subroutine set_bgfieldevol_efoldexploreosc
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!inflationary evolution: find end of inflation + store relevant quantities
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  function bg_field_evol(infIni,efoldDataNum,infObs,ptrStart,stopAtThisValue,isStopAtMax &
       ,hubbleStopAtThisValue,whichStopVelocitySign) 
!integrate the background until epsilon > epsilonStop, and returns some
!physical values (type infbgphys) for which epsilon = epsilon1EndInf  (=1)

    use fieldprec, only : transfert
    use infsolvers, only : easydverk, tunedverk, zbrent
    use infdilaton, only : conformal_factor_square
    use infpotential, only : potential
    use infio
    
    implicit none
    
    type(infbgphys), intent(in) :: infIni
!number of wanted stored efold. If accuracy is not enough, the real
!number of stored efold is modified !up to
!efoldBeforeEndObs/efoldStepDefault
    integer, optional, intent(in) :: efoldDataNum
    type(infbgphys), optional, intent(out) :: infObs
    type(infbgdata), optional, pointer :: ptrStart    

    real(kp), optional, intent(in) :: stopAtThisValue
    logical, optional, intent(in) :: isStopAtMax
    real(kp), optional, intent(in) :: hubbleStopAtThisValue
    integer, optional, intent(in) :: whichStopVelocitySign
    type(infbgphys) :: bg_field_evol
    
    
    real(kp) :: epsilon1, epsilon1JF
    real(kp) :: epsilon2, epsilon3

!if ptrBgdata input without data number, this is the default storage step
    real(kp), parameter :: efoldStepDefault = 1._kp

!we cannot discover inflation longer on this computer
    real(kp) :: efoldHuge

!how many efold after end inflation are stored (work only with
!useVelocity=T)
    real(kp) :: efoldExploreOsc
    
!observable perturbations were produced after that
    real(kp), parameter :: efoldBeforeEndObs = 120._kp
    real(kp) :: efoldObs

    real(kp) :: hubbleSquare, hubble, hubbleEndInf
    real(kp) :: hubbleSquareIni, hubbleIni
    
    real(kp) :: efold,efoldNext,efoldStepObs

    real(kp) :: efoldBeforeEndInf,efoldEndInf
    real(kp) :: efoldAfterEndInf

    real(kp), dimension(fieldNum) :: field,derivField
    real(kp), dimension(fieldNum) :: fieldEndInf,derivFieldEndInf
    real(kp), dimension(fieldNum) :: derivFieldAfterEndInf
    real(kp), dimension(2*fieldNum) :: bgVar, bgVarIni

!standard integration accuracy 
    real(kp), parameter :: tolEvol = tolkp

!backward integration accuracy (sometime instable, mayneed extra precision)
    real(kp) :: tolBackEvol

!if true derivField=Dfield/Dtphys, otherwise derivField=Dfield/Defold is
!used for the integration. In both cases, only Dfield/Defold is
!stored
    logical, parameter :: useVelocity = .true.

!Value of "epsilon1" at which inflation stops. Default to BgEvolEpsilonStopValue=1
    real(kp) :: epsilonStop
!By default, we use epsilon1 (EF), but for non-minimal coupling, you
!want epsilon1 (JF). Default to BgEvolUseOtherEpsilon=F
    logical :: useOtherEpsilon

!zbrent accuracy on efoldEnd for which epsilon=epsilonStop
    real(kp), parameter :: tolEfoldEnd = tolkp

!another tests, checked *after* epsilon1 values to stop
!integration. May be convenient for hybrid-like one field potential:
!stop integration when the min(max) value of the matter fields is
!below of above a given value, or according to the total number of
!efolds, or to the value of the hubble parameter
    real(kp) :: efoldMaxiStop
    logical :: checkMatterStop
    logical :: stopAtMax

    logical :: checkHubbleStop
    integer :: stopVelocitySign

!whether or not to determine accurately the end of inflation
    logical, parameter :: accurateEndInf = .true.
    
    real(kp) :: matterStopValue, hubbleStopValue
    integer, parameter :: stopIndexMin = 1, stopIndexMax = 1

    logical :: inflate, longEnoughObs
  
!to make old f77 routines discussing together           
    type(transfert) :: stopData, findData

!to store the data as chained list    
    type(infbgdata), pointer :: ptrCurrent
    type(infbgdata), pointer :: ptrPrevious
       
    integer :: neqs


!initialisation   !!!!!!!!!!!!!!!
    tolBackEvol = tolEvol
    efoldHuge = 1._kp/epsilon(1._kp)
    neqs = 2*fieldNum

    if (present(infObs)) then
       infObs = infIni
    endif

    efoldExploreOsc = BgEvolEfoldExploreOsc

!initialization of inflation ending check and conditional stops
    epsilonStop = BgEvolEpsilonStopValue
    useOtherEpsilon = BgEvolUseOtherEpsilon
    efoldMaxiStop = BgEvolEfoldMaxiStop
    
    if (present(stopAtThisValue)) then
       checkMatterStop = .true.
       matterStopValue = stopAtThisValue
       if (display) then
          write(*,*)'bg_field_evol: check for Matter Stop enabled'
          write(*,*)'matterStopValue= ',matterStopValue
       endif
    else
       checkMatterStop = BgEvolCheckForMatterStop
       matterStopValue = BgEvolMatterStopValue
    endif

    if (present(isStopAtMax)) then
       stopAtMax = isStopAtMax
       if (display) write(*,*)'bg_field_evol: stopAtMax is',stopAtMax
    else
       stopAtMax = BgEvolStopAtMax
    endif

    if (present(hubbleStopAtThisValue)) then
       checkHubbleStop = .true.
       hubbleStopValue = hubbleStopAtThisValue
       if (display) then
          write(*,*)'bg_field_evol: check for Hubble Stop enabled'
          write(*,*)'hubbleStopValue= ',hubbleStopValue
       endif
    else
       checkHubbleStop = BgEvolCheckForHubbleStop
       hubbleStopValue = BgEvolHubbleStopValue
    endif

    if (present(whichStopVelocitySign)) then
       stopVelocitySign = whichStopVelocitySign
       if (display) write(*,*)'bg_field_evol: stopVelocitySign is',stopVelocitySign
    else
       stopVelocitySign = BgEvolStopVelocitySign
    endif
    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    efoldHuge = min(efoldMaxiStop,efoldHuge)


!checks
    if (present(ptrStart)) then
       if (associated(ptrStart)) then
          stop 'bg_field_evol: ptr to bgdata already associated'
       endif
    endif



    if ((.not.useVelocity).and.efoldExploreOsc.ne.0) then
       write(*,*)'bg_field_evol: oscillation exploration disabled!'
       efoldExploreOsc = 0._kp
    endif

    
    
!set initial conditions   

    hubbleSquareIni = hubble_parameter_square(infIni%field,infIni%fieldDot,.false.)
    if (hubbleSquareIni.ge.0._kp) then
       hubbleIni = sqrt(hubbleSquareIni)
    else      
       stop 'bg_field_evol: hubbleSquareIni < 0'
    endif
    bgVarIni(1:fieldNum) = infIni%field(1:fieldNum)
    if (useVelocity) then
       bgVarIni(fieldNum+1:2*fieldNum) = infIni%fieldDot(1:fieldNum)*hubbleIni
    else
       bgVarIni(fieldNum+1:2*fieldNum) = infIni%fieldDot(1:fieldNum)
    endif


!localize rougthly the end of inflation: in fprime, a test stops dverk when epsilon>1
    efold = infIni%efold
    bgVar = bgVarIni

!initialize (all other subtypes may change)
    stopData%yesno1 = useOtherEpsilon
    stopData%mat = checkMatterStop
    stopData%ismax = stopAtMax
    stopData%hub = checkHubbleStop
    stopData%vsign = stopVelocitySign
    stopData%check = .false.
    stopData%update = .false.
    stopData%xend = efoldHuge
    stopData%real1 = epsilonStop
    stopData%real2 = matterStopValue
    stopData%real3 = hubbleStopValue
    stopData%int1 = stopIndexMin
    stopData%int2 = stopIndexMax

!after completion stopData%mat and stopData%hub are changed, they are
!true or false according to the outcome of the end of inflation
    
    
    if (useVelocity) then
!derivField=Dfield/Dtphys
       call easydverk(neqs,bg_field_dot_coupled,efold,bgVar,efoldHuge,tolEvol &
            ,stopData)
    else
!derivField=Dfield/Defold
       call easydverk(neqs,bg_field_dot_decoupled,efold,bgVar,efoldHuge,tolEvol &
            ,stopData)
    endif

    
!up to dverk exploration
!  efoldAfterEndInf = stopData%xend and something like bgVar =
!  stopData%ptr after allocation

!    print *,'efold bg',efold,bgVar
!    print *,'epsilon1',slowroll_first_parameter(bgVar(1:2), bgVar(3:4) &
!         , useVelocity)


    derivFieldAfterEndInf(1:fieldNum) = bgVar(fieldNum+1:2*fieldNum)

    efoldAfterEndInf = efold
    efoldBeforeEndInf = efoldAfterEndInf - efoldStepDefault
  
    
    if (efoldBeforeEndInf.le.infIni%efold) then
       if (display) write(*,*)'bg_field_evol: inflation is too short!'
       bg_field_evol = infIni       
       return
    endif

    


!precise determination of efoldEndInf up to tolEfoldEnd provided
!inflation is longer than efoldStepDefault efold


!in these cases, accurate determination of the end of inflation is not needed
    if ((.not.accurateEndInf).or.(efoldAfterEndInf.eq.efoldMaxiStop)) then
       if (display) write(*,*)'bg_field_evol: endinf not determined accurately'
       efoldEndInf = efoldAfterEndInf
       fieldEndInf = bgVar(1:fieldNum)
       derivFieldEndInf = bgVar(fieldNum+1:2*fieldNum)
       hubbleEndInf  &
            = sqrt(hubble_parameter_square(fieldEndInf,derivFieldEndInf,useVelocity))

    else


!careful localisation of the end of inflation move to
!efoldBeforeEndInf: background integration maybe unstable, precision
!is pushed to maximum accuracy


       tolBackEvol = epsilon(1._kp)

       if (useVelocity) then
          call tunedverk(neqs,bg_field_dot_coupled,efold,bgVar,efoldBeforeEndInf &
               ,tolBackEvol)        
       else
          call tunedverk(neqs,bg_field_dot_decoupled,efold,bgVar,efoldBeforeEndInf &
               ,tolBackEvol)
       endif
    
!use zbrent zero finder in [efoldBeforeEndInf, efoldAfterEndInf]
       findData%yesno1 = useVelocity
       findData%yesno2 = useOtherEpsilon       
       findData%real1 = efold

       allocate(findData%ptrvector1(2*fieldNum))
       allocate(findData%ptrvector2(2*fieldNum))

       findData%ptrvector1 = bgVar
       findData%real2 = tolEvol       

       findData%int1 = stopIndexMin
       findData%int2 = stopIndexMax

!find at tolEfoldEnd precision: the interval should be small
!(background integration of inflation is unstable)
!    print *,'go in zbrent set',efoldBeforeEndInf,efoldAfterEndInf,findData%real1    

       if (checkHubbleStop.and.stopData%hub) then
          
          findData%real3 = hubbleStopValue
          
          efoldEndInf = zbrent(find_endinf_hubble,efoldBeforeEndInf &
               ,efoldAfterEndInf,tolEfoldEnd,findData)

       elseif (checkMatterStop.and.stopData%mat) then

          findData%yesno3 = stopAtMax
          findData%real3 = matterStopValue
          
          efoldEndInf = zbrent(find_endinf_matter,efoldBeforeEndInf &
               ,efoldAfterEndInf,tolEfoldEnd,findData)
       else
          
          efoldEndInf = zbrent(find_endinf_epsilon,efoldBeforeEndInf &
               ,efoldAfterEndInf,tolEfoldEnd,findData)

       endif
       
!read the results in the findData buffer
       fieldEndInf = findData%ptrvector2(1:fieldNum)
       derivFieldEndInf = findData%ptrvector2(fieldNum+1:2*fieldNum)
       hubbleEndInf &
            = sqrt(hubble_parameter_square(fieldEndInf,derivFieldEndInf,useVelocity))
    endif



!set the output values to bg_field_evol return

    bg_field_evol%hubble = hubbleEndInf
    bg_field_evol%efold = efoldEndInf    
    bg_field_evol%field = fieldEndInf

    if (useVelocity) then
       bg_field_evol%fieldDot = derivFieldEndInf/hubbleEndInf
    else
       bg_field_evol%fieldDot = derivFieldEndInf
    endif
    bg_field_evol%epsilon1 = slowroll_first_parameter(fieldEndInf, derivFieldEndInf &
         , useVelocity)
    bg_field_evol%epsilon2 = slowroll_second_parameter(fieldEndInf, derivFieldEndInf &
         , useVelocity)
    bg_field_evol%epsilon3 = slowroll_third_parameter(fieldEndInf, derivFieldEndInf &
         , useVelocity)
    bg_field_evol%epsilon1JF = slowroll_first_parameter_JF(fieldEndInf, derivFieldEndInf &
         , useVelocity)

    if (associated(findData%ptrvector1)) deallocate(findData%ptrvector1)
    if (associated(findData%ptrvector2)) deallocate(findData%ptrvector2)

   
!save some data in memory if ptr input is present, recomputes the
!background knowing the end of inflation. Whatever the total number of
!efold, the last relevant efold (after efoldBeforeEndObs) are
!sampled according to efoldStepDefault

     efoldObs = efoldEndInf - efoldBeforeEndObs
!     print *,'efoldObs efoldEndInf',efoldObs,efoldEndInf

    if (present(efoldDataNum)) then
       if (efoldDataNum.le.1) return
    endif



    if (present(ptrStart)) then       
!for debugging       
       if (dump_file) then
          call delete_file('field.dat')
          call delete_file('derivfield.dat')
          call delete_file('geom.dat')
          call delete_file('epsilon1.dat')
          call delete_file('epsilon2.dat')
          call delete_file('epsilon3.dat')
          call delete_file('hubble.dat')
          call delete_file('potential.dat')
          call delete_file('confsquare.dat')
          call delete_file('a2chi2.dat')
       endif


       if (present(efoldDataNum)) then          
          if (efoldObs - infIni%efold.gt.0._kp) then
             longEnoughObs = .true.
!             efoldStepNoObs = 2._kp*(efoldObs - infIni%efold)/real(efoldDataNum-1)
             efoldStepObs = 2._kp*efoldBeforeEndObs/real(efoldDataNum-1,kp)
!             efoldStepNoObs = max(efoldStepNoObs,efoldStepObs)
          else
             longEnoughObs = .false.              
!             efoldStepNoObs = 0._kp
             efoldStepObs = (efoldEndInf  - infIni%efold)/real(efoldDataNum/2-1,kp)
          endif
       else
          longEnoughObs = .false.
!          efoldStepNoObs = efoldStepDefault
          efoldStepObs = efoldStepDefault
       endif

!initialisation      
       efold = infIni%efold      
       bgVar = bgVarIni

       inflate = .true.
       
       ptrCurrent => null()
       ptrPrevious => null()

!       i=0

       do 

!compute the saved physical quantities    
          field =  bgVar(1:fieldNum)
          derivField = bgVar(fieldNum+1:2*fieldNum) 

          hubbleSquare = hubble_parameter_square(field,derivField,useVelocity)
          epsilon1 = slowroll_first_parameter(field,derivField,useVelocity)
          epsilon2 =  slowroll_second_parameter(field,derivField,useVelocity)
          epsilon3 =  slowroll_third_parameter(field,derivField,useVelocity)
          epsilon1JF = slowroll_first_parameter_JF(field,derivField,useVelocity)

          if (hubbleSquare.ge.0._kp) then
             hubble = sqrt(hubbleSquare)
          else
             print *,'efold',efold
             stop 'bg_field_evol: hubbleSquare < 0'
          endif
          

!save the physical quantities as a chained list
         
          if (.not.associated(ptrStart)) then
             allocate(ptrStart); ptrStart%ptr => null(); ptrCurrent => ptrStart
          else
             allocate(ptrCurrent); ptrCurrent%ptr => null(); ptrPrevious%ptr => ptrCurrent
          endif
          ptrCurrent%bg%efold = efold
          ptrCurrent%bg%hubble = hubble
          ptrCurrent%bg%epsilon1 = epsilon1
          ptrCurrent%bg%epsilon1JF = epsilon1JF
          ptrCurrent%bg%epsilon2 = epsilon2
          ptrCurrent%bg%epsilon3 = epsilon3
          ptrCurrent%bg%field = field
!          i=i+1
!             print *,'stored efold',ptrCurrent%bg%efold,i
          if (useVelocity) then
             ptrCurrent%bg%fieldDot = derivField/hubble             
          else
             ptrCurrent%bg%fieldDot = derivField
          endif
          ptrPrevious => ptrCurrent
          
!slow down a lot the computation: for test!
          if (dump_file) then
             call livewrite('field.dat',efold,bgVar(1),bgVar(2))
             call livewrite('derivfield.dat',efold, ptrCurrent%bg%fieldDot(1))
             call livewrite('hubble.dat',efold,hubble)
             call livewrite('epsilon1.dat',efold,epsilon1,epsilon1JF)
             call livewrite('epsilon2.dat',efold,epsilon2)
             call livewrite('epsilon3.dat',efold,epsilon3)
             call livewrite('potential.dat',efold,potential(bgVar(1:fieldNum)))
             call livewrite('confsquare.dat',efold &
                  ,conformal_factor_square(bgVar(matterNum+1:fieldNum)))
             call livewrite('a2chi2.dat',efold &
                  ,conformal_factor_square(bgVar(matterNum+1:fieldNum)) &
                  *bgVar(1)*bgVar(1))
          endif
                       

          if (efold.eq.efoldObs) then
             if (present(infObs)) then
                if (infObs%efold.eq.infIni%efold) then
                   infObs%efold = efold
                   infObs%field = field
                   infObs%hubble = hubble
                   infObs%fieldDot = ptrCurrent%bg%fieldDot
                   infObs%epsilon1 = epsilon1
                   infObs%epsilon2 = epsilon2
                   infObs%epsilon1JF = epsilon1JF
!                   print *,'infObs set',infObs
                endif
             endif
          endif

          
          if (longEnoughObs.and.(efold.lt.efoldObs)) then
             efoldNext = min(efold + efoldStepObs,efoldObs)
          else
             efoldNext = min(efold + efoldStepObs,efoldAfterEndInf + efoldExploreOsc)
          endif
               
!          print *,'efoldNExt',efoldNext
!          read(*,*)


!avoid the next step
          if (.not.inflate) exit

!integration again with the stopping criteria
          if (useVelocity) then
             if (efoldNext.ge.(efoldAfterEndInf + efoldExploreOsc)) then
                inflate = .false.
             endif
!derivField=Dfield/Dtphys
             call easydverk(neqs,bg_field_dot_coupled,efold,bgVar,efoldNext,tolEvol)

          else
             if (efoldNext.ge.efoldAfterEndInf) then
                inflate = .false.
             endif
!derivField=Dfield/Defold
             call easydverk(neqs,bg_field_dot_decoupled,efold,bgVar,efoldNext,tolEvol)
          endif                      

       enddo
    
       ptrCurrent => null()
       ptrPrevious => null()

    endif


  end function bg_field_evol




  function find_endinf_epsilon(efold,findData)
    use infsolvers, only : tunedverk
    use fieldprec, only : transfert
    implicit none
    real(kp), intent(in) :: efold
    type(transfert), optional, intent(inout) :: findData
    real(kp) :: find_endinf_epsilon

    logical :: useVelocity, useOtherEpsilon
    real(kp) :: efoldStart
    real(kp), dimension(2*fieldNum) :: bgVar
    real(kp), dimension(fieldNum) :: field, derivField

    useOtherEpsilon = findData%yesno2
    useVelocity = findData%yesno1
    efoldStart = findData%real1


    if ((.not.associated(findData%ptrvector1)) &
         .or.(.not.associated(findData%ptrvector2))) then
       stop 'find_endif: ptrvector not associated'
    endif

    bgVar = findData%ptrvector1
       
!backward integration from efoldStart found in bg_evol to the wanted efold
    if (useVelocity) then
       call tunedverk(2*fieldNum,bg_field_dot_coupled,efoldStart,bgVar &
            ,efold,findData%real2)
    else
       call tunedverk(2*fieldNum,bg_field_dot_decoupled,efoldStart,bgVar &
         ,efold,findData%real2)
    endif

    findData%ptrvector2 = bgVar

    field = bgVar(1:fieldNum)
    derivField = bgVar(fieldNum+1:2*fieldNum)


!difference between the epsilon corresponding to the current efold and
!the one wanted to end inflation (=1). The efoldEnd is the zero of
!this function. Striclty speaking inflation ends when epsilon1JF=1.
        
    if (useOtherEpsilon) then
       find_endinf_epsilon &
            = ptr_alternate_stop_parameter(field,derivField,findData%yesno1) - 1._kp
    else
       find_endinf_epsilon &
            = slowroll_first_parameter(field,derivField,findData%yesno1) - 1._kp
    endif

    findData%real3 = find_endinf_epsilon + 1._kp
       

  end function find_endinf_epsilon




  function find_endinf_matter(efold,findData)
    use infsolvers, only : tunedverk
    use fieldprec, only : transfert
    implicit none
    real(kp), intent(in) :: efold
    type(transfert), optional, intent(inout) :: findData
    real(kp) :: find_endinf_matter

    integer :: ifMin,ifMax
    logical :: useVelocity, stopForMatterMax
    real(kp) :: efoldStart, matterStop
    real(kp), dimension(2*fieldNum) :: bgVar
    real(kp), dimension(fieldNum) :: field

    
    useVelocity = findData%yesno1
    stopForMatterMax = findData%yesno3
    efoldStart = findData%real1
    matterStop = findData%real3

    ifMin = findData%int1
    ifMax = findData%int2


    if ((.not.associated(findData%ptrvector1)) &
         .or.(.not.associated(findData%ptrvector2))) then
       stop 'find_endif: ptrvector not associated'
    endif

    bgVar = findData%ptrvector1
       
!backward integration from efoldStart found in bg_evol to the wanted efold
    if (useVelocity) then
       call tunedverk(2*fieldNum,bg_field_dot_coupled,efoldStart,bgVar &
            ,efold,findData%real2)
    else
       call tunedverk(2*fieldNum,bg_field_dot_decoupled,efoldStart,bgVar &
         ,efold,findData%real2)
    endif

    findData%ptrvector2 = bgVar

    field = bgVar(1:fieldNum)
   
!difference between the min matter field value corresponding to the
!current efold and the one wanted to end inflation
!(matterMinStop). The efoldEnd is the zero of this function.
     
    if (stopForMatterMax) then
       find_endinf_matter =  maxval(field(ifMin:ifMax)) -  matterStop
    else
       find_endinf_matter =  minval(field(ifMin:ifMax)) -  matterStop 
    endif


  end function find_endinf_matter

  

  function find_endinf_hubble(efold,findData)
    use infsolvers, only : tunedverk
    use fieldprec, only : transfert
    implicit none
    real(kp), intent(in) :: efold
    type(transfert), optional, intent(inout) :: findData
    real(kp) :: find_endinf_hubble
    
    logical :: useVelocity
    real(kp) :: efoldStart, hubbleStop,hubbleSquare
    real(kp), dimension(2*fieldNum) :: bgVar
    real(kp), dimension(fieldNum) :: field,derivField

    
    useVelocity = findData%yesno1
    
    efoldStart = findData%real1
    hubbleStop = findData%real3
  
    if ((.not.associated(findData%ptrvector1)) &
         .or.(.not.associated(findData%ptrvector2))) then
       stop 'find_endif: ptrvector not associated'
    endif

    bgVar = findData%ptrvector1
       
!backward integration from efoldStart found in bg_evol to the wanted efold
    if (useVelocity) then
       call tunedverk(2*fieldNum,bg_field_dot_coupled,efoldStart,bgVar &
            ,efold,findData%real2)
    else
       call tunedverk(2*fieldNum,bg_field_dot_decoupled,efoldStart,bgVar &
         ,efold,findData%real2)
    endif

    findData%ptrvector2 = bgVar

    field = bgVar(1:fieldNum)
    derivField = bgVar(fieldNum+1:2*fieldNum)

    hubbleSquare = hubble_parameter_square(field,derivField,useVelocity)


    find_endinf_hubble =  hubbleSquare -  hubbleStop**2


  end function find_endinf_hubble




  subroutine bg_field_dot_decoupled(neqs,efold,bgVar,bgVarDot,stopData)
!for derivField=Dfield/Defold the field equations decouple from the
!Hubble flow. However, this decomposition becomes singular when the
!potential vanishes and the integration fails a that point. Harmless
!for the inflationary era.

    use fieldprec, only : transfert
    implicit none
        
    integer :: neqs    
    real(kp) :: efold
    real(kp), dimension(neqs) :: bgVar
    real(kp), dimension(neqs) :: bgVarDot
    type(transfert), optional, intent(inout) :: stopData
   
    integer :: i
    real(kp), dimension(fieldNum) :: dlnPotVec
    real(kp), dimension(fieldNum) :: field
    real(kp), dimension(fieldNum) :: fieldDot, fieldDotDot
    real(kp) :: epsilon, hubbleSquare
      
    logical :: stopNow

    stopNow = .false.

    field = bgVar(1:fieldNum)
    fieldDot = bgVar(fieldNum+1:2*fieldNum)

!    fieldDotSquare = dot_product(fieldDot,matmul(metric(field),fieldDot))

    if (present(stopData)) then
!use other method or not to stop inflation
       if (stopData%check) then
          
          if (stopData%yesno1) then
             epsilon = ptr_alternate_stop_parameter(field,fieldDot,.false.)
          else
             !epsilon = fieldDotSquare/2._kp
             epsilon = slowroll_first_parameter(field,fieldDot,.false.)
          endif

          if (stopData%mat) then
             if (stopData%ismax) then
                stopNow = (maxval(field(stopData%int1:stopData%int2)) &
                     .gt.stopData%real2)
             else
                stopNow = (minval(field(stopData%int1:stopData%int2)) &
                     .lt.stopData%real2)
             endif
             if (stopNow) stopData%hub = .false.
          endif

          if (stopData%hub) then
             hubbleSquare = hubble_parameter_square(field,fieldDot,.false.)
             stopNow = (hubbleSquare.lt.stopData%real3**2)
             if (stopData%vsign > 0) then
                stopNow = stopNow .and. (fieldDot(stopData%int1).ge.0._kp)
             elseif (stopData%vsign < 0) then
                stopNow = stopNow .and. (fieldDot(stopData%int1).le.0._kp)
             endif
             if (stopNow) stopData%mat = .false.
          endif

          if (epsilon.gt.stopData%real1) then
             stopNow = .true.
             stopData%mat = .false.
             stopData%hub = .false.
          endif
          
          if (stopNow) then
             stopData%update = .true.
             stopData%xend = efold 
          endif

       endif
    endif
    
!    do i=1,fieldNum
!       christVec(i) = dot_product(fieldDot,matmul(christoffel(i,:,:),fieldDot))
!    enddo
     
!    dlnPotVec = matmul(metric_inverse(field),deriv_ln_potential(field))
          
!    fieldDotDot = -christVec - (3._kp - fieldDotSquare/2._kp)*(fieldDot + dlnPotVec)

!    bgVarDot(1:fieldNum) = fieldDot
!    bgVarDot(fieldNum+1:2*fieldNum) = fieldDotDot

    bgVarDot(1:fieldNum) = fieldDot
    bgVarDot(fieldNum+1:2*fieldNum) = field_second_derivative(field,fieldDot)
    
  end subroutine bg_field_dot_decoupled





  subroutine bg_field_dot_coupled(neqs,efold,bgVar,bgVarDot,stopData)                       
!for derivField=Dfield/Dtphys the field equations are coupled to the
!hubble flow. This avoid the singular behavior of the decoupled
!equations and allows to properly sample the oscillations of the field
!at the end of inflation.

    use fieldprec, only : transfert

    implicit none
        
    integer :: neqs    
    real(kp) :: efold
    real(kp), dimension(neqs) :: bgVar
    real(kp), dimension(neqs) :: bgVarDot
    type(transfert), optional, intent(inout) :: stopData

    integer :: i
    real(kp), dimension(fieldNum) :: dPotVec
    real(kp), dimension(fieldNum) :: field, velocity
    real(kp), dimension(fieldNum) :: fieldDot
    real(kp) :: epsilon, hubbleSquare, hubble
    
    logical :: stopNow

    stopNow = .false.

    field = bgVar(1:fieldNum)
    velocity = bgVar(fieldNum+1:2*fieldNum)

!    christoffel = connection_affine(field)
!    velocitySquare = dot_product(velocity,matmul(metric(field),velocity))    
!    do i=1,fieldNum
!       christVec(i) = dot_product(velocity(:),matmul(christoffel(i,:,:),velocity(:)))
!    enddo
   
!    dPotVec = matmul(metric_inverse(field),deriv_potential(field))
   
    hubbleSquare = hubble_parameter_square(field,velocity,.true.)
    if (hubbleSquare.ge.0.) then
       hubble = sqrt(hubbleSquare)
    else
       print *,'efold= ',efold
       print *,'field= ',field
       print *,'velocity= ',velocity 
       stop 'bg_field_dot_coupled: hubbleSquare < 0'
    endif
    
    if (present(stopData)) then
!use other method to stop inflation
       if (stopData%check) then
          if (stopData%yesno1) then
             epsilon = ptr_alternate_stop_parameter(field,velocity,.true.)
          else
             !             epsilon = velocitySquare/2._kp/hubbleSquare
             epsilon = slowroll_first_parameter(field,velocity,.true.)
          endif          

          if (stopData%mat) then
             if (stopData%ismax) then
                stopNow = (maxval(field(stopData%int1:stopData%int2)) &
                     .gt.stopData%real2)
             else
                stopNow = (minval(field(stopData%int1:stopData%int2)) &
                     .lt.stopData%real2)
             endif
             if (stopNow) stopData%hub = .false.
          endif
          
          
          if (stopData%hub) then
             stopNow = (hubble.le.stopData%real3)
             if (stopData%vsign > 0) then
                stopNow = stopNow .and. (velocity(stopData%int1).ge.0._kp)
             elseif (stopData%vsign < 0) then
                stopNow = stopNow .and. (velocity(stopData%int1).le.0._kp)
             endif
             if (stopNow) stopData%mat = .false.
          endif

          if (epsilon.gt.stopData%real1) then
             stopNow = .true.
             stopData%mat = .false.
             stopData%hub = .false.
          endif
          
          if (stopNow) then
             stopData%update = .true.
             stopData%xend = efold 
          endif         
          
!          print *,'efold field',efold,field
!          print *,'eps1 eps2',epsilon1,slowroll_second_parameter(field,velocity,.true.)
!          read(*,*)
       endif              
    endif
   
    fieldDot = velocity/hubble
!    velocityDot =  -3._kp*velocity - (christVec + dPotVec)/hubble
       
!    bgVarDot(1:fieldNum) = fieldDot
!    bgVarDot(fieldNum+1:2*fieldNum) = velocityDot

    bgVarDot(1:fieldNum) = fieldDot
    bgVarDot(fieldNum+1:2*fieldNum) = velocity_derivative(field,velocity)
    
  end subroutine bg_field_dot_coupled


end module infbg
