module infbg
  use infprec, only : kp, tolkp
  use infbgmodel, only : matterNum, dilatonNum, fieldNum
  use infbgfunc, only : hubble_parameter_square, slowroll_first_parameter_JF
  use infbgfunc, only : slowroll_first_parameter, slowroll_second_parameter
  

!background evolution in the Einstein FLRW Frame for multifields:
!scalar gravity + matter fields.

  implicit none

  private

   
!for debugging
  logical, parameter :: display = .false.
  logical, parameter :: dump_file = .false.
 

!to store snapshot (ini or end, or more)
  type infbgphys    
     sequence
     real(kp) :: efold, hubble, epsilon1, epsilon1JF
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

  public set_infbg_ini
  public rescale_potential
  public bg_field_evol, bg_field_dot_coupled
         

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

  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!inflation settings: initial conditions, rescaling, normalisation...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  function set_infbg_ini(infParam)
!to start on the attractor. This is epsilon2 = 0 for epsilon1<<1
!initialy.
    use infbgmodel, only : infbgparam
    use infsigma, only : metric_inverse
    use infpotential, only : deriv_ln_potential
    use infsric, only : slowroll_initial_matter

    implicit none
    type(infbgparam), intent(in) :: infParam
    type(infbgphys) :: set_infbg_ini
    
    real(kp) :: hubbleSquareIni
    type(infbgphys) :: infIni


    infIni%field(1:matterNum) = infParam%matters
    infIni%field(matterNum+1:fieldNum) = infParam%conforms


!if the matter fields are set to 0, use slow-roll guesses
    if (all(infParam%matters == 0._kp)) then
       infIni%field(1:matterNum) = slowroll_initial_matter(infParam)
    endif

!this is starting on slow-roll
    infIni%fieldDot &
         = - matmul(metric_inverse(infIni%field),deriv_ln_potential(infIni%field))
!overides for test
!    infIni%fieldDot = 0._kp
!    print *,'initial condition fieldDot=',infIni%fieldDot

!    infIni%fieldDot(1)=1.
!    infIni%fieldDot(2)=-0.5
!    infIni%fieldDot(3)=0.6

    hubbleSquareIni = hubble_parameter_square(infIni%field,infIni%fieldDot,.false.)

    if (hubbleSquareIni.ge.0._kp) then
       infIni%hubble = sqrt(hubbleSquareIni)       
    else     
       stop 'H^2 < 0, check initial condition'
    endif

    infIni%epsilon1 = slowroll_first_parameter(infIni%field,infIni%fieldDot, .false.)
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
    infParam%consts(1) = infParam%consts(1)*scale**0.25

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



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!inflationary evolution: find end of inflation + store relevant quantities
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  function bg_field_evol(infIni,efoldDataNum,infObs,ptrStart,stopAtThisValue,isStopAtMax) 
!integrate the background until epsilon > epsilonStop, and returns some
!physical values (type infbgphys) for which epsilon = epsilon1EndInf  (=1)

    use infprec, only : transfert
    use inftools, only : easydverk, tunedverk, zbrent
    use infdilaton, only : conformal_factor_square
    use infpotential, only : potential
    use infinout
    
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
    type(infbgphys) :: bg_field_evol
    
    
    real(kp) :: epsilon1, epsilon1JF
    real(kp) :: epsilon2

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
    
    real(kp) :: efold,efoldNext,efoldStepObs, efoldStepNoObs

    real(kp) :: efoldBeforeEndInf,efoldEndInf
    real(kp) :: efoldAfterEndInf

    real(kp), dimension(fieldNum) :: field,derivField
    real(kp), dimension(fieldNum) :: fieldEndInf,derivFieldEndInf
    real(kp), dimension(fieldNum) :: derivFieldAfterEndInf
    real(kp), dimension(2*fieldNum) :: bgVar, bgVarIni

!standard integration accuracy 
!    real(kp), parameter :: tolEvol = 1e-11
    real(kp), parameter :: tolEvol = tolkp

!backward integration accuracy (sometime instable, mayneed extra precision)
    real(kp) :: tolBackEvol

!if true derivField=Dfield/Dtphys, otherwise derivField=Dfield/Defold is
!used for the integration. In both cases, only Dfield/Defold is
!stored
    logical, parameter :: useVelocity = .true.
!end inflation when epsilon1=1 in Jordan Frame, or in Einstein Frame
!Physics says in JF, but both are the same up to 2% when the dilaton coupling are set to 1
!Today dilaton couplings are 0.01 maxi, and they are constant or null in our model.
!Integration stops when epsilon1(useJF or not) > epsilon1Stop
    real(kp), parameter :: epsilon1Stop = 1
    logical, parameter :: useEpsilon1JF = .false.

!zbrent accuracy on efoldEnd for which epsilon=epsilon1Stop
    real(kp), parameter :: tolEfoldEnd = tolkp

!another test, checked after epsilon1 values to stop integration. May
!be convenient for hybrid-like one field potential: stop integration
!when the min value of the matter fields is below matterMiniStop, or
!according to the total number of efolds
    logical, parameter :: accurateEndInf = .true.
    real(kp) :: efoldMaxiStop
    logical :: checkHubbleStop
    logical :: checkMatterStop
    logical :: stopForMax   

    real(kp) :: valueStop
    integer, parameter :: stopIndexMin = 1, stopIndexMax = 1

    logical :: inflate, longEnoughObs
  
!to make f77 routines discussing together           
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

    efoldExploreOsc = 0.

!enabled by true
   
    checkHubbleStop = .false.

    checkMatterStop = .false.

    efoldMaxiStop = 200.

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
       efoldExploreOsc = 0.
    endif

    if (present(isStopAtMax)) then
       stopForMax = isStopAtMax
    else
       stopForMax = .false.
    endif

    if (present(stopAtThisValue)) then
       valueStop = stopAtThisValue
       if ((.not.checkMatterStop).and.checkHubbleStop) then
          if (display) write(*,*)'bg_field_evol: check for Hubble stop enabled'          
       else
          if (display) write(*,*)'bg_field_evol: check for Matter stop enabled'          
          if (display) write(*,*)'bg_field_evol: stopForMax is',stopForMax
          checkHubbleStop = .false.
          checkMatterStop = .true.
       endif
    else
       valueStop = tolEvol
    endif   
    
!set initial conditions   

    hubbleSquareIni = hubble_parameter_square(infIni%field,infIni%fieldDot,.false.)
    if (hubbleSquareIni.ge.0d0) then
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
    stopData%yesno1 = useEpsilon1JF
    stopData%yesno2 = checkMatterStop
    stopData%yesno3 = stopForMax
    stopData%yesno4 = checkHubbleStop
    stopData%check = .false.    
    stopData%update = .false.
    stopData%xend = efoldHuge
    stopData%real1 = epsilon1Stop
    stopData%real2 = valueStop! - 10._kp*tolEvol
    stopData%int1 = stopIndexMin
    stopData%int2 = stopIndexMax

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
       if (display) write(*,*)'bg_field_evol: inflation too short'
       bg_field_evol = infIni       
       return
    endif

    


!precise determination of efoldEndInf up to tolEfoldEnd provided
!inflation is longer than efoldStepDefault efold

!checkMatterMini stands for cases when epsilon1=epsilon1Stop does not
!define the end of inflation, so who cares about accurate
!determination
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
       findData%yesno2 = useEpsilon1JF       
       findData%yesno3 = stopForMax
       findData%real1 = efold
       allocate(findData%ptrvector1(2*fieldNum))
       allocate(findData%ptrvector2(2*fieldNum))
       findData%ptrvector1 = bgVar
       findData%real2 = tolEvol       
       findData%real3 = valueStop

       findData%int1 = stopIndexMin
       findData%int2 = stopIndexMax

!find at tolEfoldEnd precision: the interval should be small
!(background integration of inflation is unstable)
!    print *,'go in zbrent set',efoldBeforeEndInf,efoldAfterEndInf,findData%real1    

       if (checkHubbleStop.and.(.not.stopData%yesno2)) then
          efoldEndInf = zbrent(find_endinf_hubble,efoldBeforeEndInf &
               ,efoldAfterEndInf,tolEfoldEnd,findData)
       elseif (checkMatterStop.and.(.not.stopData%yesno2)) then
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
          call delete_file('hubble.dat')
          call delete_file('potential.dat')
          call delete_file('confsquare.dat')
          call delete_file('a2chi2.dat')
       endif


       if (present(efoldDataNum)) then          
          if (efoldObs - infIni%efold.gt.0.) then
             longEnoughObs = .true.
             efoldStepNoObs = 2.*(efoldObs - infIni%efold)/real(efoldDataNum-1)
             efoldStepObs = 2.*efoldBeforeEndObs/real(efoldDataNum-1)
             efoldStepNoObs = max(efoldStepNoObs,efoldStepObs)
          else
             longEnoughObs = .false.              
             efoldStepNoObs = 0.
             efoldStepObs = (efoldEndInf  - infIni%efold)/real(efoldDataNum/2-1)
          endif
       else
          longEnoughObs = .false.
          efoldStepNoObs = efoldStepDefault
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
          epsilon1JF = slowroll_first_parameter_JF(field,derivField,useVelocity)

          if (hubbleSquare.ge.0d0) then
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
             epsilon2 = slowroll_second_parameter(field,derivField,useVelocity)
             call livewrite('epsilon2.dat',efold,epsilon2)
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
                   infObs%epsilon1JF = slowroll_first_parameter_JF(field, derivField &
                        , useVelocity)
!                   print *,'infObs set',infObs
                endif
             endif
          endif

          
          if (longEnoughObs.and.(efold.lt.efoldObs)) then
             efoldNext = min(efold + efoldStepNoObs,efoldObs)
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
    use inftools, only : tunedverk
    use infprec, only : transfert
    implicit none
    real(kp), intent(in) :: efold
    type(transfert), optional, intent(inout) :: findData
    real(kp) :: find_endinf_epsilon

    logical :: useVelocity, useEpsilon1JF
    real(kp) :: efoldStart
    real(kp), dimension(2*fieldNum) :: bgVar
    real(kp), dimension(fieldNum) :: field, derivField

    useEpsilon1JF = findData%yesno2
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
        
    if (useEpsilon1JF) then
       find_endinf_epsilon &
            = slowroll_first_parameter_JF(field,derivField,findData%yesno1) - 1._kp
    else
       find_endinf_epsilon &
            = slowroll_first_parameter(field,derivField,findData%yesno1) - 1._kp
    endif

    findData%real3 = find_endinf_epsilon + 1._kp
       

  end function find_endinf_epsilon




  function find_endinf_matter(efold,findData)
    use inftools, only : tunedverk
    use infprec, only : transfert
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
    use inftools, only : tunedverk
    use infprec, only : transfert
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

    use infprec, only : transfert
    use infsigma, only : metric, metric_inverse   
    use infsigma, only : connection_affine
    use infpotential, only : deriv_ln_potential
    implicit none
        
    integer :: neqs    
    real(kp) :: efold
    real(kp), dimension(neqs) :: bgVar
    real(kp), dimension(neqs) :: bgVarDot
    type(transfert), optional, intent(inout) :: stopData
   
    integer :: i
    real(kp), dimension(fieldNum) :: dlnPotVec
    real(kp), dimension(fieldNum) :: field, christVec
    real(kp), dimension(fieldNum) :: fieldDot, fieldDotDot 
    real(kp), dimension(fieldNum,fieldNum,fieldNum) :: christoffel 
    real(kp) :: fieldDotSquare, epsilon1, hubbleSquare
  
    logical, save :: stopNow=.false.


    field = bgVar(1:fieldNum)
    fieldDot = bgVar(fieldNum+1:2*fieldNum)
    christoffel = connection_affine(field)

    fieldDotSquare = dot_product(fieldDot,matmul(metric(field),fieldDot))

    if (present(stopData)) then
!use epsilon1JF or not to stop inflation
       if (stopData%check) then
          
          if (stopData%yesno1) then
             epsilon1 = slowroll_first_parameter_JF(field,fieldDot,.false.)
          else
             epsilon1 = fieldDotSquare/2d0
          endif

          if (stopData%yesno2) then
             if (stopData%yesno3) then
                stopNow = (maxval(field(stopData%int1:stopData%int2)) &
                     .gt.stopData%real2)
             else
                stopNow = (minval(field(stopData%int1:stopData%int2)) &
                     .lt.stopData%real2)
             endif
          endif

          if (stopData%yesno4) then
             hubbleSquare = hubble_parameter_square(field,fieldDot,.false.)
             stopNow = (hubbleSquare.lt.(stopData%real2)**2)
          endif
          
          if (stopNow.or.(epsilon1.gt.stopData%real1)) then
             stopData%update = .true.
             stopData%xend = efold 
             stopData%yesno2 = .false.
             stopData%yesno4 = .false.
          endif

       endif
    endif
    
    do i=1,fieldNum
       christVec(i) = dot_product(fieldDot,matmul(christoffel(i,:,:),fieldDot))
    enddo
 
!    dlnPotVec = deriv_ln_potential_vec(field)
    
    dlnPotVec = matmul(metric_inverse(field),deriv_ln_potential(field))
          
    fieldDotDot = -christVec - (3d0 - fieldDotSquare/2d0)*(fieldDot + dlnPotVec)

    bgVarDot(1:fieldNum) = fieldDot
    bgVarDot(fieldNum+1:2*fieldNum) = fieldDotDot
       
  end subroutine bg_field_dot_decoupled





  subroutine bg_field_dot_coupled(neqs,efold,bgVar,bgVarDot,stopData)                       
!for derivField=Dfield/Dtphys the field equations are coupled to the
!hubble flow. This avoid the singular behavior of the decoupled
!equations and allows to properly sample the oscillations of the field
!at the end of inflation.

    use infprec, only : transfert
    use infsigma, only : metric, metric_inverse, connection_affine
    use infpotential, only : deriv_potential
    implicit none
        
    integer :: neqs    
    real(kp) :: efold
    real(kp), dimension(neqs) :: bgVar
    real(kp), dimension(neqs) :: bgVarDot
    type(transfert), optional, intent(inout) :: stopData

    integer :: i
    real(kp), dimension(fieldNum) :: dPotVec
    real(kp), dimension(fieldNum) :: field, velocity
    real(kp), dimension(fieldNum) :: fieldDot, velocityDot, christVec
    real(kp), dimension(fieldNum,fieldNum,fieldNum) :: christoffel 
    real(kp) :: velocitySquare, epsilon1, hubbleSquare, hubble
    
    logical, save :: stopNow = .false.


    field = bgVar(1:fieldNum)
    velocity = bgVar(fieldNum+1:2*fieldNum)
    christoffel = connection_affine(field)
   
    velocitySquare = dot_product(velocity,matmul(metric(field),velocity))    

    do i=1,fieldNum
       christVec(i) = dot_product(velocity(:),matmul(christoffel(i,:,:),velocity(:)))
    enddo
   
!    dPotVec = deriv_potential_vec(field)
    dPotVec = matmul(metric_inverse(field),deriv_potential(field))
   
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
!use epsilon1JF or not to stop inflation
       if (stopData%check) then
          if (stopData%yesno1) then
             epsilon1 = slowroll_first_parameter_JF(field,velocity,.true.)
          else
             epsilon1 = velocitySquare/2d0/hubbleSquare
          endif          
          if (stopData%yesno2) then
             if (stopData%yesno3) then
                stopNow = (maxval(field(stopData%int1:stopData%int2)) &
                     .gt.stopData%real2)
             else
                stopNow = (minval(field(stopData%int1:stopData%int2)) &
                     .lt.stopData%real2)
             endif
          endif
             
          if (stopData%yesno4) stopNow = hubble.le.stopData%real2


          if (stopNow.or.(epsilon1.gt.stopData%real1)) then
             stopData%update = .true.
             stopData%xend = efold 
             stopData%yesno2 = .false.
             stopData%yesno4 = .false.
          endif
          
!          print *,'efold field',efold,field
!          print *,'eps1 eps2',epsilon1,slowroll_second_parameter(field,velocity,.true.)
!          read(*,*)
       endif              
    endif
   
    fieldDot = velocity/hubble
    velocityDot =  -3d0*velocity - (christVec + dPotVec)/hubble
       
    bgVarDot(1:fieldNum) = fieldDot
    bgVarDot(fieldNum+1:2*fieldNum) = velocityDot
            
  end subroutine bg_field_dot_coupled


end module infbg
