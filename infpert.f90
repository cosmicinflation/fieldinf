module infpert  
!inflationary perturbations for scalar fields and tensor modes, in the
!original field basis (might have precision issues in the
!determination of ridiculously small entropy modes). efold is the
!number of efold running from efoldIni to efoldEnd during the
!inflationary era. bfold is the number of efold before the end of
!inflation: bfold = efold - efoldEnd running from efoldIni-efoldEnd to
!0.

  use infprec, only : kp, tolkp
  use infbgmodel, only : matterNum, dilatonNum, fieldNum
  use infbg, only : infbgphys
  use inftorad, only : inftoradcosmo
  implicit none

  private

  
!for debugging
  logical, parameter :: display = .false.
  logical, parameter :: dump_spectra = .false.
  logical, parameter :: dump_modes = .false.
 
!+1 due to the bardeen potential.
  integer, parameter :: scalNum = fieldNum + 1

  real(kp), parameter :: kphysOverHubbleCreate = 10._kp
 
  integer, parameter :: entroRef = 1


  public scalNum  
  public power_spectrum_tens, power_spectrum_scal
  

contains
 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!gravitational waves
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  function power_spectrum_tens(infCosmo,kmpc)        
    use inftorad, only : infhubblexit, hubble_splinexit
    use infinout, only : livewrite
    implicit none
    real(kp) :: power_spectrum_tens
    type(inftoradcosmo), intent(in) :: infCosmo    
    real(kp), intent(in) :: kmpc
    
    type(infhubblexit) :: atHkExit
    complex(kp) :: tensEnd
    real(kp) :: pi

    pi=acos(-1._kp)

!evolution of k^3/2 x h by recomputing the background (faster)
    tensEnd = pert_tens_bgevol(infCosmo,kmpc)
        
    power_spectrum_tens = (2._kp/pi/pi)*real(conjg(tensEnd)*tensEnd)

    if (dump_spectra) then      
       call livewrite('powerhK.dat',kmpc,power_spectrum_tens)
       atHkExit = hubble_splinexit(infCosmo,kmpc)
       call livewrite('powerhN.dat',atHkExit%bfold,power_spectrum_tens)
    endif

    if (display) then
        write(*,*)'power_spectrum_tens:'
       write(*,*)'kmpc= Ph= '&
            ,kmpc,power_spectrum_tens
    endif


    if (dump_modes) then
       write(*,*)'tensor evolution dumped!'
       read(*,*)
    endif
   
  end function power_spectrum_tens




  function pert_tens_bgevol(infCosmo,kmpc)
!this evolves the bg and tens = k^3/2 * h simultaneously
    use inftools, only : easydverk
    use infprec, only : transfert
    use infbg, only : operator(/=)
    use infbg, only : bg_field_dot_coupled
    use infbgfunc, only : hubble_parameter_square
    use inftorad, only : bfold_hubble_fraction, lnMpcToKappa
    use infinout

    implicit none

    type(inftoradcosmo), intent(in) :: infCosmo    
    real(kp), intent(in) :: kmpc        
    complex(kp) :: pert_tens_bgevol

    type(transfert) :: cosmoData 

!we do not want to recompute the background from billion efolds. It is
!enough to start just before the relevant modes today are created

!workaround for old omp threadprivate directive (common for
!static + save)
    real(kp) :: efoldTune
    real(kp), dimension(fieldNum) :: fieldTune, velocityTune
    type(infbgphys) :: bgPrevious
    common /bgTensSave/ efoldTune, fieldTune, velocityTune, bgPrevious
    save /bgTensSave/
!$omp threadprivate(/bgTensSave/)

    real(kp) :: bfold, efold, efoldStart, bfoldStart, bfoldNext

!output files steps 
    integer, parameter :: bfoldDataNum = 5000
    real(kp) :: bfoldStep

!accuracy for the forward integration of background only and both
!background and perturbations (any backward (unstable) integration is
!avoided). Since sub-Hubble perturbations oscillate, funny low
!accuracy on pertevol lead to funny long integration time
    real(kp), parameter :: tolBgEvol = tolkp
    real(kp), parameter :: tolEvol = 1e-11   
    
    integer, parameter :: neqs = 4 + 2*fieldNum
    integer, parameter :: neqsbg = 2*fieldNum

    complex(kp) :: tens,tensDot,tensStart,tensDotStart    
    real(kp), dimension(2*fieldNum) :: bgVar
!tens and bg together
    real(kp), dimension(neqs) :: allVar
    real(kp), dimension(fieldNum) :: field,velocity
    real(kp) :: hubbleSquare,kphysOverHubbleStart




!for transferring kmpc + other cosmo to the differential equations
    cosmoData%real1 = kmpc
    cosmoData%real2 = infCosmo%bgEnd%efold
    cosmoData%real3 = infCosmo%efoldEndToToday

!find the bfold of creation
! solution of N + ln[kappa H(N)] = ln(kmpc) -ln(aend/a0) -Nc -ln[(k/aH)_creation]    
    kphysOverHubbleStart = kphysOverHubbleCreate
    bfoldStart = bfold_hubble_fraction(kmpc,infCosmo,kphysOverHubbleStart)

    if (display) then
       write(*,*)
       write(*,*)'tensor modes: kmpc = ',kmpc
       write(*,*)'bfold quantum = ',bfoldStart
       write(*,*)
    endif


    if (dump_modes) then
       if (bfoldDataNum.le.1) stop 'pert_tens_bgevol: 1 data point?'
       bfoldStep = (infCosmo%bfoldEnd - bfoldStart)/real(bfoldDataNum-1)
       call delete_file('tens.dat')
    else
!turbo settings, let's go to the end of inflation
       bfoldStep = 1._kp/epsilon(1._kp)
    endif


!evolve the bg from efoldIni to efoldStart of mode creation only
!once. After move from previous bgStart to the new efoldStart of mode
!creation, ONLY if the mode creation is at greater efold than the
!efoldStart of the previous mode: backward integration of the bg is
!unstable by rk methods. We also check that the bginf
!model is the same to do that. When all these conditions are not filled
!we start again from infCosmo
    
!test if the bg model change print *,'testcond',
!    ((infCosmo%bgIni%efold.ne.bgPrevious%efold) &
!    .or.(any(infCosmo%bgIni%field.ne.bgPrevious%field)) &
!    .or.(any(infCosmo%bgIni%fieldDot.ne.bgPrevious%fieldDot))),
!    (infCosmo%bgIni /= bgPrevious)

    if (infCosmo%bgIni /= bgPrevious) then
       efoldTune = infCosmo%bgIni%efold
       fieldTune = infCosmo%bgIni%field
       velocityTune = infCosmo%bgIni%fieldDot * infCosmo%bgIni%hubble
    endif

    bgPrevious = infCosmo%bgIni
   
    efoldStart = bfoldStart + infCosmo%bgEnd%efold

!only forward integration allowed
    if (efoldStart.lt.efoldTune) then      
       efoldTune = infCosmo%bgIni%efold
       fieldTune = infCosmo%bgIni%field
       velocityTune = infCosmo%bgIni%fieldDot * infCosmo%bgIni%hubble
       if (efoldStart.lt.infCosmo%bgIni%efold) then
          write(*,*) 'pert_tens_bgevol: kmpc =',kmpc
          write(*,*) 'bfoldCreate= bfoldExtrem= ',bfoldStart,infCosmo%bgIni%efold
          stop
       endif
    endif
    

    efold = efoldTune
    
    bgVar(1:fieldNum) = fieldTune(1:fieldNum)
    bgVar(fieldNum+1:2*fieldNum) = velocityTune(1:fieldNum)
    
    call easydverk(neqsbg,bg_field_dot_coupled,efold,bgVar,efoldStart,tolBgEvol)

    efoldTune = efold
    fieldTune(1:fieldNum) = bgVar(1:fieldNum)
    velocityTune(1:fieldNum) = bgVar(fieldNum+1:2*fieldNum)


!set the quantum initial conditions for the perturbations at bfoldStart    

    call pert_creation_standard(tensStart,tensDotStart,bfoldStart,infCosmo,kmpc)

   
!evolve bg + tens pert from mode creation to end of inflation

    allVar(1) = real(tensStart)
    allVar(2) = aimag(tensStart)
    allVar(3) = real(tensDotStart)    
    allVar(4) = aimag(tensDotStart)

    allVar(5:4+fieldNum) = bgVar(1:fieldNum)
    allVar(5+fieldNum:4+2*fieldNum) = bgVar(fieldNum+1:2*fieldNum)

    bfold = bfoldStart

    do while (bfold.lt.infCosmo%bfoldEnd)
    
       bfoldNext = min(bfold + bfoldStep, infCosmo%bfoldEnd)

       call easydverk(neqs,pert_tens_bgdot,bfold,allVar,bfoldNext,tolEvol,cosmoData)

       tens =  cmplx(AllVar(1),AllVar(2))
       tensDot = cmplx(AllVar(3), AllVar(4))       

!for test
       if (dump_modes) then
          field = allVar(5:4+fieldNum)
          velocity = allVar(5+fieldNum:4+2*fieldNum)
          hubbleSquare = hubble_parameter_square(field,velocity,.true.)
          call livewrite('kmpc.dat',bfold, &
               kmpc*exp(-bfold)*exp(infCosmo%efoldEndToToday) &
               *exp(-lnMpcToKappa))
          call livewrite('tens.dat',bfold,abs(tens),abs(real(tens)),abs(aimag(tens)))
          call livewrite('hubble.dat',bfold,sqrt(hubbleSquare))
       endif
    
    enddo
  
    pert_tens_bgevol = tens

  end function  pert_tens_bgevol




  subroutine pert_tens_bgdot(neqs,bfold,allVar,allVarDot,cosmoData)       
    use infprec, only : transfert
    use infbg, only : bg_field_dot_coupled
    use infbgfunc, only : slowroll_first_parameter, hubble_parameter_square
    use inftorad, only : lnMpcToKappa
    implicit none
        
    integer :: neqs    
    real(kp) :: bfold
    type(transfert), optional, intent(inout) :: cosmoData    
    real(kp), dimension(neqs) :: allVar
    real(kp), dimension(neqs) :: allVarDot

    real(kp), dimension(fieldNum) :: field, velocity
    real(kp), dimension(2*fieldNum) :: bgVar, bgVarDot
    real(kp), dimension(4) :: pertVar, pertVarDot

    real(kp) :: kmpc, efoldEnd, efoldEndToToday
    real(kp) :: efold, epsilon1, hubble, kphysOverHubble

    kmpc = cosmoData%real1
    efoldEnd = cosmoData%real2
    efoldEndToToday = cosmoData%real3
        
!equation of motion for the background
    bgVar(1:2*fieldNum) = allVar(5:4+2*fieldNum)

    efold = bfold + efoldEnd
    call bg_field_dot_coupled(2*fieldNum,efold,bgVar,bgVarDot)  

!need the background for the equations of the perturbations
    field = bgVar(1:fieldNum)
    velocity = bgVar(fieldNum+1:2*fieldNum)
    epsilon1 = slowroll_first_parameter(field,velocity,.true.)
    hubble = sqrt(hubble_parameter_square(field,velocity,.true.))

    kphysOverHubble = (kmpc/hubble)* exp(-bfold + efoldEndToToday &
         - lnMpcToKappa)
!    print *,'kphysOverHubble',kphysOverHubble
!    read(*,*)
!equations of motion for the tensor modes
    pertVar(1:4) = allVar(1:4)

    pertVarDot(1:2) = pertVar(3:4)
    pertVarDot(3:4) = - (kphysOverHubble**2) * pertVar(1:2) &
         - (3d0 - epsilon1)*pertVar(3:4)

!set the return value 
    allVarDot(1:4) = pertVarDot(1:4)
    allVarDot(5:4+2*fieldNum) = bgVarDot(1:2*fieldNum)

  end subroutine pert_tens_bgdot



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!scalar modes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  function power_spectrum_scal(infCosmo,kmpc)
!computes 1/(2pi^2) * <zeta,zeta>,  <zeta,S1>, <zeta,S2>, <S1,zeta>, <S1,S1> etc...
!the evolved variables already include the k^3/2 factor.   
    use inftorad, only : infhubblexit, hubble_splinexit
    use infinout, only : livewrite

    implicit none
    real(kp), dimension(scalNum,scalNum) :: power_spectrum_scal
    type(inftoradcosmo), intent(in) :: infCosmo
    real(kp), intent(in) :: kmpc

    type(infhubblexit) :: atHkExit
    integer :: vacuumNum
    complex(kp) :: powerFrom
  
    real(kp) :: pi
    complex(kp), dimension(scalNum,fieldNum) :: obsPertFrom
   
    integer :: i,j,k,whichInVac


    obsPertFrom = 0._kp

!    vacuumNum = fieldNum

!when the dilaton is not coupled
    vacuumNum = matterNum

    pi=acos(-1._kp)

    do whichInVac=1,vacuumNum
       obsPertFrom(:,whichInVac) = pert_scal_bgevol(infCosmo,kmpc,whichInVac)
    enddo

!the two-points correlation functions
!mind the modulus --> cross correlations are complex in Fourier Space.
    
    do j=1,scalNum
       do i=1,scalNum
          powerFrom = 0._kp
          do k=1,vacuumNum
             powerFrom = powerFrom + conjg(obsPertFrom(i,k))*obsPertFrom(j,k)             
          enddo
          power_spectrum_scal(i,j) = (0.5_kp/pi/pi)*real(powerFrom)
       enddo
    enddo

    if (dump_spectra) then      
       call livewrite('powerzetaK.dat',kmpc &
            ,power_spectrum_scal(scalNum,scalNum))
       atHkExit = hubble_splinexit(infCosmo,kmpc)
       call livewrite('powerzetaN.dat',atHkExit%bfold &
            ,power_spectrum_scal(scalNum,scalNum))
    endif

    if (display) then
       write(*,*)'power_spectrum_scal:'
       write(*,*)'kmpc= Ps= '&
            ,kmpc,power_spectrum_scal(scalNum,scalNum)
    endif

    if (dump_modes) then
       write(*,*)'scalar evolution dumped!'       
       read(*,*)
    endif
    
      
  end function power_spectrum_scal





   function pert_scal_bgevol(infCosmo,kmpc,whichInVacuum)
!this evolves the bg and scalar = k^3/2 *(fieldPert,Psi) simultaneously and
!return the final values of . zeta comoving
!curvature perturbation, Psi Bardeen potential, Q Mukhanov variable.

    use inftools, only : easydverk
    use infprec, only : transfert
    use infdilaton, only : conformal_factor_square, conformal_first_gradient
    use infsigma, only : metric, deriv_metric    
    use infpotential, only : potential, deriv_potential
    use infbg, only : operator(/=)
    use infbg, only : bg_field_dot_coupled
    use infbgfunc, only : slowroll_first_parameter, slowroll_second_parameter    
    use infbgfunc, only : hubble_parameter_square
    
    use inftorad, only : bfold_hubble_fraction
    use infinout

    implicit none

    complex(kp), dimension(scalNum) :: pert_scal_bgevol
    type(inftoradcosmo), intent(in) :: infCosmo
    real(kp), intent(in) :: kmpc
    integer, intent(in) :: whichInVacuum
   
    type(transfert) :: cosmoData
 
!not intialized
    real(kp) :: efoldTune, hubbleTune, hubble
    real(kp), dimension(fieldNum) :: fieldTune, velocityTune
    type(infbgphys) :: bgPrevious
    common /bgScalSave/ efoldTune, fieldTune, velocityTune, bgPrevious   
    save /bgScalSave/
!$omp threadprivate(/bgScalSave/)

    real(kp) :: bfold, efold, efoldStart, bfoldStart, bfoldNext
  
!output files number of steps
    integer, parameter :: bfoldDataNum = 5000
    real(kp) :: bfoldStep

!field, fieldDot real; pert, pertDot complex.

!accuracy for the forward integration of background only and both
!background and perturbations (any backward (unstable) integration is
!avoided). Since sub-Hubble perturbations oscillate, funny low
!accuracy on pertevol lead to funny long integration time
    real(kp), parameter :: tolBgEvol = tolkp
    real(kp), parameter :: tolEvol = 1e-11

    integer, parameter :: neqs = 4*scalNum + 2*fieldNum
    integer, parameter :: neqsbg = 2*fieldNum

    complex(kp) :: pertStart, pertDotStart    
    complex(kp), dimension(2) :: bardeensStart
    complex(kp), dimension(scalNum) :: scalStart, scalDotStart
    complex(kp), dimension(scalNum) :: scal,scalDot

    real(kp) :: epsilon1
    real(kp) :: kphysOverHubbleStart, kinetic, kineticDot

!observable quantites, total and partial comoving curvature perturbations
    complex(kp) :: zeta, zetaJF
    complex(kp), dimension(fieldNum) :: mukhaScal,zetaScal !,entro
    complex(kp), dimension(matterNum) :: zetaMatterDensJF

    real(kp) :: sigmaDot
    real(kp), dimension(fieldNum) :: unitDot
    real(kp), dimension(fieldNum) :: field,velocity,fieldDot
    real(kp), dimension(2*fieldNum) :: bgVar
    real(kp), dimension(fieldNum,fieldNum) :: metricVal
    real(kp), dimension(fieldNum,fieldNum,fieldNum) :: metricDeriv
!field and scal together
    real(kp), dimension(neqs) :: allVar

    integer :: i
    character(len=15) :: strgWhich, strgCount


!bugbuster
    complex(kp), dimension(2) :: bardeenTest
    logical, parameter :: doJordanFrame = .false.
   
    if ((whichInVacuum.gt.fieldNum).or.(whichInVacuum.le.0)) then
       stop 'pert_scalar_bgevol: no such scalar perturbations'
    endif


!to pass kmpc to differential equations throught dverk
    cosmoData%real1 = kmpc
    cosmoData%real2 = infCosmo%bgEnd%efold
    cosmoData%real3 = infCosmo%efoldEndToToday
   
!find the bfold of creation
! solution of N + ln[kappa H(N)] = ln(kmpc) -ln(aend/a0) -Nc -ln[(k/aH)_creation]    
    kphysOverHubbleStart = kphysOverHubbleCreate
    bfoldStart = bfold_hubble_fraction(kmpc,infCosmo,kphysOverHubbleStart)
    
    if (display) then
       write(*,*)
       write(*,*)'scalar modes: kmpc = ',kmpc
       write(*,*)'bfold quantum = ',bfoldStart
       write(*,*)'mode in vacuum: ',whichInVacuum       
       write(*,*)
    endif

    if (dump_modes) then
       if (bfoldDataNum.le.1) stop 'pert_scal_bgevol: 1 data point?'
       bfoldStep = (infCosmo%bfoldEnd - bfoldStart)/real(bfoldDataNum-1)       
       write(strgWhich,*) whichInVacuum       
       call delete_file('zeta2_'//trim(adjustl(strgWhich))//'.dat')
       call delete_file('psi_'//trim(adjustl(strgWhich))//'.dat')
       call delete_file('zeta2JF_'//trim(adjustl(strgWhich))//'.dat')
       do i=1,fieldNum
          write(strgCount,*) i
          call delete_file('scal_'//trim(adjustl(strgCount))// &
                  '_'//trim(adjustl(strgWhich))//'.dat')
!          call delete_file('entro_'//trim(adjustl(strgCount))// &
!                  '_'//trim(adjustl(strgWhich))//'.dat')         
          call delete_file('zetascal_'//trim(adjustl(strgCount))// &
               '_'//trim(adjustl(strgWhich))//'.dat')         
       enddo
    else
!turbo settings, let's go to the end of inflation
       bfoldStep = 1._kp/epsilon(1._kp)
    endif
   

!test if the bg model change
    if (infCosmo%bgIni /= bgPrevious) then
       efoldTune = infCosmo%bgIni%efold
       fieldTune = infCosmo%bgIni%field
       velocityTune = infCosmo%bgIni%fieldDot * infCosmo%bgIni%hubble
    endif

    bgPrevious = infCosmo%bgIni
   
    efoldStart = bfoldStart + infCosmo%bgEnd%efold

!only forward integration allowed
    if (efoldStart.lt.efoldTune) then
       efoldTune = infCosmo%bgIni%efold
       fieldTune = infCosmo%bgIni%field
       velocityTune = infCosmo%bgIni%fieldDot * infCosmo%bgIni%hubble
       if (efoldStart.lt.infCosmo%bgIni%efold) then
          write(*,*) 'pert_scalar_bgevol: kmpc =',kmpc
          write(*,*) 'bfoldCreate= bfoldExtrem= ',bfoldStart,infCosmo%bgIni%efold
          stop
       endif
    endif
    

!evolve the bg from 0 to bfold of mode creation
    efold = efoldTune
    bgVar(1:fieldNum) = fieldTune(1:fieldNum)
    bgVar(fieldNum+1:2*fieldNum) = velocityTune(1:fieldNum)

    call easydverk(neqsbg,bg_field_dot_coupled,efold,bgVar,efoldStart,tolBgEvol)
   
!update varTune to the bfold of mode creation (the bg will be computed
!from here next time)
    efoldTune = efold
    fieldTune(1:fieldNum) = bgVar(1:fieldNum)
    velocityTune(1:fieldNum) = bgVar(fieldNum+1:2*fieldNum)
    hubbleTune = sqrt(hubble_parameter_square(fieldTune,velocityTune,.true.))

!the field initial quantum states have to be normalised with respect
!to their kinetic term (done by calling pert_creation_kinetic)
    metricVal = metric(fieldTune)
    metricDeriv = deriv_metric(fieldTune)

    kinetic = metricVal(whichInVacuum,whichInVacuum)
    kineticDot = dot_product(metricDeriv(whichInVacuum,whichInVacuum,:) &
         ,velocityTune/hubbleTune)

!    print *,'kinteic',kinetic,kineticDot

    call pert_creation_kinetic(pertStart,pertDotStart &
         ,bfoldStart,kinetic,kineticDot,infCosmo,kmpc)
!    print *,'kinetic'
!    print *,'pertStart',pertStart
!    print *,'pertDotStart',pertDotStart

!    call pert_creation_standard(pertStart,pertDotStart,bfoldStart,infCosmo,kmpc)
!    print *,'standard'
!    print *,'pertStart',pertStart
!    print *,'pertDotStart',pertDotStart

!reset everyone
    scalStart = 0._kp
    scalDotStart = 0._kp

!there is a sqrt(2) factor between scalar and tensor modes.
    scalStart(whichInVacuum) = pertStart /sqrt(2._kp)
    scalDotStart(whichInVacuum) = pertDotStart /sqrt(2._kp)


!The initial value of the Bardeen potential is fixed by the constraint
!equuations, i.e. momemtum and energy conservation. It is a bad idea
!to use these equations during the evolution since they are singular
!at Hubble exit and epsilon=1, but are fine initially, for quantum
!wavelength well below the Hubble radius
    bardeensStart = bardeen_bardeen_dot(bfoldStart,fieldTune,velocityTune &
         ,scalStart(1:fieldNum),scalDotStart(1:fieldNum),cosmoData)

    scalStart(scalNum) = bardeensStart(1)      
    scalDotStart(scalNum) = bardeensStart(2)

    
!evolve bg + tens pert from mode creation to end of inflation
!pert variables and their derivatives
    do i=1,scalNum       
       allVar(2*i-1) = real(scalStart(i))
       allVar(2*i) = aimag(scalStart(i))
       allVar(2*scalNum + 2*i-1) = real(scalDotStart(i))    
       allVar(2*scalNum + 2*i) = aimag(scalDotStart(i))
    enddo
   
!bg field and velocity
    allVar(4*scalNum+1:4*scalNum+2*fieldNum) = bgVar(1:2*fieldNum)
    
    bfold = bfoldStart
  
    do while (bfold.lt.infCosmo%bfoldEnd)
    
       bfoldNext = min(bfold + bfoldStep, infCosmo%bfoldEnd)

!       print *
!       print *,'allVar',allVar
!       print *

       call easydverk(neqs,pert_scalar_bgdot,bfold,allVar,bfoldNext,tolEvol,cosmoData)

       field = allVar(4*scalNum+1:4*scalNum+fieldNum)
       metricVal = metric(field)       
       velocity = allVar(4*scalNum+fieldNum+1:4*scalNum+2*fieldNum)
       hubble = sqrt(hubble_parameter_square(field,velocity,.true.))
       fieldDot = velocity/hubble
       epsilon1 = slowroll_first_parameter(field,velocity,.true.)       
       sigmaDot = sqrt(2._kp*epsilon1)

       do i=1,scalNum
          scal(i) = cmplx(allVar(2*i-1),allVar(2*i))
          scalDot(i) = cmplx(allVar(2*scalNum+2*i-1),allVar(2*scalNum+2*i))          
       enddo

       unitDot = fieldDot/sigmaDot
       
!curvature perturbation for each field: Psi +
!deltaField/(Dfield/Dbfold)
       do i=1,fieldNum
          mukhaScal(i) = scal(i) + fieldDot(i)*scal(scalNum)
          if (fieldDot(i).ne.0._kp) then
             zetaScal(i) = scal(scalNum) + scal(i)/fieldDot(i)
          else
             zetaScal(i) = 0._kp
          endif
!          entro(i) = scal(i) - unitDot(i) &
!               * dot_product(matmul(metricVal,unitDot),scal(1:fieldNum))
       enddo
       

      
!zeta constant density hypersurface in Jordan Frame for matter fields
       if (doJordanFrame) then
          zetaMatterDensJF(1:matterNum) &
               = curvature_matter_density_JF(field,velocity,scal,scalDot)
          zetaJF = curvature_comoving_JF(field,velocity,scal,scalDot)
       else
!comoving curvature perturbation: Psi + deltaSigma/(DSigma/Dbfold)
          zeta = scal(scalNum) & 
               + (dot_product(matmul(metricVal,fieldDot),scal(1:fieldNum))) &
               /(2._kp*epsilon1)
       endif

!for test (fully inefficient way of writing files)
       if (dump_modes) then                    


          bardeenTest = bardeen_bardeen_dot(bfold,field,velocity &
               ,scal(1:fieldNum),scalDot(1:fieldNum),cosmoData)

!          print *,'Psi',bardeenTest(1)-scal(scalNum)
!          print *,'Psidot',bardeenTest(2)-scalDot(scalNum)
!          print *,'zeta',zeta-(bardeenTest(1)+ (bardeenTest(2) + bardeenTest(1))/epsilon1)
!          print *,'ZetaScal'
!          print *,'zeta -each',zeta &
!               - (dot_product(matmul(metricVal,velocity/hubble) &
!               ,velocity/hubble*zetaScal)/(2.*epsilon1))
!          read(*,*)

          call livewrite('zeta2_'//trim(adjustl(strgWhich))//'.dat' &
               , bfold, abs(conjg(zeta)*zeta),abs(real(zeta)),abs(aimag(zeta)))

          call livewrite('zeta2JF_'//trim(adjustl(strgWhich))//'.dat' &
                  ,bfold,real(conjg(zetaJF)*zetaJF),abs(real(zetaJF)) &
                  ,abs(aimag(zetaJF)))

          do i=1,fieldNum
             write(strgCount,*) i
             call livewrite('scal_'//trim(adjustl(strgCount))// &
                  '_'//trim(adjustl(strgWhich))//'.dat' &
                  , bfold, abs(scal(i)),abs(real(scal(i))), abs(aimag(scal(i))))             
             call livewrite('zetascal_'//trim(adjustl(strgCount))// &
                  '_'//trim(adjustl(strgWhich))//'.dat' &
                  ,bfold,abs(zetaScal(i)),abs(real(zetaScal(i)))&
                  ,abs(aimag(zetaScal(i))))             
          enddo

          call livewrite('psi_'//trim(adjustl(strgWhich))//'.dat' &
               , bfold, abs(scal(scalNum)),abs(real(scal(scalNum))) &
               ,abs(aimag(scal(scalNum))))

       endif
    enddo


!zeta is returned instead from the mukhanov variable Q
!zeta = Q/sigmaDot = Q/(sqrt(2epsilon1))
!rescaled entropic perturbations are returned: S1(2) = entro1(2)/sigmaDot

!    pert_scal_bgevol(1:fieldNum) = entro(1:fieldNum)/sigmaDot
!    pert_scal_bgevol(scalNum) = zeta

    
    if (doJordanFrame) then
       pert_scal_bgevol(1:matterNum) = zetaMatterDensJF
       pert_scal_bgevol(matterNum+1:fieldNum) = scal(matterNum+1:fieldNum)
       pert_scal_bgevol(scalNum) = zetaJF       
    else
       do i=1,fieldNum
          pert_scal_bgevol(i) = zetaScal(i)-zetaScal(entroRef)
       enddo                 
       pert_scal_bgevol(scalNum) = zeta
    endif

  end function  pert_scal_bgevol






  
  subroutine pert_scalar_bgdot(neqs,bfold,allVar,allVarDot,cosmoData) 
    use infprec, only : transfert
    use infsigma, only : metric, metric_inverse, deriv_metric
    use infsigma, only : connection_affine, deriv_connection_affine
    use infpotential, only : potential, deriv_potential, deriv_second_potential
    use infbg, only : bg_field_dot_coupled
    use infbgfunc, only : slowroll_first_parameter, hubble_parameter_square

    use inftorad, only : lnMpcToKappa

    implicit none
        
    integer :: neqs    
    real(kp) :: bfold
    type(transfert), optional, intent(inout) :: cosmoData    
    real(kp), dimension(neqs) :: allVar
    real(kp), dimension(neqs) :: allVarDot

    real(kp), dimension(fieldNum) :: field, velocity, fieldDot
    
    real(kp), dimension(2*fieldNum) :: bgVar, bgVarDot
   
    real(kp) :: kmpc, efoldEnd, efoldEndToToday
    real(kp) :: efold, kphysOverHubble, pot
    real(kp) :: hubbleSquare, hubble, epsilon1


!temp variable for writing the equations of motion...   
    complex(kp) :: bardeenPsi, bardeenPsiDot, bardeenPsiDotDot
    complex(kp), dimension(fieldNum) :: fieldPert,fieldPertDot,fieldPertDotDot

!the coefficients entering in the equations of motion


!the physical quantities entering in the above coefficients
    real(kp), dimension(fieldNum) :: potDeriv, potDerivVec
    real(kp), dimension(fieldNum,fieldNum) :: metricVal, metricInv, potSecond
    real(kp), dimension(fieldNum,fieldNum) :: connectXfieldDot, potSecondVec
    real(kp), dimension(fieldNum,fieldNum) :: connectDerivXfieldDotSquare
    real(kp), dimension(fieldNum,fieldNum) :: metricPotDerivs, Vff
    real(kp), dimension(fieldNum,fieldNum,fieldNum) :: metricDeriv,connect
    real(kp), dimension(fieldNum,fieldNum,fieldNum,fieldNum) :: connectDeriv

    integer :: i,j



!get the wave number from dverk
    kmpc = cosmoData%real1       
    efoldEnd = cosmoData%real2
    efoldEndToToday = cosmoData%real3

!equation of motion for the background (we use the physical velocity)
    bgVar(1:2*fieldNum) = allVar(4*scalNum+1:4*scalNum+2*fieldNum)

    efold = bfold + efoldEnd
    call bg_field_dot_coupled(2*fieldNum,efold,bgVar,bgVarDot)  


!need some background quantities for the perturbations
    field = bgVar(1:fieldNum)
    velocity = bgVar(fieldNum+1:2*fieldNum)
    hubbleSquare = hubble_parameter_square(field,velocity,.true.)
    hubble = sqrt(hubbleSquare)
    fieldDot = velocity/hubble    

!remind this is also: -(DHubble/Dbfold)/Hubble and sigmaDot^2/2
    epsilon1 = slowroll_first_parameter(field,velocity,.true.)    

    kphysOverHubble = (kmpc/hubble) * exp(-bfold + efoldEndToToday &
         - lnMpcToKappa)


!the scalar perturbations (real and imaginary parts)
!    fieldPert(1:fieldNum) = allVar(1:2*fieldNum)
!    bardeenPsi = allVar(2*fieldNum+1:2*scalNum)

!their derivatives
!    fieldPertDot(1:fieldNum) = allVar(2*scalNum+1:2*scalNum+2*fieldNum)
!    bardeenPsiDot = allVar(2*scalNum+2*fieldNum+1:4*scalNum)

    do i=1,fieldNum
       fieldPert(i) = cmplx(allVar(2*i-1),allVar(2*i))
       fieldPertDot(i) = cmplx(allVar(2*scalNum+2*i-1),allVar(2*scalNum+2*i))
    enddo
    bardeenPsi = cmplx(allVar(2*fieldNum+1),allVar(2*scalNum))
    bardeenPsiDot = cmplx(allVar(2*scalNum+2*fieldNum+1),allVar(4*scalNum))



!their equations of motion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! warm up !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    metricVal = metric(field)
    metricInv = metric_inverse(field)
    metricDeriv = deriv_metric(field)

    connect = connection_affine(field)
    connectDeriv = deriv_connection_affine(field)

!all potential related functions are in unit of H^2
    pot = potential(field)/hubbleSquare
    potDeriv = deriv_potential(field)/hubbleSquare
    potSecond = deriv_second_potential(field)/hubbleSquare

    potDerivVec = matmul(metricInv,potDeriv)
    potSecondVec = matmul(metricInv,potSecond)

    do i=1,fieldNum       
       connectXfieldDot(i,:) = matmul(connect(i,:,:),fieldDot)
       do j=1,fieldNum
          connectDerivXfieldDotSquare(i,j) &
               = dot_product(matmul(connectDeriv(i,:,:,j),fieldDot),fieldDot)
          metricPotDerivs(i,j) &
               = dot_product(matmul(metricDeriv(:,:,j),potDerivVec),metricInv(i,:))
          Vff(i,j) = connectDerivXfieldDotSquare(i,j) + potSecondVec(i,j) &
               - metricPotDerivs(i,j)
       enddo
    enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! equations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    fieldPertDotDot = - (3._kp - epsilon1) * fieldPertDot &
         - 2._kp * matmul(connectXfieldDot,fieldPertDot) - matmul(Vff,fieldPert) &
         - kphysOverHubble**2 * fieldPert &
         + 4._kp * fieldDot * bardeenPsiDot - 2._kp * potDerivVec * bardeenPsi

    bardeenPsiDotDot = - (7._kp - epsilon1) * bardeenPsiDot &
         - (2._kp * pot + kphysOverHubble**2) * bardeenPsi &
         - dot_product(potDeriv,fieldPert)
             
!
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! end equations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!specify the return value: glue all things together for integration.

!    allVarDot(1:2*fieldNum) = fieldPertDot(1:fieldNum)
!    allVarDot(2*fieldNum+1:2*scalNum) = bardeenPsiDot
!    allVarDot(2*scalNum+1:2*scalNum+2*fieldNum) = fieldpertDotDot(1:fieldNum)
!    allVarDot(2*scalNum+2*fieldNum+1:4*scalNum) = bardeenPsiDotDot

    do i=1,fieldNum
       allVarDot(2*i-1)=real(fieldPertDot(i))
       allVarDot(2*i) = aimag(fieldPertDot(i))
       allVarDot(2*scalNum+2*i-1) = real(fieldPertDotDot(i))
       allVarDot(2*scalNum+2*i) = aimag(fieldPertDotDot(i))
    enddo

    allVarDot(2*fieldNum+1) = real(bardeenPsiDot)
    allVarDot(2*scalNum) = aimag(bardeenPsiDot)
    allVarDot(2*scalNum+2*fieldNum+1) = real(bardeenPsiDotDot)
    allVarDot(4*scalNum) = aimag(bardeenPsiDotDot)
    
    allVarDot(4*scalNum+1:4*scalNum+2*fieldNum) = bgVarDot(1:2*fieldNum)

  end subroutine pert_scalar_bgdot
  



  
  function bardeen_bardeen_dot(bfold,field,velocity,fieldPert,fieldPertDot,cosmoData)
!return the Bardeen potential and its derivative wrt efold from the
!constraint equations
    use infprec, only : transfert
    use infsigma, only : metric, deriv_metric
    use infpotential, only : deriv_potential
    use infbgfunc, only : slowroll_first_parameter, hubble_parameter_square
    use inftorad, only : lnMpcToKappa

    implicit none
    complex(kp), dimension(2) :: bardeen_bardeen_dot
    real(kp), intent(in) :: bfold
    real(kp), dimension(fieldNum), intent(in) :: field, velocity
    complex(kp), dimension(fieldNum), intent(in) :: fieldPert,fieldPertDot    
    type(transfert), intent(in) :: cosmoData 
    
    
    real(kp) :: epsilon1, hubbleSquare, hubble
    real(kp) :: kmpc, efoldEndToToday, kphysOverHubble
    

    real(kp), dimension(fieldNum) :: potDeriv,metricDerivXFieldDotSquare, fieldDot
    real(kp), dimension(fieldNum) :: metricXFieldDot
    real(kp), dimension(fieldNum,fieldNum) :: metricVal
    real(kp), dimension(fieldNum,fieldNum,fieldNum) :: metricDeriv

    integer :: i


    kmpc = cosmoData%real1           
    efoldEndToToday = cosmoData%real3

    epsilon1 = slowroll_first_parameter(field,velocity,.true.)
    hubbleSquare = hubble_parameter_square(field,velocity,.true.)
    hubble = sqrt(hubbleSquare)

    fieldDot = velocity/hubble

    metricVal = metric(field)
    metricDeriv = deriv_metric(field)
    potDeriv = deriv_potential(field)/hubbleSquare
    
    kphysOverHubble = (kmpc/hubble) * exp(-bfold + efoldEndToToday &
         - lnMpcToKappa)
    
    metricXFieldDot = matmul(metricVal,fieldDot)

    do i=1,fieldNum
       metricDerivXFieldDotSquare(i) = dot_product(matmul(metricDeriv(:,:,i),fieldDot) &
            ,fieldDot)
    enddo

!bardeen potential
    bardeen_bardeen_dot(1) &
         = (1.5_kp * dot_product(metricXFieldDot,fieldPert) &
         + 0.25_kp * dot_product(metricDerivXFieldDotSquare,fieldPert) &
         + 0.5_kp * dot_product(potDeriv,fieldPert) &
         + 0.5_kp * dot_product(metricXFieldDot,fieldPertDot)) &
         / (epsilon1 - kphysOverHubble**2)

!bardeen potential derivative wrt bfold time
    bardeen_bardeen_dot(2) = 0.5_kp * dot_product(metricXFieldDot,fieldPert) &
         - bardeen_bardeen_dot(1)


  end function bardeen_bardeen_dot




  function curvature_matter_density_JF(field,velocity,scal,scalDot)
!this is the curvature on constant density hypersurface for the matter
!fields zetamat = Psi + (delta Rho) / (DRho/Defold), if one matter
!field only. zetamat_i = Psi + (delta Rho_i) / (DField_i/Defold)^2
!otherwise. In case of multiple matter fields, rho_i is not well
!defined due to possible couplings through the matter potential.

    use infdilaton, only : conformal_factor_square, conformal_first_gradient
    use infmatter, only : deriv_matter_potential
    use infbgfunc, only : hubble_parameter_square
    implicit none

    complex(kp), dimension(matterNum) :: curvature_matter_density_JF
    real(kp), dimension(fieldNum), intent(in) :: field, velocity
    complex(kp), dimension(scalNum), intent(in) :: scal, scalDot

    
    complex(kp), dimension(matterNum) :: matterPert, matterPertDot
    complex(kp), dimension(dilatonNum) :: dilatonPert, dilatonPertDot
    complex(kp) :: bardeenPsi, bardeenPhi

    real(kp) :: hubbleSquare, confSquare
    real(kp), dimension(matterNum) :: matter, matterDot,derivMatterPot
    real(kp), dimension(dilatonNum) :: dilaton, confFirstGrad
    real(kp), dimension(fieldNum) :: fieldDot
    
    bardeenPsi = scal(scalNum)
    bardeenPhi = bardeenPsi

    matterPert = scal(1:matterNum)
    matterPertDot = scalDot(1:matterNum)
    dilatonPert = scal(matterNum+1:fieldNum)
    dilatonPertDot = scalDot(matterNum+1:fieldNum)

    hubbleSquare = hubble_parameter_square(field,velocity,.true.)
    fieldDot = velocity/sqrt(hubbleSquare)
    matterDot(1:matterNum) = fieldDot(1:matterNum)
    
    matter = field(1:matterNum)
    dilaton = field(matterNum+1:fieldNum)
    derivMatterPot = deriv_matter_potential(matter)
    confSquare = conformal_factor_square(dilaton)
    confFirstGrad = conformal_first_gradient(dilaton)


    curvature_matter_density_JF = bardeenPsi + bardeenPhi/3._kp &
         - (2._kp/3._kp) * dot_product(confFirstGrad,dilatonPert) &
         - (1._kp/3._kp) * (matterPertDot/matterDot &
         + confSquare * derivMatterPot / hubbleSquare / matterDot**2._kp &
         * matterPert)

  end function curvature_matter_density_JF




  function curvature_comoving_JF(field,velocity,scal,scalDot)
!this is the comoiving curvature perturbation in the Jordan Frame,
!i.e.  zeta = PsiJF + (PsiJF' + HubbleConfJF x PhiJF)/(HubbleConfJF x
!epsilon1JF)

    use infdilaton, only : conformal_first_gradient, conformal_second_gradient
    use infbgfunc, only : slowroll_first_parameter_JF, hubble_parameter_square
    implicit none

    complex(kp) :: curvature_comoving_JF
    real(kp), dimension(fieldNum), intent(in) :: field, velocity
    complex(kp), dimension(scalNum), intent(in) :: scal, scalDot

        
    complex(kp), dimension(dilatonNum) :: dilatonPert, dilatonPertDot
    complex(kp) :: bardeenPsi, bardeenPsiDot, bardeenPhi
    complex(kp) :: bardeenPsiJF, bardeenPsiJFDot, bardeenPhiJF

    real(kp) :: hubbleSquare, epsilon1JF
    real(kp), dimension(fieldNum) :: fieldDot
    real(kp), dimension(dilatonNum) :: dilaton, dilatonDot,confFirstGrad
    real(kp), dimension(dilatonNum,dilatonNum) :: confSecondGrad


    hubbleSquare = hubble_parameter_square(field,velocity,.true.)
    epsilon1JF = slowroll_first_parameter_JF(field,velocity,.true.)

    fieldDot = velocity/sqrt(hubbleSquare)

    dilaton = field(matterNum+1:fieldNum)
    dilatonDot = fieldDot(matterNum+1:fieldNum)
    
    bardeenPsi = scal(scalNum)
    bardeenPhi = bardeenPsi
    bardeenPsiDot = scalDot(scalNum)

    dilatonPert = scal(matterNum+1:fieldNum)
    dilatonPertDot = scalDot(matterNum+1:fieldNum)

      
    confFirstGrad = conformal_first_gradient(dilaton)
    confSecondGrad = conformal_second_gradient(dilaton)
   
    bardeenPsiJF = bardeenPsi - dot_product(confFirstGrad,dilatonPert)
    bardeenPhiJF = bardeenPhi + dot_product(confFirstGrad,dilatonPert)

    bardeenPsiJFDot = bardeenPsiDot - dot_product(confFirstGrad,dilatonPertDot) &
         - dot_product(matmul(confSecondGrad,dilatonDot),dilatonPert)

    curvature_comoving_JF = bardeenPsiJF + (bardeenPsiJFDot &
         + (1._kp + dot_product(confFirstGrad,dilatonDot))*bardeenPhiJF) &
         / ((1._kp + dot_product(confFirstGrad,dilatonDot))*epsilon1JF)

  end function curvature_comoving_JF



 
  subroutine pert_creation_standard(pertIni,pertDotIni,bfoldCreate,infCosmo,kmpc)
!For standard constant kinetic terms.
!mode creation at different times corresponding to a fix physical
!scale with respect to the Hubble radius: k/aH = Cte
!kmpc is the wavenumber today in unit of Mpc^-1. H is the Hubble
!parameter during inflation in unit of kappa. Evolution of the perturbations
!involve the dimensionless quantity: k/(aH)
! k/aH = (k/a0)/(kappa H) x aend/a x a0/aend x kappa
!      = kmpc/kappaH x exp[-N] x a0/aend x (kappa/1Mpc)
!      = kmpc/kappaH x exp[-N] x exp[ln(a0/aend)] x exp[-Nc]
! where Nc = ln[1Mpc/sqrt(8pi)/lPl]   <---- Nc=lnMpcToKappa()
! N number of efold before the end of inflation (efold-efoldEnd)

    use inftorad, only : scaleFactorToday, lnMpcToKappa
    implicit none
   
    real(kp), intent(in) :: kmpc
    type(inftoradcosmo), intent(in) :: infCosmo
    real(kp), intent(in) :: bfoldCreate
    complex(kp), intent(out) :: pertIni, pertDotIni
     
    complex(kp) :: qMode, qModePrimeOverk
   

!set the quantum initial conditions at found bfold, in bfold time
!ex: for grav wav
! h = qmode/a
! h = a0/aend x aend/acreate x qmode/a0
! Dh/Dbfold = a H h' = a H (qmode/a)' = (k/aH) x a0/aend x aend/acreate 
!                                                   x (qmode'/k)/a0 - h

    call quantum_creation(kmpc,qmode,qmodePrimeOverk)

    pertIni = qmode/scaleFactorToday
    pertDotIni = (kphysOverHubbleCreate*qmodePrimeOverk - qmode)/scaleFactorToday

!normalisation of k in Mpc^-1 today    
    pertIni = pertIni*exp(infCosmo%efoldEndToToday - bfoldCreate - lnMpcToKappa)
    pertDotIni = pertDotIni*exp(infCosmo%efoldEndToToday - bfoldCreate - lnMpcToKappa)

  end subroutine pert_creation_standard




  subroutine pert_creation_kinetic(pertIni,pertDotIni,bfoldCreate,kinetic,kineticDot &
       ,infCosmo,kmpc)
!for non-standard kinetic terms K(eta) Dfield Dfield instead of Dfield Dfield
    use inftorad, only : scaleFactorToday, lnMpcToKappa
    implicit none
   
    real(kp), intent(in) :: kmpc
    type(inftoradcosmo), intent(in) :: infCosmo
    real(kp), intent(in) :: bfoldCreate, kinetic, kineticDot
    complex(kp), intent(out) :: pertIni, pertDotIni

    real(kp) :: sqrtKinetic, lnDotSqrtKinetic
    complex(kp) :: qMode, qModePrimeOverk
   
!set the quantum initial conditions at found bfold, in bfold time
!ex: for grav wav, qmod designs the standard normalised quantum mode
! h = qmode/a/K^1/2
! h = a0/aend x aend/acreate x qmode/a0/K^1/2
! Dh/Dbfold = a H h' = a H (qmode/a/K^1/2)'
! =(k/aH) x a0/aend x aend/acreate x (qmode'/k)/a0/K^1/2 - h x (1+1/2 DLn(K)/Dbfold)

    call quantum_creation(kmpc,qmode,qmodePrimeOverk)

    sqrtKinetic = sqrt(kinetic)
    lnDotSqrtKinetic = 0.5_kp * kineticDot/kinetic

    pertIni = qmode /scaleFactorToday /sqrtKinetic
    pertDotIni = (kphysOverHubbleCreate*qmodePrimeOverk - (1._kp + lnDotSqrtKinetic)*qmode) &
         /scaleFactorToday /sqrtKinetic

!normalisation of k in Mpc^-1 today  
    pertIni = pertIni*exp(infCosmo%efoldEndToToday - bfoldCreate - lnMpcToKappa)
    pertDotIni = pertDotIni*exp(infCosmo%efoldEndToToday - bfoldCreate - lnMpcToKappa)

  end subroutine pert_creation_kinetic




  subroutine quantum_creation(k,modeIni,modePrimeIniOverFreq)
!creates the initial quantum mode for a given k. The Bogoliubov
!coefficients are alpha and beta such as |alpha|^2 - |beta|^2 = 1. The
!creation time is not specified here and Prime denotes derivative with
!respect to the conformal time. Calling this function for all k at a
!same given conformal time with alpha=1,beta=0 is a Bunch-Davies vacuum.
!This is for rescaled modes in Fourier space: k^{3/2}*mode and
!in unit hbar=c=1, k is in unit of what you wish (cause rescaled modes).

    implicit none
    real(kp), intent(in) :: k
    complex(kp), intent(out) :: modeIni, modePrimeIniOverFreq
    real(kp), parameter :: pi = 3.141592653589793238
    complex(kp), parameter :: alpha = (1._kp,0._kp)
    complex(kp), parameter :: beta = (0._kp,0._kp)

  
    modeIni = k*(alpha + beta)
    modePrimeIniOverFreq = - k*cmplx(0._kp,1._kp)*(alpha - beta)
!    print *,'qmode qmodeDot',modeIni,modePrimeIniOverFreq
  end subroutine quantum_creation


  
end module infpert
