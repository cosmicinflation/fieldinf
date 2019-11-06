!example advanced usage of fieldinf

!integrate power law inflation; which has to end by instability; and
!set up a spline of the power spectra. Checkout fieldmain for the
!simplest usage and more pedagogical comments

program gurumain

  use fieldprec, only : kp
  use infbgmodel
  use infbg 
  use inftorad
  use infpert
  use infbounds, only : field_stopinf
  use infio, only : livewrite, delete_file
  use infpowspline, only : set_power_scal_spline, set_power_tens_spline
  use infpowspline, only : splineval_power_scal, splineval_power_tens

  implicit none


  type(infbgparam) :: infParam
  type(infbgphys) :: bgIni, bgEnd, bgObs
  type(inftoradcosmo) :: infCosmo
  type(infhubblexit) :: atHkExit

  logical :: setupDone
  real(kp) :: lnReheat

  real(kp), dimension(2) :: fieldstop

!store infbgphys data all along the trajectory
  type(infbgdata), pointer :: ptrBgData => null()
  type(infbgdata), pointer :: ptrRun => null()
!number of evaluation we want
  integer, parameter :: nBgData = 1000

  real(kp) :: efold, bfold, phi, eps1, eps2, hubble, ns

!perturbation modes in Mpc^-1
  real(kp) :: kmpc, lnkmpc, lnKmpcMin, lnKmpcMax
  real(kp), dimension(:), allocatable :: lnkmpcvec

!power spectra
  real(kp) :: powerTens, powerZeta
  real(kp), dimension(scalNum,scalNum) :: powerScal

  integer :: i, nkmpcvec, nkmpc



!initialization
  infParam%consts = 0.

!model name
  infParam%name = 'powlaw'

!power law inflation has 3 params V = c1^4 exp[-c2 F] + the field
!value at which inflation stop
!c1
  infparam%consts(1) = 1e-4
!c2
  infparam%consts(2) = 0.02

!in planck unit
  infParam%consts(matterParamNum) = 10


!initial field values (in Planck unit). If zero, triggers an automatic
!guess from slow-roll integration using libaspic; but here I want it
!to start at 1 for instance
  infParam%matters = 1._kp

!update parameters + consistency checks
  setupDone = set_infbg_param(infParam)
  if (.not.setupDone) stop 'model parameters incorrect'

!set and return all initial physical quantities in bgIni

  bgIni = set_infbg_ini(infParam)

  write(*,*)'Initial Conditions'
  call print_infbgphys(bgIni,'bgIni=')
  write(*,*)

  write(*,*)'INITIAL CONDITIONS SET!'



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!background integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!return a two dimensional array fieldstop(1) is the field value at
!which inflation stops, fieldstop(2) is either -1 or +1 if this value
!is a minimum of maximum
  fieldstop = field_stopinf(infParam)

!integrate the background till the end of inflation, here defined by
!field=fieldstop(1) and this is a maximal field value, the last entry
!is .true. (advanced usage).

!Moreover, we store all the infbgphys data in a chain list starting at
!the address pointed by ptrBgData, with nBgData points. Because we
!have set bgIni by hand above, our integration could be very long; so
!inputing bgObs here returns a infbgphys at a point which is at the
!boundary to be observable. Starting a new integration from BgObs
!instead of bgIni will yield *exactly* the same observable predictions
!bgObs is evaluated at 120 efolds before the end of inflation, you can
!change this value in infbgmodel.f90; as well the maximum number of
!e-folds explorable.

  bgEnd = bg_field_evol(bgIni,nBgData,bgObs,ptrBgData,fieldstop(1) &
       ,(fieldstop(2)==1._kp))

  write(*,*)'End of Inflation'
  call print_infbgphys(bgEnd,'bgEnd=')
  write(*,*)

  write(*,*)'Barely Observable (120 e-folds before the end)'
  call print_infbgphys(bgObs,'bgObs=')
  write(*,*)


  write(*,*)'BACKGROUND INTEGRATED'  

!let's dump the field trajectory and all that into a file

  call delete_file('guruindex.dat')
  call delete_file('guruphi.dat')
      
  ptrRun => ptrBgdata
  do while (associated(ptrRun))
     efold = ptrRun%bg%efold
     bfold = ptrRun%bg%efold - bgEnd%efold
     phi = ptrrun%bg%field(1)
     eps1 = ptrRun%bg%epsilon1
     eps2 = ptrRun%bg%epsilon2

     ns = 1._kp - 2._kp*eps1 - eps2

     call livewrite('guruindex.dat',bfold,ns,1._kp-infParam%consts(2)**2)
     call livewrite('guruphi.dat',bfold,phi)

     ptrRun => ptrRun%ptr

  enddo
     
  ptrRun => null()

  write(*,*)'BACKGROUND DATA DUMPED'




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Linear perturbations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  lnReheat = 0._kp
  
  infCosmo = set_inftorad_cosmo(infParam,bgIni,bgEnd,lnReheat)

!prepare for the perturbation
  call set_infpert_ini(infCosmo,ptrBgData)


!let's set a spline of the primordial power spectra over 10 points
  nkmpcvec = 10
  lnKmpcMin = -18._kp  
  lnKmpcMax = 0._kp

  allocate(lnkmpcvec(nkmpcvec))
  do i = 1, nkmpcvec
     lnkmpcvec(i) = lnkmpcMin + real(i-1)*(lnkmpcMax - lnkmpcMin)/real(nkmpcvec - 1)
  enddo

!set the spline
  call set_power_scal_spline(infCosmo,lnkmpcvec)
  call set_power_tens_spline(infCosmo,lnkmpcvec)

  deallocate(lnkmpcvec)

  write(*,*)'POWER SPECTRA SPLINES COMPUTED'
  write(*,*)

!check the spline, input any value of kmpc

  kmpc = 0.05
  powerScal = splineval_power_scal(kmpc)
  powerTens = splineval_power_tens(kmpc)

  write(*,*)'kmpc= Pzeta= Ph= ',kmpc,powerScal(scalNum,scalNum),powerTens
  write(*,*)

!or dumps the whole spectra

  call delete_file('gurupow.dat')

!number of modes (spline being fast, we can get a lot)
  nkmpc = 10000

  do i=1,nkmpc-1
     lnkmpc =  lnkmpcMin + real(i-1)*(lnkmpcMax - lnkmpcMin)/real(nkmpc - 1)
     kmpc = exp(lnkmpc)
     
     powerScal = splineval_power_scal(kmpc)
     powerZeta = powerScal(scalNum,scalNum)

     powerTens = splineval_power_tens(kmpc)

     call livewrite('gurupow.dat',kmpc,powerZeta,powerTens)

  enddo
 

end program gurumain
