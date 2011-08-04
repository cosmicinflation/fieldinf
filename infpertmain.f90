!simple driver to get the scalar and tensor power spectra at the end
!of inflation.

program infpertmain
  use infprec, only : kp
  use infbgmodel
  use infbg
  use inftorad
  use infpert
  use infbgspline
  use infpowspline
  use infinout

  implicit none
  
  integer :: inum,i
  real(kp) :: bfold, efold
  real(kp) :: bfoldExit, hubbleExit,kphysOverHubbleExit
  
  type(infbgparam) :: infParam
  type(infbgphys) :: bgIni, bgObs,bgEnd 
  type(infbgdata), pointer :: ptrBgdata => null()
  type(infbgdata), pointer :: ptrRun => null()
  type(infhubblexit) :: atHkExit,atHstarExit
  type(inftoradcosmo) :: infCosmo

  real(kp) :: pi,error
  real(kp) :: powerTens, powerTensSpline
  real(kp) :: kmpc, lnkmpc,lnkmpcMin,lnkmpcMax
  real(kp), dimension(:), allocatable :: lnkmpcvec
  integer :: nkmpc,nkmpcvec
  real(kp), dimension(scalNum,scalNum) :: powerScal, powerScalSpline
 
  real(kp) :: hubble, epsilon1,epsilon1JF
  real(kp), dimension(fieldNum) :: field, fieldDot
  real(kp), dimension(fieldNum,fieldNum) :: metricVal,metricInv 

  real(kp) :: kstar,const,eps1star,eps2star,Hstar,bfoldStar
  real(kp) :: powerZetaSlowRoll, powerTensSlowRoll
  real(kp), dimension(fieldNum) :: fieldStar
  real(kp), dimension(fieldNum) :: fieldDotStar

  real(kp) :: lnReheatCorr = 0._kp

  logical :: setupDone = .false.
  logical :: check_spline = .false.
  logical :: slowroll = .false.

  pi = acos(-1._kp)


!the model parameters
  infParam%name = 'largef'
  infParam%consts(1) = 1.
  infParam%consts(2) = 2.
!  infParam%consts(3) = 4.623622988896771E-004
    infParam%consts(3) = 1e-2 * sqrt(8._kp*pi)
  infParam%consts(4) = 0.
  infParam%consts(5) = -1.
  infParam%consts(6)=0.

  infParam%consts(7) = 0.1_kp
  infParam%consts(8) = 1000._kp

  infParam%consts(1) = infParam%consts(3) * sqrt(2._kp*pi)

  infParam%conforms(1) = 1.
  
  infParam%matters(1) = 0.*infParam%consts(3)



!initial conditions for the background   
  setupDone = set_infbg_param(infParam)
  print *,'Setup done',setupDone
  bgIni = set_infbg_ini(infParam)


!evolves the bg
  print *,'bgIni',bgIni
  print *
  bgEnd = bg_field_evol(bgIni,1000,bgObs,ptrBgdata)

  print *,'bgEnd',bgEnd
  print *
 
!set the pert model and spline the bgdata     
  call set_infbg_spline(bgEnd,ptrBgdata)
  infCosmo = set_inftorad_cosmo(infParam,bgObs,bgEnd,lnReheatCorr)
  print *,'infCosmo',infCosmo
  print *
  read(*,*)


!just for testing the storage and background splining
  if (associated(ptrBgdata)) then

     call delete_file('splinhubble.dat')
     call delete_file('splinepsilons.dat')
     call delete_file('splinfields.dat')
     call delete_file('splinfieldsdot.dat')
     call delete_file('kmpc.dat')
     
     error=0d0
     ptrRun => ptrBgdata
     do while (associated(ptrRun))

        efold = ptrRun%bg%efold
        bfold = ptrRun%bg%efold - bgEnd%efold

        hubble = splineval_hubble_parameter(bfold)
        epsilon1 = splineval_epsilon1(bfold)
        epsilon1JF = splineval_epsilon1JF(bfold)
        field = splineval_field(bfold)
        fieldDot= splineval_fielddot(bfold)

        metricVal = metric(field)
!        metricInv = metric_inverse(field)

        call livewrite('splinhubble.dat',bfold,hubble,ptrRun%bg%hubble)
        call livewrite('splinepsilons.dat',bfold,epsilon1,epsilon1JF)
!        call livewrite('splinfields.dat',efold,field(1),field(2),field(3))
!        call livewrite('splinfieldsdot.dat',efold,fieldDot(1),fieldDot(2),fieldDot(3))

        
        error = max(abs(1d0 - ptrRun%bg%hubble/hubble),error)
        inum = inum + 1
        ptrRun => ptrRun%ptr     
     enddo
     ptrRun => null()
     print *
     print *,'error on hubble spline = ',error
     print *
!     call free_infbg_data(ptrBgdata)
!     read(*,*)
  endif



!dump the power spectra
!-14
!-12 18

  kphysOverHubbleExit = 1d0
  lnKmpcMin = -18d0  
  lnKmpcMax = 18d0

!  lnKmpcMin=-30d0
!  lnKmpcMax=-20d0

  nkmpc = 101


  kmpc = exp(lnkmpcMin)


!setup the P(k) splines for further evaluation
  if (check_spline) then
     nkmpcvec = 12
  
     allocate(lnkmpcvec(nkmpcvec))
     do i = 1, nkmpcvec
        lnkmpcvec(i) = lnkmpcMin + real(i-1)*(lnkmpcMax - lnkmpcMin)/real(nkmpcvec - 1)
     enddo

     call set_power_scal_spline(infCosmo,lnkmpcvec)
     call set_power_tens_spline(infCosmo,lnkmpcvec)
     deallocate(lnkmpcvec)

     print *,'spline computed on nkmpc=',nkmpc,'points'
     read(*,*)
  endif

!check the first order slow-roll expansion of P(k)
  if (slowroll) then
     const= 0.577215664901 + log(2.) - 2.
     kstar=0.05
     atHstarExit = hubble_splinexit(infCosmo,kstar)
     eps1star = atHstarExit%epsilon1

     Hstar = atHstarExit%Hubble
     bfoldStar = atHstarExit%bfold
     print *,'const = ',const
     print *,'bfoldStar = ',bfoldStar
     fieldStar=splineval_field(bfoldStar)
     fieldDotStar=splineval_fielddot(bfoldStar)

     eps2star=slowroll_second_parameter(fieldStar,fieldDotStar,.false.)
     print *,'epsilon1*= epsilon2*= ',eps1star,eps2star
     print *

     do i = 1, nkmpc-1
        lnkmpc =  lnkmpcMin + real(i-1)*(lnkmpcMax - lnkmpcMin)/real(nkmpc - 1)

        kmpc = exp(lnkmpc)
     
        powerZetaSlowRoll = (Hstar**2/(8.*pi**2*eps1star))*(1.-2.*(const+1.)*eps1star &
             - const*eps2star - (2.*eps1star + eps2star)*log(kmpc/kstar))
        
        powerTensSlowRoll = (2*Hstar**2/(pi**2))*(1.-2.*(const+1)*eps1star &
             - 2.*eps1star*log(kmpc/kstar))

        atHkExit = hubble_splinexit(infCosmo,kmpc)

        call livewrite('powerzetaSRN.dat',atHkExit%bfold,powerZetaSlowRoll)
        call livewrite('powerhSRN.dat',atHkExit%bfold,powerTensSlowRoll)
        call livewrite('powerzetaSRK.dat',kmpc,powerZetaSlowRoll)
        call livewrite('powerhSRK.dat',kmpc,powerTensSlowRoll)

     enddo

     print *,'slow roll P(k) dumped with pivot kstar = ',kstar

  endif


!full numerical P(k) 
  do i = 1, nkmpc-1
     lnkmpc =  lnkmpcMin + real(i-1)*(lnkmpcMax - lnkmpcMin)/real(nkmpc - 1)

     kmpc = exp(lnkmpc)

     print *,'lnkmpc= kmpc=',lnkmpc,kmpc

     atHkExit = hubble_splinexit(infCosmo,kmpc)
    
     print *,'At kmpc exit: H= N= ',atHkExit%hubble,atHkExit%bfold
     print *
     
     powerScal = power_spectrum_scal(infCosmo,kmpc)
     powerTens = power_spectrum_tens(infCosmo,kmpc)

!spline evaluation at the wanted kmpc
     if (check_spline) then
        powerScalSpline = splineval_power_scal(kmpc)
        powerTensSpline = splineval_power_tens(kmpc)
     endif

!(de)comment as required     
     call livewrite('powerzetaN.dat',atHkExit%bfold, powerScal(scalNum,scalNum) &
          , powerScalSpline(scalNum,scalNum))
     call livewrite('powerent1N.dat',atHkExit%bfold, powerScal(1,1), powerScalSpline(1,1))
!     call livewrite('powerent2N.dat',atHkExit%bfold, powerScal(2,2), powerScalSpline(2,2))
!     call livewrite('powerent3N.dat',atHkExit%bfold, powerScal(3,3), powerScalSpline(3,3))
     call livewrite('powerhN.dat',atHkExit%bfold,powerTens)
!     call livewrite('powerzetaent2N.dat',atHkExit%bfold,powerScal(scalNum,2) &
!          ,powerScalSpline(scalNum,2))
!     call livewrite('powerzetaent3N.dat',atHkExit%bfold,powerScal(scalNum,3) &
!          ,powerScalSpline(scalNum,3))
!     call livewrite('powerent1ent2N.dat',atHkExit%bfold,powerScal(1,2) &
!          ,powerScalSpline(1,2))

     call livewrite('kmpcN.dat',kmpc,atHkExit%bfold)

     call livewrite('powerzetaK.dat',kmpc,powerScal(scalNum,scalNum) &
          , powerScalSpline(scalNum,scalNum))
     call livewrite('powerent1K.dat',kmpc,powerScal(1,1), powerScalSpline(1,1))
!     call livewrite('powerent2K.dat',kmpc,powerScal(2,2), powerScalSpline(2,2))
!     call livewrite('powerent3K.dat',kmpc,powerScal(3,3),powerScalSpline(3,3))
     call livewrite('powerhK.dat',kmpc,powerTens)
!     call livewrite('powerzetaent2K.dat',kmpc,powerScal(scalNum,2) &
!          ,powerScalSpline(scalNum,2))
!     call livewrite('powerzetaent3K.dat',kmpc,powerScal(scalNum,3) &
!          ,powerScalSpline(scalNum,3))
!     call livewrite('powerent2ent3K.dat',kmpc,powerScal(2,3) &
!          ,powerScalSpline(2,3))
  enddo

   
!free the spline data
  call free_infbg_spline()
  call free_power_scal_spline()
  call free_power_tens_spline()
 
end program infpertmain


