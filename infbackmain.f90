program infbackmain
  use infprec, only : kp
  use infsric
  use infbgmodel
  use infbg
  use infinout
  use infbounds
 
!  use rgisr

  implicit none

   
  real(kp), dimension(fieldNum) :: field, fieldDot
  type(infbgdata), pointer :: ptrToBgdata => null()
  type(infbgdata), pointer :: ptrRun => null()
  
  type(infbgparam) :: infParam
  type(infbgphys) :: infIni, infObs, infEnd
 
  integer :: inum = 0

  real(kp) :: matter,efold,efoldGo,efoldFin,hubble
  real(kp) :: epsilon1,epsilon1JF, epsilon1SR, epsilon2SR
  real(kp) :: epsmax, alpha, beta, gamma, p, lambda

  real(kp) :: ricci, ricciOverH2, x, xend, xini
  real(kp), dimension(dilatonNum) :: dilaton
  real(kp), dimension(fieldNum,fieldNum) :: metricVal,metricInv  
  real(kp), dimension(2*fieldNum) :: bgVar
  real(kp), dimension(2) :: fieldstop
  real(kp) :: tol
  logical :: paramCheck = .false.
  integer :: ind

!inflation model
  infParam%name = 'radiag'

!parameters (see infbgmodel.f90)
  infParam%consts(1) = 1e-4

  alpha = 2.
  
  infParam%consts(2) = alpha

!fieldstop
  infParam%consts(matterParamNum) = 0

  infParam%conforms = 1
!initial field value.
  infParam%matters(1) = 0

!  print *,'xvmax=', lmi2_x_max_potential(gamma,beta)

!  print *,'xini'


!set the parameters
  paramCheck =  set_infbg_param(infParam)
  print *,'setup params', paramCheck
  print *,'infPAram',infParam
   
!set initial condition  
  infIni = set_infbg_ini(infParam)
  print *,'infIni',infIni

!checkout fieldstop values
  fieldstop = field_stopinf(infParam)
  print *,'fieldstop',fieldstop

!evolves the background till the end of inflation and store the
!results with 5000 points
  if (fieldstop(1).ne.0._kp) then
     infEnd = bg_field_evol(infIni,10000,infObs,ptrToBgdata,fieldstop(1),(fieldstop(2)==1._kp))
  else
     infEnd = bg_field_evol(infIni,10000,infObs,ptrToBgdata)
  endif
     

!physical quantities at the end of inflation
  print *,'infEnd', infEnd, (infEnd==infIni)
  print *,'infObs',infObs
  print *


!test the chain list
  inum =0
  if (associated(ptrToBgdata)) then
     ptrRun => ptrToBgdata
     do while (associated(ptrRun))
!        print *,'efold hubble =',ptrRun%bg%efold, ptrRun%bg%hubble
        inum = inum + 1
        ptrRun => ptrRun%ptr        
     enddo
     ptrRun => null()
     print *,'inum',inum
     print *,'count',count_infbg_data(ptrToBgdata)
  endif


!  stop

!test rescaling
!  call rescale_potential(2._kp,infParam,infIni,infEnd,infObs,ptrToBgdata)
!  print *,'afterRescale'
!  print *,'infParam',infParam
!  print *,'infIni',infIni
!  print *,'infEnd',infEnd
!  print *,'infObs',infObs
!  print *
!  read(*,*)
!  inum=0
   
  epsmax = 0._kp

  if (associated(ptrToBgdata)) then
     ptrRun => ptrToBgdata

     call delete_file('resfield.dat')
     call delete_file('resfielDotd.dat')
     call delete_file('reshubble.dat')
     call delete_file('resepsilons.dat')
     call delete_file('resricci.dat')

!     xend = lmi1_x_endinf(gamma,beta)

     xend = infEnd%field(1)


     do while (associated(ptrRun))
!        print *,'efold hubble =',ptrRun%bg%efold, ptrRun%bg%hubble
        efold = ptrRun%bg%efold
        field = ptrRun%bg%field
        fieldDot = ptrRun%bg%fieldDot
        hubble = ptrRun%bg%hubble
        epsilon1 = ptrRun%bg%epsilon1

!        epsmax = max(epsilon1,epsmax)
!        x = rgi_x_trajectory(efold-infEnd%efold,xend,alpha)
!        epsilon1SR = rgi_epsilon_one(field(1),alpha)
!        epsilon2SR = rgi_epsilon_two(field(1),alpha)

!        epsilon1JF =  ptrRun%bg%epsilon1JF
        ricciOverH2 = 6._kp*(2._kp-epsilon1)
        ricci = ricciOverH2*hubble*hubble
        call livewrite('resfield.dat',efold,field(1),x)
        call livewrite('resslowroll.dat',efold,epsilon1SR, epsilon2SR)
        call livewrite('resfieldDot.dat',efold,fieldDot(1))
        call livewrite('reshubble.dat',efold,hubble)
        call livewrite('resepsilons.dat',efold,epsilon1)
        call livewrite('resricci.dat',efold,riccioverH2,ricci)
        inum = inum + 1
        ptrRun => ptrRun%ptr             
     enddo
     ptrRun => null()
     print *,'inum',inum
     print *,'count',count_infbg_data(ptrToBgdata)
!     print *,'espmax ',epsmax
  endif

  if (associated(ptrToBgdata)) then
     call free_infbg_data(ptrToBgdata)
  endif
 
 

end program infbackmain


 
