program infbackmain
  use infprec, only : kp
  use infsric
  use infbgmodel
  use infbg
  use infinout

  implicit none

   
  real(kp), dimension(fieldNum) :: field, fieldDot
  type(infbgdata), pointer :: ptrToBgdata => null()
  type(infbgdata), pointer :: ptrRun => null()
  
  type(infbgparam) :: infParam
  type(infbgphys) :: infIni, infObs, infEnd
 
  integer :: inum = 0

  real(kp) :: matter,efold,efoldGo,efoldFin,hubble,epsilon1,epsilon1JF
  real(kp) :: ricci, ricciOverH2
  real(kp), dimension(dilatonNum) :: dilaton
  real(kp), dimension(fieldNum,fieldNum) :: metricVal,metricInv  
  real(kp), dimension(2*fieldNum) :: bgVar
  real(kp) :: tol
  logical :: paramCheck = .false.
  integer :: ind

!inflation model
  infParam%name = 'mixlf'

!parameters (see infbgmodel.f90)
  infParam%consts(1) = 1e-5
  infParam%consts(2) = 2
  infParam%consts(3) = 0.
  infParam%consts(4) = 0.
  infParam%consts(5) = 1.
  infParam%consts(6) = 0.5
  infParam%consts(7:11) = 0.
  infParam%consts(12) = 4


  infParam%conforms = 1
!initial field value.
  infParam%matters(1) = 0



!set the parameters
  paramCheck =  set_infbg_param(infParam)
  print *,'setup params', paramCheck
  print *,'infPAram',infParam
   
!set initial condition  
  infIni = set_infbg_ini(infParam)
  print *,'infIni',infIni
 
!evolves the background till the end of inflation and store the
!results with 5000 points
  infEnd = bg_field_evol(infIni,50000,infObs,ptrToBgdata)
  

!physical quantities at the end of inflation
  print *,'infEnd', infEnd, (infEnd==infIni)
  print *,'infObs',infObs
  print *
  read(*,*)

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
  call rescale_potential(2._kp,infParam,infIni,infEnd,infObs,ptrToBgdata)

  print *,'afterRescale'
  print *,'infParam',infParam
  print *,'infIni',infIni
  print *,'infEnd',infEnd
  print *,'infObs',infObs
  print *
  read(*,*)
  inum=0
   
  
  if (associated(ptrToBgdata)) then
     ptrRun => ptrToBgdata

     call delete_file('resfield.dat')
     call delete_file('resfielDotd.dat')
     call delete_file('reshubble.dat')
     call delete_file('resepsilons.dat')
     call delete_file('resricci.dat')

     do while (associated(ptrRun))
!        print *,'efold hubble =',ptrRun%bg%efold, ptrRun%bg%hubble
        efold = ptrRun%bg%efold
        field = ptrRun%bg%field
        fieldDot = ptrRun%bg%fieldDot
        hubble = ptrRun%bg%hubble
        epsilon1 = ptrRun%bg%epsilon1
        epsilon1JF =  ptrRun%bg%epsilon1JF
        ricciOverH2 = 6._kp*(2._kp-epsilon1)
        ricci = ricciOverH2*hubble*hubble
        call livewrite('resfield.dat',efold,field(1))
        call livewrite('resfieldDot.dat',efold,fieldDot(1))
        call livewrite('reshubble.dat',efold,hubble)
        call livewrite('resepsilons.dat',efold,epsilon1,epsilon1JF)
        call livewrite('resricci.dat',efold,riccioverH2,ricci)
        inum = inum + 1
        ptrRun => ptrRun%ptr             
     enddo
     ptrRun => null()
     print *,'inum',inum
     print *,'count',count_infbg_data(ptrToBgdata)
  endif

  if (associated(ptrToBgdata)) then
     call free_infbg_data(ptrToBgdata)
  endif

 
  infIni = set_infbg_ini(infParam)
  print *,'After ini again'
  print *,'infPAram',infParam
  print *,'infIni',infIni

!evolves the background
  infEnd = bg_field_evol(infIni,100,infObs,ptrToBgdata)
  

!print the physical quantities at the end of inflation
  print *,'infEnd', infEnd, (infEnd==infIni)
  print *,'infObs',infObs
  print *
 

end program infbackmain


 
