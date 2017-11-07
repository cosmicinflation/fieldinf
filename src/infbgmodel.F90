module infbgmodel
  use infprec, only : kp, lenshort
  use infmatter, only : potParamNum, matterNum
  use infdilaton, only : confParamNum, dilatonNum
  implicit none

  private
 
  logical, parameter :: display = .false.

  
!totalnumber of background scalar fields

  integer, parameter :: fieldNum = matterNum + dilatonNum


  
  integer, parameter :: matterParamNum = potParamNum + 2
  integer, parameter :: infParamNum = matterParamNum + confParamNum



  type infbgparam
!label identifier
     character(len=lenshort) :: name
!coupling constants [mattercouplings,dilatoncouplings]
     real(kp), dimension(infParamNum) :: consts     
![dilaton degrees of freedom]
     real(kp), dimension(dilatonNum) :: conforms
!matterfield values
     real(kp), dimension(matterNum) :: matters
  end type infbgparam

  

  interface operator (==)     
     module procedure infparam_equal
  end interface


  interface operator (/=)
     module procedure infparam_unequal
  end interface


  public infbgparam
  public operator(==),operator(/=)


  public fieldNum, dilatonNum, matterNum

  public infParamNum, matterParamNum

  public set_infbg_param
   

contains


  function infparam_equal(infparamA, infparamB)
    implicit none
    type(infbgparam), intent(in) :: infparamA, infparamB
    logical :: infparam_equal

    infparam_equal = ((infparamA%name == infparamB%name) &
         .and. all(infparamA%conforms == infParamB%conforms) &
!         .and. all(infparamA%conformsDot == infParamB%conformsDot) &
         .and. all(infparamA%matters == infparamA%matters) &
!         .and. all(infparamA%mattersDot == infparamA%mattersDot) &
         .and. all(infparamA%consts == infparamB%consts))

  end function infparam_equal




  function infparam_unequal(infparamA, infparamB)
    implicit none
    type(infbgparam), intent(in) :: infparamA, infparamB
    logical :: infparam_unequal

    infparam_unequal = ((infparamA%name /= infparamB%name) &
         .or. any(infparamA%conforms /= infParamB%conforms) &
!         .or. any(infparamA%conformsDot /= infParamB%conformsDot) &
         .or. any(infparamA%matters /= infparamA%matters) &
!         .or. any(infparamA%mattersDot /= infparamA%mattersDot) &
         .or. any(infparamA%consts /= infparamB%consts))

  end function infparam_unequal


    



  function set_infbg_param(infParam)
    use infmatter, only : set_potential_param
    implicit none
    logical :: set_infbg_param
    type(infbgparam), intent(in) :: infParam    
    real(kp), dimension(matterParamNum) :: matterParam
  
    logical :: badParams = .true.
    
    real(kp) :: kpbuffer

    matterParam = 0._kp
    set_infbg_param=.false.
    
!the matter potential is parametrized as (see infmatter.f90)
!
! U = [(p1 + p4 ln F) F^p2 + p3]^p5
!   + [p6 + p7 exp(p8 F) + p9 cos(p10 F + p11)] F^p12
!   + p13 F^p14 + p15 F^p16 exp(p17 F^p18) 
!
!where the p are the potential params. In terms of the  "matterParams", they read
!
! p1 = sign(m1) m1^4 
! p2 = m2
! p3 = sign(m3) m3^4
! p4 = sign(m4) m4^4
! p5 = m5
! p6 = sign(m6) m6^4 
! p7 = sign(m7) m7^4 
! p8 = m8
! p9 = sign(m9) m9^4 
!p10 = m10
!p11 = m11
!p12 = m12
!p13 = sign(m13) m13^4
!p14 = m14
!p15 = sign(m15) m15^4
!p16 = m16
!p17 = m17
!p18 = m18


!The matterParams mi are set from the ci params according to the model
!under scrutiny. Only the ci (infparam%consts) are public.
!
! m19=c19
!
! is a field value that bounds the initial field values (ex, the throat size for kklt)
!
! m20=c20
!
!is a field value that stops inflation (ex, hybrid, kklt) instead of
! the condition epsilon1 = 1
!
!
!


!default initialization for
!fieldUv limit
    matterParam(matterParamNum-1) = infParam%consts(matterParamNum-1)
!fieldStop value
    matterParam(matterParamNum) = infParam%consts(matterParamNum)

!name passed to infmatter module    
    infname = infParam%name
    

    case ('largef','lfi')
! large field inflation: LFI

! U = c1^4 F^c2       

       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp))

       
       if (badParams) then
          write(*,*)'consts = ',infParam%consts(1:2)
          stop 'large field: improper params'
       endif


       matterParam(1) = infParam%consts(1)
       matterParam(2) = infParam%consts(2)
       matterParam(3) = 0._kp
       matterParam(4) = 0._kp
       matterParam(5) = 1._kp



    case ('smallf','sfi')
! small field inflation: SFI

! U = c1^4 [1 - (F/c3)^c2]

       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp) &
            .or.(infParam%consts(3).le.0._kp))


       if (badParams) then
          write(*,*)'model name: ',infParam%name
          write(*,*)'consts = ',infParam%consts(1:3)
          stop 'small field: improper params'
       endif

       matterParam(1) = - infParam%consts(1) &
            /(infParam%consts(3)**(infParam%consts(2)/4._kp))
       matterParam(2) = infParam%consts(2)
       matterParam(3) = infParam%consts(1)
       matterParam(4) = 0._kp
       matterParam(5) = 1._kp

            
       if (maxval(infParam%matters).gt.infParam%consts(3)) then
          write(*,*)'set_infbg_param: not small fields initially'
          write(*,*)'model name: ',infParam%name         
          write(*,*)'matterScale = ',infParam%consts(3)
          write(*,*)'matterField = ',maxval(infParam%matters)          
       endif


    case ('branei','bi')
! brane inflation: BI

! U = c1^4 [1 - (F/c3)^-c2]

       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp) &
            .or.(infParam%consts(3).le.0._kp))
       
       if (badParams) then
          write(*,*)'model name: ',infParam%name
          write(*,*)'consts = ',infParam%consts(1:3)
          stop 'brane inflation: improper params'
       endif

       matterParam(1) = - infParam%consts(1) &
            /(infParam%consts(3)**(-infParam%consts(2)/4._kp))
       matterParam(2) = -infParam%consts(2)
       matterParam(3) = infParam%consts(1)
       matterParam(4) = 0._kp
       matterParam(5) = 1._kp

      
       

    case ('hybrid','vhi')
! valley hybrid inflation: VHI

! U = c1^4 [1 + (F/c3)^c2]

!fieldstop value required (checkout infbounds.f90)

       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp) &
            .or.(infParam%consts(3).le.0._kp))            

        if  (badParams) then
          write(*,*)'model name: ',infParam%name
          write(*,*)'consts = ',infParam%consts(1:infParamNum)
          stop 'valley hybrid inflation: improper params'
       endif
          
       matterParam(1) = infParam%consts(1) &
               /(infParam%consts(3)**(infParam%consts(2)/4._kp))
       matterParam(2) = infParam%consts(2)
       matterParam(3) = infParam%consts(1)
       matterParam(4) = 0._kp
       matterParam(5) = 1._kp


    case ('dysusy','dsi')
! dynamical susy inflation: DSI

! U = c1^4 [1 + (F/c3)^(-c2)]
 
!fieldstop value required (checkout infbounds.f90)

       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp) &
            .or.(infParam%consts(3).le.0._kp))            

        if  (badParams) then
          write(*,*)'model name: ',infParam%name
          write(*,*)'consts = ',infParam%consts(1:infParamNum)
          stop 'dynamical susy inflation: improper params'
       endif
          
       matterParam(1) = infParam%consts(1) &
               /(infParam%consts(3)**(-infParam%consts(2)/4._kp))
       matterParam(2) = -infParam%consts(2)
       matterParam(3) = infParam%consts(1)
       matterParam(4) = 0._kp
       matterParam(5) = 1._kp


    case ('grunma','runma1','runma2','runma3','runma4')
! (generalized) running mass inflation: RMI

! U = c1^4 { 1 + c4[1/c2 - ln(F/c3)] F^c2 }

!fieldstop value required (checkout infbounds.f90)

       badParams = ((infParam%consts(3).le.0._kp).or.(infParam%consts(3).gt.1._kp) &
            .or.(infParam%consts(1).le.0._kp))

       badParams = badParams.or. ( &
            (infParam%name.ne.'grunma').and.(infParam%consts(2).ne.2._kp) )
       badParams = badParams .or. ( &
            (infParam%name.eq.'runma1').and.(infParam%consts(4).le.0._kp) )
       badParams = badParams .or. ( &
            (infParam%name.eq.'runma2').and.(infParam%consts(4).le.0._kp) )
       badParams = badParams .or. ( &
            (infParam%name.eq.'runma3').and.(infParam%consts(4).ge.0._kp) )
       badParams = badParams .or. ( &
            (infParam%name.eq.'runma4').and.(infParam%consts(4).ge.0._kp) )

       
       if  (badParams) then
          write(*,*)'model name: ',infParam%name
          write(*,*)'consts = ',infParam%consts(1:4)
          stop '(generalized) running mass: improper params'
       endif


       kpbuffer = infParam%consts(4)*(1._kp/infParam%consts(2) &
            + log(infParam%consts(3)))
       matterParam(1) = infParam%consts(1) * sign(abs(kpbuffer)**0.25_kp,kpbuffer)
       matterParam(2) = infParam%consts(2)
       matterParam(3) = infParam%consts(1)

       matterParam(4) = -infParam%consts(1) &
            * sign(abs(infParam%consts(4))**0.25_kp,infParam%consts(4))
     
       matterParam(5) = 1._kp
       

    case('kklmmt','kklti')
! KKLT inflation: KKTLI
! U = c1^4 / [1 + (F/c3)^(-c2)] with c2 > 0 for c5=-1
!
! + Phi_string/mu related to the flux number N
!
! + PhiUv related to the brane tension and the string coupling

!the case c5=+1 is branei

       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp) &
            .or.(infParam%consts(3).le.0._kp) &            
            .or.(infParam%consts(matterParamNum-1).lt.0._kp))

       if (badParams) then
          write(*,*)'model name: ',infParam%name
          write(*,*)'consts = ',infParam%consts(1:4)
          stop 'kklmmt: improper params'
       endif


       matterParam(1) = 1._kp/infParam%consts(1) &
            /(infParam%consts(3)**(-infParam%consts(2)/4._kp))
       matterParam(2) = -infParam%consts(2)
       matterParam(3) = 1._kp/infParam%consts(1)
       matterParam(4) = 0._kp
       matterParam(5) = -1._kp

!those are necessary model parameters: checkout infbounds.f90
!fieldUv: prevents brane out of the throat
!flux number N
       
            

    case ('rcquad','rcqi')
! radiatively corrected quartic inflation: RCQI

! U = c1^4 F^4 [1 - c2 ln(F)]
          
       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp))
       
       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:2)
          stop 'radiative corrected quartic field: improper params'
       endif

       matterParam(1) = infParam%consts(1)
       matterParam(2) = 4._kp
       matterParam(3) = 0._kp
       matterParam(4) = -infParam%consts(1) &
            *sign(abs(infParam%consts(2))**0.25_kp,infParam%consts(2))
       matterParam(5) = 1._kp


    case ('gswli','li')
! global susy loop inflation: LI

! U = c1^4 [1 + c2 ln(F)]
          
       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp))
       
       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:2)
          stop 'global susy loop inflation: improper params'
       endif

       matterParam(1) = 0._kp
       matterParam(2) = 0._kp
       matterParam(3) = infParam%consts(1)
       matterParam(4) = infParam%consts(1) * infParam%consts(2)**0.25_kp
       matterParam(5) = 1._kp                


    case ('colwei','cwi')
!Coleman-Weinberg inflation: CWI

! U = c1^4 [1 + c2 ln(F/c3) (F/c3)^c4]
       
       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp))
       badParams = badParams.or.(infParam%consts(4).ne.4._kp)
       badParams = badParams.or.(infParam%consts(2).ne.(4._kp*exp(1._kp)))

       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:4)
          stop 'Coleman-Weinberg inflation: improper params'
       endif

       matterParam(1) = -infParam%consts(1) &
            /(infParam%consts(3)**(0.25_kp*infParam%consts(4))) &
            * sign(abs(infParam%consts(2)*log(infParam%consts(3)))**0.25_kp &
            , infParam%consts(2)*log(infParam%consts(3)))
       matterParam(2) = infParam%consts(4)
       matterParam(3) = infParam%consts(1)
       matterParam(4) = infParam%consts(1) &
            * (infParam%consts(2)/infParam%consts(3)**infParam%consts(4))**0.25_kp
       matterParam(5) = 1._kp

    case ('tdwell','dwi')
!topological double well inflation: DWI

! U = c1^4 [(F/c2)^2 - 1]^2

       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(2).lt.sqrt(8._kp)))
       
       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:2)
          stop 'topological double well: improper params'
       endif

       matterParam(1) = sqrt(infParam%consts(1)/infParam%consts(2))
       matterParam(2) = 2._kp
       matterParam(3) = -sqrt(infParam%consts(1))
       matterParam(4) = 0._kp
       matterParam(5) = 2._kp


    case ('betexp','bei')
! beta exponential inflation: BEI

! U =c1^4 (1 - c2 F)^c3

       badParams = (infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp)
       badParams = badParams.or.(infParam%consts(3).le.0._kp)
       
       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:3)
          stop 'beta-exponential: improper params'
       endif


       matterParam(5) = infParam%consts(3)
       matterParam(4) = 0._kp
       matterParam(3) = infParam%consts(1)**(1._kp/infParam%consts(3))
       matterParam(2) = 1._kp
       matterParam(1) = - infParam%consts(2)**0.25_kp * matterParam(3)


    case ('radiag','rgi')
!radion assisted gauge inflation: RGI

!U = c1^4 / [1 + c2 F^-2]

       badParams = (infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp)
       
       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:2)
          stop 'radion assisted gauge: improper params'
       endif

       matterParam(5) = -1._kp
       matterParam(4) = 0._kp
       matterParam(3) = 1._kp/infParam%consts(1)
       matterParam(2) = -2._kp
       matterParam(1) = infParam%consts(2)**0.25_kp * matterParam(3)


    case ('nszero','csi')
!constant spectrum inflation: CSI

!U = c1^4 / [1 - c2 F]^2

       badParams = (infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp)
       
       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:2)
          stop 'constant spectrum inflation: improper params'
       endif

       matterParam(5) = -2._kp
       matterParam(4) = 0._kp
       matterParam(3) = 1._kp/sqrt(infParam%consts(1))
       matterParam(2) = 1._kp
       matterParam(1) = - infParam%consts(2)**0.25_kp * matterParam(3)


    case ('sugrab','sbi')
!supergravity brane inflation: SBI

!U = c1^4 [ 1 + F^4 (c3 ln F - c2) ]

       badParams = (infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp)
       badParams = badParams.or.(infParam%consts(3).le.0._kp)
       
       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:3)
          stop 'supergravity brane inflation: improper params'
       endif
       
       matterParam(5) = 1._kp
       matterParam(4) = infParam%consts(1) * infParam%consts(3)**0.25_kp
       matterParam(3) = infParam%consts(1)
       matterParam(2) = 4._kp
       matterParam(1) = - infParam%consts(1) * infParam%consts(2)**0.25_kp

    case ('logpot','logpo1','logpo2','logpo3','witorh')
!logarithimc potential inflation: LPI
!Witten O.Raifeartaigh inflation: WRI

!U = c1^4 (F/c3)^c2 [ln(F/c3)]^c4

       badParams = (infParam%consts(3).le.0._kp).or.(infParam%consts(1).le.0._kp)

       badParams = badParams.or.( &
            (infParam%name.eq.'logpo2').and.(modulo(infParam%consts(4),2._kp).ne.0._kp))
       badParams = badParams.or.( &
            (infParam%name.eq.'logpo3').and.(modulo(infParam%consts(4),2._kp).ne.0._kp))

       badParams = badParams.or.((infParam%name.eq.'witorh').and. &
            ((infParam%consts(2).ne.0._kp).or.(infParam%consts(4).ne.2._kp)))

       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:4)
          stop 'logpot / Witten O.Raifeartaigh inflation: improper params'
       endif
       
       matterParam(5) = infParam%consts(4)
       matterParam(4) = (infParam%consts(1) &
            /infParam%consts(3)**(infParam%consts(2)*0.25_kp))**(1._kp/infParam%consts(4))
       matterParam(3) = 0._kp
       matterParam(2) = infParam%consts(2)/infParam%consts(4)
       matterParam(1) = -sign(abs(log(infParam%consts(3)))**0.25_kp,log(infParam%consts(3))) &
            * (infParam%consts(1)/infParam%consts(3)**(infParam%consts(2)*0.25_kp)) &
            **(1._kp/infParam%consts(4))


    case ('ostach','osti')
!open string tachyonic inflation: OSTI

!U = -c1^4 (F/c3)^2 ln[(F/c3)^2]

       badParams = (infParam%consts(3).le.0._kp).or.(infParam%consts(1).le.0._kp)
       badParams = badParams.or.(infParam%consts(2).ne.2._kp)
       
       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:3)
          stop 'open string tachyonic inflation: improper params'
       endif
       
       matterParam(5) = 1._kp
       matterParam(4) = -infParam%consts(2)**0.25_kp * infParam%consts(1) &
            /infParam%consts(3)**(infParam%consts(2)*0.25_kp)
       matterParam(3) = 0._kp
       matterParam(2) = infParam%consts(2)
       matterParam(1) = infParam%consts(2)**0.25_kp &
            * sign(abs(log(infParam%consts(3)))**0.25_kp,log(infParam%consts(3))) &
            * (infParam%consts(1)/infParam%consts(3)**(infParam%consts(2)*0.25_kp))


     case ('invmon','imi')
!inverse monomial inflation: IMI

! U = c1^4 F^(-c2)       

       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp))
       
       if (badParams) then
          write(*,*)'consts = ',infParam%consts(1:2)
          stop 'inverse monomial inflation: improper params'
       endif

       matterParam(1) = infParam%consts(1)
       matterParam(2) = -infParam%consts(2)
       matterParam(3) = 0._kp
       matterParam(4) = 0._kp
       matterParam(5) = 1._kp

#ifndef PP5


    case ('gmixlf','gmlfi')
!generalized mixed large field inflation: GMLFI

! U = c1^4 F^c2 [1 + c3 F^c4]

       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp) &
            .or.(infParam%consts(3).le.0._kp).or.(infParam%consts(4).le.0._kp))
       
       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:infParamNum)
          stop 'generalized mixed large field inflation: improper params'
       endif


       matterParam(1) = infParam%consts(1)
       matterParam(2) = infParam%consts(2)
       matterParam(3) = 0._kp
       matterParam(4) = 0._kp
       matterParam(5) = 1._kp
       matterParam(6) = infParam%consts(1) &
            *sign(abs(infParam%consts(3))**0.25_kp,infParam%consts(3))
       matterParam(7:11) = 0._kp
       matterParam(12) = infParam%consts(2) + infParam%consts(4)



    case ('lfcorr','rclfi')
!radiatively corrected large field inflation: RCLFI

!U = c1^4 F^c2 [1 - c3 F^c4 ln(F)]       


       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp).or. &
            (infParam%consts(3).le.0._kp).or.(infParam%consts(4).le.0._kp))
       
       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:4)
          stop 'radiatively corrected large field inflation: improper params'
       endif

       matterParam(1) = 0._kp
       matterParam(2) = infParam%consts(2) + infParam%consts(4)
       matterParam(3) = 0._kp
       matterParam(4) = - infParam%consts(1) &
            * sign(abs(infParam%consts(3))**0.25_kp,infParam%consts(3))
       matterParam(5) = 1._kp
       matterParam(6) = infParam%consts(1)
       matterParam(7:11) = 0._kp
       matterParam(12) = infParam%consts(2)


    case ('rcmass','rcmi')
!radiatively corrected massive inflation: RCMI

!U = c1^4 F^2 [1 - c2 F^2 ln F]

!c2  is 2 x alpha

       badParams = ((infParam%consts(1).le.0._kp).or.infParam%consts(2).le.0._kp)

       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:2)
          stop 'radiatively corrected massive inflation: improper params'
       endif
       
       matterParam(1) = 0._kp
       matterParam(2) = 4._kp
       matterParam(3) = 0._kp
       matterParam(4) = - infParam%consts(1) &
            * sign(abs(infParam%consts(2))**0.25_kp,infParam%consts(2))
       matterParam(5) = 1._kp
       matterParam(6) = infParam%consts(1)
       matterParam(7:11) = 0._kp
       matterParam(12) = 2._kp


    case ('natinf','ni')
!natural inflation: NI

!U = c1^4 [1 + cos(F/c2)]

       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp))
       
       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:2)
          stop 'Natural inflation: improper params'
       endif
       
       matterParam(1:4) = 0._kp
       matterParam(5) = 2._kp
       matterParam(6) = infparam%consts(1)
       matterParam(7) = 0._kp
       matterParam(8) = 0._kp
       matterParam(9) = infparam%consts(1)
       matterParam(10) = 1._kp/infParam%consts(2)
       matterParam(11) = 0._kp
       matterParam(12) = 0._kp


    case ('hybnat','hni')
!hybrid natural inflation: HNI

!U = c1^4 [1 + c3 cos(F/c2)]

       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp) &
            .or. (infParam%consts(3).gt.1._kp) &
            .or. (infParam%consts(3).le.0._kp))
       
       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:3)
          stop 'hybrid natural inflation: improper params'
       endif
       
       matterParam(1:4) = 0._kp
       matterParam(5) = 2._kp
       matterParam(6) = infparam%consts(1)
       matterParam(7) = 0._kp
       matterParam(8) = 0._kp
       matterParam(9) = infparam%consts(1)*infParam%consts(3)**0.25_kp
       matterParam(10) = 1._kp/infParam%consts(2)
       matterParam(11) = 0._kp
       matterParam(12) = 0._kp
       

    case ('nrcoli','ncli')
!non-renormalisable corrected loop inflation: NCLI

!U = c1^4 [1 + c2 ln(F) + (F/c3)^(4+2c4)]

       badParams = (infParam%consts(1).le.0._kp) &
            .or. (infParam%consts(2).le.0._kp) &
            .or. (infParam%consts(2).gt.1._kp)

       badParams = badParams &
            .or. (infParam%consts(3).le.0._kp) &
            .or. (infParam%consts(4).le.0._kp)
       
       if (badParams) then
          write(*,*)'model name: ',infParam%name
          write(*,*)'consts= ',infParam%consts
          stop 'non-renormalizable corrected loop inflation: improper params'
       endif

       matterParam(1) = infParam%consts(1)
       matterParam(2:3) = 0._kp
       matterParam(4) = infParam%consts(1)*infParam%consts(2)**0.25_kp
       matterParam(5) = 1._kp
       matterParam(6) = infparam%consts(1) &
            / infParam%consts(3)**(1._kp + 0.5_kp*infParam%consts(4))
       matterParam(7:11) = 0._kp
       matterParam(12) = 4._kp + 2._kp*infParam%consts(4)
       


    case ('exsusy','esi')
!exponential susy inflation: ESI

!U = c1^4 [1 - exp(-c2 F)]


       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp))
       
       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:2)
          stop 'exponential susy inflation: improper params'
       endif

       matterParam(1:4) = 0._kp
       matterParam(5) = 2._kp
       matterParam(6) = infParam%consts(1)
       matterParam(7) = -infParam%consts(1)
       matterParam(8) = -infParam%consts(2)
       matterParam(9:12) = 0._kp


    case ('powlaw','pli')
!power law inflation: PLI

!U = c1^4 exp[-c2 F]

!fieldstop value required (checkout infbounds.f90)
       
       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp))
       
       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:2)
          stop 'power law inflation: improper params'
       endif

       matterParam(1:4) = 0._kp
       matterParam(5) = 2._kp
       matterParam(6) = 0._kp
       matterParam(7) = infParam%consts(1)
       matterParam(8) = -infParam%consts(2)
       matterParam(9:12) = 0._kp


     

    case ('hfline','hf1i')
!U = c1^4 [ (1 + c2 F)^2 - 2/3 c2^2 ]

!c2 is A1
       badParams = (infParam%consts(1).le.0._kp)

      
       matterParam(1) = sqrt(infParam%consts(1)) &
            *  sign(abs(infParam%consts(2))**0.25_kp,infParam%consts(2))
       matterParam(2) = 1._kp
       matterParam(3) = sqrt(infParam%consts(1))
       matterParam(4) = 0._kp
       matterParam(5) = 2._kp
       matterParam(6) = - infParam%consts(1) &
            * (2._kp/3._kp*infParam%consts(2)**2)**0.25_kp
       matterParam(7:12) = 0._kp


    case ('interm','ii')
!intermediate inflation: II

! U = c1^4 F^(-c2) [1 - c2^(2/6) F^(-2)]

       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp))           
       
       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:infParamNum)
          stop 'intermediate inflation: improper params'
       endif


       matterParam(1) = infParam%consts(1)
       matterParam(2) = -infParam%consts(2)
       matterParam(3) = 0._kp
       matterParam(4) = 0._kp
       matterParam(5) = 1._kp
       matterParam(6) = -infParam%consts(1) * (infParam%consts(2)**2/6._kp)**0.25_kp
       matterParam(7:11) = 0._kp
       matterParam(12) = -infParam%consts(2) - 2._kp


    case ('twisti','twi')
!twisted inflation: TWI

!U = c1^4 [1 - c2 (F/c3)^2 exp(-F/c3)]

!c2 = 32/[92 ζ(5)]~0.33183220
!fieldstop value required (checkout infbounds.f90)

       badParams = (infParam%consts(1).le.0._kp).or.(infParam%consts(3).le.0._kp)
!       badParams = badParams.or.(infParam%consts(3).gt.0.04228)
       badParams = badParams.or.(abs(infParam%consts(2)-0.33183220).gt.epsilon(1.))       
       
       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:3)
          stop 'twisted inflation: improper params'
       endif

       matterParam(1:12) = 0._kp
       matterParam(3) = infParam%consts(1)
       matterParam(5) = 1._kp
       matterParam(7) = - infParam%consts(1) * infParam%consts(2)**0.25_kp &
            /sqrt(infParam%consts(3))
       matterParam(8) = -1._kp/infParam%consts(3)
       matterParam(12) = 2._kp


    case ('nckahi','ncki')
!non-canonical Kahler inflation: NCKI

!U = c1^4 [ 1 + c2 ln F + c3 F^2 ]

       badParams = ((infParam%consts(1).le.0._kp) &
            .or.(infParam%consts(2).le.0._kp))
       
       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:3)
          stop 'non-canonical Kahler inflation: improper params'
       endif
          
       matterParam(1:2) = 0._kp
       matterParam(3) = infParam%consts(1)
       matterParam(4) = infParam%consts(1)*infParam%consts(2)**0.25_kp
       matterParam(5) = 1._kp

       matterParam(6) = infParam%consts(1)*sign(abs(infParam%consts(3))**0.25_kp,infParam%consts(3))
       matterParam(7:11) = 0._kp
       matterParam(12) = 2._kp


    case ('oifold','oi')
!orientifold inflation: OI

!U = c1^4 (F/c3)^4 [ c2 ln^2(F/c3) - 1 ]

       badParams = ((infParam%consts(1).le.0._kp).or.infParam%consts(2).le.0._kp)
       badParams = badParams.or.(infParam%consts(3).le.0._kp)
       
       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:3)
          stop 'orientifold inflation: improper params'
       endif

       matterParam(12) = 4._kp
       matterParam(7:11) = 0._kp
       matterParam(6) = -infParam%consts(1)*infParam%consts(2)**0.25_kp &
            /infParam%consts(3)
       matterParam(5) = 2._kp
       matterParam(4) = sqrt(infParam%consts(1)/infParam%consts(3))
       matterParam(3) = 0._kp
       matterParam(2) = 2._kp
       matterParam(1) = -sqrt(infParam%consts(1)/infParam%consts(3)) &
            * sign(abs(log(infParam%consts(3)))**0.25_kp,log(infParam%consts(3)))

    case ('ssbinf','spsyb1','spsyb2','spsyb3','spsyb4','spsyb5','spsyb6')
!spontaneous symmetry breaking inflation: SSBI

!U = c1^4 [1 + alpha F^2 + beta F^4 ]

       badParams = (infParam%consts(1).le.0._kp)

       badParams = badParams.or. ( (infparam%name.eq.'spsyb1').and.(.not.( &
            (infParam%consts(2).gt.0._kp).and.(infParam%consts(3).gt.0._kp))))
       badParams = badParams.or. ( (infparam%name.eq.'spsyb2').and.(.not.( &
            (infParam%consts(2).lt.0._kp).and.(infParam%consts(3).lt.0._kp))))
       badParams = badParams.or. ( (infparam%name.eq.'spsyb3').and.(.not.( &
            (infParam%consts(2).gt.0._kp).and.(infParam%consts(3).lt.0._kp))))
       badParams = badParams.or. ( (infparam%name.eq.'spsyb4').and.(.not.( &
            (infParam%consts(2).gt.0._kp).and.(infParam%consts(3).lt.0._kp))))
       badParams = badParams.or. ( (infparam%name.eq.'spsyb5').and.(.not.( &
            (infParam%consts(2).lt.0._kp).and.(infParam%consts(3).gt.0._kp))))
       badParams = badParams.or. ( (infparam%name.eq.'spsyb6').and.(.not.( &
            (infParam%consts(2).lt.0._kp).and.(infParam%consts(3).gt.0._kp))))
       

       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:3)
          stop 'spontaneous symmetry breaking inflation: improper params'
       endif


       matterParam(12) = 4._kp
       matterParam(7:11) = 0._kp
       matterParam(6) = infParam%consts(1) &
            * sign(abs(infParam%consts(3))**0.25_kp,infParam%consts(3))
       matterParam(5) = 1._kp
       matterParam(4) = 0._kp
       matterParam(3) = infParam%consts(1)
       matterParam(2) = 2._kp
       matterParam(1) = infParam%consts(1) &
            * sign(abs(infParam%consts(2))**0.25_kp,infParam%consts(2))


#ifndef PP12

    case ('kahmod','kmi')
!Kahler moduli inflation: KMI

!U = c1^4 [1 - c2 F^c3 exp(-c4 F^c5)]

       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp) &
            .or. (infParam%consts(3).ne.infParam%consts(5)) )

       badParams = ( badParams .or. &
            .not. ((infParam%consts(3).eq.1._kp) &
            .or. (infParam%consts(3).eq.4._kp/3._kp)) )
            
       badParams = ( badParams .or. &
            ((infParam%consts(3).eq.1._kp).and.(infParam%consts(4).ne.1._kp)) )
       
        if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:5)
          stop 'kahler moduli: improper params'
       endif

       matterParam(1:4) = 0._kp      
       matterParam(5) = 2._kp
       matterParam(6:12) = 0._kp

       matterParam(13) = infParam%consts(1)
       matterParam(14) = 0._kp
       matterParam(15) = -infParam%consts(1) &
            * sign(abs(infParam%consts(2))**0.25_kp,infParam%consts(2))
       matterParam(16) = infParam%consts(3)
       matterParam(17) = -infParam%consts(4)
       matterParam(18) = infParam%consts(5)


    case ('higgsi','hi','si')
!Higgs/Starobinski inflation: HI/SI

!U = c1^4 [1 - exp(-c2 F)]^2 with c2 = sqrt(2/3)


       badParams = ((infParam%consts(1).le.0._kp) &
            .or.(infParam%consts(2).ne.sqrt(2._kp/3._kp)))
       
       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:2)
          stop 'Higgs/Starobinski inflation: improper params'
       endif

       matterParam(3) = infParam%consts(1)
       matterParam(5) = 1._kp
       matterParam(7) = - 2**0.25_kp * infParam%consts(1)
       matterParam(8) = - infParam%consts(2)
       matterParam(15) = infParam%consts(1)
       matterParam(17) = -2._kp*infParam%consts(2)
       matterParam(18) = 1._kp

    case ('logmdi','logmd1','logmd2','lmi')
!logamediate inflation: LMI

!U = c1^4 F^c2 exp(-c3 F^c4) with c2 = 4*(1-c4)

       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(3).le.0._kp) &
            .or. (infParam%consts(4).lt.0._kp) .or. (infParam%consts(4).gt.1._kp ))

       badParams = ( badParams .or. &
            .not. ( (infParam%consts(2).eq.4._kp*(1._kp-infParam%consts(4)) ) ) )                  
       
       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:4)
          stop 'logamediate inflation: improper params'
       endif

       matterParam(1:4) = 0._kp      
       matterParam(5) = 2._kp
       matterParam(6:14) = 0._kp
      
       matterParam(15) = infParam%consts(1)
       matterParam(16) = infParam%consts(2)
       matterParam(17) = -infParam%consts(3)
       matterParam(18) = infParam%consts(4)


    case ('nmssmi','gmssmi','orifpt','nrifpt','grifpt')
!(generalised) minimal supersymmetry inflation: (G)MSSMI
!(generelised) renormalisable inflection point inflation: G(RIPI)

!U = c1^4 [ (F/c6)^2 -  c3 (F/c6)^c2 + c4 (F/c6)^c5] 

       badParams = (infParam%consts(1).le.0._kp) &
            .or.(infParam%consts(3).lt.0._kp) &
            .or.(infParam%consts(4).lt.0._kp) &
            .or.(infParam%consts(6).le.0._kp)

!c2=6; c5=10 and c3=2/3  and c4=1/5 c6=mu
       if (infParam%name.eq.'nmssmi') then
          badParams = badParams.or.(infParam%consts(2).ne.6._kp) &
               .or.(infParam%consts(5).ne.10._kp) &
               .or.(infParam%consts(3).ne.2._kp/3._kp) &
               .or.(infParam%consts(4).ne.1._kp/5._kp)

!c2=6; c5=10; c3=-2/3 alpha and c4=alpha/5
       elseif (infParam%name.eq.'gmssmi') then
          badParams = badParams.or.(infParam%consts(2).ne.6._kp) &
               .or.(infParam%consts(5).ne.10._kp) &
               .or.(1.5_kp*infParam%consts(3).ne.5._kp*infParam%consts(4))

!c2=3; c5=4; c3=4/3 and c4=1/2 (grifpt with alpha=1)
       elseif (infParam%name.eq.'nrifpt') then
          badParams = badParams.or.(infParam%consts(2).ne.3._kp) &
               .or. (infParam%consts(5).ne.4._kp) &
               .or. (3._kp/4._kp * infParam%consts(3).ne.1._kp) &
               .or. (2._kp*infParam%consts(4).ne.1._kp)

!c2=3; c5=4; c3=4/3 alpha and c4=alpha/2
       elseif (infParam%name.eq.'grifpt') then
          badParams = badParams.or.(infParam%consts(2).ne.3._kp) &
               .or. (infParam%consts(5).ne.4._kp) &
               .or. (3._kp/4._kp * infParam%consts(3).ne.2._kp*infParam%consts(4))

!original parametrization, for test only, use nrifpt instead which is
!grifpt with alpha=1
! c2=3; c5=4 and c4=9/32 c3^2 c6=1
       elseif (infParam%name.eq.'orifpt') then
          badParams = badParams.or.(infParam%consts(2).ne.3._kp) &
               .or.(infParam%consts(5).ne.4._kp) &
               .or.(9._kp/32._kp*infParam%consts(3)**2.ne.infParam%consts(4)) &
               .or.(infParam%consts(6).ne.1._kp)

       endif
       

       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:5)
          stop 'MSSM / RINFPT inflation: improper params'
       endif

       matterParam(1) = infParam%consts(1)/sqrt(infParam%consts(6))
       matterParam(2) = 2._kp
       matterParam(3:4) = 0._kp
       matterParam(5) = 1._kp
       
       matterParam(6) = -infParam%consts(1) &
            *sign(abs(infparam%consts(3)/infParam%consts(6)**infParam%consts(2))**0.25_kp &
            , infparam%consts(3))
       matterParam(7:10) = 0._kp
       matterParam(11) = infParam%consts(6)
       matterParam(12) = infparam%consts(2)

       matterParam(13) = infparam%consts(1) &
            * sign(abs(infparam%consts(4)/infParam%consts(6)**infParam%consts(5))**0.25_kp &
            ,infParam%consts(4))
       matterParam(14) = infParam%consts(5)
       matterParam(15:18) = 0._kp


    case ('bsusyb','bsusybi')
!brane susy breaking inflation: BSUSYBI

!U = c1^4 [ exp(c2 F) + c3 exp(c4 F) ]

!fieldstop value required (checkout infbounds.f90)

       badParams = (infParam%consts(1).lt.0._kp &
            .or.(infParam%consts(2).ne.sqrt(6._kp)) &
            .or.(infParam%consts(3).ne.1._kp))

       
       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:4)
          stop 'Brane SUSY breaking inflation: improper params'
       endif
       
       matterParam(1:4) = 0._kp
       matterParam(2) = 2._kp
       matterParam(5) = 2._kp
       
       matterParam(6) = 0._kp
       matterParam(7) = infParam%consts(1)
       matterParam(8) = infParam%consts(2)
       matterParam(9:11) = 0._kp
       matterParam(12) = 0._kp

       matterParam(13:14) = 0._kp
       matterParam(15) = infParam%consts(1)*infParam%consts(3)**0.25_kp
       matterParam(16) = 0._kp
       matterParam(17) = infParam%consts(4)
       matterParam(18) = 1._kp
       

    case ('nformi','nform1','nform2','nform3','nform4')
!N-formalism inflation: NFI

!U = c1^4 exp[-c2 F^c3]

!fieldstop required for nform2 and nform4

       badParams = (infParam%consts(1).le.0._kp)

       badParams = badParams.or. ( (infparam%name.eq.'nform1').and.(.not.( &
            (infParam%consts(2).gt.0._kp).and.(infParam%consts(3).gt.1._kp))))
       badParams = badParams.or. ( (infparam%name.eq.'nform2').and.(.not.( &
            (infParam%consts(2).lt.0._kp).and.(infParam%consts(3).gt.1._kp))))
       badParams = badParams.or. ( (infparam%name.eq.'nform3').and.(.not.(( &
            (infParam%consts(2).lt.0._kp).and.(infParam%consts(3).lt.1._kp).and.&
            (infParam%consts(3).gt.0._kp)) .or. (&
            (infParam%consts(2).gt.0._kp).and.(infParam%consts(3).lt.0._kp)))))
       badParams = badParams.or. ( (infparam%name.eq.'nform4').and.(.not.(( &
            (infParam%consts(2).gt.0._kp).and.(infParam%consts(3).lt.1._kp).and.&
            (infParam%consts(3).gt.0._kp) ) .or. (&
            (infParam%consts(2).lt.0._kp).and.(infParam%consts(3).lt.0._kp)))))
       
       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:3)
          stop 'N-formalism inflation: improper params!'
       endif

       matterParam(1:14) = 0._kp
       matterParam(5) = 2._kp
       matterParam(15)= infParam%consts(1)
       matterParam(16) = 0._kp
       matterParam(17) = -infParam%consts(2)
       matterParam(18) = infParam%consts(3)


    case ('scaaai','saai')
!superconformal alpha attractor A inflation: SAAI

!U = c1^4 {1 - exp[-sqrt(2/3/c2) F]}^2

       badParams = (infParam%consts(1).le.0._kp) &
            .or. (infParam%consts(2).le.0._kp)
       
       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:2)
          stop 'superconformal alpha attractor A inflation: improper params'
       endif
       
       matterParam(1:16) = 0._kp
       matterParam(2) = 2._kp
       matterParam(5) = 2._kp
       
       matterParam(6) = infParam%consts(1)
       matterParam(7) = -2._kp**(0.25_kp)*infParam%consts(1)
       matterParam(8) = -sqrt(2._kp/3._kp/infparam%consts(2))
       matterParam(15) = infParam%consts(1)
       matterParam(17) = -2._kp*sqrt(2._kp/3._kp/infParam%consts(2))
       matterParam(18) = 1._kp


#endif
#endif

!potentials not encompassed in the generic formula
#ifdef PPNAME
    case ('mhitop','mhi')
!Mutated hilltop inflation: MHI

!U = c1^4 [1 - 1/cosh(F/c2)]

       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp))
       
       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:2)
          stop 'mutated hilltop inflation: improper params'
       endif

       matterParam(1) = infParam%consts(1)
       matterParam(2) = infParam%consts(2)


    case ('ricci1', 'ricci2')
!R plus R^p inflation: RPI

!F -> sqrt(2/3) F
!U = c1^4 e^(-2 F) [e^F - 1]^[c2/(c2-1/2)]

       badParams = ((infParam%consts(1).le.0._kp) &
            .or.(infParam%consts(2).lt.1._kp))
       
       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:2)
          stop 'R + R^p inflation: improper params'
       endif

       matterParam(1) = infParam%consts(1)
       matterParam(2) = infParam%consts(2)


    case ('ccorsi','corsi1','corsi2','corsi3','ccsi')
!cubiquely corrected Starobinksi inflation: CCSI

!F -> sqrt(2/3) F
!U = c1^4 e^(-2 F) [e^F - 1]^2
!       x {1 + 2 Sqrt[1 + 3 c2 (e^F - 1)] + 2 c2 (e^F -1 )}
!       / {1 + Sqrt[1 + 3 c2 (e^F - 1)]}^3

       badParams = ((infParam%consts(1).le.0._kp) &
            .or.(abs(infParam%consts(2)).gt.1._kp))

       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:2)
          stop 'cubiquely corrected Starobinksi inflation: improper params'
       endif


       matterParam(1) = infParam%consts(1)
       matterParam(2) = infParam%consts(2)

!fieldstop required for corsi2

    case ('tipinf','ti')
!tip inflation: TI

!U = c1^4 [ 1 + cos(F/c3) + c2 sin^2(F/c3) ]

       badParams = ((infParam%consts(1).le.0._kp) &
            .or. (infParam%consts(2).le.0._kp) &
            .or. (infParam%consts(3).le.0._kp))

       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:3)
          stop 'tip inflation: improper params'
       endif

       matterParam(1) = infParam%consts(1)
       matterParam(2) = infParam%consts(2)
       matterParam(3) = infParam%consts(3)

    case ('psenat','psni')
!pseudo-natural inflation: PSNI

!U = c1^4 { 1 + c2 ln[cos(F/c3)] }

       badParams = ((infParam%consts(1).le.0._kp) &
            .or. (infParam%consts(2).le.0._kp) &
            .or. (infParam%consts(3).le.0._kp))

       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:3)
          stop 'pseudo natural inflation: improper params'
       endif

       matterParam(1) = infParam%consts(1)
       matterParam(2) = infParam%consts(2)
       matterParam(3) = infParam%consts(3)


    case ('arctan','ai')
!arctangent inflation: AI

!U = c1^4 [1 - 2/pi arctan(F/c2) ]
       
       badParams = ((infParam%consts(1).le.0._kp) &
            .or. (infParam%consts(2).le.0._kp))

       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:2)
          stop 'arctan inflation: improper params'
       endif

       matterParam(1) = infParam%consts(1)
       matterParam(2) = infParam%consts(2)      

    case ('fixnsa','fixnsb','fixnsc','cnai','cnbi','cnci')
!fix spectral index inflation: CN(ABC)I

!U = c1^4 {3 - (3 + 3 c2^2) tanh^2[c2 F/sqrt(2)] }
!U = c1^4 {-3 + (3 - 3 c2^2) tan^2[c2 F/sqrt(2)] }

!U = c1^4 {-3 + (3 + alpha^2) /tanh^2[c2 F/sqrt(2)] }
!fieldstop value required (checkout infbounds.f90)


       badParams = ((infParam%consts(1).le.0._kp) &
            .or. (infParam%consts(2).le.0._kp))

       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:2)
          stop 'Constant ns inflation: improper params'
       endif

       matterParam(1) = infParam%consts(1)
       matterParam(2) = infParam%consts(2)

    
    case ('fixnsd','cndi')
!fix spectral index inflation: CNDI

!U = c1^4 / [ 1 + c3 cos(c2 F) ]^2
!fieldstop value required (checkout infbounds.f90)

       badParams =  ((infParam%consts(1).le.0._kp) &
            .or. (infParam%consts(2).le.0._kp) &
            .or. (infParam%consts(3).le.0._kp))

       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:3)
          stop 'Constant ns inflation: improper params'
       endif

       matterParam(1) = infParam%consts(1)
       matterParam(2) = infParam%consts(2)
       matterParam(3) = infParam%consts(3)


    case ('dualsb','di')
!dual inflation from foftly broken N=2 super Yang-Mills theories: DI

!U(k2) = M^4 { 1 + Vo(c1) - 2(K-E)/(k2 K) - pi/(k2 K K') [nu(k2)]^2 Heaviside[nu(k2)]}
!nu(k2) = 1 - 8 sqrt(2)/(pi^2 c1) K/sqrt(k2)(E'-K')^2
!M^4 = c1^2 c2^4 / pi^2
!dF/dk2 = [4 sqrt(2)/pi] sqrt(K K')/k2^(3/2)

       badParams = ((infParam%consts(1).le.0._kp) &
            .or.(infParam%consts(1).gt.1._kp) &
            .or.(infParam%consts(2).le.0._kp) &
            .or.(infParam%consts(2).gt.1._kp))

       if (badParams) then
          write(*,*)'model name: ',infParam%name
          write(*,*)'consts = ',infParam%consts(1:3)
          stop 'dual inflation: improper params'
       endif
!f
       matterParam(1) = infParam%consts(1)
!lambda
       matterParam(2) = infParam%consts(2)


    case ('mukhai','vfmi')
!Mukhanov inflation: VFMI

!U = M^4 {1 - c2/2/[1 + (2-c1)/[2 sqrt(3 c1)] F]^(2c1/(2-c1))} *
!       exp[ 3 c2/(1-c1) ({1 + (2-c1)/[2 sqrt(3 c1)] F}^[2(1-c1)/(2-c1)] - 1)]

       badParams = ((infParam%consts(1).le.0._kp) &
            .or.(infParam%consts(2).le.0._kp) &
            .or.(infParam%consts(3).le.0._kp))
            
       if (badParams) then
          write(*,*)'model name: ',infParam%name
          write(*,*)'consts = ',infParam%consts(1:3)
          stop 'Mukhanov inflation: improper params'
       endif

!M4
       matterParam(1) = infParam%consts(1)
!alpha
       matterParam(2) = infParam%consts(2)
!beta
       matterParam(3) = infParam%consts(3)



    case ('axhtop','ahi')
!axion hilltop inflation: AHI

!U = c1^4[uplift - 2 cos(F/c2) + (pi - F/c2) sin(F/c2)]
!uplift = 4.820572476962922009964861354665487247290442049596

       badParams = ((infParam%consts(1).le.0._kp) &
            .or.(infParam%consts(2).le.0._kp))

       if (badParams) then
          write(*,*)'model name: ',infParam%name
          write(*,*)'consts = ',infParam%consts(1:2)
          stop 'axion hilltop inflation: improper params'
       endif

!M4
       matterParam(1) = infParam%consts(1)
!phi0
       matterParam(2) = infParam%consts(2)


    case ('sbkahi','sbki')
!symmetry breaking Kahler inflation: SBKI

!U = c1^4 F^2 exp[c2 F^2 + c2^2/6 F^4]

       badParams = ((infParam%consts(1).le.0._kp))

       if (badParams) then
          write(*,*)'model name: ',infParam%name
          write(*,*)'consts = ',infParam%consts(1:2)
          stop 'symmetry breaking Kahler inflation: improper params'
       endif

!M4
       matterParam(1) = infParam%consts(1)
!alpha
       matterParam(2) = infParam%consts(2)



    case ('sduali','sdi')
!S-dual inflation: SDI

!U = c1^4/cosh(F/c2)

       badParams = ((infParam%consts(1).le.0._kp) &
            .or.(infParam%consts(2).le.0._kp))

       if (badParams) then
          write(*,*)'model name: ',infParam%name
          write(*,*)'consts = ',infParam%consts(1:2)
          stop 'S-dual inflation: improper params'
       endif
       
!M4
       matterParam(1) = infParam%consts(1)
!phi0
       matterParam(2) = infParam%consts(2)

!fieldstop required for sduali


    case ('fibrei','fi')
!fiber inflation: FI

!U = c1^4 {3 - c2/c3 + [1+2/3 c2] exp[-4/sqrt(3) F] - 4(1+c2/6) exp[-F/sqrt(3)] + c2/c3 exp[2 c3/sqrt(3) F]}

       badParams = ((infParam%consts(1).le.0._kp) &
            .or.(infParam%consts(2).gt.1._kp) &
            .or.(infParam%consts(3).le.0._kp))

       if (badParams) then
          write(*,*)'model name: ',infParam%name
          write(*,*)'consts = ',infParam%consts(1:3)
          stop 'fiber inflation: improper params'
       endif

!M4
       matterParam(1) = infParam%consts(1)

!delta
       matterParam(2) = infParam%consts(2)

!n+1
       matterParam(3) = infParam%consts(3)


    case ('scaabi','sabi')
!superconformal alpha attractor B inflation: SABI

!U = c1^4 [tanh(F/c3)]^(c2)

       badParams = (infParam%consts(1).le.0._kp) &
            .or. (infParam%consts(2).le.0._kp) &
            .or. (infParam%consts(3).le.0._kp)

       if (badParams) then
          write(*,*)'model name: ',infParam%name
          write(*,*)'consts = ',infParam%consts(1:3)
          stop 'superconformal alpha attractor B inflation: improper params'
       endif
!M4
       matterParam(1) = infParam%consts(1)
!2 n
       matterParam(2) = infParam%consts(2)

!sqrt(6 alpha)
       matterParam(3) = infParam%consts(3)**0.25_kp


    case ('scaaci','saci')
!superconformal alpha attractor C inflation: SACI

!U = c1^4 {tanh(F/c3)/[1 + tanh(F/c3)]}^c2

       badParams = (infParam%consts(1).le.0._kp) &
            .or. (infParam%consts(2).le.0._kp) &
            .or. (infParam%consts(3).le.0._kp)

       if (badParams) then
          write(*,*)'model name: ',infParam%name
          write(*,*)'consts = ',infParam%consts(1:3)
          stop 'superconformal alpha attractor C inflation: improper params'
       endif

!M4
       matterParam(1) = infParam%consts(1)
!2 n
       matterParam(2) = infParam%consts(2)

!sqrt(6 alpha)
       matterParam(3) = infParam%consts(3)**0.25_kp



    case ('hyperb','hbi')
!hyperbolic inflation: HBI
       
!U = c1^4 sinh(F/c3)^c2

       badParams = (infParam%consts(1).le.0._kp) &
            .or. (infParam%consts(2).le.0._kp) &
            .or. (infParam%consts(3).le.infParam%consts(2)/sqrt(2._kp))

       if (badParams) then
          write(*,*)'model name: ',infParam%name
          write(*,*)'consts = ',infParam%consts(1:3)
          stop 'hyperbolic inflation: improper params'
       endif

!M4
       matterParam(1) = infParam%consts(1)
!n       
       matterParam(2) = infParam%consts(2)
!mu
       matterParam(3) = infParam%consts(3)
       


    case ('smearh','shi')
!smeared higgs inflation: SHI

!U = c1^4 { [1-(F/c3)^2]^2 + c2 (F/c3)^4 [ln(F/c3)-1/4] + alpha/4 }

       badParams = (infParam%consts(1).le.0._kp) &
            .or. (infParam%consts(2).le.0._kp) &
            .or. (infParam%consts(3).le.0._kp)
       
       if (badParams) then
          write(*,*)'model name: ',infParam%name
          write(*,*)'consts = ',infParam%consts(1:3)
          stop 'smeared higgs inflation: improper params'
       endif

!M4
       matterParam(1) = infParam%consts(1)
!alpha       
       matterParam(2) = infParam%consts(2)
!mu
       matterParam(3) = infParam%consts(3)
       

    case ('rcplat','rcpi')
!radiatively corrected plateau inflation

!U = c1^4 F^c2 [ 1 + c3 ln(F) + c4 ln(F)^2 ]

       badParams = (infParam%consts(1).le.0._kp) &
            .or. (infParam%consts(2).le.0._kp) &
            .or. (infParam%consts(4).le.0._kp) &
            .or. (infParam%consts(3)**2.gt.4._kp*infParam%consts(4))

       if (badParams) then
          write(*,*)'model name: ',infParam%name
          write(*,*)'consts = ',infParam%consts(1:4)
          stop 'radiatively corrected plateau inflation: improper params'
       endif
       
!M4
       matterParam(1) = infParam%consts(1)
!p
       matterParam(2) = infParam%consts(2)
!alpha^1/4
       matterParam(3) = sign(abs(infParam%consts(3))**0.25_kp,infParam%consts(3))
!beta^1/4
       matterParam(4) = infParam%consts(4)**0.25_kp

    case ('f-term')
!F-term inflation

!2 params kappa and M
       
       matterParam(1) = 1._kp
       matterParam(2) = infParam%consts(2)
       matterParam(5) = infParam%consts(3)


#endif

    case default
       write(*,*)'model name: ',infParam%name
       stop 'set_infbg_param: no such a model'

    end select

  
    if (display) then
       write(*,*)'set_infbg_param: model is ',infParam%name
!       write(*,*)'matterParam = ',matterParam
    endif


!update static potential parameters
    call set_potential_param(matterParam, infParam%name)
    set_infbg_param = .true.
    
  end function set_infbg_param


end module infbgmodel
