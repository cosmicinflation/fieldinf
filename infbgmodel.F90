module infbgmodel
  use infprec, only : kp, lenshort
  use infmatter, only : potParamNum, matterNum
  use infdilaton, only : confParamNum, dilatonNum
  implicit none

  private
 
  logical, parameter :: display = .true.

  
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
         .and. all(infparamA%matters == infparamA%matters) &
         .and. all(infparamA%consts == infparamB%consts))

  end function infparam_equal




  function infparam_unequal(infparamA, infparamB)
    implicit none
    type(infbgparam), intent(in) :: infparamA, infparamB
    logical :: infparam_unequal

    infparam_unequal = ((infparamA%name /= infparamB%name) &
         .or. any(infparamA%conforms /= infParamB%conforms) &
         .or. any(infparamA%matters /= infparamA%matters) &
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
! p3 = m3^4
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
!
!The matterParams mi are set from the ci params according to the model
!under scrutiny. Only the ci (infparam%consts) are public.
!
! m15=c15
!
! is a field value that bounds the initial field values (ex, the throat size for kklt)
!
! m16=c16
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


  select case (infParam%name)


    case ('largef')

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



    case ('smallf')

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


       

    case ('hybrid')

! U = c1^4 [1 + (F/c3)^c2]

       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp) &
            .or.(infParam%consts(3).le.0._kp))            

       

        if  (badParams) then
          write(*,*)'model name: ',infParam%name
          write(*,*)'consts = ',infParam%consts(1:infParamNum)
          stop 'hybrid models: improper params'
       endif
          
       matterParam(1) = infParam%consts(1) &
               /(infParam%consts(3)**(infParam%consts(2)/4._kp))
       matterParam(2) = infParam%consts(2)
       matterParam(3) = infParam%consts(1)
       matterParam(4) = 0._kp
       matterParam(5) = 1._kp



    case ('runmas')

! U = c1^4 { 1 + c4[1/c2 - ln(F/c3)] F^c2 }

       badParams = ((infParam%consts(3).le.0._kp).or.(infParam%consts(3).gt.1._kp) &
            .or.(infParam%consts(2).le.0._kp) &
            .or.(infParam%consts(1).le.0._kp))
       
       if  (badParams) then
          write(*,*)'model name: ',infParam%name
          write(*,*)'consts = ',infParam%consts(1:4)
          stop 'running mass: improper params'
       endif
       
       kpbuffer = infParam%consts(4)/infParam%consts(2) &
               + infParam%consts(4)*log(infParam%consts(3))
       matterParam(1) = infParam%consts(1)*kpbuffer**0.25_kp

       matterParam(2) = infParam%consts(2)
       matterParam(3) = infParam%consts(1)

       kpbuffer = infParam%consts(4)
       matterParam(4) = -infParam%consts(1)*kpbuffer**0.25_kp
       
       matterParam(1) = infParam%consts(1) * sign(abs(kpbuffer)**0.25,kpbuffer)

       kpbuffer = infParam%consts(4)
       matterParam(4) = - infParam%consts(1) * sign(abs(kpbuffer)**0.25,kpbuffer)

       matterParam(5) = 1._kp

!fieldstop value required (checkout infbounds.f90)
       


    case('kklmmt')

! U = c1^4 * [1 - (F/c3)^(-c2)] with c2 > 0 for c5=1 or
! U = c1^4 / [1 + (F/c3)^(-c2)] with c2 > 0 for c5=-1
!
! c6 is Phi_string related to the flux number N
!
! c7 is PhiUv related to the brane tension of the string coupling


       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp) &
            .or.(infParam%consts(3).le.0._kp) &
            .or.(abs(infParam%consts(5)).ne.1._kp) &
            .or.(infParam%consts(matterParamNum-1).lt.0._kp))

       if (badParams) then
          write(*,*)'model name: ',infParam%name
          write(*,*)'consts = ',infParam%consts(1:5)
          stop 'kklmmt: improper params'
       endif


       matterParam(1) = - sign(infParam%consts(1)**infParam%consts(5) &
            /(infParam%consts(3)**(-infParam%consts(2)/4._kp)) &
            ,infParam%consts(5))
       matterParam(2) = -infParam%consts(2)
       matterParam(3) = infParam%consts(1)**infParam%consts(5)
       matterParam(4) = 0._kp
       matterParam(5) = infParam%consts(5)

!those are necessary model parameters: checkout infbounds.f90
!fieldUv: prevents brane out of the throat
!flux number N
       
            


    case ('rcquad')
! U = c1^4 F^4 [1 - c2 ln(F)]
          
       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp))

       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:2)
          stop 'radiative corrected quadratic field: improper params'
       endif

       matterParam(1) = infParam%consts(1)
       matterParam(2) = 4._kp
       matterParam(3) = 0._kp
       matterParam(4) = -infParam%consts(1) &
            *sign(abs(infParam%consts(2))**0.25_kp,infParam%consts(2))
       matterParam(5) = 1._kp


    case ('gswli')
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
          

   
#ifndef PP5

    case ('mixlf')
! U = c1^4 F^c2 [1 + c3 F^c4]

       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp) &
            .or.(infParam%consts(3).le.0._kp).or.(infParam%consts(4).le.0._kp))

       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:infParamNum)
          stop 'mixed field: improper params'
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



    case ('lfcorr')
!U = c1^4 F^c2 [1 - c3 F^c4 ln(F)]       


       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp).or. &
            (infParam%consts(3).le.0._kp).or.(infParam%consts(4).le.0._kp))

       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:4)
          stop 'radiative corrected large field: improper params'
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


    case ('rcmass')
!U = c1^4 F^2 [1 - c2 F^2 ln F]

!c2  is 2 x alpha

       badParams = ((infParam%consts(1).le.0._kp).or.infParam%consts(2).le.0._kp)

       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:2)
          stop 'radiative corrected massive: improper params'
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


    case ('natinf')
!U = c1^4 [1 + c2 cos(F/c3)]

!c2 is +-1 for plus sign natural inflation or -1 for minus signa
!natural inflation

       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(3).le.0._kp) &
            .or.(abs(infParam%consts(2)).ne.1._kp))

       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:3)
          stop 'natural with plus or minus: improper params'
       endif
       
       matterParam(1:4) = 0._kp
       matterParam(5) = 2._kp
       matterParam(6) = infparam%consts(1)
       matterParam(7) = 0._kp
       matterParam(8) = 0._kp
       matterParam(9) = infparam%consts(1) &
            * sign(abs(infParam%consts(2))**0.25_kp,infParam%consts(2))
       matterParam(10) = 1._kp/infParam%consts(3)
       matterParam(11) = 0._kp
       matterParam(12) = 0._kp
       

    case ('exsusy')
!U = c1^4 [1 - exp(-c2 F)]


       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp))

       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:2)
          stop 'expoential susy: improper params'
       endif

       matterParam(1:4) = 0._kp
       matterParam(5) = 2._kp
       matterParam(6) = infParam%consts(1)
       matterParam(7) = -infParam%consts(1)
       matterParam(8) = -infParam%consts(2)
       matterParam(9:12) = 0._kp


    case ('powlaw')
!U = c1^4 exp[-c2 F]
       
       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp))

       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:2)
          stop 'power law: improper params'
       endif

       matterParam(1:4) = 0._kp
       matterParam(5) = 2._kp
       matterParam(6) = 0._kp
       matterParam(7) = infParam%consts(1)
       matterParam(8) = -infParam%consts(2)
       matterParam(9:12) = 0._kp

!fieldstop value required (checkout infbounds.f90)
     


    case ('hfline')
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


    case ('interm')
! U = c1^4 F^(-c2) [1 - c2^(2/6) F^(-2)]

       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp))           

       if (badParams) then          
          write(*,*)'model name: ',infParam%name          
          write(*,*)'consts = ',infParam%consts(1:infParamNum)
          stop 'mixed field: improper params'
       endif


       matterParam(1) = infParam%consts(1)
       matterParam(2) = -infParam%consts(2)
       matterParam(3) = 0._kp
       matterParam(4) = 0._kp
       matterParam(5) = 1._kp
       matterParam(6) = -infParam%consts(1) * (infParam%consts(2)**2/6._kp)**0.25_kp
       matterParam(7:11) = 0._kp
       matterParam(12) = -infParam%consts(2) - 2._kp


#ifndef PP12

    case ('kahmod')
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


#endif
#endif


    case default
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
