module infbgmodel
  use infprec, only : kp
  implicit none

  private
 
  logical, parameter :: display = .true.

  
  integer, parameter :: matterNum = 1
  integer, parameter :: dilatonNum = 1  
  integer, parameter :: fieldNum = matterNum + dilatonNum

#ifdef PP5
  integer, parameter :: potParamNum = 5
#elif defined(PP12)
  integer, parameter :: potParamNum = 12
#else
  integer, parameter :: potParamNum = 14
#endif



  integer, parameter :: matterParamNum = potParamNum + 2
  integer, parameter :: dilatonParamNum = 0
  integer, parameter :: infParamNum = matterParamNum + dilatonParamNum



  type infbgparam
!label identifier
     character(len=6) :: name
!coupling constants [mattercouplings,dilatoncouplings]
     real(kp), dimension(infParamNum) :: consts     
![dilaton degrees of freedom]
     real(kp), dimension(dilatonNum) :: conforms
!matterfield values
     real(kp), dimension(matterNum) :: matters
  end type infbgparam



  real(kp), dimension(matterParamNum), save :: matterParam = 0._kp
  real(kp), dimension(potParamNum), save :: potParam = 0._kp


  interface operator (==)     
     module procedure infparam_equal
  end interface


  interface operator (/=)
     module procedure infparam_unequal
  end interface


  public infbgparam
  public operator(==),operator(/=)

  public matterParam
  public fieldNum, dilatonNum, matterNum
  public infParamNum, matterParamNum, dilatonParamNum

  public set_infbg_param

  public matter_potential
  public deriv_matter_potential, deriv_second_matter_potential

  public conformal_factor_square 
  public conformal_first_gradient, conformal_second_gradient

  public metric, metric_inverse, deriv_metric


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
    implicit none
    logical :: set_infbg_param
    type(infbgparam), intent(in) :: infParam    

    logical :: badParams = .true.

    real(kp) :: kpbuffer

    set_infbg_param=.false.

!the matter potential is parametrized as
!
! U = [(p1 + p4 ln F) F^p2 + p3]^p5
!   + [p6 + p7 exp(p8 F) + p9 cos(p10 F + p11)] F^p12
!   + p13 F^p14
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
!Other notations: M=c1; mu=c3, p=c2; nu=c4

  

  select case (infParam%name)


    case ('largef')

! U = c1^4 F^c2       

       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp) &
            .or.(infParam%consts(3).ne.0._kp).or.(infParam%consts(5).ne.1._kp))

       if (badParams) then          
          write(*,*)'model name: ',infParam%name
          write(*,*)'consts = ',infParam%consts(1:infParamNum)
          stop 'large field: improper params'
       endif


       matterParam(1) = infParam%consts(1)
       matterParam(2) = infParam%consts(2)
       matterParam(3) = 0._kp
       matterParam(4) = 0._kp
       matterParam(5) = 1._kp

!fieldUv limit
       matterParam(matterParamNum-1) = infParam%consts(matterParamNum-1)
!fieldStop value
       matterParam(matterParamNum) = infParam%consts(matterParamNum)


    case ('smallf')

! U = c1^4 [1 - (F/c3)^c2]

       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp) &
            .or.(infParam%consts(3).le.0._kp).or.(infParam%consts(5).ne.1._kp))

       if (badParams) then
          write(*,*)'model name: ',infParam%name
          write(*,*)'consts = ',infParam%consts(1:infParamNum)
          stop 'small field: improper params'
       endif

       matterParam(1) = - infParam%consts(1) &
            /(infParam%consts(3)**(infParam%consts(2)/4._kp))
       matterParam(2) = infParam%consts(2)
       matterParam(3) = infParam%consts(1)
       matterParam(4) = 0._kp
       matterParam(5) = 1._kp

!fieldUv limit
       matterParam(matterParamNum-1) = infParam%consts(matterParamNum-1)
!fieldStop value
       matterParam(matterParamNum) = infParam%consts(matterParamNum)

            
       if (maxval(infParam%matters/infParam%conforms(1)) &
            .gt.infParam%consts(3)) then
          write(*,*)'set_infbg_param: not small fields initially'
          write(*,*)'model name: ',infParam%name         
          write(*,*)'matterScale = ',infParam%consts(3)
          write(*,*)'matterField = ',maxval(infParam%matters)          
       endif


       

    case ('hybrid')

! U = c1^4 [1 + (F/c3)^c2]

       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp) &
            .or.(infParam%consts(3).le.0._kp) &
            .or.(infParam%consts(5).ne.1._kp))

       

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

!fieldUv limit
       matterParam(matterParamNum-1) = infParam%consts(matterParamNum-1)
!fieldstop value
       matterParam(matterParamNum) =infParam%consts(matterParamNum)

 

    case ('runmas')

! U = c1^4 { 1 + c4[1/c2 - ln(F/c3)] F^c2 }

       badParams = ((infParam%consts(3).le.0._kp).or.(infParam%consts(3).gt.1._kp) &
            .or.(infParam%consts(2).le.0._kp) &
            .or.(infParam%consts(1).le.0._kp) &
            .or.(infParam%consts(5).ne.1._kp))


       if  (badParams) then
          write(*,*)'model name: ',infParam%name
          write(*,*)'consts = ',infParam%consts(1:infParamNum)
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

!fieldUv limit
       matterParam(matterParamNum-1) = infParam%consts(matterParamNum-1)
!fieldstop value
       matterParam(matterParamNum) = infParam%consts(matterParamNum)



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
          write(*,*)'consts = ',infParam%consts(1:infParamNum)
          stop 'kklmmt: improper params'
       endif


       matterParam(1) = - sign(infParam%consts(1)**infParam%consts(5) &
            /(infParam%consts(3)**(-infParam%consts(2)/4._kp)) &
            ,infParam%consts(5))
       matterParam(2) = -infParam%consts(2)
       matterParam(3) = infParam%consts(1)**infParam%consts(5)
       matterParam(4) = 0._kp
       matterParam(5) = infParam%consts(5)

!fieldUv: prevents brane out of the throat
       matterParam(matterParamNum-1) = infParam%consts(matterParamNum-1)
!flux number N
       matterParam(matterParamNum) = infParam%consts(matterParamNum)
            

#ifndef PP5
    case ('mixinf')
! U = c1^4 [F^c2 + c6 F^c12]      

       badParams = ((infParam%consts(1).le.0._kp).or.(infParam%consts(2).le.0._kp) &
            .or.(infParam%consts(3).ne.0._kp).or.(infParam%consts(5).ne.1._kp) &
            .or.(infParam%consts(6).le.0._kp).or.(infParam%consts(12).le.0._kp) &
            .or.(infParam%consts(7).ne.0._kp).or.(infParam%consts(8).ne.0._kp) &
            .or.(infParam%consts(9).ne.0._kp).or.(infParam%consts(10).ne.0._kp) &
            .or.(infParam%consts(11).ne.0._kp))

       if (badParams) then          
          write(*,*)'model name: ',infParam%name
          write(*,*)'infParamNum = ',infParamNum
          write(*,*)'consts = ',infParam%consts(1:infParamNum)
          stop 'mixed field: improper params'
       endif


       matterParam(1) = infParam%consts(1)
       matterParam(2) = infParam%consts(2)
       matterParam(3) = 0._kp
       matterParam(4) = 0._kp
       matterParam(5) = 1._kp
       matterParam(6) = infParam%consts(1) &
            *sign(abs(infParam%consts(6))**0.25_kp,infParam%consts(6))
       matterParam(7:11) = 0._kp
       matterParam(12) = infParam%consts(12)

!fieldUv limit
       matterParam(matterParamNum-1) = infParam%consts(matterParamNum-1)
!fieldStop value
       matterParam(matterParamNum) = infParam%consts(matterParamNum)


#endif



    case default
       stop 'set_infbg_param: no such a model'

    end select

  
    if (display) then
       write(*,*)'set_infbg_param: model is ',infParam%name
!       write(*,*)'matterParam = ',matterParam
    endif


!update static potential parameters
    call set_potential_param()
    set_infbg_param = .true.
    
  end function set_infbg_param




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!full numerical integration functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!matter sector: potentials
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine set_potential_param()
    implicit none

    potParam(1) = sign(matterParam(1)**4,matterParam(1)) !/matterParam(3)**4
    potParam(2) = matterParam(2)
    potParam(3) = matterParam(3)**4
    potParam(4) = sign(matterParam(4)**4,matterParam(4)) !/matterParam(3)**4
    potParam(5) = matterParam(5)

#ifndef PP5
    potParam(6) = sign(matterParam(6)**4,matterParam(6))
    potParam(7) = sign(matterParam(7)**4,matterParam(7))
    potParam(8) = matterParam(8)
    potParam(9) = sign(matterParam(9)**4,matterParam(9))
    potParam(10) = matterParam(10)
    potParam(11) = matterParam(11)
    potParam(12) = matterParam(12)


#ifndef PP12
    potParam(13) = sign(matterParam(13)**4,matterParam(13))
    potParam(14) = matterParam(14)
#endif
#endif

    if (display) write(*,*)'set_potential_param: potParam = ',potParam

  end subroutine set_potential_param



  
  function matter_potential(matter)
    implicit none
    real(kp) :: matter_potential
    real(kp), dimension(matterNum) :: matter

    real(kp) :: chi,lnchi

    chi = matter(1)

    if (chi.gt.0._kp) then
       lnchi = log(chi)
    else
       lnchi = -huge(1._kp)
    endif

! From p1 to p5
    matter_potential = (potParam(3) + &
         (potParam(1) + potParam(4)*lnchi)*chi**potParam(2) &
         )**potParam(5)

! From p6 to p12
#ifndef PP5
    matter_potential = matter_potential &
         + ( potParam(6) + potParam(7)*exp(potParam(8)*chi) &
         + potParam(9)*cos(potParam(10)*chi+potParam(11)) ) * chi**(potParam(12))

! From p13 to p14
#ifndef PP12
    matter_potential = matter_potential + potParam(13)*chi**potParam(14)
#endif
#endif

    if (matter_potential.lt.0._kp) then
       write(*,*)'matterField = ',chi
       write(*,*)'potParam = ',potParam
       write(*,*)'matter_potential = ',matter_potential
       stop 'infbgmodel: matter_potential < 0!'
    endif
  end function matter_potential
 

 
  function deriv_matter_potential(matter)
    implicit none
    real(kp), dimension(matterNum) :: deriv_matter_potential
    real(kp), dimension(matterNum) :: matter

    real(kp) :: chi, lnchi

    chi = matter(1)

    if (chi.gt.0._kp) then
       lnchi = log(chi)
    else
       lnchi = -huge(1._kp)
    endif

    deriv_matter_potential(1) = potParam(5)*chi**(potParam(2)-1._kp) &
         * (potParam(1)*potParam(2) + potParam(4) + potParam(4)*potParam(2)*lnchi) &
         * ( (potParam(1) + potParam(4)*lnchi)*chi**potParam(2) + potParam(3) ) &
         ** (potParam(5)-1._kp)

#ifndef PP5
    deriv_matter_potential(1) = deriv_matter_potential(1) &
         + chi**(-1 + potParam(12))*(potParam(6) + exp(chi*potParam(8))*potParam(7) &
         + Cos(chi*potParam(10)+potParam(11))*potParam(9))*potParam(12) &
         + chi**potParam(12)*(exp(chi*potParam(8))*potParam(7)*potParam(8) &
         - potParam(9)*potParam(10)*Sin(chi*potParam(10)+potParam(11)))

#ifndef PP12
    deriv_matter_potential(1) = deriv_matter_potential(1) &
         + potParam(13)*potParam(14)*chi**(potParam(14)-1._kp)
#endif
#endif

  end function deriv_matter_potential
  


  function deriv_second_matter_potential(matter)
    implicit none
    real(kp), dimension(matterNum,matterNum) :: deriv_second_matter_potential
    real(kp), dimension(matterNum) :: matter

    real(kp) :: chi, lnchi

    chi = matter(1)
     
    if (chi.gt.0._kp) then
       lnchi = log(chi)
    else
       lnchi = -huge(1._kp)
    endif

    deriv_second_matter_potential(1,1) = chi**(-2._kp + potParam(2)) * (potParam(3) &
         + chi**potParam(2)*(potParam(1) + lnchi*potParam(4)))**(-2._kp &
         + potParam(5))*((potParam(1)*(-1._kp + potParam(2))*potParam(2) &
         + lnchi*(-1._kp + potParam(2))*potParam(2)*potParam(4) &
         + (-1._kp + 2._kp*potParam(2))*potParam(4))*(potParam(3)+chi**potParam(2) &
         * (potParam(1) + lnchi*potParam(4))) + chi**potParam(2) &
         * (potParam(1)*potParam(2)+potParam(4)+lnchi*potParam(2)*potParam(4))**2._kp &
         * (-1._kp + potParam(5)))*potParam(5)

#ifndef PP5
    deriv_second_matter_potential(1,1) = deriv_second_matter_potential(1,1) &
         + chi**potParam(12)*(exp(chi*potParam(8))*potParam(7)*potParam(8)**2 &
         - Cos(chi*potParam(10)+potParam(11))*potParam(9)*potParam(10)**2) &
         + chi**(-2 + potParam(12))*(potParam(6)+exp(chi*potParam(8))*potParam(7)&
         + Cos(chi*potParam(10)+potParam(11))*potParam(9))*(-1+potParam(12))*potParam(12)&
         + 2*chi**(-1 + potParam(12))*potParam(12)*(exp(chi*potParam(8))&
         *potParam(7)*potParam(8) - potParam(9)*potParam(10) &
         *Sin(chi*potParam(10)+potParam(11)))

#ifndef PP12
    deriv_second_matter_potential(1,1) = deriv_second_matter_potential(1,1) &
         + potParam(13)*potParam(14)*(potParam(14)-1._kp)*chi**(potParam(14)-2._kp)
#endif
#endif

  end function deriv_second_matter_potential





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!gravity sector: metric on the sigma-model field manifold
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




 function metric(field)
    implicit none
    real(kp), dimension(fieldNum), intent(in) :: field
    real(kp), dimension(fieldNum,fieldNum) :: metric
       
    real(kp) :: confSquare
    real(kp), dimension(dilatonNum) :: dilaton
    integer :: i

    metric = 0._kp

    dilaton = field(matterNum+1:fieldNum)
    confSquare = conformal_factor_square(dilaton)

    do i=1,matterNum
       metric(i,i) = confSquare
    enddo

    do i=1,dilatonNum
       metric(i+matterNum,i+matterNum) = 1._kp
    enddo

  end function metric



  function metric_inverse(field)
    implicit none
    real(kp), dimension(fieldNum) :: field
    real(kp), dimension(fieldNum,fieldNum) :: metricVal, metric_inverse
    real(kp) :: det

    integer :: i

    metric_inverse = 0._kp

    metricVal = metric(field)
       
    det = 1._kp
    do i=1,fieldNum
       det = det*metricVal(i,i)
    enddo

    if (det.ne.0.) then
       do i=1,fieldNum
          metric_inverse(i,i) = 1._kp/metricVal(i,i)
       enddo
    else
       stop 'inverse_metric: singularity in the sigma-model!'
    endif

  end function metric_inverse


  
  function deriv_metric(field)
    implicit none
    real(kp), dimension(fieldNum), intent(in) :: field
    real(kp), dimension(fieldNum,fieldNum,fieldNum) :: deriv_metric

    real(kp) :: confSquare
    real(kp), dimension(dilatonNum) :: dilaton, confFirstGrad
    integer :: i

    deriv_metric = 0._kp

    dilaton = field(matterNum+1:fieldNum)
    confSquare = conformal_factor_square(dilaton)
    confFirstGrad = conformal_first_gradient(dilaton)
    
    do i=1,matterNum
       deriv_metric(i,i,matterNum+1:fieldNum) &
            = 2._kp * confFirstGrad(1:dilatonNum)* confSquare
    enddo
    


  end function deriv_metric


    
  function conformal_factor_square(dilaton)
!A^2
    implicit none
    real(kp), dimension(dilatonNum), intent(in) :: dilaton
    real(kp) :: conformal_factor_square

    conformal_factor_square = 1._kp


  end function conformal_factor_square

  

  function conformal_first_gradient(dilaton)
!alpha_field = Dln(A)/Dfield
    implicit none
    real(kp), dimension(dilatonNum) :: conformal_first_gradient
    real(kp), dimension(dilatonNum) :: dilaton

        
    conformal_first_gradient = 0._kp


  end function conformal_first_gradient


 
  function conformal_second_gradient(dilaton)
!D alpha_field / D field
    implicit none
    real(kp), dimension(dilatonNum,dilatonNum) :: conformal_second_gradient
    real(kp), dimension(dilatonNum) :: dilaton
   
    conformal_second_gradient = 0._kp   


  end function conformal_second_gradient




end module infbgmodel
