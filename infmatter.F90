module infmatter
  use infprec, only : kp, lenshort
  implicit none

  private


  logical, parameter :: display = .true.

!number of matter sector field

  integer, parameter :: matterNum = 1



!Some compilation directives intented for speed optimization. If
!nothing is defined, the matter potential has its maximal number of
!parameters and is parametrized as:
!
! U = [(p1 + p4 ln F) F^p2 + p3]^p5
!   + [p6 + p7 exp(p8 F) + p9 cos(p10 F + p11)] F^p12
!   + p13 F^p14
!
! Defining PP5 truncates this expression to the first five, PP12 to
! the first 12.  Defining PPNAME requires to explicitly code the
! potential, its first and second derivative for each case models
! triggered by potName.
!
! For matterParam vector as input the potential params are set to:
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



#ifdef PP5
  integer, parameter :: potParamNum = 5
#elif defined(PP12)
  integer, parameter :: potParamNum = 12
#else
  integer, parameter :: potParamNum = 14
#endif
  

  real(kp), dimension(potParamNum), save :: potParam = 0._kp
  character(len=lenshort), save :: potName = 'x'
!$omp threadprivate(potParam, potName)
 
  public potParamNum, matterNum

  public set_potential_param, matter_potential, deriv_matter_potential
  public deriv_ln_matter_potential, deriv_second_matter_potential

contains


  subroutine set_potential_param(matterParam, cname)
    implicit none
    real(kp), dimension(:), intent(in) :: matterParam
    character(len=*), intent(in), optional :: cname

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

#ifdef PPNAME
    potName(1:lenshort) = cname(1:lenshort)
#endif

    if (display) write(*,*)'set_potential_param: potParam = ',potParam



  end subroutine set_potential_param



  
  function matter_potential(matter)
    implicit none
    real(kp) :: matter_potential
    real(kp), dimension(matterNum) :: matter

    real(kp) :: chi,lnchi

    chi = matter(1)

#ifndef PPNAME

    if (chi.gt.0._kp) then
       lnchi = log(chi)
    else
       lnchi = -huge(1._kp)
    endif

! From p1 to p5
    matter_potential = (potParam(3) + &
         (potParam(1) + potParam(4)*lnchi)*chi**potParam(2) &
         )**potParam(5)

#ifndef PP5
! From p6 to p12
    matter_potential = matter_potential &
         + ( potParam(6) + potParam(7)*exp(potParam(8)*chi) &
         + potParam(9)*cos(potParam(10)*chi+potParam(11)) ) * chi**(potParam(12))

! From p13 to p14
#ifndef PP12
    matter_potential = matter_potential + potParam(13)*chi**potParam(14)
#endif
#endif


#else

    select case (potName)

       case ('largef')
          matter_potential = potParam(1) * chi**potParam(2)
          
       case default
          stop 'matter_potential: model not found!'

    end select

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

    real(kp) :: chi, lnchi, expchi

    chi = matter(1)

#ifndef PPNAME
       
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
    
    expchi = exp(chi*potParam(8))

    deriv_matter_potential(1) = deriv_matter_potential(1) &
         + chi**(-1 + potParam(12))*(potParam(6) + expchi*potParam(7) &
         + Cos(chi*potParam(10)+potParam(11))*potParam(9))*potParam(12) &
         + chi**potParam(12)*(expchi*potParam(7)*potParam(8) &
         - potParam(9)*potParam(10)*Sin(chi*potParam(10)+potParam(11)))

#ifndef PP12
    deriv_matter_potential(1) = deriv_matter_potential(1) &
         + potParam(13)*potParam(14)*chi**(potParam(14)-1._kp)
#endif
#endif


#else

    select case (potName)

    case ('largef')
       deriv_matter_potential(1) = potParam(1)*potParam(2)*chi**(potParam(2)-1._kp)

    case default
       stop 'deriv_matter_potential: model not found!'

    end select


#endif

  end function deriv_matter_potential
  


  function deriv_second_matter_potential(matter)
    implicit none
    real(kp), dimension(matterNum,matterNum) :: deriv_second_matter_potential
    real(kp), dimension(matterNum) :: matter

    real(kp) :: chi, lnchi, expchi, coschi

    chi = matter(1)
   
#ifndef PPNAME
     
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

    expchi = exp(chi*potParam(8))
    coschi = cos(chi*potParam(10)+potParam(11))

    deriv_second_matter_potential(1,1) = deriv_second_matter_potential(1,1) &
         + chi**potParam(12)*(expchi*potParam(7)*potParam(8)**2 &
         - coschi*potParam(9)*potParam(10)**2) &
         + chi**(-2 + potParam(12))*(potParam(6)+expchi*potParam(7)&
         + coschi*potParam(9))*(-1+potParam(12))*potParam(12)&
         + 2*chi**(-1 + potParam(12))*potParam(12)*(expchi &
         *potParam(7)*potParam(8) - potParam(9)*potParam(10) &
         *Sin(chi*potParam(10)+potParam(11)))

#ifndef PP12
    deriv_second_matter_potential(1,1) = deriv_second_matter_potential(1,1) &
         + potParam(13)*potParam(14)*(potParam(14)-1._kp)*chi**(potParam(14)-2._kp)
#endif
#endif


#else

    select case (potName)

    case ('largef')
       deriv_second_matter_potential(1,1) = potParam(1)*potParam(2)*(potParam(2)-1._kp) &
            * chi**(potParam(2)-2._kp)

    case default
       stop 'deriv_second_matter_potential: model not found!'

    end select

#endif


  end function deriv_second_matter_potential




  function deriv_ln_matter_potential(matter)    
    implicit none
    real(kp), dimension(matterNum) :: deriv_ln_matter_potential
    real(kp), dimension(matterNum) :: matter

    real(kp) :: matterPotential

    matterPotential = matter_potential(matter)

    if (matterPotential.ne.0._kp) then
       deriv_ln_matter_potential = deriv_matter_potential(matter)/matterPotential      
    else
       stop 'infmatter:deriv_ln_matter_potential: matter_potential vanishes!'
    endif

  end function deriv_ln_matter_potential

  



end module infmatter
