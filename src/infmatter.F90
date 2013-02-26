module infmatter
  use infprec, only : kp, lenshort, pi
#if defined (PPNAME) && !defined (NOASPIC)
  include 'libaspic.h'
#endif

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
!   + p13 F^p14 + p15 F^p16 exp(p17 F^p18)
!
! Defining PP5 truncates this expression to the first five, PP12 to
! the first 12.  Defining PPNAME requires to explicitly code the
! potential, its first and second derivative for each case models
! triggered by potName. But wait, that's already done in libaspic,
! so cool!
!
! For matterParam vector as input the potential params are set to:
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


#ifdef PP5
  integer, parameter :: potParamNum = 5
#elif defined(PP12)
  integer, parameter :: potParamNum = 12
#else
  integer, parameter :: potParamNum = 18
#endif
  

  real(kp), dimension(potParamNum), save :: potParam = 0._kp
  character(len=lenshort), save :: potName = 'x'
 
!some more easy to use alias for PPNAME
  real(kp), save :: alpha, beta, gam, lambda ,p, q, mu, nu, M4

  real(kp) :: M, kappa, phic
 
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
    potParam(3) = sign(matterParam(3)**4,matterParam(3))
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
    potParam(15) = sign(matterParam(15)**4,matterParam(15))
    potParam(16) = matterParam(16)
    potParam(17) = matterParam(17)
    potParam(18) = matterParam(18)
#endif
#endif

#ifdef PPNAME

    if (present(cname)) then
       potName(1:lenshort) = cname(1:lenshort)
    else
       stop 'set_potential_param: PPNAME defined but not name as an input!'
    endif

    select case (cname)

       case ('largef')
          M4 = potParam(1)
          p = potParam(2)

       case ('rcquad')
          M4 = potParam(1)
          alpha = -potParam(4)/potParam(1)

       case ('smallf')
          M4 = potParam(3)
          p = potParam(2)
          mu = (-potParam(3)/potParam(1))**(1._kp/potParam(2))

       case ('gswli')
          M4 = potParam(3)
          alpha = potParam(4)/potParam(3)

       case ('colwei')
          M4 = potParam(3)
          alpha = potParam(4)/potParam(3)*exp(-potParam(2)*potParam(1)/potParam(4))
          mu = exp(-potParam(1)/potParam(4))

       case ('tdwell')
          M4 = potParam(3)**2
          mu = sqrt(-potParam(3)/potParam(1))

       case ('betexp')
          M4 = potParam(3)**potParam(5)
          beta = 1._kp/potParam(5)
          lambda = -potParam(1)*potParam(5)/potParam(3)

       case ('radiag')
          M4 = 1._kp/potParam(3)
          alpha = potParam(1)/potParam(3)
          
       case ('hybrid')
          M4 = potParam(3)
          p = potParam(2)
          mu = (potParam(3)/potParam(1))**(1._kp/potParam(2))

       case ('dysusy')
          M4 = potParam(3)
          p = -potParam(2)
          mu = (potParam(3)/potParam(1))**(1._kp/potParam(2))

       case ('nszero')
          M4= 1._kp/potParam(3)**2
          alpha = - potParam(1)/potParam(3)
          
       case ('sugrab')
          M4 = potParam(3)
          alpha = -potParam(1)/potParam(3)
          beta = potParam(4)/potParam(3)

       case ('runmas','runma1','runma2','runma3','runma4')
          M4 = potParam(3)          
          nu = -potParam(4)/potParam(3)
          mu = exp(-potParam(1)/potParam(4) - 1._kp/2._kp)
          potName = 'runmas'


       case ('logpot','logpo1','logpo2','logpo3')
          q = potParam(5)
          p = potParam(2)*potParam(5)
          mu = exp(-potParam(1)/potParam(4))
          M4 = potParam(4)**potParam(5) &
               * exp(-potParam(1)/potParam(4)*potParam(2)*potParam(5))
          potName = 'logpot'

#ifndef PP5
       case ('gmixlf')
          M4 = potParam(1)
          p = potParam(2)
          q = potParam(12) - potParam(2)
          alpha = potParam(6)/potParam(1)

       case ('rcmass')
          M4 = potParam(6)
          alpha = -0.5_kp*potParam(4)/potParam(6)
       
       case ('natinf')
          M4 = potParam(6)
          alpha = potParam(9)/potParam(6)
          mu = 1._kp/potParam(10)
        
       case ('exsusy')
          M4 = potParam(6)
          q = -potParam(8)

       case ('powlaw')
          M4 = potParam(7)
          alpha = -potParam(8)

       case ('hfline')
          M4 = potParam(3)**2
          alpha = potParam(1)/potParam(3)

       case ('interm')
          M4 = potParam(1)
          beta = - potParam(2)

       case ('twisti')
          M4 = potParam(3)
          mu = -1._kp/potParam(8)

       case ('nckahi')
          M4= potParam(3)
          alpha = potParam(4)/potParam(3)
          beta = potParam(6)/potParam(3)

       case ('oifold')
          mu = exp(-potParam(1)/potParam(4))
          M4 = potParam(4)**2 * exp(-4*potParam(1)/potParam(4))
          alpha = -potParam(6)/potParam(4)**2

       case ('ssbinf','spsyb1','spsyb2','spsyb3','spsyb4','spsyb5','spsyb6')
          M4 = potParam(3)
          alpha = potParam(1)/potParam(3)
          beta = potParam(6)/potParam(3)
          potName='ssbinf'

#ifndef PP12
       case ('kahmod')
          M4 = potParam(13)
          alpha = -potParam(15)/potParam(13)
          beta = -potParam(17)
          p = potParam(16)
          q = potParam(18)

          if (p.eq.1._kp) potName = 'kahmo1'
          if (p.eq.4._kp/3._kp) potName = 'kahmo2'

       case ('higgsi')
          M4 = potParam(3)

       case ('logmdi','logmd1','logmd2')
          M4 = potParam(15)
          alpha = potParam(16)
          beta = -potParam(17)
          gam = potParam(18)
          potName = 'logmdi'

       case ('nmssmi','gmssmi','rinfpt')
!for convenience only
          mu = potParam(11)

          M4 = potParam(1)*potParam(11)**2
          p = potParam(12)
          q = potParam(14)
          alpha = - potParam(6)*potParam(11)**(potParam(12)-2._kp)/potParam(1)
          beta = potParam(13)*potParam(11)**(potParam(5)-2._kp)/potParam(1)

          

       case ('bsusyb')
          M4 = potParam(7)
          gam = potParam(17)/sqrt(6._kp)

#endif
#endif     

       case ('mhitop')
          M4 = potParam(1)
          mu = potParam(2)

       case ('ricci1','ricci2')
          M4 = potParam(1)
          p = potParam(2)
          q = sqrt(3._kp/2._kp)

       case ('tipinf')
          M4 = potParam(1)
          alpha = potParam(2)
          mu = potParam(3)**0.25_kp

       case ('psenat')
          M4 = potParam(1)
          alpha = potParam(2)
          mu = potParam(3)**0.25_kp

       case ('arctan')
          M4 = potParam(1)
          mu = potParam(2)

       case ('fixnsa','fixnsb','fixnsc')
          M4 = potParam(1)
          alpha = potParam(2)

       case ('fixnsd')
          M4 = potParam(1)
          alpha = potParam(2)
          beta = potParam(3)**0.25_kp

!       case ('f-term')
!          M4 = potParam(1)
!          kappa = potParam(2)
!          mu = potParam(5)
!          lambda = kappa**2  * mu**4
!          phic = sqrt(2._kp) * mu
!          M = 2._kp * mu
          

       case default

          stop 'set_potential_param: not such a model'

       end select



#endif

    if (display) write(*,*)'set_potential_param: potParam = ',potParam



  end subroutine set_potential_param



  
  function matter_potential(matter)
    implicit none
    real(kp) :: matter_potential
    real(kp), dimension(matterNum) :: matter

    real(kp) :: chi,lnchi

!    real(kp) :: phi,psi

    chi = matter(1)
    
!    phi = matter(1)
!    psi=  matter(2)


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
    matter_potential = matter_potential + potParam(13)*chi**potParam(14) &
         + potParam(15)*chi**potParam(16)*exp(potParam(17)*chi**potParam(18))
#endif
#endif

    
#elif !defined(NOASPIC)

    select case (potName)

       case ('largef')
          matter_potential = lfi_norm_potential(chi,p)
          
       case ('gmixlf')
          matter_potential = gmlfi_norm_potential(chi,p,q,alpha)
          
       case ('rcmass')
          matter_potential = rcmi_norm_potential(chi,alpha)

       case ('rcquad')
          matter_potential = rcqi_norm_potential(chi,alpha)

       case ('natinf')
          matter_potential = ni_norm_potential(chi,mu)
          
       case ('exsusy')
          matter_potential = esi_norm_potential(chi,q)

       case ('powlaw')
          matter_potential = pli_norm_potential(chi,alpha)

       case ('kahmo1')
          matter_potential = kmii_norm_potential(chi,alpha)

       case ('hfline')
          matter_potential = hf1i_norm_potential(chi,alpha)

       case ('smallf')
          matter_potential = sfi_norm_potential(chi/mu,p)

       case ('gswli')
          matter_potential = li_norm_potential(chi,alpha)

       case ('interm')
          matter_potential = ii_norm_potential(chi,beta)

       case ('kahmo2')
          matter_potential = kmiii_norm_potential(chi,alpha,beta)

       case ('colwei')
          matter_potential = cwi_norm_potential(chi/mu,alpha,mu)

       case ('higgsi')
          matter_potential = hi_norm_potential(chi)

       case ('twisti')
          matter_potential = twi_norm_potential(chi,mu)

       case ('tdwell')
          matter_potential = dwi_norm_potential(chi/mu)

       case ('mhitop')
          matter_potential = mhi_norm_potential(chi/mu)

       case ('logmdi')
          matter_potential = lmi_norm_potential(chi,gam,beta)

       case ('ricci1')
          matter_potential = rpi1_norm_potential(chi/q,p)

       case ('ricci2')
          matter_potential = rpi2_norm_potential(chi/q,p)

       case ('betexp')
          matter_potential = bei_norm_potential(chi,lambda,beta)

       case ('radiag')
          matter_potential = rgi_norm_potential(chi,alpha)

       case ('nmssmi')
          matter_potential = mssmi_norm_potential(chi/mu,mu)

       case ('rinfpt')
          matter_potential = ripi_norm_potential(chi,alpha)

       case ('gmssmi')
          matter_potential = gmssmi_norm_potential(chi/mu,1.5_kp*alpha,mu)

       case ('bsusyb')
          matter_potential = bsusybi_norm_potential(chi,gam)

       case ('tipinf')
          matter_potential = ti_norm_potential(chi/mu,alpha)

       case ('psenat')
          matter_potential = psni_norm_potential(chi/mu,alpha)

       case ('nckahi')
          matter_potential = ncki_norm_potential(chi,alpha,beta)

       case ('hybrid')
          matter_potential = vhi_norm_potential(chi/mu,p)

       case ('dysusy')
          matter_potential = dsi_norm_potential(chi/mu,p)

       case ('arctan')
          matter_potential = ai_norm_potential(chi/mu)

       case ('fixnsa')
          matter_potential = cnai_norm_potential(chi,alpha)

       case ('fixnsb')
          matter_potential = cnbi_norm_potential(chi,alpha)

       case ('fixnsc')
          matter_potential = cnci_norm_potential(chi,alpha)

       case ('fixnsd')
          matter_potential = cndi_norm_potential(chi,alpha,beta)

       case ('nszero')
          matter_potential = csi_norm_potential(chi,alpha)

       case ('oifold')
          matter_potential = oi_norm_potential(chi/mu,alpha,mu)
          
       case ('sugrab')
          matter_potential = sbi_norm_potential(chi,alpha,beta)

       case ('ssbinf')
          matter_potential = ssbi_norm_potential(chi,alpha,beta)

       case ('runmas')
          matter_potential = rmi_norm_potential(chi/mu,2._kp*nu,mu)

       case ('logpot')
          matter_potential = lpi_norm_potential(chi/mu,p,q)

!       case ('f-term')
!          matter_potential = lambda * ( ( 1._kp - psi**2 / M**2 )**2   &
!               + 2._kp * phi**2 * psi**2 / M**2 / phic**2 &
!               + 16._kp * lambda / M**4 * log(2._kp)/ (4._kp * Pi**2 * phic) * phi) 

       case default
          write(*,*)'name is ',potName
          stop 'matter_potential: model not found!'

    end select
!ensures normalistion x M^4
    matter_potential = M4 * matter_potential

#else
    stop 'matter_potential: FieldInf not built agains libsrmodels!'
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

    real(kp) :: chi, lnchi, expchi, expchi2

!    real(kp) :: phi,psi

    chi = matter(1)

!    phi = matter(1)
!    psi= matter(2)


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

    expchi2 = exp(chi**potParam(18)*potParam(17))

    deriv_matter_potential(1) = deriv_matter_potential(1) &
         + potParam(13)*potParam(14)*chi**(potParam(14)-1._kp) &
         + chi**(-1._kp+potParam(16))*expchi2*potParam(15)*potParam(16)&
         + chi**(-1 + potParam(16) + potParam(18))*expchi2*potParam(15) &
         *potParam(17)*potParam(18)
#endif
#endif


#elif !defined(NOASPIC)

    select case (potName)

       case ('largef')
          deriv_matter_potential(1) = lfi_norm_deriv_potential(chi,p)
          
       case ('gmixlf')
          deriv_matter_potential(1) = gmlfi_norm_deriv_potential(chi,p,q,alpha)
          
       case ('rcmass')
          deriv_matter_potential(1) = rcmi_norm_deriv_potential(chi,alpha)

       case ('rcquad')
          deriv_matter_potential(1) = rcqi_norm_deriv_potential(chi,alpha)

       case ('natinf')
          deriv_matter_potential(1) = ni_norm_deriv_potential(chi,mu)
                
       case ('exsusy')
          deriv_matter_potential(1) = esi_norm_deriv_potential(chi,q)

       case ('powlaw')
          deriv_matter_potential(1) = pli_norm_deriv_potential(chi,alpha)

       case ('kahmo1')
          deriv_matter_potential(1) = kmii_norm_deriv_potential(chi,alpha)

       case ('hfline')
          deriv_matter_potential(1) = hf1i_norm_deriv_potential(chi,alpha)

       case ('smallf')
          deriv_matter_potential(1) = sfi_norm_deriv_potential(chi/mu,p)/mu

       case ('gswli')
          deriv_matter_potential(1) = li_norm_deriv_potential(chi,alpha)

       case ('interm')
          deriv_matter_potential(1) = ii_norm_deriv_potential(chi,beta)

       case ('kahmo2')
          deriv_matter_potential(1) = kmiii_norm_deriv_potential(chi,alpha,beta)

       case ('colwei')
          deriv_matter_potential(1) = cwi_norm_deriv_potential(chi/mu,alpha,mu)/mu

       case ('higgsi')
          deriv_matter_potential(1) = hi_norm_deriv_potential(chi)

       case ('twisti')
          deriv_matter_potential(1) = twi_norm_deriv_potential(chi,mu)

       case ('tdwell')
          deriv_matter_potential(1) = dwi_norm_deriv_potential(chi/mu)/mu

       case ('mhitop')
          deriv_matter_potential(1) = mhi_norm_deriv_potential(chi/mu)/mu

       case ('logmdi')
          deriv_matter_potential(1) = lmi_norm_deriv_potential(chi,gam,beta)

       case ('ricci1')
          deriv_matter_potential(1) = rpi1_norm_deriv_potential(chi/q,p)/q

       case ('ricci2')
          deriv_matter_potential(1) = rpi2_norm_deriv_potential(chi/q,p)/q

       case ('betexp')
          deriv_matter_potential(1) = bei_norm_deriv_potential(chi,lambda,beta)

       case ('radiag')
          deriv_matter_potential(1) = rgi_norm_deriv_potential(chi,alpha)

       case ('nmssmi')
          deriv_matter_potential(1) = mssmi_norm_deriv_potential(chi/mu,mu)/mu

       case ('rinfpt')
          deriv_matter_potential(1) = ripi_norm_deriv_potential(chi,alpha)

       case ('gmssmi')
          deriv_matter_potential(1) = gmssmi_norm_deriv_potential(chi/mu,1.5_kp*alpha,mu)/mu

       case ('bsusyb')
          deriv_matter_potential(1) = bsusybi_norm_deriv_potential(chi,gam)

       case ('tipinf')
          deriv_matter_potential(1) = ti_norm_deriv_potential(chi/mu,alpha)/mu

       case ('psenat')
          deriv_matter_potential(1) = psni_norm_deriv_potential(chi/mu,alpha)/mu

       case ('nckahi')
          deriv_matter_potential(1) = ncki_norm_deriv_potential(chi,alpha,beta)

       case ('hybrid')
          deriv_matter_potential(1) = vhi_norm_deriv_potential(chi/mu,p)/mu

       case ('dysusy')
          deriv_matter_potential(1) = dsi_norm_deriv_potential(chi/mu,p)/mu

       case ('arctan')
          deriv_matter_potential(1) = ai_norm_deriv_potential(chi/mu)/mu

       case ('fixnsa')
          deriv_matter_potential(1) = cnai_norm_deriv_potential(chi,alpha)

       case ('fixnsb')
          deriv_matter_potential(1) = cnbi_norm_deriv_potential(chi,alpha)

       case ('fixnsc')
          deriv_matter_potential(1) = cnci_norm_deriv_potential(chi,alpha)

       case ('fixnsd')
          deriv_matter_potential(1) = cndi_norm_deriv_potential(chi,alpha,beta)

       case ('nszero')
          deriv_matter_potential(1) = csi_norm_deriv_potential(chi,alpha)

       case ('oifold')
          deriv_matter_potential(1) = oi_norm_deriv_potential(chi/mu,alpha)/mu
          
       case ('sugrab')
          deriv_matter_potential(1) = sbi_norm_deriv_potential(chi,alpha,beta)

       case ('ssbinf')
          deriv_matter_potential(1) = ssbi_norm_deriv_potential(chi,alpha,beta)

       case ('runmas')
          deriv_matter_potential(1) = rmi_norm_deriv_potential(chi/mu,2._kp*nu,mu)/mu

       case ('logpot')
          deriv_matter_potential(1) = lpi_norm_deriv_potential(chi/mu,p,q)/mu

!       case ('f-term')
!          deriv_matter_potential(1) = lambda * (16._kp * lambda / M**4 * log(2._kp) &
!               / (4._kp * Pi**2 * phic)+ 4._kp * phi *psi**2 / M**2 / phic**2 )
!          deriv_matter_potential(2) = lambda*(- 4._kp * (1._kp - psi**2 / M**2 ) &
!               * psi / M**2 +4._kp * psi * phi**2 / M**2 / phic**2 ) 

       case default
          write(*,*)'name is ',potName
          stop 'deriv_matter_potential: model not found!'

    end select

    deriv_matter_potential = M4 * deriv_matter_potential
    
#else
    stop 'deriv_matter_potential: FieldInf not built with libaspic support!'
#endif

  end function deriv_matter_potential
  


  function deriv_second_matter_potential(matter)
    implicit none
    real(kp), dimension(matterNum,matterNum) :: deriv_second_matter_potential
    real(kp), dimension(matterNum) :: matter

    real(kp) :: chi, lnchi, expchi, expchi2, coschi

!    real(kp) :: phi, psi

    chi = matter(1)
    
!    phi = matter(1)
!    psi = matter(2)


#ifndef PPNAME
     
    if (chi.gt.0._kp) then
       lnchi = log(chi)
    else
       lnchi = 0._kp
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

    expchi2 = exp(chi**potParam(18)*potParam(17))

    deriv_second_matter_potential(1,1) = deriv_second_matter_potential(1,1) &
         + potParam(13)*potParam(14)*(potParam(14)-1._kp)*chi**(potParam(14)-2._kp) &
         + chi**(-2._kp+potParam(16))*expchi2*potParam(15)*(-1._kp+potParam(16)) &
         * potParam(16)+2._kp*chi**(-2._kp+potParam(16)+potParam(18))*expchi2 &
         * potParam(15)*potParam(16)*potParam(17)*potParam(18) &
         +chi**potParam(16)*potParam(15)*(chi**(-2._kp+potParam(18))*expchi2 &
         *potParam(17)*(-1._kp+potParam(18))*potParam(18) &
         +chi**(-2._kp+2._kp*potParam(18))*expchi2*potParam(17)**2*potParam(18)**2)    
#endif
#endif

    

#elif !defined(NOASPIC)

    select case (potName)

       case ('largef')
          deriv_second_matter_potential(1,1) = lfi_norm_deriv_second_potential(chi,p)
          
       case ('gmixlf')
          deriv_second_matter_potential(1,1) = gmlfi_norm_deriv_second_potential(chi,p,q,alpha)
          
       case ('rcmass')
          deriv_second_matter_potential(1,1) = rcmi_norm_deriv_second_potential(chi,alpha)

       case ('rcquad')
          deriv_second_matter_potential(1,1) = rcqi_norm_deriv_second_potential(chi,alpha)

       case ('natinf')
          deriv_second_matter_potential(1,1) = ni_norm_deriv_second_potential(chi,mu)

       case ('exsusy')
          deriv_second_matter_potential(1,1) = esi_norm_deriv_second_potential(chi,q)

       case ('powlaw')
          deriv_second_matter_potential(1,1) = pli_norm_deriv_second_potential(chi,alpha)

       case ('kahmo1')
          deriv_second_matter_potential(1,1) = kmii_norm_deriv_second_potential(chi,alpha)

       case ('hfline')
          deriv_second_matter_potential(1,1) = hf1i_norm_deriv_second_potential(chi,alpha)

       case ('smallf')
          deriv_second_matter_potential(1,1) = sfi_norm_deriv_second_potential(chi/mu,p)/mu/mu

       case ('gswli')
          deriv_second_matter_potential(1,1) = li_norm_deriv_second_potential(chi,alpha)

       case ('interm')
          deriv_second_matter_potential(1,1) = ii_norm_deriv_second_potential(chi,beta)

       case ('kahmo2')
          deriv_second_matter_potential(1,1) = kmiii_norm_deriv_second_potential(chi,alpha,beta)

       case ('colwei')
          deriv_second_matter_potential(1,1) = cwi_norm_deriv_second_potential(chi/mu,alpha,mu)/mu/mu

       case ('higgsi')
          deriv_second_matter_potential(1,1) = hi_norm_deriv_second_potential(chi)

       case ('twisti')
          deriv_second_matter_potential(1,1) = twi_norm_deriv_second_potential(chi,mu)

       case ('tdwell')
          deriv_second_matter_potential(1,1) = dwi_norm_deriv_second_potential(chi/mu)/mu/mu

       case ('mhitop')
          deriv_second_matter_potential(1,1) = mhi_norm_deriv_second_potential(chi/mu)/mu/mu

       case ('logmdi')
          deriv_second_matter_potential(1,1) = lmi_norm_deriv_second_potential(chi,gam,beta)

       case ('ricci1')
          deriv_second_matter_potential(1,1) = rpi1_norm_deriv_second_potential(chi/q,p)/q/q

       case ('ricci2')
          deriv_second_matter_potential(1,1) = rpi2_norm_deriv_second_potential(chi/q,p)/q/q

       case ('betexp')
          deriv_second_matter_potential(1,1) = bei_norm_deriv_second_potential(chi,lambda,beta)

       case ('radiag')
          deriv_second_matter_potential(1,1) = rgi_norm_deriv_second_potential(chi,alpha)

       case ('nmssmi')
          deriv_second_matter_potential(1,1) = mssmi_norm_deriv_second_potential(chi/mu,mu)/mu/mu

       case ('rinfpt')
          deriv_second_matter_potential(1,1) = ripi_norm_deriv_second_potential(chi,alpha)

       case ('gmssmi')
          deriv_second_matter_potential(1,1) = gmssmi_norm_deriv_second_potential(chi/mu &
               ,1.5_kp*alpha,mu)/mu/mu

       case ('bsusyb')
          deriv_second_matter_potential(1,1) = bsusybi_norm_deriv_second_potential(chi,gam)

       case ('tipinf')
          deriv_second_matter_potential(1,1) = ti_norm_deriv_second_potential(chi/mu,alpha)/mu/mu

       case ('psenat')
          deriv_second_matter_potential(1,1) = psni_norm_deriv_second_potential(chi/mu &
               ,alpha)/mu/mu

       case ('nckahi')
          deriv_second_matter_potential(1,1) = ncki_norm_deriv_second_potential(chi,alpha,beta)

       case ('hybrid')
          deriv_second_matter_potential(1,1) = vhi_norm_deriv_second_potential(chi/mu,p)/mu/mu

       case ('dysusy')
          deriv_second_matter_potential(1,1) = dsi_norm_deriv_second_potential(chi/mu,p)/mu/mu

       case ('arctan')
          deriv_second_matter_potential(1,1) = ai_norm_deriv_second_potential(chi/mu)/mu/mu

       case ('fixnsa')
          deriv_second_matter_potential(1,1) = cnai_norm_deriv_second_potential(chi,alpha)

       case ('fixnsb')
          deriv_second_matter_potential(1,1) = cnbi_norm_deriv_second_potential(chi,alpha)

       case ('fixnsc')
          deriv_second_matter_potential(1,1) = cnci_norm_deriv_second_potential(chi,alpha)

       case ('fixnsd')
          deriv_second_matter_potential(1,1) = cndi_norm_deriv_second_potential(chi,alpha,beta)

       case ('nszero')
          deriv_second_matter_potential(1,1) = csi_norm_deriv_second_potential(chi,alpha)

       case ('oifold')
          deriv_second_matter_potential(1,1) = oi_norm_deriv_second_potential(chi/mu &
               ,alpha)/mu/mu

       case ('sugrab')
          deriv_second_matter_potential(1,1) = sbi_norm_deriv_second_potential(chi,alpha,beta)

       case ('ssbinf')
          deriv_second_matter_potential(1,1) = ssbi_norm_deriv_second_potential(chi,alpha,beta)

       case ('runmas')
          deriv_second_matter_potential(1,1) = rmi_norm_deriv_second_potential(chi/mu &
               ,2._kp*nu,mu)/mu/mu

       case ('logpot')
          deriv_second_matter_potential(1,1) = lpi_norm_deriv_second_potential(chi/mu,p,q)/mu/mu

!       case ('f-term')
!          deriv_second_matter_potential(1,1) = lambda * ( 4.  *psi**2 / M**2 / phic**2 ) 
!          deriv_second_matter_potential(1,2) = lambda * 8. * phi * psi / M**2 / phic**2 
!          deriv_second_matter_potential(2,1) = lambda * 8. * phi * psi / M**2 / phic**2 
!          deriv_second_matter_potential(2,2) = lambda*(- 4. * (1. - psi**2 / M**2 )  &
!               / M**2  + 8. * psi**2 / M**4  + 4. * phi**2 / M**2 / phic**2 ) 


       case default
          write(*,*)'name is ',potName
          stop 'deriv_second_matter_potential: model not found!'

    end select

    deriv_second_matter_potential = M4 * deriv_second_matter_potential   
           
#else
    stop 'deriv_second_matter_potential: FieldInf not built agains libsrmodels!'
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
