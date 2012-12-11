!provides theoretical bounds on the model parameters

module infbounds
  use infprec, only : kp  
  implicit none

  private

  logical, parameter :: display = .true.
  
  public field_stopinf, field_thbound
    
contains

  

  function field_stopinf(infParam)
!return the field stop value in (1) and if it is a maximum (+1) or
!minimum (-1) in (2)
    use infbgmodel, only : infbgparam, matterParamNum
    implicit none
    type(infbgparam), intent(in) :: infParam    
    real(kp), dimension(2) :: field_stopinf

    
    select case (infParam%name)

    case ('largef')
       field_stopinf(1) = infParam%consts(matterParamNum)
       field_stopinf(2) = -1._kp
      
    case ('hybrid')

       field_stopinf(1) = infParam%consts(matterParamNum)
       field_stopinf(2) = -1._kp
       
    case ('runmas')

       field_stopinf(1) = infParam%consts(matterParamNum)
!for nu>0
       if (infParam%consts(matterParamNum).lt.infParam%consts(3)) then
          field_stopinf(2) = -1._kp
       else
          field_stopinf(2) = +1._kp
       endif
!reverse if nu<0
       if (infParam%consts(4).lt.0._kp) then
          field_stopinf(2) = -field_stopinf(2)
       endif

    case ('kklmmt')
       field_stopinf(1) = infParam%consts(matterParamNum)*infParam%consts(3)
       field_stopinf(2) = -1._kp
       print *, 'fieldStop = infParam(matterParamNum)'

    case ('mixlf')
       field_stopinf(1) = infParam%consts(matterParamNum)
       field_stopinf(2) = -1._kp


    case ('powlaw')
       field_stopinf(1) = infParam%consts(matterParamNum)
       field_stopinf(2) = +1._kp


    case ('interm')
       field_stopinf(1) = infParam%consts(matterParamNum)
       field_stopinf(2) = +1._kp


    case ('twisti')
       field_stopinf(1) = infParam%consts(matterParamNum)
       field_stopinf(2) = -1._kp
    

    case ('logmd2')
       field_stopinf(1) = infParam%consts(matterParamNum)
       field_stopinf(2) = +1._kp


    case ('ricci2')
       field_stopinf(1) = infParam%consts(matterParamNum)
       field_stopinf(2) = +1._kp


    case ('bsusyb')
       field_stopinf(1) = infParam%consts(matterParamNum)
       field_stopinf(2) = -1._kp

    case ('dysusy')

       field_stopinf(1) = infParam%consts(matterParamNum)
       field_stopinf(2) = +1._kp

    end select

  end function field_stopinf



  function field_thbound(infParam)
!returns a theoretical bound on the allowed field values if any. May be
!from the stochastic regime or uv limit in brane setup
    use infbgmodel, only : infbgparam, matterParamNum
    implicit none
    type(infbgparam), intent(in) :: infParam    
    real(kp), dimension(2) :: field_thbound

    
    select case (infParam%name)

    case ('kklmmt')
       field_thbound(1) = infParam%consts(matterParamNum-1)
       field_thbound(2) = +1._kp
       print *, 'fieldUv = infParam%consts(matterParamNum-1)'

    case default
       stop 'no theoretical field bound implemented for this model!'

    end select

  end function field_thbound

end module infbounds
