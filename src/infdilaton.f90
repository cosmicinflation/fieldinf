module infdilaton
  use fieldprec, only : kp

  implicit none

  private

!number of scalar gravity fields
  integer, parameter :: dilatonNum = 0


!number if conformal related parameters
  integer, parameter :: confParamNum = 0



  public dilatonNum, confParamNum

  public conformal_factor_square, conformal_first_gradient, conformal_second_gradient


contains

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




end module infdilaton
