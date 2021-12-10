!   This file is part of FieldInf
!
!   Copyright (C) 2005-2021 C. Ringeval
!   
!   FieldInf is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   FieldInf is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with FieldInf.  If not, see <https://www.gnu.org/licenses/>.





module infpotential
  use fieldprec, only : kp  
  use infmatter, only : matter_potential, deriv_matter_potential
  use infmatter, only : deriv_ln_matter_potential, deriv_second_matter_potential
  use infdilaton, only : conformal_factor_square, conformal_first_gradient
  use infdilaton, only : conformal_second_gradient
  use infbgmodel, only : matterNum, fieldNum, dilatonNum

! This module gives the full potential and all its derivatives
!
! S = 1/(2kappa^2) {{ int{sqrt(-g) d^4 x} [ R - metric Dfield Dfield - 2*potential(Field)] }}  
!
    
  implicit none

  private


  public potential, deriv_potential, deriv_second_potential, deriv_ln_potential


contains

  
  function potential(field)
    implicit none
    real(kp), dimension(fieldNum), intent(in) :: field

    real(kp), dimension(matterNum) :: matter
    real(kp), dimension(dilatonNum) :: dilaton
    real(kp) :: potential  

    matter = field(1:matterNum)
    dilaton = field(matterNum+1:fieldNum)

    potential = matter_potential(matter)*(conformal_factor_square(dilaton))**2

  end function potential




  function deriv_potential(field)
!with respect to the fields   
    implicit none
    real(kp), dimension(fieldNum) :: deriv_potential
    real(kp), dimension(fieldNum), intent(in) :: field

    real(kp), dimension(matterNum) :: matter, derivMatterPot
    real(kp), dimension(dilatonNum) :: dilaton, confFirstGrad

    real(kp) :: potentialVal
    

    matter = field(1:matterNum)
    dilaton = field(matterNum+1:fieldNum)

    potentialVal = potential(field)
    derivMatterPot = deriv_matter_potential(matter)
    

    deriv_potential(1:matterNum) &
         =  derivMatterPot(1:matterNum)*(conformal_factor_square(dilaton))**2

    if (dilatonNum.eq.0) return

    confFirstGrad = conformal_first_gradient(dilaton)

    deriv_potential(matterNum+1:fieldNum) &
         = 4._kp*confFirstGrad(1:dilatonNum)*potentialVal    
       
  end function deriv_potential




  function deriv_second_potential(field)   
    implicit none
    real(kp), dimension(fieldNum,fieldNum) :: deriv_second_potential
    real(kp), dimension(fieldNum), intent(in) :: field

    real(kp), dimension(fieldNum) :: derivPot
    real(kp), dimension(matterNum) :: matter, derivMatterPot
    real(kp), dimension(matterNum,matterNum) :: derivSecondMatterPot
    real(kp), dimension(dilatonNum) :: dilaton, confFirstGrad
    real(kp), dimension(dilatonNum,dilatonNum) :: confSecondGrad

    real(kp) :: potentialVal
    integer :: i,j

    matter = field(1:matterNum)
    dilaton = field(matterNum+1:fieldNum)

    potentialVal = potential(field)

    derivPot = deriv_potential(field)
    derivMatterPot = deriv_matter_potential(matter)
    derivSecondMatterPot = deriv_second_matter_potential(matter)

    
    deriv_second_potential(1:matterNum,1:matterNum) &
         = derivSecondMatterPot(1:matterNum,1:matterNum) &
         *(conformal_factor_square(dilaton))**2

    if (dilatonNum.eq.0) return

    confFirstGrad = conformal_first_gradient(dilaton)
    confSecondGrad = conformal_second_gradient(dilaton)

    do i=1,matterNum
       deriv_second_potential(i,matterNum+1:fieldNum) &
            = derivPot(i) * 4._kp*confFirstGrad(1:dilatonNum)
       deriv_second_potential(matterNum+1:fieldNum,i) &
            = deriv_second_potential(i,matterNum+1:fieldNum)
    enddo
   
    do i=1,dilatonNum
       do j=1,dilatonNum
          deriv_second_potential(matterNum+i,matterNum+j) = potentialVal &
               * (4._kp*confSecondGrad(i,j) + 4._kp*confFirstGrad(i) &
               * 4._kp*confFirstGrad(j))
       enddo
    enddo

  end function deriv_second_potential




  function deriv_ln_potential(field)
    implicit none
    real(kp), dimension(fieldNum) :: deriv_ln_potential
    real(kp), dimension(fieldNum), intent(in) :: field

    real(kp), dimension(matterNum) :: matter,derivLnMatterPot
    real(kp), dimension(dilatonNum) :: dilaton,confFirstGrad

    matter = field(1:matterNum)
    dilaton = field(matterNum+1:fieldNum)

    derivLnMatterPot = deriv_ln_matter_potential(matter)
    
    deriv_ln_potential(1:matterNum) =  derivLnMatterPot(1:matterNum)


    if (dilatonNum.eq.0) return
    confFirstGrad = conformal_first_gradient(dilaton)

    deriv_ln_potential(matterNum+1:fieldNum) = 4._kp*confFirstGrad(1:dilatonNum)
    
  end function deriv_ln_potential



end module infpotential
