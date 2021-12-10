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






module infbgfunc
  use infbgmodel, only : matterNum, dilatonNum, fieldNum
  use fieldprec, only : kp
  use infpotential, only : potential, deriv_potential
  use infsigma, only : metric, metric_inverse

  implicit none

  private

! Hubble flow and energy density functions (exact)

  public hubble_parameter_square
  public slowroll_first_parameter, slowroll_second_parameter
  public slowroll_first_parameter_JF

  public matter_energy_density, matter_energy_density_JF

contains


  function hubble_parameter_square(field,derivField,useVelocity)
!in unit of kappa^2
    implicit none
    real(kp) :: hubble_parameter_square
    real(kp), dimension(fieldNum), intent(in) :: field, derivfield
    logical, intent(in) :: useVelocity

    real(kp) :: derivFieldSquare

    derivFieldSquare = dot_product(derivField,matmul(metric(field),derivField))

    if (useVelocity) then
       hubble_parameter_square = (derivFieldSquare + 2._kp*potential(field))/6._kp
    else
       hubble_parameter_square = 2._kp*potential(field)/(6._kp - derivFieldSquare)       
    endif
   
  end function hubble_parameter_square




  function slowroll_first_parameter(field,derivField,useVelocity)
!epsilon1 = epsilon    
    implicit none
    real(kp) :: slowroll_first_parameter
    real(kp), dimension(fieldNum), intent(in) :: field, derivField
    logical, intent(in) :: useVelocity
    
    real(kp) :: derivFieldSquare, hubbleSquare

    derivFieldSquare = dot_product(derivField,matmul(metric(field),derivField))
    
    if (useVelocity) then       
       hubbleSquare = hubble_parameter_square(field,derivField,useVelocity)
       slowroll_first_parameter = derivFieldSquare/2d0/hubbleSquare
    else
       slowroll_first_parameter = derivFieldSquare/2d0
    endif
    

  end function slowroll_first_parameter




  function slowroll_second_parameter(field,derivField,useVelocity)
!epsilon2 = 2(epsilon - delta)
    implicit none
    real(kp) :: slowroll_second_parameter
    real(kp), dimension(fieldNum), intent(in) :: field, derivField
    logical, intent(in) :: useVelocity

    real(kp), dimension(fieldNum) :: fieldDot
    real(kp) :: epsilon1
    real(kp) :: hubbleSquare, derivPotFieldDot

    hubbleSquare = hubble_parameter_square(field,derivField,useVelocity)
    
    if (useVelocity) then
       fieldDot = derivField/sqrt(hubbleSquare)
    else
       fieldDot = derivField
    endif

    epsilon1 = slowroll_first_parameter(field,derivField,useVelocity)
    
    if (epsilon1.ne.0.) then
       derivPotFieldDot = dot_product(deriv_potential(field),fieldDot)

       slowroll_second_parameter = -6d0 + 2d0*epsilon1 &
            - derivPotFieldDot/(epsilon1*hubbleSquare)
    else
       slowroll_second_parameter = 0.
    endif

  end function slowroll_second_parameter




  function slowroll_first_parameter_JF(field,derivField,useVelocity)
!this is epsilon1 EF, up to 2% when dilaton couplings = 1   
    use infdilaton, only : conformal_first_gradient, conformal_second_gradient
    use infsigma, only : connection_affine
    implicit none
    real(kp) :: slowroll_first_parameter_JF
    real(kp), dimension(fieldNum), intent(in) :: field, derivField
    logical, intent(in) :: useVelocity

    real(kp), dimension(fieldNum) :: christVec, fieldDot, potDerivVec
    real(kp), dimension(dilatonNum) :: dilaton, dilatonDot,confFirstGrad    
    real(kp), dimension(dilatonNum,dilatonNum) :: confSecondGrad
    real(kp), dimension(fieldNum,fieldNum,fieldNum) :: christoffel
    real(kp) :: derivFieldSquare, hubbleSquare, epsilon1, epsilon1Xfactor, shift
    real(kp) :: confFirstGradXdilDot

    integer :: i
    
    hubbleSquare = hubble_parameter_square(field,derivField,useVelocity)

    if (useVelocity) then
       fieldDot = derivField/sqrt(hubbleSquare)
    else
       fieldDot = derivField
    endif

    dilaton = field(matterNum+1:fieldNum)
    dilatonDot = fieldDot(matterNum+1:fieldNum)
    confFirstGrad = conformal_first_gradient(dilaton)
    confSecondGrad = conformal_second_gradient(dilaton)
!    potDerivVec = deriv_potential_vec(field)
    potDerivVec = matmul(metric_inverse(field),deriv_potential(field))
    christoffel = connection_affine(field)

    do i=1,fieldNum
       christVec(i) = dot_product(fieldDot(:),matmul(christoffel(i,:,:),fieldDot(:)))
    enddo
   
    confFirstGradXdilDot = dot_product(confFirstGrad,dilatonDot)

    epsilon1 = slowroll_first_parameter(field,derivField,useVelocity)

    
    epsilon1Xfactor = (epsilon1 + confFirstGradXdilDot)/(1._kp + confFirstGradXdilDot)
   

    shift = ((3._kp - epsilon1)*confFirstGradXdilDot &
         + dot_product(confFirstGrad,christVec(matterNum+1:fieldNum)) &
         + dot_product(confFirstGrad,potDerivVec(matterNum+1:fieldNum)/hubbleSquare) &
         - dot_product(dilatonDot, matmul(confSecondGrad,dilatonDot))) &
         / (1._kp + confFirstGradXdilDot)**2
         

    slowroll_first_parameter_JF = epsilon1Xfactor + shift


  end function slowroll_first_parameter_JF



  function matter_energy_density(field,velocity)
!energy density of the matter fields in the Einstein Frame
!A^2 (Dchi/Dtphys)^2 + A^4 U   
    use infmatter, only : matter_potential
    use infdilaton, only : conformal_factor_square
    implicit none
    real(kp) :: matter_energy_density
    real(kp), dimension(fieldNum), intent(in) :: field,velocity

    real(kp), dimension(matterNum) :: matter, matterVel
    real(kp), dimension(dilatonNum) :: dilaton
    real(kp) :: confSquare


    matter = field(1:matterNum)
    dilaton = field(matterNum+1:fieldNum)
    matterVel = velocity(1:matterNum)
    confSquare = conformal_factor_square(dilaton)

    matter_energy_density &
         = 0.5d0*confSquare * dot_product(matterVel,matterVel) &
         + matter_potential(matter) * confSquare**2

  end function matter_energy_density




  function matter_energy_density_JF(field,velocity)
!energy density of the matter fields in the Jordan Frame
!EF/A^4    
    use infdilaton, only : conformal_factor_square
    implicit none
    real(kp) :: matter_energy_density_JF
    real(kp), dimension(fieldNum), intent(in) :: field,velocity
    
    real(kp), dimension(dilatonNum) :: dilaton
    real(kp) :: confSquare
    real(kp) :: matterEnergyDensEF

    dilaton = field(matterNum+1:fieldNum)

    confSquare = conformal_factor_square(dilaton)
    matterEnergyDensEF = matter_energy_density(field,velocity)

    matter_energy_density_JF = matterEnergyDensEF/confSquare**2

  end function matter_energy_density_JF



end module infbgfunc
