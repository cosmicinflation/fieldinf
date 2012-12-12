module infsigma
  use infprec, only : kp
  use infdilaton, only : conformal_factor_square, conformal_first_gradient
  use infdilaton, only : conformal_second_gradient
  use infbgmodel, only : fieldNum, matterNum, dilatonNum

  implicit none

! This module gives the metric on the sigma model manifold and all its derivatives
!
! S = 1/(2kappa^2) {{ int{sqrt(-g) d^4 x} [ R - metric Dfield Dfield - 2*potential(Field)] }}  
!

  private
  

  public metric, metric_inverse, deriv_metric
  public connection_affine, deriv_connection_affine


contains

 
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

    metricVal = metric(field)

    metric_inverse = 0._kp
               
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

    if (dilatonNum.eq.0) return

    dilaton = field(matterNum+1:fieldNum)
    confSquare = conformal_factor_square(dilaton)
    confFirstGrad = conformal_first_gradient(dilaton)
    
    do i=1,matterNum
       deriv_metric(i,i,matterNum+1:fieldNum) &
            = 2._kp * confFirstGrad(1:dilatonNum)* confSquare
    enddo
    
  end function deriv_metric



  function connection_affine(field)
!first index contravariant, 2 others covariant and symmetric 
    implicit none
    real(kp), dimension(fieldNum) :: field    
    real(kp), dimension(fieldNum,fieldNum) :: metricInv
    real(kp), dimension(dilatonNum,dilatonNum) :: metricInvConf
    real(kp), dimension(fieldNum,fieldNum,fieldNum) :: connection_affine
    real(kp) ::  confSquare
    real(kp), dimension(dilatonNum) :: dilaton
    real(kp), dimension(dilatonNum) :: confFirstGradVec, confFirstGrad
   
    integer :: i

    connection_affine = 0._kp

    if (dilatonNum.eq.0) return

    metricInv = metric_inverse(field)
    metricInvConf(1:dilatonNum,1:dilatonNum) &
         = metricInv(matterNum+1:fieldNum,matterNum+1:fieldNum)

    dilaton=field(matterNum+1:fieldNum)
    confFirstGrad = conformal_first_gradient(dilaton)
    confFirstGradVec = matmul(metricInvConf,confFirstGrad)
    confSquare = conformal_factor_square(dilaton)
    
    do i=1,dilatonNum
       connection_affine(1:matterNum,1:matterNum,matterNum+i) &
            = confFirstGrad(i)
       connection_affine(1:matterNum,matterNum+i,1:matterNum) &
            = confFirstGrad(i)       
    enddo

    do i=1,matterNum
       connection_affine(matterNum+1:fieldNum,i,i) &
            = - confSquare * confFirstGradVec(1:dilatonNum)
    enddo

   
  end function connection_affine




  function deriv_connection_affine(field)
!first partial derivative of the christoffel: first index contrariant,
!3 others covariant   
    implicit none

    real(kp), dimension(fieldNum) :: field
    real(kp), dimension(fieldNum,fieldNum) :: metricInv
    real (kp), dimension(fieldNum,fieldNum,fieldNum) :: metricDeriv
    real(kp), dimension(fieldNum,fieldNum,fieldNum,fieldNum) :: deriv_connection_affine    

    real(kp) :: confSquare
    real(kp), dimension(dilatonNum) :: dilaton, confFirstGrad, confFirstGradVec
    real(kp), dimension(dilatonNum,dilatonNum) :: metricInvConf
    real(kp), dimension(dilatonNum,dilatonNum) :: confSecondGrad, confSecondGradVec
    integer :: i,j,k

    deriv_connection_affine = 0._kp

    if (dilatonNum.eq.0) return

    dilaton(1:dilatonNum) = field(matterNum+1:fieldNum)
    confSquare = conformal_factor_square(dilaton)

!    metricDeriv = deriv_metric(field)
    metricInv = metric_inverse(field)
    metricInvConf = metricInv(matterNum+1:fieldNum,matterNum+1:fieldNum)
    
    confFirstGrad = conformal_first_gradient(dilaton)
    confFirstGradVec = matmul(metricInvConf,confFirstGrad)
    confSecondGrad = conformal_second_gradient(dilaton)
    confSecondGradVec = matmul(metricInvConf,confSecondGrad)
    

    do j=1,dilatonNum
       do i=1,dilatonNum

          deriv_connection_affine(1:matterNum,1:matterNum,matterNum+i,matterNum+j) &
               = confSecondGrad(i,j)
          deriv_connection_affine(1:matterNum,matterNum+i,1:matterNum,matterNum+j) & 
               = confSecondGrad(i,j)

          do k=1,matterNum
             deriv_connection_affine(matterNum+i,k,k,matterNum+j) &
                  = - confSquare * (confSecondGradVec(i,j) &
                  + 2._kp * confFirstGradVec(i) * confFirstGrad(j)) 

!zero when the the metric for the dilaton is constant and diagonal
!                  + confSquare * dot_product(metricInvConf(i,:) &
!                  ,matmul(metricDeriv(matterNum+1:fieldNum,matterNum+1:fieldNum,j) &
!                  ,confFirstGradVec))
          enddo

       enddo
    enddo


  end function deriv_connection_affine



end module infsigma
