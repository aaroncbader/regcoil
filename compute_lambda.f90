subroutine compute_lambda

  use global_variables, only: nlambda, lambda_min, lambda_max, lambda
  use stel_kinds

  implicit none

  integer :: j

  allocate(lambda(nlambda))
  
  lambda(1) = 0
  do j = 1,nlambda-1
     lambda(j+1) = lambda_min * exp((log(lambda_max/lambda_min)*(j-1))/(nlambda-2))
  end do

  print *,"We will use the following values of the regularization weight lambda:"
  print *,lambda

end subroutine compute_lambda


