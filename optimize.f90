!This will eventually be a subroutine to call the winding surface optimization
module optimize
  use stel_kinds
  integer :: npop = 100
  real(dp), dimension(:,:), allocatable :: xarr !array to pass into swarm
  real(dp), dimension(:), allocatable :: ub, lb !starting bounds

contains

subroutine run_optimize

  use global_variables
  use init_surface_mod
  use pswarm

  implicit none
  integer :: i,k
  real(dp) :: r0, dtheta, dzeta, f
  real(dp) :: d2rdtheta2, d2rdthetadzeta, d2rdzeta2 !dummy variables
  real(dp), dimension(ntotal_ws*2) :: x


  !The sixth value is the major radius. We'll vary this a little to play with it
  !r0 = rmnc_ws(6)
  !print *,'original r0',r0
  !call auto_regularization_solve()
  allocate(ub(ntotal_ws*2))
  allocate(lb(ntotal_ws*2))

  x(1:ntotal_ws) = rmnc_ws
  x(ntotal_ws+1:ntotal_ws*2) = zmns_ws
 
  !make initial bounds guesses
  do i = 1,ntotal_ws*2
     if (x(i) > 0) then
        ub(i) = x(i)*2
        lb(i) = x(i)/2
     else
        ub(i) = x(i)/2
        lb(i) = x(i)*2
     end if
     !print *,lb(i),x(i),ub(i)
  end do
  call swarm_optimize(ntotal_ws, 10, lb, ub, 2, get_lambda)
  ! do i = 0,1
  !    x(6) = r0 + i/10.
  !    call get_lambda(x, f)
  !    print *,'calc of f for r',x(6),f
  ! end do
    
end subroutine run_optimize

subroutine get_lambda(x, f)
  use global_variables
  use init_surface_mod

  implicit none
  real(dp) :: r0, dtheta, dzeta, f
  real(dp), dimension(ntotal_ws*2) :: x
  real(dp) :: d2rdtheta2, d2rdthetadzeta, d2rdzeta2 !dummy variables

  
  rmnc_ws = x(1:ntotal_ws)
  zmns_ws = x(ntotal_ws+1:ntotal_ws*2)

  call calc_nescin_vars(r_coil, ntotal_ws, xm_ws, xn_ws, rmnc_ws, zmns_ws, rmns_ws, zmnc_ws, drdtheta_coil, drdzeta_coil, d2rdtheta2, d2rdthetadzeta, d2rdzeta2, ntheta_coil, nzetal_coil, theta_coil, zetal_coil, .False.)
     
  !print *,'calculating normals'
  dtheta = theta_coil(2) - theta_coil(1)
  dzeta = zeta_coil(2) - zeta_coil(1)
  call calc_normals(ntheta_coil, nzeta_coil, nzetal_coil, normal_coil, norm_normal_coil, drdtheta_coil, drdzeta_coil)
  call calc_volume(ntheta_coil, nzetal_coil, r_coil, dtheta, dzeta, sum(norm_normal_coil), area_coil)
  !print *,'rebuilding matrices'
  call build_matrices(.False.)

  call auto_regularization_solve()
  if (exit_code == 0) then 
     f = lambda(nlambda)
  else
     f = 1.0
  end if
  print *,'final f',f
  
end subroutine get_lambda
end module optimize
