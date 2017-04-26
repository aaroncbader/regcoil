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
  integer :: i,k, nproc
  real(dp) :: r0, dtheta, dzeta, f
  real(dp) :: d2rdtheta2, d2rdthetadzeta, d2rdzeta2 !dummy variables
  real(dp), dimension(ntotal_ws*2) :: x



  allocate(ub(ntotal_ws*2))
  allocate(lb(ntotal_ws*2))

  x(1:ntotal_ws) = rmnc_ws
  x(ntotal_ws+1:ntotal_ws*2) = zmns_ws

  !make initial bounds guesses
  do i = 1,ntotal_ws*2
     if (x(i) > 0) then
        ub(i) = x(i)*1.1
        lb(i) = x(i)/1.1
     else
        ub(i) = x(i)/1.1
        lb(i) = x(i)*1.1
     end if
  end do
  
  call swarm_optimize(ntotal_ws*2, opt_npop, lb, ub, opt_niter, get_lambda, &
       x, num_procs)
  if (my_pn == 0) then
     print *,'best surface'
     print *,'m, n, rmnc, zmns'
     do i = 1,ntotal_ws
        print *,xm_ws(i), xn_ws(i), x(i), x(ntotal_ws+i)
     end do
  end if
     
    
end subroutine run_optimize

subroutine faux_distance(mindist)
  use global_variables
  
  implicit none
  real(dp), dimension(ntheta_coil, nzeta_coil) :: dist
  real(dp) :: mindist

  if ((ntheta_coil .ne. ntheta_plasma) .or. (nzeta_coil .ne. nzeta_plasma)) then
     print *,'dimensions of plasma and coil do not match, fix in input'
     mindist = 0.0
     return
  end if
  
  dist = (r_plasma(1,:,:) - r_coil(1,:,:))**2 + &
         (r_plasma(2,:,:) - r_coil(2,:,:))**2 + &
         (r_plasma(3,:,:) - r_coil(3,:,:))**2
  mindist = sqrt(minval(minval(dist,1)))
  return
end subroutine faux_distance
  

subroutine get_lambda(x, f)
  use global_variables
  use init_surface_mod

  implicit none
  integer :: i
  real(dp) :: r0, dtheta, dzeta, f
  real(dp), dimension(ntotal_ws*2) :: x
  real(dp) :: d2rdtheta2, d2rdthetadzeta, d2rdzeta2 !dummy variables
  real(dp) :: mindist
  
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
  
  call faux_distance(mindist)
  

  if (exit_code == 0) then 
     !f = -1*((volume_coil)**(1.0/3.0) + log(lambda(nlambda_autoreg)))
     f = chi2_B(nlambda_autoreg) - (volume_coil)**(1.0/3.0)
  else
     f = 100.0
  end if
  if (mindist < 0.3) then
     f = f+50
  end if


end subroutine get_lambda
end module optimize
