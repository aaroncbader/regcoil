!This will eventually be a subroutine to call the winding surface optimization

subroutine optimize

  use global_variables
  use init_surface_mod
  use stel_kinds
  implicit none
  integer :: i
  real(dp) :: r0, dtheta, dzeta
  real(dp) :: d2rdtheta2, d2rdthetadzeta, d2rdzeta2 !dummy variables
  print *,'mnmax_coil',mnmax_coil
  !The sixth value is the major radius. We'll vary this a little to play with it
  
  r0 = rmnc_ws(6)
  do i = 0,1
     rmnc_ws(6) = r0 + i/10.
     print *,'calculating nescin vars'
     call calc_nescin_vars(r_coil, ntotal_ws, xm_ws, xn_ws, rmnc_ws, zmns_ws, rmns_ws, zmnc_ws, drdtheta_coil, drdzeta_coil, d2rdtheta2, d2rdthetadzeta, d2rdzeta2, ntheta_coil, nzetal_coil, theta_coil, zetal_coil, .False.)
     print *,'calculating normals'
     dtheta = theta_coil(2) - theta_coil(1)
     dzeta = zeta_coil(2) - zeta_coil(1)
     call calc_normals(ntheta_coil, nzeta_coil, nzetal_coil, normal_coil, norm_normal_coil, drdtheta_coil, drdzeta_coil)
     call calc_volume(ntheta_coil, nzetal_coil, r_coil, dtheta, dzeta, sum(norm_normal_coil), area_coil)
     print *,'rebuilding matrices'
     call build_matrices(.False.)

     call auto_regularization_solve()

     !exit_code 0 means that auto_regularization_solve completed as normal
     !exit_code < 0 means it failed
     !if it exits properly, the calculated lambda is lambda(nlambda)
     print *,'final lambda',nlambda,lambda(nlambda)
     print *,'exit code',exit_code
  end do
  
  
  
  
end subroutine optimize
  
