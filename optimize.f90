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
        ub(i) = x(i)*opt_divisor
        lb(i) = x(i)/opt_divisor
     else
        ub(i) = x(i)/opt_divisor
        lb(i) = x(i)*opt_divisor
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

!A subroutine to check whether a point is inside a polygon,
!adapted from http://www.ecse.rpi.edu/~wrf/Research/Short_Notes/pnpoly.html
subroutine inside_polygon(xp, yp, xx, yy, poly_size, is_in) 

  implicit none
  integer :: poly_size, in_polygon
  integer :: i,j
  logical :: mx, my, nx, ny, is_in
  real(dp) :: xp, yp, qq
  real(dp), dimension(poly_size) :: xx, yy, x, y

  x = xx - xp
  y = yy - yp

  in_polygon = -1

  do i=1,poly_size
     j = 1+modulo(i, poly_size)
     mx = x(i).ge.0.0
     nx = x(j).ge.0.0
     my = y(i).ge.0.0
     ny = y(j).ge.0.0
     if(.not.((my.or.ny).and.(mx.or.nx)).or.(mx.and.nx)) then
       cycle
     end if
     
     if((my.and.ny.and.(mx.or.nx).and..not.(mx.and.nx))) then
       in_polygon = -in_polygon
       cycle
     end if
     
     qq = (y(i)*x(j)-x(i)*y(j))/(x(j)-x(i))

     if (qq.lt.0) then
       cycle
     else if (qq == 0) then
        ! This condition means that it landed on the surface exactly
        in_polygon=0
        return
     else
        in_polygon=-in_polygon
     end if
  end do

  if (in_polygon.lt.0) then
    in_polygon = 0
  end if
  if (in_polygon == 0) then
     is_in = .false.
  else
     is_in = .true.
  end if
  
end subroutine inside_polygon

subroutine get_inside_poly(is_inside)
  use global_variables

  implicit none
  real(dp) :: x0, y0
  real(dp), dimension(ntheta_plasma, nzeta_plasma) :: xp, yp
  real(dp), dimension(ntheta_coil, nzeta_coil) :: xc, yc
  integer :: i, is_inside
  logical :: inside_section

  if (nzeta_coil .ne. nzeta_plasma)  then
     print *,'dimensions of plasma and coil zeta do not match, fix in input'
     is_inside = 0
     return
  end if
  
  xp = sqrt(r_plasma(1,:,:)**2 + r_plasma(2,:,:)**2)
  yp = r_plasma(3,:,:)
  xc = sqrt(r_coil(1,:,:)**2 + r_coil(2,:,:)**2)
  yc = r_coil(3,:,:)

  is_inside = 1
  do i = 1,nzeta_coil
     x0 = sum(xp(:,i))/ntheta_plasma
     y0 = sum(yp(:,i))/ntheta_plasma
     call inside_polygon(x0, y0, xc(:,i), yc(:,i), ntheta_coil, inside_section)
     if (.not. inside_section) then
        is_inside = 0
        return
     end if
  end do
     
end subroutine get_inside_poly  


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
  integer :: i, is_inside
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

  call get_inside_poly(is_inside)
  

  if (exit_code == 0) then 
     !f = -1*((volume_coil)**(1.0/3.0) + log(lambda(nlambda_autoreg)))
     f = opt_chi2b_par * chi2_B(nlambda_autoreg) - &
         opt_vol_par * (volume_coil)**(1.0/3.0) - &
         opt_lambda_par * log(lambda(nlambda_autoreg))
  else
     f = 100.0
  end if
  if (mindist < opt_min_dist) then
     f = f+50
  end if
  if (is_inside == 0) then
     f = f+200
  end if

end subroutine get_lambda
end module optimize
