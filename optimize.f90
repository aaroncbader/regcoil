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

subroutine dist_squared(x1,y1,z1,x2,y2,z2,dist)
  implicit none
  real(dp) :: x1, y1, z1, x2, y2, z2, dist

  dist = (x1-x2)**2 + (y1 - y2)**2 + (z1 - z2)**2
end subroutine dist_squared

!This calculates the minimum distance from a surface
!The method is as follows. Start with corresponding indices for the plasma
!and the coil. This is the initial guess as to the closest corresponding points
!Then search the four neighboring points in theta,zeta indices. If any are
!closer, restart the search from the closest point.
subroutine min_surf_distance(mindist)
  use global_variables

  implicit none

  integer :: i,j, is, js, iu, id, ju, jd, minind
  real(dp), dimension(4) :: d
  real(dp) :: mindist
  real(dp) :: xp0, yp0, zp0, d0, dn

  if ((ntheta_coil .ne. ntheta_plasma) .or. (nzeta_coil .ne. nzeta_plasma)) then
     print *,'dimensions of plasma and coil do not match, fix in input'
     mindist = 0.0
     return
  end if

  !set an arbitrarily high initial distance
  mindist = 1.0E20
  do i = 1,ntheta_coil
     do j = 1,nzeta_coil
        xp0 = r_plasma(1,i,j)
        yp0 = r_plasma(2,i,j)
        zp0 = r_plasma(3,i,j)
        call dist_squared(xp0, yp0, zp0, r_coil(1,i,j), r_coil(2,i,j), &
             r_coil(3,i,j), d0)
        dn = 0
        is = i
        js = j
        do while (dn < d0)
           !first and last poloidal indices are the same, so we skip them on
           !wraparound
           iu = is+1
           if (iu > ntheta_coil) iu = 2
           id = is-1
           if (id < 1) id = ntheta_coil-1
           !wrapping around toroidally is much harder, we don't bother
           ju = js+1
           if (ju > nzeta_coil) ju = nzeta_coil
           jd = js-1
           if (jd < 1) jd = 1
           ! calculate the four distances
           call dist_squared(xp0, yp0, zp0, &
                r_coil(1,iu,js), r_coil(2,iu,js), r_coil(3,iu,js), d(1))
           call dist_squared(xp0, yp0, zp0, &
                r_coil(1,id,js), r_coil(2,id,js), r_coil(3,id,js), d(2))
           call dist_squared(xp0, yp0, zp0, &
                r_coil(1,is,ju), r_coil(2,is,ju), r_coil(3,is,ju), d(3))
           call dist_squared(xp0, yp0, zp0, &
                r_coil(1,is,jd), r_coil(2,is,jd), r_coil(3,is,jd), d(4))
           !If the smallest distance of these four points is less than
           !our original guess we're done, exit the loop
           dn = minval(d)
           if (dn > d0) exit
           d0 = dn
           !we found a smaller distance, recenter around this point
           minind = minloc(d,1)
           select case (minind)
              case (1)
                 is = iu
              case (2)
                 is = id
              case (3)
                 js = ju
              case (4)
                 js = jd
            end select
         end do
         !did we find a new minimum?
         if (d0 < mindist) mindist = d0
      end do
   end do
   mindist = sqrt(mindist)   
   
end subroutine min_surf_distance

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
  
  call min_surf_distance(mindist)

  call get_inside_poly(is_inside)
  

  if (exit_code == 0) then 
     !f = -1*((volume_coil)**(1.0/3.0) + log(lambda(nlambda_autoreg)))
     f = opt_chi2b_par * chi2_B(nlambda_autoreg) - &
         opt_vol_par * (volume_coil)**(1.0/3.0) - &
         opt_lambda_par * log(lambda(nlambda_autoreg))
  else
     f = 100000.0
  end if
  if (mindist < opt_min_dist) then
     f = f + (opt_min_dist - mindist)*opt_min_par
  end if
  if (is_inside == 0) then
     f = f+100000
  end if

end subroutine get_lambda
end module optimize
