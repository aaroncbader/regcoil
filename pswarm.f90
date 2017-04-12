module pswarm
  use stel_kinds
  implicit none
  !current position, best position and velocity
  real(dp), dimension(:,:), allocatable :: ps_x, ps_v, ps_bx 
  real(dp), dimension(:), allocatable :: ps_gbx !best known position vector
  real(dp), dimension(:), allocatable :: ps_lb, ps_ub
  real(dp), dimension(:), allocatable :: ps_bf
  real(dp) :: ps_gbf
  integer :: ps_n, ps_npop, ps_gb
  integer :: ps_debug=0

contains 
!Allocate the particle swarm data, these include the actual
!data values for each particle, the global best and the upper and lower
!bound arrays.
subroutine alloc_ps()
  implicit none

  allocate(ps_x(ps_npop, ps_n))
  allocate(ps_v(ps_npop, ps_n))
  allocate(ps_bx(ps_npop, ps_n))
  allocate(ps_bf(ps_npop))
  allocate(ps_gbx(ps_n))
  allocate(ps_lb(ps_n))
  allocate(ps_ub(ps_n))
end subroutine alloc_ps

!Initialize the particles with random values
subroutine init_pop(fcn, xguess)

  use stel_kinds
  implicit none
  integer :: n, i
  real(dp) :: x, diff
  real(dp), dimension(ps_n) :: xguess
  external fcn
  ps_x(1,:) = xguess
  ps_v(1,:) = 0
  call fcn(ps_x(1,:), ps_bf(1))
  ps_gbf = ps_bf(1)
  ps_bx(1,:) = ps_x(1,:)
  ps_gbx = ps_x(1,:)
  
  call random_seed()
  do n = 2,ps_npop
     do i = 1,ps_n
        call random_number(x)
        diff = ps_ub(i) - ps_lb(i)
        ps_x(n,i) = x*diff + ps_lb(i)
        ps_bx(n,i) = ps_x(n,i)
        call random_number(x)
        ps_v(n,i) = x*2*diff - diff
     end do
     !evaluate the function
     call fcn(ps_x(n,:), ps_bf(n))
     !set the particle best to current value

     if (ps_bf(n) < ps_gbf) then
        ps_gbf = ps_bf(n)
        ps_gbx = ps_x(n,:)
     end if
     
  end do
end subroutine init_pop

subroutine iterate(fcn)
  use stel_kinds
  implicit none
  external fcn
  integer :: n,i
  real(dp) :: rp, rg, fcnval
  real(dp) :: omega, phip, phig

  !For now we set the optimization parameters
  omega = 0.5
  phip = 1.0
  phig = 1.0

  ! update the velocity vectors
  do n = 1,ps_npop
     do i = 1,ps_n
        call random_number(rp)
        call random_number(rg)
        !write some diagnostic output
        if (ps_debug == 1) then
           write (*,*) n,i,ps_v(n,i), ps_bx(n,i), ps_x(n,i), ps_gbx(i)
        end if
        ps_v(n,i) = omega*ps_v(n,i) + phip*rp*(ps_bx(n,i) - ps_x(n,i)) + &
             phig*rg*(ps_gbx(i) - ps_x(n,i))
        if (ps_debug == 1) then
           write (*,*) 'new velocity',ps_v(n,i)
        end if
     end do
  end do

  !update the position vectors
  do n = 1,ps_npop
     if (ps_debug == 1) then
        write (*,*) '------particle',n,'-------'
        write (*,*) 'old position',n,ps_x(n,:)
     end if
     ps_x(n,:) = ps_x(n,:) + ps_v(n,:)
     if (ps_debug == 1) then
        write (*,*) 'new position',n,ps_x(n,:)
     end if
     call fcn(ps_x(n,:), fcnval)
     !are we better now?
     if (fcnval < ps_bf(n)) then
        ps_bf(n) = fcnval
        ps_bx(n,:) = ps_x(n,:)
        !check if the global is better
        if (fcnval < ps_gbf) then
           ps_gbf = fcnval
           ps_gbx = ps_x(n,:)
           !write(*,*) 'New global best, from particle',n
           !write(*,*) ps_gbf, ps_gbx
        end if
     end if
  end do
end subroutine iterate

subroutine swarm_optimize(nx, pop, lb, ub, n_iter, fcn, xguess)
  use stel_kinds
  implicit none
  integer :: n, pop, n_iter, i, nx
  real(dp), dimension(nx) :: lb, ub, xguess
  external fcn

  ps_npop = pop
  ps_n = nx

  call alloc_ps()
  
  ps_lb = lb
  ps_ub = ub
  call init_pop(fcn, xguess)
  
  if (ps_debug == 1) then
     write (*,*) 'initial setup'
     do n = 1,ps_npop
        write(*,*) '-----particle',n,'-------'
        write(*,*) 'pos',ps_x(n,:)
        write(*,*) 'vel',ps_v(n,:)
        write(*,*) 'score',ps_bf(n)
     end do
     write (*,*) 'initial best', ps_gbx
  end if

  do i = 1,n_iter
     call iterate(fcn)
     if (ps_debug == 2) then
        if (mod(i,10) == 0) then
           do n = 1,ps_npop
              write(*,*) n,ps_x(n,:)
           end do
        end if
     end if
     write (*,*) 'iter ',i,'complete, minimum value: ',ps_gbf
  end do

  xguess = ps_gbx
  print *,'final minimal calculated lambda',ps_gbf

end subroutine swarm_optimize
  

end module
