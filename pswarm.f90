!A simple particle swarm implementation with mpi capabilities.
!Each processor follows a set chunk of particles only, and after
!each iteration the global best position and values are calculated
!and broadcast. Then the particle positions are updated.
!
! To do- 
! convergence checking
! randomize based off of date

module pswarm
  use stel_kinds
  implicit none
  !current position, best position and velocity
  include 'mpif.h'
  real(dp), dimension(:,:), allocatable :: ps_x, ps_v, ps_bx 
  real(dp), dimension(:), allocatable :: ps_gbx !best known position vector
  real(dp), dimension(:), allocatable :: ps_lb, ps_ub
  real(dp), dimension(:), allocatable :: ps_bf
  real(dp) :: ps_gbf, ps_gbind
  integer :: ps_n, ps_npop, ps_gb
  integer :: nprocs
  integer :: ps_debug=0
  integer :: procnum
  !MPI variables, indices for each processor to know when to start and stop
  integer, dimension(:), allocatable :: startinds, endinds 

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
  allocate(startinds(nprocs))
  allocate(endinds(nprocs))
end subroutine alloc_ps



  

  

!Initialize the particles with random values
subroutine init_pop(fcn, xguess)

  use stel_kinds
  implicit none
  integer :: n, i, ierr, k, status(MPI_STATUS_SIZE)
  real(dp) :: x, diff
  real(dp), dimension(ps_n) :: xguess 
  real(dp), dimension(nprocs) :: gbf_arr !best indices and best scores
  
  real(dp) :: gbf_temp
  integer :: send_data_tag = 2001, return_data_tag = 2002
  integer :: proc_with_min=0, root_pn

  external fcn

  root_pn = 0
  
  !do the first one and set up global best
  !we have all procs do this, because they need to wait anyway
  ps_x(1,:) = xguess
  ps_v(1,:) = 0
  call fcn(ps_x(1,:), ps_bf(1))
  ps_gbf = ps_bf(1)
  ps_bx(1,:) = ps_x(1,:)
  ps_gbx = ps_x(1,:)

  
  
  do n = startinds(procnum+1), endinds(procnum+1)
     if (n==1) cycle !skip the first one, since we've already done it
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
        ps_gbind = n
     end if 
  end do
  
  
  !recover the global best
  gbf_arr = 0
  if (procnum == 0) then
     gbf_arr(1) = ps_gbf
     do i = 1,nprocs-1
        call mpi_irecv(gbf_temp,1,MPI_REAL8,i,11,MPI_COMM_WORLD,&
             status, ierr)
        gbf_arr(i+1) = gbf_temp
     end do
  else
     gbf_temp = ps_gbf
     call mpi_isend(gbf_temp,1,MPI_REAL8,0,11, MPI_COMM_WORLD, status, ierr)
  end if
  
  !which processor found the minimum?
  if (procnum == 0) then
     proc_with_min = minloc(gbf_arr, 1)-1     
  end if
  
  
  call mpi_bcast(proc_with_min, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  
  
  !broadcast best result to all other processors
  call mpi_bcast(ps_gbf, 1, MPI_REAL8, proc_with_min, MPI_COMM_WORLD, ierr)
  call mpi_bcast(ps_gbx, ps_n, MPI_REAL8, proc_with_min, &
       MPI_COMM_WORLD, ierr)
  
  
  
  call mpi_barrier(MPI_COMM_WORLD, ierr)
  
end subroutine init_pop

subroutine iterate(fcn)
  use stel_kinds
  implicit none
  include 'mpif.h'
  external fcn
  integer :: n,i, ierr, status, proc_with_min
  real(dp) :: rp, rg, fcnval
  real(dp) :: omega, phip, phig
  real(dp), dimension(nprocs) :: gbf_arr !best indices and best scores
  real(dp) :: gbf_temp
  real(dp), dimension(ps_n) :: gbx_temp

  !For now we set the optimization parameters
  omega = 0.5
  phip = 1.0
  phig = 1.0

  
  
  ! update the velocity vectors
  do n = startinds(procnum+1), endinds(procnum+1)
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
  do n = startinds(procnum+1), endinds(procnum+1)
     
     
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

  !recover the global best
  if (procnum == 0) then
     gbf_arr(1) = ps_gbf
     do i = 1,nprocs-1
        call mpi_irecv(gbf_temp,1,MPI_REAL8,i,11,MPI_COMM_WORLD,&
             status, ierr)
        gbf_arr(i+1) = gbf_temp
     end do
  else
     gbf_temp = ps_gbf
     call mpi_isend(gbf_temp,1,MPI_REAL8,0,11, MPI_COMM_WORLD, &
          status, ierr)
  end if
  
  
  !which processor found the minimum?
  if (procnum == 0) then
     proc_with_min = minloc(gbf_arr, 1)-1     
  end if
  
  call mpi_bcast(proc_with_min, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  !broadcast best result to all other processors
  call mpi_bcast(ps_gbf, 1, MPI_REAL8, proc_with_min, &
       MPI_COMM_WORLD, ierr)
  call mpi_bcast(ps_gbx, ps_n, MPI_REAL8, proc_with_min, &
       MPI_COMM_WORLD, ierr)
  call mpi_barrier(MPI_COMM_WORLD, ierr)
      
end subroutine iterate

subroutine get_mpi_inds()
  
  integer :: minperproc, overflow, i, prevend

  !The minimum number per proc
  minperproc = ps_npop/nprocs
  !This is the amount of extra variables
  overflow = ps_npop - (nprocs * minperproc)
  prevend = 0
  
  do i = 1,nprocs
     startinds(i) = prevend + 1
     if (overflow > 0) then
        endinds(i) = startinds(i) + minperproc
        overflow = overflow - 1
     else
        endinds(i) = startinds(i) + minperproc - 1
     end if
     prevend = endinds(i)
  end do
  


end subroutine get_mpi_inds

subroutine swarm_optimize(nx, pop, lb, ub, n_iter, fcn, xguess, nproc)
  use stel_kinds
  implicit none
  integer :: i,k
  integer :: n, pop, n_iter, nx, nproc, ierr
  real(dp), dimension(nx) :: lb, ub, xguess
  integer, dimension(:), allocatable :: seed
  external fcn

  call MPI_COMM_RANK(MPI_COMM_WORLD, procnum, ierr)

  ps_npop = pop
  ps_n = nx
  nprocs = nproc !this is not the best way to code this, but was having trouble passing stuff through fortran
  call alloc_ps()
  
  ps_lb = lb
  ps_ub = ub
  call get_mpi_inds()

  !initialize random number generator
  !will use date and time later, for now we want reproducibility
  call random_seed(size=k)
  allocate(seed(1:k))
  seed(:) = procnum+1
  call random_seed(put=seed)
  deallocate(seed)

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
     if (procnum == 0) & 
          write (*,*) 'iter ',i,'complete, minimum value: ',ps_gbf
  end do

  xguess = ps_gbx
  if (procnum == 0) print *,'final minimal calculated lambda',ps_gbf

end subroutine swarm_optimize
  

end module
