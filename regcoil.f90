! Main program

program regcoil

  use global_variables, only: totalTime, outputFilename, general_option, my_pn, num_procs
  use init_plasma_mod
  use optimize

  implicit none
  include 'mpif.h'

  integer :: tic, toc, countrate, ierr

  call MPI_INIT(ierr)
  call MPI_COMM_RANK (MPI_COMM_WORLD, my_pn, ierr)
  call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)

  if (my_pn == 0) then
     print *,"This is REGCOIL,"
     print *,"a regularized least-squares method for computing stellarator coil shapes."
  end if
  call system_clock(tic,countrate)

  call read_input()
  call validate_input()
  call compute_lambda()

  ! Define the position vector and normal vector at each grid point for the surfaces:
  call init_plasma()
  call init_coil_surface()

  ! Initialize some of the vectors and matrices needed:
  call read_bnorm()
  call build_matrices(.True.)

  select case (general_option)
  case (1)
     call solve()
  case (2)
     call compute_diagnostics_for_nescout_potential()
  case (3)
     call svd_scan()
  case (4)
     call auto_regularization_solve()
  case (5)
     call run_optimize() 
  case default
     print *,"Invalid general_option:",general_option
     stop
  end select

  call system_clock(toc)
  totalTime = real(toc-tic)/countrate
 
  if (my_pn == 0) then
     call write_output()
     print *,"REGCOIL complete. Total time=",totalTime,"sec."
     print *,"You can run regcoilPlot ",trim(outputFilename)," to plot results."
  end if
  call MPI_FINALIZE(ierr)

end program regcoil
