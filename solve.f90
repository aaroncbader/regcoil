subroutine solve

  use global_variables
  use stel_constants
  use stel_kinds

  implicit none

  integer :: iflag, tic, toc, countrate
  real(dp), dimension(:,:), allocatable :: matrix, this_current_potential
  real(dp), dimension(:), allocatable :: RHS, solution
  real(dp), dimension(:), allocatable :: KDifference_x, KDifference_y, KDifference_z
  real(dp), dimension(:,:), allocatable :: this_K2_times_N
  real(dp) :: factor_theta, factor_zeta
  integer :: ilambda, itheta, izeta

  ! Variables needed by LAPACK:
  integer :: INFO, LWORK
  real(dp), dimension(:), allocatable :: WORK
  integer, dimension(:), allocatable :: IPIV

  allocate(matrix(num_basis_functions, num_basis_functions), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(RHS(num_basis_functions), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(solution(num_basis_functions), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(WORK(1), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(IPIV(num_basis_functions), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'

  allocate(chi2_B(nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(chi2_K(nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(max_Bnormal(nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(max_K(nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(current_potential(ntheta_coil,nzeta_coil,nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(single_valued_current_potential_thetazeta(ntheta_coil,nzeta_coil,nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(this_current_potential(ntheta_coil,nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(single_valued_current_potential_mn(num_basis_functions,nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(Bnormal_total(ntheta_plasma,nzeta_plasma,nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(K2(ntheta_coil,nzeta_coil,nlambda), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(KDifference_x(ntheta_coil*nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(KDifference_y(ntheta_coil*nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(KDifference_z(ntheta_coil*nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'
  allocate(this_K2_times_N(ntheta_coil,nzeta_coil), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'


  ! Call LAPACK's DSYSV in query mode to determine the optimal size of the work array
  call DSYSV('U',num_basis_functions, 1, matrix, num_basis_functions, IPIV, RHS, num_basis_functions, WORK, -1, INFO)
  LWORK = WORK(1)
  print *,"Optimal LWORK:",LWORK
  deallocate(WORK)
  allocate(WORK(LWORK), stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error!'  

  factor_zeta  = net_poloidal_current_Amperes / twopi
  factor_theta = net_toroidal_current_Amperes / twopi

  do ilambda = 1,nlambda
     print "(a,e10.3,a,i3,a,i3,a)"," Solving system for lambda=",lambda(ilambda)," (",ilambda," of ",nlambda,")"
     call system_clock(tic,countrate)

     matrix = matrix_B + lambda(ilambda) * matrix_K
     RHS    =    RHS_B + lambda(ilambda) *    RHS_K

     call system_clock(toc)
     print *,"  Additions: ",real(toc-tic)/countrate," sec."
     call system_clock(tic)

     ! Compute solution = matrix \ RHS.
     ! Use LAPACK's DSYSV since matrix is symmetric.
     ! Note: RHS will be over-written with the solution.
     call DSYSV('U',num_basis_functions, 1, matrix, num_basis_functions, IPIV, RHS, num_basis_functions, WORK, LWORK, INFO)
     if (INFO /= 0) then
        print *, "!!!!!! Error in LAPACK DSYSV: INFO = ", INFO
        !stop
     end if
     solution = RHS

     call system_clock(toc)
     print *,"  DSYSV: ",real(toc-tic)/countrate," sec."
     call system_clock(tic)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! Now compute diagnostics
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     single_valued_current_potential_mn(:,ilambda) = solution
     this_current_potential = reshape(matmul(basis_functions, solution), (/ ntheta_coil, nzeta_coil /)) ! Could use BLAS2 for this line for extra speed.
     single_valued_current_potential_thetazeta(:,:,ilambda) = this_current_potential
     do izeta = 1,nzeta_coil
        do itheta = 1,ntheta_coil
           this_current_potential(itheta,izeta) = this_current_potential(itheta,izeta) &
                + factor_zeta*zeta_coil(izeta) + factor_theta*theta_coil(itheta)
        end do
     end do
     current_potential(:,:,ilambda) = this_current_potential

     KDifference_x = d_x - matmul(f_x, solution)
     KDifference_y = d_y - matmul(f_y, solution)
     KDifference_z = d_z - matmul(f_z, solution)
     this_K2_times_N = reshape(KDifference_x*KDifference_x + KDifference_y*KDifference_y + KDifference_z*KDifference_z, (/ ntheta_coil, nzeta_coil /)) &
          / norm_normal_coil
     chi2_K(ilambda) = nfp * dtheta_coil * dzeta_coil * sum(this_K2_times_N)
     K2(:,:,ilambda) = this_K2_times_N / norm_normal_coil

     Bnormal_total(:,:,ilambda) = (reshape(matmul(g,solution),(/ ntheta_plasma, nzeta_plasma /)) / norm_normal_plasma) &
          + Bnormal_from_plasma_current + Bnormal_from_net_coil_currents

     max_Bnormal(ilambda) = maxval(abs(Bnormal_total(:,:,ilambda)))
     max_K(ilambda) = sqrt(maxval(K2(:,:,ilambda)))

     chi2_B(ilambda) = nfp * dtheta_plasma * dzeta_plasma &
          * sum(Bnormal_total(:,:,ilambda) * Bnormal_total(:,:,ilambda) * norm_normal_plasma)

     call system_clock(toc)
     print *,"  Diagnostics: ",real(toc-tic)/countrate," sec."
     print *,"  chi2_B:",chi2_B(ilambda),"chi2_K:",chi2_K(ilambda)
  end do

end subroutine solve

    ! Here is the LAPACK documentation for solving a symmetric linear system:

!!$
!!$
!!$*> \brief <b> DSYSV computes the solution to system of linear equations A * X = B for SY matrices</b>
!!$*
!!$*  =========== DOCUMENTATION ===========
!!$*
!!$* Online html documentation available at 
!!$*            http://www.netlib.org/lapack/explore-html/ 
!!$*
!!$*> \htmlonly
!!$*> Download DSYSV + dependencies 
!!$*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsysv.f"> 
!!$*> [TGZ]</a> 
!!$*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsysv.f"> 
!!$*> [ZIP]</a> 
!!$*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsysv.f"> 
!!$*> [TXT]</a>
!!$*> \endhtmlonly 
!!$*
!!$*  Definition:
!!$*  ===========
!!$*
!!$*       SUBROUTINE DSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,
!!$*                         LWORK, INFO )
!!$* 
!!$*       .. Scalar Arguments ..
!!$*       CHARACTER          UPLO
!!$*       INTEGER            INFO, LDA, LDB, LWORK, N, NRHS
!!$*       ..
!!$*       .. Array Arguments ..
!!$*       INTEGER            IPIV( * )
!!$*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * )
!!$*       ..
!!$*  
!!$*
!!$*> \par Purpose:
!!$*  =============
!!$*>
!!$*> \verbatim
!!$*>
!!$*> DSYSV computes the solution to a real system of linear equations
!!$*>    A * X = B,
!!$*> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
!!$*> matrices.
!!$*>
!!$*> The diagonal pivoting method is used to factor A as
!!$*>    A = U * D * U**T,  if UPLO = 'U', or
!!$*>    A = L * D * L**T,  if UPLO = 'L',
!!$*> where U (or L) is a product of permutation and unit upper (lower)
!!$*> triangular matrices, and D is symmetric and block diagonal with
!!$*> 1-by-1 and 2-by-2 diagonal blocks.  The factored form of A is then
!!$*> used to solve the system of equations A * X = B.
!!$*> \endverbatim
!!$*
!!$*  Arguments:
!!$*  ==========
!!$*
!!$*> \param[in] UPLO
!!$*> \verbatim
!!$*>          UPLO is CHARACTER*1
!!$*>          = 'U':  Upper triangle of A is stored;
!!$*>          = 'L':  Lower triangle of A is stored.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] N
!!$*> \verbatim
!!$*>          N is INTEGER
!!$*>          The number of linear equations, i.e., the order of the
!!$*>          matrix A.  N >= 0.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] NRHS
!!$*> \verbatim
!!$*>          NRHS is INTEGER
!!$*>          The number of right hand sides, i.e., the number of columns
!!$*>          of the matrix B.  NRHS >= 0.
!!$*> \endverbatim
!!$*>
!!$*> \param[in,out] A
!!$*> \verbatim
!!$*>          A is DOUBLE PRECISION array, dimension (LDA,N)
!!$*>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!!$*>          N-by-N upper triangular part of A contains the upper
!!$*>          triangular part of the matrix A, and the strictly lower
!!$*>          triangular part of A is not referenced.  If UPLO = 'L', the
!!$*>          leading N-by-N lower triangular part of A contains the lower
!!$*>          triangular part of the matrix A, and the strictly upper
!!$*>          triangular part of A is not referenced.
!!$*>
!!$*>          On exit, if INFO = 0, the block diagonal matrix D and the
!!$*>          multipliers used to obtain the factor U or L from the
!!$*>          factorization A = U*D*U**T or A = L*D*L**T as computed by
!!$*>          DSYTRF.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] LDA
!!$*> \verbatim
!!$*>          LDA is INTEGER
!!$*>          The leading dimension of the array A.  LDA >= max(1,N).
!!$*> \endverbatim
!!$*>
!!$*> \param[out] IPIV
!!$*> \verbatim
!!$*>          IPIV is INTEGER array, dimension (N)
!!$*>          Details of the interchanges and the block structure of D, as
!!$*>          determined by DSYTRF.  If IPIV(k) > 0, then rows and columns
!!$*>          k and IPIV(k) were interchanged, and D(k,k) is a 1-by-1
!!$*>          diagonal block.  If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0,
!!$*>          then rows and columns k-1 and -IPIV(k) were interchanged and
!!$*>          D(k-1:k,k-1:k) is a 2-by-2 diagonal block.  If UPLO = 'L' and
!!$*>          IPIV(k) = IPIV(k+1) < 0, then rows and columns k+1 and
!!$*>          -IPIV(k) were interchanged and D(k:k+1,k:k+1) is a 2-by-2
!!$*>          diagonal block.
!!$*> \endverbatim
!!$*>
!!$*> \param[in,out] B
!!$*> \verbatim
!!$*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
!!$*>          On entry, the N-by-NRHS right hand side matrix B.
!!$*>          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] LDB
!!$*> \verbatim
!!$*>          LDB is INTEGER
!!$*>          The leading dimension of the array B.  LDB >= max(1,N).
!!$*> \endverbatim
!!$*>
!!$*> \param[out] WORK
!!$*> \verbatim
!!$*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!!$*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!!$*> \endverbatim
!!$*>
!!$*> \param[in] LWORK
!!$*> \verbatim
!!$*>          LWORK is INTEGER
!!$*>          The length of WORK.  LWORK >= 1, and for best performance
!!$*>          LWORK >= max(1,N*NB), where NB is the optimal blocksize for
!!$*>          DSYTRF.
!!$*>          for LWORK < N, TRS will be done with Level BLAS 2
!!$*>          for LWORK >= N, TRS will be done with Level BLAS 3
!!$*>
!!$*>          If LWORK = -1, then a workspace query is assumed; the routine
!!$*>          only calculates the optimal size of the WORK array, returns
!!$*>          this value as the first entry of the WORK array, and no error
!!$*>          message related to LWORK is issued by XERBLA.
!!$*> \endverbatim
!!$*>
!!$*> \param[out] INFO
!!$*> \verbatim
!!$*>          INFO is INTEGER
!!$*>          = 0: successful exit
!!$*>          < 0: if INFO = -i, the i-th argument had an illegal value
!!$*>          > 0: if INFO = i, D(i,i) is exactly zero.  The factorization
!!$*>               has been completed, but the block diagonal matrix D is
!!$*>               exactly singular, so the solution could not be computed.
!!$*> \endverbatim
