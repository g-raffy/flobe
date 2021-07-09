program diag
#ifdef USE_MAGMA_GPU
   use magma
!    use magma_param
!    use magma_dfortran
!    use magmablas_dfortran
#endif
implicit none


integer :: argc, info, ndim

character(len=8) :: arg0, arg1


call get_command_argument(0,arg0)

argc = command_argument_count()
if (argc /= 1) then
    write(6,'("Usage: ",a," NDIM, where NDIM is a positive integer")') trim(arg0);
    stop
end if

call get_command_argument(1,arg1,status=info)
if (info /= 0) then
    write(6,'("Error reading argument: info = ",i2)') info
    write(6,'("Usage: ",a," NDIM, where NDIM is a positive integer")') trim(arg0);
stop
end if

read(arg1,*,iostat=info) ndim
if (info /= 0) then
    write(6,'("Error converting argument to integer: info = ",i2)') info
    write(6,'("Usage: ",a," NDIM, where NDIM is a positive integer")') trim(arg0);
stop
end if


if (ndim < 1) then
    write(6,'("Usage: ",a," NDIM, where NDIM is a positive integer")') trim(arg0);
stop
end if

    call test_dsyev(ndim)

stop
end program diag

subroutine test_dsyev(ndim)
    integer, intent(in) :: ndim
    integer, parameter :: dp = kind(1.0d0)
    real :: tstart, tstop

    real(dp), allocatable :: evalues(:)
    real(dp), allocatable :: work(:)
    real(dp), allocatable :: amat(:,:)
    integer :: i, info, j, lwork
    integer :: ldda
#ifdef USE_MAGMA_GPU
    magma_devptr_t :: d_amat
    magma_devptr_t :: queue  !! really a CPU pointer
    real(dp), allocatable :: wa(:,:)
    integer, allocatable :: iwork(:)
    integer :: nb
#endif
    real(dp) :: x

    allocate(evalues(ndim))
    evalues = 0.0_dp

#ifdef USE_MAGMA_GPU
    ! If JOBZ = MagmaVec and N > 1, LWORK >= max( 2*N + N*NB, 1 + 6*N + 2*N**2 ).
! int lwork = max(2 * DIM + DIM * magma_get_dsytrd_nb(DIM), 1 + 6 * DIM + 2 * int(pow(DIM, 2)));
    ! nb = magma_get_dsytrd_nb(ndim)
    nb = 100
    write(6,*) 'nb = ', nb
    lwork = max(2 * ndim + ndim * nb, 1 + 6 * ndim + 2 * ndim * ndim)
    ldda = ceiling(real(ndim)/32)*32
#else
    lwork = 3*ndim
#endif
    allocate(work(lwork))
    
    allocate(amat(ndim, ndim))

    !call random_number(amat)
    !amat = 0.5_dp*(amat + transpose(amat))
    do j = 1, ndim
        do i = 1, j
           call random_number(x)
           amat(i,j) = x
           amat(j,i) = x
        end do
    end do

call cpu_time(tstart)

#ifdef USE_MAGMA_GPU
    !! allocate GPU memory
    info = magmaf_dmalloc( d_amat, ldda*ndim )
    if (d_amat == 0) then
        print "(a)", "failed to allocate d_amat"
        return
    endif

    ! copy A to dA
    call magmaf_queue_create( 0, queue )
    call magmaf_dsetmatrix( ndim, ndim, amat, ndim, d_amat, ldda, queue )
    call magmaf_queue_destroy( queue )

    ! subroutine magmaf_dsyevd_gpu( jobz, uplo, n, dA, ldda, w, wA, ldwa, work, lwork, iwork,  &
    !     liwork, info )
    ! character        :: jobz
    ! character        :: uplo
    ! integer          :: n
    ! magma_devptr_t   :: dA
    ! integer          :: ldda
    ! double precision :: w(*)
    ! double precision :: wA(*)
    ! integer          :: ldwa
    ! double precision :: work(*)
    ! integer          :: lwork
    ! integer          :: iwork(*)
    ! integer          :: liwork
    ! integer          :: info
    ! end
    ! [in]	jobz	magma_vec_t

    !     = MagmaNoVec: Compute eigenvalues only;
    !     = MagmaVec: Compute eigenvalues and eigenvectors.

    ! [in]	uplo	magma_uplo_t

    !     = MagmaUpper: Upper triangle of A is stored;
    !     = MagmaLower: Lower triangle of A is stored.

    ! [in]	n	INTEGER The order of the matrix A. N >= 0.
    ! [in,out]	dA	DOUBLE PRECISION array on the GPU, dimension (LDDA, N). On entry, the symmetric matrix A. If UPLO = MagmaUpper, the leading N-by-N upper triangular part of A contains the upper triangular part of the matrix A. If UPLO = MagmaLower, the leading N-by-N lower triangular part of A contains the lower triangular part of the matrix A. On exit, if JOBZ = MagmaVec, then if INFO = 0, A contains the orthonormal eigenvectors of the matrix A. If JOBZ = MagmaNoVec, then on exit the lower triangle (if UPLO=MagmaLower) or the upper triangle (if UPLO=MagmaUpper) of A, including the diagonal, is destroyed.
    ! [in]	ldda	INTEGER The leading dimension of the array DA. LDDA >= max(1,N).
    ! [out]	w	DOUBLE PRECISION array, dimension (N) If INFO = 0, the eigenvalues in ascending order.
    ! 	wA	(workspace) DOUBLE PRECISION array, dimension (LDWA, N)
    ! [in]	ldwa	INTEGER The leading dimension of the array wA. LDWA >= max(1,N).
    ! [out]	work	(workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)) On exit, if INFO = 0, WORK[0] returns the optimal LWORK.
    ! [in]	lwork	INTEGER The length of the array WORK.

    !     If N <= 1, LWORK >= 1.
    !     If JOBZ = MagmaNoVec and N > 1, LWORK >= 2*N + N*NB.
    !     If JOBZ = MagmaVec and N > 1, LWORK >= max( 2*N + N*NB, 1 + 6*N + 2*N**2 ). NB can be obtained through magma_get_dsytrd_nb(N).
    !     If LWORK = -1, then a workspace query is assumed; the routine only calculates the optimal sizes of the WORK and IWORK arrays, returns these values as the first entries of the WORK and IWORK arrays, and no error message related to LWORK or LIWORK is issued by XERBLA.

    ! [out]	iwork	(workspace) INTEGER array, dimension (MAX(1,LIWORK)) On exit, if INFO = 0, IWORK[0] returns the optimal LIWORK.
    ! [in]	liwork	INTEGER The dimension of the array IWORK.

    !     If N <= 1, LIWORK >= 1.
    !     If JOBZ = MagmaNoVec and N > 1, LIWORK >= 1.
    !     If JOBZ = MagmaVec and N > 1, LIWORK >= 3 + 5*N.
    !     If LIWORK = -1, then a workspace query is assumed; the routine only calculates the optimal sizes of the WORK and IWORK arrays, returns these values as the first entries of the WORK and IWORK arrays, and no error message related to LWORK or LIWORK is issued by XERBLA.

    ! [out]	info	INTEGER

    !     = 0: successful exit
    !     < 0: if INFO = -i, the i-th argument had an illegal value
    !     > 0: if INFO = i and JOBZ = MagmaNoVec, then the algorithm failed to converge; i off-diagonal elements of an intermediate tridiagonal form did not converge to zero; if INFO = i and JOBZ = MagmaVec, then the algorithm failed to compute an eigenvalue while working on the submatrix lying in rows and columns INFO/(N+1) through mod(INFO,N+1).

    allocate(wa(ldda, ndim))
    allocate(iwork(1))

    call magmaf_dsyevd_gpu ('n', 'u', ndim, d_amat, ldda, evalues, wa, ldda, work, lwork, iwork, 1, info)
#else
    ! subroutine dsyev   (       character       JOBZ,
    !            character       UPLO,
    !            integer         N,
    !            double precision, dimension( lda, * )   A,
    !            integer         LDA,
    !            double precision, dimension( * )        W,
    !            double precision, dimension( * )        WORK,
    !            integer         LWORK,
    !            integer         INFO 
    !    )               
    call dsyev ('n','u',ndim,amat,ndim,evalues,work,lwork,info)
#endif
    call cpu_time(tstop)
    write(6,*) 'info = ', info
    
    write(6,'("Time taken by dsyev for matrix size ",i8," was ",f10.2," seconds")') ndim, tstop-tstart
    write(6,'(10f15.7)') evalues(1:10)
    write(6,'(10f15.7)') evalues(ndim-9:ndim)
    
    deallocate(evalues,work,amat)
    end
     