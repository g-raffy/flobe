#define MADIAG1_VERSION "1.0.0"
subroutine print_usage(prog_path)
    character(len=*), intent(in) :: prog_path
    write(6,'("madiag1 v",a,": benchmark that performs a square matrix diaginalization in double precision, using eigen vectors")') MADIAG1_VERSION;
    write(6,'()');
    write(6,'("Usage: ",a," <NDIM> <NUM_LOOPS>")') trim(prog_path);
    write(6,'("   <NDIM> positive integer representing the size of the square matrices to diagonalize")');
    write(6,'("   <NUM_LOOPS> positive integer representing the number of times the diagonalization is performed")');
end subroutine

program diag

implicit none


integer :: argc, info, ndim, num_loops

character(len=8) :: arg0, arg1, arg2


call get_command_argument(0,arg0)

argc = command_argument_count()
if (argc /= 2) then
    call print_usage(trim(arg0))
    stop
end if

call get_command_argument(1,arg1,status=info)
if (info /= 0) then
    write(6,'("Error reading argument: info = ",i2)') info
    call print_usage(trim(arg0))
stop
end if

call get_command_argument(2,arg2,status=info)
if (info /= 0) then
    write(6,'("Error reading argument: info = ",i2)') info
    call print_usage(trim(arg0))
stop
end if

read(arg1,*,iostat=info) ndim
if (info /= 0) then
    write(6,'("Error converting ndim argument to integer: info = ",i2)') info
    call print_usage(trim(arg0))
stop
end if

read(arg2,*,iostat=info) num_loops
if (info /= 0) then
    write(6,'("Error converting num_loops argument to integer: info = ",i2)') info
    call print_usage(trim(arg0))
stop
end if


if (ndim < 1) then
    call print_usage(trim(arg0))
stop
end if

    call test_dsyev(ndim, num_loops)

stop
end program diag

module mod_diag
    contains
    !---------------------------------------------------------!
    !Calls the LAPACK diagonalization subroutine DSYEV        !
    !input:  amat(ndim,ndim) = real symmetric matrix to be diagonalized!
    !            ndim  = size of amat                               !
    !output: amat(ndim,ndim) = orthonormal eigenvectors of amat           !
    !        eig(ndim) = eigenvalues of amat in ascending order     !
    !---------------------------------------------------------!
subroutine diagonalize(amat, ndim, evalues, info)
#if defined(USE_MAGMA_DSYEVD) || defined(USE_MAGMA_DSYEVD_GPU)
    use magma
 !    use magma_param
 !    use magma_dfortran
 !    use magmablas_dfortran
#endif
 
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    real(dp), allocatable :: amat(:,:)
    integer, intent(in) :: ndim
    real(dp), allocatable :: evalues(:)
    integer, intent(out) :: info

    real(dp), allocatable :: work(:)
    integer :: lwork
    integer :: lda
#ifdef USE_MAGMA_DSYEVD_GPU
    magma_devptr_t :: d_amat
    magma_devptr_t :: queue  !! really a CPU pointer
    real(dp), allocatable :: wa(:,:)
#endif

#if defined(USE_MAGMA_DSYEVD) || defined(USE_MAGMA_DSYEVD_GPU)
    integer :: liwork
    integer, allocatable :: iwork(:)
    integer :: nb
#endif
    logical :: compute_eigenvectors = .TRUE.
    character :: jobz

    if (compute_eigenvectors) then
        jobz = 'v' ! novec
    else
        jobz = 'n' ! vec
    endif

#if defined(USE_MAGMA_DSYEVD) || defined(USE_MAGMA_DSYEVD_GPU)
    call magmaf_init()
#endif
    
#if defined(USE_MAGMA_DSYEVD) || defined(USE_MAGMA_DSYEVD_GPU)
    ! NB can be obtained through magma_get_dsytrd_nb(N).
    ! nb = magma_get_dsytrd_nb(ndim)
    ! magma_get_dsytrd_nb doesn't seem to exist so I set an arbitrary value here
    nb = 100
    write(6,*) 'nb = ', nb

    if (compute_eigenvectors) then
        ! If JOBZ = MagmaVec and N > 1, LWORK >= max( 2*N + N*NB, 1 + 6*N + 2*N**2 ).
        ! int lwork = max(2 * DIM + DIM * magma_get_dsytrd_nb(DIM), 1 + 6 * DIM + 2 * int(pow(DIM, 2)));
        lwork = max(2 * ndim + ndim * nb, 1 + 6 * ndim + 2 * ndim * ndim)
    else
        ! If JOBZ = MagmaNoVec and N > 1, LWORK >= 2*N + N*NB.
        lwork = 2 * ndim + ndim * nb
    endif
#else
    lwork = 3*ndim
#endif
    allocate(work(lwork))

    lda = ceiling(real(ndim)/32)*32
    
#if defined(USE_MAGMA_DSYEVD_GPU)
    !! allocate GPU memory
    info = magmaf_dmalloc( d_amat, lda*ndim )
    if (d_amat == 0) then
        print "(a)", "failed to allocate d_amat"
        return
    endif

    ! copy A to dA
    call magmaf_queue_create( 0, queue )
    call magmaf_dsetmatrix( ndim, ndim, amat, ndim, d_amat, lda, queue )
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

    allocate(wa(lda, ndim))
#endif

#if defined(USE_MAGMA_DSYEVD) || defined(USE_MAGMA_DSYEVD_GPU)
    if (compute_eigenvectors) then
        liwork = 5 * ndim + 3
    else
        liwork = 1
    endif
    allocate(iwork(liwork))
#endif

#if defined(USE_MAGMA_DSYEVD_GPU)
    call magmaf_dsyevd_gpu (jobz, 'u', ndim, d_amat, lda, evalues, wa, lda, work, lwork, iwork, liwork, info)
#endif

#ifdef USE_MAGMA_DSYEVD
    call magmaf_dsyevd (jobz, 'u', ndim, amat, ndim, evalues, work, lwork, iwork, liwork, info)
#endif

#ifdef USE_DSYEV

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
    call dsyev (jobz,'u',ndim,amat,ndim,evalues,work,lwork,info)
#endif

    if (info .ne. 0) then
        stop 1
    end if

    deallocate(work)

#if defined(USE_MAGMA_DSYEVD) || defined(USE_MAGMA_DSYEVD_GPU)
    call magmaf_finalize()
#endif

end subroutine diagonalize

endmodule mod_diag

subroutine set_random_seed(seed)
    integer :: seed
    integer :: seed_array_size
    INTEGER, ALLOCATABLE :: seed_array (:)
    CALL RANDOM_SEED (SIZE = seed_array_size)  ! I is set to the size of
    !                              ! the seed array
    ALLOCATE (seed_array(seed_array_size))
    seed_array = seed
    CALL RANDOM_SEED (PUT=seed_array(1:seed_array_size))
end subroutine

subroutine test_dsyev(ndim, num_loops)
    use mod_diag, only: diagonalize
    implicit none
    integer, intent(in) :: ndim
    integer, intent(in) :: num_loops
    integer, parameter :: dp = kind(1.0d0)
    real :: tstart, tstop

    INTEGER :: c1,c2,cr,cm,s
    REAL :: a_diff, diff, rate

    real(dp), allocatable :: evalues(:)
    real(dp), allocatable :: amat(:,:),playmat(:,:)
    real(dp) :: x
    integer :: i, j, info

    ! First initialize the system_clock
    CALL system_clock(count_rate=cr)
    CALL system_clock(count_max=cm)
    rate = REAL(cr)
    WRITE(*,*) "system_clock rate ",rate

    diff = 0.0
    a_diff = 0.0
    s = 0

    allocate(amat(ndim, ndim))
    allocate(playmat(ndim, ndim))

    call set_random_seed(42)

    !call random_number(amat)
    !amat = 0.5_dp*(amat + transpose(amat))
    do j = 1, ndim
        do i = 1, j
           call random_number(x)
           amat(i,j) = x
           amat(j,i) = x
        end do
    end do

    allocate(evalues(ndim))
    evalues = 0.0_dp


    call cpu_time(tstart)
    call system_clock(c1)

    do j = 1, num_loops
        playmat = amat
        call diagonalize(playmat, ndim, evalues, info)
    end do

    call cpu_time(tstop)
    call system_clock(c2)
    if ( (c2 - c1)/rate < (tstop - tstart) ) s = s + 1
    diff = (c2 - c1)/rate - (tstop - tstart) + diff
    a_diff = ABS((c2 - c1)/rate - (tstop - tstart)) + a_diff

    write(6,*) 'info = ', info
    
    write(6,'("Time taken by dsyev for matrix size ",i8," was ",f10.2," seconds")') ndim, tstop-tstart
    write(6,'(10f15.7)') evalues(1:10)
    write(6,'(10f15.7)') evalues(ndim-9:ndim)
    
    WRITE(*,*) "system_clock : ",(c2 - c1)/rate
    WRITE(*,*) "cpu_time     : ",(tstop - tstart)
    WRITE(*,*) "sc < ct      : ",s
    WRITE(*,*) "mean diff    : ",diff
    WRITE(*,*) "abs mean diff: ",a_diff

    deallocate(evalues, amat, playmat)
    end
     