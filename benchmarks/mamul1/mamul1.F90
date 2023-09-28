program mamul1

implicit none


integer :: argc, info, ndim, num_loops

character(len=8) :: arg0, arg1, arg2


call get_command_argument(0,arg0)

argc = command_argument_count()
if (argc /= 2) then
    write(6,'("Usage: ",a," NDIM NUM_LOOPS, where NDIM is a positive integer")') trim(arg0);
    stop
end if

call get_command_argument(1,arg1,status=info)
if (info /= 0) then
    write(6,'("Error reading argument: info = ",i2)') info
    write(6,'("Usage: ",a," NDIM NUM_LOOPS, where NDIM is a positive integer")') trim(arg0);
stop
end if

call get_command_argument(2,arg2,status=info)
if (info /= 0) then
    write(6,'("Error reading argument: info = ",i2)') info
    write(6,'("Usage: ",a," NDIM NUM_LOOPS, where NDIM is a positive integer")') trim(arg0);
stop
end if

read(arg1,*,iostat=info) ndim
if (info /= 0) then
    write(6,'("Error converting ndim argument to integer: info = ",i2)') info
    write(6,'("Usage: ",a," NDIM NUM_LOOPS, where NDIM is a positive integer")') trim(arg0);
stop
end if

read(arg2,*,iostat=info) num_loops
if (info /= 0) then
    write(6,'("Error converting num_loops argument to integer: info = ",i2)') info
    write(6,'("Usage: ",a," NDIM NUM_LOOPS, where NDIM is a positive integer")') trim(arg0);
stop
end if


if (ndim < 1) then
    write(6,'("Usage: ",a," NDIM NUM_LOOPS, where NDIM is a positive integer")') trim(arg0);
stop
end if

    call test_dgemm(ndim, num_loops)

stop
end program mamul1

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

subroutine print_matrix(mat, ndim)
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    real(dp), intent(in) :: mat(ndim, ndim)
    integer, intent(in) :: ndim
    integer :: irow
    do irow = 1, ndim
        write(6, *) mat(irow,:)
    end do
end subroutine

subroutine test_dgemm(ndim, num_loops)
#if defined(USE_MAGMA_DGEMM) || defined(USE_MAGMA_DGEMM_GPU)
    use magma, only: magmaf_init, magmaf_finalize
    use magma, only: magmablasf_dgemm  !, magmaf_dgemm_gpu
#endif

    implicit none
    integer, intent(in) :: ndim
    integer, intent(in) :: num_loops
    integer, parameter :: dp = kind(1.0d0)
    real :: tstart, tstop
    integer(8) :: num_ops
    real :: gflops

    INTEGER :: c1,c2,cr,cm,s
    REAL :: a_diff, diff, rate

    real*8 :: amat(ndim,ndim)
    real*8 :: bmat(ndim,ndim)
    real*8 :: cmat(ndim,ndim)
    ! real*8, allocatable :: amat(:,:)
    ! real*8, allocatable :: bmat(:,:)
    ! real*8, allocatable :: cmat(:,:)
    real(dp) :: x
    integer :: i, j

#if defined(USE_MAGMA_DGEMM) || defined(USE_MAGMA_DGEMM_GPU)
    call magmaf_init()
#endif

    ! First initialize the system_clock
    CALL system_clock(count_rate=cr)
    CALL system_clock(count_max=cm)
    rate = REAL(cr)
    WRITE(*,*) "system_clock rate ",rate

    diff = 0.0
    a_diff = 0.0
    s = 0

    ! allocate(amat(ndim, ndim))
    ! allocate(bmat(ndim, ndim))
    ! allocate(cmat(ndim, ndim))

    call set_random_seed(42)

    !call random_number(amat)
    !amat = 0.5_dp*(amat + transpose(amat))
    do j = 1, ndim
        do i = 1, ndim
           call random_number(x)
           amat(i,j) = x
           call random_number(x)
           bmat(i,j) = x
        end do
    end do

    call cpu_time(tstart)
    call system_clock(c1)

    do j = 1, num_loops
        ! playmat = amat

#if defined(USE_MAGMA_DGEMM_GPU)
        ! call magmaf_dgemm_gpu ()
#endif

#ifdef USE_MAGMA_DGEMM
        call magmablasf_dgemm ('N', 'N', ndim, ndim, ndim, 1.0d0, amat, ndim, bmat, ndim, 0.0d0, cmat, ndim)
#endif

#ifdef USE_DGEMM
        ! subroutine dgemm 	( 	character  	TRANSA,
        ! 		character  	TRANSB,
        ! 		integer  	M,
        ! 		integer  	N,
        ! 		integer  	K,
        ! 		double precision  	ALPHA,
        ! 		double precision, dimension(lda,*)  	A,
        ! 		integer  	LDA,
        ! 		double precision, dimension(ldb,*)  	B,
        ! 		integer  	LDB,
        ! 		double precision  	BETA,
        ! 		double precision, dimension(ldc,*)  	C,
        ! 		integer  	LDC 
        ! 	) 	        
        call dgemm('N', 'N', ndim, ndim, ndim, 1.0d0, amat, ndim, bmat, ndim, 0.0d0, cmat, ndim)
#endif

    end do

    call cpu_time(tstop)
    call system_clock(c2)
    if ( (c2 - c1)/rate < (tstop - tstart) ) s = s + 1
    diff = (c2 - c1)/rate - (tstop - tstart) + diff
    a_diff = ABS((c2 - c1)/rate - (tstop - tstart)) + a_diff

    ! check one of the elements of cmat (the last one here: cmat(ndim, ndim))
    x = 0.0d0
    do i = 1, ndim
       x = x + amat(ndim, i) * bmat(i, ndim)
    end do

    ! write(6, *) 'amat = '
    ! call print_matrix(amat, ndim)

    ! write(6, *) 'bmat = '
    ! call print_matrix(bmat, ndim)

    ! write(6, *) 'cmat = '
    ! call print_matrix(cmat, ndim)
    write(6, '("expected cmat(", i0, ", ", i0, ")", e23.15e3)') ndim, ndim, x
    write(6, '("computed cmat(", i0, ", ", i0, ")", e23.15e3)') ndim, ndim, cmat(ndim, ndim)

    num_ops = real(ndim) * real(ndim) * real(ndim) * 2
    gflops = num_ops / (tstop-tstart) / 1.0e9


    write(6, '("Time taken by dgemm for matrix size ",i8," was ",f10.2," seconds")') ndim, tstop-tstart
    WRITE(*,*) "gflops       : ", gflops
    
    WRITE(*,*) "system_clock : ",(c2 - c1)/rate
    WRITE(*,*) "cpu_time     : ",(tstop - tstart)
    WRITE(*,*) "sc < ct      : ",s
    WRITE(*,*) "mean diff    : ",diff
    WRITE(*,*) "abs mean diff: ",a_diff

#if defined(USE_MAGMA_DGEMM) || defined(USE_MAGMA_DGEMM_GPU)
    call magmaf_finalize()
#endif


    deallocate(amat, bmat, cmat)
    end
