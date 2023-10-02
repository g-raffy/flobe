
enable_language (Fortran)

set(MADIAG1_USE_MAGMA "OFF" CACHE BOOL "if set, madiag1 build uses magma (matrix algebra on gpu)")

set(MADIAG1_MAGMA_API "CPU_MEM_API" CACHE STRING "which magma API to use when building madiag1: CPU_MEM_API for BLAS compatible API (uses matrices stored on CPU memory) or GPU_MEM_API (use matrices stored on GPU memory)")

add_executable(madiag1 diag_f.F90)

# Allow arbitrary long lines. Needed as preprocessing could generate long line lengths.
target_compile_options(madiag1 PUBLIC -ffree-line-length-none)

if (MADIAG1_USE_MAGMA)
	find_package( MAGMA REQUIRED )
	if( MADIAG1_MAGMA_API STREQUAL "CPU_MEM_API" )
		target_compile_definitions(madiag1 PUBLIC MADIAG1_USE_MAGMA_DSYEVD)
	elseif( MADIAG1_MAGMA_API STREQUAL "GPU_MEM_API" )
		target_compile_definitions(madiag1 PUBLIC MADIAG1_USE_MAGMA_DSYEVD_GPU)
	else()
		message(FATAL_ERROR "unexpected value for MADIAG1_MAGMA_API : ${MADIAG1_MAGMA_API}")
	endif()
	target_link_libraries(madiag1 "${MAGMA_LIBRARIES}")
else()
	find_package( BLAS REQUIRED )
	find_package( LAPACK REQUIRED )
	# message("BLAS_LIBRARIES=${BLAS_LIBRARIES}")
	# message("LAPACK_LIBRARIES=${LAPACK_LIBRARIES}")
	target_compile_definitions(madiag1 PUBLIC USE_DSYEV)

	# Link Blas and Lapack libraries
	target_link_libraries(madiag1 "${LAPACK_LIBRARIES}")
	target_link_libraries(madiag1 "${BLAS_LIBRARIES}")
endif()

install(TARGETS madiag1)
