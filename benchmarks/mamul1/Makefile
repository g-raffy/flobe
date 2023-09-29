#PKG_CONFIG_PATH ${PKG_CONFIG_PATH}:/usr/local/magma/lib/pkgconfig


MAGMA_CFLAGS   = $(shell export PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:/opt/magma/magma-2.7.0/lib/pkgconfig ; pkg-config --cflags magma)

MAGMA_F90FLAGS = -Dmagma_devptr_t="integer(kind=8)" \
                  $(shell export PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:/opt/magma/magma-2.7.0/lib/pkgconfig ; pkg-config --cflags-only-I magma)

MAGMA_LIBS     = $(shell export PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:/opt/magma/magma-2.7.0/lib/pkgconfig ; pkg-config --libs   magma)


.PHONY: all
all: mamul1_mkl mamul1_mkl_omp mamul1_magma_gpu mamul1_lapack

mamul1_lapack: mamul1.F90
	bash -l -c "gfortran -O3 -DUSE_DGEMM $< -o $@ -llapack"

mamul1_mkl: mamul1.F90
	bash -l -c "module load compilers/ifort/latest && ifort -DUSE_DGEMM -O3 $< -o $@ -lmkl_core -lmkl_intel_lp64 -lmkl_sequential"

mamul1_mkl_omp: mamul1.F90
	bash -l -c "module load compilers/ifort/latest && ifort -DUSE_DGEMM -O3 $< -o $@ -lmkl_intel_thread -lmkl_core -lmkl_intel_lp64 -liomp5 -lpthread -lm -ldl"

mamul1_magma_gpu: mamul1.F90
	bash -l -c 'module load compilers/ifort/latest && ifort -DUSE_MAGMA_DGEMM_GPU $(MAGMA_F90FLAGS) -O3 $< -o $@ $(MAGMA_LIBS)'

.PHONY: clean
clean:
	rm -Rf mamul1_mkl mamul1_mkl_omp mamul1_magma mamul1_magma_gpu mamul1_lapack *.mod
