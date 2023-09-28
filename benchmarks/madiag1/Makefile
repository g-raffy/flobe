#PKG_CONFIG_PATH ${PKG_CONFIG_PATH}:/usr/local/magma/lib/pkgconfig


MAGMA_CFLAGS   = $(shell export PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:/opt/magma/magma-2.5.4/lib/pkgconfig ; pkg-config --cflags magma)

MAGMA_F90FLAGS = -Dmagma_devptr_t="integer(kind=8)" \
                  $(shell export PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:/opt/magma/magma-2.5.4/lib/pkgconfig ; pkg-config --cflags-only-I magma)

MAGMA_LIBS     = $(shell export PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:/opt/magma/magma-2.5.4/lib/pkgconfig ; pkg-config --libs   magma)


.PHONY: all
all: diag_mkl diag_mkl_omp diag_magma diag_magma_gpu diag_lapack

diag_lapack: diag_f.F90
	bash -l -c "gfortran -O3 -DUSE_DSYEV $< -o $@ -llapack"

diag_mkl: diag_f.F90
	bash -l -c "module load compilers/ifort/latest && ifort -DUSE_DSYEV -O3 $< -o $@ -lmkl_core -lmkl_intel_lp64 -lmkl_sequential"

diag_mkl_omp: diag_f.F90
	bash -l -c "module load compilers/ifort/latest && ifort -DUSE_DSYEV -O3 $< -o $@ -lmkl_intel_thread -lmkl_core -lmkl_intel_lp64 -liomp5 -lpthread -lm -ldl"

diag_magma: diag_f.F90
	bash -l -c 'module load compilers/ifort/latest && ifort -DUSE_MAGMA_DSYEVD $(MAGMA_F90FLAGS) -O3 $< -o $@ $(MAGMA_LIBS)'

diag_magma_gpu: diag_f.F90
	bash -l -c 'module load compilers/ifort/latest && ifort -DUSE_MAGMA_DSYEVD_GPU $(MAGMA_F90FLAGS) -O3 $< -o $@ $(MAGMA_LIBS)'

.PHONY: clean
clean:
	rm -Rf diag_mkl diag_mkl_omp diag_magma diag_magma_gpu diag_lapack *.mod
