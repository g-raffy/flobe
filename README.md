# flobe (FLOating point BEnchmark)

flobe is a tool suite aimed at measuring floating point performance of processors

tsv files made by hand (running `madiag1_star.bash` without `madiag1_batch.bash`):
- `benchmarks/madiag1/measurements/fr.univ-rennes1.ipr.physix/diag_bench.tsv`
- `benchmarks/madiag1/measurements/fr.univ-rennes1.ipr.physix/diag_mkl_bench.tsv`
- [file://./benchmarks/madiag1/measurements/fr.univ-rennes1.ipr.physix/dsyev_bench.tsv]
	- 19106.00 measure made on 28/06/2021 with `graffy@physix90:/opt/ipr/cluster/work.local/graffy/bug3177$ time ./diag_mkl_omp 32000`
	- apparently made with 
- `benchmarks/madiag1/measurements/fr.univ-rennes1.ipr.physix/multi_diag_mkl_omp.tsv`

## how to use

### 1. run madiag1 benchmarks on each of the target machines

get a copy of flobe:
```sh
bob@physix90:/opt/ipr/work.local/bob$ git clone git@github.com:g-raffy/flobe.git
```

configure and build flobe
```sh
bob@physix90:/opt/ipr/work.local/bob$ mkdir -p ./flobe_build
bob@physix90:/opt/ipr/work.local/bob$ mkdir -p ./flobe_install
bob@physix90:/opt/ipr/work.local/bob$ cd ./flobe_build
``` 

basic configuration (uses default blas and lapack library):
```sh
bob@physix90:/opt/ipr/work.local/bob/flobe_build$ cmake -DCMAKE_INSTALL_PREFIX=~/flobe_install -DUSE_MAGMA=FALSE ../flobe
```

if you want to configure to use Intel Math Kernel Libraries:
```sh
bob@physix90:/opt/ipr/work.local/bob/flobe_build$ cmake -DCMAKE_INSTALL_PREFIX=~/flobe_install -DUSE_MAGMA=FALSE -DBLA_VENDOR=Intel10_64lp ../flobe
```

if you want to configure to use magma:
```sh
bob@physix90:/opt/ipr/work.local/bob/flobe_build$ cmake -DCMAKE_INSTALL_PREFIX=~/flobe_install -DUSE_MAGMA=TRUE -DMAGMA_API=CPU_MEM_API ../flobe
```
then build an install the executables
```sh
bob@physix90:/opt/ipr/work.local/bob/flobe_build$ make install
```


run the benchmark, here with `num_cores=72` `num_runs=10`

```sh
bob@physix90:/opt/ipr/work.local/bob/flobe$ cd ./benchmarks/madiag1
bob@physix90:/opt/ipr/work.local/bob/flobe/benchmarks/madiag1$ ./madiag_batch.bash 72 physix90 10
# machine_id matrix_size num_loops num_processes duration_min(s) duration_max(s)
physix90 128 10000 1 12.38 12.38
physix90 128 10000 2 12.38 12.45
...
physix90 8192 1 71 123.45 234.56
physix90 8192 1 72 123.45 234.56
```

once the execution of all benchmarks is complete, the results are store in `<machine-id>_<process-id>.tsv`
```sh
bob@physix90:/opt/ipr/work.local/bob/flobe/benchmarks/madiag1$ cat physix90_368318.tsv
#machine_id     matrix_size     num_loops       num_processes   duration_min(s) duration_max(s)
physix90        128     10000   1       12.38   12.38
physix90        128     10000   2       12.38   12.45
physix90        128     10000   3       12.57   12.62
...
physix90        128     10000   70      14.67   20.25
physix90        128     10000   71      14.42   19.52
physix90        128     10000   72      14.80   20.28
physix90        256     2000    1       11.98   11.98
physix90        256     2000    2       11.96   12.33
physix90        256     2000    3       12.08   12.12
...
...
physix90        256     2000    70      14.07   19.63
physix90        256     2000    71      14.53   18.22
physix90        256     2000    72      14.09   19.40
physix90        512     300     1       11.74   11.74
physix90        512     300     2       11.91   11.96
physix90        512     300     3       11.94   11.97
...
physix90        4096    1       70      585.42  824.14
physix90        4096    1       71      221.22  316.54
physix90        4096    1       72      331.30  650.13
physix90        8192    1       1       335.81  335.81
physix90        8192    1       2       203.35  332.07
physix90        8192    1       3       202.12  328.06
...
physix90        8192    1       64      1370.33 2175.10
physix90        8192    1       65      1384.70 1709.24
physix90        8192    1       66      2517.80 4155.65
```

### 2. gather results and plot them on a workstation


```sh
bob@bob-ws:~/work$ git clone git@github.com:g-raffy/flobe.git
bob@bob-ws:~/work$ cd ./flobe/benchmarks/madiag1
```

optional: gather all the log files that were generated (one per run of `madiag1_star.bash`).
```sh
bob@bob-ws:~/work/flobe/benchmarks/madiag1$ rsync --rsh=ssh -va bob@physix90.ipr.univ-rennes1.fr:/opt/ipr/work.local/bob/flobe/benchmarks/madiag1/cohab* ./measurements/fr.univ-rennes1.ipr.physix/physix90/
```

retreive the tsv file containing the results of the benchmark `physix90_368318.tsv` (`368318` here is the id of the process, so it can be anything)
```sh
bob@bob-ws:~/work/flobe/benchmarks/madiag1$ rsync --rsh=ssh -va bob@physix90.ipr.univ-rennes1.fr:/opt/ipr/work.local/bob/flobe/benchmarks/madiag1/physix90* ./measurements/fr.univ-rennes1.ipr.physix/physix90/
```

generate plots `graphics/*.svg`
```sh
bob@bob-ws:~/work/flobe/benchmarks/madiag1$ python ./graphics.py
```









