#!/bin/bash
# performs a benchmark of matrix diagonalization by running multiple matrix diagonlization code simultaneously (one per core)
# launches the given diagonalization executable one time on each core and measures the time it takes
# each run of the executable executes on one core, all runs are run simultaneously
# the results of the benchmark are in
# - <stdout> : summary on the benchmark, in human readable form
# - $LOG_FILE : stores the mixed <stdout> of the $NUM_CORES matrix diagonalization code running simultaneously
set -o errexit
NUM_CORES="$1" # the number of instances of STRESSER_COMMAND that will be run
MATRIX_SIZE="$2" # eg 128
NUM_LOOPS="$3" # eg 100
EXE_NAME="$4" # variant of the diagonalization program to allow to choose which blas library to use, eg 'diag_mkl_omp'

make $EXE_NAME  # make sure that the executable is built
case $EXE_NAME in
	'diag_mkl_omp' | 'diag_mkl')
		module load compilers/ifort/latest
		;;
esac

STRESSER_COMMAND="./$EXE_NAME $MATRIX_SIZE $NUM_LOOPS"

LOG_FILE="./cohab_$(hostname)_${NUM_CORES}_${MATRIX_SIZE}_${NUM_LOOPS}.log"
date > "$LOG_FILE"
# launch the executable $NUM_CORES times
for c in $(seq 1 $NUM_CORES)
do
	echo "num cores : $c" >> "$LOG_FILE"
	OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 $STRESSER_COMMAND >> "$LOG_FILE" &  # launches the executable in the background
done
wait  # wait for all executables to complete
date >> "$LOG_FILE"


OUT_DURATIONS_FILE_PATH="./durations.txt"  # temporary file that stores the duration of each run in seconds
MATRIX_FOOTPRINT=$(echo dummy | awk "{printf(\"%d\", $MATRIX_SIZE*$MATRIX_SIZE*8)}" )  # size of a matrix in bytes
PROCESS_FOOTPRINT=$(echo dummy | awk "{printf(\"%d\", $MATRIX_FOOTPRINT*3)}" )  # size of all the matrices usied in the diaginalization (graffy: where does the 3 come from ?), in bytes
COMPUTER_FOOTPRINT=$(echo dummy | awk "{printf(\"%d\", $PROCESS_FOOTPRINT*$NUM_CORES)}" )
strings $LOG_FILE | awk '/^Time taken/ {print $(NF-1)}' | sort -n > $OUT_DURATIONS_FILE_PATH  # the executable is expected to output a line such as "Time taken by dsyev for matrix size 256 was 1.23456 seconds"
MIN_DURATION=$(head -1 $OUT_DURATIONS_FILE_PATH)
MAX_DURATION=$(tail -1 $OUT_DURATIONS_FILE_PATH)

# output a summary of the benchmark in stdout
echo "number of cores used (1 process per core) : $NUM_CORES"
echo "size of the matrix to be diagonalized : $MATRIX_SIZE x $MATRIX_SIZE ($(numfmt --to=iec-i --suffix=B $MATRIX_FOOTPRINT))"
echo "number of identical diagonalizations performed on each core : $NUM_LOOPS"
echo "memory used : $(numfmt --to=iec-i --suffix=B $COMPUTER_FOOTPRINT) ($(numfmt --to=iec-i --suffix=B $PROCESS_FOOTPRINT) per core)"
echo "computation duration range on each core [$MIN_DURATION;$MAX_DURATION] seconds"
