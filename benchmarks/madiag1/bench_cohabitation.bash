#!/bin/bash
set -o errexit
NUM_CORES="$1" # the nember of instances of STRESSER_COMMAND that 
MATRIX_SIZE="$2" # eg 128
NUM_LOOPS="$3" # eg 100
EXE_NAME="$4" # eg 'diag_mkl_omp'

make $EXE_NAME
case $EXE_NAME in
	'diag_mkl_omp' | 'diag_mkl')
		module load compilers/ifort/latest
		;;
esac

STRESSER_COMMAND="./$EXE_NAME $MATRIX_SIZE $NUM_LOOPS"

LOG_FILE="./cohab_$(hostname)_${NUM_CORES}_${MATRIX_SIZE}_${NUM_LOOPS}.log"
date > "$LOG_FILE"
for c in $(seq 1 $NUM_CORES)
do
	echo "num cores : $c" >> "$LOG_FILE"
	OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 $STRESSER_COMMAND >> "$LOG_FILE" &
done
wait
date >> "$LOG_FILE"

MATRIX_FOOTPRINT=$(echo dummy | awk "{printf(\"%d\", $MATRIX_SIZE*$MATRIX_SIZE*8)}" )
PROCESS_FOOTPRINT=$(echo dummy | awk "{printf(\"%d\", $MATRIX_FOOTPRINT*3)}" )
COMPUTER_FOOTPRINT=$(echo dummy | awk "{printf(\"%d\", $PROCESS_FOOTPRINT*$NUM_CORES)}" )
strings $LOG_FILE | awk '/^Time taken/ {print $(NF-1)}' | sort -n > ./durations.txt
MIN_DURATION=$(head -1 ./durations.txt)
MAX_DURATION=$(tail -1 ./durations.txt)
echo "number of cores used (1 process per core) : $NUM_CORES"
echo "size of the matrix to be diagonalized : $MATRIX_SIZE x $MATRIX_SIZE ($(numfmt --to=iec-i --suffix=B $MATRIX_FOOTPRINT))"
echo "number of identical diagonalizations performed on each core : $NUM_LOOPS"
echo "memory used : $(numfmt --to=iec-i --suffix=B $COMPUTER_FOOTPRINT) ($(numfmt --to=iec-i --suffix=B $PROCESS_FOOTPRINT) per core)"
echo "computation duration range on each core [$MIN_DURATION;$MAX_DURATION] seconds"
