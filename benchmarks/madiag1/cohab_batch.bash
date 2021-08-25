#!/bin/bash
# generates tsv file by calling bench_cohabitation.bash with multiple parameters
NUM_CORES="$1" # the number of physical cores of the machine running this test
MACHINE_ID="$2" # the name of the machine running this test (eg physix12)
NUM_RUNS=$3
set -o errexit

case $EXE_NAME in
	'diag_mkl_omp' | 'diag_mkl')
		module load compilers/ifort/latest
		;;
esac

tmp_dir='/tmp/bug3177/cohab_batch/$$'
mkdir -p $tmp_dir


tsv_file_path="./${MACHINE_ID}_$$.tsv"

function add_tsv_line()
{
	line="$1"
	tsv_file_path="$2"

	echo $line
	echo $line >> $tsv_file_path
}

tsv_header_line=$(printf "# %s\t%s\t%s\t%s\t%s\t%s" 'machine_id' 'matrix_size' 'num_loops' 'num_processes' 'duration_min(s)' 'duration_max(s)')
add_tsv_line "$tsv_header_line" "$tsv_file_path"


for run_index in $(seq 1 $NUM_RUNS)
do
	num_loops=''
	for matrix_size in 128 256 512 1024 2048 4096 8192
	do
		# adjust num_loops so that the bench takes about 12 seconds
		case "$matrix_size" in
			64) # just for test
				num_loops=100
				;;
			128)
				num_loops=10000
				;;
			256)
				num_loops=2000
				;;
			512)
				num_loops=300
				;;
			1024)
				num_loops=40
				;;
			2048)
				num_loops=5
				;;
			4096)
				num_loops=1
				;;
			8192)
				num_loops=1
				;;
			*)
				echo "error : unhandled matrix size ($matrix_size)"
				exit 1
				;;
		esac
		for num_processes in $(seq 1 $NUM_CORES)
		do
			./bench_cohabitation.bash $num_processes $matrix_size $num_loops diag_mkl_omp > $tmp_dir/measure.log
			durations=$(cat $tmp_dir/measure.log | grep '^computation duration' | sed 's/computation duration range on each core \[\([0-9.]*\);\([0-9.]*\)\] seconds/\1 \2/')
			min_duration=$(echo $durations | awk '{printf("%s",$1);}')
			max_duration=$(echo $durations | awk '{printf("%s",$2);}')
			tsv_line=$(printf "%s\t%s\t%s\t%s\t%s\t%s" "$MACHINE_ID" "$matrix_size" "$num_loops" "$num_processes" "$min_duration" "$max_duration")
			add_tsv_line "$tsv_line" "$tsv_file_path"
		done
	done
done