#!/usr/bin/env python3
from typing import List, Dict, Any
from pathlib import Path
import subprocess
import argparse
import os
import re
from datetime import datetime
from shutil import rmtree

SoftwareId = str  # eg 'ifort', 'gfortran', 'mkl', 'openblas'
Version = str  # eg 19.1.0


class SoftwareDef():
    product_name: SoftwareId
    product_version: Version

    def __init__(self, product_name: SoftwareId, product_version: Version):
        self.product_name = product_name
        self.product_version = product_version

    def as_id(self) -> str:
        return f'{self.product_name}-{self.product_version}'


class Config():
    fortran_compiler: SoftwareDef
    blas_lib: SoftwareDef

    def __init__(self, fortran_compiler: SoftwareDef, blas_lib: SoftwareDef):
        self.fortran_compiler = fortran_compiler
        self.blas_lib = blas_lib


class TsvFile():
    columns: List[str]
    tsv_file_path: Path

    def __init__(self, tsv_file_path: Path, columns: List[str]):
        self.columns = columns
        self.tsv_file_path = tsv_file_path
        parent_dir = self.tsv_file_path.parent
        if not parent_dir.exists():
            os.makedirs(parent_dir, exist_ok=True)

    def add_record(self, cols: Dict[str, Any]):
        assert len(cols) == len(self.columns), f'the record you attempt to add has {len(cols)} fields {cols} while {self.tsv_file_path} is configured with {len(self.columns)} fields {self.columns}'
        # print(type(self.columns))
        add_header = not self.tsv_file_path.exists()
        col_values = [str(cols[key]) for key in self.columns]
        line = '\t'.join(col_values)
        with open(self.tsv_file_path, 'at', encoding='utf-8') as tsv_file:
            if add_header:
                header = '\t'.join(self.columns)
                tsv_file.write(f'#{header}\n')
            tsv_file.write(line + '\n')


def main():

    host_fqdn = 'alambix98.ipr.univ-rennes.fr'  # todo: remove hardcoded value
    cpu_id = 'intel-xeon-gold-6248r'  # todo: remove hardcoded value
    num_cores = 48  # todo: remove hardcoded value

    configs: List[Config] = []
    gfortran_compiler = SoftwareDef(product_name='gfortran', product_version='12.2.0')   # the gfortran version on alambix98 # todo: remove hardcoded value

    blas_lib = SoftwareDef(product_name='mkl', product_version='2023.0.0')
    configs.append(Config(fortran_compiler=gfortran_compiler, blas_lib=blas_lib))

    # blas_lib = SoftwareDef(product_name='openblas', product_version='0.3.21')   # the default blas on physix98 # todo: remove hardcoded value
    # configs.append(Config(fortran_compiler=gfortran_compiler, blas_lib=blas_lib))

    arg_parser = argparse.ArgumentParser(description='measures the performance of madiag2 in various configurations')
    arg_parser.add_argument('--flobe-root-dir', required=True)

    args = arg_parser.parse_args()

    flobe_root_dir = Path(args.flobe_root_dir)

    build_dir = Path('/tmp/flobe/madiag2')

    cols_id = [
        'measurement_time',
        'host_fqdn',
        'cpu_id',
        'num_cpu',
        'compiler_id',
        'blas_id',
        'max_cores',
        'matrix_size',
        'num_loops',
        'with_eigen_vectors',
        'duration',
        'cpu_time'
    ]
    tsv_file = TsvFile(flobe_root_dir / f'benchmarks/madiag2/measurements/{host_fqdn}/runbench.tsv', cols_id)

    max_cores = num_cores

    with_eigen_vectors = False
    for config in configs:
        if build_dir.exists():
            rmtree(build_dir)
        os.makedirs(build_dir, exist_ok=True)
        prefix = None
        cmake_defs = []
        cmake_defs.append('-DUSE_MAGMA=FALSE')
        if config.blas_lib.product_name == 'mkl':
            prefix = f'. /etc/profile; module load lib/mkl/{config.blas_lib.product_version}'
            cmake_defs.append('-DBLA_VENDOR=Intel10_64lp')
        cmake_command = f'{prefix + "; " if prefix is not None else ""}cmake {" ".join(cmake_defs)} {flobe_root_dir}/benchmarks/madiag2'
        print(cmake_command)
        subprocess.run(cmake_command, shell=True, cwd=build_dir, check=True, capture_output=True)
        subprocess.run(['make'], cwd=build_dir, check=True, capture_output=True)
        num_loops = None
        for matrix_size in [10000, 14000, 24000, 32768]:  # [128, 256, 512, 1024, 2048, 4096, 8192]:
            for with_eigen_vectors in [False, True]:
                for max_cores in [48]:  # [1, 2, 4, 8, 16, 32, num_cores]:  # todo: replace hardcoded values
                    # adjust num_loops so that the bench takes about 12 seconds
                    num_loops_for_all_cores = {
                        64: 100,  # just for test
                        128: 10000,
                        256: 2000,
                        512: 300,
                        1024: 40,
                        2048: 5,
                        4096: 1,
                        8192: 1,
                        10000: 1,
                        12000: 1,
                        14000: 1,
                        16384: 1,
                        24000: 1,
                        32768: 1,
                    }[matrix_size]
                    num_loops = max(1, int(float(num_loops_for_all_cores) * max_cores / num_cores))
                    print(f'num_loops_for_all_cores: {num_loops_for_all_cores}, max_cores: {max_cores}, num_loops: {num_loops}')
                    run_command = f'{prefix + "; " if prefix is not None else ""}OMP_NUM_THREADS={max_cores} MKL_NUM_THREADS={max_cores} ./madiag2 {matrix_size} {num_loops} {"true" if with_eigen_vectors else "false"}'
                    print(f'run_command = {run_command}')
                    completed_process = subprocess.run(run_command, cwd=build_dir, check=True, shell=True, capture_output=True)

                    duration = None
                    cpu_time = None
                    for line in completed_process.stdout.decode('utf-8').split('\n'):
                        # print(line)
                        match = re.match(r'\s*- wallclock time:\s*(?P<duration>[0-9.]+)', line)
                        if match:
                            duration = float(match['duration'])

                        match = re.match(r'\s*- cpu_time:\s*(?P<cpu_time>[0-9.]+)', line)
                        if match:
                            cpu_time = float(match['cpu_time'])

                    print(f'duration : {duration} seconds, cpu time: {cpu_time}')
                    record = {
                        'measurement_time': datetime.now().isoformat(),
                        'host_fqdn': host_fqdn,
                        'cpu_id': cpu_id,
                        'num_cpu': 2,
                        'compiler_id': config.fortran_compiler.as_id(),
                        'blas_id': config.blas_lib.as_id(),
                        'max_cores': max_cores,
                        'matrix_size': matrix_size,
                        'with_eigen_vectors': with_eigen_vectors,
                        'num_loops': num_loops,
                        'duration': duration,
                        'cpu_time': cpu_time}
                    tsv_file.add_record(record)


main()
