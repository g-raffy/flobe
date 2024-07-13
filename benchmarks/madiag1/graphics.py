#!/usr/bin/env python3
import numpy
import matplotlib.pyplot as plt
from pathlib import Path
import abc


class IFigureHandler(object):
    """
    specifies what to do with generated figures
    """

    @abc.abstractmethod
    def on_figure_ended(self, fig):
        """
        :param matplotlib.Figure fig:
        """
        pass

    @abc.abstractmethod
    def on_finalize(self):
        """
        called after all figures have been created
        """
        pass


class ScreenFigureHandler(IFigureHandler):
    """
    displays figures on screen
    """
    def __init__(self):
        pass

    def on_figure_ended(self, fig):
        pass

    def on_finalize(self):
        plt.show()


class SvgFigureHandler(IFigureHandler):
    """
    saves figures as svg files
    """

    def __init__(self, out_svg_dir_path):
        """
        :param str out_svg_dir_path: where to save the svg files
        """
        self._out_svg_dir_path = out_svg_dir_path

    def on_figure_ended(self, fig):
        file_name = None
        if hasattr(fig, 'file_name'):
            file_name = fig.file_name
        else:
            file_name = fig.axes[0].get_title()
        fig.savefig(self._out_svg_dir_path + '/' + file_name + '.svg')

    def on_finalize(self):
        pass


def merge_structured_arrays(array1, array2):
    n1 = len(array1)
    n2 = len(array2)
    array_out = array1.copy()
    array_out.resize(n1 + n2)
    array_out[n1:] = array2
    return array_out


def get_diag_perf_table():
    diag_perf_table = None
    for tsv_file in ['measurements/fr.univ-rennes1.ipr.physix/physix90/physix90_368318.tsv', 'measurements/fr.univ-rennes1.ipr.physix/physix12/physix12_9079.tsv', 'measurements/fr.univ-rennes1.ipr.physix/physix96/physix96_189748.tsv']:
        table = numpy.genfromtxt(tsv_file, dtype=("|U10", int, int, int, float, float), names=True, delimiter='\t')
        # table is what numpy calls a structured array
        if diag_perf_table is None:
            diag_perf_table = table
        else:
            diag_perf_table = merge_structured_arrays(diag_perf_table, table)
    # print(diag_perf_table.dtype)
    # print(diag_perf_table[0])
    # print(diag_perf_table[-1])
    
    return diag_perf_table


def plot_load_to_diag_perf():

    machine_to_color = {'physix12': 'red',
                        'physix90': 'blue',
                        'physix96': 'green'}

    machine_to_num_cores = {'physix12': 32,
                            'physix90': 72,
                            'physix96': 48}

    diag_perf_table = get_diag_perf_table()
    figures = []
    show_for_one_matrix = True
    mat_sizes = numpy.unique(diag_perf_table['matrix_size'])
    for matrix_size in mat_sizes:

        fig, ax = plt.subplots()
        max_duration = 0.0
        for machine_id in ['physix12', 'physix96']:
            # machine_id matrix_size num_loops num_processes duration_min(s) duration_max(s)
            # physix12 128 10000 1 23.33 23.33
            # physix12_table = numpy.genfromtxt('measurements/fr.univ-rennes1.ipr.physix/madiag1/physix12/physix12_9079.tsv', dtype=("|U10", int, int, int, float, float), names=True, delimiter='\t')
            machine_table = diag_perf_table[diag_perf_table['machine_id'] == machine_id]
            # print(machine_table.dtype)
            # print(machine_table.shape)
            perfs_for_size = machine_table[machine_table['matrix_size'] == matrix_size]
            num_diags_per_core = perfs_for_size['num_loops'][0]
            # print(perfs_for_size)
            x = perfs_for_size['num_processes']
            y = perfs_for_size['duration_maxs']
            if show_for_one_matrix:
                # scale = 1.0 / float(num_diags_per_core)
                y *= 1.0 / float(num_diags_per_core)
                num_diags_per_core = 1
            a = x * 1.0 / float(machine_to_num_cores[machine_id])
            x = a
            # print(machine_to_color[machine_id])
            print(matrix_size, machine_id, max(y))
            max_duration = max(max_duration, max(y))
            ax.scatter(x, y, color=machine_to_color[machine_id], s=2.0, label=machine_id)
        ax.set_ylim(ymin=0.0, ymax=max_duration)
        plt.xlabel(u'machine load ratio')
        plt.ylabel(u'max duration of the %d diagonalizations (s)' % num_diags_per_core)
        plt.title(u'%d matrix diagonalization %d x %d per core' % (num_diags_per_core, matrix_size, matrix_size))
        # plt.legend(bbox_to_anchor=(0.2, 1.0))
        ax.legend()
        figures.append(fig)
    return figures


def create_matsize_to_diag_perf_plot() -> plt.Figure:
    # config	num_cores	num_gpu	matrix_size	cpu(s)	gpu_mem_stage1(MiB)	gpu_mem_stage2(MiB)
    table = numpy.genfromtxt(Path('measurements/fr.univ-rennes1.ipr.physix/diag_bench.tsv'), dtype=("|U10", int, int, int, float, float), names=True, delimiter='\t', usecols=(0, 1, 2, 3, 4))
    print(table.dtype.names)

    # #measurement_time	host_fqdn	cpu_id	num_cpu	compiler_id	blas_id	max_cores	matrix_size	num_loops	with_eigen_vectors	duration	cpu_time
    table2 = numpy.genfromtxt(Path('../madiag2/measurements/alambix98.ipr.univ-rennes.fr/runbench.tsv'), dtype=("|U30", "|U30", "|U30", int, "|U30", "|U30", int, int, int, bool, float, float), names=True, delimiter='\t', usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11))
    print(table2.dtype.names)

    fig, ax = plt.subplots()
    gpu_table = table[table['config'] == 'magma']
    cpu_table = table[table['config'] == 'mkl_omp']

    # print(table2)
    cpu2_table = table2[table2['blas_id'] == 'mkl-2023.0.0']
    cpu2_table = cpu2_table[cpu2_table['max_cores'] == 48]
    cpu2_table = cpu2_table[cpu2_table['with_eigen_vectors'] == True]  # noqa
    print(cpu2_table)

    # cpu2_table = table2[table2['blas_id'] == 'mkl-2023.0.0' and table2['with_eigen_vectors'] == True and table2['max_cores'] == 48]
    ax.plot(gpu_table['matrix_size'], gpu_table['cpus'], 'o-', label='magmaf_dsyevd_gpu on nvidia p100')
    ax.plot(cpu_table['matrix_size'], cpu_table['cpus'], 'o-', label='mkl 2020.0.1 dsyevd on 4 18-core xeon 6154')
    ax.plot(cpu2_table['matrix_size'], cpu2_table['duration'] / cpu2_table['num_loops'], 'o-', label='mkl 2023.0.0 dsyevd on 2 24-core xeon 6248r')
    # ax.scatter(x=gpu_table['matrix_size'], y=gpu_table['cpus'], label='magmaf_dsyevd_gpu on nvidia p100')
    # ax.scatter(x=cpu_table['matrix_size'], y=cpu_table['cpus'], label='mkl 2020.0.1 dsyevd on 4 18-core xeon 6154')
    plt.xlabel('n: matrix size (n*n matrix)')
    plt.ylabel('duration of the dsyevd computation (s)')
    plt.title('performance of dsyevd with eigenvectors')
    ax.legend()
    return fig


def interpolate(req_x: float, vx: numpy.array, vy: numpy.array) -> float:
    '''
    req_x: requested x
    '''
    # print(f'req_x={req_x}')
    assert len(vx) == len(vy)
    xl = vx[0]  # left
    yl = vy[0]
    xr = vx[-1]  # right
    yr = vy[-1]
    for index in range(len(vx)):
        xx = vx[index]
        yy = vy[index]
        # print(f'xx={xx}, xl={xl}, xr={xr}')
        if xx < req_x and xx > xl:
            xl = xx
            yl = yy
        if xx > req_x and xx < xr:
            xr = xx
            yr = yy
    # print(f'xl={xl}, xr= {xr}')
    if xl > req_x or xr < req_x:
        raise IndexError(f'failed to find the surrounding elements for x = {req_x}')
    ratio = (req_x - xl) / (xr - xl)
    req_y = yl + (yr - yl) * ratio
    # print(f'ratio={ratio}, yl= {yl}, yr= {yr}, req_y= {req_y}')
    return req_y


def create_matsize_to_diag_gpu_gain_plot() -> plt.Figure:
    # config	num_cores	num_gpu	matrix_size	cpu(s)	gpu_mem_stage1(MiB)	gpu_mem_stage2(MiB)
    table = numpy.genfromtxt(Path('measurements/fr.univ-rennes1.ipr.physix/diag_bench.tsv'), dtype=("|U10", int, int, int, float, float), names=True, delimiter='\t', usecols=(0, 1, 2, 3, 4))
    print(table.dtype.names)

    # #measurement_time	host_fqdn	cpu_id	num_cpu	compiler_id	blas_id	max_cores	matrix_size	num_loops	with_eigen_vectors	duration	cpu_time
    table2 = numpy.genfromtxt(Path('../madiag2/measurements/alambix98.ipr.univ-rennes.fr/runbench.tsv'), dtype=("|U30", "|U30", "|U30", int, "|U30", "|U30", int, int, int, bool, float, float), names=True, delimiter='\t', usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11))
    print(table2.dtype.names)

    fig, ax = plt.subplots()
    gpu_table = table[table['config'] == 'magma']

    # print(table2)
    cpu2_table = table2[table2['blas_id'] == 'mkl-2023.0.0']
    cpu2_table = cpu2_table[cpu2_table['max_cores'] == 48]
    cpu2_table = cpu2_table[cpu2_table['with_eigen_vectors'] == True]  # noqa
    print(cpu2_table)

    gpu_matrix_sizes = gpu_table['matrix_size']
    gpu_dsyevd_durations = gpu_table['cpus']
    cpu_matrix_size = cpu2_table['matrix_size']
    cpu_dsyevd_durations = cpu2_table['duration'] / cpu2_table['num_loops']

    mat_sizes = []
    gpu_gains = []   # numpy.zeros(shape=(len(gpu_matrix_sizes),))
    print(gpu_gains)
    for gpu_measur_index in range(len(gpu_matrix_sizes)):
        mat_size = gpu_matrix_sizes[gpu_measur_index]
        if mat_size > 20000:  # ignore matrix sizes above 20000 (request by flique)
            continue  
        gpu_dsyevd_duration = gpu_dsyevd_durations[gpu_measur_index]
        try:
            cpu_dsyevd_duration = interpolate(mat_size, cpu_matrix_size, cpu_dsyevd_durations)
        except IndexError:
            continue
        mat_sizes.append(mat_size)
        gpu_gains.append(cpu_dsyevd_duration / gpu_dsyevd_duration)

    # ax.plot(mat_sizes, gpu_gains, 'o-'), label='gpu nvidia P100 speedup')
    ax.plot(mat_sizes, gpu_gains, 'o-')
    plt.xlabel('n: matrix size (n*n matrix)')
    plt.ylabel('computational speedup')
    # plt.ylabel('dsyevd computation speedup\ncompared to dual cpu intel xeon 6248r')
    # plt.title('gpu dsyevd speedup over cpu')
    fig.file_name = 'gpu-speedup'
    # ax.legend()
    return fig


if __name__ == '__main__':
    figure_handler = SvgFigureHandler(out_svg_dir_path='./graphics')
    # figure_handler = ScreenFigureHandler()
    # figures = plot_load_to_diag_perf()
    # for fig in figures:
    #     figure_handler.on_figure_ended(fig)

    fig = create_matsize_to_diag_gpu_gain_plot()
    figure_handler.on_figure_ended(fig)
    figure_handler.on_finalize()
