#!/usr/bin/env python3
import numpy
import matplotlib.pyplot as plt
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
        fig.savefig(self._out_svg_dir_path + '/' + fig.axes[0].get_title() + '.svg')

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
    for tsv_file in ['cohab/physix90/physix90_368318.tsv', 'cohab/physix12/physix12_9079.tsv', 'cohab/physix96/physix96_189748.tsv']:
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
        for machine_id in ['physix12', 'physix90', 'physix96']:
            # machine_id matrix_size num_loops num_processes duration_min(s) duration_max(s)
            # physix12 128 10000 1 23.33 23.33
            # physix12_table = numpy.genfromtxt('cohab/physix12/physix12_9079.tsv', dtype=("|U10", int, int, int, float, float), names=True, delimiter='\t')
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


if __name__ == '__main__':
    figure_handler = SvgFigureHandler(out_svg_dir_path='./graphics')
    # figure_handler = ScreenFigureHandler()
    figures = plot_load_to_diag_perf()
    for fig in figures:
        figure_handler.on_figure_ended(fig)
    figure_handler.on_finalize()
