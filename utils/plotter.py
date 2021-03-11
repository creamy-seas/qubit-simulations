from typing import Union

import numpy as np
import matplotlib.pyplot as plt


def plot_column_data(mpl_axes, output_file_name):

    array_in = np.loadtxt(output_file_name).transpose()

    mpl_axes.plot(array_in[0], array_in[1])
    plt.show()


def save_onto_white_background(
    mpl_axes, output_file_name, save_format: Union["svg", "png", "pdf"]
):
    if save_format != "svg":
        mpl_axes.set_facecolor("white")
    plt.savefig(f"{output_file_name}.{save_format}", transparent=True)


def set_geometry(left=0, c2=30, width=1280, height=1600):
    plt.get_current_fig_manager().window.setGeometry(left, c2, width, height)


def customize_figure_size(left, right, bottom, top):
    """ Set the figure size

    0 < left < right < 1
    0 < bottom < top < 1

    -1 for no change
    """

    if left != -1:
        plt.rcParams["figure.subplot.left"] = left
    if right != -1:
        plt.rcParams["figure.subplot.right"] = right
    if top != -1:
        plt.rcParams["figure.subplot.top"] = top
    if bottom != -1:
        plt.rcParams["figure.subplot.bottom"] = bottom


def compare_to_matlab(coo_matrix, output_file_name):
    """
    __ Description __
    outputs the non zeros row, col, element values of a coo matrix
    used to compare with the Matlab simulations
    """

    row_array = []
    col_array = []
    element_real_array = []
    element_imaginary_array = []

    # 1 - collect element
    for row, col, elm in zip(coo_matrix.row, coo_matrix.col, coo_matrix.data):
        row_array.append(row)
        col_array.append(col)
        element_real_array.append(elm.real)
        element_imaginary_array.append(elm.imag)

    row_array = np.array(row_array)
    col_array = np.array(col_array)
    element_real_array = np.array(element_real_array)
    element_imaginary_array = np.array(element_imaginary_array)

    # 2 - sort array by the column
    _sort_idx = np.argsort(col_array)
    row_array = np.array(row_array)[_sort_idx]
    col_array = np.array(col_array)[_sort_idx]
    element_real_array = np.array(element_real_array)[_sort_idx]
    element_imaginary_array = np.array(element_imaginary_array)[_sort_idx]

    with open(output_file_name, "w") as fout:
        for row, col, element_real, element_imaginary in zip(
            row_array, col_array, element_real_array, element_imaginary_array
        ):
            string_to_write = "{_row}\t{_col}\t{_real:.1f} {_im:.2f}\n".format(
                _row=row + 1, _col=col + 1, _real=element_real, _im=element_imaginary
            )

        fout.write(string_to_write)
