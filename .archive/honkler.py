import numpy as np
import matplotlib.pyplot as plt


def save_ree(axes, file_name, save_format):
    """
    __ Parameters __
    axes: axes to save
    save_format: png, svg, pdf

    __ Description __
    save graphic with white background
    """
    if(save_format != "svg"):
        axes.set_facecolor("white")
    plt.savefig("%s.%s" % (file_name, save_format), transparent=True)


def config_plot_size(left, right, bottom, top):
    """
    __ Parameters__
    0 < left < right < 1
    0 < bottom < top < 1

    -1 to keep the same

    __ Description __
    sets figure size
    """

    if (left != -1):
        plt.rcParams["figure.subplot.left"] = left
    if (right != -1):
        plt.rcParams["figure.subplot.right"] = right
    if (top != -1):
        plt.rcParams["figure.subplot.top"] = top
    if (bottom != -1):
        plt.rcParams["figure.subplot.bottom"] = bottom


def write_comparisson_matrix(coo_matrix, fileName):
    """
    __ Parameters __
    coomatrix - elements to print
    fileName - file to write them to

    __ Description __
    outputs the non zeros row, col, element values of a coo matrix
    used to compare with the Matlab simulations
    """

    row_array = []
    col_array = []
    elm_array = []
    elI_array = []

    # 1 - collect element
    for row, col, elm in zip(coo_matrix.row, coo_matrix.col, coo_matrix.data):
        row_array.append(row)
        col_array.append(col)
        elm_array.append(elm.real)
        elI_array.append(elm.imag)

    row_array = np.array(row_array)
    col_array = np.array(col_array)
    elm_array = np.array(elm_array)
    elI_array = np.array(elI_array)

    # 2 - sort array by the column
    sorted_index = np.argsort(col_array)
    row_array = np.array(row_array)[sorted_index]
    col_array = np.array(col_array)[sorted_index]
    elm_array = np.array(elm_array)[sorted_index]
    elI_array = np.array(elI_array)[sorted_index]

    # 3 - write
    fout = open(fileName, 'w')
    for row, col, elm, elI in zip(row_array, col_array, elm_array, elI_array):
        string_to_write = "%i\t%i\t%.1f  %.2f\n" % (row + 1, col + 1, elm, elI)
        fout.write(string_to_write)
    fout.close()


def plot_column_data(ax, fileName):
    """
    __ Parameters __
    fileName

    __ Description __
    Plots an x-y joined graph
    """

    ax_type = str(type(ax)).split('\'')[1]
    # if(type(ax) == "numpy.ndarray"):
    # if ()

    array_in = np.loadtxt(fileName).transpose()
    print(len(array_in[0]))

    ax.plot(array_in[0], array_in[1])
    plt.show()
