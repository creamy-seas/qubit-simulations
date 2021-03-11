# def sort_array(self, array_to_sort, column_to_sort_by):
#     """
#     __ Parameters __
#     array_to_sort: array to perform sorting for
#     column_to_sort_by: which column to use for sorting

#     __ Description __
#     sorts the supplied array
#     """
#     # 1 - sort array by the i-th column
#     sorted_index = np.argsort(array_to_sort[column_to_sort_by])
#     for i in range(0, len(array_to_sort)):
#         array_to_sort[i] = np.array(array_to_sort[i])[sorted_index]
#     return array_to_sort

# def import_transmission_spectrum(self, file_name, plot_axes, plot_list):
#     """Fiel should have the format

#     # no_of_x_col by no_of_y_col
#     xCol(field) yCol(freq) zCol(tranmission)

#     __ Parameters __

#     plot_list: list of the field points to plot

#     *** Format is ***
#     # no-of-xCol by no-of-yCol
#     xCol(field) yCol(freq) zCol(transmission)

#     __ Description __
#     Import the transmission array, supplied as a commented file. The comment
#     must specify the number of field points used
#     """

#     logging.debug("Importing transmission data file '{file_name}'")

#     with open(file_name) as fin:
#         first_line = fin.readline()

#     no_field_points = int(first_line.split()[1])
#     no_freq_points = int(first_line.split()[3])
#     logging.debug(
#         f"{no_field_points} field points and {no_freq_points} frequency points"
#     )

#     transmission_array = np.vsplit(np.loadtxt(file_name), no_field_points)
#     for i in range(0, len(transmission_array)):
#         transmission_array[i] = transmission_array[i].transpose()

#     if self.plot_or_not:
#         logging.debug("Plotting, cos why not")
#         for i in plot_list:
#             if (i < 0) or (i >= len(transmission_array)):
#                 raise ValueError(
#                     "{i} is outside the allowed field values [0:{no_field_points}]"
#                 )

#             plot_axes.plot(transmission_array[i][1], transmission_array[i][2])
#     logging.debug("==> 'import_transmission_spectrum' finished")

# def save_plot(self, mpl_axes, filename, dpi=100):
#     mpl_axes.set_facecolor("white")
#     plt.savefig(f"output/{filename}.png", dpi=dpi)
