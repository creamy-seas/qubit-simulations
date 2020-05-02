from PIL import Image
import numpy as np
from scipy.optimize import curve_fit


def parabolla(x, a, b, c):
    return a * x ** 2 + b * x + c


def fit_parabolla(fitting_array):
    popt, pcov = curve_fit(parabolla, fitting_array[0], fitting_array[1])
    return popt, pcov


def fit_and_plot_parabolla(fitting_array, plotting_axis_instance):
    fitting, _ = fit_parabolla(fitting_array)
    x = np.linspace(min(fitting_array[0]), max(fitting_array[0]), 1000)
    y = parabolla(x, fitting[0], fitting[1], fitting[2])

    plotting_axis_instance.plot(x, y)
    plotting_axis_instance.scatter(
        fitting_array[0], fitting_array[1], color="C3", marker="."
    )
    plotting_axis_instance.set_title(
        "%.1f" % (fitting[0] * 2), fontdict={"fontsize": 12}
    )


def delete_entries(entriesToDelete, array_in):
    # Set entires to NAN and drop them
    for i in entriesToDelete:
        array_in[0][i] = np.NAN
        array_in[1][i] = np.NAN
    array_x = array_in[0][~np.isnan(array_in[0])]
    array_y = array_in[1][~np.isnan(array_in[1])]

    # Construct new array
    array_out = np.zeros((2, len(array_x)))
    array_out[0] = array_x
    array_out[1] = array_y

    return array_out


def extract_curve_from_image(image_file, r, g, b, x_range, y_range):
    """
    Take an image and using the supplied r, g, b values get a coordinate list:
    r:          [low, high]
    g:          [low, high]
    b:          [low, high]
    x_range:    [low, high]
    y_range:    [low, high]
    """
    # Load image and convert it to a 3-colour array
    with Image.open(image_file) as raw_image:
        raw_image.show()
        raw_image_array = np.asarray(raw_image)

    x_dim = len(raw_image_array)
    y_dim = len(raw_image_array[0])
    processed_image_array = np.zeros((x_dim, y_dim, 3), dtype=np.uint8)

    x_coordinates = []
    y_coordinates = []
    # Run filtering
    for i in range(0, x_dim):
        for j in range(0, y_dim):
            if (
                (raw_image_array[i][j][0] <= r[1])
                and (raw_image_array[i][j][0] >= r[0])
                and (raw_image_array[i][j][1] <= g[1])
                and (raw_image_array[i][j][1] >= g[0])
                and (raw_image_array[i][j][2] <= b[1])
                and (raw_image_array[i][j][2] >= b[0])
            ):
                # Copy over colour
                processed_image_array[i][j][0] = raw_image_array[i][j][0]
                processed_image_array[i][j][1] = raw_image_array[i][j][1]
                processed_image_array[i][j][2] = raw_image_array[i][j][2]
                # Evaluate coordinate
                x = x_range[0] + float(x_range[1] - x_range[0]) * float(j) / float(
                    y_dim
                )
                # this is mish mash trial and _
                y = y_range[1] - float(y_range[1] - y_range[0]) * float(i) / float(
                    x_dim
                )
                x_coordinates.append(x)
                y_coordinates.append(y)
            else:
                processed_image_array[i][j][0] = 255
                processed_image_array[i][j][1] = 255
                processed_image_array[i][j][2] = 255

    processed_image = Image.fromarray(processed_image_array, "RGB")
    processed_image.show()

    return np.column_stack((x_coordinates, y_coordinates)).transpose()


def print_image_array(image_file):
    raw_image = Image.open(image_file)
    fitting_array = np.asarray(raw_image)
    print(fitting_array)
    return fitting_array
