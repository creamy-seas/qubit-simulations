from twin import twin
import numpy as np
from multiprocessing import Process
import multiprocessing
import time


def individual_run(twinInstance, EC, EJ, alpha, assymetry, file_to_write):
    # 1 - create instance, load experimental data, perform simulation and crank error
    twinInstance.override_parameters(EC, EJ, alpha, assymetry)
    twinInstance.simulate()
    error = twinInstance.experimental_data_error()

    # 2 - string file to write to file
    string_to_write = str(twinInstance.EC) + "\t" + str(twinInstance.EJ) + "\t" + \
        str(twinInstance.alpha) + "\t" + str(twinInstance.assymetry) + \
        "\t" + str(error) + "\t" + str(twinInstance.states_per_island) + "\n"

    # 3 - write file
    with open(file_to_write, 'a') as file_to_write:
        file_to_write.write(string_to_write)


def prepare_file(EC_param, EJ_param, AL_param, AS_param):
    """
    __ Parameters __
    arrays of [min, max, no_points]

    __ Description __
    creates a file that program will read inputs from
    """

    # 1 - define arrays
    EC_array = np.linspace(EC_param[0], EC_param[1], EC_param[2])
    EJ_array = np.linspace(EJ_param[0], EJ_param[1], EJ_param[2])
    AL_array = np.linspace(AL_param[0], AL_param[1], AL_param[2])
    AS_array = np.linspace(AS_param[0], AS_param[1], AS_param[2])

    # 2 - open file and fill it up
    fout = open("temp/parameters_file.txt", 'w')
    for EC in EC_array:
        for EJ in EJ_array:
            for AL in AL_array:
                for AS in AS_array:
                    string_to_write = "%.5f\t%.5f\t%.5f\t%.5f\n" % (
                        EC, EJ, AL, AS)
                    fout.write(string_to_write)

    fout.close()


def read_and_delete(file_name, no_lines_to_read):
    """
    __ Parameters __
    file_name to read from
    no_lines_to_read is the number of lines to read from begging of the file and cut out

    __ Description __
    reads in certain number of lines and deletes them from the file

    !!!!!!!!!! File must have format !!!!!!!!!!!!!!!!!!!!
    EC \tab EJ \tab alpha \tab assymetry \n
    """

    # 1 - extract the first "no_lines_to_read" as individual strings and the rest as bulk
    with open(file_name, "r+") as fin:
        split_data = fin.read().split("\n", no_lines_to_read)

    # 2 - save the split values
    stored_values = np.zeros((no_lines_to_read, 4))
    for i in range(0, no_lines_to_read):
        stored_values[i] = split_data[i].split("\t")

    # 3 - write new file, with the leftover data
    with open(file_name, "w") as fout:
        fout.write(split_data[no_lines_to_read])

    # fin.seek(0)
    # fin.write(split_data[no_lines_to_read])
    # fin.truncate()
    return stored_values


if (__name__ == "__main__"):
    # 1 - parameters
    file_to_write = "temp/temp_out.txt"
    number_of_threads = 16
    EC_param = [10, 10, 1]
    EJ_param = [73.5, 73.5, 1]
    AL_param = [1.0, 1.0, 1]
    AS_param = [1.015, 1.016, 2]

    # EC_param = [2, 42, 21]
    # EJ_param = [1, 101, 41]
    # AL_param = [1.0, 1.2, 21]
    # AS_param = [1, 1.1, 21]

    # 2 - prepare the file      <------ do not run this every time
    prepare_file(EC_param, EJ_param, AL_param, AS_param)

    # 3 - create class instances
    twinInstance = []
    for i in range(number_of_threads):
        twinInstance.append(twin(1, 1, 5, 300, False, False))
        twinInstance[i].experimental_data_load(twinInstance[i].ax, True)

    try:
        while True:
            start = time.time()
            # 4 - keep reading parameters from the data file
            parameter_array = read_and_delete(
                "temp/parameters_file.txt", number_of_threads)
            p = []

            # 5 - launch parrallel simulations using the parameter array
            for i in range(0, number_of_threads):
                p.append(Process(target=individual_run,
                                 args=(twinInstance[i],
                                       parameter_array[i][0],
                                       parameter_array[i][1],
                                       parameter_array[i][2],
                                       parameter_array[i][3],
                                       file_to_write)))
                p[i].start()

            # 6 - collect parrallel arrays together and then repeat loop
            for i in range(0, number_of_threads):
                p[i].join()

            end = time.time()
            print(end - start)

    except ValueError:
        # 7 - once the file is down to it's last parameters, launch the final set of threads

        number_of_lines_left = sum(
            1 for line in open("temp/parameters_file.txt"))
        parameter_array = read_and_delete(
            "temp/parameters_file.txt", number_of_lines_left)
        p = []

        # 8 - launch parrallel simulations using the parameter array
        for i in range(0, number_of_lines_left):
            p.append(Process(target=individual_run,
                             args=(twinInstance[i],
                                   parameter_array[i][0],
                                   parameter_array[i][1],
                                   parameter_array[i][2],
                                   parameter_array[i][3],
                                   file_to_write)))
            p[i].start()

            # 9 - collect parrallel arrays together and then repeat loop
        for i in range(0, number_of_lines_left):
            p[i].join()
