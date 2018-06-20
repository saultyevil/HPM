import csv
import numpy as np
from itertools import islice
from collections import Counter
from matplotlib import pyplot as plt


def find_dir(work_dir=True):
    """
    """
    if work_dir is False:
        # Figure out which OS is in use
        exit_flag = False
        while exit_flag is False:
            exit_flag = True
            print("Operating system in use:\nWindows(W)\nLinux(L)\nOS X(O)")
            OS = input("> ")
            if OS == "W":
                main_data_dir = "F:\\Google Drive\\RESEARCH PROJECT\\data\\"
            elif OS == "L":
                main_data_dir = "/home/saultyevil/Downloads/PH600/data/"
            elif OS == "O":
                main_data_dir = "/Users/saultyevil/Google Drive/PH600/data/"
            else:
                print("Invalid input.")
                exit_flag = False
        # Figure out if the default directory is correct
        exit_flag = False
        while exit_flag is False:
            exit_flag = True
            print("\nMain data directory: {}\nIs this correct?\n[Y/N]".format(
                main_data_dir))
            correct_dir = input("> ")
            if correct_dir == "N":
                main_data_dir = input("Please input data directory path.\n> ")
            elif correct_dir == "Y":
                break
            else:
                print("Invalid input.")
                exit_flag = False
    else:
        main_data_dir = "./"
    uwish2_name = input("Enter filename of the UWISH2 data:\n> ")
    ukidss_name = input("Enter the filename of the UKIDSS GPS data:\n> ")
    uwish2_path = main_data_dir + uwish2_name
    ukidss_path = main_data_dir + ukidss_name
    return uwish2_path, ukidss_path


def read_data(uwish2_filename, ukidss_filename, start_pos_uwish2,
              start_pos_ukidss):
    """
    """
    uwish2_file = open(uwish2_filename, "r")
    uwish2_CSV = csv.reader(uwish2_file)
    ukidss_file = open(ukidss_filename, "r")
    ukidss_CSV = csv.reader(ukidss_file)
    n_sources_uwish2 = sum(1 for line in uwish2_file) - (start_pos_uwish2 + 1)
    n_sources_ukidss = sum(1 for line in ukidss_file) - (start_pos_ukidss + 1)
    uwish2_data = np.zeros((n_sources_uwish2, 4))
    ukidss_data = np.zeros((n_sources_ukidss, 10))
    uwish2_file.seek(0)
    ukidss_file.seek(0)
    # Read uwish2 data from file
    for source in islice(uwish2_CSV, start_pos_uwish2, None):
        for sourceID, source in enumerate(uwish2_CSV):
            uwish2_data[sourceID, 0] = source[3]    # RA
            uwish2_data[sourceID, 1] = source[4]    # DEC
            uwish2_data[sourceID, 2] = source[45]   # distance
            uwish2_data[sourceID, 3] = source[28]   # H2 mag
    # Read ukidss data from file
    for source in islice(ukidss_CSV, start_pos_ukidss, None):
        for sourceID, source in enumerate(ukidss_CSV):
            ukidss_data[sourceID, 0] = source[2]    # RA
            ukidss_data[sourceID, 1] = source[3]    # DEC
            ukidss_data[sourceID, 2] = source[18]   # distance
            ukidss_data[sourceID, 3] = source[17]   # date
            ukidss_data[sourceID, 4] = source[6]    # J mag
            ukidss_data[sourceID, 5] = source[7]    # J mag error
            ukidss_data[sourceID, 6] = source[8]    # H mag
            ukidss_data[sourceID, 7] = source[9]    # H mag error
            ukidss_data[sourceID, 8] = source[10]   # K mag
            ukidss_data[sourceID, 9] = source[11]   # K mag error
    uwish2_file.close()
    ukidss_file.close()
    # Sort the data via distance using ugly, hacky method :^)
    uwish2_data = uwish2_data[uwish2_data[:, 2].argsort()]
    ukidss_data = ukidss_data[ukidss_data[:, 2].argsort()]
    print("\nUWISH2 RA[0]\tUKIDSS RA[0]\n{}\t{}".format(uwish2_data[0, 0],
          ukidss_data[0, 0]))
    print("UWISH2 DEC[0]\tUKIDSS DEC[0]\n{}\t{}".format(uwish2_data[0, 1],
          ukidss_data[0, 1]))
    return uwish2_data, ukidss_data


def source_match(uwish2, ukidss):
    """
    """
    print("\nSources will be matched with objects within a 1 arcsecond region "
          "between both images and with similar H2/K-band magnitudes.")
    uwish2_match = np.zeros((1, 4))
    ukidss_match = np.zeros((1, 10))
    n_sources_uwish2 = uwish2.shape[0]
    n_sources_ukidss = ukidss.shape[0]
    deg_limit = 1/3600
    k_match = 0.5
    cand_coords = input("Coordinates of potential HPM object from UWISH2, "
                        "separated by a comma:\n> ")
    cand_coords = np.array(cand_coords.replace(",", "").split(), dtype=float)
    cand_H2 = uwish2[uwish2[:, 0] == cand_coords[0], 3]
    # Check against the RA, DEC and K/H2 limits to see if the target HPM star
    # in UWISH2 can be matched with a star in UKIDSS. If the HPM candidnate
    # isn't matched, the program will quit
    # The candidate isn't stored just yet...
    cand_match = False
    for i in range(n_sources_ukidss):
        RA_lims = ukidss[i, 0] + np.array([-deg_limit, deg_limit])
        DEC_lims = ukidss[i, 1] + np.array([-deg_limit, deg_limit])
        k_lims = ukidss[i, 8] + np.array([-k_match, k_match])
        if RA_lims[0] < cand_coords[0] < RA_lims[1]:
            if DEC_lims[0] < cand_coords[1] < DEC_lims[1]:
                if k_lims[0] < cand_H2 < k_lims[1]:
                    # uwish2_match[0, :] = uwish2[i, :]
                    # ukidss_match[0, :] = ukidss[i, :]
                    cand_match = True
                    break
    if cand_match is True:
        print("\n---------------------------------------------")
        print("PASS: Match found for potential HPM object.")
        # Now try to find matching objects for the rest of the objects - using
        # Python lists as this was a more robust way as the number of objects
        # which will be matched it unknown at the start and numpy vstack was a
        # bit odd
        for i in range(n_sources_ukidss):
            RA_lims = ukidss[i, 0] + np.array([-1, 1]) * deg_limit
            DEC_lims = ukidss[i, 1] + np.array([-1, 1]) * deg_limit
            k_lims = ukidss[i, 8] + np.array([-1, 1]) * k_match
            for j in range(n_sources_uwish2):
                if RA_lims[0] < uwish2[j, 0] < RA_lims[1]:
                    if DEC_lims[0] < uwish2[j, 1] < DEC_lims[1]:
                        if k_lims[0] < uwish2[j, 3] < k_lims[1]:
                            ukidss_match = np.vstack([ukidss_match,
                                                      ukidss[i, :]])
                            uwish2_match = np.vstack([uwish2_match,
                                                      uwish2[j, :]])
        uwish2_match = np.delete(uwish2_match, 0, 0)
        ukidss_match = np.delete(ukidss_match, 0, 0)
        # Now we need to delete objects which have been matched with multiple
        # objects since we can't be 100% sure without manually checking if the
        # object has been matched correctly
        # Find the number of RA matches
        uwish2_no_dup = np.zeros((1, 4))
        ukidss_no_dup = np.zeros((1, 10))
        matches_dict = Counter(ukidss_match[:, 0])
        for i in range(ukidss_match.shape[0]):
            if matches_dict[ukidss_match[i, 0]] == 1:
                uwish2_no_dup = np.vstack([uwish2_no_dup, uwish2_match[i, :]])
                ukidss_no_dup = np.vstack([ukidss_no_dup, ukidss_match[i, :]])
        uwish2_match = np.delete(uwish2_no_dup, 0, 0)
        ukidss_match = np.delete(ukidss_no_dup, 0, 0)
        print("      {} objects have been matched."
              .format(uwish2_match.shape[0]))
        print("---------------------------------------------")
        return uwish2_match, ukidss_match, cand_coords
    else:
        print("\nFAIL: A match was not found for the potential HPM object.")
        restart()
    return


def calc_PM(uwish2, ukidss, colour_correct=False):
    """
    """
    uwish2_epoch = float(input("Enter the date of the UWISH2 data:\n>  "))
    ukidss_epoch = ukidss[0, 3]
    epochs = np.array([uwish2_epoch, ukidss_epoch])
    d_time = uwish2_epoch - ukidss_epoch
    pm_RA_uncorrected = ((uwish2[:, 0] - ukidss[:, 0]) * 3600) / d_time
    pm_DEC_uncorrected = ((uwish2[:, 1] - ukidss[:, 1]) * 3600) / d_time
    # Calculate the median of the data to subtract from all of the values to
    # correct for offset between the two survey images
    pm_RA_median = np.median(pm_RA_uncorrected)
    pm_DEC_median = np.median(pm_DEC_uncorrected)
    # Multiply the proper motion by 1000 for units of milliarcsecond
    pm_DEC = 1000 * (pm_DEC_uncorrected - pm_DEC_median)
    pm_RA = 1000 * np.cos(np.radians(pm_DEC)) * (pm_RA_uncorrected -
                                                 pm_RA_median)
    pm_magnitude = np.sqrt(pm_RA ** 2 + pm_DEC ** 2)
    pm = np.array([pm_RA, pm_DEC, pm_magnitude]).T
    # Remove objects due to their colour -- this is due to reddencing
    # as well as objects having crap values/a false colour as they are
    # likely noise/false detections or have been influenced by noise
    # This check will be done by considering the colour difference between
    # the H and K band magnitudes. This will only be applied to objects which
    # have a proper motion magnitude of less than 200 mas/yr
    if colour_correct is True:
        uwish2_fix = np.zeros((1, 4))
        ukidss_fix = np.zeros((1, 10))
        pm_fix = np.zeros((1, 3))
        pm_cutoff = 200
        upper_inter = 0.7
        x_coeff = 1.55
        for i in range(ukidss.shape[0]):
            if pm_magnitude[i] < pm_cutoff:
                H_K_colour = ukidss[i, 6] - ukidss[i, 8]
                if 0 < H_K_colour < 0.7:
                    reddening_inter = x_coeff * H_K_colour + upper_inter
                    reddening_no_inter = x_coeff * H_K_colour
                    if reddening_no_inter < H_K_colour < reddening_inter:
                        uwish2_fix = np.vstack([uwish2_fix, uwish2[i, :]])
                        ukidss_fix = np.vstack([ukidss_fix, ukidss[i, :]])
                        pm_fix = np.vsatck([pm_fix, pm[i, :]])
        uwish2 = np.delete(uwish2_fix, 0, 0)
        ukidss = np.delete(ukidss_fix, 0, 0)
        pm = np.delete(pm_fix, 0, 0)
    return uwish2, ukidss, pm, epochs


def output_PM(uwish2, ukidss, pm, cand_coords, epochs):
    """
    """
    # Write some stuff out about the potential HPM candidate
    cand_index = np.where(uwish2[:, 0] == cand_coords[0])
    print("\n----------------------------------------------------------")
    print("HPM Object:\nPM RA: {:3.3f} mas/yr\nPM_DEC: {:3.3f} mas/yr\n"
          "PM_MAG: {:3.3f} mas/yr".format(float(pm[cand_index, 0]),
                                          float(pm[cand_index, 1]),
                                          float(pm[cand_index, 2])))
    print("----------------------------------------------------------")
    # Graph the final output
    pm_mag_std_dev = np.std(pm[:, 2])
    sigma_lev1 = plt.Circle((0, 0), 1 * pm_mag_std_dev, color="r", fill=False)
    sigma_lev2 = plt.Circle((0, 0), 2 * pm_mag_std_dev, color="r", fill=False)
    sigma_lev3 = plt.Circle((0, 0), 3 * pm_mag_std_dev, color="r", fill=False)
    sigma_lev4 = plt.Circle((0, 0), 4 * pm_mag_std_dev, color="r", fill=False)
    sigma_lev5 = plt.Circle((0, 0), 5 * pm_mag_std_dev, color="r", fill=False)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(pm[np.where(uwish2[:, 0] != cand_coords[0]), 0],
            pm[np.where(uwish2[:, 0] != cand_coords[0]), 1], "kx")
    ax.plot(pm[cand_index, 0], pm[cand_index, 1], "rD")
    ax.set_xlabel(r"$\mu_{RA} \cos(\delta)$ (mas/yr)")
    ax.set_ylabel(r"$\mu_{\delta}$ (mas/yr)")
    ax.add_patch(sigma_lev1)
    ax.add_patch(sigma_lev2)
    ax.add_patch(sigma_lev3)
    ax.add_patch(sigma_lev4)
    ax.add_patch(sigma_lev5)
    plt.axis("equal")
    plt.grid(True)
    plt.savefig("({})({}).pdf".format(cand_coords[0], cand_coords[1]))
    plt.show()
    # Output all of the data to file
    output_file = open("({})({}).txt".format(cand_coords[0], cand_coords[1]),
                       "w")
    output_file.write("UWISH2\t{}\tUKIDSS\t{}\n".format(epochs[0], epochs[1]))
    output_file.write("CAND_RA\t{}\tCAND_DEC\t{}\n".format(cand_coords[0],
                      cand_coords[1]))
    output_file.write("UWISH2_RA\tUWISH2_DEC\tUWISH2_DIST\tUWISH2_H2"
                      "\tUKIDSS_RA\tUKIDSS_DEC\tUKIDSS_DIST\tUKIDSS_J"
                      "\tUKIDSS_J_ERR\tUKIDSS_H\tUKIDSS_H_ERR\tUKIDSS_K"
                      "\tUKIDSS_K_ERR\tPM_RA\tPM_DEC\tPM_MAG\n")
    for i in range(uwish2.shape[0]):
        output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}"
                          "\t{}\t{}\t{}\t{}\n"
                          .format(uwish2[i, 0], uwish2[i, 1],
                                  uwish2[i, 2], uwish2[i, 3], ukidss[i, 0],
                                  ukidss[i, 1], ukidss[i, 2], ukidss[i, 3],
                                  ukidss[i, 4], ukidss[i, 5], ukidss[i, 6],
                                  ukidss[i, 7], ukidss[i, 8], ukidss[i, 9],
                                  pm[i, 0], pm[i, 1], pm[i, 2]))
    output_file.close()
    return


def restart(work_dir=True):
    """
    """
    exit_flag = False
    while exit_flag is False:
        restart = input("Restart program? [Y/N]\n> ")
        print("")
        if restart == "Y":
            main(work_dir)
        elif restart == "N":
            exit_flag = True
        else:
            print("Unknown input. Try 'Y or 'N'.")
    return


def main(work_dir=True):
    """
    """
    print("----------------------------------------------------------")
    print("| Calculate  the proper motion of a candidate HPM object |")
    print("----------------------------------------------------------")
    uwish2_path, ukidss_path = find_dir()
    # Read in data from file for both surveys
    uwish2_skip = 14
    ukidss_skip = 18
    print("\nReading in data from file. Skipping {} lines for UWISH2 and"
          " {} lines for UKIDSS data files.".format(uwish2_skip, ukidss_skip))
    uwish2_data, ukidss_data = read_data(uwish2_path, ukidss_path, uwish2_skip,
                                         ukidss_skip)
    # Match the sources between the two data sets
    uwish2_data, ukidss_data, cand_coords = source_match(uwish2_data,
                                                         ukidss_data)
    # Calculate the proper motions
    uwish2_data, ukidss_data, pm, epochs = calc_PM(uwish2_data, ukidss_data)
    # Plot the data and output to file
    output_PM(uwish2_data, ukidss_data, pm, cand_coords, epochs)
    restart(work_dir)
    return


# Call the main function to run the script :-)
if __name__ == "__main__":
    work_dir = True
    main(work_dir)
