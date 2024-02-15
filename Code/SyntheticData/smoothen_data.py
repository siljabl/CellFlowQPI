import numpy as np
import scipy as sc
import matplotlib.pyplot as plt

# import rescaled and registered data
in_folder  = "../../Data/RescaledRegistered/"
out_folder = "../../Data/InitialConditions/" #RescaledRegisteredSmoothed/"
experiment = "Well1-1_resc_reg"

xlim = [500, 900]
ylim = [500, 900]

for i in range(900, 901):
    # define filenames
    in_file = in_folder + experiment + str(i) + ".tif"
    out_file = out_folder + experiment + str(i)

    data = plt.imread(in_file)
    smoothed_data = sc.ndimage.gaussian_filter(data, sigma=3)
    #np.save(out_file, smoothed_data[xlim[0]:xlim[1], ylim[0]:ylim[1]])
    np.save(out_file, smoothed_data)


