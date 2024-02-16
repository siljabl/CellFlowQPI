import sys
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt

from PIL import Image

# import rescaled and registered data
input   = sys.argv[1]
out_dir = sys.argv[2]
format  = sys.argv[3]
sigma = 3

# smoothen
im = plt.imread(input)
smoothed_data = sc.ndimage.gaussian_filter(im, sigma=sigma)

# save
filename = input.split('/')[-1]
if format == "tif":
    im_tif  = Image.fromarray(im.astype(np.uint8))
    im_tif.save(out_dir + filename)

if format == "npy":
    filename = filename.split('.')[0]
    np.save(out_dir + filename, smoothed_data)
