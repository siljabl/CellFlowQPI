import sys
import scipy as sc
import numpy as np

import matplotlib.pyplot as plt

sys.path.append("src")
from functions import smoothen
from PIL import Image

# import rescaled and registered data
input   = sys.argv[1]
out_dir = "../../Data/InitialConditions" #sys.argv[2]
sigma = 2

# smoothen
im = plt.imread(input)
smoothed_data = smoothen(im)

# save
filename = input.split('/')[-1]
im_tif  = Image.fromarray(smoothed_data.astype(np.uint8))
im_tif.save(out_dir + filename)

# endre vekt mellom u og v
# divergensfelt
