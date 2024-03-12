import sys
import numpy as np
from PIL import Image

input = sys.argv[1]
filename = sys.argv[2]

img = Image.open(input)
img_arr = np.array(img)

np.savetxt(filename + ".txt", img_arr)

