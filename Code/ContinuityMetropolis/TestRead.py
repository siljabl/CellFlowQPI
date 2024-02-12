import numpy as np
import sys
import PIL
from PIL import Image

toreadfrom = sys.argv[1] 

img = Image.open(toreadfrom)
#img.show()

img_arr = np.array(img)

#print(np.max(img_arr) + "\t" + np.min(img_arr))

maxval = 1.0 / float(np.max(img_arr))

for i in range (0, 1024): #1048576
	for j in range (0, 1024):
		printval = float(img_arr[i,j]) * maxval
		print(str(i) + "\t" + str(j) + "\t" + str(printval))

