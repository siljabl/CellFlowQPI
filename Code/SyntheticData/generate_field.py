import os
import sys
import numpy as np

sys.path.append('src')
from velocity_fields import *
import plot_functions as plot
import save_functions as save

# Folders
v_file  = "testED.txt"  # change this
field_dir = sys.argv[1]

# Creating data folder
if os.path.isdir(field_dir) == 0:
    os.mkdir(field_dir)
    os.mkdir(field_dir + "/plots/")


''' VELOCITY FIELD '''
# Meshgrid
size = [1024, 1024]
x   = np.linspace(-.5, .5, size[0]+1)
y   = np.linspace(-.5, .5, size[1]+1)
X,Y = np.meshgrid(x,y)

# Define velocity field
velocity = field_3(X, Y)
velocity = normalize(*velocity, 1)

save.velocity_field(velocity, field_dir)
plot.im_velocity(velocity, field_dir)