import os
import sys
import numpy as np

sys.path.append('src')
from functions import *
from velocity_fields import *
import plot_functions as plot
import save_functions as save

# Folders
v_file  = "testED.txt"
in_folder  = "../../Data/InitialConditions/"
out_folder = "../../Data/Synthetic/Test/"

# Creating data folder
if os.path.isdir(out_folder) == 0:
    os.mkdir(out_folder)
    os.mkdir(out_folder + "plots")

''' VELOCITY FIELD '''
# Meshgrid
size = [1024, 1024]
x   = size[0] * np.linspace(-.5, .5, size[0])
y   = size[1] * np.linspace(-.5, .5, size[1])
X,Y = np.meshgrid(x,y)

# Define velocity field
velocity = 0.2*field_3(X,Y) + field_exp(x,y, 10, v_file)
velocity = normalize(*velocity, 1)

idx, region = sub_region(f, x, y)
save.velocity_field(velocity, idx, out_folder)