import os
import sys
import numpy as np

sys.path.append('src')
from functions import *
from velocity_fields import *
import plot_functions as plot
import save_functions as save

# Folders
file    = "Well1-1_resc_reg2"
in_folder  = "../Data/InitialConditions/"
out_folder = "../Data/Synthetic/Field4/"

# Creating data folder
if os.path.isdir(out_folder) == 0:
    os.mkdir(out_folder)
    os.mkdir(out_folder + "tif")


# Experimental parameters
d_cell = 40             # approximate cell diameter in pixels
u_max  = d_cell / 4     # max displacement in pixels per frame


# Import initial conditions
init_cond = np.load(in_folder + file + ".npy")
f = 0.1                 # Defining subset

# Padd to ensure initial conditions are well behaved
pad_width = 0
init_cond = np.pad(init_cond, mode='linear_ramp', end_values=0, pad_width=pad_width)



''' VELOCITY FIELD '''
# Meshgrid
size = np.shape(init_cond)
x   = size[0] * np.linspace(-.5, .5, size[0]) / d_cell
y   = size[1] * np.linspace(-.5, .5, size[1]) / d_cell
X,Y = np.meshgrid(x,y)

# Define velocity field
velocity = field_4(X,Y)
velocity = normalize(*velocity, u_max)

idx, region = sub_region(f, x, y)
plot.velocity_field(velocity, X, Y, init_cond, region, out_folder)
save.velocity_field(velocity, idx, out_folder)



''' GENERATE DATA '''
# set time parameters
dt      = 0.1  # in frames
t_max   = 2    # in frames
t_steps = int(t_max / dt)

intensity = RK_integration(init_cond, velocity, t_steps=t_steps, dt=dt)
plot.intensity(intensity, X, t_max, t_steps, region, out_folder, file)
save.intensity(intensity, idx, t_max, t_steps, out_folder, file)
plot.mass_conservation(intensity, idx, dt, out_folder)



# SAVE PARAMETERS
config = {"initial condition": file,
          "d_cell": d_cell, 
          "u_max": u_max, 
          "dt": dt,
          "pad_width": pad_width}

save.config(config, out_folder)
