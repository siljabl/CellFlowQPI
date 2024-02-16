import os
import sys
import numpy as np

sys.path.append("src")
from functions import *
from velocity_fields import normalize
import integration
import plot_functions as plot
import save_functions as save

# Folders
im_file   = sys.argv[1]
field_dir = sys.argv[2]

# Creating data folder
tif_dir = field_dir + "tif/"
if os.path.isdir(tif_dir) == 0:
    os.mkdir(tif_dir)

# Experimental parameters
d_cell = 40             # approximate cell diameter in pixels
u_max  = d_cell / 4     # max displacement in pixels per frame


''' IMPORT IMAGE DATA AND VELOCITY FIELD '''
# Import initial conditions
# read tif directly?
# Padd to ensure initial conditions are well behaved
f = 1/4                 # Defining subset. change so remove to new exponent of 2
pad_width = 10
init_cond = np.load(im_file)
init_cond = np.pad(init_cond, mode='linear_ramp', end_values=0, pad_width=pad_width)

# Import velocity field and normalize to u_max
u, v = np.loadtxt(field_dir + "x_velocity_full.txt"), np.loadtxt(field_dir + "y_velocity_full.txt")
u, v = normalize(*[u,v], u_max)

# Ensure intensity data and velocity field have same size.
intensity, velocity = size(init_cond, np.array([u, v]))
idx, region = sub_region(f, intensity)

plot.velocity_field(velocity, intensity, region, field_dir)
save.velocity_field(velocity, tif_dir, idx)


''' GENERATE DATA '''
# set time parameters
dt      = 0.1  # in frames
t_max   = 2    # in frames
t_steps = int(t_max / dt)

intensity = integration.RK(intensity, velocity, t_steps=t_steps, dt=dt)
save.intensity(intensity, idx, t_max, t_steps, tif_dir, im_file)
plot.mass_conservation(intensity, idx, dt, field_dir)


# SAVE PARAMETERS
im_size  = np.shape(init_cond)
tif_size = (int(im_size[0] * (1 - 2*f)), int(im_size[1] * (1 - 2*f)))

config = {"im_file"         : im_file,
          "full resolution" : im_size,
          "tif resolution"  : tif_size,
          "pad_width"       : pad_width,
          "d_cell"          : d_cell, 
          "u_max"           : u_max, 
          "dt"              : dt}

save.config(config, field_dir)
